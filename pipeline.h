#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <cstdint>
#include <cstdlib>
#include <cstddef>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <new>
#include <optional>
#include <queue>
#include <vector>
#include "nr_common.h"

template <typename T>
class BoundedQueue {
public:
    explicit BoundedQueue(std::size_t max_size = 0)
        : max_size_(max_size), closed_(false) {}

    bool push(T item) {
        std::unique_lock<std::mutex> lock(mtx_);
        not_full_.wait(lock, [&] {
            return closed_ || max_size_ == 0 || queue_.size() < max_size_;
        });
        if (closed_) return false;
        queue_.push(std::move(item));
        lock.unlock();
        not_empty_.notify_one();
        return true;
    }

    // push without waiting: false (item untouched) if the queue is full or closed
    bool try_push(T&& item) {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            if (closed_ || (max_size_ != 0 && queue_.size() >= max_size_)) return false;
            queue_.push(std::move(item));
        }
        not_empty_.notify_one();
        return true;
    }

    std::optional<T> pop() {
        std::unique_lock<std::mutex> lock(mtx_);
        not_empty_.wait(lock, [&] { return !queue_.empty() || closed_; });
        if (queue_.empty()) return std::optional<T>();
        T item = std::move(queue_.front());
        queue_.pop();
        lock.unlock();
        not_full_.notify_one();
        return std::optional<T>(std::move(item));
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            closed_ = true;
        }
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    mutable std::mutex      mtx_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    std::queue<T>           queue_;
    std::size_t             max_size_;
    bool                    closed_;
};

// Raw aligned buffer, without initialization (no zero-fill).
// capacity_ = how many elements fit; size_ = how many are actually valid.
class AlignedBuffer {
public:
    AlignedBuffer(std::size_t capacity, std::size_t align_bytes)
        : data_(nullptr), capacity_(capacity), size_(0)
    {
        // posix_memalign requires an alignment that is a multiple of sizeof(void*)
        // and a power of two; 64 is fine.
        void* p = nullptr;
        std::size_t bytes = capacity_ * sizeof(cfloat);
        // Avoid size 0 (implementation-defined behavior)
        if (bytes == 0) bytes = align_bytes;
        int rc = ::posix_memalign(&p, align_bytes, bytes);
        if (rc != 0 || p == nullptr) {
            throw std::bad_alloc();
        }
        data_ = static_cast<cfloat*>(p);
        // NO initialization: cfloat is trivially constructible,
        // the producer fills the first size_ elements before use.
    }

    ~AlignedBuffer() { std::free(data_); }

    AlignedBuffer(const AlignedBuffer&)            = delete;
    AlignedBuffer& operator=(const AlignedBuffer&) = delete;

    cfloat*       data()       noexcept { return data_; }
    const cfloat* data() const noexcept { return data_; }

    std::size_t size()     const noexcept { return size_; }
    std::size_t capacity() const noexcept { return capacity_; }

    void set_size(std::size_t n) noexcept { size_ = n; }

private:
    cfloat*     data_;
    std::size_t capacity_;
    std::size_t size_;
};

class BufferPool {
public:
    typedef AlignedBuffer Buffer;

    struct Recycler {
        BufferPool* pool;
        void operator()(Buffer* b) const noexcept {
            if (!b) return;
            if (pool) pool->release(b);
            else      delete b;
        }
    };

    typedef std::unique_ptr<Buffer, Recycler> Handle;

    BufferPool(std::size_t capacity_elems,
               std::size_t align_bytes,
               std::size_t initial_count,
               std::size_t max_count)
        : capacity_(capacity_elems),
          align_(align_bytes),
          max_count_(max_count),
          allocated_(0),
          closed_(false)
    {
        for (std::size_t i = 0; i < initial_count && i < max_count; ++i) {
            free_.push(new Buffer(capacity_, align_));
            ++allocated_;
        }
    }

    ~BufferPool() {
        close();
        std::lock_guard<std::mutex> lock(mtx_);
        while (!free_.empty()) {
            delete free_.front();
            free_.pop();
        }
    }

    Handle acquire() {
        std::unique_lock<std::mutex> lock(mtx_);
        not_empty_.wait(lock, [&] {
            return closed_ || !free_.empty() || allocated_ < max_count_;
        });
        if (closed_) {
            Recycler rec;
            rec.pool = this;
            return Handle(nullptr, rec);
        }

        Buffer* b;
        if (!free_.empty()) {
            b = free_.front();
            free_.pop();
        } else {
            b = new Buffer(capacity_, align_);
            ++allocated_;
        }
        Recycler rec;
        rec.pool = this;
        return Handle(b, rec);
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            closed_ = true;
        }
        not_empty_.notify_all();
    }

private:
    void release(Buffer* b) {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            if (closed_) { delete b; --allocated_; return; }
            free_.push(b);
        }
        not_empty_.notify_one();
    }

    mutable std::mutex      mtx_;
    std::condition_variable not_empty_;
    std::queue<Buffer*>     free_;
    std::size_t             capacity_;
    std::size_t             align_;
    std::size_t             max_count_;
    std::size_t             allocated_;
    bool                    closed_;
};

// CORESET#0 in the grid of the frames (from the MIB), searched by the DCI workers
struct CORESET0Grid {
    bool valid = false;
    int k_start = 0;          // first subcarrier in the grid (N_RB*12 carriers)
    int N_RB = 0;
    int duration = 0;
    uint32_t PCI = 0;         // n_shift of the REG bundle interleaver, DM-RS and data scrambling ID
    // Type0-PDCCH (SIB1) monitoring occasions: slots n0 and n0+1 of the frames with SFN mod 2 == frame_parity,
    // from first_symbol (the other common search spaces on CORESET#0 are searched from symbol 0 of every slot)
    int n0 = 0;
    int frame_parity = 0;
    int first_symbol = 0;
};

struct Job {
    uint16_t LFN;
    long long rx_ofs;
    BufferPool::Handle Y;
    std::vector<cfloat> Y_CP_corr;
    int mu;
    std::vector<bool> is_Y_valid;
    float noise_est;
    CORESET0Grid cs0;

    Job() : LFN(0), rx_ofs(0), Y(nullptr, BufferPool::Recycler()), Y_CP_corr(0), mu(0), is_Y_valid(0), noise_est(0) {}
    Job(uint16_t LFN, long long rx_ofs, BufferPool::Handle hdl, std::vector<cfloat> &Y_CP_corr, int mu, std::vector<bool> &is_Y_valid, float noise_est)
        : LFN(LFN), rx_ofs(rx_ofs), Y(std::move(hdl)), Y_CP_corr(Y_CP_corr), mu(mu), is_Y_valid(is_Y_valid), noise_est(noise_est) {}

    // Move-only (because of unique_ptr)
    Job(Job&& other) noexcept : LFN(other.LFN), rx_ofs(other.rx_ofs), Y(std::move(other.Y)), Y_CP_corr(std::move(other.Y_CP_corr)), mu(other.mu), is_Y_valid(std::move(other.is_Y_valid)), noise_est(other.noise_est), cs0(other.cs0) {}
    Job& operator=(Job&& other) noexcept {
        LFN = other.LFN;
        rx_ofs = other.rx_ofs;
        Y  = std::move(other.Y);
        Y_CP_corr = std::move(other.Y_CP_corr);
        mu = other.mu;
        is_Y_valid = std::move(other.is_Y_valid);
        noise_est = other.noise_est;
        cs0 = other.cs0;
        return *this;
    }
    Job(const Job&)            = delete;
    Job& operator=(const Job&) = delete;
};

// what a Result carries
enum ResultKind : uint8_t {
    RES_DCI,        // decoded DCI
    RES_TRS,        // NZP-CSI-RS density 3 (TRS): slot, symbol, k_ini=k0, DCI_NID=NID
    RES_CSI,        // NZP-CSI-RS density 1, no CDM: same fields as RES_TRS
    RES_CSI_FD2,    // NZP-CSI-RS density 1, CDM2 (fd2): same fields as RES_TRS
    RES_FRAME_DONE, // a worker finished the frame starting at rx_ofs
    RES_TICK        // RX thread: one second of stream time is over, rx.* holds its stats
};

// receiver state over one second of stream time (RES_TICK)
struct RxStatus {
    long long sec;            // second index: rx_ofs/samp_rate in [sec, sec+1)
    bool      partial;        // last, incomplete second of the run
    uint8_t   state;          // 0 AGC, 1 ACQUIRE, 2 TRACK (at the end of the second)
    int16_t   PCI;            // -1 if not tracking
    double    f_SSB;          // SSB center frequency [Hz]
    uint16_t  ssb_found;      // frames with the expected SS block
    uint16_t  ssb_missed;     // frames without it
    float     snr_dB;         // average SS SNR over the frames with the SS block
    float     cfo_Hz;
    float     sfo_ppm;
    float     gain_dB;        // NaN with --file
    uint32_t  frames_shipped; // frames sent to the FFT worker during the second
    uint16_t  frames_dropped; // frames not sent because fft_q was full (live only: the workers can't keep up)
    uint16_t  overflows;
    uint16_t  acquisitions;   // ACQUIRE -> TRACK transitions
    uint16_t  track_losses;   // TRACK -> AGC transitions (ACQUIRE with --file)
    uint16_t  acq_failures;   // ACQUIRE -> AGC transitions (--ACQ_max_fails searches without a cell)
    uint16_t  pbch_ok;        // PBCH decoded (CRC ok)
    uint16_t  pbch_fail;      // PBCH decoding attempts that failed
    int16_t   SFN;            // SFN of the frame being received at the end of the second, -1 until the MIB is decoded
};

struct Result {
    ResultKind kind;
    long long rx_ofs;
    uint16_t LFN;
    uint8_t slot;
    uint16_t k_ini;
    uint16_t k_fin;
    uint8_t AL;
    uint8_t duration; // CORESET symbols
    uint16_t DCI_NID;
    uint16_t RNTI;
    uint8_t DCI_len;
    uint8_t DCI[140];
    uint8_t symbol; // CSI-RS: symbol within the slot
    bool CORESET0;  // DCI: found in CORESET#0 (k_ini..k_fin spans its CCEs, which are interleaved)
    RxStatus rx;    // RES_TICK only
};

struct FFTJob {
    uint16_t		LFN;
    long long           rx_ofs;
    std::vector<cfloat> Y;          // raw copy: time-domain for unprocessed symbols,
                                    // freq-domain for tracking symbols (PSS/SSS)
    std::vector<bool>   is_Y_valid; // true for symbols already in freq domain
    std::vector<cfloat> Y_CP_corr;
    int                 mu;
    float               noise_est;
    float               Cfo_eF;
    CORESET0Grid        cs0;
};

#endif
