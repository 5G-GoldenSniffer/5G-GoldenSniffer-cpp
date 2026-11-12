#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <atomic>
#include <thread>
#include <algorithm>
#include "gf2solve.h"
#include "polar.h"
#include "NRscrambling.h"
#include "DCIfactory.h"

void DCIfactory::encode(std::vector<uint8_t> &b_tilde, const std::vector<uint8_t> &a, uint32_t rnti_crc,
                uint32_t rnti_scr, uint32_t N_ID, uint32_t E, uint32_t crc_init, uint32_t x1_init)
{
    uint32_t k;
    
    size_t A = a.size();
    // CRC attachment 38.212 7.3.2 p.193, the last 16 parity bits scrambled with the RNTI
    const uint32_t L = 24;
    size_t K = A+L;
    uint32_t p = crc24c(a.data(), A, crc_init) ^ rnti_crc;
    for(k = 0; k < A; k++)
        DCIEnc_c[k] = a[k];
    for(k = 0; k < L; k++)
        DCIEnc_c[A+k] = (p>>(23-k)) & 0x00000001;
    
    // Channel coding 38.212 7.3.3 p.194 (n_max = 9, I_IL = 1, n_PC = 0) and
    // rate matching 38.212 7.3.4 p.194 (I_BIL = 0)
    if(!polar || polar->K != K || polar->E != E)
        polar = std::make_shared<const PolarCode>(K, E, 9, true);
    polar->encode(DCIEnc_c, DCIEnc_f);
    
    // scrambling 38.211 7.3.2.3 p.99, 38.211 5.2.1 p.18
    //std::vector<uint8_t> DCIEnc_c_DCI(E);
    scrambling(DCIEnc_c_DCI, E, ((rnti_scr<<16) ^ N_ID) & 0x7fffffff);
    if(x1_init==0){
        //std::vector<uint8_t> DCIEnc_x1(E);
        scrambling_x1(DCIEnc_x1, E);
        for(k = 0; k < E; k++)
            DCIEnc_c_DCI[k] ^= DCIEnc_x1[k];
    }
    for(k = 0; k < E; k++)
        b_tilde[k] = DCIEnc_f[k] ^ DCIEnc_c_DCI[k];
}


// the decoding matrices of the hypotheses (log2AL, A) taken from next, in decreasing E (the largest first)
static void init_tables(DCItables &T, std::atomic<uint32_t> &next, const std::vector<std::pair<uint32_t,uint32_t>> &hyp)
{
    // matrix-based DCI encoding
    const uint32_t N_RNTI_CRC = 16;
    DCIfactory enc; // encoding buffers of this thread
    std::vector<uint8_t> a;
    std::vector<uint8_t> b_tilde(108*16);

    for(uint32_t i; (i = next.fetch_add(1)) < hyp.size(); ){
        const uint32_t log2AL = hyp[i].first, A = hyp[i].second;
        const uint32_t E = 108u<<log2AL;
        const uint32_t M = A+N_RNTI_CRC;

        mzd_t *G_A = mzd_init(E, M);
        mzp_t *P_A = mzp_init(E);
        mzd_t *G_B = mzd_init(E, M);
        mzp_t *P_B = mzp_init(E);
        mzv_t *b0 = mzv_init(E);

        a.assign(A, 0);
        uint32_t RNTI_CRC = 0x0000;
        uint32_t RNTI_SCR = 0x0000;

        // both cases use the rows from a[i] and RNTI_CRC
        for(uint32_t ic = 0; ic < A; ic++){
            a[ic] = 1;
            enc.encode(b_tilde, a, RNTI_CRC, RNTI_SCR, 0, E, 0, 0);
            for(uint32_t ir = 0; ir < E; ir++){
                mzd_write_bit(G_A,ir,ic,b_tilde[ir]);
                mzd_write_bit(G_B,ir,ic,b_tilde[ir]);
            }
            a[ic] = 0;
        }
        for(uint32_t ic = A; ic < M; ic++){
            RNTI_CRC = (0x0001<<(M-1 - ic));
            enc.encode(b_tilde, a, RNTI_CRC, RNTI_SCR, 0, E, 0, 0);
            for(uint32_t ir = 0; ir < E; ir++){
                mzd_write_bit(G_A,ir,ic,b_tilde[ir]);
                mzd_write_bit(G_B,ir,ic,b_tilde[ir]);
            }
        }
        RNTI_CRC = 0x0000;

        // RNTI_SCR is combined with RNTI_CRC in Case A
        for(uint32_t ic = A+1; ic < M; ic++){
            RNTI_SCR = (0x0001<<(M-1 - ic));
            enc.encode(b_tilde, a, RNTI_CRC, RNTI_SCR, 0, E, 0, 0);
            for(uint32_t ir = 0; ir < E; ir++)
                mzd_write_bit(G_A,ir,ic,b_tilde[ir] ^ mzd_read_bit(G_A,ir,ic));
        }
        RNTI_SCR = 0x0000;

        mzd_plu(G_A, P_A);
        mzd_plu(G_B, P_B);

        // the scrambling with N_ID (linear in N_ID, the same for every A) is removed in decode()
        enc.encode(b_tilde, a, RNTI_CRC, RNTI_SCR, 0, E, 0x00ffffff, 0x01);
        for(uint32_t ir = 0; ir < E; ir++)
            b0->d[ir] = b_tilde[ir];

        T.G_A[log2AL][A-A_min] = G_A;
        T.P_A[log2AL][A-A_min] = P_A;
        T.G_B[log2AL][A-A_min] = G_B;
        T.P_B[log2AL][A-A_min] = P_B;
        T.b0 [log2AL][A-A_min] = b0;
    }
}

void DCIfactory::init()
{
    // the hypotheses, the most expensive first, built in parallel
    std::vector<std::pair<uint32_t,uint32_t>> hyp;
    for(int log2AL = log2AL_max; log2AL >= (int)log2AL_min; log2AL--)
        for(uint32_t A = A_min; A <= A_max[log2AL]; A++)
            hyp.push_back({log2AL, A});
    auto T = std::make_shared<DCItables>();
    std::atomic<uint32_t> next(0);
    std::vector<std::thread> thr(std::max(1U, std::thread::hardware_concurrency()));
    for(auto &t : thr)
        t = std::thread(init_tables, std::ref(*T), std::ref(next), std::cref(hyp));
    for(auto &t : thr)
        t.join();
    tables = T;
}

bool DCIfactory::decode(uint32_t &RNTI_CRC, std::vector<uint8_t> &a, const std::vector<uint8_t> &b_tilde, uint32_t log2AL, uint32_t A, uint32_t N_ID)
{
    const uint32_t N_RNTI_CRC = 16;
    const uint32_t E = 108u<<log2AL;
    const DCItables &T = *tables;

    // b_NID: the x2 sequence of the scrambling with c_init = N_ID (38.211 7.3.2.3), computed once per N_ID
    uint8_t *b_NID = this->b_NID[log2AL];
    if(b_NID_key[log2AL] != N_ID+1){
        scrambling(b_NID, E, N_ID);
        scrambling_x1(DCIEnc_x1, E);
        for(uint32_t ir = 0; ir < E; ir++)
            b_NID[ir] ^= DCIEnc_x1[ir];
        b_NID_key[log2AL] = N_ID+1;
    }

    // G_{A,B} * [a;RNTI_CRC] = b_tilde - b_tilde0 - b_NID
    const uint8_t *b0 = T.b0[log2AL][A-A_min]->d;
    for(uint32_t ir = 0; ir < E; ir++)
        b[ir] = b_NID[ir]^b0[ir]^b_tilde[ir];

    // case B, then case A
    uint8_t x[DCIEnc_Kmax];
    for(int tr = 0; tr < 2; tr++){
        const bool ok = tr == 0 ? mzdv_lu_solve_perm(T.G_B[log2AL][A-A_min], T.P_B[log2AL][A-A_min], b, x)
                                : mzdv_lu_solve_perm(T.G_A[log2AL][A-A_min], T.P_A[log2AL][A-A_min], b, x);
        if(ok){
            for(uint32_t ir = 0; ir < A; ir++)
                a[ir] = x[ir];
            uint32_t tmp = 0x0000;
            for(uint32_t ir = A; ir < A+N_RNTI_CRC; ir++)
                tmp |= (x[ir]<<(A+N_RNTI_CRC-1-ir));
            RNTI_CRC = tmp;
            return true;
        }
    }
    return false;
}

bool DCIfactory::selftest(bool verbose)
{
    bool is_ok = true;
    
    // matrix-based DCI encoding and decoding
    uint32_t RNTI_CRC;
    uint32_t RNTI_SCR;
    uint32_t N_ID;
    uint32_t crc_init = 0x00ffffff;
    uint32_t x1_init = 0x01;

    std::vector<uint8_t> a(140);
    std::vector<uint8_t> b_tilde(108*16);
    std::vector<uint8_t> a_est(140);
    
    for(uint32_t log2AL = log2AL_min; log2AL <= log2AL_max; log2AL++){
        uint32_t AL = 1<<log2AL;
        if(verbose)
            std::cout<<"AL: "<<int(AL)<<" A:";
        uint32_t E = 108*AL;
        for(uint32_t A = A_min; A <= A_max[log2AL]; A++){
            // the hypotheses the search never tries: with so few redundant bits, a DCI scrambled with the RNTI
            // (case A) is sometimes also a solution of case B, decoded first, with another RNTI
            if(A+16+DCI_min_redundancy > E)
                continue;
            if(verbose)
                std::cout<<" "<<int(A);
            
            // test
            RNTI_CRC = random() & 0xffff;
            N_ID = random() & 0xffff;
            a.resize(A);
            for(uint32_t ia = 0; ia < A; ia++)
                a[ia] = random() & 0x01;
            crc_init = 0x00ffffff;
            x1_init = 0x01;
            
            bool is_case_A = random() & 1;
            
            // Case A:
            // - RNTI_SCR = RNTI_CRC & 0x7fff
            // - N_ID = DMRS_scrID
            // Case B:
            // - RNTI_SCR = 0x0000
            // - N_ID = cellID (== DMRS_scrID)
            RNTI_SCR = is_case_A ? RNTI_CRC & 0x7fff : 0x0000;
            encode(b_tilde,a,RNTI_CRC,RNTI_SCR,N_ID,E,crc_init,x1_init);
            
            uint32_t RNTI_CRC_est;
            
            // a noiseless DCI must be decoded, with its own message and RNTI
            bool is_equal = decode(RNTI_CRC_est, a_est, b_tilde, log2AL, A, N_ID) && RNTI_CRC_est == RNTI_CRC;
            for(uint32_t ir = 0; ir < A && is_equal; ir++)
                if(a_est[ir] != a[ir])
                    is_equal = false;
            if(!is_equal){
                is_ok = false;
                if(verbose)
                    std::cout<<"!"<<std::flush;
            }
            if(verbose)
                std::cout <<std::flush;
        }
        if(verbose)
            std::cout<<std::endl;
    }
    return is_ok;
}
