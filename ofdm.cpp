#include <cmath>
#include "ofdm.h"

// single layer
void DMRS_demod(std::vector<uint8_t> &c_est, int rot, cfloat *C_est, int N_DMRS, int sc_spacing, int sc_offset)
{
    switch(rot){
    case 0: // C_est is not rotated with respect to the true C
        for(int i = 0; i < N_DMRS; i++){
            c_est[2*i] = C_est[i*sc_spacing+sc_offset].real()<0; // x
            c_est[2*i+1] = C_est[i*sc_spacing+sc_offset].imag()<0; // y
        }
        break;
    case 1: // C_est is rotated by pi/2
        for(int i = 0; i < N_DMRS; i++){
            c_est[2*i] = -C_est[i*sc_spacing+sc_offset].imag()<0; // -y
            c_est[2*i+1] = C_est[i*sc_spacing+sc_offset].real()<0; // x
        }
        break;
    case 2: // C_est is rotated by pi
        for(int i = 0; i < N_DMRS; i++){
            c_est[2*i] = -C_est[i*sc_spacing+sc_offset].real()<0; // -x
            c_est[2*i+1] = -C_est[i*sc_spacing+sc_offset].imag()<0; // -y
        }
        break;
    case 3: // C_est is rotated by -pi/2
        for(int i = 0; i < N_DMRS; i++){
            c_est[2*i] = C_est[i*sc_spacing+sc_offset].imag()<0; // y
            c_est[2*i+1] = -C_est[i*sc_spacing+sc_offset].real()<0; // -x
        }
        break;
    }
}

// the data subcarriers for a PDCCH with the DMRS above
void PDCCH_demod(uint8_t *b_tilde, int rot, cfloat *C_est, int N_DMRS, int sc_spacing, int sc_offset)
{
    switch(rot){
    case 0: // C_est is not rotated with respect to the true C
        for(int i = 0; i < N_DMRS; i++){
            b_tilde[6*i  ] = C_est[i*sc_spacing+sc_offset-1].real()<0; // x
            b_tilde[6*i+1] = C_est[i*sc_spacing+sc_offset-1].imag()<0; // y
            b_tilde[6*i+2] = C_est[i*sc_spacing+sc_offset+1].real()<0; // x
            b_tilde[6*i+3] = C_est[i*sc_spacing+sc_offset+1].imag()<0; // y
            b_tilde[6*i+4] = C_est[i*sc_spacing+sc_offset+2].real()<0; // x
            b_tilde[6*i+5] = C_est[i*sc_spacing+sc_offset+2].imag()<0; // y
        }
        break;
    case 1: // C_est is rotated by pi/2
        for(int i = 0; i < N_DMRS; i++){
            b_tilde[6*i  ] = -C_est[i*sc_spacing+sc_offset-1].imag()<0; // -y
            b_tilde[6*i+1] =  C_est[i*sc_spacing+sc_offset-1].real()<0; // x
            b_tilde[6*i+2] = -C_est[i*sc_spacing+sc_offset+1].imag()<0; // -y
            b_tilde[6*i+3] =  C_est[i*sc_spacing+sc_offset+1].real()<0; // x
            b_tilde[6*i+4] = -C_est[i*sc_spacing+sc_offset+2].imag()<0; // -y
            b_tilde[6*i+5] =  C_est[i*sc_spacing+sc_offset+2].real()<0; // x
        }
        break;
    case 2: // C_est is rotated by pi
        for(int i = 0; i < N_DMRS; i++){
            b_tilde[6*i  ] = -C_est[i*sc_spacing+sc_offset-1].real()<0; // -x
            b_tilde[6*i+1] = -C_est[i*sc_spacing+sc_offset-1].imag()<0; // -y
            b_tilde[6*i+2] = -C_est[i*sc_spacing+sc_offset+1].real()<0; // -x
            b_tilde[6*i+3] = -C_est[i*sc_spacing+sc_offset+1].imag()<0; // -y
            b_tilde[6*i+4] = -C_est[i*sc_spacing+sc_offset+2].real()<0; // -x
            b_tilde[6*i+5] = -C_est[i*sc_spacing+sc_offset+2].imag()<0; // -y
        }
        break;
    case 3: // C_est is rotated by -pi/2
        for(int i = 0; i < N_DMRS; i++){
            b_tilde[6*i  ] =  C_est[i*sc_spacing+sc_offset-1].imag()<0; // y
            b_tilde[6*i+1] = -C_est[i*sc_spacing+sc_offset-1].real()<0; // -x
            b_tilde[6*i+2] =  C_est[i*sc_spacing+sc_offset+1].imag()<0; // y
            b_tilde[6*i+3] = -C_est[i*sc_spacing+sc_offset+1].real()<0; // -x
            b_tilde[6*i+4] =  C_est[i*sc_spacing+sc_offset+2].imag()<0; // y
            b_tilde[6*i+5] = -C_est[i*sc_spacing+sc_offset+2].real()<0; // -x
        }
        break;
    }
}

// writes into C_est only the points Y[sym*N_FFT+(k_ini:k_step:k_fin)], unwrapped with the minimum phase increment
// - if we want all the points, as for a non-interleaved DCI, k_step=1
// - for an NZP-CSI-RS with density==3: k_step=4, density==1: k_step=12, density==0.5: k_step=24
//
void unwrap_QPSK(cfloat *C_est, const cfloat *Y, size_t sym, size_t N_FFT, int k_ini, int k_fin, int k_step)
{
    // phase unwrap of the entire DCI, as it is QPSK modulated
    cfloat Ycorr = Y[sym*N_FFT + k_ini];
    C_est[0] = cfloat(1.0,1.0);
    int k = k_ini + k_step;
    int i = 1;
    while(k < k_fin){
        cfloat tmp = Y[sym*N_FFT + k]*std::conj(Ycorr);
        if(tmp.real() <= -std::abs(tmp.imag())){ // abs(arg(tmp)) >= 0.75*M_PI
            Ycorr = -Y[sym*N_FFT + k]; // -1*Y
            C_est[i] = cfloat(-1.0,-1.0); // -1*exp(1i*pi/4)
        }
        else if(tmp.imag() < -std::abs(tmp.real())){ // arg(tmp) < -0.25*M_PI
            Ycorr = cfloat(-Y[sym*N_FFT + k].imag(),Y[sym*N_FFT + k].real()); // 1i*Y
            C_est[i] = cfloat(1.0,-1.0); // -1i*exp(1i*pi/4)
        }
        else if(tmp.imag() > std::abs(tmp.real())){ // arg(tmp) > 0.25*M_PI
            Ycorr = cfloat(Y[sym*N_FFT + k].imag(),-Y[sym*N_FFT + k].real()); // -1i*Y
            C_est[i] = cfloat(-1.0,1.0); // 1i*exp(1i*pi/4)
        }
        else{ // abs(arg(tmp)) <= 0.25*M_PI
            Ycorr = Y[sym*N_FFT + k]; // 1*Y
            C_est[i] = cfloat(1.0,1.0); // 1*exp(1i*pi/4)
        }
        k += k_step;
        i++;
    }
}

void ofdm_fft(size_t sym, cfloat *Y,float phase_comp,float Cfo_eF,size_t sym_ofs,double Ts,size_t N_FFT,fftwf_plan Y_symb_to_freq,cfloat *Y_tmpbuf,size_t CP_slack)
{
    cfloat w = std::exp(cfloat(0.0,phase_comp-2.0*M_PI*Cfo_eF*sym_ofs*Ts));
    cfloat wi = std::exp(cfloat(0.0,-2.0*M_PI*Cfo_eF*Ts));
    for(unsigned n = 0; n < N_FFT; n++){
        Y[sym*N_FFT+n] *= w;
        w *= wi;
    }
    fftwf_execute(Y_symb_to_freq);
    for(int n = 0; n < (int)N_FFT/2; n++) // positive frequencies
        Y[sym*N_FFT+N_FFT/2+n] = Y_tmpbuf[n];
    for(int n = (int)N_FFT/2; n < (int)N_FFT; n++) // negative frequencies
        Y[sym*N_FFT-N_FFT/2+n] = Y_tmpbuf[n];

    w = std::exp(cfloat(0.0,-M_PI*CP_slack)); // 2*pi*(-N_FFT/2)*CP_slack/N_FFT
    wi = std::exp(cfloat(0.0,(float)(2.0*M_PI*CP_slack)/(float)N_FFT));
    for(int n = 0; n < (int)N_FFT; n++){
        Y[sym*N_FFT + n] *= w;
        w *= wi;
    }
}
