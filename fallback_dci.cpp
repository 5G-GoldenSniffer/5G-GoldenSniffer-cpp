#include "fallback_dci.h"

// TS 38.214 Table 5.1.3.1-1 (MCS index table 1 for PDSCH): spectral efficiency Qm*R
static const double MCS_efficiency[29] = {
    0.2344, 0.3066, 0.3770, 0.4902, 0.6016, 0.7402, 0.8770, 1.0273, 1.1758, 1.3262,
    1.3281, 1.4766, 1.6953, 1.9141, 2.1602, 2.4063, 2.5703,
    2.5664, 2.7305, 3.0293, 3.3223, 3.6094, 3.9023, 4.2129, 4.5234, 4.8164, 5.1152, 5.3320, 5.5547
};

static unsigned read_msb(const uint8_t* bits, int n)
{
    unsigned v = 0;
    for(int i = 0; i < n; i++)
        v = (v<<1) | bits[i];
    return v;
}

static const uint16_t SI_RNTI = 0xffff, P_RNTI = 0xfffe;

bool DCI_is_DL(const Result& r)
{
    return r.RNTI == SI_RNTI || r.RNTI == P_RNTI || (r.DCI_len && r.DCI[0]);
}

bool decode_fallback(const Result& r, int BWP_size, Grant& g)
{
    if(BWP_size <= 0 || r.kind != RES_DCI || r.DCI_len == 0 || r.RNTI == SI_RNTI || r.RNTI == P_RNTI)
        return false;
    g.DL = r.DCI[0];
    // FDA size: ceil(log2(N(N+1)/2)) bits
    const unsigned N_RIV = BWP_size*(BWP_size+1)/2;
    int N_FDA = 0;
    while((1u<<N_FDA) < N_RIV)
        N_FDA++;
    if(g.DL){
        if(r.DCI_len != 28 + N_FDA)
            return false;
    }
    else{
        // 0_0 is zero-padded to the size of 1_0 when both are monitored in the same search space (38.212 7.3.1.0)
        if(r.DCI_len == 28 + N_FDA){
            for(int i = 20 + N_FDA; i < r.DCI_len; i++)
                if(r.DCI[i])
                    return false;
        }
        else if(r.DCI_len != 20 + N_FDA)
            return false;
    }

    // UL with frequency hopping: the N_UL_hop MSBs of the FDA are the hopping offset, the RIV is in the rest (38.214 6.1.2.2.2)
    int N_hop = 0;
    if(!g.DL && r.DCI[1+N_FDA+4])
        N_hop = BWP_size < 50 ? 1 : 2;

    // RIV (38.214 5.1.2.2.2, hDecodeRIV.m)
    const unsigned RIV = read_msb(r.DCI+1+N_hop, N_FDA-N_hop);
    if(RIV >= N_RIV)
        return false;
    int L = RIV/BWP_size + 1;
    int start = RIV - (L-1)*BWP_size;
    if(L > BWP_size - start){
        L = BWP_size - L + 2;
        start = BWP_size - 1 - start;
    }
    g.RB_start = start;
    g.L_RB = L;

    // MCS after the TDRA and VRB-to-PRB/hopping flag fields
    g.MCS = read_msb(r.DCI+1+N_FDA+5, 5);
    // RV after the NDI
    g.RV = read_msb(r.DCI+1+N_FDA+5+5+1, 2);
    // PDSCH over the rest of the slot after the CORESET, PUSCH over the whole slot
    int nDMRS = 3; // assume two DM-RS in the slot
    const int N_symb = 14 - nDMRS - (g.DL ? r.duration : 0);
    g.TBS = g.MCS <= 28 ? L*12*N_symb*MCS_efficiency[g.MCS] : 0;
    return true;
}
