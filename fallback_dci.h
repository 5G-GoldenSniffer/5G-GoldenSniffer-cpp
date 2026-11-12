#ifndef FALLBACK_DCI_H_
#define FALLBACK_DCI_H_

#include "pipeline.h"

// resource allocation of a fallback DCI and the transport block size it implies
struct Grant {
    bool DL;    // format 1_0 (PDSCH) or 0_0 (PUSCH)
    int RB_start;
    int L_RB;
    int MCS;
    int RV;     // redundancy version: 0 for new data
    double TBS; // approximate TBS [bits], 0 if MCS > 28 (retransmission: TBS of the initial transmission)
};

// downlink DCI: SI-RNTI and P-RNTI only schedule the PDSCH, with a format 1_0 that has no identifier
// field (38.212 7.3.1.2.1); for the other RNTIs the first DCI bit is the identifier (1: downlink)
bool DCI_is_DL(const Result& r);

// fallback DCI in a BWP of BWP_size RBs (38.212 7.3.1.1.1 and 7.3.1.2.1):
//   DL, format 1_0: identifier (1 bit, 1), FDA (RIV), TDRA (4), VRB-to-PRB (1), MCS (5), NDI (1), RV (2), ... 28 + FDA bits
//   UL, format 0_0: identifier (1 bit, 0), FDA (RIV), TDRA (4), hopping flag (1), MCS (5), NDI (1), RV (2), ... 20 + FDA bits,
//                   or 28 + FDA bits zero-padded to the size of 1_0
// False if r isn't a DCI of either size or the RIV is invalid, and for SI-RNTI and P-RNTI, whose format 1_0
// has another layout and whose FDA refers to CORESET#0 rather than to this BWP.
bool decode_fallback(const Result& r, int BWP_size, Grant& g);

#endif
