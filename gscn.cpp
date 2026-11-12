#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>
#include "gscn.h"

// 38.104
long long GSCN_to_freq(int GSCN)
{
    if (GSCN >=2 && GSCN <= 7498){
    int N = (GSCN+1)/3;
    int M = 3+2*(GSCN-3*N);
    return N*1200000LL+M*50000LL;
    }
    else if(GSCN>=7499 && GSCN<=22255){
        return 3000000000LL+(GSCN-7499)*1440000LL;
    }
    else if(GSCN>=22256 && GSCN<=26639){
        return 24250080000LL+(GSCN-22256)*17280000LL;
    }
    return 0LL;
}

void freq_range_to_GSCN(double fc, double BW, int mu, std::vector<int>& GSCN_range)
{
    //std::cout<<"freq_range_to_GSCN(fc:"<<fc<<",BW:"<<BW<<")"<<std::endl;
    double SCS = 15000 << mu;
    // PSS/SSS subcarriers: 56:182
    // GSCN == subcarrier 120 => -64:62
    int GSCN_min, GSCN_max;
    
    double f_min = fc - BW/2 + 64*SCS - SCS;
    if(f_min < 3e9){
        // K_min = {0,1,2}
        // N_min = ceil((f_min-50e3-100e3*K_min)/1.2e6);
        // *GSCN_min = 3*N_min+K_min-1;
        int min_val = 7499;
        for(int K_min = 0; K_min <3; K_min++)
            min_val = std::min(min_val,3*(int)ceil((f_min-50e3-100e3*K_min)/1.2e6)+K_min-1);
        GSCN_min = min_val;
    }
    else if(f_min < 24250e6)
        GSCN_min = 7499+ceil((f_min-3e9)/1.44e6);
    else
        GSCN_min = 22256+ceil((f_min-24250.08e6)/17.28e6);

    double f_max = fc + BW/2 - 62*SCS;
    if(f_max < 3e9){
        // K_max = {0,1,2}
        // N_max = floor((f_max-50e3-100e3*K_max)/1.2e6);
        // *GSCN_max = 3*N_max+K_max-1;
        int max_val = 0;
        for(int K_max = 0; K_max < 3; K_max++)
            max_val = std::max(max_val,3*(int)floor((f_max-50e3-100e3*K_max)/1.2e6)+K_max-1);
        GSCN_max = max_val;
    }
    else if(f_max < 24250e6)
        GSCN_max = 7499+floor((f_max-3e9)/1.44e6);
    else
        GSCN_max = 22256+floor((f_max-24250.08e6)/17.28e6);

    // every GSCN of the band (an earlier list of known operator GSCNs is in the git history)
    if(f_max < 3e9){
//        if((GSCN_min%3)!=0)
//            GSCN_min += 3-(GSCN_min%3);
//        if((GSCN_max%3)!=0)
//            GSCN_max -= (GSCN_max%3);
//        GSCN_range.resize((GSCN_max-GSCN_min)/3+1);
//        for(size_t n = 0; n < GSCN_range.size(); n++)
//            GSCN_range[n] = GSCN_min + 3*n;
        GSCN_range.resize(GSCN_max-GSCN_min+1);
        for(size_t n = 0; n < GSCN_range.size(); n++)
            GSCN_range[n] = GSCN_min + n;
    }
    else if(f_min > 3e9){
        GSCN_range.resize(GSCN_max-GSCN_min+1);
        for(size_t n = 0; n < GSCN_range.size(); n++)
            GSCN_range[n] = GSCN_min + n;
    }
}
