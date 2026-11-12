#ifndef GSCN_H_
#define GSCN_H_

#include <vector>

// GSCN to Hz, [fc-BW/2,fc+BW/2] to GSCN_min:GSCN_max (38.104)
long long GSCN_to_freq(int GSCN);
void freq_range_to_GSCN(double fc, double BW, int mu, std::vector<int>& GSCN_range);

#endif
