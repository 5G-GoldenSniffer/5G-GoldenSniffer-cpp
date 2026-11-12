#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include "polar.h"

// interleaving pattern Pi_IL^max of 38.212 Table 5.3.1.1-1
static const uint8_t PI_IL_max[164] = {0, 2, 4, 7, 9, 14, 19, 20, 24, 25, 26, 28, 31, 34, 42,
    45, 49,    50, 51, 53, 54, 56, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 72, 76, 77, 
    81, 82, 83, 87, 88, 89, 91, 93, 95, 98, 101, 104, 106, 108, 110, 111, 113, 115,
    118, 119, 120, 122, 123, 126, 127, 129, 132, 134, 138, 139, 140, 1, 3, 5, 8, 10,
    15, 21, 27, 29, 32, 35, 43, 46, 52, 55, 57, 60, 63, 68, 73, 78, 84, 90, 92, 94,
    96, 99, 102, 105, 107, 109, 112, 114, 116, 121, 124, 128, 130, 133, 135, 141, 6,
    11, 16, 22, 30, 33, 36, 44, 47, 64, 74, 79, 85, 97, 100, 103, 117, 125, 131,
    136, 142, 12, 17, 23, 37, 48, 75, 80, 86, 137, 143, 13, 18, 38, 144, 39, 145,
    40, 146, 41, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163};

// polar sequence Q_0^{Nmax-1} of 38.212 Table 5.3.1.2-1, in increasing reliability
static const uint16_t Q_max[1024] = {
    0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 128, 12, 33, 65, 20, 256, 34, 
    24, 36, 7, 129, 66, 512, 11, 40, 68, 130, 19, 13, 48, 14, 72, 257, 21, 132, 35, 
    258, 26, 513, 80, 37, 25, 22, 136, 260, 264, 38, 514, 96, 67, 41, 144, 28, 69, 
    42, 516, 49, 74, 272, 160, 520, 288, 528, 192, 544, 70, 44, 131, 81, 50, 73, 
    15, 320, 133, 52, 23, 134, 384, 76, 137, 82, 56, 27, 97, 39, 259, 84, 138, 145, 
    261, 29, 43, 98, 515, 88, 140, 30, 146, 71, 262, 265, 161, 576, 45, 100, 640, 
    51, 148, 46, 75, 266, 273, 517, 104, 162, 53, 193, 152, 77, 164, 768, 268, 
    274, 518, 54, 83, 57, 521, 112, 135, 78, 289, 194, 85, 276, 522, 58, 168, 139, 
    99, 86, 60, 280, 89, 290, 529, 524, 196, 141, 101, 147, 176, 142, 530, 321, 
    31, 200, 90, 545, 292, 322, 532, 263, 149, 102, 105, 304, 296, 163, 92, 47, 
    267, 385, 546, 324, 208, 386, 150, 153, 165, 106, 55, 328, 536, 577, 548, 
    113, 154, 79, 269, 108, 578, 224, 166, 519, 552, 195, 270, 641, 523, 275, 
    580, 291, 59, 169, 560, 114, 277, 156, 87, 197, 116, 170, 61, 531, 525, 642, 
    281, 278, 526, 177, 293, 388, 91, 584, 769, 198, 172, 120, 201, 336, 62, 282, 
    143, 103, 178, 294, 93, 644, 202, 592, 323, 392, 297, 770, 107, 180, 151, 
    209, 284, 648, 94, 204, 298, 400, 608, 352, 325, 533, 155, 210, 305, 547, 
    300, 109, 184, 534, 537, 115, 167, 225, 326, 306, 772, 157, 656, 329, 110, 
    117, 212, 171, 776, 330, 226, 549, 538, 387, 308, 216, 416, 271, 279, 158, 
    337, 550, 672, 118, 332, 579, 540, 389, 173, 121, 553, 199, 784, 179, 228, 
    338, 312, 704, 390, 174, 554, 581, 393, 283, 122, 448, 353, 561, 203, 63, 
    340, 394, 527, 582, 556, 181, 295, 285, 232, 124, 205, 182, 643, 562, 286, 
    585, 299, 354, 211, 401, 185, 396, 344, 586, 645, 593, 535, 240, 206, 95, 
    327, 564, 800, 402, 356, 307, 301, 417, 213, 568, 832, 588, 186, 646, 404, 
    227, 896, 594, 418, 302, 649, 771, 360, 539, 111, 331, 214, 309, 188, 449, 
    217, 408, 609, 596, 551, 650, 229, 159, 420, 310, 541, 773, 610, 657, 333, 
    119, 600, 339, 218, 368, 652, 230, 391, 313, 450, 542, 334, 233, 555, 774, 
    175, 123, 658, 612, 341, 777, 220, 314, 424, 395, 673, 583, 355, 287, 183, 
    234, 125, 557, 660, 616, 342, 316, 241, 778, 563, 345, 452, 397, 403, 207, 
    674, 558, 785, 432, 357, 187, 236, 664, 624, 587, 780, 705, 126, 242, 565, 
    398, 346, 456, 358, 405, 303, 569, 244, 595, 189, 566, 676, 361, 706, 589, 
    215, 786, 647, 348, 419, 406, 464, 680, 801, 362, 590, 409, 570, 788, 597, 
    572, 219, 311, 708, 598, 601, 651, 421, 792, 802, 611, 602, 410, 231, 688, 
    653, 248, 369, 190, 364, 654, 659, 335, 480, 315, 221, 370, 613, 422, 425, 
    451, 614, 543, 235, 412, 343, 372, 775, 317, 222, 426, 453, 237, 559, 833, 
    804, 712, 834, 661, 808, 779, 617, 604, 433, 720, 816, 836, 347, 897, 243, 
    662, 454, 318, 675, 618, 898, 781, 376, 428, 665, 736, 567, 840, 625, 238, 
    359, 457, 399, 787, 591, 678, 434, 677, 349, 245, 458, 666, 620, 363, 127, 
    191, 782, 407, 436, 626, 571, 465, 681, 246, 707, 350, 599, 668, 790, 460, 
    249, 682, 573, 411, 803, 789, 709, 365, 440, 628, 689, 374, 423, 466, 793, 
    250, 371, 481, 574, 413, 603, 366, 468, 655, 900, 805, 615, 684, 710, 429, 
    794, 252, 373, 605, 848, 690, 713, 632, 482, 806, 427, 904, 414, 223, 663, 
    692, 835, 619, 472, 455, 796, 809, 714, 721, 837, 716, 864, 810, 606, 912, 
    722, 696, 377, 435, 817, 319, 621, 812, 484, 430, 838, 667, 488, 239, 378, 
    459, 622, 627, 437, 380, 818, 461, 496, 669, 679, 724, 841, 629, 351, 467, 
    438, 737, 251, 462, 442, 441, 469, 247, 683, 842, 738, 899, 670, 783, 849, 
    820, 728, 928, 791, 367, 901, 630, 685, 844, 633, 711, 253, 691, 824, 902, 
    686, 740, 850, 375, 444, 470, 483, 415, 485, 905, 795, 473, 634, 744, 852, 
    960, 865, 693, 797, 906, 715, 807, 474, 636, 694, 254, 717, 575, 913, 798, 
    811, 379, 697, 431, 607, 489, 866, 723, 486, 908, 718, 813, 476, 856, 839, 
    725, 698, 914, 752, 868, 819, 814, 439, 929, 490, 623, 671, 739, 916, 463, 
    843, 381, 497, 930, 821, 726, 961, 872, 492, 631, 729, 700, 443, 741, 845, 
    920, 382, 822, 851, 730, 498, 880, 742, 445, 471, 635, 932, 687, 903, 825, 
    500, 846, 745, 826, 732, 446, 962, 936, 475, 853, 867, 637, 907, 487, 695, 
    746, 828, 753, 854, 857, 504, 799, 255, 964, 909, 719, 477, 915, 638, 748, 
    944, 869, 491, 699, 754, 858, 478, 968, 383, 910, 815, 976, 870, 917, 727, 
    493, 873, 701, 931, 756, 860, 499, 731, 823, 922, 874, 918, 502, 933, 743, 
    760, 881, 494, 702, 921, 501, 876, 847, 992, 447, 733, 827, 934, 882, 937, 
    963, 747, 505, 855, 924, 734, 829, 965, 938, 884, 506, 749, 945, 966, 755, 
    859, 940, 830, 911, 871, 639, 888, 479, 946, 750, 969, 508, 861, 757, 970, 
    919, 875, 862, 758, 948, 977, 923, 972, 761, 877, 952, 495, 703, 935, 978, 
    883, 762, 503, 925, 878, 735, 993, 885, 939, 994, 980, 926, 764, 941, 967, 
    886, 831, 947, 507, 889, 984, 751, 942, 996, 971, 890, 509, 949, 973, 1000, 
    892, 950, 863, 759, 1008, 510, 979, 953, 763, 974, 954, 879, 981, 982, 927, 
    995, 765, 956, 887, 985, 997, 986, 943, 891, 998, 766, 511, 988, 1001, 951, 
    1002, 893, 975, 894, 1009, 955, 1004, 1010, 957, 983, 958, 987, 1012, 999, 
    1016, 767, 989, 1003, 990, 1005, 959, 1011, 1013, 895, 1006, 1014, 1017, 
    1018, 991, 1020, 1007, 1015, 1019, 1021, 1022, 1023};

uint32_t crc24c(const uint8_t *a, uint32_t A, uint32_t init)
{
    // 38.212 5.1: remainder of a(D)*D^24 modulo g_CRC24C(D), MSB first
    const uint32_t CRC24C = 0x01b2b117;
    uint32_t p = init;
    for(uint32_t k = 0; k < A+24; k++){
        p <<= 1;
        if(p & 0x01000000)
            p ^= CRC24C;
        if(k < A)
            p ^= a[k];
    }
    return p & 0x00ffffff;
}

PolarCode::PolarCode(uint32_t K, uint32_t E, uint32_t n_max, bool I_IL) : K(K), E(E)
{
    // 38.212 5.3.1: code length
    uint32_t n1 = std::ceil(std::log2(float(E)));
    if(16*E <= (9U<<n1) && 16*K < 9*E)
        n1--;
    uint32_t n2 = 3 + std::ceil(std::log2(float(K))); // R_min = 1/8
    const uint32_t n_min = 5;
    n = std::max(n_min, std::min(std::min(n1, n2), n_max));
    N = 1U<<n;

    // 5.3.1.1 interleaving
    PI.resize(K);
    if(I_IL){
        const uint32_t K_IL_max = 164;
        uint32_t k = 0;
        for(uint32_t m = 0; m < K_IL_max; m++)
            if(PI_IL_max[m] >= K_IL_max-K)
                PI[k++] = PI_IL_max[m] - (K_IL_max-K);
    }
    else
        for(uint32_t k = 0; k < K; k++)
            PI[k] = k;

    // 5.4.1.1 sub-block interleaver
    static const uint32_t P[32] = {0, 1, 2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19, 12,
                                   20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31};
    J.resize(N);
    for(uint32_t k = 0; k < N; k++)
        J[k] = P[(32*k)/N]*(N/32) + k%(N/32);

    // 5.4.1.2 bit selection
    if(E >= N)
        rm = REPETITION;
    else if(16*K <= 7*E)
        rm = PUNCTURING;
    else
        rm = SHORTENING;

    // 5.3.1.2: the K most reliable positions of Q_0^{N-1} outside Q_F,tmp
    std::vector<uint16_t> Q;
    Q.reserve(N);
    for(uint32_t i = 0; i < 1024; i++)
        if(Q_max[i] < N)
            Q.push_back(Q_max[i]);
    std::vector<bool> Q_Ftmp(N, false);
    if(rm == PUNCTURING){
        for(uint32_t k = 0; k < N-E; k++)
            Q_Ftmp[J[k]] = true;
        uint32_t k_max;
        if(4*E >= 3*N)
            k_max = std::ceil(float(3.0*N)/4.0 - float(E)/2.0);
        else
            k_max = std::ceil(float(9.0*N)/16.0 - float(E)/4.0);
        for(uint32_t k = 0; k < k_max; k++)
            Q_Ftmp[k] = true;
    }
    else if(rm == SHORTENING)
        for(uint32_t k = E; k < N; k++)
            Q_Ftmp[J[k]] = true;
    frozen.assign(N, 1);
    uint32_t k = 0;
    for(int i = N-1; k < K; i--)
        if(!Q_Ftmp[Q[i]]){
            frozen[Q[i]] = 0;
            k++;
        }
    for(uint32_t i = 0; i < N; i++)
        if(!frozen[i])
            info.push_back(i);
}

void PolarCode::encode(const uint8_t *c, uint8_t *e) const
{
    std::vector<uint8_t> d(N, 0);
    for(uint32_t k = 0; k < K; k++)
        d[info[k]] = c[PI[k]];
    // d = u*G_N
    for(uint32_t half = 1; half < N; half <<= 1)
        for(uint32_t j = 0; j < N; j += 2*half)
            for(uint32_t k = 0; k < half; k++)
                d[j+k] ^= d[j+k+half];
    // e_k = y_{...} = d[J[...]]
    for(uint32_t k = 0; k < E; k++){
        uint32_t i = rm == REPETITION ? k%N : rm == PUNCTURING ? k+N-E : k;
        e[k] = d[J[i]];
    }
}

namespace {
// state of a decoding path: the LLRs (alpha) and the codeword of the last decoded left child (beta) of the node
// at every level s < n of the tree, stored at offset 2^s-1 (level n is the channel)
struct Path {
    std::vector<float> alpha;
    std::vector<uint8_t> beta;
    std::vector<uint8_t> u; // information bits decided so far
    float metric;
};
}

void PolarCode::decode_list(const float *llr_e, int L, std::vector<std::vector<uint8_t>> &c_list) const
{
    // rate recovery: LLRs of the codeword d
    const float known = 1e20f;
    std::vector<float> llr(N, 0.0f);
    for(uint32_t k = 0; k < E; k++){
        uint32_t i = rm == REPETITION ? k%N : rm == PUNCTURING ? k+N-E : k;
        llr[J[i]] += llr_e[k];
    }
    if(rm == SHORTENING)
        for(uint32_t i = E; i < N; i++)
            llr[J[i]] = known; // shortened bits are 0

    std::vector<Path> paths(1), next;
    paths[0].alpha.resize(N);
    paths[0].beta.resize(N);
    paths[0].metric = 0.0f;
    std::vector<uint8_t> c0(N), c1(N);

    // LLR of u_i for path p, from the root down to the leaf
    auto descend = [&](Path &p, uint32_t i){
        int start = i == 0 ? n-1 : __builtin_ctz(i);
        for(int s = start; s >= 0; s--){
            const uint32_t h = 1U<<s;
            const float *P = (uint32_t)s+1 == n ? llr.data() : &p.alpha[(2*h)-1];
            float *C = &p.alpha[h-1];
            if(((i>>s) & 1) == 0){ // left child: f
                for(uint32_t k = 0; k < h; k++){
                    float a = P[k], b = P[k+h];
                    float m = std::min(std::fabs(a), std::fabs(b));
                    C[k] = (a < 0) != (b < 0) ? -m : m;
                }
            }
            else{ // right child: g, with the codeword of the left sibling
                const uint8_t *B = &p.beta[h-1];
                for(uint32_t k = 0; k < h; k++)
                    C[k] = P[k+h] + (B[k] ? -P[k] : P[k]);
            }
        }
    };
    // partial sums: combine the decided bit u_i with the left siblings up to the first left child
    auto ascend = [&](Path &p, uint32_t i, uint8_t bit){
        uint8_t *c = c0.data(), *t = c1.data();
        c[0] = bit;
        uint32_t s = 0;
        for(; s < n && ((i>>s) & 1); s++){
            const uint32_t h = 1U<<s;
            const uint8_t *B = &p.beta[h-1];
            for(uint32_t k = 0; k < h; k++){
                t[k] = B[k] ^ c[k];
                t[k+h] = c[k];
            }
            std::swap(c, t);
        }
        if(s < n)
            std::copy(c, c+(1U<<s), &p.beta[(1U<<s)-1]);
    };

    struct Fork { int path; uint8_t bit; float metric; };
    std::vector<Fork> forks;
    for(uint32_t i = 0; i < N; i++){
        for(auto &p : paths)
            descend(p, i);
        if(frozen[i]){
            for(auto &p : paths){
                if(p.alpha[0] < 0)
                    p.metric -= p.alpha[0];
                ascend(p, i, 0);
            }
            continue;
        }
        // every path forks into u_i = 0 and u_i = 1: keep the L with the smallest metric
        forks.clear();
        for(size_t ip = 0; ip < paths.size(); ip++){
            const float a = paths[ip].alpha[0];
            forks.push_back({(int)ip, 0, paths[ip].metric + (a < 0 ? -a : 0.0f)});
            forks.push_back({(int)ip, 1, paths[ip].metric + (a > 0 ?  a : 0.0f)});
        }
        const size_t n_keep = std::min(forks.size(), (size_t)L);
        std::partial_sort(forks.begin(), forks.begin()+n_keep, forks.end(),
                          [](const Fork &a, const Fork &b){ return a.metric < b.metric; });
        std::vector<int> uses(paths.size(), 0);
        for(size_t f = 0; f < n_keep; f++)
            uses[forks[f].path]++;
        next.clear();
        for(size_t f = 0; f < n_keep; f++){
            const Fork &fk = forks[f];
            if(--uses[fk.path] == 0)
                next.push_back(std::move(paths[fk.path])); // last user of the path
            else
                next.push_back(paths[fk.path]);
            Path &p = next.back();
            p.metric = fk.metric;
            p.u.push_back(fk.bit);
            ascend(p, i, fk.bit);
        }
        paths.swap(next);
    }

    std::sort(paths.begin(), paths.end(), [](const Path &a, const Path &b){ return a.metric < b.metric; });
    c_list.assign(paths.size(), std::vector<uint8_t>(K));
    for(size_t ip = 0; ip < paths.size(); ip++)
        for(uint32_t k = 0; k < K; k++)
            c_list[ip][PI[k]] = paths[ip].u[k];
}
