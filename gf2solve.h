#ifndef _GF2SOLVE_H
#define _GF2SOLVE_H 1

#include <cstdint>

typedef struct{
    int n;
    int *d;
} mzp_t;

mzp_t *mzp_init(int n);

typedef struct{
    int n;
    uint8_t *d;
} mzv_t;

mzv_t *mzv_init(int n);

typedef struct{ 
    int nr; // rows
    int nc; // columns
    int nw; // 64 bit words in a row
    uint64_t *d; // 2 * nr * nw words, the first nr rows are the matrix, the next nr hold the L of the LU factorization
} mzd_t;

mzd_t *mzd_init(int nr, int nc);
bool mzd_read_bit(mzd_t const *A, int ir, int ic);
void mzd_write_bit(mzd_t *A, int ir, int ic, bool b);
void mzd_swap_rows(mzd_t *A, int ir1, int ir2, int num_col);
int mzd_plu(mzd_t *A, mzp_t *P);
// most 64 bit words in a row of the matrices solved by mzdv_lu_solve_perm (up to 256 columns)
#define MZD_MAX_NW 4
bool mzdv_lu_solve_perm(mzd_t const *A, mzp_t const *P, uint8_t const *y, uint8_t *x);

#endif
