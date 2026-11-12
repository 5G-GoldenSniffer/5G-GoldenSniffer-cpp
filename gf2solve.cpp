#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "gf2solve.h"

mzp_t *mzp_init(int n)
{
    mzp_t *p = (mzp_t *)malloc(sizeof(mzp_t));
    p->n = n;
    p->d = (int *)malloc(n*sizeof(int));
    for(int i = 0; i < n; i++)
        p->d[i] = i;
    return p;
}

mzv_t *mzv_init(int n){
    mzv_t *v = (mzv_t *)malloc(sizeof(mzv_t));
    v->n = n;
    v->d = (uint8_t *)calloc(n,sizeof(uint8_t));
    return v;
}

mzd_t *mzd_init(int nr, int nc){
    int nw = nc>>6;
    if((nw<<6) < nc) nw++;
    mzd_t *m = (mzd_t *)malloc(sizeof(mzd_t));
    m->nr = nr;
    m->nc = nc;
    m->nw = nw;
    m->d = (uint64_t *)calloc(2*nr*nw,sizeof(uint64_t));
    return m;
}

bool mzd_read_bit(mzd_t const *A, int ir, int ic)
{
    assert(ir >= 0 && ir < 2*A->nr && ic >= 0 && ic < A->nc);
    int iw = ic >> 6;
    return (A->d[ir*A->nw + iw] >> (ic - (iw << 6))) & 0x01;
}

void mzd_write_bit(mzd_t *A, int ir, int ic, bool b)
{
    assert(ir >= 0 && ir < 2*A->nr && ic >= 0 && ic < A->nc);
    int iw = ic >> 6;
    if (b)
        A->d[ir*A->nw + iw] |=  (1ull << (ic - (iw << 6)));
    else
        A->d[ir*A->nw + iw] &= ~(1ull << (ic - (iw << 6)));
}

void mzd_swap_rows(mzd_t *A, int ir1, int ir2, int num_col)
{
    assert(ir1 >= 0 && ir1 < 2*A->nr && ir2 >= 0 && ir2 < 2*A->nr && num_col >= 0 && num_col <= A->nc);
    uint64_t tmp_uint64;
    int nw, nb_last_word;
    if(num_col == A->nc){
        nw = A->nw;
        nb_last_word = 0;
    }
    else{
        nw = num_col >> 6; 
        nb_last_word = num_col - (nw<<6);
    }
    for(int iw = 0; iw < nw; iw++){
        tmp_uint64 = A->d[ir1*A->nw + iw];
        A->d[ir1*A->nw + iw] = A->d[ir2*A->nw + iw];
        A->d[ir2*A->nw + iw] = tmp_uint64;
    }
    if(nb_last_word){
        uint64_t mask = (1ull << nb_last_word) - 1ull;
        tmp_uint64 = A->d[ir1*A->nw + nw] & mask;
        A->d[ir1*A->nw + nw] = (A->d[ir1*A->nw + nw] & ~mask) | (A->d[ir2*A->nw + nw] & mask);
        A->d[ir2*A->nw + nw] = (A->d[ir2*A->nw + nw] & ~mask) | tmp_uint64;
    }
}

int mzd_plu(mzd_t *A, mzp_t *P)
{
    int ir,ic,iw;
    // set the L part to the identity matrix
    for(ir = 0; ir < A->nr; ir++){
        for(iw = 0; iw < A->nw; iw++)
            A->d[(ir + A->nr)*A->nw + iw] = 0ull;
        if(ir < A->nc)
            mzd_write_bit(A,ir + A->nr,ir,1);
    }
    for(ic = 0; ic < A->nc; ic++){
        bool pivot_found = mzd_read_bit(A,ic,ic);
        if(!pivot_found){
            // find the pivot and swap the rows in A and in P
            for(ir = ic+1; ir < A->nr; ir++){
                if(mzd_read_bit(A,ir,ic)){
                    mzd_swap_rows(A,ir,ic,A->nc);
                    mzd_swap_rows(A,ir+A->nr,ic+A->nr,ic);
                    int tmp_int = P->d[ir];
                    P->d[ir] = P->d[ic];
                    P->d[ic] = tmp_int;
                    pivot_found = true;
                    break;
                }
            }
        }
        if(!pivot_found)
            break;
        for(ir = ic+1; ir < A->nr; ir++){
            if(mzd_read_bit(A,ir,ic)){
                // add row ic of A to row ir
                for(iw = 0; iw < A->nw; iw++)
                    A->d[ir*A->nw + iw] ^= A->d[ic*A->nw + iw];
                // L(ir,ic) = 1
                mzd_write_bit(A,ir+A->nr,ic,1);
            }
        }
    }
    assert(ic == A->nc);
    return ic;
}

// solves A x = P y with the factorization of mzd_plu, a 64 bit word of the rows at a time; y is not changed.
// Returns false at the first row of P y inconsistent with the solution of the first nc ones, otherwise x[0..nc-1] = x
bool mzdv_lu_solve_perm(mzd_t const *A, mzp_t const *P, uint8_t const *y, uint8_t *x)
{
    assert(A->nw <= MZD_MAX_NW && P->n == A->nr);
    const int nw = A->nw;
    const uint64_t *L = A->d + (size_t)A->nr*nw;
    // forward substitution, packed in z: z(ir) = y(P(ir)) + L(ir,1:ir-1)*z(1:ir-1)
    uint64_t z[MZD_MAX_NW] = {0};
    for(int ir = 0; ir < A->nc; ir++){
        const uint64_t *l = L + (size_t)ir*nw;
        uint64_t s = 0;
        for(int iw = 0; iw <= (ir>>6); iw++)
            s ^= l[iw] & z[iw];
        z[ir>>6] |= (uint64_t)((y[P->d[ir]] ^ __builtin_parityll(s)) & 1) << (ir & 63);
    }
    // the other rows must be consistent: y(P(ir)) + L(ir,1:nc)*z(1:nc) = 0
    for(int ir = A->nc; ir < A->nr; ir++){
        const uint64_t *l = L + (size_t)ir*nw;
        uint64_t s = 0;
        for(int iw = 0; iw < nw; iw++)
            s ^= l[iw] & z[iw];
        if((y[P->d[ir]] ^ __builtin_parityll(s)) & 1)
            return false;
    }
    // back substitution, packed in xw: x(ir) = z(ir) + U(ir,ir+1:nc)*x(ir+1:nc)
    uint64_t xw[MZD_MAX_NW] = {0};
    for(int ir = A->nc-1; ir >= 0; ir--){
        const uint64_t *u = A->d + (size_t)ir*nw;
        uint64_t s = 0;
        for(int iw = ir>>6; iw < nw; iw++)
            s ^= u[iw] & xw[iw];
        const uint64_t bit = ((z[ir>>6] >> (ir & 63)) ^ __builtin_parityll(s)) & 1;
        xw[ir>>6] |= bit << (ir & 63);
        x[ir] = bit;
    }
    return true;
}
