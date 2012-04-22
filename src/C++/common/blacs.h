// -*- mode: c++ -*-
#ifndef _BLACS_H_
#define _BLACS_H_

#include "index.h"
#include "import.h"

#define BLOCK_CYCLICAL_2D 1
#define DLEN_ 9
#define DTYPE_ 0
#define CTXT_ 1
#define M_ 2
#define N_ 3
#define MB_ 4
#define NB_ 5
#define RSRC_ 6
#define CSRC_ 7
#define LLD_ 8
#define ITHVAL_ 9

#ifdef _WIN32
#define blacs_pinfo_ BLACS_PINFO
#define blacs_setup_ BLACS_SETUP
#define blacs_gridinit_ BLACS_GRIDINIT
#define blacs_gridmap_ BLACS_GRIDMAP
#define blacs_abort_ BLACS_ABORT
#define blacs_gridexit_ BLACS_GRIDEXIT
#define blacs_barrier_ BLACS_BARRIER
#define blacs_gridinfo_ BLACS_GRIDINFO
#define blacs_pcoord_ BLACS_PCOORD
#define blacs_pnum_ BLACS_PNUM
#define blacs_get_ BLACS_GET
#define blacs_set_ BLACS_SET
#define blacs_exit_ BLACS_EXIT
#define numroc_ NUMROC
#endif

#ifdef __cplusplus
extern "C" 
{
#endif
    DLLIMPORT void blacs_pinfo_ (blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_setup_ (blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_gridinit_ (blas_idx_t &, const char *, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_gridmap_ (blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_abort_ (blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_gridexit_ (blas_idx_t &);
    DLLIMPORT void blacs_barrier_ (blas_idx_t &, char *);
    DLLIMPORT void blacs_gridinfo_ (blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_pcoord_ (blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_pnum_ (blas_idx_t &, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_get_ (blas_idx_t &, blas_idx_t &, blas_idx_t &);
    DLLIMPORT void blacs_set_ (blas_idx_t &, blas_idx_t &, blas_idx_t *);
    DLLIMPORT void blacs_exit_ (blas_idx_t &);
    DLLIMPORT blas_idx_t numroc_ (blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &, blas_idx_t &);

#ifdef __cplusplus
}
#endif


#endif
