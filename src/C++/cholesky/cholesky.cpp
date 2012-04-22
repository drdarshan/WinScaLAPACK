#include <mpi.h>
#include <cstdio>
#include <cassert>
#include "block_cyclic_mat.h"
#include "scalapack.h"

static std::shared_ptr<block_cyclic_mat_t> make_tridiagonal(std::shared_ptr<blacs_grid_t> grid, blas_idx_t n_global)
{
    // First create a matrix with 2 on the diagonal
    auto a = block_cyclic_mat_t::diagonal(grid, n_global, n_global, 2.0);
    
    // Then set the off-diagonal entries to -1
    // See: http://icl.cs.utk.edu/lapack-forum/archives/scalapack/msg00055.html
    char uplo        = 'L';
    blas_idx_t n     = n_global - 1;
    blas_idx_t ia    = 2;
    blas_idx_t ja    = 1;
    double zero      =  0.0;
    double minus_one = -1.0;
    pdlaset_(uplo, n, n, zero, minus_one, a -> local_data(), ia, ja, a -> descriptor());

    uplo = 'U';
    ia = 1;
    ja = 2;
    pdlaset_(uplo, n, n, zero, minus_one, a -> local_data(), ia, ja, a -> descriptor());

    return a;
}

static double potrf_flops(blas_idx_t N)
{
    // From: 
    // https://icl.cs.utk.edu/svn/scalapack-dev/scalapack/trunk/TESTING/LIN/pdlltdriver.f
    return (1.0/3.0 * N * N * N + 1.0/2.0 * N * N)/1024.0/1024.0/1024.0;
}

static void chol_driver(blas_idx_t n_global)
{
    auto grid = std::make_shared<blacs_grid_t>();    
    auto a    = make_tridiagonal(grid, n_global);    

    // Compute Cholesky factorization of A in-place
    char       uplo     ='U';
    blas_idx_t ia       = 1, ja = 1, info;

    MPI_Barrier (MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    pdpotrf_ (uplo, n_global, a->local_data(), ia, ja, a->descriptor(), info);
    assert(info == 0);

    double t1 = MPI_Wtime() - t0;
  
    double t_glob;
    MPI_Reduce(&t1, &t_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (grid->iam() == 0) 
    {
        double gflops = potrf_flops(n_global)/t_glob/grid->nprocs();
        printf("\n"
            "MATRIX CHOLESKY FACTORIZATION BENCHMARK SUMMARY\n"
            "===============================================\n"
            "N = %d\tNP = %d\tNP_ROW = %d\tNP_COL = %d\n"
            "Time for PxPOTRF = %10.7f seconds\tGflops/Proc = %10.7f\n",
            n_global, grid->nprocs(), grid->nprows(), grid->npcols(), 
            t_glob, gflops);fflush(stdout);
    }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  blas_idx_t n_global = 4096;
  
  if (argc > 1)
  {
    n_global = blas_idx_t(atol(argv[1]));
  }

  chol_driver(n_global);
  MPI_Finalize();
}
