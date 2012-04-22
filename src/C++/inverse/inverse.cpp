#include <mpi.h>
#include <cstdio>
#include <cassert>
#include "block_cyclic_mat.h"
#include "scalapack.h"

static double getri_flops(blas_idx_t N)
{
    // From:
    // https://icl.cs.utk.edu/svn/scalapack-dev/scalapack/trunk/TESTING/LIN/pdinvdriver.f    
    // Factorization: 2/3 N^3 - 1/2 N^2
    // Inverse      : 4/3 N^3 - N^2

    return ((2.0/3.0 * N * N * N) - (1.0/2.0 * N * N) + 
        (4.0/3.0 * N * N * N) - (N * N))/1024.0/1024.0/1024.0;
}

static void inv_driver(blas_idx_t n_global)
{
    auto grid = std::make_shared<blacs_grid_t>();

    // Create a NxN random matrix A
    auto a = block_cyclic_mat_t::random(grid, n_global, n_global);        

    // Create a NxN matrix to hold A^{-1}
    auto ai = block_cyclic_mat_t::constant(grid, n_global, n_global);

    // Copy A to A^{-1} since it will be overwritten during factorization
    std::copy_n(a->local_data(), a->local_size(), ai->local_data());

    MPI_Barrier (MPI_COMM_WORLD);

    double t0 = MPI_Wtime();
    
    // Factorize A 
    blas_idx_t ia = 1, ja = 1;
    std::vector<blas_idx_t> ipiv(a->local_rows() + a->row_block_size() + 100);
    blas_idx_t info;
    pdgetrf_(n_global, n_global, 
        ai->local_data(), ia, ja, ai->descriptor(), 
        ipiv.data(), 
        info);
    assert(info == 0);
    double t_factor = MPI_Wtime() - t0;

    // Compute A^{-1} based on the LU factorization

    // Compute workspace for double and integer work arrays on each process
    blas_idx_t lwork  = 10;
    blas_idx_t liwork = 10;
    std::vector<double>     work (lwork); 
    std::vector<blas_idx_t> iwork(liwork);

    lwork = liwork = -1;    
    pdgetri_(n_global, 
        ai->local_data(), ia, ja, ai->descriptor(), 
        ipiv.data(), 
        work.data(), lwork, iwork.data(), liwork, info);
    assert(info == 0);
    lwork  = static_cast<blas_idx_t>(work[0]);
    liwork = static_cast<size_t>(iwork[0]);
    work.resize(lwork);
    iwork.resize(liwork);

    // Now compute the inverse
    t0 = MPI_Wtime();
    pdgetri_(n_global, 
        ai->local_data(), ia, ja, ai->descriptor(), 
        ipiv.data(), 
        work.data(), lwork, iwork.data(), liwork, info);
    assert(info == 0);
    double t_solve = MPI_Wtime() - t0;

    // Verify that the inverse is correct using A*A^{-1} = I
    auto identity = block_cyclic_mat_t::diagonal(grid, n_global, n_global);

    // Compute I = A * A^{-1} - I and verify that the ||I|| is small    
    char nein = 'N';
    double alpha = 1.0, beta = -1.0;
    pdgemm_(nein, nein, n_global, n_global, n_global, alpha, 
        a->local_data() , ia, ja, a->descriptor(),
        ai->local_data(), ia, ja, ai->descriptor(),
        beta,
        identity->local_data(), ia, ja, identity->descriptor());

    // Compute 1-norm of the result
    char norm='1';
    work.resize(identity->local_cols());
    double err = pdlange_(norm, n_global, n_global, 
        identity->local_data(), ia, ja, identity->descriptor(), work.data());

    double t_total = t_factor + t_solve;
    double t_glob;
    MPI_Reduce(&t_total, &t_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (grid->iam() == 0) 
    {
        double gflops = getri_flops(n_global)/t_glob/grid->nprocs();
        printf("\n"
            "MATRIX INVERSE BENCHMARK SUMMARY\n"
            "================================\n"
            "N = %d\tNP = %d\tNP_ROW = %d\tNP_COL = %d\n"
            "Time for PxGETRF + PxGETRI = %10.7f seconds\tGflops/Proc = %10.7f, Error = %f\n",
            n_global, grid->nprocs(), grid->nprows(), grid->npcols(), 
            t_glob, gflops, err);fflush(stdout);
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

    inv_driver(n_global);
    MPI_Finalize();
}
