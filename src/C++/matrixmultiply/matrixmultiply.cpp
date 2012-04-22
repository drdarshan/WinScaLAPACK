#include <mpi.h>
#include "block_cyclic_mat.h"
#include "scalapack.h"

static double gemm_flops(blas_idx_t M, blas_idx_t N, blas_idx_t K)
{
    return (2.0 * M * N * K)/(1024.0 * 1024.0 * 1024.0);
}

static void dgemm_driver(blas_idx_t m_global, blas_idx_t n_global, blas_idx_t k_global)
{
    auto grid = std::make_shared<blacs_grid_t>();

    auto a = block_cyclic_mat_t::random(grid, m_global, k_global);
    auto b = block_cyclic_mat_t::random(grid, k_global, n_global);
    auto c = block_cyclic_mat_t::random(grid, m_global, n_global);

    MPI_Barrier(MPI_COMM_WORLD);

    double alpha = 1.0, beta = 0.0;

    double t0 = MPI_Wtime();
    char NEIN = 'N';
    blas_idx_t ia = 1, ja = 1, ib = 1, jb = 1, ic = 1, jc = 1;

    pdgemm_ (NEIN, NEIN, m_global, n_global, k_global, 
        alpha, 
        a->local_data(), ia, ja, a->descriptor(), 
        b->local_data(), ib, jb, b->descriptor(),
        beta,
        c->local_data(), ic, jc, c->descriptor());
    double t1 = MPI_Wtime() - t0;

    double t_glob;
    MPI_Reduce(&t1, &t_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

    if (grid->iam() == 0) 
    { 
        double gflops = gemm_flops(m_global, n_global, k_global)/t_glob/grid->nprocs();

        printf("\n"
            "MATRIX MULTIPLY BENCHMARK SUMMARY\n"
            "=================================\n"
            "M = %d\tN = %d\tK = %d\tNP = %d\tNP_ROW = %d\tNP_COL = %d\n"
            "Time for PxGEMM = %10.7f seconds\tGFlops/Proc = %10.7f\n", 
            m_global, n_global, k_global, grid->nprocs(), grid->nprows(), grid->npcols(),
            t_glob, gflops); fflush(stdout);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    blas_idx_t m_global = 4096;
    blas_idx_t n_global = 4096;
    blas_idx_t k_global = 4096;


    if (argc > 1)
    {
        m_global = blas_idx_t(atol(argv[1]));
    }

    if (argc > 2)
    {
        n_global = blas_idx_t(atol(argv[2]));
    }

    if (argc > 3)
    {
        k_global = blas_idx_t(atol(argv[3]));
    }
    
    dgemm_driver(m_global, n_global, k_global);
    MPI_Finalize();
}