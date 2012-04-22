#include <mpi.h>
#include <cstdio>
#include <cassert>
#include "block_cyclic_mat.h"
#include "scalapack.h"

static double gesv_flops(blas_idx_t N, blas_idx_t NR)
{
    // From:
    // https://icl.cs.utk.edu/svn/scalapack-dev/scalapack/trunk/TESTING/LIN/pdludriver.f
    // Factorization: 2/3 N^3 - 1/2 N^2
    // Back solve   : NR * 2 N^2

    return ((2.0/3.0 * N * N * N) - (1.0/2.0 * N * N) + (NR * N * N))/1024.0/1024.0/1024.0;
}

static void lu_driver(blas_idx_t m_global, blas_idx_t n_global = 1)
{
    auto grid = std::make_shared<blacs_grid_t>();

    // Create a MxM random matrix A
    auto a = block_cyclic_mat_t::random(grid, m_global, m_global);    

    // Save the local data for A since it is overwritten during factorization
    std::vector<double> a_save(a->local_data(), a->local_data() + a->local_size());

    // Compute the 1-norm of A
    double norm_a = 0;
    {        
        blas_idx_t ia = 1, ja = 1;
        char norm = '1';
        std::vector<double> work(a->local_cols());
        norm_a = pdlange_(norm, m_global, m_global, a->local_data(), ia, ja, a->descriptor(),
            work.data());
    }

    // Create a MxN right-hand-side matrix filled with the value 42
    // This is overwritten with the solution of Ax = b
    auto x = block_cyclic_mat_t::constant(grid, m_global, n_global, 42.0);

    std::vector<blas_idx_t> ipiv(a->local_rows() + a->row_block_size() + 100);
    blas_idx_t ia = 1, ja = 1;
    blas_idx_t ib = 1, jb = 1;
    blas_idx_t info;

    MPI_Barrier (MPI_COMM_WORLD);
    
    // First compute Ax = b
    double t0 = MPI_Wtime();    
    pdgesv_ (m_global, n_global, 
        a->local_data(), ia, ja, a->descriptor(), 
        ipiv.data(), 
        x->local_data(), ib, jb, x->descriptor(), info);
    assert(info == 0);
    double t1 = MPI_Wtime() - t0;

    // Then form r = Ax - b
    // This is done by first letting r = b
    // followed by r = alpha * Ax + beta * r where alpha = 1.0, beta = -1.0    
    auto r = block_cyclic_mat_t::constant(grid, m_global, n_global, 42.0);
    char nein='N';
    double alpha = 1.0;
    double beta  = -1.0;
    pdgemm_(nein, nein, m_global, n_global, m_global, 
        alpha,        
        a_save.data()  , ia, ja, a->descriptor(),
        x->local_data(), ib, jb, x->descriptor(),
        beta,
        r->local_data(), ib, jb, r->descriptor());
        
    // Now compute Infinity norm of r
    double norm_r = -2.0;
    {
        char norm='I';
        std::vector<double> work(r->local_rows());
        norm_r = pdlange_(norm, 
            m_global, n_global, 
            r->local_data(), ib, jb, r->descriptor(), 
            work.data());
    }    
    
    // Compute the error
    // ||Ax - b||_oo/ (M x ||A||_1)

    double err = norm_r/m_global/norm_a;
    
    double t_glob;
    MPI_Reduce(&t1, &t_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (grid->iam() == 0) 
    {
        double gflops = gesv_flops(m_global, n_global)/t_glob/grid->nprocs();
        printf("\n"
            "MATRIX SOLVE BENCHMARK SUMMARY\n"
            "==============================\n"
            "N = %d\tNRHS = %d\tNP = %d\tNP_ROW = %d\tNP_COL = %d\n"
            "Time for PxGESV = %10.7f seconds\tGflops/Proc = %10.7f, Error = %f\n",
            m_global, n_global, grid->nprocs(), grid->nprows(), grid->npcols(), 
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

  lu_driver(n_global);
  MPI_Finalize();
}
