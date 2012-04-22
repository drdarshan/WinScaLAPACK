#ifndef _BLACS_GRID_H_
#define _BLACS_GRID_H_
#include "index.h"

/// <summary>
///   A class representing a two-dimensional BLACS
///   process grid.
/// </summary>
class blacs_grid_t 
{
public:
    /// <summary>
    ///   Creates a two-dimensional process grid such the
    ///   number of process rows is roughly equal to the 
    ///   number of process columns.
    /// </summary>
    /// <remark>
    ///   The constructor calls the BLACS_GRIDINIT subroutine
    ///   to create a two-dimensional process grid with processes
    ///   numbered in row-major order and then calls BLACS_GRIDINFO
    ///   to populate information such as the number of process rows
    ///   and process columns
    /// </remark>
    blacs_grid_t();
    
    /// <summary>
    ///   Returns the number of rows in the process grid.
    /// </summary>
    blas_idx_t nprows() const;

    /// <summary>
    ///   Returns the number of rows in the process grid.
    /// </summary>
    blas_idx_t npcols() const;

    /// <summary>
    ///   Returns the total number of processes in the 
    ///   process grid, which is the number of process rows
    ///   times the number of process columns.
    /// </summary>
    blas_idx_t nprocs() const;

    /// <summary>
    ///   Returns the rank of the calling process.
    /// </summary>
    blas_idx_t iam() const;

    /// <summary>
    ///   Returns the process row of the calling process.
    /// </summary>
    blas_idx_t myprow() const;

    /// <summary>
    ///   Returns the process column of the calling process.
    /// </summary>
    blas_idx_t mypcol() const;

    /// <summary>
    ///   Returns the underlying BLACS context object.
    /// </summary>
    blas_idx_t context() const;

    /// <summary>
    ///   Given a global number of rows and a row block size
    ///   returns the local number of rows in the calling process.
    /// </summary>
    /// <param name="global_rows">
    ///   The global number of rows in the matrix.
    /// </param>
    /// <param name="row_block_size">
    ///   The blocking factor for the rows. This corresponds
    ///   to the factor MB in ScaLAPACK.
    /// </param>    
    /// <param name="row_offset">
    ///   The starting global row offset of the matrix. This
    ///   corresponds to the factor IA in ScaLAPACK and defaults
    ///   to zero.
    /// </param>
    /// <remark>
    ///   This method internally calls the NUMROC subroutine in BLACS.
    /// </remark>
    blas_idx_t local_rows(blas_idx_t global_rows, blas_idx_t row_block_size, blas_idx_t row_offset = 0);

    /// <summary>
    ///   Given a global number of columns and a column block size
    ///   returns the local number of columns in the calling process.
    /// </summary>
    /// <param name="global_cols">
    ///   The global number of columns in the matrix.
    /// </param>
    /// <param name="col_block_size">
    ///   The blocking factor for the columns. This corresponds
    ///   to the factor NB in ScaLAPACK.
    /// </param>    
    /// <param name="col_offset">
    ///   The starting global column offset of the matrix. This
    ///   corresponds to the factor JA in ScaLAPACK and defaults
    ///   to zero.
    /// </param>
    /// <remark>
    ///   This method internally calls the NUMROC subroutine in BLACS.
    /// </remark>
    blas_idx_t local_cols(blas_idx_t global_cols, blas_idx_t col_block_size, blas_idx_t col_offset = 0);
    
    virtual ~blacs_grid_t();

private:
    blas_idx_t m_ictxt;
    blas_idx_t m_iam;
    blas_idx_t m_nprocs;
    blas_idx_t m_nprows;
    blas_idx_t m_npcols;
    blas_idx_t m_myrow;
    blas_idx_t m_mycol;

    // Mark this class as non-copyable
    blacs_grid_t(const blacs_grid_t&);
    const blacs_grid_t& operator=(const blacs_grid_t&);
};

#endif // _BLACS_GRID_H_
