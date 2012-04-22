// -*- mode: c++ -*-
#ifndef _BLOCK_CYCLIC_MAT_H_
#define _BLOCK_CYCLIC_MAT_H_

#include <memory>
#include <vector>
#include "blacs.h"
#include "blacs_grid.h"

/// <summary>
///   A class that represents a two-dimensional block-cyclically distributed 
///   matrix in ScaLAPACK
/// </summary>
class block_cyclic_mat_t 
{
public:
    enum fill_t {ZERO, CONSTANT, DIAGONAL, RANDOM};

    /// <summary>
    ///   Constructs a new block-cyclically distributed matrix.
    /// </summary>
    /// <param name="grid">
    ///   The BLACS grid on which the matrix must be distributed.    
    /// </param>
    /// <param name="global_rows">
    ///   The global number of rows in the matrix, corresponds to the 
    ///   the parameter M_A in ScaLAPACK.
    /// </param>
    /// <param name="global_cols">
    ///   The global number of columns in the matrix, corresponds to the 
    ///   the parameter N_A in ScaLAPACK.
    /// </param>
    /// <param name="row_block_size">
    ///   The row block size, corresponds to the parameter MB_A in ScaLAPACK.
    /// </param>
    /// <param name="col_block_size">
    ///   The column block size, corresponds to the parameter NB_A in ScaLAPACK.
    /// </param>
    /// <param name="fill">
    ///   An enumeration describing how the elements of the matrix must be 
    ///   populated. It can take the following values:
    ///     ZERO: Fill the matrix with zeros (DEFAULT)
    ///     CONSTANT: Fill the matrix with a constant value.
    ///     DIAGONAL: Set the diagonals of the matrix to a constant value.
    ///     RANDOM: Fill the matrix with random values in (0,1) drawn from 
    ///         a uniform distribution.
    /// </param>
    /// <param name="alpha">
    ///   The constant value used for populating the elements or the diagonal
    ///   of the matrix, defaults to 0.
    ///   If fill is CONSTANT, then alpha is the value to which all elements
    ///   of the matrix must be set to.
    ///   If fill is DIAGONAL, then alpha represents the diagonal value.
    /// </param>
    block_cyclic_mat_t (std::shared_ptr<blacs_grid_t> grid, 
        blas_idx_t global_rows, blas_idx_t global_cols, 
        blas_idx_t row_block_size = s_block_size, blas_idx_t col_block_size = s_block_size,
        fill_t fill = ZERO, double alpha = 0.0);
    
    /// <summary>
    ///   Utility function for constructing a distributed matrix with random entries.
    /// </summary>
    static std::shared_ptr<block_cyclic_mat_t>  random   (std::shared_ptr<blacs_grid_t> grid, blas_idx_t global_rows, blas_idx_t global_cols);

    /// <summary>
    ///   Utility function for constructing a distributed matrix with a constant value.
    /// </summary>
    static std::shared_ptr<block_cyclic_mat_t>  constant (std::shared_ptr<blacs_grid_t> grid, blas_idx_t global_rows, blas_idx_t global_cols, double alpha = 0.0);
    
    /// <summary>
    ///   Utility function for constructing a distributed matrix with a constant diagonal value.
    /// </summary>
    static std::shared_ptr<block_cyclic_mat_t>  diagonal (std::shared_ptr<blacs_grid_t> grid, blas_idx_t global_rows, blas_idx_t global_cols, double alpha = 1.0);
    
    /// <summary>
    ///   Returns the total number of elements in the local part of the matrix
    ///   in the calling rank.
    /// </summary>
    blas_idx_t local_size() const;

    /// <summary>
    ///   Returns the number of rows in the local part of the matrix
    ///   in the calling rank, LOCr(M_A) in ScaLAPACK convention.
    /// </summary>
    
    blas_idx_t local_rows() const;

    /// <summary>
    ///   Returns the number of columns in the local part of the matrix
    ///   in the calling rank, LOCc(N_A) in ScaLAPACK convention.
    /// </summary>

    blas_idx_t local_cols() const;

    /// <summary>
    ///   Returns the row block size, MB_A.
    /// </summary>
    blas_idx_t row_block_size() const;
    
    /// <summary>
    ///   Returns the column block size NB_A.
    /// </summary>
    blas_idx_t col_block_size() const;

    /// <summary>
    ///   Returns the global number of rows of the matrix, M_A.
    /// </summary>
    blas_idx_t global_rows() const;

    /// <summary>
    ///   Returns the global number of columns of the matrix, N_A.
    /// </summary>
    blas_idx_t global_cols() const;

    /// <summary>
    ///   Returns the local data for the matrix in the calling rank.
    /// </summary>
    double* local_data();    

    /// <summary>
    ///   Returns the ScaLAPACK matrix descriptor, DESC_A.
    /// </summary>
    blas_idx_t* descriptor();

    /// <summary>
    ///   Returns the BLACS grid on this this matrix is distributed.
    /// </summary>
    std::shared_ptr<blacs_grid_t> grid();

    /// <summary>
    ///   Prints out the local portion of a distributed matrix.
    /// </summary>
    void print() const; 

private:
    std::vector<double> m_local_data;
    blas_idx_t   m_local_size;
    blas_idx_t   m_local_rows;
    blas_idx_t   m_local_cols;
    blas_idx_t   m_mb;
    blas_idx_t   m_nb;
    blas_idx_t   m_global_rows;
    blas_idx_t   m_global_cols;
    blas_idx_t   m_desc[DLEN_];
    std::shared_ptr<blacs_grid_t> m_grid;    

    block_cyclic_mat_t(const block_cyclic_mat_t&);
    block_cyclic_mat_t operator=(const block_cyclic_mat_t&);

    static const blas_idx_t s_block_size = 64;
};

#endif
