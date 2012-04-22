#include <cmath>
#include <cassert>

#include "blacs.h"
#include "blacs_grid.h"


blacs_grid_t::blacs_grid_t()
{    
    blacs_pinfo_ (m_iam, m_nprocs);

    blas_idx_t negone = -1, zero = 0, one = 1;    
    blacs_get_ (negone, zero, m_ictxt);
    
    m_nprows = blas_idx_t(sqrt(double(m_nprocs)));

    while(m_nprocs % m_nprows)
        m_nprows --;

    m_npcols = m_nprocs/m_nprows;

    assert(m_nprows * m_npcols == m_nprocs);

    const char* row_major = "Row";    

    blacs_gridinit_ (m_ictxt, row_major, m_nprows, m_npcols);    
    blacs_gridinfo_ (m_ictxt, m_nprows, m_npcols, m_myrow, m_mycol);
}

blas_idx_t blacs_grid_t::local_rows(blas_idx_t global_rows, blas_idx_t row_block_size, blas_idx_t row_offset /*= 0*/)
{
    return numroc_ (global_rows, row_block_size, m_myrow, row_offset, m_nprows);    
}

blas_idx_t blacs_grid_t::local_cols(blas_idx_t global_cols, blas_idx_t col_block_size, blas_idx_t col_offset /*= 0*/)
{
    return numroc_ (global_cols, col_block_size, m_mycol, col_offset, m_npcols); 
}

blas_idx_t blacs_grid_t::nprows() const
{
    return m_nprows;
}

blas_idx_t blacs_grid_t::npcols() const
{
    return m_npcols;
}

blas_idx_t blacs_grid_t::myprow() const
{
    return m_myrow;
}

blas_idx_t blacs_grid_t::mypcol() const
{
    return m_mycol;
}

blacs_grid_t::~blacs_grid_t()
{
    blacs_gridexit_(m_ictxt);
}

blas_idx_t blacs_grid_t::iam() const
{
    return m_iam;
}

blas_idx_t blacs_grid_t::nprocs() const
{
    return m_nprocs;
}

blas_idx_t blacs_grid_t::context() const
{
    return m_ictxt;
}
