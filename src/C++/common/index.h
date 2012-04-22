// -*- mode: c++ -*-
#ifndef _INDEX_H_
#define _INDEX_H_

#include <cstdint>
#ifndef BLASINDEX64
typedef int32_t blas_idx_t;
#else
typedef int64_t blas_idx_t;
#endif

#endif // _INDEX_H_