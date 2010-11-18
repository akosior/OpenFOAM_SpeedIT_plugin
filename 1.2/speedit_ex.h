/*
 * speedit_ex.h - part of the SpeedIT Extreme toolkit
 * Copyright 2010 (C) Vratis Ltd
 * email: support@vratis.com
 * 
 */

#ifndef __SPEEDIT_EX_H__
#define __SPEEDIT_EX_H__

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32

  #ifdef WINDLL
    //#error You should not compile this file as a part of dll library!
    #define DLLAPI __declspec(dllexport)
  #else
    #define DLLAPI __declspec(dllimport)
  #endif

#else
 #define DLLAPI
#endif

/*
  Error handling
*/

#define CUBLAS_ERR_BASE 10000

enum {
  ERR_OK                    = 0,

  ERR_GPU_MEM_ALLOC         = 20000,
  ERR_BAD_GPU_POINTER       = 20001,
  ERR_WRONG_COPY_DIRECTION  = 20002,
  ERR_DOUBLE_UNSUPPORTED    = 20003,
  ERR_BICGSTAB_FAILED       = 20004,
  ERR_OMEGA_VANISHED        = 20005,
  ERR_TOO_LITTLE_ITER       = 20006, 
  ERR_ZERO_ON_DIAGONAL      = 20007, 

  ERR_UNKNOWN               = 20999
} ;

DLLAPI
const char* si_errstr(int err_code) ;

typedef enum {
  P_NONE = 0,
  P_DIAG
} PRECOND_TYPE ;

/*
  Initialization
*/
DLLAPI
int si_init     (void) ;
DLLAPI
int si_shutdown (void) ;

/*
  Memory management functions.
*/

// Allocate GPU buffer
//
DLLAPI
int si_gsmalloc(int size, float**  out_ptr) ;
DLLAPI
int si_gdmalloc(int size, double** out_ptr) ;
DLLAPI
int si_gimalloc(int size, int**    out_ptr) ;
DLLAPI
int si_gvmalloc(int size, void**   out_ptr) ;

// Free GPU buffer
//
DLLAPI
int si_gsfree(float**  out_ptr) ;
DLLAPI
int si_gdfree(double** out_ptr) ;
DLLAPI
int si_gifree(int**    out_ptr) ;
DLLAPI
int si_gvfree(void**   out_ptr) ;

// Copy data from CPU to GPU memory
//
DLLAPI
int si_c2gvcopy(int size, const void*   in_ptr, void*   out_ptr) ;
DLLAPI
int si_c2gscopy(int size, const float*  in_ptr, float*  out_ptr) ;
DLLAPI
int si_c2gdcopy(int size, const double* in_ptr, double* out_ptr) ;
DLLAPI
int si_c2gicopy(int size, const int*    in_ptr, int*    out_ptr) ;

// Copy data from GPU to CPU memory
//
DLLAPI
int si_g2cvcopy(int size, const void*   in_ptr, void*   out_ptr) ;
DLLAPI
int si_g2cscopy(int size, const float*  in_ptr, float*  out_ptr) ;
DLLAPI
int si_g2cdcopy(int size, const double* in_ptr, double* out_ptr) ;
DLLAPI
int si_g2cicopy(int size, const int*    in_ptr, int*    out_ptr) ;

// Allocate buffer in GPU memory and copy data fro CPU memory
//
DLLAPI
int si_c2gvmcopy(int size, const void*   in_ptr, void**   out_ptr) ;
DLLAPI
int si_c2gsmcopy(int size, const float*  in_ptr, float**  out_ptr) ;
DLLAPI
int si_c2gdmcopy(int size, const double* in_ptr, double** out_ptr) ;
DLLAPI
int si_c2gimcopy(int size, const int*    in_ptr, int**    out_ptr) ;


/*
  Sparse BLAS Level 3 routines
*/

// All pointers have to be addresses of buffers in CPU memory

DLLAPI
int si_cscsrmv(int n_rows, const float*  vals, const int* c_idx, const int* r_idx, 
                           const float*  x,    float*  y) ;
DLLAPI
int si_cdcsrmv(int n_rows, const double* vals, const int* c_idx, const int* r_idx,
                           const double* x,    double* y) ;

// All pointers have to be addresses of buffers in GPU memory

DLLAPI
int si_gscsrmv(int n_rows, const float*  vals, const int* c_idx, const int* r_idx,
                           const float*  x,    float*  y) ;
DLLAPI
int si_gdcsrmv(int n_rows, const double* vals, const int* c_idx, const int* r_idx,
                           const double* x,    double* y) ;


/*
  Linear equation system solvers
*/

// All pointers have to be addresses of buffers in CPU memory

DLLAPI
int si_cscsrbicgstab(      int    n_rows, 
                     const float* vals, const int* c_idx, const int *r_idx, 
                           float* x, 
                     const float* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           float* eps) ;
DLLAPI
int si_cdcsrbicgstab(      int    n_rows, 
                     const double* vals, const int* c_idx, const int *r_idx, 
                           double* x, 
                     const double* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           double* eps) ;
DLLAPI
int si_cscsrcg(      int    n_rows, 
                     const float* vals, const int* c_idx, const int *r_idx, 
                           float* x, 
                     const float* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           float* eps) ;
DLLAPI
int si_cdcsrcg(      int    n_rows, 
                     const double* vals, const int* c_idx, const int *r_idx, 
                           double* x, 
                     const double* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           double* eps) ;

// All pointers except n_iter and eps have to be addresses of buffers in GPU 
// memory.

DLLAPI
int si_gscsrbicgstab(      int    n_rows, 
                     const float* vals, const int* c_idx, const int *r_idx, 
                           float* x, 
                     const float* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           float* eps) ;
DLLAPI
int si_gdcsrbicgstab(      int    n_rows, 
                     const double* vals, const int* c_idx, const int *r_idx, 
                           double* x, 
                     const double* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           double* eps) ;
DLLAPI
int si_gscsrcg(      int    n_rows, 
                     const float* vals, const int* c_idx, const int *r_idx, 
                           float* x, 
                     const float* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           float* eps) ;
DLLAPI
int si_gdcsrcg(      int    n_rows, 
                     const double* vals, const int* c_idx, const int *r_idx, 
                           double* x, 
                     const double* b, 
                     PRECOND_TYPE precond, 
                           int*   n_iter, 
                           double* eps) ;

#ifdef __cplusplus
}
#endif

#endif
