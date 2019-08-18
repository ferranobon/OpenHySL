#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#if _FLOAT_
#define hysl_float_t  float

#if _HDF5_
#define H5T_NATIVE_HYSL_FLOAT H5T_NATIVE_FLOAT
#endif

#if _MPI_
#define MPI_HYSL_FLOAT MPI_FLOAT
#endif

#define hysl_abs fabsf

#define hysl_copy scopy_
#define hysl_scal sscal_
#define hysl_axpy saxpy_
#define hysl_gemv sgemv_
#define hysl_symv ssymv_
#define hysl_spmv sspmv_

#define hysl_lacpy slacpy_
#define hysl_lascl slascl_
#define hysl_sygv  ssygv_
#define hysl_potrf spotrf_
#define hysl_potri spotri_
#define hysl_pptrf spptrf_
#define hysl_pptri spptri_

#define hysl_pcopy  pscopy_
#define hysl_pscal  psscal_
#define hysl_paxpy  psaxpy_
#define hysl_placpy pslacpy_
#define hysl_plascl pslascl_
#define hysl_psymv  pssymv_
#define hysl_pgemv  psgemv_
#define hysl_pgeadd psgeadd_
#define hysl_ptradd pstradd_

#define hysl_ppotrf pspotrf_
#define hysl_ppotri pspotri_

#if _MKL_
#define hysl_mkl_csrmv  mkl_scsrmv
#define hysl_mkl_dnscsr mkl_sdnscsr
#define hysl_mkl_csradd mkl_scsradd
#endif /* _MKL_ */


#else /* FLOAT */
#define hysl_float_t double

#if _HDF5_
#define H5T_NATIVE_HYSL_FLOAT H5T_NATIVE_DOUBLE
#endif

#if _MPI_
#define MPI_HYSL_FLOAT MPI_DOUBLE
#endif

#define hysl_abs fabs

#define hysl_copy dcopy_
#define hysl_scal dscal_
#define hysl_axpy daxpy_
#define hysl_gemv dgemv_
#define hysl_symv dsymv_
#define hysl_spmv dspmv_


#define hysl_lacpy dlacpy_
#define hysl_lascl dlascl_
#define hysl_sygv  dsygv_
#define hysl_potrf dpotrf_
#define hysl_potri dpotri_
#define hysl_pptrf dpptrf_
#define hysl_pptri dpptri_

#define hysl_pcopy  pdcopy_
#define hysl_pscal  pdscal_
#define hysl_paxpy  pdaxpy_
#define hysl_placpy pdlacpy_
#define hysl_plascl pdlascl_
#define hysl_psymv  pdsymv_
#define hysl_pgemv  pdgemv_
#define hysl_pgeadd pdgeadd_
#define hysl_ptradd pdtradd_

#define hysl_ppotrf pdpotrf_
#define hysl_ppotri pdpotri_

#if _MKL_
#define hysl_mkl_csrmv  mkl_dcsrmv
#define hysl_mkl_dnscsr mkl_ddnscsr
#define hysl_mkl_csradd mkl_dcsradd
#endif /* _MKL_ */

#endif /* _FLOAT_ */

#if _MKL_
#define HYSL_INT    MKL_INT


#define scopy_ scopy
#define sscal_ sscal
#define saxpy_ saxpy
#define sgemv_ sgemv
#define ssymv_ ssymv
#define sspmv_ sspmv

#define dcopy_ dcopy
#define dscal_ dscal
#define daxpy_ daxpy
#define dgemv_ dgemv
#define dsymv_ dsymv
#define dspmv_ dspmv


#else
#define HYSL_INT    int

#endif

#endif  /* DEFINITIONS_H_ */ 
