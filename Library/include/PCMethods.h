/**
 * \file PCMethods.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 8th of September 2013
 * 
 * \todo Single precision routines.
 *
 * \brief Routines to deal the correction of the acceleration in Predictor-Corrector methods.
 *
 * Routines for calculating correction of the acceleration in Predictor-Corrector methods. The
 * routines are available for general, packed and sparse storage and they make use of the BLAS library to
 * perform the linear algebra operations and they support both single and double precision. Sparse BLAS
 * operations are supported through the Intel MKL library.
 */

#ifndef PCMETHODS_H_
#define PCMETHODS_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

void PC_PredictorStep_Displacement ( const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
				     const MatrixVector_t *const AccT, const hysl_float_t a9, const hysl_float_t a10,
				     MatrixVector_t *const DispTdT_Pred );
void PC_PredictorStep_Velocity ( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 const hysl_float_t a6, MatrixVector_t *const VelTdT_Pred );
void PC_CorrectorStep_Displacement ( const MatrixVector_t *const DispTdT_Pred, const MatrixVector_t *const AccTdT,
				     const hysl_float_t a8, MatrixVector_t *const DispTdT);
void PC_CorrectorStep_Velocity ( const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const AccTdT,
				 const hysl_float_t a7, MatrixVector_t *const VelTdT);

void PC_ReactionForces_Numerical (const MatrixVector_t *const DispTdT_Pred, MatrixVector_t *const K,
				  MatrixVector_t *const RForceTdT);

void PC_ReactionForces_Numerical_PS (const MatrixVector_t *const DispTdT_Pred, MatrixVector_t *const K,
				     MatrixVector_t *const RForceTdT);

void PC_Calculate_Acceleration (const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT,
				const MatrixVector_t *const RForceTdT, const MatrixVector_t *const RForceT,
				const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const VelT_Pred,
				const MatrixVector_t *const AccT, const MatrixVector_t *const K,
				const MatrixVector_t *const C, const MatrixVector_t *const Meinv,
				const hysl_float_t alpha_H, const hysl_float_t a7, const hysl_float_t a8,
				MatrixVector_t *const AccTdT );

void PC_Calculate_Acceleration_PS (const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT,
				   const MatrixVector_t *const RForceTdT, const MatrixVector_t *const RForceT,
				   const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const VelT_Pred,
				   const MatrixVector_t *const AccT, const MatrixVector_t *const K,
				   const MatrixVector_t *const C, const MatrixVector_t *const Meinv,
				   const hysl_float_t alpha_H, const hysl_float_t a7, const hysl_float_t a8,
				   MatrixVector_t *const AccTdT );

/**
 * \brief Corrects the acceleration for Predictor-Corrector methods. General storage version.
 *void PC_CorrectorStep_Velocity ( const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const AccTdT,
				 const HYSL_FLOAT a7, MatrixVector_t *const VelTdT)
{
 * This routine calculates the corrected acceleration required by predictor-corrector methods through:
 *
 * \f[\ddot{\vec u}_{corr}= \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * where:
 * - \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$ is corrected acceleration at time \f$t+\Delta t\f$,
 * - \f$\mathcal{M}^{-1}\f$ is inverse of the mass matrix,
 * - \f$\vec l_i^{t + \Delta t}\f$ is the input load vector at time \f$t+\Delta t\f$,
 * - \f$\vec f_c= \vec f_r + \vec f_s\f$ is the vector with the simulated (\f$\vec f_s\f$) and measured
 *   (\f$\vec f_r\f$) non linear restoring forces,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the packed storage or the
 * sparse version the routines PC_Correct_Acceleration_PS() or PC_Correct_Acceleration_Sp() should be used instead.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The matrix must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 *
 * \param[in]     MInv        The inverse of the mass matrix \f$\mathcal M^{-1}\f$.
 * \param[in]     In_LoadT    The input load vector \f$l_i^{t+\Delta t}\f$.
 * \param[in]     fc          Vector containing the non-linear restoring forces \f$\vec f_c= \vec f_r + \vec
 *                            f_s\f$.
 * \param[in,out] Tempvec     Auxiliary vector. As an input, only the size of the vector is referenced, not
 *                            its elements.
 * \param[out]    AccTdT_Corr Corrected acceleration vector \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$.
 *
 * \post \c AccTdT_Corr is the corrected acceleration required by predictor-corrector methods:
 *
 * \f[\ddot{\vec u_c} = \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * \sa MatrixVector_t.
 */
void PC_Correct_Acceleration( const MatrixVector_t *const MInv, const MatrixVector_t *const In_LoadT,
			      const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
			      MatrixVector_t *const AccTdT_Corr );


/**
 * \brief Corrects the acceleration for Predictor-Corrector methods. Packed storage version.
 *
 * This routine calculates the corrected acceleration required by predictor-corrector methods through:
 *
 * \f[\ddot{\vec u}_{corr}= \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * where:
 * - \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$ is corrected acceleration at time \f$t+\Delta t\f$,
 * - \f$\mathcal{M}^{-1}\f$ is inverse of the mass matrix,
 * - \f$\vec l_i^{t + \Delta t}\f$ is the input load vector at time \f$t+\Delta t\f$,
 * - \f$\vec f_c= \vec f_r + \vec f_s\f$ is the vector with the simulated (\f$\vec f_s\f$) and measured
 *   (\f$\vec f_r\f$) non linear restoring forces,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the general storage or the
 * sparse version the routines PC_Correct_Acceleration() or PC_Correct_Acceleration_Sp() should be used
 * instead.
 *
 * \pre 
 * - \c MInv must be properly initialised through the MatrixVector_Create_PS() routine.
 * - The rest of the elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - The matrix must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 * routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 *
 * \param[in]     MInv        The inverse of the mass matrix \f$\mathcal M^{-1}\f$ in packed storage.
 * \param[in]     In_LoadT    The input load vector \f$l_i^{t+\Delta t}\f$.
 * \param[in]     fc          Vector containing the non-linear restoring forces \f$\vec f_c= \vec f_r + \vec
 *                            f_s\f$.
 * \param[in,out] Tempvec     Auxiliary vector. As an input, only the size of the vector is referenced, not
 *                            its elements.
 * \param[out]    AccTdT_Corr Corrected acceleration vector \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$.
 *
 * \post \c AccTdT_Corr is the corrected acceleration required by predictor-corrector methods:
 *
 * \f[\ddot{\vec u_c} = \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * \sa MatrixVector_t.
 */
void PC_Correct_Acceleration_PS( const MatrixVector_t *const MInv, const MatrixVector_t *const In_LoadT,
				 const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
				 MatrixVector_t *const AccTdT_Corr );

/**
 * \brief Corrects the acceleration for Predictor-Corrector methods. Sparse version.
 *
 * This routine calculates the corrected acceleration required by predictor-corrector methods through:
 *
 * \f[\ddot{\vec u}_{corr}= \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * where:
 * - \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$ is corrected acceleration at time \f$t+\Delta t\f$,
 * - \f$\mathcal{M}^{-1}\f$ is inverse of the mass matrix,
 * - \f$\vec l_i^{t + \Delta t}\f$ is the input load vector at time \f$t+\Delta t\f$,
 * - \f$\vec f_c= \vec f_r + \vec f_s\f$ is the vector with the simulated (\f$\vec f_s\f$) and measured
 *   (\f$\vec f_r\f$) non linear restoring forces,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For the general or packed storage version the routines
 * PC_Correct_Acceleration() or PC_Correct_Acceleration_PS() should be used instead.
 *
 * \pre 
 * - \c MInv must be properly initialised through the MatrixVector_Create_Sp() routine.
 * - \c MInv must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - \c MInv must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 *
 * \param[in]     MInv        The inverse of the mass matrix \f$\mathcal M^{-1}\f$ in packed storage.
 * \param[in]     In_LoadT    The input load vector \f$l_i^{t+\Delta t}\f$.
 * \param[in]     fc          Vector containing the non-linear restoring forces \f$\vec f_c= \vec f_r + \vec
 *                            f_s\f$.
 * \param[in,out] Tempvec     Auxiliary vector. As an input, only the size of the vector is referenced, not
 *                            its elements.
 * \param[out]    AccTdT_Corr Corrected acceleration vector \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$.
 *
 * \post \c AccTdT_Corr is the corrected acceleration required by predictor-corrector methods:
 *
 * \f[\ddot{\vec u_c} = \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void PC_Correct_Acceleration_Sp( const MatrixVector_Sp_t *const MInv, const MatrixVector_t *const In_LoadT,
				 const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
				 MatrixVector_t *const AccTdT_Corr );

/**
 * \brief Corrects the acceleration for Predictor-Corrector methods. Sparse version.
 *
 * This routine calculates the corrected acceleration required by predictor-corrector methods through:
 *
 * \f[\ddot{\vec u}_{corr}= \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * where:
 * - \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$ is corrected acceleration at time \f$t+\Delta t\f$,
 * - \f$\mathcal{M}^{-1}\f$ is inverse of the mass matrix,
 * - \f$\vec l_i^{t + \Delta t}\f$ is the input load vector at time \f$t+\Delta t\f$,
 * - \f$\vec f_c= \vec f_r + \vec f_s\f$ is the vector with the simulated (\f$\vec f_s\f$) and measured
 *   (\f$\vec f_r\f$) non linear restoring forces,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For the general, packed or sparse formats, the routines
 * PC_Correct_Acceleration(), PC_Correct_Acceleration_PS() or PC_Correct_Acceleration_Sp() should be used
 * instead.
 *
 * \pre 
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The matrix \c MInv must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 *
 * \param[in]     MInv        The inverse of the mass matrix \f$\mathcal M^{-1}\f$ in packed storage.
 * \param[in]     In_LoadT    The input load vector \f$l_i^{t+\Delta t}\f$.
 * \param[in]     fc          Vector containing the non-linear restoring forces \f$\vec f_c= \vec f_r + \vec
 *                            f_s\f$.
 * \param[in,out] Tempvec     Auxiliary vector. As an input, only the size of the vector is referenced, not
 *                            its elements.
 * \param[out]    AccTdT_Corr Corrected acceleration vector \f$\ddot{\vec u}_{cor}^{t + \Delta t}\f$.
 *
 * \post \c AccTdT_Corr is the corrected acceleration required by predictor-corrector methods:
 *
 * \f[\ddot{\vec u_c} = \mathcal{M}^{-1}(\vec l_i + \vec f_c)\f]
 *
 * \sa PMatrixVector_t.
 */
void PC_Correct_Acceleration_MPI( const PMatrixVector_t *const MInv, const PMatrixVector_t *const In_LoadT,
				  const PMatrixVector_t *const fc, PMatrixVector_t *const Tempvec,
				  PMatrixVector_t *const AccTdT_Corr );

#endif
