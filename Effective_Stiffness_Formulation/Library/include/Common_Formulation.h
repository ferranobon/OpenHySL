/**
 * \file Common_Formulation.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 9th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 * \todo Parallelise the loop in InputLoad_Generate_LoadVectorForm().
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines
 * make use of the BLAS library to perform the linear algebra operations and they support both single and
 * double precision. Sparse BLAS operations are supported through the Intel MKL library.
 */

#ifndef COMMON_FORMULATION_H_
#define COMMON_FORMULATION_H_

#include "MatrixVector.h"

/**
 * \brief Calculates the new state.
 *
 * The explicit part of the new state is computed. This routine does not depend on a specific formulation
 * and therefore it can be used in displacement, velocity and acceleration control formulations.
 * 
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t - \vec f_{err}^t)\f]
 *
 * where:
 * - \f$\vec n_0^{t + \Delta t}\f$ is the explicit part of the new state vector at time \f$t + \Delta t\f$,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - and \f$\vec f_{err}^t\f$ is the error compensation force at time \f$t\f$.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations.
 * 
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - \c Gain must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower part in FORTRAN).
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * 
 * \param[in] IGain the inverted gain matrix \f$\mathcal{G}^{1}\f$.
 * \param[in] Eff_ForceT the effective force vector \f$\vec f_{eff}^t\f$.
 * \param[in] In_LoadT the input load vector \f$\vec l_i^t\f$.
 * \param[in] Err_ForceT the error force vector \f$\vec f_{err}^t\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid
 * allocating and deallocating memory each step. On input only the number of rows is used.
 * \param[out] VecTdT_0 is new explicit state (displacement, velocity or acceleration depending on the selected formulation)
 * \f$\vec n_0^{t + \Delta t}\f$.
 *
 * \post
 * - \c VecTdT_0 is the result of the operation:
 *
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t + \vec f_{err}^t)\f]
 *
 * \sa MatrixVector_t and EffK_Calc_Effective_Force().
 */
void Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT, const MatrixVector_t *const In_LoadT,
		       const MatrixVector_t *const Err_ForceT, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

/**
 * \brief Calculates the input load as absolute values.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in] Stiff The stiffness matrix \f$\mathcal K\f$.
 * \param[in] Damp The viscous damping matrix \f$\mathcal C\f$.
 * \param[in] GDisp Vector containing the ground displacement of the earthquake at a certain step \f$\vec u_g\f$.
 * \param[in] GVel Vector containing the ground velocity of the earthquake at a certain step \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the size
 * of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering absolute values.
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * \sa MatrixVector_t.
 */
void InputLoad_AbsValues( const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp,
			  const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			  MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as relative value.
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground acceleration vector should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in] Mass The mass matrix \f$\mathcal M\f$.
 * \param[in] GAcc Vector containing the ground acceleration of the earthquake at a certain step \f$\ddot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the size
 * of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering relative values.
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * \sa MatrixVector_t.
 */
void InputLoad_RelValues( const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc,
			  MatrixVector_t *const InLoad );

/**
 * \brief Generates an input load vector form from a given pattern.
 *
 * Through a given pattern, an input load vector form is generated. This load vector form indicates which
 * equations are going to have an external force applied and it is generated given a pattern that repeats
 * itself until the end of the vector is reached. For example given the pattern:
 *
 * <tt>6 1 1 1 0 0 -2</tt>
 *
 * where the first number indicates its length, the resulting output vector form would be:
 *
 * <tt>LoadVector = {1.0 1.0 1.0 0.0 0.0 -2.0
 *                   1.0 1.0 1.0 0.0 0.0 -2.0
 *                             ....
 *                   1.0 1.0 1.0 0.0 0.0 -2.0}</tt>
 *
 * \pre
 * - \c LoadVector must be properly initialised through MatrixVector_Create().
 * - The first element of the pattern must be indicative of its length.
 * - The length of \c LoadVector must be a multiple of the length of the specified pattern.
 *
 * \param[in] DOF The pattern to be applied.
 * \param[in,out] LoadVectorForm The load vector form. On input, only the dimensions are referenced.
 *
 * \post \c LoadVectorForm contains the specified pattern \f$N = \frac{Order}{n}\f$ times, where \f$Order\f$
 * is the length of this vector and \f$n\f$ the length of the load pattern.
 *
 * \sa MatrixVector_t and Algorithm_GetExcitedDOF().
 */
void InputLoad_Generate_LoadVectorForm( int *DOF, MatrixVector_t *const LoadVectorForm );


/**
 * \brief Given load vector form and a ground measurement, a excitation vector is generated.
 *
 * This routine applies a ground motion measurement (displacement, velocity or acceleration) into a load
 * vector form. The output of this routine will be a vector with the specified measurement value applied to
 * those parts of the load vector form that are different to 0.0. In this way, only the desired degrees of
 * freedom (or system equations) will have an be external load applied. The output is ready to be used in
 * calculating the input load as relative (InputLoad_RelValues()) or absolute values
 * (InputLoad_AbsValues()). Therefore if the following \c LoadForm is considered:
 *
 * <tt>LoadForm = {1.0 1.0 1.0 0.0 0.0 -2.0
 *                   1.0 1.0 1.0 0.0 0.0 -2.0
 *                             ....
 *                   1.0 1.0 1.0 0.0 0.0 -2.0}</tt>
 *
 * the resulting \c LoadVector if the ground measurement is \c A would be:
 *
 * <tt>LoadVector = {A A A 0.0 0.0 -2.0*A
 *                   A A A 0.0 0.0 -2.0*A
 *                             ....
 *                   A A A 0.0 0.0 -2.0*A}</tt>
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 * routine.
 * - The dimensions of the vectors  must be the identical.
 * - In \c LoadForm the equivalent elements of the degrees of freedom to be excited must have a value
 * different than 0.
 *
 * \param[in] LoadForm A load vector form.
 * \param[in] Value The desired ground measurement to be applied.
 * \param[out] LoadVector The load vector.
 *
 * \post \c LoadVector is the result of multiplying the \c LoadForm vector by the given ground measurement. It
 * is ready to be used in calculating the input load as relative (InputLoad_RelValues()) or absolute values
 * (InputLoad_AbsValues())
 *
 * \sa MatrixVector_t and InputLoad_Generate_LoadVectorForm().
 */
void InputLoad_Apply_LoadVectorForm ( const MatrixVector_t *const LoadForm, const double Value, MatrixVector_t *const LoadVector );

#endif /* COMMON_FORMULATION_H_ */
