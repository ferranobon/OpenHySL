/**
 * \file Input_Load.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 17th of March 2013
 * 
 * \todo Parallelise the loop in InputLoad_Generate_LoadVectorForm().
 * \todo Single precision routines.
 *
 * \brief Routines to deal with input loads.
 *
 * Routines for calculating the excitation that enters the structure in relative or absolute values. The
 * routines are available for general, packed and sparse storage and they make use of the BLAS library to
 * perform the linear algebra operations and they support both single and double precision. Sparse BLAS
 * operations are supported through the Intel MKL library.
 */

#ifndef INPUT_LOAD_H_
#define INPUT_LOAD_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

#include "Definitions.h"

/**
 * \brief Calculates the input load as absolute values. General storage version.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the packed storage or the
 * sparse version the routines InputLoad_AbsValues_PS() or InputLoad_AbsValues_Sp() should be used instead.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Stiff  The stiffness matrix \f$\mathcal K\f$.
 * \param[in]     Damp   The viscous damping matrix \f$\mathcal C\f$.
 * \param[in]     GDisp  Vector containing the ground displacement of the earthquake at a certain step \f$\vec
 *                       u_g\f$.
 * \param[in]     GVel   Vector containing the ground velocity of the earthquake at a certain step
 *                       \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the
 *                       size of the vector is referenced, not its elements.
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
 * \brief Calculates the input load as absolute values. Packed storage version.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the general storage or the
 * sparse version the routines InputLoad_AbsValues() or InputLoad_AbsValues_Sp() should be used instead.
 *
 * \pre
 * - \c Stiff and \c Damp should be symmetrical matrices in packed storage format with the upper triangular
 *   part referenced (lower part in FORTRAN). They should also be properly initialised through the
 *   MatrixVector_Create_PS() routine.
 * - The rest of the elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - The dimensions of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Stiff  The stiffness matrix \f$\mathcal K\f$.
 * \param[in]     Damp   The viscous damping matrix \f$\mathcal C\f$.
 * \param[in]     GDisp  Vector containing the ground displacement of the earthquake at a certain step \f$\vec
 *                       u_g\f$.
 * \param[in]     GVel   Vector containing the ground velocity of the earthquake at a certain step
 *                       \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering absolute values.
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * \sa MatrixVector_t.
 */
void InputLoad_AbsValues_PS( const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp,
			     const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			     MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as absolute values. Sparse version.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For the general or packed storage version the routines
 * InputLoad_AbsValues() or InputLoad_AbsValues_PS() should be used instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Stiff  The stiffness matrix \f$\mathcal K\f$.
 * \param[in]     Damp   The viscous damping matrix \f$\mathcal C\f$.
 * \param[in]     GDisp  Vector containing the ground displacement of the earthquake at a certain step \f$\vec
 *                       u_g\f$.
 * \param[in]     GVel   Vector containing the ground velocity of the earthquake at a certain step
 *                       \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering absolute values.
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void InputLoad_AbsValues_Sp( const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			     MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as absolute values. MPI version.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of PBLAS routines to perform the linear algebra operations. For the, general, packed storage
 * or  the sparse version the routines InputLoad_AbsValues(), InputLoad_AbsValues_PS() or
 * InputLoad_AbsValues_Sp() should be used instead.
 *
 * \pre 
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Stiff  The stiffness matrix \f$\mathcal K\f$.
 * \param[in]     Damp   The viscous damping matrix \f$\mathcal C\f$.
 * \param[in]     GDisp  Vector containing the ground displacement of the earthquake at a certain step \f$\vec
 *                       u_g\f$.
 * \param[in]     GVel   Vector containing the ground velocity of the earthquake at a certain step
 *                       \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering absolute values.
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * \sa PMatrixVector_t.
 */
void InputLoad_AbsValues_MPI( PMatrixVector_t *const Stiff, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const GDisp, PMatrixVector_t *const GVel,
			      PMatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as relative value. General storage version.
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the packed storage or the
 * sparse version the routines InputLoad_RelValues_PS() or InputLoad_RelValues_Sp() should be used instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the vectors and the matrix must be coherent since it will not be checked in the
 *   routine.
 * - The ground acceleration vector should have already the right values on the right position since the
 *   routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Mass   The mass matrix \f$\mathcal M\f$.
 * \param[in]     GAcc   Vector containing the ground acceleration of the earthquake at a certain step
 *                       \f$\ddot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the
 *                       size of the vector is referenced, not its elements.
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
 * \brief Calculates the input load as relative value. Packed storage version.
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For the general storage or the
 * sparse version the routines InputLoad_RelValues() or InputLoad_RelValues_Sp() should be used instead.
 *
 * \pre
 * - \c Mass should be a symmetrical matrix in packed storage format with the upper triangular part referenced
 *   (lower part in FORTRAN). It should also be properly initialised through the MatrixVector_Create_PS()
 *   routine.
 * - The rest of the elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - The dimensions of the vectors and the matrix must be coherent since it will not be checked in the
 *   routine.
 * - The ground acceleration vector should have already the right values on the right position since the
 *   routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Mass   The mass matrix \f$\mathcal M\f$.
 * \param[in]     GAcc   Vector containing the ground acceleration of the earthquake at a certain step
 *                       \f$\ddot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering relative values.
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * \sa MatrixVector_t.
 */
void InputLoad_RelValues_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc,
			     MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as relative value. Sparse version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For the general or packed storage version the routines
 * InputLoad_RelValues() or InputLoad_RelValues_PS() should be used instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - \c Mass must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - \c Mass must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - \c Mass must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground acceleration vector should have already the right values on the right position since the
 *   routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Mass   The mass matrix \f$\mathcal M\f$.
 * \param[in]     GAcc   Vector containing the ground acceleration of the earthquake at a certain step \f$\ddot{\vec
 *                       u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering relative values.
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void InputLoad_RelValues_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_t *const GAcc, 
			     MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as relative value. MPI version.
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step
 * through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of PBLAS routines to perform the linear algebra operations. For the general, packed or sparse
 * formats, the routines InputLoad_RelValues(), InputLoad_RelValues_PS() or InputLoad_RelValues_Sp() should be
 * used instead.
 *
 * \pre
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the vectors and the matrix must be coherent since it will not be checked in the
 *   routine.
 * - The ground acceleration vector should have already the right values on the right position since the
 *   routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in]     Mass   The mass matrix \f$\mathcal M\f$.
 * \param[in]     GAcc   Vector containing the ground acceleration of the earthquake at a certain step
 *                       \f$\ddot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the
 *                       size of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering relative values.
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * \sa PMatrixVector_t.
 */
void InputLoad_RelValues_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const GAcc,
			      PMatrixVector_t *const InLoad );

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
 * \param[in]     DOF            The pattern to be applied.
 * \param[in,out] LoadVectorForm The load vector form. On input, only the dimensions are referenced.
 *
 * \post \c LoadVectorForm contains the specified pattern \f$N = \frac{Order}{n}\f$ times, where \f$Order\f$
 *       is the length of this vector and \f$n\f$ the length of the load pattern.
 *
 * \sa MatrixVector_t and Algorithm_GetExcitedDOF().
 */
void InputLoad_Generate_LoadVectorForm( int *DOF, MatrixVector_t *const LoadVectorForm );

/**
 * \brief Generates an input load vector form from a given pattern. MPI version.
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
 * It makes use of ScaLAPACK auxiliary routines in order to indentify the local indeces.
 *
 * \pre
 * - \c LoadVector must be properly initialised through PMatrixVector_Create().
 * - The first element of the pattern must be indicative of its length.
 * - The length of \c LoadVector must be a multiple of the length of the specified pattern.
 *
 * \param[in]     DOF            The pattern to be applied.
 * \param[in,out] LoadVectorForm The load vector form. On input, only the dimensions are referenced.
 *
 * \post \c LoadVectorForm contains the specified pattern \f$N = \frac{Order}{n}\f$ times, where \f$Order\f$
 *       is the length of this vector and \f$n\f$ the length of the load pattern.
 *
 * \sa PMatrixVector_t and Algorithm_GetExcitedDOF().
 */
void InputLoad_Generate_LoadVectorForm_MPI( int *DOF, PMatrixVector_t *const LoadVectorForm );

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
 *                 1.0 1.0 1.0 0.0 0.0 -2.0
 *                             ....
 *                 1.0 1.0 1.0 0.0 0.0 -2.0}</tt>
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
 *   routine.
 * - The dimensions of the vectors  must be the identical.
 * - In \c LoadForm the equivalent elements of the degrees of freedom to be excited must have a value
 *   different than 0.
 *
 * \param[in]  LoadForm   A load vector form.
 * \param[in]  Value      The desired ground measurement to be applied.
 * \param[out] LoadVector The load vector.
 *
 * \post \c LoadVector is the result of multiplying the \c LoadForm vector by the given ground measurement. It
 *       is ready to be used in calculating the input load as relative (InputLoad_RelValues()) or absolute
 *       values (InputLoad_AbsValues()).
 *
 * \sa MatrixVector_t and InputLoad_Generate_LoadVectorForm().
 */
void InputLoad_Apply_LoadVectorForm( const MatrixVector_t *const LoadForm, const HYSL_FLOAT Value,
				     MatrixVector_t *const LoadVector );

/**
 * \brief Given load vector form and a ground measurement, a excitation vector is generated. MPI version.
 *
 * This routine applies a ground motion measurement (displacement, velocity or acceleration) into a load
 * vector form. The output of this routine will be a vector with the specified measurement value applied to
 * those parts of the load vector form that are different to 0.0. In this way, only the desired degrees of
 * freedom (or system equations) will have an be external load applied. The output is ready to be used in
 * calculating the input load as relative (InputLoad_RelValues_MPI()) or absolute values
 * (InputLoad_AbsValues_MPI()). Therefore if the following \c LoadForm is considered:
 *
 * <tt>LoadForm = {1.0 1.0 1.0 0.0 0.0 -2.0
 *                 1.0 1.0 1.0 0.0 0.0 -2.0
 *                             ....
 *                 1.0 1.0 1.0 0.0 0.0 -2.0}</tt>
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
 *   routine.
 * - The dimensions of the vectors  must be the identical.
 * - In \c LoadForm the equivalent elements of the degrees of freedom to be excited must have a value
 *   different than 0.
 *
 * \param[in]  LoadForm   A load vector form.
 * \param[in]  Value      The desired ground measurement to be applied.
 * \param[out] LoadVector The load vector.
 *
 * \post \c LoadVector is the result of multiplying the \c LoadForm vector by the given ground measurement. It
 *       is ready to be used in calculating the input load as relative (InputLoad_RelValues_MPI()) or absolute
 *       values (InputLoad_AbsValues_MPI()).
 *
 * \sa PMatrixVector_t and InputLoad_Generate_LoadVectorForm_MPI().
 */
void InputLoad_Apply_LoadVectorForm_MPI( PMatrixVector_t *const LoadForm, const HYSL_FLOAT Value,
					 PMatrixVector_t *const LoadVector );

#endif /* INPUT_LOAD_H_ */
