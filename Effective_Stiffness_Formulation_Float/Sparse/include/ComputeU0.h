/**
 * \file ComputeU0.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the functions used during the Compute U0 phase.
 *
 * This file contains the prototypes of the functions that are called during the Compute U0 phase of the substructure
 * algorithm. This includes a function to calculate the input load and functions to decouple vectors into their non-coupling
 * and coupling part.
 *
 */

#ifndef COMPUTEU0_H_
#define COMPUTEU0_H_

#include "MatrixVector.h"
#include "Initiation.h"

/**
 * \brief Calculates the effective force vector \f$\overrightarrow{f_{eff}}\f$ according to the formulation using the effective stiffness matrix.
 * 
 * The effective force vector is calculated according to the formulation using the effective stiffness matrix in page 53 (\cite Dorka_1998). It is used to compute the implicit
 * part of the displacement vector.
 *
 * \f[\overrightarrow{f_{eff}} = \mathcal{M}(a_0\vec{u} + a_2\vec{\dot{u}} + a_3\vec{\ddot{u}}) + \mathcal{C}(a_1\vec{u} + a_4\vec{\dot{u}} + a_5\vec{\ddot{u}})\f]
 *
 * This routine assumes that the matrices are symmetrical and in general storage. Only the upper part will be referenced (lower part in FORTRAN).
 *
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - The size of the vectors and matrices must be coeherent, since it will not be checked in the routine.
 *
 * \param[in] Mass is the Mass matrix \f$\mathcal{M}\f$.
 * \param[in] Damp is the Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in] Disp The displacement vector \f$\vec{u}\f$.
 * \param[in] Vel is the velocity vector \f$\dot{\vec{u}}\f$.
 * \param[in] Acc is the acceleration vector \f$\ddot{\vec{u}}\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid allocating and deallocating memory each step. On input
 * only the number of rows is used
 * \param[in] a0 A constant of value \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$.
 * \param[in] a1 A constant of value \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$.
 * \param[in] a2 A constant of value \f$a_2 = \frac{1}{\beta_N\Delta t}\f$.
 * \param[in] a3 A constant of value \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$.
 * \param[in] a4 A constant of value \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$.
 * \param[in] a5 A constant of value \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$.
 * \param[out] Eff_Force is the effective force vector \f$\overrightarrow{f_{eff}}\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation \f$\overrightarrow{f_{eff}} = \mathcal{M}(a_0\vec{u} + a_2\vec{\dot{u}} + a_3\vec{\ddot{u}}) + \mathcal{C}(a_1\vec{u} + a_4\vec{\dot{u}} + a_5\vec{\ddot{u}})\f$.
 *
 * \sa MatrixVector and AlgConst.
 */
void EffK_Calc_Effective_Force( const MatrixVector *const Mass, const MatrixVector *const Damp,
				const MatrixVector *const Disp, const MatrixVector *const Vel,
				const MatrixVector *const Acc, MatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
			       MatrixVector *const Eff_Force );
/**
 * \brief Calculates the implicit displacement according to the formulation using the effective stiffness matrix.
 *
 * The displacement is computed according to the formulation using the effective stiffness matrix in page 53 of (\cite Dorka_1998). This is the implicit displacement.
 * 
 * \f[\overrightarrow{u_0} = \mathcal{K}^{-1}_e(\overrightarrow{f_{eff}} + \overrightarrow{l_i} + \overrightarrow{f(e)})\f]
 * 
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - \c Keinv must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower part in FORTRAN).
 * - The size of the vectors and matrices must be coherent, since it will not be checked in the routine.
 * 
 * \param[in] Eff_Force the effective force vector \f$\overrightarrow{f_{eff}}\f$.
 * \param[in] In_Load the input load vector \f$\overrightarrow{l_i}\f$.
 * \param[in] Err_Force the error force vector coming from the PID controller \f$\overrightarrow{f(e)}\f$.
 * \param[in] PID_P the proportional part of the PID controller.
 * \param[in] Keinv the inverted effective stiffness matrix \f$\mathcal{K}^{1}_e\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid allocating and deallocating memory each step. On input
 * only the number of rows is used.
 * \param[out] Disp0 is new implicit displacement \f$\overrightarrow{u_0}\f$.
 *
 * \post
 * - \c Disp0 is the result of the operation \f$\overrightarrow{u_0} = \mathcal{K}^{-1}_e(\overrightarrow{f_{eff}} + \overrightarrow{l_i} + \overrightarrow{f(e)})\f$.
 *
 * \sa MatrixVector and EffK_Calc_Effective_Force().
 */
void EffK_ComputeU0( const MatrixVector *const Eff_Force, const MatrixVector *const In_Load,
		     const MatrixVector *const Err_Force, const float PID_P, const MatrixVector *const Keinv,
		     MatrixVector *const Tempvec, MatrixVector *const Disp0 );

/**
 * \brief Copies the non-coupling part a vector.
 *
 * The non-coupling part a vector (\f$Order - Order_C\f$) is copied. This routine makes use of the level 1 BLAS dcopy_().
 *
 * \pre
 * - The global vector (length \f$Order\f$) must be properly initialised through the Init_MatrixVector() routine.
 * - The non-coupling vector must be of length \f$Order- OrderC\f$ and properly initialised through the Init_MatrixVector() routine.
 * - The number of rows of the vectors must be indicative of their length.
 * - The coupling nodes are assumed to be consecutive (the first node is assumed to be PosCouple).
 * - \e Order is the number of rows of the input vector.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * \param[in] VectorX The global vector.
 * \param[in,out] VectorXm The vector that will contain the non-coupling elements of \c VectorX. As an input,
 * only the size of the vector is referenced, not its elements.
 * \param[in] CNodes Structure containing the coupling nodes.
 *
 * \post
 * - \c VectorXm contains only the non-coupling nodes of \c Vec.
 *
 * \sa MatrixVector.
 */
void CreateVectorXm( const MatrixVector *const VectorX, MatrixVector *const VectorXm, const Coupling_Node *const CNodes );

/**
 * \brief Copies the coupling nodes of a vector.
 *
 * The coupling positions of a vector are copied. This routine makes use of the level 1 BLAS dcopy_().
 *
 * \pre
 * - The global vector (length \f$Order\f$) must be properly initialised through the Init_MatrixVector() routine.
 * - \f$Order > Order_C\f$.
 * - The coupling vector must be at least of length \f$Order_C\f$ and properly initialised.
 * - The coupling nodes are assumed to be consecutive (the first node is assumed to be PosCouple).
 * - The coupling nodes are assumed to be in increasing order of rows and in one based index.
 *
 * \param[in] VecX The global vector.
 * \param[out] VecXc The vector that will contain the coupling elements of \c VectorX.
 * \param[in] CNodes Structure containing the coupling nodes in increasing order of rows.
 *
 * \post
 * - \c VectorXc is a vector of length \$Order_C$ and contains only the coupling elements of \c Vec.
 *
 * \sa MatrixVector.
 */
void CreateVectorXc( const MatrixVector *const VecX, float *VecXc, const Coupling_Node *const CNodes );

/**
 * \brief Calculates the effective force vector \f$\overrightarrow{f_{eff}}\f$ according to the formulation using the effective stiffness matrix. This is the sparse version
 * of the routine and it requires Intel's MKL library.
 * 
 * The effective force vector is calculated according to the formilation using the effective stiffness matrix in page 53 (\cite Dorka_1998). It is used to compute the implicit
 * part of the displacement vector .
 *
 * \f[\overrightarrow{f_{eff}} = \mathcal{M}(a_0\vec{u} + a_2\vec{\dot{u}} + a_3\vec{\ddot{u}}) + \mathcal{C}(a_1\vec{u} + a_4\vec{\dot{u}} + a_5\vec{\ddot{u}})\f]
 *
 * The mass and viscous damping matrices are assumed to be in Intel's MKL CSR-\em three \em array \em variation format and only the upper part (lower part in FORTRAN) is
 * referenced. This routine therefore requires the program to be linked to the Intel's MKL library.
 *
 * \pre
 * - All vectors must be properly initialised through the Init_MatrixVector() routine.
 * - The mass and stiffness matrices must be in Intel's MKL CSR-\em three \em array \em variation format in one based index. Only the upper part (lower part in FORTRAN)
 *   is referenced.
 * - The size of the vectors and matrices must be coeherent, since it will not be checked in the routine.
 *
 * \param[in] Mass is the Mass matrix \f$\mathcal{M}\f$.
 * \param[in] Damp is the Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in] Disp The displacement vector \f$\vec{u}\f$.
 * \param[in] Vel is the velocity vector \f$\dot{\vec{u}}\f$.
 * \param[in] Acc is the acceleration vector \f$\ddot{\vec{u}}\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid allocating and deallocating memory each step. On input
 * only the number of rows is used
 * \param[in] a0 A constant of value \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$.
 * \param[in] a1 A constant of value \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$.
 * \param[in] a2 A constant of value \f$a_2 = \frac{1}{\beta_N\Delta t}\f$.
 * \param[in] a3 A constant of value \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$.
 * \param[in] a4 A constant of value \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$.
 * \param[in] a5 A constant of value \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$.
 * \param[out] Eff_Force is the effective force vector \f$\overrightarrow{f_{eff}}\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation \f$\overrightarrow{f_{eff}} = \mathcal{M}(a_0\vec{u} + a_2\vec{\dot{u}} + a_3\vec{\ddot{u}}) + \mathcal{C}(a_1\vec{u} + a_4\vec{\dot{u}} + a_5\vec{\ddot{u}})\f$.
 *
 * \sa MatrixVector and AlgConst.
 */
void EffK_Calc_Effective_Force_Sparse( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp,
				const MatrixVector *const Disp, const MatrixVector *const Vel,
				const MatrixVector *const Acc, MatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
				       MatrixVector *const Eff_Force );

#endif /* COMPUTEU0_H_ */
