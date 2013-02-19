/**
 * \file Error_Compensation.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 18th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the error force compensation. The routines make use of the BLAS library to
 * perform the linear algebra operations and they support both single and double precision. Sparse BLAS
 * operations are supported through the Intel MKL library.
 */

#ifndef ERROR_COMPENSATION_H_
#define ERROR_COMPENSATION_H_

#include "MatrixVector.h"

/**
 * \brief Stores the proportional, integral and derivative constant of the PID compensator.
 */
typedef struct {
     double P;  /*!< \brief Proportional constant.*/
     double I;  /*!< \brief Integral constant.*/
     double D;  /*!< \brief Derivative constant.*/
} PID_t;

/**
 * \brief Error force compensation using a PID controller.
 *
 * This routine implements the PID controller for calculating the error force. The operation performed is:
 *
 * \f[\vec e^{t+\Delta t} = (\vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t} + \vec l_i^{t} - (\mathcal M \ddot{\vec u}^{t + \Delta t} + \mathcal C \dot{\vec u}^{t+\Delta t} + \mathcal K u^{t + \Delta t})\f]
 * \f[\vec f_e^{t + \Delta t} = \biggl[\vec e^t + I\Delta t\sum_i^t \vec e^t + \frac{D}{\Delta t} (\vec e^t - \vec e^{t-\Delta t})\biggr]\f]
 *
 * where:
 * - \f$\vec e^t\f$ is the equilibrium error at the end of the time step,
 * - \f$\vec f_e^t\f$ is the compensation force,
 * - \f$f_r^{t+\Delta t}\f$ and \f$f_s^{t+\Delta t}\f$ are the calculated and measured force vectors at time \f$t + \Delta t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$.
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\dot{\vec u}^{t+\Delta t}\f$ is the velocity vector at time \f$t + \Delta t\f$,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u^{t+\Delta t}\f$ is the displacement vector at time \f$t + \Delta t\f$,
 * - and \f$P\f$, \f$I\f$ and \f$D\f$ are the proportional, integral and derivative constants of the PIR error compensator.
 * 
 * It makes use of BLAS routines to perform the lineal algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The sizes of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The PID constants must be properly initialised.
 *
 * \param[in] Mass The mass matrix.
 * \param[in] Damp The proportional viscous damping matrix.
 * \param[in] Stiff The stiffness matrix.
 * \param[in] AccTdT is the acceleration vector \f$\ddot{\vec u}^{t+\Delta t}\f$.
 * \param[in] VelTdT is the velocity vector \f$\dot{\vec u}^{t+\Delta t}\f$.
 * \param[in] DispTdT is the displacement vector \f$\vec u^{t+\Delta t}\f$.
 * \param[in] fc is the coupling force vector \f$\vec f_c^{t+\Delta t} = \vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t}\f$.
 * \param[in] LoadTdT is the input load vector \f$l_i^t\f$.
 * \param[in] PID Values of the PID controller and it contains a proportional \f$P\f$, integral \f$I\f$ and derivative \f$D\f$ constants.
 * \param[in,out] fe is the compensation force vector \f$\vec fe^t\f$. As an input, only the size of the matrix is referenced, not its elements.
 *
 * \post \c fe is the compensation force vector and contains the result of:
 *
 * \f[\vec e^{t+\Delta t} = (\vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t} + \vec l_i^{t} - (\mathcal M \ddot{\vec u}^{t + \Delta t} + \mathcal C \dot{\vec u}^{t+\Delta t} + \mathcal K u^{t + \Delta t})\f]
 * \f[\vec f_e^{t + \Delta t} = \biggl[\vec e^t + I\Delta t\sum_i^t \vec e^t + \frac{D}{\Delta t} (\vec e^t - \vec e^{t-\Delta t})\biggr]\f]
 *
 * \sa MatrixVector and PID_t.
 */
void ErrorForce_PID( const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *Stiff,
		     const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
		     const MatrixVector *const fc, const MatrixVector *const LoadTdT, const PID_t *const PID, MatrixVector *const fe );

/**
 * \brief Error force compensation using a PID controller. Sparse version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
 *
 * This routine implements the PID controller for calculating the error force. The operation performed is:
 *
 * \f[\vec e^{t+\Delta t} = (\vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t} + \vec l_i^{t} - (\mathcal M \ddot{\vec u}^{t + \Delta t} + \mathcal C \dot{\vec u}^{t+\Delta t} + \mathcal K u^{t + \Delta t})\f]
 * \f[\vec f_e^{t + \Delta t} = \biggl[\vec e^t + I\Delta t\sum_i^t \vec e^t + \frac{D}{\Delta t} (\vec e^t - \vec e^{t-\Delta t})\biggr]\f]
 *
 * where:
 * - \f$\vec e^t\f$ is the equilibrium error at the end of the time step,
 * - \f$\vec f_e^t\f$ is the compensation force,
 * - \f$f_r^{t+\Delta t}\f$ and \f$f_s^{t+\Delta t}\f$ are the calculated and measured force vectors at time \f$t + \Delta t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$.
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\dot{\vec u}^{t+\Delta t}\f$ is the velocity vector at time \f$t + \Delta t\f$,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u^{t+\Delta t}\f$ is the displacement vector at time \f$t + \Delta t\f$,
 * - and \f$P\f$, \f$I\f$ and \f$D\f$ are the proportional, integral and derivative constants of the PIR error compensator.
 * 
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - All elements of type \c MatrixVector_Sp must be properly intialised through the Init_MatrixVector_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The sizes of the matrices must be the identical.
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The PID constants must be properly initialised.
 *
 * \param[in] Mass The mass matrix.
 * \param[in] Damp The proportional viscous damping matrix.
 * \param[in] Stiff The stiffness matrix.
 * \param[in] AccTdT is the acceleration vector \f$\ddot{\vec u}^{t+\Delta t}\f$.
 * \param[in] VelTdT is the velocity vector \f$\dot{\vec u}^{t+\Delta t}\f$.
 * \param[in] DispTdT is the displacement vector \f$\vec u^{t+\Delta t}\f$.
 * \param[in] fc is the coupling force vector \f$\vec f_c^{t+\Delta t} = \vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t}\f$.
 * \param[in] LoadTdT is the input load vector \f$l_i^t\f$.
 * \param[in] PID Values of the PID controller and it contains a proportional \f$P\f$, integral \f$I\f$ and derivative \f$D\f$ constants.
 * \param[in,out] fe is the compensation force vector \f$\vec fe^t\f$. As an input, only the size of the matrix is referenced, not its elements.
 *
 * \post \c fe is the compensation force vector and contains the result of:
 *
 * \f[\vec e^{t+\Delta t} = (\vec f_r^{t + \Delta t} + \vec f_s^{t+\Delta t} + \vec l_i^{t} - (\mathcal M \ddot{\vec u}^{t + \Delta t} + \mathcal C \dot{\vec u}^{t+\Delta t} + \mathcal K u^{t + \Delta t})\f]
 * \f[\vec f_e^{t + \Delta t} = \biggl[\vec e^t + I\Delta t\sum_i^t \vec e^t + \frac{D}{\Delta t} (\vec e^t - \vec e^{t-\Delta t})\biggr]\f]
 *
 * \sa MatrixVector, MatrixVector_Sp and PID_t.
 */
void ErrorForce_PID_Sp( const MatrixVector_Sp *const Mass, const MatrixVector_Sp *const Damp, const MatrixVector_Sp *Stiff,
			const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			const MatrixVector *const fc, const MatrixVector *const LoadTdT, const PID_t *const PID, MatrixVector *const fe );
#endif /* ERRORCOMPENSATION_H_*/

