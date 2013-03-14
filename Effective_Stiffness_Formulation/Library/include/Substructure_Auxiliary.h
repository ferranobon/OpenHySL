#ifndef SUBSTRUCTURE_AUXILIARY_H_
#define SUBSTRUCTURE_AUXILIARY_H_

#include "MatrixVector.h"
#include "Substructure.h"

/**
 * \brief Joins the non-coupling of a vector.
 *
 * The non-coupling (\f$Order-OrderC\f$) part of the vector is added to the global vector of size (\f$Order\f$). It performs the
 * operation:
 *
 * \f[\hat{\vec n}^{t + \Delta t} = \mathcal{\hat G}^{-1} (\hat{\vec f}_r^{n_{sub}-1} + \hat{\vec f}_s^{n_{sub}-1}) = \mathcal{\hat G}^{-1} \hat{\vec f}_c^{n_{sub}-1}\f]
 * \f[\hat{\vec n}^{t + \Delta t} \longrightarrow \vec n^{t + \Delta t}\f]
 *
 * where:
 * - \f$\hat{\vec n}^{t + \Delta t}\f$ is the new displacement of the non-coupling nodes after applying the restoring forces,
 * - \f$\mathcal{G}\f$ is a part of the gain matrix with the non-coupling elements of a column with a node marked as coupling node.
 * - \f$f_r^{i-1}\f$ and \f$f_s^{i-1}\f$ are the calculated and measured force vectors at sub-step \f$n_{sub}-1\f$, with $n_{sub} being the number of sub-steps,
 * - and \f$n^{t + \Delta t}\f$ is the new displacement vector including both parts: the coupling and non-coupling part.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - \c Gain_m must be a matrix of size \f$Order\f$x\f$OrderC\f$ and have the output of Build_MatrixXm().
 * - \c VecTdT_m must be of length \f$Order - OrderC\f$.
 * - \c VecTdT must be of length \f$Order\f$.
 * - \c fcprevsub must be of length \f$OrderC\f$.
 * - The coupling nodes are assumed to be consecutive and in increasing equation number.
 * - \e Order is the number of rows and columns of the global matrix.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \param[in] VecTdT_m The non-coupling part of the vector (displacement, velocity or acceleration depending on the formulation).
 * \param[in] Gain_m The non-coupling part of the gain matrix \f$\mathcal G\f$.
 * \param[in] fcprevsub The non-coupling part of the vector.
 * \param[in] CNodes Structure containing the coupling nodes.
 * \param[in,out] VecTdT The non-coupling part of the vector. Only the non-coupling positions are modified.
 *
 * \post
 * - \c VecTdT contains the non-coupling part of the vector, leaving the coupling nodes untouched. This is accomplished through: 
 *
 * \f[\hat{\vec n}^{t + \Delta t} = \mathcal{\hat G}^{-1} (\hat{\vec f}_r^{n_{sub}-1} + \hat{\vec f}_s^{n_{sub}-1}) = \mathcal{\hat G}^{-1} \hat{\vec f}_c^{n_{sub}-1}\f]
 * \f[\hat{\vec n}^{t + \Delta t} \longrightarrow \vec n^{t + \Delta t}\f]
 *
 *
 * \sa MatrixVector_t.
 */
void Substructure_JoinNonCouplingPart( MatrixVector_t *const VecTdT_m, const MatrixVector_t *const Gain_m,
			   const MatrixVector_t *const fcprevsub, const CouplingNode_t *const CNodes,
			   MatrixVector_t *const VecTdT );

/**
 * \brief Construction of the coupling matrix. General storage version.
 *
 * This routine copies the values in the coupling positions of a symmetric matrix
 * constructing a sub-matrix of \f$ Size = Number~of~coupling~nodes^2\f$. Since this
 * matrix may have to be sent to the experimental facility in order to perform the
 * sub-stepping process, it is treated as a full matrix (although it is
 * symmetrical). Therefore, all elements are copied if \f$Number~of~coupling~nodes \geq
 * 1\f$. For example, given the symmetric matrix \f$\mathcal{A}\f$ and coupling positions
 * in 2 and 5, the resulting \f$\mathcal{A}_{Couple}\f$ would be as follows:
 *
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & -1 & 3  & 4  & 5\\
 *   * & \mathbf{5}  & 4  & -3 & \mathbf{2}\\
 *   * & *  & 3  & 6  & 7\\
 *   * & *  & *  & 11 & 4\\
 *   * & *  & *  & *  & \mathbf{8}\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{Couple} = \begin{pmatrix}
 * 5 & 2\\
 * 2 & 8\\
 * \end{pmatrix}\f]
 *
 * If \c Mat is in packed storage, the routine Substructure_MatrixXc_PS() should be used
 * instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - \c Mat has to be a symmetrical matrix containing at least the upper part (lower part
 *   in FORTRAN routines) in general storage.
 * - The size of \c MatCouple should be \f$Size \geq Number~of~coupling~nodes^2\f$.
 * - The coupling nodes have to be properly initialised through the
 *   Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based
 *   index.
 *
 * \param[in]  Mat       The matrix that will be decoupled.
 * \param[in]  CNodes    Structure containing the coupling nodes in increasing order of
 *                       rows.
 * \param[out] MatCouple The matrix where the coupling nodes are saved.
 *
 * \post \c MatCouple is a symmetrical matrix \f$Size = Number~of~coupling~nodes^2\f$ in
 * general storage that contains the values in the coupling position of the specified
 * matrix. All the values are referenced.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 *
 */
void Substructure_MatrixXc( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes, MatrixVector_t *const MatCouple );

/**
 * \brief Construction of the coupling matrix. General storage version.
 *
 * This routine copies the values in the coupling positions of a symmetric matrix in
 * packed storage constructing a sub-matrix of \f$ Size =
 * Number~of~coupling~nodes^2\f$. Since this matrix may have to be sent to the
 * experimental facility in order to perform the sub-stepping process, it is treated as a
 * full matrix (although it is symmetrical). Therefore, all elements are copied if
 * \f$Number~of~coupling~nodes \geq 1\f$. For example, given the symmetric matrix
 * \f$\mathcal{A}\f$ and coupling positions in 2 and 5, the resulting
 * \f$\mathcal{A}_{Couple}\f$ would be as follows:
 *
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & -1 & 3  & 4  & 5\\
 *     & \mathbf{5}  & 4  & -3 & \mathbf{2}\\
 *     &    & 3  & 6  & 7\\
 *     &    &    & 11 & 4\\
 *     &    &    &    & \mathbf{8}\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{Couple} = \begin{pmatrix}
 * 5 & 2\\
 * 2 & 8\\
 * \end{pmatrix}\f]
 *
 * If \c Mat is in general storage, the routine Substructure_MatrixXc() should be used
 * instead.
 *
 * \pre
 * - \c Mat must be initialised through the MatrixVector_Create_PS() routine.
 * - \c Mat must be symmetrical and in packed storage. The upper triangular part packed in
 *   rows (lower triangular part packed in columns in FORTRAN) must be present.
 * - \c Matcouple must be in general storage and properly initialised through the
 *   MatrixVector_Create() routine.
 * - The size of \c MatCouple should be \f$Size \geq Number~of~coupling~nodes^2\f$.
 * - The coupling nodes have to be properly initialised through the
 *   Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based
 *   index.
 *
 * \param[in]  Mat       The matrix that will be decoupled in packed storage.
 * \param[in]  CNodes    Structure containing the coupling nodes in increasing order of
 *                       rows.
 * \param[out] MatCouple The matrix where the coupling nodes are saved.
 *
 * \post \c MatCouple is a symmetrical matrix \f$Size = Number~of~coupling~nodes^2\f$ in
 * general storage that contains the values in the coupling position of the specified
 * matrix. All the values are referenced.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 *
 */
void Substructure_MatrixXc_PS( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes, MatrixVector_t *const MatCouple );

/**
 * \brief Construction of the non-coupling part of a given matrix. General storage
 * version.
 *
 * This routine copies the non-coupling values of a column with coupling degrees of
 * freedom of the symmetric matrix \c Mat, constructing a matrix of \f$Size = (Order -
 * Order_C)\cdot Order_C\f$. Where:
 *
 * - \f$Order\f$ is the number of rows and columns of the input matrix.
 * - \f$Order_C\f$ is the number of coupling degrees of freedom.
 *
 * For example, given the given the symmetric matrix \f$\mathcal{A}\f$ and coupling
 * positions in 2 and 5, the resulting \f$\mathcal{A}_{cm}\f$ would be as follows:
 * 
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & \mathbf{-1} & 3  & 4  & \mathbf{5}\\
 *   * & 5  & \mathbf{4}  & \mathbf{-3} & 2\\
 *   * & *  & 3  & 6  & \mathbf{7}\\
 *   * & *  & *  & 11 & \mathbf{4}\\
 *   * & *  & *  & *  & 8\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{cm} = \begin{pmatrix}
 * -1 & 5\\
 *  4 & 7\\
 * -3 & 4\\
 * \end{pmatrix}\f]
 *
 * It makes use of BLAS routines to perform the linear algebra operations. If \c Mat is in
 * packed storage, the routine Substructure_MatrixXcm_PS() should be used instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - \c Mat has to be a symmetrical matrix containing at least the upper part (lower part
 *   in FORTRAN routines) in general storage.
 * - \f$Order > Order_C\f$.
 * - \c Mat must be of \f$Size = Order\cdot Order\f$.
 * - \c MatXcm must be of \f$Size = (Order - Order_C)\cdot Order_C\f$.
 * - The coupling nodes have to be properly initialised through the
 *   Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based
 *   index.
 *
 * \param[in]     Mat     The matrix that will be decoupled in general storage.
 * \param[in]     CNodes  Structure containing the coupling nodes in increasing order of
 *                        rows.
 * \param[in,out] MatXcm  The matrix where the non-coupling elemets of a column with a
 *                        coupling node are stored.
 *
 * \post \c MatXcm is a general matrix of \f$Size = (Order - Order_C)\cdot Order_C\f$ with
 * the non-coupling elements of the columns with coupling nodes.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 */
void Substructure_MatrixXcm( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes, MatrixVector_t *const MatXcm );

/**
 * \brief Construction of the non-coupling part of a given matrix. Packed storage version.
 *
 * This routine copies the non-coupling values of a column with coupling degrees of
 * freedom of the symmetric matrix \c Mat in packed storage, constructing a matrix of
 * \f$Size = (Order - Order_C)\cdot Order_C\f$. Where:
 *
 * - \f$Order\f$ is the number of rows and columns of the input matrix.
 * - \f$Order_C\f$ is the number of coupling degrees of freedom.
 *
 * For example, given the given the symmetric matrix \f$\mathcal{A}\f$ in packed storage
 * and coupling positions in 2 and 5, the resulting \f$\mathcal{A}_{cm}\f$ would be as
 * follows:
 * 
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & \mathbf{-1} & 3  & 4  & \mathbf{5}\\
 *     & 5  & \mathbf{4}  & \mathbf{-3} & 2\\
 *     &    & 3  & 6  & \mathbf{7}\\
 *     &    &    & 11 & \mathbf{4}\\
 *     &    &    &    & 8\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{cm} = \begin{pmatrix}
 * -1 & 5\\
 *  4 & 7\\
 * -3 & 4\\
 * \end{pmatrix}\f]
 *
 * If \c Mat is in general storage, the routine Substructure_MatrixXcm() should be used
 * instead.
 *
 * \pre
 * - \c Mat must be initialised through the MatrixVector_Create_PS() routine.
 * - \c Mat must be symmetrical and in packed storage. The upper triangular part packed in
 *   rows (lower triangular part packed in columns in FORTRAN) must be present.
 * - \c MatXcm must be in general storage and properly initialised through the
 *   MatrixVector_Create() routine.
 * - \f$Order > Order_C\f$.
 * - \c Mat must be of \f$Size = Order\cdot Order\f$.
 * - \c MatXcm must be of \f$Size = (Order - Order_C)\cdot Order_C\f$.
 * - The coupling nodes have to be properly initialised through the
 *   Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based
 *   index.
 *
 * \param[in]     Mat     The matrix that will be decoupled in packed storage.
 * \param[in]     CNodes  Structure containing the coupling nodes in increasing order of
 *                        rows.
 * \param[in,out] MatXcm  The matrix where the non-coupling elemets of a column with a
 *                        coupling node are stored.
 *
 * \post \c MatXcm is a general matrix of \f$Size = (Order - Order_C)\cdot Order_C\f$ with
 * the non-coupling elements of the columns with coupling nodes.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 */
void Substructure_MatrixXcm_PS( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes, MatrixVector_t *const MatXcm );

/**
 * \brief Copies the non-coupling part a vector.
 *
 * The non-coupling part a vector (\f$Order - Order_C\f$) is copied. It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - The global vector \c VectorX (length \f$Order\f$) must be properly initialised through the MatrixVector_Create() routine.
 * - The non-coupling vector \c VectorXm must be of length \f$Order- OrderC\f$ and properly initialised through the MatrixVector_Create() routine.
 * - The number of rows of the vectors must be indicative of their length.
 * - The coupling nodes have to be properly initialised through the Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based index.
 * - \e Order is the number of rows of the input vector.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * \param[in] VectorX The global vector.
 * \param[in] CNodes Structure containing the coupling nodes.
 * \param[in,out] VectorXm The vector that will contain the non-coupling elements of \c VectorX. As an input,
 * only the size of the vector is referenced, not its elements.
 *
 * \post
 * - \c VectorXm contains only the non-coupling nodes of \c Vec.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 */
void Substructure_VectorXm( const MatrixVector_t *const VectorX, const CouplingNode_t *const CNodes, MatrixVector_t *const VectorXm );

/**
 * \brief Copies the coupling nodes of a vector.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - \f$Order > Order_C\f$.
 * - The coupling vector \c VecXc must be at least of length \f$Order_C\f$.
 * - The coupling nodes have to be properly initialised through the Substructure_ReadCouplingNodes() routine.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based index.
 *
 * \param[in] VecX The global vector.
 * \param[in] CNodes Structure containing the coupling nodes in increasing order of rows.
 * \param[out] VecXc The vector that will contain the coupling elements of \c VectorX.
 *
 * \post
 * - \c VecXc is a vector of length \f$Order_C\f$ and contains only the coupling elements of \c VecX.
 *
 * \sa MatrixVector_t and CouplingNode_t.
 */
void Substructure_VectorXc( const MatrixVector_t *const VecX, const CouplingNode_t *const CNodes, MatrixVector_t *const VecXc );

#endif /* SUBSTRUCTURE_AUXILIARY_H_ */
