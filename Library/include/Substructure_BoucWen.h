#ifndef _SUBSTRUCTURE_BOUCWEN_H_
#define _SUBSTRUCTURE_BOUCWEN_H_

#include "Definitions.h"

#define BOUCWEN_NUMPARAM_INIT           6   /*!< \brief Number of required parameters in order to initialise a
					     * substructure of type Bouc-Wen. */
#define BOUCWENDEG_NUMPARAM_INIT        12   /*!< \brief Number of required parameters in order to initialise a
					     * substructure of type Bouc-Wen with material degradation. */
#define BOUCWENBABERNOORI_NUMPARAM_INIT 18  /*!< \brief Number of required parameters in order to initialise a
					     * substructure of type Bouc-Wen-Baber-Noori. */

enum BoucWen_Id { BOUC_WEN,            /*!< \brief Classic Bouc-Wen model. */
		  BOUC_WEN_DEG,        /*!< \brief Bouc-Wen model with material degradation. */
		  BOUC_WEN_BABER_NOORI /*!< \brief Bouc-Wen-Baber-Noori model. */
};

typedef struct BoucWen {

     char *Description;   /*!< \brief Optional description of the substructure. */
     int BoucWen_Type;    /*!< \brief Controls wich type of Bouc Wen Model will be computed. */
     hysl_float_t BW_Force; /*!< \brief Force in the Bouc-Wen hysteresis law. */
     
     /*******************************************************************************************************/
     /************************************ Basic Bouc-Wen parameters ****************************************/
     /*******************************************************************************************************/
     
     hysl_float_t alpha;    /*!< \brief Post-yield stiffness ratio \f$\alpha = k_y/k_e\f$ with \f$k_y\$ the post
			   * yield stiffness and \f$ke\f$ de pre-yield stiffness. */
     hysl_float_t Fy;       /*!< \brief Yield force */
     hysl_float_t ko;       /* ko is the elasttic stiffness \f$k_0 = F_y/u_y\f$ where \f$F_y\f$ is the post
			   * yield stiffness and \f$u_y\f$ the yield displacement. */
     hysl_float_t beta;     /*!< \brief Bouc-Wen model coefficient. */
     hysl_float_t gamma;    /*!< \brief Bouc-Wen model coefficient. */
     hysl_float_t n;        /*!< \brief Hardening - Softening parameter. Controls the transition from linear to
			   * non-linear range (as n increases the transition becomes sharper; n is usually
			   * grater or equal to 1). */

     /*******************************************************************************************************/
     /******************************** Parameters for material degradation **********************************/
     /*******************************************************************************************************/
     hysl_float_t A0;       /*!< \brief Hysteresis amplitude. */
     hysl_float_t deltaA;   /*!< \brief Control parameter of the hysteresis amplitude with respect to the
			   * energy. */
     hysl_float_t nu0;      /*!< \brief Strength degradation. */
     hysl_float_t deltaNu;  /*!< \brief Strength degradation parameter. With \f$\delta_\nu = 0\f$ no strength
			   * degradation is included in the model. */
     hysl_float_t eta0;     /*!< \brief Stiffness degradation. */
     hysl_float_t deltaEta; /*!< \brief Stiffness degradation parameter. With \f$\delta_\eta = 0\f$ no stiffness
			   * degradation is included in the model.*/

     /*******************************************************************************************************/
     /********************************* Parameters for pinching modelling ***********************************/
     /*******************************************************************************************************/
     hysl_float_t vs0;      /*!< \brief Pinching severity. With \f$\zeta_s = 0\f$ there is no pinching effect
			   * included in the model. */
     hysl_float_t p;        /*!< \brief Initial pinching parameter. With \f$p = 0\f$ there is no pinching effect
			   * included in the model. */
     hysl_float_t q;        /*!< \brief Pinching parameter. */
     hysl_float_t psi0;     /*!< \brief Pinching parameter. */
     hysl_float_t deltaPsi; /*!< \brief Controls the change of pinching in the model */
     hysl_float_t lambda;   /*!< \brief Pinching parameter. */

     
     /*******************************************************************************************************/
     /************************************* Newton-Rhapson variables ****************************************/
     /*******************************************************************************************************/
     int maxIter;
     hysl_float_t tolerance, startPoint;
     hysl_float_t z_old, e_old;
     hysl_float_t DispT;
} BoucWen_t;

typedef struct BoucWenSurface {
     char *Description;  /*!< \brief Optional description of the substructure. */
     int BoucWen_Type;   /*!< \brief Controls wich type of Bouc Wen Model will be computed. */
     hysl_float_t BW_Force[2];
     
     /*******************************************************************************************************/
     /************************************ Basic Bouc-Wen parameters ****************************************/
     /*******************************************************************************************************/
     
     hysl_float_t alpha;    /*!< \brief Post-yield stiffness ratio \f$\alpha = k_y/k_e\f$ with \f$k_y\$ the post
			   * yield stiffness and \f$ke\f$ de pre-yield stiffness. */
     hysl_float_t Fy;       /*!< \brief Yield force */
     hysl_float_t ko;       /* ko is the elasttic stiffness \f$k_0 = F_y/u_y\f$ where \f$F_y\f$ is the post
			   * yield stiffness and \f$u_y\f$ the yield displacement. */
     hysl_float_t beta;     /*!< \brief Bouc-Wen model coefficient. */
     hysl_float_t gamma;    /*!< \brief Bouc-Wen model coefficient. */
     hysl_float_t n;        /*!< \brief Hardening - Softening parameter. Controls the transition from linear to
			   * non-linear range (as n increases the transition becomes sharper; n is usually
			   * grater or equal to 1). */

     /*******************************************************************************************************/
     /******************************** Parameters for material degradation **********************************/
     /*******************************************************************************************************/
     hysl_float_t A0;       /*!< \brief Hysteresis amplitude. */
     hysl_float_t deltaA;   /*!< \brief Control parameter of the hysteresis amplitude with respect to the
			   * energy. */
     hysl_float_t nu0;      /*!< \brief Strength degradation. */
     hysl_float_t deltaNu;  /*!< \brief Strength degradation parameter. With \f$\delta_\nu = 0\f$ no strength
			   * degradation is included in the model. */
     hysl_float_t eta0;     /*!< \brief Stiffness degradation. */
     hysl_float_t deltaEta; /*!< \brief Stiffness degradation parameter. With \f$\delta_\eta = 0\f$ no stiffness
			   * degradation is included in the model.*/

     /*******************************************************************************************************/
     /********************************* Parameters for pinching modelling ***********************************/
     /*******************************************************************************************************/
     hysl_float_t vs0;      /*!< \brief Pinching severity. With \f$\zeta_s = 0\f$ there is no pinching effect
			   * included in the model. */
     hysl_float_t p;        /*!< \brief Initial pinching parameter. With \f$p = 0\f$ there is no pinching effect
			   * included in the model. */
     hysl_float_t q;        /*!< \brief Pinching parameter. */
     hysl_float_t psi0;     /*!< \brief Pinching parameter. */
     hysl_float_t deltaPsi; /*!< \brief Controls the change of pinching in the model */
     hysl_float_t lambda;   /*!< \brief Pinching parameter. */

     
     /*******************************************************************************************************/
     /************************************* Newton-Rhapson variables ****************************************/
     /*******************************************************************************************************/
     int maxIter;
     hysl_float_t tolerance;
     hysl_float_t z[2], z_old[2];
     hysl_float_t DispT[2];
} BoucWenSurface_t;

void Substructure_BoucWen_Init (const hysl_float_t alpha, const hysl_float_t ko, const hysl_float_t Fy,
				const hysl_float_t beta, const hysl_float_t gamma, const hysl_float_t n,
				const hysl_float_t A0, const hysl_float_t deltaA, const hysl_float_t nu0,
				const hysl_float_t deltaNu, const hysl_float_t eta0, const hysl_float_t deltaEta,
				const hysl_float_t vs0, const hysl_float_t p, const hysl_float_t q,
			        const hysl_float_t lambda, const hysl_float_t psi0, const hysl_float_t deltaPsi,
				const int BoucWen_Type, const char *Description, BoucWen_t *const Sub );

void Substructure_BoucWenSurface_Init (const hysl_float_t alpha, const hysl_float_t ko,
				       const hysl_float_t Fy, const hysl_float_t beta, const hysl_float_t gamma, const hysl_float_t n,
				       const hysl_float_t A0, const hysl_float_t deltaA, const hysl_float_t nu0,
				       const hysl_float_t deltaNu, const hysl_float_t eta0, const hysl_float_t deltaEta,
				       const int BoucWen_Type, const char *Description, BoucWenSurface_t *const Sub );

void Substructure_BoucWen ( const hysl_float_t DispTdT, BoucWen_t *const Sub, hysl_float_t *const force );
void Substructure_BoucWenSurface ( const hysl_float_t DispTdT1, const hysl_float_t DispTdT2, BoucWenSurface_t *const Sub, hysl_float_t *const force1, hysl_float_t *const force2 );

void right_matrix_division (const hysl_float_t *const vector, const hysl_float_t *const matrix, hysl_float_t *const output);

void Substructure_BoucWen_Destroy ( BoucWen_t *const Sub );
void Substructure_BoucWenSurface_Destroy( BoucWenSurface_t *const Sub );

#endif /* _SUBSTRUCTURE_BOUCWEN_H_ */
