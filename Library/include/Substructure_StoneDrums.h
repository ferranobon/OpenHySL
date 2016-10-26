#ifndef _SUBSTRUCTURE_STONEDRUMS_H_
#define _SUBSTRUCTURE_STONEDRUMS_H_

#include "Substructure.h"
#include "Definitions.h"

#define STONEDRUMS_NUMPARAM_INIT 12

typedef struct StoneDrums {
    
     char *Description;  /*!< \brief Optional description of the substructure. */
     int BoucWen_Type;   /*!< \brief Controls wich type of Bouc Wen Model will be computed. */
     int PrevDOF;
     HYSL_FLOAT Result_Fc;

     /*******************************************************************************************************/
     /************************************ Basic Bouc-Wen parameters ****************************************/
     /*******************************************************************************************************/
     
     HYSL_FLOAT alpha;    /*!< \brief Post-yield stiffness ratio \f$\alpha = k_y/k_e\f$ with \f$k_y\$ the post
			   * yield stiffness and \f$ke\f$ de pre-yield stiffness. */
     double ko;           /* ko is the elasttic stiffness \f$k_0 = F_y/u_y\f$ where \f$F_y\f$ is the post
			   * yield stiffness and \f$u_y\f$ the yield displacement. */
     HYSL_FLOAT beta;     /*!< \brief Bouc-Wen model coefficient. */
     HYSL_FLOAT gamma;    /*!< \brief Bouc-Wen model coefficient. */
     HYSL_FLOAT n;        /*!< \brief Hardening - Softening parameter. Controls the transition from linear to
			   * non-linear range (as n increases the transition becomes sharper; n is usually
			   * grater or equal to 1). */

     /*******************************************************************************************************/
     /******************************** Parameters for material degradation **********************************/
     /*******************************************************************************************************/
     HYSL_FLOAT A0;       /*!< \brief Hysteresis amplitude. */
     HYSL_FLOAT deltaA;   /*!< \brief Control parameter of the hysteresis amplitude with respect to the
			   * energy. */
     HYSL_FLOAT nu0;      /*!< \brief Strength degradation. */
     HYSL_FLOAT deltaNu;  /*!< \brief Strength degradation parameter. With \f$\delta_\nu = 0\f$ no strength
			   * degradation is included in the model. */
     HYSL_FLOAT eta0;     /*!< \brief Stiffness degradation. */
     HYSL_FLOAT deltaEta; /*!< \brief Stiffness degradation parameter. With \f$\delta_\eta = 0\f$ no stiffness
			   * degradation is included in the model.*/

     /*******************************************************************************************************/
     /********************************* Parameters for pinching modelling ***********************************/
     /*******************************************************************************************************/
     HYSL_FLOAT vs0;      /*!< \brief Pinching severity. With \f$\zeta_s = 0\f$ there is no pinching effect
			   * included in the model. */
     HYSL_FLOAT p;        /*!< \brief Initial pinching parameter. With \f$p = 0\f$ there is no pinching effect
			   * included in the model. */
     HYSL_FLOAT q;        /*!< \brief Pinching parameter. */
     HYSL_FLOAT psi0;     /*!< \brief Pinching parameter. */
     HYSL_FLOAT deltaPsi; /*!< \brief Controls the change of pinching in the model */
     HYSL_FLOAT lambda;   /*!< \brief Pinching parameter. */

     
     /*******************************************************************************************************/
     /************************************* Newton-Rhapson variables ****************************************/
     /*******************************************************************************************************/
     int maxIter;
     HYSL_FLOAT tolerance, startPoint;
     HYSL_FLOAT z_old, e_old;
     HYSL_FLOAT DispT;
     
} StoneDrums_t;

void Substructure_StoneDrums_Init (const int PrevDOF, const HYSL_FLOAT alpha, const HYSL_FLOAT ko,
				   const HYSL_FLOAT beta, const HYSL_FLOAT gamma, const HYSL_FLOAT n,
				   const HYSL_FLOAT A0, const HYSL_FLOAT deltaA, const HYSL_FLOAT nu0,
				   const HYSL_FLOAT deltaNu, const HYSL_FLOAT eta0, const HYSL_FLOAT deltaEta,
				   const int BoucWen_Type, const char *Description, StoneDrums_t *const Sub );

void Substructure_StoneDrums ( HYSL_FLOAT DispTdT, StoneDrums_t *const Sub, HYSL_FLOAT *force );

void Substructure_StoneDrums_Destroy( StoneDrums_t *const Sub );

#endif /* _SUBSTRUCTURE_STONEDRUMS_H_ */
