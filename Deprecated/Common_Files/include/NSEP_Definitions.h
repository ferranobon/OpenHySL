/**
 * \file NSEP_Definitions.h
 * \author Ferran Obón Santacana
 * \author Kung-Juin Wang
 * \version 1.0
 * \date 4th of November 2011
 *
 * \brief Definition of various NSEP and PNSE constants
 *
 * This file contains the definition of the constants of the NSEP protocol. This
 * definitions are used to identify the type of message that is sent between NSEP
 * compliant facilities.
 *
 */

#ifndef NSEP_DEFINITIONS_H
#define NSEP_DEFINITIONS_H

/* login information related constants */
#define __AP_ERROR			        0
#define __AP_LOGIN_REQUEST	                1
#define __AP_LOGIN_REPLY	                2
#define __AP_LOGIN_REPLY_NOTLOGINPACKET		0
#define __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT	1
#define __AP_LOGIN_REPLY_WRONGPASSWORD		2
#define __AP_LOGIN_REPLY_EXCEEDTRIALTIME	3
#define __AP_LOGIN_REPLY_ALREADYLOGGEDIN	4
#define __AP_LOGIN_REPLY_OK			5

#define __AP_END				2

/* NSEP related constants */
#define NSEP_SETPUSH	(__AP_END+1)
#define NSEP_QUERY	(__AP_END+2)
#define NSEP_EXPINFO	(__AP_END+3)
#define NSEP_EXPSTATE	(__AP_END+4)
#define NSEP_CLNSTATE	(__AP_END+5)
#define NSEP_CMD	(__AP_END+6)
#define NSEP_CSIG	(__AP_END+7)
#define NSEP_IM		(__AP_END+8)

/* NSEP client state constants */
#define NSEP_CS_DISCONNECTED  0
#define NSEP_CS_NOTREADY      1
#define NSEP_CS_RUNNING       2
#define NSEP_CS_FINISHED      3

/* Experimental State */
#define NSEP_ES_NOTREADY     0
#define NSEP_ES_RUNNING      1
#define NSEP_ES_INTERRUPTED  2
#define NSEP_ES_FINISHED     3

#endif /* NSEP_DEFINITIONS_H */
