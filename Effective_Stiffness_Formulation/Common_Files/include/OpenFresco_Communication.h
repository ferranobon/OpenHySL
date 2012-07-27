/**
 * \file OpenFresco_Communication.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 21st of November 2011
 *
 * \brief Communication routines compliant with OpenFresco
 *
 * Contains the prototypes of the functions that deal with exchanging data using the
 * OpenFresco protocol.
 */

#ifndef OPENFRESCO_COMMUNICATION_H
#define OPENFRESCO_COMMUNICATION_H

int Communicate_With_OpenFresco( const float *const Data_To_Send, float *const Data_To_Receive, unsigned int Size, int WhatToDo );

#endif /* OPENFRESCO_COMMUNICATION_H */
