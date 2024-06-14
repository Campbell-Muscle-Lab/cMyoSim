#pragma once

/**
 * @file    protocol.h
 * @brief   header file for the protocol
 * @author  Ken Campbell
 */

#include "global_definitions.h"

#include "gsl_vector.h"

class protocol
{
public:
    /**
    * Constructor - creates a protocol from a tab-delimited text file
    * param protocol_file_string a char array defining the protocol
    */
    protocol(std::string protocol_file_string);

    /**
     * Destructor
     */
    ~protocol(void);


    // Variables
    int no_of_time_points;          /**< integer number of time-points */

    gsl_vector* dt;                 /**< gsl_vector holding dt for
                                         each time-point */
    gsl_vector* delta_hsl;          /**< gsl_vector holding delta-hsl for
                                         each time-point */
    gsl_vector* sim_mode;           /**< gsl_vector holding sim_mode for
                                         each time-point */
    gsl_vector* pCa;                /**< gsl_vector holding pCa for
                                         each time-point */
};
