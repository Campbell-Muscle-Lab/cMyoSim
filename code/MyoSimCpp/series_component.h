#pragma once

/**
/* @file		series_commponent.h
/* @brief		Header file for a muscle object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declarations
class muscle;

class series_component
{
public:
	/**
	* Constructor
	*/
	series_component(muscle* set_p_parent_m);

	/**
	* Destructor
	*/
	~series_component(void);

	// Variables
	muscle* p_parent_m;						/**< Pointer to a parent_muscle */

	double sc_k_stiff;						/**< double defining the series component stiffness in N m^-2 nm^-1 */

	double sc_sigma;
	double sc_L;

	double sc_extension;					/**< double defining the length of the series component in nm */

	double sc_force;						/**< double defining the stress in the series component in N m^-2 */

	// Other functions
	double return_series_force(double series_extension);

	double return_series_extension(double muscle_force);

};
