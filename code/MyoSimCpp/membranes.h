#pragma once

#pragma once

/**
/* @file		muscle.h
/* @brief		Header file for a muscle object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declarations
class model;
class options;
class half_sarcomere;
class series_component;

class membranes
{
public:
	/**
	* Constructor
	*/
	membranes(half_sarcomere* p_parent_hs);

	/**
	* Destructor
	*/
	~membranes(void);

	// Variables
	half_sarcomere* p_parent_hs;					/**< pointer to the parent half-sarcomere */

	double memb_Ca_cytosol;							/**< double with cytosolic Ca concentration in M */

	// Other functions

};


