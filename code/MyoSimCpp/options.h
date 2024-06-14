#pragma once

/**
/* @file		options.h
/* @brief		Header file for an options object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

class options
{
public:
	/**
	* Constructor
	*/
	options(std::string options_file_string);

	/**
	* Destructor
	*/
	~options(void);

	// Variables
	std::string options_file_string;	/**< string holding the options file string */

	double max_rate;					/**< double holding max rate for a transition */

	double bin_x_min;					/**< double holding the min bin position (nm) */

	double bin_x_max;					/**< double holding the max bin position (nm) */

	double bin_x_width;					/**< double holding the bin width (nm) */

	double myofibril_force_tolerance;	/**< double defining the force tolerance for
												myofibril calculations */

	int myofibril_max_iterations;		/**< double defiining the max number of iterations
												for myofibril calculations */

	double myofibril_max_delta_hs_length;
										/**< double defiining the max length change to test
												during myofibril calculations */

	double movement_enhancement_width;	//*< double holding movement enhancement width */

	char rate_file_string[_MAX_PATH];	/**< char array holding rate file string */

	// Other functions

	void initialise_options_from_JSON_file_string(std::string JSON_options_file_string);
};
