#pragma once

/**
/* @file		muscle.h
/* @brief		Header file for a muscle object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

#include "gsl_vector.h"

// Forward declarations
class model;
class options;
class half_sarcomere;
class series_component;

class protocol;
class MyoSim_data;

//class heart_rate;
//class kinetic_scheme;

class muscle
{
public:
	/**
	* Constructor
	*/
	muscle(std::string JSON_model_file_string, std::string options_file_string);

	/**
	* Destructor
	*/
	~muscle(void);

	// Variables
	model* p_model [MAX_NO_OF_HALF_SARCOMERES];
											/**< array of pointers to model objects, one for each
													half-sarcomere */

	options* p_options;						/**< Pointer to an options object */

	protocol* p_protocol;					/**< pointer to a protocol object */

	MyoSim_data* p_data;					/**< pointer to a data structure */

	double cum_time_s;						/**< doubel with cumulative time in s */

	double m_length;						/**< double defining the length of the muscle in nm */

	double m_force;							/**< double defining the stress in the muscle in N m^-2 */

	series_component* p_sc;					/**< pointer to a series component */

	half_sarcomere* p_hs [MAX_NO_OF_HALF_SARCOMERES];
											/**< array of pointers to half-sarcomere objects */

	// Other functions
	void implement_protocol(std::string protocol_file_string, std::string results_file_string);

	void implement_time_step(int protocol_index);

	void implement_single_half_sarcomere(double time_step_s, double pCa, double sim_mode, double delta_hsl);

	void implement_myofibril(double time_step_s, double pCa, double sim_mode, double delta_hsl);

	double return_muscle_length(void);

	size_t length_control_myofibril_with_series_compliance(double time_step_s, double pCa, double delta_hsl);

	size_t worker_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* p, gsl_vector* f);

	void write_rates_file(void);
};

