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
class muscle;
class model;
class myofilaments;
class membranes;

class half_sarcomere
{
public:
	/**
	* Constructor
	*/
	half_sarcomere(muscle* set_p_parent_m, model* set_p_model, int set_hs_index);

	/**
	* Destructor
	*/
	~half_sarcomere(void);

	// Variables
	muscle* p_parent_m;						/**< pointer to a parent_muscle */

	model* p_model;							/**< pointer to a model for the half-sarcomere */

	int hs_index;							/**< integer holding the index number for the half-sarcomere */

	double time_s;							/**< time for this sarcomere */

	double pCa;								/**< pCa in the half-sarcomere */

	double hs_length;						/**< double defining the actual length of the half_sarcomere in nm */

	double hs_command_length;				/**< double defining the command length of the half_sarcomere in nm */

	double hs_slack_length;					/**< double defining the slack length of the half_sarcomere in nm */

	double hs_force;						/**< double defining the stress in the half_sarcomere in N m^-2 */

	double hs_pas_int_force;				/**< double defining the internal passive force in N m^-2 */

	double hs_pas_ext_force;				/**< double defining the external passive force in N m^-2 */

	double hs_viscous_force;				/**< double defining the viscous force in N m^-2 */

	double hs_velocity;						/**< double defining the velocity of the half-sarcomere in nm s^-1 */

	double hs_prop_fibrosis;				/**< double defining the proportion of the half-sarcomere
													cross-sectional area that is fibrosis */

	double hs_prop_myofilaments;			/**< double defining the proportion of the non-fibrotis
													cross-sectional area that is myofilaments */

	double hs_m_head_density;				/**< double defining the number of heads in a half-sarcomere
													with a cross-sectional area of 1 m^2 */

	myofilaments* p_myof;					/**< pointer to a myofilament oject */
	
	membranes* p_memb;						/**< pointer to a membranes object */

	// Other functions
	void sarcomere_kinetics(double time_step_s, double pCa);
											// Implements kinetics

	void change_hs_length(double delta_hsl, double time_step_s);

	double return_hs_length_for_stress(double target_stress, double time_step_s);

	double return_stress_after_delta_hsl(double delta_hsl);
};

