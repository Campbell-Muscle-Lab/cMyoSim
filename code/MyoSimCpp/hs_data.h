#pragma once

/**
* @file		hs_data.h
* @brief	header file for the hs_data class
* @uathor	Ken Campbell
*/

#include "gsl_vector.h"
#include "gsl_matrix.h"

class MyoSim_data;
class options;
class model;

class hs_data
{
public:
	/**
	* Constructor
	* param integer number of time-points
	*/
	hs_data(int no_of_time_points,
		MyoSim_data* set_p_data,
		options* set_p_options,
		model* set_p_model);

	/**
	* Destructor
	*/
	~hs_data(void);

	// Variables
	int no_of_time_points;		/**< integer number of time-points in the simulation */

	int hs_id;					/**< integer holding the id of the half-sarcomere */

	MyoSim_data* p_data;				/**< pointer to a data object */

	options* p_options;			/**< pointer to an options object */

	model* p_model;				/**< pointer to a model object used to set the size of the matrices
										holding the bs and cb state variables */

	gsl_vector* hs_time;		/**< gsl_vector holding time in s for each time-point */

	gsl_vector* hs_pCa;			/**< gsl_vector holding pCa for each time-point */

	gsl_vector* hs_length;		/**< gsl_vector holding hs_length (nm) for each time-point */

	gsl_vector* hs_command_length;
								/**< gsl_vector holding command_length for each time point */

	gsl_vector* hs_slack_length;
								/**< gsl_vector holding slack length for each time-point */

	gsl_vector* hs_force;		/**< gsl_vector holding hs_force for each time-point */

	gsl_vector* hs_pas_int_force;
								/**< gsl_vector holding passive internal force for each time-point */

	gsl_vector* hs_viscous_force;
								/**< gsl_vector holding viscous force within half-sarcomere */

	gsl_vector* hs_pas_ext_force;
								/**< gsl_vector holding hs_extracellular_force for each time-point */

	gsl_matrix* hs_a_pops;		/**< gsl_matrix holding the proportion of binding sites
									 in each state at each time-point */

	gsl_matrix* hs_m_pops;		/**< gsl_matrix holding the proportion of cross-bridges
									 in each state at each time-point */
};
