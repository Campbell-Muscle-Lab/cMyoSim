#pragma once

/**
* @file		MyoSim_data.h
* @brief	header file for the MyoSim_data class
* @uathor	Ken Campbell
*/

#include "gsl_vector.h"
#include "gsl_matrix.h"

#include "global_definitions.h"

class options;
class model;
class hs_data;

class MyoSim_data
{
public:
	/**
	* Constructor
	* param integer number of time-points
	*/
	MyoSim_data(int no_of_time_points,
		options* set_p_options,
		model* set_p_model);

	/**
	* Destructor
	*/
	~MyoSim_data(void);

	/**
	* void write_data_to_delimited_file(char output_file_string[], char delimiter)

	* @param output_file_string[] a character array holding the output file name
	* @param delimiter a character holding the delimiter
	* @return void
	*/
	void write_data_to_delimited_file(std::string output_file_string, char delimiter = '\t');

	// Variables
	int no_of_time_points;		/**< integer number of time-points in the simulation */

	options* p_options;
								/**< pointer to an options object */

	model* p_model;
								/**< pointer to a model object
									 used to set the size of the matrices holding
									 the bs and cb state variables */

	gsl_vector* time;			/**< gsl_vector holding time in second at
									 each time-point */

	gsl_vector* m_length;		/**< gsl_vector holding muscle length (nm) for each time-point */

	gsl_vector* m_force;		/**< gsl_vector holding muscle force (N m^-2) for each time-point */

	gsl_vector* sc_extension;
								/**< gsl_vector holding series component extension (nm)
										for each time-step */

	gsl_vector* sc_force;
								/**< gsl_vector holding series component force
										in N m^-2 for each time-step */

	hs_data* p_hs_data[MAX_NO_OF_HALF_SARCOMERES];
								/**< pointer to an array of data objects for  individual half-sarcomeres */
};