/**
* @file		hs_data.cpp
* @brief	Source file for the hs_data class
* @author	ken Campbell
*/

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "MyoSim_data.h"
#include "hs_data.h"
#include "options.h"
#include "model.h"
#include "kinetic_scheme.h"

#include "gsl_vector.h"

using namespace std::filesystem;

// Constructor
hs_data::hs_data(int set_no_of_time_points,
	MyoSim_data* set_p_data = NULL,
	options* set_p_options = NULL,
	model* set_p_model = NULL)
{
	// Initialise

	p_data = set_p_data;

	no_of_time_points = set_no_of_time_points;
	p_options = set_p_options;
	p_model = set_p_model;

	// Allocate space for data vectors
	hs_time = gsl_vector_alloc(no_of_time_points);
	hs_pCa = gsl_vector_alloc(no_of_time_points);
	hs_length = gsl_vector_alloc(no_of_time_points);
	hs_command_length = gsl_vector_alloc(no_of_time_points);
	hs_slack_length = gsl_vector_alloc(no_of_time_points);
	hs_force = gsl_vector_alloc(no_of_time_points);
	hs_pas_int_force = gsl_vector_alloc(no_of_time_points);
	hs_viscous_force = gsl_vector_alloc(no_of_time_points);
	hs_pas_ext_force = gsl_vector_alloc(no_of_time_points);

	// Allocate space for data matrices
	hs_a_pops = gsl_matrix_alloc(no_of_time_points, 2);
	hs_m_pops = gsl_matrix_alloc(no_of_time_points, p_model->p_m_scheme[0]->no_of_states);

	// Set to zero
	gsl_vector_set_zero(hs_time);
	gsl_vector_set_zero(hs_pCa);
	gsl_vector_set_zero(hs_length);
	gsl_vector_set_zero(hs_command_length);
	gsl_vector_set_zero(hs_slack_length);
	gsl_vector_set_zero(hs_force);
	gsl_vector_set_zero(hs_pas_int_force);
	gsl_vector_set_zero(hs_viscous_force);
	gsl_vector_set_zero(hs_pas_ext_force);

	gsl_matrix_set_zero(hs_a_pops);
	gsl_matrix_set_zero(hs_m_pops);
}

// Destructor
hs_data::~hs_data(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(hs_time);
	gsl_vector_free(hs_pCa);
	gsl_vector_free(hs_length);
	gsl_vector_free(hs_command_length);
	gsl_vector_free(hs_slack_length);
	gsl_vector_free(hs_force);
	gsl_vector_free(hs_pas_int_force);
	gsl_vector_free(hs_viscous_force);
	gsl_vector_free(hs_pas_ext_force);

	// Then the matrices
	gsl_matrix_free(hs_a_pops);
	gsl_matrix_free(hs_m_pops);
}
