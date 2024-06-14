/**
* @file		MyoSim_data.cpp
* @brief	Source file for the MyoSim_data class
* @author	ken Campbell
*/

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "MyoSim_data.h"
#include "options.h"
#include "model.h"
#include "hs_data.h"

//#include "hs_data.h"
#include "kinetic_scheme.h"

#include "gsl_vector.h"
#include "gsl_math.h"

using namespace std::filesystem;

// Constructor
MyoSim_data::MyoSim_data(int set_no_of_time_points,
	options* set_p_options = NULL,
	model* set_p_model = NULL)
{
	// Initialise

	no_of_time_points = set_no_of_time_points;
	p_options = set_p_options;
	p_model = set_p_model;

	// Allocate space for data vectors
	time = gsl_vector_alloc(no_of_time_points);
	m_length = gsl_vector_alloc(no_of_time_points);
	m_force = gsl_vector_alloc(no_of_time_points);

	sc_extension = gsl_vector_alloc(no_of_time_points);
	sc_force = gsl_vector_alloc(no_of_time_points);

	// Set to zero
	gsl_vector_set_zero(time);
	gsl_vector_set_zero(m_length);
	gsl_vector_set_zero(m_force);

	gsl_vector_set_all(sc_extension, GSL_NAN);
	gsl_vector_set_all(sc_force, GSL_NAN);


	// Now create the data objects for each half-sarcomere
	for (int hs_counter = 0; hs_counter < p_model->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs_data[hs_counter] = new hs_data(no_of_time_points,
			this,
			p_options,
			p_model);
	}
}

// Destructor
MyoSim_data::~MyoSim_data(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(time);
	gsl_vector_free(m_length);
	gsl_vector_free(m_force);
	gsl_vector_free(sc_extension);
	gsl_vector_free(sc_force);

	// Now the half-sarcomere data structures
	for (int hs_counter = 0; hs_counter < p_model->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_hs_data[hs_counter];
	}
}

// Functions

void MyoSim_data::write_data_to_delimited_file(std::string output_file_string, char delimiter)
{
	//! Writes data to a delimited file

	// Variables
	char file_string[_MAX_PATH];

	FILE* output_file;

	// Code

	// Make sure results directory exists
	sprintf_s(file_string, _MAX_PATH, "%s", output_file_string.c_str());
	path output_file_path(file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Results folder could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, file_string, "w");
	if (err != 0)
	{
		printf("Results file: %s\ncould not be opened\n",
			output_file_string.c_str());
		exit(1);
	}

	// Write header
	fprintf_s(output_file, "time%c", delimiter);
	fprintf_s(output_file, "m_length%c", delimiter);
	fprintf_s(output_file, "m_force%c", delimiter);
	fprintf_s(output_file, "sc_extension%c", delimiter);
	fprintf_s(output_file, "sc_force%c", delimiter);

	// Now add in the headers for each half-sarcomere
	for (int hs_counter = 0; hs_counter < p_model->no_of_half_sarcomeres; hs_counter++)
	{
		fprintf_s(output_file, "hs_%i_pCa%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_command_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_slack_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_pas_int_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_viscous_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_pas_ext_force%c", hs_counter + 1, delimiter);

		// Build pops as loops
		for (int j = 0; j < 2; j++)
		{
			fprintf_s(output_file, "hs_%i_a_pop_%i%c", hs_counter + 1, j + 1, delimiter);
		}

		for (int j = 0; j < p_model->p_m_scheme[0]->no_of_states; j++)
		{
			fprintf_s(output_file, "hs_%i_m_pop_%i", hs_counter + 1, j + 1);

			if ((j == (p_model->p_m_scheme[0]->no_of_states - 1)) &&
				(hs_counter == (p_model->no_of_half_sarcomeres - 1)))
			{
				fprintf_s(output_file, "\n");
			}
			else
			{
				fprintf_s(output_file, "%c", delimiter);
			}
		}
	}

	// Loop through points
	for (int i = 0; i < no_of_time_points; i++)
	{
		fprintf_s(output_file, "%g%c%.3f%c%g%c%g%c%g%c",
			gsl_vector_get(time, i), delimiter,
			gsl_vector_get(m_length, i), delimiter,
			gsl_vector_get(m_force, i), delimiter,
			gsl_vector_get(sc_extension, i), delimiter,
			gsl_vector_get(sc_force, i), delimiter);

		// Now add in the data for each  half-sarcomere
		for (int hs_counter = 0; hs_counter < p_model->no_of_half_sarcomeres; hs_counter++)
		{
			fprintf_s(output_file, "%.3g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_pCa, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_command_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_slack_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_pas_int_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_viscous_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_pas_ext_force, i), delimiter);

			// Build pops as loops
			for (int j = 0; j < 2 ; j++)
			{
				fprintf_s(output_file, "%g%c",
					gsl_matrix_get(p_hs_data[hs_counter]->hs_a_pops, i, j), delimiter);
			}

			for (int j = 0; j < p_model->p_m_scheme[0]->no_of_states; j++)
			{
				fprintf_s(output_file, "%g",
					gsl_matrix_get(p_hs_data[hs_counter]->hs_m_pops, i, j));

				if ((j == (p_model->p_m_scheme[0]->no_of_states - 1)) &&
					(hs_counter == (p_model->no_of_half_sarcomeres - 1)))
				{
					fprintf_s(output_file, "\n");
				}
				else
				{
					fprintf_s(output_file, "%c", delimiter);
				}
			}
		}
	}

	// Tidy up
	fclose(output_file);

	std::cout << "Finished writing data to: " << output_file_string << "\n";
}
