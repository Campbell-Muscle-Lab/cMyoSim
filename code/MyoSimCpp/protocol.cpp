/**
* @file		protocol.cpp
* @brief	Source file for the protocol class
* @author	ken Campbell
*/

#include <stdio.h>
#include <string>

#include "protocol.h"

#include "gsl_vector.h"

// Constructor
protocol::protocol(std::string protocol_file_string)
{
	// Creates a protocol from the experimental file

	// Variables
	FILE* protocol_file;

	int counter;
	int no_of_columns = 4;

	double temp_float;

	char temp_string[_MAX_PATH];

	char file_string[_MAX_PATH];

	// Code

	// Open file
	sprintf_s(file_string, _MAX_PATH, "%s", protocol_file_string.c_str());
	fopen_s(&protocol_file, file_string, "r");

	if (!protocol_file)
	{
		printf("Error: protocol could not open from file: %s\n", protocol_file_string.c_str());
		exit(1);
	}

	// Scan through file to get number of points
	counter = 0;
	while (!feof(protocol_file))
	{
		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
		counter = counter + 1;
	}

	no_of_time_points = ((counter - 1) / no_of_columns) - 1;

	// Allocate space
	dt = gsl_vector_alloc(no_of_time_points);
	pCa = gsl_vector_alloc(no_of_time_points);
	delta_hsl = gsl_vector_alloc(no_of_time_points);
	sim_mode = gsl_vector_alloc(no_of_time_points);

	// Jump back to begininning of file
	rewind(protocol_file);

	// Skip the header
	for (int i = 0; i < no_of_columns; i++)
	{
		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
	}

	// Now read data
	for (int i = 0; i < no_of_time_points; i++)
	{
		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
		temp_float = atof(temp_string);
		gsl_vector_set(dt, i, temp_float);

		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
		temp_float = atof(temp_string);
		gsl_vector_set(pCa, i, temp_float);

		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
		temp_float = atof(temp_string);
		gsl_vector_set(delta_hsl, i, temp_float);

		fscanf_s(protocol_file, "%s", temp_string, _MAX_PATH);
		temp_float = atof(temp_string);
		gsl_vector_set(sim_mode, i, temp_float);
	}

	// Tidy up
	fclose(protocol_file);
}

// Destructor
protocol::~protocol(void)
{
	// Recover space
	gsl_vector_free(dt);
	gsl_vector_free(pCa);
	gsl_vector_free(delta_hsl);
	gsl_vector_free(sim_mode);
}