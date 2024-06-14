/**
/* @file		series_component.cpp
/* @brief		Source file for a series_component object
/* @author		Ken Campbell
*/

#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>

//#include "global_definitions.h"
#include "muscle.h"
#include "model.h"
#include "series_component.h"

#include "gsl_math.h"

// Constructor
series_component::series_component(muscle* set_p_parent_m)
{
    // Initialise

    // Code
    std::cout << "In series_component constructor\n";

    // Set defining variables
    p_parent_m = set_p_parent_m;

	// Set the stiffness
	sc_k_stiff = p_parent_m->p_model[0]->sc_k_stiff;
	sc_sigma = p_parent_m->p_model[0]->sc_sigma;
	sc_L = p_parent_m->p_model[0]->sc_L;

	printf("sc constructor: k_stiff: %g\n", sc_k_stiff);

	// Set the initial length to 0
	sc_extension = 0.0;

	// And the force to 0
	sc_force = 0.0;
}

// Destructor
series_component::~series_component(void)
{
    // Code
    std::cout << "In series_component destructor\n";

    // Tidy up
}

// Other functions
double series_component::return_series_extension(double muscle_force)
{
	//! Returns the extension of the series component for a given force

	// Variables
	double ext;

	// Code
	if (gsl_isnan(sc_k_stiff))
	{
		printf("Error in series_component::return_series_extension, sc_k_stiff is NAN");
		exit(1);
	}

	if (sc_k_stiff != 0.0)
		ext = muscle_force / sc_k_stiff;
	else
		ext = sc_L * log((muscle_force / sc_sigma) + 1.0);


	//printf("sc_L: %g  sc_sigma: %g  muscle_force: %g  ext: %g\n",
	//	sc_L, sc_sigma, muscle_force, ext);

	return ext;
}

double series_component::return_series_force(double series_extension)
{
	//! Returns the force in the series component for a given extension

	// Variables
	double series_force;

	// Code
	if (sc_k_stiff != 0.0)
		series_force = series_extension * sc_k_stiff;
	else
		series_force = sc_sigma * (exp(series_extension / sc_L) - 1);

	/*printf("sc_L: %g  sc_sigma: %g  muscle_force: %g  ext: %g\n",
		sc_L, sc_sigma, series_force, series_extension);*/


	return series_force;
}