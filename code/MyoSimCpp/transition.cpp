/**
* @file		transition.cpp
* @brief	Source file for the transition class
* @author	Ken Campbell
*/

#include <cstdio>
#include <math.h>

#include "transition.h"
#include "m_state.h"
#include "kinetic_scheme.h"
#include "myofilaments.h"
#include "half_sarcomere.h"
#include "muscle.h"
#include "global_definitions.h"
#include "JSON_functions.h"

#include "model.h"
#include "options.h"

#include "rapidjson\document.h"

#include "gsl_vector.h"
#include "gsl_math.h"
#include "gsl_const_mksa.h"


// Constructor
transition::transition(const rapidjson::Value& tr, m_state* set_p_parent_m_state)
{
	// Set p_parent_m_state
	p_parent_m_state = set_p_parent_m_state;

	// Set transition_type to unknown - will be set later on 
	transition_type = 'x';

	JSON_functions::check_JSON_member_int(tr, "new_state");
	new_state = tr["new_state"].GetInt();

	JSON_functions::check_JSON_member_string(tr, "rate_type");
	sprintf_s(rate_type, _MAX_PATH, tr["rate_type"].GetString());

	if (JSON_functions::check_JSON_member_exists(tr, "ATP_required"))
	{
		JSON_functions::check_JSON_member_string(tr, "ATP_required");
		char temp_string[_MAX_PATH];
		sprintf_s(temp_string, _MAX_PATH, tr["ATP_required"].GetString());
		ATP_required = (char)tolower(temp_string[0]);
	}

	// Read in parameters
	JSON_functions::check_JSON_member_array(tr, "rate_parameters");
	const rapidjson::Value& rp = tr["rate_parameters"];

	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);

	for (int i = 0; i < (int)rp.Size(); i++)
	{
		gsl_vector_set(rate_parameters, i, rp[i].GetDouble());
	}
}

transition::transition()
{
	// Default constructor - used if there is no defined transition
	new_state = 0;
	transition_type = 'x';
	sprintf_s(rate_type, _MAX_PATH, "");
	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);
}

// Destructor
transition::~transition(void)
{
	// Tidy up
	gsl_vector_free(rate_parameters);
}

// Functions

double transition::calculate_rate(double x, double x_ext, double force, half_sarcomere* p_hs)
{
	//! Returns the rate for a transition with a given x
	//! 

	// Variables
	double rate = 0.0;

	options* p_options = p_parent_m_state->p_parent_scheme->p_options;

	// Code

	// Constant
	if (!strcmp(rate_type, "constant"))
	{
		rate = gsl_vector_get(rate_parameters, 0);
	}

	// Force-dependent
	if (!strcmp(rate_type, "force_dependent"))
	{
		rate = gsl_vector_get(rate_parameters, 0) *
			(1.0 + (gsl_max(force, 0.0) * gsl_vector_get(rate_parameters, 1)));
	}

	//// Force and adjacent hs dependent
	//if (!strcmp(rate_type, "force_and_adjacent_hs_dependent"))
	//{
	//	// Assume that first half-sarcomere starts with it's Z-line on the left hand side

	//	int hs_ind;
	//	int hs_across_Z;
	//	int hs_across_M;

	//	half_sarcomere* p_hs_across_Z = NULL;
	//	half_sarcomere* p_hs_across_M = NULL;

	//	muscle* p_parent_m;

	//	model* p_model;

	//	double factor_across_Z = 0;
	//	double factor_across_M = 0;

	//	// Calculate node_force factor
	//	double force_factor = gsl_max(force, 0.0) * gsl_vector_get(rate_parameters, 1);

	//	// Set the pointers
	//	hs_ind = p_hs->hs_index;
	//	p_parent_m = p_hs->p_parent_m;
	//	p_model = p_hs->p_model;

	//	// Need to decide if this is a half-sarcomere with Z-disk or M-line on left
	//	if (GSL_IS_EVEN(hs_ind))
	//	{
	//		// Half-sarcomere has Z-disk on low (left) side
	//		if (hs_ind > 1)
	//		{
	//			hs_across_Z = hs_ind - 1;
	//			p_hs_across_Z = p_parent_m->p_hs[hs_across_Z];
	//			factor_across_Z = gsl_vector_get(rate_parameters, 2) *
	//				(p_hs->hs_titin_force - p_hs_across_Z->hs_titin_force);
	//		}

	//		if (hs_ind < (p_fs_model->no_of_half_sarcomeres - 1))
	//		{
	//			hs_across_M = hs_ind + 1;
	//			p_hs_across_M = p_parent_m->p_hs[hs_across_M];
	//			factor_across_M = gsl_vector_get(rate_parameters, 3) *
	//				(p_hs->hs_titin_force - p_hs_across_M->hs_titin_force);
	//		}
	//	}
	//	else
	//	{
	//		// Half-sarcomere has M-line on low (left) side
	//		if (hs_ind > 0)
	//		{
	//			hs_across_M = hs_ind - 1;
	//			p_hs_across_M = p_parent_m->p_hs[hs_across_M];
	//			factor_across_M = gsl_vector_get(rate_parameters, 3) *
	//				(p_hs->hs_titin_force - p_hs_across_M->hs_titin_force);
	//		}

	//		if (hs_ind < (p_fs_model->no_of_half_sarcomeres - 1))
	//		{
	//			hs_across_Z = hs_ind + 1;
	//			p_hs_across_Z = p_parent_m->p_hs[hs_across_Z];
	//			factor_across_Z = gsl_vector_get(rate_parameters, 2) *
	//				(p_hs->hs_titin_force - p_hs_across_Z->hs_titin_force);
	//		}
	//	}

	//	rate = gsl_vector_get(rate_parameters, 0) *
	//		(1.0 + (node_force_factor + gsl_max(factor_across_Z, 0.0) + gsl_max(factor_across_M, 0.0)));

	//	//printf("hs_index: %i  titin_f: %g\t\tnode_f: %g\t\tf_across_Z: %g\t\tf_across_M: %g\n",
	//		//hs_ind, p_hs->hs_titin_force, node_force_factor, factor_across_Z, factor_across_M);
	//}

	// Gaussian
	if (!strcmp(rate_type, "gaussian"))
	{
		model* p_model = p_parent_m_state->p_parent_scheme->p_model;
		double k_cb = p_model->m_k_cb;

		// Set parameters
		double amp = gsl_vector_get(rate_parameters, 0);
		double k_modifier = 1.0;

		if (!gsl_isnan(gsl_vector_get(rate_parameters, 1)))
			k_modifier = gsl_vector_get(rate_parameters, 1);

		// Calculate
		rate = amp *
			exp(-(0.5 * k_modifier * k_cb * gsl_pow_int(x, 2)) /
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));
	}

	// Gaussian_hsl
	if (!strcmp(rate_type, "gaussian_hsl"))
	{
		// Distance between surface of thick and thin filaments is
		// (2/3)*d_1,0 - r_thin - r_thick
		// Assume d_1_0 at hsl = 1100 nm is 37 nm, r_thin = 5.5 nm, t_thick = 7.5 nm
		// d at hsl = x is (2/3) * (37 / sqrt(x/1100)) - 5/5 - 7.5
		// first passage time to position y is t = y^2 / (2*D)
		// rate is proportional to 1/t
		// rate at hsl == x is ref_rate * (y_ref / y_x)^2
		// See PMID 35450825 and first passage in
		// Mechanics of motor proteins and the cytoskeleton, Joe Howard book

		model* p_model = p_parent_m_state->p_parent_scheme->p_model;
		double k_cb = p_model->m_k_cb;

		double hs_length;
		double y_ref;		// distance between filaments at 1100 nm
		double y_actual;	// distance between filaments at current hsl
		double r_thick = 7.5;
		double r_thin = 5.5;

		// Set parameters
		double amp = gsl_vector_get(rate_parameters, 0);
		double k_modifier = 1.0;

		if (!gsl_isnan(gsl_vector_get(rate_parameters, 1)))
			k_modifier = gsl_vector_get(rate_parameters, 1);

		// Deduce filament separation
		y_ref = ((2.0 / 3.0) * 37.0) - r_thick - r_thin;

		if (p_hs == NULL)
			hs_length = 1100.0;
		else
			hs_length = p_hs->hs_length;

		y_actual = (2.0 / 3.0) * (37.0 / sqrt(hs_length / 1100.0)) - r_thick - r_thin;

		rate = amp *
			exp(-(0.5 * k_modifier * k_cb * gsl_pow_int(x, 2)) /
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));

		rate = rate * gsl_pow_2(y_ref / y_actual);
	}
	// Force-dependent Gaussian

	if (!strcmp(rate_type, "force_dependent_gaussian"))
	{

		model* p_model = p_parent_m_state->p_parent_scheme->p_model;
		double k_cb = p_model->m_k_cb;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
					(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature)) *
			(1.0 + gsl_max(force, 0.0) * gsl_vector_get(rate_parameters, 1));
	}

	// Poly
	if (!strcmp(rate_type, "poly"))
	{
		double x_center = gsl_vector_get(rate_parameters, 3); // optional parameter defining the zero of the polynomial

		if (gsl_isnan(x_center)) { // optional parameter is not specified, use the state extension instead
			x_center = x_ext;
		}	

		rate = gsl_vector_get(rate_parameters, 0) +
				(gsl_vector_get(rate_parameters, 1) *
					gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 2)));

	}

	// Poly_asymmetric
	if (!strcmp(rate_type, "poly_asym"))
	{
		double x_center = gsl_vector_get(rate_parameters, 5); // optional parameter defining the zero of the polynomial

		if (gsl_isnan(x_center)) { // optional parameter is not specified, use the state extension instead
			x_center = x_ext;
		}

		if (x > -x_center)
			rate = gsl_vector_get(rate_parameters, 0) +
				(gsl_vector_get(rate_parameters, 1) *
					gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 3)));
		else
			rate = gsl_vector_get(rate_parameters, 0) +
			(gsl_vector_get(rate_parameters, 2) *
				gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 4)));
	}

	// Decreasing load-dependent exponential
	if (!strcmp(rate_type, "exp"))
	{
		double A = gsl_vector_get(rate_parameters, 0);
		double B = gsl_vector_get(rate_parameters, 1);
		double C = gsl_vector_get(rate_parameters, 2);
		double x_center = gsl_vector_get(rate_parameters, 3);
		double x_wall = gsl_vector_get(rate_parameters, 4);

		if (x < x_wall)
			rate = A + B * exp(-C * (x + x_center));
		else
			rate = p_options->max_rate;
	}

	if (!strcmp(rate_type, "sigmoid"))
	{
		// Variables
		double y_left = gsl_vector_get(rate_parameters, 0);
		double y_amp = gsl_vector_get(rate_parameters, 1);
		double k = gsl_vector_get(rate_parameters, 2);
		double x_mid = gsl_vector_get(rate_parameters, 3);

		if (gsl_isnan(x_mid))
		{
			// optional parameter is not specified, use the state extension instead
			model* p_model = p_parent_m_state->p_parent_scheme->p_model;
			kinetic_scheme* p_scheme = p_parent_m_state->p_parent_scheme;
			m_state* p_new_state = p_scheme->p_m_states[new_state - 1];

			double x_new_ext = p_new_state->extension;
			double x_current_ext = p_parent_m_state->extension;
			double x_mid = 0.5 * (x_new_ext + x_current_ext);
		}

		rate = y_left + y_amp * (1 /
			(1 + exp(-k * (x + x_mid))));
	}

	if (!strcmp(rate_type, "exp_wall"))
	{
		// Variables

		model* p_model = p_parent_m_state->p_parent_scheme->p_model;
		double k0 = gsl_vector_get(rate_parameters, 0);
		double F = p_model->m_k_cb * (x + x_ext);
		double d = gsl_vector_get(rate_parameters, 1);
		double x_wall = gsl_vector_get(rate_parameters, 2);
		double x_smooth = gsl_vector_get(rate_parameters, 3);

		// Code
		rate = k0 * exp(-(F * d) /
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));

		rate = rate + p_options->max_rate * (1 /
			(1 + exp(-x_smooth * (x - x_wall))));
	}

	if (!strcmp(rate_type, "bi_wall"))
	{
		// Variables
		double base = gsl_vector_get(rate_parameters, 0);
		double x_pos = gsl_vector_get(rate_parameters, 1);
		double k_pos = gsl_vector_get(rate_parameters, 2);
		double x_neg = gsl_vector_get(rate_parameters, 3);
		double k_neg = gsl_vector_get(rate_parameters, 4);

		rate = base +
			p_options->max_rate *
			((1 / (1 + exp(-k_pos * (x - x_pos)))) +
				(1 / (1 + exp(k_neg * (x - x_neg)))));

	}

	// Curtail at max rate

	if (rate > (p_options->max_rate))
		rate = p_options->max_rate;

	if (rate < 0.0)
		rate = 0.0;

	if (gsl_isnan(rate))
		rate = 0.0;
	
	// Return
	return rate;
}

