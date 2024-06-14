/**
/* @file		myofilaments.cpp
/* @brief		Source file for a myofilaments object
/* @author		Ken Campbell
*/

#include <iostream>
#include <filesystem>

#include "myofilaments.h"
#include "half_sarcomere.h"
#include "muscle.h"
#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"
#include "membranes.h"
#include "model.h"
#include "options.h"

#include "gsl_math.h"
#include "gsl_errno.h"
#include "gsl_odeiv2.h"
#include "gsl_interp.h"
#include "gsl_spline.h"

// Constructor
myofilaments::myofilaments(half_sarcomere* set_p_parent_hs)
{
	//! Constructor

	// Code

	// Set the pointers to the appropriate places
	p_parent_hs = set_p_parent_hs;

	p_model = p_parent_hs->p_model;
	p_options = p_parent_hs->p_parent_m->p_options;

	// Set pointers safe
	x = NULL;
	y = NULL;
	m_y_indices = NULL;
	m_state_pops = NULL;
	m_state_stresses = NULL;

	// Initialize
	myof_head_density = p_model->m_head_density;
	myof_k_cb = p_model->m_k_cb;

	myof_int_pas_sigma = p_model->int_pas_sigma;
	myof_int_pas_L = p_model->int_pas_L;
	myof_int_pas_slack_hsl = p_model->int_pas_slack_hsl;

	myof_ext_pas_sigma = p_model->ext_pas_sigma;
	myof_ext_pas_L = p_model->ext_pas_L;
	myof_ext_pas_slack_hsl = p_model->ext_pas_slack_hsl;

	myof_fil_compliance_factor = p_model->m_fil_compliance_factor;

	myof_thick_fil_length = p_model->m_thick_fil_length;
	myof_bare_zone_length = p_model->m_bare_zone_length;
	myof_thin_fil_length = p_model->m_thin_fil_length;

	myof_a_k_on = p_model->a_k_on;
	myof_a_k_off = p_model->a_k_off;
	myof_a_k_coop = p_model->a_k_coop;

	myof_a_off = 0.0;
	myof_a_on = 0.0;

	myof_m_bound = 0.0;
	myof_f_overlap = 1.0;

	myof_m_no_of_isotypes = p_model->m_no_of_isotypes;


	//for (int i = 0; i < myof_m_no_of_isotypes; i++)
	//{
	//	printf("i: %i\n", i);
	//	int s = p_model->p_m_scheme[i]->no_of_states;
	//	printf("s: %i\n", s);
	//	m_pops_array[i] = (double*)malloc(p_model->p_m_scheme[i]->no_of_states * sizeof(double));
	//	m_stresses_array[i] = (double*)malloc(p_model->p_m_scheme[i]->no_of_states * sizeof(double));
	//}

	myof_stress_cb = 0.0;
	myof_stress_myof = 0.0;
	myof_stress_int_pas = 0.0;
	myof_stress_ext_pas = 0.0;
	myof_stress_viscous = 0.0;
	myof_stress_total = 0.0;

	no_of_bin_positions = 0;
	y_length = 0;

	max_shift = 0.0;
	n_max_sub_steps = 0;

	cb_dump_file_string;
	cb_dump_file_defined = false;

	myof_ATP_flux = 0.0;
}

// Destructor
myofilaments::~myofilaments(void)
{
	//! Destructor

	// Code

	// Tidy up
	if (x != NULL)
	{
		gsl_vector_free(x);
	}
	if (y != NULL)
	{
		gsl_matrix_free(y);
	}
	if (m_y_indices != NULL)
	{
		gsl_matrix_int_free(m_y_indices);
	}
	if (m_state_pops != NULL)
	{
		gsl_matrix_free(m_state_pops);
	}
	if (m_state_stresses != NULL)
	{
		gsl_matrix_free(m_state_stresses);
	}

	//for (int i = 0; i < myof_m_no_of_isotypes; i++)
	//{
	//	std::free(m_pops_array[i]);
	//		std::free(m_stresses_array[i]);
	//}

	cout << "max_shift: " << max_shift << " n_max_sub_steps: " << n_max_sub_steps << "\n";
}

// Other functions
void myofilaments::initialise_simulation(void)
{
	//! Function adds data fields to main results object

	// Variables

	// Code

	// Now do lots of stuff specific to this class
	// Code
	no_of_bin_positions = 1 + (int)((p_options->bin_x_max - p_options->bin_x_min) /
		p_options->bin_x_width);

	// Assign x
	x = gsl_vector_alloc(no_of_bin_positions);

	for (int i = 0; i < no_of_bin_positions; i++)
	{
		gsl_vector_set(x, i, p_options->bin_x_min + ((double)i * p_options->bin_x_width));
	}

	// Now set the length of the system from the kinetic scheme + 2 for thin filament
	y_length = (size_t)(p_model->p_m_scheme[0]->no_of_detached_states) +
			(size_t)(p_model->p_m_scheme[0]->no_of_attached_states * no_of_bin_positions) +
			2;

	// Allocate space
	y = gsl_matrix_alloc(myof_m_no_of_isotypes, y_length);

	// Set indices
	a_off_index = y_length - 2;
	a_on_index = y_length - 1;

	// Set the calc length
	calc_length = ((size_t)myof_m_no_of_isotypes * (y_length - 2)) + 2;

	// Initialise with all myosins in y[0] and actins in y[end-2]
	gsl_matrix_set_zero(y);
	for (int i = 0; i < myof_m_no_of_isotypes; i++)
	{
		gsl_matrix_set(y, i, 0, gsl_vector_get(p_model->m_isotype_props, i));
		gsl_matrix_set(y, i, a_off_index, 1.0);
	}

	// Initialise and zero the m_bin_indices
	m_y_indices = gsl_matrix_int_alloc(p_model->p_m_scheme[0]->no_of_states, 2);
	gsl_matrix_int_set_zero(m_y_indices);

	int y_index = 0;
	for (int state_counter = 0; state_counter < p_model->p_m_scheme[0]->no_of_states; state_counter++)
	{
		// Check whether state is attached
		if (p_model->p_m_scheme[0]->p_m_states[state_counter]->state_type == 'A')
		{
			gsl_matrix_int_set(m_y_indices, state_counter, 0, y_index);
			y_index = y_index + no_of_bin_positions - 1;
			gsl_matrix_int_set(m_y_indices, state_counter, 1, y_index);
		}
		else
		{
			gsl_matrix_int_set(m_y_indices, state_counter, 0, y_index);
			gsl_matrix_int_set(m_y_indices, state_counter, 1, y_index);
		}

		y_index = y_index + 1;
	}

	cout << "bin indices \n";
	for (int i = 0; i < p_model->p_m_scheme[0]->no_of_states; i++)
	{
		cout << gsl_matrix_int_get(m_y_indices, i, 0) << "   " << gsl_matrix_int_get(m_y_indices, i, 1) << "\n";
	}

	// Initialise and set the bin populations
	m_state_pops = gsl_matrix_alloc(myof_m_no_of_isotypes, p_model->p_m_scheme[0]->no_of_states);
	gsl_matrix_set_zero(m_state_pops);
	for (int i = 0; i < myof_m_no_of_isotypes ; i++)
	{
		gsl_matrix_set(m_state_pops, i, 0, gsl_matrix_get(y, i, 0));
	}

	// And the state forces
	m_state_stresses = gsl_matrix_alloc(myof_m_no_of_isotypes, p_model->p_m_scheme[0]->no_of_states);
	gsl_matrix_set_zero(m_state_stresses);

	// Update class variables
	myof_a_off = gsl_matrix_get(y, 0, a_off_index);
	myof_a_on = gsl_matrix_get(y, 0, a_on_index);

	//// Set the pointer to the results object
	//p_cmv_results_beat = p_parent_hs->p_cmv_results_beat;

	//// Now add the results fields
	//p_cmv_results_beat->add_results_field("myof_a_off", &myof_a_off);
	//p_cmv_results_beat->add_results_field("myof_a_on", &myof_a_on);
	//p_cmv_results_beat->add_results_field("myof_f_overlap", &myof_f_overlap);
	//p_cmv_results_beat->add_results_field("myof_m_bound", &myof_m_bound);
	//p_cmv_results_beat->add_results_field("myof_ATP_flux", &myof_ATP_flux);

	/*for (int i = 0; i < p_m_scheme->no_of_states; i++)
	{
		string label = string("myof_m_pop_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &m_pops_array[i]);
	}

	for (int i = 0; i < p_m_scheme->no_of_states; i++)
	{
		string label = string("myof_m_stress_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &m_stresses_array[i]);
	}

	p_cmv_results_beat->add_results_field("myof_stress_cb", &myof_stress_cb);
	p_cmv_results_beat->add_results_field("myof_stress_int_pas", &myof_stress_int_pas);
	p_cmv_results_beat->add_results_field("myof_stress_ext_pas", &myof_stress_ext_pas);
	p_cmv_results_beat->add_results_field("myof_stress_myof", &myof_stress_myof);
	p_cmv_results_beat->add_results_field("myof_stress_total", &myof_stress_total);*/
}

void myofilaments::implement_time_step(double time_step_s)
{
	//! Function implements a time ste
	
	// Code
	run_kinetics(time_step_s);
}

// This function is not a member of the myofilaments class but is used to interace
// with the GSL ODE system. It must appear before the myofilaments class members
// that calls it, and communicates with the myofilament class through a pointer to
// the class object

int myof_calculate_derivs(double t, const double y[], double f[], void* params)
{
	// Function sets derivs

	// Variables
	(void)(t);

	myofilaments* p_myof = (myofilaments*)params;

	double rate;

	double hs_stress;
	double hs_length;

	double x_pos;
	double x_ext;

	double f_overlap;
	double m_bound;

	double J_on;
	double J_off;

	double flux;

	int current_ind;
	int new_ind;

	// Code

	// Calculat f_overlap
	p_myof->calculate_f_overlap();
	f_overlap = p_myof->myof_f_overlap;

	// Calculate state populations
	p_myof->calculate_m_state_pops(y);
	m_bound = p_myof->myof_m_bound;

	// Set state variables
	hs_stress = p_myof->myof_stress_myof;
	hs_length = p_myof->p_parent_hs->hs_length;

	// Initalise derivs
	for (int i = 0; i < p_myof->calc_length; i++)
	{
		f[i] = 0.0;
	}

	// Zero the flux
	p_myof->myof_ATP_flux = 0.0;

	// Start with myosin

	// Work through the states
	for (int iso_counter = 0; iso_counter < p_myof->myof_m_no_of_isotypes; iso_counter++)
	{
		for (int state_counter = 0; state_counter < p_myof->p_model->p_m_scheme[iso_counter]->no_of_states;
				state_counter++)
		{
			char current_state_type = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->state_type;

			// Now through the transitions
			for (int t_counter = 0; t_counter < p_myof->p_model->p_m_scheme[iso_counter]->max_no_of_transitions;
				t_counter++)
			{
				int new_state = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->
					p_transitions[t_counter]->new_state;

				if (new_state == 0)
				{
					// Transition is not allowed out - skip out
					continue;
				}

				char new_state_type = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[new_state - 1]->state_type;

				if ((current_state_type == 'S') || (current_state_type == 'D'))
				{
					// Deatched state
					if ((new_state_type == 'S') || (new_state_type == 'D'))
					{
						// Detached to detached
						rate = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->p_transitions[t_counter]->
							calculate_rate(0, 0, hs_stress, p_myof->p_parent_hs);

						// Find current index
						current_ind = (iso_counter* (p_myof->y_length - 2)) +
							gsl_matrix_int_get(p_myof->m_y_indices, state_counter, 0);

						// Find flux
						flux = rate * y[current_ind];

						// Cross-bridges leaving current state
						f[current_ind] = f[current_ind] - flux;

						// Cross-bridges arriving at new state
						new_ind = (iso_counter * (p_myof->y_length - 2)) +
							gsl_matrix_int_get(p_myof->m_y_indices, new_state - 1, 0);

						f[new_ind] = f[new_ind] + flux;
					}
					else
					{
						// Detached to attached
						current_ind = (iso_counter * (p_myof->y_length - 2)) + 
							gsl_matrix_int_get(p_myof->m_y_indices, state_counter, 0);

						// Cycle through the bins
						for (int bin_index = 0; bin_index < p_myof->no_of_bin_positions;
							bin_index++)
						{
							// Get x position
							x_pos = gsl_vector_get(p_myof->x, bin_index);
							x_ext = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->extension;

							// Calculate rate
							rate = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->
								p_transitions[t_counter]->calculate_rate(x_pos, x_ext, hs_stress, p_myof->p_parent_hs);

							// Adjust for movement enhancement
							if (!gsl_isnan(p_myof->p_options->movement_enhancement_width))
							{
								rate = rate * (1.0 + fabs(p_myof->p_parent_hs->hs_velocity *
														p_myof->obj_time_step_s /
														p_myof->p_options->movement_enhancement_width));
							}

							flux = p_myof->p_options->bin_x_width *
								rate * y[current_ind] * (y[p_myof->a_on_index] - m_bound);

							new_ind = (iso_counter * (p_myof->y_length - 2)) + 
										gsl_matrix_int_get(p_myof->m_y_indices, new_state - 1, 0) +
										bin_index;

							// Cross-bridges leaving this state
							f[current_ind] = f[current_ind] - flux;

							// Cross-bridges arriving at new state
							f[new_ind] = f[new_ind] + flux;
						}
					}
				}
				else
				{
					// We are in an attached state
					if ((new_state_type == 'S') || (new_state_type == 'D'))
					{
						// Attached to detached
						new_ind = (iso_counter * (p_myof->y_length - 2)) +
							gsl_matrix_int_get(p_myof->m_y_indices, new_state - 1, 0);

						// Cycle through the bins
						for (int bin_index = 0; bin_index < p_myof->no_of_bin_positions;
							bin_index++)
						{
							// Get x position
							x_pos = gsl_vector_get(p_myof->x, bin_index);
							x_ext = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->extension;

							current_ind = (iso_counter * (p_myof->y_length - 2)) +
								gsl_matrix_int_get(p_myof->m_y_indices, state_counter, 0) + bin_index;

							// Calculate rate
							rate = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->
								p_transitions[t_counter]->calculate_rate(x_pos, x_ext, hs_stress, p_myof->p_parent_hs);

							flux = rate * y[current_ind];

							// Cross-bridges leaving this state
							f[current_ind] = f[current_ind] - flux;

							// Cross-bridges arriving at new state
							f[new_ind] = f[new_ind] + flux;

							// Check for ATP
							if (p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->p_transitions[t_counter]->ATP_required == 'y')
							{
								p_myof->myof_ATP_flux = p_myof->myof_ATP_flux + flux;
							}
						}
					}
					else
					{
						// Cross-bridges transitioning between bound states

						// Cycle through the bins
						for (int bin_index = 0; bin_index < p_myof->no_of_bin_positions;
							bin_index++)
						{
							// Get x position
							x_pos = gsl_vector_get(p_myof->x, bin_index);
							x_ext = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->extension;

							current_ind = (iso_counter * (p_myof->y_length - 2)) + 
											gsl_matrix_int_get(p_myof->m_y_indices, state_counter, 0) +
											bin_index;

							// Calculate rate
							rate = p_myof->p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->
								p_transitions[t_counter]->calculate_rate(x_pos, x_ext, hs_stress, p_myof->p_parent_hs);

							new_ind = (iso_counter * (p_myof->y_length - 2)) + 
										gsl_matrix_int_get(p_myof->m_y_indices, new_state - 1, 0) +
										bin_index;

							flux = rate * y[current_ind];

							// Cross-bridges leaving this state
							f[current_ind] = f[current_ind] - flux;

							// Cross-bridgse arriving at new state
							f[new_ind] = f[new_ind] + flux;
						}
					}
				}
			}
		}
	}

	// Now handle the actin
	int a_on_ind = ((p_myof->myof_m_no_of_isotypes - 1) * (p_myof->y_length - 2)) +
						p_myof->a_on_index;
	int a_off_ind = ((p_myof->myof_m_no_of_isotypes - 1) * (p_myof->y_length - 2)) +
						p_myof->a_off_index;

	if (f_overlap > 0.0)
	{
		J_on = p_myof->myof_a_k_on *
			(p_myof->p_parent_hs->p_memb->memb_Ca_cytosol) *
				(f_overlap - y[a_on_ind]) *
			(1.0 + (p_myof->myof_a_k_coop * (y[a_on_ind] / f_overlap)));

		J_off = p_myof->myof_a_k_off *
			(y[a_on_ind] - m_bound) *
			(1.0 + (p_myof->myof_a_k_coop * ((f_overlap - y[a_on_ind]) / f_overlap)));
	}
	else
	{
		J_on = 0.0;
		J_off = p_myof->myof_a_k_off * (y[a_on_ind] - m_bound);
	}

	f[a_off_ind] = -J_on + J_off;
	f[a_on_ind] = -f[a_off_ind];

	return GSL_SUCCESS;
}

void myofilaments::run_kinetics(double time_step_s)
{
	//! Code advances the simulation by time_step

	// Variables
	double eps_abs = 1e-6;
	double eps_rel = 1e-6;

	int status;

	int calc_length;
	int calc_counter;

	double t_start_s = 0.0;
	double t_stop_s = time_step_s;

	double* y_calc;

	double holder;
	double adjustment;

	// Code

	// Allocate memory for y
	// calc_length = y_length - 2 for each isotype + 2 for actin
	calc_length = 2 + myof_m_no_of_isotypes * (y_length - 2);
	
	y_calc = (double*)malloc(calc_length * sizeof(double));
	
	// Fill y_calc
	calc_counter = 0;
	for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
	{
		for (size_t i = 0; i < (y_length-2); i++)
		{
			y_calc[calc_counter] = gsl_matrix_get(y, iso_counter, i);

			calc_counter = calc_counter + 1;
		}
	}
	// Now the thin states
	y_calc[calc_length - 2] = gsl_matrix_get(y, 0, y_length - 2);
	y_calc[calc_length - 1] = gsl_matrix_get(y, 0, y_length - 1);

	//// Check y_calc
	//printf("Before\n");
	//for (int i = 0; i < y_length; i++)
	//{
	//	printf("[%i]: %g\n", i, y_calc[i]);
	//}
	//printf("\n");

	obj_time_step_s = time_step_s;

	gsl_odeiv2_system sys = { myof_calculate_derivs, NULL, calc_length, this };

	gsl_odeiv2_driver* d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
			0.5 * time_step_s, eps_abs, eps_rel);

	status = gsl_odeiv2_driver_apply(d, &t_start_s, t_stop_s, y_calc);

	gsl_odeiv2_driver_free(d);

	if (status != GSL_SUCCESS)
	{
		std::cout << "Integration problem in myofilaments::implement_time_step\n";
		exit(1);
	}
	else
	{
		//Unpack, noting how many bridges we have
		calc_counter = 0;
		for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
		{
			holder = 0.0;

			for (int i = 0; i < (y_length - 2); i++)
			{
				if (y_calc[calc_counter] < 0.0)
					y_calc[calc_counter] = 0.0;

				gsl_matrix_set(y, iso_counter, i, y_calc[calc_counter]);

				holder = holder + y_calc[calc_counter];

				calc_counter = calc_counter + 1;
			}

			// Add m problems back into the first state
			adjustment = gsl_vector_get(p_model->m_isotype_props, iso_counter) - holder;

			gsl_matrix_set(y, iso_counter, 0, gsl_matrix_get(y, iso_counter, 0) + adjustment);

			if (fabs(adjustment) > 0.001)
			{
				double t = p_parent_hs->p_parent_m->cum_time_s;
				cout << "t: " << t << " fast sliding : " << adjustment << "\n";
			}

			// Now update the actins
			gsl_matrix_set(y, iso_counter, y_length - 2, y_calc[calc_length - 2]);
			gsl_matrix_set(y, iso_counter, y_length - 1, y_calc[calc_length - 1]);
		}
	}

	//// Check y_calc
	//printf("After\n");
	//for (int i = 0; i < y_length; i++)
	//{
	//	printf("[%i]: %g\n", i, y_calc[i]);
	//}
	//printf("\n");

	// Update class variables

	// Populations
	calculate_m_state_pops(y_calc);
	/*for (int i = 0; i < p_m_scheme->no_of_states; i++)
	{
		m_pops_array[i] = gsl_vector_get(m_state_pops, i);
	}*/

	myof_a_off = gsl_matrix_get(y, 0, a_off_index);
	myof_a_on = gsl_matrix_get(y, 0, a_on_index);

	//// Forces
	//calculate_m_state_stresses();
	//for (int i = 0; i < p_m_scheme->no_of_states; i++)
	//{
	//	m_stresses_array[i] = gsl_vector_get(m_state_stresses, i);
	//}

	//calculate_stresses();

	// Tidy up
	free(y_calc);
}

void myofilaments::calculate_m_state_pops(const double y_calc[])
{
	//! Function returns m_bound

	// Variables
	double holder;
	double bound_holder;

	// Code

	bound_holder = 0.0;

	for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
	{
		for (int state_counter = 0; state_counter < p_model->p_m_scheme[iso_counter]->no_of_states;
			state_counter++)
		{
			holder = 0.0;

			int state_start = (iso_counter * (y_length - 2)) +
				gsl_matrix_int_get(m_y_indices, state_counter, 0);
			int state_stop = (iso_counter * (y_length - 2)) +
				gsl_matrix_int_get(m_y_indices, state_counter, 1);

			for (int i = state_start ; i <= state_stop ; i++)
			{
				holder = holder + y_calc[i];
			}
			gsl_matrix_set(m_state_pops, iso_counter, state_counter, holder);

			if (p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->state_type == 'A')
			{
				bound_holder = bound_holder + holder;
			}
		}
	}

	myof_m_bound = bound_holder;
}

void myofilaments::calculate_f_overlap(void)
{
	//! Calculate f_overlap
	//! 
	//!     <-- Thin--->
	//!    |------------         |
	//!    |                     |
	//!    |      --------|------|
	//!           <--Thick------>|
	//!                    <Bare>|

	// Variables
	double x_no_overlap;
	double x_overlap;
	double max_x_overlap;
	double protrusion;

	// Code

	x_no_overlap = p_parent_hs->hs_length - myof_thick_fil_length;
	x_overlap = myof_thin_fil_length - x_no_overlap;
	max_x_overlap = myof_thick_fil_length - myof_bare_zone_length;
	protrusion = myof_thin_fil_length - p_parent_hs->hs_length;

	if (x_overlap == 0.0)
		myof_f_overlap = 0.0;

	if ((x_overlap > 0.0) && (x_overlap <= max_x_overlap))
		myof_f_overlap = x_overlap / max_x_overlap;

	if (x_overlap > max_x_overlap)
		myof_f_overlap = 1.0;

	if (protrusion > 0.0)
	{
		x_overlap = max_x_overlap - protrusion;
		myof_f_overlap = x_overlap / max_x_overlap;
	}
}

void myofilaments::calculate_m_state_stresses(void)
{
	//| Function calculates cross-bridge stresses

	// Variables

	int bin_index;

	double x_bin;
	double x_ext;
	double bin_pop;

	double holder;

	// Code

	// Cycle through through the states
	for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
	{
		for (int state_counter = 0; state_counter < p_model->p_m_scheme[iso_counter]->no_of_states;
			state_counter++)
		{
			char state_type = p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->state_type;

			holder = 0.0;

			if (state_type == 'A')
			{
				x_ext = p_model->p_m_scheme[iso_counter]->p_m_states[state_counter]->extension;

				// Set the bin index
				bin_index = gsl_matrix_int_get(m_y_indices, state_counter, 0);

				for (int i = 0; i < no_of_bin_positions; i++)
				{
					// Get bin position
					x_bin = gsl_vector_get(x, i);

					/// Get population
					bin_pop = gsl_matrix_get(y, iso_counter, bin_index);

					// Add force to holder
					holder = holder + (bin_pop * (x_bin + x_ext));

					bin_index = bin_index + 1;
				}
			}

			// Adjust for fibrosis, myofilament area, and units
			holder = (1.0 - p_parent_hs->hs_prop_fibrosis) *
				p_parent_hs->hs_prop_myofilaments *
				p_parent_hs->hs_m_head_density * 1e-9 * myof_k_cb * holder;

			gsl_matrix_set(m_state_stresses, iso_counter, state_counter, holder);
		}
	}
}

void myofilaments::calculate_stresses(bool check_only)
{
	//! Code calculates forces

	// Code
	calculate_cb_stress();
	calculate_int_pas_stress();
	calculate_ext_pas_stress();

	myof_stress_myof = myof_stress_cb + myof_stress_int_pas;
	myof_stress_total = myof_stress_myof + myof_stress_ext_pas;
}

double myofilaments::calculate_cb_stress(bool check_only)
{
	//! Code calculates cb force

	// Variables
	double holder = 0.0;

	// Code

	// First update the matrix
	calculate_m_state_stresses();

	// Then sum it up
	for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
	{
		for (int i = 0; i < p_model->p_m_scheme[iso_counter]->no_of_states; i++)
		{
			holder = holder + gsl_matrix_get(m_state_stresses, iso_counter, i);
		}
	}

	if (check_only == false)
		myof_stress_cb = holder;

	return holder;
}

double myofilaments::calculate_int_pas_stress(bool check_only, double delta_hsl)
{
	//! Code calculates internal passive force for a given delta_hsl
	//! If check_only is true, does not set class variable

	// Variables

	double x;			// length relative to slack in nm

	double pas_stress;	// intracellular passive stress

	// Code

	x = (p_parent_hs->hs_length + delta_hsl) - myof_int_pas_slack_hsl;

	if (x > 0.0)
	{
		pas_stress = myof_int_pas_sigma * (exp(x / myof_int_pas_L) - 1.0);
	}
	else
	{
		pas_stress = -myof_int_pas_sigma * (exp(fabs(x) / myof_int_pas_L) - 1.0);
	}

	if (check_only == false)
		myof_stress_int_pas = pas_stress;

	return pas_stress;
}

double myofilaments::calculate_ext_pas_stress(bool check_only, double delta_hsl)
{
	//! Code calculates external passive force for a given delta_hsl
	//! If check_only is true, does not set class variable

	// Variables

	double x;			// length relative to slack in nm

	double pas_stress;	// intracellular passive stress

	// Code

	x = (p_parent_hs->hs_length + delta_hsl) - myof_ext_pas_slack_hsl;

	if (x > 0.0)
	{
		pas_stress = myof_ext_pas_sigma * (exp(x / myof_ext_pas_L) - 1.0);
	}
	else
	{
		pas_stress = -myof_ext_pas_sigma * (exp(fabs(x) / myof_ext_pas_L) - 1.0);
	}

	if (check_only == false)
		myof_stress_ext_pas = pas_stress;

	return pas_stress;
}

double myofilaments::return_stress_after_delta_hsl(double delta_hsl)
{
	//! Function returns stress after given delta_hsl

	// Variables
	double delta_int_pas_stress;
	double delta_ext_pas_stress;
	double delta_cb_stress;
	double new_stress;

	// Code
	delta_int_pas_stress = calculate_int_pas_stress(true, delta_hsl) -
		myof_stress_int_pas;

	delta_ext_pas_stress = calculate_ext_pas_stress(true, delta_hsl) -
		myof_stress_ext_pas;

	delta_cb_stress = (1.0 - p_parent_hs->hs_prop_fibrosis) *
		p_parent_hs->hs_prop_myofilaments *
		p_parent_hs->hs_m_head_density * 1e-9 * myof_k_cb * myof_m_bound *
		myof_fil_compliance_factor * delta_hsl;

	new_stress = myof_stress_total +
		delta_int_pas_stress +
		delta_ext_pas_stress +
		delta_cb_stress;

	return new_stress;
}

void myofilaments::move_cb_populations(double delta_hsl)
{
	//! Code displaces cross-bridge populations

	// Variables
	bool x_defined = false;

	double* x_calc = NULL;
	double* y_calc = NULL;
	double* y_temp;

	double x_shift;

	int n_sub_steps;
	double s;
	bool keep_going;

	const gsl_interp_type* interp_type = gsl_interp_linear;
	gsl_interp_accel* acc = NULL;
	gsl_spline* spline = NULL;

	// Code

	// Skip out if delta_hsl == 0
	if (delta_hsl == 0.0)
	{
		return;
	}

	// Allocate memory for y_calc
	x_calc = (double*)malloc(no_of_bin_positions * sizeof(double));
	y_calc = (double*)malloc(no_of_bin_positions * sizeof(double));
	y_temp = (double*)malloc(no_of_bin_positions * sizeof(double));

	// Cycle through states
	for (int iso_counter = 0; iso_counter < myof_m_no_of_isotypes; iso_counter++)
	{
		for (size_t state_counter = 0; state_counter < (size_t)(p_model->p_m_scheme[0]->no_of_states); state_counter++)
		{
			if (p_model->p_m_scheme[0]->p_m_states[state_counter]->state_type == 'A')
			{
				// We need to move populations

				// Set x_calc
				if (x_defined == false)
				{
					// First time
					// Set x
					for (size_t ind = 0; ind < no_of_bin_positions; ind++)
					{
						x_calc[ind] = gsl_vector_get(x, ind);
					}

					// Initialise for interpolation
					acc = gsl_interp_accel_alloc();
					spline = gsl_spline_alloc(interp_type, no_of_bin_positions);

					// Work out the shift
					x_shift = myof_fil_compliance_factor * delta_hsl;

					// Note we are set up
					x_defined = true;
				}

				// Subdivide if necessary
				n_sub_steps = 1;
				s = x_shift;
				keep_going = true;

				while (keep_going)
				{
					if (fabs(s) < 1.0)
					{
						keep_going = false;
					}
					else
					{
						n_sub_steps = n_sub_steps + 1;
						s = x_shift / (double)(n_sub_steps);
					}
				}

				for (int repeat = 1; repeat <= n_sub_steps; repeat++)
				{
					if ((repeat == 1) && (n_sub_steps > n_max_sub_steps))
					{
						cout << "n_sub_steps: " << n_sub_steps << "  x_shift: " << x_shift << "   s: " << s << "\n";
					}
					for (size_t ind = 0; ind < no_of_bin_positions; ind++)
					{
						if (repeat == 1)
						{
							y_calc[ind] = gsl_matrix_get(y, iso_counter, ind + gsl_matrix_int_get(m_y_indices, state_counter, 0));
						}
						else
						{
							y_calc[ind] = y_temp[ind];
						}
					}

					gsl_spline_init(spline, x_calc, y_calc, no_of_bin_positions);

					// Calculate at the new positions
					for (size_t ind = 0; ind < no_of_bin_positions; ind++)
					{
						double new_pos = x_calc[ind] - s;

						if ((new_pos < p_options->bin_x_min) || (new_pos > p_options->bin_x_max))
							y_temp[ind] = 0.0;
						else
						{
							y_temp[ind] = gsl_spline_eval(spline, new_pos, acc);
						}
					}
				}

				for (size_t ind = 0; ind < no_of_bin_positions; ind++)
				{
					// Assign 
					gsl_matrix_set(y, iso_counter, gsl_matrix_int_get(m_y_indices, state_counter, 0) + ind,
						y_temp[ind]);
				}
			}
		}
	}

	if ((fabs(x_shift) > fabs(max_shift)) && (p_parent_hs->time_s > 0.001))
	{
		max_shift = x_shift;
		n_max_sub_steps = n_sub_steps;
	}

	// Tidy up
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	free(x_calc);
	free(y_calc);
	free(y_temp);

}
//
//void myofilaments::dump_cb_distributions(void)
//{
//	//! Code dumps cb distributions
//
//	// Variables
//
//	FILE* output_file;
//
//	path options_file_path;		// path for options file
//	path output_file_path;		// path for output file
//
//	string open_mode;
//
//	int attached_state_counter;
//
//	// Code
//
//	// Adjust the rates file
//	if (cb_dump_file_defined == false)
//	{
//		// First time through, set the class member value
//		path base_dir;
//
//		if (p_cmv_options->cb_dump_relative_to == "this_file")
//		{
//			options_file_path = path(p_cmv_options->options_file_string);
//			base_dir = options_file_path.parent_path();
//		}
//		else
//		{
//			base_dir = path(p_cmv_options->rates_dump_relative_to);
//		}
//
//		output_file_path = base_dir / p_cmv_options->cb_dump_file_string;
//		cb_dump_file_string = output_file_path.string();
//	}
//
//	// Make sure directory exists
//	output_file_path = absolute(path(cb_dump_file_string));
//
//	if (!(is_directory(output_file_path.parent_path())))
//	{
//		if (create_directories(output_file_path.parent_path()))
//		{
//			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
//		}
//		else
//		{
//			std::cout << "\nError: Folder for cb dump file could not be created: " <<
//				output_file_path.parent_path().string() << "\n";
//			exit(1);
//		}
//	}
//
//	if (cb_dump_file_defined == false)
//		open_mode = "w";
//	else
//		open_mode = "a";
//
//	// Check file can be opened, abort if not
//	errno_t err = fopen_s(&output_file, cb_dump_file_string.c_str(), open_mode.c_str());
//
//	if (err != 0)
//	{
//		printf("dump_cb_distributions_to_file(): %s\ncould not be opened\n",
//			cb_dump_file_string.c_str());
//		exit(1);
//	}
//	else
//	{
//		if (cb_dump_file_defined == false)
//			cout << "Writing cb distributions to: " << cb_dump_file_string << "\n";
//	}
//
//	if (cb_dump_file_defined == false)
//	{
//		// Write the header
//		attached_state_counter = 0;
//		for (int state_counter = 0; state_counter < p_m_scheme->no_of_states; state_counter++)
//		{
//			if (p_m_scheme->p_m_states[state_counter]->state_type == 'A')
//			{
//				attached_state_counter = attached_state_counter + 1;
//
//				for (int ind = 0; ind < no_of_bin_positions; ind++)
//				{
//					fprintf(output_file, "M%i_%.0f", attached_state_counter, 10 * gsl_vector_get(x, ind));
//					if ((ind == (no_of_bin_positions - 1)) &&
//						(attached_state_counter == p_m_scheme->no_of_attached_states))
//					{
//						fprintf(output_file, "\n");
//					}
//					else
//					{
//						fprintf(output_file, "\t");
//					}
//				}
//			}
//		}
//	}
//
//	// Now write data
//	// Write the header
//	attached_state_counter = 0;
//	for (int state_counter = 0; state_counter < p_m_scheme->no_of_states; state_counter++)
//	{
//		if (p_m_scheme->p_m_states[state_counter]->state_type == 'A')
//		{
//			attached_state_counter = attached_state_counter + 1;
//
//			for (int ind = 0; ind < no_of_bin_positions; ind++)
//			{
//				fprintf(output_file, "%.5f", gsl_vector_get(y,
//					(gsl_matrix_int_get(m_y_indices, state_counter, 0) + ind)));
//
//				if ((ind == (no_of_bin_positions - 1)) &&
//					(attached_state_counter == p_m_scheme->no_of_attached_states))
//				{
//					fprintf(output_file, "\n");
//				}
//				else
//				{
//					fprintf(output_file, "\t");
//				}
//			}
//		}
//	}
//
//	// Tidy up
//	fclose(output_file);
//
//	// Update
//	cb_dump_file_defined = true;
//}
