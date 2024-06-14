/**
/* @file		half_sarcomere.cpp
/* @brief		Source file for a half_sarcomere object
/* @author		Ken Campbell
*/

#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>

#include "global_definitions.h"
#include "muscle.h"
#include "model.h"
#include "half_sarcomere.h"
#include "myofilaments.h"
#include "membranes.h"

#include "gsl_math.h"
#include "gsl_errno.h"
#include "gsl_math.h"
#include "gsl_roots.h"
#include "gsl_const_num.h"

////#include "kinetic_scheme.h"
//#include "m_state.h"
//#include "transition.h"
//#include "JSON_functions.h"
//#include "kinetic_scheme.h"

// Structure used for root finding for force control mode
struct force_control_params
{
    double target_force;
    double time_step;
    half_sarcomere* p_hs;
};

// Constructor
half_sarcomere::half_sarcomere(muscle* set_p_parent_m, model* set_p_model, int set_hs_index)
{
    // Initialise

    // Code
    std::cout << "In half-sarcomere constructor\n";

    // Set defining variables
    
    p_parent_m = set_p_parent_m;

    p_model = set_p_model;

    hs_index = set_hs_index;

    time_s = 0;

    // Set some values
    hs_length = p_model->initial_hs_length;

    hs_prop_fibrosis = p_model->prop_fibrosis;

    hs_prop_myofilaments = p_model->prop_myofilaments;

    hs_m_head_density = p_model->m_head_density;

    hs_velocity = 0;

    // Create the membranes
    p_memb = new membranes(this);

    // Create the myofilaments
    p_myof = new myofilaments(this);

    p_myof->initialise_simulation();

    // Set some values
    hs_command_length = hs_length;
    hs_slack_length = GSL_NAN;
    hs_force = p_myof->myof_stress_total;
    hs_pas_int_force = p_myof->myof_stress_int_pas;
    hs_pas_ext_force = p_myof->myof_stress_ext_pas;
    hs_viscous_force = p_myof->myof_stress_viscous;


    printf("myof_a_off: %g\n", p_myof->myof_a_off);
}

// Destructor
half_sarcomere::~half_sarcomere(void)
{
    // Code
    std::cout << "In half_sarcomere destructor: " << hs_index << "\n";

    // Tidy up
    delete p_memb;
    delete p_myof;
}

// Other functions
void half_sarcomere::sarcomere_kinetics(double time_step_s, double pCa)
{
    //! Function implements kinetics
    
    // Update Ca
    p_memb->memb_Ca_cytosol = pow(10, -pCa);
    
    p_myof->implement_time_step(time_step_s);

    p_myof->calculate_stresses(true);

    hs_force = p_myof->myof_stress_total;
}

void half_sarcomere::change_hs_length(double delta_hsl, double time_step_s)
{
    //! Changes half-sarcomere length by delta_hsl

    // Code

    // Adjust this half-sarcomere
    hs_length = hs_length + delta_hsl;

    // And the cross-bridges
    p_myof->move_cb_populations(delta_hsl);

    // Now update wall stress
    p_myof->calculate_stresses();

    // Now set half-sarcomere velocity
    hs_velocity = delta_hsl / time_step_s;

    hs_force = p_myof->myof_stress_total;
}

double half_sarcomere::return_stress_after_delta_hsl(double delta_hsl)
{
    //! Function returns wall stress after given delta_hsl

    // Variables
    double stress;

    // Code
    stress = p_myof->return_stress_after_delta_hsl(delta_hsl);

    if (!gsl_finite(stress))
    {
        if (delta_hsl > 0)
            stress = DBL_MAX;
        else
            stress = -DBL_MAX;
    }

    return stress;
}

struct gsl_root_params
{
    half_sarcomere* p_hs_temp;
    double stress_target;
};

double hsl_stress_root_finder(double x, void* params)
{
    //! Function for the root finder

    // Variables
    double new_stress;
    double stress_difference;

    // Code

    // Unpack the pointer
    struct gsl_root_params* p = (struct gsl_root_params*)params;

    // Calculate the new stress after a length change
    new_stress = p->p_hs_temp->return_stress_after_delta_hsl(x);

    // Calculate the difference between the new stress and the target
    stress_difference = new_stress - p->stress_target;

    // Return the difference
    return stress_difference;
}

double half_sarcomere::return_hs_length_for_stress(double target_stress, double time_step_s)
{
    //! Code returns the hs_length at which force is equal to the target

    // Variables
    double x_lo = -500;
    double x_hi = 500;
    double r;

    const gsl_root_fsolver_type* T;
    gsl_root_fsolver* s;

    gsl_function F;
    struct gsl_root_params params = { this, target_stress };

    int status;
    int iter = 0;
    int max_iter = 100;

    double epsabs = 0.01;
    double epsrel = 0.01;

    // Code

    F.function = &hsl_stress_root_finder;
    F.params = &params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(s);

        r = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, epsabs, epsrel);

    } while ((status == GSL_CONTINUE) && (iter < max_iter));

    gsl_root_fsolver_free(s);

    // Calculate the length

    return (hs_length + r);
}
