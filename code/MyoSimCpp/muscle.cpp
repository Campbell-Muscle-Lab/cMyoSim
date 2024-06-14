/**
/* @file		muscle.cpp
/* @brief		Source file for a muscle object
/* @author		Ken Campbell
*/

#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>

#include "muscle.h"
#include "model.h"
#include "options.h"
#include "half_sarcomere.h"
#include "myofilaments.h"
#include "membranes.h"
#include "series_component.h"
#include "protocol.h"
#include "MyoSim_data.h"
#include "hs_data.h"
#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"

#include "gsl_math.h"
#include "gsl_vector.h"
#include "gsl_multiroots.h"

// Structure used for root-finding for myofibril in length-or force control mode
struct m_control_params
{
    double time_step;
    muscle* p_m;
    double target_force;
};

// This is a function used by the root finding algorithm that handles the recasting of pointers
int wrapper_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* params, gsl_vector* f);

struct force_control_params
{
    double target_force;
    double time_step;
    half_sarcomere* p_hs;
    double delta_hsl;
};

// Constructor
muscle::muscle(std::string JSON_model_file_string, std::string JSON_options_file_string)
{
    // Initialise

    // Code
    std::cout << "In muscle constructor\n";

    // Create a new options objects
    p_options = new options(JSON_options_file_string);

    // Create a new model object - we need this to deduce the number of half-sarcomeres
    p_model[0] = new model(JSON_model_file_string, p_options);

    // Set some safe pointers
    p_protocol = NULL;

    p_data = NULL;

    // Dump rate_functions to file
    if (strlen(p_options->rate_file_string) > 0)
        write_rates_file();

    // Create the half-sarcomeres, building up m_length as we go

    m_length = 0.0;

    for (int i = 0; i < p_model[0]->no_of_half_sarcomeres; i++)
    {
        // Load another model if required
        if (i > 0)
        {
            p_model[i] = new model(JSON_model_file_string, p_options);
        }

        p_hs[i] = new half_sarcomere(this, p_model[i], i);

        m_length = m_length + p_hs[i]->hs_length;
    }

    // Check if we need to make a series component
    if (p_model[0]->sc_exists == false)
    {
        p_sc = NULL;
    }
    else
    {
        p_sc = new series_component(this);

        m_length = m_length + p_sc->sc_extension;
    }
}

// Destructor
muscle::~muscle(void)
{
    // Code
    std::cout << "In muscle destructor\n";

    // Tidy up

    // Delete the series component if it is there
    if (p_sc != NULL)
        delete p_sc;

    // Delete the half-sarcomeres
    int no_of_half_sarcomeres = p_model[0]->no_of_half_sarcomeres;
    for (int i = 0; i < no_of_half_sarcomeres; i++)
    {
        delete p_hs[i];
        delete p_model[i];
    }

    delete p_options;

    delete p_protocol;

    delete p_data;
}

// Other functions
void muscle::implement_protocol(std::string protocol_file_string, std::string results_file_string)
{
    //! Function runs a simulation
    
    // Variables
    

    // Code

    // Make a protocol
    p_protocol = new protocol(protocol_file_string);

    // Make a data object
    p_data = new MyoSim_data(p_protocol->no_of_time_points, p_options, p_model[0]);

    // Reset the time
    cum_time_s = 0;

    // Now implement the protocol
    for (int i = 0 ; i < p_protocol->no_of_time_points ; i++)
    {
        implement_time_step(i);
    }

    // Output main results file
    cout << "Muscle: Attempting to write results to: " << results_file_string << "\n";
    p_data->write_data_to_delimited_file(results_file_string);
}

void muscle::implement_time_step(int protocol_index)
{
    //! Code runs a single time_step
    
    // Variables
    double time_step_s;
    double pCa;
    double sim_mode;
    double delta_hsl;

    double new_length;
    double adjustment;

    int hs_counter;

    // Code

    // Pull values from the protocol

    //// Update the hs_command_length
    //p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_command_length +
    //    gsl_vector_get(p_protocol->delta_hsl, protocol_index);

    // Pull values from protocol
    time_step_s = gsl_vector_get(p_protocol->dt, protocol_index);
    pCa = gsl_vector_get(p_protocol->pCa, protocol_index);
    sim_mode = gsl_vector_get(p_protocol->sim_mode, protocol_index);
    delta_hsl = gsl_vector_get(p_protocol->delta_hsl, protocol_index);

        // Update time
    cum_time_s = cum_time_s + time_step_s;

    // Check for simple model
    if ((p_model[0]->no_of_half_sarcomeres == 1) && (p_sc == NULL))
    {
        // single half sarcomere with no series compliance
        implement_single_half_sarcomere(time_step_s, pCa, sim_mode, delta_hsl);
    }
    else
    {
        // We have a myofibril
        if (sim_mode < 0)
        {
            length_control_myofibril_with_series_compliance(time_step_s, pCa, delta_hsl);
        }
    }     

    // Update the muscle length
    m_length = return_muscle_length();

    // Update the muscle force
    m_force = p_hs[0]->hs_force;

    // Update data
    gsl_vector_set(p_data->time, protocol_index, cum_time_s);
    gsl_vector_set(p_data->m_length, protocol_index, m_length);
    gsl_vector_set(p_data->m_force, protocol_index, m_force);

    if (p_sc != NULL)
    {
        gsl_vector_set(p_data->sc_extension, protocol_index, p_sc->sc_extension);
        gsl_vector_set(p_data->sc_force, protocol_index, p_sc->sc_force);
    }

    // Now the half-sarcomeres
    for (int hs_counter = 0; hs_counter < p_model[0]->no_of_half_sarcomeres; hs_counter++)
    {
        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_pCa,
            protocol_index,
            p_hs[hs_counter]->pCa);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_length,
            protocol_index,
            p_hs[hs_counter]->hs_length);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_command_length,
            protocol_index,
            p_hs[hs_counter]->hs_command_length);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_slack_length,
            protocol_index,
            p_hs[hs_counter]->hs_slack_length);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_force,
            protocol_index,
            p_hs[hs_counter]->hs_force);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_pas_int_force,
            protocol_index,
            p_hs[hs_counter]->hs_pas_int_force);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_viscous_force,
            protocol_index,
            p_hs[hs_counter]->hs_viscous_force);

        gsl_vector_set(p_data->p_hs_data[hs_counter]->hs_pas_ext_force,
            protocol_index,
            p_hs[hs_counter]->hs_pas_ext_force);

        // Update pops
        
        // Actin
        gsl_matrix_set(p_data->p_hs_data[hs_counter]->hs_a_pops,
            protocol_index, 0, p_hs[hs_counter]->p_myof->myof_a_off);
        gsl_matrix_set(p_data->p_hs_data[hs_counter]->hs_a_pops,
            protocol_index, 1, p_hs[hs_counter]->p_myof->myof_a_on);

        for (int i = 0; i < p_model[hs_counter]->p_m_scheme[0]->no_of_states; i++)
        {
            gsl_matrix_set(p_data->p_hs_data[hs_counter]->hs_m_pops,
                protocol_index, i, gsl_matrix_get(p_hs[hs_counter]->p_myof->m_state_pops, 0, i));
        }
    }
}

void muscle::implement_single_half_sarcomere(double time_step_s, double pCa, double sim_mode, double delta_hsl)
{
    //! Code implements a time-step for a single half-sarcomere

    // Variables

    double new_length;
    double adjustment;

    int hs_counter;

    // Simple model
    hs_counter = 0;

    // Update the kinetics
    p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;
    p_hs[hs_counter]->sarcomere_kinetics(time_step_s, pCa);

    // Semi-clever check for comparing sim_mode to -1.0
    if (gsl_fcmp(sim_mode, -1.0, 1e-3) == 0)
    {
        // Check slack length mode for ktr
        p_hs[hs_counter]->hs_slack_length =
            p_hs[hs_counter]->return_hs_length_for_stress(0.0, time_step_s);

        // The hs_length cannot be shorter than its slack length
        new_length = GSL_MAX(p_hs[hs_counter]->hs_slack_length,
            p_hs[hs_counter]->hs_command_length);

        adjustment = new_length - p_hs[hs_counter]->hs_length;
    }
    else
    {
        // Over-write slack length
        p_hs[hs_counter]->hs_slack_length = GSL_NAN;

        // Are we in force control
        if (sim_mode >= 0.0)
        {
            // Force control
            new_length = p_hs[hs_counter]->return_hs_length_for_stress(sim_mode, time_step_s);

            adjustment = new_length - p_hs[hs_counter]->hs_length;

            // Update the command length which changes depending on force control
            p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_length;
        }
        else
        {
            // Length control
            adjustment = delta_hsl;
        }
    }

    // Apply the adjustment
    p_hs[hs_counter]->change_hs_length(adjustment, time_step_s);
}

double muscle::return_muscle_length(void)
{
    //! Code returns muscle length

    // Variables

    double holder;

    // Code

    // Sum up half-sarcomere length
    holder = 0.0;
    for (int i = 0; i < p_model[0]->no_of_half_sarcomeres; i++)
    {
        holder = holder + p_hs[i]->hs_length;
    }

    // Check for series component
    if (p_sc != NULL)
    {
        holder = holder + p_sc->sc_extension;
    }

    // Return
    return holder;
}

int wrapper_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* p, gsl_vector* f)
{
    //! This is a wrapper around muscle::check_residuals_for_myofibril_length_control()
    //! that handles the re-casting of pointers

    // Variables
    int f_return_value;

    struct m_control_params* params =
        (struct m_control_params*)p;

    // Code

    muscle* p_m = params->p_m;

    f_return_value = (int)p_m->worker_length_control_myofibril_with_series_compliance(x, params, f);

    return f_return_value;
}

size_t muscle::length_control_myofibril_with_series_compliance(double time_step_s, double pCa, double delta_hsl)
{
    //! Tries to impose length control on a system with a series compliance
    //! and/or 1 or more half-sarcomeres

    // Try to find a vector x such that the forces in each half-sarcomere
    // and the force in the series elastic element (which is the length
    // of the muscle - the combined length of the half-sarcomeres)
    // are all equal. x is filled with an initial guess.

    // Variables

    // Stuff to do with the root finding
    int status;
    int myofibril_iterations;
    size_t x_length;
    gsl_vector* x;

    // Other stuff
    double holder_hs_length;

    size_t max_lattice_iterations;

    // Code

    // First adjust the muscle length
    m_length = m_length + ((double)p_model[0]->no_of_half_sarcomeres * delta_hsl);

    // Now run kinetics
    for (int hs_counter = 0; hs_counter < p_model[0]->no_of_half_sarcomeres; hs_counter++)
    {
        p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;

        p_hs[hs_counter]->sarcomere_kinetics(time_step_s, pCa);
    }

    // Now deduce the length of x
    x_length = p_model[0]->no_of_half_sarcomeres + 1;

    // Allocate the vector
    x = gsl_vector_alloc(x_length);

    // The x-vector has all the half-sarcomere lengths followed by the force in the first_half_sarcomere
    for (int hs_counter = 0; hs_counter < p_model[0]->no_of_half_sarcomeres; hs_counter++)
    {
        gsl_vector_set(x, hs_counter, p_hs[hs_counter]->hs_length);
    }
    gsl_vector_set(x, x_length - 1, p_sc->sc_force);

    // Do the root finding
    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;
    const size_t calculation_size = x_length;

    m_control_params* par = new m_control_params;
    par->p_m = this;
    par->time_step = time_step_s;

    gsl_multiroot_function f = { &wrapper_length_control_myofibril_with_series_compliance, calculation_size, par };

    //T = gsl_multiroot_fsolver_hybrid;
    T = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc(T, calculation_size);
    gsl_multiroot_fsolver_set(s, &f, x);

    myofibril_iterations = 0;

    do
    {
        gsl_vector* y = gsl_vector_alloc(x_length);

        status = gsl_multiroot_fsolver_iterate(s);

        myofibril_iterations++;

        if (status)
        {
            printf("Myofibril multiroot solver break - Status: %i\t", status);

            if (status == GSL_EBADFUNC)
            {
                printf("Bad function value\n");
            }

            if (status == GSL_ENOPROG)
            {
                printf("Not making progress\n");
            }

            if (status == GSL_ENOPROGJ)
            {
                printf("Jacobian evaluations are not helping\n");
            }
        }

        status = gsl_multiroot_test_delta(s->dx, s->x, p_options->myofibril_force_tolerance, 0);

        gsl_vector_free(y);
    } while ((status == GSL_CONTINUE) && (myofibril_iterations < p_options->myofibril_max_iterations));

    // At this point, the s->x vector contains the lengths of the n half-sarcomeres
    // followed by the force in the series element

    // Implement the change
    holder_hs_length = 0.0;

    for (int hs_counter = 0; hs_counter < p_model[0]->no_of_half_sarcomeres; hs_counter++)
    {
        double new_hs_length = gsl_vector_get(s->x, hs_counter);
        double delta_hsl = new_hs_length - p_hs[hs_counter]->hs_length;

        p_hs[hs_counter]->change_hs_length(delta_hsl, time_step_s);

        holder_hs_length = holder_hs_length + new_hs_length;
    }

    p_sc->sc_extension = (m_length - holder_hs_length);
    p_sc->sc_force = p_sc->return_series_force(p_sc->sc_extension);

    // Update muscle force
    m_force = p_sc->sc_force;

    // Tidy up
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    delete par;

    // Return the max number of lattice iterations
    return 1.0;
}

size_t muscle::worker_length_control_myofibril_with_series_compliance(
    const gsl_vector* x, void* p, gsl_vector* f)
{
    //! This code calculates the f_vector that is minimized towards 0 to
    //! ensure that the force in the myofibril and its length are constrained

    // Variables
    struct m_control_params* params =
        (struct m_control_params*)p;

    double delta_hsl;
    double cum_hs_length;
    double target_force;
    double test_se_length;
    double force_diff;

    // Code

    // The f-vector is the difference between the force in each half-sarcomere and
    // the force in the series component
    // The x vector has a series of lengths followed by a muscle force
    // Calculate the force in each half-sarcomere and compare it to the force
    // Store up the half-sarcomere lengths as you go, and use that to calculate the length
    // of the series component

    // We need the force-control params for the calculation

    // Serial operation

    // Extract the target force
    target_force = gsl_vector_get(x, x->size - 1);

    cum_hs_length = 0;

    for (int hs_counter = 0; hs_counter < p_model[0]->no_of_half_sarcomeres; hs_counter++)
    {
        // Update the length
        cum_hs_length = cum_hs_length + gsl_vector_get(x, hs_counter);

        delta_hsl = gsl_vector_get(x, hs_counter) - p_hs[hs_counter]->hs_length;

        // Constrain delta_hsl to a plausible range
        if (delta_hsl > p_options->myofibril_max_delta_hs_length)
            delta_hsl = p_options->myofibril_max_delta_hs_length;
        if (delta_hsl < -p_options->myofibril_max_delta_hs_length)
            delta_hsl = -p_options->myofibril_max_delta_hs_length;

        force_diff = p_hs[hs_counter]->return_stress_after_delta_hsl(delta_hsl) - target_force;

        gsl_vector_set(f, hs_counter, force_diff);
    }

    // Now deduce the series elastic force
    test_se_length = m_length - cum_hs_length;
    force_diff = p_sc->return_series_force(test_se_length) - gsl_vector_get(x, x->size - 1);

    gsl_vector_set(f, f->size - 1, force_diff);

    return GSL_SUCCESS;
}

void muscle::write_rates_file(void)
{
    //! Function writes the m and c rate functions to file in JSON format

    // Variables
    int isotype_counter;					// isotype counter

    char file_write_mode[_MAX_PATH];		// mode for opening file

    char JSON_append_string[_MAX_PATH];		// written after scheme to keep JSON
    // structure, should be , if other entries follow
    // otherwise ""

    FILE* output_file;						// pointer for output file

    // Make sure directory exists
    std::filesystem::path output_file_path(p_options->rate_file_string);

    if (!(is_directory(output_file_path.parent_path())))
    {
        if (create_directories(output_file_path.parent_path()))
        {
            std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
        }
        else
        {
            std::cout << "\nError: folder for rates file could not be created: " <<
                output_file_path.parent_path().string() << "\n";
            exit(1);
        }
    }

    // Check file can be opened in write mode, abort if not
    errno_t err = fopen_s(&output_file, p_options->rate_file_string, "w");

    if (err != 0)
    {
        printf("muscle::write_rates_file(): %s\ncould not be opened\n",
            p_options->rate_file_string);
        exit(1);
    }

    // Start JSON structure
    fprintf_s(output_file, "{\n\t\"MyoSim_rates\":\n\t{\n");
    fprintf_s(output_file, "\t\t\"myosin\":\n");
    fprintf_s(output_file, "\t\t[\n");
    fclose(output_file);

    // Set the file write mode
    sprintf_s(file_write_mode, _MAX_PATH, "a");

    // Now cycle through the m isotypes
    for (isotype_counter = 0; isotype_counter < p_model[0]->m_no_of_isotypes; isotype_counter++)
    {
        // Set the append string
        if (isotype_counter < (p_model[0]->m_no_of_isotypes - 1))
        {
            sprintf_s(JSON_append_string, _MAX_PATH, ",");
        }
        else
        {
            sprintf_s(JSON_append_string, _MAX_PATH, "");
        }

        p_model[0]->p_m_scheme[isotype_counter]->write_rate_functions_to_file(
            p_options->rate_file_string, file_write_mode,
            JSON_append_string, p_hs[0]);
    }

    // Tidy up
    fopen_s(&output_file, p_options->rate_file_string, "a");
    fprintf(output_file, "\t\t]\n");
    fprintf(output_file, "\t}\n");
    fprintf(output_file, "}\n");
    fclose(output_file);
}

