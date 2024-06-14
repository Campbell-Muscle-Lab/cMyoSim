/**
/* @file		model.cpp
/* @brief		Source file for a model object
/* @author		Ken Campbell
*/

#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>

#include "global_definitions.h"
#include "model.h"
#include "options.h"
#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"
#include "JSON_functions.h"

#include "gsl_math.h"
#include "gsl_vector.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

// Constructor
model::model(std::string JSON_model_file_string, options* set_p_options)
{
	// Initialise

	// Code
	std::cout << "In model constructor\n";

    p_options = set_p_options;

    initialise_model_from_JSON_file_string(JSON_model_file_string);
}

// Destructor
model::~model(void)
{
	// Code
	std::cout << "In model destructor\n";

    // Tidy up
    for (int i = 0; i < m_no_of_isotypes; i++)
    {
        delete p_m_scheme[i];
    }

    // And the vectors
    gsl_vector_free(m_isotype_props);
}

void model::initialise_model_from_JSON_file_string(std::string JSON_model_file_string)
{
	//! Code initialises a model from a JSON formatted file
	

    errno_t file_error;
    FILE* fp;
    char readBuffer[65536];

    // Code
    file_error = fopen_s(&fp, JSON_model_file_string.c_str(), "rb");
    if (file_error != 0)
    {
        std::cout << "Error opening JSON model file: " << JSON_model_file_string;
        exit(1);
    }

    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

    rapidjson::Document doc;
    doc.ParseStream(is);

    fclose(fp);

    // Now work through the JSON file
    JSON_functions::check_JSON_member_object(doc, "muscle");
    const rapidjson::Value& mus = doc["muscle"];

    JSON_functions::check_JSON_member_int(mus, "no_of_half_sarcomeres");
    no_of_half_sarcomeres = mus["no_of_half_sarcomeres"].GetInt();

    JSON_functions::check_JSON_member_int(mus, "no_of_myofibrils");
    no_of_myofibrils = mus["no_of_myofibrils"].GetInt();

    JSON_functions::check_JSON_member_number(mus, "initial_hs_length");
    initial_hs_length = mus["initial_hs_length"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "prop_fibrosis");
    prop_fibrosis = mus["prop_fibrosis"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "prop_myofilaments");
    prop_myofilaments = mus["prop_myofilaments"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "m_filament_density");
    m_filament_density = mus["m_filament_density"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "m_heads_per_filament");
    m_heads_per_filament = mus["m_heads_per_filament"].GetDouble();

    if (JSON_functions::check_JSON_member_exists(mus, "temperature"))
    {
        JSON_functions::check_JSON_member_number(mus, "temperature");
        temperature = mus["temperature"].GetDouble();
    }
    else
        temperature = 310;

    // Check if the series elastic component is specified
    sc_exists = false;
    sc_k_stiff = GSL_NAN;
    sc_sigma = GSL_NAN;
    sc_L = GSL_NAN;

    if (JSON_functions::check_JSON_member_exists(doc, "series_component"))
    {
        sc_exists = true;

        const rapidjson::Value& sc = doc["series_component"];

        // Check components
        if (JSON_functions::valid_JSON_member_number(sc, "sc_k_stiff"))
            sc_k_stiff = sc["sc_k_stiff"].GetDouble();
    
        if (JSON_functions::valid_JSON_member_number(sc, "sc_sigma"))
            sc_sigma = sc["sc_sigma"].GetDouble();
        
        if (JSON_functions::valid_JSON_member_number(sc, "sc_L"))
            sc_L = sc["sc_L"].GetDouble();
    }

    // Intracellular passive properties
    JSON_functions::check_JSON_member_object(doc, "intra_pas_parameters");
    const rapidjson::Value& intra_pas = doc["intra_pas_parameters"];

    JSON_functions::check_JSON_member_number(intra_pas, "int_pas_sigma");
    int_pas_sigma = intra_pas["int_pas_sigma"].GetDouble();

    JSON_functions::check_JSON_member_number(intra_pas, "int_pas_L");
    int_pas_L = intra_pas["int_pas_L"].GetDouble();

    JSON_functions::check_JSON_member_number(intra_pas, "int_pas_slack_hsl");
    int_pas_slack_hsl = intra_pas["int_pas_slack_hsl"].GetDouble();

    // Extracellular passive properties
    JSON_functions::check_JSON_member_object(doc, "extra_pas_parameters");
    const rapidjson::Value& extra_pas = doc["extra_pas_parameters"];

    JSON_functions::check_JSON_member_number(extra_pas, "ext_pas_sigma");
    ext_pas_sigma = extra_pas["ext_pas_sigma"].GetDouble();

    JSON_functions::check_JSON_member_number(extra_pas, "ext_pas_L");
    ext_pas_L = extra_pas["ext_pas_L"].GetDouble();

    JSON_functions::check_JSON_member_number(extra_pas, "ext_pas_slack_hsl");
    ext_pas_slack_hsl = extra_pas["ext_pas_slack_hsl"].GetDouble();

    // Thin filament properties
    JSON_functions::check_JSON_member_object(doc, "thin_parameters");
    const rapidjson::Value& thin_parameters = doc["thin_parameters"];

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_on");
    a_k_on = thin_parameters["a_k_on"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_off");
    a_k_off = thin_parameters["a_k_off"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_coop");
    a_k_coop = thin_parameters["a_k_coop"].GetDouble();

    // Now try to load myosin system
    JSON_functions::check_JSON_member_object(doc, "m_parameters");
    const rapidjson::Value& m_parameters = doc["m_parameters"];

    JSON_functions::check_JSON_member_number(m_parameters, "m_k_cb");
    m_k_cb = m_parameters["m_k_cb"].GetDouble();

    JSON_functions::check_JSON_member_number(m_parameters, "m_fil_compliance_factor");
    m_fil_compliance_factor = m_parameters["m_fil_compliance_factor"].GetDouble();

    JSON_functions::check_JSON_member_array(m_parameters, "m_isotype_proportions");
    const rapidjson::Value& mip = m_parameters["m_isotype_proportions"];

    m_no_of_isotypes = mip.Size();
    m_isotype_props = gsl_vector_alloc(MAX_NO_OF_ISOTYPES);
    gsl_vector_set_zero(m_isotype_props);

    for (int i = 0; i < (int)mip.Size(); i++)
    {
        gsl_vector_set(m_isotype_props, i, mip[i].GetDouble());
    }

    // Kinetic scheme for myosin - this is complicated so it's done in a different file
    JSON_functions::check_JSON_member_array(doc, "m_kinetics");
    const rapidjson::Value& m_ks = doc["m_kinetics"].GetArray();

    // Set the kinetic scheme for each isotype
    m_no_of_cb_states = 0;
    for (rapidjson::SizeType i = 0 ; i < (rapidjson::SizeType)m_no_of_isotypes ; i++)
    {
        p_m_scheme[i] = create_kinetic_scheme(m_ks[i], this, p_options);
        m_no_of_cb_states = GSL_MAX(m_no_of_cb_states, p_m_scheme[i]->no_of_states);
    }

    // Some calculations
    m_head_density = m_filament_density * m_heads_per_filament;
}

kinetic_scheme* model::create_kinetic_scheme(const rapidjson::Value& ks, model* p_model, options* p_options)
{
    //! Loads kinetic scheme

    // Variables
    kinetic_scheme* p_scheme;

    // Create the kinetic scheme
    p_scheme = new kinetic_scheme(ks, this, p_options);

    // Return the pointer
    return p_scheme;
}
