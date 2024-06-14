#pragma once

/**
/* @file		model.h
/* @brief		Header file for a model object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>
#include <string>

#include "gsl_vector.h"

// Definitions for JSON parsing
#ifndef _RAPIDJSON_DOCUMENT
#define _RAPIDJSON_DOCUMENT
#include "rapidjson/document.h"
#endif

#include "global_definitions.h"

// Forward declararations
class half_sarcomere;
class options;
class membranes;
class heart_rate;
class kinetic_scheme;

//
//class FiberSim_model;
//
//struct cmv_model_valve_structure;
//struct cmv_model_rc_structure;
//struct cmv_model_gc_structure;

class model
{
public:
	/**
	* Constructor
	*/
	model(std::string JSON_model_file_string, options* set_p_options);

	/**
	* Destructor
	*/
	~model(void);

	// Other variables
	options* p_options;					/**< Pointer to an options object */

	int no_of_half_sarcomeres;			/**< Number of half-sarcomeres in model */
	int no_of_myofibrils;				/**< Number of myofibrils in model */

	double initial_hs_length;			/**< double defining the initial hs_length in nm */

	double prop_fibrosis;				/**< double defining the proportion of the
                                             cross-section occupied by fibrosis. This
                                             contributes extracellular passive tension */

	double prop_myofilaments;           /**< double defining the proportion of the
											 non-fibrotic cross-section occupied by
											 myofilaments. This contributes titin and
											 cross-bridged mediated force */

	double m_filament_density;          /**< double defining the number of thick filaments
											 per square meter of cross-section */

	double m_heads_per_filament;		/**< double defining the number of heads per thick
											 filament */

	double m_head_density;				/**< double defining the number of heads in a half-sarcomere
											 per square meter of cross-section */

	bool sc_exists;						/**< bool defining whether a series component is specified */

	double sc_k_stiff;                  /**< double definiing the stiffness of the
											 series elastic component  in N m^-1 */

	double sc_sigma;					/**< double defining a multiplier for an exponential
											 series elastic component */

	double sc_L;						/**< double defining a curvature for an exponential
											 series elastic component */


	// Thin parameters
	double a_k_on;						/**< double defining the Ca2+-dependent on rate for thin
											 filament activation */

	double a_k_off;						/**< double defining the off rate for thin filament units */

	double a_k_coop;					/**< double defining the cooperativity of the thin filament */

	// Myofilament parameters
	double m_k_cb;						/**< double defining the stiffness of a cross-bridge link */

	double m_fil_compliance_factor;		/**< double defining the proportion of an applied length change
												that is imposed on cross-bridge links */

	double m_thick_fil_length;			/**< double defining the thick filament length in nm */

	double m_thin_fil_length;			/**< double defining the thin filament length in nm */

	double m_bare_zone_length;			/**< double defining the bare zone length in nm */

	double int_pas_sigma;				/**< double with sigma for intracellular passive stress in
												N m^-2 */

	double int_pas_L;					/**< double with curvature for intracellular passive stress
												exponential in nm */

	double int_pas_slack_hsl;			/**< double with slack length for intracellular passive stress
												exponential in nm */

	double ext_pas_sigma;				/**< double with sigma for extracellular passive stress in
												N m^-2 */

	double ext_pas_L;					/**< double with curvature for extracellular passive stress
												exponential in nm */

	double ext_pas_slack_hsl;			/**< double with slack length for extracellular passive stress
												exponential in nm */


	// Kinetics
	kinetic_scheme* p_m_scheme[MAX_NO_OF_ISOTYPES];
										/**< pointer to an array of kinetic schemes */

	int m_no_of_cb_states;				/**< integer with the maximum number of m states in a scheme */

	int m_no_of_isotypes;				/**< integer defining the number of isotypes */

	gsl_vector* m_isotype_props;		/**< gsl_vector holding the myosin isotype proportions */

	// Other stuff
	double temperature;					/**< double defining temperature in K */


	// Other functions

	/**
	/* Function initialises a model object from file
	*/
	void initialise_model_from_JSON_file_string(std::string JSON_model_file_string);

	/**
	/* Function returns a pointer to a myosin scheme
	*/
	kinetic_scheme* create_kinetic_scheme(const rapidjson::Value& ks, model* p_model, options* p_options);
};