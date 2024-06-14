#pragma once

/**
/* @file		myofilaments.h
/* @brief		Header file for a myofilaments object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include <iostream>

#include "gsl_vector.h"
#include "gsl_matrix.h"

#include "global_definitions.h"

// Forward declararations
class half_sarcomere;
class model;
class options;
class kinetic_scheme;

using namespace std;

class myofilaments
{
public:
	/**
	 * Constructor
	 */
	myofilaments(half_sarcomere* set_p_parent_hs);

	/**
	* Destructor
	*/
	~myofilaments(void);

	// Variables
	half_sarcomere* p_parent_hs;			/**< Pointer to parent half-sarcomere */

	model* p_model;							/**< Pointer to model */

	options* p_options;						/**< Pointer to options */

	double myof_head_density;				/**< double with number of cross-bridges in
													a half-sarcomere with an area of
													1 m^-2 */

	double myof_k_cb;						/**< double with stiffness of cross-bridge
													link in N m^-1 */

	double myof_a_k_on;						/**< double with a_k_on, rate constant
													for thin filament activation
													M^-1 s^-1 */

	double myof_a_k_off;					/**< double with a_k_off, rate constant
													for thin filament deactivation
													s^-1 */

	double myof_a_k_coop;					/**< double for thin filament
													cooperativity, dimensionless */

	int no_of_bin_positions;				/**< integer with the number of
													positions cross-bridge
													distributions are evaluated
													at */

	int myof_m_no_of_isotypes;				/**< integer with the number of myosin isotypes */

	gsl_vector* x;							/**< gsl_vector with bin positions */

	size_t y_length;						/**< number of elements for a single m_isotype, caclulated as
												 no_of_detached_states + (no_of_attached_states * no_of_bins) +
												 2 for thin filament */
	
	size_t calc_length;						/**<number of elements in the array calculated as
												 m_no_of_isotypes *
													(no_of_detached_states + (no_of_attached_states * no_of_bins)) +
												2 for thin filament */

	gsl_matrix* y;							/**< gsl_matrix with the system */

	size_t a_off_index;						// indices
	size_t a_on_index;

	gsl_matrix_int* m_y_indices;			/**< gsl_matrix_int holding the start and
													stop indices in y for
													myosin state */

	gsl_matrix* m_state_pops;				/**< gsl_vector holding m_state populations */

	double* m_pops_array;					/**< array of doubles holding state populations
													used to communicate with gsl_results */

	double myof_a_off;						/**< double with proportion of off thin sites */

	double myof_a_on;						/**< double with proportion of on thin sites */

	double myof_m_bound;					/**< double with proportion of myosins
													that are bound */

	double myof_f_overlap;					/**< double defining the proportion
													of filaments that are in overlap */

	gsl_matrix* m_state_stresses;			/**< gsl_vector holding m_state stresses */

	double* m_stresses_array;				/**< array of doubles holding cross-bridge
													state stresses, used to communicate
													with cmv_results */

	double myof_stress_cb;					/**< double holding the cross-bridge force */

	double myof_stress_int_pas;				/**< double holding the int pas stress */

	double myof_stress_ext_pas;				/**< double holding the ext pas stress */

	double myof_stress_viscous;				/**< double holding the viscous stress */

	double myof_stress_myof;				/**< double holding the cb and int pas stress */

	double myof_stress_total;				/**< double holding the total stress */

	double myof_fil_compliance_factor;		/**< double holding the myofilament compliance
													factor. When the half-sarcomere changes
													length by x, cross-bridges are displaced
													by x * this variable */

	double myof_thick_fil_length;		/**< double describing thick filament length */

	double myof_bare_zone_length;		/**< double describing bare zone length */

	double myof_thin_fil_length;		/**< double describing thin filament length */

	double myof_ATP_flux;				/**< double decribing flux of ATP requiring
												transitions */

	double myof_int_pas_sigma;			/**< double with sigma for intracellular passive stress in
												N m^-2 */

	double myof_int_pas_L;				/**< double with curvature for intracellular passive stress
												exponential in nm */

	double myof_int_pas_slack_hsl;		/**< double with slack length for intracellular passive stress
												exponential in nm */

	double myof_ext_pas_sigma;			/**< double with sigma for extracellular passive stress in
												N m^-2 */

	double myof_ext_pas_L;				/**< double with curvature for extracellular passive stress
												exponential in nm */

	double myof_ext_pas_slack_hsl;		/**< double with slack length for extracellular passive stress
												exponential in nm */

	double max_shift;
	int n_max_sub_steps;

	double obj_time_step_s;

	double myof_mean_stress_int_pas;	/**< double holding the mean pas int
												stress over a cardiac cycle */

	string cb_dump_file_string;			/**< string hold cb dump file */
	bool cb_dump_file_defined;

	// Functions

	void initialise_simulation(void);

	void implement_time_step(double time_step_s);

	void run_kinetics(double time_step_s);

	void calculate_f_overlap(void);

	void calculate_m_state_pops(const double y[]);

	void calculate_m_state_stresses(void);

	void calculate_stresses(bool check_only = false);

	double calculate_cb_stress(bool check_only = false);

	double calculate_int_pas_stress(bool check_only = false, double delta_hsl = 0.0);

	double calculate_ext_pas_stress(bool check_only = false, double delta_hsl = 0.0);

	double return_stress_after_delta_hsl(double delta_hsl);

	void move_cb_populations(double delta_hsl);

	//void dump_cb_distributions(void);
};
