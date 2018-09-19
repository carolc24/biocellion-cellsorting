/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

#ifndef __MY_DEFINE_H__
#define __MY_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

#define SYSTEM_DIMENSION 3

/* MODEL START */

/* Two cell types to be sorted */
typedef enum _cell_type_e {
        CELL_TYPE_A,
	CELL_TYPE_B,
        NUM_CELL_TYPES
} cell_type_e;

/* Real values associated with cells */
typedef enum _agent_state_real_e {
  AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_A,
  AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_B,
  AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_C,
  AGENT_STATE_REAL_ELLIPSOID_A,
  AGENT_STATE_REAL_ELLIPSOID_B,
  AGENT_STATE_REAL_ELLIPSOID_C,
  AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A,
  AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B,
  AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C,
  AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D,
  AGENT_STATE_REAL_MASS,
  AGENT_STATE_REAL_YOUNGS_MODULUS,
  AGENT_STATE_REAL_POISSONS_RATIO,
  NUM_AGENT_STATE_REALS
} agent_state_real_e;

/* Int values associated with cells */
typedef enum _agent_state_int_e {
  NUM_AGENT_STATE_INTS
} agent_state_int_e;

/* Internal reals */
typedef enum _agent_state_internal_real_e {
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_X,
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Y,
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Z,
  AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_X,
  AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Y,
  AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Z,
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_X,
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Y,
  AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Z,
  NUM_AGENT_STATE_INTERNAL_REALS
} agent_state_internal_real_e;

/* Internal ints */
typedef enum _agent_state_internal_int_e {
 NUM_AGENT_STATE_INTERNAL_INTS
} agent_state_internal_int_e;

/* Mechanical stuff */
typedef enum _agent_mech_real_e {
  AGENT_MECH_REAL_FORCE_X,
  AGENT_MECH_REAL_FORCE_Y,
  AGENT_MECH_REAL_FORCE_Z,
  AGENT_MECH_REAL_MOMENT_X,
  AGENT_MECH_REAL_MOMENT_Y,
  AGENT_MECH_REAL_MOMENT_Z,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_X,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_X,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Y,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Y,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Z,
  AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Z,
  NUM_AGENT_MECH_REALS
} agent_mech_real_e;

/* Two diffusible elements to act as chemoattractants */
typedef enum _diffusible_elem_e {
        DIFFUSIBLE_ELEM_CHEMOATTRACTANT,
        NUM_DIFFUSIBLE_ELEMS
} diffusible_elem_e;

/* Number of cells used in this model */
const S32 A_INI_N_CELLS[NUM_CELL_TYPES] = { 1000, 1000 };

/* Fixed radius of the cells used in the model */
const REAL A_CELL_RADIUS[NUM_CELL_TYPES] = { 1.0, 1.0 };

/* Interface Grid Variables to output in the Summary Ouput
 GRID_SUMMARY_REAL_LIVE_CELLS is the output variable for total number of cells that is available at each time output interval. */
typedef enum _grid_summary_real_e {
        GRID_SUMMARY_REAL_LIVE_CELLS,
        NUM_GRID_SUMMARY_REALS
} grid_summary_real_e;

/* A Uniform Random Number Generator. One should not use c++ function, rand(), as it is not thread-safe
 Instead, users should use Biocellion's built-in RNG */
typedef enum _model_rng_type_e {
        MODEL_RNG_UNIFORM,
        MODEL_RNG_GAUSSIAN,
        NUM_MODEL_RNGS
} model_rng_type_e;

typedef enum _extra_scalar_e {
  PARTICLE_EXTRA_OUTPUT_RADIUS,
  NUM_PARTICLE_EXTRA_OUTPUT_SCALAR_VARS
} extra_scalar_e;

/* Infomation needed to render spheres as ellipsoids  */
typedef enum _extra_vector_e {
  PARTICLE_EXTRA_OUTPUT_SCALE,
  PARTICLE_EXTRA_OUTPUT_ORIENT,
  NUM_PARTICLE_EXTRA_OUTPUT_VECTOR_VARS
} extra_vector_e;

/* IF_GRID_SPACING is the unit length of each voxel in the Simulatoion Domain
 The Simulation Domain size is set in the model XML file
 The Grid spacing can not be less than maximun cell agent diameter */
const REAL IF_GRID_SPACING = 4.0; // Micrometers

/* A baseline time step is the largest time step used in Biocellion
 Users can split a baseline time step into one or more state-and-grid time steps */
const REAL BASELINE_TIME_STEP_DURATION = 1.0;
const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1;

/* Required variables for chemotaxis */
const REAL DIFFUSIBLE_ELEM_DECAY_RATE[NUM_DIFFUSIBLE_ELEMS] = { 0.1 };
const REAL DIFFUSIBLE_ELEM_DIFFUSION_COEFFICIENT[NUM_DIFFUSIBLE_ELEMS] = { 0.5 };
const REAL A_CELL_CHEMOATTRACTANT_SECRETION_RATE[NUM_CELL_TYPES] = { 0.2, 0 };
const REAL A_CELL_CHEMOTAXIS_FORCE_STRENGTH[NUM_CELL_TYPES] = { 0.4, -0.1 };

/* Mechanical Interaction Force Constant */
typedef enum _cell_mech_real_e {
  CELL_MECH_REAL_FORCE_X,/* shoving & adhesion */
  CELL_MECH_REAL_FORCE_Y,/* shoving & adhesion */
  CELL_MECH_REAL_FORCE_Z,/* shoving & adhesion */
  NUM_CELL_MECH_REALS
} cell_mech_real_e;

/* Cell Adhesion Constant */
const REAL A_CELL_ADHESION = 0.5;

/* Cell Shoving Constant */
const REAL A_CELL_SPRING_CONSTANT = 0.05;

/* Cell diffusion coefficient */
const REAL A_CELL_DIFFUSION_COEFF[NUM_CELL_TYPES] = { 0.002, 0.002 };

/* Ellipsoid mechanics */
const S32 MECH_INTRCT_ELLIPSOID_MAX_ITERS = 100;
const REAL MECH_INTRCT_ELLIPSOID_EPSILON = 1e-10;

const REAL AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION = 0.001;
const S32 AGENT_TRANSLATION_ROTATION_INTEGRATION_STEPS_PER_BASELINE_TIME_STEP = (S32) ROUND( BASELINE_TIME_STEP_DURATION / 0.001);
const REAL AGENT_ARTIFICIAL_LINEAR_DRAG_COEFF_SCALE_FACTOR = 200;
const REAL AGENT_TRANSLATION_MAX_DISPLACEMENT_PER_BASELINE_TIME_STEP = 0.1;
const REAL AGENT_ARTIFICIAL_ANGULAR_DRAG_COEFF_SCALE_FACTOR = 10.0;
const REAL AGENT_ROTATION_MAX_ANGULAR_DISPLACEMENT_PER_BASELINE_TIME_STEP = (MY_PI / 180.0 ) * 5.0;

const REAL AGENT_DEFORMATION_NORMAL_STRESS_LARGE_DIFF_RATIO = 0.3;
const REAL AGENT_DEFORMATION_NORMAL_STRESS_TINY_DIFF_RATIO = 1e-5;
const REAL AGENT_DEFORMATION_NORMAL_STRESS_SMOOTHING_RATE = 0.01;
const REAL AGENT_DEFORMATION_MAX_STRETCH_RATIO_CHANGE_PER_BASELINE_TIME_STEP = 0.01;

const BOOL A_PERIODIC_DOMAIN[DIMENSION] = { true, true, true };

/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

