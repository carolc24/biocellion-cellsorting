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

/* One diffusible element to act as chemoattractant */
typedef enum _diffusible_elem_e {
        DIFFUSIBLE_ELEM_CHEMOATTRACTANT,
        NUM_DIFFUSIBLE_ELEMS
} diffusible_elem_e;

/* Number of cells used in this model */
const S32 A_INI_N_CELLS[NUM_CELL_TYPES] = { 100, 100 };

/* Fixed radius of the cells used in the model */
const REAL A_CELL_RADIUS[NUM_CELL_TYPES] = { 2.0, 2.0 };

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
const REAL DIFFUSIBLE_ELEM_DIFFUSION_COEFFICIENT[NUM_DIFFUSIBLE_ELEMS] = { 1.0 };
const REAL A_CELL_CHEMOATTRACTANT_SECRETION_RATE[NUM_CELL_TYPES] = { 0.1, 0 };
const REAL A_CELL_CHEMOTAXIS_FORCE_STRENGTH[NUM_CELL_TYPES] = { 0.5, 0 };
/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

