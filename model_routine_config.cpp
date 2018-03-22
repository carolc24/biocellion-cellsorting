/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
	/* MODEL START */

        /* Initialize Interface Grid Spacing (which was declared in model_define.h) */
        ifGridSpacing = IF_GRID_SPACING;
	
	/* MODEL END */

	return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
	/* MODEL START */

  callInfo.numComputeMechIntrctIters = 0;
  callInfo.numUpdateIfGridVarPreStateAndGridStepIters = 1;
  callInfo.numUpdateIfGridVarPostStateAndGridStepIters = 0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
	/* MODEL START */

        /* We set the Boundary Types of the Simulation Domain 
         There are only two Boundary Types: DOMAIN_BDRY_TYPE_PERIODIC, DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL */
        a_domainBdryType[0] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-x direction walls
        a_domainBdryType[1] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-y direction walls
        a_domainBdryType[2] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-z direction walls

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBdryType( pde_buffer_bdry_type_e& pdeBufferBdryType ) {
	/* MODEL START */

        /* Set the Boundary Type at between Buffer space, and Simulation space
	 We don't use buffers in this model, so set as PDE_BUFFER_BDRY_TYPE_HARD_WALL */
        pdeBufferBdryType = PDE_BUFFER_BDRY_TYPE_HARD_WALL;

	/* MODEL END */

	return;
}

void ModelRoutine::updateTimeStepInfo( TimeStepInfo& timeStepInfo ) {
	/* MODEL START */

        /* There are two different TimeStep sizes in Biocellion: 
	Baseline time step and the finer, num state and grid time step */
        timeStepInfo.durationBaselineTimeStep = BASELINE_TIME_STEP_DURATION;
        timeStepInfo.numStateAndGridTimeStepsPerBaseline = NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

	/* MODEL END */

	return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& mechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
	/* MODEL START */

	/* Update the sync method in model_routine_mech_intrct.cpp and model_routine_grid.cpp
	 Since mechanical interactions and PDEs are not used in this model, these are dummy */ 
        mechIntrctSyncMethod = SYNC_METHOD_PER_ATTR;
        updateIfGridVarSyncMethod = SYNC_METHOD_PER_ATTR;

	/* MODEL END */

	return;
}

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {/* set the mechanical interaction range & the numbers of model specific variables */
	/* MODEL START */

        /* Provide information about the discrete agent types in the user model */ 
        v_spAgentInfo.resize( NUM_CELL_TYPES );

        for( S32 i = 0 ; i < NUM_CELL_TYPES ; i++ ) {
                SpAgentInfo info;

                info.dMax = IF_GRID_SPACING;
                info.numBoolVars = 0;
                info.numStateModelReals = 0;
                info.numStateModelInts = 0;
                info.v_mechIntrctModelRealInfo.clear();
                info.v_mechIntrctModelIntInfo.clear();
                info.v_odeNetInfo.clear();

                v_spAgentInfo[i] = info;
        }

	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateJunctionEndInfo( Vector<JunctionEndInfo>& v_junctionEndInfo ) {/* set the numbers of model specific variables */
	/* MODEL START */

	/* No junctions in this model */
	v_junctionEndInfo.clear();
	
	/* MODEL END */

	return;
}

void ModelRoutine::updatePhiPDEInfo( Vector<PDEInfo>& v_phiPDEInfo ) {
	/* MODEL START */

  CHECK(NUM_DIFFUSIBLE_ELEMS == 1);
  v_phiPDEInfo.resize( NUM_DIFFUSIBLE_ELEMS );

  PDEInfo pdeInfo;
  GridPhiInfo gridPhiInfo;

  pdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_TIME_DEPENDENT_LINEAR;
  pdeInfo.numLevels = 3;
  pdeInfo.v_tagExpansionSize.assign( 3, 0 );
  pdeInfo.numTimeSteps = 1;
  pdeInfo.callAdjustRHSTimeDependentLinear = false;

  pdeInfo.mgSolveInfo.numPre = 3;
  pdeInfo.mgSolveInfo.numPost = 3;
  pdeInfo.mgSolveInfo.numBottom = 3;
  pdeInfo.mgSolveInfo.vCycle = true;
  pdeInfo.mgSolveInfo.maxIters = 30;
  pdeInfo.mgSolveInfo.epsilon = 1e-8;
  pdeInfo.mgSolveInfo.hang = 1e-8;
  pdeInfo.mgSolveInfo.normThreshold = 1e-18;
  pdeInfo.pdeIdx = 0;

  pdeInfo.advectionInfo.courantNumber = 0.5;  //dummy

  gridPhiInfo.elemIdx = DIFFUSIBLE_ELEM_CHEMOATTRACTANT;
  gridPhiInfo.name = "chemoattractant";

  /* PDE Boundary conditions */
  gridPhiInfo.aa_bcType[0][0] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcType[0][1] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcType[1][0] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcType[1][1] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcType[2][0] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcType[2][1] = BC_TYPE_NEUMANN_CONST;
  gridPhiInfo.aa_bcVal[0][0] = 0.0;
  gridPhiInfo.aa_bcVal[0][1] = 0.0;
  gridPhiInfo.aa_bcVal[1][0] = 0.0;
  gridPhiInfo.aa_bcVal[1][1] = 0.0;
  gridPhiInfo.aa_bcVal[2][0] = 0.0;
  gridPhiInfo.aa_bcVal[2][1] = 0.0;

  /* PDE error/warning and corrections */
  gridPhiInfo.errorThresholdVal = -1.0e-15;
  gridPhiInfo.warningThresholdVal = -1.0e-18;
  gridPhiInfo.setNegToZero = true;

  pdeInfo.v_gridPhiInfo.assign( 1, gridPhiInfo );

  v_phiPDEInfo[DIFFUSIBLE_ELEM_CHEMOATTRACTANT] = pdeInfo;

	
	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {
	/* MODEL START */

	/* No model_routine_grid.cpp in this model */
        v_ifGridModelRealInfo.clear();
        v_ifGridModelIntInfo.clear();

	/* MODEL END */

	return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
	/* MODEL START */

        /* We use two different RGN's in this model */
        CHECK( NUM_MODEL_RNGS == 2 );

        v_rngInfo.resize( NUM_MODEL_RNGS );

        RNGInfo rngInfo;
        /* Uniform distribution (min=0 and max=1) */
        rngInfo.type = RNG_TYPE_UNIFORM;
        rngInfo.param0 = 0.0;
        rngInfo.param1 = 1.0;
        rngInfo.param2 = 0.0;/* dummy */
        v_rngInfo[MODEL_RNG_UNIFORM] = rngInfo;

        /* Gaussian distribution (mean=0 and std=1) */
        rngInfo.type = RNG_TYPE_GAUSSIAN;
        rngInfo.param0 = 0.0;
        rngInfo.param1 = 1.0;
        rngInfo.param2 = 0.0;/* dummy */
        v_rngInfo[MODEL_RNG_GAUSSIAN] = rngInfo;

	/* MODEL END */

	return;
}

void ModelRoutine::updateFileOutputInfo( FileOutputInfo& fileOutputInfo ) {
	/* MODEL START */

	/* FileOutputInfo class holds the information related to file output of simulation results. */
        fileOutputInfo.particleOutput = true;                          
        //fileOutputInfo.particleNumExtraOutputVars = 0;

	fileOutputInfo.v_gridPhiOutput.assign( NUM_DIFFUSIBLE_ELEMS, true );
	fileOutputInfo.v_gridPhiOutputDivideByKappa.assign( NUM_DIFFUSIBLE_ELEMS, false);

	/* MODEL END */

	return;
}

void ModelRoutine::updateSummaryOutputInfo( Vector<SummaryOutputInfo>& v_summaryOutputRealInfo, Vector<SummaryOutputInfo>& v_summaryOutputIntInfo ) {
	/* MODEL START */

        /* Declare information you want to output when you run the simulation
         Biocellion supports a summary (reduction) mechanism for interface grid
         variables. Users add a fixed number of grid state variables to every unit box in the interface
         grid (or assign a fixed number of attributes to the interface grid), and for each attribute, users
         can ask Biocellion to visit every unit box to reduce the variables for the attribute to a single
         value. summary type e is used to set the reduction method. Choose one of these summary types:
         {SUMMARY_TYPE_SUM, SUMMARY_TYPE_AVG, SUMMARY_TYPE_MIN, SUMMARY_TYPE_MAX} */
	SummaryOutputInfo info;        
	v_summaryOutputIntInfo.clear();
	v_summaryOutputRealInfo.resize( NUM_GRID_SUMMARY_REALS );
	info.name = "Number of Live Cells";
	info.type = SUMMARY_TYPE_SUM;
	v_summaryOutputRealInfo[GRID_SUMMARY_REAL_LIVE_CELLS] = info;

	/* MODEL END */

	return;
}

void ModelRoutine::initGlobal( Vector<U8>& v_globalData ) {
	/* MODEL START */

	/* nothing to do */
	
	/* MODEL END */

	return;
}

void ModelRoutine::init( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::term( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::setPDEBuffer( const VIdx& startVIdx, const VIdx& regionSize, BOOL& isPDEBuffer ) {
	/* MODEL START */

	/* No buffers in this model */
	isPDEBuffer = false;

	/* MODEL END */

	return;
}

void ModelRoutine::setHabitable( const VIdx& vIdx, BOOL& isHabitable ) {
	/* MODEL START */

	/* Used to check if the initialized/simulated cells are placed in habitable space */
	isHabitable = true;

	/* MODEL END */

	return;
}

