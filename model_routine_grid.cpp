/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, UBEnv& ubEnv ) {
	/* MODEL START */

  CHECK( ubEnv.getIfGridPhiArray().size() == NUM_DIFFUSIBLE_ELEMS);
  //CHECK( ubEnv.getModelRealArray().size() == NUM_GRID_SUMMARY_REALS);
  //CHECK( ubEnv.getModelIntArray().size() == 0);

  ubEnv.setIfGridPhi(DIFFUSIBLE_ELEM_CHEMOATTRACTANT, 0.0);

  /* for( S32 elemIdx = 0 ; elemIdx < NUM_GRID_SUMMARY_REALS ; elemIdx++ ) {
    ubEnv.setModelReal(elemIdx, 0.0);
    }*/
	/* MODEL END */

	return;
}

void ModelRoutine::initIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridKappa ) {/* relevant only if v_gridPhiOutputDivideByKappa[pdeIdx] is set to true in updateFileOutputInfo() */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 iter, const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, NbrUBEnv& nbrUBEnv/* [INOUT] */ ) {
	/* MODEL START */

  if (pre) {
    for (S32 elemIdx = 0; elemIdx < NUM_DIFFUSIBLE_ELEMS; elemIdx++ ) {
	REAL phi = nbrUBEnv.getIfGridPhi(0, 0, 0, elemIdx);
	const UBAgentData& ubAgentData = *(nbrUBAgentData.getConstPtr(0, 0, 0));
	for(U16 l = 0; l < (U16)ubAgentData.v_spAgent.size(); l++) {
	  const SpAgent& spAgent = ubAgentData.v_spAgent[l];
	  S32 agentType = spAgent.state.getType();
	  phi += A_CELL_CHEMOATTRACTANT_SECRETION_RATE[agentType];
	}
	nbrUBEnv.setIfGridPhi(0, 0, 0, elemIdx, phi);
    }
  }

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridKappa ) {
	/* MODEL START */

  gridKappa = 1.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridKappaDomainBdry( const S32 pdeIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& kappa ) {

  kappa = 1.0;

  return;

}

void ModelRoutine::updateIfSubgridAlpha( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

  gridAlpha = -DIFFUSIBLE_ELEM_DECAY_RATE[elemIdx];

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& subgridVOffset0, const UBAgentData& ubAgentData0, const UBEnv& ubEnv0, const VIdx& vIdx1, const VIdx& subgridVOffset1, const UBAgentData& ubAgentData1, const UBEnv& ubEnv1, REAL& gridBeta ) {
	/* MODEL START */

  gridBeta = DIFFUSIBLE_ELEM_DIFFUSION_COEFFICIENT[elemIdx];

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
	/* MODEL START */

  gridBeta = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

  gridRHS = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::adjustIfSubgridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const REAL gridPhi, REAL& gridRHS/* INOUT */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {/* for Wnt & SFRP */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
	/* MODEL START */

  for(S32 pdeIdx = 0; pdeIdx < NUM_DIFFUSIBLE_ELEMS; pdeIdx++ ) {
    v_finestLevel[pdeIdx] = 2;
  }

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridDirichletBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdryPhi) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridNeumannBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdrySlope) {
	/* MODEL START */

  bdrySlope = 0.0;

	/* MODEL END */

	return;
}
