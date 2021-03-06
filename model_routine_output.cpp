/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* UESR START */

#include "model_define.h"

/* UESR END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentOutput( const VIdx& vIdx, const SpAgent& spAgent, REAL& color, Vector<REAL>& v_extraScalar, Vector<VReal>& v_extraVector ) {
	/* MODEL START */

        color = spAgent.state.getType();
	v_extraScalar[PARTICLE_EXTRA_OUTPUT_RADIUS] = 1.0;
	v_extraVector[PARTICLE_EXTRA_OUTPUT_SCALE][0] = spAgent.state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_A);
	v_extraVector[PARTICLE_EXTRA_OUTPUT_SCALE][1] = spAgent.state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_B);
	v_extraVector[PARTICLE_EXTRA_OUTPUT_SCALE][2] = spAgent.state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_C);
	v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][0] = 0.0;
	v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][1] = 0.0;
	v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][2] = 0.0;

	REAL real;
	VReal vImg;
	Quaternion q;
	REAL thetaX, thetaY, thetaZ;
	
	real = spAgent.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A );
        vImg[0] = spAgent.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B );
        vImg[1] = spAgent.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C );
        vImg[2] = spAgent.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D );
        q.set( real, vImg );
        Quaternion::toTaitBryanIntrinsicZ1stX2ndY3rd( q, thetaZ, thetaX, thetaY );

        v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][0] = thetaX * ( 180.0 / MY_PI );/* paraview rotation: 1) rotateZ( orient[2] ), 2) rotateX( orient[0] ), 3) rotateY( orient[1] ) */
        v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][1] = thetaY * ( 180.0 / MY_PI );
        v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][2] = thetaZ * ( 180.0 / MY_PI );

        for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {
                CHECK( ( v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] >= -180.0 ) && ( v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] <= 180.0 ) );
                if( v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] < 0.0 ) {
                        v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] += 360.0;
                }
                CHECK( ( v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] >= 0.0 ) || ( v_extraVector[PARTICLE_EXTRA_OUTPUT_ORIENT][dim] < 360.0 ) );
        }


	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateSummaryVar( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<REAL>& v_realVal/* [elemIdx] */, Vector<S32>& v_intVal/* [elemIdx] */ ) {
	/* MODEL START */

        CHECK( v_realVal.size() == NUM_GRID_SUMMARY_REALS );
        CHECK( v_intVal.size() == 0 );

        const UBAgentData& ubAgentData = *( nbrUBAgentData.getConstPtr( 0, 0, 0 ) );

        REAL count = 0.0;

        /* Count the number of cells placed in the Simulation Domain */
        for (S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
                count += 1.0 ;
        }

        /* GRID_SUMMARY_REAL_LIVE_CELLS is set in model_routine_config.cpp */
        v_realVal[GRID_SUMMARY_REAL_LIVE_CELLS] = count;

	/* MODEL END */

	return;
}

