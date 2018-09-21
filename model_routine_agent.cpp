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

#if HAS_SPAGENT

static void computeAgentTranslationRotationAndDeformation( const VReal& vPos, const MechIntrctData& mechIntrctData, const JunctionData& junctionData, SpAgentState& state, VReal fwdDir, REAL chemForce, VReal& vDisp);
static void integrateTranslation( const VReal& vForce, VReal& oldStaggeredVLinear, VReal& vDisp, VReal& oldVDisp, Ellipsoid e, REAL m, REAL t);
static void integrateRotation( Quaternion& qMoment, const VReal& vMoment, VReal& oldStaggeredVAngular, Quaternion& qRot, double* a_inertia, REAL artificialAngularDragCoeff, REAL t);
static void computeStretchRatio( const JunctionData& junctionData, const REAL E, const REAL nu, const VReal& normalVStress, REAL a_stretchRatio[DIMENSION]);
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {/* initialization */
	/* MODEL START */

     	if( init == true ) {

      		/* Place initial cells on the Simulation Domain. 
                 When init is True, the placement occurs at the beginning of the Simulation */

	  //S32 numACells = A_INI_N_CELLS[0];
	  //S32 numBCells = A_INI_N_CELLS[1];
	  //S32 N = numACells + numBCells;

	  //for( S32 i = 0 ; i < NUM_CELL_TYPES ; i++ ){

                        /* Loop over all the initialized cells */
	  for( S32 j = 0 ; j < 2; j++ ){
                                VIdx vIdx;
				VReal vOffset;
				
				vIdx[0] = startVIdx[0];
				vIdx[1] = startVIdx[1];
				vIdx[2] = startVIdx[2] + (regionSize[2] - 1) * j;

				vOffset[0] = 0.0;
				vOffset[1] = 0.0;
				vOffset[2] = 0.0;
								
                                /* Initialize states of each cell */
				VReal vScale;
				Quaternion qRot;
				REAL E;
				REAL nu;
				REAL density;
				SpAgentState state;
				
				E = 1e1;
				nu = 0.5;  // incompressible
				density = 0.25;
				
				vScale = { A_CELL_RADIUS[j], A_CELL_RADIUS[j], 0.5 * A_CELL_RADIUS[j] };
				//qRot.set(1.0, Util::getModelRand(MODEL_RNG_UNIFORM), Util::getModelRand(MODEL_RNG_UNIFORM), Util::getModelRand(MODEL_RNG_UNIFORM));
				qRot.set(1.0, 0.0, 0.0, 0.0);
				qRot = Quaternion::normalize(qRot);
				
                                state.setType(j);                                                     
                                state.setModelReal(AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_A, vScale[0]);
				state.setModelReal(AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_B, vScale[1]);
				state.setModelReal(AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_C, vScale[2]);
				state.setModelReal(AGENT_STATE_REAL_ELLIPSOID_A, vScale[0]);
				state.setModelReal(AGENT_STATE_REAL_ELLIPSOID_B, vScale[1]);
				state.setModelReal(AGENT_STATE_REAL_ELLIPSOID_C, vScale[2]);
				state.setModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A, qRot.getReal() );
				state.setModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B, qRot.getImgI() );
				state.setModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C, qRot.getImgJ() );
				state.setModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D, qRot.getImgK() );
				state.setModelReal(AGENT_STATE_REAL_MASS, (4 * MY_PI / 3 ) * vScale[0] * vScale[1] * vScale[2] * density );
				state.setModelReal(AGENT_STATE_REAL_YOUNGS_MODULUS, E);
				state.setModelReal(AGENT_STATE_REAL_POISSONS_RATIO, nu);

				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_X, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Y, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Z, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_X, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Y, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Z, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_X, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Y, 0.0);
				state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Z, 0.0);
				state.setMechIntrctBdryEllipsoid(vScale, qRot);

                                /* Initialize return vectors {v_spAgentState, v_spAgentVIdx, v_spAgentOffset} by using .push_back() */
                                CHECK( ifGridHabitableBoxData.get( vIdx ) == true );
                                v_spAgentVIdx.push_back( vIdx );
				v_spAgentOffset.push_back( vOffset );
                                v_spAgentState.push_back( state );

                        }
                }
	/* MODEL END */

	return;
}

void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const NbrUBEnv& nbrUBEnv, const Vector<double>& v_y, Vector<double>& v_f ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
	/* MODEL START */

	/* nothing to do */	

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, BOOL& divide, BOOL& disappear ) {
	/* MODEL START */

 	divide = false;
        disappear = false;

	/* MODEL END */

	return;
}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& vDisp ) {/* if not dividing or disappearing */
	/* MODEL START */

  agentType_t type = state.getType();

  VReal vPos;

  for(S32 dim = 0; dim < DIMENSION; dim++) {
    vPos[dim] = ((REAL)vIdx[dim] + 0.5) * IF_GRID_SPACING + vOffset[dim];
  }
  printf("Position: %f, %f, %f\n", vPos[0], vPos[1], vPos[2]);
  
        /* Find a random unit vector */
  VReal fwdDir;
  VReal bckDir;
  REAL mag = 0.0;

  for (S32 dim = 0; dim < SYSTEM_DIMENSION; dim++) {
    fwdDir[dim] = -0.5 + Util::getModelRand(MODEL_RNG_UNIFORM);
    mag += fwdDir[dim] * fwdDir[dim];
  }
  mag = sqrt(mag);

  for (S32 dim = 0; dim < SYSTEM_DIMENSION; dim++) {
    fwdDir[dim] /= mag; // Normalize
    bckDir[dim] = -fwdDir[dim]; // Opposite vector
  }

  /* Find neighbor box that includes endpoint of fwdDir and bckDir */
  VIdx fwdOffset;
  VIdx bckOffset;

  for (S32 dim = 0; dim < SYSTEM_DIMENSION; dim++) {
    fwdOffset[dim] = fwdDir[dim] > 0.5 ? 1 : (fwdDir[dim] < 0.5 ? -1 : 0);
    bckOffset[dim] = bckDir[dim] > 0.5 ? 1 : (bckDir[dim] < 0.5 ? -1 : 0);
  }

  /* Read chemoattractant values in each of the two neighbor boxes */
  REAL fwdVal = nbrUBEnv.getValidFlag(fwdOffset) ? nbrUBEnv.getIfGridPhi(fwdOffset[0], fwdOffset[1], fwdOffset[2], DIFFUSIBLE_ELEM_CHEMOATTRACTANT) : 0.0;
  REAL bckVal = nbrUBEnv.getValidFlag(bckOffset) ? nbrUBEnv.getIfGridPhi(bckOffset[0], bckOffset[1], bckOffset[2], DIFFUSIBLE_ELEM_CHEMOATTRACTANT) : 0.0;

  //REAL fwdVal = nbrUBEnv.getIfGridPhi(fwdOffset[0], fwdOffset[1], fwdOffset[2], 0);
  //REAL bckVal = nbrUBEnv.getIfGridPhi(bckOffset[0], bckOffset[1], bckOffset[2], 0);

  //REAL tmpVal = nbrUBEnv.getIfGridPhi(fwdOffset, DIFFUSIBLE_ELEM_CHEMOATTRACTANT);
  //if (tmpVal > 0.0) printf("Chemoattractant: %f\n", tmpVal);

  //printf("Chemoattractant neighbors: %f %f %f\nVector: %d %d %d\nValue: %f\n", nbrUBEnv.getIfGridPhi(-1, -1, -1, 0), nbrUBEnv.getIfGridPhi(0,0,0,0), nbrUBEnv.getIfGridPhi(1,1,1,0), fwdOffset[0], fwdOffset[1], fwdOffset[2], fwdVal);
  
  /* Calculate chemotactic force vector and apply displacement if any */
  REAL chemForce = A_CELL_CHEMOTAXIS_FORCE_STRENGTH[state.getType()] * (fwdVal - bckVal); 

  computeAgentTranslationRotationAndDeformation(vPos, mechIntrctData, junctionData, state, fwdDir, chemForce, vDisp);

  VReal vScale;
  Quaternion qRot;

  qRot.set( state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D));
  vScale[0] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_A);
  vScale[1] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_B);
  vScale[2] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_C);

  state.setMechIntrctBdryEllipsoid(vScale, qRot);

  Ellipsoid e;
  e.vScale = vScale;
  e.qRot = qRot;

  for(S32 dim = 0; dim < DIMENSION - 1; dim++) {
    REAL lMin = 0.0;
    REAL rMax = Info::getDomainSize(dim) * IF_GRID_SPACING;

    if(A_PERIODIC_DOMAIN[dim] == false) {
      VReal vRet = e.getAxisAlignedBoundingBox();
      REAL lBound;
      REAL rBound;

      lBound = ((REAL) vIdx[dim] + 0.5) * IF_GRID_SPACING + vOffset[dim] + vDisp[dim] - vRet[dim];
      if (lBound < lMin)
	vDisp[dim] += lMin - lBound;

      rBound = ((REAL) vIdx[dim] + 0.5) * IF_GRID_SPACING + vOffset[dim] + vDisp[dim] - vRet[dim];
      if (rBound > rMax)
	vDisp -= rBound - rMax;
    }
  }

  for (int i=0; i<SYSTEM_DIMENSION; i++) {
    if (FABS(vDisp[i]) > IF_GRID_SPACING) {
      vDisp[i] = vDisp[i] * IF_GRID_SPACING / FABS(vDisp[i]);
    }
    CHECK(vDisp[i] <= IF_GRID_SPACING);
  }

  //printf("displacement after adjust: %f, %f, %f\n", vDisp[0], vDisp[1], vDisp[2]);
	/* MODEL END */

	return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

static void integrateTranslation(const VReal& vForce, VReal& oldStaggeredVLinear, VReal& vDisp, VReal& oldVDisp, Ellipsoid e, REAL m, REAL t) {
    REAL frontalArea;
    REAL artificialLinearDragCoeff;
    VReal newStaggeredVLinear;
    VReal vLinear;
    VReal newVDisp;
    VReal vScale = e.vScale;

    if ( oldStaggeredVLinear.lengthSquare() > 0.0) {
      frontalArea = e.compute2DProjectionArea( VReal::normalize( oldStaggeredVLinear ));
    } else if (vForce.lengthSquare() > 0.0 ) {
      CHECK( oldVDisp.lengthSquare() <= 0.0 );
      frontalArea = e.compute2DProjectionArea ( VReal::normalize(vForce));
    } else {
      REAL avgRadius = CBRT( vScale[0] * vScale[1] * vScale[2] );
      frontalArea = MY_PI * avgRadius * avgRadius;
    }

    artificialLinearDragCoeff = frontalArea * AGENT_ARTIFICIAL_LINEAR_DRAG_COEFF_SCALE_FACTOR;

    if (oldStaggeredVLinear.length() < ( artificialLinearDragCoeff * oldStaggeredVLinear.lengthSquare() ) * (t / m)) {
      WARNING("oldStaggeredVLinear: Time step too large, t=" << t);
    }

    if (oldStaggeredVLinear.lengthSquare() > 0.0) {
      newStaggeredVLinear = oldStaggeredVLinear + (vForce - VReal::normalize( oldStaggeredVLinear) * (artificialLinearDragCoeff * oldStaggeredVLinear.lengthSquare() ) ) * ( t / m );
    } else {
      newStaggeredVLinear = oldStaggeredVLinear + vForce * ( t / m);
    }

    newVDisp = vDisp + newStaggeredVLinear * t;
    vLinear = (newVDisp - oldVDisp) / (t * 2.0);

    if (vLinear.length() < (artificialLinearDragCoeff * vLinear.lengthSquare()) * (t / m)) {
      WARNING("vLinear: Time step too large, t=" << t);
    }

    if (vLinear.lengthSquare() > 0.0) {
      newStaggeredVLinear = oldStaggeredVLinear + (vForce - VReal::normalize( vLinear ) * (artificialLinearDragCoeff * vLinear.lengthSquare())) * (t / m);
    } else {
      newStaggeredVLinear = oldStaggeredVLinear + vForce * (t / m );
    }

    newVDisp = vDisp + newStaggeredVLinear * t;

    oldStaggeredVLinear = newStaggeredVLinear;
    oldVDisp = vDisp;
    vDisp = newVDisp;
}

static void integrateRotation(Quaternion& qMoment, const VReal& vMoment, VReal& oldStaggeredVAngular, Quaternion& qRot, double* a_inertia, REAL artificialAngularDragCoeff, REAL t) {
  
    VReal vTmp;
    Quaternion qTmp;

    VReal newStaggeredVAngular;
    VReal vAngular;

    VReal vAngularDot;
    VReal newStaggeredVAngularDot;

    SquareMatrix4 matQ;
    Quaternion qDot;
    Quaternion newStaggeredQDot;

    Quaternion newStaggeredQRot;
    Quaternion newQRot;

    vTmp[0] = ( -artificialAngularDragCoeff * oldStaggeredVAngular[0] ) / a_inertia[0];
    vTmp[1] = ( -artificialAngularDragCoeff * oldStaggeredVAngular[1] ) / a_inertia[1];
    vTmp[2] = ( -artificialAngularDragCoeff * oldStaggeredVAngular[2] ) / a_inertia[2];

    for ( S32 dim=0; dim<DIMENSION; dim++ ){
      if (FABS(oldStaggeredVAngular[dim]) < FABS( vTmp[dim] * t ) ) {
	WARNING("Time step too large, t=" << t);
      }
    }

    vAngularDot[0] = ((qMoment.getImgI() - artificialAngularDragCoeff * oldStaggeredVAngular[0] ) + (a_inertia[1] - a_inertia[2]) * oldStaggeredVAngular[1] * oldStaggeredVAngular[2] ) / a_inertia[0];
    vAngularDot[1] = ((qMoment.getImgJ() - artificialAngularDragCoeff * oldStaggeredVAngular[1] ) + (a_inertia[2] - a_inertia[0]) * oldStaggeredVAngular[2] * oldStaggeredVAngular[0] ) / a_inertia[1];
    vAngularDot[2] = ((qMoment.getImgK() - artificialAngularDragCoeff * oldStaggeredVAngular[2] ) + (a_inertia[0] - a_inertia[1]) * oldStaggeredVAngular[0] * oldStaggeredVAngular[1] ) / a_inertia[2];

    vAngular[0] = oldStaggeredVAngular[0] + vAngularDot[0] * (t * 0.5);
    vAngular[1] = oldStaggeredVAngular[1] + vAngularDot[1] * (t * 0.5);
    vAngular[2] = oldStaggeredVAngular[2] + vAngularDot[2] * (t * 0.5);

    vTmp[0] = ( -artificialAngularDragCoeff * vAngular[0] ) / a_inertia[0];
    vTmp[1] = ( -artificialAngularDragCoeff * vAngular[1] ) / a_inertia[1];
    vTmp[2] = ( -artificialAngularDragCoeff * vAngular[2] ) / a_inertia[2];

    for ( S32 dim=0; dim<DIMENSION; dim++) {
      if ( FABS(vAngular[dim]) < FABS( vTmp[dim] * t )) {
	WARNING("Time step too large, t=" << t);
      }
    }

    vAngularDot[0] = ((qMoment.getImgI() - artificialAngularDragCoeff * vAngular[0] ) + (a_inertia[1] - a_inertia[2]) * vAngular[1] * vAngular[2] ) / a_inertia[0];
    vAngularDot[1] = ((qMoment.getImgJ() - artificialAngularDragCoeff * vAngular[1] ) + (a_inertia[2] - a_inertia[0]) * vAngular[2] * vAngular[0] ) / a_inertia[1];
    vAngularDot[2] = ((qMoment.getImgK() - artificialAngularDragCoeff * vAngular[2] ) + (a_inertia[0] - a_inertia[1]) * vAngular[0] * vAngular[1] ) / a_inertia[2];

    vAngular[0] = oldStaggeredVAngular[0] + vAngularDot[0] * (t * 0.5);
    vAngular[1] = oldStaggeredVAngular[1] + vAngularDot[1] * (t * 0.5);
    vAngular[2] = oldStaggeredVAngular[2] + vAngularDot[2] * (t * 0.5);

    newStaggeredVAngular[0] = oldStaggeredVAngular[0] + vAngularDot[0] * t;
    newStaggeredVAngular[1] = oldStaggeredVAngular[1] + vAngularDot[1] * t;
    newStaggeredVAngular[2] = oldStaggeredVAngular[2] + vAngularDot[2] * t;

    matQ.aa_val[0][0] = qRot.getReal() * 0.5;
    matQ.aa_val[0][1] = -qRot.getImgI() * 0.5;
    matQ.aa_val[0][2] = -qRot.getImgJ() * 0.5;
    matQ.aa_val[0][3] = -qRot.getImgK() * 0.5;

    matQ.aa_val[1][0] = qRot.getImgI() * 0.5;
    matQ.aa_val[1][1] = qRot.getReal() * 0.5;
    matQ.aa_val[1][2] = -qRot.getImgK() * 0.5;
    matQ.aa_val[1][3] = -qRot.getImgJ() * 0.5;

    matQ.aa_val[2][0] = qRot.getImgJ() * 0.5;
    matQ.aa_val[2][1] = qRot.getImgK() * 0.5;
    matQ.aa_val[2][2] = qRot.getReal() * 0.5;
    matQ.aa_val[2][3] = -qRot.getImgI() * 0.5;

    matQ.aa_val[3][0] = qRot.getImgK() * 0.5;
    matQ.aa_val[3][1] = -qRot.getImgJ() * 0.5;
    matQ.aa_val[3][2] = qRot.getImgI() * 0.5;
    matQ.aa_val[3][3] = qRot.getReal() * 0.5;

    qDot = matQ.mult( Quaternion( 0.0, vAngular ) );

    newStaggeredQRot = qRot + qDot * (t * 0.5);

    matQ.aa_val[0][0] = newStaggeredQRot.getReal() * 0.5;
    matQ.aa_val[0][1] = -newStaggeredQRot.getImgI() * 0.5;
    matQ.aa_val[0][2] = -newStaggeredQRot.getImgJ() * 0.5;
    matQ.aa_val[0][3] = -newStaggeredQRot.getImgK() * 0.5;

    matQ.aa_val[1][0] = newStaggeredQRot.getImgI() * 0.5;
    matQ.aa_val[1][1] = newStaggeredQRot.getReal() * 0.5;
    matQ.aa_val[1][2] = -newStaggeredQRot.getImgK() * 0.5;
    matQ.aa_val[1][3] = -newStaggeredQRot.getImgJ() * 0.5;

    matQ.aa_val[2][0] = newStaggeredQRot.getImgJ() * 0.5;
    matQ.aa_val[2][1] = newStaggeredQRot.getImgK() * 0.5;
    matQ.aa_val[2][2] = newStaggeredQRot.getReal() * 0.5;
    matQ.aa_val[2][3] = -newStaggeredQRot.getImgI() * 0.5;

    matQ.aa_val[3][0] = newStaggeredQRot.getImgK() * 0.5;
    matQ.aa_val[3][1] = -newStaggeredQRot.getImgJ() * 0.5;
    matQ.aa_val[3][2] = newStaggeredQRot.getImgI() * 0.5;
    matQ.aa_val[3][3] = newStaggeredQRot.getReal() * 0.5;

    newStaggeredQDot = matQ.mult(Quaternion(0.0, newStaggeredVAngular));
    newQRot = qRot + newStaggeredQDot * t;

    newQRot = Quaternion::normalize(newQRot);

    CHECK(FABS(1.0 - newQRot.norm()) < REAL_EPSILON * 1e1);

    qTmp = Quaternion::qpqStar( Quaternion( 0.0, newStaggeredVAngular), qRot);
    qTmp = Quaternion::qStarpq( qTmp, newQRot );

    oldStaggeredVAngular = qTmp.getImg();

    qMoment = Quaternion::qStarpq( Quaternion( 0.0, vMoment), newQRot );

    qRot = newQRot;
}

static void computeAgentTranslationRotationAndDeformation(const VReal& vPos, const MechIntrctData& mechIntrctData, const JunctionData& junctionData, SpAgentState& state, VReal fwdDir, REAL chemForce, VReal& vDisp) {

  agentType_t type = state.getType();

  VReal vScale;
  Quaternion qRot;

  VReal vForce;
  VReal vMoment;
  REAL m;
  REAL theta;

  VReal oldStaggeredVLinear;
  VReal oldVDisp;

  Ellipsoid e;

  Quaternion qMoment;

  VReal oldStaggeredVAngular;
  

  REAL a_inertia[DIMENSION];
  REAL approxSurfaceArea;
  REAL artificialAngularDragCoeff;

  Quaternion orgQRot;
  REAL dotProdSquare;

  qRot.set( state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C), state.getModelReal(AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D));

  vScale[0] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_A);
  vScale[1] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_B);
  vScale[2] = state.getModelReal(AGENT_STATE_REAL_ELLIPSOID_C);

  vForce[0] = mechIntrctData.getModelReal(AGENT_MECH_REAL_FORCE_X) + chemForce * fwdDir[0];
  vForce[1] = mechIntrctData.getModelReal(AGENT_MECH_REAL_FORCE_Y) + chemForce * fwdDir[1];
  vForce[2] = mechIntrctData.getModelReal(AGENT_MECH_REAL_FORCE_Z) + chemForce * fwdDir[2];

  if (state.getType() == 0) vForce[2] += 5.0;
  if (state.getType() == 1) vForce[2] -= 5.0;

  vMoment[0] = mechIntrctData.getModelReal(AGENT_MECH_REAL_MOMENT_X);
  vMoment[1] = mechIntrctData.getModelReal(AGENT_MECH_REAL_MOMENT_Y);
  vMoment[2] = mechIntrctData.getModelReal(AGENT_MECH_REAL_MOMENT_Z);

  m = state.getModelReal(AGENT_STATE_REAL_MASS);

  //printf("Force: %f, %f, %f\n", vForce[0], vForce[1], vForce[2]);
  //vForce[2] -= m * 9.8;  // gravity

  /* Translation */

  e.vScale = vScale;
  e.qRot = qRot;

  oldStaggeredVLinear[0] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_X );
  oldStaggeredVLinear[1] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Y );
  oldStaggeredVLinear[2] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Z );

  for (int i=0; i<SYSTEM_DIMENSION; i++) {   
    oldVDisp[i] = oldStaggeredVLinear[i] * AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION * -1.0;
  }
  //printf("Old displacement: %f %f %f\n", oldVDisp[0], oldVDisp[1], oldVDisp[2]);

  if (vForce.dotProduct(oldStaggeredVLinear) < 0.0) {
    oldStaggeredVLinear = VReal::ZERO;
    oldVDisp = VReal::ZERO;
  }

  vDisp = VReal::ZERO;

  REAL timeElapsed = 0.0;
  REAL curTimeStep = AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION;
  REAL hMax = 0.0;
  VReal tmpVLinear = oldStaggeredVLinear;
  VReal tmpVDisp = vDisp;
  VReal tmpOldVDisp = oldVDisp;

  REAL error = 0.0;
  REAL maxError = 1e-8;

  while( timeElapsed < BASELINE_TIME_STEP_DURATION ) {
    timeElapsed += curTimeStep;
    
    integrateTranslation(vForce, oldStaggeredVLinear, vDisp, oldVDisp, e, m, curTimeStep);  // One whole step (vDisp = r_h)
    //integrateTranslation(vForce, tmpVLinear, tmpVDisp, tmpOldVDisp, e, m, curTimeStep / 2); 
    // integrateTranslation(vForce, tmpVLinear, tmpVDisp, tmpOldVDisp, e, m, curTimeStep / 2);  // Two half steps (tmpVDisp = r_h/2)

    /* error = 0.0;
    for (int i=0; i<DIMENSION; i++) {
      error += POW(vDisp[i] - tmpVDisp[i], 2);
    }
    error = (8.0/7.0) * POW(error, 0.5);
    
    hMax = POW(maxError / error, 0.25) * curTimeStep;
    if (hMax < 0.5 * curTimeStep) curTimeStep = 1.8 * hMax;  // correction step
    if (curTimeStep < 0.0001) curTimeStep = 0.0001;

    if (curTimeStep + timeElapsed > BASELINE_TIME_STEP_DURATION) curTimeStep = BASELINE_TIME_STEP_DURATION - timeElapsed;

    tmpVDisp = vDisp;
    tmpVLinear = oldStaggeredVLinear;
    tmpOldVDisp = oldVDisp;*/
  }

  if (vDisp.length() > 0.0) {
    VReal vPointOn = e.getIntersectingPointOn( VReal::normalize( vDisp ) );
    if (vDisp.length() > vPointOn.length() ) {
      WARNING("vDisp: Time step too large");
    }
  }

  //printf("displacement: %f, %f, %f\n", vDisp[0], vDisp[1], vDisp[2]);

  state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_X, oldStaggeredVLinear[0]);
  state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Y, oldStaggeredVLinear[1]);
  state.setInternalModelReal(AGENT_STATE_INTERNAL_REAL_STAGGERED_VELOCITY_Z, oldStaggeredVLinear[2]);

  /* rotation */

  a_inertia[0] = (vScale[1] * vScale[1] + vScale[2] * vScale[2] ) * (m / 5.0);
  a_inertia[1] = (vScale[0] * vScale[0] + vScale[2] * vScale[2] ) * (m / 5.0);
  a_inertia[2] = (vScale[1] * vScale[1] + vScale[0] * vScale[0] ) * (m / 5.0);

  approxSurfaceArea = POW((POW(vScale[0] * vScale[1], 1.6075) + POW(vScale[1] * vScale[2], 1.6075) + POW(vScale[2] * vScale[0], 1.6075) ) /3.0, 1/1.6075); // Wikipedia2018:Ellipsoid
  artificialAngularDragCoeff = approxSurfaceArea * AGENT_ARTIFICIAL_ANGULAR_DRAG_COEFF_SCALE_FACTOR;

  orgQRot = qRot;

  oldStaggeredVAngular[0] = state.getInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_X);
  oldStaggeredVAngular[1] = state.getInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Y);
  oldStaggeredVAngular[2] = state.getInternalModelReal(AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Z);

  CHECK(FABS(1.0 - qRot.norm()) < REAL_EPSILON * 1e1);

  qMoment = Quaternion::qStarpq( Quaternion(0.0, vMoment), qRot );

  timeElapsed = 0.0;
  curTimeStep = AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION;
  error = 0.0;

  Quaternion tmpQRot = qRot;
  VReal tmpVAngular = oldStaggeredVAngular;
  Quaternion tmpQMoment = qMoment;

  while (timeElapsed < BASELINE_TIME_STEP_DURATION) {
    timeElapsed += curTimeStep;
    
    integrateRotation(qMoment, vMoment, oldStaggeredVAngular, qRot, a_inertia, artificialAngularDragCoeff, curTimeStep);
    //integrateRotation(tmpQMoment, vMoment, tmpVAngular, tmpQRot, a_inertia, artificialAngularDragCoeff, curTimeStep / 2);
    //integrateRotation(tmpQMoment, vMoment, tmpVAngular, tmpQRot, a_inertia, artificialAngularDragCoeff, curTimeStep / 2);

    /*error = 0.0;
    error += POW(qRot.getReal() - tmpQRot.getReal(), 2);
    error += POW(qRot.getImgI() - tmpQRot.getImgI(), 2);
    error += POW(qRot.getImgJ() - tmpQRot.getImgJ(), 2);
    error += POW(qRot.getImgK() - tmpQRot.getImgK(), 2);

    error = (8.0/7.0) * POW(error, 0.5);

    hMax = POW(maxError/error, 0.25) * curTimeStep;
    if (hMax < 0.5 * curTimeStep) curTimeStep = 1.8 * hMax;
    if (curTimeStep + timeElapsed > BASELINE_TIME_STEP_DURATION) curTimeStep = BASELINE_TIME_STEP_DURATION - timeElapsed;

    tmpQRot = qRot;
    tmpVAngular = oldStaggeredVAngular;
    tmpQMoment = qMoment;*/
  }

  dotProdSquare = orgQRot.dotProduct(qRot);
  dotProdSquare *= dotProdSquare;

  if (dotProdSquare > 1.0) {
    CHECK(dotProdSquare - 1.0 < REAL_EPSILON * 1e1);
    dotProdSquare = 1.0;
  }
  theta = ACOS(2.0 * dotProdSquare - 1.0);
  if (theta > AGENT_ROTATION_MAX_ANGULAR_DISPLACEMENT_PER_BASELINE_TIME_STEP) {
    WARNING("time step too large");
  }

  state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_X, oldStaggeredVAngular[0]);
  state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Y, oldStaggeredVAngular[1]);
  state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_STAGGERED_ANGULAR_VELOCITY_Z, oldStaggeredVAngular[2]);

  state.setModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A, qRot.getReal());
  state.setModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B, qRot.getImgI());
  state.setModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C, qRot.getImgJ());
  state.setModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D, qRot.getImgK());

    VReal undeformedVScale;
    VReal oldNormalVStress;

    VReal lowNormalVForce;
    VReal highNormalVForce;

    REAL a_stretchRatio[DIMENSION];
    REAL E;
    REAL nu;

    VReal normalVForce;
    VReal normalVStress;
    VReal vDiff;

    undeformedVScale[0] = state.getModelReal( AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_A );
    undeformedVScale[1] = state.getModelReal( AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_B );
    undeformedVScale[2] = state.getModelReal( AGENT_STATE_REAL_UNDEFORMED_ELLIPSOID_C );

    oldNormalVStress[0] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_X );
    oldNormalVStress[1] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Y );
    oldNormalVStress[2] = state.getInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Z );

    lowNormalVForce[0] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_X );
    lowNormalVForce[1] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Y );
    lowNormalVForce[2] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Z );
    highNormalVForce[0] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_X );
    highNormalVForce[1] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Y );
    highNormalVForce[2] = mechIntrctData.getModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Z );

    if (state.getType() == 0) lowNormalVForce[2] += 5.0;
    if (state.getType() == 1) highNormalVForce[2] -= 5.0;
    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
      a_stretchRatio[dim] = vScale[dim] / undeformedVScale[dim];
      CHECK( a_stretchRatio[dim] > 0.0 );
    }

    E = state.getModelReal( AGENT_STATE_REAL_YOUNGS_MODULUS );
    nu = state.getModelReal( AGENT_STATE_REAL_POISSONS_RATIO );
    if( nu < 1e-2 ) {
      WARNING( "id=" << junctionData.getCurId() << " type=" << type << ", nu (" << nu << ") is too small, and this routine may not be numerically stable." );
    }

    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
      REAL low = lowNormalVForce[dim] * -1.0;/* -:compression, +:tension */
      REAL high = highNormalVForce[dim];/* -:compression, +:tension */
      if( low * high > 0.0 ) {
	//printf("Deforming.....");
	if( low < 0.0 ) {/* compression */
	  normalVForce[dim] = FMAX( low, high );
	}
	else {/* tension */
	  normalVForce[dim] = FMIN( low, high );
	}
      }
      else {
	normalVForce[dim] = 0.0;
      }
    }

    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
      normalVStress[dim] = normalVForce[dim] / ( vScale[( dim + 1 ) % DIMENSION] * vScale[( dim + 2 ) % DIMENSION] );/* this should be close to true stress assuming that vScale changes only slowly per baseline time-step */
    }
    vDiff = normalVStress - oldNormalVStress;

    if( vDiff.length() > E * AGENT_DEFORMATION_NORMAL_STRESS_LARGE_DIFF_RATIO ) {
      REAL scale = ( E * AGENT_DEFORMATION_NORMAL_STRESS_LARGE_DIFF_RATIO ) / vDiff.length();
      vDiff *= scale;
      normalVStress = oldNormalVStress + vDiff * AGENT_DEFORMATION_NORMAL_STRESS_SMOOTHING_RATE;
    }
    else if( vDiff.length() < E * AGENT_DEFORMATION_NORMAL_STRESS_TINY_DIFF_RATIO ) {
      normalVStress = oldNormalVStress + vDiff;
    }
    else {
      normalVStress = oldNormalVStress + vDiff * AGENT_DEFORMATION_NORMAL_STRESS_SMOOTHING_RATE;
    }
    printf("Stress: %f\n", normalVStress[2]);
    computeStretchRatio( junctionData, E, nu, normalVStress, a_stretchRatio );

    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
      REAL oldScale = vScale[dim];
      REAL newScale = undeformedVScale[dim] * a_stretchRatio[dim];
      printf("Stretch ratio: %f, %f, %f\n", a_stretchRatio[0], a_stretchRatio[1], a_stretchRatio[2]);
      CHECK( a_stretchRatio[dim] > 0.0 );
      CHECK( oldScale > 0.0 );
      CHECK( newScale > 0.0 );
      if( FABS( newScale / oldScale - 1.0 ) > AGENT_DEFORMATION_MAX_STRETCH_RATIO_CHANGE_PER_BASELINE_TIME_STEP ) {
	WARNING( "id=" << junctionData.getCurId() << " type=" << type << ", too large smoothing rate (" << AGENT_DEFORMATION_NORMAL_STRESS_SMOOTHING_RATE << "), new strain / old strain (" << ( newScale / oldScale ) << ", dim=" << dim << ") is too large." );
      }
      vScale[dim] = newScale;
    }

    state.setModelReal( AGENT_STATE_REAL_ELLIPSOID_A, vScale[0] );
    state.setModelReal( AGENT_STATE_REAL_ELLIPSOID_B, vScale[1] );
    state.setModelReal( AGENT_STATE_REAL_ELLIPSOID_C, vScale[2] );
    state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_X, normalVStress[0] );
    state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Y, normalVStress[1] );
    state.setInternalModelReal( AGENT_STATE_INTERNAL_REAL_BODY_FIXED_NORMAL_STRESS_Z, normalVStress[2] );
    
}
#endif

static void computeStretchRatio( const JunctionData& junctionData, const REAL E, const REAL nu, const VReal& normalVStress, REAL a_stretchRatio[DIMENSION]/* [INOUT] */ ) {
#if INCLUDE_DEFORMATION_HOOKS_LAW
  for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
    REAL trueStrain = ( 1 / E ) * ( normalVStress[dim] - nu * ( normalVStress[( dim + 1 ) % DIMENSION] + normalVStress[( dim + 2 ) % DIMENSION] ) );
    a_stretchRatio[dim] = EXP( trueStrain );
  }
#elif INCLUDE_DEFORMATION_NEO_HOOKEAN
  REAL C1 = E / ( 4.0 * ( 1.0 + nu ) );
  REAL initNorm = REAL_MAX;
  REAL norm = REAL_MAX;
  S32 iters = 0;
  SquareMatrix3 matJInv;/* inverse matrix of Jacobian */
  CHECK( nu <= 0.5 );
  if( nu >= AGENT_DEFORMATION_INCOMPRESSIBLE_MIN_POISSONS_RATIO ) {/* incompressible */
    while( iters < AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS ) {/* Newton's method to solve a system of non-linear equations (unknowns: a_stretchRatio[0], a_stretchRatio[1], a_stretchRatio[2] */
      VReal vF;
      VReal vTmp;

      /* update vF */

      vF[0] = normalVStress[0] - normalVStress[2] - 2.0 * C1 * ( a_stretchRatio[0] * a_stretchRatio[0] - a_stretchRatio[2] * a_stretchRatio[2] ); 
      vF[1] = normalVStress[1] - normalVStress[2] - 2.0 * C1 * ( a_stretchRatio[1] * a_stretchRatio[1] - a_stretchRatio[2] * a_stretchRatio[2] ); 
      vF[2] = a_stretchRatio[0] * a_stretchRatio[1] * a_stretchRatio[2] - 1.0;

      norm = vF.length();
      if( iters == 0 ) {
	initNorm = norm;
      }

      if( ( iters >= AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS ) || ( norm < initNorm * AGENT_DEFORMATION_NEWTONS_METHOD_EPSILON_INIT_NORM ) || ( norm < E * AGENT_DEFORMATION_NEWTONS_METHOD_EPSILON_E ) || ( norm < AGENT_DEFORMATION_NEWTONS_METHOD_NORM_THRESHOLD ) ) {/* convergence */
	break;
      }
      else {
	/* update matJInv */

	matJInv.aa_val[0][0] = -4.0 * C1 * a_stretchRatio[0];
	matJInv.aa_val[0][1] = 0.0;
	matJInv.aa_val[0][2] = 4.0 * C1 * a_stretchRatio[2];

	matJInv.aa_val[1][0] = 0.0;
	matJInv.aa_val[1][1] = -4.0 * C1 * a_stretchRatio[1];
	matJInv.aa_val[1][2] = 4.0 * C1 * a_stretchRatio[2];

	matJInv.aa_val[2][0] = a_stretchRatio[1] * a_stretchRatio[2];
	matJInv.aa_val[2][1] = a_stretchRatio[2] * a_stretchRatio[0];
	matJInv.aa_val[2][2] = a_stretchRatio[0] * a_stretchRatio[1];

	CHECK( FABS( matJInv.determinant() ) > 0.0 );
	matJInv.invert();

	/* update a_stretchRatio */

	vTmp = matJInv.mult( vF );
	for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
	  a_stretchRatio[dim] -= vTmp[dim];
	  CHECK( a_stretchRatio[dim] > 0.0 );
	}
      }

      iters++;
    }
  }
  else {
    REAL D1 = ( E * nu ) / ( 2.0 * ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) ) );
    while( iters < AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS ) {/* Newton's method to solve a system of non-linear equations (unknowns: a_stretchRatio[0], a_stretchRatio[1], a_stretchRatio[2] */
      VReal vF;
      VReal vTmp;
      REAL norm;
      REAL J;
      REAL I1;

      /* update vF */

      J = a_stretchRatio[0] * a_stretchRatio[1] * a_stretchRatio[2];
      I1 = a_stretchRatio[0] * a_stretchRatio[0] + a_stretchRatio[1] * a_stretchRatio[1] + a_stretchRatio[2] * a_stretchRatio[2];

      vF[0] = normalVStress[0] - normalVStress[2] - 2.0 * C1 * POW( J, -5.0 / 3.0 ) * ( a_stretchRatio[0] * a_stretchRatio[0] - a_stretchRatio[2] * a_stretchRatio[2] );
      vF[1] = normalVStress[1] - normalVStress[2] - 2.0 * C1 * POW( J, -5.0 / 3.0 ) * ( a_stretchRatio[1] * a_stretchRatio[1] - a_stretchRatio[2] * a_stretchRatio[2] );
      vF[2] = normalVStress[2] - 2.0 * C1 * POW( J, -5.0 / 3.0 ) * ( a_stretchRatio[2] * a_stretchRatio[2] - I1 / 3.0 ) - 2.0 * D1 * ( J - 1.0 );

      norm = vF.length();
      if( iters == 0 ) {
	initNorm = norm;
      }

      if( ( iters >= AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS ) || ( norm < initNorm * AGENT_DEFORMATION_NEWTONS_METHOD_EPSILON_INIT_NORM ) || ( norm < E * AGENT_DEFORMATION_NEWTONS_METHOD_EPSILON_E ) || ( norm < AGENT_DEFORMATION_NEWTONS_METHOD_NORM_THRESHOLD ) ) {/* convergence */
	break;
      }
      else {
	/* update matJInv */

	matJInv.aa_val[0][0] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[1] * a_stretchRatio[2] ) * ( a_stretchRatio[0] * a_stretchRatio[0] - a_stretchRatio[2] * a_stretchRatio[2] ) + POW( J, -5.0 / 3.0 ) * ( 2.0 * a_stretchRatio[0] ) );
	matJInv.aa_val[0][1] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[2] * a_stretchRatio[0] ) * ( a_stretchRatio[0] * a_stretchRatio[0] - a_stretchRatio[2] * a_stretchRatio[2] ) + POW( J, -5.0 / 3.0 ) * ( -2.0 * a_stretchRatio[1] ) );
	matJInv.aa_val[0][2] = -2.0 * C1 * ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[0] * a_stretchRatio[1] ) * ( a_stretchRatio[0] * a_stretchRatio[0] - a_stretchRatio[2] * a_stretchRatio[2] );

	matJInv.aa_val[1][0] = -2.0 * C1 * ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[1] * a_stretchRatio[2] ) * ( a_stretchRatio[1] * a_stretchRatio[1] - a_stretchRatio[2] * a_stretchRatio[2] );
	matJInv.aa_val[1][1] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[2] * a_stretchRatio[0] ) * ( a_stretchRatio[1] * a_stretchRatio[1] - a_stretchRatio[2] * a_stretchRatio[2] ) + POW( J, -5.0 / 3.0 ) * ( 2.0 * a_stretchRatio[1] ) );
	matJInv.aa_val[1][2] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[0] * a_stretchRatio[1] ) * ( a_stretchRatio[1] * a_stretchRatio[1] - a_stretchRatio[2] * a_stretchRatio[2] ) + POW( J, -5.0 / 3.0 ) * ( -2.0 * a_stretchRatio[2] ) );

	matJInv.aa_val[2][0] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[1] * a_stretchRatio[2] ) * ( a_stretchRatio[2] * a_stretchRatio[2] - I1 / 3.0 ) + POW( J, -5.0 / 3.0 ) * ( ( -2.0 / 3.0 ) * a_stretchRatio[0] ) ) - 2.0 * D1 * ( a_stretchRatio[1] * a_stretchRatio[2] );
	matJInv.aa_val[2][1] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[2] * a_stretchRatio[0] ) * ( a_stretchRatio[2] * a_stretchRatio[2] - I1 / 3.0 ) + POW( J, -5.0 / 3.0 ) * ( ( -2.0 / 3.0 ) * a_stretchRatio[1] ) ) - 2.0 * D1 * ( a_stretchRatio[2] * a_stretchRatio[0] );
	matJInv.aa_val[2][2] = -2.0 * C1 * ( ( -5.0 / 3.0 ) * POW( J, -8.0 / 3.0 ) * ( a_stretchRatio[0] * a_stretchRatio[1] ) * ( a_stretchRatio[2] * a_stretchRatio[2] - I1 / 3.0 ) + POW( J, -5.0 / 3.0 ) * ( ( 4.0 / 3.0 ) * a_stretchRatio[2] ) ) - 2.0 * D1 * ( a_stretchRatio[0] * a_stretchRatio[1] );

	CHECK( FABS( matJInv.determinant() ) > 0.0 );
	matJInv.invert();

	/* update a_stretchRatio */

	vTmp = matJInv.mult( vF );
	for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
	  a_stretchRatio[dim] -= vTmp[dim];
	  CHECK( a_stretchRatio[dim] > 0.0 );
	}
      }

      iters++;
    }
  }

  if( iters >= AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS ) {
    WARNING( "id=" << junctionData.getCurId() << ", max iteration(" << AGENT_DEFORMATION_NEWTONS_METHOD_MAX_ITERS << ") count reached, initNorm=" << initNorm << " norm=" << norm );
  } 
#endif
  return;
}
