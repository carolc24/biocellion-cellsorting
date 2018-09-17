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

VReal computeOvlpForce( const SpAgentState& state0, const SpAgentState& state1, const Ellipsoid& e0, const Ellipsoid& e1, VReal& minPotentialVPointOnE0, VReal& minPotentialVPointOnE0FromE1, VReal& minPotentialVPointOnE1, VReal& minPotentialVPointOnE1FromE0);

void computeMinPotentialPoints( const Ellipsoid& e0, const Ellipsoid& e1, const VReal& vDir, const REAL& dist, VReal& minPotentialVPointOnE0, VReal& minPotentialVPointOnE1);

void ModelRoutine::initJunctionSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const UBEnv& ubEnv0, const VIdx& vIdx1, const SpAgent& spAgent1, const UBEnv& ubEnv1, const VReal& vDir/* unit direction vector from spAgent1 to spAgent0 */, const REAL dist, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */
) {
	/* MODEL START */

	link = false;

	/* MODEL END */

	return;
}

void ModelRoutine::computeMechIntrctSpAgent( const S32 iter, const VIdx& vIdx0, const SpAgent& spAgent0, const UBEnv& ubEnv0, const VIdx& vIdx1, const SpAgent& spAgent1, const UBEnv& ubEnv1, const VReal& vDir/* unit direction vector from spAgent1 to spAgent0 */, const REAL dist, MechIntrctData& mechIntrctData0, MechIntrctData& mechIntrctData1, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */, BOOL& unlink ) {
  /* MODEL START */

  agentType_t type0 = spAgent0.state.getType();
  agentType_t type1 = spAgent1.state.getType();
  
  link = false;
  unlink = false;

  VReal vScale0;
  VReal vScale1;
  Ellipsoid e0;
  Ellipsoid e1;

  VReal vForce;
  VReal vMoment0;
  VReal vMoment1;

  VReal equivMinPotentialVPointOnE0;
  VReal equivMinPotentialVPointOnE0FromE1;
  VReal equivMinPotentialVPointOnE1;
  VReal equivMinPotentialVPointOnE1FromE0;
  S32 idx0;
  S32 idx1;

  BOOL closerThanEquilibriumDist;

  vScale0[0] = spAgent0.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_A );
  vScale0[1] = spAgent0.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_B );
  vScale0[2] = spAgent0.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_C );
  e0.qRot.set( spAgent0.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A ), spAgent0.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B ), spAgent0.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C ), spAgent0.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D ));

  vScale1[0] = spAgent1.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_A );
  vScale1[1] = spAgent1.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_B );
  vScale1[2] = spAgent1.state.getModelReal( AGENT_STATE_REAL_ELLIPSOID_C );
  e1.qRot.set( spAgent1.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_A ), spAgent1.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_B ), spAgent1.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_C ), spAgent1.state.getModelReal( AGENT_STATE_REAL_ROTATIONAL_QUATERNION_D ));

  vForce = VReal::ZERO;
  vMoment0 = VReal::ZERO;
  vMoment1 = VReal::ZERO;

  e0.vScale = vScale0;
  e1.vScale = vScale1;

  computeMinPotentialPoints( e0, e1, vDir, dist, equivMinPotentialVPointOnE0, equivMinPotentialVPointOnE1);
  equivMinPotentialVPointOnE0FromE1 = equivMinPotentialVPointOnE0 + vDir * dist;
  equivMinPotentialVPointOnE1FromE0 = equivMinPotentialVPointOnE1 - vDir * dist;

  if (e0.overlaps( e1, vDir * -1.0, dist) == true) {
    //printf("Overlap detected\n");
    VReal vPos0;
    VReal vPos1;

    vForce = computeOvlpForce( spAgent0.state, spAgent1.state, e0, e1, equivMinPotentialVPointOnE0, equivMinPotentialVPointOnE0FromE1, equivMinPotentialVPointOnE1, equivMinPotentialVPointOnE1FromE0);

     vPos0 = ( equivMinPotentialVPointOnE0 + equivMinPotentialVPointOnE1FromE0) * 0.5;
     vMoment0 = VReal::crossProduct( vPos0, vForce);

     vPos1 = ( equivMinPotentialVPointOnE1 + equivMinPotentialVPointOnE0FromE1) * 0.5;
     vMoment1 = VReal::crossProduct( vPos1, vForce*-1.0);

     closerThanEquilibriumDist = true;
  } else {
    closerThanEquilibriumDist = false;
  }

  mechIntrctData0.setModelReal( AGENT_MECH_REAL_FORCE_X, vForce[0] );
  mechIntrctData0.setModelReal( AGENT_MECH_REAL_FORCE_Y, vForce[1] );
  mechIntrctData0.setModelReal( AGENT_MECH_REAL_FORCE_Z, vForce[2] );
  mechIntrctData0.setModelReal( AGENT_MECH_REAL_MOMENT_X, vMoment0[0] );
  mechIntrctData0.setModelReal( AGENT_MECH_REAL_MOMENT_Y, vMoment0[1] );
  mechIntrctData0.setModelReal( AGENT_MECH_REAL_MOMENT_Z, vMoment0[2] );

  mechIntrctData1.setModelReal( AGENT_MECH_REAL_FORCE_X, vForce[0]*-1.0 );
  mechIntrctData1.setModelReal( AGENT_MECH_REAL_FORCE_Y, vForce[1]*-1.0 );
  mechIntrctData1.setModelReal( AGENT_MECH_REAL_FORCE_Z, vForce[2]*-1.0 );
  mechIntrctData1.setModelReal( AGENT_MECH_REAL_MOMENT_X, vMoment1[0] );
  mechIntrctData1.setModelReal( AGENT_MECH_REAL_MOMENT_Y, vMoment1[1] );
  mechIntrctData1.setModelReal( AGENT_MECH_REAL_MOMENT_Z, vMoment1[2] );
  
  if( vForce.lengthSquare() > 0.0 ) {
    e0.vScale = vScale0;
    e1.vScale = vScale1;

    //if( A_AGENT_DEFORMABLE[type0] == true ) {/* deformation */
      VReal normalVForce;
      VReal normalVDir;
      VReal intersectingVPos;

      Quaternion qTmp;

      if( vPos0.lengthSquare() > 0.0 ) {
	intersectingVPos = e0.getIntersectingPointOn( VReal::normalize( vPos0 ) );

	CHECK( FABS( 1.0 - e0.qRot.norm() ) < REAL_EPSILON * 1e1 );/* to assure q^-1 == q^* */
	qTmp = Quaternion::qStarpq( Quaternion( 0.0, intersectingVPos ), e0.qRot );/* now in the body-fixed frame */
	intersectingVPos = qTmp.getImg();

	normalVDir[0] = ( 2.0 * intersectingVPos[0] ) / ( e0.vScale[0] * e0.vScale[0] );
	normalVDir[1] = ( 2.0 * intersectingVPos[1] ) / ( e0.vScale[1] * e0.vScale[1] );
	normalVDir[2] = ( 2.0 * intersectingVPos[2] ) / ( e0.vScale[2] * e0.vScale[2] );
	normalVDir = VReal::normalize( normalVDir );
	CHECK( FABS( 1.0 - normalVDir.length() ) < REAL_EPSILON * 1e1 );

	normalVForce = vForce;
	CHECK( FABS( 1.0 - e0.qRot.norm() ) < REAL_EPSILON * 1e1 );/* to assure q^-1 == q^* */
	qTmp = Quaternion::qStarpq( Quaternion( 0.0, normalVForce ), e0.qRot );/* now in the body-fixed frame */
	normalVForce = qTmp.getImg();

	normalVForce = normalVDir * normalVDir.dotProduct( normalVForce );

	if( intersectingVPos[0] < 0.0 ) {/* low side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_X, normalVForce[0] );
	}
	else if( intersectingVPos[0] > 0.0 ) {/* high side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_X, normalVForce[0] );
	}

	if( intersectingVPos[1] < 0.0 ) {/* low side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Y, normalVForce[1] );
	}
	else if( intersectingVPos[1] > 0.0 ) {/* high side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Y, normalVForce[1] );
	}

	if( intersectingVPos[2] < 0.0 ) {/* low side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Z, normalVForce[2] );
	}
	else if( intersectingVPos[2] > 0.0 ) {/* high side */
	  mechIntrctData0.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Z, normalVForce[2] );
	}
      }
      //}

      //if( A_AGENT_DEFORMABLE[type1] == true ) {/* deformation */
      //VReal normalVForce;
      //VReal normalVDir;
      //VReal intersectingVPos;

      //Quaternion qTmp;

      if( vPos1.lengthSquare() > 0.0 ) {
	intersectingVPos = e1.getIntersectingPointOn( VReal::normalize( vPos1 ) );

	CHECK( FABS( 1.0 - e1.qRot.norm() ) < REAL_EPSILON * 1e1 );/* to assure q^-1 == q^* */
	qTmp = Quaternion::qStarpq( Quaternion( 0.0, intersectingVPos ), e1.qRot );/* now in the body-fixed frame */
	intersectingVPos = qTmp.getImg();

	normalVDir[0] = ( 2.0 * intersectingVPos[0] ) / ( e1.vScale[0] * e1.vScale[0] );
	normalVDir[1] = ( 2.0 * intersectingVPos[1] ) / ( e1.vScale[1] * e1.vScale[1] );
	normalVDir[2] = ( 2.0 * intersectingVPos[2] ) / ( e1.vScale[2] * e1.vScale[2] );
	normalVDir = VReal::normalize( normalVDir );
	CHECK( FABS( 1.0 - normalVDir.length() ) < REAL_EPSILON * 1e1 );

	normalVForce = vForce * -1.0;
	CHECK( FABS( 1.0 - e1.qRot.norm() ) < REAL_EPSILON * 1e1 );/* to assure q^-1 == q^* */
	qTmp = Quaternion::qStarpq( Quaternion( 0.0, normalVForce ), e1.qRot );/* now in the body-fixed frame */
	normalVForce = qTmp.getImg();

	normalVForce = normalVDir * normalVDir.dotProduct( normalVForce );

	if( intersectingVPos[0] < 0.0 ) {/* low side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_X, normalVForce[0] );
	}
	else if( intersectingVPos[0] > 0.0 ) {/* high side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_X, normalVForce[0] );
	}

	if( intersectingVPos[1] < 0.0 ) {/* low side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Y, normalVForce[1] );
	}
	else if( intersectingVPos[1] > 0.0 ) {/* high side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Y, normalVForce[1] );
	}

	if( intersectingVPos[2] < 0.0 ) {/* low side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_LOW_Z, normalVForce[2] );
	}
	else if( intersectingVPos[2] > 0.0 ) {/* high side */
	  mechIntrctData1.incModelReal( AGENT_MECH_REAL_BODY_FIXED_NORMAL_FORCE_HIGH_Z, normalVForce[2] );
	}
      }

      //}
  }
/* MODEL END */

  return;
  
}

void computeMinPotentialPoints( const Ellipsoid& e0, const Ellipsoid& e1, const VReal& vDir, const REAL& dist, VReal& minPotentialVPointOnE0, VReal& minPotentialVPointOnE1) {
  minPotentialVPointOnE0 = e0.getMinPotentialPointOn( e1, vDir * -dist );
  minPotentialVPointOnE1 = e1.getMinPotentialPointOn( e0, vDir * dist );

  return;
}

VReal computeOvlpForce( const SpAgentState& state0, const SpAgentState& state1, const Ellipsoid& e0, const Ellipsoid& e1, VReal& minPotentialVPointOnE0, VReal& minPotentialVPointOnE0FromE1, VReal& minPotentialVPointOnE1, VReal& minPotentialVPointOnE1FromE0) {
  VReal a_principalAxisVDir0[2];
  REAL a_principalCurvature0[2];

  VReal a_principalAxisVDir1[2];
  REAL a_principalCurvature1[2];

  REAL oneOverR0Prime;
  REAL oneOverR0DoublePrime;
  REAL oneOverR1Prime;
  REAL oneOverR1DoublePrime;
  VReal principalVDir0;
  VReal principalVDir1;
  REAL cosAlpha;
  REAL cos2Alpha;
  REAL APlusB;
  REAL BMinusA;
  REAL A;
  REAL B;
  REAL Re;

  REAL E0;
  REAL E1;
  REAL nu0;
  REAL nu1;
  REAL EStar;

  VReal vN;
  REAL sigmaN;
  REAL correction;
  REAL f;

  e0.computePrincipalCurvatures( minPotentialVPointOnE0, a_principalAxisVDir0, a_principalCurvature0 );
  e1.computePrincipalCurvatures( minPotentialVPointOnE1, a_principalAxisVDir1, a_principalCurvature1 );

  oneOverR0Prime = FABS( a_principalCurvature0[0] );
  oneOverR0DoublePrime = FABS( a_principalCurvature0[1] );
  oneOverR1Prime = FABS( a_principalCurvature1[0] );
  oneOverR1DoublePrime = FABS( a_principalCurvature1[1] );

  principalVDir0 = a_principalAxisVDir0[0] + a_principalAxisVDir0[1];
  principalVDir1 = a_principalAxisVDir1[0] + a_principalAxisVDir1[1];
  principalVDir0 = VReal::normalize(principalVDir0);
  principalVDir1 = VReal::normalize(principalVDir1);

  APlusB = (oneOverR0Prime + oneOverR0DoublePrime + oneOverR1Prime + oneOverR1DoublePrime) * 0.5;
  cosAlpha = principalVDir0.dotProduct(principalVDir1);
  cos2Alpha = 2.0 * cosAlpha * cosAlpha - 1.0;
  if (cos2Alpha > 1.0 - REAL_EPSILON * 1e1) {
    BMinusA = FABS((oneOverR0Prime - oneOverR0DoublePrime) + (oneOverR1Prime - oneOverR1DoublePrime))*0.5;
  } else if (cos2Alpha > (1.0 - REAL_EPSILON * 1e1) * -1.0 ) {
    BMinusA = FABS((oneOverR0Prime - oneOverR0DoublePrime) - (oneOverR1Prime - oneOverR1DoublePrime))*0.5;
  } else {
    BMinusA = SQRT((oneOverR0Prime - oneOverR0DoublePrime) * (oneOverR0Prime - oneOverR0DoublePrime) + (oneOverR1Prime - oneOverR1DoublePrime) * (oneOverR1Prime - oneOverR1DoublePrime) + 2.0 * (oneOverR0Prime - oneOverR0DoublePrime) * (oneOverR1Prime - oneOverR1DoublePrime) * cos2Alpha) * 0.5;
  }

  CHECK(APlusB >= BMinusA);
  B = (APlusB + BMinusA) * 0.5;
  A = APlusB - B;
  CHECK(A*B > 0.0);
  Re = (1.0 / SQRT(A*B)) * 0.5;

  E0 = state0.getModelReal(AGENT_STATE_REAL_YOUNGS_MODULUS);
  E1 = state1.getModelReal(AGENT_STATE_REAL_YOUNGS_MODULUS);
  nu0 = state0.getModelReal(AGENT_STATE_REAL_POISSONS_RATIO);
  nu1 = state1.getModelReal(AGENT_STATE_REAL_POISSONS_RATIO);

  EStar = 1.0 / (( 1.0 - nu0 * nu0) / E0 + (1.0 - nu1 * nu1) / E1);

  vN = (((minPotentialVPointOnE1FromE0 - minPotentialVPointOnE0) +( minPotentialVPointOnE1 - minPotentialVPointOnE0FromE1))*0.5);
  sigmaN=vN.length();
  vN = VReal::normalize(vN);
  correction = POW(1.0 - POW(POW(B/A, 0.0684) - 1.0, 1.531), -3.0/2.0);

  f = (4.0/3.0) * EStar * SQRT(Re) * POW(sigmaN, 3.0/2.0) * correction;

  return vN * f;
}
#endif

