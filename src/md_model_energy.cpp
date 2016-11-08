/***************************************************************************
* Copyright (C) 2013-2014 Alexander V. Popov.
* 
* This file is part of Amateur Modeling of Biopolymers: Restoration, 
* Optimization, Solvation & Initial Analysis (AMBROSIA) source code.
* 
* AMBROSIA source code is free software; you can redistribute it and/or 
* modify it under the terms of the GNU General Public License as 
* published by the Free Software Foundation; either version 2 of 
* the License, or (at your option) any later version.
* 
* AMBROSIA source code is distributed in the hope that it will be 
* useful, but WITHOUT ANY WARRANTY; without even the implied 
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
* See the GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software 
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
***************************************************************************/
#include "md_main.h"
#include "md_model.h"
#include "md_math.h"

cModel *cModel :: g_pModel = NULL;

//-----------------------------------
// Energy calculation & Optimization
//-----------------------------------
vec_t cModel :: CalcBondEnergy( const cBond &source, const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flBondLength;
	vec3_t vecDelta;

	// Get pointers to positions
	assert( pPositions != NULL );
	const vec_t *pPosI = *(pPositions + source.mAtomIndices[0]);
	const vec_t *pPosJ = *(pPositions + source.mAtomIndices[1]);

	// Calculate bond length
	Vec3Sub( pPosI, pPosJ, vecDelta );
	flBondLength = Vec3Len( vecDelta );

	// Calculate energy term
	// E = Kh*(b-b0)^2
	flEnergy = source.mHarmonicConstant * SQR( flBondLength - source.mLength );

	// Calculate forces on atoms
	// Fi = -dE/dRi = -2*Kh*(b-b0)*(Rij/|Rij|)
	// Fj = -dE/dRj = 2*Kh*(b-b0)*(Rij/|Rij|)
	if ( pForces ) {
		vec3_t vecForce;

		vec_t *pForceI = *(pForces + source.mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + source.mAtomIndices[1]);

		Vec3Scale( vecDelta, 2 * source.mHarmonicConstant * ( 1 - source.mLength / flBondLength ), vecForce );
		Vec3Sub( pForceI, vecForce, pForceI );
		Vec3Add( pForceJ, vecForce, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
	}

	return flEnergy;
}

vec_t cModel :: CalcAngleEnergy( const cAngle &source, const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flTheta, flDot;
	vec3_t vecDelta1, vecDelta2;
	vec_t flLenIJInv, flLenKJInv;

	// Get pointers to positions
	assert( pPositions != NULL );
	const vec_t *pPosI = *(pPositions + source.mAtomIndices[0]);
	const vec_t *pPosJ = *(pPositions + source.mAtomIndices[1]);
	const vec_t *pPosK = *(pPositions + source.mAtomIndices[2]);

	// Calculate angle
	Vec3Sub( pPosI, pPosJ, vecDelta1 );
	Vec3Sub( pPosK, pPosJ, vecDelta2 );
	Vec3Nrm2( vecDelta1, flLenIJInv );
	Vec3Nrm2( vecDelta2, flLenKJInv );
	flDot = Vec3Dot( vecDelta1, vecDelta2 );
	CLAMP( flDot, -1, 1 );

	flTheta = acos( flDot );

	// Calculate energy term
	// E = Kh*(theta-theta0)^2
	vec_t angleDiff = AngleDiff( flTheta, source.mTheta );
	flEnergy = source.mHarmonicConstant * SQR( angleDiff );

	// Calculate forces on atoms
	// S = 2*Kh*(theta-theta0)*(1/sqrt(1-flDot*flDot))
	// Fi = -dE/dRi = (S/|Rij|)*[R(n)kj - D*R(n)ij]
	// Fk = -dE/dRk =(S/|Rkj|)*[R(n)ij - D*R(n)kj]
	// Fj = -dE/dRj = -Fi-Fk
	if ( pForces ) {
		vec3_t vecForceI, vecForceK;
		vec_t scalar;

		vec_t *pForceI = *(pForces + source.mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + source.mAtomIndices[1]);
		vec_t *pForceK = *(pForces + source.mAtomIndices[2]);

		scalar = 2 * source.mHarmonicConstant * angleDiff / ( sqrt( 1 - flDot*flDot ) );

		Vec3MA( vecDelta2, -flDot, vecDelta1, vecForceI );
		Vec3Scale( vecForceI, scalar * flLenIJInv, vecForceI );
		Vec3MA( vecDelta1, -flDot, vecDelta2, vecForceK );
		Vec3Scale( vecForceK, scalar * flLenKJInv, vecForceK );

		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Add( pForceK, vecForceK, pForceK );
		Vec3Sub( pForceJ, vecForceI, pForceJ );
		Vec3Sub( pForceJ, vecForceK, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
		CHECK_NAN3( pForceK );
	}

	return flEnergy;
}


vec_t cModel :: CalcImproperEnergy( const cImproper &source, const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flXi, flDot;
	vec3_t vecIJ, vecKJ, vecKL, vecMJ, vecNK;
	vec_t flLenMJ, flLenNK;

	// Get pointers to positions
	// NOTE: improper center in config is atom J, but here we assume center to be atom I
	assert( pPositions != NULL );
	const vec_t *pPosI = *(pPositions + source.mAtomIndices[1]);
	const vec_t *pPosJ = *(pPositions + source.mAtomIndices[0]);
	const vec_t *pPosK = *(pPositions + source.mAtomIndices[2]);
	const vec_t *pPosL = *(pPositions + source.mAtomIndices[3]);

	// Calculate torsion
	Vec3Sub( pPosI, pPosJ, vecIJ );
	Vec3Sub( pPosK, pPosJ, vecKJ );
	Vec3Sub( pPosK, pPosL, vecKL );
	Vec3Cross( vecIJ, vecKJ, vecMJ );
	Vec3Cross( vecKJ, vecKL, vecNK );
	flLenMJ = Vec3Nrm( vecMJ );
	flLenNK = Vec3Nrm( vecNK );
	flDot = Vec3Dot( vecMJ, vecNK );
	CLAMP( flDot, -1, 1 );

	flXi = acos( flDot );

	if ( Vec3Dot( vecIJ, vecNK ) < 0 )
		flXi = -flXi;

	// Calculate energy term
	// E = Kh*(xi-xi0)^2
	vec_t angleDiff = AngleDiff( flXi, source.mXi );
	flEnergy = source.mHarmonicConstant * SQR( angleDiff );

	// Calculate forces on atoms
	// Fi = -dE/dRi = -2*Kh*(xi-xi0)*(|Rkj|/|Rmj|^2)*Rmj
	// Fj = -dE/dRj = [(Rij.Rkj)/|Rkj|^2 - 1]*Fi - [(Rkl.Rkj)/|Rkj|^2]*Fl
	// Fk = -Fi-Fj-Fl
	// Fl = -dE/dRl = 2*Kh*(xi-xi0)*(|Rkj|/|Rnk|^2)*Rnk
	if ( pForces ) {
		vec3_t vecForceI, vecForceJ, vecForceL;
		vec_t flLenKJ, scalar;

		vec_t *pForceI = *(pForces + source.mAtomIndices[1]);
		vec_t *pForceJ = *(pForces + source.mAtomIndices[0]);
		vec_t *pForceK = *(pForces + source.mAtomIndices[2]);
		vec_t *pForceL = *(pForces + source.mAtomIndices[3]);

		flLenKJ = Vec3Len( vecKJ );
		scalar = 2 * source.mHarmonicConstant * angleDiff;
		flLenMJ *= flLenMJ;
		flLenNK *= flLenNK;

		Vec3Scale( vecMJ, -scalar * ( flLenKJ / flLenMJ ), vecForceI );
		Vec3Scale( vecNK, scalar * ( flLenKJ / flLenNK ), vecForceL );

		flLenKJ *= flLenKJ;
		Vec3Scale( vecForceI, Vec3Dot( vecIJ, vecKJ ) / flLenKJ - 1, vecForceJ );
		Vec3MA( vecForceJ, -Vec3Dot( vecKL, vecKJ ) / flLenKJ, vecForceL, vecForceJ );

		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Add( pForceL, vecForceL, pForceL );
		Vec3Add( pForceJ, vecForceJ, pForceJ );
		Vec3Sub( pForceK, vecForceI, pForceK );
		Vec3Sub( pForceK, vecForceL, pForceK );
		Vec3Sub( pForceK, vecForceJ, pForceK );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
		CHECK_NAN3( pForceK );
	}

	return flEnergy;
}

vec_t cModel :: CalcTorsionEnergy( const cTorsion &source, const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flScale, flDot, flCosN;
	vec3_t vecIJ, vecKJ, vecKL, vecIM, vecNL;
	vec_t flLenKJSq, flLenIM, flLenNL;
	vec_t flSine;
	cTorsionHarm *pHarm;
	
	// Get pointers to positions
	assert( pPositions != NULL );
	const vec_t *pPosI = *(pPositions + source.mAtomIndices[0]);
	const vec_t *pPosJ = *(pPositions + source.mAtomIndices[1]);
	const vec_t *pPosK = *(pPositions + source.mAtomIndices[2]);
	const vec_t *pPosL = *(pPositions + source.mAtomIndices[3]);

	Vec3Sub( pPosI, pPosJ, vecIJ );
	Vec3Sub( pPosK, pPosJ, vecKJ );
	Vec3Sub( pPosK, pPosL, vecKL );
	flLenKJSq = Vec3Dot( vecKJ, vecKJ );
	Vec3MA( vecIJ, -Vec3Dot( vecIJ, vecKJ ) / flLenKJSq, vecKJ, vecIM );
	Vec3MA( vecKL, -Vec3Dot( vecKL, vecKJ ) / flLenKJSq, vecKJ, vecNL );

	flLenIM = Vec3Nrm( vecIM );
	flLenNL = Vec3Nrm( vecNL );
	flDot = -Vec3Dot( vecIM, vecNL );
	CLAMP( flDot, -1, 1 );

	flEnergy = 0;

	// Calculate harmonics
	pHarm = source.mpHarmonic;
	while ( pHarm ) {
		if ( pHarm->mHarmonicConstant != 0 ) {
			switch ( pHarm->mPeriodicity ) {
			case 1: flCosN = flDot; break;
			case 2: flCosN = flDot*flDot*2-1; break;
			case 3: flCosN = flDot*(flDot*flDot*4-3); break;
			case 4: flCosN = flDot*flDot*8*(flDot*flDot-1)+1; break;
			case 5: flCosN = 4*flDot*flDot*flDot*(flDot*flDot*4-5)+flDot*5; break;
			case 6: flCosN = 16*flDot*flDot*flDot*flDot*(flDot*flDot*2-3)+flDot*flDot*18-1; break;
			default: assume_0;
			}

			// Get scale factor
			flScale = source.mBarrierScale ? pHarm->mHarmonicScale : 1;

			// Get sine
			if ( pHarm->mPhaseSin != 0 )
				flSine = sqrt( 1 - flCosN*flCosN ) * pHarm->mPhaseSin;
			else
				flSine = 0;

			// Calculate energy term
			// E = Kh*(1 + cos(n*omega - phase)) = Kh*(1 + cos(n*omega) * cos(phase) + sin(n*omega) * sin(phase)))
			flEnergy += pHarm->mHarmonicConstant * flScale * ( 1 + pHarm->mPhaseCos * flCosN + flSine );
		}

		// Advance to the next harmonic
		pHarm = pHarm->mpNext;
	}

	// Calculate forces on atoms
	// Fi = -dE/dRi = Kh*dcos(n*omega-phi)*[Rnl(n)+Rim(n)*cos(omega)]/|Rim|
	// Fj = -dE/dRj = [(Rij.Rkj)/(|Rkj|^2)-1]*Fi - [(Rkl.Rkj)/|Rkj|^2]*Fl
	// Fk = -Fi-Fj-Fl
	// Fl = -dE/dRl = -Kh*dcos(n*omega-phi)[Rim(n)+Rnl(n)*cos(omega)]/|Rnl|
	if ( pForces ) {
		vec3_t vecForceI, vecForceJ, vecForceL;
		vec_t flSinN, flDCosN, scalar = 0;

		vec_t *pForceI = *(pForces + source.mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + source.mAtomIndices[1]);
		vec_t *pForceK = *(pForces + source.mAtomIndices[2]);
		vec_t *pForceL = *(pForces + source.mAtomIndices[3]);

		// Calculate harmonics
		pHarm = source.mpHarmonic;
		while ( pHarm ) {
			if ( pHarm->mHarmonicConstant != 0 ) {
				switch ( pHarm->mPeriodicity ) {
				case 1: flDCosN = 1; break;
				case 2: flDCosN = flDot*4; break;
				case 3: flDCosN = flDot*flDot*12-3; break;
				case 4: flDCosN = 16*flDot*(flDot*flDot*2-1); break;
				case 5: flDCosN = 20*flDot*flDot*(flDot*flDot*4-3)+5; break;
				case 6: flDCosN = 192*flDot*flDot*flDot*(flDot*flDot-1)+flDot*36; break;
				default: assume_0;
				}
				flScale = source.mBarrierScale ? pHarm->mHarmonicScale : 1;
				flSine = 0;
				if ( pHarm->mPhaseSin != 0 ) {
					switch ( pHarm->mPeriodicity ) {
					case 1: flCosN = flDot; break;
					case 2: flCosN = flDot*flDot*2-1; break;
					case 3: flCosN = flDot*(flDot*flDot*4-3); break;
					case 4: flCosN = flDot*flDot*8*(flDot*flDot-1)+1; break;
					case 5: flCosN = 4*flDot*flDot*flDot*(flDot*flDot*4-5)+flDot*5; break;
					case 6: flCosN = 16*flDot*flDot*flDot*flDot*(flDot*flDot*2-3)+flDot*flDot*18-1; break;
					default: assume_0;
					}
					flSinN = sqrt( 1 - flCosN*flCosN );
					if ( flSinN != 0 ) flSine = flCosN * pHarm->mPhaseSin / flSinN;
				}
				scalar += pHarm->mHarmonicConstant * flScale * flDCosN * ( pHarm->mPhaseCos - flSine );
				CHECK_NAN( scalar );
			}

			// Advance to the next harmonic
			pHarm = pHarm->mpNext;
		}

		Vec3MA( vecNL, flDot, vecIM, vecForceI );
		Vec3Scale( vecForceI, scalar / flLenIM, vecForceI );
		Vec3MA( vecIM, flDot, vecNL, vecForceL );
		Vec3Scale( vecForceL, -scalar / flLenNL, vecForceL );

		Vec3Scale( vecForceI, Vec3Dot( vecIJ, vecKJ ) / flLenKJSq - 1, vecForceJ );
		Vec3MA( vecForceJ, -Vec3Dot( vecKL, vecKJ ) / flLenKJSq, vecForceL, vecForceJ );

		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Add( pForceL, vecForceL, pForceL );
		Vec3Add( pForceJ, vecForceJ, pForceJ );
		Vec3Sub( pForceK, vecForceI, pForceK );
		Vec3Sub( pForceK, vecForceL, pForceK );
		Vec3Sub( pForceK, vecForceJ, pForceK );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
		CHECK_NAN3( pForceK );
		CHECK_NAN3( pForceL );
	}

	return flEnergy;
}

vec_t cModel :: CalcSpecialEnergy( const vec4_t *pPositions, vec4_t *pForces, vec_t *pPRPart, vec_t *pDRPart )
{
	vec_t flPosRestrEnergy = 0;
	vec_t flDistRestrEnergy = 0;
	vec3_t vecDelta;
	vec_t flValue;

	assert( pPositions != NULL );

	const vec4_t *pPos = pPositions;
	vec4_t *pForce = pForces;

	// position restraining
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++pPos, ++pForce ) {
		if ( it->mFlags & AF_RESTRAINED ) {
			// Calculate square of offset
			Vec3Sub( *pPos, it->mOriginalPosition, vecDelta );
			flValue = Vec3LenSq( vecDelta );

			// Calculate energy term
			// E = Kh*(r-r0)^2
			flPosRestrEnergy += it->mPosRestraintHarmConst * flValue;

			// Calculate force on atom
			// F = -dE/dR = -2*Kh*(r-r0)
			if ( pForces ) {
				Vec3MA( *pForce, -2 * it->mPosRestraintHarmConst, vecDelta, *pForce );
				CHECK_NAN3( *pForce );
			}
		}
	}

	// distance restraining
	for ( DistRestArray::iterator it = m_distRestraints.begin(); it != m_distRestraints.end(); ++it ) {
		const vec_t *pPosI = *(pPositions + it->mAtomIndices[0]);
		const vec_t *pPosJ = *(pPositions + it->mAtomIndices[1]);

		// Calculate bond length
		Vec3Sub( pPosI, pPosJ, vecDelta );
		flValue = Vec3LenSq( vecDelta );

		// Calculate energy term
		// E = Kh*(b^2-b0^2)^2
		flDistRestrEnergy += it->mHarmonicConstant * SQR( flValue - it->mLengthSquared );

		// Calculate forces on atoms
		// Fi = -dE/dRi = -4*Kh*(b^2-b0^2)*Rij
		// Fj = -dE/dRj = 4*Kh*(b^2-b0^2)*Rij
		if ( pForces ) {
			vec3_t vecForce;

			vec_t *pForceI = *(pForces + it->mAtomIndices[0]);
			vec_t *pForceJ = *(pForces + it->mAtomIndices[1]);

			Vec3Scale( vecDelta, 4 * it->mHarmonicConstant * ( flValue - it->mLengthSquared ), vecForce );
			Vec3Sub( pForceI, vecForce, pForceI );
			Vec3Add( pForceJ, vecForce, pForceJ );
			CHECK_NAN3( pForceI );
			CHECK_NAN3( pForceJ );
		}
	}

	if ( pPRPart )	*pPRPart	+= flPosRestrEnergy;
	if ( pDRPart )	*pDRPart	+= flDistRestrEnergy;

	return flPosRestrEnergy + flDistRestrEnergy;
}

vec_t cModel :: CalcVdWEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flSwitch_F = 1, flSwitch_dF = 0;

	// 1-4 scale was taken from AMBER (J Am Chem Soc 1995, V. 117, N 19)
#if !defined(BIOPASED_COMPAT_ENERGY)
	const vec_t vdw14Scale = vec_t( 0.5 );
#else
	const vec_t vdw14Scale = vec_t( 0.25 );
#endif

	// Calculate square of distance
	vec_t flDistanceSquared = SQR( flDistance );

	// Soft cut-off potential using switch function
	// See http://www.ks.uiuc.edu/Research/namd/1.5/ug/node52.html
	if ( flDistanceSquared > m_flSwitchR1Squared ) {
		vec_t invDenom = vec_t( 1.0 ) / CUBE( m_flCutoffR1Squared - m_flSwitchR1Squared );

		// Calculate switch function
		// F = [(cutoff^2 - dist^2)^2 * (cutoff^2 + 2*dist^2 - 3*switch^2)] / (cutoff^2 - switch^2)^3
		flSwitch_F = SQR( m_flCutoffR1Squared - flDistanceSquared ) * ( m_flCutoffR1Squared + 2 * flDistanceSquared - 3 * m_flSwitchR1Squared ) * invDenom;

		// Calculate switch function derivative
		// dF = 12*dist*(cutoff^2 - dist^2)*(switch^2 - dist^2) / (cutoff^2 - switch^2)^3
		flSwitch_dF = 12 * flDistance * ( m_flCutoffR1Squared - flDistanceSquared ) * ( m_flSwitchR1Squared - flDistanceSquared ) * invDenom;
	}

	// Calculate van der Waals potential in kcal/mol
	// E = (C12/r^6 - C6)/r^6
	vec_t r6i = SIXTH( flDistanceInv );
	vec_t C12_r6i = pPair->mpVdWPair->mC12 * r6i;
	vec_t C6 = pPair->mpVdWPair->mC6;

	// Scale coefficients for 1-4 neighbours
	if ( pPair->mFlags & NBPF_NEIGHBOURS_14 ) {
		C12_r6i *= vdw14Scale;
		C6 *= vdw14Scale;
	}

	flEnergy = ( C12_r6i - C6 ) * r6i;

	if ( pForces ) {
		vec3_t vecForceI;

		vec_t *pForceI = *(pForces + pPair->mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + pPair->mAtomIndices[1]);

		// Calculate van der Waals forces on atoms
		// Fi = -dE/dRi = [2*C12/r^6 - C6]*6*(Rij/|Rij|^8)
		// Fj = -Fi
		if ( flDistanceSquared > m_flSwitchR1Squared )
			Vec3Scale( pvecDelta, ( 2 * C12_r6i - C6 ) * 6 * r6i * flDistanceInv * flSwitch_F - flEnergy * flSwitch_dF, vecForceI );
		else
			Vec3Scale( pvecDelta, ( 2 * C12_r6i - C6 ) * 6 * r6i * flDistanceInv, vecForceI );

		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Sub( pForceJ, vecForceI, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
	}

	return flEnergy * flSwitch_F;
}

vec_t cModel :: CalcElectrostaticEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces )
{
	vec_t flEnergy;
	vec_t flDistScale = 1;

	// Electrostatic energy constant is used to convert internal units of Coulomb energy to kcal/mol
	// eec = [charge of electron, CGSE]^2 * (calories in Joule) * NAvogadro * 10.0
	const vec_t eec = vec_t( SQR( 4.80320427 ) * 0.238846 * 6.02214129 * 10.0 );

	// 1-4 scale was taken from AMBER (J Am Chem Soc 1995, V. 117, N 19)
#if !defined(BIOPASED_COMPAT_ENERGY)
	const vec_t es14Scale = vec_t( 1.0 / 1.2 );
#else
	const vec_t es14Scale = vec_t( 0.125 );
#endif

	// Calculate square of distance
	vec_t flDistanceSquared = SQR( flDistance );

	// Get charges
	vec_t qI = *(m_charges + pPair->mAtomIndices[0]);
	vec_t qJ = *(m_charges + pPair->mAtomIndices[1]);
	vec_t qIqJ = eec * qI * qJ;

	// Scale coefficients for 1-4 neighbours
	if ( pPair->mFlags & NBPF_NEIGHBOURS_14 ) {
		qIqJ *= es14Scale;
	}

	// Calculate Coulomb potential in kcal/mol in reaction field terms (J Chem Phys 1995, V. 102, N 13)
	// E = eec*qi*qj*[1/r + A*r^2 + B]
	flEnergy = qIqJ * ( flDistanceInv + m_flCRF_A * flDistanceSquared + m_flCRF_B );

	if ( m_iRDIE ) {
		// Use distance dependence
		flDistScale = flDistanceInv * m_flRDIEInv;
		flEnergy *= flDistScale;
	}

	if ( pForces ) {
		vec3_t vecForceI;
		vec_t scalar;

		vec_t *pForceI = *(pForces + pPair->mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + pPair->mAtomIndices[1]);

		// Calculate Coulomb forces on atoms in reaction field terms
		// Fi = -dE/dRi = eec*qi*qj*(1/|Rij|^2 - 2*A*|Rij|)*R(n)ij
		// Fj = -Fi
		scalar = qIqJ * ( flDistanceInv * flDistanceInv - 2 * m_flCRF_A * flDistance );

		if ( m_iRDIE != 0 ) {
			// Use distance dependence
			scalar = flDistScale * scalar - flEnergy * flDistanceInv * flDistanceInv;
		}

		vec3_t testI, testJ;
		Vec3Copy( pForceI, testI );
		Vec3Copy( pForceJ, testJ );

		Vec3Scale( pvecDelta, scalar, vecForceI );
		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Sub( pForceJ, vecForceI, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );	
	}

	return flEnergy;
}

vec_t cModel :: CalcHBondEnergy( const cHBTriplet *pTriplet, const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flEnergy;
	vec3_t vecDH, vecAH;
	vec_t flLenAH, flLenDHInv, flLenAHInv, flCosTheta;
	vec_t flAngleDep, flRA = 0, flRB = 0, flCB = 0;

	// inverse square of sigma
	// sigmaSq = SQR( cos( 5*pi/6 ) - cos( pi ) ) = 0.01795
	const vec_t sigmaInvSq = vec_t( 1.0 / 0.01795 );

	assert( m_iHBModel == 126 || m_iHBModel == 128 );

	// Get atomic positions
	const vec_t *pPosD = *(pPositions + pTriplet->mAtomIndices[0]);
	const vec_t *pPosH = *(pPositions + pTriplet->mAtomIndices[1]);
	const vec_t *pPosA = *(pPositions + pTriplet->mAtomIndices[2]);

	// Calculate geometry data
	Vec3Sub( pPosD, pPosH, vecDH );
	Vec3Sub( pPosA, pPosH, vecAH );
	Vec3Nrm2( vecDH, flLenDHInv );
	flLenAH = Vec3Nrm2( vecAH, flLenAHInv );
	flCosTheta = Vec3Dot( vecDH, vecAH );
	CLAMP( flCosTheta, -1, 1 );

	// Calculate angle-dependent term
	flAngleDep = exp( -SQR( flCosTheta + 1 ) * sigmaInvSq );

	// Calculate distance-dependent terms
	if ( flLenAH > pTriplet->mpHBPair->mR ) {
		if ( m_iHBModel == 126 ) {
			flRA = SIXTH( flLenAHInv );
			flRB = flRA;
			flCB = 6;
		} else {
			flRA = SQR( flLenAHInv * flLenAHInv );
			flRB = SQR( flRA );
			flCB = 8;
		}
		// Calculate energy
		// E = (RA*A - B)*RB*angleDep
		flEnergy = (flRA*pTriplet->mpHBPair->mA - pTriplet->mpHBPair->mB)*flRB;
	} else {
		flEnergy = -pTriplet->mpHBPair->mE;
	}

	// Scale by angle-dependent term and global h-bond scale
	flEnergy *= flAngleDep * m_flHBScale;
/*
	Log().DPrintf( "HB: triplet [%s-%i]%s-[%s-%i]%s..[%s-%i]%s (%s..%s)\tdyh = %f\tang = %f\tE = %f\n",
					m_atoms.at(pTriplet->mAtomIndices[0]).mResidueTitle, m_atoms.at(pTriplet->mAtomIndices[0]).mResidueNumber, m_atoms.at(pTriplet->mAtomIndices[0]).mAtomTitle, 
					m_atoms.at(pTriplet->mAtomIndices[1]).mResidueTitle, m_atoms.at(pTriplet->mAtomIndices[1]).mResidueNumber, m_atoms.at(pTriplet->mAtomIndices[1]).mAtomTitle, 
					m_atoms.at(pTriplet->mAtomIndices[2]).mResidueTitle, m_atoms.at(pTriplet->mAtomIndices[2]).mResidueNumber, m_atoms.at(pTriplet->mAtomIndices[2]).mAtomTitle, 
					Config().LookupFFTitle( m_atoms.at(pTriplet->mAtomIndices[1]).mFFCode ), 
					Config().LookupFFTitle( m_atoms.at(pTriplet->mAtomIndices[2]).mFFCode ),
				 	flLenAH, RAD2DEG( acos( flCosTheta ) ), flEnergy );
*/

	if ( pForces ) {
		// Calculate forces on atoms (i=A, j=H, k=D)
		vec3_t vecForceI;
		vec3_t vecForceK;

		vec_t *pForceI = *(pForces + pTriplet->mAtomIndices[2]);
		vec_t *pForceJ = *(pForces + pTriplet->mAtomIndices[1]);
		vec_t *pForceK = *(pForces + pTriplet->mAtomIndices[0]);

		// Fi = -dE/dRi = (2*A*RA-B)*6*RB*R(n)ah/(|Rah|)*angleDep - (E*angleDepD/|Rah|)*[R(n)dh - cos*R(n)ah]
		// Fk = -dE/dRk = -(E*angleDepD/|Rkj|)*[R(n)ah - cos*R(n)dh]
		// Fj = -Fi-Fk

		vec_t flAngleDepD = 2 * ( flCosTheta + 1 ) * sigmaInvSq;
		vec_t scalar2 = flEnergy*flAngleDepD;

		Vec3MA( vecDH, -flCosTheta, vecAH, vecForceI );
		Vec3Scale( vecForceI, scalar2 * flLenAHInv, vecForceI );
		Vec3MA( vecAH, -flCosTheta, vecDH, vecForceK );
		Vec3Scale( vecForceK, scalar2 * flLenDHInv, vecForceK );

		if ( flLenAH > pTriplet->mpHBPair->mR ) {
			vec_t scalar1 = (12*pTriplet->mpHBPair->mA*flRA-flCB*pTriplet->mpHBPair->mB)*flRB*flLenAHInv*flAngleDep*m_flHBScale;
			Vec3MA( vecForceI, scalar1, vecAH, vecForceI );
		}
	
		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Add( pForceK, vecForceK, pForceK );
		Vec3Sub( pForceJ, vecForceI, pForceJ );
		Vec3Sub( pForceJ, vecForceK, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
		CHECK_NAN3( pForceK );
	}

	return flEnergy;
}

vec_t cModel :: CalcSolvationEnergy_GS( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces )
{
	// Gaussian solvent-exclusion model
	// Only for heavy atom pairs
	// Lazaridis, Karplus. Proteins 1999: 35, 133-152
	vec_t flEnergy;
	vec_t flDistanceInvSq, flXij2, flXji2, flAVExpI, flAVExpJ;

	cSolvParms *pSolvI = m_atoms.at( pPair->mAtomIndices[0] ).mpSolvation;
	cSolvParms *pSolvJ = m_atoms.at( pPair->mAtomIndices[1] ).mpSolvation;

	assert( pSolvI != NULL );
	assert( pSolvJ != NULL );

	flDistanceInvSq = flDistanceInv * flDistanceInv;
	flXij2 = SQR( flDistance - pSolvI->u.GS.mR ) * pSolvI->u.GS.mB;
	flXji2 = SQR( flDistance - pSolvJ->u.GS.mR ) * pSolvJ->u.GS.mB;

	// E = (1/|Rij|^2) * (Ai * exp(-Xij2) * Vj + Aj * exp(-Xji2) * Vi)
	flAVExpI = pSolvI->u.GS.mA * exp( flXij2 ) * pSolvJ->u.GS.mV;
	flAVExpJ = pSolvJ->u.GS.mA * exp( flXji2 ) * pSolvI->u.GS.mV;
	flEnergy = flDistanceInvSq * ( flAVExpI + flAVExpJ );

	if ( pForces ) {
		vec3_t vecForceI;
		vec_t scalar;

		vec_t *pForceI = *(pForces + pPair->mAtomIndices[0]);
		vec_t *pForceJ = *(pForces + pPair->mAtomIndices[1]);

		// Calculate solvent forces on atoms
		// Fi = -dE/dRi = R(n)ij*(2/|Rij|^2)*[Ai*Vj*exp(-Xij2)*((Ri-|Rij|)*Bi+1/|Rij|)+Aj*Vi*exp(-Xji2)*((Rj-|Rij|)*Bj+1/|Rij|)]
		// Fj = -Fi
		scalar = 2 * flDistanceInvSq * ( flAVExpI * ( ( pSolvI->u.GS.mR - flDistance ) * pSolvI->u.GS.mB + flDistanceInv ) + flAVExpJ * ( ( pSolvJ->u.GS.mR - flDistance ) * pSolvJ->u.GS.mB + flDistanceInv ) );

		Vec3Scale( pvecDelta, scalar, vecForceI );
		Vec3Add( pForceI, vecForceI, pForceI );
		Vec3Sub( pForceJ, vecForceI, pForceJ );
		CHECK_NAN3( pForceI );
		CHECK_NAN3( pForceJ );
	}

	return flEnergy;
}

vec_t cModel :: CalcSolvationEnergy_GB( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces )
{
	// Generalized Born model
	UNREFERENCED_PARAMETER( pPair );
	UNREFERENCED_PARAMETER( pvecDelta );
	UNREFERENCED_PARAMETER( flDistance );
	UNREFERENCED_PARAMETER( flDistanceInv );
	UNREFERENCED_PARAMETER( pForces );

	//!TODO: SolvGBorn
	return 0;
}

vec_t cModel :: CalcSolvationEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces )
{
	switch ( m_iSolvModel ) {
	case SOLV_GAUSS:
		if ( pPair->mFlags & NBPF_BOTH_HEAVY )
			return CalcSolvationEnergy_GS( pPair, pvecDelta, flDistance, flDistanceInv, pForces );
		break;
	case SOLV_GBORN:
		if ( pPair->mFlags & NBPF_BOTH_HEAVY )
			return CalcSolvationEnergy_GB( pPair, pvecDelta, flDistance, flDistanceInv, pForces );
		break;
	default:
		assume_0;
	}
	return 0;
}

void cModel :: CalcNonBondedEnergyIndividualR1( int num, int threadNum )
{
	vec3_t vecDelta;
	vec_t flDist, flDistInv;
	const cNonBondedPair *pPair = &m_nbPairsR1[num];

	// Get pointers to positions
	const vec_t *pPosI = *(m_pMTPositions + pPair->mAtomIndices[0]);
	const vec_t *pPosJ = *(m_pMTPositions + pPair->mAtomIndices[1]);

	// Calculate geometry data
	Vec3Sub( pPosI, pPosJ, vecDelta );
	flDist = Vec3Nrm2( vecDelta, flDistInv );

	// Calculate VdW interaction
	m_flMTVdWEnergy[threadNum] += CalcVdWEnergy( pPair, vecDelta, flDist, flDistInv, m_pMTForces[threadNum] );
	// Calculate electrostatic interaction
	m_flMTESEnergy[threadNum] += CalcElectrostaticEnergy( pPair, vecDelta, flDist, flDistInv, m_pMTForces[threadNum] );
	// Calculate solvation
	if ( m_iSolvModel )
		m_flMTSolvEnergy[threadNum] += CalcSolvationEnergy( pPair, vecDelta, flDist, flDistInv, m_pMTForces[threadNum] );
}

void cModel :: CalcNonBondedEnergyIndividualR2( int num, int threadNum )
{
	vec3_t vecDelta;
	vec_t flDist, flDistInv;
	const cNonBondedPair *pPair = &m_nbPairsR2[num];

	// Get pointers to positions
	const vec_t *pPosI = *(m_pMTPositions + pPair->mAtomIndices[0]);
	const vec_t *pPosJ = *(m_pMTPositions + pPair->mAtomIndices[1]);
	
	// Calculate geometry data
	Vec3Sub( pPosI, pPosJ, vecDelta );
	flDist = Vec3Nrm2( vecDelta, flDistInv );

	// Calculate electrostatic interaction
	m_flMTESEnergy[threadNum] += CalcElectrostaticEnergy( pPair, vecDelta, flDist, flDistInv, m_pMTForces[threadNum] );
}

void cModel :: CalcHBondEnergyIndividual( int num, int threadNum )
{
	const cHBTriplet *pTriplet = &m_hbTriplets[num];

	// Calculate H-bond interaction
	m_flMTHBondEnergy[threadNum] += CalcHBondEnergy( pTriplet, m_pMTPositions, m_pMTForces[threadNum] );
}

vec_t cModel :: CalcNonBondedEnergy( const vec4_t *pPositions, vec4_t *pForces, vec_t *pVdWPart, vec_t *pHBPart, vec_t *pSolvPart, vec_t *pESPart )
{
	size_t nbpR1size = m_nbPairsR1.size();
	size_t nbpR2size = m_nbPairsR2.size();

	assert( pPositions != NULL );

	m_pMTPositions = pPositions;
	m_pMTForces[0] = pForces;
	for ( int i = 1; i < m_iNumThreads; ++i )
		m_pMTForces[i] = pForces ? m_physForcesMT + m_header.mNumAtoms * (i-1) : NULL;

	memset( m_flMTVdWEnergy, 0, sizeof(vec_t)*m_iNumThreads );
	memset( m_flMTESEnergy, 0, sizeof(vec_t)*m_iNumThreads );
	memset( m_flMTHBondEnergy, 0, sizeof(vec_t)*m_iNumThreads );
	memset( m_flMTSolvEnergy, 0, sizeof(vec_t)*m_iNumThreads );

	if ( pForces && ( m_iNumThreads > 1 ) )
		memset( m_physForcesMT, 0, sizeof(vec4_t)*m_header.mNumAtoms*(m_iNumThreads-1) );

	// Loop through all NB pairs in R1
	if ( nbpR1size )
		ThreadManager().RunThreadsOn( this, nbpR1size, &cModel::CalcNonBondedEnergyIndividualR1 );

	// Loop through all NB pairs in R2
	if ( nbpR2size )
		ThreadManager().RunThreadsOn( this, nbpR2size, &cModel::CalcNonBondedEnergyIndividualR2 );

	if ( m_iHBModel ) {
		// Loop through all HB triplets
		size_t hbTsize = m_hbTriplets.size();
		if ( hbTsize ) ThreadManager().RunThreadsOn( this, hbTsize, &cModel::CalcHBondEnergyIndividual );
	}

	// Reduction
	vec4_t *partialForces = m_physForcesMT;
	for ( int i = 1; i < m_iNumThreads; ++i, partialForces += m_header.mNumAtoms ) {
		m_flMTVdWEnergy[0] += m_flMTVdWEnergy[i];
		m_flMTESEnergy[0] += m_flMTESEnergy[i];
		m_flMTHBondEnergy[0] += m_flMTHBondEnergy[i];
		m_flMTSolvEnergy[0] += m_flMTSolvEnergy[i];

		// Accumulate forces
		if ( pForces ) {
			for ( int j = 0; j < m_header.mNumAtoms; ++j ) {
				Vec3Add( pForces[j], partialForces[j], pForces[j] );
			}
		}
	}

	if ( m_iSolvModel == SOLV_GAUSS ) {
		// Esolv = Eref - E
		m_flMTSolvEnergy[0] = m_flSolvGref - m_flMTSolvEnergy[0];
	}

	// Add to partial energies, if needed
	if ( pVdWPart )	*pVdWPart	+= m_flMTVdWEnergy[0];
	if ( pHBPart )	*pHBPart	+= m_flMTHBondEnergy[0];
	if ( pESPart )	*pESPart	+= m_flMTESEnergy[0];
	if ( pSolvPart )*pSolvPart	+= m_flMTSolvEnergy[0];

	return m_flMTVdWEnergy[0] + m_flMTHBondEnergy[0] + m_flMTESEnergy[0] + m_flMTSolvEnergy[0];
}

vec_t cModel :: CalcEnergyAndForces( const vec4_t *pPositions, vec4_t *pForces )
{
	vec_t flTotalEnergy = 0;

	// increment profile counter
	m_prof.mECounter++;

	// Clear forces
	if ( pForces )
		memset( pForces, 0, sizeof(vec4_t) * m_header.mNumAtoms );

	// Build non-bonded lists
	ProfileStart();
	BuildNonBondedPairLists( pPositions, true, true );
	ProfileEnd( &m_prof.mNBPLTime );

	// Collect bond energy and forces
	ProfileStart();
	for ( BondArray::const_iterator it = m_bonds.begin(); it != m_bonds.end(); ++it )
		flTotalEnergy += CalcBondEnergy( *it, pPositions, pForces );
	ProfileEnd( &m_prof.mEBondTime );

	// Collect angle energy and forces
	ProfileStart();
	for ( AngleArray::const_iterator it = m_angles.begin(); it != m_angles.end(); ++it )
		flTotalEnergy += CalcAngleEnergy( *it, pPositions, pForces );
	ProfileEnd( &m_prof.mEAngleTime );

	// Collect improper energy and forces
	ProfileStart();
	for ( ImproperArray::const_iterator it = m_impropers.begin(); it != m_impropers.end(); ++it )
		flTotalEnergy += CalcImproperEnergy( *it, pPositions, pForces );
	ProfileEnd( &m_prof.mEImproperTime );

	// Collect torsion energy and forces
	ProfileStart();
	for ( TorsionArray::const_iterator it = m_torsions.begin(); it != m_torsions.end(); ++it )
		flTotalEnergy += CalcTorsionEnergy( *it, pPositions, pForces );
	ProfileEnd( &m_prof.mETorsionTime );

	// Collect non-bonded energy
	ProfileStart();
	flTotalEnergy += CalcNonBondedEnergy( pPositions, pForces, NULL, NULL, NULL, NULL );
	ProfileEnd( &m_prof.mENBTime );

	// Collect special energy and forces
	ProfileStart();
	flTotalEnergy += CalcSpecialEnergy( pPositions, pForces, NULL, NULL );
	ProfileEnd( &m_prof.mESpecialTime );

	return flTotalEnergy;
}

void cModel :: CalcEnergy( bool initial )
{
	vec_t	flBondEnergy = 0;
	vec_t	flAngleEnergy = 0;
	vec_t	flTorsionEnergy = 0;
	vec_t	flImproperEnergy = 0;
	vec_t	flVdWEnergy = 0;
	vec_t	flHBondEnergy = 0;
	vec_t	flSolvEnergy = 0;
	vec_t	flESEnergy = 0;
	vec_t	flNonBondedEnergy = 0;
	vec_t	flPosRestrEnergy = 0;
	vec_t	flDistRestrEnergy = 0;
	vec_t	flTotalEnergy = 0;
	char	szCommentLine[81];

	// increment profile counter
	m_prof.mECounter++;

	// Build non-bonded lists
	ProfileStart();
	BuildNonBondedPairLists( m_physPositions, true, true );
	ProfileEnd( &m_prof.mNBPLTime );

	// Collect bond energy
	ProfileStart();
	for ( BondArray::const_iterator it = m_bonds.begin(); it != m_bonds.end(); ++it )
		flBondEnergy += CalcBondEnergy( *it, m_physPositions, NULL );
	ProfileEnd( &m_prof.mEBondTime );

	// Collect angle energy
	ProfileStart();
	for ( AngleArray::const_iterator it = m_angles.begin(); it != m_angles.end(); ++it )
		flAngleEnergy += CalcAngleEnergy( *it, m_physPositions, NULL );
	ProfileEnd( &m_prof.mEAngleTime );

	// Collect improper energy
	ProfileStart();
	for ( ImproperArray::const_iterator it = m_impropers.begin(); it != m_impropers.end(); ++it )
		flImproperEnergy += CalcImproperEnergy( *it, m_physPositions, NULL );
	ProfileEnd( &m_prof.mEImproperTime );

	// Collect torsion energy
	ProfileStart();
	for ( TorsionArray::const_iterator it = m_torsions.begin(); it != m_torsions.end(); ++it )
		flTorsionEnergy += CalcTorsionEnergy( *it, m_physPositions, NULL );
	ProfileEnd( &m_prof.mETorsionTime );

	// Calculate partial charges
	if ( initial && m_iAutoQ )
		CalcPartialCharges( m_physPositions );

	// Collect non-bonded energy
	ProfileStart();
	flNonBondedEnergy += CalcNonBondedEnergy( m_physPositions, NULL, &flVdWEnergy, &flHBondEnergy, &flSolvEnergy, &flESEnergy );
	ProfileEnd( &m_prof.mENBTime );

	// Collect special energy and forces
	ProfileStart();
	flTotalEnergy += CalcSpecialEnergy( m_physPositions, NULL, &flPosRestrEnergy, &flDistRestrEnergy );
	ProfileEnd( &m_prof.mESpecialTime );

	// Total energy sum
	flTotalEnergy = flBondEnergy + flAngleEnergy + flTorsionEnergy + flImproperEnergy + flNonBondedEnergy;

	// Generate comment
	const char *hbTitle = ( m_iHBModel == 126 ) ? "engHb126" : "engHb128";
	const char *solvTitle = "eSolvGS";

	memset( m_header.mComment, 0, sizeof(m_header.mComment) );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "Generated by %s %i.%i\n\n", PROGRAM_TITLE, PROGRAM_VERSION_MAJOR, PROGRAM_VERSION_MINOR ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eVbondDef", flBondEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eVangDef", flAngleEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eImpDef", flImproperEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eTorsDef", flTorsionEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eVdW612", flVdWEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", hbTitle, flHBondEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eCoulomb", flESEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", solvTitle, flSolvEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "ePosRst", flPosRestrEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eDistRst", flDistRestrEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eGeoRst", flPosRestrEnergy +  + flDistRestrEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "eGeoDef", flBondEnergy + flAngleEnergy + flImproperEnergy + flTorsionEnergy ); strcat_s( m_header.mComment, szCommentLine );
	sprintf_s( szCommentLine, sizeof(szCommentLine), "%-9s: %13.5f\n", "ePOTENT", flTotalEnergy ); strcat_s( m_header.mComment, szCommentLine );

	Save( "molEnOpt.pdb" );

	m_header.mComment[0] = 0;
}

vec_t cModel :: CalcDeltaEnergy( vec_t delta )
{
	for ( int i = 0; i < m_header.mNumAtoms; ++i )
		Vec3MA( m_enOptPositions[i], delta, m_enOptForces[i], m_enOptTemp[i] );

	return CalcEnergyAndForces( m_enOptTemp, NULL );
}

vec_t cModel :: EnergyMinimizationFunction( vec_t delta )
{
	return g_pModel->CalcDeltaEnergy( delta );
}

void cModel :: Optimize( int nSteps )
{
	// Hard limit of iterations (can be increased)
	const int maxIterations = 10000;

	vec_t flEnergy;
	vec3_t *g, *h, vecShift;
	vec_t gradSquared, shift;
	int iteration;

	if ( !nSteps )
		return;

	// Allocate memory for the temporary positions and forces
	m_enOptPositions = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms );
	assert( m_enOptPositions != NULL );
	memset( m_enOptPositions, 0, sizeof(vec4_t) * m_header.mNumAtoms );
	m_enOptForces = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms );
	assert( m_enOptForces != NULL );
	memset( m_enOptForces, 0, sizeof(vec4_t) * m_header.mNumAtoms );
	m_enOptTemp = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms );
	assert( m_enOptTemp != NULL );
	memset( m_enOptTemp, 0, sizeof(vec4_t) * m_header.mNumAtoms );

	g = new vec3_t[m_header.mNumAtoms];
	h = new vec3_t[m_header.mNumAtoms];

	if ( m_iAutoQ )	CalcPartialCharges( m_physPositions );
	flEnergy = CalcEnergyAndForces( m_physPositions, m_physForces );

	gradSquared = 0;
	for ( int i = 0; i < m_header.mNumAtoms; ++i ) {
		gradSquared += Vec3LenSq( m_physForces[i] );
		Vec3Negate( m_physForces[i], g[i] );
		Vec3Copy( g[i], h[i] );
		Vec3Copy( h[i], m_physForces[i] );
	}

	// Print initial value
	Log().DPrintf( "Maximum number of iterations = %i (hard limit %i)\n", nSteps, maxIterations );
	Log().DPrintf( "Initial energy: E = %f\n", flEnergy );

	// Optimize iteratively
	// NB: since we use static member, the following code is NOT thread-safe!
	g_pModel = this;
	for ( iteration = 0; iteration < maxIterations; ++iteration ) {
		// Check if we have done all steps specified by config
		if ( iteration == nSteps )
			break;

		// Perform line minimization
		for ( int i = 0; i < m_header.mNumAtoms; ++i ) {
			Vec3Copy( m_physPositions[i], m_enOptPositions[i] );
			Vec3Copy( m_physForces[i], m_enOptForces[i] );
		}
		LineMinimization( m_header.mNumAtoms, m_physPositions, m_physForces, &cModel::EnergyMinimizationFunction );

		// Calculate shift
		shift = 0;
		for ( int i = 0; i < m_header.mNumAtoms; ++i ) {
			Vec3Sub( m_physPositions[i], m_enOptPositions[i], vecShift );
			shift += Vec3Len( vecShift );
		}

		// Calculate current energy
		if ( m_iAutoQ )	CalcPartialCharges( m_physPositions );
		flEnergy = CalcEnergyAndForces( m_physPositions, m_physForces );

		// Update gradient
		vec_t gg = 0;
		vec_t dgg = 0;

		gradSquared = 0;
		for ( int i = 0; i < m_header.mNumAtoms; ++i ) {
			gg += Vec3LenSq( g[i] );

			gradSquared += Vec3LenSq( m_physForces[i] );	// Fletcher-Reeves
			dgg += Vec3Dot( m_physForces[i], g[i] );		// add: Polak-Ribiere
		}

		dgg += gradSquared;

		if ( gg == 0 ) {
			Log().DPrintf( "Reached zero gradient!\n" );
			break;
		}

		vec_t gam = dgg / gg;

		for ( int i = 0; i < m_header.mNumAtoms; ++i ) {
			Vec3Negate( m_physForces[i], g[i] );
			Vec3MA( g[i], gam, h[i], h[i] );
			Vec3Copy( h[i], m_physForces[i] );
		}

		// Save snapshot
		CalcEnergy( false );
		
		Log().Printf( "Energy optimization step %3i is done...\n", iteration + 1 );
		Log().DPrintf( "Step %i: E = %f [dgg = %f, shift = %f, gradsq = %f]\n", iteration + 1, flEnergy, dgg, shift, gradSquared );
	}
	// NB: end of NOT thread-safe code

	// Free memory
	COM_AlignedFree( m_enOptPositions );
	COM_AlignedFree( m_enOptForces );
	COM_AlignedFree( m_enOptTemp );

	delete [] g;
	delete [] h;

	Log().Printf( "Energy optimization is done!\n" );
}
