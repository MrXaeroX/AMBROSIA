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

//-------------------------------
// Molecular topology functions
//-------------------------------

void cModel :: BuildTopology( void )
{
	// Vector holding information about residues (was something restored in it, or not)
	std::vector<bool> m_residueRestored;
	m_residueRestored.reserve( m_header.mNumResidues );

	// 1st pass: restore heavy atoms
	RestoreHeavyAtoms( m_residueRestored );

	// 2nd pass: search for S-S bonds
	FindSSBonds();

	// 3rd pass: restore hydrogens
	RestoreHydrogens( m_residueRestored );

	// 4th pass: build lists of bonded interactions (bonds, angles, torsions)
	InitializeBonds();

	// 5th pass: apply restraining parameters
	InitializeRestraints();

	// 6th pass: init physics
	InitAtomPhysics();
}

int cModel :: LookupAtomIndex( const char *atName, int residueSequenceNum, int chainId )
{
	int atIndex = 0;
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++atIndex ) {
		if ( it->mChainNumber == chainId && it->mResidueSequenceNumber == residueSequenceNum && !_stricmp( it->mAtomTitle, atName ) )
			return atIndex;
	}
	return -1;
}

vec_t cModel :: CalcDihedral( const vec_t *v1, const vec_t *v2, const vec_t *v3, const vec_t *v4 )
{
	vec3_t d1, d2, n1, n2;
	vec_t dot, angle;

	Vec3Sub( v3, v2, d1 );
	Vec3Sub( v4, v2, d2 );
	Vec3Cross( d2, d1, n1 );
	Vec3Nrm( n1 );
	Vec3Sub( v2, v1, d1 );
	Vec3Sub( v3, v1, d2 );
	Vec3Cross( d2, d1, n2 );
	Vec3Nrm( n2 );
	Vec3Sub( v2, v3, d1 );
	dot = Vec3Dot( n1, n2 );

	if ( dot < -1 ) dot = -1;
	if ( dot > 1 ) dot = 1;

	angle = acos( dot );

	Vec3Cross( n1, n2, d2 );
	dot = Vec3Dot( d2, d1 );

	if ( dot < 0 )
		angle = -angle;

	return angle;
}

void cModel :: RodriguesGibbs( const vec3_t *prevCoords, const vec_t r, const vec_t theta, const vec_t phi, vec_t *outPosition )
{
	// Reconstruct atom origin using Rodrigues-Gibbs Formulation.
	// Jerod Parsons, J. Bradley Holmes, J. Maurice Rojas, Jerry Tsai, Charlie E. M. Strauss (2005). 
	// "Practical conversion from torsion space to Cartesian space for in silico protein synthesis". 
	// Journal of Computational Chemistry 26 (10): 1063–1068.
	vec3_t axisU, axisV, vc, vs, cr, K, K2;
	vec_t rcos, rsin;

	Vec3Sub( prevCoords[0], prevCoords[2], vc );
	Vec3Sub( prevCoords[0], prevCoords[1], axisV );
	Vec3Cross( vc, axisV, axisU );
	Vec3Nrm( axisU );
	Vec3Nrm( axisV );

	// rotate V-axis around U-axis (normal) by bond angle complement theta
	rsin = sin( theta );
	rcos = cos( theta );
	
	Vec3Scale( axisV, rsin, vs );
	Vec3Scale( axisV, rcos, vc );
	Vec3Cross( axisU, vs, cr );
	Vec3Add( vc, cr, cr );
	Vec3MA( cr, Vec3Dot( axisV, axisU ) * ( 1 - rcos ), axisU, K );

	// rotate K around V-axis (BC) by dihedral phi
	rsin = sin( phi );
	rcos = cos( phi );
	
	Vec3Scale( K, rsin, vs );
	Vec3Scale( K, rcos, vc );
	Vec3Cross( axisV, vs, cr );
	Vec3Add( vc, cr, cr );
	Vec3MA( cr, Vec3Dot( K, axisV ) * ( 1 - rcos ), axisV, K2 );

	Vec3MA( prevCoords[0], r, K2, outPosition );
}

void cModel :: GetDummyCoords( const AtomArray &residueAtoms, vec3_t *duCoords )
{
	vec3_t backbone[3];
	vec3_t prevCoords[3];
	const cResidueAtom *backboneP[3];
	const cResidueAtom *duP[3];
	int i, j, numDU = 0;
	int backboneSize = 0;

	// Lookup info
	const cResidueLocation *pCurrentInfo = Config().LookupResidueInfo( m_pszCurrentResidueTitle, m_eCurrentResidueLocation );
	if ( !pCurrentInfo ) {
		// Not in database? That's really bad.
		Log().Fatal( "Residue %s (%s): not found in the database!\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	}

	// Iterate through all atoms
	for ( ResidueAtomArray::const_iterator rit = pCurrentInfo->mAtoms.begin(); rit != pCurrentInfo->mAtoms.end(); ++rit ) {
		// Get DU info
		if ( rit->mType == 'D' ) {
			if ( numDU < 3 ) duP[numDU++] = &(*rit);
			continue;
		}
		// Check the backbone only
		if ( rit->mType == 'B' ) {
			// Check if the backbone atom is connected to DU in a proper way:
			// 0: 3 2 1
			// 1: X 3 2
			// 2: X X 3
			for ( i = backboneSize, j = 3; i < 3; ++i, --j ) {
				if ( rit->mPrevious[i] != j )
					break;
			}
			if ( i != 3 )
				continue;
			// Valid backbone atom
			AtomArray::const_iterator it = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit->mTitle ) );
			if ( it == residueAtoms.end() )
				continue;

			backboneP[backboneSize] = &(*rit);
			Vec3Copy( it->mOriginalPosition, backbone[backboneSize] );
			if ( ++backboneSize == 3 )
				break;
		}
	}

	if ( numDU != 3 )
		Log().Fatal( "Residue %s (%s): corrupt Z-matrix (missing DU triple)\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	if ( backboneSize != 3 ) {
		Log().Fatal( "Residue %s (%s): corrupt Z-matrix (there is no triple of backbone atoms connected sequentially via DU)\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	}
	
	// Ready to restore DU coords
	RodriguesGibbs( backbone, backboneP[0]->mR, backboneP[1]->mTheta, backboneP[2]->mPhi, duCoords[2] );

	Vec3Copy( duCoords[2], prevCoords[0] );
	Vec3Copy( backbone[0], prevCoords[1] );
	Vec3Copy( backbone[1], prevCoords[2] );
	RodriguesGibbs( prevCoords, duP[2]->mR, backboneP[0]->mTheta, backboneP[1]->mPhi, duCoords[1] );

	Vec3Copy( duCoords[1], prevCoords[0] );
	Vec3Copy( duCoords[2], prevCoords[1] );
	Vec3Copy( backbone[0], prevCoords[2] );
	RodriguesGibbs( prevCoords, duP[1]->mR, duP[2]->mTheta, backboneP[0]->mPhi, duCoords[0] );

	Log().Verbose( " Chain %i, residue %s-%i (%s): DU coords:\n   DU1 (%8.4f, %8.4f, %8.4f)\n   DU2 (%8.4f, %8.4f, %8.4f)\n   DU3 (%8.4f, %8.4f, %8.4f)\n",
		m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, cConfig::LocationToString( m_eCurrentResidueLocation ),
		duCoords[0][0], duCoords[0][1], duCoords[0][2], duCoords[1][0], duCoords[1][1], duCoords[1][2], duCoords[2][0], duCoords[2][1], duCoords[2][2] );

#if 0
	// Check restored DU coords
	vec3_t test;

	Vec3Copy( duCoords[2], prevCoords[0] );
	Vec3Copy( duCoords[1], prevCoords[1] );
	Vec3Copy( duCoords[0], prevCoords[2] );
	RodriguesGibbs( prevCoords, backboneP[0]->mR, backboneP[0]->mTheta, backboneP[0]->mPhi, test );
	Log().DPrintf( "DU Test for `%4s': expected ( %8.4f, %8.4f, %8.4f ), got ( %8.4f, %8.4f, %8.4f )\n",
		backboneP[0]->mTitle, backbone[0][0], backbone[0][1], backbone[0][2], test[0], test[1], test[2] );

	Vec3Copy( backbone[0], prevCoords[0] );
	Vec3Copy( duCoords[2], prevCoords[1] );
	Vec3Copy( duCoords[1], prevCoords[2] );
	RodriguesGibbs( prevCoords, backboneP[1]->mR, backboneP[1]->mTheta, backboneP[1]->mPhi, test );
	Log().DPrintf( "DU Test for `%4s': expected ( %8.4f, %8.4f, %8.4f ), got ( %8.4f, %8.4f, %8.4f )\n",
		backboneP[1]->mTitle, backbone[1][0], backbone[1][1], backbone[1][2], test[0], test[1], test[2] );

	Vec3Copy( backbone[1], prevCoords[0] );
	Vec3Copy( backbone[0], prevCoords[1] );
	Vec3Copy( duCoords[2], prevCoords[2] );
	RodriguesGibbs( prevCoords, backboneP[2]->mR, backboneP[2]->mTheta, backboneP[2]->mPhi, test );
	Log().DPrintf( "DU Test for `%4s': expected ( %8.4f, %8.4f, %8.4f ), got ( %8.4f, %8.4f, %8.4f )\n",
		backboneP[2]->mTitle, backbone[2][0], backbone[2][1], backbone[2][2], test[0], test[1], test[2] );
#endif
}

void cModel :: GetBackboneCoords( const AtomArray &residueAtoms, vec3_t *backboneCoords )
{
	const cResidueAtom *backboneP[3];
	int backboneSize = 0;

	// Lookup info
	const cResidueLocation *pCurrentInfo = Config().LookupResidueInfo( m_pszCurrentResidueTitle, m_eCurrentResidueLocation );
	if ( !pCurrentInfo ) {
		// Not in database? That's really bad.
		Log().Fatal( "Residue %s (%s): not found in the database!\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	}

	// Get coords of last 3 backbone atoms
	// Iterate through all atoms
	for ( ResidueAtomArray::const_reverse_iterator rit = pCurrentInfo->mAtoms.rbegin(); rit != pCurrentInfo->mAtoms.rend(); ++rit ) {
		// Check the backbone only
		if ( rit->mType == 'B' ) {
			// Valid backbone atom
			AtomArray::const_iterator it = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit->mTitle ) );
			if ( it == residueAtoms.end() )
				continue;

			backboneP[2-backboneSize] = &(*rit);
			Vec3Copy( it->mOriginalPosition, backboneCoords[2-backboneSize] );
			if ( ++backboneSize == 3 )
				break;
		}
	}

	if ( backboneSize != 3 )
		Log().Fatal( "Residue %s (%s): corrupt Z-matrix (there is no triple of backbone atoms connected sequentially via DU)\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );

	Log().Verbose( " Chain %i, residue %s-%i (%s): backbone coords:\n  %4s (%8.4f, %8.4f, %8.4f)\n  %4s (%8.4f, %8.4f, %8.4f)\n  %4s (%8.4f, %8.4f, %8.4f)\n",
		m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, cConfig::LocationToString( m_eCurrentResidueLocation ),
		backboneP[0]->mTitle, backboneCoords[0][0], backboneCoords[0][1], backboneCoords[0][2], 
		backboneP[1]->mTitle, backboneCoords[1][0], backboneCoords[1][1], backboneCoords[1][2], 
		backboneP[2]->mTitle, backboneCoords[2][0], backboneCoords[2][1], backboneCoords[2][2] );
}

cAtom *cModel :: RestoreAtom( const cResidueLocation *pResidueInfo, const cResidueAtom *pResidueAtom, const vec3_t *backboneCoords, AtomArray &residueAtoms )
{
	cAtom newAtom;
	vec_t phiOffset = 0;

	assert( backboneCoords != NULL );

	memset( &newAtom, 0, sizeof(newAtom) );

	strncpy_s( newAtom.mAtomSymbol, pResidueAtom->mSymbol, sizeof(newAtom.mAtomSymbol)-1 );
	strncpy_s( newAtom.mAtomTitle, pResidueAtom->mTitle, sizeof(newAtom.mAtomTitle)-1 );
	strncpy_s( newAtom.mResidueTitle, m_pszCurrentResidueTitle, sizeof(newAtom.mResidueTitle)-1 );

	newAtom.mChainNumber = m_iCurrentChain;
	newAtom.mResidueNumber = m_iCurrentResidueNumber;
	newAtom.mResidueSequenceNumber = m_iCurrentResidueSequence;

	// Calculate new position from the Z-matrix
	vec3_t prevCoords[3];
	
	for ( int i = 0; i < 3; ++i ) {
		// Get previous residue atom
		ResidueAtomArray::const_iterator itRPrev = std::find_if( pResidueInfo->mAtoms.begin(), pResidueInfo->mAtoms.end(), ResidueAtomSearchByIndexFunctor( pResidueAtom->mPrevious[i] ) );
		if ( itRPrev == pResidueInfo->mAtoms.end() )
			Log().Fatal( "Residue %s (%s): corrupt Z-matrix for atom `%s' (%i)\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ), pResidueAtom->mTitle, pResidueAtom->mIndex );

		// Get atom corresponding the residue atom
		if ( itRPrev->mType == 'D' ) {
			// Dummy atom denotes previous residue
			Vec3Copy( backboneCoords[itRPrev->mIndex-1], prevCoords[i] );
		} else {
			AtomArray::const_iterator itPrev = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( itRPrev->mTitle ) );
			if ( itPrev == residueAtoms.end() )
				Log().Fatal( "Residue %s (%s): corrupt Z-matrix for atom `%s' (%i)\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ), pResidueAtom->mTitle, pResidueAtom->mIndex );
			Vec3Copy( itPrev->mOriginalPosition, prevCoords[i] );
		}
	}

	// If we are restoring an atom, check if the reference frame is rotated, and fix dihedral angle (phi).
	// Rotation of the reference frame is a dihedral angle R-P0-R'-P1, where Pn is a previous heavy atom (prevCoords[n]),
	// R is a reference heavy successor (HS) position (from topology), and R' is an actual HS position.
	// HS is defined as any heavy atom, which has the same parents as the current atom.
	const cAtom *pHeavySuccessor = NULL;
	const cResidueAtom *pResHeavySuccessor = NULL;

	for ( ResidueAtomArray::const_iterator rit = pResidueInfo->mAtoms.begin(); rit != pResidueInfo->mAtoms.end(); ++rit ) {
		if ( rit->mType == 'D' || rit->mType == 'H' )
			continue;
		if ( rit->mPrevious[0] != pResidueAtom->mPrevious[0] || 
			 rit->mPrevious[1] != pResidueAtom->mPrevious[1] || 
			 rit->mPrevious[2] != pResidueAtom->mPrevious[2] )
			continue;
		AtomArray::const_iterator itSuccessor = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit->mTitle ) );
		if ( itSuccessor == residueAtoms.end() )
			continue;
		pResHeavySuccessor = &(*rit);
		pHeavySuccessor = &(*itSuccessor);
		break;
	}
	if ( pHeavySuccessor && pResHeavySuccessor ) {
		vec3_t succCoords;
		RodriguesGibbs( prevCoords, pResHeavySuccessor->mR, pResHeavySuccessor->mTheta, pResHeavySuccessor->mPhi, succCoords );
		phiOffset = CalcDihedral( pHeavySuccessor->mOriginalPosition, prevCoords[0], prevCoords[1], succCoords );
		Log().Verbose( " Reference frame corrected by %.3f degrees, based on heavy successor `%s'\n", RAD2DEG( phiOffset ), pHeavySuccessor->mAtomTitle );
	}

	// Restore coords from Z-matrix
	RodriguesGibbs( prevCoords, pResidueAtom->mR, pResidueAtom->mTheta, pResidueAtom->mPhi + phiOffset, newAtom.mOriginalPosition );

	// Add the atom to the residue
	residueAtoms.push_back( newAtom );
	m_header.mNumAtoms++;

	return &residueAtoms.at( residueAtoms.size() - 1 );
}

int cModel :: RestoreResidueHeavyAtoms( AtomArray &residueAtoms, AtomArray &newAtoms, const vec3_t *backboneCoords )
{
	vec3_t dummyCoords[3];
	bool dummyCalculated = false;
	int cRestored = 0;
	
	// If backbone coords are NULL, pass dummy
	if ( !backboneCoords )
		backboneCoords = dummyCoords;

	// Lookup info
	const cResidueLocation *pCurrentInfo = Config().LookupResidueInfo( m_pszCurrentResidueTitle, m_eCurrentResidueLocation );
	if ( !pCurrentInfo ) {
		// Not in database? That's really bad.
		Log().Fatal( "Residue %s (%s): not found in the database!\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	}

	Log().DPrintf( " Chain %i, residue %s-%i (%s)\n",
		m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, cConfig::LocationToString( m_eCurrentResidueLocation ) );

	for ( ResidueAtomArray::const_iterator rit = pCurrentInfo->mAtoms.begin(); rit != pCurrentInfo->mAtoms.end(); ++rit ) {
		// Ignore dummy atoms
		if ( rit->mType == 'D' )
			continue;

		// Ignore hydrogens, if automatic restoration will take place
		// Otherwise, push existing hydrogens to array
		if ( rit->mType == 'H' && !Config().Parameters().mReadHydrogens ) {
			// Restored automatically later
			continue;
		}

		// Get FF code
		int ffCode = Config().LookupFFCode( rit->mFFTitle );
		if ( ffCode < 0 )
			Log().Fatal( "Undefined force field atom type: `%s'\n", rit->mFFTitle );
		
		// Find this atom in the residue
		AtomArray::iterator found = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit->mTitle ) );
		if ( found != residueAtoms.end() ) {
			if ( m_eCurrentResidueLocation == RL_BEG || m_eCurrentResidueLocation == RL_ISO )
				found->mFlags |= AF_RL_BEGIN;
			if ( m_eCurrentResidueLocation == RL_END || m_eCurrentResidueLocation == RL_ISO )
				found->mFlags |= AF_RL_END;
			found->mSortNumber = rit->mIndex;
			found->mFFCode = ffCode;
			found->mpResidueAtom = &(*rit);
			newAtoms.push_back( *found );
		} else {
			// If backbone atom is missing, this is a fatal error
			if ( rit->mType == 'B' ) {
				Log().Fatal( "Chain %i, residue %s-%i: missing backbone atom `%s'!\n", 
					m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, rit->mTitle);
			}
			// If hydrogen atom is missing, it's OK for now
			else if ( rit->mType == 'H' ) {
				Log().DPrintf( " Chain %i, residue %s-%i: missing hydrogen `%s'!\n", 
					m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, rit->mTitle );
			}
			// Sidechain atoms will be restored
			else {
				Log().DPrintf( " Chain %i, residue %s-%i: restoring missing heavy atom `%s'...\n", 
					m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, rit->mTitle );

				// Also restore all heavy atoms that are connected through this atom
				// To ensure this, strip title for such atoms, so they will be also treated as missing
				for ( ResidueAtomArray::const_iterator rit2 = rit + 1; rit2 != pCurrentInfo->mAtoms.end(); ++rit2 ) {
					if ( rit2->mPrevious[0] == rit->mIndex ) {
						AtomArray::iterator strip = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit2->mTitle ) );
						if ( strip != residueAtoms.end() ) {
							Log().DPrintf( " Chain %i, residue %s-%i: stripping atom `%s' connected to `%s'\n", 
								m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, rit2->mTitle, rit->mTitle );
							strip->mAtomTitle[0] = '\0';
						}
					}
				}

				if ( backboneCoords == dummyCoords && !dummyCalculated ) {
					GetDummyCoords( residueAtoms, dummyCoords );
					dummyCalculated = true;
				}

				cAtom *restored = RestoreAtom( pCurrentInfo, &(*rit), backboneCoords, residueAtoms );
				if ( restored ) {
					if ( m_eCurrentResidueLocation == RL_BEG || m_eCurrentResidueLocation == RL_ISO )
						restored->mFlags |= AF_RL_BEGIN;
					if ( m_eCurrentResidueLocation == RL_END || m_eCurrentResidueLocation == RL_ISO )
						restored->mFlags |= AF_RL_END;
					restored->mSortNumber = rit->mIndex;
					restored->mFFCode = ffCode;
					restored->mpResidueAtom = &(*rit);
					newAtoms.push_back( *restored );
					++cRestored;
				}
			}
		}
	}

	return cRestored;
}

int cModel :: RestoreResidueHydrogens( AtomArray &residueAtoms, AtomArray &newAtoms, const vec3_t *backboneCoords )
{
	vec3_t dummyCoords[3];
	bool dummyCalculated = false;
	int cRestored = 0;
	
	// If backbone coords are NULL, pass dummy
	if ( !backboneCoords )
		backboneCoords = dummyCoords;

	// Lookup info
	const cResidueLocation *pCurrentInfo = Config().LookupResidueInfo( m_pszCurrentResidueTitle, m_eCurrentResidueLocation );
	if ( !pCurrentInfo ) {
		// Not in database? That's really bad.
		Log().Fatal( "Residue %s (%s): not found in the database!\n", m_pszCurrentResidueTitle, cConfig::LocationToString( m_eCurrentResidueLocation ) );
	}

	Log().DPrintf( " Chain %i, residue %s-%i (%s)\n",
		m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, cConfig::LocationToString( m_eCurrentResidueLocation ) );

	for ( ResidueAtomArray::const_iterator rit = pCurrentInfo->mAtoms.begin(); rit != pCurrentInfo->mAtoms.end(); ++rit ) {
		// Ignore dummy atoms
		if ( rit->mType == 'D' )
			continue;
		// Find this atom in the residue
		AtomArray::iterator found = std::find_if( residueAtoms.begin(), residueAtoms.end(), AtomSearchFunctor( rit->mTitle ) );
		if ( found != residueAtoms.end() ) {
			// This atom is already restored
			newAtoms.push_back( *found );
			continue;
		}

		if ( rit->mType != 'H' )
			Log().Fatal( "Internal error: heavy atom `%s' was not restored!\n", rit->mTitle );

		// Get FF code
		int ffCode = Config().LookupFFCode( rit->mFFTitle );
		if ( ffCode < 0 )
			Log().Fatal( "Undefined force field atom type: `%s'\n", rit->mFFTitle );
		
		Log().DPrintf( " Chain %i, residue %s-%i: restoring missing hydrogen `%s'...\n", 
						m_iCurrentChain, m_pszCurrentResidueTitle, m_iCurrentResidueSequence, rit->mTitle );

		if ( backboneCoords == dummyCoords && !dummyCalculated ) {
			GetDummyCoords( residueAtoms, dummyCoords );
			dummyCalculated = true;
		}

		cAtom *restored = RestoreAtom( pCurrentInfo, &(*rit), backboneCoords, residueAtoms );
		if ( restored ) {
			if ( m_eCurrentResidueLocation == RL_BEG || m_eCurrentResidueLocation == RL_ISO )
				restored->mFlags |= AF_RL_BEGIN;
			if ( m_eCurrentResidueLocation == RL_END || m_eCurrentResidueLocation == RL_ISO )
				restored->mFlags |= AF_RL_END;
			restored->mSortNumber = rit->mIndex;
			restored->mFFCode = ffCode;
			restored->mpResidueAtom = &(*rit);
			newAtoms.push_back( *restored );
			++cRestored;
		}
	}

	return cRestored;
}

void cModel :: RestoreHeavyAtoms( std::vector<bool> &resInfo )
{
	AtomArray newAtoms;
	AtomArray residueAtoms;
	vec3_t backboneCoords[3];
	bool backboneValid = false;
	bool isolated = true;
	int cRestored = 0, cRestoredTotal = 0;

	m_iCurrentChain = -1;
	m_iCurrentResidueNumber = -1;
	m_iCurrentResidueSequence = 0;
	m_pszCurrentResidueTitle = NULL;
	m_eCurrentResidueLocation = RL_BEG;

	Log().NewLine();
	Log().TPrintf( "Restoring heavy atoms...\n" );

	newAtoms.reserve( m_atoms.size() );
	residueAtoms.reserve( 128 );

	// Replace an array of atoms read from the source structure, with a new array, with all heavy atoms in place.
	for ( AtomArray::const_iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		// Check if we started a new chain
		if ( m_iCurrentChain != it->mChainNumber ) {
			if ( m_eCurrentResidueLocation == RL_BEG )
				m_eCurrentResidueLocation = RL_ISO;
			else
				m_eCurrentResidueLocation = RL_END;
		}
		// Check if we are in the previous residue
		else if ( m_iCurrentResidueNumber == it->mResidueNumber ) {
			// Append atom to the residue
			residueAtoms.push_back( *it );
			continue;
		}

		if ( m_iCurrentChain >= 0 ) {
			// The residue has been finished
			// Check if it misses any heavy atoms
			// Also append atoms to the overall array
			cRestored = RestoreResidueHeavyAtoms( residueAtoms, newAtoms, backboneValid ? backboneCoords : NULL );
			cRestoredTotal += cRestored;
			resInfo.push_back( cRestored > 0 );

			// Get backbone coords
			if ( residueAtoms.size() >= 3 ) {
				GetBackboneCoords( residueAtoms, backboneCoords );
				backboneValid = true;
			}
			isolated = false;
			m_eCurrentResidueLocation = RL_INT;
		}
		
		// Clear the residue
		residueAtoms.clear();

		// Check if we started a new chain
		if ( m_iCurrentChain != it->mChainNumber ) {
			m_iCurrentChain = it->mChainNumber;
			m_eCurrentResidueLocation = RL_BEG;
			backboneValid = false;
			isolated = true;
		} 

		m_iCurrentResidueNumber = it->mResidueNumber;
		m_iCurrentResidueSequence = it->mResidueSequenceNumber;
		m_pszCurrentResidueTitle = it->mResidueTitle;

		// Push this atom
		residueAtoms.push_back( *it );
	}

	// The last residue has been finished
	// Check if it misses any heavy atoms
	// Also append atoms to the overall array
	m_eCurrentResidueLocation = isolated ? RL_ISO : RL_END;
	cRestored = RestoreResidueHeavyAtoms( residueAtoms, newAtoms, backboneValid ? backboneCoords : NULL );
	cRestoredTotal += cRestored;
	resInfo.push_back( cRestored > 0 );

	// Assign list of heavy atoms to the current list of atoms
	m_atoms.swap( newAtoms );
	m_header.mNumAtoms = (int)m_atoms.size();

	// Sort atoms
	std::sort( m_atoms.begin(), m_atoms.end(), AtomSortFunctor() );

	// Report count of restored atoms
	Log().DPrintf( "%i heavy atom(s) restored\n", cRestored );

	// Write PDB file
	Save( "molAddHvyAt.pdb" );	
}

void cModel :: RenameResidue( int residueNumber, const char *newName )
{
	const cResidueLocation *pCurrentInfo = NULL;

	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		if ( it->mResidueNumber == residueNumber ) {
			// get pointer to new info
			if ( !pCurrentInfo ) {
				// determine location
				eResidueLocation loc = RL_INT;
				if ( ( it->mFlags & (AF_RL_BEGIN|AF_RL_END) ) == (AF_RL_BEGIN|AF_RL_END) )
					loc = RL_ISO;
				else if ( it->mFlags & AF_RL_BEGIN )
					loc = RL_BEG;
				else if ( it->mFlags & AF_RL_END )
					loc = RL_END;
				// lookup info
				pCurrentInfo = Config().LookupResidueInfo( newName, loc );
				if ( !pCurrentInfo ) {
					// Not in database? That's really bad.
					Log().Fatal( "Residue %s (%s): not found in the database!\n", newName, cConfig::LocationToString( loc ) );
				}
			}
			// copy new name
			strcpy_s( it->mResidueTitle, newName );
			
			// change atom properties
			ResidueAtomArray::const_iterator rit = std::find_if( pCurrentInfo->mAtoms.begin(), pCurrentInfo->mAtoms.end(), ResidueAtomSearchByTitleFunctor( it->mAtomTitle ) );
			if ( rit != pCurrentInfo->mAtoms.end() ) {
				it->mpResidueAtom = &(*rit);
				it->mFFCode = Config().LookupFFCode( rit->mFFTitle );
				it->mSortNumber = rit->mIndex;
			}
		}
	}
}

void cModel :: FindSSBonds( void )
{
	std::vector<SSBondSource> ssBondSrc;
	cSSBond newSSBond;
	int numSSbonds = 0;

	Log().NewLine();
	Log().TPrintf( "Searching for S-S bonds...\n" );

	const vec_t maxSSDistSq = SQR( Config().Parameters().mSSBondDist );

	// Generate list of all atoms that can potentially form an S-S bond (SG in CYS, CYX)
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		if ( _stricmp( it->mAtomTitle, "SG" ) )
			continue;

		if (( _stricmp( it->mResidueTitle, "CYX" ) ) && 
			( _stricmp( it->mResidueTitle, "CYS" ) || !Config().Parameters().mAutoSSBonds )  ) 
			 continue;

		SSBondSource src;
		src.pAtom1 = &(*it);
		src.pAtom2 = NULL;
		ssBondSrc.push_back( src );
	}

	// Build list of S-S bonds
	newSSBond.mAtomIndices[0] = -1;
	newSSBond.mAtomIndices[1] = -1;

	for ( std::vector<SSBondSource>::iterator it = ssBondSrc.begin(); it != ssBondSrc.end(); ++it ) {
		if ( it->pAtom2 != NULL )
			continue;

		vec_t flMinDistSq = maxSSDistSq;

		std::vector<SSBondSource>::iterator bestIt = ssBondSrc.end();

		for ( std::vector<SSBondSource>::iterator it2 = it+1; it2 != ssBondSrc.end(); ++it2 ) {
			vec3_t delta;
			vec_t deltaLenSq;

			if ( it2->pAtom2 != NULL )
				continue;
			
			Vec3Sub( it->pAtom1->mOriginalPosition, it2->pAtom1->mOriginalPosition, delta );
			deltaLenSq = Vec3LenSq( delta );
			if ( deltaLenSq < flMinDistSq ) {
				flMinDistSq = deltaLenSq;
				bestIt = it2;
			}
		}

		if ( bestIt != ssBondSrc.end() ) {
			// remember pointers to atoms
			it->pAtom2 = bestIt->pAtom1;
			bestIt->pAtom2 = it->pAtom1;

			// SS-bond found!
			++numSSbonds;
			Log().DPrintf( "Found S-S bond: %s-%i `%s' - %s-%i `%s' (d = %.3f)\n",
				it->pAtom1->mResidueTitle, it->pAtom1->mResidueSequenceNumber, it->pAtom1->mAtomTitle,
				it->pAtom2->mResidueTitle, it->pAtom2->mResidueSequenceNumber, it->pAtom2->mAtomTitle,
				sqrt( flMinDistSq ) );

			// rename residues to CYX, if they are CYS
			if ( !_stricmp( it->pAtom1->mResidueTitle, "CYS" ) )
				RenameResidue( it->pAtom1->mResidueNumber, "CYX" );
			if ( !_stricmp( it->pAtom2->mResidueTitle, "CYS" ) )
				RenameResidue( it->pAtom2->mResidueNumber, "CYX" );

			// remember the bond
			newSSBond.mResidueNumbers[0] = it->pAtom1->mResidueNumber;
			newSSBond.mResidueNumbers[1] = it->pAtom2->mResidueNumber;
			strcpy_s( newSSBond.mAtomTitles[0], it->pAtom1->mAtomTitle );
			strcpy_s( newSSBond.mAtomTitles[1], it->pAtom2->mAtomTitle );
			m_SSbonds.push_back( newSSBond );
		}
	}

	Log().DPrintf( "%i S-S bond(s) found\n", numSSbonds );
}

void cModel :: RestoreHydrogens( const std::vector<bool> &resInfo )
{
	AtomArray newAtoms;
	AtomArray residueAtoms;
	int residueCounter = 0;
	vec3_t backboneCoords[3];
	bool backboneValid = false;
	bool isolated = true;
	int cRestored = 0;

	m_iCurrentChain = -1;
	m_iCurrentResidueNumber = -1;
	m_iCurrentResidueSequence = 0;
	m_pszCurrentResidueTitle = NULL;
	m_eCurrentResidueLocation = RL_BEG;

	newAtoms.reserve( m_atoms.size() );
	residueAtoms.reserve( 128 );

	// Restore all missing hydrogen atoms
	// If Hread is not set, there are no hydrogens at all
	// Otherwise, there are some hydrogens read from the source structure (but probably not all)

	Log().NewLine();
	Log().TPrintf( "Restoring hydrogen atoms...\n" );

	// Replace an array of atoms read from the source structure, with a new array, with ALL atoms in place.
	for ( AtomArray::const_iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		// Check if we started a new chain
		if ( m_iCurrentChain != it->mChainNumber ) {
			m_eCurrentResidueLocation = RL_END;
		}
		// Check if we are in the previous residue
		else if ( m_iCurrentResidueNumber == it->mResidueNumber ) {
			// Append atom to the residue
			// Don't append, if the residue was somehow modified
			if ( strcmp( it->mAtomSymbol, "H" ) || !resInfo[residueCounter] )
				residueAtoms.push_back( *it );
			continue;
		}

		if ( m_iCurrentChain >= 0 ) {
			// The residue has been finished
			// Check if it misses any heavy atoms
			// Also append atoms to the overall array
			cRestored += RestoreResidueHydrogens( residueAtoms, newAtoms, backboneValid ? backboneCoords : NULL );
			++residueCounter;

			// Get backbone coords
			if ( residueAtoms.size() >= 3 ) {
				GetBackboneCoords( residueAtoms, backboneCoords );
				backboneValid = true;
			}
			isolated = false;
			m_eCurrentResidueLocation = RL_INT;
		}
		
		// Clear the residue
		residueAtoms.clear();

		// Check if we started a new chain
		if ( m_iCurrentChain != it->mChainNumber ) {
			m_iCurrentChain = it->mChainNumber;
			m_eCurrentResidueLocation = RL_BEG;
			backboneValid = false;
			isolated = true;
		} 

		m_iCurrentResidueNumber = it->mResidueNumber;
		m_iCurrentResidueSequence = it->mResidueSequenceNumber;
		m_pszCurrentResidueTitle = it->mResidueTitle;

		// Push this atom
		// Don't push, if the residue was somehow modified
		if ( strcmp( it->mAtomSymbol, "H" ) || !resInfo[residueCounter] )
			residueAtoms.push_back( *it );
	}

	// The last residue has been finished
	// Check if it misses any hydrogen atoms
	// Also append atoms to the overall array
	m_eCurrentResidueLocation = isolated ? RL_ISO : RL_END;
	cRestored += RestoreResidueHydrogens( residueAtoms, newAtoms, backboneValid ? backboneCoords : NULL );
	++residueCounter;

	// Assign list of heavy atoms to the current list of atoms
	m_atoms.swap( newAtoms );
	m_header.mNumAtoms = (int)m_atoms.size();

	// Sort atoms, if needed
	if ( cRestored )
		std::sort( m_atoms.begin(), m_atoms.end(), AtomSortFunctor() );

	// Report count of restored atoms
	Log().DPrintf( "%i hydrogen atom(s) restored\n", cRestored );

	// Write PDB file
	Save( "molAllAtXYZ.pdb" );
}

bool cModel :: IsBondedConnection( const AtomArray::const_iterator &atom1, const AtomArray::const_iterator &atom2 )
{
	assert( atom1->mpResidueAtom != NULL );
	assert( atom2->mpResidueAtom != NULL );

	// check if atoms are on the same chain
	if ( atom1->mChainNumber != atom2->mChainNumber )
		return false;

	// check if atoms are on the same residue
	// the check is simpler in that case
	if ( atom1->mResidueNumber == atom2->mResidueNumber ) {
		// atom2 is always prior to atom1
		return ( atom1->mpResidueAtom->mPrevious[0] == atom2->mpResidueAtom->mIndex );
	} else {
		// if atom1 is at the beginning of the topology, and atom2 is in the previous residue 
		// and at the end of its topology, then they are definitely connected.
		if ( atom1->mpResidueAtom->mFlags & TF_TOPOLOGY_BEGIN ) {
			if ( ( atom2->mpResidueAtom->mFlags & TF_TOPOLOGY_END ) &&
				   atom2->mResidueNumber == atom1->mResidueNumber - 1 )
				return true;
		}
	}

	return false;
}

bool cModel :: IsBondedConnection( const cAtom &atom1, const cAtom &atom2 )
{
	assert( atom1.mpResidueAtom != NULL );
	assert( atom2.mpResidueAtom != NULL );

	// check if atoms are on the same chain
	if ( atom1.mChainNumber != atom2.mChainNumber )
		return false;

	// check if atoms are on the same residue
	// the check is simpler in that case
	if ( atom1.mResidueNumber == atom2.mResidueNumber ) {
		// atom2 is always prior to atom1
		return ( atom1.mpResidueAtom->mPrevious[0] == atom2.mpResidueAtom->mIndex );
	} else {
		// if atom1 is at the beginning of the topology, and atom2 is in the previous residue 
		// and at the end of its topology, then they are definitely connected.
		if ( atom1.mpResidueAtom->mFlags & TF_TOPOLOGY_BEGIN ) {
			if ( ( atom2.mpResidueAtom->mFlags & TF_TOPOLOGY_END ) &&
				   atom2.mResidueNumber == atom1.mResidueNumber - 1 )
				return true;
		}
	}

	return false;
}

bool cModel :: IsLoopConnection( const AtomArray::const_iterator &atom1, const AtomArray::const_iterator &atom2 )
{
	assert( atom1->mpResidueAtom != NULL );
	assert( atom2->mpResidueAtom != NULL );

	// check if atoms are on the same chain
	if ( atom1->mChainNumber != atom2->mChainNumber )
		return false;

	// check if atoms are on the same residue
	if ( atom1->mResidueNumber != atom2->mResidueNumber )
		return false;

	// check loop bits
	if ( atom1->mpResidueAtom->mLoopFlags & atom2->mpResidueAtom->mLoopFlags )
		return true;

	return false;
}

void cModel :: AddRingImproper( AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3, AtomArray::iterator &atom4 )
{
	vec_t flXi, flDot;
	vec3_t vecIJ, vecKJ, vecKL, vecMJ, vecNK;
	vec_t flLenMJ, flLenNK;

	cImproper newImproper;
	newImproper.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
	newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
	newImproper.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
	newImproper.mAtomIndices[3] = (int)(atom4 - m_atoms.begin());
	newImproper.mHarmonicConstant = 40;

	// Get pointers to positions
	// NOTE: improper center in config is atom J, but here we assume center to be atom I
	const vec_t *pPosI = atom2->mOriginalPosition;
	const vec_t *pPosJ = atom1->mOriginalPosition;
	const vec_t *pPosK = atom3->mOriginalPosition;
	const vec_t *pPosL = atom4->mOriginalPosition;

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

	if ( flXi < vec_t( M_PI * -0.5 ) )
		newImproper.mXi = -M_PI;
	else if ( flXi > vec_t( M_PI * 0.5 ) )
		newImproper.mXi = M_PI;
	else
		newImproper.mXi = 0;

	Log().DPrintf( " Ring impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s (ideal xi = %f)\n", 
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle,
					newImproper.mXi );

	m_impropers.push_back( newImproper );
}

void cModel :: FindTorsions( const AtomArray::iterator &start, const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3, bool checkSort )
{
	// Loop through all previous atoms and find neighbours (assume atoms are sorted by topology means)
	for ( AtomArray::iterator atom4 = start; atom4 != end; ++atom4 ) {
		if ( atom4 == atom1 || atom4 == atom2 || atom4 == atom3 )
			continue;
		if ( IsBondedConnection( atom3, atom4 ) || IsBondedConnection( atom4, atom3 ) ) {
			if ( checkSort && ( atom4 > atom1 ) )
				continue;

			// Check if all atoms are in an aromatic ring
			bool bRing = ( atom1->mpResidueAtom->mFlags & TF_TOPOLOGY_RING ) &&
						 ( atom2->mpResidueAtom->mFlags & TF_TOPOLOGY_RING ) &&
						 ( atom3->mpResidueAtom->mFlags & TF_TOPOLOGY_RING ) &&
						 ( atom4->mpResidueAtom->mFlags & TF_TOPOLOGY_RING );

			// Add connectivity information
			// If atoms are 1-4 in a ring, add them to 1-2-3 list to skip non-bonded pairing
			assert( atom1->mpConnectivity != NULL );
			assert( atom4->mpConnectivity != NULL );
			int atIndex1 = (int)(atom4 - m_atoms.begin());
			int atIndex2 = (int)(atom1 - m_atoms.begin());
			if ( bRing ) {
				if ( std::find( atom1->mpConnectivity->mNeighbours_12_13.begin(), atom1->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12_13.end() )
					atom1->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
				if ( std::find( atom4->mpConnectivity->mNeighbours_12_13.begin(), atom4->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom4->mpConnectivity->mNeighbours_12_13.end() )
					atom4->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );
			} else {
				if ( std::find( atom1->mpConnectivity->mNeighbours_14.begin(), atom1->mpConnectivity->mNeighbours_14.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_14.end() )
					atom1->mpConnectivity->mNeighbours_14.push_back( atIndex1 );
				if ( std::find( atom4->mpConnectivity->mNeighbours_14.begin(), atom4->mpConnectivity->mNeighbours_14.end(), atIndex2 ) == atom4->mpConnectivity->mNeighbours_14.end() )
					atom4->mpConnectivity->mNeighbours_14.push_back( atIndex2 );
			}

			// Add ring improper, if any
			if ( bRing )
				AddRingImproper( atom1, atom2, atom3, atom4 );

			// Find torsion info
			const cTorsionInfo *pInfo = Config().LookupTorsionInfo( atom1->mFFCode, atom2->mFFCode, atom3->mFFCode, atom4->mFFCode );
			if ( !pInfo ) {
				Log().Warning( "Missing force field parameters for torsion `%s'-`%s'-`%s'-`%s'!\n", 
					Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ), Config().LookupFFTitle( atom3->mFFCode ), Config().LookupFFTitle( atom4->mFFCode ) );
				continue;
			}
			
			// Check for Kharm = 0
			vec_t totalKharm = 1;
			cTorsionHarm *pHarm = pInfo->mpHarmonic;
			while ( pHarm ) {
				totalKharm *= pHarm->mHarmonicConstant;
				pHarm = pHarm->mpNext;
			}
			if ( totalKharm == 0 ) {
				Log().Verbose( "Torsion `%s'-`%s'-`%s'-`%s' ignored due to zero barrier height\n",
					Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ), Config().LookupFFTitle( atom3->mFFCode ), Config().LookupFFTitle( atom4->mFFCode ) );
				continue;
			}
			
			// New torsion found
			cTorsion newTorsion;
			newTorsion.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
			newTorsion.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
			newTorsion.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
			newTorsion.mAtomIndices[3] = (int)(atom4 - m_atoms.begin());
			newTorsion.mpHarmonic = pInfo->mpHarmonic;

			// If all atoms are parts of planar ring structure, don't divide barrier (make dihedral more ridig)
			newTorsion.mBarrierScale = !bRing;

			Log().DPrintf( " Found tors:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
				atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
				atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
				atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
				atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle );
			m_torsions.push_back( newTorsion );
		}
	}
}

void cModel :: FindImpropers( const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3 )
{
	// Loop through all previous atoms and find neighbours (assume atoms are sorted by topology means)
	for ( AtomArray::iterator atom4 = atom2; atom4 != end; ++atom4 ) {
		if ( atom4 == atom1 || atom4 == atom2 || atom4 == atom3 )
			continue;
		if ( IsBondedConnection( atom4, atom2 ) ) {
			// Improper:
			// atom1 -> atom2 <- atom3
			//			  ^
			//			  |
			//			atom4

			// Find improper info
			const cImproperInfo *pInfo = Config().LookupImproperInfo( atom1->mFFCode, atom2->mFFCode, atom3->mFFCode, atom4->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom4 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
			pInfo = Config().LookupImproperInfo( atom1->mFFCode, atom2->mFFCode, atom4->mFFCode, atom3->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom4 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom3 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle,
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
			pInfo = Config().LookupImproperInfo( atom3->mFFCode, atom2->mFFCode, atom4->mFFCode, atom1->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom3 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom4 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom1 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle,
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
			pInfo = Config().LookupImproperInfo( atom3->mFFCode, atom2->mFFCode, atom1->mFFCode, atom4->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom3 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom1 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom4 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
			pInfo = Config().LookupImproperInfo( atom4->mFFCode, atom2->mFFCode, atom1->mFFCode, atom3->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom4 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom1 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom3 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
			pInfo = Config().LookupImproperInfo( atom4->mFFCode, atom2->mFFCode, atom3->mFFCode, atom1->mFFCode );
			if ( pInfo ) {
				// New improper found
				cImproper newImproper;
				newImproper.mAtomIndices[0] = (int)(atom4 - m_atoms.begin());
				newImproper.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newImproper.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
				newImproper.mAtomIndices[3] = (int)(atom1 - m_atoms.begin());
				newImproper.mHarmonicConstant = pInfo->mHarmonicConstant;
				newImproper.mXi = pInfo->mXi;

				Log().DPrintf( " Found impr:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom4->mResidueTitle, atom4->mResidueSequenceNumber, atom4->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
					atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle,
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle );

				m_impropers.push_back( newImproper );
				continue;
			}
		}
	}
}

void cModel :: FindAngles( const AtomArray::iterator &start, AtomArray::iterator &atom1, AtomArray::iterator &atom2 )
{
	// Loop through all previous atoms and find neighbours (assume atoms are sorted by topology means)
	for ( AtomArray::iterator atom3 = start; atom3 != atom1; ++atom3 ) {
		if ( atom3 == atom2 )
			continue;

		int bondType = 0;
		if ( IsBondedConnection( atom2, atom3 ) )
			bondType = 1;
		else if ( IsBondedConnection( atom3, atom2 ) )
			bondType = 2;

		if ( bondType ) {
			// New angle found
			cAngle newAngle;
			newAngle.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
			newAngle.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
			newAngle.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
			newAngle.mHarmonicConstant = 63;			// Default for CA-CA-CA angle
			newAngle.mTheta = vec_t( DEG2RAD( 120.0 ) );// Default for CA-CA-CA angle

			// Find angle info
			const cAngleInfo *pInfo = Config().LookupAngleInfo( atom1->mFFCode, atom2->mFFCode, atom3->mFFCode );
			if ( !pInfo ) {
				Log().Warning( "Missing force field parameters for angle `%s'-`%s'-`%s', using default parameters\n", 
					Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ), Config().LookupFFTitle( atom3->mFFCode ) );
			} else {
				newAngle.mHarmonicConstant = pInfo->mHarmonicConstant;
				newAngle.mTheta = pInfo->mTheta;
			}

			Log().DPrintf( " Found angle: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
				atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
				atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
				atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle );
			m_angles.push_back( newAngle );

			// Add connectivity information
			assert( atom1->mpConnectivity != NULL );
			assert( atom3->mpConnectivity != NULL );
			int atIndex1 = newAngle.mAtomIndices[2];
			int atIndex2 = newAngle.mAtomIndices[0];
			if ( std::find( atom1->mpConnectivity->mNeighbours_12_13.begin(), atom1->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12_13.end() )
				atom1->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
			if ( std::find( atom3->mpConnectivity->mNeighbours_12_13.begin(), atom3->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom3->mpConnectivity->mNeighbours_12_13.end() )
				atom3->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

			// Collect improper torsions for the angle
			FindImpropers( atom1, atom1, atom2, atom3 );

			// Collect proper torsions for the angle
			FindTorsions( start, atom1, atom1, atom2, atom3, bondType == 2 );
		}
	}
}

void cModel :: FindLoopAngles( const AtomArray::iterator &start, const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2 )
{
	// Loop through all previous atoms and find neighbours (assume atoms are sorted by topology means)
	for ( AtomArray::iterator atom3 = start; atom3 != end; ++atom3 ) {
		if ( atom3 == atom1 || atom3 == atom2 )
			continue;

		// Check if atom1 is top of angle
		if ( IsBondedConnection( atom1, atom3 ) || IsBondedConnection( atom3, atom1 ) ) {
			// New angle found
			cAngle newAngle;
			newAngle.mAtomIndices[0] = (int)(atom2 - m_atoms.begin());
			newAngle.mAtomIndices[1] = (int)(atom1 - m_atoms.begin());
			newAngle.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
			newAngle.mHarmonicConstant = 63;			// Default for CA-CA-CA angle
			newAngle.mTheta = vec_t( DEG2RAD( 120.0 ) );// Default for CA-CA-CA angle

			// Find angle info
			const cAngleInfo *pInfo = Config().LookupAngleInfo( atom2->mFFCode, atom1->mFFCode, atom3->mFFCode );
			if ( !pInfo ) {
				Log().Warning( "Missing force field parameters for angle `%s'-`%s'-`%s', using default parameters\n", 
					Config().LookupFFTitle( atom2->mFFCode ), Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom3->mFFCode ) );
			} else {
				newAngle.mHarmonicConstant = pInfo->mHarmonicConstant;
				newAngle.mTheta = pInfo->mTheta;
			}

			Log().DPrintf( " Found loop angle: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
				atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
				atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
				atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle );
			m_angles.push_back( newAngle );

			// Add connectivity information
			assert( atom2->mpConnectivity != NULL );
			assert( atom3->mpConnectivity != NULL );
			int atIndex1 = newAngle.mAtomIndices[2];
			int atIndex2 = newAngle.mAtomIndices[0];
			if ( std::find( atom2->mpConnectivity->mNeighbours_12_13.begin(), atom2->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom2->mpConnectivity->mNeighbours_12_13.end() )
				atom2->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
			if ( std::find( atom3->mpConnectivity->mNeighbours_12_13.begin(), atom3->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom3->mpConnectivity->mNeighbours_12_13.end() )
				atom3->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

			// Collect improper torsions for the angle
			FindImpropers( end, atom2, atom1, atom3 );

			// Collect proper torsions for the angle
			FindTorsions( start, end, atom2, atom1, atom3, false );
			FindTorsions( start, end, atom3, atom1, atom2, false );
		}
		// Check if atom2 is top of angle
		else if ( IsBondedConnection( atom2, atom3 ) || IsBondedConnection( atom3, atom2 ) ) {
			// New angle found
			cAngle newAngle;
			newAngle.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
			newAngle.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
			newAngle.mAtomIndices[2] = (int)(atom3 - m_atoms.begin());
			newAngle.mHarmonicConstant = 63;			// Default for CA-CA-CA angle
			newAngle.mTheta = vec_t( DEG2RAD( 120.0 ) );// Default for CA-CA-CA angle

			// Find angle info
			const cAngleInfo *pInfo = Config().LookupAngleInfo( atom1->mFFCode, atom2->mFFCode, atom3->mFFCode );
			if ( !pInfo ) {
				Log().Warning( "Missing force field parameters for angle `%s'-`%s'-`%s', using default parameters\n", 
					Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ), Config().LookupFFTitle( atom3->mFFCode ) );
			} else {
				newAngle.mHarmonicConstant = pInfo->mHarmonicConstant;
				newAngle.mTheta = pInfo->mTheta;
			}

			Log().DPrintf( " Found loop angle: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
				atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
				atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle,
				atom3->mResidueTitle, atom3->mResidueSequenceNumber, atom3->mAtomTitle );
			m_angles.push_back( newAngle );

			// Add connectivity information
			assert( atom1->mpConnectivity != NULL );
			assert( atom3->mpConnectivity != NULL );
			int atIndex1 = newAngle.mAtomIndices[2];
			int atIndex2 = newAngle.mAtomIndices[0];
			if ( std::find( atom1->mpConnectivity->mNeighbours_12_13.begin(), atom1->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12_13.end() )
				atom1->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
			if ( std::find( atom3->mpConnectivity->mNeighbours_12_13.begin(), atom3->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom3->mpConnectivity->mNeighbours_12_13.end() )
				atom3->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

			// Collect improper torsions for the angle
			FindImpropers( atom1, atom1, atom2, atom3 );

			// Collect proper torsions for the angle
			FindTorsions( start, atom1, atom1, atom2, atom3, false );
			FindTorsions( start, atom1, atom3, atom2, atom1, false );
		}
	}
}

void cModel :: FindBonds( const AtomArray::iterator &start, AtomArray::iterator &atom1 )
{
	// Loop through all previous atoms and find neighbours (assume atoms are sorted by topology means)
	// Pass 1: search for regular bonds
	for ( AtomArray::iterator atom2 = start; atom2 != atom1; ++atom2 ) {
		if ( IsBondedConnection( atom1, atom2 ) ) {
			// New bond found
			cBond newBond;
			newBond.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
			newBond.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
			newBond.mHarmonicConstant = 469;	// Default for CA-CA bond
			newBond.mLength = vec_t( 1.400 );	// Default for CA-CA bond

			// Find bond info
			const cBondInfo *pInfo = Config().LookupBondInfo( atom1->mFFCode, atom2->mFFCode );
			if ( !pInfo ) {
				Log().Warning( "Missing force field parameters for bond `%s'-`%s', using default parameters\n", 
					Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ) );
			} else {
				newBond.mHarmonicConstant = pInfo->mHarmonicConstant;
				newBond.mLength = pInfo->mLength;
			}

			Log().DPrintf( " Found bond:  %3s-%-3i %4s :: %3s-%-3i %4s\n", 
				atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
				atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle );
			m_bonds.push_back( newBond );

			// Add connectivity information
			assert( atom1->mpConnectivity != NULL );
			assert( atom2->mpConnectivity != NULL );
			int atIndex1 = newBond.mAtomIndices[1];
			int atIndex2 = newBond.mAtomIndices[0];
			if ( std::find( atom1->mpConnectivity->mNeighbours_12_13.begin(), atom1->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12_13.end() )
				atom1->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
			if ( std::find( atom2->mpConnectivity->mNeighbours_12_13.begin(), atom2->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom2->mpConnectivity->mNeighbours_12_13.end() )
				atom2->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

			// Check for hydrogen-bonding hydrogen
			if ( atom1->mAtomSymbol[0] == 'H' && Config().IsHBHydrogen( atom1->mFFCode ) ) {
				if ( std::find( atom2->mpConnectivity->mNeighbours_12H.begin(), atom2->mpConnectivity->mNeighbours_12H.end(), atIndex2 ) == atom2->mpConnectivity->mNeighbours_12H.end() )
					atom2->mpConnectivity->mNeighbours_12H.push_back( atIndex2 );
			}
			if ( atom2->mAtomSymbol[0] == 'H' && Config().IsHBHydrogen( atom2->mFFCode ) ) {
				if ( std::find( atom1->mpConnectivity->mNeighbours_12H.begin(), atom1->mpConnectivity->mNeighbours_12H.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12H.end() )
					atom1->mpConnectivity->mNeighbours_12H.push_back( atIndex1 );
			}

			// Collect angles for this bond
			FindAngles( start, atom1, atom2 );
		}
	}

	// Only for atoms involved in loops
	if ( atom1->mpResidueAtom->mLoopFlags ) {
		// Pass 2: search for loops
		// First, get the final atom of the residue
		AtomArray::iterator end;
		for ( end = atom1 + 1; end != m_atoms.end(); ++end ) {
			if ( end->mChainNumber != atom1->mChainNumber ||
				 end->mResidueNumber != atom1->mResidueNumber ) {
				--end;
				break;
			}
		}
		if ( end == m_atoms.end() )
			--end;

		Log().Verbose( "Residue end is %s\n", end->mAtomTitle );

		for ( AtomArray::iterator atom2 = start; atom2 != atom1; ++atom2 ) {
			// Check whether the atom is involved in loop
			if ( !atom2->mpResidueAtom->mLoopFlags )
				continue;

			if ( IsLoopConnection( atom1, atom2 ) ) {
				// New bond found
				cBond newBond;
				newBond.mAtomIndices[0] = (int)(atom1 - m_atoms.begin());
				newBond.mAtomIndices[1] = (int)(atom2 - m_atoms.begin());
				newBond.mHarmonicConstant = 469;	// Default for CA-CA bond
				newBond.mLength = vec_t( 1.400 );	// Default for CA-CA bond

				// Find bond info
				const cBondInfo *pInfo = Config().LookupBondInfo( atom1->mFFCode, atom2->mFFCode );
				if ( !pInfo ) {
					Log().Warning( "Missing force field parameters for bond `%s'-`%s', using default parameters\n", 
						Config().LookupFFTitle( atom1->mFFCode ), Config().LookupFFTitle( atom2->mFFCode ) );
				} else {
					newBond.mHarmonicConstant = pInfo->mHarmonicConstant;
					newBond.mLength = pInfo->mLength;
				}

				Log().DPrintf( " Found loop:  %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom1->mResidueTitle, atom1->mResidueSequenceNumber, atom1->mAtomTitle,
					atom2->mResidueTitle, atom2->mResidueSequenceNumber, atom2->mAtomTitle );
				m_bonds.push_back( newBond );

				// Add connectivity information
				assert( atom1->mpConnectivity != NULL );
				assert( atom2->mpConnectivity != NULL );
				int atIndex1 = newBond.mAtomIndices[1];
				int atIndex2 = newBond.mAtomIndices[0];
				if ( std::find( atom1->mpConnectivity->mNeighbours_12_13.begin(), atom1->mpConnectivity->mNeighbours_12_13.end(), atIndex1 ) == atom1->mpConnectivity->mNeighbours_12_13.end() )
					atom1->mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
				if ( std::find( atom2->mpConnectivity->mNeighbours_12_13.begin(), atom2->mpConnectivity->mNeighbours_12_13.end(), atIndex2 ) == atom2->mpConnectivity->mNeighbours_12_13.end() )
					atom2->mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

				// Collect angles for this loop
				FindLoopAngles( start, end, atom1, atom2 );
			}
		}
	}
}

void cModel :: InitializeSSBonds( void )
{
	// Register covalent bonds and angles for S-S bonds
	for ( SSBondArray::const_iterator it = m_SSbonds.begin(); it != m_SSbonds.end(); ++it ) {
		const cAtom &atom1 = m_atoms.at( it->mAtomIndices[0] );
		const cAtom &atom2 = m_atoms.at( it->mAtomIndices[1] );

		cBond newBond;
		newBond.mAtomIndices[0] = it->mAtomIndices[0];
		newBond.mAtomIndices[1] = it->mAtomIndices[1];
		newBond.mHarmonicConstant = 166;	// Default for S-S bond
		newBond.mLength = vec_t( 2.038 );	// Default for S-S bond

		// Find bond info
		const cBondInfo *pInfo = Config().LookupBondInfo( atom1.mFFCode, atom2.mFFCode );
		if ( !pInfo ) {
			Log().Warning( "Missing force field parameters for bond `%s'-`%s', using default parameters\n", 
				Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ) );
		} else {
			newBond.mHarmonicConstant = pInfo->mHarmonicConstant;
			newBond.mLength = pInfo->mLength;
		}

		m_bonds.push_back( newBond );

		// Add connectivity information
		// No need to check for duplicates, because atoms belong to different residues
		assert( atom1.mpConnectivity != NULL );
		assert( atom2.mpConnectivity != NULL );
		int atIndex1 = it->mAtomIndices[1];
		int atIndex2 = it->mAtomIndices[0];
		atom1.mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
		atom2.mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

		// Add angles and torsions
		for ( IntArray::const_iterator ni = atom1.mpConnectivity->mNeighbours_12_13.begin(); 
				ni != atom1.mpConnectivity->mNeighbours_12_13.end(); ++ni ) {
			const cAtom &atom3 = m_atoms.at( *ni );
			if ( IsBondedConnection( atom1, atom3 ) ) {
				// New angle found
				cAngle newAngle;
				newAngle.mAtomIndices[0] = it->mAtomIndices[1];
				newAngle.mAtomIndices[1] = it->mAtomIndices[0];
				newAngle.mAtomIndices[2] = *ni;
				newAngle.mHarmonicConstant = 68;				// Default for S-S-CT angle
				newAngle.mTheta = vec_t( DEG2RAD( 103.70 ) );	// Default for S-S-CT angle

				// Find angle info
				const cAngleInfo *pInfo = Config().LookupAngleInfo( atom2.mFFCode, atom1.mFFCode, atom3.mFFCode );
				if ( !pInfo ) {
					Log().Warning( "Missing force field parameters for angle `%s'-`%s'-`%s', using default parameters\n", 
						Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ) );
				} else {
					newAngle.mHarmonicConstant = pInfo->mHarmonicConstant;
					newAngle.mTheta = pInfo->mTheta;
				}

				Log().DPrintf( " Found S-S angle: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom2.mResidueTitle, atom2.mResidueSequenceNumber, atom2.mAtomTitle,
					atom1.mResidueTitle, atom1.mResidueSequenceNumber, atom1.mAtomTitle,
					atom3.mResidueTitle, atom3.mResidueSequenceNumber, atom3.mAtomTitle );
				m_angles.push_back( newAngle );

				// Add connectivity information
				assert( atom3.mpConnectivity != NULL );
				int atIndex1 = newAngle.mAtomIndices[2];
				int atIndex2 = newAngle.mAtomIndices[0];
				atom2.mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
				atom3.mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

				// Add dihedrals X-X-S-S
				for ( IntArray::const_iterator ni2 = atom3.mpConnectivity->mNeighbours_12_13.begin(); 
					ni2 != atom3.mpConnectivity->mNeighbours_12_13.end(); ++ni2 ) {
					const cAtom &atom4 = m_atoms.at( *ni2 );
					if ( IsBondedConnection( atom3, atom4 ) ) {
						// Add connectivity information
						assert( atom4.mpConnectivity != NULL );
						int atIndex1 = it->mAtomIndices[1];
						int atIndex2 = *ni2;
						atom2.mpConnectivity->mNeighbours_14.push_back( atIndex2 );
						atom4.mpConnectivity->mNeighbours_14.push_back( atIndex1 );

						// Find torsion info
						const cTorsionInfo *pInfo = Config().LookupTorsionInfo( atom4.mFFCode, atom3.mFFCode, atom1.mFFCode, atom2.mFFCode );
						if ( !pInfo ) {
							Log().Warning( "Missing force field parameters for torsion `%s'-`%s'-`%s'-`%s'!\n", 
								Config().LookupFFTitle( atom4.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ) );
							continue;
						}

						// Check for Kharm = 0
						vec_t totalKharm = 1;
						cTorsionHarm *pHarm = pInfo->mpHarmonic;
						while ( pHarm ) {
							totalKharm *= pHarm->mHarmonicConstant;
							pHarm = pHarm->mpNext;
						}
						if ( totalKharm == 0 ) {
							Log().Verbose( "Torsion `%s'-`%s'-`%s'-`%s' ignored due to zero barrier height\n",
								Config().LookupFFTitle( atom4.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ) );
							continue;
						}

						// New torsion found
						cTorsion newTorsion;
						newTorsion.mAtomIndices[0] = *ni2;
						newTorsion.mAtomIndices[1] = *ni;
						newTorsion.mAtomIndices[2] = it->mAtomIndices[0];
						newTorsion.mAtomIndices[3] = it->mAtomIndices[1];
						newTorsion.mpHarmonic = pInfo->mpHarmonic;
						newTorsion.mBarrierScale = true;

						Log().DPrintf( " Found S-S tors:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
							atom4.mResidueTitle, atom4.mResidueSequenceNumber, atom4.mAtomTitle,
							atom3.mResidueTitle, atom3.mResidueSequenceNumber, atom3.mAtomTitle,
							atom1.mResidueTitle, atom1.mResidueSequenceNumber, atom1.mAtomTitle,
							atom2.mResidueTitle, atom2.mResidueSequenceNumber, atom2.mAtomTitle );
						m_torsions.push_back( newTorsion );
					}
				}

				// Add dihedrals X-S-S-X
				for ( IntArray::const_iterator ni2 = atom2.mpConnectivity->mNeighbours_12_13.begin(); 
					ni2 != atom2.mpConnectivity->mNeighbours_12_13.end(); ++ni2 ) {
					const cAtom &atom4 = m_atoms.at( *ni2 );
					if ( IsBondedConnection( atom2, atom4 ) ) {
						// Add connectivity information
						assert( atom4.mpConnectivity != NULL );
						int atIndex1 = *ni;
						int atIndex2 = *ni2;
						atom3.mpConnectivity->mNeighbours_14.push_back( atIndex2 );
						atom4.mpConnectivity->mNeighbours_14.push_back( atIndex1 );

						// Find torsion info
						const cTorsionInfo *pInfo = Config().LookupTorsionInfo( atom3.mFFCode, atom1.mFFCode, atom2.mFFCode, atom4.mFFCode );
						if ( !pInfo ) {
							Log().Warning( "Missing force field parameters for torsion `%s'-`%s'-`%s'-`%s'!\n", 
								Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom4.mFFCode ) );
							continue;
						}
			
						// Check for Kharm = 0
						vec_t totalKharm = 1;
						cTorsionHarm *pHarm = pInfo->mpHarmonic;
						while ( pHarm ) {
							totalKharm *= pHarm->mHarmonicConstant;
							pHarm = pHarm->mpNext;
						}
						if ( totalKharm == 0 ) {
							Log().Verbose( "Torsion `%s'-`%s'-`%s'-`%s' ignored due to zero barrier height\n",
								Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom4.mFFCode ) );
							continue;
						}

						// New torsion found
						cTorsion newTorsion;
						newTorsion.mAtomIndices[0] = *ni;
						newTorsion.mAtomIndices[1] = it->mAtomIndices[0];
						newTorsion.mAtomIndices[2] = it->mAtomIndices[1];
						newTorsion.mAtomIndices[3] = *ni2;
						newTorsion.mpHarmonic = pInfo->mpHarmonic;
						newTorsion.mBarrierScale = true;

						Log().DPrintf( " Found S-S tors:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
							atom3.mResidueTitle, atom3.mResidueSequenceNumber, atom3.mAtomTitle,
							atom1.mResidueTitle, atom1.mResidueSequenceNumber, atom1.mAtomTitle,
							atom2.mResidueTitle, atom2.mResidueSequenceNumber, atom2.mAtomTitle,
							atom4.mResidueTitle, atom4.mResidueSequenceNumber, atom4.mAtomTitle );
						m_torsions.push_back( newTorsion );
					}
				}
			}
		}
		for ( IntArray::const_iterator ni = atom2.mpConnectivity->mNeighbours_12_13.begin(); 
				ni != atom2.mpConnectivity->mNeighbours_12_13.end(); ++ni ) {
			const cAtom &atom3 = m_atoms.at( *ni );
			if ( IsBondedConnection( atom2, atom3 ) ) {
				// New angle found
				cAngle newAngle;
				newAngle.mAtomIndices[0] = it->mAtomIndices[0];
				newAngle.mAtomIndices[1] = it->mAtomIndices[1];
				newAngle.mAtomIndices[2] = *ni;
				newAngle.mHarmonicConstant = 68;				// Default for S-S-CT angle
				newAngle.mTheta = vec_t( DEG2RAD( 103.70 ) );	// Default for S-S-CT angle

				// Find angle info
				const cAngleInfo *pInfo = Config().LookupAngleInfo( atom1.mFFCode, atom2.mFFCode, atom3.mFFCode );
				if ( !pInfo ) {
					Log().Warning( "Missing force field parameters for angle `%s'-`%s'-`%s', using default parameters\n", 
						Config().LookupFFTitle( atom1.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ) );
				} else {
					newAngle.mHarmonicConstant = pInfo->mHarmonicConstant;
					newAngle.mTheta = pInfo->mTheta;
				}

				Log().DPrintf( " Found S-S angle: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
					atom1.mResidueTitle, atom1.mResidueSequenceNumber, atom1.mAtomTitle,
					atom2.mResidueTitle, atom2.mResidueSequenceNumber, atom2.mAtomTitle,
					atom3.mResidueTitle, atom3.mResidueSequenceNumber, atom3.mAtomTitle );
				m_angles.push_back( newAngle );

				// Add connectivity information
				assert( atom3.mpConnectivity != NULL );
				int atIndex1 = newAngle.mAtomIndices[2];
				int atIndex2 = newAngle.mAtomIndices[1];
				atom1.mpConnectivity->mNeighbours_12_13.push_back( atIndex1 );
				atom3.mpConnectivity->mNeighbours_12_13.push_back( atIndex2 );

				// Add dihedrals X-X-S-S
				for ( IntArray::const_iterator ni2 = atom3.mpConnectivity->mNeighbours_12_13.begin(); 
					ni2 != atom3.mpConnectivity->mNeighbours_12_13.end(); ++ni2 ) {
					const cAtom &atom4 = m_atoms.at( *ni2 );
					if ( IsBondedConnection( atom3, atom4 ) ) {
						// Add connectivity information
						assert( atom4.mpConnectivity != NULL );
						int atIndex1 = it->mAtomIndices[0];
						int atIndex2 = *ni2;
						atom1.mpConnectivity->mNeighbours_14.push_back( atIndex2 );
						atom4.mpConnectivity->mNeighbours_14.push_back( atIndex1 );

						// Find torsion info
						const cTorsionInfo *pInfo = Config().LookupTorsionInfo( atom4.mFFCode, atom3.mFFCode, atom2.mFFCode, atom1.mFFCode );
						if ( !pInfo ) {
							Log().Warning( "Missing force field parameters for torsion `%s'-`%s'-`%s'-`%s'!\n", 
								Config().LookupFFTitle( atom4.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ) );
							continue;
						}

						// Check for Kharm = 0
						vec_t totalKharm = 1;
						cTorsionHarm *pHarm = pInfo->mpHarmonic;
						while ( pHarm ) {
							totalKharm *= pHarm->mHarmonicConstant;
							pHarm = pHarm->mpNext;
						}
						if ( totalKharm == 0 ) {
							Log().Verbose( "Torsion `%s'-`%s'-`%s'-`%s' ignored due to zero barrier height\n",
								Config().LookupFFTitle( atom4.mFFCode ), Config().LookupFFTitle( atom3.mFFCode ), Config().LookupFFTitle( atom2.mFFCode ), Config().LookupFFTitle( atom1.mFFCode ) );
							continue;
						}

						// New torsion found
						cTorsion newTorsion;
						newTorsion.mAtomIndices[0] = *ni2;
						newTorsion.mAtomIndices[1] = *ni;
						newTorsion.mAtomIndices[2] = it->mAtomIndices[1];
						newTorsion.mAtomIndices[3] = it->mAtomIndices[0];
						newTorsion.mpHarmonic = pInfo->mpHarmonic;
						newTorsion.mBarrierScale = true;

						Log().DPrintf( " Found S-S tors:  %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s :: %3s-%-3i %4s\n", 
							atom4.mResidueTitle, atom4.mResidueSequenceNumber, atom4.mAtomTitle,
							atom3.mResidueTitle, atom3.mResidueSequenceNumber, atom3.mAtomTitle,
							atom2.mResidueTitle, atom2.mResidueSequenceNumber, atom2.mAtomTitle,
							atom1.mResidueTitle, atom1.mResidueSequenceNumber, atom1.mAtomTitle );
						m_torsions.push_back( newTorsion );
					}
				}
			}
		}
	}
}

void cModel :: InitializeBonds( void )
{
	// The topology must be complete at this point
	// Any missing atom leads to FATAL error!
	// Now walk the structure and collect all bonded interactions

	// Since we put atomic indices from m_atoms array, into other arrays,
	// don't change the m_atoms array anymore! The only you can do is to
	// append new atoms to the end of the array (e.g. explicit water).
	int atomCounter = 0;
	AtomArray::iterator itPrevResidue;
	AtomArray::iterator itCurResidue;
	std::map<const char*,int,StringCompareFunctor> residueStats;

	Log().NewLine();
	Log().TPrintf( "Initializing bonded interactions...\n" );

	m_iCurrentChain = -1;
	m_iCurrentResidueNumber = -1;

	// Allocate and setup connectivity
	cConnectivity *pConnectivity = m_connectivity = new cConnectivity[m_atoms.size()];
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++pConnectivity ) {
		// Setup connectivity pointer
		it->mpConnectivity = pConnectivity;
	}

	// Get atomic indices for S-S bonds
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++atomCounter ) {
		for ( SSBondArray::iterator it2 = m_SSbonds.begin(); it2 != m_SSbonds.end(); ++it2 ) {
			if ( ( it2->mAtomIndices[0] < 0 ) && 
				 ( it2->mResidueNumbers[0] == it->mResidueNumber ) && 
				 !_stricmp( it2->mAtomTitles[0], it->mAtomTitle ) ) {
				it2->mAtomIndices[0] = atomCounter;
				Log().Verbose( " S-S bond: atom %s-%i resolved to %i\n", 
					it2->mAtomTitles[0], it2->mResidueNumbers[0], it2->mAtomIndices[0] );
			} else if ( ( it2->mAtomIndices[1] < 0 ) && 
				 ( it2->mResidueNumbers[1] == it->mResidueNumber ) && 
				 !_stricmp( it2->mAtomTitles[1], it->mAtomTitle ) ) {
				it2->mAtomIndices[1] = atomCounter;
				Log().Verbose( " S-S bond: atom %s-%i resolved to %i\n", 
					it2->mAtomTitles[1], it2->mResidueNumbers[1], it2->mAtomIndices[1] );
			}
		}
	}

	// Search for bonds
	m_iMaxAtomsInResidue = 0;
	int curAtoms = 0;
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		// Check for changing a chain
		if ( it->mChainNumber != m_iCurrentChain ) {
			m_iCurrentResidueNumber = -1;
			m_iCurrentChain = it->mChainNumber;
		}
		// Check for changing a residue
		if ( it->mResidueNumber != m_iCurrentResidueNumber ) {
			itPrevResidue = (m_iCurrentResidueNumber == -1) ? it : itCurResidue;
			itCurResidue = it;
			m_iCurrentResidueNumber = it->mResidueNumber;
			if ( curAtoms > m_iMaxAtomsInResidue ) 
				m_iMaxAtomsInResidue = curAtoms;
			curAtoms = 0;
			// Add to stats
			std::map<const char*,int,StringCompareFunctor>::iterator stati = residueStats.find( it->mResidueTitle );
			if ( stati == residueStats.end() ) 
				residueStats[it->mResidueTitle] = 1;
			else
				stati->second++;
		}
		++curAtoms;

		// Find bonds
		FindBonds( itPrevResidue, it );
	}
	if ( curAtoms > m_iMaxAtomsInResidue ) 
		m_iMaxAtomsInResidue = curAtoms;

	// Initialize bonds, angles and torsions for S-S bonds
	if ( m_SSbonds.size() > 0 )
		InitializeSSBonds();

	Log().DPrintf( "%i bond(s) found\n", m_bonds.size() );
	Log().DPrintf( "%i angle(s) found\n", m_angles.size() );
	Log().DPrintf( "%i improper(s) found\n", m_impropers.size() );
	Log().DPrintf( "%i torsion(s) found\n", m_torsions.size() );

	// Sort connectivity for binary search
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		assert( it->mpConnectivity != NULL );
		std::sort( it->mpConnectivity->mNeighbours_12_13.begin(), it->mpConnectivity->mNeighbours_12_13.end() );
		std::sort( it->mpConnectivity->mNeighbours_14.begin(), it->mpConnectivity->mNeighbours_14.end() );
	}

	// Verbose connectivity
	if ( Log().VerboseMode() ) {
		size_t sizeofConnectivityInfo = 0;
		Log().DPrintf( "--- Connectivity information ---\n" );

		for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
			Log().DPrintf( " %4s-%-4i %4s: (123) ", it->mResidueTitle, it->mResidueSequenceNumber, it->mAtomTitle );
			if ( it->mpConnectivity->mNeighbours_12_13.size() ) {
				sizeofConnectivityInfo += it->mpConnectivity->mNeighbours_12_13.size() * sizeof(IntArray::value_type);
				for ( IntArray::iterator itc = it->mpConnectivity->mNeighbours_12_13.begin(); itc != it->mpConnectivity->mNeighbours_12_13.end(); ++itc ) {
					const cAtom &pAt = m_atoms.at( *itc );
					Log().DPrintf( "[%s-%i %s] ", pAt.mResidueTitle, pAt.mResidueSequenceNumber, pAt.mAtomTitle );
				}
			} else {
				Log().DPrintf( "---" );
			}
			Log().DPrintf( "(1-4) " );
			if ( it->mpConnectivity->mNeighbours_14.size() ) {
				sizeofConnectivityInfo += it->mpConnectivity->mNeighbours_14.size() * sizeof(IntArray::value_type);
				for ( IntArray::iterator itc = it->mpConnectivity->mNeighbours_14.begin(); itc != it->mpConnectivity->mNeighbours_14.end(); ++itc ) {
					const cAtom &pAt = m_atoms.at( *itc );
					Log().DPrintf( "[%s-%i %s] ", pAt.mResidueTitle, pAt.mResidueSequenceNumber, pAt.mAtomTitle );
				}
			} else {
				Log().DPrintf( "---" );
			}
			Log().DPrintf( "\n" );
		}
		Log().DPrintf( "--- %.1f kb connectivity size ---\n", sizeofConnectivityInfo / 1024.0f );
	}
	
	// Report residue statistics
	Log().NewLine();
	Log().TPrintf( "Residue statistics:\n" );
	int counter = 1;
	for ( std::map<const char*,int,StringCompareFunctor>::iterator stati = residueStats.begin(); stati != residueStats.end(); ++stati, ++counter ) {
		Log().DPrintf( " %3i: %-4s %5i (%5.2f%%)\n", counter, stati->first, stati->second, (float)stati->second * 100.0f / m_header.mNumResidues );
	}
	Log().DPrintf( " Total %i residues\n", m_header.mNumResidues );
	Log().DPrintf( " Max. atoms in residue: %i\n", m_iMaxAtomsInResidue );
}

void cModel :: InitializeRestraints( void )
{
	int numPosResAt = 0;
	int numDistResAt = 0;
	int counter;

	const cConfig::PosRestrVector &posR = Config().PositionalRestraints();
	const cConfig::DistRestrVector &distR = Config().DistantRestraints();

	// check if there are some...
	if ( !posR.size() && !distR.size() )
		return;

	Log().NewLine();
	Log().TPrintf( "Initializing restraints...\n" );

	// apply positional restraints
	counter = 1;
	for ( cConfig::PosRestrVector::const_iterator rit = posR.begin(); rit != posR.end(); ++rit, ++counter ) {
		for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
			if ( it->mChainNumber == rit->mChain && it->mResidueSequenceNumber >= rit->mFirstResidue && it->mResidueSequenceNumber <= rit->mLastResidue ) {
				if ( !rit->mBackboneOnly || ( it->mpResidueAtom && it->mpResidueAtom->mType == 'B' ) ) {
					it->mFlags |= AF_RESTRAINED;
					it->mPosRestraintHarmConst = rit->mHarmConst;
					Log().DPrintf( "POSR  %5i  %-4s%-3s %c%4i    %8.3f\n", counter, it->mAtomTitle, it->mResidueTitle, it->mChainNumber + 'A' - 1, it->mResidueSequenceNumber, it->mPosRestraintHarmConst );
					++numPosResAt;
				}
			}
		}
	}

	// apply distant restraints
	counter = 1;
	for ( cConfig::DistRestrVector::const_iterator rit = distR.begin(); rit != distR.end(); ++rit, ++counter ) {
		cDistRestrain r;
		bool atomsFound = true;
		for ( int i = 0; i < 2; ++i ) {
			r.mAtomIndices[i] = LookupAtomIndex( rit->mAtomTitle[i], rit->mResidue[i], rit->mChain[i] );
			if ( r.mAtomIndices[i] < 0 ) {
				Log().Warning( "Couldn't lookup atom index for `%s'-%i-%i (distant restraint #%i)\n", 
								rit->mAtomTitle[i], rit->mResidue[i], rit->mChain[i], counter );
				atomsFound = false;
				break;
			}
		}
		if ( !atomsFound )
			continue;
		if ( r.mAtomIndices[0] == r.mAtomIndices[1] ) {
			Log().Warning( "Atoms are the same (distant restraint #%i)\n" );
			continue;
		}
		r.mLengthSquared = SQR( rit->mDistance );
		r.mHarmonicConstant = rit->mHarmConst / ( 8 * r.mLengthSquared );

		Log().DPrintf( "DISTR %5i  %-4s%-3s %c%4i   %5i  %-4s%-3s %c%4i   %8.3f %8.3f\n", 
			r.mAtomIndices[0], m_atoms.at( r.mAtomIndices[0] ).mAtomTitle, m_atoms.at( r.mAtomIndices[0] ).mResidueTitle, m_atoms.at( r.mAtomIndices[0] ).mChainNumber + 'A' - 1, m_atoms.at( r.mAtomIndices[0] ).mResidueSequenceNumber, 
			r.mAtomIndices[1], m_atoms.at( r.mAtomIndices[1] ).mAtomTitle, m_atoms.at( r.mAtomIndices[1] ).mResidueTitle, m_atoms.at( r.mAtomIndices[1] ).mChainNumber + 'A' - 1, m_atoms.at( r.mAtomIndices[1] ).mResidueSequenceNumber, 
			r.mHarmonicConstant, r.mLengthSquared );
	
		++numDistResAt;
		m_distRestraints.push_back( r );
	}

	Log().DPrintf( "%i atom(s) with positional restraints\n", numPosResAt );
	Log().DPrintf( "%i bond(s) with distant restraints\n", numDistResAt );
}

bool cModel :: InitAtomSolvation( cAtom *pAtom )
{
	if (!( pAtom->mFlags & AF_HEAVY ))
		return false;

	cSolvParms *pSolvation = pAtom->mpSolvation;

	switch ( m_iSolvModel ) {
	case SOLV_GAUSS:
		{
			const cVdWInfo *pVdWInfo = Config().LookupVdWInfo( pAtom->mFFCode );
			if ( !pVdWInfo ) {
				// Set defaults
				Log().Warning( "No VdW parameters for atom FF code: `%s'\n", Config().LookupFFTitle( pAtom->mFFCode ) );
				pSolvation->u.GS.mR = 0;
			} else {
				pSolvation->u.GS.mR = pVdWInfo->mR;
			}
			int sffCode = Config().LookupSFFCode( pAtom->mpResidueAtom->mSFFTitle );
			const cSolvGSInfo *pSolvInfo = Config().LookupSolvGSInfo( sffCode );
			if ( !pSolvInfo ) {
				// Set defaults
				Log().Warning( "No SolvGS parameters for atom SFF code: `%s'\n", pAtom->mpResidueAtom->mSFFTitle );
				pSolvation->u.GS.mA = 0;
				pSolvation->u.GS.mB = 0;
				pSolvation->u.GS.mV = 0;
			} else {
				pSolvation->u.GS.mA = vec_t( pSolvInfo->mGfree / ( 2.0 * M_PI * sqrt( M_PI ) * pSolvInfo->mLambda ) );
				pSolvation->u.GS.mB = vec_t( -1.0 / ( pSolvInfo->mLambda * pSolvInfo->mLambda ) );
				pSolvation->u.GS.mV = pSolvInfo->mVolume;
				m_flSolvGref += pSolvInfo->mGref;
			}
		}
		break;
	case SOLV_GBORN:
		Log().Fatal( "SOLV_GBORN is not implemented!\n" );
		break;
	default: assume_0;
	}
	return true;
}

void cModel :: InitAtomPhysics( void )
{
	std::set<int> usedFFCodes;
	vec3_t vecMins = { 99999, 99999, 99999 };
	vec3_t vecMaxs = { -99999, -99999, -99999 };

	m_flSolvGref = 0;

	assert( m_physAtoms == NULL );
	assert( m_physPositions == NULL );
	assert( m_physForces == NULL );
	assert( m_physForcesMT == NULL );
	assert( m_charges == NULL );
	assert( m_solvation == NULL );
	assert( m_qeqparms == NULL );
	assert( m_vdWpairs == NULL );
	assert( m_HBpairs == NULL );
	assert( m_header.mNumAtoms == m_atoms.size() );

	// Allocate aligned memory
	m_physAtoms = (cPhysAtom*)COM_AlignedMalloc( sizeof(cPhysAtom) * m_header.mNumAtoms );
	assert( m_physAtoms != NULL );
	memset( m_physAtoms, 0, sizeof(cPhysAtom) * m_header.mNumAtoms );
	m_physPositions = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms );
	assert( m_physPositions != NULL );
	memset( m_physPositions, 0, sizeof(vec4_t) * m_header.mNumAtoms );
	m_physForces = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms );
	assert( m_physForces != NULL );
	memset( m_physForces, 0, sizeof(vec4_t) * m_header.mNumAtoms );

	// Allocate additional memory for multithreaded calculations
	if ( m_iNumThreads > 1 ) {
		m_physForcesMT = (vec4_t*)COM_AlignedMalloc( sizeof(vec4_t) * m_header.mNumAtoms * ( m_iNumThreads - 1 ) );
		assert( m_physForcesMT != NULL );
		memset( m_physForcesMT, 0, sizeof(vec4_t) * m_header.mNumAtoms * ( m_iNumThreads - 1 ) );
	}

	// Allocate charges (unaligned, because they are scalars and we don't care)
	m_charges = new vec_t[m_header.mNumAtoms];
	assert( m_charges != NULL );

	// Clear memory
	memset( m_physAtoms, 0, sizeof(cPhysAtom) * m_header.mNumAtoms );
	memset( m_charges, 0, sizeof(vec_t) * m_header.mNumAtoms );

	if ( m_iSolvModel ) {
		// Allocate solvation parameters (unaligned)
		m_solvation = new cSolvParms[m_header.mNumHeavyAtoms];
		assert( m_solvation != NULL );
		// Clear memory
		memset( m_solvation, 0, sizeof(cSolvParms) * m_header.mNumHeavyAtoms );
	}
	if ( m_iAutoQ ) {
		// Allocate QEq parameters (unaligned)
		m_qeqparms = new cQEqParms[m_header.mNumAtoms];
		assert( m_qeqparms != NULL );
		// Clear memory
		memset( m_qeqparms, 0, sizeof(cQEqParms) * m_header.mNumAtoms );
	}

	// Fill properties
	cPhysAtom *pCurrentPhysAtom = m_physAtoms;
	cSolvParms *pSolvation = m_solvation;
	cQEqParms *pQEq = m_qeqparms;
	vec4_t *pCurrentPosition = m_physPositions;
	vec4_t *pCurrentForce = m_physForces;
	vec_t *pCurrentCharge = m_charges;
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++pCurrentPhysAtom, ++pCurrentPosition, ++pCurrentForce, ++pCurrentCharge ) {
		// Assign pointers
		it->mCurrentPosition = &(*pCurrentPosition)[0];
		it->mCurrentForce = &(*pCurrentForce)[0];

		// Assign flags
		if ( it->mpConnectivity->mNeighbours_12H.size() > 0 )
			it->mFlags |= AF_HB_DONOR;
		if ( Config().IsHBAcceptor( it->mFFCode ) )
			it->mFlags |= AF_HB_ACCEPTOR;

		// Copy position
		Vec3Copy( it->mOriginalPosition, it->mCurrentPosition );

		// Add to bounds
		for ( int i = 0; i < 3; ++i ) {
			if ( it->mOriginalPosition[i] < vecMins[i] ) vecMins[i] = it->mOriginalPosition[i];
			if ( it->mOriginalPosition[i] > vecMaxs[i] ) vecMaxs[i] = it->mOriginalPosition[i];
		}

		// Get atomic properties
		const cAtomInfo *pAtomInfo = Config().LookupAtomInfo( it->mAtomSymbol );
		if ( !pAtomInfo )
			Log().Fatal( "Unknown atom type: `%s'\n", it->mAtomSymbol );

		// Setup mass and radius
		pCurrentPhysAtom->mMassReciprocal = pAtomInfo->mMassReciprocal;
		pCurrentPhysAtom->mRadius = pAtomInfo->mRadius;

		// Setup charge
		*pCurrentCharge = it->mpResidueAtom->mCharge;

		// Setup solvation
		if ( m_iSolvModel != SOLV_NONE ) {
			it->mpSolvation = pSolvation;
			if ( InitAtomSolvation( &(*it) ) ) {
				++pSolvation;
			} else {
				it->mpSolvation = NULL;
			}
		}

		// Setup QEq
		if ( m_iAutoQ ) {
			const cAtomInfo *pAtInfo = Config().LookupAtomInfo( it->mAtomSymbol );
			if ( pAtInfo ) {
				pQEq->mXi = pAtInfo->mXi;
				pQEq->mJ0 = pAtInfo->mJ0;
			}
			eResidueLocation loc = RL_INT;
			if ( ( it->mFlags & (AF_RL_BEGIN|AF_RL_END) ) == (AF_RL_BEGIN|AF_RL_END) )
				loc = RL_ISO;
			else if ( it->mFlags & AF_RL_BEGIN )
				loc = RL_BEG;
			else if ( it->mFlags & AF_RL_END )
				loc = RL_END;

			const cResidueLocation *pResInfo = Config().LookupResidueInfo( it->mResidueTitle, loc );
			if ( pResInfo )
				pQEq->mQTot = pResInfo->mTotalCharge;

			it->mpQEq = pQEq++;
		}

		// Remember FF code and clear VdW code (will be set later)
		usedFFCodes.insert( it->mFFCode );
		it->mVdWCode = -1;

		// clear HB code (will be set later)
		it->mHBCode = -1;

		// Remember pointer
		it->mpPhysAtom = pCurrentPhysAtom;
	}

	// Report bounds, density, and lambda
	// lambda is a length of free particle trace: 1/(sqrt(2)*pi*d^2*n)
	// is may be utilized as a max. direct electrostatic distance (with no screening)
	// it is usually to short, however
	Log().NewLine();
	Log().DPrintf( "System AABB: (%.2f, %.2f, %.2f) - (%.2f, %.2f, %.2f)\n",
					vecMins[0], vecMins[1], vecMins[2], vecMaxs[0], vecMaxs[1], vecMaxs[2] );

	vec_t volume = (vecMaxs[0] - vecMins[0])*(vecMaxs[1] - vecMins[1])*(vecMaxs[2] - vecMins[2]);
	vec_t n = m_atoms.size() / volume;
	vec_t d2 = vec_t( 1.40 * 1.40 * 4 );
	vec_t lambda = vec_t( 1.0 ) / ( sqrt( vec_t( 2.0 ) ) * M_PI * d2 * n );

	Log().DPrintf( "Average density: %.2f molecules per cubic Angstrom\n", n );
	Log().DPrintf( "Lambda for R = 1.4: %.4f Angstroms\n", lambda );

	if ( m_iSolvModel == SOLV_GAUSS )
		Log().DPrintf( "SolvGS: Total reference free energy: %.4f\n", m_flSolvGref );

	// Allocate VdW pairs
	size_t numvdWpairs = ( usedFFCodes.size() * ( usedFFCodes.size() + 1 ) ) >> 1;
	m_vdWpairOffset = 2 * (int)usedFFCodes.size() - 1;
	Log().Verbose( "VdW: %i pair table size, %i lookup offset\n", numvdWpairs, m_vdWpairOffset );
	Log().Verbose( "VdW: %i unique FF codes in set\n", usedFFCodes.size() );

	m_vdWpairs = new cVdWPair[numvdWpairs];
	assert( m_vdWpairs != NULL );
	memset( m_vdWpairs, 0, sizeof(cVdWPair) * numvdWpairs );

	// Calculate pair parameters
	int pi = 0;
	int vdw_i = 0;
	for ( std::set<int>::iterator it = usedFFCodes.begin(); it != usedFFCodes.end(); ++it, ++vdw_i ) {
		vec_t Ri, Ei;
		const cVdWInfo *pVdWInfoI = Config().LookupVdWInfo( *it );
		if ( !pVdWInfoI ) {
			Log().Warning( "No VdW parameters for atom FF code: `%s'\n", Config().LookupFFTitle( *it ) );
			Ri = 0;
			Ei = 0;
		} else {
			Ri = pVdWInfoI->mR;
			Ei = pVdWInfoI->mE;
		}

		int vdw_j = vdw_i;
		for ( std::set<int>::iterator it2 = it; it2 != usedFFCodes.end(); ++it2, ++vdw_j ) {
			Log().Verbose( "[%i] Pairing: %i - %i (%i - %i)\n", pi++, *it, *it2, vdw_i, vdw_j );
			vec_t Rj, Ej;
			const cVdWInfo *pVdWInfoJ = Config().LookupVdWInfo( *it2 );
			if ( !pVdWInfoJ ) {
				Log().Warning( "No VdW parameters for atom FF code: `%s'\n", Config().LookupFFTitle( *it2 ) );
				Rj = 0;
				Ej = 0;
			} else {
				Rj = pVdWInfoJ->mR;
				Ej = pVdWInfoJ->mE;
			}

			assert( vdw_i <= vdw_j );
			int cellIndex = ( ( vdw_i * ( m_vdWpairOffset - vdw_i ) ) >> 1 ) + vdw_j;
			cVdWPair *pPair = m_vdWpairs + cellIndex;

			assert( cellIndex >= 0 && cellIndex < (int)numvdWpairs );
			assert( pPair->mC12 == 0 );
			assert( pPair->mC6 == 0 );

			vec_t Rij = Ri + Rj;
			vec_t Eij = sqrt( Ei * Ej );

			pPair->mC6 = 2 * Eij * SIXTH( Rij );
			pPair->mC12 = Eij * TWELFTH( Rij );
		}
	}

	// Allocate HB pairs
	const cConfig::HBCodeVector &hbHydrogens = Config().GetHBHydrogens();
	const cConfig::HBCodeVector &hbAcceptors = Config().GetHBAcceptors();

	size_t numHBpairs = hbHydrogens.size() * hbAcceptors.size();
	m_HBpairOffset = (int)hbAcceptors.size();
	Log().Verbose( "HB: %i pair table size, %i lookup offset\n", numHBpairs, m_HBpairOffset );

	m_HBpairs = new cHBPair[numHBpairs];
	assert( m_HBpairs != NULL );
	memset( m_HBpairs, 0, sizeof(cHBPair) * numHBpairs );

	// Calculate pair parameters
	int h_i = 0;
	for ( cConfig::HBCodeVector::const_iterator it = hbHydrogens.begin(); it != hbHydrogens.end(); ++it, h_i += m_HBpairOffset ) {
		int a_i = h_i;
		for ( cConfig::HBCodeVector::const_iterator it2 = hbAcceptors.begin(); it2 != hbAcceptors.end(); ++it2, ++a_i ) {
			const cHBInfo *pHBInfo = Config().LookupHBInfo( *it, *it2 );
			if ( !pHBInfo ) {
				Log().Warning( "Missing hydrogen bond parameters for `%s..%s'!\n",
					Config().LookupFFTitle( *it ), Config().LookupFFTitle( *it2 ) );
				continue;
			}
			cHBPair *pPair = m_HBpairs + a_i;
			assert( pPair->mA == 0 );
			assert( pPair->mB == 0 );

			pPair->mR = pHBInfo->mR;
			pPair->mE = pHBInfo->mE;

			switch ( m_iHBModel ) {
			case 126:
				pPair->mA = pHBInfo->mE * TWELFTH( pHBInfo->mR );
				pPair->mB = 2 * pHBInfo->mE * SIXTH( pHBInfo->mR );
				break;
			case 128:
				pPair->mA = pHBInfo->mE * TWELFTH( pHBInfo->mR );
				pPair->mB = vec_t( 1.5 ) * pHBInfo->mE * EIGHTH( pHBInfo->mR );
				break;
			default: assume_0;
			}
		}
	}

	// Assign VdW and HB codes
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it ) {
		int ffcode = it->mFFCode;
		int code = 0;
		for ( std::set<int>::iterator it2 = usedFFCodes.begin(); it2 != usedFFCodes.end(); ++it2, ++code ) {
			if ( *it2 == ffcode ) {
				it->mVdWCode = code;
				break;
			}
		}
		code = 0;
		for ( cConfig::HBCodeVector::const_iterator it2 = hbHydrogens.begin(); it2 != hbHydrogens.end(); ++it2, ++code ) {
			if ( *it2 == ffcode ) {
				it->mHBCode = code;
			}
		}
		code = 0;
		for ( cConfig::HBCodeVector::const_iterator it2 = hbAcceptors.begin(); it2 != hbAcceptors.end(); ++it2, ++code ) {
			if ( *it2 == ffcode ) {
				it->mHBCode = code;
			}
		}
		assert( it->mVdWCode >= 0 );
	}
}
