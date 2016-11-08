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

//--------------------------------------------------------------------
// Generation of non-bonded pair lists and hydrogen bond triplet list
//--------------------------------------------------------------------
// BuildNonBondedPairListsIndividual: must be thread-safe!
//--------------------------------------------------------------------
void cModel :: BuildNonBondedPairListsIndividual( int num, int threadNum )
{
	const vec4_t *pPosI, *pPosJ, *pPosH;
	const cAtom *pAtI, *pAtJ, *pAtH;
	vec3_t vecDelta, vecDeltaH;
	vec_t vecDistSquared, vecDistHSquared;
	cNonBondedPair nbPair;
	cHBTriplet hbTriplet;

	if ( m_iNumThreads > 1 ) {
		if ( m_NBPR1 ) m_nbTempPairsR1[threadNum].clear();
		if ( m_NBPR2 ) m_nbTempPairsR2[threadNum].clear();
		if ( m_HBTUp ) m_hbTempTriplets[threadNum].clear();
	}

	pPosI = m_pMTPositions + num;
	pAtI = &m_atoms.at( num );
	pPosJ = pPosI + 1;
	pAtJ = pAtI + 1;
	nbPair.mAtomIndices[0] = num;

//	Log().DPrintf( "BuildNonBondedPairListsIndividual: %i %i\n", num, threadNum );

	for ( int j = num + 1; j < m_header.mNumAtoms; ++j, ++pPosJ, ++pAtJ ) {
		// Calculate squared distance
		Vec3Sub( *pPosI, *pPosJ, vecDelta );
		vecDistSquared = Vec3LenSq( vecDelta );

		// Check if too far
		if ( vecDistSquared >= m_flCutoffR2Squared )
			continue;

		// Check if atoms have bonded 1-2 or 1-3 connection
		assert( pAtI->mpConnectivity != NULL );
		if ( std::binary_search( pAtI->mpConnectivity->mNeighbours_12_13.begin(), pAtI->mpConnectivity->mNeighbours_12_13.end(), j ) )
			continue;

		nbPair.mAtomIndices[1] = j;
		nbPair.mFlags = 0;
		nbPair.mpVdWPair = NULL;

		// Check if atoms have bonded 1-4 connection
		if ( std::binary_search( pAtI->mpConnectivity->mNeighbours_14.begin(), pAtI->mpConnectivity->mNeighbours_14.end(), j ) )
			nbPair.mFlags |= NBPF_NEIGHBOURS_14;

		if ( ( pAtI->mFlags & AF_HEAVY ) && ( pAtJ->mFlags & AF_HEAVY ) )
			nbPair.mFlags |= NBPF_BOTH_HEAVY;
			
		// Add to NBPR1 list
		if ( vecDistSquared < m_flCutoffR1Squared ) {
			if ( m_NBPR1 ) {
				// Get pointer to VdW pair
				int cellIndex;
				if ( pAtI->mVdWCode <= pAtJ->mVdWCode )
					cellIndex = ( ( pAtI->mVdWCode * ( m_vdWpairOffset - pAtI->mVdWCode ) ) >> 1 ) + pAtJ->mVdWCode;
				else
					cellIndex = ( ( pAtJ->mVdWCode * ( m_vdWpairOffset - pAtJ->mVdWCode ) ) >> 1 ) + pAtI->mVdWCode;
				nbPair.mpVdWPair = m_vdWpairs + cellIndex;
				if ( m_iNumThreads > 1 )
					m_nbTempPairsR1[threadNum].push_back( nbPair );
				else
					m_nbPairsR1.push_back( nbPair );
			}
			// Add to HBT list
			if ( m_HBTUp ) {
				// One atom must be a donor, and another - an acceptor
				// This may form a valid HB triplet, loop through all hydrogens connected to the donor
				if ( ( pAtI->mFlags & AF_HB_DONOR ) && ( pAtJ->mFlags & AF_HB_ACCEPTOR ) ) {
					// Atom I is donor, build a triplet I-H..J
					hbTriplet.mAtomIndices[0] = num;
					hbTriplet.mAtomIndices[2] = j;

					// Loop through all donor's hydrogens
					for ( IntArray::const_iterator ih = pAtI->mpConnectivity->mNeighbours_12H.begin(); 
						  ih != pAtI->mpConnectivity->mNeighbours_12H.end(); ++ih ) {
						pPosH = m_pMTPositions + *ih;
						// Calculate squared distance H..A
						Vec3Sub( *pPosH, *pPosJ, vecDeltaH );
						vecDistHSquared = Vec3LenSq( vecDeltaH );
						// Check cut-off
						if ( vecDistHSquared < m_flCutoffHBSquared ) {
							pAtH = &m_atoms.at( *ih );
							hbTriplet.mAtomIndices[1] = *ih;
							assert( pAtH->mHBCode >= 0 );
							assert( pAtJ->mHBCode >= 0 );
							hbTriplet.mpHBPair = m_HBpairs + pAtH->mHBCode * m_HBpairOffset + pAtJ->mHBCode;
							if ( m_iNumThreads > 1 )
								m_hbTempTriplets[threadNum].push_back( hbTriplet );
							else
								m_hbTriplets.push_back( hbTriplet );
						/*	Log().DPrintf( "HB1: found triplet (%s%i)%s-(%s%i)%s..(%s%i)%s (%s..%s) [i = %i, j = %i]\n",
								pAtI->mResidueTitle, pAtI->mResidueNumber, pAtI->mAtomTitle,
								pAtH->mResidueTitle, pAtH->mResidueNumber, pAtH->mAtomTitle,
								pAtJ->mResidueTitle, pAtJ->mResidueNumber, pAtJ->mAtomTitle,
								Config().LookupFFTitle( pAtH->mFFCode ), Config().LookupFFTitle( pAtJ->mFFCode ), num, j );*/
						}
					}
				}
				if ( ( pAtI->mFlags & AF_HB_ACCEPTOR ) && ( pAtJ->mFlags & AF_HB_DONOR ) ) {
					// Atom J is donor, build a triplet J-H..I
					hbTriplet.mAtomIndices[0] = j;
					hbTriplet.mAtomIndices[2] = num;

					// Loop through all donor's hydrogens
					for ( IntArray::const_iterator ih = pAtJ->mpConnectivity->mNeighbours_12H.begin(); 
						  ih != pAtJ->mpConnectivity->mNeighbours_12H.end(); ++ih ) {
						pPosH = m_pMTPositions + *ih;
						// Calculate squared distance H..A
						Vec3Sub( *pPosH, *pPosI, vecDeltaH );
						vecDistHSquared = Vec3LenSq( vecDeltaH );
						// Check cut-off
						if ( vecDistHSquared < m_flCutoffHBSquared ) {
							pAtH = &m_atoms.at( *ih );
							hbTriplet.mAtomIndices[1] = *ih;
							assert( pAtH->mHBCode >= 0 );
							assert( pAtI->mHBCode >= 0 );
							hbTriplet.mpHBPair = m_HBpairs + pAtH->mHBCode * m_HBpairOffset + pAtI->mHBCode;
							if ( m_iNumThreads > 1 )
								m_hbTempTriplets[threadNum].push_back( hbTriplet );
							else
								m_hbTriplets.push_back( hbTriplet );
							/*Log().DPrintf( "HB2: found triplet (%s%i)%s-(%s%i)%s..(%s%i)%s (%s..%s) [i = %i, j = %i]\n",
								pAtJ->mResidueTitle, pAtJ->mResidueNumber, pAtJ->mAtomTitle,
								pAtH->mResidueTitle, pAtH->mResidueNumber, pAtH->mAtomTitle,
								pAtI->mResidueTitle, pAtI->mResidueNumber, pAtI->mAtomTitle,
								Config().LookupFFTitle( pAtH->mFFCode ), Config().LookupFFTitle( pAtI->mFFCode ), num, j ); */
						}
					}
				}
			}
		}
		// Add to NBPR2 list
		else if ( m_NBPR2 ) {
			if ( m_iNumThreads > 1 )
				m_nbTempPairsR2[threadNum].push_back( nbPair );
			else
				m_nbPairsR2.push_back( nbPair );
		}
	}

	if ( m_iNumThreads > 1 ) {
		// merge lists
		ThreadManager().EnterCriticalSection();
			if ( m_NBPR1 ) m_nbPairsR1.insert( m_nbPairsR1.end(), m_nbTempPairsR1[threadNum].begin(), m_nbTempPairsR1[threadNum].end() );
			if ( m_NBPR2 ) m_nbPairsR2.insert( m_nbPairsR2.end(), m_nbTempPairsR2[threadNum].begin(), m_nbTempPairsR2[threadNum].end() );
			if ( m_HBTUp ) m_hbTriplets.insert( m_hbTriplets.end(), m_hbTempTriplets[threadNum].begin(), m_hbTempTriplets[threadNum].end() );
		ThreadManager().LeaveCriticalSection();
	}
}

void cModel :: BuildNonBondedPairLists( const vec4_t *pPositions, const bool nbpR1, const bool nbpR2 )
{
	// Setup local data
	m_pMTPositions = pPositions;
	m_NBPR1 = nbpR1;
	m_NBPR2 = nbpR2;

	if ( !m_NBPInit ) {
		m_NBPInit = true;
		for ( int i = 0; i < m_iNumThreads; ++i ) {
			m_nbTempPairsR1[i].reserve( m_header.mNumAtoms );
			m_nbTempPairsR2[i].reserve( m_header.mNumAtoms );
			m_hbTempTriplets[i].reserve( m_header.mNumAtoms << 1 );
		}
	}

	assert( m_NBPR1 || m_NBPR2 );

	// Collect HB triplets, if updating R1 interactions, and HB model is not null
	m_HBTUp = m_NBPR1 && ( m_iHBModel != 0 );

	// Clear current lists
	if ( m_NBPR1 ) m_nbPairsR1.clear();
	if ( m_NBPR2 ) m_nbPairsR2.clear();
	if ( m_HBTUp ) m_hbTriplets.clear();

	// Collect non-bonded pairs (NBP) and hydrogen bonding triplets (HBT)
	ThreadManager().RunThreadsOn( this, m_header.mNumAtoms - 1, &cModel::BuildNonBondedPairListsIndividual );
	
	// Report stats
	Log().Verbose( "BuildNonBondedPairLists:\n" );
	if ( m_NBPR1 ) Log().Verbose( " R1 = %f: %i pairs (%i kb)\n", m_flCutoffR1, m_nbPairsR1.size(), ( m_nbPairsR1.size() * sizeof( cNonBondedPair ) >> 10 ) );
	if ( m_NBPR2 ) Log().Verbose( " R2 = %f: %i pairs (%i kb)\n", m_flCutoffR2, m_nbPairsR2.size(), ( m_nbPairsR2.size() * sizeof( cNonBondedPair ) >> 10 ) );
	if ( m_HBTUp ) Log().Verbose( " HB = %f: %i triplets (%i kb)\n", m_flCutoffHB, m_hbTriplets.size(), ( m_hbTriplets.size() * sizeof( cHBTriplet ) >> 10 ) );
}
