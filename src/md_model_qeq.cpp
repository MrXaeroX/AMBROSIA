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
// Initialization of qEq data (will be performed once)
//--------------------------------------------------------------------
void cModel :: InitQEq( void )
{
	assert( m_jValues == NULL );
	assert( m_cValues == NULL );
	assert( m_iValues == NULL );
	assert( m_iMaxAtomsInResidue > 0 );

	m_jValues = new vec_t[m_iMaxAtomsInResidue*m_iMaxAtomsInResidue];
	assert( m_jValues != NULL );
	memset( m_jValues, 0, sizeof(vec_t)*m_iMaxAtomsInResidue*m_iMaxAtomsInResidue );
	m_cValues = new vec_t[m_iMaxAtomsInResidue*m_iMaxAtomsInResidue];
	assert( m_cValues != NULL );
	memset( m_cValues, 0, sizeof(vec_t)*m_iMaxAtomsInResidue*m_iMaxAtomsInResidue );
	m_iValues = new int[m_iMaxAtomsInResidue];
	assert( m_iValues != NULL );
	memset( m_iValues, 0, sizeof(int)*m_iMaxAtomsInResidue );
}

//--------------------------------------------------------------------
// Calculation of qEq data for a residue
//--------------------------------------------------------------------
void cModel :: CalcQEq( const vec4_t *pPositions, const int firstAtom, const int numAtoms, vec_t *pCharges )
{
	vec3_t delta;
	vec_t dist;

	AtomArray::const_iterator itf = m_atoms.begin() + firstAtom;
	AtomArray::const_iterator it = itf;

	for ( int i = 0; i < numAtoms; ++i, ++it ) {
		m_jValues[i*numAtoms+i] = it->mpQEq->mJ0;
		m_cValues[i] = vec_t( 1.0 );
	}

	for ( int i = 0; i < numAtoms; ++i ) {
		const vec_t *pPosI = *(pPositions + i);
		for ( int j = 0; j < numAtoms; ++j ) {
			const vec_t *pPosJ = *(pPositions + j);
			if ( i != j ) {
				Vec3Sub( pPosI, pPosJ, delta );
				dist = Vec3Len( delta );
				m_jValues[i*numAtoms+j] = vec_t( 1.0 ) / ((dist / vec_t( 14.4 )) + (vec_t( 1.0 )/(m_jValues[i*numAtoms+i]+m_jValues[j*numAtoms+j])));
			}
		}
		if ( i > 0 ) {
			for ( int j = 0; j < numAtoms; ++j ) {
				m_cValues[i*numAtoms+j] = m_jValues[j] - m_jValues[i*numAtoms+j];
			}
		}
	}

	// jValues are not needed anymore, use for temp values for ludcmp
	ludcmp( m_cValues, numAtoms, m_iValues, m_jValues );

	// total charge of the residue
	pCharges[0] = itf->mpQEq->mQTot;
	it = itf + 1;
	for ( int i = 1; i < numAtoms; ++i, ++it )
		pCharges[i] = it->mpQEq->mXi - itf->mpQEq->mXi;
	lubksb( m_cValues, numAtoms, m_iValues, pCharges );

	// dump results
	Log().TPrintf( "--- Charge equilibration results for %s-%i (QTot = %f) ---\n", itf->mResidueTitle, itf->mResidueSequenceNumber, itf->mpQEq->mQTot );
	it = itf;
	vec_t test = 0;
	for ( int i = 0; i < numAtoms; ++i, ++it ) {
		Log().DPrintf( "[%4s] %5s %.5f\n", it->mResidueTitle, it->mAtomTitle, pCharges[i] );
		test += pCharges[i];
	}
	Log().DPrintf( " Total = %.2f\n", test );
	Log().NewLine();
}

//--------------------------------------------------------------------
// Charge equilibration for the calculation of partial charges
//--------------------------------------------------------------------
void cModel :: CalcPartialCharges( const vec4_t *pPositions )
{
	int currentAtom = 0;
	int firstAtom = -1;
	vec_t *pCharges = m_charges;

	if ( !m_jValues )
		InitQEq();

	// Calculate partial charges per residue
	m_iCurrentChain = -1;
	m_iCurrentResidueNumber = -1;
	for ( AtomArray::iterator it = m_atoms.begin(); it != m_atoms.end(); ++it, ++currentAtom ) {
		// Check for changing a chain
		if ( it->mChainNumber != m_iCurrentChain ) {
			m_iCurrentResidueNumber = -1;
			m_iCurrentChain = it->mChainNumber;
		}
		// Check for changing a residue
		if ( it->mResidueNumber != m_iCurrentResidueNumber ) {
			if ( firstAtom >= 0 )
				CalcQEq( pPositions + firstAtom, firstAtom, currentAtom - firstAtom, pCharges + firstAtom );
			m_iCurrentResidueNumber = it->mResidueNumber;
			firstAtom = currentAtom;
		}
	}
	if ( currentAtom > firstAtom )
		CalcQEq( pPositions + firstAtom, firstAtom, currentAtom - firstAtom, pCharges + firstAtom );
}
