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
#include "md_pdb_format.h"

char cModel :: m_fileName[MAX_OSPATH] = DEFAULT_MODEL_FILENAME;
char cModel :: m_modelName[MAX_OSPATH] = "";
cModel :: IOMap cModel :: m_ioMap;
bool cModel :: m_bProfiling = false;

cModel :: cModel()
{
	memset( &m_prof, 0, sizeof(m_prof) );
	memset( &m_profTime, 0, sizeof(m_profTime) );
	m_profDepth = 0;
	m_iNumThreads = 0;

	memset( &m_header, 0, sizeof(m_header) );
	m_atoms.reserve( 8192 );
	m_bonds.reserve( 8192 );
	m_SSbonds.reserve( 16 );
	m_distRestraints.reserve( 16 );
	m_angles.reserve( 16384 );
	m_torsions.reserve( 16384 );
	m_impropers.reserve( 4096 );
	m_nbPairsR1.reserve( 2097152 );
	m_nbPairsR2.reserve( 2097152 );
	m_hbTriplets.reserve( 8192 );

	m_physAtoms = NULL;
	m_physPositions = NULL;
	m_physForces = NULL;
	m_physForcesMT = NULL;
	m_charges = NULL;
	m_vdWpairs = NULL;
	m_HBpairs = NULL;
	m_connectivity = NULL;
	m_solvation = NULL;
	m_qeqparms = NULL;
	m_jValues = NULL;
	m_cValues = NULL;
	m_iValues = NULL;
	m_NBPInit = false;
}

cModel :: ~cModel()
{
	Clear();
}

void cModel :: SetFileName( const char *name )
{
	strncpy_s( m_fileName, name, sizeof(m_fileName)-1 );
}

void cModel :: SetModelName( const char *name )
{
	strncpy_s( m_modelName, name, sizeof(m_modelName)-1 );
}

void cModel :: InitializeIO( void )
{
	// Register all your formats here
	m_ioMap[".pdb"] = new cPDBFormat;
}

void cModel :: InitializeParams( void )
{
	// Setup number of threads
	m_iNumThreads = ThreadManager().CountThreads();
	assert( m_iNumThreads > 0 );

	// Setup cut-off radii
	m_flCutoffR1 = Config().Parameters().mShortCutoffDist;
	m_flCutoffR1Squared = m_flCutoffR1 * m_flCutoffR1;
	m_flSwitchR1 = Config().Parameters().mShortSwitchDist;
	m_flSwitchR1Squared = m_flSwitchR1 * m_flSwitchR1;
	m_flCutoffR2 = Config().Parameters().mLongCutoffDist;
	m_flCutoffR2Squared = m_flCutoffR2 * m_flCutoffR2;
	m_flCutoffHB = Config().Parameters().mHBondCutoffDist;
	m_flCutoffHBSquared = m_flCutoffHB * m_flCutoffHB;

	// Setup reaction field constants A and B
	// Calculate Crf
	const vec_t e = Config().Parameters().mSolventDielectricConstant;
	const vec_t k = vec_t( 1.0 ) / Config().Parameters().mDebyeRadius;
	vec_t Crf = ( ( 2 - 2 * e ) * ( 1 + k * m_flCutoffR2 ) - e * ( k * m_flCutoffR2 ) * ( k * m_flCutoffR2 ) ) /
				( ( 1 + 2 * e ) * ( 1 + k * m_flCutoffR2 ) + e * ( k * m_flCutoffR2 ) * ( k * m_flCutoffR2 ) );

	// A = -0.5*Crf / Rc^3
	m_flCRF_A = vec_t( -0.5 ) * Crf / ( m_flCutoffR2Squared * m_flCutoffR2 );
	// B = -(1-0.5*Crf)/Rc
	m_flCRF_B = ( vec_t( 0.5 ) * Crf - 1 ) / m_flCutoffR2;

	// Setup AutoQ flag
	m_iAutoQ = Config().Parameters().mAutoQ;

	// Setup RDIE
	m_iRDIE = 0;
	if ( Config().Parameters().mDistanceDependentDielectricConstant != 0 ) {
		m_iRDIE = 1;
		m_flRDIEInv = vec_t( 1.0 / Config().Parameters().mDistanceDependentDielectricConstant );
	}

	// Setup H-bond model
	m_iHBModel = Config().Parameters().mHBondModel;
	m_flHBScale = Config().Parameters().mHBondScale;

	// Setup solvation model
	if ( Config().Parameters().mSolvationModel == NULL ) {
		m_iSolvModel = SOLV_NONE;
	} else if ( !_stricmp( Config().Parameters().mSolvationModel, "None" ) ) {
		m_iSolvModel = SOLV_NONE;
	} else if ( !_stricmp( Config().Parameters().mSolvationModel, "GS" ) ) {
		m_iSolvModel = SOLV_GAUSS;
	} else if ( !_stricmp( Config().Parameters().mSolvationModel, "GB" ) ) {
		m_iSolvModel = SOLV_GBORN;
	} else {
		Log().Fatal( "Unknown solvation model: `%s'\n", Config().Parameters().mSolvationModel );
	}
}

void cModel :: Clear( void )
{
	memset( &m_header, 0, sizeof(m_header) );
	m_atoms.clear();
	m_bonds.clear();
	m_SSbonds.clear();
	m_distRestraints.clear();
	m_angles.clear();
	m_impropers.clear();
	m_torsions.clear();
	m_nbPairsR1.clear();
	m_nbPairsR2.clear();
	m_hbTriplets.clear();
	
	if ( m_physAtoms ) {
		COM_AlignedFree( m_physAtoms );
		m_physAtoms = NULL;
	}
	if ( m_physPositions ) {
		COM_AlignedFree( m_physPositions );
		m_physPositions = NULL;
	}
	if ( m_physForces ) {
		COM_AlignedFree( m_physForces );
		m_physForces = NULL;
	}
	if ( m_physForcesMT ) {
		COM_AlignedFree( m_physForcesMT );
		m_physForcesMT = NULL;
	}

	delete [] m_charges;
	delete [] m_vdWpairs;
	delete [] m_HBpairs;
	delete[] m_connectivity;
	delete[] m_solvation;
	delete[] m_qeqparms;
	delete[] m_jValues;
	delete[] m_cValues;
	delete[] m_iValues;
	
	m_charges = NULL;
	m_vdWpairs = NULL;
	m_HBpairs = NULL;
	m_connectivity = NULL;
	m_solvation = NULL;
	m_qeqparms = NULL;
	m_jValues = NULL;
	m_cValues = NULL;
	m_iValues = NULL;
}

void cModel :: Load( void )
{
	const char *ext = COM_FileExtension( m_fileName );

	IOMap::iterator pReader = m_ioMap.find( ext );
	if ( pReader == m_ioMap.end() || !pReader->second )
		Log().Fatal( "Unsupported input structure file format: %s\n", ext );

	Clear();

	pReader->second->ReadFile( m_fileName, &m_header, &m_atoms );

	Log().TPrintf( "Model loaded: %i atom(s), %i residue(s), %i chain(s)\n", 
		m_header.mNumAtoms, m_header.mNumResidues, m_header.mNumChains );
}

void cModel :: Save( const char *const filename )
{
	char localfilename[MAX_OSPATH];
	const char *pfilename = filename;
	const char *ext = COM_FileExtension( filename );

	if ( m_modelName[0] ) {
		sprintf_s( localfilename, sizeof(localfilename), "%s.%s", m_modelName, filename );
		pfilename = localfilename;
	}

	IOMap::iterator pWriter = m_ioMap.find( ext );
	if ( pWriter == m_ioMap.end() || !pWriter->second )
		Log().Fatal( "Unsupported output structure file format: %s\n", ext );

	pWriter->second->WriteFile( pfilename, &m_header, &m_atoms );
}

void cModel :: ProfileStart( void )
{
	if ( m_bProfiling ) {
		assert( m_profDepth < MAX_PROFILER_DEPTH );
		m_profTime[m_profDepth++] = COM_Milliseconds();
	}
}

void cModel :: ProfileEnd( unsigned int *pCounter )
{
	if ( m_bProfiling ) {
		--m_profDepth;
		assert( m_profDepth >= 0 );
		*pCounter += ( COM_Milliseconds() - m_profTime[m_profDepth] );
	}
}

void cModel :: ProfileReport( void )
{
	if ( m_bProfiling ) {
		Log().TPrintf( "--- PERFORMANCE PROFILING RESULTS ---\n" );
		Log().DPrintf( " NBPL building:\t\t%8.3f ms\n", (vec_t)m_prof.mNBPLTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Bond energy:\t\t%8.3f ms\n", (vec_t)m_prof.mEBondTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Angle energy:\t\t%8.3f ms\n", (vec_t)m_prof.mEAngleTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Improper energy:\t%8.3f ms\n", (vec_t)m_prof.mEImproperTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Torsion energy:\t%8.3f ms\n", (vec_t)m_prof.mETorsionTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Non-bonded energy:\t%8.3f ms\n", (vec_t)m_prof.mENBTime / (vec_t)m_prof.mECounter );
		Log().DPrintf( " Special energy:\t%8.3f ms\n", (vec_t)m_prof.mESpecialTime / (vec_t)m_prof.mECounter );
	}
}
