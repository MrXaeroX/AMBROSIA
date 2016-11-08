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
#include "md_math.h"

char cConfig :: m_userConfigName[MAX_OSPATH] = "";
char cConfig :: m_parmsName[MAX_OSPATH] = DEFAULT_PARMS_FILENAME;

TYPEDESCRIPTION cConfig :: m_parmData[] = {
	DEFINE_PARAM_FIELD( "LoadModel",		mLoadModel,				FIELD_INTEGER ),
	DEFINE_PARAM_FIELD( "Hread",			mReadHydrogens,			FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "AutoSSBonds",		mAutoSSBonds,			FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "SSBondDist",		mSSBondDist,			FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "AutoQ",			mAutoQ,					FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "AutoRSN",			mAutoRSN,				FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "EngCalc",			mEnergyCalc,			FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "EngOptim",			mEnergyOptimize,		FIELD_FLAG ),
	DEFINE_PARAM_FIELD( "nOptStep",			mEnergyOptimizeSteps,	FIELD_INTEGER ),
	DEFINE_PARAM_FIELD( "ShortCutoffDist",	mShortCutoffDist,		FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "ShortSwitchDist",	mShortSwitchDist,		FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "LongCutoffDist",	mLongCutoffDist,		FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "SolventEConst",	mSolventDielectricConstant, FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "RDIE",				mDistanceDependentDielectricConstant, FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "DebyeRadius",		mDebyeRadius,			FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "HBondModel",		mHBondModel,			FIELD_INTEGER ),
	DEFINE_PARAM_FIELD( "HBondScale",		mHBondScale,			FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "HBondCutoffDist",	mHBondCutoffDist,		FIELD_FLOAT ),
	DEFINE_PARAM_FIELD( "SolvModel",		mSolvationModel,		FIELD_STRING ),
};

cConfig :: cConfig()
{
	memset( &m_parms, 0, sizeof(m_parms) );
	m_posRestraints.reserve( 16 );
	m_distRestraints.reserve( 16 );
	InitDefaults();
}

cConfig :: ~cConfig()
{
	FreeParams();
	Clear();
}

void cConfig :: InitDefaults( void )
{
	m_parms.mLoadModel = 0;
	m_parms.mSSBondDist = DEFAULT_SS_BOND_DISTANCE;
	m_parms.mShortCutoffDist = DEFAULT_SHORT_CUTOFF_DISTANCE;
	m_parms.mShortSwitchDist = DEFAULT_SHORT_SWITCH_DISTANCE;
	m_parms.mLongCutoffDist = DEFAULT_LONG_CUTOFF_DISTANCE;
	m_parms.mSolventDielectricConstant = DEFAULT_SOLVENT_DIELECTRIC_CONST;
	m_parms.mDistanceDependentDielectricConstant = DEFAULT_RDIE_CONST;
	m_parms.mDebyeRadius = DEFAULT_SOLVENT_DEBYE_RADIUS;
	m_parms.mHBondModel = DEFAULT_HBOND_MODEL;
	m_parms.mHBondScale = DEFAULT_HBOND_SCALE;
	m_parms.mHBondCutoffDist = DEFAULT_HBOND_CUTOFF_DISTANCE;
	m_parms.mSolvationModel = NULL;
	m_posRestraints.clear();
	m_distRestraints.clear();
}

void cConfig :: ValidateParams( void )
{
	if ( m_parms.mShortSwitchDist > m_parms.mShortCutoffDist ) {
		Log().Warning( "ShortSwitchDist > ShortCutoffDist!\n" );
		m_parms.mShortSwitchDist = m_parms.mShortCutoffDist;
	}
	if ( m_parms.mLongCutoffDist <= m_parms.mShortCutoffDist ) {
		Log().Warning( "LongCutoffDist <= ShortCutoffDist!\n" );
		m_parms.mLongCutoffDist = m_parms.mShortCutoffDist * 2;
	}
	if ( m_parms.mHBondModel != 0 && m_parms.mHBondModel != 126 && m_parms.mHBondModel != 128 ) {
		Log().Warning( "Unknown H-bond model (%i), using default (%i)\n", m_parms.mHBondModel, DEFAULT_HBOND_MODEL );
		m_parms.mHBondModel = DEFAULT_HBOND_MODEL;
	}
	if ( m_parms.mHBondScale < vec_t( 0.00001 ) ) {
		Log().Warning( "Bad H-bond energy scale (%f), using default (%f)\n", m_parms.mHBondScale, vec_t(DEFAULT_HBOND_SCALE) );
		m_parms.mHBondScale = DEFAULT_HBOND_SCALE;
	}
}

void cConfig :: SetUserConfigName( const char *name )
{
	strncpy_s( m_userConfigName, name, sizeof(m_userConfigName)-1 );
}

void cConfig :: SetParmsFileName( const char *name )
{
	strncpy_s( m_parmsName, name, sizeof(m_parmsName)-1 );
}

const char *cConfig :: LocationToString( eResidueLocation loc )
{
	switch ( loc ) {
	case RL_INT: return "int";
	case RL_BEG: return "beg";
	case RL_END: return "end";
	case RL_ISO: return "iso";
	default: return "???";
	}
}

void cConfig :: Clear( void )
{
	// Clear FF code info
	for ( FFCodeMap::iterator it = m_ffCodeInfo.begin(); it != m_ffCodeInfo.end(); ++it )
		delete [] it->first;
	m_ffCodeInfo.clear();

	// Clear SFF code info
	for ( FFCodeMap::iterator it = m_sffCodeInfo.begin(); it != m_sffCodeInfo.end(); ++it )
		delete [] it->first;
	m_sffCodeInfo.clear();

	// Clear atom info
	for ( AtomInfoMap::iterator it = m_atomInfo.begin(); it != m_atomInfo.end(); ++it )
		delete [] it->first;
	m_atomInfo.clear();

	// Clear residue info
	for ( ResidueInfoMap::iterator it = m_residueInfo.begin(); it != m_residueInfo.end(); ++it )
		delete [] it->first;
	m_residueInfo.clear();

	// Clear alias info
	for ( ResidueAliasMap::iterator it = m_aliasInfo.begin(); it != m_aliasInfo.end(); ++it )
		delete [] it->first;
	m_aliasInfo.clear();

	// Clear torsion info
	for ( TorsionInfoList::iterator it = m_torsionInfo.begin(); it != m_torsionInfo.end(); ++it ) {
		cTorsionHarm *pHarm = it->mpHarmonic;
		cTorsionHarm *pNextHarm;
		while ( pHarm ) {
			pNextHarm = pHarm->mpNext;
			delete pHarm;
			pHarm = pNextHarm;
		}
	}
	m_torsionInfo.clear();

	// Clear info
	m_bondInfo.clear();
	m_angleInfo.clear();
	m_improperInfo.clear();
	m_vdwInfo.clear();
	m_hbInfo.clear();
	m_hbHydrogens.clear();
	m_hbAcceptors.clear();
	m_solvGSInfo.clear();
}

void cConfig :: LoadConfig( const char *const name )
{
	FILE *fp;

	// Open the file
	if ( COM_fopen_local( &fp, name, "rb" ) )
		Log().Fatal( "Couldn't load config \"%s\"!\n", name );

	Log().TPrintf( "Loading config: \"%s\"\n", name );

	// Load the whole file into memory
	size_t dataSize = COM_FileLength( fp );
	if ( dataSize ) {
		byte *fileData = new byte[dataSize+1];
		memset( fileData, 0, dataSize+1 );
		size_t readBytes = fread( fileData, 1, dataSize, fp );
		if ( readBytes != dataSize )
			Log().Fatal( "File read error! (read %i bytes, expected %i bytes)\n", readBytes, dataSize );

		fileData[dataSize] = '\0';
		ParseConfig( fileData );
		delete [] fileData;
	}

	// Close the file
	fclose( fp );
}

void cConfig :: ParseConfig( byte *fileData )
{
	byte *current = fileData;
	bool hbParms = false;
	
	m_lineNumber = 1;

	// Parse config tokens
 	while ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) != NULL ) {
 		if ( !_stricmp( m_token, "atom" ) ) {
			current = ParseAtom( current );
		} else if ( !_stricmp( m_token, "qeqparm" ) ) {
			current = ParseQEq( current );
		} else if ( !_stricmp( m_token, "residue" ) ) {
			current = ParseResidue( current );
		} else if ( !_stricmp( m_token, "bond" ) ) {
			current = ParseBond( current );
		} else if ( !_stricmp( m_token, "angle" ) ) {
			current = ParseAngle( current );
		} else if ( !_stricmp( m_token, "improper" ) ) {
			current = ParseImproper( current );
		} else if ( !_stricmp( m_token, "torsion" ) ) {
			current = ParseTorsion( current );
		} else if ( !_stricmp( m_token, "vdwparm" ) ) {
			current = ParseVdWParms( current );
		} else if ( !_stricmp( m_token, "hbparm" ) ) {
			hbParms = true;
			current = ParseHBParms( current );
		} else if ( !_stricmp( m_token, "solvparm" ) ) {
			current = ParseSolvationParms( current );
		} else {
			// Undefined token
			Log().Warning( "Unexpected token at line #%i: `%s'\n", m_lineNumber, m_token );
		}
	}

	// Sort lists of known hydrogens and acceptors for a binary search
	if ( hbParms ) {
		std::sort( m_hbHydrogens.begin(), m_hbHydrogens.end() );
		std::sort( m_hbAcceptors.begin(), m_hbAcceptors.end() );
	}
}

void cConfig :: LoadUserConfig( void )
{
	// Load if file name is not empty
	if ( strlen( m_userConfigName ) )
		LoadConfig( m_userConfigName );
}

void cConfig :: FreeParams( void )
{
	int numParms = sizeof(m_parmData) / sizeof(m_parmData[0]);
	for ( int i = 0; i < numParms; ++i ) {
		if ( m_parmData[i].mFieldType == FIELD_STRING )
			delete [] (*(char**)((byte*)&m_parms + m_parmData[i].mFieldOffset ));
	}

	memset( &m_parms, 0, sizeof(m_parms) );
}

void cConfig :: LoadParams( void )
{
	FILE *fp;

	int numParms = sizeof(m_parmData) / sizeof(m_parmData[0]);
	for ( int i = 0; i < numParms; ++i ) {
		m_parmData[i].mDefault = true;
	}

	// Open the file
	if ( fopen_s( &fp, m_parmsName, "rb" ) )
		Log().Fatal( "Couldn't load parameters from \"%s\"!\n", m_parmsName );

	Log().TPrintf( "Loading parameters: \"%s\"\n", m_parmsName );

	// Load the whole file into memory
	size_t dataSize = COM_FileLength( fp );
	if ( dataSize ) {
		byte *fileData = new byte[dataSize+1];
		memset( fileData, 0, dataSize+1 );
		size_t readBytes = fread( fileData, 1, dataSize, fp );
		if ( readBytes != dataSize )
			Log().Fatal( "File read error! (read %i bytes, expected %i bytes)\n", readBytes, dataSize );

		fileData[dataSize] = '\0';
		ParseParams( fileData );
		delete [] fileData;
	}

	// Close the file
	fclose( fp );

	// Validate params
	ValidateParams();

	// Print parameters to the log
	PrintParams();
}

void cConfig :: PrintParams( void )
{
	int numParms = sizeof(m_parmData) / sizeof(m_parmData[0]);
	int counter;

	Log().DPrintf( "=== Current parameters: ===\n" );

	for ( int i = 0; i < numParms; ++i ) {
		switch ( m_parmData[i].mFieldType ) {
		default:
		case FIELD_INTEGER:
			Log().DPrintf( " %s = %i%s\n", m_parmData[i].mFieldName, *(int*)((byte*)&m_parms + m_parmData[i].mFieldOffset ), m_parmData[i].mDefault ? " (D)" : "" );
			break;
		case FIELD_FLOAT:
			Log().DPrintf( " %s = %f%s\n", m_parmData[i].mFieldName, *(vec_t*)((byte*)&m_parms + m_parmData[i].mFieldOffset ), m_parmData[i].mDefault ? " (D)" : "" );
			break;
		case FIELD_FLAG:
			Log().DPrintf( " %s = %i%s\n", m_parmData[i].mFieldName, *(int*)((byte*)&m_parms + m_parmData[i].mFieldOffset ), m_parmData[i].mDefault ? " (D)" : "" );
			break;
		case FIELD_STRING:
			Log().DPrintf( " %s = '%s'%s\n", m_parmData[i].mFieldName, *(char**)((byte*)&m_parms + m_parmData[i].mFieldOffset ), m_parmData[i].mDefault ? " (D)" : "" );
			break;
		}
	}

	if ( m_posRestraints.size() ) {
		Log().DPrintf( "=== Positional restraints: ===\n" );
		counter = 1;
		for ( PosRestrVector::const_iterator it = m_posRestraints.begin(); it != m_posRestraints.end(); ++it, ++counter )
			Log().DPrintf( " %i) Chain %c: residues %i-%i (Kh = %f, type = %s)\n", counter, 'A' + it->mChain - 1, it->mFirstResidue, it->mLastResidue, it->mHarmConst, it->mBackboneOnly ? "backbone" : "all" );
	}
	if ( m_distRestraints.size() ) {
		Log().DPrintf( "=== Distant restraints: ===\n" );
		counter = 1;
		for ( DistRestrVector::const_iterator it = m_distRestraints.begin(); it != m_distRestraints.end(); ++it, ++counter )
			Log().DPrintf( " %i) %c-%i-%-4s : %c-%i-%-4s (Kh = %f, dist = %f)\n", counter,
							'A' + it->mChain[0] - 1, it->mResidue[0], it->mAtomTitle[0],
							'A' + it->mChain[1] - 1, it->mResidue[1], it->mAtomTitle[1],
							it->mHarmConst, it->mDistance );
	}
}

int cConfig :: ChainIndexFromToken( const char *token )
{
	if ( !token || !*token )
		return 0;

	if ( *token >= 'A' && *token <= 'Z' )
		return *token - 'A' + 1;
	else if ( *token >= 'a' && *token <= 'z' )
		return *token - 'a' + 1;
	else if ( *token >= '0' && *token <= '9' )
		return atoi( token );
	else {
		Log().Warning( "ChainIndexFromToken: bad token `%s'\n", token );
		return 1;
	}
}

void cConfig :: ParseParams( byte *fileData )
{
	char parmName[128];
	byte *current = fileData;
	int numParms = sizeof(m_parmData) / sizeof(m_parmData[0]);
	int i;
	
	m_lineNumber = 1;

	FreeParams();
	InitDefaults();

	// Parse config tokens
 	while ( ( current = COM_Parse( current, true, parmName, sizeof(parmName), m_lineNumber ) ) != NULL ) {
		if ( parmName[0] != '$' ) {
			Log().Error( "Syntax error at line #%i: parameters must begin with `$'\n", m_lineNumber );
			break;
		}

		// Check for restraint blocks
		if ( !_stricmp( parmName + 1, "PosRestraint" ) ) {
			// parse positional restraints
			current = COM_MatchToken( current, true, "{", m_lineNumber );
			while ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) != NULL ) {
				if ( !_stricmp( m_token, "}" ) )
					break;
				cPosRestraintInfo ri;
				ri.mChain = ChainIndexFromToken( m_token );
				// parse residue indices
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				ri.mFirstResidue = atoi( m_token );
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				ri.mLastResidue = atoi( m_token );
				// parse harmonic constant
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				ri.mHarmConst = atof( m_token );
				// parse type
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				if ( !_stricmp( m_token, "all" ) )
					ri.mBackboneOnly = false;
				else if ( !_stricmp( m_token, "backbone" ) )
					ri.mBackboneOnly = true;
				else {
					ri.mBackboneOnly = false;
					Log().Warning( "Unknown positional restraint type at line #%i: `%s'\n", m_lineNumber, m_token );
				}
				m_posRestraints.push_back( ri );
			}
			continue;
		}
		else if ( !_stricmp( parmName + 1, "DistRestraint" ) ) {
			// parse distant restraints
			current = COM_MatchToken( current, true, "{", m_lineNumber );
			while ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) != NULL ) {
				if ( !_stricmp( m_token, "}" ) )
					break;
				cDistRestraintInfo ri;
				for ( i = 0; i < 2; ++i ) { 
					// parse chain identifier
					if ( i > 0 ) {
						if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
							Log().Error( "Line #%i is incomplete\n", m_lineNumber );
							break;
						}
					}
					ri.mChain[i] = ChainIndexFromToken( m_token );
					// parse residue index
					if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
						Log().Error( "Line #%i is incomplete\n", m_lineNumber );
						break;
					}
					ri.mResidue[i] = atoi( m_token );
					// parse atom title
					if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
						Log().Error( "Line #%i is incomplete\n", m_lineNumber );
						break;
					}
					strncpy_s( ri.mAtomTitle[i], sizeof(ri.mAtomTitle[i]), m_token, sizeof(ri.mAtomTitle[i])-1 );
				}
				// parse harmonic constant
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				ri.mHarmConst = atof( m_token );
				// parse distance
				if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
					Log().Error( "Line #%i is incomplete\n", m_lineNumber );
					break;
				}
				ri.mDistance = atof( m_token );
				if ( ri.mDistance < DEFAULT_MIN_DIST_RESTR_LENGTH ) {
					Log().Warning( "Distant restraint length is too low at line #%i: `%f'\n", m_lineNumber, ri.mDistance );
					ri.mDistance = DEFAULT_MIN_DIST_RESTR_LENGTH;
				}
				m_distRestraints.push_back( ri );
			}
			continue;
		}

		// Parse value, if any
		if ( COM_TokenAvailable( current, false ) ) {
			current = COM_MatchToken( current, true, "=", m_lineNumber );
			if ( !current )
				break;
			if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
				Log().Error( "Line #%i is incomplete\n", m_lineNumber );
				break;
			}
		}

		// Find parm and assign the value
		for ( i = 0; i < numParms; ++i ) {
			if ( !_stricmp( parmName + 1, m_parmData[i].mParmName ) ) {
				switch ( m_parmData[i].mFieldType ) {
				default:
				case FIELD_INTEGER:
					*(int*)((byte*)&m_parms + m_parmData[i].mFieldOffset ) = atoi( m_token );
					break;
				case FIELD_FLOAT:
					*(vec_t*)((byte*)&m_parms + m_parmData[i].mFieldOffset ) = COM_Atof( m_token );
					break;
				case FIELD_FLAG:
					*(int*)((byte*)&m_parms + m_parmData[i].mFieldOffset ) = 1;
					break;
				case FIELD_STRING:
					*(char**)((byte*)&m_parms + m_parmData[i].mFieldOffset ) = COM_AllocString( m_token );
					break;
				}
				m_parmData[i].mDefault = false;
				break;
			}
		}

		// Check if parm was not found
		if ( i == numParms ) {
			Log().Warning( "Unknown parameter at line #%i: `%s'\n", m_lineNumber, parmName );
		}
	}
}

byte *cConfig :: ParseAtom( byte *fileData )
{
	cAtomInfo atomInfo;
	char key[2][16];
	byte *current = fileData;

	// Parse atom definition
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}

	atomInfo.mMassReciprocal = COM_Atof( m_token );
	if ( atomInfo.mMassReciprocal != 0 )
		atomInfo.mMassReciprocal = vec_t(1) / atomInfo.mMassReciprocal;
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
		
	atomInfo.mRadius = COM_Atof( m_token );

	// clear qEq parms
	atomInfo.mXi = 0;
	atomInfo.mJ0 = 0;

	Log().Verbose( "Atom:\t'%-2s'\t\tinvmass = %f\tradius = %f\n", key, atomInfo.mMassReciprocal, atomInfo.mRadius );
	m_atomInfo[COM_AllocString(key[0])] = atomInfo;
	return current;
}

byte *cConfig :: ParseQEq( byte *fileData )
{
	cAtomInfo *pAtomInfo = NULL;
	char key[16];
	byte *current = fileData;

	// Parse atom definition
	if ( ( current = COM_Parse( current, false, key, sizeof(key), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key );
	pAtomInfo = const_cast<cAtomInfo*>(LookupAtomInfo( key ));
	if ( !pAtomInfo ) {
		Log().Error( "Line #%i: undefined atom `%s'\n", m_lineNumber, key );
		return NULL;
	}

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	pAtomInfo->mXi = COM_Atof( m_token );

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	pAtomInfo->mJ0 = COM_Atof( m_token );

	Log().Verbose( "Atom:\t'-%2s'\t\txi = %f\t\tj0 = %f\n", key, pAtomInfo->mXi, pAtomInfo->mJ0 );
	return current;
}

byte *cConfig :: ParseResidueLocation( cResidueLocation *pLoc, byte *fileData )
{
	byte *current = fileData;
	vec_t totalChargeCalc = 0.0f;
	int topoFlags = 0;
	int numHeavyAtoms = 0;

	// Read location parameters
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	pLoc->mNumAtoms = atoi( m_token );
	if ( !pLoc->mNumAtoms ) {
		Log().Error( "invalid number of atoms in the location: %i\n", pLoc->mNumAtoms );
		return NULL;
	}
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	pLoc->mTotalCharge = COM_Atof( m_token );

	pLoc->mNumLoops = 0;

	// Start reading the block
	current = COM_MatchToken( current, true, "{", m_lineNumber );
	if ( !current )
		return NULL;

	// Read atoms
	for ( int i = 0; i < pLoc->mNumAtoms; ++i ) {
		cResidueAtom rAtom;
		memset( &rAtom, 0, sizeof(rAtom) );

		if ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		rAtom.mIndex = atoi( m_token );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		strncpy_s( rAtom.mTitle, m_token, sizeof(rAtom.mTitle)-1 );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		strncpy_s( rAtom.mFFTitle, m_token, sizeof(rAtom.mFFTitle)-1 );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		if ( m_token[0] != 'D' && m_token[0] != 'B' && m_token[0] != 'S' && m_token[0] != 'R' && m_token[0] != 'H' )
			Log().Warning( "Line #%i: unknown atom type `%c'\n", m_lineNumber, m_token[0] );
		rAtom.mType = m_token[0];

		if ( m_token[0] == 'B' || m_token[0] == 'S' || m_token[0] == 'R' )
			++numHeavyAtoms;

		// First backbone atom is start of the topology
		if ( m_token[0] == 'B' && !( topoFlags & TF_TOPOLOGY_BEGIN ) ) {
			rAtom.mFlags |= TF_TOPOLOGY_BEGIN;
			topoFlags |= TF_TOPOLOGY_BEGIN;
		}
		if ( m_token[0] == 'R' ) {
			rAtom.mFlags |= TF_TOPOLOGY_RING;
			topoFlags |= TF_TOPOLOGY_RING;
		}

		for ( int j = 0; j < 3; ++j ) {
			if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
				Log().Error( "Line #%i is incomplete\n", m_lineNumber );
				return NULL;
			}
			rAtom.mPrevious[j] = atoi( m_token );
			if ( rAtom.mPrevious[j] >= rAtom.mIndex )
				Log().Warning( "Error in topology at line #%i: previous atom %i >= current!\n", m_lineNumber, j );
		}

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		rAtom.mR = COM_Atof( m_token );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		// In config, there is a bond angle
		// We need its complement to 180 deg
		rAtom.mTheta = DEG2RAD( 180 - COM_Atof( m_token ) );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		rAtom.mPhi = DEG2RAD( COM_Atof( m_token ) );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		rAtom.mCharge = COM_Atof( m_token );
		totalChargeCalc += rAtom.mCharge;

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		strncpy_s( rAtom.mSymbol, m_token, sizeof(rAtom.mSymbol)-1 );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		strncpy_s( rAtom.mSFFTitle, m_token, sizeof(rAtom.mSFFTitle)-1 );

		if ( strcmp( rAtom.mSymbol, "H" ) && strcmp( rAtom.mSymbol, "-" ) && !strcmp( rAtom.mSFFTitle, "-" ) )
			Log().Warning( "Heavy atom `%s' without solvation force field code!\n", rAtom.mTitle );

		pLoc->mAtoms.push_back( rAtom );
	}

	// Read loops, if any
	while ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) != NULL ) {
		if ( m_token[0] == '}' )
			break;

		if ( !_stricmp( m_token, "loop" ) ) {
			// Parse loop definition
			// Atom1
			if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
				Log().Error( "Line #%i is incomplete\n", m_lineNumber );
				break;
			}
			ResidueAtomArray::iterator at1 = std::find_if( pLoc->mAtoms.begin(), pLoc->mAtoms.end(), ResidueAtomSearchByTitleFunctor( m_token ) );
			if ( at1 == pLoc->mAtoms.end() ) {
				Log().Error( "Loop error at line #%i: atom `%s' not found\n", m_lineNumber, m_token );
				break;
			}
			if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
				Log().Error( "Line #%i is incomplete\n", m_lineNumber );
				break;
			}
			ResidueAtomArray::iterator at2 = std::find_if( pLoc->mAtoms.begin(), pLoc->mAtoms.end(), ResidueAtomSearchByTitleFunctor( m_token ) );
			if ( at2 == pLoc->mAtoms.end() ) {
				Log().Error( "Loop error at line #%i: atom `%s' not found\n", m_lineNumber, m_token );
				break;
			}
			if ( at1 == at2 ) {
				Log().Error( "Loop error at line #%i: atom is connected to itself\n", m_lineNumber );
				break;
			}

			// Check if too many loops
			if ( pLoc->mNumLoops == 32 ) {
				Log().Error( "Loop error at line #%i: too many loops (limit is 32)!\n", m_lineNumber );
				break;
			}

			// Set loop bits
			int loopBit = ( 1 << pLoc->mNumLoops );
			at1->mLoopFlags |= loopBit;
			at2->mLoopFlags |= loopBit;
			pLoc->mNumLoops++;
		} else {
			// Undefined token
			Log().Warning( "Unexpected token at line #%i: `%s'\n", m_lineNumber, m_token );
		}
	}

	if ( (vec_t)fabs( totalChargeCalc - pLoc->mTotalCharge ) > ON_EPSILON )
		Log().Warning( "Total charge (%f) doesn't match sum of partial charges (%f)\n", pLoc->mTotalCharge, totalChargeCalc );

	// Mark last atom at the end of the topology
	for ( ResidueAtomArray::reverse_iterator it = pLoc->mAtoms.rbegin(); it != pLoc->mAtoms.rend(); ++it ) {
		if ( it->mType == 'B' ) {
			if ( ( numHeavyAtoms > 1 ) && ( it->mFlags & TF_TOPOLOGY_BEGIN ) )
				Log().Warning( "Only one backbone atom in the topology!\n" );
			it->mFlags |= TF_TOPOLOGY_END;
			break;
		}
	}

	return current;
}

void cConfig :: DumpResidue( const char *title, const cResidueInfo *pInfo )
{
	static FILE *fp = NULL;
	static FILE *fp2 = NULL;
	static const char *locNames[] = { "INT", "BEG", "END", "ISO" };
	std::vector<std::pair<int,int> > loopInfo;

	if ( !fp ) {
		fopen_s( &fp, "topology.dat", "w" );
		if ( !fp ) return;
	}
	if ( !fp2 ) {
		fopen_s( &fp2, "solvGS.dat", "w" );
		if ( !fp2 ) return;
	}

	// Export topology in BioPASED format
	for ( int i = 0; i < RL_MAX; ++i ) {
		if ( !pInfo->mLocation[i].mNumAtoms )
			continue;

		loopInfo.clear();

		fprintf( fp, "$MTRES\n" );
		fprintf( fp, "%-4s %3s  NA %5i   %.1f\n", title, locNames[i], pInfo->mLocation[i].mNumAtoms, pInfo->mLocation[i].mTotalCharge );

		if ( i == RL_INT )
			fprintf( fp2, "!%s\n", title );

		for ( int j = 0; j < pInfo->mLocation[i].mNumAtoms; ++j ) {
			const cResidueAtom *at = &pInfo->mLocation[i].mAtoms.at( j );
			fprintf( fp, "%4i  %-4s  %-2s    %c %4i%4i%4i%10.3f%10.3f%10.3f%10.5f\n",
				j + 1,at->mTitle, !_stricmp( at->mFFTitle, "CI" ) ? "CT" : at->mFFTitle, 
				at->mType, at->mPrevious[0], at->mPrevious[1], at->mPrevious[2],
				at->mR, 180 - RAD2DEG(at->mTheta), RAD2DEG(at->mPhi), at->mCharge );
			if ( at->mLoopFlags ) {
				for ( int k = 0; k < 32; ++k ) {
					if ( at->mLoopFlags & ( 1 << k ) ) {
						for ( int j2 = j + 1; j2 < pInfo->mLocation[i].mNumAtoms; ++j2 ) {
							const cResidueAtom *at2 = &pInfo->mLocation[i].mAtoms.at( j2 );
							if ( at2->mLoopFlags & ( 1 << k ) ) {
								loopInfo.push_back( std::make_pair( j + 1, j2 + 1 ) );
								break;
							}
						}
					}
				}
			}
			if ( i == RL_INT ) {
				if ( at->mSFFTitle[0] && at->mSFFTitle[0] != '-' && _stricmp( at->mSymbol, "P" ) ) {
					int solvCode = 0;
					if ( !_stricmp( at->mSFFTitle, "CC" ) || !_stricmp( at->mSFFTitle, "CT" ) )
						solvCode = 1;
					else if ( !_stricmp( at->mSFFTitle, "CR" ) || !_stricmp( at->mSFFTitle, "CY" ) )
						solvCode = 2;
					else if ( !_stricmp( at->mSFFTitle, "C1" ) )
						solvCode = 3;
					else if ( !_stricmp( at->mSFFTitle, "C2" ) )
						solvCode = 4;
					else if ( !_stricmp( at->mSFFTitle, "C3" ) )
						solvCode = 5;
					else if ( !_stricmp( at->mSFFTitle, "CA" ) || !_stricmp( at->mSFFTitle, "CX" ) )
						solvCode = 6;
					else if ( !_stricmp( at->mSFFTitle, "N1" )  || !_stricmp( at->mSFFTitle, "NX" ) )
						solvCode = 7;
					else if ( !_stricmp( at->mSFFTitle, "NR" ) )
						solvCode = 8;
					else if ( !_stricmp( at->mSFFTitle, "N2" ) )
						solvCode = 9;
					else if ( !_stricmp( at->mSFFTitle, "N3" ) )
						solvCode = 10;
					else if ( !_stricmp( at->mSFFTitle, "NG" ) )
						solvCode = 11;
					else if ( !_stricmp( at->mSFFTitle, "NP" ) )
						solvCode = 12;
					else if ( !_stricmp( at->mSFFTitle, "OH" ) )
						solvCode = 13;
					else if ( !_stricmp( at->mSFFTitle, "O1" ) )
						solvCode = 14;
					else if ( !_stricmp( at->mSFFTitle, "O2" ) )
						solvCode = 15;
					else if ( !_stricmp( at->mSFFTitle, "S0" ) || !_stricmp( at->mSFFTitle, "SM" ) )
						solvCode = 16;
					else if ( !_stricmp( at->mSFFTitle, "S1" ) )
						solvCode = 17;
					else
						Log().Error( "Funny SFF code `%s' in residue `%s'\n", at->mSFFTitle, title );

					fprintf( fp2, "%-8s %4s%8i\n", at->mTitle, title, solvCode );
				}
			}
		}

		if ( loopInfo.size() ) {
			fprintf( fp, "LOOP%4i\n",  loopInfo.size() );
			for ( std::vector<std::pair<int,int> >::const_iterator it = loopInfo.begin(); it != loopInfo.end(); ++it ) {
				fprintf( fp, "D:%4i%4i\n", it->first, it->second );
			}
		}

		fprintf( fp, "DONE\n" );
	}

	fflush( fp );
	fflush( fp2 );
}

byte *cConfig :: ParseResidue( byte *fileData )
{
	cResidueInfo residueInfo;
	char title[16];
	char alias[16];
	char *allocTitle = NULL;
	byte *current = fileData;
	int numLocations = 0;

	// Parse residue topology
	for ( int i = 0; i < RL_MAX; ++i ) {
		residueInfo.mLocation[i].mNumAtoms = 0;
		residueInfo.mLocation[i].mTotalCharge = 0;
	}

	// Read title and check its length
	if ( ( current = COM_Parse( current, false, title, sizeof(title), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( title );
	if ( strlen( title ) > 4 )
		Log().Warning( "title of the residue is too long: `%s'\n", title );

	allocTitle = COM_AllocString(title);

	// Start reading the block
	current = COM_MatchToken( current, true, "{", m_lineNumber );
	if ( !current )
		return NULL;

	// Read the next token
	while ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) != NULL ) {
		// Check if residue is complete
		if ( m_token[0] == '}' ) {
			if ( !numLocations ) {
				Log().Warning( "'%s': residue topology defines no locations!\n", title );
			} else {
				if ( Log().VerboseMode() ) {
					Log().DPrintf( "Residue: %s (%i locations)\n", title, numLocations );
					for ( int i = 0; i < RL_MAX; ++i ) {
						Log().DPrintf( " location '%s': %i atoms (%f charge)\n", LocationToString( (eResidueLocation)i ), residueInfo.mLocation[i].mNumAtoms, residueInfo.mLocation[i].mTotalCharge );
						for ( int j = 0; j < residueInfo.mLocation[i].mNumAtoms; ++j ) {
							const cResidueAtom *at = &residueInfo.mLocation[i].mAtoms.at( j );
							Log().DPrintf( "  %5i %4s %5s %c %4i %4i %4i %8.3f %8.3f %8.3f %8.5f %2s %5s %i %i", 
								at->mIndex, at->mTitle, at->mFFTitle, at->mType,
								at->mPrevious[0], at->mPrevious[1], at->mPrevious[2],
								at->mR, at->mTheta, at->mPhi, at->mCharge, at->mSymbol, at->mSFFTitle, at->mFlags, at->mLoopFlags );
							Log().DPrintf( "\n" );
						}
					}
				}
				//DumpResidue( allocTitle, &residueInfo );
				m_residueInfo[allocTitle] = residueInfo;
			}
			return current;
		}
		// Check for new location
		else if ( !_stricmp( m_token, "int" ) ) {
			current = ParseResidueLocation( &residueInfo.mLocation[RL_INT], current );
			++numLocations;
		} 
		else if ( !_stricmp( m_token, "beg" ) ) {
			current = ParseResidueLocation( &residueInfo.mLocation[RL_BEG], current );
			++numLocations;
		}
		else if ( !_stricmp( m_token, "end" ) ) {
			current = ParseResidueLocation( &residueInfo.mLocation[RL_END], current );
			++numLocations;
		}
		else if ( !_stricmp( m_token, "iso" ) ) {
			current = ParseResidueLocation( &residueInfo.mLocation[RL_ISO], current );
			++numLocations;
		}
		else if ( !_stricmp( m_token, "alias" ) ) {
			// Read alias and check its length
			if ( ( current = COM_Parse( current, false, alias, sizeof(alias), m_lineNumber ) ) == NULL ) {
				Log().Error( "Line #%i is incomplete\n", m_lineNumber );
				return NULL;
			}
			_strupr_s( alias );
			if ( strlen( alias ) > 4 )
				Log().Warning( "alias of the residue is too long: `%s'\n", alias );

			// Check if there is already an alias
			ResidueAliasMap::const_iterator ita = m_aliasInfo.find( alias );
			if ( ita != m_aliasInfo.end() ) {
				Log().Error( "alias `%s' already in use!\n", alias );
			} else {
				m_aliasInfo[COM_AllocString(alias)] = allocTitle;
				if ( Log().VerboseMode() )
					Log().DPrintf( "Alias: `%s' -> `%s'\n", alias, allocTitle );
			}
		}
		else {
			Log().Warning( "Unexpected token at line #%i: `%s'\n", m_lineNumber, m_token );
		}
	}
		
	Log().Error( "Line #%i is incomplete\n", m_lineNumber );
	return NULL;
}

byte *cConfig :: ParseBond( byte *fileData )
{
	cBondInfo bondInfo;
	char key[2][16];
	int ffCode1, ffCode2, ffTemp = 0;
	byte *current = fileData;

	// Parse bond definition
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	ffCode1 = FetchFFCode( key[0] );

	if ( ( current = COM_Parse( current, false, key[1], sizeof(key[1]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[1] );
	ffCode2 = FetchFFCode( key[1] );
	
	if ( ffCode1 > ffCode2 ) {
		ffTemp = ffCode2;
		ffCode2 = ffCode1;
		ffCode1 = ffTemp;
		ffTemp = 1;
	}
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	
	bondInfo.mHarmonicConstant = COM_Atof( m_token );
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	
	bondInfo.mLength = COM_Atof( m_token );
	
	Log().Verbose( "Bond:\t%s(%i)-%s(%i)\tkharm = %f\tlength = %f\n", key[ffTemp], ffCode1, key[ffTemp^1], ffCode2, bondInfo.mHarmonicConstant, bondInfo.mLength );
	m_bondInfo[std::make_pair(ffCode1, ffCode2)] = bondInfo;
	return current;
}

byte *cConfig :: ParseAngle( byte *fileData )
{
	cAngleInfo angleInfo;
	char key[3][16];
	int ffCode1, ffCode2, ffCode3, ffTemp = 0;
	byte *current = fileData;

	// Parse angle definition
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	ffCode1 = FetchFFCode( key[0] );

	if ( ( current = COM_Parse( current, false, key[2], sizeof(key[2]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[2] );
	ffCode3 = FetchFFCode( key[2] );

	if ( ( current = COM_Parse( current, false, key[1], sizeof(key[1]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[1] );
	ffCode2 = FetchFFCode( key[1] );
	
	if ( ffCode1 > ffCode2 ) {
		ffTemp = ffCode2;
		ffCode2 = ffCode1;
		ffCode1 = ffTemp;
		ffTemp = 1;
	}
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	angleInfo.mHarmonicConstant = COM_Atof( m_token );
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	angleInfo.mTheta = DEG2RAD( COM_Atof( m_token ) );
	
	Log().Verbose( "Angle:\t%s(%i)-%s(%i)-%s(%i)\tkharm = %f\ttheta = %f\n", key[ffTemp], ffCode1, key[2], ffCode3, key[ffTemp^1], ffCode2, angleInfo.mHarmonicConstant, angleInfo.mTheta );

	m_angleInfo[std::make_pair(ffCode3, std::make_pair(ffCode1,ffCode2))] = angleInfo;
	return current;
}

byte *cConfig :: ParseImproper( byte *fileData )
{
	cImproperInfo improperInfo;
	char key[4][16];
	int ffCodeI, ffCodeJ, ffCodeK, ffCodeL = 0;
	byte *current = fileData;

	// Parse improper definition
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	ffCodeI = FetchFFCode( key[0] );

	if ( ( current = COM_Parse( current, false, key[1], sizeof(key[1]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[1] );
	ffCodeJ = FetchFFCode( key[1] );

	if ( ( current = COM_Parse( current, false, key[2], sizeof(key[2]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[2] );
	ffCodeK = FetchFFCode( key[2] );

	if ( ( current = COM_Parse( current, false, key[3], sizeof(key[3]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[3] );
	ffCodeL = FetchFFCode( key[3] );
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	improperInfo.mHarmonicConstant = COM_Atof( m_token ) * vec_t( 0.5 );
	
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	
	improperInfo.mXi = DEG2RAD( COM_Atof( m_token ) );
	
	Log().Verbose( "Improper:\t%s(%i)-%s(%i)-%s(%i)-%s(%i)\tkharm = %f\txi = %f\n", 
		key[0], ffCodeI, key[1], ffCodeJ, key[2], ffCodeK, key[3], ffCodeL, 
		improperInfo.mHarmonicConstant, improperInfo.mXi );

	m_improperInfo[std::make_pair(std::make_pair(ffCodeI,ffCodeJ), std::make_pair(ffCodeK,ffCodeL))] = improperInfo;
	return current;
}

byte *cConfig :: ParseTorsion( byte *fileData )
{
	cTorsionInfo torsionInfo;
	char key[4][16];
	int numHarmonics = 0;
	byte *current = fileData;
	byte *store;
	vec_t phase;

	// Parse torsion definition
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	torsionInfo.mFFCode[0] = FetchFFCode( key[0] );

	if ( ( current = COM_Parse( current, false, key[1], sizeof(key[1]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[1] );
	torsionInfo.mFFCode[1] = FetchFFCode( key[1] );

	if ( ( current = COM_Parse( current, false, key[2], sizeof(key[2]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[2] );
	torsionInfo.mFFCode[2] = FetchFFCode( key[2] );

	if ( ( current = COM_Parse( current, false, key[3], sizeof(key[3]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[3] );
	torsionInfo.mFFCode[3] = FetchFFCode( key[3] );

	torsionInfo.mpHarmonic = NULL;
	cTorsionHarm *pLastHarm = NULL;

	while ( current != NULL ) {
		// Parse all harmonics
		store = current;
	
		if ( ( current = COM_Parse( current, true, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}

		// if not a digit, then no more harmonics
		if ( m_token[0] < '0' || m_token[0] > '9' ) {
			current = store;
			break;
		}

		// allocate harmonic
		cTorsionHarm *pHarm = new cTorsionHarm;
		pHarm->mpNext = NULL;
		if ( !pLastHarm )
			torsionInfo.mpHarmonic = pHarm;
		else
			pLastHarm->mpNext = pHarm;
		pLastHarm = pHarm;

		pHarm->mHarmonicScale = vec_t( 1.0 ) / (vec_t)atoi( m_token );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		pHarm->mHarmonicConstant = COM_Atof( m_token );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		phase = DEG2RAD( COM_Atof( m_token ) );
		pHarm->mPhaseCos = cos( phase );
		pHarm->mPhaseSin = sqrt( 1 - pHarm->mPhaseCos*pHarm->mPhaseCos);

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		pHarm->mPeriodicity = atoi( m_token );
		if ( pHarm->mPeriodicity < 1 || pHarm->mPeriodicity > 6 ) {
			Log().Error( "Line #%i: periodicity must be in range [1,6]!\n", m_lineNumber );
			pHarm->mPeriodicity = 1;
		}

		++numHarmonics;
	}

	Log().Verbose( "Torsion:\t%s(%i)-%s(%i)-%s(%i)-%s(%i) [%i harmonics]", 
		key[0], torsionInfo.mFFCode[0], key[1], torsionInfo.mFFCode[1], key[2], torsionInfo.mFFCode[2], key[3], torsionInfo.mFFCode[3],
		numHarmonics );
	pLastHarm = torsionInfo.mpHarmonic;
	while ( pLastHarm ) {
		Log().Verbose( "\tscale = %f\tkharm = %f\tphaseCos = %f\tphaseSin = %f\tperiodicity = %i", 
						pLastHarm->mHarmonicScale, pLastHarm->mHarmonicConstant, pLastHarm->mPhaseCos, pLastHarm->mPhaseSin, pLastHarm->mPeriodicity );
		pLastHarm = pLastHarm->mpNext;
	}
	Log().Verbose( "\n" );

	// Torsions with "any" atoms must be at the end, so we first return a specialized torsion
	if ( torsionInfo.mFFCode[0] >= 0 && torsionInfo.mFFCode[1] >= 0 && torsionInfo.mFFCode[2] >= 0 && torsionInfo.mFFCode[3] >= 0 ) 
		m_torsionInfo.push_front( torsionInfo );
	else
		m_torsionInfo.push_back( torsionInfo );
	return current;
}

byte *cConfig :: ParseVdWParms( byte *fileData )
{
	cVdWInfo vdwInfo;
	char key[16];
	int ffCode;
	byte *current = fileData;

	// Parse FF code
	if ( ( current = COM_Parse( current, false, key, sizeof(key), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key );
	ffCode = FetchFFCode( key );

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	vdwInfo.mR = COM_Atof( m_token );

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	vdwInfo.mE = COM_Atof( m_token );

	Log().Verbose( "VdWParm:\t%s\tR = %f\tE = %f\n", 
		key, vdwInfo.mR, vdwInfo.mE );

	m_vdwInfo[ffCode] = vdwInfo;
	return current;
}

byte *cConfig :: ParseHBParms( byte *fileData )
{
	cHBInfo hbInfo;
	char key[2][16];
	int ffCode1, ffCode2;
	byte *current = fileData;
	HBCodeVector::const_iterator vi;

	// Parse FF codes
	if ( ( current = COM_Parse( current, false, key[0], sizeof(key[0]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[0] );
	ffCode1 = FetchFFCode( key[0] );
	if ( ( current = COM_Parse( current, false, key[1], sizeof(key[1]), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( key[1] );
	ffCode2 = FetchFFCode( key[1] );

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	hbInfo.mR = COM_Atof( m_token );

	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	hbInfo.mE = COM_Atof( m_token );

	Log().Verbose( "HBParm:\t%s-%s\tR = %f\tE = %f\n", 
		key[0], key[1], hbInfo.mR, hbInfo.mE );

	vi = std::find( m_hbHydrogens.begin(), m_hbHydrogens.end(), ffCode1 );
	if ( vi == m_hbHydrogens.end() ) m_hbHydrogens.push_back( ffCode1 );

	vi = std::find( m_hbAcceptors.begin(), m_hbAcceptors.end(), ffCode2 );
	if ( vi == m_hbAcceptors.end() ) m_hbAcceptors.push_back( ffCode2 );

	m_hbInfo[std::make_pair( ffCode1, ffCode2 )] = hbInfo;
	return current;
}

byte* cConfig :: ParseSolvationParms( byte *fileData )
{
	char key[16];
	byte *current = fileData;

	// Parse solvation model name
	if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
		Log().Error( "Line #%i is incomplete\n", m_lineNumber );
		return NULL;
	}
	_strupr_s( m_token );

	if ( !strcmp( m_token, "GS" ) ) {
		// Parse Gaussian solvation parameters
		int sffCode;
		cSolvGSInfo sInfo;

		// Parse FF code
		if ( ( current = COM_Parse( current, false, key, sizeof(key), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		_strupr_s( key );
		sffCode = FetchSFFCode( key );

		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		sInfo.mVolume = COM_Atof( m_token );
		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		sInfo.mGref = COM_Atof( m_token );
		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		sInfo.mGfree = COM_Atof( m_token );
		if ( ( current = COM_Parse( current, false, m_token, sizeof(m_token), m_lineNumber ) ) == NULL ) {
			Log().Error( "Line #%i is incomplete\n", m_lineNumber );
			return NULL;
		}
		sInfo.mLambda = COM_Atof( m_token );

		Log().Verbose( "SolvGSParm:\t%s\tV = %f\tGref = %f\tGfree = %f\tLambda = %f\n", 
			key, sInfo.mVolume, sInfo.mGref, sInfo.mGfree, sInfo.mLambda );

		m_solvGSInfo[sffCode] = sInfo;
	} else {
		Log().Error( "Unknown solvation model name `%s' as line #%i\n", m_token, m_lineNumber );
		return NULL;
	}

	return current;
}

int cConfig :: FetchFFCode( const char *const title )
{
	// Check for "any" atom
	if ( !strcmp( title, "X" ) )
		return -1;

	FFCodeMap::const_iterator it = m_ffCodeInfo.find( title );
	if ( it == m_ffCodeInfo.end() ) {
		// Register a new FF code
		int newCode = (int)m_ffCodeInfo.size();
		m_ffCodeInfo[COM_AllocString(title)] = newCode;
		return newCode;
	}
	return it->second;
}

int cConfig :: FetchSFFCode( const char *const title )
{
	if ( !strcmp( title, "-" ) )
		return -1;

	FFCodeMap::const_iterator it = m_sffCodeInfo.find( title );
	if ( it == m_sffCodeInfo.end() ) {
		// Register a new SFF code
		int newCode = (int)m_sffCodeInfo.size();
		m_sffCodeInfo[COM_AllocString(title)] = newCode;
		return newCode;
	}
	return it->second;
}

int cConfig :: LookupFFCode( const char *const title )
{
	FFCodeMap::const_iterator it = m_ffCodeInfo.find( title );
	if ( it == m_ffCodeInfo.end() )
		return -1;
	return it->second;
}

int cConfig :: LookupSFFCode( const char *const title )
{
	FFCodeMap::const_iterator it = m_sffCodeInfo.find( title );
	if ( it == m_sffCodeInfo.end() )
		return -1;
	return it->second;
}

const char *cConfig :: LookupFFTitle( const int code )
{
	if ( code == -1 )
		return "X";

	FFCodeMap::const_iterator it = std::find_if( m_ffCodeInfo.begin(), m_ffCodeInfo.end(), FindPairByValueFunctor<const char*,int>( code ) );
	if ( it == m_ffCodeInfo.end() )
		return "???";
	return it->first;
}

const cAtomInfo *cConfig :: LookupAtomInfo( const char *const symbol )
{
	AtomInfoMap::const_iterator it = m_atomInfo.find( symbol );
	if ( it == m_atomInfo.end() )
		return NULL;
	return &it->second;
}

const cBondInfo *cConfig :: LookupBondInfo( int ffCode1, int ffCode2 )
{
	std::pair<int,int> pair;
	if ( ffCode1 <= ffCode2 )
		pair = std::make_pair( ffCode1, ffCode2 );
	else
		pair = std::make_pair( ffCode2, ffCode1 );

	BondInfoMap::const_iterator it = m_bondInfo.find( pair );
	if ( it == m_bondInfo.end() )
		return NULL;
	return &it->second;
}

const cAngleInfo *cConfig :: LookupAngleInfo( int ffCode1, int ffCode2, int ffCode3 )
{
	std::pair<int,std::pair<int,int> > pair;

	if ( ffCode1 <= ffCode3 )
		pair = std::make_pair( ffCode2, std::make_pair( ffCode1, ffCode3 ) );
	else
		pair = std::make_pair( ffCode2, std::make_pair( ffCode3, ffCode1 ) );

	AngleInfoMap::const_iterator it = m_angleInfo.find( pair );
	if ( it == m_angleInfo.end() )
		return NULL;
	return &it->second;
}

const cHBInfo *cConfig :: LookupHBInfo( int ffCodeH, int ffCodeA )
{
	HBInfoMap::const_iterator it = m_hbInfo.find( std::make_pair( ffCodeH, ffCodeA ) );
	if ( it == m_hbInfo.end() )
		return NULL;
	return &it->second;
}

bool cConfig :: IsHBHydrogen( int ffCode )
{
	return std::binary_search( m_hbHydrogens.begin(), m_hbHydrogens.end(), ffCode );
}

bool cConfig :: IsHBAcceptor( int ffCode )
{
	return std::binary_search( m_hbAcceptors.begin(), m_hbAcceptors.end(), ffCode );
}

const cImproperInfo *cConfig :: LookupImproperInfo( int ffCodeI, int ffCodeJ, int ffCodeK, int ffCodeL )
{
	std::pair<std::pair<int,int>,std::pair<int,int> > pair = std::make_pair( std::make_pair( ffCodeI, ffCodeJ ), std::make_pair( ffCodeK, ffCodeL ) );
	ImproperInfoMap::const_iterator it = m_improperInfo.find( pair );
	if ( it == m_improperInfo.end() )
		return NULL;
	return &it->second;
}

const cTorsionInfo *cConfig :: LookupTorsionInfo( int ffCodeI, int ffCodeJ, int ffCodeK, int ffCodeL )
{
	TorsionInfoList::const_iterator it = std::find_if( m_torsionInfo.begin(), m_torsionInfo.end(), TorsionInfoSearchFunctor( ffCodeI, ffCodeJ, ffCodeK, ffCodeL ) );

	if ( it == m_torsionInfo.end() )
		return NULL;
	return &(*it);
}

const cVdWInfo *cConfig :: LookupVdWInfo( int ffCode )
{
	VdWInfoMap::const_iterator it = m_vdwInfo.find( ffCode );
	if ( it == m_vdwInfo.end() )
		return NULL;
	return &it->second;
}

const cSolvGSInfo *cConfig :: LookupSolvGSInfo( int sffCode )
{
	if ( sffCode < 0 )
		return NULL;
	SolvGSInfoMap::const_iterator it = m_solvGSInfo.find( sffCode );
	if ( it == m_solvGSInfo.end() )
		return NULL;
	return &it->second;
}

bool cConfig :: CheckResidueInfo( const char *const name )
{
	// Find the residue
	ResidueInfoMap::const_iterator it = m_residueInfo.find( name );
	if ( it == m_residueInfo.end() ) {
		// Check aliases
		ResidueAliasMap::const_iterator ita = m_aliasInfo.find( name );
		if ( ita == m_aliasInfo.end() )
			return false;
		it = m_residueInfo.find( ita->second );
		if ( it == m_residueInfo.end() )
			return false;
	}
	return true;
}

const cResidueLocation *cConfig :: LookupResidueInfo( const char *const name, eResidueLocation loc )
{
	const cResidueLocation *pLocation = NULL;

	if ( loc < 0 || loc >= RL_MAX )
		Log().Fatal( "LookupResidueInfo: bad location %i\n", loc );

	// Find the residue
	ResidueInfoMap::const_iterator it = m_residueInfo.find( name );
	if ( it == m_residueInfo.end() ) {
		// Check aliases
		ResidueAliasMap::const_iterator ita = m_aliasInfo.find( name );
		if ( ita == m_aliasInfo.end() )
			return NULL;
		it = m_residueInfo.find( ita->second );
		if ( it == m_residueInfo.end() )
			return NULL;
	}

	pLocation = &it->second.mLocation[loc];

	// Check if the location was defined
	if ( !pLocation || pLocation->mNumAtoms <= 0)
		return NULL;

	return pLocation;
}
