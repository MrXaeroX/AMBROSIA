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
#ifndef MD_CONFIG_H
#define MD_CONFIG_H
/**
* @file
* @brief	Declaration of a config class.
* @details	The file defines the config class.
*/

/** @brief	Default short cut-off distance (R1), in angstroms. */
#define DEFAULT_SHORT_CUTOFF_DISTANCE		9
/** @brief	Default short switch distance, in angstroms. */
#define DEFAULT_SHORT_SWITCH_DISTANCE		7
/** @brief	Default long cut-off distance (R2), in angstroms. */
#define DEFAULT_LONG_CUTOFF_DISTANCE		15
/** @brief	Default maximum length of a recognizable S-S bond. */
#define DEFAULT_SS_BOND_DISTANCE			vec_t( 4.2 )
/** @brief	Default dielectric constant for a solvent (assume pure water). */
#define DEFAULT_SOLVENT_DIELECTRIC_CONST	80
/** @brief	Default dielectric constant for a solvent (assume pure water). */
#define DEFAULT_RDIE_CONST					vec_t( 1.0 )
/** @brief	Default debye radius for a solvent, in angstroms (assume pure water). */
#define DEFAULT_SOLVENT_DEBYE_RADIUS		10000
/** @brief	Default hydrogen bond model. */
#define DEFAULT_HBOND_MODEL					128
/** @brief	Default hydrogen bond energy scale. */
#define DEFAULT_HBOND_SCALE					1
/** @brief	Default hydrogen bond cut-off distance, in angstroms. */
#define DEFAULT_HBOND_CUTOFF_DISTANCE		vec_t( 3.5 )
/** @brief	Minimum distant restraint length, in angstroms. */
#define DEFAULT_MIN_DIST_RESTR_LENGTH		vec_t( 0.1 )

/**
* @brief	Enumeration of possible residue locations.
*/
typedef enum {
	RL_INT = 0,
	RL_BEG,
	RL_END,
	RL_ISO,
	RL_MAX
} eResidueLocation;

/**
* @brief	Enumeration of solvation models.
*/
typedef enum {
	SOLV_NONE = 0,
	SOLV_GAUSS,
	SOLV_GBORN
} eSolvModel;

/**
* @brief	Description of an atom.
*/
typedef struct stAtomInfo {
	/** @brief	Resiprocal mass (in a.e.m.). */
	vec_t	mMassReciprocal;
	/** @brief	Atomic radius (in angstroms). */
	vec_t	mRadius;
	/** @brief	Electronegativity (qEq). */
	vec_t	mXi;
	/** @brief	Self-repulsion (qEq). */
	vec_t	mJ0;
} cAtomInfo;

/** @brief	This atom is a backbone atom at the beginning of the topology. */
#define TF_TOPOLOGY_BEGIN	( 1 << 0 )
/** @brief	This atom is a backbone atom at the end of the topology. */
#define TF_TOPOLOGY_END		( 1 << 1 )
/** @brief	This atom is a part of planar ring structure. */
#define TF_TOPOLOGY_RING	( 1 << 2 )

/**
* @brief	Description of a single atom of a residue.
*/
typedef struct stResidueAtom {
	/** @brief	Index of the atom in the Z-matrix. */
	int		mIndex;
	/** @brief	Previous atoms in the Z-matrix */
	int		mPrevious[3];
	/** @brief	Z-matrix R-coordinate */
	vec_t	mR;
	/** @brief	Z-matrix Theta-coordinate */
	vec_t	mTheta;
	/** @brief	Z-matrix Phi-coordinate */
	vec_t	mPhi;
	/** @brief	Partial charge on the atom */
	vec_t	mCharge;
	/** @brief	PDB title of the atom */
	char	mTitle[8];
	/** @brief	ForceField title of the atom */
	char	mFFTitle[8];
	/** @brief	Solvent ForceField title of the atom */
	char	mSFFTitle[8];
	/** @brief	Symbol of the atom */
	char	mSymbol[4];
	/** @brief	Type of the atom (B, S, H, etc.) */
	int		mType;
	/** @brief	Topology flags */
	int		mFlags;
	/** @brief	Loop flags (bits set are loop indices) */
	int		mLoopFlags;
} cResidueAtom;

/** @brief	Type for an array of #cResidueAtom objects. */
typedef std::vector<cResidueAtom> ResidueAtomArray;

/** @brief Residue atom search by index functor. */
class ResidueAtomSearchByIndexFunctor {
public:
	/** 
	* @brief Constructor.
	* @param index : index of the atom in the Z-matrix.
	*/
	ResidueAtomSearchByIndexFunctor( int index ) : m_index( index ) { assert( index > 0 ); }
	/** @brief Compare operator. */
	bool operator() ( const cResidueAtom &at ) { return ( at.mIndex == m_index ); }
private:
	int m_index;
};

/** @brief Residue atom search by title functor. */
class ResidueAtomSearchByTitleFunctor {
public:
	/** 
	* @brief Constructor.
	* @param title : title of the atom in the Z-matrix.
	*/
	ResidueAtomSearchByTitleFunctor( const char *title ) : m_title( title ) { assert( title != NULL ); }
	/** @brief Compare operator. */
	bool operator() ( const cResidueAtom &at ) { return !strcmp( at.mTitle, m_title ); }
private:
	const char *m_title;
};

/**
* @brief	Description of a residue in a certain location.
*/
typedef struct stResidueLocation {
	/** @brief	Number of atoms in the residue. */
	int		mNumAtoms;
	/** @brief	Number of loops in the residue. */
	int		mNumLoops;
	/** @brief	Total charge of the residue. */
	vec_t	mTotalCharge;
	/** @brief	Array of definitions of atoms, containing #cResidueAtom records. */
	ResidueAtomArray mAtoms;
} cResidueLocation;

/**
* @brief	Description of a residue.
*/
typedef struct stResidueInfo {
	/** @brief	Array of all possible locations, containing #cResidueLocation records. */
	cResidueLocation mLocation[RL_MAX];
} cResidueInfo;

/**
* @brief	Description of a covalent bond.
*/
typedef struct stBondInfo {
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*angstrom^2)@f$). */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium bond length (in angstroms). */
	vec_t	mLength;
} cBondInfo;

/**
* @brief	Description of a covalent angle.
*/
typedef struct stAngleInfo {
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*rad^2)@f$). */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium bond angle (in radians). */
	vec_t	mTheta;
} cAngleInfo;

/**
* @brief	Description of an improper torsion angle.
*/
typedef struct stImproperInfo {
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*rad^2)@f$). */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium improper angle (in radians). */
	vec_t	mXi;
} cImproperInfo;

/**
* @brief	Description of a single harmonic of a torsion angle.
*/
typedef struct stTorsionHarm {
	/** @brief	Barrier force constant (in kcal/mol). */
	vec_t	mHarmonicConstant;
	/** @brief	Barrier scale (for non-ring structures). */
	vec_t	mHarmonicScale;
	/** @brief	Cosine of phase shift. */
	vec_t	mPhaseCos;
	/** @brief	Sine of phase shift. */
	vec_t	mPhaseSin;
	/** @brief	Periodicity of a barrier. */
	int		mPeriodicity;
	/** @brief	Pointer to the next harmonic. */
	struct stTorsionHarm *mpNext;
} cTorsionHarm;

/**
* @brief	Description of a proper torsion angle.
*/
typedef struct stTorsionInfo {
	/** @brief	Force field codes of atoms I-J-K-L (-1 = any). */
	int		mFFCode[4];
	/** @brief	Fourier series representing the torsion. */
	cTorsionHarm *mpHarmonic;
} cTorsionInfo;


/** @brief Torsion info search by atomic indices. */
class TorsionInfoSearchFunctor {
public:
	/** 
	* @brief Constructor.
	* @param i : index of atom I.
	* @param j : index of atom J.
	* @param k : index of atom K.
	* @param l : index of atom L.
	*/
	TorsionInfoSearchFunctor( int i, int j, int k, int l ) : m_i( i ), m_j( j ), m_k( k ), m_l( l ) {}
	/** @brief	Compare operator. */
	bool operator() ( const std::list<cTorsionInfo>::const_reference check ) { 
		return ( ( compare( check.mFFCode[0], m_i ) && compare( check.mFFCode[1], m_j ) && compare( check.mFFCode[2], m_k ) && compare( check.mFFCode[3], m_l ) ) ||
					 ( compare( check.mFFCode[3], m_i ) && compare( check.mFFCode[2], m_j ) && compare( check.mFFCode[1], m_k ) && compare( check.mFFCode[0], m_l ) ) );
	}
private:
	bool compare( int isrc, int iref ) { return ( isrc < 0 || isrc == iref ); }
private:
	int m_i, m_j, m_k, m_l;
};

/**
* @brief	Description of van der Waals parameters of an atom.
*/
typedef struct stVdWInfo {
	/** @brief	Rmin/2, in angstroms. */
	vec_t	mR;
	/** @brief	Well depth, in kcal/mol. */
	vec_t	mE;
} cVdWInfo;

/**
* @brief	Description of H-bond parameters of an atom.
*/
typedef struct stHBInfo {
	/** @brief	Rmin/2, in angstroms. */
	vec_t	mR;
	/** @brief	Well depth, in kcal/mol. */
	vec_t	mE;
} cHBInfo;

/**
* @brief	Description of Gaussian solvation parameters of an atom.
*/
typedef struct stSolvGSInfo {
	/** @brief	Atomic volume, in cubic angstroms. */
	vec_t	mVolume;
	/** @brief	Reference deltaG, in kcal/mol. */
	vec_t	mGref;
	/** @brief	Free deltaG, in kcal/mol. */
	vec_t	mGfree;
	/** @brief	Correlation length, in angstroms. */
	vec_t	mLambda;
} cSolvGSInfo;

/**
* @brief	Description of position restraint.
*/
typedef struct stPosRestraintInfo {
	/** @brief	Index of the chain. */
	int		mChain;
	/** @brief	Index of the first residue. */
	int		mFirstResidue;
	/** @brief	Index of the last residue. */
	int		mLastResidue;
	/** @brief	Harmonic constant (in @f$ kcal/(mol*angstrom^2)@f$). */
	vec_t	mHarmConst;
	/** @brief	True if restraint should be applied to backbone only. */
	bool	mBackboneOnly;
} cPosRestraintInfo;

/**
* @brief	Description of distant restraint.
*/
typedef struct stDistRestraintInfo {
	/** @brief	Index of the chain. */
	int		mChain[2];
	/** @brief	Index of the residue. */
	int		mResidue[2];
	/** @brief	Title of the atom. */
	char	mAtomTitle[2][8];
	/** @brief	Harmonic constant (in @f$ kcal/(mol*angstrom^2)@f$). */
	vec_t	mHarmConst;
	/** @brief	True if restraint should be applied to backbone only. */
	vec_t	mDistance;
} cDistRestraintInfo;

/**
* @brief	Simulation parameters that are set up via parameters' file.
*/
typedef struct stParameters {
	/** @brief	Defines which model to read (from files that support multiple models, e.g. PDB). */
	int		mLoadModel;
	/** @brief	Defines whether to read H-atoms from the source structure, or restore them automatically. */
	int		mReadHydrogens;
	/** @brief	Automatic detection of S-S bonds. */
	int		mAutoSSBonds;
	/** @brief	Maximum length of a recognizable S-S bond. */
	vec_t	mSSBondDist;
	/** @brief	Automatic (re)calculation of partial charges. */
	int		mAutoQ;
	/** @brief	Automatic rebase of residue serial numbers. */
	int		mAutoRSN;
	/** @brief	Calculate energy for the source structure. */
	int		mEnergyCalc;
	/** @brief	Optimize energy for the source structure. */
	int		mEnergyOptimize;
	/** @brief	Perform N steps of optimization procedure. */
	int		mEnergyOptimizeSteps;
	/** @brief	Van der Waals cut-off distance, in angstroms. */
	vec_t	mShortCutoffDist;
	/** @brief	Van der Waals switch distance, in angstroms. */
	vec_t	mShortSwitchDist;
	/** @brief	Electrostatic cut-off distance, in angstroms. */
	vec_t	mLongCutoffDist;
	/** @brief	Scale of distance-dependent dielectric constant (0 = no distance dependence). */
	vec_t	mDistanceDependentDielectricConstant;
	/** @brief	Dielectric constant of the solvent. */
	vec_t	mSolventDielectricConstant;
	/** @brief	Debye radius of the solvent, in angstroms. */
	vec_t	mDebyeRadius;
	/** @brief	Hydrogen bonding model (either 126 or 128). */
	int		mHBondModel;
	/** @brief	Hydrogen bonding energy scale. */
	vec_t	mHBondScale;
	/** @brief	Hydrogen bonding cut-off distance, in angstroms. */
	vec_t	mHBondCutoffDist;
	/** @brief	Solvation model: NULL (none), GS (gaussian), GB (generalized Born). */
	char	*mSolvationModel;
} cParameters;

/**
* @brief	Enumeration of field types used in description of parameters' fields.
*/
typedef enum {
	FIELD_INTEGER = 0,
	FIELD_FLOAT,
	FIELD_STRING,
	FIELD_FLAG
} eFieldType;

/**
* @brief	Description of parametric fields from #cParameters, that are parsed.
*/
typedef struct stParametersDesc {
	/** @brief	Parameter name. */
	const char	*mParmName;
	/** @brief	Field name. */
	const char	*mFieldName;
	/** @brief	Field offset in the structure. */
	intptr_t	mFieldOffset;
	/** @brief	Field type. */
	eFieldType	mFieldType;
	/** @brief	Use default. */
	bool		mDefault;
} TYPEDESCRIPTION;

/**
* @brief	Helper macro to define fields in a #TYPEDESCRIPTION.
*/
#define DEFINE_PARAM_FIELD(x,y,z)	{ x, #y, offsetof(cParameters,y), z, true }

/**
* @brief	This object incapsulates persistant molecular dynamics configuration.
* @details	The following parameters are defined in config:\n
* -# atomic properties (mass and radius);
* -# covalent bond properties.
*/
class cConfig
{
	DECLARE_SINGLETON( cConfig );
	~cConfig();

public:
	/** @brief	Type for a map of force field atom name and code. */
	typedef std::map<const char*,int,StringCompareFunctor> FFCodeMap;
	/** @brief	Type for a map of #cAtomInfo objects. */
	typedef std::map<const char*,cAtomInfo,StringCompareFunctor> AtomInfoMap;
	/** @brief	Type for a map of #cResidueInfo objects. */
	typedef std::map<const char*,cResidueInfo,StringCompareFunctor> ResidueInfoMap;
	/** @brief	Type for a map of #cResidueInfo objects. */
	typedef std::map<const char*,const char*,StringCompareFunctor> ResidueAliasMap;
	/** @brief	Type for a map of #cBondInfo objects. */
	typedef std::map<std::pair<int, int>,cBondInfo> BondInfoMap;
	/** @brief	Type for a map of #cAngleInfo objects. */
	typedef std::map<std::pair<int,std::pair<int, int> >,cAngleInfo> AngleInfoMap;
	/** @brief	Type for a map of #cImproperInfo objects. */
	typedef std::map<std::pair<std::pair<int,int>,std::pair<int,int> >,cImproperInfo> ImproperInfoMap;
	/** @brief	Type for a list of #cTorsionInfo objects (list is used instead of map, because we support "any" type of atom). */
	typedef std::list<cTorsionInfo> TorsionInfoList;
	/** @brief	Type for a map of #cVdWInfo objects. */
	typedef std::map<int,cVdWInfo> VdWInfoMap;
	/** @brief	Type for a map of #cHBInfo objects. */
	typedef std::map<std::pair<int, int>,cHBInfo> HBInfoMap;
	/** @brief	Type for a list of hydrogen/acceptor codes. */
	typedef std::vector<int> HBCodeVector;
	/** @brief	Type for a map of #cSolvGSInfo objects. */
	typedef std::map<int,cSolvGSInfo> SolvGSInfoMap;
	/** @brief	Type for a list of positional restraints. */
	typedef std::vector<cPosRestraintInfo> PosRestrVector;
	/** @brief	Type for a list of distant restraints. */
	typedef std::vector<cDistRestraintInfo> DistRestrVector;

public:
	/**
	* @brief	Sets name of a user-defined configuration file.
	* @details	Default is empty (no user-defined file).
	* @param name : new configuration file name.
	*/
	static void SetUserConfigName( const char *name );
	/**
	* @brief	Returns a name of a user-defined configuration file.
	* @return The name of a user-defined configuration file.
	*/
	static const char *UserConfigName( void ) { return m_userConfigName; }
	/**
	* @brief	Sets name of a file with simulation parameters.
	* @details	Default is #DEFAULT_PARMS_FILENAME.
	* @param name : new file name.
	*/
	static void SetParmsFileName( const char *name );
	/**
	* @brief	Returns a name of a file with simulation parameters.
	* @return The name of a file with simulation parameters.
	*/
	static const char *ParmsFileName( void ) { return m_parmsName; }
	/**
	* @brief	Returns a string for a given location index (useful for logging).
	* @param loc : location index.
	* @return The name of the location.
	*/
	static const char *LocationToString( eResidueLocation loc );
	/**
	* @brief	Loads and parses configuration file.
	* @details	If any config data already exists, new data is appended, overriding old data when necessary.
	* @param name : config file name.
	*/
	void LoadConfig( const char *const name );
	/**
	* @brief	Loads and parses user-defined configuration file.
	* @details	File name is defined with #SetUserConfigName. If the name is empty, nothing is loaded.
	*/
	void LoadUserConfig( void );
	/**
	* @brief	Loads and parses file with simulation parameters.
	* @details	File name is defined with #SetParmsFileName.
	*/
	void LoadParams( void );
	/**
	* @brief	Clears the config, frees all data allocated.
	*/
	void Clear( void );
	/**
	* @brief	Get force field code for a given force field title.
	* @param title : force field title of the atom ("HW", "OW", "N*", etc.), must be uppercase.
	* @return Force field code identifier.
	*/
	int LookupFFCode( const char *const title );
	/**
	* @brief	Get force field title for a given force field code.
	* @param code : force field code identifier.
	* @return Force field title of the atom ("HW", "OW", "N*", etc.), uppercase.
	*/
	const char *LookupFFTitle( const int code );
	/**
	* @brief	Get solvent force field code for a given solvent force field title.
	* @param title : solvent force field title of the atom, must be uppercase.
	* @return Solvent force field code identifier.
	*/
	int LookupSFFCode( const char *const title );
	/**
	* @brief	Get solvent force field title for a given solvent force field code.
	* @param code : solvent force field code identifier.
	* @return Solvent force field title of the atom, uppercase.
	*/
	const char *LookupSFFTitle( const int code );
	/**
	* @brief	Get description of an atom for the atomic symbol.
	* @param symbol : atomic symbol ("H", "C", "N", etc.), must be uppercase.
	* @return Constant pointer to structure describing the atom.
	*/
	const cAtomInfo *LookupAtomInfo( const char *const symbol );
	/**
	* @brief	Get description of a bond for given atoms.
	* @param ffCode1 : force field code for the first atom.
	* @param ffCode2 : force field code for the second atom.
	* @return Constant pointer to structure describing the bond.
	*/
	const cBondInfo *LookupBondInfo( int ffCode1, int ffCode2 );
	/**
	* @brief	Get description of an angle for given atoms.
	* @param ffCode1 : force field code for the first atom.
	* @param ffCode2 : force field code for the second atom.
	* @param ffCode3 : force field code for the third atom.
	* @return Constant pointer to structure describing the angle.
	*/
	const cAngleInfo *LookupAngleInfo( int ffCode1, int ffCode2, int ffCode3 );
	/**
	* @brief	Get description of an improper torsion angle for given atoms.
	* @param ffCodeI : force field code for the atom I.
	* @param ffCodeJ : force field code for the atom J (center).
	* @param ffCodeK : force field code for the atom K.
	* @param ffCodeL : force field code for the atom L (terminator).
	* @return Constant pointer to structure describing the improper.
	*/
	const cImproperInfo *LookupImproperInfo( int ffCodeI, int ffCodeJ, int ffCodeK, int ffCodeL );
	/**
	* @brief	Get description of a proper torsion angle for given atoms.
	* @param ffCodeI : force field code for the atom I.
	* @param ffCodeJ : force field code for the atom J.
	* @param ffCodeK : force field code for the atom K.
	* @param ffCodeL : force field code for the atom L.
	* @return Constant pointer to structure describing the torsion.
	*/
	const cTorsionInfo *LookupTorsionInfo( int ffCodeI, int ffCodeJ, int ffCodeK, int ffCodeL );
	/**
	* @brief	Get description of van der Waals parameters for a given atom.
	* @param ffCode : force field code for the atom.
	* @return Constant pointer to structure describing van der Waals parameters.
	*/
	const cVdWInfo *LookupVdWInfo( int ffCode );
	/**
	* @brief	Get description of hydrogen bond parameters for a given atom pair.
	* @param ffCodeH : force field code for the hydrogen atom.
	* @param ffCodeA : force field code for the acceptor atom.
	* @return Constant pointer to structure describing hydrogen bond parameters.
	*/
	const cHBInfo *LookupHBInfo( int ffCodeH, int ffCodeA );
	/**
	* @brief	Return whether an atom is a hydrogen that can form a hydrogen bond.
	* @param ffCode : force field code for the atom.
	* @return True or false.
	*/
	bool IsHBHydrogen( int ffCode );
	/**
	* @brief	Return whether an atom is an acceptor that can form a hydrogen bond.
	* @param ffCode : force field code for the atom.
	* @return True or false.
	*/
	bool IsHBAcceptor( int ffCode );
	/**
	* @brief	Return list of known hydrogen-bonding hydrogen FF codes.
	* @return #HBCodeVector.
	*/
	const HBCodeVector &GetHBHydrogens( void ) { return m_hbHydrogens; }
	/**
	* @brief	Return list of known hydrogen-bonding acceptor FF codes.
	* @return #HBCodeVector.
	*/
	const HBCodeVector &GetHBAcceptors( void ) { return m_hbAcceptors; }
	/**
	* @brief	Get description of Gaussian solvation parameters for a given atom pair.
	* @param sffCode : force field code for the atom.
	* @return Constant pointer to structure describing solvation parameters.
	*/
	const cSolvGSInfo *LookupSolvGSInfo( int sffCode );
	/**
	* @brief	Get pointer to a structure describing the residue in a given location.
	* @param name : name of the residue ("GLY", "ALA", "G", etc.), must be uppercase.
	* @param loc : location of the residue.
	* @return Constant pointer to a structure describing the residue in a given location.
	*/
	const cResidueLocation *LookupResidueInfo( const char *const name, eResidueLocation loc );
	/**
	* @brief	Check whether a named residue is defined in the topology library.
	* @param name : name of the residue ("GLY", "ALA", "G", etc.), must be uppercase.
	* @return True (exists) or false (doesn't exist).
	*/
	bool CheckResidueInfo( const char *const name );
	/**
	* @brief	Helper function to access simulation parameters.
	* @return Constant reference to structure describing simulation parameters.
	*/
	const cParameters &Parameters( void ) { return m_parms; }
	/**
	* @brief	Helper function to access positional restraints.
	* @return Constant reference to vector of positional restraints.
	*/
	const PosRestrVector &PositionalRestraints( void ) { return m_posRestraints; }
	/**
	* @brief	Helper function to access distant restraints.
	* @return Constant reference to vector of distant restraints.
	*/
	const DistRestrVector &DistantRestraints( void ) { return m_distRestraints; }

private:
	void	InitDefaults( void );
	void	ValidateParams( void );
	int		FetchFFCode( const char *const title );
	int		FetchSFFCode( const char *const title );
	int		ChainIndexFromToken( const char *token );
	void	ParseConfig( byte *fileData );
	byte*	ParseAtom( byte *fileData );
	byte*	ParseQEq( byte *fileData );
	byte*	ParseResidue( byte *fileData );
	byte*	ParseResidueLocation( cResidueLocation *pLoc, byte *fileData );
	byte*	ParseBond( byte *fileData );
	byte*	ParseAngle( byte *fileData );
	byte*	ParseImproper( byte *fileData );
	byte*	ParseTorsion( byte *fileData );
	byte*	ParseVdWParms( byte *fileData );
	byte*	ParseHBParms( byte *fileData );
	byte*	ParseSolvationParms( byte *fileData );
	void	ParseParams( byte *fileData );
	void	FreeParams( void );
	void	PrintParams( void );
	void	DumpResidue( const char *title, const cResidueInfo *pInfo );

private:
	static char		m_userConfigName[MAX_OSPATH];
	static char		m_parmsName[MAX_OSPATH];
	static TYPEDESCRIPTION m_parmData[];

private:
	int				m_lineNumber;
	cParameters		m_parms;
	PosRestrVector	m_posRestraints;
	DistRestrVector	m_distRestraints;
	FFCodeMap		m_ffCodeInfo;
	FFCodeMap		m_sffCodeInfo;
	AtomInfoMap		m_atomInfo;
	ResidueInfoMap	m_residueInfo;
	ResidueAliasMap	m_aliasInfo;
	BondInfoMap		m_bondInfo;
	AngleInfoMap	m_angleInfo;
	ImproperInfoMap m_improperInfo;
	TorsionInfoList	m_torsionInfo;
	VdWInfoMap		m_vdwInfo;
	HBInfoMap		m_hbInfo;
	HBCodeVector	m_hbHydrogens;
	HBCodeVector	m_hbAcceptors;
	SolvGSInfoMap	m_solvGSInfo;
	char			m_token[1024];
};

/**
* @brief	Helper function to get global config singleton.
* @return Reference to the global config object.
*/
inline cConfig& Config( void )
{
	return cConfig::Instance();
}

#endif /*MD_CONFIG_H*/
