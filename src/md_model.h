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
#ifndef MD_MODEL_H
#define MD_MODEL_H
/**
* @file
* @brief	Declaration of a model class.
* @details	The file defines the model class.
*/

/** @brief	Maximum depth of nested functions to profile. */
#define MAX_PROFILER_DEPTH	8

/**
* @brief	Model header (global data).
* @details	This is a group of persistant data that can be saved or restored.
*/
typedef struct stModelHeader {
	/** @brief	Number of atoms in a model. */
	int			mNumAtoms;
	/** @brief	Number of heavy atoms (i.e. not hydrogens) in a model. */
	int			mNumHeavyAtoms;
	/** @brief	Number of residues in a model. */
	int			mNumResidues;
	/** @brief	Number of chains in a model. */
	int			mNumChains;
	/** @brief	Header comment. */
	char		mComment[8192];
} cModelHeader;

/**
* @brief	Connectivity information for an atom.
*/
typedef struct stConnectivity {
	/** @brief  List of first and second neighbours (1-2 and 1-..-3). */
	IntArray mNeighbours_12_13;
	/** @brief  List of third neighbours (1-..-..-4). */
	IntArray mNeighbours_14;
	/** @brief  List of first neighbour hydrogens. */
	IntArray mNeighbours_12H;
} cConnectivity;

/**
* @brief	Van der Waals pair parameters.
*/
typedef struct stVdWPair {
	/** @brief  Value of C12 for the pair. */
	vec_t	mC12;
	/** @brief  Value of C6 for the pair. */
	vec_t	mC6;
} cVdWPair;

/**
* @brief	Hydrogen bonding pair parameters.
*/
typedef struct stHBPair {
	/** @brief  Value of Rmin for the pair. */
	vec_t	mR;
	/** @brief  Value of E for the pair. */
	vec_t	mE;
	/** @brief  Value of A (C12) for the pair. */
	vec_t	mA;
	/** @brief  Value of B (C6 or C8) for the pair. */
	vec_t	mB;
} cHBPair;

/**
* @brief	Solvation parameters.
*/
typedef struct stSolvParms {
	union {
		struct {
			/** @brief	Van der Waals radius, in angstroms. */
			vec_t	mR;
			/** @brief	Volume, in cubic angstroms. */
			vec_t	mV;
			/** @brief	A = 2*dGfree/[4*pi*sqrt(pi)*lambda]. */
			vec_t	mA;
			/** @brief	B = -1/lambda^2. */
			vec_t	mB;
		} GS;
	} 
	/** @brief	Union for different solvent model parameters. */
	u;
} cSolvParms;

/**
* @brief	QEq parameters.
*/
typedef struct stQEqParms {
	/** @brief	Electronegativity. */
	vec_t	mXi;
	/** @brief	Self-repulsion. */
	vec_t	mJ0;
	/** @brief	Total charge of the residue. */
	vec_t	mQTot;
} cQEqParms;

/** @brief	Atom is heavy (i.e. not hydrogen). */
#define AF_HEAVY			( 1 << 0 )
/** @brief	This atom is a potential hydrogen bond donor. */
#define AF_HB_DONOR			( 1 << 1 )
/** @brief	This atom is a potential hydrogen bond acceptor. */
#define AF_HB_ACCEPTOR		( 1 << 2 )
/** @brief	Residue is at the beginning of a chain. */
#define AF_RL_BEGIN			( 1 << 3 )
/** @brief	Residue is at the end of a chain. */
#define AF_RL_END			( 1 << 4 )
/** @brief	Atom's position is harmonically restrained. */
#define AF_RESTRAINED		( 1 << 5 )

/**
* @brief	Single atom of a model.
*/
typedef struct stAtom {
	/** @brief	Pointer to physics properties (NULL until topology is complete). */
	struct stPhysAtom *mpPhysAtom;
	/** @brief	Pointer to connectivity information (NULL until topology is complete). */
	cConnectivity *mpConnectivity;
	/** @brief	Pointer to solvation parameters (NULL until topology is complete and some solvation model is enabled). */
	cSolvParms	*mpSolvation;
	/** @brief	Pointer to qEq parameters (NULL until topology is complete and autoQ is enabled). */
	cQEqParms	*mpQEq;
	/** @brief	Atomic flags (see AF_xxx defines). */
	int		mFlags;
	/** @brief	Pointer to physics position (NULL until topology is complete). */
	vec_t	*mCurrentPosition;
	/** @brief	Pointer to physics force (NULL until topology is complete). */
	vec_t	*mCurrentForce;
	/** @brief	Original position of the atom (in angstroms). */
	vec3_t	mOriginalPosition;
	/** @brief	Harmonic constant for position restraining (in @f$ kcal/(mol*angstrom^2)@f$). */
	vec_t	mPosRestraintHarmConst;
	/** @brief	Force field code for the atom. */
	int		mFFCode;
	/** @brief	Van der Waals code for the atom (for pair table lookup). */
	int		mVdWCode;
	/** @brief	Hydrogen bond code for the atom (for pair table lookup). */
	int		mHBCode;
	/** @brief	Sort number of the atom (used by sorting procedures). */
	int		mSortNumber;
	/** @brief	Serial number of the residue (one-based index). */
	int		mResidueNumber;
	/** @brief	Serial number of the chain (one-based index). */
	int		mChainNumber;
	/** @brief	Sequence number of the residue (matches PDB data). */
	int		mResidueSequenceNumber;
	/** @brief	Pointer to residue atom description. */
	const cResidueAtom *mpResidueAtom;
	/** @brief	Atom symbol. */
	char	mAtomSymbol[4];
	/** @brief	Atom title. */
	char	mAtomTitle[8];
	/** @brief	Residue title. */
	char	mResidueTitle[8];
} cAtom;

/**
* @brief	Physical properties of an atom.
*/
typedef struct stPhysAtom {
	/** @brief	Velocity of the atom in homogeneous space (in angstroms per picosecond). */
	vec4_t	mVelocity;
	/** @brief	Acceleration of the atom in homogeneous space (in angstroms per sq. picosecond). */
	vec4_t	mAccel;
	/** @brief	Resiprocal mass (in mol/g), derived from #cAtomInfo. */
	vec_t	mMassReciprocal;
	/** @brief	Atomic radius (in angstroms), derived from #cAtomInfo. */
	vec_t	mRadius;
	/** @brief	Align to 16-bytes. */
	vec_t	mUnused[2];
} cPhysAtom;

/**
* @brief	Single covalent bond of a model.
*/
typedef struct stBond {
	/** @brief	Indices of atoms comprising the bond. */
	int		mAtomIndices[2];
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*angstrom^2)@f$), derived from #cBondInfo. */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium bond length (in angstroms), derived from #cBondInfo. */
	vec_t	mLength;
} cBond;

/**
* @brief	Single covalent disulfide (S-S) bond of a model.
*/
typedef struct stSSBond {
	/** @brief	Indices of atoms comprising the bond. */
	int		mAtomIndices[2];
	/** @brief	Residue numbers. */
	int		mResidueNumbers[2];
	/** @brief	Titles of atoms. */
	char	mAtomTitles[2][8];
} cSSBond;

/**
* @brief	Single distant restraint of a model.
*/
typedef struct stDistRestrain {
	/** @brief	Indices of atoms comprising the restraint. */
	int		mAtomIndices[2];
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*angstrom^2)@f$). */
	vec_t	mHarmonicConstant;
	/** @brief	Square of equilibrium interatomic distance (in square angstroms). */
	vec_t	mLengthSquared;
} cDistRestrain;

/**
* @brief	Single covalent angle of a model.
*/
typedef struct stAngle {
	/** @brief	Indices of atoms comprising the angle. */
	int		mAtomIndices[3];
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*rad^2)@f$), derived from #cAngleInfo. */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium angle value (in radians), derived from #cAngleInfo. */
	vec_t	mTheta;
} cAngle;

/**
* @brief	Single improper torsion angle of a model.
*/
typedef struct stImproper {
	/** @brief	Indices of atoms comprising the improper (I-J-K-L). */
	int		mAtomIndices[4];
	/** @brief	Harmonic force constant (in @f$ kcal/(mol*rad^2)@f$), derived from #cImproperInfo. */
	vec_t	mHarmonicConstant;
	/** @brief	Equilibrium angle value (in radians), derived from #cImproperInfo. */
	vec_t	mXi;
} cImproper;

/**
* @brief	Single proper torsion angle of a model.
*/
typedef struct stTorsion {
	/** @brief	Indices of atoms comprising the torsion (I-J-K-L). */
	int		mAtomIndices[4];
	/** @brief	Boolean flag defines whether to scale barrier heights in harmonics. */
	bool	mBarrierScale;
	/** @brief	Pointer to linked list of Fourier series describing harmonics. */
	cTorsionHarm *mpHarmonic;
} cTorsion;

/** @brief	Non-bonded pair flag: 1-4 connection. */
#define NBPF_NEIGHBOURS_14		( 1 << 0 )
/** @brief	Non-bonded pair flag: both atoms are heavy. */
#define NBPF_BOTH_HEAVY			( 1 << 1 )

/**
* @brief	Non-bonded pair.
*/
typedef struct stNonBondedPair {
	/** @brief	Indices of atoms comprising the pair. */
	int		mAtomIndices[2];
	/** @brief	Pair flags. */
	int		mFlags;
	/** @brief	Pointer to VdW pair (can be NULL). */
	const cVdWPair *mpVdWPair;
} cNonBondedPair;

/**
* @brief	Hydrogen bond triplet (D = donor, H = hydrogen, A = acceptor).
*/
typedef struct stHBTriplet {
	/** @brief	Indices of atoms comprising the triplet (D-H..A). */
	int		mAtomIndices[3];
	/** @brief	Pointer to hydrogen bond pair (H..A). */
	const cHBPair *mpHBPair;
} cHBTriplet;

/** @brief Atom search functor. */
class AtomSearchFunctor {
public:
	/** 
	* @brief Constructor.
	* @param title : title of the atom to search.
	*/
	AtomSearchFunctor( const char *const title ) : m_title( title ) {}
	/** @brief Compare operator. */
	bool operator() ( const cAtom &at ) { return !strcmp( at.mAtomTitle, m_title ); }
private:
	const char *m_title;
};

/**
* @brief	Atom sort functor.
* @details	Sort atoms by chain, then by residuenumber, then by sortnumber.
*/
class AtomSortFunctor {
public:
	/** @brief Compare operator. */
	bool operator() ( const cAtom &at1, const cAtom &at2 ) { 
		if ( at1.mChainNumber != at2.mChainNumber )
			return at1.mChainNumber < at2.mChainNumber;
		if ( at1.mResidueNumber != at2.mResidueNumber )
			return at1.mResidueNumber < at2.mResidueNumber;
		return at1.mSortNumber < at2.mSortNumber;
	}
};

/**
* @brief	Performance profiling data.
*/
typedef struct stProfileData {
	/** @brief	Profile energy counter. */
	unsigned int	mECounter;
	/** @brief	Non-bonded pair list build profile time (ms). */
	unsigned int	mNBPLTime;
	/** @brief	Bond energy profile time (ms). */
	unsigned int	mEBondTime;
	/** @brief	Angle energy profile time (ms). */
	unsigned int	mEAngleTime;
	/** @brief	Improper energy profile time (ms). */
	unsigned int	mEImproperTime;
	/** @brief	Torsion energy profile time (ms). */
	unsigned int	mETorsionTime;
	/** @brief	Non-bonded energy profile time (ms). */
	unsigned int	mENBTime;
	/** @brief	Special energy profile time (ms). */
	unsigned int	mESpecialTime;
} cProfileData;

/**
* @brief	Object representing a 3D-model of biopolymer.
* @details	This object is used for the following tasks:\n
* -# loading the model from a structure file;
* -# saving current model state to a structure file.
*/
class cModel
{
	DECLARE_SINGLETON( cModel );
	~cModel();

	friend class cThreadManager;
public:
	/**
	* @brief	Type for an array of atoms.
	*/
	typedef std::vector<cAtom> AtomArray;
	/**
	* @brief	Type for an array of bonds.
	*/
	typedef std::vector<cBond> BondArray;
	/**
	* @brief	Type for an array of S-S bonds.
	*/
	typedef std::vector<cSSBond> SSBondArray;
	/**
	* @brief	Type for an array of distant restraints.
	*/
	typedef std::vector<cDistRestrain> DistRestArray;
	/**
	* @brief	Type for an array of angles.
	*/
	typedef std::vector<cAngle> AngleArray;
	/**
	* @brief	Type for an array of impropers.
	*/
	typedef std::vector<cImproper> ImproperArray;
	/**
	* @brief	Type for an array of torsions.
	*/
	typedef std::vector<cTorsion> TorsionArray;
	/**
	* @brief	Type for an array of non-bonded pairs (NBP).
	*/
	typedef std::vector<cNonBondedPair> NBPArray;
	/**
	* @brief	Type for an array of hydrogen bond triplets (HBT).
	*/
	typedef std::vector<cHBTriplet> HBTArray;
	/**
	* @brief	Type for a map of IO objects.
	*/
	typedef std::map<const char*,cDataFormat<cAtom,cModelHeader>*,StringCompareFunctor> IOMap;

private:
	typedef struct {
		cAtom *pAtom1;
		cAtom *pAtom2;
	} SSBondSource;

public:
	/**
	* @brief	Sets name of a source structure file.
	* @details	Default is #DEFAULT_MODEL_FILENAME.
	* @param name : new source structure file name.
	*/
	static void SetFileName( const char *name );
	/**
	* @brief	Returns a name of a source structure file.
	* @return The name of a source structure file.
	*/
	static const char *FileName( void ) { return m_fileName; }
	/**
	* @brief	Sets name of a model.
	* @details	This is a prefix for all output structure files. Default is empty, that means no prefix is added.
	* @param name : new name of a model.
	*/
	static void SetModelName( const char *name );
	/**
	* @brief	Returns a name of a model.
	* @return The name of a model.
	*/
	static const char *ModelName( void ) { return m_modelName; }
	/**
	* @brief	Enable performance profiling (debug feature).
	* @details	When profiling is enabled, most performance-critical subroutines accumulate their execution time. The statistics is displayed when program terminates.
	*/
	static void EnableProfiling( void ) { m_bProfiling = true; }
	/**
	* @brief	Initialize IO objects for different data formats.
	* @details	This function should be modified upon implementation of a new data format.
	*/
	static void InitializeIO( void );
	/**
	* @brief	Static function to access "CalcDeltaEnergy" member of model object from mathlib.
	* @details	This function is not thread-safe.
	* @param delta : coordinate offset to apply temporarily to the model.
	* @return The energy of the structure with coordinate offset applied.
	*/
	static vec_t EnergyMinimizationFunction( vec_t delta );
	/**
	* @brief	Loads a source structure file into memory.
	* @details	Structure file type is determined by file name extension, default is @b PDB.
	* This also prepares the structure for further operation.
	*/
	void Load( void );
	/**
	* @brief	Saves current snapshot to the structure file.
	* @details	Structure file type is determined by file name extension, default is @b PDB.
	* @param filename : name of the structure file to write snapshot to.
	*/
	void Save( const char *const filename );
	/**
	* @brief	Clears the model, frees all data allocated.
	*/
	void Clear( void );
	/**
	* @brief	Initialize persistant modelling parameters.
	*/
	void InitializeParams( void );
	/**
	* @brief	Builds model's topology.
	* @details	The function performs the following tasks:
	*	-# restores missing heavy atoms;
	*	-# restores missing hydrogen atoms;
	*	-# initializes arrays of bonds, angles, impropers and torsions.
	* This must be called after the model is loaded.
	*/
	void BuildTopology( void );
	/**
	* @brief	Calculates energy for the structure.
	* @details	The function performs the following tasks:
	*	-# calculates total energy;
	*	-# writes PDB file with remarks.
	* @param initial : true if this is an initial energy calculation.
	*/
	void CalcEnergy( bool initial );
	/**
	* @brief	Optimizes the structure using local energy minimization.
	* @details	The function performs optimization using conjugate gradients method.
	* @param nSteps : maximum number of optimization steps.
	*/
	void Optimize( int nSteps );
	/**
	* @brief	Print performance profiling results to the log file.
	*/
	void ProfileReport( void );

private:
	void ProfileStart( void );
	void ProfileEnd( unsigned int *pCounter );
	int LookupAtomIndex( const char *atName, int residueSequenceNum, int chainId );
	vec_t CalcDihedral( const vec_t *v1, const vec_t *v2, const vec_t *v3, const vec_t *v4 );
	void RodriguesGibbs( const vec3_t *prevCoords, const vec_t r, const vec_t theta, const vec_t phi, vec_t *outPosition );
	void GetDummyCoords( const AtomArray &residueAtoms, vec3_t *duCoords );
	void GetBackboneCoords( const AtomArray &residueAtoms, vec3_t *backboneCoords );
	cAtom *RestoreAtom( const cResidueLocation *pResidueInfo, const cResidueAtom *pResidueAtom, const vec3_t *backboneCoords, AtomArray &residueAtoms );
	void RenameResidue( int residueNumber, const char *newName );
	int RestoreResidueHeavyAtoms( AtomArray &residueAtoms, AtomArray &newAtoms, const vec3_t *backboneCoords );
	int RestoreResidueHydrogens( AtomArray &residueAtoms, AtomArray &newAtoms, const vec3_t *backboneCoords );
	void RestoreHeavyAtoms( std::vector<bool> &resInfo );
	void RestoreHydrogens( const std::vector<bool> &resInfo );
	void InitializeSSBonds( void );
	void InitializeBonds( void );
	void InitializeRestraints( void );
	bool IsBondedConnection( const AtomArray::const_iterator &atom1, const AtomArray::const_iterator &atom2 );
	bool IsBondedConnection( const cAtom &atom1, const cAtom &atom2 );
	bool IsLoopConnection( const AtomArray::const_iterator &atom1, const AtomArray::const_iterator &atom2 );
	void AddRingImproper( AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3, AtomArray::iterator &atom4 );
	void FindBonds( const AtomArray::iterator &start, AtomArray::iterator &atom1 );
	void FindAngles( const AtomArray::iterator &start, AtomArray::iterator &atom1, AtomArray::iterator &atom2 );
	void FindLoopAngles( const AtomArray::iterator &start, const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2 );
	void FindImpropers( const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3 );
	void FindTorsions( const AtomArray::iterator &start, const AtomArray::iterator &end, AtomArray::iterator &atom1, AtomArray::iterator &atom2, AtomArray::iterator &atom3, bool checkSort );
	void FindSSBonds( void );
	bool InitAtomSolvation( cAtom *pAtom );
	void InitAtomPhysics( void );
	void BuildNonBondedPairLists( const vec4_t *pPositions, const bool nbpR1, const bool nbpR2 );
	void BuildNonBondedPairListsIndividual( int num, int threadNum );
	void InitQEq( void );
	void CalcQEq( const vec4_t *pPositions, const int firstAtom, const int numAtoms, vec_t *pCharges );
	void CalcPartialCharges( const vec4_t *pPositions );
	void CalcNonBondedEnergyIndividualR1( int num, int threadNum );
	void CalcNonBondedEnergyIndividualR2( int num, int threadNum );
	void CalcHBondEnergyIndividual( int num, int threadNum );

	vec_t CalcBondEnergy( const cBond &source, const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcAngleEnergy( const cAngle &source, const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcImproperEnergy( const cImproper &source, const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcTorsionEnergy( const cTorsion &source, const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcSpecialEnergy( const vec4_t *pPositions, vec4_t *pForces, vec_t *pPRPart, vec_t *pDRPart );
	vec_t CalcVdWEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces );
	vec_t CalcElectrostaticEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces );
	vec_t CalcHBondEnergy( const cHBTriplet *pTriplet, const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcSolvationEnergy( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces );
	vec_t CalcSolvationEnergy_GS( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces );
	vec_t CalcSolvationEnergy_GB( const cNonBondedPair *pPair, const vec_t *pvecDelta, const vec_t flDistance, const vec_t flDistanceInv, vec4_t *pForces );
	vec_t CalcNonBondedEnergy( const vec4_t *pPositions, vec4_t *pForces, vec_t *pVdWPart, vec_t *pHBPart, vec_t *pSolvPart, vec_t *pESPart );
	vec_t CalcEnergyAndForces( const vec4_t *pPositions, vec4_t *pForces );
	vec_t CalcDeltaEnergy( vec_t delta );

private:
	static char		m_fileName[MAX_OSPATH];
	static char		m_modelName[MAX_OSPATH];
	static IOMap	m_ioMap;
	static bool		m_bProfiling;

private:
	cModelHeader	m_header;
	AtomArray		m_atoms;
	BondArray		m_bonds;
	SSBondArray		m_SSbonds;
	DistRestArray	m_distRestraints;
	AngleArray		m_angles;
	ImproperArray	m_impropers;
	TorsionArray	m_torsions;
	NBPArray		m_nbPairsR1;
	NBPArray		m_nbPairsR2;
	HBTArray		m_hbTriplets;
	cPhysAtom		*m_physAtoms;
	vec4_t			*m_physPositions;
	vec4_t			*m_physForces;
	vec4_t			*m_physForcesMT;
	cConnectivity	*m_connectivity;
	cSolvParms		*m_solvation;
	cQEqParms		*m_qeqparms;
	vec_t			*m_charges;
	cVdWPair		*m_vdWpairs;
	int				m_vdWpairOffset;
	cHBPair			*m_HBpairs;
	int				m_HBpairOffset;
	cProfileData	m_prof;
	int				m_profDepth;
	int				m_profTime[MAX_PROFILER_DEPTH];
	int				m_iNumThreads;

	// Persistant modelling parameters
	vec_t			m_flCutoffR1;
	vec_t			m_flCutoffR1Squared;
	vec_t			m_flSwitchR1;
	vec_t			m_flSwitchR1Squared;
	vec_t			m_flCutoffR2;
	vec_t			m_flCutoffR2Squared;
	vec_t			m_flCRF_A;
	vec_t			m_flCRF_B;
	vec_t			m_flRDIEInv;
	vec_t			m_flCutoffHB;
	vec_t			m_flCutoffHBSquared;
	vec_t			m_flHBScale;
	vec_t			m_flSolvGref;
	int				m_iRDIE;
	int				m_iHBModel;
	int				m_iSolvModel;
	int				m_iAutoQ;

	// Topology private data
	int				m_iCurrentChain;
	int				m_iCurrentResidueNumber;
	int				m_iCurrentResidueSequence;
	eResidueLocation m_eCurrentResidueLocation;
	const char		*m_pszCurrentResidueTitle;
	int				m_iMaxAtomsInResidue;

	// Generic multithreaded calculation private data
	const vec4_t	*m_pMTPositions;
	vec4_t			*m_pMTForces[MAXTHREADS];
	vec_t			m_flMTVdWEnergy[MAXTHREADS];
	vec_t			m_flMTESEnergy[MAXTHREADS];
	vec_t			m_flMTHBondEnergy[MAXTHREADS];
	vec_t			m_flMTSolvEnergy[MAXTHREADS];

	// NBPL private data
	NBPArray		m_nbTempPairsR1[MAXTHREADS];
	NBPArray		m_nbTempPairsR2[MAXTHREADS];
	HBTArray		m_hbTempTriplets[MAXTHREADS];
	bool			m_NBPR1;
	bool			m_NBPR2;
	bool			m_HBTUp;
	bool			m_NBPInit;

	// Energy optimization private data
	static cModel	*g_pModel;
	vec4_t			*m_enOptPositions;
	vec4_t			*m_enOptForces;
	vec4_t			*m_enOptTemp;

	// Charge equilibration private data
	vec_t			*m_jValues;
	vec_t			*m_cValues;
	int				*m_iValues;
};

/**
* @brief	Helper function to get global model singleton.
* @return Reference to the global model object.
*/
inline cModel& Model( void )
{
	return cModel::Instance();
}

#endif /*MD_MODEL_H*/
