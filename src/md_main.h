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
#ifndef MD_MAIN_H
#define MD_MAIN_H
/** @mainpage AMBROSIA Software
*
* @author Alexander V. Popov <apopov@niboch.nsc.ru>
*
* [TOC]
*
* @section intro Introduction
* AMBROSIA (<b>A</b>mateur <b>M</b>odeling of <b>B</b>iopolymers: <b>R</b>estoration, <b>O</b>ptimization, <b>S</b>olvation & <b>I</b>nitial <b>A</b>nalysis) is a molecular dynamics modeling tool. It is cross-platform (Win32, Linux) software supporting both 32-bit and 64-bit architectures, and is optimized for modern multi-core CPUs. AMBROSIA performs the following tasks:
* - biopolymer structure refinement (restoration of hydrogens and heavy sidechain atoms);
* - disulfide bond detection;
* - charge equilibration;
* - structure optimization using energy minimization.
*
* The following tasks are under development:
* - explicit solvation;
* - solvent-accessible surface calculation;
* - structure optimization using molecular dynamics simulation;
* - global conformational analsys (RMSD/RMSF, free energy, etc.).
*
* @section ff Force Field
* The program is based on the [AMBER] <b>ff99</b> force field with some parametric (ff99SB, parmbsc0) and functional form (improper, hydrogen bonding) modifications. AMBROSIA force field has the following general form:
* > \f$E = E_{bond} + E_{angle} + E_{torsion} + E_{improper} + E_{nonb} + E_{special}\f$
* Non-bonded energy comprises the hydrophobic (van der Waals), electrostatic, hydrogen bonding and implicit solvation interactions:
* > \f$E_{nonb} = E_{vdW} + E_{coulomb} + E_{H-bond} + E_{solv}\f$
* Special [energy] term is actually an artificial term used to restrain atomic positions or interatomic distances to given values.
* @subsection bond Bond energy
* Covalent bonds are modeled as rigid harmonic oscillators with the following formulation:
* > \f$E_{bond} = \sum K_{hb}(b_{ij}-b_{ij0})^2\f$
* where \f$b_{ij}\f$ is a distance between atoms \f$i\f$ and \f$j\f$ connected with a covalent bond, and \f$b_{ij0}\f$ is an equilibrium bond length. The summation runs over all covalent bonds in the model.
*
* The harmonic constant is defined in a configuration file (amber99.cfg). Please note that [AMBER] defines it as the force constant scaled by two. [AMBER] bond parameters are fully compatible with AMBROSIA.
* @subsection angle Bond angle energy
* Covalent bond angles are modeled as rigid harmonic oscillators with the following formulation:
* > \f$E_{angle} = \sum K_{ha}(\theta_{ijk}-\theta_{ijk0})^2\f$
* where \f$\theta_{ijk}\f$ is an angle between atoms \f$i\f$, \f$j\f$ and \f$k\f$ connected with covalent bonds \f$j-i\f$ and \f$j-k\f$, and \f$\theta_{ijk0}\f$ is an equilibrium angle value. The summation runs over all covalent bond angles in the model.
*
* The harmonic constant is defined in a configuration file (amber99.cfg). Please note that [AMBER] defines it as the force constant scaled by two. [AMBER] bond angle parameters are fully compatible with AMBROSIA.
* @subsection torsion Torsion energy
* Covalent bond torsion energy is usually decomposed into fourier series:
* > \f$E_{torsion} = \sum \sum_{n} \frac{K_{ht}}{s}(1+cos(n\omega-\phi))\f$
* where \f$\omega\f$ is a dihedral angle defined by atoms \f$i\f$, \f$j\f$, \f$k\f$ and \f$l\f$. Number of fourier series (\f$n\f$), barrier scale (\f$s\f$), and phase angle (\f$\phi\f$) are defined in a configuration file. The summation runs over all fourier series defined for the torsion angle, and for all torsion angles in the model.
* If all atoms comprised by the torsion, belong to a planar ring structure, the harmonic constant is not divided by the barrier scale \f$s\f$.
*
* The harmonic constant is also defined in a configuration file (amber99.cfg). Please note that [AMBER] defines it as the force constant scaled by two. [AMBER] torsion angle parameters (including non-degenerate phase angles) are fully compatible with AMBROSIA.
* @subsection improper Improper torsion energy
* Impropers are used to keep certain dihedral structures in a well-known conformation (e.g. planar aromatic rings). Alike bond angles, impropers are modeled as rigid harmonic oscillators with the following formulation:
* > \f$E_{improper} = \sum K_{hi}(\xi_{ijkl}-\xi_{ijkl0})^2\f$
* where \f$\xi_{ijkl}\f$ is a dihedral angle in a tetrahedron defined by atoms \f$i\f$, \f$j\f$, \f$k\f$ and \f$l\f$ with the top atom \f$j\f$, and \f$\xi_{ijkl0}\f$ is an equilibrium dihedral angle value. The summation runs over all improper angles in the model.
*
* The harmonic constant is defined in a configuration file (amber99.cfg). Since [AMBER] doesn't treat improper dihedrals separately, the functional form of the improper energy is derived from [GROMOS].
* Please note that ring impropers are created automatically using information from the topology, with the default harmonic constant of 40.0 kcal/mol*A<sup>2</sup>. There is no need to define any impropers for aromatic rings, since they are implicitly described in the topology.
* @subsection nonb Non-bonded intraction energy
* Non-bonded interactions comprise two primary and two secondary interactions. The primary interactions are:
* - hydrophobic (van der Waals) interaction;
* - electrostatic (Coulomb) interaction.
*
* The secondary (and actually optional) interactions include:
* - hydrogen bonding interaction;
* - interaction with an implicit solvent.
*
* The summation runs in principle over all the non-bonded atomic pairs within certain cut-off radii (except hydrogen bonding interactions, that are summed over triplets). 1-4 neighbours have special scales applied to van der Waals and electrostatic energy and forces: 0.5 for van der Waals, and 1/1.2 for electrostatics. Such scales were documented for the [AMBER] force field.
* 
* <b>Van der Waals</b> interaction has a functional form of 12-6 potential with a soft switch function:
* > \f$E_{vdW} = \sum E_{ij}\left(\frac{R_{ij}^{12}}{r_{ij}^{12}} - 2\frac{R_{ij}^6}{r_{ij}^6}\right)\f$
* where \f$r_{ij}\f$ is a pairwise distance between atoms \f$i\f$ and \f$j\f$, \f$E_{ij} = \sqrt{E_iE_j}\f$, \f$R_{ij} = R_i+R_j\f$, \f$E_i\f$ is a potential well depth for atom \f$i\f$, and \f$R_i\f$ is a van der Waals radius of atom \f$i\f$. Both van der Waals parameters \f$E_i\f$ and \f$R_i\f$ were derived from [AMBER].
* Switch function has the functional form of that utilized in [NAMD][NAMDSF]:
* > \f$F_{switch} = \frac{(R_c^2 - r_{ij}^2)^2(R_c^2 + 2r_{ij}^2 - 3R_s^2)}{(R_c^2 - R_s^2)^3}\f$
* where \f$R_c\f$ is a short cut-off radius, and \f$R_s\f$ is a switch radius (\f$R_s<=R_c\f$; note that if \f$R_s=R_c\f$, the switch function is not applied).
*
* <b>Electrostatic</b> interaction has a form if Coulomb law with reaction field correction:
* > \f$E_{coulomb} = \frac{1}{Dr_{ij}}q_iq_j\left(\frac{1}{r_{ij}}-\frac{C_{rf}r_{ij}^2}{R_c^3}-\frac{1-0.5C_{rf}}{R_c}\right)\f$\n
* > \f$C_{rf} = \frac{(2-2\epsilon)(1+{\kappa}R_c)-\epsilon({\kappa}R_c)^2}{(1+2\epsilon)(1+{\kappa}R_c)+\epsilon({\kappa}R_c)^2}\f$
* where \f$r_{ij}\f$ is a pairwise distance between atoms \f$i\f$ and \f$j\f$, \f$R_c\f$ is a long cut-off radius, \f$\epsilon\f$ is a relative dielectic permittivity (80.0 for water), \f$\kappa\f$ is a Debye radius, and \f$D\f$ is a distance-dependent dielectric constant scale.
* By default, it is equal to 1.0, that corresponds a simple linear distant-dependent constant. When \f$D = 0\f$, the dielectric constant doesn't depend on the distance.
*
* Reaction field correction is a common practice for taking into account screening effects; for example, it is utilized in [GROMOS].
*
* AMBROSIA offers two slightly different models of <b>hydrogen bonding</b> energy: 12-6 and 12-8. The first form has a longer range than the second. The default model is 12-8.
* > \f$E_{H-bond-126} = S_H\sum E_H\left(\frac{R_H^{12}}{r_{ij}^{12}} - 2\frac{R_H^6}{r_{ij}^6}\right)e^{-\frac{(cos\theta-cos\theta_0)^2}{\sigma^2}}\f$\n
* > \f$E_{H-bond-128} = S_H\sum E_H\left(\frac{R_H^{12}}{r_{ij}^{12}} - 1.5\frac{R_H^8}{r_{ij}^8}\right)e^{-\frac{(cos\theta-cos\theta_0)^2}{\sigma^2}}\f$
* where \f$r_{ij}\f$ is a pairwise distance between hydrogen atom and acceptor atom, \f$S_H\f$ is a global hydrogen bonding energy scale (1.0 by default; can be increased when hydrogen bonding is essential for stability of the structure), \f$R_H\f$ and \f$E_H\f$ are parametrized values for equilibrium distance and well depth, \f$\theta\f$ is a triplet angle formed by donor, hydrogen and acceptor atoms, \f$\theta_0\f$ is an ideal angle (\f$\pi\f$ radians), and \f$\sigma\f$ is a constant (\f$cos\frac{5\pi}{6}-cos\pi\f$).
*
* These models use common parameters derived from [BioPASED] software package.
*
* Implicit <b>solvation</b> energy term depends on a solvation model selected and is described [below](@ref solvation).
*
* @subsection special Special energy
* AMBROSIA supports two types of restraining functions:
* - positional restraints;
* - distant restraints.
* 
* Positional restraints allow to keep atoms softly approaching their initial positions. The potential has the following functional form:
* > \f$E_{rp} = \sum K_{hr}(r_i-r_{i0})^2\f$
* where \f$r_i\f$ is the current position of an atom \f$i\f$, and \f$r_{i0}\f$ is an initial position of the atom (before any optimizations). The summation runs over all positionally restrained atoms of the model.
*
* Distant restraints allow to keep two different atoms within the distance \f$d\f$. AMBROSIA utilizes a soft form of the potential (unlike the covalent-bond potential):
* > \f$E_{rd} = \sum \frac{K_{hr}}{2d_{ij0}^2}(d_{ij}^2-d_{ij0}^2)^2\f$
* where \f$d_{ij}\f$ is the distance between atoms \f$i\f$ and \f$j\f$, and \f$d_{ij0}\f$ is the equilibrium distance. The summation runs over all pairs of distantly restrained atoms of the model.
*
* The harmonic constant \f$K_{hr}\f$ of a restraint is, as usual, a force constant scaled by two.
* @section solvation Solvation
* Currently, AMBROSIA supports only the Gaussian solvent-exclusion implicit solvation model, also known as EEF1, described by T. Lazaridis and M. Karplus (Proteins 1999: 35, 133-152) and implemented in [CHARMM] and [BioPASED]. Other implicit (e.g. Generalized Born) and explicit solvation models are under development.
* @subsection impsolv Implicit solvation models
* <b>Gaussian solvent-exclusion model</b> supplements an electrostatic energy term with an additional solvation energy term:
* > \f$E_{solv} = \sum\Delta{G_i}^{ref}-\sum\frac{1}{2\pi\sqrt{\pi}r_{ij}^2}\left\{\frac{\Delta{G_i}^{free}}{\lambda_i}e^{-\left({\frac{r_{ij}-R_i}{\lambda_i}}\right)^2}V_j+\frac{\Delta{G_j}^{free}}{\lambda_j}e^{-\left({\frac{r_{ij}-R_j}{\lambda_j}}\right)^2}V_i\right\}\f$
* where \f$r_{ij}\f$ is a pairwise distance between two heavy atoms, \f$\Delta{G_i}^{ref}\f$ is the reference solvation energy of atom \f$i\f$, \f$\Delta{G_i}^{free}\f$ is the solvation free energy of the free (isolated) atom \f$i\f$, \f$R_i\f$ is a van der Waals radius of atom \f$i\f$ (derived from [AMBER]), and \f$\lambda_i\f$ is a correlation length of atom \f$i\f$. The specific parameters were taken from the reference paper.
*
* The summation runs over the pairs of heavy atoms, for which the solvation parameters are defined within the topology, within the short cut-off radius.
* This implicit solvation model is purely analytical and extremely fast.
* @subsection expsolv Explicit solvation models
* Under development.
*
*
* [AMBER]:		http://ambermd.org/											"AMBER Software"
* [GROMOS]:		http://www.gromos.net/										"GROMOS Software"
* [NAMDSF]:		http://www.ks.uiuc.edu/Research/namd/1.5/ug/node52.html		"Van der Waals switch function in NAMD"
* [CHARMM]:		http://www.charmm.org/										"CHARMM Software"
* [BioPASED]:	http://biopased.niboch.nsc.ru/								"GUI-BioPASED"
*/
/**
* @file
* @brief	Main program file, referencing all necessary includes.
* @details	The file defines common constants and macros, includes headers, etc.
*/
/** @brief	Program short name. */
#define PROGRAM_NAME				"Ambrosia"
/** @brief	Program capitalized short name. */
#define PROGRAM_NAME_CAP			"AMBROSIA"
/** @brief	Program full name. */
#define PROGRAM_TITLE				"AMBROSIA"
/** @brief	Program build string. */
#define PROGRAM_BUILDSTRING			__DATE__ " " __TIME__
/** @brief	Program major version number. */
#define PROGRAM_VERSION_MAJOR		1
/** @brief	Program minor version number. */
#define PROGRAM_VERSION_MINOR		0
/** @brief	Program build configuration string. */
#if defined(_DEBUG)
#define PROGRAM_CONFIGSTRING		"Debug"
#else
#define PROGRAM_CONFIGSTRING		"Release"
#endif
/** @brief	Program OS name. */
#if defined(WIN32)
#define PROGRAM_OSNAME				"Win32"
#if !defined(_M_X64)
#define	PROGRAM_CPUSTRING			"x86"
#else
#define	PROGRAM_CPUSTRING			"x64"
#endif
#elif defined(LINUX)
#define PROGRAM_OSNAME				"Linux"
#if defined(__i386__)
#define	PROGRAM_CPUSTRING			"x86"
#else
#define	PROGRAM_CPUSTRING			"x64"
#endif
#endif

/** @brief	Default log file name (can be overriden with "-o" command-line argument). */
#define DEFAULT_LOG_FILENAME		PROGRAM_NAME ".log"
/** @brief	Default source structure file name (can be overriden with "-c" command-line argument). */
#define DEFAULT_MODEL_FILENAME		"input.pdb"
/** @brief	Default file with simulation parameters name (can be overriden with "-i" command-line argument). */
#define DEFAULT_PARMS_FILENAME		"input.par"

/** @brief	Maximum size of OS file path. */
#define MAX_OSPATH					260

/** @brief	Enable some compatibilities with BioPASED. */
#define BIOPASED_COMPAT_PDB
//#define BIOPASED_COMPAT_ENERGY

#include <cstdio>
#include <cstring>
#include <cassert>
#include <ctime>
#include <cstdarg>
#include <cmath>
#include <cctype>

#include <vector>
#include <list>
#include <set>
#include <map>
#include <algorithm>

#include <setjmp.h>
#include <errno.h>

#if defined(WIN32)
#define NOMINMAX
#include <windows.h>
#include <shlobj.h>
#include <io.h>
#else
#include <unistd.h>
#include <dirent.h>
#include <signal.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#endif

/** @brief	Helper macro to declare a class as a singleton. */
#define DECLARE_SINGLETON(x)	\
	\
	private: x(); \
	static x& _Instance() { static x singleton_instance; return singleton_instance; } \
	public: static x& Instance() { 	typedef x& (*pfn_instance)(); static pfn_instance pf = &_Instance; return pf(); }

/** @brief	Exit program if fatal error has been occured. */
extern void FatalExit( void );

#include "md_types.h"
#include "md_scrti.h"
#include "md_common.h"
#include "md_log.h"
#include "md_config.h"
#include "md_format.h"
#include "md_threads.h"

#endif /*MD_MAIN_H*/
