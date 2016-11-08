# AMBROSIA
## Amateur Modeling of Biopolymers: Restoration, Optimization, Solvation & Initial Analysis
AMBROSIA is a molecular dynamics modeling tool. It is cross-platform (Win32, Linux) software supporting both 32-bit and 64-bit architectures, and is optimized for modern multi-core CPUs. AMBROSIA performs the following tasks:
- biopolymer structure refinement (restoration of hydrogens and heavy sidechain atoms);
- disulfide bond detection;
- charge equilibration;
- structure optimization using energy minimization.

The following tasks are under development:
- explicit solvation;
- solvent-accessible surface calculation;
- structure optimization using molecular dynamics simulation;
- global conformational analsys (RMSD/RMSF, free energy, etc.).

## Force Field
The program is based on the AMBER ff99 force field with some parametric (ff99SB, parmbsc0) and functional form (improper, hydrogen bonding) modifications. AMBROSIA force field has the following general form:

E=E<sub>bond</sub>+E<sub>angle</sub>+E<sub>torsion</sub>+E<sub>improper</sub>+E<sub>nonb</sub>+E<sub>special</sub>

Non-bonded energy comprises the van der Waals, electrostatic, hydrogen bonding and implicit solvation interactions:

E<sub>nonb</sub>=E<sub>vdW</sub>+E<sub>coulomb</sub>+E<sub>H−bond</sub>+E<sub>solv</sub>

Optional *special energy* term is actually an artificial term used to restrain atomic positions or interatomic distances to given values.

Currently, AMBROSIA supports only the Gaussian solvent-exclusion implicit solvation model, also known as EEF1, described by T. Lazaridis and M. Karplus (Proteins 1999: 35, 133-152) and implemented in CHARMM and BioPASED. Other implicit (e.g. Generalized Born) and explicit solvation models are under development.

## More Information
- [PDF Manual](../blob/master/AMBROSIA.pdf)
