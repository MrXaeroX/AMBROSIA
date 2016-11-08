# AMBROSIA
Tightly Associated Solvent Shell Extractor
Version 1.0

TASSE is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. MDTRA is open-source software; the source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Features:
TODO

Usage:
tasse [arguments]

Arguments:
 -tf  : source topology filename
 -tn  : source topology nature (0 = autodetect)
 -cf  : source coordinate filename
 -cn  : source coordinate nature (0 = autodetect)
 -trf : coordinate trajectory filename
 -trn : coordinate trajectory nature (0 = autodetect)
 -n   : starting snapshot of the trajectory (skip N previous snapshots)
 -o   : output PDB filename
 -i   : output information filename
 -e   : relative dielectric permittivity of the medium (default 80)
 -r   : electrostatic cut-off radius, in angstroms (default 15.0)
 -ce  : electrostatic scale coefficient (default 0.25)
 -d   : h-bond maximum length, in angstroms (default 4.0)
 -h   : h-bond cut-off absolute energy, in kcal/mol (default 1)
 -ch  : h-bond scale coefficient (default 0.75)
 -p   : probability (trajectory occurence) cut-off (default 0.9)
 -ng  : don't group similar donor/acceptor atoms
 -q   : read atomic charges from source (not applicable to PDB nature)
 -s   : solvent residue title (default HOH)
 -t   : number of threads (default is autodetect)
 -low : low thread priority (yield resources to other programs)
 -est : show progress pacifier (estimate completion time)
 -v   : verbose mode (print log messages to the console)
 -x   : don't process trajectory, only convert source to output PDB
 -?   : print help for arguments (this message)
 
 Natures:
 0    : TCT : autodetect
 1    : xxT : PDB list (*.lst) : list of PDB files, each line names a single snapshot file
 2    : TCx : PDB file (*.pdb) : a single shapshot file
 3    : Txx : AMBER prmtop (*.prmtop) : LEaP topology
 4    : xCx : AMBER coords (*.inpcrd, *.rst) : LEaP coordinates or restart file
 5    : xxT : AMBER mdcrd (*.mdcrd) : MD trajectory in text format

TCT means applicable to (T)opology, (C)oordinate, or (T)rajectory natures; x means not applicable.

Notice:
The GUI version of the Program performs hardware accelerated rendering using OpenGL rendering API. If you experience any problems related to plot rendering, make sure you have latest video card drivers installed.

Tools used:
Qt 4.8.7 (LGPL)

Copyright:
The author of the Program is Alexander V. Popov. 
Copyright © 2011-2016. All rights reserved.
