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
#ifndef MD_PDB_FORMAT_H
#define MD_PDB_FORMAT_H
/**
* @file
* @brief	Declaration of a PDB file IO class.
* @details	The file defines PDB data format IO operations.
*/
/**
* @brief	PDB file format.
* @details	PDB file IO subroutines for #cModelHeader and #cAtom array used by #cModel class.
*/
class cPDBFormat : public cDataFormat<cAtom,cModelHeader>
{
public:
	virtual void ReadFile( const char *const filename, DataHeader *const header, DataArray *const data );
	virtual void WriteFile( const char *const filename, const DataHeader *const header, const DataArray *const data );

private:
	bool IsRecognizableGroup( const char *const resName );
	bool IsRecognizableHG( const char *const resName );
	void AutoAtomSymbol( DataElement &atom );
	void ReadAtom( const char *const buffer, size_t bufferSize, DataArray *const data );
	void WriteAtom( const DataElement &atom, int serialNumber, FILE *fp );

private:
	static const char *m_sRecognizableAA[];
	static const char *m_sRecognizableNA[];
	static const char *m_sRecognizableHG[];

private:
	int		m_atomCounter;
	int		m_heavyAtomCounter;
	int		m_residueCounter;
	int		m_chainCounter;
	int		m_residueChainCounter;
	int		m_oldResidue;
	int		m_lastChainIdentifier;
	int		m_lastResidueSequenceNumber;
	char	m_lastResidueTitle[8];
};

#endif /*MD_PDB_FORMAT_H*/
