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
#ifndef MD_FORMAT_H
#define MD_FORMAT_H
/**
* @file
* @brief	Declaration of a data format interface.
* @details	The file defines data format pure abstract class.
*/
/**
* @brief	Abstract data format.
* @tparam T : data structure of element to work with.
* @tparam U : data structure of header to work with.
*/
template<typename T, typename U>
struct cDataFormat
{
	/**
	* @brief	Data type alias.
	*/
	typedef T	DataElement;
	/**
	* @brief	Header type alias.
	*/
	typedef U	DataHeader;
	/**
	* @brief	Array of data type (STL vector).
	*/
	typedef std::vector<T> DataArray;
	/**
	* @brief	Reads array of data from file.
	* @param filename : file name to read.
	* @param header : pointer to header to read in.
	* @param data : pointer to data to read in.
	*/
	virtual void ReadFile( const char *const filename, DataHeader *const header, DataArray *const data ) = 0;
	/**
	* @brief	Writes array of data to file.
	* @param filename : file name to write.
	* @param header : pointer to header to write out.
	* @param data : pointer to data to write out.
	*/
	virtual void WriteFile( const char *const filename, const DataHeader *const header, const DataArray *const data ) = 0;
};

#endif /*MD_FORMAT_H*/
