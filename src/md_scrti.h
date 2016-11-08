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
#ifndef MD_SECURE_CRT_IMPL_H
#define MD_SECURE_CRT_IMPL_H
/**
* @file
* @brief	Secure CRT implementation.
* @details	The file declares secure CRT function prototypes with a custom implementation for compilers that don't support them (e.g. GCC).
*/
#if !defined(WIN32)

/** @brief	Macro replacement for sprintf_s. */
#define sprintf_s(buffer, buffer_size, stringbuffer, ...)	sprintf(buffer, stringbuffer, __VA_ARGS__)
/** @brief	Macro replacement for _vsnprintf_s. */
#define _vsnprintf_s(buffer, buffer_size, n, format, arg )	vsnprintf(buffer, n, format, arg)
/** @brief	Macro replacement for _stricmp. */
#define _stricmp strcasecmp
/** @brief	Macro replacement for _strnicmp. */
#define _strnicmp strncasecmp

/** @brief	Implementation of fopen_s. */
extern int fopen_s( FILE **f, const char *name, const char *mode );
/** @brief	Implementation of strcpy_s. */
extern int strcpy_s( char *strDestination, size_t numberOfElements, const char *strSource );
/** @brief	Implementation of strcpy_s. */
extern int strcpy_s( char *strDestination, const char *strSource );
/** @brief	Implementation of strncpy_s. */
extern int strncpy_s( char *strDestination, size_t numberOfElements, const char *strSource, size_t count );
/** @brief	Implementation of strncpy_s. */
extern int strncpy_s( char *strDestination, const char *strSource, size_t count );
/** @brief	Implementation of strcat_s. */
extern int strcat_s( char *strDestination, size_t numberOfElements, const char *strSource );
/** @brief	Implementation of strcat_s. */
extern int strcat_s( char *strDestination, const char *strSource );
/** @brief	Implementation of strncat_s. */
extern int strncat_s( char *strDestination, size_t numberOfElements, const char *strSource, size_t count );
/** @brief	Implementation of strncat_s. */
extern int strncat_s( char *strDestination, const char *strSource, size_t count );
/** @brief	Implementation of _strupr_s. */
extern int _strupr_s( char *str );
/** @brief	Implementation of _strupr_s. */
extern int _strupr_s( char *str, size_t numberOfElements );
/** @brief	Implementation of getenv_s. */
extern int getenv_s( size_t *pReturnValue, char* buffer, size_t numberOfElements, const char *varname );

#endif //!WIN32
#endif //MD_SECURE_CRT_IMPL_H
