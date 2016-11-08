/***************************************************************************
* Copyright (C) 2011-2013 Alexander V. Popov.
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
#if !defined(WIN32)
#include "md_main.h"

int fopen_s( FILE **f, const char *name, const char *mode ) 
{
    int ret = 0;
    *f = fopen(name, mode);
    if (!*f) ret = errno;
    return ret;
}

int strcpy_s( char *strDestination, size_t numberOfElements, const char *strSource )
{
	UNREFERENCED_PARAMETER( numberOfElements );
	strcpy( strDestination, strSource );
	return 0;
}

int strcpy_s( char *strDestination, const char *strSource )
{
	strcpy( strDestination, strSource );
	return 0;
}

int strncpy_s( char *strDestination, size_t numberOfElements, const char *strSource, size_t count )
{
	UNREFERENCED_PARAMETER( numberOfElements );
	strncpy( strDestination, strSource, count );
	return 0;
}

int strncpy_s( char *strDestination, const char *strSource, size_t count )
{
	strncpy( strDestination, strSource, count );
	return 0;
}

int strcat_s( char *strDestination, size_t numberOfElements, const char *strSource )
{
	UNREFERENCED_PARAMETER( numberOfElements );
	strcat( strDestination, strSource );
	return 0;
}

int strcat_s( char *strDestination, const char *strSource )
{
	strcat( strDestination, strSource );
	return 0;
}

int strncat_s( char *strDestination, size_t numberOfElements, const char *strSource, size_t count )
{
	UNREFERENCED_PARAMETER( numberOfElements );
	strncat( strDestination, strSource, count );
	return 0;
}

int strncat_s( char *strDestination, const char *strSource, size_t count )
{
	strncat( strDestination, strSource, count );
	return 0;
}

int _strupr_s( char *str )
{
	for ( char *p = str; *p; ++p )
		*p = toupper( *p );
	return 0;
}

int _strupr_s( char *str, size_t numberOfElements )
{
	UNREFERENCED_PARAMETER( numberOfElements );
	return _strupr_s( str );
}

int getenv_s( size_t *pReturnValue, char* buffer, size_t numberOfElements, const char *varname )
{
	const char *env = getenv( varname );
	if ( !env ) return 1;
	strncpy( buffer, varname, numberOfElements-1 );
	if ( pReturnValue ) *pReturnValue = strlen( buffer );
	return 0;
}

#endif //!WIN32
