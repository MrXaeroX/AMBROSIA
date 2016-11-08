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
#ifndef MD_COMMON_H
#define MD_COMMON_H
/**
* @file
* @brief	Common utility functions.
*/

/** @brief Generic "offsetof" macro. */
#ifndef offsetof
#define offsetof(s,m)	(size_t)&(((s *)0)->m)
#endif

/** @brief Generic "UNREFERENCED_PARAMETER" macro. */
#if !defined(UNREFERENCED_PARAMETER)
#define UNREFERENCED_PARAMETER(x) (void)(x)
#endif

/** @brief "assume_0" macro (compiler optimization hint). */
#if defined(_DEBUG)
#define assume_0	assert( 1 == 0 )
#else
#if defined(_MSC_VER)
#define assume_0	__assume( 0 )
#elif defined(__GNUC__)
#define assume_0	__builtin_unreachable()
#else
#define assume_0
#endif
#endif

/** @brief	Type for an array of integers. */
typedef std::vector<int> IntArray;

/** @brief String compare functor for IO object map */
class StringCompareFunctor {
public:
	/** @brief Compare operator. */
	bool operator()( char const *a, char const *b ) {
		return std::strcmp( a, b ) < 0;
	}
};

/** 
* @brief Find pair by key. 
* @tparam T: first pair element type (key).
* @tparam U: first pair element type (value).
*/
template<typename T, typename U>
class FindPairByKeyFunctor {
public:
	/** 
	* @brief Constructor.
	* @param check : key to search for.
	*/
	FindPairByKeyFunctor( T check ) : m_check( check ) {}
	/** @brief Compare operator. */
	bool operator() ( const std::pair<T,U> &check ) { return ( check.first == m_check ); }
private:
	T m_check;
};
/** 
* @brief Find pair by value. 
* @tparam T: first pair element type (key).
* @tparam U: first pair element type (value).
*/
template<typename T, typename U>
class FindPairByValueFunctor {
public:
	/** 
	* @brief Constructor.
	* @param check : value to search for.
	*/
	FindPairByValueFunctor( U check ) : m_check( check ) {}
	/** @brief Compare operator. */
	bool operator() ( const std::pair<T,U> &check ) { return ( check.second == m_check ); }
private:
	U m_check;
};

/** 
* @brief Find data for file search funtions. 
*/
typedef struct stFindData
{
	/** @brief	Attributes. */
	int		mAttrib;
	/** @brief	Creation time. */
	qword	mTimeCreate;
	/** @brief	Modification time. */
	qword	mTimeWrite;
	/** @brief	File name (relative to search path). */
	char	mFileName[MAX_OSPATH];
} cFindData;

extern void *COM_AlignedMalloc( size_t size );
extern void COM_AlignedFree( void *baseptr );
extern const char *COM_FileExtension( const char *const name );
extern void COM_Trim( char *string );
extern void COM_Substr( const char *const string, int start, int count, char *buffer, size_t bufferSize );
extern char *COM_AllocString( const char *const src );
extern size_t COM_FileLength( FILE *fp );
extern const char *COM_GetHomePath( void );
extern int COM_fopen_local( FILE **fp, const char *filename, const char *mode );
extern byte *COM_Parse( byte *data, bool crossLine, char *token, int tokenSize, int &lineCounter );
extern byte *COM_MatchToken( byte *data, bool crossLine, const char *const token, int &lineCounter );
extern bool COM_TokenAvailable( byte *data, bool crossLine );
extern int COM_Milliseconds( void );
extern void COM_LogElapsedTime( void );
extern intptr_t COM_FindFirst( const char *dir, const char *mask, cFindData *pdata );
extern int COM_FindNext( const char *dir, const char *mask, intptr_t handle, cFindData *pdata );
extern int COM_FindClose( intptr_t handle );

/** @brief Helper function to avoid type cast when we use single-precision vec_t. */
inline vec_t COM_Atof( const char *string )
{
	return (vec_t)atof( string );
}

#endif /*MD_COMMON_H*/
