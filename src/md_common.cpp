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
#include "md_main.h"

/**
* @brief	Allocate 16-byte aligned memory.
* @param size : number of bytes to allocate.
*/
void *COM_AlignedMalloc( size_t size )
{
	char *ptr, *ptr2, *aligned_ptr;

	ptr = (char *)malloc(size + 16 + sizeof(size_t));
	if (!ptr) return NULL;

	ptr2 = ptr + sizeof(size_t);
	aligned_ptr = ptr2 + (16 - ((size_t)ptr2 & 15));

	ptr2 = aligned_ptr - sizeof(size_t);
	*((size_t *)ptr2) = (size_t)(aligned_ptr - ptr);

	return aligned_ptr;
}

/**
* @brief	Free 16-byte aligned memory.
* @param baseptr : pointer to data.
*/
void COM_AlignedFree( void *baseptr ) 
{
	char *ptr = (char*)baseptr - *((size_t*)baseptr - 1);
	free( ptr );
}

/**
* @brief	Get file extension.
* @param name : file name with extension.
* @return File extension.
*/
const char *COM_FileExtension( const char *const name )
{
	static char ext[8];
	int i;

	if ( !name )
		return "";

	const char *s = name + strlen(name) - 1;

	while ( s != name && *s != '.' )
		s--;
	if ( s == name )
		return "";
	
	for ( i = 0; i < 7 && *s; ++i, ++s )
		ext[i] = *s;

	ext[i] = 0;
	return ext;
}

/**
* @brief	Trim whitespaces from the string.
* @details	The function operates in place, and the original string becomes no more valid.
* @param string : zero-terminated string to trim.
*/
void COM_Trim( char *string )
{
	size_t counter = 0;
	int len = (int)strlen( string );
	char *start = string;
	char *end = string + len - 1;

	// trim left
	for ( int i = 0; i < len; ++i, ++start ) {
		if ( (byte)(*start) > 32 )
			break;
	}

	// trim right
	for ( int i = len-1; i >= 0; --i, --end ) {
		if ( (byte)(*end) > 32 )
			break;
	}

	// check if string was trimmed away
	if ( end < start ) {
		string[0] = '\0';
		return;
	}

	// check if trimming is unneeded
	if ( start == string && end == string + len - 1 )
		return;

	// build new string
	do {
		string[counter++] = *start++;
	} while ( start <= end );
	string[counter] = '\0';
}

/**
* @brief	Extract a substring.
* @details	The function operates in place, and the original string becomes no more valid.
* @param string : zero-terminated source string.
* @param start : first character to extract.
* @param count : number of characters to extract.
* @param buffer : buffer to place the extracted substring.
* @param bufferSize : size of the buffer.
*/
void COM_Substr( const char *const string, int start, int count, char *buffer, size_t bufferSize )
{
	memset( buffer, 0, bufferSize );

	int len = (int)strlen( string );
	if ( len <= start )
		return;

	count = std::min( count, (int)bufferSize - 1 );
	count = std::min( len - start, count );

	for ( int i = 0; i < count; ++i, ++buffer )
		*buffer = string[start+i];
}

/**
* @brief	Returns a string allocated on heap.
* @param src : source string.
* @return Pointer to the new string.
*/
char *COM_AllocString( const char *const src )
{
	if ( !src )
		return NULL;

	size_t slen = strlen(src);
	char *s = new char[slen+1];
	if ( slen ) memcpy( s, src, slen );
	s[slen] = '\0';
	return s;
}

/**
* @brief	Returns length of a file in bytes.
* @param fp : file handle (FILE).
* @return File length (in bytes).
*/
size_t COM_FileLength( FILE *fp )
{
	long pos, end;
	
	pos = ftell( fp );
	fseek( fp, 0, SEEK_END );
	end = ftell( fp );
	fseek( fp, pos, SEEK_SET );
	return end;
}

/**
* @brief	Return home path name, using xxxHOME environment variable.
* @return pointer to string containing home path (or default home path).
*/
const char *COM_GetHomePath( void )
{
	static char *homepath = NULL;

	if ( !homepath ) {
		size_t iReturnSize = 0;
		homepath = new char[MAX_OSPATH];
		if ( getenv_s( &iReturnSize, homepath, MAX_OSPATH, PROGRAM_NAME_CAP "HOME" ) ) {
			memset( homepath, 0, MAX_OSPATH );
		}
		if ( !homepath[0] ) {
#if defined(WIN32)
			if ( SHGetSpecialFolderPath( 0, homepath, CSIDL_PROGRAM_FILES, FALSE ) ) {
				strcat_s( homepath, MAX_OSPATH, "\\" PROGRAM_NAME );
			} else {
				strcpy_s( homepath, MAX_OSPATH, "C:\\" PROGRAM_NAME );
			}
#else
			strcpy_s( homepath, MAX_OSPATH, "/usr/local/bin/" PROGRAM_NAME );
#endif
		}
		Log().DPrintf( "Home path:\t\t%s\n", homepath );
	}
	return homepath;
}

/**
* @brief	Open file from the program's local directory.
* @param fp : pointer to file handle.
* @param filename : file name.
* @param mode : file access mode ("r", "w", etc.).
* @return error code.
*/
int COM_fopen_local( FILE **fp, const char *filename, const char *mode )
{
	char fullpath[MAX_OSPATH];

	if ( fopen_s( fp, filename, mode ) ) {
		strcpy_s( fullpath, COM_GetHomePath() );
		strcat_s( fullpath, "/" );
		strcat_s( fullpath, filename );
		return fopen_s( fp, fullpath, mode );
	}
	return 0;
}

/**
* @brief	Parses a token out of a string.
* @param data : pointer to current file data.
* @param crossLine : set whether to allow reading a token from the next line.
* @param token : pointer to token buffer to be filled (can be NULL, to skip the token).
* @param tokenSize : maximum size of the token buffer.
* @param lineCounter : holds current line number of the parser.
* @return pointer to the next file data to parse, NULL if nothing left unparsed.
*/
byte *COM_Parse( byte *data, bool crossLine, char *token, int tokenSize, int &lineCounter )
{
	int c;

	if ( !data )
		return NULL;
	
	int len = 0;

	if ( token )
		token[0] = 0;
			
	// skip whitespace
skipwhite:
	while ( ( c = *data ) <= ' ' ) {
		if ( c == 0 ) {
			if ( token ) token[len] = 0;
			return NULL;
		}
		if ( c == '\n' ) {
			if ( !crossLine ) {
				if ( token ) token[len] = 0;
				return NULL;
			}
			++lineCounter;
		}
		++data;
	}
	
	// skip C++ style comments
	if ( c == '/' && data[1] == '/' ) {
		while ( *data && *data != '\n' ) ++data;
		goto skipwhite;
	}
	// skip backslash comments
	if ( c == '\\' && data[1] == '\\' ) {
		while ( *data && *data != '\n' ) ++data;
		goto skipwhite;
	}
	// skip C style comments
	if ( c == '/' && data[1] == '*' ) {
		// Skip "/*"
		data += 2;
		while ( *data  ) {
			if ( *data == '*' && data[1] == '/' ) {
				data += 2;
				break;
			}
			if ( *data == '\n' )
				++lineCounter;
			++data;
		}
		goto skipwhite;
	}
	
	// handle quoted strings specially
	if ( c == '\"' ) {
		++data;
		while ( 1 ) {
			c = *data++;
			if ( c == '\"' || !c ) {
				if ( token ) token[len] = 0;
				return data;
			}
			if ( c == '\n' )
				++lineCounter;
			if ( token && len < tokenSize-1 ) {
				token[len] = c;
				++len;
			}
		}
	}

	// parse single characters
	if ( c == '{' || c == '}' || c == ')' || c == '(' || c == ',' || c == '=' ) {
		if ( token && len < tokenSize-1 ) {
			token[len] = c;
			++len;
		}
		if ( token ) token[len] = 0;
		return data+1;
	}

	// parse a regular word
	do {
		if ( token && len < tokenSize-1 ) {
			token[len] = c;
			++len;
		}
		data++;
		c = *data;
		if ( c == '{' || c == '}' || c == ')' || c == '(' || c == ',' || c == '=' )
			break;
	} while ( c > ' ' );
	
	if ( token ) 
		token[len] = 0;

	return data;
}

/**
* @brief	Gets the next token and checks whether it matches the token passed.
* @param data : pointer to current file data.
* @param crossLine : set whether to allow reading a token from the next line.
* @param token : token to match.
* @param lineCounter : holds current line number of the parser.
* @return pointer to the next file data to parse, NULL if token was not matched.
*/
byte *COM_MatchToken( byte *data, bool crossLine, const char *const token, int &lineCounter )
{
	char localtoken[1024];
	byte *newdata = COM_Parse( data, crossLine, localtoken, sizeof(localtoken), lineCounter );

	if ( !newdata || !newdata[0] ) {
		Log().Error( "end of file found before the token `%s' was matched\n", token );
		return NULL;
	}

	if ( strcmp( localtoken, token ) != 0 ) {
		Log().Error( "expected `%s', found `%s' at line #%i\n", token, localtoken, lineCounter );
		return NULL;
	}

	return newdata;
}

/**
* @brief	Checks whether a token is available in the file data.
* @param data : pointer to current file data.
* @param crossLine : set whether to allow reading a token from the next line.
* @return true if token is available, false otherwise.
*/
bool COM_TokenAvailable( byte *data, bool crossLine )
{
	int temp;
	byte *newdata = COM_Parse( data, crossLine, NULL, 0, temp );
	return ( newdata != NULL );
}

/**
* @brief	Return a number of milliseconds passed from program initialization.
* @details	Initialization in performed on the first call to this function.
* @return number of milliseconds.
*/
int COM_Milliseconds( void )
{
	static unsigned int sys_timeBase = 0;

#if defined(WIN32)
	if (!sys_timeBase) {
		sys_timeBase = timeGetTime();
		return sys_timeBase;
	}

	return timeGetTime() - sys_timeBase;
#else
	struct timeval tp;
	gettimeofday(&tp, NULL);

	if (!sys_timeBase) {
		sys_timeBase = tp.tv_sec;
		return tp.tv_usec/1000;
	}

	int curtime = (tp.tv_sec - sys_timeBase)*1000 + tp.tv_usec/1000;
	return curtime;
#endif
}

static void COM_SecondsToDHMS( unsigned int elapsed_time, unsigned int& days, unsigned int& hours, unsigned int& minutes, unsigned int& seconds )
{
	seconds = elapsed_time % 60;
	elapsed_time /= 60;

	minutes = elapsed_time % 60;
	elapsed_time /= 60;

	hours = elapsed_time % 24;
	elapsed_time /= 24;

	days = elapsed_time;
}

/**
* @brief	Print time elapsed from program startup to log file in human-readable form.
*/
void COM_LogElapsedTime( void )
{
    unsigned int days = 0;
    unsigned int hours = 0;
    unsigned int minutes = 0;
    unsigned int seconds = 0;
	vec_t totalseconds = vec_t( COM_Milliseconds() * 0.001 );

    COM_SecondsToDHMS( (unsigned int)totalseconds, days, hours, minutes, seconds );

    if (days) {
		Log().DPrintf( "%.2f seconds elapsed [%ud %uh %um %us]\n", totalseconds, days, hours, minutes, seconds );
    } else if (hours) {
        Log().DPrintf( "%.2f seconds elapsed [%uh %um %us]\n", totalseconds, hours, minutes, seconds );
    } else if (minutes) {
        Log().DPrintf( "%.2f seconds elapsed [%um %us]\n", totalseconds, minutes, seconds );
    } else {
        Log().DPrintf( "%.2f seconds elapsed\n", totalseconds );
    }
}


static const char *string_contains( const char *str1, const char *str2, int casesensitive ) 
{
	int len = (int)(strlen(str1) - strlen(str2));
	for ( int i = 0; i <= len; ++i, ++str1 ) {
		int j = 0;
		for ( ; str2[j]; ++j ) {
			if ( casesensitive ) {
				if ( str1[j] != str2[j] )
					break;
			} else {
				if ( toupper( str1[j] ) != toupper( str2[j] ) )
					break;
			}
		}
		if ( !str2[j] )
			return str1;
	}
	return NULL;
}

static bool str_filter( const char *filter, const char *name, int casesensitive )
{
	char buf[1024];
	const char *ptr;
	int i, found;

	while(*filter) {
		if (*filter == '*') {
			filter++;
			for (i = 0; *filter; i++) {
				if (*filter == '*' || *filter == '?') break;
				buf[i] = *filter;
				filter++;
			}
			buf[i] = '\0';
			if (strlen(buf)) {
				ptr = string_contains(name, buf, casesensitive);
				if (!ptr) return false;
				name = ptr + strlen(buf);
			}
		}
		else if (*filter == '?') {
			filter++;
			name++;
		}
		else if (*filter == '[' && *(filter+1) == '[') {
			filter++;
		}
		else if (*filter == '[') {
			filter++;
			found = false;
			while(*filter && !found) {
				if (*filter == ']' && *(filter+1) != ']') break;
				if (*(filter+1) == '-' && *(filter+2) && (*(filter+2) != ']' || *(filter+3) == ']')) {
					if (casesensitive) {
						if (*name >= *filter && *name <= *(filter+2)) found = true;
					}
					else {
						if (toupper(*name) >= toupper(*filter) &&
							toupper(*name) <= toupper(*(filter+2))) found = true;
					}
					filter += 3;
				}
				else {
					if (casesensitive) {
						if (*filter == *name) found = true;
					}
					else {
						if (toupper(*filter) == toupper(*name)) found = true;
					}
					filter++;
				}
			}
			if (!found) return false;
			while(*filter) {
				if (*filter == ']' && *(filter+1) != ']') break;
				filter++;
			}
			filter++;
			name++;
		}
		else {
			if (casesensitive) {
				if (*filter != *name) return false;
			}
			else {
				if (toupper(*filter) != toupper(*name)) return false;
			}
			filter++;
			name++;
		}
	}
	return true;
}

/**
* @brief	Get the information on a matching first file, or -1 if there is no match.
* @param dir : directory to search files in.
* @param mask : search mask.
* @param pdata : pointer to output file information.
* @return search handle, -1 if no files were found.
*/
intptr_t COM_FindFirst( const char *dir, const char *mask, cFindData *pdata )
{
#if defined(WIN32)
	struct _finddata_t findinfo;
	char search[MAX_OSPATH];

	strncpy_s(search, dir, MAX_OSPATH-1);
	strncat_s(search, "/", MAX_OSPATH-1);
	strncat_s(search, mask, MAX_OSPATH-1);

	intptr_t h = _findfirst( search, &findinfo );
	if ( h == -1 )
		return -1;

	strncpy_s( pdata->mFileName, findinfo.name, MAX_OSPATH-1 );
	pdata->mAttrib = findinfo.attrib;
	pdata->mTimeCreate = findinfo.time_create;
	pdata->mTimeWrite = findinfo.time_write;
	return h;
#else
	DIR	*fdir;
	struct dirent *d;
	struct stat st;
	char filename[MAX_OSPATH];
	int iAttribs;
	bool bFound = false;

	if ( ( fdir = opendir(dir) ) == NULL )
		return -1;

	while ( ( d = readdir( fdir ) ) != NULL ) {
		iAttribs = 0;
		sprintf_s( filename, sizeof(filename), "%s/%s", dir, d->d_name );
		if ( stat( filename, &st ) == -1 )
			continue;

		if ( st.st_mode & S_IFDIR ) {
			if ( !_stricmp(d->d_name, ".") || !_stricmp(d->d_name, "..") ) 
				continue;
			iAttribs |= 0x10;
		}

		if ( !str_filter( mask, d->d_name, false ) )
			continue;

		strncpy_s( pdata->mFileName, d->d_name, MAX_OSPATH-1 );
		pdata->mAttrib = iAttribs;
		pdata->mTimeCreate = st.st_ctime;
		pdata->mTimeWrite = st.st_mtime;
		bFound = true;
		break;
	}
	
	if ( !bFound )
		return -1;

	return (intptr_t)fdir;
#endif
}

/**
* @brief	Get the information on the next matching file, or -1 if there is no more matches.
* @param dir : directory to search files in.
* @param mask : search mask.
* @param handle : search handle returned by #COM_FindFirst.
* @param pdata : pointer to output file information.
* @return 0 (success) or -1 (fail).
*/
int COM_FindNext( const char *dir, const char *mask, intptr_t handle, cFindData *pdata )
{
#if defined(WIN32)
	struct _finddata_t findinfo;

	UNREFERENCED_PARAMETER( dir );
	UNREFERENCED_PARAMETER( mask );

	if ( handle == -1 )
		return -1;

	int h = _findnext( handle, &findinfo );
	if ( h == -1 )
		return -1;

	strncpy_s( pdata->mFileName, findinfo.name, MAX_OSPATH-1 );
	pdata->mAttrib = findinfo.attrib;
	pdata->mTimeCreate = findinfo.time_create;
	pdata->mTimeWrite = findinfo.time_write;
	return 0;
#else
	DIR	*fdir;
	struct dirent *d;
	struct stat st;
	char filename[MAX_OSPATH];
	int iAttribs;
	bool bFound = false;
	
	fdir = (DIR*)handle;
	if ( !fdir )
		return -1;

	while ( ( d = readdir(fdir) ) != NULL ) {
		iAttribs = 0;
		sprintf_s( filename, sizeof(filename), "%s/%s", dir, d->d_name );
		if ( stat( filename, &st ) == -1 )
			continue;

		if ( st.st_mode & S_IFDIR ) {
			if ( !_stricmp(d->d_name, ".") || !_stricmp(d->d_name, "..") ) 
				continue;
			iAttribs |= 0x10;
		}

		if ( !str_filter( mask, d->d_name, false ) )
			continue;

		strncpy_s( pdata->mFileName, d->d_name, MAX_OSPATH-1 );
		pdata->mAttrib = iAttribs;
		pdata->mTimeCreate = st.st_ctime;
		pdata->mTimeWrite = st.st_mtime;
		bFound = true;
		break;
	}
	
	if ( !bFound )
		return -1;

	return 0;
#endif
}

/**
* @brief	Close the search object.
* @param handle : search handle returned by #COM_FindFirst.
* @return 0 (success) or -1 (fail).
*/
int COM_FindClose( intptr_t handle )
{
#if defined(WIN32)
	return _findclose( handle );
#else
	DIR	*fdir = (DIR*)handle;
	if ( !fdir )
		return -1;
	return closedir( fdir );
#endif
}