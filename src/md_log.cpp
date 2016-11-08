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

char cLog :: m_fileName[MAX_OSPATH] = DEFAULT_LOG_FILENAME;
bool cLog :: m_bVerbose = false;

cLog :: cLog()
{
	m_fp = NULL;
	m_fatalError = false;
	m_errorCount = 0;
	m_warningCount = 0;
	memset( m_timestr, 0, sizeof(m_timestr) );
	memset( m_logstr, 0, sizeof(m_logstr) );
	Open();
}

cLog :: ~cLog()
{
	Close();
}

void cLog :: BuildTimeString( void )
{
	time_t rawtime;
	time( &rawtime );

#if defined(WIN32)
	struct tm timeinfo;
	struct tm *ptimeinfo = &timeinfo;
	localtime_s( &timeinfo, &rawtime );
#else
	struct tm *ptimeinfo = localtime( &rawtime );
#endif

	int y = 1900 + ptimeinfo->tm_year;
	int m = ptimeinfo->tm_mon+1;
	int d = ptimeinfo->tm_mday;
	int hh = ptimeinfo->tm_hour;
	int mm = ptimeinfo->tm_min;
	int ss = ptimeinfo->tm_sec;
 
	sprintf_s( m_timestr, sizeof(m_timestr), "%02i.%02i.%04i %02i:%02i:%02i", d, m, y, hh, mm, ss );
}

void cLog :: Open( void )
{
	fopen_s( &m_fp, m_fileName, "w" );
	Print( "", PROGRAM_TITLE " initialized\n", true );
	Print( "", "--------------------------------------------------------\n", false );
	DPrintf( PROGRAM_TITLE " %i.%i " PROGRAM_CPUSTRING " " PROGRAM_CONFIGSTRING " " PROGRAM_OSNAME " (build " PROGRAM_BUILDSTRING ")\n", PROGRAM_VERSION_MAJOR, PROGRAM_VERSION_MINOR );
	Print( "", "--------------------------------------------------------\n", false );
}

void cLog :: Close( void )
{
	if ( m_fp ) {
		Print( "", "\n", false );
		if ( m_fatalError )
			Print( "", "*** " PROGRAM_TITLE " terminated due to fatal error! ***\n", true );
		else
			Print( "", PROGRAM_TITLE " successfully finished!\n", true );

		DPrintf( "%i error(s), %i warning(s)\n", m_errorCount, m_warningCount );
		Print( "", "--------------------------------------------------------\n", false );
		COM_LogElapsedTime();

		fclose( m_fp );
		m_fp = NULL;
	}
}

void cLog :: SetFileName( const char *name )
{
	strncpy_s( m_fileName, name, sizeof(m_fileName)-1 );
}

void cLog :: Print( const char *const prefix, const char *const str, bool addTimeStamp )
{
	if ( !m_fp )
		return;

#if defined(WIN32) && defined(_DEBUG)
	OutputDebugString( str );
#endif

	if ( addTimeStamp ) {
		BuildTimeString();
		fprintf( m_fp, "[%s] %s%s", m_timestr, prefix, str );
	} else {
		fprintf( m_fp, "%s%s", prefix, str );
	}
	fflush( m_fp );
}

void cLog :: DPrintf( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	Print( "", m_logstr, false );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: CPrintf( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	fputs( m_logstr, stdout );
	fflush( stdout );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: Verbose( const char *fmt, ... )
{
	if ( m_bVerbose ) {
		va_list	argptr;

		ThreadManager().EnterCriticalSection();

		va_start( argptr, fmt );
		_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
		va_end( argptr );

		Print( "", m_logstr, false );

		ThreadManager().LeaveCriticalSection();
	}
}

void cLog :: TPrintf( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	Print( "", m_logstr, true );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: Printf( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	Print( "", m_logstr, true );
	fputs( m_logstr, stdout );
	fflush( stdout );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: Warning( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	++m_warningCount;
	Print( "WARNING: ", m_logstr, true );

	fputs( "WARNING: ", stdout );
	fputs( m_logstr, stdout );
	fflush( stdout );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: Error( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	++m_errorCount;
	Print( "ERROR: ", m_logstr, true );

	fputs( "ERROR: ", stdout );
	fputs( m_logstr, stdout );
	fflush( stdout );

	ThreadManager().LeaveCriticalSection();
}

void cLog :: Fatal( const char *fmt, ... )
{
	va_list	argptr;

	ThreadManager().EnterCriticalSection();

	va_start( argptr, fmt );
	_vsnprintf_s( m_logstr, sizeof(m_logstr), sizeof(m_logstr)-1, fmt, argptr );
	va_end( argptr );

	++m_errorCount;
	m_fatalError = true;
	Print( "", "\n", false );
	Print( "FATAL ERROR: ", m_logstr, true );

	fputs( "\nFATAL ERROR: ", stdout );
	fputs( m_logstr, stdout );
	fflush( stdout );

	ThreadManager().LeaveCriticalSection();

	FatalExit();
}

void cLog :: NewLine( void )
{
	Print( "", "\n", false );
}


