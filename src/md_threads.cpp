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
#include "md_model.h"

int cThreadManager :: m_iMaxThreads = 0;
int cThreadManager :: m_iCritEnter = 0;

cThreadManager :: cThreadManager()
{
	m_iNumThreads = 1;
	m_pThreadModel = NULL;
	m_ThreadFunc = NULL;
	m_ModelThreadFunc = NULL;
	m_bThreaded = false;
	memset( m_ThreadInfo, 0, sizeof(m_ThreadInfo) );
#if !defined(SINGLETHREADED)
#if defined(WIN32)
	InitializeCriticalSection( &m_hCrit );
	memset( m_hThread, 0, sizeof(m_hThread) );
	memset( m_hEvent, 0, sizeof(m_hEvent) );
	m_hSemaphore = NULL;
#else
	pthread_mutex_init( &m_hMutex, NULL );
	memset( m_hThread, 0, sizeof(m_hThread) );
#endif
#endif //SINGLETHREADED
}

cThreadManager :: ~cThreadManager()
{
	Cleanup();
}

void cThreadManager :: Cleanup( void )
{
#if !defined(SINGLETHREADED)
#if defined(WIN32)
	for ( int i = 0; i < MAXTHREADS; ++i ) {
		if ( m_hThread[i] )
			TerminateThread( m_hThread[i], 0 );
		if ( m_hEvent[i] )
			CloseHandle( m_hEvent[i] );
	}
	if ( m_hSemaphore )
		CloseHandle( m_hSemaphore );
	DeleteCriticalSection( &m_hCrit );
#else
	for ( int i = 0; i < MAXTHREADS; ++i ) {
		if ( m_hThread[i] )
			pthread_kill( m_hThread[i], SIGTERM );
	}
	pthread_cond_destroy( &m_hSemaphore.condBlock );
	pthread_cond_destroy( &m_hSemaphore.condWait );
	pthread_mutex_destroy( &m_hThreadMutex );
	pthread_mutex_destroy( &m_hSemaphore.mutexBlock );
	pthread_mutex_destroy( &m_hSemaphore.mutexWait );
	pthread_mutex_destroy( &m_hMutex );
#endif
#endif //SINGLETHREADED
}

void cThreadManager :: SetMaxThreadCount( int iMaxThreads )
{
	m_iMaxThreads = iMaxThreads;
}

#if defined(WIN32)
DWORD WINAPI cThreadManager :: ThreadEntryStub( LPVOID pParam )
{
	ThreadManager().ThreadFunction( (int)pParam );
	return 0;
}
#else
void *cThreadManager :: ThreadEntryStub( void *pParam )
{
	ThreadManager().ThreadFunction( (int)pParam );
	return NULL;
}
#endif

void cThreadManager :: InitThreads( void )
{
#if defined(SINGLETHREADED)
	m_iNumThreads = 1;
	Log().TPrintf( "Multithreading disabled at compile time (-DSINGLETHREADED)\n" );
#else
#if defined(DEBUGTHREADS)
	Log().TPrintf( "Thread debugging enabled (-DDEBUGTHREADS), expect performance impact!\n" );
#endif
	m_iNumThreads = m_iMaxThreads;

	if ( !m_iNumThreads ) {
#if defined(WIN32)
		SYSTEM_INFO info;
        GetSystemInfo( &info );
        m_iNumThreads = info.dwNumberOfProcessors;
#else
		// poll /proc/cpuinfo
		FILE *fp = NULL;
		if ( fopen_s( &fp, "/proc/cpuinfo", "r" ) ) {
			m_iNumThreads = 1;
		} else {
			char buf[1024];
			memset( buf, 0, sizeof(buf) );
			m_iNumThreads = 0;
			while ( !feof( fp ) ) {
				if ( !fgets( buf, 1023, fp ) )
					break;
				if ( !_strnicmp( buf, "processor", 9 ) )
					++m_iNumThreads;
			}
			fclose( fp );
		}
#endif
	}

	if ( m_iNumThreads < 1 )
		m_iNumThreads = 1;
	else if ( m_iNumThreads > MAXTHREADS )
		m_iNumThreads = MAXTHREADS;

	Log().TPrintf( "Initializing multithreading...\n" );
	Log().DPrintf( "...using %i thread(s)\n", m_iNumThreads );

	if ( m_iNumThreads == 1 ) {
		Log().Verbose( "...single-threaded mode\n" );
		return;
	}

	// create thread sync object
	Log().Verbose( "...creating semaphore\n" );
#if defined(WIN32)
	m_hSemaphore = CreateSemaphore( NULL, 0, m_iNumThreads, NULL );
	if ( !m_hSemaphore ) {
		Log().Error( "Unable to create semaphore\n" );
		Log().Verbose( "...single-threaded mode\n" );
		m_iNumThreads = 1;
		return;
	}
#else
	if ( pthread_mutex_init( &m_hSemaphore.mutexBlock, NULL ) != 0 ||
		 pthread_mutex_init( &m_hSemaphore.mutexWait, NULL ) != 0 ||
		 pthread_cond_init( &m_hSemaphore.condBlock, NULL ) != 0 ||
		 pthread_cond_init( &m_hSemaphore.condWait, NULL ) != 0 ) {
		Log().Error( "Unable to create semaphore\n" );
		Log().Verbose( "...single-threaded mode\n" );
		m_iNumThreads = 1;
		return;
	}
	m_hSemaphore.count = m_iNumThreads;
	m_hSemaphore.maxCount = m_iNumThreads;
#endif

	// initialize all threads and events
	Log().Verbose( "...creating threads\n" );
#if defined(WIN32)
	for ( int i = 0; i < m_iNumThreads; ++i ) {
		DWORD threadid;
		m_hEvent[i] = CreateEvent( NULL, TRUE, FALSE, NULL );
		if ( !m_hEvent[i] ) {
			Log().Error( "Unable to create event %i of %i\n", i + 1, m_iNumThreads );
			Log().Verbose( "...single-threaded mode\n" );
			m_iNumThreads = 1;
			break;
        }
        m_hThread[i] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)&cThreadManager::ThreadEntryStub, (LPVOID) i, CREATE_SUSPENDED, &threadid );
        if ( !m_hThread[i] ) {
			Log().Error( "Unable to create thread %i of %i\n", i + 1, m_iNumThreads );
			Log().Verbose( "...single-threaded mode\n" );
			m_iNumThreads = 1;
			break;
        }
    }
#else
	if ( pthread_mutex_init( &m_hThreadMutex, NULL ) != 0 ) {
		Log().Error( "Unable to create thread mutex\n" );
		Log().Verbose( "...single-threaded mode\n" );
		m_iNumThreads = 1;
		return;
	}
	pthread_attr_t threadAttr;
	pthread_attr_init( &threadAttr );
	for ( int i = 0; i < m_iNumThreads; ++i ) {
		if ( pthread_create( &m_hThread[i], &threadAttr, &cThreadManager::ThreadEntryStub, (void*)i ) == -1 ) {
			Log().Error( "Unable to create thread %i of %i\n", i + 1, m_iNumThreads );
			Log().Verbose( "...single-threaded mode\n" );
			m_iNumThreads = 1;
			break;
        }
	}
#endif
	if ( m_iNumThreads <= 1 )
		return;

	// start all threads
	Log().Verbose( "...starting up\n" );
#if defined(WIN32)
    for ( int i = 0; i < m_iNumThreads; ++i ) {
		if ( ResumeThread( m_hThread[i] ) == 0xFFFFFFFF ) {
			Log().Error( "Unable to start thread %i of %i\n", i + 1, m_iNumThreads );
			Log().Verbose( "...single-threaded mode\n" );
			m_iNumThreads = 1;
			break;
        }
    }
#else
	for ( int i = 0; i < m_iNumThreads; ++i ) {
		if ( pthread_detach( m_hThread[i] ) != 0 ) {
			Log().Error( "Unable to detach thread %i of %i\n", i + 1, m_iNumThreads );
			Log().Verbose( "...single-threaded mode\n" );
			m_iNumThreads = 1;
			break;
        }
    }
#endif
	if ( m_iNumThreads <= 1 )
		return;

	// initial sync
	Log().Verbose( "...initial synchronization\n" );
#if defined(WIN32)
	WaitForMultipleObjects( m_iNumThreads, m_hEvent, TRUE, INFINITE );
	for ( int i = 0; i < m_iNumThreads; ++i )
		ResetEvent( m_hEvent[i] );
#else
	m_bThreaded = true;
	pthread_mutex_lock( &m_hSemaphore.mutexWait );
	if ( m_hSemaphore.count > 0 )
		pthread_cond_wait( &m_hSemaphore.condWait, &m_hSemaphore.mutexWait );
	pthread_mutex_lock( &m_hThreadMutex );
	pthread_mutex_unlock( &m_hSemaphore.mutexWait );
	pthread_mutex_unlock( &m_hThreadMutex );
	m_bThreaded = false;
#endif
	Log().TPrintf( "Multithreading initialized\n" );
#endif //SINGLETHREADED
}

void cThreadManager :: StartThreads( size_t workSize )
{
#if !defined(SINGLETHREADED)
	int partial, residual;

	// build threadinfo
	partial = (int)workSize / m_iNumThreads;
	assert( partial * m_iNumThreads <= (int)workSize );

	for ( int i = 0, start = 0; i < m_iNumThreads; ++i, start += partial ) {
		m_ThreadInfo[i].mStart = start;
		m_ThreadInfo[i].mEnd = start + partial;
	}

	// last thread does the residual
	residual = (int)( workSize - partial * m_iNumThreads );
	m_ThreadInfo[m_iNumThreads-1].mEnd += residual;

#if defined(DEBUGTHREADS)
	Log().DPrintf( "%i threads started (worksize = %i)\n", m_iNumThreads, workSize );
	Log().DPrintf( "Thread work distribution:\n" );
	for ( int i = 0; i < m_iNumThreads; ++i )
		Log().DPrintf( " Thread #%i: %i to %i\n", i, m_ThreadInfo[i].mStart, m_ThreadInfo[i].mEnd );
#endif

	m_bThreaded = true;
#if defined(WIN32)
	ReleaseSemaphore( m_hSemaphore, m_iNumThreads, NULL );
	WaitForMultipleObjects( m_iNumThreads, m_hEvent, TRUE, INFINITE );
	for ( int i = 0; i < m_iNumThreads; ++i )
		ResetEvent( m_hEvent[i] );
#else
	m_hSemaphore.count = m_hSemaphore.maxCount;
#if defined(DEBUGTHREADS)
	Log().DPrintf( "...broadcasting block cond\n" );
#endif
	pthread_cond_broadcast( &m_hSemaphore.condBlock );
#if defined(DEBUGTHREADS)
	Log().DPrintf( "...waiting until semaphore reaches 0\n" );
#endif
	pthread_mutex_lock( &m_hSemaphore.mutexWait );
	if ( m_hSemaphore.count > 0 ) {
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...semaphore is %i\n", m_hSemaphore.count );
#endif
		pthread_cond_wait( &m_hSemaphore.condWait, &m_hSemaphore.mutexWait );
	}
	pthread_mutex_lock( &m_hThreadMutex );
	pthread_mutex_unlock( &m_hSemaphore.mutexWait );
#if defined(DEBUGTHREADS)
	Log().DPrintf( "...semaphore reached 0\n" );
#endif
	pthread_mutex_unlock( &m_hThreadMutex );
#endif
	m_bThreaded = false;
#endif
}

void cThreadManager :: ThreadFunction( int threadIndex )
{
#if !defined(SINGLETHREADED)
	const cThreadInfo *ti = &m_ThreadInfo[threadIndex];
#if defined(WIN32)
	const HANDLE *pHandle = &m_hEvent[threadIndex];
#endif

	do {
#if !defined(WIN32) && defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : do thread work\n", threadIndex );
#endif
		// do thread work
		if ( m_pThreadModel && m_ModelThreadFunc ) {
			assert( m_bThreaded == true );
			for ( int i = ti->mStart; i < ti->mEnd; ++i ) {
				(m_pThreadModel->*m_ModelThreadFunc)( i, threadIndex );
			}
		} else if ( m_ThreadFunc ) {
			assert( m_bThreaded == true );
			for ( int i = ti->mStart; i < ti->mEnd; ++i ) {
				m_ThreadFunc( i, threadIndex );
			}
		}
		// signal that we are done
#if defined(WIN32)
		SignalObjectAndWait( *pHandle, m_hSemaphore, INFINITE, FALSE );
#else
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : [waiting] decrementing semaphore\n", threadIndex );
#endif
		pthread_mutex_lock( &m_hSemaphore.mutexBlock );
		pthread_mutex_lock( &m_hThreadMutex );
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : decrementing semaphore\n", threadIndex );
		if ( m_hSemaphore.count <= 0 )
			Log().DPrintf( "...thread( %i ) : ERROR! semaphore is %i\n", threadIndex, m_hSemaphore.count );
#endif
		m_hSemaphore.count--;
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : semaphore is now %i\n", threadIndex, m_hSemaphore.count );
#endif
		if ( !m_hSemaphore.count )
			pthread_cond_signal( &m_hSemaphore.condWait );
		pthread_mutex_unlock( &m_hSemaphore.mutexBlock );
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : waiting for block cond\n", threadIndex );
#endif
		pthread_cond_wait( &m_hSemaphore.condBlock, &m_hThreadMutex );
		pthread_mutex_unlock( &m_hThreadMutex );
#if defined(DEBUGTHREADS)
		Log().DPrintf( "...thread( %i ) : wait done!\n", threadIndex );
#endif
#endif
	} while ( true );
#endif
}

void cThreadManager :: EnterCriticalSection( void )
{
#if !defined(SINGLETHREADED)
	if ( m_bThreaded ) {
#if defined(WIN32)
		::EnterCriticalSection( &m_hCrit );
#else
		pthread_mutex_lock( &m_hMutex );
#endif
#if defined(DEBUGTHREADS)
		if ( m_iCritEnter )
			printf( "EnterCriticalSection: recursive call\n" );
		++m_iCritEnter;
#endif
	}
#endif
}

void cThreadManager :: LeaveCriticalSection( void )
{
#if !defined(SINGLETHREADED)
	if ( m_bThreaded ) {
#if defined(DEBUGTHREADS)
		if ( !m_iCritEnter )
			printf( "LeaveCriticalSection: called without EnterCriticalSection\n" );
		else
			--m_iCritEnter;
#endif
#if defined(WIN32)
		::LeaveCriticalSection( &m_hCrit );
#else
		pthread_mutex_unlock( &m_hMutex );
#endif
	}
#endif
}
