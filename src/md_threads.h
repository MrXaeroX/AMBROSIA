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
#ifndef MD_THREADS_H
#define MD_THREADS_H
/**
* @file
* @brief	Declaration of a thread class.
* @details	The file defines the thread class.
*/

class cModel;

/** @brief	Maximum number of threads supported. */
#define MAXTHREADS	64

/** @brief	Define this macro to ignore multithreading capabilities (everything will be executed in the main thread). */
//#define SINGLETHREADED

/** @brief	Define this macro to debug multithreading. */
//#define DEBUGTHREADS

/**
* @brief	Thread function pointer.
*/
typedef void (*cThreadFunc)(int, int);
/**
* @brief	Thread class-member function pointer.
*/
typedef void (cModel::*cModelThreadFunc)(int, int);

/**
* @brief	Local thread information.
*/
typedef struct stThreadInfo {
	/** @brief	Work start. */
	int	mStart;
	/** @brief	Work end. */
	int	mEnd;
} cThreadInfo;

#if !defined(WIN32)
/**
* @brief	POSIX semaphore implementation.
*/
typedef struct stSemaphore {
	/** @brief	Blocking mutex for threads. */
	pthread_mutex_t	mutexBlock;
	/** @brief	Blocking mutex for the main thread. */
	pthread_mutex_t	mutexWait;
	/** @brief	Conditional event for threads. */
	pthread_cond_t	condBlock;
	/** @brief	Conditional event for the main thread. */
	pthread_cond_t	condWait;
	/** @brief	Semaphore initial count. */
	int		count;
	/** @brief	Semaphore maximum count. */
	int		maxCount;
} cSemaphore;
#endif

/**
* @brief	Thread manager object incapsulates multithreading execution capabilities.
*/
class cThreadManager
{
	DECLARE_SINGLETON( cThreadManager );
	~cThreadManager();

public:
	/**
	* @brief	Sets maximum number of threads to use.
	* @details	Default is 0 (autodetect processor capabilities).
	* @param iMaxThreads : maximum number of threads ( 0 = autodetect ).
	*/
	static void SetMaxThreadCount( int iMaxThreads );
	/**
	* @brief	Returns maximum number of threads to use.
	* @return The maximum number of threads to use (set in #SetMaxThreadCount).
	*/
	static int MaxThreadCount( void ) { return m_iMaxThreads; }

#if defined(WIN32)
	/**
	* @brief	Win32 thread entry function.
	* @param pParam : thread parameter (thread index).
	*/
	static DWORD WINAPI ThreadEntryStub( LPVOID pParam );
#else
	/**
	* @brief	POSIX thread entry function.
	* @param pParam : thread parameter (thread index).
	*/
	static void *ThreadEntryStub( void *pParam );
#endif

public:
	/**
	* @brief	Stop threads, free objects, etc.
	*/
	void Cleanup( void );
	/**
	* @brief	Initialize threads (determine an actual thread count, create thread and sync objects, etc.).
	*/
	void InitThreads( void );
	/**
	* @brief	Return number of threads in use.
	*/
	int CountThreads( void ) { return m_iNumThreads; }
	/**
	* @brief	Threaded function enters a critical section.
	*/
	void EnterCriticalSection( void );
	/**
	* @brief	Threaded function leaves a critical section.
	*/
	void LeaveCriticalSection( void );
	/**
	* @brief	Wrapper for global (or static class member) threaded function calls.
	* @param	workSize : number of iterations for the worker.
	* @param	threadWorkerFn : worker function.
	*/
	void RunThreadsOn( size_t workSize, void (*threadWorkerFn)(int, int) )
	{
#if defined(SINGLETHREADED)
		for ( int i = 0; i < (int)workSize; ++i ) 
			(*threadWorkerFn)( i , 0 );
#else
		if ( m_iNumThreads <= 1 ) {
			for ( int i = 0; i < (int)workSize; ++i ) 
				(*threadWorkerFn)( i , 0 );
			return;
		}
		m_pThreadModel = NULL;
		m_ThreadFunc = threadWorkerFn;
		StartThreads( workSize );
#endif
	}
	/**
	* @brief	Wrapper for class-member threaded function calls.
	* @param	pObject : class instance.
	* @param	workSize : number of iterations for the worker.
	* @param	threadWorkerFn : worker function (non-static member function of pObject).
	*/
	void RunThreadsOn( cModel *pObject, size_t workSize, void (cModel::*threadWorkerFn)(int, int) )
	{
#if defined(SINGLETHREADED)
		for ( int i = 0; i < (int)workSize; ++i ) 
			(pObject->*threadWorkerFn)( i , 0 );
#else
		if ( m_iNumThreads <= 1 ) {
			for ( int i = 0; i < (int)workSize; ++i ) 
				(pObject->*threadWorkerFn)( i , 0 );
			return;
		}
		m_pThreadModel = pObject;
		m_ModelThreadFunc = threadWorkerFn;
		StartThreads( workSize );
#endif
	}
	/*
	void FakeRunThreadsOn( cModel *pObject, size_t workSize, void (cModel::*threadWorkerFn)(int, int) )
	{
		for ( int i = 0; i < (int)workSize; ++i ) 
			(pObject->*threadWorkerFn)( i , 0 );
	}
	*/
private:
	void ThreadFunction( int threadIndex );
	void StartThreads( size_t workSize );

private:
	static int			m_iMaxThreads;
	static int			m_iCritEnter;

private:
	int					m_iNumThreads;
	bool				m_bThreaded;
	cModel				*m_pThreadModel;
	cThreadFunc			m_ThreadFunc;
	cModelThreadFunc	m_ModelThreadFunc;
	cThreadInfo			m_ThreadInfo[MAXTHREADS];
#if defined(WIN32)
	HANDLE				m_hThread[MAXTHREADS];
	HANDLE				m_hEvent[MAXTHREADS];
	HANDLE				m_hSemaphore;
	CRITICAL_SECTION	m_hCrit;
#else
	pthread_t			m_hThread[MAXTHREADS];
	pthread_mutex_t		m_hThreadMutex;
	cSemaphore			m_hSemaphore;
	pthread_mutex_t		m_hMutex;
#endif
};

/**
* @brief	Helper function to get global thread manager singleton.
* @return Reference to the global thread manager object.
*/
inline cThreadManager& ThreadManager( void )
{
	return cThreadManager::Instance();
}

#endif /*MD_THREADS_H*/
