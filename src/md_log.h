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
#ifndef MD_LOG_H
#define MD_LOG_H
/**
* @file
* @brief	Declaration of a log class.
* @details	The file defines the log class.
*/
/**
* @brief	Object responsible for logging and error reporting.
* @details	This object is used for the following tasks:\n
* -# printing messages to log file with or without timestamps;
* -# printing messages to standard output (console);
* -# printing warnings and errors;
* -# printing @b fatal @b errors with a termination of the whole program.
*/
class cLog
{
	DECLARE_SINGLETON( cLog );
	~cLog();

public:
	/**
	* @brief	Sets name of a log file. 
	* @details	Default is #DEFAULT_LOG_FILENAME.
	* @param name : new log file name.
	*/
	static void SetFileName( const char *name );
	/**
	* @brief	Returns a name of a log file.
	* @return The name of a log file.
	*/
	static const char *FileName( void ) { return m_fileName; }
	/**
	* @brief	Enables verbose output. 
	*/
	static void EnableVerboseMode( void ) { m_bVerbose = true; }
	/**
	* @brief	Returns whether verbose output is enabled.
	*/
	static bool VerboseMode( void ) { return m_bVerbose; }
	/**
	* @brief	Manually close logging.
	*/
	void Close( void );
	/**
	* @brief	Print formatted text both to the log file and stdout.
	* @param fmt : message / format string.
	*/
	void Printf( const char *fmt, ... );
	/**
	* @brief	Print formatted text to the log file.
	* @param fmt : message / format string.
	*/
	void DPrintf( const char *fmt, ... );
	/**
	* @brief	Print formatted text to the console ONLY.
	* @param fmt : message / format string.
	*/
	void CPrintf( const char *fmt, ... );
	/**
	* @brief	Print formatted text to the log file if verbose mode is set.
	* @param fmt : message / format string.
	*/
	void Verbose( const char *fmt, ... );
	/**
	* @brief	Print formatted text to the log file without timestamp prefix.
	* @param fmt : message / format string.
	*/
	void TPrintf( const char *fmt, ... );
	/**
	* @brief	Print formatted warning message both to the log file and stdout.
	* @details	This also increments warning counter.
	* @param fmt : message / format string.
	*/
	void Warning( const char *fmt, ... );
	/**
	* @brief	Print formatted error message both to the log file and stdout.
	* @details	This also increments error counter.
	* @param fmt : message / format string.
	*/
	void Error( const char *fmt, ... );
	/**
	* @brief	Print formatted error message both to the log file and stdout.
	* @details	This also increments error counter.
	* @details	Program execution is immediately terminated after this call.
	* @param fmt : message / format string.
	*/
	void Fatal( const char *fmt, ... );
	/**
	* @brief	Print newline to log file.
	* @details	Just for formatting in a more readable way.
	*/
	void NewLine( void );

private:
	void Open( void );
	void Print( const char *const prefix, const char *const str, bool addTimeStamp );
	void BuildTimeString( void );

private:
	static char m_fileName[MAX_OSPATH];
	static bool m_bVerbose;
	FILE *m_fp;
	bool m_fatalError;
	int	m_warningCount;
	int	m_errorCount;
	char m_timestr[64];
	char m_logstr[2048];
};

/**
* @brief	Helper function to get global log singleton.
* @return Reference to the global log object.
*/
inline cLog& Log( void )
{
	return cLog::Instance();
}

#endif /*MD_LOG_H*/
