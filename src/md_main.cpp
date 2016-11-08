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
#include "md_math.h"

static jmp_buf g_fatal_error;

void FatalExit( void )
{
	Log().Close();
	ThreadManager().Cleanup();
	longjmp( g_fatal_error, -1 );
}

static void PrintLogo( void )
{
	printf( "%s %i.%i (build %s)\n", PROGRAM_TITLE, PROGRAM_VERSION_MAJOR, PROGRAM_VERSION_MINOR, PROGRAM_BUILDSTRING );
}

static void PrintSettings( void )
{
	Log().DPrintf( "Floating-point precision: %s\n", sizeof(vec_t) == 8 ? "DBL" : "SNG" );
	Log().DPrintf( "Log file name:\t\t%s\n", cLog::FileName() );
	Log().DPrintf( "Input structure:\t%s\n", cModel::FileName() );
	Log().DPrintf( "Input parameters:\t%s\n", cConfig::ParmsFileName() );
	Log().DPrintf( "User configuration:\t%s\n", cConfig::UserConfigName() );
	Log().DPrintf( "Model name:\t\t%s\n", cModel::ModelName() );
	if ( cThreadManager::MaxThreadCount() ) 
		Log().DPrintf( "Max. # of threads:\t%i\n", cThreadManager::MaxThreadCount() );
	else
		Log().DPrintf( "Max. # of threads:\t(auto)\n" );
}

static void ParseCommandLine( int argc, char *argv[] )
{
	for ( int i = 1; i < argc; ++i ) {
		if ( !strcmp( argv[i], "-o" ) ) {
			// -o <name> : specify log file name
			if ( i >= argc - 1 ) Log().Fatal( "-o: missing log file name\n" );
			cLog::SetFileName( argv[++i] );
		} else if ( !strcmp( argv[i], "-c" ) ) {
			// -c <name> : specify model file name
			if ( i >= argc - 1 ) Log().Fatal( "-c: missing input structure file name\n" );
			cModel::SetFileName( argv[++i] );
		} else if ( !strcmp( argv[i], "-mn" ) ) {
			// -mn <name> : specify model name
			if ( i >= argc - 1 ) Log().Fatal( "-mn: missing model name\n" );
			cModel::SetModelName( argv[++i] );
		} else if ( !strcmp( argv[i], "-i" ) ) {
			// -i <name> : specify parms file name
			if ( i >= argc - 1 ) Log().Fatal( "-i: missing input parameters' file name\n" );
			cConfig::SetParmsFileName( argv[++i] );
		} else if ( !strcmp( argv[i], "-u" ) ) {
			// -u <name> : specify user config file name
			if ( i >= argc - 1 ) Log().Fatal( "-u: missing user config file name\n" );
			cConfig::SetUserConfigName( argv[++i] );
		} else if ( !strcmp( argv[i], "-v" ) ) {
			// -v : enable verbose output to log file
			cLog::EnableVerboseMode();
		} else if ( !strcmp( argv[i], "-prof" ) ) {
			// -prof : enable performance profiling
			cModel::EnableProfiling();
		} else if ( !strcmp( argv[i], "-threads" ) ) {
			// -threads <number> : specify maximum number of threads
			if ( i >= argc - 1 ) Log().Fatal( "-threads: missing a number\n" );
			cThreadManager::SetMaxThreadCount( atoi( argv[++i] ) );
		} else if ( argv[i][0] == '-' ) {
			Log().Warning( "%s: unknown command-line parameter\n", argv[i] );
		} else {
			Log().Warning( "%s: unknown command-line argument\n", argv[i] );
		}
	}
}

static void LoadConfigFiles( void )
{
	cFindData config;
	char configPath[MAX_OSPATH];
	char configName[MAX_OSPATH];
	int result;

	const char *cpszConfigExt = "*.cfg";
	strcpy_s( configPath, "../conf" );

	// search in current path
	intptr_t configSearch = COM_FindFirst( configPath, cpszConfigExt, &config );
	if ( configSearch == -1 ) {
		// search in home path
		strcpy_s( configPath, COM_GetHomePath() );
		strcat_s( configPath, "/conf" );
		configSearch = COM_FindFirst( configPath, cpszConfigExt, &config );
		if ( configSearch == -1 )
			return;
	}

	do {
		sprintf_s( configName, sizeof(configName), "%s/%s", configPath, config.mFileName );
		Config().LoadConfig( configName );
	} while ( ( result = COM_FindNext( configPath, cpszConfigExt, configSearch, &config ) ) == 0 );

	COM_FindClose( configSearch );
	Config().LoadUserConfig();
}

int main( int argc, char *argv[] )
{
	int stepCounter = 0;

	PrintLogo();

	if ( setjmp( g_fatal_error ) )
		return 1;

	// parse command line and print settings
	ParseCommandLine( argc, argv );
	PrintSettings();
	COM_GetHomePath();

	// alignment checks
	if ( sizeof( cPhysAtom ) % 16 != 0 )
		Log().Fatal( "sizeof(cPhysAtom) is not multiple of 16!\n" );

	// initialize time
	COM_Milliseconds();

	// Initialize threads
	ThreadManager().InitThreads();

	// load config files
	LoadConfigFiles();

	// load parameters
	Config().LoadParams();
	Model().InitializeParams();

	// initialize data formats
	cModel::InitializeIO();

	// load model structure and build the topology
	Log().NewLine();
	Log().Printf( "Step %2i: Loading model from file...\n", ++stepCounter );
	Model().Load();

	Log().NewLine();
	Log().Printf( "Step %2i: Building molecular topology...\n", ++stepCounter );
	Model().BuildTopology();

	// calculate energy
	if ( Config().Parameters().mEnergyCalc ) {
		Log().NewLine();
		Log().Printf( "Step %2i: Calculating energy...\n", ++stepCounter );
		Model().CalcEnergy( true );
	}

	// minimize energy
	if ( Config().Parameters().mEnergyOptimize ) {
		Log().NewLine();
		Log().Printf( "Step %2i: Start energy optimization...\n", ++stepCounter );
		Model().Optimize( Config().Parameters().mEnergyOptimizeSteps );
	}

	Model().ProfileReport();
	return 0;
}
