////////////////////////////////////////////////////////////////////////////////
//   wmc4.cc for wmc4                                                         //
//                                                                            //
//   This is main function of the simulation program.                         //
//                                                                            //
//                    - 20. Oct. 2017. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <unistd.h>

#include "WMC4DetCon.hh"
#include "WMC4PriGenAct.hh"
#include "WMC4ActIni.hh"
#include "WMC4ConMan.hh"

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4String.hh"
#include "Randomize.hh"

#include "QGSP_BERT.hh"

#include "TString.h"

// Declaration of PrintHelp()
void PrintHelp();

//////////////////////////////////////////////////
//   Main function                              //
//////////////////////////////////////////////////
int main(int argc, char** argv)
{
	// Read options
	int flag_b = 0, flag_g = 0, flag_h = 0, flag_m = 0, flag_c = 0;
	const char* optDic = "bghm:c:"; // Option dictionary
	int option;
	char* macro;
	char* config;
	while ( (option = getopt(argc, argv, optDic)) != -1 ) // -1 means getopt() parses all options.
	{
		switch ( option )
		{
			case 'b' :
				flag_b = 1;
				break;
			case 'g' :
				flag_g = 1;
				break;
			case 'h' :
				flag_h = 1;
				break;
			case 'm' :
				flag_m = 1;
				macro = optarg;
				break;
			case 'c' :
				flag_c = 1;
				config = optarg;
				break;
			case '?' :
				flag_h = 1;
				break;
		}
	}

	// In case of '-h' option activated
	if ( flag_h )
	{
		PrintHelp();
		return 0;
	}

	// Randomizer
	CLHEP::RanluxEngine defaultEngine(1234567, 4);
	G4Random::setTheEngine(&defaultEngine);
	G4int seed = time(NULL);
	G4Random::setTheSeed(seed);

	// Detect interactive mode (if flag_g) and define UI session
	G4UIExecutive* UI = 0;
	if ( flag_g ) UI = new G4UIExecutive(argc, argv);

	// Get configuration from file
	// File is in the same path with binary
	char buff[128];
	getcwd(buff, 128);
	TString cfgFileName;
	if ( flag_c )
	{
		cfgFileName.Form("%s/%s", buff, config);
	}
	else
	{
		cfgFileName.Form("%s/config.cfg", buff);
	}
	WMC4ConMan* CM = new WMC4ConMan(cfgFileName.Data());
	
	// Run manager
//	G4MTRunManager* RM = new G4MTRunManager();
	G4RunManager* RM = new G4RunManager();
//	G4int nofTRD = CM -> GetNofTRD();
//	RM -> SetNumberOfThreads(nofTRD);

	// Detector construction from configuration (Geometry)
	// We define everything about geomtrical setup in this class.
	WMC4DetCon* DC = new WMC4DetCon(CM);
	RM -> SetUserInitialization(DC);

	// Physics list to be used
	G4VModularPhysicsList* PL = new QGSP_BERT;
	PL -> SetVerboseLevel(0);
	RM -> SetUserInitialization(PL);

	// User actions
	RM -> SetUserInitialization(new WMC4ActIni(CM));

	// Initialize
	RM -> Initialize();

	// Visualization manger
	G4VisManager* VM = new G4VisExecutive();
	VM -> Initialize();

	// Get the pointer to the user interface manager
	G4UImanager* UM = G4UImanager::GetUIpointer();

	// Process macro or start UI session
	G4String command = "/control/execute ";
	G4String command2;
	if ( !UI )
	{
		// terminal mode
		G4UIsession* US = new G4UIterminal(new G4UItcsh);
		if ( flag_m )
		{
			command2 = macro;
			UM -> ApplyCommand(command + command2);
		}
		if ( !flag_b ) US -> SessionStart();

		delete US;
	}
	else
	{
		// Interactive mode
		command2 = "init_vis.mac";
		UM -> ApplyCommand(command + command2);
		if ( flag_m )
		{
			command2 = macro;
			UM -> ApplyCommand(command + command2);
		}
		UI -> SessionStart();

		delete UI;
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program.
	delete VM;
	delete RM;

	std::cout << "bye bye :)" << std::endl;

	return 0;
}

//////////////////////////////////////////////////
//   Print help message                         //
//////////////////////////////////////////////////
void PrintHelp()
{
	std::cout << "usage: wmc4 [-b] [-g] [-m macrofile] [-c configfile]" << std::endl;
	std::cout << std::endl;
	std::cout << "Examples:" << std::endl;
	std::cout << "  wmc4 -b -m myRun.mac -c myConfig.cfg  # Run in batch mode with macro and config." << std::endl;
	std::cout << "  wmc4 -g                               # Run in graphical mode." << std::endl;
	std::cout << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "  -b  Execute in batch mode"         << std::endl;
	std::cout << "  -g  Execute in graphical mode"     << std::endl;
	std::cout << "      Note: Default is command mode" << std::endl;
	std::cout << "  -h  Show help message"             << std::endl;
	std::cout << "  -m  Run with macro"                << std::endl;
	std::cout << "  -c  Use specified config"          << std::endl;
	std::cout << "      Default is config.cfg"         << std::endl;
	std::cout << std::endl;
	std::cout << "bye bye :)" << std::endl;
	std::cout << std::endl;
}
