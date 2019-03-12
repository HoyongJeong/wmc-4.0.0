////////////////////////////////////////////////////////////////////////////////
//   WMC4RunAct.cc                                                            //
//                                                                            //
//   Definitions of WMC4RunAct class's member functions. Details of user      //
// actions are here.                                                          //
//                                                                            //
//                    - 21. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <ctime>

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

#include "WMC4RunAct.hh"
#include "WMC4Ana.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4RunAct::WMC4RunAct(WMC4ConMan* CM): G4UserRunAction(), m_CM(CM)
{
	// Create analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in WMC4Ana.hh
	G4AnalysisManager* AM = G4AnalysisManager::Instance();
	G4cout << "Using " << AM -> GetType() << G4endl;

	// Default settings
	// I don't like chatterbox...
	AM -> SetVerboseLevel(0);

	// Directories
	// Creating histograms
	AM -> CreateH1("FPCEDep", "Energy deposition in a tube", 2000, 0., 0.05*MeV); // HID =  0;
	AM -> CreateH1("FPC1", "Hit count on FPC1", 261, -130.5, 130.5);              // HID =  1;
	AM -> CreateH1("FPC2", "Hit count on FPC2", 261, -130.5, 130.5);              // HID =  2;
	AM -> CreateH1("FPC3", "Hit count on FPC3", 261, -130.5, 130.5);              // HID =  3;
	AM -> CreateH1("FPC4", "Hit count on FPC4", 261, -130.5, 130.5);              // HID =  4;
	AM -> CreateH1("FPC5", "Hit count on FPC5", 261, -130.5, 130.5);              // HID =  5;
	AM -> CreateH1("FPC6", "Hit count on FPC6", 261, -130.5, 130.5);              // HID =  6;
	AM -> CreateH1("FPC7", "Hit count on FPC7", 261, -130.5, 130.5);              // HID =  7;
	AM -> CreateH1("FPC8", "Hit count on FPC8", 261, -130.5, 130.5);              // HID =  8;
	AM -> CreateH1("FPC1CM", "Cluster multiplicity of FPC1", 11, -0.5, 10.5);     // HID =  9;
	AM -> CreateH1("FPC2CM", "Cluster multiplicity of FPC2", 11, -0.5, 10.5);     // HID = 10;
	AM -> CreateH1("FPC3CM", "Cluster multiplicity of FPC3", 11, -0.5, 10.5);     // HID = 11;
	AM -> CreateH1("FPC4CM", "Cluster multiplicity of FPC4", 11, -0.5, 10.5);     // HID = 12;
	AM -> CreateH1("FPC5CM", "Cluster multiplicity of FPC5", 11, -0.5, 10.5);     // HID = 13;
	AM -> CreateH1("FPC6CM", "Cluster multiplicity of FPC6", 11, -0.5, 10.5);     // HID = 14;
	AM -> CreateH1("FPC7CM", "Cluster multiplicity of FPC7", 11, -0.5, 10.5);     // HID = 15;
	AM -> CreateH1("FPC8CM", "Cluster multiplicity of FPC8", 11, -0.5, 10.5);     // HID = 16;
	AM -> CreateH1("FPCCS", "Cluster size", 11, -0.5, 10.5);                      // HID = 17;

	// Creating ntuple
	AM -> CreateNtuple("data", "data");
	AM -> CreateNtupleIColumn("event_id"); // Column ID = 0
	AM -> CreateNtupleDColumn("theta");    // Column ID = 1
	AM -> CreateNtupleDColumn("phi");      // Column ID = 2
	G4String columnName;
	for ( int i = 0; i < 216; i++ )        // Column ID = 3 ~ 218
	{
		columnName  = "DetID";
		columnName += G4UIcommand::ConvertToString(i);
		AM -> CreateNtupleDColumn(columnName.data());
	}
	AM -> FinishNtuple();
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4RunAct::~WMC4RunAct()
{
	delete G4AnalysisManager::Instance();
}

//////////////////////////////////////////////////
//   Begin of run action                        //
//////////////////////////////////////////////////
void WMC4RunAct::BeginOfRunAction(const G4Run* /*run*/)
{
	// All actions defined here will be excuted at the beginning of every run.
	// What is a run? You may type "/run/beamOn [someNumber]".
	// Whenever you do this, "one run" runs.

	// Get analysis manager
	G4AnalysisManager* AM = G4AnalysisManager::Instance();

	// Get current time to include it to file name
	// This time info is going to be used to generate output file name.
	time_t rawTime;
	struct tm* timeInfo;
	char buffer[80];
	time(&rawTime);
	timeInfo = localtime(&rawTime);
	strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", timeInfo);
	std::string sTime(buffer);

	// Get kinetic energy to include it to file name
	G4double dKE = m_CM -> GetBeamKE();
	G4int iKE = dKE / MeV;

	// Open an output file
	// The default file name is set in WMC4RunAct::WMC4RunAct(),
	// and it can be overwritten in a macro
	G4String fileName;
	if ( m_CM -> GetUseInput() )
	{
		fileName = "pluto";
	}
	else
	{
		if ( m_CM -> GetUseTar() ) fileName = m_CM -> GetTarMat();
		else fileName = "empty";
		fileName += "-";
		fileName += std::to_string(iKE);
	}
	fileName += "-";
	fileName += sTime;

	AM -> OpenFile(fileName);
}

//////////////////////////////////////////////////
//   End of run action                          //
//////////////////////////////////////////////////
void WMC4RunAct::EndOfRunAction(const G4Run* run)
{
	// save histograms & ntuple
	G4AnalysisManager* AM = G4AnalysisManager::Instance();
	// You must save. Otherwise, file will be just empty...
	AM -> Write();
	// You must close the file. Otherwise, file will be crahsed...
	AM -> CloseFile();
}
