////////////////////////////////////////////////////////////////////////////////
//   WMC4SteAct.cc                                                            //
//                                                                            //
//   Definitions of WMC4SteAct class's member functions. Details of user      //
// actions are here.                                                          //
//                                                                            //
//                    - 21. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VProcess.hh"

#include "WMC4SteAct.hh"
#include "WMC4EveAct.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4SteAct::WMC4SteAct(WMC4EveAct* EA): G4UserSteppingAction(), m_EA(EA)
{
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4SteAct::~WMC4SteAct()
{
}

//////////////////////////////////////////////////
//   User stepping action                       //
//////////////////////////////////////////////////
void WMC4SteAct::UserSteppingAction(const G4Step* step)
{
	// I wrote examples of some information which can be extracted from a step.
	// Uncomment whatever you want to use.

	// Who am I? Where am I? What am I undergoing?
	// // Track ID
//	G4int trackID = step -> GetTrack() -> GetTrackID();
	// // Particle name
//	G4String parName = step -> GetTrack() -> GetDefinition() -> GetParticleName();
	// // Particle ID
//	G4int parID = step -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
	// // Particle charge
//	G4int parCharge = step -> GetTrack() -> GetDefinition() -> GetPDGCharge();
	// // Process name
//	const G4VProcess* creProc = step -> GetTrack() -> GetCreatorProcess();
//	G4String procName = step -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();
	// // Physical volume
	G4String namePrePV = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName();
//	G4String namePostPV;
//	G4VPhysicalVolume* postPV = step -> GetPostStepPoint() -> GetPhysicalVolume();
//	if ( postPV != 0 ) namePostPV = postPV -> GetName();
//	else namePostPV = "outside";
	// // Status
//	G4StepStatus preStat = step -> GetPreStepPoint() -> GetStepStatus();
//	G4StepStatus postStat = step -> GetPostStepPoint() -> GetStepStatus();
	// // Position
//	G4ThreeVector prePos = step -> GetPreStepPoint() -> GetPosition();
//	G4ThreeVector postPos = step -> GetPostStepPoint() -> GetPosition();
	// // Momentum
//	G4ThreeVector preMom = step -> GetPreStepPoint() -> GetMomentum();
//	G4ThreeVector postMom = step -> GetPostStepPoint() -> GetMomentum();
	// // Kinetic energy
//	G4double preKinEgy = step -> GetPreStepPoint() -> GetKineticEnergy();
//	G4double postKinEgy = step -> GetPostStepPoint() -> GetKineticEnergy();
	// // Energy deposition
	G4double eDep = step -> GetTotalEnergyDeposit();

	// Send information to the event action object.
	if      ( namePrePV.contains("FWC1PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))-   1, eDep); //    0 ~   23 ( 24)
	else if ( namePrePV.contains("FWC2PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+  23, eDep); //   24 ~   47 ( 24)
	else if ( namePrePV.contains("FTHPV_"      ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  6))+  47, eDep); //   48 ~   95 ( 48)
	else if ( namePrePV.contains("FRH1PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+  95, eDep); //   96 ~  119 ( 24)
	else if ( namePrePV.contains("FRH2PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+ 119, eDep); //  120 ~  143 ( 24)
	else if ( namePrePV.contains("FRH3PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+ 143, eDep); //  144 ~  167 ( 24)
	else if ( namePrePV.contains("FRH4PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+ 167, eDep); //  168 ~  191 ( 24)
	else if ( namePrePV.contains("FRH5PV_"     ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0,  7))+ 191, eDep); //  192 ~  215 ( 24)
	else if ( namePrePV.contains("FPC1TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+ 346, eDep); //  216 ~  345 (130)
	else if ( namePrePV.contains("FPC1TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+ 345, eDep); //  346 ~  475 (130)
	else if ( namePrePV.contains("FPC2TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+ 606, eDep); //  476 ~  605 (130)
	else if ( namePrePV.contains("FPC2TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+ 605, eDep); //  606 ~  735 (130)
	else if ( namePrePV.contains("FPC3TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+ 866, eDep); //  736 ~  865 (130)
	else if ( namePrePV.contains("FPC3TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+ 865, eDep); //  866 ~  995 (130)
	else if ( namePrePV.contains("FPC4TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+1126, eDep); //  996 ~ 1125 (130)
	else if ( namePrePV.contains("FPC4TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+1125, eDep); // 1126 ~ 1255 (130)
	else if ( namePrePV.contains("FPC5TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+1386, eDep); // 1256 ~ 1385 (130)
	else if ( namePrePV.contains("FPC5TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+1385, eDep); // 1386 ~ 1515 (130)
	else if ( namePrePV.contains("FPC6TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+1646, eDep); // 1516 ~ 1645 (130)
	else if ( namePrePV.contains("FPC6TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+1645, eDep); // 1646 ~ 1775 (130)
	else if ( namePrePV.contains("FPC7TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+1906, eDep); // 1776 ~ 1905 (130)
	else if ( namePrePV.contains("FPC7TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+1905, eDep); // 1906 ~ 2035 (130)
	else if ( namePrePV.contains("FPC8TubePV_-") ) m_EA -> AddEDep(-std::stoi(namePrePV.remove(0, 12))+2166, eDep); // 2036 ~ 2165 (130)
	else if ( namePrePV.contains("FPC8TubePV_" ) ) m_EA -> AddEDep( std::stoi(namePrePV.remove(0, 11))+2165, eDep); // 2166 ~ 2295 (130)
}
