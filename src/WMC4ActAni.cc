////////////////////////////////////////////////////////////////////////////////
//   WMC4ActIni.cc                                                            //
//                                                                            //
//   Definitions of WMC4ActIni class's member functions.                      //
// All actions must be initialized here in order to use multi thread.         //
//                                                                            //
//                    - 28. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "WMC4ActIni.hh"
#include "WMC4PriGenAct.hh"
#include "WMC4RunAct.hh"
#include "WMC4EveAct.hh"
#include "WMC4SteAct.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4ActIni::WMC4ActIni(WMC4ConMan* CM): G4VUserActionInitialization(), m_CM(CM)
{
	// Pointer of configuration manager class is initialized
	// when this class is constructed. (CM object)
	// This pointer will be delivered to user actions.
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4ActIni::~WMC4ActIni()
{
}

//////////////////////////////////////////////////
//   Build For Master                           //
//////////////////////////////////////////////////
void WMC4ActIni::BuildForMaster() const
{
	// So, this part is for master. This program is possible to do multithread.
	// A thread will care things as a master.

	// Master thread must have this.
	SetUserAction(new WMC4RunAct(m_CM));
}

//////////////////////////////////////////////////
//   Build                                      //
//////////////////////////////////////////////////
void WMC4ActIni::Build() const
{
	// All user actions are here.
	SetUserAction(new WMC4PriGenAct(m_CM));
	SetUserAction(new WMC4RunAct(m_CM));
  
	WMC4EveAct* WMC4EA = new WMC4EveAct(m_CM);
	SetUserAction(WMC4EA);
  
	SetUserAction(new WMC4SteAct(WMC4EA));
}
