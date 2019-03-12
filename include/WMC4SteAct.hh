#ifndef WMC4STEACT_h
#define WMC4STEACT_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4SetAct.hh                                                            //
//                                                                            //
//   This file is a header for WMC4SteAct class. User can add user-defined    //
// stepping action in this class. So this class works at every step.          //
// The most busiest class.                                                    //
//                                                                            //
//                    - 23. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4UserSteppingAction.hh"
#include "WMC4EveAct.hh"

class WMC4EveAct;

class WMC4SteAct: public G4UserSteppingAction
{
  public:
	WMC4SteAct(WMC4EveAct* EA);
	virtual ~WMC4SteAct();

	virtual void UserSteppingAction(const G4Step*);

  private:
	WMC4EveAct* m_EA;
};

#endif
