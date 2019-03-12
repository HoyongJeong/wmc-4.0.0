#ifndef WMC4RUNACT_h
#define WMC4RUNACT_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4RunAct.hh                                                            //
//                                                                            //
//   This file is a header for WMC4RunAct class. In addition to the basic     //
// actions, user can add additional user-defined actions.                     //
//                                                                            //
//                    - 21. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4UserRunAction.hh"

#include "WMC4ConMan.hh"

class G4Run;

class WMC4ConMan;

class WMC4RunAct: public G4UserRunAction
{
  public:
	WMC4RunAct(WMC4ConMan* CM);
	virtual ~WMC4RunAct();

	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);

  private:
	// Configuration
	WMC4ConMan* m_CM;
};

#endif
