#ifndef WMC4EVEACT_h
#define WMC4EVEACT_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4EveAct.hh                                                            //
//                                                                            //
//   This file is a header for WMC4EveAct class. User can add user-defined    //
// event action in this class.                                                //
//                                                                            //
//                    - 22. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <map>

#include "G4UserEventAction.hh"
#include "G4String.hh"

#include "WMC4ConMan.hh"
#include "WMC4FPCEve.hh"
#include "WMC4DetMap.hh"

class G4Event;

class WMC4ConMan;

class WMC4EveAct: public G4UserEventAction
{
  public:
	WMC4EveAct(WMC4ConMan* CM);
	virtual ~WMC4EveAct();

	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);

	void AddEDep(G4int, G4double);

  private:
	WMC4ConMan* m_CM;

	std::map<G4int, G4double> m_EDep;

	// Tracking
	std::vector<G4int> m_DetIDByPosOrder;
	G4double m_FPCThreshold;
	WMC4FPCEve* m_FPCEvent;
	WMC4DetMap* m_DetMap;
};

#endif
