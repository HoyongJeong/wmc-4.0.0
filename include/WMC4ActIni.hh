#ifndef WMC4ACTINI_h
#define WMC4ACTINI_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4ActIni.hh                                                            //
//                                                                            //
//   This file is a header for WMC4ActIni class. Every actions are            //
// initialized through this class.                                            //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VUserActionInitialization.hh"
#include "G4String.hh"

#include "WMC4ConMan.hh"

class WMC4ConMan;

class WMC4ActIni: public G4VUserActionInitialization
{
  public:
	WMC4ActIni(WMC4ConMan* CM);
	virtual ~WMC4ActIni();

	virtual void BuildForMaster() const;
	virtual void Build() const;

  private:
	// Configuration manager
	WMC4ConMan* m_CM;
};

#endif
