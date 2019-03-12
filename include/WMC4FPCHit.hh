#ifndef WMC4FPCHIT_h
#define WMC4FPCHIT_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCHit.hh                                                            //
//                                                                            //
//   This file is a header for WMC4FPCHit class. This class stores a FPC hit. //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VHit.hh"

class WMC4FPCHit: public G4VHit
{
  public:
	WMC4FPCHit();
	WMC4FPCHit(G4int, G4int, G4double);
	virtual ~WMC4FPCHit();

	// Set & get
	void     SetLayerID(G4int);
	G4int    GetLayerID();
	void     SetTubeID(G4int);
	G4int    GetTubeID();
	void     SetPosID(G4int);
	G4int    GetPosID();
	void     SetEDep(G4double);
	G4double GetEDep();

  private:
	G4int m_LayerID;
	G4int m_TubeID;
	G4int m_PosID;
	G4double m_EDep;
};

#endif
