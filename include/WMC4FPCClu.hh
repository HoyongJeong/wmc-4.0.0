#ifndef WMC4FPCCLU_h
#define WMC4FPCCLU_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCClu.hh                                                            //
//                                                                            //
//   This file is a header for WMC4FPCClu class. This class stores a FPC      //
// cluster.                                                                   //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "WMC4ConMan.hh"
#include "WMC4FPCHit.hh"

class WMC4FPCClu
{
  public:
	WMC4FPCClu();
	~WMC4FPCClu();

	void Clear();
	void AddHit(WMC4FPCHit* fpcHit);

	// Set & get
	G4int GetClusterSize();
	G4int GetLayerID();
	G4double GetMean();
	G4double GetDev();

  private:
	std::vector<WMC4FPCHit*> m_HitsList;
	G4int m_ClusterSize;
};

#endif
