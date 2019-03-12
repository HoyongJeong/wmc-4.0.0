#ifndef WMC4DETMAP_h
#define WMC4DETMAP_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4DetMap.hh                                                            //
//                                                                            //
//   This file is a header for WMC4DetMap class. I made this class to get     //
// mapping between IDs or names                                               //
//                                                                            //
//                    - 22. Dec. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <map>
#include "G4String.hh"

class WMC4DetMap
{
  public:
	// Constructor and destructor
	WMC4DetMap();
	~WMC4DetMap();

	// Initializaion
	void Init();

	// Get
	const G4int GetFPCLayerIDFromDetID(G4int);
	const G4int GetFPCTubeIDFromDetID(G4int);
};

#endif
