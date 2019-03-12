////////////////////////////////////////////////////////////////////////////////
//   WMC4DetMap.cc                                                            //
//                                                                            //
//   Definitions of WMC4DetMap class's member functions                       //
//                                                                            //
//                    - 22. Dev. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4UIsession.hh"

#include "WMC4DetMap.hh"

//////////////////////////////////////////////////
//   Constructors and destructor                //
//////////////////////////////////////////////////
WMC4DetMap::WMC4DetMap()
{
	Init();
}

WMC4DetMap::~WMC4DetMap()
{
}

//////////////////////////////////////////////////
//   Initialization                             //
//////////////////////////////////////////////////
void WMC4DetMap::Init()
{
}

//////////////////////////////////////////////////
//   Get                                        //
//////////////////////////////////////////////////
const G4int WMC4DetMap::GetFPCLayerIDFromDetID(G4int detId)
{
	G4int i;
	i = detId - 216;            // Translation by 216 plastic Scint. elements
	i = i / 260;                // Each plane has 260 tubes

	return i + 1;
}

const G4int WMC4DetMap::GetFPCTubeIDFromDetID(G4int detId)
{
	G4int i;
	i = detId - 216;            // Translation by 216 plastic Scint. elements
	i = i % 260;                // Each plane has 260 tubes
	if ( i < 130 ) i = i - 130; // From -130 to  -1
	else           i = i - 129; // From    1 to 130

//	G4cout << "From " << detId << " to " << i << G4endl;
	return i;
}
