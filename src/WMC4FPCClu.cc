////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCClu.cc                                                            //
//                                                                            //
//   Definitions of WMC4FPCClu class's member functions.                      //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "WMC4FPCClu.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4FPCClu::WMC4FPCClu()
{
	Clear();
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4FPCClu::~WMC4FPCClu()
{
}

//////////////////////////////////////////////////
//   Clear                                      //
//////////////////////////////////////////////////
void WMC4FPCClu::Clear()
{
	// Hit pointers will be deleted by FPC event class.
	// So I don't delete them here...
	m_HitsList.clear();
	m_ClusterSize = 0;
}

//////////////////////////////////////////////////
//   Add FPC hit                                //
//////////////////////////////////////////////////
void WMC4FPCClu::AddHit(WMC4FPCHit* fpcHit)
{
	m_HitsList.push_back(fpcHit);
	m_ClusterSize++;
}

//////////////////////////////////////////////////
//   Number of total hits                       //
//////////////////////////////////////////////////
G4int WMC4FPCClu::GetClusterSize()
{
	return m_ClusterSize;
}

//////////////////////////////////////////////////
//   Layer ID                                   //
//////////////////////////////////////////////////
G4int WMC4FPCClu::GetLayerID()
{
	if ( m_ClusterSize == 0 ) return 0;
	else return m_HitsList[0] -> GetLayerID();
}

//////////////////////////////////////////////////
//   Mean position                              //
//////////////////////////////////////////////////
G4double WMC4FPCClu::GetMean()
{
	if ( m_ClusterSize == 0 ) return 0.;

	G4double mean = 0.;
	std::vector<WMC4FPCHit*>::iterator hitIt    = m_HitsList.begin();
	std::vector<WMC4FPCHit*>::iterator endHitIt = m_HitsList.end();
	for ( ; hitIt != endHitIt; hitIt++ ) mean += 1. * (*hitIt) -> GetPosID() / m_ClusterSize;

	return mean;
}

//////////////////////////////////////////////////
//   Standard deviation                         //
//////////////////////////////////////////////////
G4double WMC4FPCClu::GetDev()
{
	if ( m_ClusterSize == 0 ) return 0.;

	G4double dev = 0.;
	std::vector<WMC4FPCHit*>::iterator hitIt    = m_HitsList.begin();
	std::vector<WMC4FPCHit*>::iterator endHitIt = m_HitsList.end();
	for ( ; hitIt != endHitIt; hitIt++ )
	{
		dev += 1. * ((*hitIt) -> GetPosID() - GetMean())
		          * ((*hitIt) -> GetPosID() - GetMean()) / m_ClusterSize;
	}
	dev = std::pow(dev, 0.5);

	return dev;
}
