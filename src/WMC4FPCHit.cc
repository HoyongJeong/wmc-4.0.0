////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCHit.cc                                                            //
//                                                                            //
//   Definitions of WMC4FPCHit class's member functions.                      //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "WMC4FPCHit.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4FPCHit::WMC4FPCHit()
{
	m_LayerID = 0;
	m_TubeID  = 0;
	m_PosID   = 0;
	m_EDep    = 0;
}

WMC4FPCHit::WMC4FPCHit(G4int layerID, G4int tubeID, G4double eDep)
{
	m_LayerID = layerID;
	m_TubeID  = tubeID;
	m_EDep    = eDep;
	if      ( m_TubeID < - 120 ) m_PosID = m_TubeID + 120;
	else if ( m_TubeID >   120 ) m_PosID = m_TubeID - 120;
	else                         m_PosID = m_TubeID;
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4FPCHit::~WMC4FPCHit()
{
}

//////////////////////////////////////////////////
//   Set & get                                  //
//////////////////////////////////////////////////
void WMC4FPCHit::SetLayerID(G4int layerID)
{
	m_LayerID = layerID;
}
G4int WMC4FPCHit::GetLayerID()
{
	return m_LayerID;
}

void WMC4FPCHit::SetTubeID(G4int tubeID)
{
	m_TubeID = tubeID;
}
G4int WMC4FPCHit::GetTubeID()
{
	return m_TubeID;
}

void WMC4FPCHit::SetPosID(G4int posID)
{
	m_PosID = posID;
}
G4int WMC4FPCHit::GetPosID()
{
	return m_PosID;
}

void WMC4FPCHit::SetEDep(G4double eDep)
{
	m_EDep = eDep;
}
G4double WMC4FPCHit::GetEDep()
{
	return m_EDep;
}
