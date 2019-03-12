////////////////////////////////////////////////////////////////////////////////
//   WMC4EveAct.cc                                                            //
//                                                                            //
//   Definitions of WMC4EveAct class's member functions. Details of user      //
// actions are here.                                                          //
//                                                                            //
//                    - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include "WMC4EveAct.hh"
#include "WMC4Ana.hh"
#include "WMC4FPCHit.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4EveAct::WMC4EveAct(WMC4ConMan* CM): m_CM(CM)
{
	// Initialize
	m_FPCThreshold = m_CM -> GetFPCThr();

	m_FPCEvent = new WMC4FPCEve(m_CM);
	m_DetMap   = new WMC4DetMap();

	// Sorted Det ID list of the tubes by position order
	// By doing this, it becomes easy to design clustering algorithm
	m_DetIDByPosOrder.clear();
	for ( G4int i = 0; i < 8; i++ )
	{
		for ( G4int j = 226 + 260*i; j < 336 + 260*i; j++ ) m_DetIDByPosOrder.push_back(j);
		for ( G4int j = 336 + 260*i; j < 346 + 260*i; j++ )
		{
			m_DetIDByPosOrder.push_back(j      );
			m_DetIDByPosOrder.push_back(j - 120);
		}
		for ( G4int j = 346 + 260*i; j < 356 + 260*i; j++ )
		{
			m_DetIDByPosOrder.push_back(j      );
			m_DetIDByPosOrder.push_back(j + 120);
		}
		for ( G4int j = 356 + 260*i; j < 466 + 260*i; j++ ) m_DetIDByPosOrder.push_back(j);
	}
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4EveAct::~WMC4EveAct()
{
	delete m_DetMap;
	delete m_FPCEvent;
}

//////////////////////////////////////////////////
//   Begin of event action                      //
//////////////////////////////////////////////////
void WMC4EveAct::BeginOfEventAction(const G4Event* /* anEvent */)
{
	// Initialize
	m_EDep.clear();
	m_FPCEvent -> Clear();
	// Number of detectors are 24 + 24 + 48 + 24 + 24 + 24 + 24 + 24 + 260*8 = 2296
	for ( int i = 0; i < 2296; i++ ) m_EDep[i] = 0.0 * MeV;
}

//////////////////////////////////////////////////
//   End of event action                        //
//////////////////////////////////////////////////
void WMC4EveAct::EndOfEventAction(const G4Event* anEvent)
{
	// Get analysis manager
	G4AnalysisManager* AM = G4AnalysisManager::Instance();

	// Tracking
	G4int detID, layerID, tubeID;
	std::vector<G4int>::iterator it    = m_DetIDByPosOrder.begin();
	std::vector<G4int>::iterator endIt = m_DetIDByPosOrder.end();
	for ( ; it != endIt; it++ )
	{
		detID = *it;
		if ( m_EDep[detID] < m_FPCThreshold ) continue;

		WMC4FPCHit* hit = new WMC4FPCHit();
		layerID = m_DetMap -> GetFPCLayerIDFromDetID(detID);
		tubeID  = m_DetMap -> GetFPCTubeIDFromDetID(detID);

		hit -> SetLayerID(layerID);
		hit -> SetTubeID(tubeID);
		hit -> SetEDep(m_EDep[detID]);
		if      ( tubeID < - 120 ) hit -> SetPosID(tubeID + 120);
		else if ( tubeID >   120 ) hit -> SetPosID(tubeID - 120);
		else                       hit -> SetPosID(tubeID      );

		m_FPCEvent -> AddHit(hit);

		// Fill EDep hist and count hists
		AM -> FillH1(0, m_EDep[detID]);
		AM -> FillH1(layerID, tubeID);
	}

	// Calculation
	m_FPCEvent -> Clustering();
	if ( m_FPCEvent -> GetNClusters() )
		for ( G4int i = 0; i < m_FPCEvent -> GetNClusters(); i++ )
			AM -> FillH1(17, m_FPCEvent -> GetCluster(i) -> GetClusterSize());

	// Fill CM hists
	for ( G4int i = 0; i < 8; i++ ) AM -> FillH1(i + 9, m_FPCEvent -> GetLayerMultiplicity(i + 1));

	if ( m_FPCEvent -> GetLayerMultiplicity(1) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(2) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(3) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(4) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(5) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(6) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(7) == 1 &&
	     m_FPCEvent -> GetLayerMultiplicity(8) == 1 )
	{
		m_FPCEvent -> CalculateAngle();
	}

	// Fill ntuple
	AM -> FillNtupleIColumn(0, anEvent -> GetEventID() );
	AM -> FillNtupleDColumn(1, m_FPCEvent -> GetTheta());
	AM -> FillNtupleDColumn(2, m_FPCEvent -> GetPhi()  );
	for ( int i = 0; i < 216; i++ )
	{
		if      ( i <  24 && m_CM -> GetFWC1Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i <  48 && m_CM -> GetFWC2Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i <  96 && m_CM -> GetFTHThr()  <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i < 120 && m_CM -> GetFRH1Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i < 144 && m_CM -> GetFRH2Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i < 168 && m_CM -> GetFRH3Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i < 192 && m_CM -> GetFRH4Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else if ( i < 216 && m_CM -> GetFRH5Thr() <= m_EDep[i] ) AM -> FillNtupleDColumn(i+3, m_EDep[i]);
		else                                                     AM -> FillNtupleDColumn(i+3, -1       );
	}
	AM -> AddNtupleRow();
}

//////////////////////////////////////////////////
//   Add energy depositon                       //
//////////////////////////////////////////////////
void WMC4EveAct::AddEDep(G4int id, G4double eDep)
{
	m_EDep[id] += eDep;
}
