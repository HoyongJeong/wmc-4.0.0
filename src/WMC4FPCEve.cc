////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCEve.cc                                                            //
//                                                                            //
//   Definitions of WMC4FPCEve class's member functions.                      //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"

#include "WMC4FPCEve.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4FPCEve::WMC4FPCEve(WMC4ConMan* CM)
{
	m_CM = CM;
	Init();
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4FPCEve::~WMC4FPCEve()
{
	Clear();
}

//////////////////////////////////////////////////
//   Initialize                                 //
//////////////////////////////////////////////////
void WMC4FPCEve::Init()
{
	Clear();
}

//////////////////////////////////////////////////
//   Clear                                      //
//////////////////////////////////////////////////
void WMC4FPCEve::Clear()
{
	std::vector<WMC4FPCHit*>::iterator hitIt    = m_HitsList.begin();
	std::vector<WMC4FPCHit*>::iterator endHitIt = m_HitsList.end();
	for ( ; hitIt != endHitIt; hitIt++ ) delete *hitIt;
	m_HitsList.clear();
	m_NTotalHits = 0;

	std::vector<WMC4FPCClu*>::iterator cluIt    = m_ClustersList.begin();
	std::vector<WMC4FPCClu*>::iterator endCluIt = m_ClustersList.end();
	for ( ; cluIt != endCluIt; cluIt++ ) delete *cluIt;
	m_ClustersList.clear();
	m_NClusters  = 0;

	for ( G4int i = 0; i < 8; i++ ) m_LayerMultiplicity[i] = 0;

	m_Theta   =     0;
	m_Phi     =     0;
	m_Tracked = false;
}

//////////////////////////////////////////////////
//   Add FPC hit                                //
//////////////////////////////////////////////////
void WMC4FPCEve::AddHit(WMC4FPCHit* fpcHit)
{
	m_HitsList.push_back(fpcHit);
	m_NTotalHits++;
}

//////////////////////////////////////////////////
//   Clustering                                 //
//////////////////////////////////////////////////
void WMC4FPCEve::Clustering()
{
	// Don't have to do it if there is no hit
	if ( m_NTotalHits == 0 ) return;

	// Prepare
	m_FlagClusterEnd = false;
	WMC4FPCClu* cluster;

	// No condition for first hit
	m_LastHitLayer = m_HitsList[0] -> GetLayerID();
	m_LastHitPos   = m_HitsList[0] -> GetPosID();
	cluster = new WMC4FPCClu();
	cluster -> AddHit(m_HitsList[0]);
	for ( G4int i = 1; i < m_NTotalHits; i++ )
	{
		// Where are we now?
		m_CurrentHitLayer = m_HitsList[i] -> GetLayerID();
		m_CurrentHitPos   = m_HitsList[i] -> GetPosID();

		// Check whether cluster is being continued
		if ( (m_CurrentHitLayer != m_LastHitLayer                                   )
		  || (m_CurrentHitPos != m_LastHitPos + 1 && m_CurrentHitPos != m_LastHitPos) )
			m_FlagClusterEnd = true;

		// If cluster ends, add it and initialize
		if ( m_FlagClusterEnd )
		{
			m_ClustersList.push_back(cluster);
			m_NClusters++;
			m_LayerMultiplicity[m_LastHitLayer-1]++;
			m_FlagClusterEnd = false;
			cluster = new WMC4FPCClu();
		}
		
		// The show must go on
		m_LastHitLayer = m_CurrentHitLayer;
		m_LastHitPos  = m_CurrentHitPos;
		cluster -> AddHit(m_HitsList[i]);
	}

	// The last cluster
	m_ClustersList.push_back(cluster);
	m_NClusters++;
	m_LayerMultiplicity[m_LastHitLayer-1]++;
}

//////////////////////////////////////////////////
//   Calculate angle                            //
//////////////////////////////////////////////////
void WMC4FPCEve::CalculateAngle()
{
	std::vector<WMC4FPCClu*>::iterator cluIt    = m_ClustersList.begin();
	std::vector<WMC4FPCClu*>::iterator endCluIt = m_ClustersList.end();
	for ( ; cluIt != endCluIt; cluIt++ )
	{
		if      ( (*cluIt) -> GetLayerID() == 1 )
		{
			m_U1  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 2 )
		{
			m_U2  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 3 )
		{
			m_V1  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 4 )
		{
			m_V2  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 5 )
		{
			m_X1  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 6 )
		{
			m_X2  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 7 )
		{
			m_Y1  = (*cluIt) -> GetMean();
		}
		else if ( (*cluIt) -> GetLayerID() == 8 )
		{
			m_Y2  = (*cluIt) -> GetMean();
		}
	}

	// Real position scaling
	G4double distance = m_CM -> GetTubeDia() + 0.1 * mm;
	m_U1  *= distance * 0.5;
	m_U2  *= distance * 0.5;
	m_V1  *= distance * 0.5;
	m_V2  *= distance * 0.5;
	m_X1  *= distance * 0.5;
	m_X2  *= distance * 0.5;
	m_Y1  *= distance * 0.5;
	m_Y2  *= distance * 0.5;

	// Merge coordinates
	m_U  = 0.5 * ( m_U1 + m_U2 );
	m_V  = 0.5 * ( m_V1 + m_V2 );
	m_X  = 0.5 * ( m_X1 + m_X2 );
	m_Y  = 0.5 * ( m_Y1 + m_Y2 );

	// Projection
	G4double fpc12PZ = m_CM -> GetFPC12PZ();
	G4double fpc34PZ = m_CM -> GetFPC34PZ();
	G4double fpc56PZ = m_CM -> GetFPC56PZ();
	G4double fpc78PZ = m_CM -> GetFPC78PZ();
	m_V  *= fpc12PZ / fpc34PZ;
	m_X  *= fpc12PZ / fpc56PZ;
	m_Y  *= fpc12PZ / fpc78PZ;

	// Position vector from (U, V) and (X, Y)
	m_PUV.setX(m_U);
	m_PUV.setY(m_V);
	m_PUV.setZ(0);
	m_PXY.setX(m_X);
	m_PXY.setY(m_Y);
	m_PXY.setZ(0);
	m_PUV.rotateZ(-135*degree);
	m_PXY.rotateZ(- 90*degree);

	// Distance between two points
	m_Delta = std::sqrt((m_PUV.x()-m_PXY.x())*(m_PUV.x()-m_PXY.x())
	                  + (m_PUV.y()-m_PXY.y())*(m_PUV.y()-m_PXY.y()));
	
	// Do two points match?
	if ( m_Delta <= m_CM -> GetTubeDia() )
	{
		m_Tracked = true;
		m_x     = 0.5 * (m_PUV.x() + m_PXY.x());
		m_y     = 0.5 * (m_PUV.y() + m_PXY.y());
		m_rho   = std::sqrt(m_x*m_x + m_y*m_y);
		m_Theta = std::atan(m_rho / fpc12PZ);
		m_Phi   = std::atan2(m_y, m_x);
	}
}


//////////////////////////////////////////////////
//   Number of total hits                       //
//////////////////////////////////////////////////
G4int WMC4FPCEve::GetNTotalHits()
{
	return m_NTotalHits;
}

//////////////////////////////////////////////////
//   Number of clusters                         //
//////////////////////////////////////////////////
G4int WMC4FPCEve::GetNClusters()
{
	return m_NClusters;
}

//////////////////////////////////////////////////
//   Get cluster                                //
//////////////////////////////////////////////////
WMC4FPCClu* WMC4FPCEve::GetCluster(G4int iClu)
{
	return m_ClustersList[iClu];
}

//////////////////////////////////////////////////
//   Get layer multiplicity                     //
//////////////////////////////////////////////////
G4int WMC4FPCEve::GetLayerMultiplicity(G4int iLayer)
{
	if ( iLayer > 0 && iLayer < 9 ) return m_LayerMultiplicity[iLayer-1];
	else                            return -1;
}

//////////////////////////////////////////////////
//   Is tracked?                                //
//////////////////////////////////////////////////
G4bool WMC4FPCEve::GetTracked()
{
	return m_Tracked;
}

//////////////////////////////////////////////////
//   Theta and phi                              //
//////////////////////////////////////////////////
G4double WMC4FPCEve::GetTheta()
{
	if ( m_Tracked ) return m_Theta;
	else             return -1000;
}
G4double WMC4FPCEve::GetPhi()
{
	if ( m_Tracked ) return m_Phi;
	else             return -1000;
}
