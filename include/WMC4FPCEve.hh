#ifndef WMC4FPCEVE_h
#define WMC4FPCEVE_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4FPCEve.hh                                                            //
//                                                                            //
//   This file is a header for WMC4FPCEve class. The class collects FPC hits, //
// organizes clusters from the hits, and finds tracks.                        //
//                                                                            //
//                    - 22. Nov. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"

#include "WMC4ConMan.hh"
#include "WMC4FPCHit.hh"
#include "WMC4FPCClu.hh"

class WMC4FPCEve
{
  public:
	WMC4FPCEve(WMC4ConMan* CM);
	~WMC4FPCEve();

	void Init();
	void Clear();
	void AddHit(WMC4FPCHit* fpcHit);

	// Calculation
	void Clustering();
	void CalculateAngle();

	// Set & get
	G4int GetNTotalHits();
	G4int GetNClusters();
	WMC4FPCClu* GetCluster(G4int);
	G4int GetLayerMultiplicity(G4int);
	G4bool GetTracked();
	G4double GetTheta();
	G4double GetPhi();

  private:
	WMC4ConMan* m_CM;

	std::vector<WMC4FPCHit*> m_HitsList;
	G4int m_NTotalHits;
	G4int m_NClusters;

	// Clustering
	std::vector<WMC4FPCClu*> m_ClustersList;
	G4bool m_FlagClusterEnd;
	G4int m_CurrentHitLayer;
	G4int m_CurrentHitPos;
	G4int m_LastHitLayer;
	G4int m_LastHitPos;
	G4int m_LayerMultiplicity[8];

	// Tracking info
	G4double m_U1, m_U2, m_V1, m_V2, m_X1, m_X2, m_Y1, m_Y2;
	G4double m_U, m_V, m_X, m_Y;
	G4ThreeVector m_PUV;
	G4ThreeVector m_PXY;
	G4double m_Delta;
	G4double m_x, m_y, m_rho;
	G4double m_Theta, m_Phi;
	G4bool m_Tracked;
};

#endif
