#ifndef WMC4PRIGENACT_h
#define WMC4PRIGENACT_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4PriGenAct.hh                                                         //
//                                                                            //
//   This file is a header for WMC4PriGenAct class. You can set primary       //
// beam options in this class.                                                //
//                                                                            //
//                    - 21. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"

#include "WMC4ConMan.hh"
#include "WMC4PluInp.hh"

class G4ParticleGun;

class WMC4ConMan;

class WMC4PriGenAct: public G4VUserPrimaryGeneratorAction
{
  public:
	WMC4PriGenAct(WMC4ConMan* CM);
	~WMC4PriGenAct();

	virtual void GeneratePrimaries(G4Event* anEvent);

  private:
	WMC4ConMan* m_CM;

	G4ParticleGun*   m_PG;
	G4ParticleTable* m_PT;

	G4double m_BeamPX, m_BeamPY;
	G4double m_BeamDX, m_BeamDY;
	G4double m_WorldZ;
	G4ThreeVector m_GunPos;
	G4ParticleDefinition* m_Par;
	G4ThreeVector m_MomDir;
	G4double m_Mom;
	G4double m_KinEgy;
	G4ThreeVector m_Pol;

	G4bool m_UseInput;
	WMC4PluInp* m_PlutoInput;
	int m_NPar;
	vector<int> m_ID;
	vector<double> m_PX;
	vector<double> m_PY;
	vector<double> m_PZ;
	vector<double> m_Theta;
	vector<double> m_Phi;
	vector<double> m_KE;
	vector<double> m_X;
	vector<double> m_Y;
	vector<double> m_Z;
};

#endif
