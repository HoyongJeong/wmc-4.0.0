#ifndef WMC4PLUINP_h
#define WMC4PLUINP_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4PluInt.hh                                                            //
//                                                                            //
//   This file is a header for WMC4PluInp class. WMC4PluInp class manages     //
// Pluto input.                                                               //
//                                                                            //
//                    - 10. Dec. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "PParticle.h"

class WMC4PluInp
{
  public:
	WMC4PluInp();
	~WMC4PluInp();

	bool IsGood();
	void Clear();
	void ReadEntry(long iEntry);

	int GetNPar();
	vector<int>    GetID();
	vector<double> GetPX();
	vector<double> GetPY();
	vector<double> GetPZ();
	vector<double> GetTheta();
	vector<double> GetPhi();
	vector<double> GetKE();
	vector<double> GetX();
	vector<double> GetY();
	vector<double> GetZ();

  private:
	TFile* m_File;
	TTree* m_Tree;

	bool m_Good; 
	// Declaration of leave types
	int m_NPar;
	TClonesArray* m_ParArr;
	PParticle* m_Par;

	// Particle properties
	vector<int>    m_ID;
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
