////////////////////////////////////////////////////////////////////////////////
//   WMC4PluInp.cc                                                            //
//                                                                            //
//   Definitions of WMC4PluInp class's member functions.                      //
//                                                                            //
//                    - 10. Dec. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "WMC4PluInp.hh"

//////////////////////////////////////////////////
//   Constructor                                //
//////////////////////////////////////////////////
WMC4PluInp::WMC4PluInp()
{
	m_Good = true;
	m_File = 0;
	m_Tree = 0;
	m_ParArr = new TClonesArray("PParticle", 20);

	m_File = new TFile("pluto.root");

	if ( !m_File -> IsOpen() )
	{
		m_Good = false;
		return;
	}

	m_Tree = (TTree*) m_File -> Get("data");
	if ( !m_Tree )
	{
		m_Good = false;
		return;
	}

	m_Tree -> SetBranchAddress("Npart"    , &m_NPar  );
	m_Tree -> SetBranchAddress("Particles", &m_ParArr);

	// Initialize
	Clear();
}

//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
WMC4PluInp::~WMC4PluInp()
{
	if ( m_File ) delete m_File;
}

//////////////////////////////////////////////////
//   Is file open okay?                         //
//////////////////////////////////////////////////
bool WMC4PluInp::IsGood()
{
	return m_Good;
}

//////////////////////////////////////////////////
//   Clear particle info                        //
//////////////////////////////////////////////////
void WMC4PluInp::Clear()
{
	m_ID.clear();
	m_PX.clear();
	m_PY.clear();
	m_PZ.clear();
	m_KE.clear();
	m_X.clear();
	m_Y.clear();
	m_Z.clear();
	m_Theta.clear();
	m_Phi.clear();
}

//////////////////////////////////////////////////
//   Get entry                                  //
//////////////////////////////////////////////////
void WMC4PluInp::ReadEntry(long iEntry)
{
	if ( !m_Tree ) return;
	m_Tree -> GetEntry(iEntry);

	Clear();
	for ( int i = 0; i < m_NPar; i++ )
	{
		m_Par = (PParticle*) m_ParArr -> At(i);
		m_ID.push_back(m_Par -> ID());
		m_PX.push_back(m_Par -> Vect().X());
		m_PY.push_back(m_Par -> Vect().Y());
		m_PZ.push_back(m_Par -> Vect().Z());
		m_Theta.push_back(m_Par -> Vect().Theta());
		m_Phi.push_back(m_Par -> Vect().Phi());
		m_KE.push_back(m_Par -> KE());
		m_X.push_back(m_Par -> X());
		m_Y.push_back(m_Par -> Y());
		m_Z.push_back(m_Par -> Z());
	}
}

//////////////////////////////////////////////////
//   How many particles in an event             //
//////////////////////////////////////////////////
int WMC4PluInp::GetNPar()
{
	return m_NPar;
}

//////////////////////////////////////////////////
//   Get particle ID                            //
//////////////////////////////////////////////////
vector<int> WMC4PluInp::GetID()
{
	return m_ID;
}

//////////////////////////////////////////////////
//   Get momentum                               //
//////////////////////////////////////////////////
vector<double> WMC4PluInp::GetPX()
{
	return m_PX;
}

vector<double> WMC4PluInp::GetPY()
{
	return m_PY;
}

vector<double> WMC4PluInp::GetPZ()
{
	return m_PZ;
}

vector<double> WMC4PluInp::GetTheta()
{
	return m_Theta;
}

vector<double> WMC4PluInp::GetPhi()
{
	return m_Phi;
}

//////////////////////////////////////////////////
//   Get kinetic energy                         //
//////////////////////////////////////////////////
vector<double> WMC4PluInp::GetKE()
{
	return m_KE;
}

//////////////////////////////////////////////////
//   Get vertex                                 //
//////////////////////////////////////////////////
vector<double> WMC4PluInp::GetX()
{
	return m_X;
}

vector<double> WMC4PluInp::GetY()
{
	return m_Y;
}

vector<double> WMC4PluInp::GetZ()
{
	return m_Y;
}
