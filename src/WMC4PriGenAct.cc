////////////////////////////////////////////////////////////////////////////////
//   WMC4PriGenAct.cc                                                         //
//                                                                            //
//   Definitions of WMC4PriGenAct class's member functions.                   //
//                                                                            //
//                    - 21. Nov. 2017. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////


#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include "WMC4PriGenAct.hh"

//////////////////////////////////////////////////
//   Constructor and destructor                 //
//////////////////////////////////////////////////
WMC4PriGenAct::WMC4PriGenAct(WMC4ConMan* CM): m_CM(CM)
{
	m_PG = new G4ParticleGun();

	// Gun position
	m_WorldZ = m_CM -> GetWorldZ();
	m_BeamPX = m_CM -> GetBeamPX();
	m_BeamPY = m_CM -> GetBeamPY();
	m_BeamDX = m_CM -> GetBeamDX();
	m_BeamDY = m_CM -> GetBeamDY();

	// Set particle definition
	m_PT = G4ParticleTable::GetParticleTable();
	m_Par = m_PT -> FindParticle("proton");
	m_PG -> SetParticleDefinition(m_Par);

	// Momentum
	m_MomDir = G4ThreeVector(0., 0., 1.);
	m_PG -> SetParticleMomentumDirection(m_MomDir);

	// Kinetic energy
	m_KinEgy = m_CM -> GetBeamKE();
	m_PG -> SetParticleEnergy(m_KinEgy);

	// Polarization
	m_Pol = G4ThreeVector(0., 1., 0.);
	m_PG -> SetParticlePolarization(m_Pol);

	// Read Pluto input
	m_UseInput = m_CM -> GetUseInput();
	if ( m_UseInput )
	{
		m_PlutoInput = new WMC4PluInp();
		if ( m_PlutoInput -> IsGood() ) G4cout << "[Notice] pluto.root has been loaded." << G4endl;
		else                            G4cout << "[Error] Cannot open pluto.root!"      << G4endl;
	}
}

WMC4PriGenAct::~WMC4PriGenAct()
{
	delete m_PG;
}

//////////////////////////////////////////////////
//   Shoot!                                     //
//////////////////////////////////////////////////
void WMC4PriGenAct::GeneratePrimaries(G4Event* anEvent)
{
	int EID = anEvent -> GetEventID();

	// Read PLUTO
	if ( m_UseInput && m_PlutoInput -> IsGood() )
	{
		m_PlutoInput -> ReadEntry(EID);
		m_NPar  = m_PlutoInput -> GetNPar();
		m_ID    = m_PlutoInput -> GetID();
		m_PX    = m_PlutoInput -> GetPX();
		m_PY    = m_PlutoInput -> GetPY();
		m_PZ    = m_PlutoInput -> GetPZ();
		m_Theta = m_PlutoInput -> GetTheta();
		m_Phi   = m_PlutoInput -> GetPhi();
		m_KE    = m_PlutoInput -> GetKE();
		m_X     = m_PlutoInput -> GetX();
		m_Y     = m_PlutoInput -> GetY();
		m_Z     = m_PlutoInput -> GetZ();

		// Shoot particles
		for ( int nP = 0; nP < m_NPar; nP++ )
		{
			// Set particle definition
			if      ( m_ID[nP] ==  14 ) m_Par = m_PT -> FindParticle("proton");
//			else if ( m_ID[nP] == 614 ) m_Par = m_PT -> FindParticle("carbon");
			else continue;
			m_PG -> SetParticleDefinition(m_Par);

			// Momentum
			m_MomDir = G4ThreeVector(m_PX[nP], m_PY[nP], m_PZ[nP]);
			m_PG -> SetParticleMomentumDirection(m_MomDir);

//			G4cout << m_Theta[nP] << "\t" << m_Phi[nP] << G4endl;

			// Kinetic energy
			m_PG -> SetParticleEnergy(m_KE[nP] * GeV);

			// Vertex position
			m_GunPos = G4ThreeVector(m_X[nP] * mm, m_Y[nP] * mm, m_Z[nP] * mm);
			m_PG -> SetParticlePosition(m_GunPos);

			m_PG -> GeneratePrimaryVertex(anEvent);
		}
	}
	else
	{
		// 2D Gaussian
		G4double dX = G4RandGauss::shoot(0., m_BeamDX / mm); // 1 sigma
		G4double dY = G4RandGauss::shoot(0., m_BeamDY / mm); // 1 sigma
		m_GunPos = G4ThreeVector(m_BeamPX + dX * mm, m_BeamPY + dY * mm, - m_WorldZ/2);
//		G4cout << m_GunPos << G4endl;
		m_PG -> SetParticlePosition(m_GunPos);

		m_PG -> GeneratePrimaryVertex(anEvent);
	}
}
