////////////////////////////////////////////////////////////////////////////////
//   WMC4ConMan.cc                                                            //
//                                                                            //
//   Definitions of WMC4ConMan class's member functions                       //
//                                                                            //
//                    - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"

#include "WMC4ConMan.hh"

//////////////////////////////////////////////////
//   Constructors and destructor                //
//////////////////////////////////////////////////
WMC4ConMan::WMC4ConMan()
{
	SetDefault();
	Init();
}

WMC4ConMan::WMC4ConMan(const char* fileName)
{
	SetDefault();
	Init(fileName);
}

WMC4ConMan::~WMC4ConMan()
{
}

//////////////////////////////////////////////////
//   Initialization function                    //
//////////////////////////////////////////////////
void WMC4ConMan::Init(const char* fileName)
{
	if ( !fileName ) SetDefault();
	else if ( !Load(fileName) ) SetDefault();
//	SetDefault();
}

void WMC4ConMan::SetDefault()
{
	// When the class fails to find configuration file, this setup is applied.

	// Computing
	m_NofTRD   =           4;
	// World
	m_WorldX   = 6000.0 * mm;
	m_WorldY   = 6000.0 * mm;
	m_WorldZ   = 6000.0 * mm;
	// Input/output
	m_UseInput =       false;
	m_OutForm  =           0;
	// Beam
	m_BeamPX   =    0.0 * mm;
	m_BeamPY   =    0.0 * mm;
	m_BeamDX   =    1.0 * mm;
	m_BeamDY   =    1.0 * mm;
	m_BeamKE   =  232.8 * MeV;
	// Pipe & chamber
	m_PipeDia  =   80.0 * mm;
	m_VacThi   =    1.0 * mm;
	// Target
	m_UseTar   =        true;
	m_TarMat   =         "C";
	m_TarThi   =    1.0 * mm;
	// FWC1
	m_UseFWC1  =        true;
	m_FWC1Mask.assign(24, true);
	m_FWC1DX1  =    9.5 * mm;
	m_FWC1DX2  =   97.1 * mm;
	m_FWC1DY   =  356.3 * mm;
	m_FWC1Thi  =    3.0 * mm;
	m_FWC1RMin =   39.0 * mm;
	m_FWC1RMax =  334.8 * mm;
	m_FWC1PR   =  206.4 * mm;
	m_FWC1PZ   = 1075.0 * mm;
	m_FWC1Thr  =    0.1 * MeV;
	// FWC2
	m_UseFWC2  =        true;
	m_FWC2Mask.assign(24, true);
	m_FWC2DX1  =    9.5 * mm;
	m_FWC2DX2  =   98.2 * mm;
	m_FWC2DY   =  339.6 * mm;
	m_FWC2Thi  =    3.0 * mm;
	m_FWC2RMin =   39.0 * mm;
	m_FWC2RMax =  339.6 * mm;
	m_FWC2PR   =  208.8 * mm;
	m_FWC2PZ   = 1162.0 * mm;
	m_FWC2Thr  =    0.1 * MeV;
	// FTH
	m_UseFTH   =        true;
	m_FTHMask.assign(48, true);
	m_FTHDX1   =    6.3 * mm;
	m_FTHDX2   =   74.7 * mm;
	m_FTHDY    =  531.7 * mm;
	m_FTHThi   =    5.0 * mm;
	m_FTHRMin  =   48.3 * mm;
	m_FTHRMax  =  580.0 * mm;
	m_FTHPR    =  314.2 * mm;
	m_FTHPZ    = 1577.0 * mm;
	m_FTHThr   =    0.1 * MeV;
	// FRH1
	m_UseFRH1  =        true;
	m_FRH1Mask.assign(24, true);
	m_FRH1DX1  =   11.0 * mm;
	m_FRH1DX2  =  155.0 * mm;
	m_FRH1DY   =  548.0 * mm;
	m_FRH1Thi  =  110.0 * mm;
	m_FRH1RMin =   42.0 * mm;
	m_FRH1RMax =  590.0 * mm;
	m_FRH1PR   =  316.0 * mm;
	m_FRH1PZ   = 1766.5 * mm;
	m_FRH1Thr  =    0.1 * MeV;
	// FRH2
	m_UseFRH2  =        true;
	m_FRH2Mask.assign(24, true);
	m_FRH2DX1  =   11.0 * mm;
	m_FRH2DX2  =  171.0 * mm;
	m_FRH2DY   =  608.0 * mm;
	m_FRH2Thi  =  110.0 * mm;
	m_FRH2RMin =   42.0 * mm;
	m_FRH2RMax =  650.0 * mm;
	m_FRH2PR   =  346.0 * mm;
	m_FRH2PZ   = 1901.5 * mm;
	m_FRH2Thr  =    0.1 * MeV;
	// FRH3
	m_UseFRH3  =        true;
	m_FRH3Mask.assign(24, true);
	m_FRH3DX1  =   11.0 * mm;
	m_FRH3DX2  =  186.9 * mm;
	m_FRH3DY   =  668.0 * mm;
	m_FRH3Thi  =  110.0 * mm;
	m_FRH3RMin =   42.0 * mm;
	m_FRH3RMax =  710.0 * mm;
	m_FRH3PR   =  376.0 * mm;
	m_FRH3PZ   = 2036.5 * mm;
	m_FRH3Thr  =    0.1 * MeV;
	// FRH4
	m_UseFRH4  =        true;
	m_FRH4Mask.assign(24, true);
	m_FRH4DX1  =   11.0 * mm;
	m_FRH4DX2  =  208.0 * mm;
	m_FRH4DY   =  788.0 * mm;
	m_FRH4Thi  =  150.0 * mm;
	m_FRH4RMin =   42.0 * mm;
	m_FRH4RMax =  790.0 * mm;
	m_FRH4PR   =  436.0 * mm;
	m_FRH4PZ   = 2194.0 * mm;
	m_FRH4Thr  =    0.1 * MeV;
	// FRH5
	m_UseFRH5  =        true;
	m_FRH5Mask.assign(24, true);
	m_FRH5DX1  =   11.0 * mm;
	m_FRH5DX2  =  208.0 * mm;
	m_FRH5DY   =  788.0 * mm;
	m_FRH5Thi  =  150.0 * mm;
	m_FRH5RMin =   42.0 * mm;
	m_FRH5RMax =  790.0 * mm;
	m_FRH5PR   =  436.0 * mm;
	m_FRH5PZ   = 2374.0 * mm;
	m_FRH5Thr  =    0.1 * MeV;
	// Wrapping
	m_MylThi   =    0.02 * mm;
	m_TedThi   =    0.08 * mm;
	// FPC
	m_FPC12Phi = 225.0 * degree;
	m_FPC34Phi = 315.0 * degree;
	m_FPC56Phi = 270.0 * degree;
	m_FPC78Phi =   0.0 * degree;
	m_FPC12PZ  = 1219.0 * mm;
	m_FPC34PZ  = 1279.0 * mm;
	m_FPC56PZ  = 1339.0 * mm;
	m_FPC78PZ  = 1399.0 * mm;
	m_FPCAlThi =   0.02 * mm;
	m_TubeDia  =    8.0 * mm;
	m_WallThi  =  0.026 * mm;
	m_WireDia  =   0.02 * mm;
	m_TubeL    =  950.0 * mm;
	m_FPCGas   =           0;
	m_FPCThr   =    0.1 * MeV;
	m_FPCVis   =        true;
	// Trigger
	m_Trigger  =          11;
}

//////////////////////////////////////////////////
//   Load configuration file                    //
//////////////////////////////////////////////////
bool WMC4ConMan::Load(const char* fileName)
{
	// Open file
	std::ifstream file(fileName, std::ifstream::in);
	if ( !file.is_open() ) return false;

	// Read line by line
	std::string line;
	while ( std::getline(file, line) )
	{
		if ( StartsWith(line, "#") ) continue;

		if ( StartsWith(line, "NTHREADS") )
		{
			char hfile[1000];
			sscanf(line.data(), "NTHREADS %s", hfile);
			m_SNofTRD = hfile;
			m_NofTRD = std::stoi(m_SNofTRD);
		}
		// Read configurations: World
		if ( StartsWith(line, "WORLDX") )
		{
			char hfile[1000];
			sscanf(line.data(), "WORLDX %s", hfile);
			m_SWorldX = hfile;
			m_WorldX = std::stod(m_SWorldX) * mm;
		}
		if ( StartsWith(line, "WORLDY") )
		{
			char hfile[1000];
			sscanf(line.data(), "WORLDY %s", hfile);
			m_SWorldY = hfile;
			m_WorldY = std::stod(m_SWorldY) * mm;
		}
		if ( StartsWith(line, "WORLDZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "WORLDZ %s", hfile);
			m_SWorldZ = hfile;
			m_WorldZ = std::stod(m_SWorldZ) * mm;
		}
		// Read configurations: Input/output
		if ( StartsWith(line, "DATAINPUT") )
		{
			char hfile[1000];
			sscanf(line.data(), "DATAINPUT %s", hfile);
			m_SUseInput = hfile;
			if      ( m_SUseInput == "ON"  ) m_UseInput = true;
			else if ( m_SUseInput == "OFF" ) m_UseInput = false;
		}
		if ( StartsWith(line, "OUTPUTFORMAT") )
		{
			char hfile[1000];
			sscanf(line.data(), "OUTPUTFORMAT %s", hfile);
			m_SOutForm = hfile;
			m_OutForm = std::stoi(m_SOutForm);
		}
		// Read configurations: Beam
		if ( StartsWith(line, "BEAMPOSITIONX") )
		{
			char hfile[1000];
			sscanf(line.data(), "BEAMPOSITIONX %s", hfile);
			m_SBeamPX = hfile;
			m_BeamPX = std::stod(m_SBeamPX) * mm;
		}
		if ( StartsWith(line, "BEAMPOSITIONY") )
		{
			char hfile[1000];
			sscanf(line.data(), "BEAMPOSITIONY %s", hfile);
			m_SBeamPY = hfile;
			m_BeamPY = std::stod(m_SBeamPY) * mm;
		}
		if ( StartsWith(line, "BEAMDELTAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "BEAMDELTAX %s", hfile);
			m_SBeamDX = hfile;
			m_BeamDX = std::stod(m_SBeamDX) * mm;
		}
		if ( StartsWith(line, "BEAMDELTAY") )
		{
			char hfile[1000];
			sscanf(line.data(), "BEAMDELTAY %s", hfile);
			m_SBeamDY = hfile;
			m_BeamDY = std::stod(m_SBeamDY) * mm;
		}
		if ( StartsWith(line, "BEAMKINETICENERGY") )
		{
			char hfile[1000];
			sscanf(line.data(), "BEAMKINETICENERGY %s", hfile);
			m_SBeamKE = hfile;
			m_BeamKE = std::stod(m_SBeamKE) * MeV;
		}
		// Read configurations: Pipe & chamber
		if ( StartsWith(line, "PIPEDIAMETER") )
		{
			char hfile[1000];
			sscanf(line.data(), "PIPEDIAMETER %s", hfile);
			m_SPipeDia = hfile;
			m_PipeDia = std::stod(m_SPipeDia) * mm;
		}
		if ( StartsWith(line, "VACUUMWALLTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "VACUUMWALLTHICKNESS %s", hfile);
			m_SVacThi = hfile;
			m_VacThi = std::stod(m_SVacThi) * mm;
		}
		// Read configurations: Target
		if ( StartsWith(line, "USETARGET") )
		{
			char hfile[1000];
			sscanf(line.data(), "USETARGET %s", hfile);
			m_SUseTar = hfile;
			if      ( m_SUseTar == "ON"  ) m_UseTar = true;
			else if ( m_SUseTar == "OFF" ) m_UseTar = false;
		}
		if ( StartsWith(line, "TARGETMATERIAL") )
		{
			char hfile[1000];
			sscanf(line.data(), "TARGETMATERIAL %s", hfile);
			m_STarMat = hfile;
			m_TarMat = m_STarMat;
		}
		if ( StartsWith(line, "TARGETTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "TARGETTHICKNESS %s", hfile);
			m_STarThi = hfile;
			m_TarThi = std::stod(m_STarThi) * mm;
		}
		// Read configurations: FWC1
		if ( StartsWith(line, "USEFWC1") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFWC1 %s", hfile);
			m_SUseFWC1 = hfile;
			if      ( m_SUseFWC1 == "ON"  ) m_UseFWC1 = true;
			else if ( m_SUseFWC1 == "OFF" ) m_UseFWC1 = false;
		}
		if ( StartsWith(line, "FWC1ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1ELEMENTMASK %s", hfile);
			m_SFWC1Mask = hfile;
			if ( m_SFWC1Mask.length() != 24 ) m_FWC1Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFWC1Mask.data()[i];
					if ( temp == 48 ) m_FWC1Mask[i] = false;
					else m_FWC1Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FWC1DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1DX1 %s", hfile);
			m_SFWC1DX1 = hfile;
			m_FWC1DX1 = std::stod(m_SFWC1DX1) * mm;
		}
		if ( StartsWith(line, "FWC1DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1DX2 %s", hfile);
			m_SFWC1DX2 = hfile;
			m_FWC1DX2 = std::stod(m_SFWC1DX2) * mm;
		}
		if ( StartsWith(line, "FWC1DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1DY %s", hfile);
			m_SFWC1DY = hfile;
			m_FWC1DY = std::stod(m_SFWC1DY) * mm;
		}
		if ( StartsWith(line, "FWC1THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1THICKNESS %s", hfile);
			m_SFWC1Thi = hfile;
			m_FWC1Thi = std::stod(m_SFWC1Thi) * mm;
		}
		if ( StartsWith(line, "FWC1RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1RMIN %s", hfile);
			m_SFWC1RMin = hfile;
			m_FWC1RMin = std::stod(m_SFWC1RMin) * mm;
		}
		if ( StartsWith(line, "FWC1RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1RMAX %s", hfile);
			m_SFWC1RMax = hfile;
			m_FWC1RMax = std::stod(m_SFWC1RMax) * mm;
		}
		if ( StartsWith(line, "FWC1POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1POSITIONR %s", hfile);
			m_SFWC1PR = hfile;
			m_FWC1PR = std::stod(m_SFWC1PR) * mm;
		}
		if ( StartsWith(line, "FWC1POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1POSITIONZ %s", hfile);
			m_SFWC1PZ = hfile;
			m_FWC1PZ = std::stod(m_SFWC1PZ) * mm;
		}
		if ( StartsWith(line, "FWC1THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC1THRESHOLD %s", hfile);
			m_SFWC1Thr = hfile;
			m_FWC1Thr = std::stod(m_SFWC1Thr) * MeV;
		}
		// Read configurations: FWC2
		if ( StartsWith(line, "USEFWC2") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFWC2 %s", hfile);
			m_SUseFWC2 = hfile;
			if      ( m_SUseFWC2 == "ON"  ) m_UseFWC2 = true;
			else if ( m_SUseFWC2 == "OFF" ) m_UseFWC2 = false;
		}
		if ( StartsWith(line, "FWC2ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2ELEMENTMASK %s", hfile);
			m_SFWC2Mask = hfile;
			if ( m_SFWC2Mask.length() != 24 ) m_FWC2Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFWC2Mask.data()[i];
					if ( temp == 48 ) m_FWC2Mask[i] = false;
					else m_FWC2Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FWC2DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2DX1 %s", hfile);
			m_SFWC2DX1 = hfile;
			m_FWC2DX1 = std::stod(m_SFWC2DX1) * mm;
		}
		if ( StartsWith(line, "FWC2DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2DX2 %s", hfile);
			m_SFWC2DX2 = hfile;
			m_FWC2DX2 = std::stod(m_SFWC2DX2) * mm;
		}
		if ( StartsWith(line, "FWC2DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2DY %s", hfile);
			m_SFWC2DY = hfile;
			m_FWC2DY = std::stod(m_SFWC2DY) * mm;
		}
		if ( StartsWith(line, "FWC2THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2THICKNESS %s", hfile);
			m_SFWC2Thi = hfile;
			m_FWC2Thi = std::stod(m_SFWC2Thi) * mm;
		}
		if ( StartsWith(line, "FWC2RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2RMIN %s", hfile);
			m_SFWC2RMin = hfile;
			m_FWC2RMin = std::stod(m_SFWC2RMin) * mm;
		}
		if ( StartsWith(line, "FWC2RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2RMAX %s", hfile);
			m_SFWC2RMax = hfile;
			m_FWC2RMax = std::stod(m_SFWC2RMax) * mm;
		}
		if ( StartsWith(line, "FWC2POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2POSITIONR %s", hfile);
			m_SFWC2PR = hfile;
			m_FWC2PR = std::stod(m_SFWC2PR) * mm;
		}
		if ( StartsWith(line, "FWC2POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2POSITIONZ %s", hfile);
			m_SFWC2PZ = hfile;
			m_FWC2PZ = std::stod(m_SFWC2PZ) * mm;
		}
		if ( StartsWith(line, "FWC2THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FWC2THRESHOLD %s", hfile);
			m_SFWC2Thr = hfile;
			m_FWC2Thr = std::stod(m_SFWC2Thr) * MeV;
		}
		// Read configurations: FTH
		if ( StartsWith(line, "USEFTH") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFTH %s", hfile);
			m_SUseFTH = hfile;
			if      ( m_SUseFTH == "ON"  ) m_UseFTH = true;
			else if ( m_SUseFTH == "OFF" ) m_UseFTH = false;
		}
		if ( StartsWith(line, "FTHELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHELEMENTMASK %s", hfile);
			m_SFTHMask = hfile;
			if ( m_SFTHMask.length() != 48 ) m_FTHMask.assign(48, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFTHMask.data()[i];
					if ( temp == 48 ) m_FTHMask[i] = false;
					else m_FTHMask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FTHDX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHDX1 %s", hfile);
			m_SFTHDX1 = hfile;
			m_FTHDX1 = std::stod(m_SFTHDX1) * mm;
		}
		if ( StartsWith(line, "FTHDX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHDX2 %s", hfile);
			m_SFTHDX2 = hfile;
			m_FTHDX2 = std::stod(m_SFTHDX2) * mm;
		}
		if ( StartsWith(line, "FTHDY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHDY %s", hfile);
			m_SFTHDY = hfile;
			m_FTHDY = std::stod(m_SFTHDY) * mm;
		}
		if ( StartsWith(line, "FTHTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHTHICKNESS %s", hfile);
			m_SFTHThi = hfile;
			m_FTHThi = std::stod(m_SFTHThi) * mm;
		}
		if ( StartsWith(line, "FTHRMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHRMIN %s", hfile);
			m_SFTHRMin = hfile;
			m_FTHRMin = std::stod(m_SFTHRMin) * mm;
		}
		if ( StartsWith(line, "FTHRMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHRMAX %s", hfile);
			m_SFTHRMax = hfile;
			m_FTHRMax = std::stod(m_SFTHRMax) * mm;
		}
		if ( StartsWith(line, "FTHPOSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHPOSITIONR %s", hfile);
			m_SFTHPR = hfile;
			m_FTHPR = std::stod(m_SFTHPR) * mm;
		}
		if ( StartsWith(line, "FTHPOSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHPOSITIONZ %s", hfile);
			m_SFTHPZ = hfile;
			m_FTHPZ = std::stod(m_SFTHPZ) * mm;
		}
		if ( StartsWith(line, "FTHTHRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FTHTHRESHOLD %s", hfile);
			m_SFTHThr = hfile;
			m_FTHThr = std::stod(m_SFTHThr) * MeV;
		}
		// Read configurations: FRH1
		if ( StartsWith(line, "USEFRH1") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFRH1 %s", hfile);
			m_SUseFRH1 = hfile;
			if      ( m_SUseFRH1 == "ON"  ) m_UseFRH1 = true;
			else if ( m_SUseFRH1 == "OFF" ) m_UseFRH1 = false;
		}
		if ( StartsWith(line, "FRH1ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1ELEMENTMASK %s", hfile);
			m_SFRH1Mask = hfile;
			if ( m_SFRH1Mask.length() != 24 ) m_FRH1Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFRH1Mask.data()[i];
					if ( temp == 48 ) m_FRH1Mask[i] = false;
					else m_FRH1Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FRH1DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1DX1 %s", hfile);
			m_SFRH1DX1 = hfile;
			m_FRH1DX1 = std::stod(m_SFRH1DX1) * mm;
		}
		if ( StartsWith(line, "FRH1DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1DX2 %s", hfile);
			m_SFRH1DX2 = hfile;
			m_FRH1DX2 = std::stod(m_SFRH1DX2) * mm;
		}
		if ( StartsWith(line, "FRH1DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1DY %s", hfile);
			m_SFRH1DY = hfile;
			m_FRH1DY = std::stod(m_SFRH1DY) * mm;
		}
		if ( StartsWith(line, "FRH1THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1THICKNESS %s", hfile);
			m_SFRH1Thi = hfile;
			m_FRH1Thi = std::stod(m_SFRH1Thi) * mm;
		}
		if ( StartsWith(line, "FRH1RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1RMIN %s", hfile);
			m_SFRH1RMin = hfile;
			m_FRH1RMin = std::stod(m_SFRH1RMin) * mm;
		}
		if ( StartsWith(line, "FRH1RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1RMAX %s", hfile);
			m_SFRH1RMax = hfile;
			m_FRH1RMax = std::stod(m_SFRH1RMax) * mm;
		}
		if ( StartsWith(line, "FRH1POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1POSITIONR %s", hfile);
			m_SFRH1PR = hfile;
			m_FRH1PR = std::stod(m_SFRH1PR) * mm;
		}
		if ( StartsWith(line, "FRH1POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1POSITIONZ %s", hfile);
			m_SFRH1PZ = hfile;
			m_FRH1PZ = std::stod(m_SFRH1PZ) * mm;
		}
		if ( StartsWith(line, "FRH1THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH1THRESHOLD %s", hfile);
			m_SFRH1Thr = hfile;
			m_FRH1Thr = std::stod(m_SFRH1Thr) * MeV;
		}
		// Read configurations: FRH2
		if ( StartsWith(line, "USEFRH2") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFRH2 %s", hfile);
			m_SUseFRH2 = hfile;
			if      ( m_SUseFRH2 == "ON"  ) m_UseFRH2 = true;
			else if ( m_SUseFRH2 == "OFF" ) m_UseFRH2 = false;
		}
		if ( StartsWith(line, "FRH2ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2ELEMENTMASK %s", hfile);
			m_SFRH2Mask = hfile;
			if ( m_SFRH2Mask.length() != 24 ) m_FRH2Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFRH2Mask.data()[i];
					if ( temp == 48 ) m_FRH2Mask[i] = false;
					else m_FRH2Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FRH2DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2DX1 %s", hfile);
			m_SFRH2DX1 = hfile;
			m_FRH2DX1 = std::stod(m_SFRH2DX1) * mm;
		}
		if ( StartsWith(line, "FRH2DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2DX2 %s", hfile);
			m_SFRH2DX2 = hfile;
			m_FRH2DX2 = std::stod(m_SFRH2DX2) * mm;
		}
		if ( StartsWith(line, "FRH2DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2DY %s", hfile);
			m_SFRH2DY = hfile;
			m_FRH2DY = std::stod(m_SFRH2DY) * mm;
		}
		if ( StartsWith(line, "FRH2THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2THICKNESS %s", hfile);
			m_SFRH2Thi = hfile;
			m_FRH2Thi = std::stod(m_SFRH2Thi) * mm;
		}
		if ( StartsWith(line, "FRH2RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2RMIN %s", hfile);
			m_SFRH2RMin = hfile;
			m_FRH2RMin = std::stod(m_SFRH2RMin) * mm;
		}
		if ( StartsWith(line, "FRH2RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2RMAX %s", hfile);
			m_SFRH2RMax = hfile;
			m_FRH2RMax = std::stod(m_SFRH2RMax) * mm;
		}
		if ( StartsWith(line, "FRH2POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2POSITIONR %s", hfile);
			m_SFRH2PR = hfile;
			m_FRH2PR = std::stod(m_SFRH2PR) * mm;
		}
		if ( StartsWith(line, "FRH2POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2POSITIONZ %s", hfile);
			m_SFRH2PZ = hfile;
			m_FRH2PZ = std::stod(m_SFRH2PZ) * mm;
		}
		if ( StartsWith(line, "FRH2THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH2THRESHOLD %s", hfile);
			m_SFRH2Thr = hfile;
			m_FRH2Thr = std::stod(m_SFRH2Thr) * MeV;
		}
		// Read configurations: FRH3
		if ( StartsWith(line, "USEFRH3") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFRH3 %s", hfile);
			m_SUseFRH3 = hfile;
			if      ( m_SUseFRH3 == "ON"  ) m_UseFRH3 = true;
			else if ( m_SUseFRH3 == "OFF" ) m_UseFRH3 = false;
		}
		if ( StartsWith(line, "FRH3ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3ELEMENTMASK %s", hfile);
			m_SFRH3Mask = hfile;
			if ( m_SFRH3Mask.length() != 24 ) m_FRH3Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFRH3Mask.data()[i];
					if ( temp == 48 ) m_FRH3Mask[i] = false;
					else m_FRH3Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FRH3DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3DX1 %s", hfile);
			m_SFRH3DX1 = hfile;
			m_FRH3DX1 = std::stod(m_SFRH3DX1) * mm;
		}
		if ( StartsWith(line, "FRH3DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3DX2 %s", hfile);
			m_SFRH3DX2 = hfile;
			m_FRH3DX2 = std::stod(m_SFRH3DX2) * mm;
		}
		if ( StartsWith(line, "FRH3DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3DY %s", hfile);
			m_SFRH3DY = hfile;
			m_FRH3DY = std::stod(m_SFRH3DY) * mm;
		}
		if ( StartsWith(line, "FRH3THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3THICKNESS %s", hfile);
			m_SFRH3Thi = hfile;
			m_FRH3Thi = std::stod(m_SFRH3Thi) * mm;
		}
		if ( StartsWith(line, "FRH3RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3RMIN %s", hfile);
			m_SFRH3RMin = hfile;
			m_FRH3RMin = std::stod(m_SFRH3RMin) * mm;
		}
		if ( StartsWith(line, "FRH3RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3RMAX %s", hfile);
			m_SFRH3RMax = hfile;
			m_FRH3RMax = std::stod(m_SFRH3RMax) * mm;
		}
		if ( StartsWith(line, "FRH3POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3POSITIONR %s", hfile);
			m_SFRH3PR = hfile;
			m_FRH3PR = std::stod(m_SFRH3PR) * mm;
		}
		if ( StartsWith(line, "FRH3POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3POSITIONZ %s", hfile);
			m_SFRH3PZ = hfile;
			m_FRH3PZ = std::stod(m_SFRH3PZ) * mm;
		}
		if ( StartsWith(line, "FRH3THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH3THRESHOLD %s", hfile);
			m_SFRH3Thr = hfile;
			m_FRH3Thr = std::stod(m_SFRH3Thr) * MeV;
		}
		// Read configurations: FRH4
		if ( StartsWith(line, "USEFRH4") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFRH4 %s", hfile);
			m_SUseFRH4 = hfile;
			if      ( m_SUseFRH4 == "ON"  ) m_UseFRH4 = true;
			else if ( m_SUseFRH4 == "OFF" ) m_UseFRH4 = false;
		}
		if ( StartsWith(line, "FRH4ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4ELEMENTMASK %s", hfile);
			m_SFRH4Mask = hfile;
			if ( m_SFRH4Mask.length() != 24 ) m_FRH4Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFRH4Mask.data()[i];
					if ( temp == 48 ) m_FRH4Mask[i] = false;
					else m_FRH4Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FRH4DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4DX1 %s", hfile);
			m_SFRH4DX1 = hfile;
			m_FRH4DX1 = std::stod(m_SFRH4DX1) * mm;
		}
		if ( StartsWith(line, "FRH4DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4DX2 %s", hfile);
			m_SFRH4DX2 = hfile;
			m_FRH4DX2 = std::stod(m_SFRH4DX2) * mm;
		}
		if ( StartsWith(line, "FRH4DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4DY %s", hfile);
			m_SFRH4DY = hfile;
			m_FRH4DY = std::stod(m_SFRH4DY) * mm;
		}
		if ( StartsWith(line, "FRH4THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4THICKNESS %s", hfile);
			m_SFRH4Thi = hfile;
			m_FRH4Thi = std::stod(m_SFRH4Thi) * mm;
		}
		if ( StartsWith(line, "FRH4RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4RMIN %s", hfile);
			m_SFRH4RMin = hfile;
			m_FRH4RMin = std::stod(m_SFRH4RMin) * mm;
		}
		if ( StartsWith(line, "FRH4RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4RMAX %s", hfile);
			m_SFRH4RMax = hfile;
			m_FRH4RMax = std::stod(m_SFRH4RMax) * mm;
		}
		if ( StartsWith(line, "FRH4POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4POSITIONR %s", hfile);
			m_SFRH4PR = hfile;
			m_FRH4PR = std::stod(m_SFRH4PR) * mm;
		}
		if ( StartsWith(line, "FRH4POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4POSITIONZ %s", hfile);
			m_SFRH4PZ = hfile;
			m_FRH4PZ = std::stod(m_SFRH4PZ) * mm;
		}
		if ( StartsWith(line, "FRH4THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH4THRESHOLD %s", hfile);
			m_SFRH4Thr = hfile;
			m_FRH4Thr = std::stod(m_SFRH4Thr) * MeV;
		}
		// Read configurations: FRH5
		if ( StartsWith(line, "USEFRH5") )
		{
			char hfile[1000];
			sscanf(line.data(), "USEFRH5 %s", hfile);
			m_SUseFRH5 = hfile;
			if      ( m_SUseFRH5 == "ON"  ) m_UseFRH5 = true;
			else if ( m_SUseFRH5 == "OFF" ) m_UseFRH5 = false;
		}
		if ( StartsWith(line, "FRH5ELEMENTMASK") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5ELEMENTMASK %s", hfile);
			m_SFRH5Mask = hfile;
			if ( m_SFRH5Mask.length() != 24 ) m_FRH5Mask.assign(24, false);
			else
			{
				int temp;
				for ( int i = 0; i < 24; i++ )
				{
					temp = m_SFRH5Mask.data()[i];
					if ( temp == 48 ) m_FRH5Mask[i] = false;
					else m_FRH5Mask[i] = true;
				}
			}
		}
		if ( StartsWith(line, "FRH5DX1") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5DX1 %s", hfile);
			m_SFRH5DX1 = hfile;
			m_FRH5DX1 = std::stod(m_SFRH5DX1) * mm;
		}
		if ( StartsWith(line, "FRH5DX2") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5DX2 %s", hfile);
			m_SFRH5DX2 = hfile;
			m_FRH5DX2 = std::stod(m_SFRH5DX2) * mm;
		}
		if ( StartsWith(line, "FRH5DY") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5DY %s", hfile);
			m_SFRH5DY = hfile;
			m_FRH5DY = std::stod(m_SFRH5DY) * mm;
		}
		if ( StartsWith(line, "FRH5THICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5THICKNESS %s", hfile);
			m_SFRH5Thi = hfile;
			m_FRH5Thi = std::stod(m_SFRH5Thi) * mm;
		}
		if ( StartsWith(line, "FRH5RMIN") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5RMIN %s", hfile);
			m_SFRH5RMin = hfile;
			m_FRH5RMin = std::stod(m_SFRH5RMin) * mm;
		}
		if ( StartsWith(line, "FRH5RMAX") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5RMAX %s", hfile);
			m_SFRH5RMax = hfile;
			m_FRH5RMax = std::stod(m_SFRH5RMax) * mm;
		}
		if ( StartsWith(line, "FRH5POSITIONR") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5POSITIONR %s", hfile);
			m_SFRH5PR = hfile;
			m_FRH5PR = std::stod(m_SFRH5PR) * mm;
		}
		if ( StartsWith(line, "FRH5POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5POSITIONZ %s", hfile);
			m_SFRH5PZ = hfile;
			m_FRH5PZ = std::stod(m_SFRH5PZ) * mm;
		}
		if ( StartsWith(line, "FRH5THRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FRH5THRESHOLD %s", hfile);
			m_SFRH5Thr = hfile;
			m_FRH5Thr = std::stod(m_SFRH5Thr) * MeV;
		}
		// Read configuration: Wrapping
		if ( StartsWith(line, "MYLARTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "MYLARTHICKNESS %s", hfile);
			m_SMylThi = hfile;
			m_MylThi = std::stod(m_SMylThi) * mm;
		}
		if ( StartsWith(line, "TEDLARTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "TEDLARTHICKNESS %s", hfile);
			m_STedThi = hfile;
			m_TedThi = std::stod(m_STedThi) * mm;
		}
		// Read configurations: FPC
		if ( StartsWith(line, "FPC12PHI") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC12PHI %s", hfile);
			m_SFPC12Phi = hfile;
			m_FPC12Phi = std::stod(m_SFPC12Phi) * degree;
		}
		if ( StartsWith(line, "FPC34PHI") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC34PHI %s", hfile);
			m_SFPC34Phi = hfile;
			m_FPC34Phi = std::stod(m_SFPC34Phi) * degree;
		}
		if ( StartsWith(line, "FPC56PHI") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC56PHI %s", hfile);
			m_SFPC56Phi = hfile;
			m_FPC56Phi = std::stod(m_SFPC56Phi) * degree;
		}
		if ( StartsWith(line, "FPC78PHI") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC78PHI %s", hfile);
			m_SFPC78Phi = hfile;
			m_FPC78Phi = std::stod(m_SFPC78Phi) * degree;
		}
		if ( StartsWith(line, "FPC12POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC12POSITIONZ %s", hfile);
			m_SFPC12PZ = hfile;
			m_FPC12PZ = std::stod(m_SFPC12PZ) * mm;
		}
		if ( StartsWith(line, "FPC34POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC34POSITIONZ %s", hfile);
			m_SFPC34PZ = hfile;
			m_FPC34PZ = std::stod(m_SFPC34PZ) * mm;
		}
		if ( StartsWith(line, "FPC56POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC56POSITIONZ %s", hfile);
			m_SFPC56PZ = hfile;
			m_FPC56PZ = std::stod(m_SFPC56PZ) * mm;
		}
		if ( StartsWith(line, "FPC78POSITIONZ") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPC78POSITIONZ %s", hfile);
			m_SFPC78PZ = hfile;
			m_FPC78PZ = std::stod(m_SFPC78PZ) * mm;
		}
		if ( StartsWith(line, "FPCSHIELDTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPCSHIELDTHICKNESS %s", hfile);
			m_SFPCAlThi = hfile;
			m_FPCAlThi = std::stod(m_SFPCAlThi) * mm;
		}
		if ( StartsWith(line, "TUBEDIAMETER") )
		{
			char hfile[1000];
			sscanf(line.data(), "TUBEDIAMETER %s", hfile);
			m_STubeDia = hfile;
			m_TubeDia = std::stod(m_STubeDia) * mm;
		}
		if ( StartsWith(line, "MYLARWALLTHICKNESS") )
		{
			char hfile[1000];
			sscanf(line.data(), "MYLARWALLTHICKNESS %s", hfile);
			m_SWallThi = hfile;
			m_WallThi = std::stod(m_SWallThi) * mm;
		}
		if ( StartsWith(line, "WIREDIAMETER") )
		{
			char hfile[1000];
			sscanf(line.data(), "WIREDIAMETER %s", hfile);
			m_SWireDia = hfile;
			m_WireDia = std::stod(m_SWireDia) * mm;
		}
		if ( StartsWith(line, "TUBELENGTH") )
		{
			char hfile[1000];
			sscanf(line.data(), "TUBELENGTH %s", hfile);
			m_STubeL = hfile;
			m_TubeL = std::stod(m_STubeL) * mm;
		}
		if ( StartsWith(line, "FPCGAS") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPCGAS %s", hfile);
			m_SFPCGas = hfile;
			m_FPCGas = std::stoi(m_SFPCGas);
		}
		if ( StartsWith(line, "FPCTHRESHOLD") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPCTHRESHOLD %s", hfile);
			m_SFPCThr = hfile;
			m_FPCThr = std::stod(m_SFPCThr) * MeV;
		}
		if ( StartsWith(line, "FPCVISUALIZE") )
		{
			char hfile[1000];
			sscanf(line.data(), "FPCVISUALIZE %s", hfile);
			m_SFPCVis = hfile;
			if      ( m_SFPCVis == "ON"  ) m_FPCVis = true;
			else if ( m_SFPCVis == "OFF" ) m_FPCVis = false;
		}
		// Read configurations: Trigger
		if ( StartsWith(line, "TRIGGER") )
		{
			char hfile[1000];
			sscanf(line.data(), "TRIGGER %s", hfile);
			m_STrigger = hfile;
			m_Trigger = std::stoi(m_STrigger);
		}
	}

	PrintConfiguration();

	return true;
}

void WMC4ConMan::PrintConfiguration()
{
	// Computing
	printf(" WMC4ConMan::Load() => NTHREADS            %s\n"    , m_SNofTRD.data()  );
	// World                                             
	printf(" WMC4ConMan::Load() => WORLDX              %s mm\n" , m_SWorldX.data()  );
	printf(" WMC4ConMan::Load() => WORLDY              %s mm\n" , m_SWorldY.data()  );
	printf(" WMC4ConMan::Load() => WORLDZ              %s mm\n" , m_SWorldZ.data()  );
	// Input/output
	printf(" WMC4ConMan::Load() => DATAINPUT           %s\n"    , m_SUseInput.data());
	printf(" WMC4ConMan::Load() => OUTPUTFORMAT        %s\n"    , m_SOutForm.data() );
	// Beam                                              
	printf(" WMC4ConMan::Load() => BEAMPOSITIONX       %s mm\n" , m_SBeamPX.data()  );
	printf(" WMC4ConMan::Load() => BEAMPOSITIONY       %s mm\n" , m_SBeamPY.data()  );
	printf(" WMC4ConMan::Load() => BEAMDELTAX          %s mm\n" , m_SBeamDX.data()  );
	printf(" WMC4ConMan::Load() => BEAMDELTAY          %s mm\n" , m_SBeamDY.data()  );
	printf(" WMC4ConMan::Load() => BEAMKINETICENERGY   %s MeV\n", m_SBeamKE.data()  );
	// Pipe & Chamber
	printf(" WMC4ConMan::Load() => PIPEDIAMETER        %s mm\n" , m_SPipeDia.data() );
	printf(" WMC4ConMan::Load() => VACUUMWALLTHICKNESS %s mm\n" , m_SVacThi.data()  );
	// Target                                            
	printf(" WMC4ConMan::Load() => USETARGET           %s\n"    , m_SUseTar.data()  );
	printf(" WMC4ConMan::Load() => TARGETMATERIAL      %s\n"    , m_STarMat.data()  );
	printf(" WMC4ConMan::Load() => TARGETTHICKNESS     %s mm\n" , m_STarThi.data()  );
	// FWC1
	printf(" WMC4ConMan::Load() => USEFWC1             %s\n"    , m_SUseFWC1.data() );
	printf(" WMC4ConMan::Load() => FWC1ELEMENTMASK     %s\n"    , m_SFWC1Mask.data());
	printf(" WMC4ConMan::Load() => FWC1DX1             %s mm\n" , m_SFWC1DX1.data() );
	printf(" WMC4ConMan::Load() => FWC1DX2             %s mm\n" , m_SFWC1DX2.data() );
	printf(" WMC4ConMan::Load() => FWC1DY              %s mm\n" , m_SFWC1DY.data()  );
	printf(" WMC4ConMan::Load() => FWC1THICKNESS       %s mm\n" , m_SFWC1Thi.data() );
	printf(" WMC4ConMan::Load() => FWC1RMIN            %s mm\n" , m_SFWC1RMin.data());
	printf(" WMC4ConMan::Load() => FWC1RMAX            %s mm\n" , m_SFWC1RMax.data());
	printf(" WMC4ConMan::Load() => FWC1POSITIONR       %s mm\n" , m_SFWC1PR.data()  );
	printf(" WMC4ConMan::Load() => FWC1POSITIONZ       %s mm\n" , m_SFWC1PZ.data()  );
	printf(" WMC4ConMan::Load() => FWC1THRESHOLD       %s MeV\n", m_SFWC1Thr.data() );
	// FWC2
	printf(" WMC4ConMan::Load() => USEFWC2             %s\n"    , m_SUseFWC2.data() );
	printf(" WMC4ConMan::Load() => FWC2ELEMENTMASK     %s\n"    , m_SFWC2Mask.data());
	printf(" WMC4ConMan::Load() => FWC2DX1             %s mm\n" , m_SFWC2DX1.data() );
	printf(" WMC4ConMan::Load() => FWC2DX2             %s mm\n" , m_SFWC2DX2.data() );
	printf(" WMC4ConMan::Load() => FWC2DY              %s mm\n" , m_SFWC2DY.data()  );
	printf(" WMC4ConMan::Load() => FWC2THICKNESS       %s mm\n" , m_SFWC2Thi.data() );
	printf(" WMC4ConMan::Load() => FWC2RMIN            %s mm\n" , m_SFWC2RMin.data());
	printf(" WMC4ConMan::Load() => FWC2RMAX            %s mm\n" , m_SFWC2RMax.data());
	printf(" WMC4ConMan::Load() => FWC2POSITIONR       %s mm\n" , m_SFWC2PR.data()  );
	printf(" WMC4ConMan::Load() => FWC2POSITIONZ       %s mm\n" , m_SFWC2PZ.data()  );
	printf(" WMC4ConMan::Load() => FWC2THRESHOLD       %s MeV\n", m_SFWC2Thr.data() );
	// FTH
	printf(" WMC4ConMan::Load() => USEFTH              %s\n"    , m_SUseFTH.data()  );
	printf(" WMC4ConMan::Load() => FTHELEMENTMASK      %s\n"    , m_SFTHMask.data() );
	printf(" WMC4ConMan::Load() => FTHDX1              %s mm\n" , m_SFTHDX1.data()  );
	printf(" WMC4ConMan::Load() => FTHDX2              %s mm\n" , m_SFTHDX2.data()  );
	printf(" WMC4ConMan::Load() => FTHDY               %s mm\n" , m_SFTHDY.data()   );
	printf(" WMC4ConMan::Load() => FTHTHICKNESS        %s mm\n" , m_SFTHThi.data()  );
	printf(" WMC4ConMan::Load() => FTHRMIN             %s mm\n" , m_SFTHRMin.data() );
	printf(" WMC4ConMan::Load() => FTHRMAX             %s mm\n" , m_SFTHRMax.data() );
	printf(" WMC4ConMan::Load() => FTHPOSITIONR        %s mm\n" , m_SFTHPR.data()   );
	printf(" WMC4ConMan::Load() => FTHPOSITIONZ        %s mm\n" , m_SFTHPZ.data()   );
	printf(" WMC4ConMan::Load() => FTHTHRESHOLD        %s MeV\n", m_SFTHThr.data()  );
	// FRH1
	printf(" WMC4ConMan::Load() => USEFRH1             %s\n"    , m_SUseFRH1.data() );
	printf(" WMC4ConMan::Load() => FRH1ELEMENTMASK     %s\n"    , m_SFRH1Mask.data());
	printf(" WMC4ConMan::Load() => FRH1DX1             %s mm\n" , m_SFRH1DX1.data() );
	printf(" WMC4ConMan::Load() => FRH1DX2             %s mm\n" , m_SFRH1DX2.data() );
	printf(" WMC4ConMan::Load() => FRH1DY              %s mm\n" , m_SFRH1DY.data()  );
	printf(" WMC4ConMan::Load() => FRH1THICKNESS       %s mm\n" , m_SFRH1Thi.data() );
	printf(" WMC4ConMan::Load() => FRH1RMIN            %s mm\n" , m_SFRH1RMin.data());
	printf(" WMC4ConMan::Load() => FRH1RMAX            %s mm\n" , m_SFRH1RMax.data());
	printf(" WMC4ConMan::Load() => FRH1POSITIONR       %s mm\n" , m_SFRH1PR.data()  );
	printf(" WMC4ConMan::Load() => FRH1POSITIONZ       %s mm\n" , m_SFRH1PZ.data()  );
	printf(" WMC4ConMan::Load() => FRH1THRESHOLD       %s MeV\n", m_SFRH1Thr.data() );
	// FRH2
	printf(" WMC4ConMan::Load() => USEFRH2             %s\n"    , m_SUseFRH2.data() );
	printf(" WMC4ConMan::Load() => FRH2ELEMENTMASK     %s\n"    , m_SFRH2Mask.data());
	printf(" WMC4ConMan::Load() => FRH2DX1             %s mm\n" , m_SFRH2DX1.data() );
	printf(" WMC4ConMan::Load() => FRH2DX2             %s mm\n" , m_SFRH2DX2.data() );
	printf(" WMC4ConMan::Load() => FRH2DY              %s mm\n" , m_SFRH2DY.data()  );
	printf(" WMC4ConMan::Load() => FRH2THICKNESS       %s mm\n" , m_SFRH2Thi.data() );
	printf(" WMC4ConMan::Load() => FRH2RMIN            %s mm\n" , m_SFRH2RMin.data());
	printf(" WMC4ConMan::Load() => FRH2RMAX            %s mm\n" , m_SFRH2RMax.data());
	printf(" WMC4ConMan::Load() => FRH2POSITIONR       %s mm\n" , m_SFRH2PR.data()  );
	printf(" WMC4ConMan::Load() => FRH2POSITIONZ       %s mm\n" , m_SFRH2PZ.data()  );
	printf(" WMC4ConMan::Load() => FRH2THRESHOLD       %s MeV\n", m_SFRH2Thr.data() );
	// FRH3
	printf(" WMC4ConMan::Load() => USEFRH3             %s\n"    , m_SUseFRH3.data() );
	printf(" WMC4ConMan::Load() => FRH3ELEMENTMASK     %s\n"    , m_SFRH3Mask.data());
	printf(" WMC4ConMan::Load() => FRH3DX1             %s mm\n" , m_SFRH3DX1.data() );
	printf(" WMC4ConMan::Load() => FRH3DX2             %s mm\n" , m_SFRH3DX2.data() );
	printf(" WMC4ConMan::Load() => FRH3DY              %s mm\n" , m_SFRH3DY.data()  );
	printf(" WMC4ConMan::Load() => FRH3THICKNESS       %s mm\n" , m_SFRH3Thi.data() );
	printf(" WMC4ConMan::Load() => FRH3RMIN            %s mm\n" , m_SFRH3RMin.data());
	printf(" WMC4ConMan::Load() => FRH3RMAX            %s mm\n" , m_SFRH3RMax.data());
	printf(" WMC4ConMan::Load() => FRH3POSITIONR       %s mm\n" , m_SFRH3PR.data()  );
	printf(" WMC4ConMan::Load() => FRH3POSITIONZ       %s mm\n" , m_SFRH3PZ.data()  );
	printf(" WMC4ConMan::Load() => FRH3THRESHOLD       %s MeV\n", m_SFRH3Thr.data() );
	// FRH4
	printf(" WMC4ConMan::Load() => USEFRH4             %s\n"    , m_SUseFRH4.data() );
	printf(" WMC4ConMan::Load() => FRH4ELEMENTMASK     %s\n"    , m_SFRH4Mask.data());
	printf(" WMC4ConMan::Load() => FRH4DX1             %s mm\n" , m_SFRH4DX1.data() );
	printf(" WMC4ConMan::Load() => FRH4DX2             %s mm\n" , m_SFRH4DX2.data() );
	printf(" WMC4ConMan::Load() => FRH4DY              %s mm\n" , m_SFRH4DY.data()  );
	printf(" WMC4ConMan::Load() => FRH4THICKNESS       %s mm\n" , m_SFRH4Thi.data() );
	printf(" WMC4ConMan::Load() => FRH4RMIN            %s mm\n" , m_SFRH4RMin.data());
	printf(" WMC4ConMan::Load() => FRH4RMAX            %s mm\n" , m_SFRH4RMax.data());
	printf(" WMC4ConMan::Load() => FRH4POSITIONR       %s mm\n" , m_SFRH4PR.data()  );
	printf(" WMC4ConMan::Load() => FRH4POSITIONZ       %s mm\n" , m_SFRH4PZ.data()  );
	printf(" WMC4ConMan::Load() => FRH4THRESHOLD       %s MeV\n", m_SFRH4Thr.data() );
	// FRH5
	printf(" WMC4ConMan::Load() => USEFRH5             %s\n"    , m_SUseFRH5.data() );
	printf(" WMC4ConMan::Load() => FRH5ELEMENTMASK     %s\n"    , m_SFRH5Mask.data());
	printf(" WMC4ConMan::Load() => FRH5DX1             %s mm\n" , m_SFRH5DX1.data() );
	printf(" WMC4ConMan::Load() => FRH5DX2             %s mm\n" , m_SFRH5DX2.data() );
	printf(" WMC4ConMan::Load() => FRH5DY              %s mm\n" , m_SFRH5DY.data()  );
	printf(" WMC4ConMan::Load() => FRH5THICKNESS       %s mm\n" , m_SFRH5Thi.data() );
	printf(" WMC4ConMan::Load() => FRH5RMIN            %s mm\n" , m_SFRH5RMin.data());
	printf(" WMC4ConMan::Load() => FRH5RMAX            %s mm\n" , m_SFRH5RMax.data());
	printf(" WMC4ConMan::Load() => FRH5POSITIONR       %s mm\n" , m_SFRH5PR.data()  );
	printf(" WMC4ConMan::Load() => FRH5POSITIONZ       %s mm\n" , m_SFRH5PZ.data()  );
	printf(" WMC4ConMan::Load() => FRH5THRESHOLD       %s MeV\n", m_SFRH5Thr.data() );
	// Wrapping
	printf(" WMC4ConMan::Load() => MYLARTHICKNESS      %s mm\n" , m_SMylThi.data()  );
	printf(" WMC4ConMan::Load() => TEDLARTHICKNESS     %s mm\n" , m_STedThi.data()  );
	// FPC
	printf(" WMC4ConMan::Load() => FPC12PHI            %s deg\n", m_SFPC12Phi.data());
	printf(" WMC4ConMan::Load() => FPC34PHI            %s deg\n", m_SFPC34Phi.data());
	printf(" WMC4ConMan::Load() => FPC56PHI            %s deg\n", m_SFPC56Phi.data());
	printf(" WMC4ConMan::Load() => FPC78PHI            %s deg\n", m_SFPC78Phi.data());
	printf(" WMC4ConMan::Load() => FPC12POSITIONZ      %s mm\n" , m_SFPC12PZ.data() );
	printf(" WMC4ConMan::Load() => FPC34POSITIONZ      %s mm\n" , m_SFPC34PZ.data() );
	printf(" WMC4ConMan::Load() => FPC56POSITIONZ      %s mm\n" , m_SFPC56PZ.data() );
	printf(" WMC4ConMan::Load() => FPC78POSITIONZ      %s mm\n" , m_SFPC78PZ.data() );
	printf(" WMC4ConMan::Load() => FPCSHIELDTHICKNESS  %s mm\n" , m_SFPCAlThi.data());
	printf(" WMC4ConMan::Load() => TUBEDIAMETER        %s mm\n" , m_STubeDia.data() );
	printf(" WMC4ConMan::Load() => MYLARWALLTHICKNESS  %s mm\n" , m_SWallThi.data() );
	printf(" WMC4ConMan::Load() => WIREDIAMETER        %s mm\n" , m_SWireDia.data() );
	printf(" WMC4ConMan::Load() => TUBELENGTH          %s mm\n" , m_STubeL.data()   );
	printf(" WMC4ConMan::Load() => FPCGAS              %s\n"    , m_SFPCGas.data()  );
	printf(" WMC4ConMan::Load() => FPCTHRESHOLD        %s MeV\n", m_SFPCThr.data()  );
	printf(" WMC4ConMan::Load() => FPCVISUALIZE        %s\n"    , m_SFPCVis.data()  );
	// Trigger
	printf(" WMC4ConMan::Load() => TRIGGER             %s\n"    , m_STrigger.data() );
}

//////////////////////////////////////////////////
//   Set & get configurations                   //
//////////////////////////////////////////////////
// Computing
void WMC4ConMan::SetNofTRD(const G4int nofTRD)
{
	m_NofTRD = nofTRD;
}
G4int WMC4ConMan::GetNofTRD() const
{	
	return m_NofTRD;
}

// World
void WMC4ConMan::SetWorldX(const G4double worldX)
{
	m_WorldX = worldX;
}
G4double WMC4ConMan::GetWorldX() const
{
	return m_WorldX;
}

void WMC4ConMan::SetWorldY(const G4double worldY)
{
	m_WorldY = worldY;
}
G4double WMC4ConMan::GetWorldY() const
{
	return m_WorldY;
}

void WMC4ConMan::SetWorldZ(const G4double worldZ)
{
	m_WorldZ = worldZ;
}
G4double WMC4ConMan::GetWorldZ() const
{
	return m_WorldZ;
}

// Input
void WMC4ConMan::SetUseInput(const G4bool useInput)
{
	m_UseInput = useInput;
}
G4bool WMC4ConMan::GetUseInput() const
{
	return m_UseInput;
}
void WMC4ConMan::SetOutForm(const G4int outForm)
{
	m_OutForm = outForm;
}
G4int WMC4ConMan::GetOutForm() const
{
	return m_OutForm;
}

// Beam
void WMC4ConMan::SetBeamPX(const G4double beamPX)
{
	m_BeamPX = beamPX;
}
G4double WMC4ConMan::GetBeamPX() const
{
	return m_BeamPX;
}

void WMC4ConMan::SetBeamPY(const G4double beamPY)
{
	m_BeamPY = beamPY;
}
G4double WMC4ConMan::GetBeamPY() const
{
	return m_BeamPY;
}

void WMC4ConMan::SetBeamDX(const G4double beamDX)
{
	m_BeamDX = beamDX;
}
G4double WMC4ConMan::GetBeamDX() const
{
	return m_BeamDX;
}

void WMC4ConMan::SetBeamDY(const G4double beamDY)
{
	m_BeamDY = beamDY;
}
G4double WMC4ConMan::GetBeamDY() const
{
	return m_BeamDY;
}

void WMC4ConMan::SetBeamKE(const G4double beamKE)
{
	m_BeamKE = beamKE;
}
G4double WMC4ConMan::GetBeamKE() const
{
	return m_BeamKE;
}

// Pipe & chamber
void WMC4ConMan::SetPipeDia(const G4double pipeDia)
{
	m_PipeDia = pipeDia;
}
G4double WMC4ConMan::GetPipeDia() const
{
	return m_PipeDia;
}

void WMC4ConMan::SetVacThi(const G4double vacThi)
{
	m_VacThi = vacThi;
}
G4double WMC4ConMan::GetVacThi() const
{
	return m_VacThi;
}

// Target
void WMC4ConMan::SetUseTar(const G4bool useTar)
{
	m_UseTar = useTar;
}
G4bool WMC4ConMan::GetUseTar() const
{
	return m_UseTar;
}

void WMC4ConMan::SetTarMat(const G4String tarMat)
{
	m_TarMat = tarMat;
}
G4String WMC4ConMan::GetTarMat() const
{
	return m_TarMat;
}

void WMC4ConMan::SetTarThi(const G4double tarThi)
{
	m_TarMat = tarThi;
}
G4double WMC4ConMan::GetTarThi() const
{
	return m_TarThi;
}

// FWC1
void WMC4ConMan::SetUseFWC1(const G4bool useFWC1)
{
	m_UseFWC1 = useFWC1;
}
G4bool WMC4ConMan::GetUseFWC1() const
{
	return m_UseFWC1;
}

void WMC4ConMan::SetFWC1Mask(const std::vector<G4bool> fwc1Mask)
{
	m_FWC1Mask = fwc1Mask;
}
std::vector<G4bool> WMC4ConMan::GetFWC1Mask() const
{
	return m_FWC1Mask;
}

void WMC4ConMan::SetFWC1DX1(const G4double fwc1DX1)
{
	m_FWC1DX1 = fwc1DX1;
}
G4double WMC4ConMan::GetFWC1DX1() const
{
	return m_FWC1DX1;
}

void WMC4ConMan::SetFWC1DX2(const G4double fwc1DX2)
{
	m_FWC1DX2 = fwc1DX2;
}
G4double WMC4ConMan::GetFWC1DX2() const
{
	return m_FWC1DX2;
}

void WMC4ConMan::SetFWC1DY(const G4double fwc1DY)
{
	m_FWC1DY = fwc1DY;
}
G4double WMC4ConMan::GetFWC1DY() const
{
	return m_FWC1DY;
}

void WMC4ConMan::SetFWC1Thi(const G4double fwc1Thi)
{
	m_FWC1Thi = fwc1Thi;
}
G4double WMC4ConMan::GetFWC1Thi() const
{
	return m_FWC1Thi;
}

void WMC4ConMan::SetFWC1RMin(const G4double fwc1RMin)
{
	m_FWC1RMin = fwc1RMin;
}
G4double WMC4ConMan::GetFWC1RMin() const
{
	return m_FWC1RMin;
}

void WMC4ConMan::SetFWC1RMax(const G4double fwc1RMax)
{
	m_FWC1RMax = fwc1RMax;
}
G4double WMC4ConMan::GetFWC1RMax() const
{
	return m_FWC1RMax;
}

void WMC4ConMan::SetFWC1PR(const G4double fwc1PR)
{
	m_FWC1PR = fwc1PR;
}
G4double WMC4ConMan::GetFWC1PR() const
{
	return m_FWC1PR;
}

void WMC4ConMan::SetFWC1PZ(const G4double fwc1PZ)
{
	m_FWC1PZ = fwc1PZ;
}
G4double WMC4ConMan::GetFWC1PZ() const
{
	return m_FWC1PZ;
}

void WMC4ConMan::SetFWC1Thr(const G4double fwc1Thr)
{
	m_FWC1Thr = fwc1Thr;
}
G4double WMC4ConMan::GetFWC1Thr() const
{
	return m_FWC1Thr;
}

// FWC2
void WMC4ConMan::SetUseFWC2(const G4bool useFWC2)
{
	m_UseFWC2 = useFWC2;
}
G4bool WMC4ConMan::GetUseFWC2() const
{
	return m_UseFWC2;
}

void WMC4ConMan::SetFWC2Mask(const std::vector<G4bool> fwc2Mask)
{
	m_FWC2Mask = fwc2Mask;
}
std::vector<G4bool> WMC4ConMan::GetFWC2Mask() const
{
	return m_FWC2Mask;
}

void WMC4ConMan::SetFWC2DX1(const G4double fwc2DX1)
{
	m_FWC2DX1 = fwc2DX1;
}
G4double WMC4ConMan::GetFWC2DX1() const
{
	return m_FWC2DX1;
}

void WMC4ConMan::SetFWC2DX2(const G4double fwc2DX2)
{
	m_FWC2DX2 = fwc2DX2;
}
G4double WMC4ConMan::GetFWC2DX2() const
{
	return m_FWC2DX2;
}

void WMC4ConMan::SetFWC2DY(const G4double fwc2DY)
{
	m_FWC2DY = fwc2DY;
}
G4double WMC4ConMan::GetFWC2DY() const
{
	return m_FWC2DY;
}

void WMC4ConMan::SetFWC2Thi(const G4double fwc2Thi)
{
	m_FWC2Thi = fwc2Thi;
}
G4double WMC4ConMan::GetFWC2Thi() const
{
	return m_FWC2Thi;
}

void WMC4ConMan::SetFWC2RMin(const G4double fwc2RMin)
{
	m_FWC2RMin = fwc2RMin;
}
G4double WMC4ConMan::GetFWC2RMin() const
{
	return m_FWC2RMin;
}

void WMC4ConMan::SetFWC2RMax(const G4double fwc2RMax)
{
	m_FWC2RMax = fwc2RMax;
}
G4double WMC4ConMan::GetFWC2RMax() const
{
	return m_FWC2RMax;
}

void WMC4ConMan::SetFWC2PR(const G4double fwc2PR)
{
	m_FWC2PR = fwc2PR;
}
G4double WMC4ConMan::GetFWC2PR() const
{
	return m_FWC2PR;
}

void WMC4ConMan::SetFWC2PZ(const G4double fwc2PZ)
{
	m_FWC2PZ = fwc2PZ;
}
G4double WMC4ConMan::GetFWC2PZ() const
{
	return m_FWC2PZ;
}

void WMC4ConMan::SetFWC2Thr(const G4double fwc2Thr)
{
	m_FWC2Thr = fwc2Thr;
}
G4double WMC4ConMan::GetFWC2Thr() const
{
	return m_FWC2Thr;
}

// FTH
void WMC4ConMan::SetUseFTH(const G4bool useFTH)
{
	m_UseFTH = useFTH;
}
G4bool WMC4ConMan::GetUseFTH() const
{
	return m_UseFTH;
}

void WMC4ConMan::SetFTHMask(const std::vector<G4bool> fthMask)
{
	m_FTHMask = fthMask;
}
std::vector<G4bool> WMC4ConMan::GetFTHMask() const
{
	return m_FTHMask;
}

void WMC4ConMan::SetFTHDX1(const G4double fthDX1)
{
	m_FTHDX1 = fthDX1;
}
G4double WMC4ConMan::GetFTHDX1() const
{
	return m_FTHDX1;
}

void WMC4ConMan::SetFTHDX2(const G4double fthDX2)
{
	m_FTHDX2 = fthDX2;
}
G4double WMC4ConMan::GetFTHDX2() const
{
	return m_FTHDX2;
}

void WMC4ConMan::SetFTHDY(const G4double fthDY)
{
	m_FTHDY = fthDY;
}
G4double WMC4ConMan::GetFTHDY() const
{
	return m_FTHDY;
}

void WMC4ConMan::SetFTHThi(const G4double fthThi)
{
	m_FTHThi = fthThi;
}
G4double WMC4ConMan::GetFTHThi() const
{
	return m_FTHThi;
}

void WMC4ConMan::SetFTHRMin(const G4double fthRMin)
{
	m_FTHRMin = fthRMin;
}
G4double WMC4ConMan::GetFTHRMin() const
{
	return m_FTHRMin;
}

void WMC4ConMan::SetFTHRMax(const G4double fthRMax)
{
	m_FTHRMax = fthRMax;
}
G4double WMC4ConMan::GetFTHRMax() const
{
	return m_FTHRMax;
}

void WMC4ConMan::SetFTHPR(const G4double fthPR)
{
	m_FTHPR = fthPR;
}
G4double WMC4ConMan::GetFTHPR() const
{
	return m_FTHPR;
}

void WMC4ConMan::SetFTHPZ(const G4double fthPZ)
{
	m_FTHPZ = fthPZ;
}
G4double WMC4ConMan::GetFTHPZ() const
{
	return m_FTHPZ;
}

void WMC4ConMan::SetFTHThr(const G4double fthThr)
{
	m_FTHThr = fthThr;
}
G4double WMC4ConMan::GetFTHThr() const
{
	return m_FTHThr;
}

// FRH1
void WMC4ConMan::SetUseFRH1(const G4bool useFRH1)
{
	m_UseFRH1 = useFRH1;
}
G4bool WMC4ConMan::GetUseFRH1() const
{
	return m_UseFRH1;
}

void WMC4ConMan::SetFRH1Mask(const std::vector<G4bool> frh1Mask)
{
	m_FRH1Mask = frh1Mask;
}
std::vector<G4bool> WMC4ConMan::GetFRH1Mask() const
{
	return m_FRH1Mask;
}

void WMC4ConMan::SetFRH1DX1(const G4double frh1DX1)
{
	m_FRH1DX1 = frh1DX1;
}
G4double WMC4ConMan::GetFRH1DX1() const
{
	return m_FRH1DX1;
}

void WMC4ConMan::SetFRH1DX2(const G4double frh1DX2)
{
	m_FRH1DX2 = frh1DX2;
}
G4double WMC4ConMan::GetFRH1DX2() const
{
	return m_FRH1DX2;
}

void WMC4ConMan::SetFRH1DY(const G4double frh1DY)
{
	m_FRH1DY = frh1DY;
}
G4double WMC4ConMan::GetFRH1DY() const
{
	return m_FRH1DY;
}

void WMC4ConMan::SetFRH1Thi(const G4double frh1Thi)
{
	m_FRH1Thi = frh1Thi;
}
G4double WMC4ConMan::GetFRH1Thi() const
{
	return m_FRH1Thi;
}

void WMC4ConMan::SetFRH1RMin(const G4double frh1RMin)
{
	m_FRH1RMin = frh1RMin;
}
G4double WMC4ConMan::GetFRH1RMin() const
{
	return m_FRH1RMin;
}

void WMC4ConMan::SetFRH1RMax(const G4double frh1RMax)
{
	m_FRH1RMax = frh1RMax;
}
G4double WMC4ConMan::GetFRH1RMax() const
{
	return m_FRH1RMax;
}

void WMC4ConMan::SetFRH1PR(const G4double frh1PR)
{
	m_FRH1PR = frh1PR;
}
G4double WMC4ConMan::GetFRH1PR() const
{
	return m_FRH1PR;
}

void WMC4ConMan::SetFRH1PZ(const G4double frh1PZ)
{
	m_FRH1PZ = frh1PZ;
}
G4double WMC4ConMan::GetFRH1PZ() const
{
	return m_FRH1PZ;
}

void WMC4ConMan::SetFRH1Thr(const G4double frh1Thr)
{
	m_FRH1Thr = frh1Thr;
}
G4double WMC4ConMan::GetFRH1Thr() const
{
	return m_FRH1Thr;
}

// FRH2
void WMC4ConMan::SetUseFRH2(const G4bool useFRH2)
{
	m_UseFRH2 = useFRH2;
}
G4bool WMC4ConMan::GetUseFRH2() const
{
	return m_UseFRH2;
}

void WMC4ConMan::SetFRH2Mask(const std::vector<G4bool> frh2Mask)
{
	m_FRH2Mask = frh2Mask;
}
std::vector<G4bool> WMC4ConMan::GetFRH2Mask() const
{
	return m_FRH2Mask;
}

void WMC4ConMan::SetFRH2DX1(const G4double frh2DX1)
{
	m_FRH2DX1 = frh2DX1;
}
G4double WMC4ConMan::GetFRH2DX1() const
{
	return m_FRH2DX1;
}

void WMC4ConMan::SetFRH2DX2(const G4double frh2DX2)
{
	m_FRH2DX2 = frh2DX2;
}
G4double WMC4ConMan::GetFRH2DX2() const
{
	return m_FRH2DX2;
}

void WMC4ConMan::SetFRH2DY(const G4double frh2DY)
{
	m_FRH2DY = frh2DY;
}
G4double WMC4ConMan::GetFRH2DY() const
{
	return m_FRH2DY;
}

void WMC4ConMan::SetFRH2Thi(const G4double frh2Thi)
{
	m_FRH2Thi = frh2Thi;
}
G4double WMC4ConMan::GetFRH2Thi() const
{
	return m_FRH2Thi;
}

void WMC4ConMan::SetFRH2RMin(const G4double frh2RMin)
{
	m_FRH2RMin = frh2RMin;
}
G4double WMC4ConMan::GetFRH2RMin() const
{
	return m_FRH2RMin;
}

void WMC4ConMan::SetFRH2RMax(const G4double frh2RMax)
{
	m_FRH2RMax = frh2RMax;
}
G4double WMC4ConMan::GetFRH2RMax() const
{
	return m_FRH2RMax;
}

void WMC4ConMan::SetFRH2PR(const G4double frh2PR)
{
	m_FRH2PR = frh2PR;
}
G4double WMC4ConMan::GetFRH2PR() const
{
	return m_FRH2PR;
}

void WMC4ConMan::SetFRH2PZ(const G4double frh2PZ)
{
	m_FRH2PZ = frh2PZ;
}
G4double WMC4ConMan::GetFRH2PZ() const
{
	return m_FRH2PZ;
}

void WMC4ConMan::SetFRH2Thr(const G4double frh2Thr)
{
	m_FRH2Thr = frh2Thr;
}
G4double WMC4ConMan::GetFRH2Thr() const
{
	return m_FRH2Thr;
}

// FRH3
void WMC4ConMan::SetUseFRH3(const G4bool useFRH3)
{
	m_UseFRH3 = useFRH3;
}
G4bool WMC4ConMan::GetUseFRH3() const
{
	return m_UseFRH3;
}

void WMC4ConMan::SetFRH3Mask(const std::vector<G4bool> frh3Mask)
{
	m_FRH3Mask = frh3Mask;
}
std::vector<G4bool> WMC4ConMan::GetFRH3Mask() const
{
	return m_FRH3Mask;
}

void WMC4ConMan::SetFRH3DX1(const G4double frh3DX1)
{
	m_FRH3DX1 = frh3DX1;
}
G4double WMC4ConMan::GetFRH3DX1() const
{
	return m_FRH3DX1;
}

void WMC4ConMan::SetFRH3DX2(const G4double frh3DX2)
{
	m_FRH3DX2 = frh3DX2;
}
G4double WMC4ConMan::GetFRH3DX2() const
{
	return m_FRH3DX2;
}

void WMC4ConMan::SetFRH3DY(const G4double frh3DY)
{
	m_FRH3DY = frh3DY;
}
G4double WMC4ConMan::GetFRH3DY() const
{
	return m_FRH3DY;
}

void WMC4ConMan::SetFRH3Thi(const G4double frh3Thi)
{
	m_FRH3Thi = frh3Thi;
}
G4double WMC4ConMan::GetFRH3Thi() const
{
	return m_FRH3Thi;
}

void WMC4ConMan::SetFRH3RMin(const G4double frh3RMin)
{
	m_FRH3RMin = frh3RMin;
}
G4double WMC4ConMan::GetFRH3RMin() const
{
	return m_FRH3RMin;
}

void WMC4ConMan::SetFRH3RMax(const G4double frh3RMax)
{
	m_FRH3RMax = frh3RMax;
}
G4double WMC4ConMan::GetFRH3RMax() const
{
	return m_FRH3RMax;
}

void WMC4ConMan::SetFRH3PR(const G4double frh3PR)
{
	m_FRH3PR = frh3PR;
}
G4double WMC4ConMan::GetFRH3PR() const
{
	return m_FRH3PR;
}

void WMC4ConMan::SetFRH3PZ(const G4double frh3PZ)
{
	m_FRH3PZ = frh3PZ;
}
G4double WMC4ConMan::GetFRH3PZ() const
{
	return m_FRH3PZ;
}

void WMC4ConMan::SetFRH3Thr(const G4double frh3Thr)
{
	m_FRH3Thr = frh3Thr;
}
G4double WMC4ConMan::GetFRH3Thr() const
{
	return m_FRH3Thr;
}

// FRH4
void WMC4ConMan::SetUseFRH4(const G4bool useFRH4)
{
	m_UseFRH4 = useFRH4;
}
G4bool WMC4ConMan::GetUseFRH4() const
{
	return m_UseFRH4;
}

void WMC4ConMan::SetFRH4Mask(const std::vector<G4bool> frh4Mask)
{
	m_FRH4Mask = frh4Mask;
}
std::vector<G4bool> WMC4ConMan::GetFRH4Mask() const
{
	return m_FRH4Mask;
}

void WMC4ConMan::SetFRH4DX1(const G4double frh4DX1)
{
	m_FRH4DX1 = frh4DX1;
}
G4double WMC4ConMan::GetFRH4DX1() const
{
	return m_FRH4DX1;
}

void WMC4ConMan::SetFRH4DX2(const G4double frh4DX2)
{
	m_FRH4DX2 = frh4DX2;
}
G4double WMC4ConMan::GetFRH4DX2() const
{
	return m_FRH4DX2;
}

void WMC4ConMan::SetFRH4DY(const G4double frh4DY)
{
	m_FRH4DY = frh4DY;
}
G4double WMC4ConMan::GetFRH4DY() const
{
	return m_FRH4DY;
}

void WMC4ConMan::SetFRH4Thi(const G4double frh4Thi)
{
	m_FRH4Thi = frh4Thi;
}
G4double WMC4ConMan::GetFRH4Thi() const
{
	return m_FRH4Thi;
}

void WMC4ConMan::SetFRH4RMin(const G4double frh4RMin)
{
	m_FRH4RMin = frh4RMin;
}
G4double WMC4ConMan::GetFRH4RMin() const
{
	return m_FRH4RMin;
}

void WMC4ConMan::SetFRH4RMax(const G4double frh4RMax)
{
	m_FRH4RMax = frh4RMax;
}
G4double WMC4ConMan::GetFRH4RMax() const
{
	return m_FRH4RMax;
}

void WMC4ConMan::SetFRH4PR(const G4double frh4PR)
{
	m_FRH4PR = frh4PR;
}
G4double WMC4ConMan::GetFRH4PR() const
{
	return m_FRH4PR;
}

void WMC4ConMan::SetFRH4PZ(const G4double frh4PZ)
{
	m_FRH4PZ = frh4PZ;
}
G4double WMC4ConMan::GetFRH4PZ() const
{
	return m_FRH4PZ;
}

void WMC4ConMan::SetFRH4Thr(const G4double frh4Thr)
{
	m_FRH4Thr = frh4Thr;
}
G4double WMC4ConMan::GetFRH4Thr() const
{
	return m_FRH4Thr;
}

// FRH5
void WMC4ConMan::SetUseFRH5(const G4bool useFRH5)
{
	m_UseFRH5 = useFRH5;
}
G4bool WMC4ConMan::GetUseFRH5() const
{
	return m_UseFRH5;
}

void WMC4ConMan::SetFRH5Mask(const std::vector<G4bool> frh5Mask)
{
	m_FRH5Mask = frh5Mask;
}
std::vector<G4bool> WMC4ConMan::GetFRH5Mask() const
{
	return m_FRH5Mask;
}

void WMC4ConMan::SetFRH5DX1(const G4double frh5DX1)
{
	m_FRH5DX1 = frh5DX1;
}
G4double WMC4ConMan::GetFRH5DX1() const
{
	return m_FRH5DX1;
}

void WMC4ConMan::SetFRH5DX2(const G4double frh5DX2)
{
	m_FRH5DX2 = frh5DX2;
}
G4double WMC4ConMan::GetFRH5DX2() const
{
	return m_FRH5DX2;
}

void WMC4ConMan::SetFRH5DY(const G4double frh5DY)
{
	m_FRH5DY = frh5DY;
}
G4double WMC4ConMan::GetFRH5DY() const
{
	return m_FRH5DY;
}

void WMC4ConMan::SetFRH5Thi(const G4double frh5Thi)
{
	m_FRH5Thi = frh5Thi;
}
G4double WMC4ConMan::GetFRH5Thi() const
{
	return m_FRH5Thi;
}

void WMC4ConMan::SetFRH5RMin(const G4double frh5RMin)
{
	m_FRH5RMin = frh5RMin;
}
G4double WMC4ConMan::GetFRH5RMin() const
{
	return m_FRH5RMin;
}

void WMC4ConMan::SetFRH5RMax(const G4double frh5RMax)
{
	m_FRH5RMax = frh5RMax;
}
G4double WMC4ConMan::GetFRH5RMax() const
{
	return m_FRH5RMax;
}

void WMC4ConMan::SetFRH5PR(const G4double frh5PR)
{
	m_FRH5PR = frh5PR;
}
G4double WMC4ConMan::GetFRH5PR() const
{
	return m_FRH5PR;
}

void WMC4ConMan::SetFRH5PZ(const G4double frh5PZ)
{
	m_FRH5PZ = frh5PZ;
}
G4double WMC4ConMan::GetFRH5PZ() const
{
	return m_FRH5PZ;
}

void WMC4ConMan::SetFRH5Thr(const G4double frh5Thr)
{
	m_FRH5Thr = frh5Thr;
}
G4double WMC4ConMan::GetFRH5Thr() const
{
	return m_FRH5Thr;
}

// Wrapping
void WMC4ConMan::SetMylThi(const G4double mylThi)
{
	m_MylThi = mylThi;
}
G4double WMC4ConMan::GetMylThi() const
{
	return m_MylThi;
}

void WMC4ConMan::SetTedThi(const G4double tedThi)
{
	m_TedThi = tedThi;
}
G4double WMC4ConMan::GetTedThi() const
{
	return m_TedThi;
}

// FPC
void WMC4ConMan::SetFPC12Phi(const G4double fpc12Phi)
{
	m_FPC12Phi = fpc12Phi;
}
G4double WMC4ConMan::GetFPC12Phi() const
{
	return m_FPC12Phi;
}

void WMC4ConMan::SetFPC34Phi(const G4double fpc34Phi)
{
	m_FPC34Phi = fpc34Phi;
}
G4double WMC4ConMan::GetFPC34Phi() const
{
	return m_FPC34Phi;
}

void WMC4ConMan::SetFPC56Phi(const G4double fpc56Phi)
{
	m_FPC56Phi = fpc56Phi;
}
G4double WMC4ConMan::GetFPC56Phi() const
{
	return m_FPC56Phi;
}

void WMC4ConMan::SetFPC78Phi(const G4double fpc78Phi)
{
	m_FPC78Phi = fpc78Phi;
}
G4double WMC4ConMan::GetFPC78Phi() const
{
	return m_FPC78Phi;
}

void WMC4ConMan::SetFPC12PZ(const G4double fpc12PZ)
{
	m_FPC12PZ = fpc12PZ;
}
G4double WMC4ConMan::GetFPC12PZ() const
{
	return m_FPC12PZ;
}

void WMC4ConMan::SetFPC34PZ(const G4double fpc34PZ)
{
	m_FPC34PZ = fpc34PZ;
}
G4double WMC4ConMan::GetFPC34PZ() const
{
	return m_FPC34PZ;
}

void WMC4ConMan::SetFPC56PZ(const G4double fpc56PZ)
{
	m_FPC56PZ = fpc56PZ;
}
G4double WMC4ConMan::GetFPC56PZ() const
{
	return m_FPC56PZ;
}

void WMC4ConMan::SetFPC78PZ(const G4double fpc78PZ)
{
	m_FPC78PZ = fpc78PZ;
}
G4double WMC4ConMan::GetFPC78PZ() const
{
	return m_FPC78PZ;
}

void WMC4ConMan::SetFPCAlThi(const G4double fpcAlThi)
{
	m_FPCAlThi = fpcAlThi;
}
G4double WMC4ConMan::GetFPCAlThi() const
{
	return m_FPCAlThi;
}

void WMC4ConMan::SetTubeDia(const G4double tubeDia)
{
	m_TubeDia = tubeDia;
}
G4double WMC4ConMan::GetTubeDia() const
{
	return m_TubeDia;
}

void WMC4ConMan::SetWireDia(const G4double wireDia)
{
	m_WireDia = wireDia;
}
G4double WMC4ConMan::GetWireDia() const
{
	return m_WireDia;
}

void WMC4ConMan::SetWallThi(const G4double wallThi)
{
	m_WallThi = wallThi;
}
G4double WMC4ConMan::GetWallThi() const
{
	return m_WallThi;
}

void WMC4ConMan::SetTubeL(const G4double tubeL)
{
	m_TubeL = tubeL;
}
G4double WMC4ConMan::GetTubeL() const
{
	return m_TubeL;
}

void WMC4ConMan::SetFPCGas(const G4int fpcGas)
{
	m_FPCGas = fpcGas;
}
G4int WMC4ConMan::GetFPCGas() const
{
	return m_FPCGas;
}

void WMC4ConMan::SetFPCThr(const G4double fpcThr)
{
	m_FPCThr = fpcThr;
}
G4double WMC4ConMan::GetFPCThr() const
{
	return m_FPCThr;
}

void WMC4ConMan::SetFPCVis(const G4bool fpcVis)
{
	m_FPCVis = fpcVis;
}
G4bool WMC4ConMan::GetFPCVis() const
{
	return m_FPCVis;
}

// Trigger
void WMC4ConMan::SetTrigger(const G4int trigger)
{
	m_Trigger = trigger;
}
G4int WMC4ConMan::GetTrigger() const
{
	return m_Trigger;
}

//////////////////////////////////////////////////
//   Test whether a line starts with ...        //
//////////////////////////////////////////////////
bool WMC4ConMan::StartsWith(const std::string& text, const std::string& token)
{
	if ( text.length() < token.length() ) return false;
	return ( text.compare(0, token.length(), token) == 0 );
}
