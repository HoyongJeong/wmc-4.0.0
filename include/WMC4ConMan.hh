#ifndef WMC4CONMAN_h
#define WMC4CONMAN_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4ConMan.hh                                                            //
//                                                                            //
//   This file is a header for WMC4ConMan class. It's a configuration         //
// manager for simulation. The class contains all configuration parameters    //
// for simulation configuration from a file.                                  //
//                                                                            //
//                    - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>
#include "globals.hh"

class WMC4ConMan
{
  public:
	// Constructors and Destructor
	WMC4ConMan();
	WMC4ConMan(const char* fileName);
	~WMC4ConMan();

	// Initialize
	void Init(const char* fileName = 0);

	// Load configuration file
	bool Load(const char* fileName);
	void PrintConfiguration();

	// Set & get configurations: Computing
	void SetNofTRD(  const G4int    nofTRD  );
	G4int    GetNofTRD()   const;
	// Set & get configurations: World
	void SetWorldX(  const G4double worldX  );
	G4double GetWorldX()   const;
	void SetWorldY(  const G4double worldY  );
	G4double GetWorldY()   const;
	void SetWorldZ(  const G4double worldZ  );
	G4double GetWorldZ()   const;
	// Set & get configurations: Input/output
	void SetUseInput(const G4bool   useInput);
	G4bool GetUseInput()   const;
	void SetOutForm( const G4int    outForm );
	G4int GetOutForm()     const;
	// Set & get configurations: Beam
	void SetBeamPX(  const G4double beamPX  );
	G4double GetBeamPX()   const;
	void SetBeamPY(  const G4double beamPY  );
	G4double GetBeamPY()   const;
	void SetBeamDX(  const G4double beamDX  );
	G4double GetBeamDX()   const;
	void SetBeamDY(  const G4double beamDY  );
	G4double GetBeamDY()   const;
	void SetBeamKE(  const G4double beamKE  );
	G4double GetBeamKE()   const;
	// Set & get configurations: Pipe & chamber
	void SetPipeDia( const G4double pipeDia );
	G4double GetPipeDia()  const;
	void SetVacThi(  const G4double vacThi  );
	G4double GetVacThi()   const;
	// Set & get configurations: Target
	void SetUseTar(  const G4bool   useTar  );
	G4bool   GetUseTar()   const;
	void SetTarMat(  const G4String tarMat  );
	G4String GetTarMat()   const;
	void SetTarThi(  const G4double tarThi  );
	G4double GetTarThi()   const;
	// Set & get configurations: FWC1
	void SetUseFWC1( const G4bool   useFWC1 );
	G4bool   GetUseFWC1()  const;
	void SetFWC1Mask(const std::vector<G4bool> fwc1Mask);
	std::vector<G4bool> GetFWC1Mask() const;
	void SetFWC1DX1( const G4double fwc1DX1 );
	G4double GetFWC1DX1()  const;
	void SetFWC1DX2( const G4double fwc1DX2 );
	G4double GetFWC1DX2()  const;
	void SetFWC1DY(  const G4double fwc1DY  );
	G4double GetFWC1DY()   const;
	void SetFWC1Thi( const G4double fwc1Thi );
	G4double GetFWC1Thi()  const;
	void SetFWC1RMin(const G4double fwc1RMin);
	G4double GetFWC1RMin() const;
	void SetFWC1RMax(const G4double fwc1RMax);
	G4double GetFWC1RMax() const;
	void SetFWC1PR(  const G4double fwc1PR  );
	G4double GetFWC1PR()   const;
	void SetFWC1PZ(  const G4double fwc1PZ  );
	G4double GetFWC1PZ()   const;
	void SetFWC1Thr( const G4double fwc1Thr );
	G4double GetFWC1Thr()  const;
	// Set & get configurations: FWC2
	void SetUseFWC2( const G4bool   useFWC2 );
	G4bool   GetUseFWC2()  const;
	void SetFWC2Mask(const std::vector<G4bool> fwc2Mask);
	std::vector<G4bool> GetFWC2Mask() const;
	void SetFWC2DX1( const G4double fwc2DX1 );
	G4double GetFWC2DX1()  const;
	void SetFWC2DX2( const G4double fwc2DX2 );
	G4double GetFWC2DX2()  const;
	void SetFWC2DY(  const G4double fwc2DY  );
	G4double GetFWC2DY()   const;
	void SetFWC2Thi( const G4double fwc2Thi );
	G4double GetFWC2Thi()  const;
	void SetFWC2RMin(const G4double fwc2RMin);
	G4double GetFWC2RMin() const;
	void SetFWC2RMax(const G4double fwc2RMax);
	G4double GetFWC2RMax() const;
	void SetFWC2PR(  const G4double fwc2PR  );
	G4double GetFWC2PR()   const;
	void SetFWC2PZ(  const G4double fwc2PZ  );
	G4double GetFWC2PZ()   const;
	void SetFWC2Thr( const G4double fwc2Thr );
	G4double GetFWC2Thr()  const;
	// Set & get configurations: FTH
	void SetUseFTH(  const G4bool   useFTH  );
	G4bool   GetUseFTH()   const;
	void SetFTHMask( const std::vector<G4bool> fthMask);
	std::vector<G4bool> GetFTHMask() const;
	void SetFTHDX1(  const G4double fthDX1  );
	G4double GetFTHDX1()   const;
	void SetFTHDX2(  const G4double fthDX2  );
	G4double GetFTHDX2()   const;
	void SetFTHDY(   const G4double fthDY   );
	G4double GetFTHDY()    const;
	void SetFTHThi(  const G4double fthThi  );
	G4double GetFTHThi()   const;
	void SetFTHRMin( const G4double fthRMin );
	G4double GetFTHRMin()  const;
	void SetFTHRMax( const G4double fthRMax );
	G4double GetFTHRMax()  const;
	void SetFTHPR(   const G4double fthPR   );
	G4double GetFTHPR()    const;
	void SetFTHPZ(   const G4double fthPZ   );
	G4double GetFTHPZ()    const;
	void SetFTHThr(  const G4double fthThr  );
	G4double GetFTHThr()   const;
	// Set & get configurations: FRH1
	void SetUseFRH1( const G4bool   useFRH1 );
	G4bool   GetUseFRH1()  const;
	void SetFRH1Mask(const std::vector<G4bool> frh1Mask);
	std::vector<G4bool> GetFRH1Mask() const;
	void SetFRH1DX1( const G4double frh1DX1 );
	G4double GetFRH1DX1()  const;
	void SetFRH1DX2( const G4double frh1DX2 );
	G4double GetFRH1DX2()  const;
	void SetFRH1DY(  const G4double frh1DY  );
	G4double GetFRH1DY()   const;
	void SetFRH1Thi( const G4double frh1Thi );
	G4double GetFRH1Thi()  const;
	void SetFRH1RMin(const G4double frh1RMin);
	G4double GetFRH1RMin() const;
	void SetFRH1RMax(const G4double frh1RMax);
	G4double GetFRH1RMax() const;
	void SetFRH1PR(  const G4double frh1PR  );
	G4double GetFRH1PR()   const;
	void SetFRH1PZ(  const G4double frh1PZ  );
	G4double GetFRH1PZ()   const;
	void SetFRH1Thr( const G4double frh1Thr );
	G4double GetFRH1Thr()  const;
	// Set & get configurations: FRH2
	void SetUseFRH2( const G4bool   useFRH2 );
	G4bool   GetUseFRH2()  const;
	void SetFRH2Mask(const std::vector<G4bool> frh2Mask);
	std::vector<G4bool> GetFRH2Mask() const;
	void SetFRH2DX1( const G4double frh2DX1 );
	G4double GetFRH2DX1()  const;
	void SetFRH2DX2( const G4double frh2DX2 );
	G4double GetFRH2DX2()  const;
	void SetFRH2DY(  const G4double frh2DY  );
	G4double GetFRH2DY()   const;
	void SetFRH2Thi( const G4double frh2Thi );
	G4double GetFRH2Thi()  const;
	void SetFRH2RMin(const G4double frh2RMin);
	G4double GetFRH2RMin() const;
	void SetFRH2RMax(const G4double frh2RMax);
	G4double GetFRH2RMax() const;
	void SetFRH2PR(  const G4double frh2PR  );
	G4double GetFRH2PR()   const;
	void SetFRH2PZ(  const G4double frh2PZ  );
	G4double GetFRH2PZ()   const;
	void SetFRH2Thr( const G4double frh2Thr );
	G4double GetFRH2Thr()  const;
	// Set & get configurations: FRH3
	void SetUseFRH3( const G4bool   useFRH3 );
	G4bool   GetUseFRH3()  const;
	void SetFRH3Mask(const std::vector<G4bool> frh3Mask);
	std::vector<G4bool> GetFRH3Mask() const;
	void SetFRH3DX1( const G4double frh3DX1 );
	G4double GetFRH3DX1()  const;
	void SetFRH3DX2( const G4double frh3DX2 );
	G4double GetFRH3DX2()  const;
	void SetFRH3DY(  const G4double frh3DY  );
	G4double GetFRH3DY()   const;
	void SetFRH3Thi( const G4double frh3Thi );
	G4double GetFRH3Thi()  const;
	void SetFRH3RMin(const G4double frh3RMin);
	G4double GetFRH3RMin() const;
	void SetFRH3RMax(const G4double frh3RMax);
	G4double GetFRH3RMax() const;
	void SetFRH3PR(  const G4double frh3PR  );
	G4double GetFRH3PR()   const;
	void SetFRH3PZ(  const G4double frh3PZ  );
	G4double GetFRH3PZ()   const;
	void SetFRH3Thr( const G4double frh3Thr );
	G4double GetFRH3Thr()  const;
	// Set & get configurations: FRH4
	void SetUseFRH4( const G4bool   useFRH4 );
	G4bool   GetUseFRH4()  const;
	void SetFRH4Mask(const std::vector<G4bool> frh4Mask);
	std::vector<G4bool> GetFRH4Mask() const;
	void SetFRH4DX1( const G4double frh4DX1 );
	G4double GetFRH4DX1()  const;
	void SetFRH4DX2( const G4double frh4DX2 );
	G4double GetFRH4DX2()  const;
	void SetFRH4DY(  const G4double frh4DY  );
	G4double GetFRH4DY()   const;
	void SetFRH4Thi( const G4double frh4Thi );
	G4double GetFRH4Thi()  const;
	void SetFRH4RMin(const G4double frh4RMin);
	G4double GetFRH4RMin() const;
	void SetFRH4RMax(const G4double frh4RMax);
	G4double GetFRH4RMax() const;
	void SetFRH4PR(  const G4double frh4PR  );
	G4double GetFRH4PR()   const;
	void SetFRH4PZ(  const G4double frh4PZ  );
	G4double GetFRH4PZ()   const;
	void SetFRH4Thr( const G4double frh4Thr );
	G4double GetFRH4Thr()  const;
	// Set & get configurations: FRH5
	void SetUseFRH5( const G4bool   useFRH5 );
	G4bool   GetUseFRH5()  const;
	void SetFRH5Mask(const std::vector<G4bool> frh5Mask);
	std::vector<G4bool> GetFRH5Mask() const;
	void SetFRH5DX1( const G4double frh5DX1 );
	G4double GetFRH5DX1()  const;
	void SetFRH5DX2( const G4double frh5DX2 );
	G4double GetFRH5DX2()  const;
	void SetFRH5DY(  const G4double frh5DY  );
	G4double GetFRH5DY()   const;
	void SetFRH5Thi( const G4double frh5Thi );
	G4double GetFRH5Thi()  const;
	void SetFRH5RMin(const G4double frh5RMin);
	G4double GetFRH5RMin() const;
	void SetFRH5RMax(const G4double frh5RMax);
	G4double GetFRH5RMax() const;
	void SetFRH5PR(  const G4double frh5PR  );
	G4double GetFRH5PR()   const;
	void SetFRH5PZ(  const G4double frh5PZ  );
	G4double GetFRH5PZ()   const;
	void SetFRH5Thr( const G4double frh5Thr );
	G4double GetFRH5Thr()  const;
	// Set & get configurations: Wrapping scintillators
	void SetMylThi(  const G4double mylThi  );
	G4double GetMylThi()   const;
	void SetTedThi(  const G4double tedThi  );
	G4double GetTedThi()   const;
	// Set & get configurations: FPC
	void SetFPC12Phi(const G4double fpc12Phi);
	G4double GetFPC12Phi() const;
	void SetFPC34Phi(const G4double fpc34Phi);
	G4double GetFPC34Phi() const;
	void SetFPC56Phi(const G4double fpc56Phi);
	G4double GetFPC56Phi() const;
	void SetFPC78Phi(const G4double fpc78Phi);
	G4double GetFPC78Phi() const;
	void SetFPC12PZ( const G4double fpc12PZ );
	G4double GetFPC12PZ()  const;
	void SetFPC34PZ( const G4double fpc34PZ );
	G4double GetFPC34PZ()  const;
	void SetFPC56PZ( const G4double fpc56PZ );
	G4double GetFPC56PZ()  const;
	void SetFPC78PZ( const G4double fpc78PZ );
	G4double GetFPC78PZ()  const;
	void SetFPCAlThi(const G4double fpcAlThi);
	G4double GetFPCAlThi() const;
	void SetTubeDia( const G4double tubeDia );
	G4double GetTubeDia()  const;
	void SetWallThi( const G4double wallThi );
	G4double GetWallThi()  const;
	void SetWireDia( const G4double wireDia );
	G4double GetWireDia()  const;
	void SetTubeL(   const G4double tubeL   );
	G4double GetTubeL()    const;
	void SetFPCGas(  const G4int    fpcGas  );
	G4int    GetFPCGas()   const;
	void SetFPCThr(  const G4double fpcThr  );
	G4double GetFPCThr()   const;
	void SetFPCVis(  const G4bool   fpcVis  );
	G4bool   GetFPCVis()   const;
	// Set & get configurations: Tigger
	void SetTrigger( const G4int    trigger );
	G4int    GetTrigger()  const;

	// Does the line start with specific word?
	bool StartsWith(const std::string& line, const std::string& token);

  private:
	void SetDefault();

	// Computing resource
	std::string m_SNofTRD;   G4int    m_NofTRD;
	// World
	std::string m_SWorldX;   G4double m_WorldX;
	std::string m_SWorldY;   G4double m_WorldY;
	std::string m_SWorldZ;   G4double m_WorldZ;
	// Input/output
	std::string m_SUseInput; G4bool   m_UseInput;
	std::string m_SOutForm;  G4int    m_OutForm;
	// Beam
	std::string m_SBeamPX;   G4double m_BeamPX;
	std::string m_SBeamPY;   G4double m_BeamPY;
	std::string m_SBeamDX;   G4double m_BeamDX;
	std::string m_SBeamDY;   G4double m_BeamDY;
	std::string m_SBeamKE;   G4double m_BeamKE;
	// Pipe and chamber
	std::string m_SPipeDia;  G4double m_PipeDia;
	std::string m_SVacThi;   G4double m_VacThi;
	// Target
	std::string m_SUseTar;   G4bool   m_UseTar;
	std::string m_STarMat;   G4String m_TarMat;
	std::string m_STarThi;   G4double m_TarThi;
	// FWC1
	std::string m_SUseFWC1;  G4bool   m_UseFWC1;
	std::string m_SFWC1Mask; std::vector<G4bool> m_FWC1Mask;
	std::string m_SFWC1DX1;  G4double m_FWC1DX1;
	std::string m_SFWC1DX2;  G4double m_FWC1DX2;
	std::string m_SFWC1DY;   G4double m_FWC1DY;
	std::string m_SFWC1Thi;  G4double m_FWC1Thi;
	std::string m_SFWC1RMin; G4double m_FWC1RMin;
	std::string m_SFWC1RMax; G4double m_FWC1RMax;
	std::string m_SFWC1PR;   G4double m_FWC1PR;
	std::string m_SFWC1PZ;   G4double m_FWC1PZ;
	std::string m_SFWC1Thr;  G4double m_FWC1Thr;
	// FWC2
	std::string m_SUseFWC2;  G4bool   m_UseFWC2;
	std::string m_SFWC2Mask; std::vector<G4bool> m_FWC2Mask;
	std::string m_SFWC2DX1;  G4double m_FWC2DX1;
	std::string m_SFWC2DX2;  G4double m_FWC2DX2;
	std::string m_SFWC2DY;   G4double m_FWC2DY;
	std::string m_SFWC2Thi;  G4double m_FWC2Thi;
	std::string m_SFWC2RMin; G4double m_FWC2RMin;
	std::string m_SFWC2RMax; G4double m_FWC2RMax;
	std::string m_SFWC2PR;   G4double m_FWC2PR;
	std::string m_SFWC2PZ;   G4double m_FWC2PZ;
	std::string m_SFWC2Thr;  G4double m_FWC2Thr;
	// FTH
	std::string m_SUseFTH;   G4bool   m_UseFTH;
	std::string m_SFTHMask;  std::vector<G4bool> m_FTHMask;
	std::string m_SFTHDX1;   G4double m_FTHDX1;
	std::string m_SFTHDX2;   G4double m_FTHDX2;
	std::string m_SFTHDY;    G4double m_FTHDY;
	std::string m_SFTHThi;   G4double m_FTHThi;
	std::string m_SFTHRMin;  G4double m_FTHRMin;
	std::string m_SFTHRMax;  G4double m_FTHRMax;
	std::string m_SFTHPR;    G4double m_FTHPR;
	std::string m_SFTHPZ;    G4double m_FTHPZ;
	std::string m_SFTHThr;   G4double m_FTHThr;
	// FRH1
	std::string m_SUseFRH1;  G4bool   m_UseFRH1;
	std::string m_SFRH1Mask; std::vector<G4bool> m_FRH1Mask;
	std::string m_SFRH1DX1;  G4double m_FRH1DX1;
	std::string m_SFRH1DX2;  G4double m_FRH1DX2;
	std::string m_SFRH1DY;   G4double m_FRH1DY;
	std::string m_SFRH1Thi;  G4double m_FRH1Thi;
	std::string m_SFRH1RMin; G4double m_FRH1RMin;
	std::string m_SFRH1RMax; G4double m_FRH1RMax;
	std::string m_SFRH1PR;   G4double m_FRH1PR;
	std::string m_SFRH1PZ;   G4double m_FRH1PZ;
	std::string m_SFRH1Thr;  G4double m_FRH1Thr;
	// FRH2
	std::string m_SUseFRH2;  G4bool   m_UseFRH2;
	std::string m_SFRH2Mask; std::vector<G4bool> m_FRH2Mask;
	std::string m_SFRH2DX1;  G4double m_FRH2DX1;
	std::string m_SFRH2DX2;  G4double m_FRH2DX2;
	std::string m_SFRH2DY;   G4double m_FRH2DY;
	std::string m_SFRH2Thi;  G4double m_FRH2Thi;
	std::string m_SFRH2RMin; G4double m_FRH2RMin;
	std::string m_SFRH2RMax; G4double m_FRH2RMax;
	std::string m_SFRH2PR;   G4double m_FRH2PR;
	std::string m_SFRH2PZ;   G4double m_FRH2PZ;
	std::string m_SFRH2Thr;  G4double m_FRH2Thr;
	// FRH3
	std::string m_SUseFRH3;  G4bool   m_UseFRH3;
	std::string m_SFRH3Mask; std::vector<G4bool> m_FRH3Mask;
	std::string m_SFRH3DX1;  G4double m_FRH3DX1;
	std::string m_SFRH3DX2;  G4double m_FRH3DX2;
	std::string m_SFRH3DY;   G4double m_FRH3DY;
	std::string m_SFRH3Thi;  G4double m_FRH3Thi;
	std::string m_SFRH3RMin; G4double m_FRH3RMin;
	std::string m_SFRH3RMax; G4double m_FRH3RMax;
	std::string m_SFRH3PR;   G4double m_FRH3PR;
	std::string m_SFRH3PZ;   G4double m_FRH3PZ;
	std::string m_SFRH3Thr;  G4double m_FRH3Thr;
	// FRH4
	std::string m_SUseFRH4;  G4bool   m_UseFRH4;
	std::string m_SFRH4Mask; std::vector<G4bool> m_FRH4Mask;
	std::string m_SFRH4DX1;  G4double m_FRH4DX1;
	std::string m_SFRH4DX2;  G4double m_FRH4DX2;
	std::string m_SFRH4DY;   G4double m_FRH4DY;
	std::string m_SFRH4Thi;  G4double m_FRH4Thi;
	std::string m_SFRH4RMin; G4double m_FRH4RMin;
	std::string m_SFRH4RMax; G4double m_FRH4RMax;
	std::string m_SFRH4PR;   G4double m_FRH4PR;
	std::string m_SFRH4PZ;   G4double m_FRH4PZ;
	std::string m_SFRH4Thr;  G4double m_FRH4Thr;
	// FRH5
	std::string m_SUseFRH5;  G4bool   m_UseFRH5;
	std::string m_SFRH5Mask; std::vector<G4bool> m_FRH5Mask;
	std::string m_SFRH5DX1;  G4double m_FRH5DX1;
	std::string m_SFRH5DX2;  G4double m_FRH5DX2;
	std::string m_SFRH5DY;   G4double m_FRH5DY;
	std::string m_SFRH5Thi;  G4double m_FRH5Thi;
	std::string m_SFRH5RMin; G4double m_FRH5RMin;
	std::string m_SFRH5RMax; G4double m_FRH5RMax;
	std::string m_SFRH5PR;   G4double m_FRH5PR;
	std::string m_SFRH5PZ;   G4double m_FRH5PZ;
	std::string m_SFRH5Thr;  G4double m_FRH5Thr;
	// Wrapping
	std::string m_SMylThi;   G4double m_MylThi;
	std::string m_STedThi;   G4double m_TedThi;
	// FPC
	std::string m_SFPC12Phi; G4double m_FPC12Phi;
	std::string m_SFPC34Phi; G4double m_FPC34Phi;
	std::string m_SFPC56Phi; G4double m_FPC56Phi;
	std::string m_SFPC78Phi; G4double m_FPC78Phi;
	std::string m_SFPC12PZ;  G4double m_FPC12PZ;
	std::string m_SFPC34PZ;  G4double m_FPC34PZ;
	std::string m_SFPC56PZ;  G4double m_FPC56PZ;
	std::string m_SFPC78PZ;  G4double m_FPC78PZ;
	std::string m_SFPCAlThi; G4double m_FPCAlThi;
	std::string m_STubeDia;  G4double m_TubeDia;
	std::string m_SWireDia;  G4double m_WireDia;
	std::string m_SWallThi;  G4double m_WallThi;
	std::string m_STubeL;    G4double m_TubeL;
	std::string m_SFPCGas;   G4int    m_FPCGas;
	std::string m_SFPCThr;   G4double m_FPCThr;
	std::string m_SFPCVis;   G4bool   m_FPCVis;
	// Trigger
	std::string m_STrigger;  G4int    m_Trigger;
};

#endif
