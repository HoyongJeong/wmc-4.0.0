////////////////////////////////////////////////////////////////////////////////
//   WMC4DetCon.cc                                                            //
//                                                                            //
//   Definitions of WMC4DetCon class's member functions.                      //
// And it describes geometry of simulation.                                   //
//                                                                            //
//                    - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "WMC4DetCon.hh"

//////////////////////////////////////////////////
//   Constructor and destructor                 //
//////////////////////////////////////////////////
WMC4DetCon::WMC4DetCon(WMC4ConMan* CM): m_CM(CM)
{
	ConstructMaterials();
	DefineDimensions();
}

WMC4DetCon::~WMC4DetCon()
{
	DestructMaterials();
}

//////////////////////////////////////////////////
//   Define dimensions                          //
//////////////////////////////////////////////////
void WMC4DetCon::DefineDimensions()
{

	// World dimensions: Laboratory size
	m_LabX     = m_CM -> GetWorldX();   // World x dimension
	m_LabY     = m_CM -> GetWorldY();   // World y dimension
	m_LabZ     = m_CM -> GetWorldZ();   // World z dimension
	// Pipe & chamber
	m_PipeDia  = m_CM -> GetPipeDia();  // Diameter of pipe (not including wall)
	m_VacThi   = m_CM -> GetVacThi();   // Thickness of vacuum chamber wall
	// Target dimensions
	m_TarD     = 50.0 * mm;             // Target diameter
	m_TarT     = m_CM -> GetTarThi();   // Thickness
	// Geometry selectors
	m_UseTar   = m_CM -> GetUseTar();
	m_UseFWC1  = m_CM -> GetUseFWC1();
	m_UseFWC2  = m_CM -> GetUseFWC2();
	m_UseFTH   = m_CM -> GetUseFTH();
	m_UseFRH1  = m_CM -> GetUseFRH1();
	m_UseFRH2  = m_CM -> GetUseFRH2();
	m_UseFRH3  = m_CM -> GetUseFRH3();
	m_UseFRH4  = m_CM -> GetUseFRH4();
	m_UseFRH5  = m_CM -> GetUseFRH5();
	m_FWC1Mask = m_CM -> GetFWC1Mask();
	m_FWC2Mask = m_CM -> GetFWC2Mask();
	m_FTHMask  = m_CM -> GetFTHMask();
	m_FRH1Mask = m_CM -> GetFRH1Mask();
	m_FRH2Mask = m_CM -> GetFRH2Mask();
	m_FRH3Mask = m_CM -> GetFRH3Mask();
	m_FRH4Mask = m_CM -> GetFRH4Mask();
	m_FRH5Mask = m_CM -> GetFRH5Mask();
	// FWC1 dimensions
	m_FWC1DX1  = m_CM -> GetFWC1DX1();
	m_FWC1DY   = m_CM -> GetFWC1DY();
	m_FWC1Thi  = m_CM -> GetFWC1Thi();
	m_FWC1RMin = m_CM -> GetFWC1RMin();
	m_FWC1RMax = m_CM -> GetFWC1RMax();
	m_FWC1PZ   = m_CM -> GetFWC1PZ();
	// FWC2 dimensions
	m_FWC2DX1  = m_CM -> GetFWC2DX1();
	m_FWC2DY   = m_CM -> GetFWC2DY();
	m_FWC2Thi  = m_CM -> GetFWC2Thi();
	m_FWC2RMin = m_CM -> GetFWC2RMin();
	m_FWC2RMax = m_CM -> GetFWC2RMax();
	m_FWC2PZ   = m_CM -> GetFWC2PZ();
	// FTH dimensions
	m_FTHDX1   = m_CM -> GetFTHDX1();
	m_FTHDY    = m_CM -> GetFTHDY();
	m_FTHThi   = m_CM -> GetFTHThi();
	m_FTHRMin  = m_CM -> GetFTHRMin();
	m_FTHRMax  = m_CM -> GetFTHRMax();
	m_FTHPZ    = m_CM -> GetFTHPZ();
	// FRH1 dimensions
	m_FRH1DX1  = m_CM -> GetFRH1DX1();
	m_FRH1DY   = m_CM -> GetFRH1DY();
	m_FRH1Thi  = m_CM -> GetFRH1Thi();
	m_FRH1RMin = m_CM -> GetFRH1RMin();
	m_FRH1RMax = m_CM -> GetFRH1RMax();
	m_FRH1PZ   = m_CM -> GetFRH1PZ();
	// FRH2 dimensions
	m_FRH2DX1  = m_CM -> GetFRH2DX1();
	m_FRH2DY   = m_CM -> GetFRH2DY();
	m_FRH2Thi  = m_CM -> GetFRH2Thi();
	m_FRH2RMin = m_CM -> GetFRH2RMin();
	m_FRH2RMax = m_CM -> GetFRH2RMax();
	m_FRH2PZ   = m_CM -> GetFRH2PZ();
	// FRH3 dimensions
	m_FRH3DX1  = m_CM -> GetFRH3DX1();
	m_FRH3DY   = m_CM -> GetFRH3DY();
	m_FRH3Thi  = m_CM -> GetFRH3Thi();
	m_FRH3RMin = m_CM -> GetFRH3RMin();
	m_FRH3RMax = m_CM -> GetFRH3RMax();	
	m_FRH3RMax = m_CM -> GetFRH3RMax();
	m_FRH3PZ   = m_CM -> GetFRH3PZ();
	// FRH4 dimensions
	m_FRH4DX1  = m_CM -> GetFRH4DX1();
	m_FRH4DY   = m_CM -> GetFRH4DY();
	m_FRH4Thi  = m_CM -> GetFRH4Thi();
	m_FRH4RMin = m_CM -> GetFRH4RMin();
	m_FRH4RMax = m_CM -> GetFRH4RMax();
	m_FRH4PZ   = m_CM -> GetFRH4PZ();
	// FRH5 dimensions
	m_FRH5DX1  = m_CM -> GetFRH5DX1();
	m_FRH5DY   = m_CM -> GetFRH5DY();
	m_FRH5Thi  = m_CM -> GetFRH5Thi();
	m_FRH5RMin = m_CM -> GetFRH5RMin();
	m_FRH5RMax = m_CM -> GetFRH5RMax();
	m_FRH5PZ   = m_CM -> GetFRH5PZ();
	// Wrapping
	m_MylThi   = m_CM -> GetMylThi();
	m_TedThi   = m_CM -> GetTedThi();
	// FPC
	m_FPC12Phi = m_CM -> GetFPC12Phi();
	m_FPC34Phi = m_CM -> GetFPC34Phi();
	m_FPC56Phi = m_CM -> GetFPC56Phi();
	m_FPC78Phi = m_CM -> GetFPC78Phi();
	m_FPC12PZ  = m_CM -> GetFPC12PZ();
	m_FPC34PZ  = m_CM -> GetFPC34PZ();
	m_FPC56PZ  = m_CM -> GetFPC56PZ();
	m_FPC78PZ  = m_CM -> GetFPC78PZ();
	m_FPCAlThi = m_CM -> GetFPCAlThi();
	m_TubeDia  = m_CM -> GetTubeDia();
	m_WallThi  = m_CM -> GetWallThi();
	m_WireDia  = m_CM -> GetWireDia();
	m_TubeL    = m_CM -> GetTubeL();
	m_FPCGas   = m_CM -> GetFPCGas();
	m_FPCVis   = m_CM -> GetFPCVis();
}

//////////////////////////////////////////////////
//   Construct everything                       //
//////////////////////////////////////////////////
G4VPhysicalVolume* WMC4DetCon::Construct()
{
	// World
	m_LabSolid = new G4Box("LabSolid", m_LabX/2, m_LabY/2, m_LabZ/2);
	m_LabLV = new G4LogicalVolume(m_LabSolid, m_AirMat, "LabLV");
	m_LabPV = new G4PVPlacement(0, G4ThreeVector(), "LabPV", m_LabLV, 0, false, 0);

	// Pipe & chamber
	// NOTE: This geometry just has been imported from WMC3.
	// Someday, it is needed to be updated to reflect current setup exactly.
	m_PCFSolid    = new G4Cons("PCFSolid"   , 30.*mm, 33.*mm, 325.*mm, 328.*mm, 455.*mm, 360.*degree, 360.*degree);
	m_PCFVacSolid = new G4Cons("PCFVacSolid",  0.   , 30.*mm,   0.   , 325.*mm, 455.*mm, 360.*degree, 360.*degree);
	m_VCF1Solid    = new G4Tubs("VCF1Solid"   ,  34. * mm,  42. * mm,              5. * mm, 0., 360. * degree);
	m_VCF2Solid    = new G4Tubs("VCF2Solid"   ,  34. * mm,  42. * mm,             13. * mm, 0., 360. * degree);
	m_VCF3Solid    = new G4Tubs("VCF3Solid"   , 351. * mm, 400. * mm,             30. * mm, 0., 360. * degree);
	m_BTFSolid     = new G4Tubs("BTFSolid"    ,  34. * mm,  52. * mm,             14. * mm, 0., 360. * degree);
	m_BTCDSolid    = new G4Tubs("BTCDSolid"   ,  30. * mm,  34. * mm,            96.5 * mm, 0., 360. * degree);
	m_BTCDVacSolid = new G4Tubs("BTCDVacSolid",   0.     ,  30. * mm,            96.5 * mm, 0., 360. * degree);
	m_BTBSolid     = new G4Tubs("BTBSolid"    ,  30. * mm,  34. * mm, m_LabZ/4 - 599. * mm, 0., 360. * degree);
	m_BTBVacSolid  = new G4Tubs("BTBVacSolid" ,   0.     ,  30. * mm, m_LabZ/4 - 599. * mm, 0., 360. * degree);
	m_BBBSolid     = new G4Tubs("BBBSolid"    ,  30. * mm,  34. * mm, m_LabZ/4 +  45. * mm, 0., 360. * degree);
	m_BBBVacSolid  = new G4Tubs("BBBVacSolid" ,   0.     ,  30. * mm, m_LabZ/4 +  45. * mm, 0., 360. * degree);
	m_VCBSolid = new G4Box("VCBSolid", 180.*mm, 2.5*mm, 13.*mm);
	G4double ang1 = 3.391   * degree;
	G4double ang2 = 27.2391 * degree;
	G4double radi = 710.064 * mm;
	m_VCWSolid = new G4Sphere("VCWSolid", radi-.2*mm, radi+.2*mm, 0., 360.*degree, ang1, ang2);
	m_PCFLV     = new G4LogicalVolume(m_PCFSolid    , m_SSMat , "PCFLV"    );
	m_PCFVacLV  = new G4LogicalVolume(m_PCFVacSolid , m_VacMat, "PCFVacLV" );
	m_VCF1LV    = new G4LogicalVolume(m_VCF1Solid   , m_SSMat , "VCF1LV"   );
	m_VCF2LV    = new G4LogicalVolume(m_VCF2Solid   , m_SSMat , "VCF2LV"   );
	m_VCF3LV    = new G4LogicalVolume(m_VCF3Solid   , m_SSMat , "VCF3LV"   );
	m_VCB1LV    = new G4LogicalVolume(m_VCBSolid    , m_SSMat , "VCB1LV"   );
	m_VCB2LV    = new G4LogicalVolume(m_VCBSolid    , m_SSMat , "VCB2LV"   );
	m_VCB3LV    = new G4LogicalVolume(m_VCBSolid    , m_SSMat , "VCB3LV"   );
	m_VCB4LV    = new G4LogicalVolume(m_VCBSolid    , m_SSMat , "VCB4LV"   );
	m_VCWLV     = new G4LogicalVolume(m_VCWSolid    , m_SSMat , "VCWLV"    );
	m_BTFLV     = new G4LogicalVolume(m_BTFSolid    , m_SSMat , "BTFLV"    );
	m_BTCDLV    = new G4LogicalVolume(m_BTCDSolid   , m_SSMat , "BTCDLV"   );
	m_BTCDVacLV = new G4LogicalVolume(m_BTCDVacSolid, m_VacMat, "BTCDVacLV");
	m_BTBLV     = new G4LogicalVolume(m_BTBSolid    , m_SSMat , "BTBLV"    );
	m_BTBVacLV  = new G4LogicalVolume(m_BTBVacSolid , m_VacMat, "BTBVacLV" );
	m_BBBLV     = new G4LogicalVolume(m_BBBSolid    , m_SSMat , "BBBLV"    );
	m_BBBVacLV  = new G4LogicalVolume(m_BBBVacSolid , m_VacMat, "BBBVacLV" );
	m_PCFLV     -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_PCFVacLV  -> SetVisAttributes(new G4VisAttributes(0               ));
	m_VCF1LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCF2LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCF3LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCB1LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCB2LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCB3LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCB4LV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_VCWLV     -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_BTFLV     -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_BTCDLV    -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_BTCDVacLV -> SetVisAttributes(new G4VisAttributes(0               ));
	m_BTBLV     -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_BTBVacLV  -> SetVisAttributes(new G4VisAttributes(0               ));
	m_BBBLV     -> SetVisAttributes(new G4VisAttributes(G4Colour::Gray()));
	m_BBBVacLV  -> SetVisAttributes(new G4VisAttributes(0               ));
	new G4PVPlacement(0, G4ThreeVector(0., 0., 545.*mm         ), "PCFPV"    , m_PCFLV    , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 545.*mm         ), "PCFVacPV" , m_PCFVacLV , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1000.*mm        ), "VCF1PV"   , m_VCF1LV   , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1113.*mm        ), "VCF2PV"   , m_VCF2LV   , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1065.*mm        ), "VCF3PV"   , m_VCF3LV   , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1449.*mm        ), "BTFPV"    , m_BTFLV    , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1101.5*mm       ), "BTCDPV"   , m_BTCDLV   , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1101.5*mm       ), "BTCDVacPV", m_BTCDVacLV, m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., m_LabZ/4+599.*mm), "BTBPV"    , m_BTBLV    , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., m_LabZ/4+599.*mm), "BTBVacPV" , m_BTBVacLV , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 45.*mm-m_LabZ/4 ), "BBBPV"    , m_BBBLV    , m_LabPV, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 45.*mm-m_LabZ/4 ), "BBBVacPV" , m_BBBVacLV , m_LabPV, false, 0);
	// Support bars
	G4ThreeVector VCB1Pos(  222. * mm,     0.     , 1113. * mm);
	G4ThreeVector VCB2Pos(    0.     ,   222. * mm, 1113. * mm);
	G4ThreeVector VCB3Pos(- 222. * mm,     0.     , 1113. * mm);
	G4ThreeVector VCB4Pos(    0.     , - 222. * mm, 1113. * mm);
	m_VCBRot = new G4RotationMatrix();
	m_VCBRot -> rotateZ(90. * degree);
	new G4PVPlacement(0       , VCB1Pos , "VCB1PV", m_VCB1LV, m_LabPV, false, 0);
	new G4PVPlacement(m_VCBRot, VCB2Pos , "VCB2PV", m_VCB2LV, m_LabPV, false, 0);
	new G4PVPlacement(0       , VCB3Pos , "VCB3PV", m_VCB3LV, m_LabPV, false, 0);
	new G4PVPlacement(m_VCBRot, VCB4Pos , "VCB4PV", m_VCB4LV, m_LabPV, false, 0);
	// Window
	m_VCWRot = new G4RotationMatrix();
	m_VCWRot -> rotateX(180. * degree);
	new G4PVPlacement(m_VCWRot, G4ThreeVector(0., 0., 997.5*mm+radi) , "VCWPV", m_VCWLV, m_LabPV, false, 0);


	// Cylindrical target
	if ( m_UseTar )
	{
		m_TarSolid = new G4Tubs("TarSolid", 0., m_TarD/2, m_TarT/2, 0. * degree, 360. * degree);
		m_TarLV = new G4LogicalVolume(m_TarSolid, m_TarMat, "TarLV");
		new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), "TarPV", m_TarLV, m_LabPV, false, 0);
	}

	// Wrapping assistant: Going to be used for Mylar and Tedlar
	// "N" means "Narrow". Due to 5 mm support, some FWC1 elements are narrower
	m_MylAttachRot  = new G4RotationMatrix();
	m_MylAttachRotN = new G4RotationMatrix();
	m_MylAttachRot2 = new G4RotationMatrix();
	m_TedAttachRot  = new G4RotationMatrix();
	m_TedAttachRotN = new G4RotationMatrix();
	m_TedAttachRot2 = new G4RotationMatrix();
	m_MylAttachRot  -> rotateZ(- 15.  * degree);
	m_MylAttachRotN -> rotateZ(- 14.2 * degree);
	m_MylAttachRot2 -> rotateZ(-  7.5 * degree);
	m_TedAttachRot  -> rotateZ(- 15.  * degree);
	m_TedAttachRotN -> rotateZ(- 14.2 * degree);
	m_TedAttachRot2 -> rotateZ(-  7.5 * degree);

	// Ready for hodoscope layers
	G4double rDMin, rDMax;

	// FWC1
	if ( m_UseFWC1 )
	{
		rDMin = 0.5 * m_FWC1DX1 / std::sin(7.5*degree);
		rDMax = rDMin + m_FWC1DY;
		m_FWC1Solid         = new G4Tubs("FWC1Solid"        , rDMin                  , rDMax                  , m_FWC1Thi/2                  , 0., 15. *degree);
		m_FWC1NSolid        = new G4Tubs("FWC1NSolid"       , rDMin                  , rDMax                  , m_FWC1Thi/2                  , 0., 14.2*degree);
		m_FWC1MylBaseSolid  = new G4Tubs("FWC1MylBaseSolid" , rDMin-m_MylThi         , rDMax+m_MylThi         , m_FWC1Thi/2+m_MylThi         , 0., 15. *degree);
		m_FWC1MylBaseNSolid = new G4Tubs("FWC1MylBaseNSolid", rDMin-m_MylThi         , rDMax+m_MylThi         , m_FWC1Thi/2+m_MylThi         , 0., 14.2*degree);
		m_FWC1TedBaseSolid  = new G4Tubs("FWC1TedBaseSolid" , rDMin-m_MylThi-m_TedThi, rDMax+m_MylThi+m_TedThi, m_FWC1Thi/2+m_MylThi+m_TedThi, 0., 15. *degree);
		m_FWC1TedBaseNSolid = new G4Tubs("FWC1TedBaseNSolid", rDMin-m_MylThi-m_TedThi, rDMax+m_MylThi+m_TedThi, m_FWC1Thi/2+m_MylThi+m_TedThi, 0., 14.2*degree);
		m_FWC1MylAttachSolid = new G4Box("FWC1MylAttachSolid", (rDMax-rDMin)/2 + m_MylThi           , m_MylThi/2, m_FWC1Thi/2 + m_MylThi           );
		m_FWC1TedAttachSolid = new G4Box("FWC1TedAttachSolid", (rDMax-rDMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FWC1Thi/2 + m_MylThi + m_TedThi);
		m_FWC1MylSolidTemp1  = new G4SubtractionSolid("FWC1MylSolidTemp1" , m_FWC1MylBaseSolid , m_FWC1Solid        , 0, G4ThreeVector());
		m_FWC1MylNSolidTemp1 = new G4SubtractionSolid("FWC1MylNSolidTemp1", m_FWC1MylBaseNSolid, m_FWC1NSolid       , 0, G4ThreeVector());
		m_FWC1TedSolidTemp1  = new G4SubtractionSolid("FWC1TedSolidTemp1" , m_FWC1TedBaseSolid , m_FWC1MylBaseSolid , 0, G4ThreeVector());
		m_FWC1TedNSolidTemp1 = new G4SubtractionSolid("FWC1TedNSolidTemp1", m_FWC1TedBaseNSolid, m_FWC1MylBaseNSolid, 0, G4ThreeVector());
		G4ThreeVector mylAttach( (rDMax+rDMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector mylAttachN((rDMax+rDMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach( (rDMax+rDMin)/2, - 0.5*m_TedThi, 0.);
		G4ThreeVector tedAttachN((rDMax+rDMin)/2, - 0.5*m_TedThi, 0.);
		m_FWC1MylSolidTemp2  = new G4UnionSolid("FWC1MylSolidTemp2" , m_FWC1MylSolidTemp1 , m_FWC1MylAttachSolid , 0, mylAttach );
		m_FWC1MylNSolidTemp2 = new G4UnionSolid("FWC1MylNSolidTemp2", m_FWC1MylNSolidTemp1, m_FWC1MylAttachSolid , 0, mylAttachN);
		m_FWC1TedSolidTemp2  = new G4UnionSolid("FWC1TedSolidTemp2" , m_FWC1TedSolidTemp1 , m_FWC1TedAttachSolid , 0, tedAttach );
		m_FWC1TedNSolidTemp2 = new G4UnionSolid("FWC1TedNSolidTemp2", m_FWC1TedNSolidTemp1, m_FWC1TedAttachSolid , 0, tedAttachN);
		mylAttach  += G4ThreeVector(0., m_MylThi, 0.);
		mylAttachN += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach  += G4ThreeVector(0., m_TedThi, 0.);
		tedAttachN += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ( 15.  * degree);
		mylAttachN.rotateZ(14.2 * degree);
		tedAttach.rotateZ( 15.  * degree);
		tedAttachN.rotateZ(14.2 * degree);
		m_FWC1MylSolid  = new G4UnionSolid("FWC1MylSolid" , m_FWC1MylSolidTemp2 , m_FWC1MylAttachSolid , m_MylAttachRot , mylAttach );
		m_FWC1MylNSolid = new G4UnionSolid("FWC1MylNSolid", m_FWC1MylNSolidTemp2, m_FWC1MylAttachSolid , m_MylAttachRotN, mylAttachN);
		m_FWC1TedSolid  = new G4UnionSolid("FWC1TedSolid" , m_FWC1TedSolidTemp2 , m_FWC1TedAttachSolid , m_TedAttachRot , tedAttach );
		m_FWC1TedNSolid = new G4UnionSolid("FWC1TedNSolid", m_FWC1TedNSolidTemp2, m_FWC1TedAttachSolid , m_TedAttachRotN, tedAttachN);
		// Due to 5 mm support, these elements are narrower by 2.6 mm
		m_SmallFWC1El[ 0] = 1;
		m_SmallFWC1El[ 5] = 1;
		m_SmallFWC1El[ 6] = 1;
		m_SmallFWC1El[11] = 1;
		m_SmallFWC1El[12] = 1;
		m_SmallFWC1El[17] = 1;
		m_SmallFWC1El[18] = 1;
		m_SmallFWC1El[23] = 1;
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FWC1Mask[i] ) continue;

			// Logical volume
			m_FWC1LVName[i]    = "FWC1LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FWC1MylLVName[i] = "FWC1MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FWC1TedLVName[i] = "FWC1TedLV_" + G4UIcommand::ConvertToString(i+1);
			if ( m_SmallFWC1El[i] )
			{
				m_FWC1LV[i]    = new G4LogicalVolume(m_FWC1NSolid   , m_PlaMat, m_FWC1LVName[i].data()   );
				m_FWC1MylLV[i] = new G4LogicalVolume(m_FWC1MylNSolid, m_MylMat, m_FWC1MylLVName[i].data());
				m_FWC1TedLV[i] = new G4LogicalVolume(m_FWC1TedNSolid, m_TedMat, m_FWC1TedLVName[i].data());
			}
			else
			{
				m_FWC1LV[i]    = new G4LogicalVolume(m_FWC1Solid    , m_PlaMat, m_FWC1LVName[i].data()   );
				m_FWC1MylLV[i] = new G4LogicalVolume(m_FWC1MylSolid , m_MylMat, m_FWC1MylLVName[i].data());
				m_FWC1TedLV[i] = new G4LogicalVolume(m_FWC1TedSolid , m_TedMat, m_FWC1TedLVName[i].data());
			}
//			m_FWC1LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
//			m_FWC1MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
			m_FWC1TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));

			// Positioning of FWC1 element
			G4ThreeVector dr(m_FWC1RMin - rDMin, 0., 0.);
			dr += G4ThreeVector((m_FWC1RMax-m_FWC1RMin) * (1. - std::cos(20. * degree)), 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FWC1Pos[i] = G4ThreeVector(0., 0., m_FWC1PZ - 0.5*(m_FWC1DY+m_FWC1DX1)*std::sin(20. * degree)) + dr;
			m_FWC1Rot[i] = new G4RotationMatrix();
			m_FWC1Rot[i] -> rotateZ(- 15. * i * degree);
			if ( i % 6 == 0 ) m_FWC1Rot[i] -> rotateZ(- 0.8 * degree);
			m_FWC1Rot[i] -> rotateZ(-  7.5    * degree);
			m_FWC1Rot[i] -> rotateY(   20.    * degree);
			m_FWC1Rot[i] -> rotateZ(    7.5   * degree);

			// Place it
			m_FWC1PVName[i]    = "FWC1PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FWC1MylPVName[i] = "FWC1MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FWC1TedPVName[i] = "FWC1TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FWC1Rot[i], m_FWC1Pos[i], m_FWC1PVName[i]   , m_FWC1LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FWC1Rot[i], m_FWC1Pos[i], m_FWC1MylPVName[i], m_FWC1MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FWC1Rot[i], m_FWC1Pos[i], m_FWC1TedPVName[i], m_FWC1TedLV[i], m_LabPV, false, 0);
		}
	}

	// FWC2
	if ( m_UseFWC2 )
	{
		rDMin = 0.5 * m_FWC2DX1 / std::sin(7.5*degree);
		rDMax = rDMin + m_FWC2DY;
		m_FWC2Solid        = new G4Tubs("FWC2Solid"       , rDMin                  , rDMax                  , m_FWC2Thi/2                  , 0., 15.*degree);
		m_FWC2MylBaseSolid = new G4Tubs("FWC2MylBaseSolid", rDMin-m_MylThi         , rDMax+m_MylThi         , m_FWC2Thi/2+m_MylThi         , 0., 15.*degree);
		m_FWC2TedBaseSolid = new G4Tubs("FWC2TedBaseSolid", rDMin-m_MylThi-m_TedThi, rDMax+m_MylThi+m_TedThi, m_FWC2Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FWC2MylAttachSolid = new G4Box("FWC2MylAttachSolid", (rDMax-rDMin)/2 + m_MylThi           , m_MylThi/2, m_FWC2Thi/2 + m_MylThi           );
		m_FWC2TedAttachSolid = new G4Box("FWC2TedAttachSolid", (rDMax-rDMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FWC2Thi/2 + m_MylThi + m_TedThi);
		m_FWC2MylSolidTemp1 = new G4SubtractionSolid("FWC2MylSolidTemp1", m_FWC2MylBaseSolid, m_FWC2Solid       , 0, G4ThreeVector(0., 0., 0.));
		m_FWC2TedSolidTemp1 = new G4SubtractionSolid("FWC2TedSolidTemp1", m_FWC2TedBaseSolid, m_FWC2MylBaseSolid, 0, G4ThreeVector(0., 0., 0.));
		G4ThreeVector mylAttach((rDMax+rDMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((rDMax+rDMin)/2, - 0.5*m_TedThi, 0.);
		m_FWC2MylSolidTemp2 = new G4UnionSolid("FWC2MylSolidTemp2", m_FWC2MylSolidTemp1, m_FWC2MylAttachSolid , 0, mylAttach);
		m_FWC2TedSolidTemp2 = new G4UnionSolid("FWC2TedSolidTemp2", m_FWC2TedSolidTemp1, m_FWC2TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FWC2MylSolid = new G4UnionSolid("FWC2MylSolid", m_FWC2MylSolidTemp2, m_FWC2MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FWC2TedSolid = new G4UnionSolid("FWC2TedSolid", m_FWC2TedSolidTemp2, m_FWC2TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FWC2Mask[i] ) continue;

			// Logical volume
			m_FWC2LVName[i]    = "FWC2LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FWC2MylLVName[i] = "FWC2MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FWC2TedLVName[i] = "FWC2TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FWC2LV[i]    = new G4LogicalVolume(m_FWC2Solid   , m_PlaMat, m_FWC2LVName[i].data()   );
			m_FWC2MylLV[i] = new G4LogicalVolume(m_FWC2MylSolid, m_MylMat, m_FWC2MylLVName[i].data());
			m_FWC2TedLV[i] = new G4LogicalVolume(m_FWC2TedSolid, m_TedMat, m_FWC2TedLVName[i].data());
//			m_FWC2LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
//			m_FWC2MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
			m_FWC2TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));

			// Positioning of FWC2 element
			G4ThreeVector dr(m_FWC2RMin - rDMin, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FWC2Pos[i] = G4ThreeVector(0., 0., m_FWC2PZ) + dr;
			m_FWC2Rot[i] = new G4RotationMatrix();
			m_FWC2Rot[i] -> rotateZ((- 15.*i + 7.5) * degree);

			// Place it
			m_FWC2PVName[i]    = "FWC2PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FWC2MylPVName[i] = "FWC2MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FWC2TedPVName[i] = "FWC2TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FWC2Rot[i], m_FWC2Pos[i], m_FWC2PVName[i]   , m_FWC2LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FWC2Rot[i], m_FWC2Pos[i], m_FWC2MylPVName[i], m_FWC2MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FWC2Rot[i], m_FWC2Pos[i], m_FWC2TedPVName[i], m_FWC2TedLV[i], m_LabPV, false, 0);
		}
	}

	// FTH
	if ( m_UseFTH )
	{
		rDMin = 0.5 * m_FTHDX1 / std::sin(3.75*degree);
		rDMax = rDMin + m_FTHDY;
		m_FTHSolid        = new G4Tubs("FTHSolid"       , rDMin                  , rDMax                  , m_FTHThi/2                  , 0., 7.5*degree);
		m_FTHMylBaseSolid = new G4Tubs("FTHMylBaseSolid", rDMin-m_MylThi         , rDMax+m_MylThi         , m_FTHThi/2+m_MylThi         , 0., 7.5*degree);
		m_FTHTedBaseSolid = new G4Tubs("FTHTedBaseSolid", rDMin-m_MylThi-m_TedThi, rDMax+m_MylThi+m_TedThi, m_FTHThi/2+m_MylThi+m_TedThi, 0., 7.5*degree);
		m_FTHMylAttachSolid = new G4Box("FTHMylAttachSolid", (rDMax-rDMin)/2+m_MylThi         , m_MylThi/2, m_FTHThi/2+m_MylThi         );
		m_FTHTedAttachSolid = new G4Box("FTHTedAttachSolid", (rDMax-rDMin)/2+m_MylThi+m_TedThi, m_TedThi/2, m_FTHThi/2+m_MylThi+m_TedThi);
		m_FTHMylSolidTemp1 = new G4SubtractionSolid("FTHMylSolidTemp1", m_FTHMylBaseSolid, m_FTHSolid       , 0, G4ThreeVector());
		m_FTHTedSolidTemp1 = new G4SubtractionSolid("FTHTedSolidTemp1", m_FTHTedBaseSolid, m_FTHMylBaseSolid, 0, G4ThreeVector());
		G4ThreeVector mylAttach((rDMax+rDMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((rDMax+rDMin)/2, - 0.5*m_TedThi, 0.);
		m_FTHMylSolidTemp2 = new G4UnionSolid("FTHMylSolidTemp2", m_FTHMylSolidTemp1, m_FTHMylAttachSolid , 0, mylAttach);
		m_FTHTedSolidTemp2 = new G4UnionSolid("FTHTedSolidTemp2", m_FTHTedSolidTemp1, m_FTHTedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(7.5 * degree);
		tedAttach.rotateZ(7.5 * degree);
		m_FTHMylSolid = new G4UnionSolid("FTHMylSolid", m_FTHMylSolidTemp2, m_FTHMylAttachSolid , m_MylAttachRot2, mylAttach);
		m_FTHTedSolid = new G4UnionSolid("FTHTedSolid", m_FTHTedSolidTemp2, m_FTHTedAttachSolid , m_TedAttachRot2, tedAttach);
		for ( int i = 0; i < 48; i++ )
		{
			if ( !m_FTHMask[i] ) continue;

			// Logical volume
			m_FTHLVName[i]    = "FTHLV_"    + G4UIcommand::ConvertToString(i+1);
			m_FTHMylLVName[i] = "FTHMylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FTHTedLVName[i] = "FTHTedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FTHLV[i]    = new G4LogicalVolume(m_FTHSolid   , m_PlaMat, m_FTHLVName[i].data()   );
			m_FTHMylLV[i] = new G4LogicalVolume(m_FTHMylSolid, m_MylMat, m_FTHMylLVName[i].data());
			m_FTHTedLV[i] = new G4LogicalVolume(m_FTHTedSolid, m_TedMat, m_FTHTedLVName[i].data());
//			m_FTHLV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
//			m_FTHMylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
			m_FTHTedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));

			// Positioning of FTH element
//			G4ThreeVector dr(m_FTHRMin - rDMin, 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((3.75+7.5*i) * degree);
			m_FTHPos[i] = G4ThreeVector(0., 0., m_FTHPZ) + dr;
			m_FTHRot[i] = new G4RotationMatrix();
			m_FTHRot[i] -> rotateZ(- 7.5 * i * degree);

			// Place it
			m_FTHPVName[i]    = "FTHPV_"    + G4UIcommand::ConvertToString(i+1);
			m_FTHMylPVName[i] = "FTHMylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FTHTedPVName[i] = "FTHTedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FTHRot[i], m_FTHPos[i], m_FTHPVName[i]   , m_FTHLV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FTHRot[i], m_FTHPos[i], m_FTHMylPVName[i], m_FTHMylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FTHRot[i], m_FTHPos[i], m_FTHTedPVName[i], m_FTHTedLV[i], m_LabPV, false, 0);
		}
	}

	// FRH1
	if ( m_UseFRH1 )
	{
		rDMin = 0.5 * m_FRH1DX1 / std::sin(7.5*degree);
		rDMax = rDMin + m_FRH1DY;
		m_FRH1Solid        = new G4Tubs("FRH1Solid"       , m_FRH1RMin                  , m_FRH1RMax                  , m_FRH1Thi/2                  , 0., 15.*degree);
		m_FRH1MylBaseSolid = new G4Tubs("FRH1MylBaseSolid", m_FRH1RMin-m_MylThi         , m_FRH1RMax+m_MylThi         , m_FRH1Thi/2+m_MylThi         , 0., 15.*degree);
		m_FRH1TedBaseSolid = new G4Tubs("FRH1TedBaseSolid", m_FRH1RMin-m_MylThi-m_TedThi, m_FRH1RMax+m_MylThi+m_TedThi, m_FRH1Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FRH1MylAttachSolid = new G4Box("FRH1MylAttachSolid", (m_FRH1RMax - m_FRH1RMin)/2 + m_MylThi           , m_MylThi/2, m_FRH1Thi/2 + m_MylThi           );
		m_FRH1TedAttachSolid = new G4Box("FRH1TedAttachSolid", (m_FRH1RMax - m_FRH1RMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FRH1Thi/2 + m_MylThi + m_TedThi);
		m_FRH1MylSolidTemp1 = new G4SubtractionSolid("FRH1MylSolidTemp1", m_FRH1MylBaseSolid, m_FRH1Solid       , 0, G4ThreeVector(0., 0., 0.));
		m_FRH1TedSolidTemp1 = new G4SubtractionSolid("FRH1TedSolidTemp1", m_FRH1TedBaseSolid, m_FRH1MylBaseSolid, 0, G4ThreeVector(0., 0., 0.));
		G4ThreeVector mylAttach((m_FRH1RMax + m_FRH1RMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((m_FRH1RMax + m_FRH1RMin)/2, - 0.5*m_TedThi, 0.);
		m_FRH1MylSolidTemp2 = new G4UnionSolid("FRH1MylSolidTemp2", m_FRH1MylSolidTemp1, m_FRH1MylAttachSolid , 0, mylAttach);
		m_FRH1TedSolidTemp2 = new G4UnionSolid("FRH1TedSolidTemp2", m_FRH1TedSolidTemp1, m_FRH1TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FRH1MylSolid = new G4UnionSolid("FRH1MylSolid", m_FRH1MylSolidTemp2, m_FRH1MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FRH1TedSolid = new G4UnionSolid("FRH1TedSolid", m_FRH1TedSolidTemp2, m_FRH1TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FRH1Mask[i] ) continue;

			// Logical volume
			m_FRH1LVName[i]    = "FRH1LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH1MylLVName[i] = "FRH1MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH1TedLVName[i] = "FRH1TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH1LV[i]    = new G4LogicalVolume(m_FRH1Solid   , m_PlaMat, m_FRH1LVName[i].data()   );
			m_FRH1MylLV[i] = new G4LogicalVolume(m_FRH1MylSolid, m_MylMat, m_FRH1MylLVName[i].data());
			m_FRH1TedLV[i] = new G4LogicalVolume(m_FRH1TedSolid, m_TedMat, m_FRH1TedLVName[i].data());
//			m_FRH1LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
//			m_FRH1MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
			m_FRH1TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

			// Positioning of FRH1 element
//			G4ThreeVector dr(m_FRH1RMin - m_FRH1DX1 / 2.0 / std::sin(7.5 * degree), 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FRH1Pos[i] = G4ThreeVector(0., 0., m_FRH1PZ) + dr;
			m_FRH1Rot[i] = new G4RotationMatrix();
			m_FRH1Rot[i] -> rotateZ(- 15. * i * degree);

			// Place it
			m_FRH1PVName[i]    = "FRH1PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH1MylPVName[i] = "FRH1MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH1TedPVName[i] = "FRH1TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FRH1Rot[i], m_FRH1Pos[i], m_FRH1PVName[i]   , m_FRH1LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FRH1Rot[i], m_FRH1Pos[i], m_FRH1MylPVName[i], m_FRH1MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FRH1Rot[i], m_FRH1Pos[i], m_FRH1TedPVName[i], m_FRH1TedLV[i], m_LabPV, false, 0);
		}
	}

	// FRH2
	if ( m_UseFRH2 )
	{
		m_FRH2Solid        = new G4Tubs("FRH2Solid"       , m_FRH2RMin                  , m_FRH2RMax                  , m_FRH2Thi/2                  , 0., 15.*degree);
		m_FRH2MylBaseSolid = new G4Tubs("FRH2MylBaseSolid", m_FRH2RMin-m_MylThi         , m_FRH2RMax+m_MylThi         , m_FRH2Thi/2+m_MylThi         , 0., 15.*degree);
		m_FRH2TedBaseSolid = new G4Tubs("FRH2TedBaseSolid", m_FRH2RMin-m_MylThi-m_TedThi, m_FRH2RMax+m_MylThi+m_TedThi, m_FRH2Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FRH2MylAttachSolid = new G4Box("FRH2MylAttachSolid", (m_FRH2RMax - m_FRH2RMin)/2 + m_MylThi           , m_MylThi/2, m_FRH2Thi/2 + m_MylThi           );
		m_FRH2TedAttachSolid = new G4Box("FRH2TedAttachSolid", (m_FRH2RMax - m_FRH2RMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FRH2Thi/2 + m_MylThi + m_TedThi);
		m_FRH2MylSolidTemp1 = new G4SubtractionSolid("FRH2MylSolidTemp1", m_FRH2MylBaseSolid, m_FRH2Solid       , 0, G4ThreeVector(0., 0., 0.));
		m_FRH2TedSolidTemp1 = new G4SubtractionSolid("FRH2TedSolidTemp1", m_FRH2TedBaseSolid, m_FRH2MylBaseSolid, 0, G4ThreeVector(0., 0., 0.));
		G4ThreeVector mylAttach((m_FRH2RMax + m_FRH2RMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((m_FRH2RMax + m_FRH2RMin)/2, - 0.5*m_TedThi, 0.);
		m_FRH2MylSolidTemp2 = new G4UnionSolid("FRH2MylSolidTemp2", m_FRH2MylSolidTemp1, m_FRH2MylAttachSolid , 0, mylAttach);
		m_FRH2TedSolidTemp2 = new G4UnionSolid("FRH2TedSolidTemp2", m_FRH2TedSolidTemp1, m_FRH2TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FRH2MylSolid = new G4UnionSolid("FRH2MylSolid", m_FRH2MylSolidTemp2, m_FRH2MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FRH2TedSolid = new G4UnionSolid("FRH2TedSolid", m_FRH2TedSolidTemp2, m_FRH2TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FRH2Mask[i] ) continue;

			// Logical volume
			m_FRH2LVName[i]    = "FRH2LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH2MylLVName[i] = "FRH2MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH2TedLVName[i] = "FRH2TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH2LV[i]    = new G4LogicalVolume(m_FRH2Solid   , m_PlaMat, m_FRH2LVName[i].data()   );
			m_FRH2MylLV[i] = new G4LogicalVolume(m_FRH2MylSolid, m_MylMat, m_FRH2MylLVName[i].data());
			m_FRH2TedLV[i] = new G4LogicalVolume(m_FRH2TedSolid, m_TedMat, m_FRH2TedLVName[i].data());
//			m_FRH2LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
//			m_FRH2MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
			m_FRH2TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

			// Positioning of FRH2 element
//			G4ThreeVector dr(m_FRH2RMin - m_FRH2DX1 / 2.0 / std::sin(7.5 * degree), 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FRH2Pos[i] = G4ThreeVector(0., 0., m_FRH2PZ) + dr;
			m_FRH2Rot[i] = new G4RotationMatrix();
			m_FRH2Rot[i] -> rotateZ(- 15. * i * degree);

			// Place it
			m_FRH2PVName[i]    = "FRH2PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH2MylPVName[i] = "FRH2MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH2TedPVName[i] = "FRH2TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FRH2Rot[i], m_FRH2Pos[i], m_FRH2PVName[i]   , m_FRH2LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FRH2Rot[i], m_FRH2Pos[i], m_FRH2MylPVName[i], m_FRH2MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FRH2Rot[i], m_FRH2Pos[i], m_FRH2TedPVName[i], m_FRH2TedLV[i], m_LabPV, false, 0);
		}
	}

	// FRH3
	if ( m_UseFRH3 )
	{
		m_FRH3Solid        = new G4Tubs("FRH3Solid"       , m_FRH3RMin                  , m_FRH3RMax                  , m_FRH3Thi/2                  , 0., 15.*degree);
		m_FRH3MylBaseSolid = new G4Tubs("FRH3MylBaseSolid", m_FRH3RMin-m_MylThi         , m_FRH3RMax+m_MylThi         , m_FRH3Thi/2+m_MylThi         , 0., 15.*degree);
		m_FRH3TedBaseSolid = new G4Tubs("FRH3TedBaseSolid", m_FRH3RMin-m_MylThi-m_TedThi, m_FRH3RMax+m_MylThi+m_TedThi, m_FRH3Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FRH3MylAttachSolid = new G4Box("FRH3MylAttachSolid", (m_FRH3RMax - m_FRH3RMin)/2 + m_MylThi           , m_MylThi/2, m_FRH3Thi/2 + m_MylThi           );
		m_FRH3TedAttachSolid = new G4Box("FRH3TedAttachSolid", (m_FRH3RMax - m_FRH3RMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FRH3Thi/2 + m_MylThi + m_TedThi);
		m_FRH3MylSolidTemp1 = new G4SubtractionSolid("FRH3MylSolidTemp1", m_FRH3MylBaseSolid, m_FRH3Solid       , 0, G4ThreeVector(0., 0., 0.));
		m_FRH3TedSolidTemp1 = new G4SubtractionSolid("FRH3TedSolidTemp1", m_FRH3TedBaseSolid, m_FRH3MylBaseSolid, 0, G4ThreeVector(0., 0., 0.));
		G4ThreeVector mylAttach((m_FRH3RMax + m_FRH3RMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((m_FRH3RMax + m_FRH3RMin)/2, - 0.5*m_TedThi, 0.);
		m_FRH3MylSolidTemp2 = new G4UnionSolid("FRH3MylSolidTemp2", m_FRH3MylSolidTemp1, m_FRH3MylAttachSolid , 0, mylAttach);
		m_FRH3TedSolidTemp2 = new G4UnionSolid("FRH3TedSolidTemp2", m_FRH3TedSolidTemp1, m_FRH3TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FRH3MylSolid = new G4UnionSolid("FRH3MylSolid", m_FRH3MylSolidTemp2, m_FRH3MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FRH3TedSolid = new G4UnionSolid("FRH3TedSolid", m_FRH3TedSolidTemp2, m_FRH3TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FRH3Mask[i] ) continue;

			// Logical volume
			m_FRH3LVName[i]    = "FRH3LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH3MylLVName[i] = "FRH3MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH3TedLVName[i] = "FRH3TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH3LV[i]    = new G4LogicalVolume(m_FRH3Solid   , m_PlaMat, m_FRH3LVName[i].data()   );
			m_FRH3MylLV[i] = new G4LogicalVolume(m_FRH3MylSolid, m_MylMat, m_FRH3MylLVName[i].data());
			m_FRH3TedLV[i] = new G4LogicalVolume(m_FRH3TedSolid, m_TedMat, m_FRH3TedLVName[i].data());
//			m_FRH3LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
//			m_FRH3MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
			m_FRH3TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

			// Positioning of FRH3 element
//			G4ThreeVector dr(m_FRH3RMin - m_FRH3DX1 / 2.0 / std::sin(7.5 * degree), 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FRH3Pos[i] = G4ThreeVector(0., 0., m_FRH3PZ) + dr;
			m_FRH3Rot[i] = new G4RotationMatrix();
			m_FRH3Rot[i] -> rotateZ(- 15. * i * degree);

			// Place it
			m_FRH3PVName[i]    = "FRH3PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH3MylPVName[i] = "FRH3MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH3TedPVName[i] = "FRH3TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FRH3Rot[i], m_FRH3Pos[i], m_FRH3PVName[i]   , m_FRH3LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FRH3Rot[i], m_FRH3Pos[i], m_FRH3MylPVName[i], m_FRH3MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FRH3Rot[i], m_FRH3Pos[i], m_FRH3TedPVName[i], m_FRH3TedLV[i], m_LabPV, false, 0);
		}
	}

	// FRH4
	if ( m_UseFRH4 )
	{
		m_FRH4Solid        = new G4Tubs("FRH4Solid"       , m_FRH4RMin                  , m_FRH4RMax                  , m_FRH4Thi/2                  , 0., 15.*degree);
		m_FRH4MylBaseSolid = new G4Tubs("FRH4MylBaseSolid", m_FRH4RMin-m_MylThi         , m_FRH4RMax+m_MylThi         , m_FRH4Thi/2+m_MylThi         , 0., 15.*degree);
		m_FRH4TedBaseSolid = new G4Tubs("FRH4TedBaseSolid", m_FRH4RMin-m_MylThi-m_TedThi, m_FRH4RMax+m_MylThi+m_TedThi, m_FRH4Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FRH4MylAttachSolid = new G4Box("FRH4MylAttachSolid", (m_FRH4RMax - m_FRH4RMin)/2 + m_MylThi           , m_MylThi/2, m_FRH4Thi/2 + m_MylThi           );
		m_FRH4TedAttachSolid = new G4Box("FRH4TedAttachSolid", (m_FRH4RMax - m_FRH4RMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FRH4Thi/2 + m_MylThi + m_TedThi);
		m_FRH4MylSolidTemp1 = new G4SubtractionSolid("FRH4MylSolidTemp1", m_FRH4MylBaseSolid, m_FRH4Solid       , 0, G4ThreeVector(0., 0., 0.));
		m_FRH4TedSolidTemp1 = new G4SubtractionSolid("FRH4TedSolidTemp1", m_FRH4TedBaseSolid, m_FRH4MylBaseSolid, 0, G4ThreeVector(0., 0., 0.));
		G4ThreeVector mylAttach((m_FRH4RMax + m_FRH4RMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((m_FRH4RMax + m_FRH4RMin)/2, - 0.5*m_TedThi, 0.);
		m_FRH4MylSolidTemp2 = new G4UnionSolid("FRH4MylSolidTemp2", m_FRH4MylSolidTemp1, m_FRH4MylAttachSolid , 0, mylAttach);
		m_FRH4TedSolidTemp2 = new G4UnionSolid("FRH4TedSolidTemp2", m_FRH4TedSolidTemp1, m_FRH4TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FRH4MylSolid = new G4UnionSolid("FRH4MylSolid", m_FRH4MylSolidTemp2, m_FRH4MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FRH4TedSolid = new G4UnionSolid("FRH4TedSolid", m_FRH4TedSolidTemp2, m_FRH4TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FRH4Mask[i] ) continue;

			// Logical volume
			m_FRH4LVName[i]    = "FRH4LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH4MylLVName[i] = "FRH4MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH4TedLVName[i] = "FRH4TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH4LV[i]    = new G4LogicalVolume(m_FRH4Solid   , m_PlaMat, m_FRH4LVName[i].data()   );
			m_FRH4MylLV[i] = new G4LogicalVolume(m_FRH4MylSolid, m_MylMat, m_FRH4MylLVName[i].data());
			m_FRH4TedLV[i] = new G4LogicalVolume(m_FRH4TedSolid, m_TedMat, m_FRH4TedLVName[i].data());
//			m_FRH4LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
//			m_FRH4MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
			m_FRH4TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

			// Positioning of FRH4 element
//			G4ThreeVector dr(m_FRH4RMin - m_FRH4DX1 / 2.0 / std::sin(7.5 * degree), 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FRH4Pos[i] = G4ThreeVector(0., 0., m_FRH4PZ) + dr;
			m_FRH4Rot[i] = new G4RotationMatrix();
			m_FRH4Rot[i] -> rotateZ(- 15. * i * degree);

			// Place it
			m_FRH4PVName[i]    = "FRH4PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH4MylPVName[i] = "FRH4MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH4TedPVName[i] = "FRH4TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FRH4Rot[i], m_FRH4Pos[i], m_FRH4PVName[i]   , m_FRH4LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FRH4Rot[i], m_FRH4Pos[i], m_FRH4MylPVName[i], m_FRH4MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FRH4Rot[i], m_FRH4Pos[i], m_FRH4TedPVName[i], m_FRH4TedLV[i], m_LabPV, false, 0);
		}
	}

	// FRH5
	if ( m_UseFRH5 )
	{
		m_FRH5Solid        = new G4Tubs("FRH5Solid"       , m_FRH5RMin                  , m_FRH5RMax                  , m_FRH5Thi/2                  , 0., 15.*degree);
		m_FRH5MylBaseSolid = new G4Tubs("FRH5MylBaseSolid", m_FRH5RMin-m_MylThi         , m_FRH5RMax+m_MylThi         , m_FRH5Thi/2+m_MylThi         , 0., 15.*degree);
		m_FRH5TedBaseSolid = new G4Tubs("FRH5TedBaseSolid", m_FRH5RMin-m_MylThi-m_TedThi, m_FRH5RMax+m_MylThi+m_TedThi, m_FRH5Thi/2+m_MylThi+m_TedThi, 0., 15.*degree);
		m_FRH5MylAttachSolid = new G4Box("FRH5MylAttachSolid", (m_FRH5RMax - m_FRH5RMin)/2 + m_MylThi           , m_MylThi/2, m_FRH5Thi/2 + m_MylThi           );
		m_FRH5TedAttachSolid = new G4Box("FRH5TedAttachSolid", (m_FRH5RMax - m_FRH5RMin)/2 + m_MylThi + m_TedThi, m_TedThi/2, m_FRH5Thi/2 + m_MylThi + m_TedThi);
		m_FRH5MylSolidTemp1 = new G4SubtractionSolid("FRH5MylSolidTemp1", m_FRH5MylBaseSolid, m_FRH5Solid       , 0, G4ThreeVector());
		m_FRH5TedSolidTemp1 = new G4SubtractionSolid("FRH5TedSolidTemp1", m_FRH5TedBaseSolid, m_FRH5MylBaseSolid, 0, G4ThreeVector());
		G4ThreeVector mylAttach((m_FRH5RMax + m_FRH5RMin)/2, - 0.5*m_MylThi, 0.);
		G4ThreeVector tedAttach((m_FRH5RMax + m_FRH5RMin)/2, - 0.5*m_TedThi, 0.);
		m_FRH5MylSolidTemp2 = new G4UnionSolid("FRH5MylSolidTemp2", m_FRH5MylSolidTemp1, m_FRH5MylAttachSolid , 0, mylAttach);
		m_FRH5TedSolidTemp2 = new G4UnionSolid("FRH5TedSolidTemp2", m_FRH5TedSolidTemp1, m_FRH5TedAttachSolid , 0, tedAttach);
		mylAttach += G4ThreeVector(0., m_MylThi, 0.);
		tedAttach += G4ThreeVector(0., m_TedThi, 0.);
		mylAttach.rotateZ(15. * degree);
		tedAttach.rotateZ(15. * degree);
		m_FRH5MylSolid = new G4UnionSolid("FRH5MylSolid", m_FRH5MylSolidTemp2, m_FRH5MylAttachSolid , m_MylAttachRot, mylAttach);
		m_FRH5TedSolid = new G4UnionSolid("FRH5TedSolid", m_FRH5TedSolidTemp2, m_FRH5TedAttachSolid , m_TedAttachRot, tedAttach);
		for ( int i = 0; i < 24; i++ )
		{
			if ( !m_FRH5Mask[i] ) continue;

			// Logical volume
			m_FRH5LVName[i]    = "FRH5LV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH5MylLVName[i] = "FRH5MylLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH5TedLVName[i] = "FRH5TedLV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH5LV[i]    = new G4LogicalVolume(m_FRH5Solid   , m_PlaMat, m_FRH5LVName[i].data()   );
			m_FRH5MylLV[i] = new G4LogicalVolume(m_FRH5MylSolid, m_MylMat, m_FRH5MylLVName[i].data());
			m_FRH5TedLV[i] = new G4LogicalVolume(m_FRH5TedSolid, m_TedMat, m_FRH5TedLVName[i].data());
//			m_FRH5LV[i] i  -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
//			m_FRH5MylLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
			m_FRH5TedLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));

			// Positioning of FRH5 element
//			G4ThreeVector dr(m_FRH5RMin - m_FRH5DX1 / 2.0 / std::sin(7.5 * degree), 0., 0.);
			G4ThreeVector dr(2. * mm, 0., 0.);
			dr.rotateZ((7.5+15.*i) * degree);
			m_FRH5Pos[i] = G4ThreeVector(0., 0., m_FRH5PZ) + dr;
			m_FRH5Rot[i] = new G4RotationMatrix();
			m_FRH5Rot[i] -> rotateZ(- 15. * i * degree);

			// Place it
			m_FRH5PVName[i]    = "FRH5PV_"    + G4UIcommand::ConvertToString(i+1);
			m_FRH5MylPVName[i] = "FRH5MylPV_" + G4UIcommand::ConvertToString(i+1);
			m_FRH5TedPVName[i] = "FRH5TedPV_" + G4UIcommand::ConvertToString(i+1);
			new G4PVPlacement(m_FRH5Rot[i], m_FRH5Pos[i], m_FRH5PVName[i]   , m_FRH5LV[i]   , m_LabPV, false, 0);
			new G4PVPlacement(m_FRH5Rot[i], m_FRH5Pos[i], m_FRH5MylPVName[i], m_FRH5MylLV[i], m_LabPV, false, 0);
			new G4PVPlacement(m_FRH5Rot[i], m_FRH5Pos[i], m_FRH5TedPVName[i], m_FRH5TedLV[i], m_LabPV, false, 0);
		}
	}

	// FPC solids
	m_WireSolid      = new G4Tubs("WireSolid"     , 0., m_WireDia/2            ,  m_TubeL/2 + .2*mm                      , 0., 360.*degree);
	m_WireSSolid     = new G4Tubs("WireSSolid"    , 0., m_WireDia/2            , (m_TubeL-m_PipeDia-2.*m_VacThi-0.2*mm)/4, 0., 360.*degree);
	m_TubeBaseSolid  = new G4Tubs("TubeBaseSolid" , 0., m_TubeDia/2            ,  m_TubeL/2 + .1*mm                      , 0., 360.*degree);
	m_TubeSBaseSolid = new G4Tubs("TubeSBaseSolid", 0., m_TubeDia/2            , (m_TubeL-m_PipeDia-2.*m_VacThi-0.2*mm)/4+.1*mm, 0., 360.*degree);
	m_WallBaseSolid  = new G4Tubs("WallBaseSolid" , 0., m_TubeDia/2 + m_WallThi,  m_TubeL/2                              , 0., 360.*degree);
	m_WallSBaseSolid = new G4Tubs("WallSBaseSolid", 0., m_TubeDia/2 + m_WallThi, (m_TubeL-m_PipeDia-2.*m_VacThi-0.2*mm)/4, 0., 360.*degree);
	m_ShieldSubSolid  = new G4Box("ShieldBaseSolid", 62. * (m_TubeDia+0.1*mm),
	                                                 m_TubeL / 2 + 0.1 * mm,
	                                                 1.5 * std::sqrt(3.) * (m_TubeDia+0.1*mm)                        );
	m_ShieldBaseSolid = new G4Box("ShieldBaseSolid", 62. * (m_TubeDia+0.1*mm) + m_FPCAlThi,
	                                                 m_TubeL / 2,
	                                                 1.5 * std::sqrt(3.) * (m_TubeDia+0.1*mm) + m_FPCAlThi + 0.1 * mm);
	m_TubeSolid   = new G4SubtractionSolid("TubeSolid"   , m_TubeBaseSolid   , m_WireSolid      , 0, G4ThreeVector());
	m_TubeSSolid  = new G4SubtractionSolid("TubeSSolid"  , m_TubeSBaseSolid  , m_WireSSolid     , 0, G4ThreeVector());
	m_WallSolid   = new G4SubtractionSolid("WallSolid"   , m_WallBaseSolid   , m_TubeBaseSolid  , 0, G4ThreeVector());
	m_WallSSolid  = new G4SubtractionSolid("WallSSolid"  , m_WallSBaseSolid  , m_TubeSBaseSolid , 0, G4ThreeVector());
	m_ShieldSolid = new G4SubtractionSolid("ShieldSolid" , m_ShieldBaseSolid , m_ShieldSubSolid , 0, G4ThreeVector());

	// FPC12
	for ( int i = 0; i < 260; i++ )
	{
		int iEl = i - 130;
		if ( iEl >= 0 ) iEl += 1;

		// Logical volume
		m_FPC1WireLVName[i] = "FPC1WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC1TubeLVName[i] = "FPC1TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC1WallLVName[i] = "FPC1WallLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2WireLVName[i] = "FPC2WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2TubeLVName[i] = "FPC2TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2WallLVName[i] = "FPC2WallLV_" + G4UIcommand::ConvertToString(iEl);
		if     ( iEl <= - 121 )
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC12PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC12PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i+1)));
		}
		else if ( iEl <= -  11 )
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC12PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC12PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
		}
		else if ( iEl <= -   1 )
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC12PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC12PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
		}
		else if ( iEl <=    10 )
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC12PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC12PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
		}
		else if ( iEl <=   120 )
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC12PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC12PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
		}
		else
		{
			m_FPC1WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC1WireLVName[i].data());
			m_FPC1TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC1TubeLVName[i].data());
			m_FPC1WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC1TubeLVName[i].data());
			m_FPC2WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC2WireLVName[i].data());
			m_FPC2TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC2TubeLVName[i].data());
			m_FPC2WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC2TubeLVName[i].data());

			m_FPC1Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC12PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i+1)));
			m_FPC2Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC12PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i)));
		}
		m_FPC1WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC1TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC1WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC2WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC2TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC2WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC1Pos[i].rotateZ(225. * degree);
		m_FPC2Pos[i].rotateZ(225. * degree);

		// Place it
		m_FPC12Rot = new G4RotationMatrix();
		m_FPC12Rot -> rotateX( 90. * degree);
		m_FPC12Rot -> rotateY(225. * degree);

		m_FPC1WirePVName[i] = "FPC1WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC1TubePVName[i] = "FPC1TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC1WallPVName[i] = "FPC1WallPV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2WirePVName[i] = "FPC2WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2TubePVName[i] = "FPC2TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC2WallPVName[i] = "FPC2WallPV_" + G4UIcommand::ConvertToString(iEl);
		new G4PVPlacement(m_FPC12Rot, m_FPC1Pos[i], m_FPC1WirePVName[i], m_FPC1WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC12Rot, m_FPC1Pos[i], m_FPC1TubePVName[i], m_FPC1TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC12Rot, m_FPC1Pos[i], m_FPC1WallPVName[i], m_FPC1WallLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC12Rot, m_FPC2Pos[i], m_FPC2WirePVName[i], m_FPC2WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC12Rot, m_FPC2Pos[i], m_FPC2TubePVName[i], m_FPC2TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC12Rot, m_FPC2Pos[i], m_FPC2WallPVName[i], m_FPC2WallLV[i], m_LabPV, false, 0);
	}

	// FPC34
	for ( int i = 0; i < 260; i++ )
	{
		int iEl = i - 130;
		if ( iEl >= 0 ) iEl += 1;

		// Logical volume
		m_FPC3WireLVName[i] = "FPC3WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC3TubeLVName[i] = "FPC3TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC3WallLVName[i] = "FPC3WallLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4WireLVName[i] = "FPC4WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4TubeLVName[i] = "FPC4TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4WallLVName[i] = "FPC4WallLV_" + G4UIcommand::ConvertToString(iEl);
		if     ( iEl <= - 121 )
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC34PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC34PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i+1)));
		}
		else if ( iEl <= -  11 )
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC34PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC34PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
		}
		else if ( iEl <= -   1 )
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC34PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC34PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
		}
		else if ( iEl <=    10 )
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC34PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC34PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
		}
		else if ( iEl <=   120 )
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC34PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC34PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
		}
		else
		{
			m_FPC3WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC3WireLVName[i].data());
			m_FPC3TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC3TubeLVName[i].data());
			m_FPC3WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC3TubeLVName[i].data());
			m_FPC4WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC4WireLVName[i].data());
			m_FPC4TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC4TubeLVName[i].data());
			m_FPC4WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC4TubeLVName[i].data());

			m_FPC3Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC34PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i+1)));
			m_FPC4Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC34PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i)));
		}
		m_FPC3WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC3TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC3WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC4WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC4TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC4WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC3Pos[i].rotateZ(315. * degree);
		m_FPC4Pos[i].rotateZ(315. * degree);

		// Place it
		m_FPC34Rot = new G4RotationMatrix();
		m_FPC34Rot -> rotateX( 90. * degree);
		m_FPC34Rot -> rotateY(315. * degree);

		m_FPC3WirePVName[i] = "FPC3WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC3TubePVName[i] = "FPC3TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC3WallPVName[i] = "FPC3WallPV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4WirePVName[i] = "FPC4WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4TubePVName[i] = "FPC4TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC4WallPVName[i] = "FPC4WallPV_" + G4UIcommand::ConvertToString(iEl);
		new G4PVPlacement(m_FPC34Rot, m_FPC3Pos[i], m_FPC3WirePVName[i], m_FPC3WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC34Rot, m_FPC3Pos[i], m_FPC3TubePVName[i], m_FPC3TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC34Rot, m_FPC3Pos[i], m_FPC3WallPVName[i], m_FPC3WallLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC34Rot, m_FPC4Pos[i], m_FPC4WirePVName[i], m_FPC4WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC34Rot, m_FPC4Pos[i], m_FPC4TubePVName[i], m_FPC4TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC34Rot, m_FPC4Pos[i], m_FPC4WallPVName[i], m_FPC4WallLV[i], m_LabPV, false, 0);
	}

	// FPC56
	for ( int i = 0; i < 260; i++ )
	{
		int iEl = i - 130;
		if ( iEl >= 0 ) iEl += 1;

		// Logical volume
		m_FPC5WireLVName[i] = "FPC5WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC5TubeLVName[i] = "FPC5TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC5WallLVName[i] = "FPC5WallLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6WireLVName[i] = "FPC6WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6TubeLVName[i] = "FPC6TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6WallLVName[i] = "FPC6WallLV_" + G4UIcommand::ConvertToString(iEl);
		if     ( iEl <= - 121 )
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC56PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC56PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i+1)));
		}
		else if ( iEl <= -  11 )
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC56PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC56PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
		}
		else if ( iEl <= -   1 )
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC56PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC56PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
		}
		else if ( iEl <=    10 )
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC56PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC56PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
		}
		else if ( iEl <=   120 )
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC56PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC56PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
		}
		else
		{
			m_FPC5WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC5WireLVName[i].data());
			m_FPC5TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC5TubeLVName[i].data());
			m_FPC5WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC5TubeLVName[i].data());
			m_FPC6WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC6WireLVName[i].data());
			m_FPC6TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC6TubeLVName[i].data());
			m_FPC6WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC6TubeLVName[i].data());

			m_FPC5Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC56PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i+1)));
			m_FPC6Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC56PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i)));
		}
		m_FPC5WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC5TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC5WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC6WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC6TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC6WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC5Pos[i].rotateZ(270. * degree);
		m_FPC6Pos[i].rotateZ(270. * degree);

		// Place it
		m_FPC56Rot = new G4RotationMatrix();
		m_FPC56Rot -> rotateX( 90. * degree);
		m_FPC56Rot -> rotateY(270. * degree);

		m_FPC5WirePVName[i] = "FPC5WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC5TubePVName[i] = "FPC5TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC5WallPVName[i] = "FPC5WallPV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6WirePVName[i] = "FPC6WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6TubePVName[i] = "FPC6TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC6WallPVName[i] = "FPC6WallPV_" + G4UIcommand::ConvertToString(iEl);
		new G4PVPlacement(m_FPC56Rot, m_FPC5Pos[i], m_FPC5WirePVName[i], m_FPC5WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC56Rot, m_FPC5Pos[i], m_FPC5TubePVName[i], m_FPC5TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC56Rot, m_FPC5Pos[i], m_FPC5WallPVName[i], m_FPC5WallLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC56Rot, m_FPC6Pos[i], m_FPC6WirePVName[i], m_FPC6WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC56Rot, m_FPC6Pos[i], m_FPC6TubePVName[i], m_FPC6TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC56Rot, m_FPC6Pos[i], m_FPC6WallPVName[i], m_FPC6WallLV[i], m_LabPV, false, 0);
	}

	// FPC78
	for ( int i = 0; i < 260; i++ )
	{
		int iEl = i - 130;
		if ( iEl >= 0 ) iEl += 1;

		// Logical volume
		m_FPC7WireLVName[i] = "FPC7WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC7TubeLVName[i] = "FPC7TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC7WallLVName[i] = "FPC7WallLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8WireLVName[i] = "FPC8WireLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8TubeLVName[i] = "FPC8TubeLV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8WallLVName[i] = "FPC8WallLV_" + G4UIcommand::ConvertToString(iEl);
		if     ( iEl <= - 121 )
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC78PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl+120)/2,-(m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4,m_FPC78PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2-std::pow(-1,i+1)));
		}
		else if ( iEl <= -  11 )
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC78PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC78PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
		}
		else if ( iEl <= -   1 )
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC78PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC78PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
		}
		else if ( iEl <=    10 )
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC78PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i+1)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*iEl/2,  (m_TubeL+m_PipeDia+2.*m_VacThi+.2*mm)/4, m_FPC78PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1, i)));
		}
		else if ( iEl <=   120 )
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSolid , m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSolid , m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSolid , m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC78PZ + (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i+1)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm) * iEl / 2, 0., m_FPC78PZ - (m_TubeDia+.1*mm) * std::sqrt(3.) / 4 * (2.-std::pow(-1, i)));
		}
		else
		{
			m_FPC7WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC7WireLVName[i].data());
			m_FPC7TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC7TubeLVName[i].data());
			m_FPC7WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC7TubeLVName[i].data());
			m_FPC8WireLV[i] = new G4LogicalVolume(m_WireSSolid, m_SSMat , m_FPC8WireLVName[i].data());
			m_FPC8TubeLV[i] = new G4LogicalVolume(m_TubeSSolid, m_GasMat, m_FPC8TubeLVName[i].data());
			m_FPC8WallLV[i] = new G4LogicalVolume(m_WallSSolid, m_MylMat, m_FPC8TubeLVName[i].data());

			m_FPC7Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC78PZ+(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i+1)));
			m_FPC8Pos[i] = G4ThreeVector((m_TubeDia+.1*mm)*(iEl-120)/2,-(m_TubeL+m_PipeDia+2*m_VacThi+.2*mm)/4,m_FPC78PZ-(m_TubeDia+.1*mm)*std::sqrt(3)/4*(2.-std::pow(-1,i)));
		}
		m_FPC7WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC7TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC7WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
		m_FPC8WireLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC8TubeLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC8WallLV[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
//		m_FPC7Pos[i].rotateZ(  0. * degree);
//		m_FPC8Pos[i].rotateZ(  0. * degree);

		// Place it
		m_FPC78Rot = new G4RotationMatrix();
		m_FPC78Rot -> rotateX( 90. * degree);
//		m_FPC78Rot -> rotateY(  0. * degree);

		m_FPC7WirePVName[i] = "FPC7WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC7TubePVName[i] = "FPC7TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC7WallPVName[i] = "FPC7WallPV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8WirePVName[i] = "FPC8WirePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8TubePVName[i] = "FPC8TubePV_" + G4UIcommand::ConvertToString(iEl);
		m_FPC8WallPVName[i] = "FPC8WallPV_" + G4UIcommand::ConvertToString(iEl);
		new G4PVPlacement(m_FPC78Rot, m_FPC7Pos[i], m_FPC7WirePVName[i], m_FPC7WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC78Rot, m_FPC7Pos[i], m_FPC7TubePVName[i], m_FPC7TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC78Rot, m_FPC7Pos[i], m_FPC7WallPVName[i], m_FPC7WallLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC78Rot, m_FPC8Pos[i], m_FPC8WirePVName[i], m_FPC8WireLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC78Rot, m_FPC8Pos[i], m_FPC8TubePVName[i], m_FPC8TubeLV[i], m_LabPV, false, 0);
		new G4PVPlacement(m_FPC78Rot, m_FPC8Pos[i], m_FPC8WallPVName[i], m_FPC8WallLV[i], m_LabPV, false, 0);
	}

	// FPC shield
	m_FPC12ShieldLVName = "FPC12ShieldLV";
	m_FPC34ShieldLVName = "FPC34ShieldLV";
	m_FPC56ShieldLVName = "FPC56ShieldLV";
	m_FPC78ShieldLVName = "FPC78ShieldLV";
	m_FPC12ShieldLV = new G4LogicalVolume(m_ShieldSolid, m_AluMat, m_FPC12ShieldLVName.data());
	m_FPC34ShieldLV = new G4LogicalVolume(m_ShieldSolid, m_AluMat, m_FPC34ShieldLVName.data());
	m_FPC56ShieldLV = new G4LogicalVolume(m_ShieldSolid, m_AluMat, m_FPC56ShieldLVName.data());
	m_FPC78ShieldLV = new G4LogicalVolume(m_ShieldSolid, m_AluMat, m_FPC78ShieldLVName.data());
	m_FPC12ShieldRot = new G4RotationMatrix();
	m_FPC34ShieldRot = new G4RotationMatrix();
	m_FPC56ShieldRot = new G4RotationMatrix();
	m_FPC78ShieldRot = new G4RotationMatrix();
	m_FPC12ShieldRot -> rotateZ(135. * degree);
	m_FPC34ShieldRot -> rotateZ(225. * degree);
	m_FPC56ShieldRot -> rotateZ( 90. * degree);
	m_FPC78ShieldRot -> rotateZ(270. * degree);
	m_FPC12ShieldPVName = "FPC12ShieldPV";
	m_FPC34ShieldPVName = "FPC34ShieldPV";
	m_FPC56ShieldPVName = "FPC56ShieldPV";
	m_FPC78ShieldPVName = "FPC78ShieldPV";
	new G4PVPlacement(m_FPC12ShieldRot, G4ThreeVector(0., 0., m_FPC12PZ), m_FPC12ShieldPVName, m_FPC12ShieldLV, m_LabPV, false, 0);
	new G4PVPlacement(m_FPC34ShieldRot, G4ThreeVector(0., 0., m_FPC34PZ), m_FPC34ShieldPVName, m_FPC34ShieldLV, m_LabPV, false, 0);
	new G4PVPlacement(m_FPC56ShieldRot, G4ThreeVector(0., 0., m_FPC56PZ), m_FPC56ShieldPVName, m_FPC56ShieldLV, m_LabPV, false, 0);
	new G4PVPlacement(m_FPC78ShieldRot, G4ThreeVector(0., 0., m_FPC78PZ), m_FPC78ShieldPVName, m_FPC78ShieldLV, m_LabPV, false, 0);

	// FPC makes graphics too slow...
	if ( !m_FPCVis )
	{
		for ( int i = 0; i < 260; i++ )
		{
			m_FPC1WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC2WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC3WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC4WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC5WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC6WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC7WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
			m_FPC8WallLV[i] -> SetVisAttributes(new G4VisAttributes(0));
		}
		m_FPC12ShieldLV -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC34ShieldLV -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC56ShieldLV -> SetVisAttributes(new G4VisAttributes(0));
		m_FPC78ShieldLV -> SetVisAttributes(new G4VisAttributes(0));
	}

	return m_LabPV;
}

void WMC4DetCon::ConstructMaterials()
{
	const G4double labTemp = 300.0 * kelvin;

	// Elements to be used to construct materials
	m_ElH  = new G4Element( "Hydrogen",  "H",  1,   1.00794*g/mole);
	m_ElC  = new G4Element(   "Carbon",  "C",  6,  12.011  *g/mole);
	m_ElN  = new G4Element( "Nitrogen",  "N",  7,  14.00674*g/mole);
	m_ElO  = new G4Element(   "Oxygen",  "O",  8,  15.9994 *g/mole);
	m_ElF  = new G4Element( "Fluorine",  "F",  9,  18.9984 *g/mole);
	m_ElMg = new G4Element("Magnesium", "Mg", 12,  24.3050 *g/mole);
	m_ElAl = new G4Element( "Aluminum", "Al", 13,  26.98154*g/mole);
	m_ElSi = new G4Element(  "Silicon", "Si", 14,  28.0855 *g/mole);
	m_ElAr = new G4Element(    "Argon", "Ar", 18,  39.948  *g/mole);
	m_ElCr = new G4Element("Chromuium", "Cr", 24,  51.9961 *g/mole);
	m_ElMn = new G4Element("Manganese", "Mn", 25,  54.93804*g/mole);
	m_ElFe = new G4Element(     "Iron", "Fe", 26,  55.845  *g/mole);
	m_ElNi = new G4Element(   "Nickel", "Ni", 28,  58.9634 *g/mole);
	m_ElCu = new G4Element(   "Copper", "Cu", 29,  63.546  *g/mole);
	m_ElZn = new G4Element(     "Zinc", "Zn", 30,  65.38   *g/mole);
	m_ElY  = new G4Element(  "Yttrium",  "Y", 39,  88.90585*g/mole);
	m_ElSn = new G4Element(      "Tin", "Sn", 50, 118.71   *g/mole);
	m_ElCe = new G4Element(   "Cerium", "Ce", 58, 140.116  *g/mole);
	m_ElLu = new G4Element( "Lutetium", "Lu", 71, 176.944  *g/mole);

	// Vacuum material
	m_VacMat = new G4Material("Vacuum", 0.1205e-5*g/cm3, 3, kStateGas, labTemp);
	m_VacMat -> AddElement(m_ElN, 75.47/99.95);
	m_VacMat -> AddElement(m_ElO, 23.20/99.95);
	m_VacMat -> AddElement(m_ElAr, 1.28/99.95);

	// Air material
	m_AirMat = new G4Material("Air", 1.1839e-3*g/cm3, 3, kStateGas, labTemp); // Density at 1 atm
	m_AirMat -> AddElement(m_ElN, 75.47/99.95);
	m_AirMat -> AddElement(m_ElO, 23.20/99.95);
	m_AirMat -> AddElement(m_ElAr, 1.28/99.95);

	// Target material
	G4String targetMaterial = m_CM -> GetTarMat();
	if ( targetMaterial == "C" )
	{
		m_TarMat = new G4Material("Graphite", 1.7*g/cm3, 3, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElC, 99.0/100.0);
		m_TarMat -> AddElement(m_ElN,  0.7/100.0);
		m_TarMat -> AddElement(m_ElO,  0.3/100.0);
	}
	else if ( targetMaterial == "Mg" )
	{
		m_TarMat = new G4Material("Magnesium", 1.738*g/cm3, 1, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElMg, 1.0);
	}
	else if ( targetMaterial == "Al" )
	{
		m_TarMat = new G4Material("Aluminum", 2.70*g/cm3, 1, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElAl, 1.0);
	}
	else if ( targetMaterial == "Sn" )
	{
		m_TarMat = new G4Material("Tin", 5.769*g/cm3, 1, kStateSolid, labTemp); // Gray tin?
		m_TarMat -> AddElement(m_ElSn, 1.0);
	}
	else if ( targetMaterial == "Si" )
	{
		m_TarMat = new G4Material("Silicon", 2.329*g/cm3, 1, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElSi, 1.0);
	}
	else if ( targetMaterial == "Ni" )
	{
		m_TarMat = new G4Material("Nickel", 8.908*g/cm3, 1, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElNi, 1.0);
	}
	else if ( targetMaterial == "Zn" )
	{
		m_TarMat = new G4Material("Zinc", 7.14*g/cm3, 1, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElZn, 1.0);
	}
	else if ( targetMaterial == "CH2" )
	{
		m_TarMat = new G4Material("Polyethelene", 0.92*g/cm3, 2, kStateSolid, labTemp);
		m_TarMat -> AddElement(m_ElC, 1);
		m_TarMat -> AddElement(m_ElH, 2);
	}
	else
	{
		std::stringstream message;
		message << "Could not find target material with name " << targetMaterial;
		G4Exception("WMC4DetCon::ConstructMaterials()", "", FatalException, message.str().c_str());
	}

	// Plastic scintillator material
	m_PlaMat = new G4Material("PlasticScint", 1.032*g/cm3, 2, kStateSolid, labTemp);
	m_PlaMat -> AddElement(m_ElC, 10);
	m_PlaMat -> AddElement(m_ElH, 11);

	// Mylar (PET) material
	m_MylMat = new G4Material("Mylar", 1.38*g/cm3, 3, kStateSolid, labTemp);
	m_MylMat -> AddElement(m_ElC, 5);
	m_MylMat -> AddElement(m_ElH, 4);
	m_MylMat -> AddElement(m_ElO, 2);

	// Tedlar material
	m_TedMat = new G4Material("Tedlar", 1.545*g/cm3, 3, kStateSolid, labTemp);
	m_TedMat -> AddElement(m_ElC, 2);
	m_TedMat -> AddElement(m_ElH, 3);
	m_TedMat -> AddElement(m_ElF, 1);

	// ALuminum material
	m_AluMat = new G4Material("Aluminum", 2.70*g/cm3, 1, kStateSolid, labTemp);
	m_AluMat -> AddElement(m_ElAl, 1.0);

	// Stainless steel material
	m_SSMat = new G4Material("StainlessSteel", 8.06*g/cm3, 6, kStateSolid, labTemp);
	m_SSMat -> AddElement(m_ElC , 0.001);
	m_SSMat -> AddElement(m_ElSi, 0.007);
	m_SSMat -> AddElement(m_ElCr, 0.18 );
	m_SSMat -> AddElement(m_ElMn, 0.01 );
	m_SSMat -> AddElement(m_ElFe, 0.712);
	m_SSMat -> AddElement(m_ElNi, 0.09 );

	// Gas meterial
//	G4int gasSpecies = m_CM -> GetFPCGas();
	m_GasMat = new G4Material("Gas", 1.69844e-3*g/cm3, 3, kStateGas, labTemp);
	m_GasMat -> AddElement(m_ElC , 1);
	m_GasMat -> AddElement(m_ElAr, 2);
	m_GasMat -> AddElement(m_ElH , 3);
}

void WMC4DetCon::DestructMaterials()
{
	delete m_GasMat;
	delete m_SSMat;
	delete m_AluMat;
	delete m_TedMat;
	delete m_MylMat;
	delete m_PlaMat;
	delete m_TarMat;
	delete m_AirMat;
	delete m_VacMat;

	delete m_ElLu;
	delete m_ElCe;
	delete m_ElSn;
	delete m_ElY;
	delete m_ElZn;
	delete m_ElCu;
	delete m_ElNi;
	delete m_ElFe;
	delete m_ElMn;
	delete m_ElCr;
	delete m_ElAr;
	delete m_ElSi;
	delete m_ElAl;
	delete m_ElMg;
	delete m_ElF;
	delete m_ElO;
	delete m_ElN;
	delete m_ElC;
	delete m_ElH;
}
