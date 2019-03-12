#ifndef WMC4DETCON_h
#define WMC4DETCON_h 1

////////////////////////////////////////////////////////////////////////////////
//   WMC4DetCon.hh                                                            //
//                                                                            //
//   This file is a header for WMC4DetCon class. It's for construction of     //
// whole geometry of simulation, which includes target geometry and detector  //
// geometry. Details of configuration can be set at config.cfg file. The      //
// configuration is tranfered at main function and used here to generate      //
// world of simulation.                                                       //
//                                                                            //
//                    - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"

#include "WMC4ConMan.hh"

class G4VPhysicalVolume;

class WMC4ConMan;

class WMC4DetCon: public G4VUserDetectorConstruction
{
  public:
	WMC4DetCon(WMC4ConMan* CM);
	virtual ~WMC4DetCon();
	virtual G4VPhysicalVolume* Construct();

  private:
	void DefineDimensions();
	void ConstructMaterials();
	void DestructMaterials();

  private:
	// Configuration manager
	WMC4ConMan* m_CM;

	// Elements
	G4Element* m_ElH;
	G4Element* m_ElC;
	G4Element* m_ElN;
	G4Element* m_ElO;
	G4Element* m_ElF;
	G4Element* m_ElMg;
	G4Element* m_ElAl;
	G4Element* m_ElSi;
	G4Element* m_ElAr;
	G4Element* m_ElCr;
	G4Element* m_ElMn;
	G4Element* m_ElFe;
	G4Element* m_ElNi;
	G4Element* m_ElCu;
	G4Element* m_ElZn;
	G4Element* m_ElY;
	G4Element* m_ElSn;
	G4Element* m_ElCe;
	G4Element* m_ElLu;

	// Materials
	G4Material* m_VacMat;
	G4Material* m_AirMat;
	G4Material* m_TarMat;
	G4Material* m_PlaMat;
	G4Material* m_MylMat;
	G4Material* m_TedMat;
	G4Material* m_AluMat;
	G4Material* m_SSMat;
	G4Material* m_GasMat;

	// Dimensions and detector setup
	G4double m_LabX, m_LabY, m_LabZ;
	G4double m_PipeDia, m_VacThi;
	G4bool   m_UseTar;
	G4double m_TarD, m_TarT;
	G4bool   m_UseFWC1, m_UseFWC2, m_UseFTH, m_UseFRH1, m_UseFRH2, m_UseFRH3, m_UseFRH4, m_UseFRH5;
	std::vector<G4bool> m_FWC1Mask, m_FWC2Mask, m_FTHMask, m_FRH1Mask, m_FRH2Mask, m_FRH3Mask, m_FRH4Mask, m_FRH5Mask;
	G4double m_FWC1DX1, m_FWC1DY,  m_FWC1Thi, m_FWC1RMin, m_FWC1RMax, m_FWC1PZ;
	G4double m_FWC2DX1, m_FWC2DY,  m_FWC2Thi, m_FWC2RMin, m_FWC2RMax, m_FWC2PZ;
	G4double m_FTHDX1 , m_FTHDY ,  m_FTHThi , m_FTHRMin , m_FTHRMax , m_FTHPZ;
	G4double m_FRH1DX1, m_FRH1DY,  m_FRH1Thi, m_FRH1RMin, m_FRH1RMax, m_FRH1PZ;
	G4double m_FRH2DX1, m_FRH2DY,  m_FRH2Thi, m_FRH2RMin, m_FRH2RMax, m_FRH2PZ;
	G4double m_FRH3DX1, m_FRH3DY,  m_FRH3Thi, m_FRH3RMin, m_FRH3RMax, m_FRH3PZ;
	G4double m_FRH4DX1, m_FRH4DY,  m_FRH4Thi, m_FRH4RMin, m_FRH4RMax, m_FRH4PZ;
	G4double m_FRH5DX1, m_FRH5DY,  m_FRH5Thi, m_FRH5RMin, m_FRH5RMax, m_FRH5PZ;
	G4double m_MylThi , m_TedThi;
	G4double m_FPC12Phi, m_FPC34Phi, m_FPC56Phi, m_FPC78Phi;
	G4double m_FPC12PZ , m_FPC34PZ , m_FPC56PZ , m_FPC78PZ;
	G4double m_FPCAlThi, m_TubeDia , m_WallThi , m_WireDia, m_TubeL;
	G4int    m_FPCGas;
	G4bool   m_FPCVis;

	// Geometry objects: Lab
	G4Box* m_LabSolid;
	G4LogicalVolume* m_LabLV;
	G4VPhysicalVolume* m_LabPV;

	// Pipe & chamber
	G4Cons* m_PCFSolid;
	G4Cons* m_PCFVacSolid;
	G4Tubs* m_VCF1Solid;
	G4Tubs* m_VCF2Solid;
	G4Tubs* m_VCF3Solid;
	G4Tubs* m_BTFSolid;
	G4Tubs* m_BTCDSolid;
	G4Tubs* m_BTCDVacSolid;
	G4Tubs* m_BTBSolid;
	G4Tubs* m_BTBVacSolid;
	G4Tubs* m_BBBSolid;
	G4Tubs* m_BBBVacSolid;
	G4Box* m_VCBSolid;
	G4Sphere* m_VCWSolid; 
	G4LogicalVolume* m_PCFLV;
	G4LogicalVolume* m_PCFVacLV;
	G4LogicalVolume* m_VCF1LV;
	G4LogicalVolume* m_VCF2LV;
	G4LogicalVolume* m_VCF3LV;
	G4LogicalVolume* m_VCB1LV;	
	G4LogicalVolume* m_VCB2LV;
	G4LogicalVolume* m_VCB3LV;
	G4LogicalVolume* m_VCB4LV;
	G4LogicalVolume* m_VCWLV;
	G4LogicalVolume* m_BTFLV;
	G4LogicalVolume* m_BTCDLV;
	G4LogicalVolume* m_BTCDVacLV;
	G4LogicalVolume* m_BTBLV;
	G4LogicalVolume* m_BTBVacLV;
	G4LogicalVolume* m_BBBLV;
	G4LogicalVolume* m_BBBVacLV;
	G4RotationMatrix* m_VCBRot;
	G4RotationMatrix* m_VCWRot;

	// Geometry objects: Target
	G4Tubs* m_TarSolid;
	G4LogicalVolume* m_TarLV;

	// Geometry objects: Wrapping assistant
	G4RotationMatrix* m_MylAttachRot;
	G4RotationMatrix* m_MylAttachRotN;
	G4RotationMatrix* m_MylAttachRot2;
	G4RotationMatrix* m_TedAttachRot;
	G4RotationMatrix* m_TedAttachRotN;
	G4RotationMatrix* m_TedAttachRot2;

	// Geometry objects: FWC1
	G4Tubs* m_FWC1Solid;
	G4Tubs* m_FWC1NSolid;
	G4Tubs* m_FWC1MylBaseSolid;
	G4Tubs* m_FWC1MylBaseNSolid;
	G4Tubs* m_FWC1TedBaseSolid;
	G4Tubs* m_FWC1TedBaseNSolid;
	G4Box* m_FWC1MylAttachSolid;
	G4Box* m_FWC1TedAttachSolid;
	G4SubtractionSolid* m_FWC1MylSolidTemp1;
	G4SubtractionSolid* m_FWC1MylNSolidTemp1;
	G4SubtractionSolid* m_FWC1TedSolidTemp1;
	G4SubtractionSolid* m_FWC1TedNSolidTemp1;
	G4UnionSolid* m_FWC1MylSolidTemp2;
	G4UnionSolid* m_FWC1MylNSolidTemp2;
	G4UnionSolid* m_FWC1TedSolidTemp2;
	G4UnionSolid* m_FWC1TedNSolidTemp2;
	G4UnionSolid* m_FWC1MylSolid;
	G4UnionSolid* m_FWC1MylNSolid;
	G4UnionSolid* m_FWC1TedSolid;
	G4UnionSolid* m_FWC1TedNSolid;
	G4String m_FWC1LVName[24];
	G4String m_FWC1MylLVName[24];
	G4String m_FWC1TedLVName[24];
	G4String m_FWC1PVName[24];
	G4String m_FWC1MylPVName[24];
	G4String m_FWC1TedPVName[24];
	G4LogicalVolume* m_FWC1LV[24];
	G4LogicalVolume* m_FWC1MylLV[24];
	G4LogicalVolume* m_FWC1TedLV[24];
	G4ThreeVector m_FWC1Pos[24];
	G4RotationMatrix* m_FWC1Rot[24];
	G4int m_SmallFWC1El[24] = {0};

	// Geometry objects: FWC2
	G4Tubs* m_FWC2Solid;
	G4Tubs* m_FWC2MylBaseSolid;
	G4Tubs* m_FWC2TedBaseSolid;
	G4Box* m_FWC2MylAttachSolid;
	G4Box* m_FWC2TedAttachSolid;
	G4SubtractionSolid* m_FWC2MylSolidTemp1;
	G4SubtractionSolid* m_FWC2TedSolidTemp1;
	G4UnionSolid* m_FWC2MylSolidTemp2;
	G4UnionSolid* m_FWC2TedSolidTemp2;
	G4UnionSolid* m_FWC2MylSolid;
	G4UnionSolid* m_FWC2TedSolid;
	G4String m_FWC2LVName[24];
	G4String m_FWC2MylLVName[24];
	G4String m_FWC2TedLVName[24];
	G4String m_FWC2PVName[24];
	G4String m_FWC2MylPVName[24];
	G4String m_FWC2TedPVName[24];
	G4LogicalVolume* m_FWC2LV[24];
	G4LogicalVolume* m_FWC2MylLV[24];
	G4LogicalVolume* m_FWC2TedLV[24];
	G4ThreeVector m_FWC2Pos[24];
	G4RotationMatrix* m_FWC2Rot[24];

	// Geometry objects: FTH
	G4Tubs* m_FTHSolid;
	G4Tubs* m_FTHMylBaseSolid;
	G4Tubs* m_FTHTedBaseSolid;
	G4Box* m_FTHMylAttachSolid;
	G4Box* m_FTHTedAttachSolid;
	G4SubtractionSolid* m_FTHMylSolidTemp1;
	G4SubtractionSolid* m_FTHTedSolidTemp1;
	G4UnionSolid* m_FTHMylSolidTemp2;
	G4UnionSolid* m_FTHTedSolidTemp2;
	G4UnionSolid* m_FTHMylSolid;
	G4UnionSolid* m_FTHTedSolid;
	G4String m_FTHLVName[48];
	G4String m_FTHMylLVName[48];
	G4String m_FTHTedLVName[48];
	G4String m_FTHPVName[48];
	G4String m_FTHMylPVName[48];
	G4String m_FTHTedPVName[48];
	G4LogicalVolume* m_FTHLV[48];
	G4LogicalVolume* m_FTHMylLV[48];
	G4LogicalVolume* m_FTHTedLV[48];
	G4ThreeVector m_FTHPos[48];
	G4RotationMatrix* m_FTHRot[48];

	// Geometry objects: FRH1
	G4Tubs* m_FRH1Solid;
	G4Tubs* m_FRH1MylBaseSolid;
	G4Tubs* m_FRH1TedBaseSolid;
	G4Box* m_FRH1MylAttachSolid;
	G4Box* m_FRH1TedAttachSolid;
	G4SubtractionSolid* m_FRH1MylSolidTemp1;
	G4SubtractionSolid* m_FRH1TedSolidTemp1;
	G4UnionSolid* m_FRH1MylSolidTemp2;
	G4UnionSolid* m_FRH1TedSolidTemp2;
	G4UnionSolid* m_FRH1MylSolid;
	G4UnionSolid* m_FRH1TedSolid;
	G4String m_FRH1LVName[24];
	G4String m_FRH1MylLVName[24];
	G4String m_FRH1TedLVName[24];
	G4String m_FRH1PVName[24];
	G4String m_FRH1MylPVName[24];
	G4String m_FRH1TedPVName[24];
	G4LogicalVolume* m_FRH1LV[24];
	G4LogicalVolume* m_FRH1MylLV[24];
	G4LogicalVolume* m_FRH1TedLV[24];
	G4ThreeVector m_FRH1Pos[24];
	G4RotationMatrix* m_FRH1Rot[24];

	// Geometry objects: FRH2
	G4Tubs* m_FRH2Solid;
	G4Tubs* m_FRH2MylBaseSolid;
	G4Tubs* m_FRH2TedBaseSolid;
	G4Box* m_FRH2MylAttachSolid;
	G4Box* m_FRH2TedAttachSolid;
	G4SubtractionSolid* m_FRH2MylSolidTemp1;
	G4SubtractionSolid* m_FRH2TedSolidTemp1;
	G4UnionSolid* m_FRH2MylSolidTemp2;
	G4UnionSolid* m_FRH2TedSolidTemp2;
	G4UnionSolid* m_FRH2MylSolid;
	G4UnionSolid* m_FRH2TedSolid;
	G4String m_FRH2LVName[24];
	G4String m_FRH2MylLVName[24];
	G4String m_FRH2TedLVName[24];
	G4String m_FRH2PVName[24];
	G4String m_FRH2MylPVName[24];
	G4String m_FRH2TedPVName[24];
	G4LogicalVolume* m_FRH2LV[24];
	G4LogicalVolume* m_FRH2MylLV[24];
	G4LogicalVolume* m_FRH2TedLV[24];
	G4ThreeVector m_FRH2Pos[24];
	G4RotationMatrix* m_FRH2Rot[24];

	// Geometry objects: FRH3
	G4Tubs* m_FRH3Solid;
	G4Tubs* m_FRH3MylBaseSolid;
	G4Tubs* m_FRH3TedBaseSolid;
	G4Box* m_FRH3MylAttachSolid;
	G4Box* m_FRH3TedAttachSolid;
	G4SubtractionSolid* m_FRH3MylSolidTemp1;
	G4SubtractionSolid* m_FRH3TedSolidTemp1;
	G4UnionSolid* m_FRH3MylSolidTemp2;
	G4UnionSolid* m_FRH3TedSolidTemp2;
	G4UnionSolid* m_FRH3MylSolid;
	G4UnionSolid* m_FRH3TedSolid;
	G4String m_FRH3LVName[24];
	G4String m_FRH3MylLVName[24];
	G4String m_FRH3TedLVName[24];
	G4String m_FRH3PVName[24];
	G4String m_FRH3MylPVName[24];
	G4String m_FRH3TedPVName[24];
	G4LogicalVolume* m_FRH3LV[24];
	G4LogicalVolume* m_FRH3MylLV[24];
	G4LogicalVolume* m_FRH3TedLV[24];
	G4ThreeVector m_FRH3Pos[24];
	G4RotationMatrix* m_FRH3Rot[24];

	// Geometry objects: FRH4
	G4Tubs* m_FRH4Solid;
	G4Tubs* m_FRH4MylBaseSolid;
	G4Tubs* m_FRH4TedBaseSolid;
	G4Box* m_FRH4MylAttachSolid;
	G4Box* m_FRH4TedAttachSolid;
	G4SubtractionSolid* m_FRH4MylSolidTemp1;
	G4SubtractionSolid* m_FRH4TedSolidTemp1;
	G4UnionSolid* m_FRH4MylSolidTemp2;
	G4UnionSolid* m_FRH4TedSolidTemp2;
	G4UnionSolid* m_FRH4MylSolid;
	G4UnionSolid* m_FRH4TedSolid;
	G4String m_FRH4LVName[24];
	G4String m_FRH4MylLVName[24];
	G4String m_FRH4TedLVName[24];
	G4String m_FRH4PVName[24];
	G4String m_FRH4MylPVName[24];
	G4String m_FRH4TedPVName[24];
	G4LogicalVolume* m_FRH4LV[24];
	G4LogicalVolume* m_FRH4MylLV[24];
	G4LogicalVolume* m_FRH4TedLV[24];
	G4ThreeVector m_FRH4Pos[24];
	G4RotationMatrix* m_FRH4Rot[24];

	// Geometry objects: FRH5
	G4Tubs* m_FRH5Solid;
	G4Tubs* m_FRH5MylBaseSolid;
	G4Tubs* m_FRH5TedBaseSolid;
	G4Box* m_FRH5MylAttachSolid;
	G4Box* m_FRH5TedAttachSolid;
	G4SubtractionSolid* m_FRH5MylSolidTemp1;
	G4SubtractionSolid* m_FRH5TedSolidTemp1;
	G4UnionSolid* m_FRH5MylSolidTemp2;
	G4UnionSolid* m_FRH5TedSolidTemp2;
	G4UnionSolid* m_FRH5MylSolid;
	G4UnionSolid* m_FRH5TedSolid;
	G4String m_FRH5LVName[24];
	G4String m_FRH5MylLVName[24];
	G4String m_FRH5TedLVName[24];
	G4String m_FRH5PVName[24];
	G4String m_FRH5MylPVName[24];
	G4String m_FRH5TedPVName[24];
	G4LogicalVolume* m_FRH5LV[24];
	G4LogicalVolume* m_FRH5MylLV[24];
	G4LogicalVolume* m_FRH5TedLV[24];
	G4ThreeVector m_FRH5Pos[24];
	G4RotationMatrix* m_FRH5Rot[24];

	// Geometry objects: FPC general
	G4Tubs* m_WireSolid;
	G4Tubs* m_WireSSolid;
	G4Tubs* m_TubeBaseSolid;
	G4Tubs* m_TubeSBaseSolid;
	G4Tubs* m_WallBaseSolid;
	G4Tubs* m_WallSBaseSolid;
	G4Box*  m_ShieldBaseSolid;
	G4Box*  m_ShieldSubSolid;
	G4SubtractionSolid* m_TubeSolid;
	G4SubtractionSolid* m_TubeSSolid;
	G4SubtractionSolid* m_WallSolid;
	G4SubtractionSolid* m_WallSSolid;
	G4SubtractionSolid* m_ShieldSolid;

	// Geometry objects: FPC12
	G4String m_FPC1WireLVName[260];
	G4String m_FPC1TubeLVName[260];
	G4String m_FPC1WallLVName[260];
	G4String m_FPC2WireLVName[260];
	G4String m_FPC2TubeLVName[260];
	G4String m_FPC2WallLVName[260];
	G4String m_FPC1WirePVName[260];
	G4String m_FPC1TubePVName[260];
	G4String m_FPC1WallPVName[260];
	G4String m_FPC2WirePVName[260];
	G4String m_FPC2TubePVName[260];
	G4String m_FPC2WallPVName[260];
	G4LogicalVolume* m_FPC1WireLV[260];
	G4LogicalVolume* m_FPC1TubeLV[260];
	G4LogicalVolume* m_FPC1WallLV[260];
	G4LogicalVolume* m_FPC2WireLV[260];
	G4LogicalVolume* m_FPC2TubeLV[260];
	G4LogicalVolume* m_FPC2WallLV[260];
	G4ThreeVector m_FPC1Pos[260];
	G4ThreeVector m_FPC2Pos[260];
	G4ThreeVector m_FPC12Pos;
	G4RotationMatrix* m_FPC12Rot;

	// Geometry objects: FPC34
	G4String m_FPC3WireLVName[260];
	G4String m_FPC3TubeLVName[260];
	G4String m_FPC3WallLVName[260];
	G4String m_FPC4WireLVName[260];
	G4String m_FPC4TubeLVName[260];
	G4String m_FPC4WallLVName[260];
	G4String m_FPC3WirePVName[260];
	G4String m_FPC3TubePVName[260];
	G4String m_FPC3WallPVName[260];
	G4String m_FPC4WirePVName[260];
	G4String m_FPC4TubePVName[260];
	G4String m_FPC4WallPVName[260];
	G4LogicalVolume* m_FPC3WireLV[260];
	G4LogicalVolume* m_FPC3TubeLV[260];
	G4LogicalVolume* m_FPC3WallLV[260];
	G4LogicalVolume* m_FPC4WireLV[260];
	G4LogicalVolume* m_FPC4TubeLV[260];
	G4LogicalVolume* m_FPC4WallLV[260];
	G4ThreeVector m_FPC3Pos[260];
	G4ThreeVector m_FPC4Pos[260];
	G4RotationMatrix* m_FPC34Rot;

	// Geometry objects: FPC56
	G4String m_FPC5WireLVName[260];
	G4String m_FPC5TubeLVName[260];
	G4String m_FPC5WallLVName[260];
	G4String m_FPC6WireLVName[260];
	G4String m_FPC6TubeLVName[260];
	G4String m_FPC6WallLVName[260];
	G4String m_FPC5WirePVName[260];
	G4String m_FPC5TubePVName[260];
	G4String m_FPC5WallPVName[260];
	G4String m_FPC6WirePVName[260];
	G4String m_FPC6TubePVName[260];
	G4String m_FPC6WallPVName[260];
	G4LogicalVolume* m_FPC5WireLV[260];
	G4LogicalVolume* m_FPC5TubeLV[260];
	G4LogicalVolume* m_FPC5WallLV[260];
	G4LogicalVolume* m_FPC6WireLV[260];
	G4LogicalVolume* m_FPC6TubeLV[260];
	G4LogicalVolume* m_FPC6WallLV[260];
	G4ThreeVector m_FPC5Pos[260];
	G4ThreeVector m_FPC6Pos[260];
	G4RotationMatrix* m_FPC56Rot;

	// Geometry objects: FPC78
	G4String m_FPC7WireLVName[260];
	G4String m_FPC7TubeLVName[260];
	G4String m_FPC7WallLVName[260];
	G4String m_FPC8WireLVName[260];
	G4String m_FPC8TubeLVName[260];
	G4String m_FPC8WallLVName[260];
	G4String m_FPC7WirePVName[260];
	G4String m_FPC7TubePVName[260];
	G4String m_FPC7WallPVName[260];
	G4String m_FPC8WirePVName[260];
	G4String m_FPC8TubePVName[260];
	G4String m_FPC8WallPVName[260];
	G4LogicalVolume* m_FPC7WireLV[260];
	G4LogicalVolume* m_FPC7TubeLV[260];
	G4LogicalVolume* m_FPC7WallLV[260];
	G4LogicalVolume* m_FPC8WireLV[260];
	G4LogicalVolume* m_FPC8TubeLV[260];
	G4LogicalVolume* m_FPC8WallLV[260];
	G4ThreeVector m_FPC7Pos[260];
	G4ThreeVector m_FPC8Pos[260];
	G4RotationMatrix* m_FPC78Rot;

	// Geometry objects: FPC shield
	G4String m_FPC12ShieldLVName;
	G4String m_FPC34ShieldLVName;
	G4String m_FPC56ShieldLVName;
	G4String m_FPC78ShieldLVName;
	G4String m_FPC12ShieldPVName;
	G4String m_FPC34ShieldPVName;
	G4String m_FPC56ShieldPVName;
	G4String m_FPC78ShieldPVName;
	G4LogicalVolume* m_FPC12ShieldLV;
	G4LogicalVolume* m_FPC34ShieldLV;
	G4LogicalVolume* m_FPC56ShieldLV;
	G4LogicalVolume* m_FPC78ShieldLV;
	G4RotationMatrix* m_FPC12ShieldRot;
	G4RotationMatrix* m_FPC34ShieldRot;
	G4RotationMatrix* m_FPC56ShieldRot;
	G4RotationMatrix* m_FPC78ShieldRot;
};
#endif
