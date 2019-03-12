/******************************************************************************
 *
 * Environment:
 *    Software development for WASA detector system at COSY
 * (c) 2016 The JEDI Collaboration
 *
 * Author List:
 *    Aleksandra Wronska     JU Krakow
 *
 *****************************************************************************/


#ifndef __TreeTrack_H_
#define __TreeTrack_H_ 1
#include "TObject.h"

const int gNlayers = 8; 

class TreeTrack : public TObject{
  
public:
  TreeTrack();
  ~TreeTrack();
  void Clear(Option_t *option = "");
  void SetPhi(Double_t phi){fPhi=phi;};  
  void SetTheta(Double_t theta){fTheta=theta;};  
  void SetPhiTrue(Double_t phi){fPhiTrue=phi;};  
  void SetThetaTrue(Double_t theta){fThetaTrue=theta;};  
  void SetETrue(Double_t etrue){fETrue=etrue;};  
  void SetADCs(Double_t* adc);  
  void SetTDCs(Double_t* tdc);  
  void SetElements(Int_t* el);  
  //  void SetPIDs(Bool_t isP, Bool_t isD);

  Double_t fTdc[gNlayers];
  Double_t fAdc[gNlayers];
  Int_t fEl[gNlayers];
  Double_t fTheta;
  Double_t fPhi;
  Double_t fThetaTrue;
  Double_t fPhiTrue;
  Double_t fETrue;
  //  Bool_t fIsP, fIsD;

private:  

  ClassDef(TreeTrack,1)
  
};

#endif
