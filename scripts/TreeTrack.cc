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


#include "TreeTrack.hh"

ClassImp(TreeTrack);


//------------------------------------------------------------------
TreeTrack::TreeTrack(){
}
//------------------------------------------------------------------                   
TreeTrack::~TreeTrack(){
}
//------------------------------------------------------------------
void TreeTrack::Clear(Option_t *option){
  for(int i=0; i<gNlayers; i++){
    fTdc[i]=-1;
    fAdc[i]=-1;
    fEl[i]=-1;

    fTheta = -1000;
    fPhi = -1000;
  }
}
//------------------------------------------------------------------
void TreeTrack::SetADCs(Double_t* adc){
  fAdc[0] = adc[0]; //fwc1
  fAdc[1] = adc[1]; //fwc2
  fAdc[2] = adc[2];  //fth
  fAdc[3] = adc[3];  //trh1
  fAdc[4] = adc[4];  //frh2
  fAdc[5] = adc[5];  //frh3
  fAdc[6] = adc[6];  //frh4
  fAdc[7] = adc[7];  //frh5
}  
//------------------------------------------------------------------
void TreeTrack::SetTDCs(Double_t* tdc){

  fTdc[0] = tdc[0]; //fwc1
  fTdc[1] = tdc[1]; //fwc2
  fTdc[2] = tdc[2];  //fth
  fTdc[3] = tdc[3];  //trh1
  fTdc[4] = tdc[4];  //frh2
  fTdc[5] = tdc[5];  //frh3
  fTdc[6] = tdc[6];  //frh4
  fTdc[7] = tdc[7];  //frh5
}
//------------------------------------------------------------------
void TreeTrack::SetElements(Int_t* el){

  fEl[0] = el[0]; //fwc1
  fEl[1] = el[1]; //fwc2
  fEl[2] = el[2];  //fth
  fEl[3] = el[3];  //trh1
  fEl[4] = el[4];  //frh2
  fEl[5] = el[5];  //frh3
  fEl[6] = el[6];  //frh4
  fEl[7] = el[7];  //frh5
}
//------------------------------------------------------------------
//void TreeTrack::SetPIDs(Bool_t isP, Bool_t isD){
//fIsP = isP;
//fIsD = isD;
//}
