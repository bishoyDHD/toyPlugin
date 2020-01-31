#include <algorithm>

#include "../../../include/track.h"

track::track(){
  fNHits=0;
  fChi2PerNDF=0;
  fFlagTargetToAC=0;
  fFlagToFMatch=kTRUE;

  fMwpc4PathLocation=0;
  fMwpc4PathTime=0;
  fMwpc4EnergyLoss=0;
  fMwpc4MX=0;
  fMwpc4MY=0;
  fMwpc4MZ=0;
  fMwpc4MXErr=0;
  fMwpc4MYErr=0;
  fMwpc4MZErr=0;
  fMwpc4SX=0;
  fMwpc4SY=0;
  fMwpc4SZ=0;
  fMwpc4SNx=0;
  fMwpc4SNy=0;
  fMwpc4SNz=0;
  fMwpc4SP=0;

  fMwpc3PathLocation=0;
  fMwpc3PathTime=0;
  fMwpc3EnergyLoss=0;
  fMwpc3MX=0;
  fMwpc3MY=0;
  fMwpc3MZ=0;
  fMwpc3MXErr=0;
  fMwpc3MYErr=0;
  fMwpc3MZErr=0;
  fMwpc3SX=0;
  fMwpc3SY=0;
  fMwpc3SZ=0;
  fMwpc3SNx=0;
  fMwpc3SNy=0;
  fMwpc3SNz=0;
  fMwpc3SP=0;

  fMwpc2PathLocation=0;
  fMwpc2PathTime=0;
  fMwpc2EnergyLoss=0;
  fMwpc2MX=0;
  fMwpc2MY=0;
  fMwpc2MZ=0;
  fMwpc2MXErr=0;
  fMwpc2MYErr=0;
  fMwpc2MZErr=0;
  fMwpc2SX=0;
  fMwpc2SY=0;
  fMwpc2SZ=0;
  fMwpc2SNx=0;
  fMwpc2SNy=0;
  fMwpc2SNz=0;
  fMwpc2SP=0;

  fTof2GapNum=0;
  fTof2PathLocation=0;
  fTof2PathTime=0;
  fTof2Part='Z';
  fTof2TimeI=0;
  fTof2TimeO=0;
  fTof2TimeMean=0;
  fTof2SX=0;
  fTof2SY=0;
  fTof2SZ=0;
  fTof2SNx=0;
  fTof2SNy=0;
  fTof2SNz=0;
  fTof2SP=0;

  fBLowEdgeSX=0;
  fBLowEdgeSY=0;
  fBLowEdgeSZ=0;
  fBLowEdgeSNx=0;
  fBLowEdgeSNy=0;
  fBLowEdgeSNz=0;
  fBLowEdgeSP=0;

  fTof1GapNum=0;
  fTof1PathLocation=0;
  fTof1PathTime=0;
  fTof1TimeU=0;
  fTof1TimeD=0;
  fTof1TimeMean=0;
  fTof1SX=0;
  fTof1SY=0;
  fTof1SZ=0;
  fTof1SNx=0;
  fTof1SNy=0;
  fTof1SNz=0;
  fTof1SP=0;

  fSftSX=0;
  fSftSY=0;
  fSftSZ=0;
  fSftSNx=0;
  fSftSNy=0;
  fSftSNz=0;
  fSftSP=0;

  fVertPathLocation=0;
  fVertPathTime=0;
  fVertMX=0;
  fVertMY=0;
  fVertMPhi=0;
  fVertMXErr=0;
  fVertMYErr=0;
  fVertMPhiErr=0;
  fVertSX=0;
  fVertSY=0;
  fVertSZ=0;
  fVertSPhi=0;

  fTof1PathTimePionPlus=0;
  fTof1SPPionPlus=0;
  fSftSPPionPlus=0;
  fVertPathTimePionPlus=0;
  fVertSPPionPlus=0;

  fTof1PathTimePositron=0;
  fTof1SPPositron=0;
  fSftSPPositron=0;
  fVertPathTimePositron=0;
  fVertSPPositron=0;
}

ClassImp(track)

trackArray::trackArray(){}

trackArray::~trackArray(){
  vectTrack.clear();
}

ClassImp(trackArray)


