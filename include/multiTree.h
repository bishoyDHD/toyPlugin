/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __MULTITREE_H_
#define __MULTITREE_H_

#include "cookerrawtree.h" // for CRTBase
#include <vector>
static const double E_kpi2=0.10854566040017; // in GeV
class CRTtrackVar:public CRTBase{
 public:
  UInt_t eventNo;
  CRTtrackVar();
  virtual ~CRTtrackVar();
  ClassDef(CRTtrackVar,1);
};
class targetE36: public CRTBase {
public:
  UInt_t run;
  UInt_t eventNumber;
  UInt_t eventFlag;
  Int_t TOF1Gap;
  Int_t TOF2Gap;
  /*
  Int_t extraTOF1;
  Float_t phiAngle;
  Float_t deltaPhiAngle;
  Float_t chiS;
  Int_t ndf;
  Float_t reducedChiS;
  Int_t leptonCounter;
  Float_t intersectTargetX;
  Float_t intersectTargetY;
  Float_t deltaX;
  Float_t deltaY;
  Float_t intersectSFTX;
  Float_t intersectSFTY;
  Float_t kaonCentroidX;
  Float_t kaonCentroidY;
  Float_t kaonCentroidXErr;
  Float_t kaonCentroidYErr;
  Float_t kaonStopIntersectX;
  Float_t kaonStopIntersectY;
  Float_t kaonStopX;
  Float_t kaonStopY;
  Float_t kaonStopXEr;
  Float_t kaonStopYEr;
  Float_t RKstop;
  Int_t kaonBarStop;
  Int_t kaonClusterSize;
  Float_t kaonReducedChiS;
  Int_t ckMultiplicity;
  Int_t cpiMultiplicity;
  Float_t trackLength;
  Float_t C2XCentroid;
  Float_t TDCDiff;
  Int_t sumADCHGLeptons;
  Float_t leptonAverageTDC;
  Float_t kaonAverageTDC;
  Float_t m_kaon;
  Float_t Dm_kaon;
  Float_t b_kaon;
  Float_t Db_kaon;
  Float_t rho_kaon;
  Float_t C_kaon;
  Float_t m_lepton;
  Float_t Dm_lepton;
  Float_t b_lepton;
  Float_t Db_lepton;
  Float_t rho_lepton;
  Float_t C_lepton;
  Bool_t hasEdgeBars;
  Bool_t isGoodTOF1;
  Int_t pruningMethod;
  Int_t kStopType;
  Int_t caseNum;
  Int_t kstopErrFlag;
  Int_t badEventFlag;*/
  targetE36();
  virtual ~targetE36();
  ClassDef(targetE36,2);
};
class mwpcE36: public CRTBase{
public:
  UInt_t run ;
  UInt_t event;
  Int_t nTracks;
  Int_t fgapNumTof2;
  Int_t nHits;
  Double_t fTof2SP;
  Double_t fVertSP;
  mwpcE36();
  virtual ~mwpcE36();
  ClassDef(mwpcE36,2);
};
class trackingE36 : public CRTBase{
public:
  Double_t nxVert, nyVert, nzVert;
  Double_t xVert, yVert, zVert;
  Double_t pVertpi0;
  Int_t tof2Gap, evtNum;
  trackingE36();
  virtual ~trackingE36();
  ClassDef(trackingE36, 2);
};

#endif

