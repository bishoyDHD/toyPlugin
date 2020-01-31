#ifndef __TRACKTREE_H_
#define __TRACKTREE_H_
#include "cookerrawtree.h" // for CRTBase

#include <vector>
#include <iostream>
#include "TClonesArray.h"
using namespace std;

#include "TObject.h"


class track : public CRTBase{
 public:
  track();
  ~track(){}


  inline void setNHits(Int_t n) { fNHits = n; }
  inline Int_t  getNHits() const { return fNHits; } 

  inline void setChi2PerNDF(Double_t chi2PerNDF) { fChi2PerNDF = chi2PerNDF; }
  inline Double_t  getChi2PerNDF() const { return fChi2PerNDF; } 

  inline void setFlagTargetToAC(Int_t flagTargetToAC) { fFlagTargetToAC = flagTargetToAC; }
  inline Int_t  getFlagTargetToAC() const { return fFlagTargetToAC; } 

  inline void setFlagToFMatch(Int_t flagToFMatch) { fFlagToFMatch = flagToFMatch; }
  inline Bool_t  getFlagToFMatch() const { return fFlagToFMatch; } 

  inline void setMwpc4PathLocation(Double_t s) { fMwpc4PathLocation = s; }
  inline Double_t  getMwpc4PathLocation() const { return fMwpc4PathLocation; }
  inline void setMwpc4PathTime(Double_t t) { fMwpc4PathTime = t; }
  inline Double_t  getMwpc4PathTime() const { return fMwpc4PathTime; }
  inline void setMwpc4EnergyLoss(Double_t energyLoss) { fMwpc4EnergyLoss = energyLoss; }
  inline Double_t  getMwpc4EnergyLoss() const { return fMwpc4EnergyLoss; }

  inline void setMwpc4MX(Double_t x) { fMwpc4MX = x; }
  inline Double_t  getMwpc4MX() const { return fMwpc4MX; }
  inline void setMwpc4MY(Double_t y) { fMwpc4MY = y; }
  inline Double_t  getMwpc4MY() const { return fMwpc4MY; }
  inline void setMwpc4MZ(Double_t z) { fMwpc4MZ = z; }
  inline Double_t  getMwpc4MZ() const { return fMwpc4MZ; }

  inline void setMwpc4MXErr(Double_t xErr) { fMwpc4MXErr = xErr; }
  inline Double_t  getMwpc4MXErr() const { return fMwpc4MXErr; }
  inline void setMwpc4MYErr(Double_t yErr) { fMwpc4MYErr = yErr; }
  inline Double_t  getMwpc4MYErr() const { return fMwpc4MYErr; }
  inline void setMwpc4MZErr(Double_t zErr) { fMwpc4MZErr = zErr; }
  inline Double_t  getMwpc4MZErr() const { return fMwpc4MZErr; }

  inline void setMwpc4SX(Double_t x) { fMwpc4SX = x; }
  inline Double_t  getMwpc4SX() const { return fMwpc4SX; }
  inline void setMwpc4SY(Double_t y) { fMwpc4SY = y; }
  inline Double_t  getMwpc4SY() const { return fMwpc4SY; }
  inline void setMwpc4SZ(Double_t z) { fMwpc4SZ = z; }
  inline Double_t  getMwpc4SZ() const { return fMwpc4SZ; }
  inline void setMwpc4SNx(Double_t nx) { fMwpc4SNx = nx; }
  inline Double_t  getMwpc4SNx() const { return fMwpc4SNx; }
  inline void setMwpc4SNy(Double_t ny) { fMwpc4SNy = ny; }
  inline Double_t  getMwpc4SNy() const { return fMwpc4SNy; }
  inline void setMwpc4SNz(Double_t nz) { fMwpc4SNz = nz; }
  inline Double_t  getMwpc4SNz() const { return fMwpc4SNz; }
  inline void setMwpc4SP(Double_t p) { fMwpc4SP = p; }
  inline Double_t  getMwpc4SP() const { return fMwpc4SP; }

  inline void setMwpc3PathLocation(Double_t s) { fMwpc3PathLocation = s; }
  inline Double_t  getMwpc3PathLocation() const { return fMwpc3PathLocation; }
  inline void setMwpc3PathTime(Double_t t) { fMwpc3PathTime = t; }
  inline Double_t  getMwpc3PathTime() const { return fMwpc3PathTime; }
  inline void setMwpc3EnergyLoss(Double_t energyLoss) { fMwpc3EnergyLoss = energyLoss; }
  inline Double_t  getMwpc3EnergyLoss() const { return fMwpc3EnergyLoss; }

  inline void setMwpc3MX(Double_t x) { fMwpc3MX = x; }
  inline Double_t  getMwpc3MX() const { return fMwpc3MX; }
  inline void setMwpc3MY(Double_t y) { fMwpc3MY = y; }
  inline Double_t  getMwpc3MY() const { return fMwpc3MY; }
  inline void setMwpc3MZ(Double_t z) { fMwpc3MZ = z; }
  inline Double_t  getMwpc3MZ() const { return fMwpc3MZ; }

  inline void setMwpc3MXErr(Double_t xErr) { fMwpc3MXErr = xErr; }
  inline Double_t  getMwpc3MXErr() const { return fMwpc3MXErr; }
  inline void setMwpc3MYErr(Double_t yErr) { fMwpc3MYErr = yErr; }
  inline Double_t  getMwpc3MYErr() const { return fMwpc3MYErr; }
  inline void setMwpc3MZErr(Double_t zErr) { fMwpc3MZErr = zErr; }
  inline Double_t  getMwpc3MZErr() const { return fMwpc3MZErr; }

  inline void setMwpc3SX(Double_t x) { fMwpc3SX = x; }
  inline Double_t  getMwpc3SX() const { return fMwpc3SX; }
  inline void setMwpc3SY(Double_t y) { fMwpc3SY = y; }
  inline Double_t  getMwpc3SY() const { return fMwpc3SY; }
  inline void setMwpc3SZ(Double_t z) { fMwpc3SZ = z; }
  inline Double_t  getMwpc3SZ() const { return fMwpc3SZ; }
  inline void setMwpc3SNx(Double_t nx) { fMwpc3SNx = nx; }
  inline Double_t  getMwpc3SNx() const { return fMwpc3SNx; }
  inline void setMwpc3SNy(Double_t ny) { fMwpc3SNy = ny; }
  inline Double_t  getMwpc3SNy() const { return fMwpc3SNy; }
  inline void setMwpc3SNz(Double_t nz) { fMwpc3SNz = nz; }
  inline Double_t  getMwpc3SNz() const { return fMwpc3SNz; }
  inline void setMwpc3SP(Double_t p) { fMwpc3SP = p; }
  inline Double_t  getMwpc3SP() const { return fMwpc3SP; }

  inline void setMwpc2PathLocation(Double_t s) { fMwpc2PathLocation = s; }
  inline Double_t  getMwpc2PathLocation() const { return fMwpc2PathLocation; }
  inline void setMwpc2PathTime(Double_t t) { fMwpc2PathTime = t; }
  inline Double_t  getMwpc2PathTime() const { return fMwpc2PathTime; }
  inline void setMwpc2EnergyLoss(Double_t energyLoss) { fMwpc2EnergyLoss = energyLoss; }
  inline Double_t  getMwpc2EnergyLoss() const { return fMwpc2EnergyLoss; }

  inline void setMwpc2MX(Double_t x) { fMwpc2MX = x; }
  inline Double_t  getMwpc2MX() const { return fMwpc2MX; }
  inline void setMwpc2MY(Double_t y) { fMwpc2MY = y; }
  inline Double_t  getMwpc2MY() const { return fMwpc2MY; }
  inline void setMwpc2MZ(Double_t z) { fMwpc2MZ = z; }
  inline Double_t  getMwpc2MZ() const { return fMwpc2MZ; }

  inline void setMwpc2MXErr(Double_t xErr) { fMwpc2MXErr = xErr; }
  inline Double_t  getMwpc2MXErr() const { return fMwpc2MXErr; }
  inline void setMwpc2MYErr(Double_t yErr) { fMwpc2MYErr = yErr; }
  inline Double_t  getMwpc2MYErr() const { return fMwpc2MYErr; }
  inline void setMwpc2MZErr(Double_t zErr) { fMwpc2MZErr = zErr; }
  inline Double_t  getMwpc2MZErr() const { return fMwpc2MZErr; }

  inline void setMwpc2SX(Double_t x) { fMwpc2SX = x; }
  inline Double_t  getMwpc2SX() const { return fMwpc2SX; }
  inline void setMwpc2SY(Double_t y) { fMwpc2SY = y; }
  inline Double_t  getMwpc2SY() const { return fMwpc2SY; }
  inline void setMwpc2SZ(Double_t z) { fMwpc2SZ = z; }
  inline Double_t  getMwpc2SZ() const { return fMwpc2SZ; }
  inline void setMwpc2SNx(Double_t nx) { fMwpc2SNx = nx; }
  inline Double_t  getMwpc2SNx() const { return fMwpc2SNx; }
  inline void setMwpc2SNy(Double_t ny) { fMwpc2SNy = ny; }
  inline Double_t  getMwpc2SNy() const { return fMwpc2SNy; }
  inline void setMwpc2SNz(Double_t nz) { fMwpc2SNz = nz; }
  inline Double_t  getMwpc2SNz() const { return fMwpc2SNz; }
  inline void setMwpc2SP(Double_t p) { fMwpc2SP = p; }
  inline Double_t  getMwpc2SP() const { return fMwpc2SP; }

  inline void setTof2GapNum(Int_t gapNum) { fTof2GapNum = gapNum; }
  inline Int_t  getTof2GapNum() const { return fTof2GapNum; }
  inline void setTof2PathLocation(Double_t s) { fTof2PathLocation = s; }
  inline Double_t  getTof2PathLocation() const { return fTof2PathLocation; }
  inline void setTof2PathTime(Double_t t) { fTof2PathTime = t; }
  inline Double_t  getTof2PathTime() const { return fTof2PathTime; }

  inline void setTof2Part(Char_t part) { fTof2Part = part; }
  inline Char_t  getTof2Part() const { return fTof2Part; }
  inline void setTof2TimeI(Double_t timeI) { fTof2TimeI = timeI; }
  inline Double_t  getTof2TimeI() const { return fTof2TimeI; }
  inline void setTof2TimeO(Double_t timeO) { fTof2TimeO = timeO; }
  inline Double_t  getTof2TimeO() const { return fTof2TimeO; }
  inline void setTof2TimeMean(Double_t timeMean) { fTof2TimeMean = timeMean; }
  inline Double_t  getTof2TimeMean() const { return fTof2TimeMean; }

  inline void setTof2SX(Double_t x) { fTof2SX = x; }
  inline Double_t  getTof2SX() const { return fTof2SX; }
  inline void setTof2SY(Double_t y) { fTof2SY = y; }
  inline Double_t  getTof2SY() const { return fTof2SY; }
  inline void setTof2SZ(Double_t z) { fTof2SZ = z; }
  inline Double_t  getTof2SZ() const { return fTof2SZ; }
  inline void setTof2SNx(Double_t nx) { fTof2SNx = nx; }
  inline Double_t  getTof2SNx() const { return fTof2SNx; }
  inline void setTof2SNy(Double_t ny) { fTof2SNy = ny; }
  inline Double_t  getTof2SNy() const { return fTof2SNy; }
  inline void setTof2SNz(Double_t nz) { fTof2SNz = nz; }
  inline Double_t  getTof2SNz() const { return fTof2SNz; }
  inline void setTof2SP(Double_t p) { fTof2SP = p; }
  inline Double_t  getTof2SP() const { return fTof2SP; }

  inline void setTof1GapNum(Int_t gapNum) { fTof1GapNum = gapNum; }
  inline Int_t  getTof1GapNum() const { return fTof1GapNum; }
  inline void setTof1PathLocation(Double_t s) { fTof1PathLocation = s; }
  inline Double_t  getTof1PathLocation() const { return fTof1PathLocation; }
  inline void setTof1PathTime(Double_t t) { fTof1PathTime = t; }
  inline Double_t  getTof1PathTime() const { return fTof1PathTime; }

  inline void setTof1TimeU(Double_t timeU) { fTof1TimeU = timeU; }
  inline Double_t  getTof1TimeU() const { return fTof1TimeU; }
  inline void setTof1TimeD(Double_t timeD) { fTof1TimeD = timeD; }
  inline Double_t  getTof1TimeD() const { return fTof1TimeD; }
  inline void setTof1TimeMean(Double_t timeMean) { fTof1TimeMean = timeMean; }
  inline Double_t  getTof1TimeMean() const { return fTof1TimeMean; }

  inline void setTof1SX(Double_t x) { fTof1SX = x; }
  inline Double_t  getTof1SX() const { return fTof1SX; }
  inline void setTof1SY(Double_t y) { fTof1SY = y; }
  inline Double_t  getTof1SY() const { return fTof1SY; }
  inline void setTof1SZ(Double_t z) { fTof1SZ = z; }
  inline Double_t  getTof1SZ() const { return fTof1SZ; }
  inline void setTof1SNx(Double_t nx) { fTof1SNx = nx; }
  inline Double_t  getTof1SNx() const { return fTof1SNx; }
  inline void setTof1SNy(Double_t ny) { fTof1SNy = ny; }
  inline Double_t  getTof1SNy() const { return fTof1SNy; }
  inline void setTof1SNz(Double_t nz) { fTof1SNz = nz; }
  inline Double_t  getTof1SNz() const { return fTof1SNz; }
  inline void setTof1SP(Double_t p) { fTof1SP = p; }
  inline Double_t  getTof1SP() const { return fTof1SP; }

  inline void setSftSX(Double_t x) { fSftSX = x; }
  inline Double_t  getSftSX() const { return fSftSX; }
  inline void setSftSY(Double_t y) { fSftSY = y; }
  inline Double_t  getSftSY() const { return fSftSY; }
  inline void setSftSZ(Double_t z) { fSftSZ = z; }
  inline Double_t  getSftSZ() const { return fSftSZ; }
  inline void setSftSNx(Double_t nx) { fSftSNx = nx; }
  inline Double_t  getSftSNx() const { return fSftSNx; }
  inline void setSftSNy(Double_t ny) { fSftSNy = ny; }
  inline Double_t  getSftSNy() const { return fSftSNy; }
  inline void setSftSNz(Double_t nz) { fSftSNz = nz; }
  inline Double_t  getSftSNz() const { return fSftSNz; }
  inline void setSftSP(Double_t p) { fSftSP = p; }
  inline Double_t  getSftSP() const { return fSftSP; }


  inline void setVertPathLocation(Double_t s) { fVertPathLocation = s; }
  inline Double_t  getVertPathLocation() const { return fVertPathLocation; }
  inline void setVertPathTime(Double_t t) { fVertPathTime = t; }
  inline Double_t  getVertPathTime() const { return fVertPathTime;}


  inline void setVertMX(Double_t x) { fVertMX = x; }
  inline Double_t  getVertMX() const { return fVertMX; }
  inline void setVertMY(Double_t y) { fVertMY = y; }
  inline Double_t  getVertMY() const { return fVertMY; }
  inline void setVertMPhi(Double_t phi) { fVertMPhi = phi; }
  inline Double_t  getVertMPhi() const { return fVertMPhi; }
  inline void setVertMXErr(Double_t xErr) { fVertMXErr = xErr; }
  inline Double_t  getVertMXErr() const { return fVertMXErr; }
  inline void setVertMYErr(Double_t yErr) { fVertMYErr = yErr; }
  inline Double_t  getVertMYErr() const { return fVertMYErr; }
  inline void setVertMPhiErr(Double_t phiErr) { fVertMPhiErr = phiErr; }
  inline Double_t  getVertMPhiErr() const { return fVertMPhiErr; }

  inline void setVertSX(Double_t x) { fVertSX = x; }
  inline Double_t  getVertSX() const { return fVertSX; }
  inline void setVertSY(Double_t y) { fVertSY = y; }
  inline Double_t  getVertSY() const { return fVertSY; }
  inline void setVertSZ(Double_t z) { fVertSZ = z; }
  inline Double_t  getVertSZ() const { return fVertSZ; }
  inline void setVertSNx(Double_t nx) { fVertSNx = nx; }
  inline Double_t  getVertSNx() const { return fVertSNx; }
  inline void setVertSNy(Double_t ny) { fVertSNy = ny; }
  inline Double_t  getVertSNy() const { return fVertSNy; }
  inline void setVertSNz(Double_t nz) { fVertSNz = nz; }
  inline Double_t  getVertSNz() const { return fVertSNz; }
  inline void setVertSP(Double_t p) { fVertSP = p; }
  inline Double_t  getVertSP() const { return fVertSP; }
  inline void setVertSPhi(Double_t phi) { fVertSPhi = phi; }
  inline Double_t  getVertSPhi() const { return fVertSPhi; }

  inline void setBLowEdgeSX(Double_t x) { fBLowEdgeSX = x; }
  inline Double_t  getBLowEdgeSX() const { return fBLowEdgeSX; }
  inline void setBLowEdgeSY(Double_t y) { fBLowEdgeSY = y; }
  inline Double_t  getBLowEdgeSY() const { return fBLowEdgeSY; }
  inline void setBLowEdgeSZ(Double_t z) { fBLowEdgeSZ = z; }
  inline Double_t  getBLowEdgeSZ() const { return fBLowEdgeSZ; }
  inline void setBLowEdgeSNx(Double_t nx) { fBLowEdgeSNx = nx; }
  inline Double_t  getBLowEdgeSNx() const { return fBLowEdgeSNx; }
  inline void setBLowEdgeSNy(Double_t ny) { fBLowEdgeSNy = ny; }
  inline Double_t  getBLowEdgeSNy() const { return fBLowEdgeSNy; }
  inline void setBLowEdgeSNz(Double_t nz) { fBLowEdgeSNz = nz; }
  inline Double_t  getBLowEdgeSNz() const { return fBLowEdgeSNz; }
  inline void setBLowEdgeSP(Double_t p) { fBLowEdgeSP = p; }
  inline Double_t  getBLowEdgeSP() const { return fBLowEdgeSP; }

  inline void setTof1PathTimePionPlus(Double_t t) { fTof1PathTimePionPlus = t; }
  inline Double_t  getTof1PathTimePionPlus() const { return fTof1PathTimePionPlus; }
  inline void setTof1SPPionPlus(Double_t p) { fTof1SPPionPlus = p; }
  inline Double_t  getTof1SPPionPlus() const { return fTof1SPPionPlus; }
  inline void setSftSPPionPlus(Double_t p) { fSftSPPionPlus = p; }
  inline Double_t  getSftSPPionPlus() const { return fSftSPPionPlus; }
  inline void setVertPathTimePionPlus(Double_t t) { fVertPathTimePionPlus = t; }
  inline Double_t  getVertPathTimePionPlus() const { return fVertPathTimePionPlus;}
  inline void setVertSPPionPlus(Double_t p) { fVertSPPionPlus = p; }
  inline Double_t  getVertSPPionPlus() const { return fVertSPPionPlus; }

  inline void setTof1PathTimePositron(Double_t t) { fTof1PathTimePositron = t; }
  inline Double_t  getTof1PathTimePositron() const { return fTof1PathTimePositron; }
  inline void setTof1SPPositron(Double_t p) { fTof1SPPositron = p; }
  inline Double_t  getTof1SPPositron() const { return fTof1SPPositron; }
  inline void setSftSPPositron(Double_t p) { fSftSPPositron = p; }
  inline Double_t  getSftSPPositron() const { return fSftSPPositron; }
  inline void setVertPathTimePositron(Double_t t) { fVertPathTimePositron = t; }
  inline Double_t  getVertPathTimePositron() const { return fVertPathTimePositron;}
  inline void setVertSPPositron(Double_t p) { fVertSPPositron = p; }
  inline Double_t  getVertSPPositron() const { return fVertSPPositron; }
  
 
private:
  Int_t fNHits;
  Double_t fChi2PerNDF;
  Int_t fFlagTargetToAC;
  Bool_t fFlagToFMatch;

  Double_t fMwpc4PathLocation;
  Double_t fMwpc4PathTime;
  Double_t fMwpc4EnergyLoss;

  Double_t fMwpc4MX;
  Double_t fMwpc4MY;
  Double_t fMwpc4MZ;
  Double_t fMwpc4MXErr;
  Double_t fMwpc4MYErr;
  Double_t fMwpc4MZErr;

  Double_t fMwpc4SX;
  Double_t fMwpc4SY;
  Double_t fMwpc4SZ;
  Double_t fMwpc4SNx;
  Double_t fMwpc4SNy;
  Double_t fMwpc4SNz;
  Double_t fMwpc4SP;

  Double_t fMwpc3PathLocation;
  Double_t fMwpc3PathTime;
  Double_t fMwpc3EnergyLoss;

  Double_t fMwpc3MX;
  Double_t fMwpc3MY;
  Double_t fMwpc3MZ;
  Double_t fMwpc3MXErr;
  Double_t fMwpc3MYErr;
  Double_t fMwpc3MZErr;

  Double_t fMwpc3SX;
  Double_t fMwpc3SY;
  Double_t fMwpc3SZ;
  Double_t fMwpc3SNx;
  Double_t fMwpc3SNy;
  Double_t fMwpc3SNz;
  Double_t fMwpc3SP;

  Double_t fMwpc2PathLocation;
  Double_t fMwpc2PathTime;
  Double_t fMwpc2EnergyLoss;

  Double_t fMwpc2MX;
  Double_t fMwpc2MY;
  Double_t fMwpc2MZ;
  Double_t fMwpc2MXErr;
  Double_t fMwpc2MYErr;
  Double_t fMwpc2MZErr;

  Double_t fMwpc2SX;
  Double_t fMwpc2SY;
  Double_t fMwpc2SZ;
  Double_t fMwpc2SNx;
  Double_t fMwpc2SNy;
  Double_t fMwpc2SNz;
  Double_t fMwpc2SP;

  Int_t fTof2GapNum;
  Double_t fTof2PathLocation;
  Double_t fTof2PathTime;

  Char_t fTof2Part;
  Double_t fTof2TimeI;
  Double_t fTof2TimeO;
  Double_t fTof2TimeMean;

  Double_t fTof2SX;
  Double_t fTof2SY;
  Double_t fTof2SZ;
  Double_t fTof2SNx;
  Double_t fTof2SNy;
  Double_t fTof2SNz;
  Double_t fTof2SP;

  Double_t fBLowEdgeSX;
  Double_t fBLowEdgeSY;
  Double_t fBLowEdgeSZ;
  Double_t fBLowEdgeSNx;
  Double_t fBLowEdgeSNy;
  Double_t fBLowEdgeSNz;
  Double_t fBLowEdgeSP;

  Int_t fTof1GapNum;
  Double_t fTof1PathLocation;
  Double_t fTof1PathTime;

  Double_t fTof1TimeU;
  Double_t fTof1TimeD;
  Double_t fTof1TimeMean;

  Double_t fTof1SX;
  Double_t fTof1SY;
  Double_t fTof1SZ;
  Double_t fTof1SNx;
  Double_t fTof1SNy;
  Double_t fTof1SNz;
  Double_t fTof1SP;

  Double_t fSftSX;
  Double_t fSftSY;
  Double_t fSftSZ;
  Double_t fSftSNx;
  Double_t fSftSNy;
  Double_t fSftSNz;
  Double_t fSftSP;

  Double_t fVertPathLocation;
  Double_t fVertPathTime;

  Double_t fVertMX;
  Double_t fVertMY;
  Double_t fVertMPhi;
  Double_t fVertMXErr;
  Double_t fVertMYErr;
  Double_t fVertMPhiErr;

  Double_t fVertSX;
  Double_t fVertSY;
  Double_t fVertSZ;
  Double_t fVertSNx;
  Double_t fVertSNy;
  Double_t fVertSNz;
  Double_t fVertSP;
  Double_t fVertSPhi;

  Double_t fTof1PathTimePionPlus;
  Double_t fTof1SPPionPlus;
  Double_t fSftSPPionPlus;
  Double_t fVertPathTimePionPlus;
  Double_t fVertSPPionPlus;

  Double_t fTof1PathTimePositron;
  Double_t fTof1SPPositron;
  Double_t fSftSPPositron;
  Double_t fVertPathTimePositron;
  Double_t fVertSPPositron;

  ClassDef(track, 1)
};


class trackArray : public CRTBase{
 public:
  trackArray();
  ~trackArray();

  vector<track> vectTrack;
  ClassDef(trackArray, 2)
};





#endif

