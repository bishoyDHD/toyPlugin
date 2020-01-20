#ifndef __COVFEFE__
#define __COVFEFE__

#include <marinateCsI.h>
#include "evalClusters.h"
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLine.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1D.h>
#include <iostream>
#include <map>
#include <iterator>

class covfefe:public Plugin{
 private:
  CATSingleCsI* treeCalib;
  CATClusterCsI* clsmar;
  CATCaliCsI* calibcsi;
  const double dE=143.5;
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  double sigdE=3.25;
  double T_0[3], tcorr[3], trefCorr[3], cdf50[3], refgaus[3];
  double phi[12]={30.,60.,0.,90.,120.,150.,180.,210.,240.,270.,300.,330.};
  double csiphi;
  double lowRange, upRange, apcsi;
  evalClusters* clustEval;
 public:
  covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~covfefe();
  int iclock, iModule;
  Int_t iUD, iFB, nbins;
  Double_t adcVal, intVal;
  TH1D* calibHist;
  TH1D* integHist;
  TH1D* tpeak;
  TH1D* Ecorr;
  TH1D* intEn;
  TH1D* phdis;
  TH1D* hkmu2;
  TH1D* timing, *phdistr;
  TH1D* tof1ang;
  TH2D* csiAng;
  Int_t index;
  TH1D* h1time[12][2][2][16];
  TH1D* h1cali[12][2][2][16];
  // add funtions with return value Long_t here:
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  // functions for cluster analysis:
  Long_t hist_clust();
  Long_t startup_clust();
  Long_t process_clust();
  Long_t finalize_clust();
  // histograms for cluster analysis
  TH1D* E_pi0[3], *M_pi0[3]; // pi0 total energy
  TH1D* waveID[2], *clustM[3], *id1;
  TH1D* g1px[2], *g1py[2], *g1pz[2];
  TH1D* g2px[2], *g2py[2], *g2pz[2];
  TH1D* pi0px[2], *pi0py[2], *pi0pz[2];
  TH1D* vertpx[2], *vertpy[2], *vertpz[2];
  TH2D* kmass[2], *h2Angle[2];
  // angles
  TH1D *h1theta[2], *h1phi[2];
  TH1D *angPP[3], *piPang[3], *pi0ang[3];
  TH2D *h2corrAng;
  // timining histograms
  TH1D* h1t_0[3], *h1tcorr[3],*h1trefCorr[3],*h1cdf50[3],*h1refgaus[3];
  TLorentzVector prim1lv,prim2lv;

  virtual Long_t cmdline(char * cmd);

  ClassDef(covfefe,1);
};
#endif
