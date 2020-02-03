#ifndef __COVFEFE__
#define __COVFEFE__

#include <marinateCsI.h>
#include "evalClusters.h"
#include "track.h"
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
  CATTimingCsI* treeTime;
  CATClusterCsI* clsmar;
  CATCaliCsI* calibcsi;
  trackArray *trackArr;
  const double dE=143.5;
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  double sigdE=3.25;
  double T_0[3], tcorr[3], trefCorr[3], cdf50[3], refgaus[3];
  double phi[12]={30.,60.,0.,90.,120.,150.,180.,210.,240.,270.,300.,330.};
  double csiphi;
  double lowRange, upRange, apcsi;
  evalClusters* clustEval;
  double t_ref1,t_ref2,t_ref3;
  double t_corr, tCalc;
  double pVert, piPpx, piPpy, piPpz, P, primCosTheta;
 public:
  covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~covfefe();
  int iclock, iModule;
  Int_t iUD, iFB, nbins;
  Double_t adcVal, intVal;
  TH1D* calibHist;
  TH1D* integHist;
  TH1D* tpeak;
  TH1D* Ecorr,*Eclust,*InvM;
  TH1D* intEn;
  TH1D* phdis;
  TH1D* hkmu2,*h1Pkpi2;
  TH1D* timing, *phdistr;
  TH1D* tof1ang,*h1angDiff,*h1score;
  TH1D* angScore,*h1tof1,*tof1ID;
  TH2D* csiAng;
  Int_t index;
  TH1D* h1time[12][2][2][16];
  TH1D* h1kmu2[12][2][2][16];
  TH1D* h1cali[12][2][2][16];
  // add funtions for calibration analysis
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  // functions for sanity analysis:
  Long_t histosCsI();
  Long_t startupCsI();
  Long_t processCsI();
  Long_t finalizeCsI();
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
  double Tconst[3]={6.60068,6.58848,6.54471};
  // angles
  TH1D *h1theta[2], *h1phi[2];
  TH1D *angClust[3],*prAng[3],*h1Mass[3],*h1Eclust[3],*h1P[3];
  TH2D *h2corrAng;
  TH1D*h1prmAng;
  // timining histograms
  TH1D* h1t_0[3], *h1tcorr[3],*h1trefCorr[3],*h1cdf50[3],*h1refgaus[3];
  TH1D* h1Tcorr,*h1diffTof2;
  TLorentzVector prim1lv,prim2lv;
  TVector3 prim1vec3,prim2vec3;

  virtual Long_t cmdline(char * cmd);

  ClassDef(covfefe,1);
};
#endif
