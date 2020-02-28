#ifndef __DET_CSI__
#define __DET_CSI__
#include <TObject.h>
#include <TStyle.h>
#include <Plugin.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <singleCsI.h>
#include "findClusters.h"
#include "clusterScore.h"
#include "csitree.h"
#include "multiTree.h"
#include "track.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
typedef std::pair<UInt_t,UInt_t> IdCsI;
class Det_CsI:public Plugin{
private:
  CRTSingleCsI *treeFit;			/// Output tree for CSI data
  CRTRawCsI *treeRaw;			/// Input tree with CSI raw data
  CRTClusterCsI *treeClus;			/// Input tree with CSI raw data
  mwpcE36* mwpcTree;
  trackArray* trackTree;   // Input tree from tracks (MWPCs)
  targetE36* tgtTree;
  //
  // Detector parameters set in init file
  //
  // constants 
  findClusters* fclusters;
  clusterScore* scoring;
  std::map<IdCsI,UInt_t> mapCsI;
  std::ifstream parfile;
  UInt_t iClock,iFB,iUD,iModule, event=0;
  double calibpar[12][2][2][16];
  double mapPhi,pcal;
  int moduleNo, clusCrys,clustM;;
  int multiCrys, singleCrys;
  bool clus_csi;
  std::vector<Double_t> crysE,phval;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiph;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiR;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiZ;
  std::map<std::pair<Double_t,Double_t>,bool> crysChk;
  std::vector<Double_t> clusEne,singleEne,singZ,singR,singTheta,singPhi;
  std::vector<Double_t> clusThetaE, clusPhiE, clusEz, clusEr;
  std::vector<Double_t> csiEdep, csiPhi, csiTheta;
  std::vector<UInt_t> channel;
  std::vector<Double_t> phdiff; // difference of pulse-height
  // indices start at zero now
  int thetaCsI[16][48];
  int phiCsI[16][48];
  // angles to be used by the clusterFinder table
  Double_t otheta, ophi, wtheta, wphi, wz, wr;
  Double_t pr2px, pr2py, pr2pz;
  Double_t cl1px, cl1py, cl1pz;
  Double_t cl2px, cl2py, cl2pz;
  Double_t cl1x,  cl1y,  cl1z, cl1r;
  Double_t cl2x,  cl2y,  cl2z, cl2r;
  Double_t cl1E,cl2E,cl1theta, cl2theta;
  Double_t cl1phi, cl2phi;
  Double_t pr1p,pr1px, pr1py, pr1pz;
  Double_t piPp,piPpx, piPpy, piPpz;
  Double_t pr2x, pr2y, pr2z;
  Double_t pr1Etot;
  Double_t E2clust=0;
  TLorentzVector prim1lv,prim2lv;
  TLorentzVector kaon;
  TVector3 prim1vec3,prim2vec3,gv1;
  double opAngle,prim2px,prim2py,prim2pz;
  double phiMWPC[12]={60,30,0,330,300,270,240,210,180,150,120,90};
  double trackPhi;//=((phiMap[gapNumTof2[j]])*M_PI/180);
  //std::cout<<"  --- Gap number check "<<gapNumTof1[j]<<" --- \n";
  double nxVert,nyVert,nzVert;
  double nXVert,nYVert,nZVert;

public:
  Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_CsI();

  clusters clusVar;
  // Main CsI analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  // CsI fit
  Long_t histos_fit();
  Long_t startup_fit();
  Long_t process_fit();
  Long_t finalize_fit();
  Long_t setIdCsI(std::map<IdCsI,UInt_t>& mapCsI);
  Long_t set_goodEvents(int, int);
  std::vector<std::pair<int,int> > listGoodEvent;
  void readFiles();
  void initVar();
  void initSingleVar();
  Long_t angleCsI(int id, int module, int channel, int yy, int zz);
  Long_t singleAng(int id, int module, int channel, int yy, int zz);
  //histograms for fit
  TH2D* h2TimeVSCsI,*h2ang;
  TH2D* h2ChargeVSCsI;
  TH1D* h1MaxDiff;
  TH1D* h1Ch;
  TH2D* h2DiffVSCsI;
  
  virtual Long_t cmdline(char * cmd);
  
  ClassDef(Det_CsI,1);
};
#endif
