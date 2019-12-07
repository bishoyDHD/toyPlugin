#ifndef __DET_CSI__
#define __DET_CSI__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <singleCsI.h>
#include "findClusters.h"
#include "clusterScore.h"
#include "csitree.h"
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
  //CRTFitCsI *treeFit;			/// Output tree for CSI data
  CRTRawCsI *treeRaw;			/// Input tree with CSI raw data
  CRTClusterCsI *treeClus;			/// Input tree with CSI raw data
  //
  // Detector parameters set in init file
  //
  // constants 
  findClusters* fclusters;
  clusterScore* scoring;
  std::map<IdCsI,UInt_t> mapCsI;
  std::ifstream parfile;
  UInt_t iClock,iFB,iUD,iModule;
  double calibpar[12][2][2][16];
  double mapPhi,pcal;
  int moduleNo, clusCrys,clustM;;
  int multiCrys, singleCrys;
  bool clus_csi;
  std::vector<double> crysE,phval;
  std::map<std::pair<double,double>,double> csiph;
  std::map<std::pair<double,double>,double> csiR;
  std::map<std::pair<double,double>,double> csiZ;
  std::map<std::pair<double,double>,bool> crysChk;
  std::vector<double> clusEne,singleEne,singZ,singR,singTheta,singPhi;
  std::vector<double> clusThetaE, clusPhiE, clusEz, clusEr;
  std::vector<double> csiEdep, csiPhi, csiTheta;
  // indices start at zero now
  int thetaCsI[16][48];
  int phiCsI[16][48];
  // angles to be used by the clusterFinder table
  double otheta, ophi, wtheta, wphi, wz, wr;
  double pr2px, pr2py, pr2pz;
  double cl1px, cl1py, cl1pz;
  double cl2px, cl2py, cl2pz;
  double cl1x,  cl1y,  cl1z, cl1r;
  double cl2x,  cl2y,  cl2z, cl2r;
  double cl1E,cl2E,cl1theta, cl2theta;
  double cl1phi, cl2phi;
  double pr1px, pr1py, pr1pz;
  double pr2x, pr2y, pr2z;
  double pr1Etot;
  double E2clust=0;

public:
  Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_CsI();

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
  Long_t angleCsI(int id, int module, int channel, int yy, int zz);
  //histograms for fit

  TH2D* h2TimeVSCsI,*h2ang;
  TH2D* h2ChargeVSCsI;
  TH1D* h1MaxDiff;
  TH2D* h2DiffVSCsI;
  

  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_CsI,1);
    };

#endif
