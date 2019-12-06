#ifndef __DET_CSI__
#define __DET_CSI__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <singleCsI.h>
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
  ClusterCsI myClust;
  std::map<IdCsI,UInt_t> mapCsI;
  std::ifstream parfile;
  UInt_t iClock,iFB,iUD,iModule;
  double calibpar[12][2][2][16];
  double mapPhi,pcal;
  int moduleNo;
  // indices start at zero now
  int thetaCsI[16][48];
  int phiCsI[16][48];
  // angles to be used by the clusterFinder table
  double otheta, ophi, wtheta, wphi, wz, wr;
  // crystal center Z:
  double crysZ[20]={-48.3449, -42.0302, -36.5676, -31.5834, -26.851,
                    -20.9203, -15.7210, -10.9616, -6.46940, -2.1341,
                    2.1341, 6.46940, 10.9616, 15.7210, 20.9203,
                    26.851, 31.5834, 36.5676, 42.0302, 48.3449};
  // crystal center r:
  double crysr[20]={16.4109, 20.727, 24.4337, 27.6979, 30.6177,
                    31.3094, 31.879, 32.2917, 32.5240, 32.5595,
                    32.5595, 32.5240, 32.2917, 31.879, 31.3094,
                    30.6177, 27.6979, 24.4337, 20.727, 16.4109};

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

  TH2D* h2TimeVSCsI;
  TH2D* h2ChargeVSCsI;
  TH1D* h1MaxDiff;
  TH2D* h2DiffVSCsI;
  

  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_CsI,1);
    };

#endif
