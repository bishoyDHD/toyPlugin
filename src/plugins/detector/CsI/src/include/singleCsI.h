#ifndef __SingleCsI_H
#define __SingleCsI_H 1
#include <memory>
#include <map>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include <TMath.h>
#define _USE_MATH_DEFINES
using std::vector;
using std::cout;
using std::endl;
using std::shared_ptr;
using std::unique_ptr;
class WaveformCsI{
public:
  UInt_t mNWave;
public:
  WaveformCsI(UInt_t nWave){mNWave=nWave;}
  Double_t waveformSingle(Double_t *x,Double_t *par);//waveform of single wave no background
  Double_t waveformDouble(Double_t *x,Double_t *par);//waveform with pile-up signals
  Double_t waveformTriple(Double_t *x,Double_t *par);//waveform with pile-up signals
  Double_t waveformQuad(Double_t *x,Double_t *par);//waveform with pile-up signals
  Double_t waveformOverrange(Double_t *x,Double_t *par);//waveform of overrange signals
  Double_t waveformCut(Double_t *x,Double_t *par);//waveform of multiple wave w/ background
};
class SingleCsI{
private:
  // cluster map
  std::map<std::pair<double,double>,double> csiClust;
  std::map<std::pair<double,double>,bool> csiCheck;
  UInt_t mRunNo;
  UInt_t mEventNo;
  UInt_t mIndexCsI;
  UInt_t mNWave;
  Double_t mPar[36];
  Int_t x1,x2,diffMax;
  vector<Double_t> mListData;
  vector<Double_t> mListTime;
  vector<Double_t> mListEnergy;
  char mName[7];
  vector<Double_t> mListLocalMax;
private:
  Double_t const dummy=-1000;
  Double_t findChi2(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>,double* xval,double yped,double ymax);
  void tryFit(TH1D*,double* xval,double yped,double ymax);
  void findLocalMax(shared_ptr<TH1D>);
  void drawWaves(shared_ptr<TH1D> h1);
  void drawWaves(TH1D* h1);
  void drawWaves(shared_ptr<TH1D> h1,shared_ptr<TF1> f1);
  Double_t chi2,ndf;
  Double_t mapPhi,pcal,Theta,acos,energy;
  Double_t lmax,lmin,ymax2,ymax3;
  Double_t ptime,rtime,cdf50;
  Double_t phdiff,tcalc;
  Double_t otheta, ophi, wtheta, wphi, wz, wr;
  int moduleNo,thetaIndex,phiIndex;
  // crystal center Z:
  Double_t crysZ[20]={-48.3449, -42.0302, -36.5676, -31.5834, -26.851,
                    -20.9203, -15.7210, -10.9616, -6.46940, -2.1341,
                    2.1341, 6.46940, 10.9616, 15.7210, 20.9203,
                    26.851, 31.5834, 36.5676, 42.0302, 48.3449};
  // crystal center r:
  Double_t crysr[20]={16.4109, 20.727, 24.4337, 27.6979, 30.6177,
                    31.3094, 31.879, 32.2917, 32.5240, 32.5595,
                    32.5595, 32.5240, 32.2917, 31.879, 31.3094,
                    30.6177, 27.6979, 24.4337, 20.727, 16.4109};
  // Parameters for the waveform fitting function
  double param[31]={1000, 35.76, 26.68, 19.85, 15.83, 0.065, 2.255, 31.21,120,
                    120.5, 800, 700., 200.1, 17.1, 0.065, 2.255, 31.21,
                    200.5, 600, 28.9, 18.0, 17.1, 0.065, 2.255, 31.21,
                    25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1};
  double parUplim[31]={1e5, 70.1, 60, 45, 20, 1.07, 4, 90, 350,
                    250, 1e5, 80, 70, 50, 1.07, 1000, 250};
  double parLowlim[31]={80.0, 15.1, 1, 5, 0., 1e-4, 0., 10,
                    70.1, 15, 10, 10, 10, 1e-4, 10, 20};
public:
  SingleCsI(){}
  ~SingleCsI(){}
  SingleCsI(UInt_t runNo){mRunNo=runNo;}
  SingleCsI(UInt_t runNo,UInt_t eventNo){mRunNo=runNo;eventNo=mEventNo;}
  SingleCsI(UInt_t runNo,UInt_t eventNo,UInt_t index){mRunNo=runNo;mEventNo=eventNo;mIndexCsI=index;}
  inline void setRunNo(UInt_t runNo){mRunNo=runNo;}
  inline void setEventNo(UInt_t eventNo){mEventNo=eventNo;}
  inline void setIndex(UInt_t index){mIndexCsI=index;}
  inline void setCsI(UInt_t index){mIndexCsI=index;}
  inline void addData(UShort_t sample){mListData.push_back(sample);}
  void setData(const vector<UShort_t>&);
  inline void setIndexTheta(int iTheta){thetaIndex=iTheta;}
  inline void setIndexPhi(UInt_t iPhi){phiIndex=iPhi;}
  void calcThetaPhi(double);
  void calTime(shared_ptr<TF1> f1);
  void initVar(); //initialize
  inline void dumbFn(){std::cout<<" -----> Hola 7 phezulu ===== "<<wtheta<<", "<<wphi<<std::endl;}
  void setAngles(int module,int channel,int yy,int zz);
  Double_t getTheta(){return wtheta;}
  Double_t getPhi(){return wphi;}
  Double_t getEdep(){return energy;}
  Double_t getphDiff(){return phdiff;}
  Double_t getR(){return wr;}
  Double_t getZ(){return wz;}
  Double_t getTime(){return rtime;}
  Double_t getpTime(){return ptime;}
  Double_t getCDF50(){return cdf50;}
  Double_t getChi2(){return chi2;}
  Double_t getNDF(){return ndf;}
  UInt_t numberWave() const{
    return mNWave;
  }
  UInt_t runNo() const{
    return mRunNo;
  }
  UInt_t eventNo() const{
    return mEventNo;
  }
  bool fit();
  char* nameCsI(unsigned int index);
  char* nameCsI(const UInt_t& iClock, const UInt_t& iFB,const UInt_t& iUD,const UInt_t& iModule);
};
class ClusterCsI{
private:
  UInt_t mRunNo;
  UInt_t mEventNo;
  UInt_t mIndexCsI;
  UInt_t ichan;
  UInt_t mNWave;
  Double_t mPar[36];
  vector<Double_t> mListData;
  vector<Double_t> mListTime;
  vector<Double_t> mListEnergy;
  char mName[7];
  vector<Double_t> mListLocalMax;
private:
  int clusCrys;
  bool clus_csi;
  Double_t phdiff,tcalc;
  std::vector<Double_t> crysE;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiph;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiR;
  std::map<std::pair<Double_t,Double_t>,Double_t> csiZ;
  std::map<std::pair<Double_t,Double_t>,bool> csiClus;
  Double_t findChi2(shared_ptr<TH1D>);
  Double_t round(double val);
  void tryFit(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>,double* xval,double yped,double ymax);
  void tryFit(TH1D*,double* xval,double yped,double ymax);
  void findLocalMax(shared_ptr<TH1D>);
  void drawWaves(shared_ptr<TH1D> h1);
  void drawWaves(TH1D* h1);
  void drawWaves(shared_ptr<TH1D> h1,shared_ptr<TF1> f1);
  Double_t mapPhi,pcal,Theta,acos,energy;
  Double_t lmax,lmin;
  Double_t ptime,rtime,cdf50;
  int moduleNo,thetaIndex,phiIndex;
  std::pair<Double_t,Double_t> angles;
  TH2D* h2ang;
  // indices start at zero now
  int thetaCsI[16][48];
  int phiCsI[16][48];
  // angles to be used by the clusterFinder table
  Double_t otheta, ophi, wtheta, wphi, wz, wr;
  // crystal center Z:
  Double_t crysZ[20]={-48.3449, -42.0302, -36.5676, -31.5834, -26.851,
                    -20.9203, -15.7210, -10.9616, -6.46940, -2.1341,
                    2.1341, 6.46940, 10.9616, 15.7210, 20.9203,
                    26.851, 31.5834, 36.5676, 42.0302, 48.3449};
  // crystal center r:
  Double_t crysr[20]={16.4109, 20.727, 24.4337, 27.6979, 30.6177,
                    31.3094, 31.879, 32.2917, 32.5240, 32.5595,
                    32.5595, 32.5240, 32.2917, 31.879, 31.3094,
                    30.6177, 27.6979, 24.4337, 20.727, 16.4109};

  // Parameters for the waveform fitting function
  Double_t param[31]={1000, 35.76, 26.68, 19.85, 15.83, 0.065, 2.255, 31.21,120,
                    120.5, 800, 700., 200.1, 17.1, 0.065, 2.255, 31.21,
                    200.5, 600, 28.9, 18.0, 17.1, 0.065, 2.255, 31.21,
                    25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1};
  Double_t parUplim[31]={1e5, 70.1, 60, 45, 20, 1.07, 4, 90, 350,
                    250, 1e5, 80, 70, 50, 1.07, 1000, 250};
  Double_t parLowlim[31]={80.0, 15.1, 1, 5, 0., 1e-4, 0., 10,
                    70.1, 15, 10, 10, 10, 1e-4, 10, 20};
public:
  ClusterCsI(){}
  ~ClusterCsI(){}
  ClusterCsI(UInt_t runNo){mRunNo=runNo;}
  ClusterCsI(UInt_t runNo,UInt_t eventNo){mRunNo=runNo;eventNo=mEventNo;}
  ClusterCsI(UInt_t runNo,UInt_t eventNo,UInt_t index){mRunNo=runNo;mEventNo=eventNo;mIndexCsI=index;}
  inline void setRunNo(UInt_t runNo){mRunNo=runNo;}
  inline void setEventNo(UInt_t eventNo){mEventNo=eventNo;}
  inline void setIndex(UInt_t index){mIndexCsI=index;}
  inline void setCsI(UInt_t index){mIndexCsI=index;}
  inline void addData(UShort_t sample){mListData.push_back(sample);}
  inline void setIndexTheta(int iTheta){thetaIndex=iTheta;}
  inline void setIndexPhi(UInt_t iPhi){phiIndex=iPhi;}
  void calcThetaPhi(double);
  void calTime(shared_ptr<TF1> f1);
  inline void dumbFn(){std::cout<<" -----> Hola 7 phezulu ===== "<<wtheta<<", "<<wphi<<std::endl;}
  void setAngles(int module,int channel,int yy,int zz);
  inline void setpCalib(double calib){pcal=calib/1000.;}
  void setData(const vector<UShort_t>&);
  // cluster vars
  bool getClustCsI(){return clus_csi;}
  inline void setClustCsI(bool val){clus_csi=val;}
  Double_t getTheta(){return wtheta;}
  Double_t getPhi(){return wphi;}
  Double_t getEdep(){return energy;}
  Double_t getphDiff(){return phdiff;}
  Double_t getR(){return wr;}
  Double_t getZ(){return wz;}
  Double_t getTime(){return rtime;}
  Double_t getpTime(){return ptime;}
  Double_t getCDF50(){return cdf50;}
  UInt_t numberWave() const{
    return mNWave;
  }
  UInt_t runNo() const{
    return mRunNo;
  }
  UInt_t eventNo() const{
    return mEventNo;
  }
  bool fit();
  char* nameCsI(unsigned int index);
  char* nameCsI(const UInt_t& iClock, const UInt_t& iFB,const UInt_t& iUD,const UInt_t& iModule);
};
#endif
