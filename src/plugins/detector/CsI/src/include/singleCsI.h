#ifndef __SingleCsI_H
#define __SingleCsI_H 1
#include <memory>
#include <map>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
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
  vector<Double_t> mListData;
  vector<Double_t> mListTime;
  vector<Double_t> mListEnergy;
  char mName[7];
  vector<Double_t> mListLocalMax;
private:
  Double_t findChi2(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>,double* xval,double yped,double ymax);
  void tryFit(TH1D*,double* xval,double yped,double ymax);
  void findLocalMax(shared_ptr<TH1D>);
  void drawWaves(shared_ptr<TH1D> h1);
  void drawWaves(TH1D* h1);
  void drawWaves(shared_ptr<TH1D> h1,shared_ptr<TF1> f1);
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
  UInt_t mNWave;
  Double_t mPar[36];
  vector<Double_t> mListData;
  vector<Double_t> mListTime;
  vector<Double_t> mListEnergy;
  char mName[7];
  vector<Double_t> mListLocalMax;
private:
  Double_t findChi2(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>);
  void tryFit(shared_ptr<TH1D>,double* xval,double yped,double ymax);
  void tryFit(TH1D*,double* xval,double yped,double ymax);
  void findLocalMax(shared_ptr<TH1D>);
  void drawWaves(shared_ptr<TH1D> h1);
  void drawWaves(TH1D* h1);
  void drawWaves(shared_ptr<TH1D> h1,shared_ptr<TF1> f1);
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
  inline void setModuleNo(int mNum){moduleNo=mNum;}
  inline void setAngles(int module,int channel,int yy,int zz){
    std::cout<<" cluster crystal init var--> \n";
    phiCsI[module][channel]=yy;
    thetaCsI[module][channel]=zz;
    std::cout<<"\n --- thetaCsI["<<module<<"]["<<channel<<"] = "<<yy<<" <---> "<<thetaCsI[module][channel]<<" ---\n";
  }
  inline void setpCalib(double calib){pcal=calib;}
  void setData(const vector<UShort_t>&);
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
