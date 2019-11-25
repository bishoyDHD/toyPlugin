#include <vector>
#include "TH1D.h"
#include "TF1.h"
using std::vector;
using std::cout;
using std::endl;
class WaveformCsI{
public:
  UInt_t mNWave;
public:
  WaveformCsI(UInt_t nWave){mNWave=nWave;}
  Double_t waveformSingle(Double_t *x,Double_t *par);//waveform of single wave no background
  Double_t waveformDouble(Double_t *x,Double_t *par);//waveform with pile-up signals
  Double_t waveformTriple(Double_t *x,Double_t *par);//waveform with pile-up signals
  Double_t waveformOverrange(Double_t *x,Double_t *par);//waveform of overrange signals
  Double_t waveformCut(Double_t *x,Double_t *par);//waveform of multiple wave w/ background
};
class SingleCsI{
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
  Double_t findChi2(TH1D*);
  void tryFit(TH1D*);
  void findLocalMax(TH1D*);
  void drawWaves(TH1D* h1);
public:
  SingleCsI(){}
  ~SingleCsI(){}
  SingleCsI(UInt_t runNo){mRunNo=runNo;}
  SingleCsI(UInt_t runNo,UInt_t eventNo){mRunNo=runNo;eventNo=mEventNo;}
  SingleCsI(UInt_t runNo,UInt_t eventNo,UInt_t index){mRunNo=runNo;mEventNo=eventNo;mIndexCsI=index;}
  void setRunNo(UInt_t runNo){mRunNo=runNo;}
  void setEventNo(UInt_t eventNo){mEventNo=eventNo;}
  void setIndex(UInt_t index){mIndexCsI=index;}
  void setCsI(UInt_t index){mIndexCsI=index;}
  void addData(UShort_t sample){mListData.push_back(sample);}
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
};
