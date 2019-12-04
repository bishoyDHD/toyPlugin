#include <Det_CsI.h>
#include <singleCsI.h>
#include <cmath>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TF1.h"
#include <cstdlib>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
using namespace std;
Long_t Det_CsI::set_goodEvents(int run, int event){
  cout<<"set good events..."<<endl;
  return 0;
};

Long_t Det_CsI::histos_fit(){
  h2TimeVSCsI=dH2("TimeVSCsI","fitted time;CsI;time",768,0.5,768.5,250,0.0,250.0);
  h2ChargeVSCsI=dH2("TimeVSCsI","fitted charge;CsI;time",768,0.5,768.5,2500,0.0,250000.0);
  h1MaxDiff=dH1("MaxDiff","max difference in a fit;diff",100,0.0,1000.0);
  h2DiffVSCsI=dH2("DiffVSCsI","max difference;CsI;diff",768,0.5,768.5,100,0.0,1000.0); 
  return 0;
}
Long_t Det_CsI::startup_fit(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  std::cout<<" -- Hell Ya! Let's get it started right here!!\n";
  return 0;
}
Double_t firstDerive(Double_t *x,Double_t *par){

  return 0;
}
Double_t waveform(Double_t *x,Double_t *par){
  //  if(x[0]<1 && x[0]>250) return 0;
  Double_t termFirst=par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[7])));
  Double_t termSecond=(x[0]-par[1])/(par[3]*par[3]);
  Double_t termThird=TMath::Exp(-(x[0]-par[1])/par[3])+par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t value=termFirst*termSecond*termThird+par[6];
  if(value>1023.0) value=1023.0;
  return value;
}
Double_t findSteepestSlope(Double_t *par){

  return 0;
}
Long_t Det_CsI::process_fit(){
  if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
  }
  for(UInt_t iCh=0;iCh<treeRaw->nChannel;iCh++){
    UInt_t myEvent=treeRaw->eventNo;
    std::stringstream ss;
    char* p=(char*)&(treeRaw->nameModule[iCh]);
    int moduleName=(p[3]-'0')*10+(p[2]-'0')-1;
    std::string nameModule;
    nameModule+=(*p);
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    std::string nameCsI;
    p=(char*)&(treeRaw->nameCsI[iCh]);
    int indexClock=(p[3]-'0')*10+(p[2]-'0')-1;
    //std::cout<< "\n Index clock 11: "<<indexClock<<endl;
    p+=3;
    nameCsI+=(*p);
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;
    iFB=0;
    if(p[1]=='b' || p[1]=='B') iFB=1;
    iUD=0;
    if(p[0]=='d' || p[0]=='D') iUD=1;
    if(p[1]=='0' || p[0]=='0'){
      std::cout<<"..... Got one!!\n";
    }
    #pragma omp parallel num_threads(8)
    iClock=indexClock; iModule=treeRaw->indexCsI[iCh]-1;
    //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
    //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;

    IdCsI myIndex(treeRaw->nameCsI[iCh],treeRaw->indexCsI[iCh]);
    UInt_t iCsI=mapCsI[myIndex];
    //    if(iCsI==144 || iCsI==400||iCsI==656||iCsI==166||iCsI==128) continue;
    SingleCsI myCsI(treeRaw->runNo,myEvent,iModule);
    //std::cout<<" top U/D || F/B :"<<p[0]<<", "<<p[1]<<std::endl;
    myCsI.nameCsI(iClock,iFB,iUD,iModule);
    myCsI.setData(treeRaw->data[iCh]);
    myCsI.fit();
  }
  return 0;
}
Long_t Det_CsI::finalize_fit(){
  cout<<"finalize_fit"<<endl;
  return 0;
}
