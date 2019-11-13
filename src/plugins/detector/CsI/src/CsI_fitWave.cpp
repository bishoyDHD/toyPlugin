#include <Det_CsI.h>
#include <singleCsI.h>
#include<cmath>
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
Long_t Det_CsI::set_goodEvents(int run, int event)
{
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
  return 0;
}
Double_t firstDerive(Double_t *x,Double_t *par){


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


}
Long_t Det_CsI::process_fit(){
  if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
  }
  for(UInt_t iCh=0;iCh<treeRaw->nChannel;iCh++){
    UInt_t myEvent=treeRaw->eventNo;
    std::stringstream ss;

    IdCsI myIndex(treeRaw->nameCsI[iCh],treeRaw->indexCsI[iCh]);
    UInt_t iCsI=mapCsI[myIndex];
    //    if(iCsI==144 || iCsI==400||iCsI==656||iCsI==166||iCsI==128) continue;
    SingleCsI myCsI(treeRaw->runNo,myEvent,iCsI);
    myCsI.setData(treeRaw->data[iCh]);
    myCsI.fit();
  }
  return 0;
}
Long_t Det_CsI::finalize_fit(){
  cout<<"finalize_fit"<<endl;
  return 0;
}
