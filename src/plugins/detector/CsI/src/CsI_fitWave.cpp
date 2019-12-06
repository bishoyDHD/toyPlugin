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
  treeClus=new CRTClusterCsI();
  makeBranch("treeClus",(TObject **) &treeClus);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  readFiles();
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
void Det_CsI::initVar(){
  int dummy=-1000;
  // init 2 gamma variables
  treeClus->cpid1Px=dummy;   treeClus->cpid2Px=dummy;  treeClus->cpid1x=dummy;   treeClus->cpid2x=dummy;
  treeClus->cpid1Py=dummy;   treeClus->cpid2Py=dummy;  treeClus->cpid1y=dummy;   treeClus->cpid2y=dummy;
  treeClus->cpid1Pz=dummy;   treeClus->cpid2Pz=dummy;  treeClus->cpid1z=dummy;   treeClus->cpid2z=dummy;
  // init pion vars
  treeClus->prim1px=dummy;       treeClus->prim1px=dummy;
  treeClus->prim1py=dummy;       treeClus->prim1py=dummy;
  treeClus->prim1pz=dummy;       treeClus->prim1pz=dummy;
  treeClus->clCosTheta=dummy;    treeClus->prCosTheta=dummy;
  // WaveID and cluster var
  treeClus->waveID=dummy;
  treeClus->dubP_1=dummy;
  treeClus->channel=dummy;
  treeClus->clusterM=dummy;      treeClus->E_prim2=dummy;
  treeClus->ClustCrys=dummy;     treeClus->M_prim2=dummy;
  treeClus->Ncrys=dummy;         treeClus->prim2M2=dummy;
  treeClus->cpid1thetaE=dummy;   treeClus->M_k=dummy;
  treeClus->cpid1phiE=dummy;     treeClus->kM2=dummy;
  treeClus->cpid2thetaE=dummy;   treeClus->cpid1E=dummy;
  treeClus->cpid2phiE=dummy;     treeClus->cpid2E=dummy;
  //E2clust=-1000;
}
// Method to read in external files such as calibration par's
void Det_CsI::readFiles(){
  string fileName="calibPar.txt";
  double calib;
  int iclock, ifb, iud, imodule;
  parfile.open(fileName);
  if(parfile){
    for(int iClock=0;iClock<12;iClock++){
      for(int iFB=0;iFB<2;iFB++){
        for(int iUD=0;iUD<2;iUD++){
          for(int iModule=0;iModule<16;iModule++){
            parfile>>iclock>>ifb>>iud>>imodule>>calib;
            calibpar[iclock][ifb][iud][imodule]=calib;
          }
        }
      }
    }
    parfile.close();
  }else{
    std::cerr<<" ****** Error opening file ''"<<fileName<<".'' Please make sure you have the file \n";
    std::abort();
  }
}
Long_t Det_CsI::angleCsI(int id, int module, int channel, int yy, int zz){
  //std::cout<<"\n --- "<<id<<"-"<<module<<"-"<<channel<<"-"<<yy<<"-"<<zz<<" ---\n";
  ClusterCsI myCsI;
  myCsI.setAngles(module,channel,yy,zz);
  //phiCsI[module][channel]=yy;
  //thetaCsI[module][channel]=zz;
  //std::cout<<"\n --- thetaCsI["<<module<<"]["<<channel<<"] = "<<yy<<" <---> "<<thetaCsI[module][channel]<<" ---\n";

  return 0;
}
Long_t Det_CsI::process_fit(){
  if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
  }
  initVar();
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
      return 0;
    }
    #pragma omp parallel num_threads(8)
    iClock=indexClock; iModule=treeRaw->indexCsI[iCh]-1;
    //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
    //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;

    if(!(treeRaw->indexCsI[iCh]==16 && iFB==0 && iUD==0) ||
       !(indexClock==0 || indexClock==2 || indexClock==4 ||
         indexClock==6 || indexClock==8 || indexClock==10)){
      if(treeRaw->indexCsI[iCh]==16) pcal=.22/1000.;
      IdCsI myIndex(treeRaw->nameCsI[iCh],treeRaw->indexCsI[iCh]);
      UInt_t iCsI=mapCsI[myIndex];
      //    if(iCsI==144 || iCsI==400||iCsI==656||iCsI==166||iCsI==128) continue;
      //SingleCsI myCsI(treeRaw->runNo,myEvent,iModule);
      // convert TKO module to module number here
      if(nameModule=="TK34") moduleNo=1;
      if(nameModule=="TK32") moduleNo=2;
      if(nameModule=="TK36") moduleNo=3;
      if(nameModule=="TK37") moduleNo=4;
      if(nameModule=="TK38") moduleNo=5;
      if(nameModule=="TK08") moduleNo=6;
      if(nameModule=="TK50") moduleNo=7;
      if(nameModule=="TK09") moduleNo=8;
      if(nameModule=="TK54") moduleNo=9;
      if(nameModule=="TK31") moduleNo=10;
      if(nameModule=="TK04") moduleNo=11;
      if(nameModule=="TK45") moduleNo=12;
      if(nameModule=="TK33") moduleNo=13;
      if(nameModule=="TK39") moduleNo=14;
      if(nameModule=="TK40") moduleNo=15;
      if(nameModule=="TK41") moduleNo=16;
      std::cout<<" naming TKO module:  "<<nameModule<<std::endl;
      ClusterCsI myCsI(treeRaw->runNo,myEvent,iModule);
      //std::cout<<" top U/D || F/B :"<<p[0]<<", "<<p[1]<<std::endl;
      myCsI.nameCsI(iClock,iFB,iUD,iModule);
      myCsI.setData(treeRaw->data[iCh]);
      myCsI.fit();
    }// end of signal CsI if-loop
  }
  return 0;
}
Long_t Det_CsI::finalize_fit(){
  cout<<"finalize_fit"<<endl;
  return 0;
}
