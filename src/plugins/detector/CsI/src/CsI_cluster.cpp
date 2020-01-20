#include <Det_CsI.h>
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
  h2ang=dH2("h2deg", "Angular distribution #phi vs #theta deg", 24.0,0.,180., 48.0,0.,360.);
  return 0;
}
Long_t Det_CsI::startup_fit(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("mwpcInfo",(TObject **) &mwpcTree);
  getBranchObject("tgtInfo",(TObject **) &tgtTree);
  treeClus=new CRTClusterCsI();
  makeBranch("treeClus",(TObject **) &treeClus);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  fclusters=new findClusters();
  scoring=new clusterScore();
  readFiles();
  scoring->init();
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
  treeClus->isBad=dummy;
  //treeClus->extraTOF1=dummy;
  channel.clear();               phdiff.clear();
  csiph.clear();   phval.clear();   crysChk.clear();
  csiR.clear();    csiZ.clear();
  clusEne.clear(); singleEne.clear();
  singZ.clear();   singR.clear();
  singTheta.clear(); singPhi.clear();
  clusThetaE.clear(); clusPhiE.clear();
  clusEz.clear();     clusEr.clear();
  csiEdep.clear();    csiPhi.clear();
  csiTheta.clear();
  clus_csi=false;
  E2clust=-1000;
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
  phiCsI[module][channel]=yy;
  thetaCsI[module][channel]=zz;
  std::cout<<"\n --- thetaCsI["<<module<<"]["<<channel<<"] = "<<yy<<" <---> "<<thetaCsI[module][channel]<<" ---\n";

  return 0;
}
Long_t Det_CsI::process_fit(){
  // Save relavant run/event info
  treeClus->eventNo=treeRaw->eventNo;
  treeClus->runNo=treeRaw->runNo;
  if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
  }
  // need to make sure this called for every event
  // avoid memory leaks
  initVar();
  if(mwpcTree->run != treeRaw->runNo && tgtTree->run != treeRaw->runNo){
    std::cout<<" Oops you are comparing to different runs \n";
    std::cout<<" ***Bailing*** Bailing*** \n";
    std::abort();
  }
  // make sure that this is a good gap event: isBad>0 (=1)
  if(mwpcTree->nTracks==0 || mwpcTree->fVertSP==-10000) return 0;
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
    if(p[1]=='\0' || p[0]=='\0'){
      std::cout<<"..... Got one!!\n";
      return 0;
    }
    #pragma omp parallel num_threads(16)
    iClock=indexClock; iModule=treeRaw->indexCsI[iCh]-1;
    //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
    //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
    // start of signal CsI if-loop
    if(!(treeRaw->indexCsI[iCh]==16 && iFB==0 && iUD==0) ||
       !(indexClock==0 || indexClock==2 || indexClock==4 ||
         indexClock==6 || indexClock==8 || indexClock==10)){
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
      //std::cout<<" checking the hell outta this shit:"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[iCh]]<<std::endl;
      int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      myCsI.nameCsI(iClock,iFB,iUD,iModule);
      myCsI.setClustCsI(false);
      myCsI.setData(treeRaw->data[iCh]);
      myCsI.setIndexTheta(thetaIndex);
      myCsI.setIndexPhi(phiIndex);
      pcal=calibpar[indexClock][iFB][iUD][iModule];
      if(treeRaw->indexCsI[iCh]==16) pcal=.22;
      myCsI.setpCalib(pcal);
      myCsI.fit();
      if(myCsI.getClustCsI()){
	clus_csi=true;
        std::cout<<" ********** "<<myCsI.getR()<<" || "<<myCsI.getTheta()<<std::endl;
	auto angles=std::make_pair(myCsI.getTheta(),myCsI.getPhi());
	csiph[angles]=myCsI.getEdep();
        crysChk[angles]=true;
        csiR[angles]=myCsI.getR();
        csiZ[angles]=myCsI.getZ();
	csiEdep.push_back(myCsI.getEdep());
	csiTheta.push_back(myCsI.getTheta());
	csiPhi.push_back(myCsI.getPhi());
	//if(myCsI.getTheta()<16.2) std::abort();
	h2ang->Fill(myCsI.getTheta(),myCsI.getPhi());
	// Fill temporary CsI information
	// Needed for future re-analysis
	//treeClus->tsig[angles]=myCsI.getCDF50();
	//treeClus->csiph[angles]=myCsI.getEdep();
        //treeClus->crysChk[angles]=true;
        //treeClus->csiR[angles]=myCsI.getR();
        //treeClus->csiZ[angles]=myCsI.getZ();
	//treeClus->csiEdep.push_back(myCsI.getphDiff());
	//treeClus->csiTheta.push_back(myCsI.getTheta());
	//treeClus->csiPhi.push_back(myCsI.getPhi());
      }
      // not sure this is needed... just in case
      myCsI.setClustCsI(false);
    }// end of signal CsI if-loop
  }
  if(clus_csi){
    // find clusters and perform cluster scoring and analysis
    // find clusters
    fclusters->clusters(csiph,csiEdep,csiTheta,csiPhi,csiR,csiZ,crysChk);
    multiCrys=fclusters->getMultiCrysClust();
    singleCrys=fclusters->getSingleCrysClust();
    // obtain Multi-cluster E, theta & phi
    clusEne=fclusters->getMultiE();
    clusThetaE=fclusters->getMultiTheta();
    clusPhiE=fclusters->getMultiPhi();
    clusEz=fclusters->getMultiZ();
    clusEr=fclusters->getMultiR();
    // obtain single cluster E, theta & phi
    singleEne=fclusters->getSingleE();
    singTheta=fclusters->getSingleTheta();
    singPhi=fclusters->getSinglePhi();
    singZ=fclusters->getSingleZ();
    singR=fclusters->getSingleR();
    TLorentzVector prim1lv,prim2lv;
    TLorentzVector kaon;
    TVector3 prim1vec3,prim2vec3,gv1;
    double opAngle,prim2px,prim2py,prim2pz;
    scoring->setScoreMass(0.09);
    if((multiCrys>=4 || singleCrys>=4)) goto exitFill;
    //if((multiCrys+singleCrys>=4)) goto exitFill;
    if(multiCrys==2 && singleCrys==0){
      std::cout<<"  This is only 2 multiCrys ---|\n";
      scoring->init();
      //scoring->clusterEval(clusEne,clusThetaE,clusPhiE);
      scoring->clusterEval(clusEne,clusEr,clusEz,clusThetaE,clusPhiE);
      pr2px=scoring->getprPx();   pr2py=scoring->getprPy();  pr2pz=scoring->getprPz();
      E2clust=scoring->getE();
      // cluster PID 1:
      scoring->setCpid(1); //NOTE: must always set this
      cl1px=scoring->getclPx();
      cl1py=scoring->getclPy();
      cl1pz=scoring->getclPz();
      cl1E=scoring->getclE();
      cl1theta=scoring->getclTheta();
      cl1phi=scoring->getclPhi();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      // cluster PID 1:
      scoring->setCpid(2);
      cl2px=scoring->getclPx();
      cl2py=scoring->getclPy();
      cl2pz=scoring->getclPz();
      cl2E=scoring->getclE();
      cl2theta=scoring->getclTheta();
      cl2phi=scoring->getclPhi();
      clustM=scoring->getClustM();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      opAngle=scoring->getOpAngleClust();
      prim2lv=scoring->getprimLV();
      //ThreeVector for angular analysis
      prim1vec3.SetXYZ(pr1px,pr1py,pr1pz);
      prim2vec3.SetXYZ(pr2px,pr2py,pr2pz);
    }else
    if(multiCrys==0 && singleCrys==2){
      std::cout<<"  This is only 2 singleCrys ---|\n";
      scoring->init();
      //scoring->clusterEval(singleEne,singTheta,singPhi);
      scoring->clusterEval(singleEne,singR,singZ,singTheta,singPhi);
      pr2px=scoring->getprPx();   pr2py=scoring->getprPy();  pr2pz=scoring->getprPz();
      E2clust=scoring->getE();
      // cluster PID 1:
      scoring->setCpid(1); //NOTE: must always set this
      cl1px=scoring->getclPx();
      cl1py=scoring->getclPy();
      cl1pz=scoring->getclPz();
      cl1E=scoring->getclE();
      cl1theta=scoring->getclTheta();
      cl1phi=scoring->getclPhi();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      // cluster PID 1:
      scoring->setCpid(2);
      cl2px=scoring->getclPx();
      cl2py=scoring->getclPy();
      cl2pz=scoring->getclPz();
      cl2E=scoring->getclE();
      cl2theta=scoring->getclTheta();
      cl2phi=scoring->getclPhi();
      clustM=scoring->getClustM();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      opAngle=scoring->getOpAngleClust();
      prim2lv=scoring->getprimLV();
      //ThreeVector for angular analysis
      prim1vec3.SetXYZ(pr1px,pr1py,pr1pz);
      prim2vec3.SetXYZ(pr2px,pr2py,pr2pz);
      //if(std::cos(prim1vec3.Angle(prim2vec3))>-.5) goto exitFilltree;
      //if(prim2lv.M()<0.04 || prim2lv.M()>.160) goto exitFilltree;
    }else
    if(multiCrys<=3 && singleCrys<=3){
      if(((multiCrys==1 && singleCrys==0)||(multiCrys==0 && singleCrys==1))){
        scoring->init(); 
        goto exitFill;
      }else
      if(multiCrys==3 && singleCrys==0){
        scoring->init();
        std::cout<<"  This is only 3 multi-Crys ---|\n";
        //scoring->clusterEval(clusEne,clusThetaE,clusPhiE);
	scoring->clusterEval(clusEne,clusEr,clusEz,clusThetaE,clusPhiE);
      }else
      if(multiCrys==0 && singleCrys==3){
        scoring->init();
        std::cout<<"  This is only 3 singleCrys ---|\n";
        //scoring->clusterEval(singleEne,singTheta,singPhi);
	scoring->clusterEval(singleEne,singR,singZ,singTheta,singPhi);
      }else{
        scoring->init();
        std::cout<<"  This is only combined ---|\n";
        //scoring->clusterEval(clusEne,singleEne,clusThetaE,clusPhiE,singTheta,singPhi);
	scoring->clusterEval(clusEne,singleEne,clusEr,clusEz,clusThetaE,clusPhiE,singR,singZ,singTheta,singPhi);
      }
      pr2px=scoring->getprPx();   pr2py=scoring->getprPy();  pr2pz=scoring->getprPz();
      E2clust=scoring->getE();
      // cluster PID 1:
      scoring->setCpid(1); //NOTE: must always set this
      cl1px=scoring->getclPx();
      cl1py=scoring->getclPy();
      cl1pz=scoring->getclPz();
      cl1E=scoring->getclE();
      cl1theta=scoring->getclTheta();
      cl1phi=scoring->getclPhi();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      // cluster PID 1:
      scoring->setCpid(2);
      cl2px=scoring->getclPx();
      cl2py=scoring->getclPy();
      cl2pz=scoring->getclPz();
      cl2E=scoring->getclE();
      cl2theta=scoring->getclTheta();
      cl2phi=scoring->getclPhi();
      clustM=scoring->getClustM();
      std::cout<<" **** g1px is => "<<scoring->getclPx()<<"\n";
      opAngle=scoring->getOpAngleClust();
      prim2lv=scoring->getprimLV();
      //ThreeVector for angular analysis
      prim1vec3.SetXYZ(pr1px,pr1py,pr1pz);
      prim2vec3.SetXYZ(pr2px,pr2py,pr2pz);
      //if(std::cos(prim1vec3.Angle(prim2vec3))>-.5) goto exitFilltree;
      // set variables for multiple clusters
      //if(prim2lv.M()<0.09 || prim2lv.M()>.180) goto exitFilltree;
      //if(E2clust<0.100 || E2clust>.300) goto exitFilltree;
    }// End of multiCrys clusters
    /****************************************************
     *     Fill tree Variables
     ****************************************************/
    // K+ Lorentz vector info. from pi+ and pi0
    prim1lv.SetPxPyPzE(pr1px, pr1py, pr1pz,pr1Etot);
    kaon=prim1lv+prim2lv;
    // Fill tree var
    //treeClus->extraTOF1=tgtTree->extraTOF1;
    treeClus->E_prim2=E2clust;
    treeClus->cpid1Px=cl1px;       treeClus->cpid2Px=cl2px;      treeClus->prim2px=prim2lv.Px();
    treeClus->cpid1Py=cl1py;       treeClus->cpid2Py=cl2py;      treeClus->prim2py=prim2lv.Py();
    treeClus->cpid1Pz=cl1pz;       treeClus->cpid2Pz=cl2pz;      treeClus->prim2pz=prim2lv.Pz();
    // position variables:
    treeClus->cpid1x=cl1x;       treeClus->cpid2x=cl2x;
    treeClus->cpid1y=cl1y;       treeClus->cpid2y=cl2y;
    treeClus->cpid1z=cl1z;       treeClus->cpid2z=cl2z;
    treeClus->cpid1r=cl1r;       treeClus->cpid2r=cl2r;
    // pi+ pi0
    treeClus->prim1px=pr1px;
    treeClus->prim1py=pr1py;
    treeClus->prim1pz=pr1pz;
    treeClus->clCosTheta=opAngle;
    treeClus->prCosTheta=std::cos(prim1vec3.Angle(prim2vec3));
    treeClus->M_prim2=prim2lv.M();
    treeClus->prim2M2=prim2lv.M2();
    treeClus->M_k=kaon.M();
    treeClus->kM2=kaon.M2();
    treeClus->cpid1E=cl1E;
    treeClus->cpid2E=cl2E;
    treeClus->cpid1thetaE=cl1theta;
    treeClus->cpid2thetaE=cl2theta;
    treeClus->cpid1phiE=cl1phi;
    treeClus->cpid2phiE=cl2phi;
    treeClus->clusterM=clustM;
    std::cout<<"\n  piPecking total Cluster Energy:  "<<E2clust<<endl;
    std::cout<<"\n  Angular1 checking (centriod)   ("<<cl1theta<<", "<<cl1phi<<")\n";
    std::cout<<"\n  Angular2 checking (centriod)   ("<<cl2theta<<", "<<cl2phi<<")\n";
    std::cout<<"\n  Checking pi0 InvMass:      "<<prim2lv.M()<<endl;
    std::cout<<"\n  Checking cos(theta):       "<<opAngle<<std::endl;
    std::cout<<"\n  Checking vertex opening    "<<std::cos(prim1vec3.Angle(prim2vec3))<<endl;
    std::cout<<"\n  Cluster multiplicity:      "<<clustM<<endl;
    treeClus->channel=(singleCrys+multiCrys);
    exitFill:
    std::cout<<" ........ singleCrys "<<singleCrys<<std::endl;
    std::cout<<" ........ Multi-Crys "<<multiCrys<<std::endl;
    std::cout<<"... End reached! Outta here! \n";
    std::cout<<" ***************************************************************************\n";
  }
  treeClus->isBad=1; // event register, this is the case for good event
  return 0;
}
Long_t Det_CsI::finalize_fit(){
  cout<<"finalize_fit"<<endl;
  return 0;
}
