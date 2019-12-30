#include <Det_CsI.h>
#include <fstream>
#include <iostream>
using namespace std;

Det_CsI::Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_){
  // Set defaults for various options
  //treeFit=0;
  //setIdCsI(mapCsI);
  std::cout<<" -- Checking the hell outta this thing!\n";
};


Det_CsI::~Det_CsI(){
};

Long_t Det_CsI::setIdCsI(map<IdCsI,UInt_t> & map){
  char fb[2]={'f','b'};
  char ud[2]={'u','d'};
  UInt_t index=1;
  for(UInt_t iHole=0;iHole<12;iHole++){
    for(int ifb=0;ifb<2;ifb++){
      for(int iud=0;iud<2;iud++){
	for(UInt_t iCsI=0;iCsI<16;iCsI++){
	  UInt_t name;
	  name=0x30000000+(((iHole+1)/10)<<24);
	  name+=0x00300000+(((iHole+1)%10)<<16);
	  name+=fb[ifb]<<8;
	  name+=ud[iud];
	  IdCsI first(name,iCsI+1);
	  map[first]=index;
	  index++;
	}
      }
    }
  }
}

Long_t Det_CsI::histos(){

  return 0;
}

Long_t Det_CsI::startup(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  treeFit=new CRTSingleCsI();
  makeBranch("singleCsI",(TObject **) &treeFit);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t Det_CsI::singleAng(int id, int module, int channel, int yy, int zz){
  //std::cout<<"\n --- "<<id<<"-"<<module<<"-"<<channel<<"-"<<yy<<"-"<<zz<<" ---\n";
  phiCsI[module][channel]=yy;
  thetaCsI[module][channel]=zz;
  std::cout<<"\n --- thetaCsI["<<module<<"]["<<channel<<"] = "<<yy<<" <---> "<<thetaCsI[module][channel]<<" ---\n";

  return 0;
}
//Initialize storage variables here
void Det_CsI::initSingleVar(){
  int dummy=-1000;
  treeFit->thSing=dummy;
  treeFit->indexCsI=dummy;      treeFit->phiSing=dummy;
  treeFit->tpeak=dummy; 
  treeFit->trise=dummy;         treeFit->typeAB=dummy;
  treeFit->calInt=dummy;        treeFit->crysID=dummy;
  treeFit->csiArrange[0]=dummy; treeFit->fb=dummy;
  treeFit->csiArrange[1]=dummy; treeFit->ped=dummy;
  treeFit->clock=dummy;
  treeFit->ovrped=dummy;
  treeFit->ovrpH=dummy;
  treeFit->ud=dummy;           
  treeFit->phei=dummy;         
  treeFit->phdstr=dummy;
  //Single pulse                 //double pulse
  treeFit->sphei=dummy;         treeFit->kmu2=dummy;
  treeFit->sptime=dummy;        treeFit->dubPed=dummy;
  treeFit->sped=dummy;          treeFit->dubphei=dummy;
  treeFit->waveID=dummy;        treeFit->intKmu2=dummy;
  for(int i=0;i<3;i++){
    treeFit->tref[i]=dummy;
    treeFit->refpk[i]=dummy;
    treeFit->tcorr[i]=dummy;
    treeFit->refmn[i]=dummy;
    treeFit->rgaus[i]=dummy;
  }
}

Long_t Det_CsI::process(){
  // Save relavant run/event info
  treeFit->eventNo=treeRaw->eventNo;
  treeFit->runNo=treeRaw->runNo;
  if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
  }
  // need to make sure this called for every event
  // avoid memory leaks
  initSingleVar();
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
      SingleCsI myCsI(treeRaw->runNo,myEvent,iModule);
      //std::cout<<" top U/D || F/B :"<<p[0]<<", "<<p[1]<<std::endl;
      //std::cout<<" checking the hell outta this shit:"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[iCh]]<<std::endl;
      int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      myCsI.nameCsI(iClock,iFB,iUD,iModule);
      myCsI.setData(treeRaw->data[iCh]);
      myCsI.setIndexTheta(thetaIndex);
      myCsI.setIndexPhi(phiIndex);
      myCsI.fit();
      // not sure this is needed... just in case
    }// end of signal CsI if-loop
  }
  /****************************************************
   *     Fill tree Variables
   ****************************************************/
  // K+ Lorentz vector info. from pi+ and pi0
  /*
  prim1lv.SetPxPyPzE(pr1px, pr1py, pr1pz,pr1Etot);
  kaon=prim1lv+prim2lv;
  // Fill tree var
  treeClus->extraTOF1=tgtTree->extraTOF1;
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
  treeClus->isBad=1; // event register, this is the case for good event*/
  return 0;
}

Long_t Det_CsI::done(){
  cout<<"***************| Finalize fit |***************\n";
  return 0;
};

Long_t Det_CsI::cmdline(char *cmd){
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_){
    return (Plugin *) new Det_CsI(in_,out_,inf_,outf_,p_);
  }
}

ClassImp(Det_CsI);
