#include <Det_CsI.h>
#include <fstream>
#include <iostream>

Det_CsI::Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_){
  // Set defaults for various options
  //treeFit=0;
  //setIdCsI(mapCsI);
  std::cout<<" -- Checking the hell outta this thing!\n";
};

Det_CsI::~Det_CsI(){
};
/*
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
*/
Long_t Det_CsI::histos(){
  h1Ch=new TH1D("h1Ch","channel7",17,-.5,16.5);
  return 0;
}

Long_t Det_CsI::startup(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  //getBranchObject("tgtInfo",(TObject **) &tgtTree);
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
  treeFit->ovrX2=dummy;
  treeFit->ovrTime=dummy;
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
  /*
  if(treeRaw->runNo && tgtTree->run != treeRaw->runNo){
    std::cout<<" Oops you are comparing to different runs \n";
    std::cout<<" ***Bailing*** Bailing*** \n";
    std::abort();
  }*/
  // since this is a calibration plugin,
  // we should only consider single crystal hits
  if(treeRaw->nChannel>7) goto exitFill;
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
    // reference timing from 3 modules
    // timing from all 3 modules will considered
    if((treeRaw->indexCsI[iCh]==16 && iFB==0 && iUD==0) &&
                    (indexClock==0 || indexClock==4 || indexClock==8)){
      int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[iCh]];
      SingleCsI myCsI(treeRaw->runNo,myEvent,iModule);
      myCsI.nameCsI(iClock,iFB,iUD,iModule);
      myCsI.setData(treeRaw->data[iCh]);
      myCsI.setIndexTheta(thetaIndex);
      myCsI.setIndexPhi(phiIndex);
      myCsI.initVar();
      myCsI.fit();
      // Three timing modules
      switch(indexClock){
        case 0:
	  std::cout<<"---------| This is reference Module: 1 \n";
          treeFit->rgaus[0]=myCsI.getTime();
          treeFit->tref[0]=myCsI.getCDF50();
	  treeFit->refpk[0]=myCsI.getpTime();
	  break;
	case 4:
	  std::cout<<"---------| This is reference Module: 2 \n";
          treeFit->rgaus[1]=myCsI.getTime();
          treeFit->tref[1]=myCsI.getCDF50();
	  treeFit->refpk[1]=myCsI.getpTime();
	  break;
	case 8:
	  std::cout<<"---------| This is reference Module: 3 \n";
          treeFit->rgaus[2]=myCsI.getTime();
          treeFit->tref[2]=myCsI.getCDF50();
	  treeFit->refpk[2]=myCsI.getpTime();
	  break;
      }// end of switch statement
    }// end of if-timing loop
    //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
    //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
    // start of signal CsI if-loop
    if(!(treeRaw->indexCsI[iCh]==16 && iFB==0 && iUD==0) ||
       !(indexClock==0 || indexClock==2 || indexClock==4 ||
         indexClock==6 || indexClock==8 || indexClock==10)){
      IdCsI myIndex(treeRaw->nameCsI[iCh],treeRaw->indexCsI[iCh]);
      UInt_t iCsI=mapCsI[myIndex];
      if(iClock==7 && iFB==1 && iUD==0){
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
	//h1Ch->Fill(iClock);
      }
      if(iClock==7 && (iFB==1 || iUD==0)){
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
        std::cout<<"============> Hell Ya this is an empty module (?) <==============\n";
	if(iUD==0)
	  h1Ch->Fill((1+iModule));
      }
      //    if(iCsI==144 || iCsI==400||iCsI==656||iCsI==166||iCsI==128) continue;
      //SingleCsI myCsI(treeRaw->runNo,myEvent,iModule);
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
      //Filling the TTree here:
      treeFit->indexCsI=treeRaw->indexCsI[iCh];
      treeFit->clock=indexClock+1;
      treeFit->ud=iUD;
      treeFit->fb=iFB;
      treeFit->tpeak=myCsI.getpTime();
      treeFit->trise=myCsI.getTime();
      if(myCsI.numberWave()==2){
        treeFit->kmu2=myCsI.getphDiff();
        treeFit->dubPed=myCsI.getPedestal();
      }
      if(myCsI.numberWave()==5 || myCsI.numberWave()==6){
        treeFit->ovrpH=myCsI.getphDiff();
        treeFit->ovrped=myCsI.getPedestal();
	treeFit->ovrX2=myCsI.getChi2();
	treeFit->ovrTime=myCsI.getTime();
      }
      treeFit->phei=myCsI.getphDiff();
      treeFit->ped=myCsI.getPedestal();
      treeFit->thSing=myCsI.getTheta();
      treeFit->phiSing=myCsI.getPhi();
      treeFit->chi2=myCsI.getChi2();
      treeFit->ndf=myCsI.getNDF();
      treeFit->waveID=myCsI.numberWave();
      // fill target variables:
      //treeFit->extraTOF1_size=tgtTree->extraTOF1_size;
      //treeFit->extraTOF1=tgtTree->extraTOF1;
      //treeFit->vec_extraTOF1=tgtTree->vec_extraTOF1;
      //treeFit->tof1Gap=tgtTree->TOF1Gap;
      //treeFit->tof2Gap=tgtTree->TOF2Gap;
      std::cout<<"---> rise time: "<<myCsI.getTime()<<"\n";
      std::cout<<"---> WaveID: "<<myCsI.numberWave()<<"\n";
      //std::cout<<"---> extraTOF1: "<<tgtTree->vec_extraTOF1->size()<<"\n";
    }// end of signal CsI if-loop
  }// End of Channel No. for-loop
  exitFill:
  std::cout<<"... End reached! Outta here! \n";
  std::cout<<" ***************************************************************************\n";
  treeFit->isBad=1; // event register, this is the case for good event
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
