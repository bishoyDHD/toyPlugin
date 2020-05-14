#include <covfefe.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TVector2.h>
#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>

covfefe::covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
  treeTime=NULL;
  calibcsi=0;
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          h1time[iClock][iFB][iUD][iModule]=NULL;
          h1cali[iClock][iFB][iUD][iModule]=NULL;
          h1kmu2[iClock][iFB][iUD][iModule]=NULL;
        }
      }
    }
  }
  calibHist=NULL;
  Ecorr=NULL;
  integHist=NULL;
  intEn=NULL;
  phdis=NULL;
  hkmu2=NULL;
  rand=new TRandom();
  std::cout<<" checking this shit \n";
};

covfefe::~covfefe(){
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<15;iModule++){
          delete h1time[iClock][iFB][iUD][iModule];
          delete h1cali[iClock][iFB][iUD][iModule];
        }
      }
    }
  }
  for(int i=0; i<3; i++){
    delete h1t_0[i];
    delete h1tcorr[i];
    delete h1trefCorr[i];
    delete h1cdf50[i];
    delete h1refgaus[i];
  }
  delete calibHist;
  delete Ecorr;
  delete integHist;
  delete intEn;
  delete phdis;
  delete hkmu2;
  delete timing;
};

Long_t covfefe::histos(){
  std::ostringstream cfdvE,rtvE,cfdvpt;
  cfdvE<<"cfdVE";rtvE<<"rtVE";cfdvpt<<"cfdVpt";
  h1sigRT=new TH1D("h1sigRT","stats",244,-0.5,60.5);
  h1sigPT=new TH1D("h1sigPT","stats",44,57.5,67.5);
  h1sigCDF50=new TH1D("h1sigCDF50","stats",84,30.5,50.5);
  h2cfd50Vph=new TH2D(cfdvE.str().c_str(),"stats",100,-.5,999.5,180,5.0,50.0);
  h2rTvph=new TH2D(rtvE.str().c_str(),"stats",100,-.5,999.5,244,-.5,60.5);
  h2cfd50Vpt=new TH2D(cfdvpt.str().c_str(),"stats",20,60.0,65.0,180,5.0,50.0);
  for(int i=0; i<3; i++){
    std::ostringstream rtime,ptime,cdf50time,csiRef,refcd50;
    std::ostringstream rtCorr,ptCorr,cdf50Corr;
    std::ostringstream diffrtime,diffptime,diffcdf50;
    std::ostringstream cfdvEcorr,rtvEcorr,cfdvptcorr;
    rtime<<"riseT_";ptime<<"pkTime_";cdf50time<<"cdf50Time_";csiRef<<"csiRefpT_";
    rtCorr<<"refRTcorr_";ptCorr<<"refpkTcorr_";cdf50Corr<<"refcdf50corr_";refcd50<<"cd50refCsI_";
    diffrtime<<"refDiffRT_";diffptime<<"refDiffpkT_";diffcdf50<<"refDiffcdf50_";
    cfdvEcorr<<"cfdVEcorr_";rtvEcorr<<"rtVEcorr_";cfdvptcorr<<"cfdVptcorr_";
    rtime<<i;ptime<<i;cdf50time<<i;
    rtCorr<<i;ptCorr<<i;cdf50Corr<<i;csiRef<<i;
    diffrtime<<i;diffptime<<i;diffcdf50<<i;refcd50<<i;
    cfdvEcorr<<i;rtvEcorr<<i;cfdvptcorr<<i;
    h1pkCsI[i]=new TH1D(csiRef.str().c_str(),"stats",90,59.95,65.05);
    h1rTime[i]=new TH1D(rtime.str().c_str(),"stats",150,-900,300);
    h1pkTime[i]=new TH1D(ptime.str().c_str(),"stats",20,-100.,100.);
    h1cdf50[i]=new TH1D(cdf50time.str().c_str(),"stats",50,-100.5,99.5);
    h1rTimeCorr[i]=new TH1D(rtCorr.str().c_str(),"stats",183,-.5,60.5);
    h1pTimeCorr[i]=new TH1D(ptCorr.str().c_str(),"stats",30,59.5,65.5);
    h1cdf50Corr[i]=new TH1D(cdf50Corr.str().c_str(),"stats",93,29.5,60.5);
    h1refrTimeDiff[i]=new TH1D(diffrtime.str().c_str(),"stats",180,-45,45.);
    h1refpTimeDiff[i]=new TH1D(diffptime.str().c_str(),"stats",30,-15,15);
    h1refcdf50Diff[i]=new TH1D(diffcdf50.str().c_str(),"stats",180,-45,45);
    h1refCDF50[i]=new TH1D(refcd50.str().c_str(),"stats",44,39.5,50.5);
    h2cfd50VphCorr[i]=new TH2D(cfdvEcorr.str().c_str(),"stats",100,-.5,999.5,84,30.5,50.5);
    h2rTvphCorr[i]=new TH2D(rtvEcorr.str().c_str(),"stats",100,-.5,999.5,244,-.5,60.5);
    h2cfd50VptCorr[i]=new TH2D(cfdvptcorr.str().c_str(),"stats",44,57.5,67.5,84,30.5,50.5);
  }
  //Other timing histograms--> will be deprecated soon
  for(int i=0; i<3; i++){
    std::ostringstream name1,name2,name3,name4,name5;
    name1<<"t_0_"; name2<<"tcorr_"; name3<<"trefCorr_"; name4<<"cdf50_"; name5<<"refgaus_";
    name1<<i;name2<<i;name3<<i;name4<<i;name5<<i;
    h1t_0[i]=new TH1D(name1.str().c_str(),"stats",50,-100,100);
    h1tcorr[i]=new TH1D(name2.str().c_str(),"stats",50,-10.,20.);
    h1trefCorr[i]=new TH1D(name3.str().c_str(),"stats",750,-1000,1000);
    h1refgaus[i]=new TH1D(name5.str().c_str(),"stats",44,34.5,45.5);
  }
  std::ostringstream nameCal, ename, nameInt, inEne;
  nameCal<<"CalibCsI";
  ename<<"E_corr";
  nameInt<<"integCsI";
  inEne<<"IntEnergy";
  calibHist=new TH1D(nameCal.str().c_str(),"stat",62.0,0,250);
  Ecorr=new TH1D(ename.str().c_str(),"stat",63.0,0,250);
  integHist=new TH1D(nameInt.str().c_str(),"stat",63.0,0,250);
  intEn=new TH1D(inEne.str().c_str(),"stat",63.0,0,250);
  phdis=new TH1D("pheight","stat",150.,0,1000);
  timing=new TH1D("timing","stat",75.,-200.,200.);
  h1Tcorr=new TH1D("tCorr","stat",40.,-100.,100.);
  phdistr=new TH1D("ptdist","stat",62.5,1,249);
  hkmu2=new TH1D("kmu2Ds","stat",150.,0,1000);
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          std::ostringstream name, name2, name3, name4, name5, kname, tname,name6;
          name<<"IntEne_"; name2<<"Mnfit_"; name3<<"F1fit_"; kname<<"kmu2Calib_"; name5<<"pHeight", name6<<"Diff";
          tname<<"time";
          name<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          tname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          kname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1time[iClock][iFB][iUD][iModule]=new TH1D(tname.str().c_str(),"stat",90.,-10.,20.);
          h1cali[iClock][iFB][iUD][iModule]=new TH1D(name.str().c_str(),"stat",187.5,0,1200);
          h1kmu2[iClock][iFB][iUD][iModule]=new TH1D(kname.str().c_str(),"stat",187.5,0,1200);
        }
      }
    }
  }
  return 0;
}

Long_t covfefe::startup(){
  getBranchObject("singleCsI",(TObject **) &treeTime); 
  calibcsi=new CATCaliCsI();
  makeBranch("marinCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::process(){
  iclock=treeTime->clock-1;
  //std::cout<<"-----> "<<treeTime->extraTOF1->size()<<std::endl;
  iModule=treeTime->indexCsI-1;
  iUD=treeTime->ud; iFB=treeTime->fb;
  adcVal=treeTime->phei;
  //std::cout<<" ... so called adc val: "<<adcVal<<std::endl;
  intVal=treeTime->intKmu2;
  phdis->Fill(treeTime->phei);
  hkmu2->Fill(treeTime->kmu2);
  // timing variables
  t_ref1=treeTime->rgaus[0];
  t_ref2=treeTime->rgaus[1];
  t_ref3=treeTime->rgaus[2];
  if(adcVal>0){
    std::cout<<" -----------------------------------" <<std::endl;
    std::cout<<" value of clock:  "<<iclock<<std::endl;
    std::cout<<" value of Module: "<<iModule<<std::endl;
    std::cout<<" value of iUD:    "<<iUD<<std::endl;
    std::cout<<" value of iFB:    "<<iFB<<std::endl;
    std::cout<<" value of adcVal: "<<adcVal<<std::endl;
    std::cout<<" Kmu2 adc value1: "<<treeTime->kmu2<<std::endl;
    std::cout<<" Kmu2 adc value2: "<<treeTime->ovrpH<<std::endl;
    std::cout<<" timing of Cryst: "<<treeTime->trise<<std::endl;
    std::cout<<" timing of Ref 1: "<<treeTime->rgaus[0]<<std::endl;
    std::cout<<" timing of Ref 2: "<<treeTime->rgaus[1]<<std::endl;
    std::cout<<" timing of Ref 3: "<<treeTime->rgaus[2]<<std::endl;
    csiphi=treeTime->phiSing;
    // reference Crystal times
    refrT[0]=treeTime->rgaus[0];
    refrT[1]=treeTime->rgaus[1];
    refrT[2]=treeTime->rgaus[2];
    refpkT[0]=treeTime->refpk[0];
    refpkT[1]=treeTime->refpk[1];
    refpkT[2]=treeTime->refpk[2];
    refcdf50[0]=treeTime->tref[0];
    refcdf50[1]=treeTime->tref[1];
    refcdf50[2]=treeTime->tref[2];
    h2cfd50Vph->Fill(adcVal,treeTime->cdf50);
    h2rTvph->Fill(adcVal,treeTime->trise);
    h2cfd50Vpt->Fill(treeTime->tpeak,treeTime->cdf50);
    for(int n=0; n<3; n++){
      h1pkCsI[n]->Fill(refpkT[n]);
    }
    h1cali[iclock][iFB][iUD][iModule]->Fill(adcVal);
    // Fill Kmu2  and timing histos
    if(treeTime->kmu2>0){
      h1sigCDF50->Fill(treeTime->cdf50);
      h1sigPT->Fill(treeTime->tpeak);
      h1sigRT->Fill(treeTime->trise);
      timeKmu2[0]=std::abs(treeTime->tpeak-refpkT[0]);
      timeKmu2[1]=std::abs(treeTime->tpeak-refpkT[1]);
      timeKmu2[2]=std::abs(treeTime->tpeak-refpkT[2]);
      h1kmu2[iclock][iFB][iUD][iModule]->Fill(treeTime->kmu2);
      h1time[iclock][iFB][iUD][iModule]->Fill(treeTime->rgaus[0]-treeTime->trise);
    }
    double tcalc=(treeTime->rgaus[0]-treeTime->cdf50)*40;
    if(treeTime->kmu2>500 && treeTime->kmu2<800){
      if(treeTime->tpeak>62.5 && treeTime->tpeak<64.5){
        double tpeakCorr=(refcdf50[0]-treeTime->cdf50)*40-1.56785e+02;
        timing->Fill(tcalc+2.102287e+02);
        if(tpeakCorr>=-43.322 && tpeakCorr<=43.322){
          h1pTimeCorr[0]->Fill(treeTime->tpeak);
          h1cdf50Corr[0]->Fill(treeTime->cdf50);
          h1rTimeCorr[0]->Fill(treeTime->trise);
        }
        // cdf50
        h1refcdf50Diff[0]->Fill((refcdf50[0]-refcdf50[1])*40);
        h1refcdf50Diff[1]->Fill((refcdf50[0]-refcdf50[2])*40);
        h1refcdf50Diff[2]->Fill((refcdf50[1]-refcdf50[2])*40);
        // riseTime
        h1refrTimeDiff[0]->Fill((refrT[0]-refrT[1])*40);
        h1refrTimeDiff[1]->Fill((refrT[0]-refrT[2])*40);
        h1refrTimeDiff[2]->Fill((refrT[1]-refrT[2])*40);
        // peakTime
        h1refpTimeDiff[0]->Fill((refpkT[0]-refpkT[1])*40);
        h1refpTimeDiff[1]->Fill((refpkT[0]-refpkT[2])*40);
        h1refpTimeDiff[2]->Fill((refpkT[1]-refpkT[2])*40);
        for(int n=0; n<3; n++){
          h1cdf50[n]->Fill((refcdf50[n]-treeTime->cdf50)*40-1.56785e+02);
          h1pkTime[n]->Fill((refpkT[n]-treeTime->tpeak)*40);
          h1rTime[n]->Fill((refrT[n]-treeTime->trise)*40-2.46849e+02);
          h1refCDF50[n]->Fill(refcdf50[n]);
          h1refgaus[n]->Fill(refgaus[n]);
          h2cfd50VphCorr[n]->Fill(treeTime->kmu2,(refcdf50[n]-treeTime->cdf50)*40-1.56785e+02);
          h2rTvphCorr[n]->Fill(treeTime->kmu2,(refrT[n]-treeTime->trise)*40-2.46849e+02);
          h2cfd50VptCorr[n]->Fill(refpkT[n],(refcdf50[n]-treeTime->cdf50)*40-1.56785e+02);
        }
      }
    }
    // additional timing histos
    if(t_ref1>0 && t_ref2>0 && t_ref3){
      for(int n=0; n<3; n++){
        T_0[n]=treeTime->tref[n]; tcorr[n]=treeTime->tcorr[n];
        refgaus[n]=treeTime->rgaus[n];
        h1t_0[n]->Fill(T_0[n]); 
	tCalc=refgaus[n]-treeTime->trise;
        h1trefCorr[n]->Fill((refgaus[n]-treeTime->trise)*40.);
	h1tcorr[n]->Fill(refgaus[n]-treeTime->trise); 
      }
    }
  }
/*
  if(t_ref1> -900 || treeTime->rgaus[1]>0 || treeTime->rgaus[2]>0){
    double t_rise=treeTime->rgaus[0];
    std::cout<<"\n timing check 1. ------> "<<t_rise<<std::endl;
    std::cout<<"\n timing check 2. ------> "<<treeTime->rgaus[1]<<std::endl;
    std::cout<<"\n timing check 3. ------> "<<treeTime->rgaus[2]<<std::endl;
    //timing->Fill(treeTime->tcsi);
  }
*/
  double pktime=treeTime->phdstr;
  if(pktime>0 && pktime<250)
    phdistr->Fill(treeTime->phdstr);

  return 0; // 0 = all ok
};

Long_t covfefe::finalize(){
  TSpectrum *s;
  gStyle->SetOptStat(0);
  //std::ofstream calib;
  //calib.open("newcalibPar.txt");
  double xmax, xx, calpar;
  TCanvas* c0=new TCanvas("c0"," Peak time ",808,700);
  c0->cd();
  timing->SetTitle("CsI Timing Study");
  timing->GetXaxis()->SetTitle("time [ns]");
  timing->GetXaxis()->SetTitleSize(0.045);
  timing->GetYaxis()->SetTitle("counts/bin");
  timing->SetLineWidth(2);
  //timing->Draw("hist");
  //timing->Write("hist");
  c0->Write();
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          //h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->SetRangeUser(400,900);
          xmax =h1kmu2[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx   =h1kmu2[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1kmu2[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
          //h1kmu2[iClock][iFB][iUD][iModule]->GetXaxis()->SetRangeUser(0.0,1200);
	  lowRange=xx-90;
          upRange=xx+70;
	  TF1* f1=new TF1("f1","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1kmu2[iClock][iFB][iUD][iModule]->Fit(f1,"QR");
	  apcsi=f1->GetMaximumX();
	  calpar=dE/apcsi;
          if(apcsi<300. || calpar<0.){
	    calpar=0.232;
          }
	  //calib<<iClock<<"\t"<<iFB<<"\t"<<iUD<<"\t"<<iModule<<"\t"<<calpar<<"\n";
	  for(int n=0; n<h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();n++){
	    double yy=h1kmu2[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x =h1kmu2[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    integHist->Fill(xnew,yy);
	    Ecorr->Fill(xnew+8.9,yy);
	  }
	  delete f1;
        }
      }
    }
  }
  //calib.close();
  gStyle->SetOptStat(0);
  //TCanvas* c1=new TCanvas("MarinateCsI","Pulse-height distribution",3508,2480);
  //TCanvas* c2=new TCanvas("MarinateCsI2","Integrated pulse-height distribution",3508,2480);
  //TCanvas* c3=new TCanvas("E_CsI","Energy CsI",808,700);
  //TCanvas* c4=new TCanvas("ECsI","Energy CsI comparison",808,700);
  //TCanvas* c5=new TCanvas("Phei","Pulse height Distribution",808,700);
  //TCanvas* c6=new TCanvas("c6"," CsI timing",808,700);
  //TCanvas* c7=new TCanvas("c7"," Peak time ",808,700);
  //c1->Divide(3,4);
  //c2->Divide(3,4);
  for(int iClock=0;iClock<1;iClock++){
    for(int iFB=0;iFB<1;iFB++){
      for(int iUD=0;iUD<1;iUD++){
        for(int iModule=0;iModule<12;iModule++){
          //c1->cd(iModule+1);
          h1time[iClock][iFB][iUD][iModule]->Write();
          //c2->cd(iModule+1);
          h1kmu2[iClock][iFB][iUD][iModule]->Write();
	}
      }
    }
  }
/*
  c1->Write();
  c2->Write();
  c3->cd();
  Ecorr->SetTitle("CsI: reconstructed energy for K_{#mu2} ");
  Ecorr->GetXaxis()->SetTitle("Energy (T_{#mu}) [MeV]");
  Ecorr->GetYaxis()->SetTitle("counts/bin");
  //Ecorr->GetYaxis()->SetRangeUser(0,1050e3);
  Ecorr->SetLineColor(kCyan);
  Ecorr->SetLineWidth(2);
  Ecorr->Write("hist");
  integHist->Write("hist same");
  auto leg0=new TLegend(0.1,0.7,0.48,0.9);
  leg0->SetHeader("Key:","C");
  leg0->AddEntry(Ecorr, "CD E_{loss} correction (#mu=153.67, #sigma=9.78)");
  leg0->AddEntry(intEn, "dE (#mu=145.8, #sigma=12.9)");
  leg0->Write();
  c3->Write();
  c4->cd();
  Ecorr->SetTitle("Comparing CsI Energy for K_{#mu2} for 2 methods ");
  Ecorr->GetXaxis()->SetTitle("Deposited energy (T_{#mu}) [MeV]");
  Ecorr->GetYaxis()->SetTitle("counts");
  Ecorr->SetLineWidth(2);
  Ecorr->Draw("hist");
  intEn->SetLineColor(kCyan+3);
  intEn->SetLineWidth(3);
  intEn->Write("hist same");
  intEn->Write("hist same");
  auto leg=new TLegend(0.1,0.7,0.48,0.9);
  leg->SetHeader("Key:","C");
  leg->AddEntry(Ecorr, "Pulse-height: E_{loss} applied (#mu=155.3, #sigma=11.7)");
  leg->AddEntry(intEn, "Integrated waveform: E_{loss} applied (#mu=153.3, #sigma=19.5)");
  leg->Write();
  c4->Write();
  c5->cd();
  phdis->SetTitle("Pulse-height distribution");
  phdis->GetXaxis()->SetTitle("pulse-height");
  phdis->GetYaxis()->SetTitle("counts/bin");
  phdis->SetLineWidth(2);
  phdis->Write("hist");
  hkmu2->SetFillColor(kMagenta);
  hkmu2->SetFillStyle(3244);
  hkmu2->SetLineWidth(2);
  hkmu2->Write("hist same");
  auto leg2=new TLegend(0.1,0.7,0.48,0.9);
  leg2->SetHeader("Key:","C");
  leg2->AddEntry(phdis, "Pulse-height distribution");
  leg2->AddEntry(hkmu2, "K_{#mu2} pulse-height distribtion");
  leg2->Write();
  c5->Write();
  c6->cd();
  h1Tcorr->SetTitle("CsI timing");
  h1Tcorr->GetXaxis()->SetTitle("time");
  h1Tcorr->GetYaxis()->SetTitle("counts/bin");
  h1Tcorr->SetLineWidth(2);
  h1Tcorr->Write("hist");
  c6->Write();
  c7->cd();
  phdistr->SetTitle("Peak time distribution");
  phdistr->GetXaxis()->SetTitle("time [ch]");
  phdistr->GetYaxis()->SetTitle("counts/bin");
  phdistr->SetLineWidth(2);
  phdistr->Write("hist");
  c7->Write();
*/
  return 0; // 0 = all ok
};

Long_t covfefe::cmdline(char* cmd){

  return 0;
}

extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new covfefe(in,out,inf_,outf_,p);
  }
}
