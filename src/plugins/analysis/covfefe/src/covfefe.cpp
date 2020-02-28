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
  h1time[12][2][2][16]=NULL;
  h1cali[12][2][2][16]=NULL;
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
  for(int i=0; i<3; i++){
    std::ostringstream name1, name2,name3,name4,name5;
    name1<<"t_0_"; name2<<"tcorr_"; name3<<"trefCorr_"; name4<<"cdf50_"; name5<<"refgaus_";
    name1<<i;name2<<i;name3<<i;name4<<i;name5<<i;
    h1t_0[i]=new TH1D(name1.str().c_str(),"stats",50,-100,100);
    h1tcorr[i]=new TH1D(name2.str().c_str(),"stats",50,-10.,20.);
    h1trefCorr[i]=new TH1D(name3.str().c_str(),"stats",750,-1000,1000);
    h1cdf50[i]=new TH1D(name4.str().c_str(),"stats",550,-1000.,1000.);
    h1refgaus[i]=new TH1D(name5.str().c_str(),"stats",50,-100,100);
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
  timing=new TH1D("timing","stat",125,-250.,250.);
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
    h1cali[iclock][iFB][iUD][iModule]->Fill(adcVal);
    // Fill Kmu2  and timing histos
    if(treeTime->kmu2>0){
      h1kmu2[iclock][iFB][iUD][iModule]->Fill(treeTime->kmu2);
      h1time[iclock][iFB][iUD][iModule]->Fill(treeTime->rgaus[0]-treeTime->trise);
      for(int n=0; n<3; n++)
        h1cdf50[n]->Fill((treeTime->rgaus[n]-treeTime->trise)*40);
      double tcalc=(treeTime->rgaus[0]-treeTime->trise)*40;
      //if(tcalc>0)
        timing->Fill(tcalc-(rand->Gaus(22.9148,40.2762)));
    }
    // additional timing histos
    if(t_ref1>0 && t_ref2>0 && t_ref3){
      for(int n=0; n<3; n++){
        T_0[n]=treeTime->tref[n]; tcorr[n]=treeTime->tcorr[n];
        refgaus[n]=treeTime->rgaus[n];
        h1t_0[n]->Fill(T_0[n]); 
        h1refgaus[n]->Fill(refgaus[n]);
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
  std::ofstream calib;
  calib.open("newcalibPar.txt");
  double xmax, xx, calpar;
  TCanvas* c0=new TCanvas("c0"," Peak time ",808,700);
  c0->cd();
  timing->SetTitle("CsI Timing Study");
  timing->GetXaxis()->SetTitle("time [ns]");
  timing->GetXaxis()->SetTitleSize(0.045);
  timing->GetYaxis()->SetTitle("counts/bin");
  timing->SetLineWidth(2);
  timing->Draw("hist");
  c0->Update();
  c0->Write();
  /*
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          xmax=h1time[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-2.5;
          upRange=xx+2.5;
	  TF1* f1=new TF1("f1","gaus",lowRange,upRange);
	  h1time[iClock][iFB][iUD][iModule]->Fit(f1,"QR");
	  t_corr=f1->GetMaximumX();
	  calib<<iClock<<"\t"<<iFB<<"\t"<<iUD<<"\t"<<iModule<<"\t"<<t_corr<<"\n";
	  for(int n=0; n<nbins;n++){
	    double yy=h1time[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x=h1time[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    if(!(iClock==7 && iFB==1 && iUD==0)){
              double xnew=(t_corr-x);//*40;
              h1Tcorr->Fill(xnew); //,yy);
	    }
	    if(!(iClock==6 && iFB==1 && (iUD==0 || iUD==1))){
              double xnew=(t_corr-x);//*40;
              h1Tcorr->Fill(xnew); //,yy);
	    }
	  }
          xmax =h1kmu2[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx   =h1kmu2[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1kmu2[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-85;
          upRange=xx+77;
	  TF1* f2=new TF1("f2","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1kmu2[iClock][iFB][iUD][iModule]->Fit(f2,"QR");
	  apcsi=f2->GetMaximumX();
	  calpar=dE/apcsi;
	  for(int n=0; n<h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();n++){
	    double yy=h1kmu2[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x =h1kmu2[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    integHist->Fill(xnew,yy);
	    Ecorr->Fill(xnew+8.9,yy);
	  }
	  delete f1;
	  delete f2;
        }
      }
    }
  }
  calib.close();
  gStyle->SetOptStat(0);
  TCanvas* c1=new TCanvas("MarinateCsI","Pulse-height distribution",3508,2480);
  TCanvas* c2=new TCanvas("MarinateCsI2","Integrated pulse-height distribution",3508,2480);
  TCanvas* c3=new TCanvas("E_CsI","Energy CsI",808,700);
  TCanvas* c4=new TCanvas("ECsI","Energy CsI comparison",808,700);
  TCanvas* c5=new TCanvas("Phei","Pulse height Distribution",808,700);
  TCanvas* c6=new TCanvas("c6"," CsI timing",808,700);
  TCanvas* c7=new TCanvas("c7"," Peak time ",808,700);
  c1->Divide(3,4);
  c2->Divide(3,4);
  for(int iClock=0;iClock<1;iClock++){
    for(int iFB=0;iFB<1;iFB++){
      for(int iUD=0;iUD<1;iUD++){
        for(int iModule=0;iModule<12;iModule++){
          c1->cd(iModule+1);
          h1time[iClock][iFB][iUD][iModule]->Draw();
          c2->cd(iModule+1);
          h1kmu2[iClock][iFB][iUD][iModule]->Draw();
	}
      }
    }
  }
  c3->cd();
  Ecorr->SetTitle("CsI: reconstructed energy for K_{#mu2} ");
  Ecorr->GetXaxis()->SetTitle("Energy (T_{#mu}) [MeV]");
  Ecorr->GetYaxis()->SetTitle("counts/bin");
  //Ecorr->GetYaxis()->SetRangeUser(0,1050e3);
  Ecorr->SetLineColor(kCyan);
  Ecorr->SetLineWidth(2);
  Ecorr->Draw("hist");
  integHist->Draw("hist same");
  auto leg0=new TLegend(0.1,0.7,0.48,0.9);
  leg0->SetHeader("Key:","C");
  leg0->AddEntry(Ecorr, "CD E_{loss} correction (#mu=153.67, #sigma=9.78)");
  leg0->AddEntry(intEn, "dE (#mu=145.8, #sigma=12.9)");
  leg0->Draw();
  c4->cd();
  Ecorr->SetTitle("Comparing CsI Energy for K_{#mu2} for 2 methods ");
  Ecorr->GetXaxis()->SetTitle("Deposited energy (T_{#mu}) [MeV]");
  Ecorr->GetYaxis()->SetTitle("counts");
  Ecorr->SetLineWidth(2);
  Ecorr->Draw("hist");
  intEn->SetLineColor(kCyan+3);
  intEn->SetLineWidth(3);
  intEn->Draw("hist same");
  auto leg=new TLegend(0.1,0.7,0.48,0.9);
  leg->SetHeader("Key:","C");
  leg->AddEntry(Ecorr, "Pulse-height: E_{loss} applied (#mu=155.3, #sigma=11.7)");
  leg->AddEntry(intEn, "Integrated waveform: E_{loss} applied (#mu=153.3, #sigma=19.5)");
  leg->Draw();
  c5->cd();
  phdis->SetTitle("Pulse-height distribution");
  phdis->GetXaxis()->SetTitle("pulse-height");
  phdis->GetYaxis()->SetTitle("counts/bin");
  phdis->SetLineWidth(2);
  phdis->Draw("hist");
  hkmu2->SetFillColor(kMagenta);
  hkmu2->SetFillStyle(3244);
  hkmu2->SetLineWidth(2);
  hkmu2->Draw("hist same");
  auto leg2=new TLegend(0.1,0.7,0.48,0.9);
  leg2->SetHeader("Key:","C");
  leg2->AddEntry(phdis, "Pulse-height distribution");
  leg2->AddEntry(hkmu2, "K_{#mu2} pulse-height distribtion");
  leg2->Draw();
  c6->cd();
  h1Tcorr->SetTitle("CsI timing");
  h1Tcorr->GetXaxis()->SetTitle("time");
  h1Tcorr->GetYaxis()->SetTitle("counts/bin");
  h1Tcorr->SetLineWidth(2);
  h1Tcorr->Draw("hist");
  c7->cd();
  phdistr->SetTitle("Peak time distribution");
  phdistr->GetXaxis()->SetTitle("time [ch]");
  phdistr->GetYaxis()->SetTitle("counts/bin");
  phdistr->SetLineWidth(2);
  phdistr->Draw("hist");

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
