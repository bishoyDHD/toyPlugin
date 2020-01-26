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

Long_t covfefe::histosCsI(){
  tof1ang=new TH1D("tof1Angles","stat",24,0,360.);
  h1angDiff=new TH1D("angDiff","stat",48.,-400.5,360.);
  angScore=new TH1D("angScore","stat",48.,-400.5,360.);
  h1score=new TH1D("h1score","stat",16.,-0.5,15.5);
  h1tof1=new TH1D("tofhits","stat",13,-0.5,12.5);
  tof1ID=new TH1D("tofID","stat",13,-0.5,12.5);
  csiAng=new TH2D("csiAngles","stat",24,0.,360.,24,0.,360.);

  return 0;
}

Long_t covfefe::startupCsI(){
  getBranchObject("singleCsI",(TObject **) &treeCalib); 
  calibcsi=new CATCaliCsI();
  makeBranch("sanityCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::processCsI(){
  iclock=treeCalib->clock-1;
  //std::cout<<"-----> "<<treeCalib->extraTOF1->size()<<std::endl;
  iModule=treeCalib->indexCsI-1;
  iUD=treeCalib->ud; iFB=treeCalib->fb;
  adcVal=treeCalib->phei;
/*
  std::cout<<" ... so called adc val: "<<adcVal<<std::endl;
  intVal=treeCalib->intKmu2;
  // timing variables
  //double t_ref=treeCalib->tref;
  if(adcVal > 0){
    std::cout<<" -----------------------------------" <<std::endl;
    std::cout<<" value of clock:  "<<iclock<<std::endl;
    std::cout<<" value of Module: "<<iModule<<std::endl;
    std::cout<<" value of iUD:    "<<iUD<<std::endl;
    std::cout<<" value of iFB:    "<<iFB<<std::endl;
    std::cout<<" value of adcVal: "<<adcVal<<std::endl;
    csiphi=treeCalib->phiSing;
    // timing histos
    if(treeCalib->extraTOF1->size()>0){
      for(UInt_t i=0; i <treeCalib->extraTOF1->size(); i++){
	index=treeCalib->extraTOF1->at(i);
        tof1ang->Fill(phi[index-1]);
        h1angDiff->Fill(phi[index-1]-csiphi);
	csiAng->Fill(csiphi,phi[index-1]);
	h1score->Fill((*(treeCalib->vec_extraTOF1))[i][1]);
	h1tof1->Fill(treeCalib->vec_extraTOF1->size());
	tof1ID->Fill(index);
	if((*(treeCalib->vec_extraTOF1))[i][1]<3){
          angScore->Fill(phi[index-1]-csiphi);
	}
	//std::cout<<".... Checking name of TOF1 & score: "<<treeCalib->extraTOF1->at(i);
	//std::cout<<" "<<(*(treeCalib->vec_extraTOF1))[i][1]<<"\n";
      }
    }
  }
*/
  return 0; // 0 = all ok
};

Long_t covfefe::finalizeCsI(){
/*
  TSpectrum *s;
  gStyle->SetOptStat(0);
  std::ofstream calib;
  //calib.open("calibPar.txt");
  double xmax, xx, calpar;
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<15;iModule++){
          xmax=h1time[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-110;
          upRange=xx+100;
	  TF1* f1=new TF1("f1","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1time[iClock][iFB][iUD][iModule]->Fit(f1,"QR");
	  apcsi=f1->GetMaximumX();
	  calpar=dE/apcsi;
	  //calib<<iClock<<"\t"<<iFB<<"\t"<<iUD<<"\t"<<iModule<<"\t"<<calpar<<"\n";
	  for(int n=0; n<nbins;n++){
	    double yy=h1time[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x=h1time[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    calibHist->Fill(xnew,yy);
	    Ecorr->Fill(xnew+8.9,yy);
	  }
          xmax=h1cali[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx=h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-11000;
          upRange=xx+10000;
	  TF1* f2=new TF1("f2","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1cali[iClock][iFB][iUD][iModule]->Fit(f2,"QR");
	  apcsi=f2->GetMaximumX();
	  calpar=dE/apcsi;
	  for(int n=0; n<h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();n++){
	    double yy=h1cali[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x=h1cali[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    integHist->Fill(xnew,yy);
	    intEn->Fill(xnew+8.9,yy);
	  }
	  delete f1;
	  delete f2;
        }
      }
    }
  }
  //calib.close();
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
          h1cali[iClock][iFB][iUD][iModule]->Draw();
	}
      }
    }
  }
  c1->Write();
  c2->Write();
  c3->cd();
  Ecorr->SetTitle("CsI: reconstructed energy for K_{#mu2} ");
  Ecorr->GetXaxis()->SetTitle("Energy (T_{#mu}) [MeV]");
  Ecorr->GetYaxis()->SetTitle("counts/bin");
  Ecorr->GetYaxis()->SetRangeUser(0,1050e3);
  Ecorr->SetLineColor(kCyan);
  Ecorr->SetLineWidth(2);
  Ecorr->Draw("hist");
  calibHist->Draw("hist same");
  auto leg0=new TLegend(0.1,0.7,0.48,0.9);
  leg0->SetHeader("Key:","C");
  leg0->AddEntry(Ecorr, "CD E_{loss} correction (#mu=153.67, #sigma=9.78)");
  leg0->AddEntry(intEn, "dE (#mu=145.8, #sigma=12.9)");
  leg0->Draw();
  c3->Write();
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
  c4->Write();
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
  c5->Write();
  c6->cd();
  timing->SetTitle("CsI timing");
  timing->GetXaxis()->SetTitle("time");
  timing->GetYaxis()->SetTitle("counts/bin");
  timing->SetLineWidth(2);
  timing->Draw("hist");
  c6->Write();
  c7->cd();
  phdistr->SetTitle("Peak time distribution");
  phdistr->GetXaxis()->SetTitle("time [ch]");
  phdistr->GetYaxis()->SetTitle("counts/bin");
  phdistr->SetLineWidth(2);
  phdistr->Draw("hist");
  c7->Write();
*/
  return 0; // 0 = all ok
};
