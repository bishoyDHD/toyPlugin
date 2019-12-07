#include <singleCsI.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TList.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TFrame.h>
#include <TSpectrum.h>
using std::exp;
void ClusterCsI::setData(const vector<UShort_t>& listD){
  for(unsigned int iS=0,nS=listD.size();iS<nS;iS++){
    mListData.push_back(listD[iS]);
  }
}
char* ClusterCsI::nameCsI(const UInt_t& iClock, const UInt_t& iFB,const UInt_t& iUD,const UInt_t& iModule){
  char ifb[3]="fb";
  char iud[3]="ud";
  sprintf(mName,"%i%c%c%i",iClock,ifb[iFB],iud[iUD],iModule);
  //std::cout<<" checking the name: "<<mName<<std::endl;
  return mName;
}
char* ClusterCsI::nameCsI(unsigned int index){
  if(index<1){
    sprintf(mName,"0000");
  }
  unsigned int iClock=(index-1)/(2*2*16)+1;
  unsigned int iInside=index-1-(iClock-1)*2*2*16;
  unsigned int iFB=iInside/(2*16);
  unsigned int iUD=(iInside-iFB*2*16)/16;
  unsigned int iPaddle=iInside%16;
  char NameFB[3]="fb";
  char NameUD[3]="ud";
  sprintf(mName,"%i%c%c%i",iClock,NameFB[iFB],NameUD[iUD],iPaddle+1);
  return mName;
}
Double_t ClusterCsI::round(double var){
  // type cast to int in order to obtain integer var
  // then divide by 100 or the value is set to correct decimal points
  // Useful upto 1e4 higher than this is basically a no-no
  double value=(int)(var*100+.5);
  return (double)value/100;
}
bool ClusterCsI::fit(){
  //std::cout<<"|| entering into fitting function..... success! \n";
  static unsigned int count=0;
  char name[256];
  //cout<<"processing "<<mIndexCsI<<" "<<nameCsI(mIndexCsI)<<endl;
  cout<<"processing "<<mIndexCsI<<" "<<mName<<endl;
  sprintf(name,"wave_run%d_%dCsI_%s_%u",mRunNo,mEventNo,mName,count);
  std::cout<<"  name of histogram... checking here: "<<name<<std::endl;
  //sprintf(name,"wave_run%d_%dCsI_%s_%u",mRunNo,mEventNo,nameCsI(mIndexCsI),count);
  char title[256];
  sprintf(title,"fAdc waveform");
  shared_ptr<TH1D> h1(new TH1D(name,title,250,0.5,250.5));
  //TH1D* h1=new TH1D(name,title,250,0.5,250.5);
  for(UInt_t iData=0,nData=mListData.size();iData<nData;iData++){
    h1->SetBinContent(iData+1,mListData[iData]);
    h1->SetBinError(iData+1,1.1);
    //h1->SetBinContent(iData,mListData[iData]);
  }
  Double_t xmax=h1->GetBinLowEdge(h1->GetMaximumBin());
  if(xmax<=57 || xmax>=68) return false;
  Double_t xmin=h1->GetBinLowEdge(h1->GetMinimumBin());
  Double_t ymax=h1->GetBinContent(h1->FindBin(xmax));
  Double_t yped=h1->GetBinContent(h1->FindBin(xmin));
  // check for type of waveform (single, double, triple etc)
  unique_ptr<TSpectrum> s(new TSpectrum(4));
  Int_t nfound=s->Search(h1.get(),2,"",0.10);
  //Int_t nfound=s->Search(h1,2,"",0.10);
  double* xpeaks=new double();
  xpeaks=s->GetPositionX();
  std::sort(xpeaks,xpeaks+nfound);
  mNWave=nfound;
  if(nfound==3){
    if(std::abs(xpeaks[2]-xpeaks[1]<35)) mNWave=2;
  }
  tryFit(h1,xpeaks,yped,ymax);

  return true;
}
//void ClusterCsI::findLocalMax(TH1D* h1){
void ClusterCsI::findLocalMax(shared_ptr<TH1D>h1){
  TF1* f1=h1->GetFunction("waveCut");
  int nBin=h1->GetNbinsX();
  double maxValue=0;
  double newPeak;
  for(int iBin=1;iBin<=nBin;iBin++){
    float vBin=h1->GetBinContent(iBin);
    float value=f1->Eval(h1->GetBinCenter(iBin));
    if(vBin-value>maxValue){
      maxValue=vBin-value;
      newPeak=h1->GetBinCenter(iBin);
    }
  }
  mListLocalMax.push_back(newPeak);
}
//void ClusterCsI::tryFit(TH1D* h1,double* xval,double yped,double ymax){
void ClusterCsI::tryFit(shared_ptr<TH1D> h1,double* xval,double yped,double ymax){
  clus_csi=false; // really just in case
  shared_ptr<WaveformCsI> myWave(new WaveformCsI(mNWave));
  double ymax2, ymax3;
  unsigned int NPar=9;
  shared_ptr<TF1> f1(new TF1("waveCut",myWave,&WaveformCsI::waveformSingle,1,250,NPar));//"WaveformCsI","waveformCut"));
  switch(mNWave){
    case 11:
      //TF1* f1=new TF1("waveCut",singlemodel().c_str(),1,250);//"WaveformCsI","waveformCut"));
      for(int n=0; n<9; n+=1){
        f1->SetParameter(n,param[n]);
        f1->SetParLimits(n,parLowlim[n],parUplim[n]);
      }
      std::cout<<" ---- fitting max Bin:  "<<xval[0]<<std::endl;
      f1->SetParameter(0,ymax);
      f1->SetParLimits(0,ymax-61.7,ymax+971.7);
      f1->SetParameter(1,xval[0]+.1);
      f1->SetParLimits(1,xval[0]-261.7,xval[0]+571.7);
      f1->SetParameter(8,yped);
      f1->SetParLimits(8,yped-161.7,yped+171.7);
      f1->SetLineStyle(6);
      f1->SetLineColor(1);
      f1->SetLineWidth(3);
      h1->Fit("waveCut","Q");
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      if(mEventNo % 100==0)
        drawWaves(h1);
      std::cout<<" cluster Fitting |--> "<<f1->GetMaximumX(xval[0]-10,xval[0]+13)<<std::endl;
      std::cout<<" cluster Fitting |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
      break;
    case 2:
      NPar=11;
      f1.reset(new TF1("waveCut",myWave,&WaveformCsI::waveformDouble,1,250,NPar));//"WaveformCsI","waveformCut"));
      //TF1* f1=new TF1("waveCut",singlemodel().c_str(),1,250);//"WaveformCsI","waveformCut"));
      for(UInt_t n=0; n<NPar+1; n+=1){
        f1->SetParameter(n,param[n]);
        f1->SetParLimits(n,parLowlim[n],parUplim[n]);
      }
      std::cout<<" ---- Waveform 2 max Bin:  "<<xval[0]<<" "<<xval[1]<<std::endl;
      ymax2=h1->GetBinContent(h1->FindBin(xval[1]));
      f1->SetParameter(0,ymax);
      f1->SetParLimits(0,ymax-61.7,ymax+971.7);
      f1->SetParameter(1,xval[0]+.1);
      f1->SetParLimits(1,xval[0]-261.7,xval[0]+571.7);
      f1->SetParameter(8,yped);
      f1->SetParLimits(8,yped-161.7,yped+11.7);
      f1->SetParameter(9,xval[1]-15.1);
      f1->SetParLimits(9,xval[1]-61.7,xval[1]+17.7);
      f1->SetParameter(10,ymax2);
      f1->SetParLimits(10,ymax2-7.7,ymax2+7.7);
      f1->SetLineStyle(6);
      f1->SetLineColor(1);
      f1->SetLineWidth(3);
      h1->Fit("waveCut","Q");
      if(f1->GetMaximumX()>=60 && f1->GetMaximumX()<=65){
	if(!clus_csi)
          clus_csi=true;
        lmax=f1->GetMaximum(),lmin=f1->GetMinimum();
        energy=(lmax-lmin)*pcal;
	calcThetaPhi(energy);
        //if(mEventNo % 100==0)
          drawWaves(h1);
        std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX(xval[0]-10,xval[0]+13)<<std::endl;
        std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      break;
    case 13:
      NPar=11;
      f1.reset(new TF1("waveCut",myWave,&WaveformCsI::waveformDouble,1,250,NPar));//"WaveformCsI","waveformCut"));
      //TF1* f1=new TF1("waveCut",singlemodel().c_str(),1,250);//"WaveformCsI","waveformCut"));
      for(UInt_t n=0; n<NPar+1; n+=1){
        f1->SetParameter(n,param[n]);
        f1->SetParLimits(n,parLowlim[n],parUplim[n]);
      }
      std::cout<<" ---- Waveform 2 max Bin:  "<<xval[0]<<" "<<xval[1]<<std::endl;
      ymax2=h1->GetBinContent(h1->FindBin(xval[1]));
      ymax3=h1->GetBinContent(h1->FindBin(xval[2]));
      f1->SetParameter(0,ymax);
      f1->SetParLimits(0,ymax-61.7,ymax+971.7);
      f1->SetParameter(1,xval[0]+.1);
      f1->SetParLimits(1,xval[0]-261.7,xval[0]+571.7);
      f1->SetParameter(8,yped);
      f1->SetParLimits(8,yped-161.7,yped+11.7);
      f1->SetParameter(9,xval[1]-15.1);
      f1->SetParLimits(9,xval[1]-61.7,xval[1]+17.7);
      f1->SetParameter(10,ymax2);
      f1->SetParLimits(10,ymax2-7.7,ymax2+7.7);
      f1->FixParameter(12,xval[2]-15.1);
      //f1->SetParLimits(12,xpeaks[2]-61.7,xpeaks[2]+71.7);
      f1->SetParameter(11,ymax3);
      f1->SetParLimits(11,ymax3-41.7,ymax3+25.77);
      f1->SetLineStyle(6);
      f1->SetLineColor(1);
      f1->SetLineWidth(3);
      h1->Fit("waveCut","Q");
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      //if(mEventNo % 100==0)
        drawWaves(h1);
      std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX(xval[0]-10,xval[0]+13)<<std::endl;
      std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
      break;
  }// end of switch statement
}
void ClusterCsI::calcThetaPhi(double Edep){
  mapPhi=2*180*(phiIndex-0.5)/48.; // from mapping init file
  wz=crysZ[thetaIndex-1];
  wr=crysr[thetaIndex-1];
  Theta=crysZ[thetaIndex-1]/std::sqrt(std::pow(crysZ[thetaIndex-1],2)+
        std::pow(crysr[thetaIndex-1],2));
  acos=TMath::ACos(Theta);
  wtheta=TMath::RadToDeg()*acos; // convert to deg.
  wtheta=round(wtheta); // make sure to obtain 2 dp
  wphi=90.-mapPhi; // world phi
  if(wphi<0) wphi=360.+wphi;
  angles=std::make_pair(wtheta,wphi);
  //h2ang->Fill(wtheta,wphi);
  cout<< " *** World Angles  "<<wtheta<<", "<<wphi<<endl;
}
//void ClusterCsI::drawWaves(TH1D* h1){
void ClusterCsI::drawWaves(shared_ptr<TH1D> h1){
  char name[256];
  //sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,mName);
  TCanvas* c1=new TCanvas(name,name,1200,900);

  if(h1!=0){
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();
    h1->SetMarkerStyle(2);
    h1->SetMarkerSize(1.2);
    h1->Draw();
    /*
    WaveformCsI* myWave=new WaveformCsI(mNWave);
    TF1* f2=new TF1("waveform",myWave,&WaveformCsI::waveformSingle,1,250,mNWave*7+1,"WaveformCsI","waveform");
    Double_t par[mNWave*7+1];
    par[mNWave*7]=h1->GetFunction("waveCut")->GetParameter(mNWave*7);
    for(unsigned int i=0;i<mNWave;i++){
      //      TF1* f3=new TF1(name,derivativeSingle,1,250,7);
      for(int iPar=0;iPar<7;iPar++){
	double myPar=h1->GetFunction("waveCut")->GetParameter(i*7+iPar);
	f2->FixParameter(i*7+iPar,myPar);
	//	f3->FixParameter(iPar,h1->GetFunction("waveCut")->GetParameter(i*7+iPar));
	
      }
      //      f3->SetLineColor(kRed);
      //      h1->GetListOfFunctions()->Add(f3);

      for(unsigned int i=0;i<7*mNWave;i++){
	par[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      Double_t time=findTime(&par[i*7]);
      TLine* line=new TLine(time,0,time,1000);
      line->Draw("same");
    }
    f2->FixParameter(mNWave*7,par[mNWave*7]);
    f2->SetLineStyle(9);
    f2->SetLineColor(kBlue);
    h1->GetListOfFunctions()->Add(f2);
    pad1->Update();
    c1->cd();

    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad2->Draw();
    pad2->cd();
    Double_t ymin = pad1->GetFrame()->GetY1();
    Double_t ymax = pad1->GetFrame()->GetY2();
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom 
    Double_t xmin = pad1->GetFrame()->GetX1();
    Double_t xmax = pad1->GetFrame()->GetX2();
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy-par[mNWave*7],xmax+0.1*dx,ymax+0.1*dy-par[mNWave*7]);
    for(unsigned int i=0;i<mNWave;i++){
      Double_t par[7];
      for(int j=0;j<7;j++){
	par[i*7+j]=h1->GetFunction("waveCut")->GetParameter(i*7+j);
      }
      Double_t energy=findEnergy(&par[i*7]);
      sprintf(name,"energy_run%d_%dCsI_%s_%i",mRunNo,mEventNo,nameCsI(mIndexCsI),i);
      TH1D* h1Energy=new TH1D(name,"energy",100,par[i*7+1],findZero(&par[i*7]));
      Double_t mySum=fillEnergy(&par[i*7],int(1000000),h1Energy);
      h1Energy->Scale(energy/mySum);
      h1Energy->SetFillColorAlpha(i+1, 0.35);
      h1Energy->SetFillStyle(3010);
      
      h1Energy->SetLineColorAlpha(0,0.2);
      //      h1Energy->Draw("][same");
      pad2->Update();
    }*/
  }
  c1->Write();
  //sprintf(name,"wave_run%d_%dCsI_%s.png",mRunNo,mEventNo,nameCsI(mIndexCsI));
  //c1->SaveAs(name);
}
//Double_t ClusterCsI::findChi2(TH1D* h1){
Double_t ClusterCsI::findChi2(shared_ptr<TH1D> h1){
  unsigned int iLow=1;
  unsigned int iUp=h1->GetNbinsX();
  TF1* funFit=h1->GetFunction("waveCut");
  Double_t chi2=0;
  for(unsigned int i=iLow;i<=iUp;i++){
    Double_t x=h1->GetBinCenter(i);
    Double_t y=h1->GetBinContent(i);
    Double_t error=h1->GetBinError(i);
    Double_t value=funFit->Eval(x);
    Double_t thisChi2=(value-y)*(value-y)/error/error;
    chi2+=thisChi2;
  }
  chi2/=(iUp-1-h1->GetFunction("waveCut")->GetNpar());
  return chi2;
}
