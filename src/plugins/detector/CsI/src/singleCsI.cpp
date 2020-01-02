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
Double_t WaveformCsI::waveformSingle(Double_t *x,Double_t *par){
  //if(x[0]<par[1]) return 0;
  //std::cout<<"  ----> entering single waveform fit function \n";
  double den=1-exp(-(x[0]-par[1])/par[6]);
  double num=(par[0]/den)*TMath::Freq((x[0]-par[1]-par[2])/par[3])*
             ((x[0]-par[1])/par[4]*exp(1-(x[0]-par[1])/par[4])+
             par[5]*(x[0]-par[1])/(par[7])*exp(1-(x[0]-par[1])/(par[7])));
  return num+par[8];
}
Double_t WaveformCsI::waveformDouble(Double_t *x,Double_t* par){
  //if(x[0]<par[1]) return 0;
  double den1=1-exp(-1*(x[0]-par[1])/par[6]);
  double den2=1-exp(-1*(x[0]-par[9])/par[6]);
  double frac=par[0]*2*TMath::Freq((x[0]-par[1]-par[2])/par[3])*
              ((x[0]-par[1])/par[4]*exp(1-(x[0]-par[1])/par[4])+
              par[5]*(x[0]-par[1])/(par[7])*exp(1-(x[0]-par[1])/(par[7])))/(den1)+
	      par[10]*2*TMath::Freq((x[0]-par[9]-par[2])/par[3])*((x[0]-par[9])/par[4]*
              exp(1-(x[0]-par[9])/par[7])+par[5]*(x[0]-par[9])/(par[7])*
              exp(1-(x[0]-par[9])/par[7]))/(den2);
  return frac+par[8];
}
Double_t WaveformCsI::waveformTriple(Double_t *x,Double_t* par){
  //if(x[0]<par[1]) return 0;
  double den1=1-exp(-1*(x[0]-par[1])/par[6]);
  double den2=1-exp(-1*(x[0]-par[9])/par[6]);
  double den3=1-exp(-1*(x[0]-par[12])/par[6]);
  double frac=par[0]*2*TMath::Freq((x[0]-par[1]-par[2])/par[3])*
              ((x[0]-par[1])/par[4]*exp(1-(x[0]-par[1])/par[4])+
              par[5]*(x[0]-par[1])/(par[7])*exp(1-(x[0]-par[1])/(par[7])))/(den1)+
              par[10]*2*TMath::Freq((x[0]-par[9]-par[2])/par[3])*
              ((x[0]-par[9])/par[4]*exp(1-(x[0]-par[9])/par[7]) +
              par[5]*(x[0]-par[9])/(par[7])*exp(1-(x[0]-par[9])/(par[7])))/(den2) +
              par[11]*2*TMath::Freq((x[0]-par[12]-par[2])/par[3])*
              ((x[0]-par[12])/par[4]*exp(1-(x[0]-par[12])/par[7]) +
              par[5]*(x[0]-par[12])/(par[7])*exp(1-(x[0]-par[12])/(par[7])))/(den3);
  return frac+par[8];
}
Double_t WaveformCsI::waveformQuad(Double_t *x,Double_t* par){
  //if(x[0]<par[1]) return 0;
  double den1=1-exp(-1*(x[0]-par[1])/par[6]);
  double den2=1-exp(-1*(x[0]-par[9])/par[6]);
  double den3=1-exp(-1*(x[0]-par[12])/par[6]);
  double den4=1-exp(-1*(x[0]-par[14])/par[6]);
  double frac=par[0]*2*TMath::Freq((x[0]-par[1]-par[2])/par[3])*
              ((x[0]-par[1])/par[4]*exp(1-(x[0]-par[1])/par[4])+
              par[5]*(x[0]-par[1])/(par[7])*exp(1-(x[0]-par[1])/(par[7])))/(den1)+
              par[10]*2*TMath::Freq((x[0]-par[9]-par[2])/par[3])*
              ((x[0]-par[9])/par[4]*exp(1-(x[0]-par[9])/par[7])+
              par[5]*(x[0]-par[9])/(par[7])*exp(1-(x[0]-par[9])/(par[7])))/(den2)+
              par[11]*2*TMath::Freq((x[0]-par[12]-par[2])/par[3])*
              ((x[0]-par[12])/par[4]*exp(1-(x[0]-par[12])/par[7])+
              par[5]*(x[0]-par[12])/(par[7])*exp(1-(x[0]-par[12])/(par[7])))/(den3)+
              par[13]*2*TMath::Freq((x[0]-par[14]-par[2])/par[3])*
              ((x[0]-par[14])/par[4]*exp(1-(x[0]-par[14])/par[7])+
              par[5]*(x[0]-par[14])/(par[7])*exp(1-(x[0]-par[14])/(par[7])))/(den4);    
  return frac+par[8];
}
Double_t WaveformCsI::waveformOverrange(Double_t *x,Double_t* par){
  //if(x[0]<par[1]) return 0;
  double den=1-exp(-(x[0]-par[1])/par[6]);
  double num=par[0]*2*TMath::Freq((x[0]-par[1]-par[2])/par[3])*
             ((x[0]-par[1])/par[4]*exp(1-(x[0]-par[1])/par[4])+
             par[5]*(x[0]-par[1])/(par[7])*exp(1-(x[0]-par[1])/(par[7])));
  return (num/(den)+par[8]);
}
Double_t WaveformCsI::waveformCut(Double_t *x,Double_t* par){
  //Double_t value=waveform(x,par);
  //if(value>1023.0) value=1023.0;
  //return value;
  return 0;
}
void SingleCsI::setData(const vector<UShort_t>& listD){
  for(unsigned int iS=0,nS=listD.size();iS<nS;iS++){
    mListData.push_back(listD[iS]);
  }
}
static Double_t waveformPure(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return 0;
  Double_t termFirst=par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[6])));
  Double_t termSecond=(x[0]-par[1])/(par[3]*par[3]);
  Double_t termThird=TMath::Exp(-(x[0]-par[1])/par[3])+par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t value=termFirst*termSecond*termThird;
  return value;
}
static Double_t findZero(Double_t *par){
  Double_t x=par[6];
  Double_t stepSize=0.01;
  while(true){
    Double_t value=waveformPure(&x,par);
    if(value<1) return x;
    x+=stepSize;
  }
}
static Double_t derivativeSingle(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return 0;
  Double_t firstExp=par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t secondExp=TMath::Exp(-(x[0]-par[1])/par[3]);
  Double_t thirdExp=TMath::Exp(par[2]*par[6]-par[2]*x[0]);
  Double_t termFirst=par[0]*(x[0]-par[1])*(-firstExp/par[5]-secondExp/par[3])/(par[3]*par[3]*(thirdExp+1));
  Double_t termSecond=par[0]*(firstExp+secondExp)/(par[3]*par[3]*(thirdExp+1));
  Double_t termThird=par[0]*par[2]*(x[0]-par[1])*(firstExp+secondExp)*thirdExp/(par[3]*par[3]*(thirdExp+1)*(thirdExp+1));
  Double_t value=termFirst+termSecond+termThird;
  return value;
}
static Double_t findTime(Double_t *par){
  Double_t x=par[1];
  Double_t slopePre=0;
  const Double_t stepSize=0.001;
  while(true){
    Double_t value=derivativeSingle(&x,par);
    //    cout<<value<<" "<<slopePre<<endl;
    if(value<slopePre) return x-stepSize/2.0;
    slopePre=value;
    x+=stepSize;
    if(x>250) return 250;
  }
}
static Double_t findEnergy(Double_t *par){
  Double_t xLow=par[1];
  Double_t xHigh=findZero(par);

  TF1* f1=new TF1("pure",waveformPure,1,250,7);
  Double_t epsilon=0.01;
  Double_t value=f1->Integral(xLow,xHigh,epsilon); // temp correction
  //Double_t value=f1->Integral(xLow,xHigh,par,epsilon); //Needs fixing!!
  return value;
}
static Double_t findPeak(Double_t *par){
  Double_t x=par[1];
  bool wasPositive=false;
  const Double_t stepSize=0.001;
  while(true){
    Double_t value=derivativeSingle(&x,par);
    if(!wasPositive){
      if(value>0) wasPositive=true;
    }
    else{
      if(value<0) return x-stepSize/2.0;
    }
    x+=stepSize;
  }
}
static Double_t fillEnergy(Double_t *par,unsigned int nEvent,TH1D* h1){
  Double_t xLow=par[1];
  Double_t xHigh=findZero(par);
  Double_t xPeak=findPeak(par);
  TRandom3* poolRandom=new TRandom3();
  for(unsigned int i=0;i<nEvent;i++){
    Double_t x=poolRandom->Uniform(xLow,xHigh);
    h1->Fill(x,waveformPure(&x,par));
  }
  return h1->Integral("width");
}
char* SingleCsI::nameCsI(const UInt_t& iClock, const UInt_t& iFB,const UInt_t& iUD,const UInt_t& iModule){
  char ifb[3]="fb";
  char iud[3]="ud";
  sprintf(mName,"%i%c%c%i",iClock,ifb[iFB],iud[iUD],iModule);
  //std::cout<<" checking the name: "<<mName<<std::endl;
  return mName;
}
char* SingleCsI::nameCsI(unsigned int index){
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
void SingleCsI::initVar(){
  wtheta=dummy;
  wphi=dummy;
  energy=dummy;
  phdiff=dummy;
  wr=dummy;
  wz=dummy;
  rtime=dummy;
  ptime=dummy;
  cdf50=dummy;
  chi2=dummy;
  ndf=dummy;
  mNWave=dummy;
}
bool SingleCsI::fit(){
  //std::cout<<"|| entering into fitting function..... success! \n";
  initVar(); // will be called at begining of every event
  static unsigned int count=0;
  char name[256];
  //cout<<"processing "<<mIndexCsI<<" "<<nameCsI(mIndexCsI)<<endl;
  cout<<"processing "<<mIndexCsI<<" "<<mName<<endl;
  std::cout<<" This is the add data method, size: "<<mListData.size()<<std::endl;
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
  if(ymax==1023){
    mNWave=5;
    // We need to make sure that the max extends for at least 3 bins
    // If not, treat as non-max wave pulse
    x1=h1->FindFirstBinAbove(1022);
    x2=h1->FindLastBinAbove(1022);
    diffMax=std::abs(x2-x1);
    if(diffMax<=3) mNWave=1;
    if(nfound==2){
      mNWave=6;
      if(diffMax<=3) mNWave=2;
      //if(diffMax>3 && diffMax<35) mNWave=5;
    }
    //return false;
  }
  if(nfound==3){
    if(std::abs(xpeaks[2]-xpeaks[1]<35)) mNWave=2;
  }
  tryFit(h1,xpeaks,yped,ymax);
  //delete xpeaks;
  /*
  Double_t chi2Reduced=findChi2(h1);
  cout<<chi2Reduced<<endl;
  if(chi2Reduced>80){
    mNWave=2;
    tryFit(h1);
  Double_t chi2Reduced=findChi2(h1);
  cout<<chi2Reduced<<endl;

  }
  for(mNWave=1;mNWave<2;mNWave++){
    tryFit(h1);
    Double_t chi2Reduced=findChi2(h1);
    if(chi2Reduced<100){
      cout<<chi2Reduced<<endl;
      break;
    }
  }
  */

  return true;
  /*

  Double_t par[8];
  for(int i=0;i<8;i++){
    par[i]=h1->GetFunction("wave")->GetParameter(i);
  }
  Double_t time=findTime(par);
  TLine* line=new TLine(time,0,time,1000);
  Double_t peak=findPeak(par);
  Double_t valuePeak=waveformSingle(&peak,par);
  TLine* line2=new TLine(peak,0,peak,valuePeak);
  line2->SetLineColor(kRed);
  line->Draw("same");
  line2->Draw("same");
  Double_t energy=findEnergy(par);
  sprintf(name,"energy_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  TH1D* h1Energy=new TH1D(name,"energy",5000,par[1],findZero(par));
  Double_t mySum=fillEnergy(par,int(10000000),h1Energy);
  cout<<mySum<<" "<<energy<<endl;
  h1Energy->Scale(energy/mySum);
  h1Energy->SetFillColor(1);
  h1Energy->SetFillStyle(3003);
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
  pad2->Range(xmin-0.1*dx,ymin-0.1*dy-par[7],xmax+0.1*dx,ymax+0.1*dy-par[7]);

  h1Energy->SetLineColorAlpha(0,0.2);
  h1Energy->Draw("][same");
  pad2->Update();
  */
}
void SingleCsI::findLocalMax(shared_ptr<TH1D>h1){
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
std::string singlemodel(){
  char norm[1000],modelname[1000];
  sprintf(norm,"1-exp(-(x-[1])/[6])");
  
  sprintf(modelname,"[0]*TMath::Freq((x-[1]-[2])/[3])/(%s)*((x-[1])/[4]*exp(1-(x-[1])/[4])+\
        [5]*((x-[1])/([7]))*exp(1-(x-[1])/[7]))+[8]",norm);
  std::string namestring=modelname;
  return namestring;
}
void SingleCsI::tryFit(shared_ptr<TH1D> h1,double* xval,double yped,double ymax){
  shared_ptr<WaveformCsI> myWave(new WaveformCsI(mNWave));
  unsigned int NPar=9;
  shared_ptr<TF1> f1(new TF1("waveCut",myWave,&WaveformCsI::waveformSingle,1,250,NPar));//"WaveformCsI","waveformCut"));
  switch(mNWave){
    case 1:
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
      if(f1->GetMaximumX()>=60 && f1->GetMaximumX()<=65){
        lmax=f1->GetMaximum();lmin=f1->GetMinimum();
        energy=(lmax-lmin);
	calcThetaPhi();
	phdiff=(lmax-lmin);
        tcalc=.5*(lmax-lmin);
        rtime=f1->GetParameter(1);
        cdf50=f1->GetX(tcalc);
	ptime=f1->GetMaximumX();
	chi2=f1->GetChisquare();
	ndf=f1->GetNDF();
        std::cout<<" -------------------- "<<cdf50<<std::endl;
	//calTime(f1);
        if(mEventNo % 10000==0)
          drawWaves(h1);
        std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX(xval[0]-10,xval[0]+13)<<std::endl;
        std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      break;
    case 2:
      NPar=11;
      f1.reset(new TF1("waveCut",myWave,&WaveformCsI::waveformDouble,1,250,NPar));
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
        lmax=f1->GetMaximum();lmin=f1->GetMinimum();
        energy=(lmax-lmin);
	calcThetaPhi();
	phdiff=(lmax-lmin);
        tcalc=.5*(lmax+lmin);
        rtime=f1->GetParameter(1);
        cdf50=f1->GetX(tcalc);
	ptime=f1->GetMaximumX();
	chi2=f1->GetChisquare();
	ndf=f1->GetNDF();
        std::cout<<" -------------------- "<<cdf50<<std::endl;
	//calTime(f1);
        if(mEventNo % 10000==0)
          drawWaves(h1);
        //std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      break;
    case 3:
      NPar=13;
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
      if(f1->GetMaximumX()>=60 && f1->GetMaximumX()<=65){
        lmax=f1->GetMaximum();lmin=f1->GetMinimum();
        energy=(lmax-lmin);
	phdiff=(lmax-lmin);
	calcThetaPhi();
        tcalc=.5*(lmax+lmin);
        rtime=f1->GetParameter(1);
        cdf50=f1->GetX(tcalc);
	ptime=f1->GetMaximumX();
	chi2=f1->GetChisquare();
	ndf=f1->GetNDF();
        std::cout<<" -------------------- "<<cdf50<<std::endl;
	//calTime(f1);
        if(mEventNo % 10000==0)
          drawWaves(h1);
        //std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      break;
    case 5:
      NPar=10;
      f1.reset(new TF1("waveCut",myWave,&WaveformCsI::waveformOverrange,1,250,NPar));
      for(UInt_t n=0; n<NPar+1; n+=1){
        f1->SetParameter(n,param[n]);
        f1->SetParLimits(n,parLowlim[n],parUplim[n]);
      }
      ymax2=h1->GetBinContent(h1->FindBin(xval[1]));
      f1->SetParameter(0,ymax*10.5);
      f1->SetParLimits(0,ymax-61.7,1e6);
      f1->SetParameter(1,xval[0]+40.1);
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
        lmax=f1->GetMaximum();lmin=f1->GetMinimum();
        energy=(lmax-lmin);
	phdiff=(lmax-lmin);
	calcThetaPhi();
        tcalc=.5*(lmax+lmin);
        rtime=f1->GetParameter(1);
        cdf50=f1->GetX(tcalc);
	ptime=f1->GetMaximumX();
	chi2=f1->GetChisquare();
	ndf=f1->GetNDF();
        std::cout<<" ------| 1023 |------ "<<cdf50<<std::endl;
	//calTime(f1);
        if(mEventNo % 1000==0)
          drawWaves(h1);
        std::cout<<"*|First bin: "<<h1->FindFirstBinAbove(1022)<<" || "<<h1->FindLastBinAbove(1022)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      break;
    case 6:
      NPar=11;
      f1.reset(new TF1("waveCut",myWave,&WaveformCsI::waveformDouble,1,250,NPar));
      for(UInt_t n=0; n<NPar+1; n+=1){
        f1->SetParameter(n,param[n]);
        f1->SetParLimits(n,parLowlim[n],parUplim[n]);
      }
      std::cout<<" ---- Waveform 2 max Bin:  "<<xval[0]<<" "<<xval[1]<<std::endl;
      ymax2=h1->GetBinContent(h1->FindBin(xval[1]));
      f1->SetParameter(0,ymax*10.5);
      f1->SetParLimits(0,ymax-61.7,ymax*5.4);
      f1->SetParameter(1,xval[0]+40.1);
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
        lmax=f1->GetMaximum();lmin=f1->GetMinimum();
        energy=(lmax-lmin);
	//calcThetaPhi(energy);
	phdiff=(lmax-lmin);
        tcalc=.5*(lmax+lmin);
	calcThetaPhi();
        rtime=f1->GetParameter(1);
        cdf50=f1->GetX(tcalc);
	ptime=f1->GetMaximumX();
	chi2=f1->GetChisquare();
        std::cout<<" -------------------- "<<cdf50<<std::endl;
	//calTime(f1);
        if(mEventNo % 1000==0)
          drawWaves(h1);
        //std::cout<<" cluster From CsI |--> "<<f1->GetMaximumX()<<" | baseline "<<f1->GetParameter(8)<<std::endl;
        //std::cout<<" ******  Checking the energy   -->"<<energy<<" || "<<(f1->GetParameter(0)-f1->GetParameter(8))*pcal<<"\n";
      }
      for(unsigned int i=0;i<NPar;i++){
        mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
      }
      break;
  }// end of switch statement
}
void SingleCsI::calcThetaPhi(){
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
  //angles=std::make_pair(wtheta,wphi);
  //h2ang->Fill(wtheta,wphi);
  cout<< " *** World Angles  "<<wtheta<<", "<<wphi<<endl;
}
void SingleCsI::drawWaves(shared_ptr<TH1D> h1){
  char name[256];
  //sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,mName);
  //TCanvas* c1=new TCanvas(name,name,1200,900);

  if(h1!=0){
    //TPad *pad1 = new TPad("pad1","",0,0,1,1);
    //pad1->Draw();
    //pad1->cd();
    h1->SetMarkerStyle(2);
    h1->SetMarkerSize(1.2);
    h1->Write();
    //h1->Draw();
  }
  //c1->Write();
  //sprintf(name,"wave_run%d_%dCsI_%s.png",mRunNo,mEventNo,nameCsI(mIndexCsI));
  //c1->SaveAs(name);
}
/*
Double_t SingleCsI::findChi2(shared_ptr<TH1D> h1){
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
}*/
