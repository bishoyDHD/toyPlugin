#include "evalClusters.h"

evalClusters::evalClusters():channelNo(7){

}
evalClusters::~evalClusters(){

}
void evalClusters::setChannel(int val=7){
  channelNo=val;
}
void evalClusters::defHistos(){
  switch(channelNo){
    case 7:
      title1="Opening angle of 2#gamma";
      title2="Invariant Mass of #pi^{0}";
      title3="Energy of 2#gamma";
      title4="Opening angle of  #pi^{+}#pi^{0}";
      prname1="#pi^{+}";
      prname2="#pi^{0}";
      clname1="#gamma";
      clname2="#gamma";
      break;
    case 71: // should really be 7.1 for pi0->e+e-gamma. switch only takes int
      title1="Opening angle of 2#gamma";
      title2="Invariant Mass of #pi^{0}";
      title3="Energy of 2#gamma";
      title4="Opening angle of  #pi^{+}#pi^{0}";
      prname1="#pi^{+}";
      prname2="#pi^{0}";
      clname1="e^{+}";
      clname2="e^{-}";
      clname3="#gamma";
      break;
    case 14: // for A' (channel 14)
      title1="Opening angle of e^{+}e^{-}";
      title2="Invariant Mass of A'";
      title3="Energy of e^{+}e^{-}";
      title4="Opening angle of #mu^{+}A'";
      prname1="#mu^{+}";
      prname2="A'";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    case 15: // for A' (channel 15)
      title1="Opening angle of e^{+}e^{-}";
      title2="Invariant Mass of A'";
      title3="Energy of e^{+}e^{-}";
      title4="Opening angle of #pi^{+}A'";
      prname1="#pi^{+}";
      prname2="A'";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    case 16: // for SM (channel 16)
      title1="Opening angle of e^{+}e^{-}";
      //title2="Invariant Mass of ";
      title3="Energy of e^{+}e^{-}";
      //title4="Opening angle of ";
      prname1="#mu^{+}";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    default:
      std::cout<<" *** Warning!\n --- Channel definition is out of range!!\n";
      std::cout<<" *** Warning!\n";
      break;
  }// end of switch statement
  // Define various histograms
  for(int i=0; i<3; i++){
    invm="invMass_", clang="clusAng_", clE="clusE_", prang="primAngle_", M2="mass2_";
    invm+=(i+1)/10+'0', clang+=(i+1)/10+'0', clE+=(i+1)/10+'0', prang+=(i+1)/10+'0', M2+=(i+1)/10+'0';
    invm+=(i+1)%10+'0', clang+=(i+1)%10+'0', clE+=(i+1)%10+'0', prang+=(i+1)%10+'0', M2+=(i+1)%10+'0';
    clustAng[i]=new TH1D(clang.c_str(),title1.c_str(),75.0,-1.1,1.1);
    invM[i]=new TH1D(invm.c_str(),title2.c_str(),84.50,0.0,.210);
    Eclust[i]=new TH1D(clE.c_str(),title3.c_str(),100.,0.0,.300);
    primAng[i]=new TH1D(prang.c_str(),title4.c_str(),75.0,-1.1,1.1);
    h1M2[i]=new TH1D(M2.c_str(),"M_{ee}^{2}",170.0,-0.210,.210);
  }
  // Cut histograms for 2 Cluster events
  h1clstAng=new TH1D("h1Ang",title1.c_str(),75.0,-1.1,1.1);
  h1invM=new TH1D("h1invM",title2.c_str(),84.50,0.0,.210);
  h1Eclust=new TH1D("h1Eclst",title3.c_str(),100.,0.0,.300);
  h1prmAng=new TH1D("h1prmAng",title4.c_str(),75.0,-1.1,1.1);
  // needed regardless of channel no.
  thetaPhi=new TH2D("h2ang", "#theta Vs. #phi",24.0,0.,M_PI,48.0,0.,2.*M_PI);
  clMltp=new TH1D("clusM"," Cluster multiplicity",10,-.5,9.5);
}
// function to plot actual histos
void evalClusters::drawHistos(){
  switch(channelNo){
    case 7:
      drawCanvas(Eclust[0],invM[0],clustAng[0],primAng[0],0);
      drawCanvas(Eclust[1],invM[1],clustAng[1],primAng[1],1);
      //drawCanvas(h1Eclust,h1invM,h1clstAng,h1prmAng,1);
      drawCanvas(thetaPhi,1);
      break;
    case 14:
      std::cout<<" will get to this soon!! \n";
      break;
  }// end of switch statement
}
void evalClusters::fillHistos(double theta1,double phi1,double theta2,double phi2,double theta3,double phi3){
  switch(channelNo){
    case 71:
      thetaPhi->Fill(theta1,phi1);
      thetaPhi->Fill(theta2,phi2);
      thetaPhi->Fill(theta3,phi3);
      break;
    case 7:
      thetaPhi->Fill(theta1,phi1);
      thetaPhi->Fill(theta2,phi2);
      break;
  }
}
// plot cluster multiplicity regardless of channel
void evalClusters::fillMltp(int val){
  clMltp->Fill(val);
}
void evalClusters::fillHistos(double InvMass,double opAng1,double Etot,double primOpAng, int n=0){
  switch(n){
  case 0: // Any number of cluster eval
    invM[0]->Fill(InvMass);
    Eclust[0]->Fill(Etot);
    clustAng[0]->Fill(opAng1);
    primAng[0]->Fill(primOpAng);
    break;
  case 1: // apply cut and eval N_clus==2
    invM[1]->Fill(InvMass);
    Eclust[1]->Fill(Etot);
    clustAng[1]->Fill(opAng1);
    primAng[1]->Fill(primOpAng);
    break;
  case 2: // apply cut and eval N_clus==2
    invM[2]->Fill(InvMass);
    Eclust[2]->Fill(Etot);
    clustAng[2]->Fill(opAng1);
    primAng[2]->Fill(primOpAng);
    break;
  }// end of switch statement
}
void evalClusters::fillHistos(double InvMass2,int n=0){
  switch(n){
  case 0: // Any number of cluster eval
    h1M2[0]->Fill(InvMass2);
    break;
  case 1: // apply cut and eval N_clus==2
    h1M2[1]->Fill(InvMass2);
    break;
  case 2: // apply cut and eval N_clus==2
    h1M2[2]->Fill(InvMass2);
    break;
  }// end of switch statement
}
void evalClusters::drawCanvas(TH2D* h2,int val){
  std::string cname="c_";
  cname+=(val+1)/10+'0';
  cname+=(val+1)%10+'0';
  TCanvas* c1=new TCanvas(cname.c_str(),"Theta vs Phi",900,800);
  c1->cd();
  std::string xtitle="#theta [rad]";
  std::string ytitle="#phi [rad]";
  h2->GetXaxis()->SetTitle(xtitle.c_str());
  h2->GetYaxis()->SetTitle(ytitle.c_str());
  h2->Draw("colz");
  c1->Write();
}
void evalClusters::singleCanvas(){
  TCanvas* c2=new TCanvas("c2","cluster multiplicity",900,800);
  c2->cd();
  clMltp->GetXaxis()->SetTitle("number of #gamma cluster");
  clMltp->GetXaxis()->SetTitleSize(0.045);
  clMltp->GetYaxis()->SetTitle("counts/bin");
  clMltp->Draw("hist");
  c2->Write();
}
void evalClusters::drawCanvas(TH1D* hist1,TH1D* hist2,TH1D* hist3,TH1D* hist4,int val=1){
  std::string cname="c";
  std::string xtitle;
  cname+=(val+1)/10+'0';
  cname+=(val+1)%10+'0';
  TCanvas* c1=new TCanvas(cname.c_str(),"Sanity plots",900,800);
  c1->Divide(2,2);
  c1->cd(1);
  xtitle="E_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="}";
  //hist3->GetXaxis()->SetLabelFont(0.35);
  hist1->GetXaxis()->SetTitle(xtitle.c_str());
  hist1->GetXaxis()->SetTitleSize(0.045);
  hist1->GetYaxis()->SetTitle("counts/bin");
  hist1->Draw("hist");
  c1->cd(2);
  xtitle="M_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="}";
  hist2->GetXaxis()->SetTitle(xtitle.c_str());
  hist2->GetXaxis()->SetTitleSize(0.045);
  hist2->GetYaxis()->SetTitle("counts/bin");
  hist2->Draw("hist");
  c1->cd(3);
  xtitle="cos(#theta_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="})";
  hist3->GetXaxis()->SetTitle(xtitle.c_str());
  hist3->GetXaxis()->SetTitleSize(0.045);
  hist3->GetYaxis()->SetTitle("counts/bin");
  hist3->Draw("hist");
  /*
  // Draw text
  auto tex=new TLatex(0.,3.," Very preliminary");
  tex->SetTextColor(17);
  tex->SetTextSize(1.22);
  tex->SetTextAngle(30.22);
  tex->SetLineWidth(2);
  tex->Draw();*/
  c1->cd(4);
  xtitle="cos(#theta_{";
  xtitle+=prname1;
  xtitle+=prname2;
  xtitle+="})";
  hist4->GetXaxis()->SetTitle(xtitle.c_str());
  hist4->GetXaxis()->SetTitleSize(0.045);
  hist4->GetYaxis()->SetTitle("counts/bin");
  hist4->Draw("hist");
  c1->Write();
}
