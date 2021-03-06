#ifndef evalClusters_H
#define evalClusters_H

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <TLine.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TH2D.h>
#include "Plugin.h"

class evalClusters{
public:
  evalClusters();
  ~evalClusters();
  void defHistos();
  // fill histos for sanity checks
  //void fillCutHistos(double InvMass,double OpAng1,double Etot,double primOpAng);
  void fillHistos(double InvMass,double OpAng1,double Etot,double primOpAng,int n);
  void fillHistos(double InvMass,double Etot,double primOpAng);
  void fillHistos(double InvMass2,int n);
  // fill histos for centroid angles
  //void fillHistos(double theta1,double phi1,double theta2,double phi2);
  void fillHistos(double theta1,double phi1,double theta2,double phi2,double theta3,double phi3);
  void drawHistos();
  void drawCanvas(TH1D* h,int val);
  void drawCanvas(TH2D* h,int val);
  void drawCanvas(TH1D* hist1,TH1D* hist2,TH1D* hist3,TH1D* hist4,int val);
  void setChannel(int val); // select channel no. for eval
  void fillMltp(int val);
  void singleCanvas();
  void singleCanvas(TH1D*);
  // functions to define and fill 1D histograms
  void defHistos(std::string name,int bin,double xlow,double xhigh);
  Long_t fillHistos(double data);
  void writeHistos();
private:
  int channelNo; //used for iterating TCanvas and selecting No. of cluster to analyze.
  int Ncrys, Nclust;
  std::string title1, title2, title3, title4, title5;
  std::string prname1, prname2, clname1, clname2, clname3;
  std::string invm, clang, clE, prang, M2;
  double invMass, cpidOpAng, primOpAng, E2clust, E3clust;
  double cpidtheta1, cpidtheta2, cpidtheta3;
  double cpidphi1, cpidphi2, cpidphi3;
  TH1D* clustAng[3], *primAng[3], *Eclust[3], *invM[3];
  TH1D* h1M2[3];
  TH1D* h1clstAng, *h1prmAng, *h1Eclust, *h1invM;
  TH1D* cl1E, *cl2E, *cl3E;
  TH1D* clMltp;
  TH1D* h1Hist;
  TH2D* thetaPhi;
  // just in case I need this later down the line
  double cpid1x, cpid1y, cpid1z;
  double cpid2x, cpid2y, cpid2z;
  double cpid3x, cpid3y, cpid3z;
};
#endif
