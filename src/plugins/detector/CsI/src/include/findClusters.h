#ifndef findClusters_h
#define findClusters_h 1

#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <TMath.h>

class findClusters{
public:
  findClusters();
  ~findClusters();
  void clusters(std::map<std::pair<double,double>,double> csiClust,const std::vector<double> &csiE,std::vector<double> csiTheta,std::vector<double> csiPhi,std::map<std::pair<double,double>,bool> csiCheck);
  void clusters(std::map<std::pair<double,double>,double> csiClust,const std::vector<double> &csiE,std::vector<double> csiTheta,std::vector<double> csiPhi,std::map<std::pair<double,double>,double> csiR,std::map<std::pair<double,double>,double>csiZ,std::map<std::pair<double,double>,bool> csiCheck);
  inline int getSingleCrysClust(){return numOfsingleClus;};
  inline int getMultiCrysClust(){return numOfClus;};
  // retrieval function for single cluster variables
  inline std::vector<double> getSingleE(){return singleEne;};
  inline std::vector<double> getSingleTheta(){return singTheta;};
  inline std::vector<double> getSinglePhi(){return singPhi;};
  inline std::vector<double> getSingleR(){return singR;};
  inline std::vector<double> getSingleZ(){return singZ;};
  // retrieval function for single cluster variables
  inline std::vector<double> getMultiE(){return clusEne;};
  inline std::vector<double> getMultiTheta(){return clusThetaE;};
  inline std::vector<double> getMultiPhi(){return clusPhiE;};
  inline std::vector<double> getMultiR(){return clusEr;};
  inline std::vector<double> getMultiZ(){return clusEz;};
private:
  typedef std::vector<double> ve;
  //extern ve indexph;
  ve indexph;
  // function that returns highest value of std::vector
  // then iterate in decending order from max->min
  std::size_t get_nthIndex(const std::vector<double> &vect, std::size_t k);
  std::vector<double> *phval, *clusth, *clusphi;
  std::pair<double,double> angP1, angP2, angP3, angP4;
  std::pair<double,double> angP5, angP6, angP7, angP8;
  std::pair<double,double> angE1, angE2, angE3, angE4, angE5;
  std::pair<double,double> angE6, angE7, angE8, angE9, angE10;
  std::pair<double,double> tppair;
  double posTheta, negTheta, posPhi, negPhi; // adding or subtracting 7.5*deg
  double rtheta,rphi,z_w,r_w;
  std::vector<double> csThet,csPhi;
  std::vector<double> singTheta, singPhi, singleEne,singR,singZ;
  std::vector<double> clusEne,clusThetaE,clusPhiE,clusEr,clusEz;
  std::map<std::pair<double,double>,bool> csiCheck;
  double ntheta,nphi,Eclus,thetaE,phiE,clusZ,clusR;
  int clusCrys;
  int numOfsingleClus,numOfClus;
};
#endif
