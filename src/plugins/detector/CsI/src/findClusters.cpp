#include "findClusters.h"

findClusters::findClusters(){
}
findClusters::~findClusters(){
}
// private function so cannot be accessed outside of class
// since only intended for this class
// ---> All good!
std::size_t findClusters::get_nthIndex(const std::vector<double> &vect,std::size_t k){
  std::vector<std::size_t> indexes(vect.size());
  std::iota(indexes.begin(), indexes.end(), 0);
  std::nth_element(indexes.begin(), indexes.begin() + k, indexes.end(),
    [&](int lhs, int rhs){
      return vect[lhs] > vect[rhs];
    }
  );
  return indexes[k];
}
void findClusters::clusters(std::map<std::pair<double,double>,double> csiClust,const std::vector<double> &csiE,std::vector<double> csiTheta,std::vector<double> csiPhi,std::map<std::pair<double,double>,double> csiR,std::map<std::pair<double,double>,double> csiZ,std::map<std::pair<double,double>,bool> check){
  std::cout<<" size of Ecsi is "<<csiClust.size()<<"\n";
  csiCheck.insert(check.begin(),check.end());
  numOfsingleClus=0; numOfClus=0;
  singleEne.clear();       clusEne.clear();
  singTheta.clear();       clusThetaE.clear();
  singPhi.clear();         clusPhiE.clear();
  clusEr.clear();          clusEz.clear();
  singR.clear();           singZ.clear();
  for(std::size_t mm=0; mm !=csiE.size(); mm++){
    const std::size_t index=get_nthIndex(csiE, mm);
    std::cout<<"  the greater index --> "<<index;// <<std::endl;
    std::cout<<" with value "<<csiE[index]<<std::endl;
    ntheta=csiTheta[index], nphi=csiPhi[index];
    tppair=std::make_pair(ntheta,nphi);
    std::cout<<"   ==>  theta and phi "<<ntheta<<"  "<<nphi<<std::endl;
    negTheta=ntheta-7.5;    posTheta=ntheta+7.5;
    negPhi=nphi-7.5;        posPhi=nphi+7.5;
    // Revolution case
    if(nphi==3.75) negPhi=356.25;
    // Performing search for non-edge case:
    angP1=std::make_pair(posTheta,nphi);        angP2=std::make_pair(negTheta,nphi);
    angP3=std::make_pair(ntheta,posPhi);        angP4=std::make_pair(ntheta,negPhi);
    angP5=std::make_pair(posTheta,posPhi);      angP6=std::make_pair(negTheta,negPhi);
    angP7=std::make_pair(negTheta,posPhi);      angP8=std::make_pair(posTheta,negPhi);
    // accounting for the complicated case of Crystal No. 10 (index=16)
    // This is strange because they are 2X larger in phi
    // from chan16 -> chan15: Back
    if(ntheta==161.25){
      angE1=std::make_pair(ntheta-7.5,nphi-7.50);  angE2=std::make_pair(ntheta-7.5,nphi+15.0);
      angE3=std::make_pair(ntheta-7.5,nphi);       angE4=std::make_pair(ntheta-7.5,nphi+7.50);
      // only phi angle of chan16 changes
      angE9=std::make_pair(ntheta, nphi+15);       angE10=std::make_pair(ntheta,nphi-15);
    }
    if(ntheta==18.75){// from chan16 -> chan15: Front
      angE5=std::make_pair(ntheta+7.5,nphi-7.50);  angE6=std::make_pair(ntheta+7.5,nphi+15.0);
      angE7=std::make_pair(ntheta+7.5,nphi);       angE8=std::make_pair(ntheta+7.5,nphi+7.5);
      // only phi angle of chan16 changes
      angE9=std::make_pair(ntheta, nphi+15);       angE10=std::make_pair(ntheta,nphi-15);
    }
    // handeling of going from chan15->chan16: Back
    if((ntheta+7.5)==161.25){
      angE1=std::make_pair(ntheta+7.5,nphi-7.50);  angE2=std::make_pair(ntheta+7.5,nphi+7.50);
      angE3=std::make_pair(ntheta+7.5,nphi);       //angE4=std::make_pair(ntheta+7.5,nphi+3.75);
    }
    // from chan15 -> chan16: Front
    if((ntheta-7.5)==18.75){
      angE5=std::make_pair(ntheta-7.5,nphi-7.50);  angE6=std::make_pair(ntheta-7.5,nphi+7.50);
      angE7=std::make_pair(ntheta-7.5,nphi);       //angE8=std::make_pair(ntheta-7.5,nphi+3.75);
    }
    /******************************************************** 
     * 
     * Dealing with the revolution case
     *
     * *******************************************************/
    std::pair<double,double> angR1, angR2, angR3, angR4, angR5, angR6;
    if(nphi==3.75){
      // only phi angle of chan16 changes
      angR1=std::make_pair(ntheta, 348.75);       angR2=std::make_pair(ntheta,356.25);
      angR3=std::make_pair(ntheta-7.5, 348.75);   angR4=std::make_pair(ntheta+7.5,356.25);
      angR5=std::make_pair(ntheta-7.5, 356.25);   angR6=std::make_pair(ntheta+7.5,348.75);
    }
    if(nphi==348.75 || nphi==356.25){
      // only phi angle of chan16 changes
      angR1=std::make_pair(ntheta, 3.75);       angR2=std::make_pair(ntheta+7.5,3.75);
      angR3=std::make_pair(ntheta-7.5, 3.75);   //angR4=std::make_pair(ntheta+7.5,3.75);
    }
    // need to set critial variables to zero at the start of each iteration
    clusCrys=0;
    Eclus=0.;
    clusZ=0.; clusR=0.;
    thetaE=0; phiE=0;
    z_w=0; r_w=0;
    //rtheta=0; rphi=0;
    // For backward double counting elimination
    if(!csiCheck[tppair]){
      std::cout<<" Already checked this Crystal... moving on \n";
      goto checkedCrys;
    }
    if(csiClust[angP1] > 0 &&  csiClust[angP2] > 0 &&  csiClust[angP3] > 0 && csiClust[angP4] > 0 &&
      csiClust[angP5] > 0 && csiClust[angP6] > 0 && csiClust[angP7] > 0 &&
      csiClust[angP8] > 0){
      std::cout<<"  Total number of cluster crystals is 8 \n";
    }else if(csiClust[angP1]>0 || csiClust[angP2]>0 || csiClust[angP3]>0 || csiClust[angP4]>0 || 
      csiClust[angP5]>0 || csiClust[angP6]>0 || csiClust[angP7]>0 || csiClust[angP8]>0 ||
      csiClust[angE1]>0 || csiClust[angE2]>0 || csiClust[angE3]>0 || csiClust[angE4]>0 ||
      csiClust[angE5]>0 || csiClust[angE6]>0 || csiClust[angE7]>0 || csiClust[angE8]>0 ||
      csiClust[angE9]>0 || csiClust[angE10]>0 || csiClust[angR1]>0 || csiClust[angR2]>0 ||
      csiClust[angR3]>0 || csiClust[angR4]>0 || csiClust[angR5]>0 || csiClust[angR6]>0){
      //clusCrys=clusCrys+1;
      if(csiClust[angP1]>0){
        std::cout<<" This crystal Cluster pulse-height P1: "<<csiClust[angP1];
        std::cout<<" ["<<std::get<0>(angP1)<<", "<<std::get<1>(angP1)<<"] \n";
        if(csiCheck[angP1]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP1];
          thetaE=thetaE+csiClust[angP1]*(std::get<0>(angP1));
          phiE  =phiE  +csiClust[angP1]*(std::get<1>(angP1));
          clusZ=clusZ  +csiClust[angP1]*csiZ[angP1];
          clusR=clusR  +csiClust[angP1]*csiR[angP1];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP1]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP2]>0){
        std::cout<<" This crystal Cluster finder in pair loop P2: "<<csiClust[angP2];
        std::cout<<" ["<<std::get<0>(angP2)<<", "<<std::get<1>(angP2)<<"] \n";
        if(csiCheck[angP2]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP2];
          thetaE=thetaE+csiClust[angP2]*(std::get<0>(angP2));
          phiE  =phiE  +csiClust[angP2]*(std::get<1>(angP2));
          clusZ=clusZ  +csiClust[angP2]*csiZ[angP2];
          clusR=clusR  +csiClust[angP2]*csiR[angP2];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP2]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP3]>0){
        std::cout<<" This crystal Cluster finder in pair loop P3: "<<csiClust[angP3];
        std::cout<<" ["<<std::get<0>(angP3)<<", "<<std::get<1>(angP3)<<"] \n";
        if(csiCheck[angP3]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP3];
          thetaE=thetaE+csiClust[angP3]*(std::get<0>(angP3));
          phiE  =phiE  +csiClust[angP3]*(std::get<1>(angP3));
          clusZ=clusZ  +csiClust[angP3]*csiZ[angP3];
          clusR=clusR  +csiClust[angP3]*csiR[angP3];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP3]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP4]>0){
        std::cout<<" This crystal Cluster finder in pair loop P4: "<<csiClust[angP4];
        std::cout<<" ["<<std::get<0>(angP4)<<", "<<std::get<1>(angP4)<<"] \n";
        if(csiCheck[angP4]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP4];
          thetaE=thetaE+csiClust[angP4]*(std::get<0>(angP4));
          phiE  =phiE  +csiClust[angP4]*(std::get<1>(angP4));
          clusZ=clusZ  +csiClust[angP4]*csiZ[angP4];
          clusR=clusR  +csiClust[angP4]*csiR[angP4];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP4]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP5]>0){
        std::cout<<" This crystal Cluster finder in pair loop P5: "<<csiClust[angP5];
        std::cout<<" ["<<std::get<0>(angP5)<<", "<<std::get<1>(angP5)<<"] \n";
        if(csiCheck[angP5]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP5];
          thetaE=thetaE+csiClust[angP5]*(std::get<0>(angP5));
          phiE  =phiE  +csiClust[angP5]*(std::get<1>(angP5));
          clusZ=clusZ  +csiClust[angP5]*csiZ[angP5];
          clusR=clusR  +csiClust[angP5]*csiR[angP5];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP5]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP6]>0){
        std::cout<<" This crystal Cluster finder in pair loop P6: "<<csiClust[angP6];
        std::cout<<" ["<<std::get<0>(angP6)<<", "<<std::get<1>(angP6)<<"] \n";
        if(csiCheck[angP6]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP6];
          thetaE=thetaE+csiClust[angP6]*(std::get<0>(angP6));
          phiE  =phiE  +csiClust[angP6]*(std::get<1>(angP6));
          clusZ=clusZ  +csiClust[angP6]*csiZ[angP6];
          clusR=clusR  +csiClust[angP6]*csiR[angP6];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP6]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP7]>0){
        std::cout<<" This crystal Cluster finder in pair loop P7: "<<csiClust[angP7];
        std::cout<<" ["<<std::get<0>(angP7)<<", "<<std::get<1>(angP7)<<"] \n";
        if(csiCheck[angP7]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP7];
          thetaE=thetaE+csiClust[angP7]*(std::get<0>(angP7));
          phiE  =phiE  +csiClust[angP7]*(std::get<1>(angP7));
          clusZ=clusZ  +csiClust[angP7]*csiZ[angP7];
          clusR=clusR  +csiClust[angP7]*csiR[angP7];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP7]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angP8]>0){
        std::cout<<" This crystal Cluster finder in pair loop P8: "<<csiClust[angP8];
        std::cout<<" ["<<std::get<0>(angP8)<<", "<<std::get<1>(angP8)<<"] \n";
        if(csiCheck[angP8]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angP8];
          thetaE=thetaE+csiClust[angP8]*(std::get<0>(angP8));
          phiE  =phiE  +csiClust[angP8]*(std::get<1>(angP8));
          clusZ=clusZ  +csiClust[angP8]*csiZ[angP8];
          clusR=clusR  +csiClust[angP8]*csiR[angP8];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angP8]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE1]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E1: "<<csiClust[angE1];
        std::cout<<" ["<<std::get<0>(angE1)<<", "<<std::get<1>(angE1)<<"] \n";
        if(csiCheck[angE1]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE1];
          thetaE=thetaE+csiClust[angE1]*(std::get<0>(angE1));
          phiE  =phiE  +csiClust[angE1]*(std::get<1>(angE1));
          clusZ=clusZ  +csiClust[angE1]*csiZ[angE1];
          clusR=clusR  +csiClust[angE1]*csiR[angE1];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE1]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE2]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E2: "<<csiClust[angE2];
        std::cout<<" ["<<std::get<0>(angE2)<<", "<<std::get<1>(angE2)<<"] \n";
        if(csiCheck[angE2]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE2];
          thetaE=thetaE+csiClust[angE2]*(std::get<0>(angE2));
          phiE  =phiE  +csiClust[angE2]*(std::get<1>(angE2));
          clusZ=clusZ  +csiClust[angE2]*csiZ[angE2];
          clusR=clusR  +csiClust[angE2]*csiR[angE2];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE2]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE3]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E3: "<<csiClust[angE3];
        std::cout<<" ["<<std::get<0>(angE3)<<", "<<std::get<1>(angE3)<<"] \n";
        if(csiCheck[angE3]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE3];
          thetaE=thetaE+csiClust[angE3]*(std::get<0>(angE3));
          phiE  =phiE  +csiClust[angE3]*(std::get<1>(angE3));
          clusZ=clusZ  +csiClust[angE3]*csiZ[angE3];
          clusR=clusR  +csiClust[angE3]*csiR[angE3];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE3]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      // Got rid of angE4 and angE8 pairs because they are not needed at this time
      // Dongwi: 07.14.2019
      if(csiClust[angE5]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E5: "<<csiClust[angE5];
        std::cout<<" ["<<std::get<0>(angE5)<<", "<<std::get<1>(angE5)<<"] \n";
        if(csiCheck[angE5]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE5];
          thetaE=thetaE+csiClust[angE5]*(std::get<0>(angE5));
          phiE  =phiE  +csiClust[angE5]*(std::get<1>(angE5));
          clusZ=clusZ  +csiClust[angE5]*csiZ[angE5];
          clusR=clusR  +csiClust[angE5]*csiR[angE5];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE5]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE6]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E6: "<<csiClust[angE6];
        std::cout<<" ["<<std::get<0>(angE6)<<", "<<std::get<1>(angE6)<<"] \n";
        if(csiCheck[angE6]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE6];
          thetaE=thetaE+csiClust[angE6]*(std::get<0>(angE6));
          phiE  =phiE  +csiClust[angE6]*(std::get<1>(angE6));
          clusZ=clusZ  +csiClust[angE6]*csiZ[angE6];
          clusR=clusR  +csiClust[angE6]*csiR[angE6];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE6]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE7]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E7: "<<csiClust[angE7];
        std::cout<<" ["<<std::get<0>(angE7)<<", "<<std::get<1>(angE7)<<"] \n";
        if(csiCheck[angE7]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE7];
          thetaE=thetaE+csiClust[angE7]*(std::get<0>(angE7));
          phiE  =phiE  +csiClust[angE7]*(std::get<1>(angE7));
          clusZ=clusZ  +csiClust[angE7]*csiZ[angE7];
          clusR=clusR  +csiClust[angE7]*csiR[angE7];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE7]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE9]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E9: "<<csiClust[angE9];
        std::cout<<" ["<<std::get<0>(angE9)<<", "<<std::get<1>(angE9)<<"] \n";
        if(csiCheck[angE9]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE9];
          thetaE=thetaE+csiClust[angE9]*(std::get<0>(angE9));
          phiE  =phiE  +csiClust[angE9]*(std::get<1>(angE9));
          clusZ=clusZ  +csiClust[angE9]*csiZ[angE9];
          clusR=clusR  +csiClust[angE9]*csiR[angE9];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE9]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angE10]>0){
        std::cout<<" Edge effect Cluster finder in pair loop E10: "<<csiClust[angE10];
        std::cout<<" ["<<std::get<0>(angE10)<<", "<<std::get<1>(angE10)<<"] \n";
        if(csiCheck[angE10]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angE10];
          thetaE=thetaE+csiClust[angE10]*(std::get<0>(angE10));
          phiE  =phiE  +csiClust[angE10]*(std::get<1>(angE10));
          clusZ=clusZ  +csiClust[angE10]*csiZ[angE10];
          clusR=clusR  +csiClust[angE10]*csiR[angE10];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angE10]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      // revolution around CsI(Tl) barrel
      if(csiClust[angR1]>0){
        std::cout<<" Revolution case Cluster pulse-height R1: "<<csiClust[angR1];
        std::cout<<" ["<<std::get<0>(angR1)<<", "<<std::get<1>(angR1)<<"] \n";
        if(csiCheck[angR1]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR1];
          thetaE=thetaE+csiClust[angR1]*(std::get<0>(angR1));
          phiE  =phiE  +csiClust[angR1]*(std::get<1>(angR1));
          clusZ=clusZ  +csiClust[angR1]*csiZ[angR1];
          clusR=clusR  +csiClust[angR1]*csiR[angR1];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR1]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angR2]>0){
        std::cout<<" Revolution case Cluster finder in pair loop R2: "<<csiClust[angR2];
        std::cout<<" ["<<std::get<0>(angR2)<<", "<<std::get<1>(angR2)<<"] \n";
        if(csiCheck[angR2]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR2];
          thetaE=thetaE+csiClust[angR2]*(std::get<0>(angR2));
          phiE  =phiE  +csiClust[angR2]*(std::get<1>(angR2));
          clusZ=clusZ  +csiClust[angR2]*csiZ[angR2];
          clusR=clusR  +csiClust[angR2]*csiR[angR2];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR2]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angR3]>0){
        std::cout<<" Edge effect Cluster finder in pair loop R3: "<<csiClust[angR3];
        std::cout<<" ["<<std::get<0>(angR3)<<", "<<std::get<1>(angR3)<<"] \n";
        if(csiCheck[angR3]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR3];
          thetaE=thetaE+csiClust[angR3]*(std::get<0>(angR3));
          phiE  =phiE  +csiClust[angR3]*(std::get<1>(angR3));
          clusZ=clusZ  +csiClust[angR3]*csiZ[angR3];
          clusR=clusR  +csiClust[angR3]*csiR[angR3];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR3]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angR4]>0){
        std::cout<<" Edge effect Cluster finder in pair loop R4: "<<csiClust[angR4];
        std::cout<<" ["<<std::get<0>(angR4)<<", "<<std::get<1>(angR4)<<"] \n";
        if(csiCheck[angR4]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR4];
          thetaE=thetaE+csiClust[angR4]*(std::get<0>(angR4));
          phiE  =phiE  +csiClust[angR4]*(std::get<1>(angR4));
          clusZ=clusZ  +csiClust[angR4]*csiZ[angR4];
          clusR=clusR  +csiClust[angR4]*csiR[angR4];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR4]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angR5]>0){
        std::cout<<" Edge effect Cluster finder in pair loop R5: "<<csiClust[angR5];
        std::cout<<" ["<<std::get<0>(angR5)<<", "<<std::get<1>(angR5)<<"] \n";
        if(csiCheck[angR5]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR5];
          thetaE=thetaE+csiClust[angR5]*(std::get<0>(angR5));
          phiE  =phiE  +csiClust[angR5]*(std::get<1>(angR5));
          clusZ=clusZ  +csiClust[angR5]*csiZ[angR5];
          clusR=clusR  +csiClust[angR5]*csiR[angR5];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR5]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      if(csiClust[angR6]>0){
        std::cout<<" Edge effect Cluster finder in pair loop R6: "<<csiClust[angR6];
        std::cout<<" ["<<std::get<0>(angR6)<<", "<<std::get<1>(angR6)<<"] \n";
        if(csiCheck[angR6]){
          clusCrys=clusCrys+1;
          Eclus=Eclus+csiClust[angR6];
          thetaE=thetaE+csiClust[angR6]*(std::get<0>(angR6));
          phiE  =phiE  +csiClust[angR6]*(std::get<1>(angR6));
          clusZ=clusZ  +csiClust[angR6]*csiZ[angR6];
          clusR=clusR  +csiClust[angR6]*csiR[angR6];
          std::cout<<" This crystal is now removed from the list: "<<std::endl;
          csiCheck[angR6]=false;
        }else{
          std::cout<<" Already checked this crystal \n";
        }
      }
      std::cout<<" some crystals actually have hits "<<clusCrys<<std::endl;
    }else{
      clusCrys=0;//clusCrys+1;
      std::cout<<" Single cluster crystals here \n"; //<<clusCrys<<std::endl;
    }
    if(csiCheck[tppair]){
      clusCrys=clusCrys+1;
      Eclus=Eclus+csiClust[tppair];
      thetaE=thetaE+csiClust[tppair]*(std::get<0>(tppair));
      phiE  =phiE  +csiClust[tppair]*(std::get<1>(tppair));
      clusZ=clusZ  +csiClust[tppair]*csiZ[tppair];
      clusR=clusR  +csiClust[tppair]*csiR[tppair];
      std::cout<<" -->  pulse-heignt for central crystal: "<<csiClust[tppair];
      std::cout<<" ["<<std::get<0>(tppair)<<", "<<std::get<1>(tppair)<<"] \n";
      std::cout<<" >>>  Cluster energy is ------------->: "<<Eclus<<" [GeV]";
    }
    if(clusCrys>=2){
      numOfClus++;
      // perform energy-weighting and convert from deg-->rad
      rtheta=TMath::DegToRad()*(thetaE/Eclus);
      rphi=TMath::DegToRad()*(phiE/Eclus);
      clusEne.push_back(Eclus);
      clusThetaE.push_back(rtheta);
      clusPhiE.push_back(rphi);
      z_w=clusZ/Eclus;
      r_w=clusR/Eclus;
      clusEz.push_back(z_w);
      clusEr.push_back(r_w);
      //h2ang->Fill(rtheta,rphi);
      //h2deg->Fill(rtheta*180./M_PI,rphi*180/M_PI);
    }
    if(clusCrys==1){
      numOfsingleClus++;
      rtheta=TMath::DegToRad()*(thetaE/Eclus);
      rphi=TMath::DegToRad()*(phiE/Eclus);
      singTheta.push_back(rtheta);
      singPhi.push_back(rphi);
      singleEne.push_back(Eclus);
      z_w=clusZ/Eclus;
      r_w=clusR/Eclus;
      singZ.push_back(z_w);
      singR.push_back(r_w);
      //h2ang->Fill(rtheta,rphi);
      //h2deg->Fill(rtheta*180./M_PI,rphi*180/M_PI);
      //h1sclus->Fill(1);
    }
    csiCheck[tppair]=false; // mute central crystal
    std::cout<<"  Number of crystals is:  "<<clusCrys<<std::endl;
    checkedCrys:
    std::cout<<"\n --------------------------------------------------------------------->\n\n";
  } // end of cluster finding routine
  csiCheck.clear();
  std::cout<<"   Checking the cluster size: "<<clusEne.size()<<", "<<singleEne.size()<<std::endl;
  std::cout<<"\n\n  Number of clusters is   :  "<<numOfClus<<std::endl;
  std::cout<<"  Number of single clusters is:  "<<numOfsingleClus<<std::endl;
  std::cout<<" ***************************************************************************\n";
}
