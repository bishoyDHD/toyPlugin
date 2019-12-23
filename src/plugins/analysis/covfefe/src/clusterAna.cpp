#include <covfefe.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "TVector2.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>

Long_t covfefe::hist_clust(){
  gStyle->SetOptStat(0);
  clustEval=new evalClusters();
  clustEval->setChannel(7); // must be called before defHistos (defaul chan7)
  clustEval->defHistos();

  return 0;
}

Long_t covfefe::startup_clust(){
  getBranchObject("treeClus",(TObject **) &clsmar); 
  //calibcsi=new CATCaliCsI();
  //makeBranch("marinCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::process_clust(){
  if(clsmar->E_prim2>0){
      // Construct CM Lorentz vectors for 2gammas
      //std::cout<<" making sure this thing works "<<clsmar->clusterM<<" "<<clsmar->E_prim2<<" "<<clsmar->Ncrys<<"\n";
      prim2lv.SetPxPyPzE(clsmar->prim2px,clsmar->prim2py,clsmar->prim2pz,clsmar->E_prim2);
      double piPpx=-1*clsmar->prim1px;
      double piPpy=-1*clsmar->prim1py;
      double piPpz=-1*clsmar->prim1pz;
      double P=std::sqrt(piPpx*piPpx+piPpy*piPpy+piPpz*piPpz);
      //double E=std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      double E=0.493677-std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      double bx=piPpx/E;
      double by=piPpy/E;
      double bz=piPpz/E;
      double b2=bx*bx+by*by+bz*bz;
      double gamma=1.0/std::sqrt(1.0-b2);
      //clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta,0);
      clustEval->fillHistos(prim2lv.M2(),0);
      //clustEval->fillHistos(clsmar->cpid1theta,clsmar->cpid1phi,clsmar->cpid2theta,clsmar->cpid2phi,0,0);
      clustEval->fillMltp(clsmar->clusterM);
      // eval only 2 clusters
      if(clsmar->clusterM==2 && clsmar->E_prim2>0.02){
        clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta,0);
        clustEval->fillHistos(prim2lv.M2(),0);
	//std::cout<<" Checking the momentum of pi+ =>"<<P<<"\n";
	// check that TOF1 multiplicity is greater than or equal to 2
	if(clsmar->extraTOF1>=2){
          clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta,1);
          clustEval->fillHistos(prim2lv.M2(),1);
	  //std::cout<<" ... Mass2: "<<prim2lv.M2()<<std::endl;
	  //std::cout<<" Px, Py, Pz: ["<<clsmar->prim2px<<", "<<clsmar->prim2py<<", "<<clsmar->prim2pz<<"]\n";
	}
      }
    //}
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize_clust(){
  clustEval->drawHistos();
  clustEval->singleCanvas();

  return 0; // 0 = all ok
};
