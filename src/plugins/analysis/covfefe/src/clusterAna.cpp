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
  tof1ID=new TH1D("tofID","stat",13,-0.5,12.5);
  h1score=new TH1D("h1score","stat",16.,-0.5,15.5);
  h1Pkpi2=new TH1D("kpi2Mom","stat",55.,0.15,.245);
  h1prmAng=new TH1D("prmAng","stat",75.0,-1.1,1.1);
  gStyle->SetOptStat(0);
  clustEval=new evalClusters();
  clustEval->setChannel(7); // must be called before defHistos (defaul chan7)
  clustEval->defHistos();

  return 0;
}

Long_t covfefe::startup_clust(){
  getBranchObject("treeClus",(TObject **) &clsmar); 
  getBranchObject("tracks",(TObject **) &trackArr); 
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
      if(trackArr->vectTrack.size()>0){
        pVert=trackArr->vectTrack[0].getVertSPPionPlus()/1000.;
        piPpx=pVert*trackArr->vectTrack[0].getVertSNx();
        piPpy=pVert*trackArr->vectTrack[0].getVertSNy();
        piPpz=pVert*trackArr->vectTrack[0].getVertSNz();
        P=std::sqrt(piPpx*piPpx+piPpy*piPpy+piPpz*piPpz);
      }
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
      // eval only 2 clusters
      if(clsmar->clusterM==2 && clsmar->E_prim2>0.02){
        prim1vec3.SetXYZ(piPpx,piPpy,piPpz);
        prim2vec3.SetXYZ(clsmar->prim2px,clsmar->prim2py,clsmar->prim2pz);
	primCosTheta=std::cos(prim1vec3.Angle(prim2vec3));
        clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,primCosTheta,0);
        clustEval->fillMltp(clsmar->extraTOF1_size);
        clustEval->fillHistos(prim2lv.M2(),0);
	//std::cout<<" Checking the momentum of pi+ =>"<<P<<"\n";
	if(clsmar->fVertSP>=19. && clsmar->fVertSP<=215.){
          clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,primCosTheta,1);
          clustEval->fillHistos(prim2lv.M2(),1);
	  // check that TOF1 multiplicity is greater than or equal to 2
	  clustEval->fillHistos(clsmar->extraTOF1->size(),3);
	  for(UInt_t id=0; id<clsmar->vec_extraTOF1->size(); id++){
            // tof1ID
	    tof1ID->Fill((*(clsmar->vec_extraTOF1))[id][0]);
            // tof1 score
	    h1score->Fill((*(clsmar->vec_extraTOF1))[id][1]);
	  }
	  h1Pkpi2->Fill(pVert);
	  h1prmAng->Fill(primCosTheta);
	  std::cout<<" ... Mass2: "<<prim2lv.M2()<<std::endl;
	  std::cout<<" Px, Py, Pz: ["<<clsmar->prim2px<<", "<<clsmar->prim2py<<", "<<clsmar->prim2pz<<"]\n";
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
