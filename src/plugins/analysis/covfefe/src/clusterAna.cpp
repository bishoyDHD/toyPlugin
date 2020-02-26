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
  std::string pname[]={"pNoCut","pKpi2Cut","ptof1Cut","pAngCut"};
  std::string Ename[]={"ENoCut","EKpi2Cut","Etof1Cut","EAngCut"};
  std::string Mname[]={"MNoCut","MKpi2Cut","Mtof1Cut","MAngCut"};
  std::string tf1IDname[]={"tof1IDNoCut","tof1IDKpi2Cut","tof1IDtof1Cut","tof1IDAngCut"};
  std::string tf1Mname[]={"tof1MNoCut","tof1MKpi2Cut","tof1Mtof1Cut","tof1AngCut"};
  std::string tf1Sname[]={"tof1SNoCut","tof1SKpi2Cut","tof1Stof1Cut","tof1AngCut"};
  std::string tf1vMname[]={"tof1vMNoCut","tof1vMKpi2Cut","tof1vMtof1Cut","tof1vMAngCut"};
  std::string prAngname[]={"prAngNoCut","prAngKpi2Cut","prAngtof1Cut","prAngAngCut"};
  std::string clustAngname[]={"clustAngNoCut","clustAngKpi2Cut","clustAngtof1Cut","clustAngAngCut"};
  std::string crysNname[]={"crysNumNoCut","crysNumKpi2Cut","crysNumtof1Cut","crysNumAngCut"};
  for(int i=0; i<4; i++){
    angClust[i]=new TH1D(clustAngname[i].c_str(),"stats",75,-1.1,1.1);
    prAng[i]=new TH1D(prAngname[i].c_str(),"stats",75,-1.1,1.1);
    h1Mass[i]=new TH1D(Mname[i].c_str(),"stats",84.5,0.0,.210);
    h1Eclust[i]=new TH1D(Ename[i].c_str(),"stats",100.,0.0,.30);
    h1P[i]=new TH1D(pname[i].c_str(),"stats",84.5,0.1,.300);
    tof1ID[i]=new TH1D(tf1IDname[i].c_str(),"stat",13,-0.5,12.5);
    tof1M[i]=new TH1D(tf1Mname[i].c_str(),"stat",13,-0.5,12.5);
    tof1S[i]=new TH1D(tf1Sname[i].c_str(),"stat",20,-0.5,19.5);
    //clustM[i]=new TH1D(crysNname[i].c_str(),"stat",13,-0.5,12.5);
    tof1IDvM[i]=new TH2D(tf1vMname[i].c_str(),"stat",13,-.5,12.5,20,-.5,19.5);
  }
  h1TofCorr=new TH1D("angCorr","stat",90.,-3.75,362.25);
  h1diffTof2=new TH1D("tof2diff","stat",25,-12.5,12.5);
  h1score=new TH1D("h1score","stat",16.,-0.5,15.5);
  h1Pkpi2=new TH1D("GapP","stat",55.,0.15,.245);
  Eclust=new TH1D("Tof1Eclust","stat",100,0.,.3);
  InvM=new TH1D("Tof1InvMass","stat",84,0.,.210);
  h1prmAng=new TH1D("prmAng","stat",75.0,-1.1,1.1);
  gStyle->SetOptStat(0);
  clustEval=new evalClusters();
  clustEval->setChannel(7); // must be called before defHistos (defaul chan7)
  clustEval->defHistos();

  return 0;
}

Long_t covfefe::startup_clust(){
  getBranchObject("treeClus",(TObject **) &clsmar); 
  //getBranchObject("tracks",(TObject **) &trackArr); 
  //calibcsi=new CATCaliCsI();
  //makeBranch("marinCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::process_clust(){
  if(clsmar->E_prim2>0.045){
      // Construct CM Lorentz vectors for 2gammas
      //std::cout<<" making sure this thing works "<<clsmar->clusterM<<" "<<clsmar->E_prim2<<" "<<clsmar->Ncrys<<"\n";
      prim2lv.SetPxPyPzE(clsmar->prim2px,clsmar->prim2py,clsmar->prim2pz,clsmar->E_prim2);
      h1diffTof2->Fill(clsmar->fgapNumTof2-clsmar->tof2Gap);
      //double E=std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      /*
      double E=0.493677-std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      double bx=piPpx/E;
      double by=piPpy/E;
      double bz=piPpz/E;
      double b2=bx*bx+by*by+bz*bz;
      double gamma=1.0/std::sqrt(1.0-b2);*/
      //clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta,0);
      clustEval->fillHistos(prim2lv.M2(),0);
      //clustEval->fillHistos(clsmar->cpid1theta,clsmar->cpid1phi,clsmar->cpid2theta,clsmar->cpid2phi,0,0);
      // eval only 2 clusters
      if(clsmar->clusterM==2 && clsmar->E_prim2>0.045){
        prim1vec3.SetXYZ(piPpx,piPpy,piPpz);
        prim2vec3.SetXYZ(clsmar->prim2px,clsmar->prim2py,clsmar->prim2pz);
	primCosTheta=std::cos(prim1vec3.Angle(prim2vec3));
        clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,primCosTheta,0);
        clustEval->fillMltp(clsmar->extraTOF1_size);
        clustEval->fillHistos(prim2lv.M2(),0);
	// fill histos to study various cut conditions
	angClust[0]->Fill(clsmar->clCosTheta);
	prAng[0]->Fill(clsmar->prCosTheta);
	h1Mass[0]->Fill(clsmar->M_prim2);
	h1Eclust[0]->Fill(clsmar->E_prim2);
	h1P[0]->Fill(clsmar->fVertSP/1000.);
	//clustM[0]->Fill(clsmar->ClustCrys);
        // Fill various tof1 histograms
        tof1M[0]->Fill(clsmar->vec_extraTOF1->size());
        for(UInt_t id=0; id<clsmar->vec_extraTOF1->size(); id++){
          // tof1ID
          tof1ID[0]->Fill((*(clsmar->vec_extraTOF1))[id][0]);
          // tof1 score
          tof1S[0]->Fill((*(clsmar->vec_extraTOF1))[id][1]);
          tof1IDvM[0]->Fill((*(clsmar->vec_extraTOF1))[id][0],(*(clsmar->vec_extraTOF1))[id][1]);
        }
	//std::cout<<" Checking the momentum of pi+ =>"<<P<<"\n";
	if(clsmar->fVertSP>=190. && clsmar->fVertSP<=215.){
          clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta,1);
          clustEval->fillHistos(prim2lv.M2(),1);
	  // check that TOF1 multiplicity is greater than or equal to 2
	  clustEval->fillHistos(clsmar->extraTOF1->size(),3);
	  // Fill various tof1 histograms
          tof1M[1]->Fill(clsmar->vec_extraTOF1->size());
	  for(UInt_t id=0; id<clsmar->vec_extraTOF1->size(); id++){
            // tof1ID
	    tof1ID[1]->Fill((*(clsmar->vec_extraTOF1))[id][0]);
            // tof1 score
	    h1score->Fill((*(clsmar->vec_extraTOF1))[id][1]);
	    tof1S[1]->Fill((*(clsmar->vec_extraTOF1))[id][1]);
            tof1IDvM[1]->Fill((*(clsmar->vec_extraTOF1))[id][0],(*(clsmar->vec_extraTOF1))[id][1]);
	  }
          // fill histos to study various cut conditions
          angClust[1]->Fill(clsmar->clCosTheta);
          prAng[1]->Fill(clsmar->prCosTheta);
          h1Mass[1]->Fill(clsmar->M_prim2);
          h1Eclust[1]->Fill(clsmar->E_prim2);
          h1P[1]->Fill(clsmar->fVertSP/1000.);
	  h1prmAng->Fill(primCosTheta);
	  // apply opening angle cut pi+pi0
	  if(clsmar->prCosTheta<=-.8){
            // Fill various tof1 histograms
            tof1M[3]->Fill(clsmar->vec_extraTOF1->size());
            for(UInt_t id=0; id<clsmar->vec_extraTOF1->size(); id++){
              // tof1ID
              tof1ID[3]->Fill((*(clsmar->vec_extraTOF1))[id][0]);
              // tof1 score
              tof1S[3]->Fill((*(clsmar->vec_extraTOF1))[id][1]);
              tof1IDvM[3]->Fill((*(clsmar->vec_extraTOF1))[id][0],(*(clsmar->vec_extraTOF1))[id][1]);
            }
            // fill histos to study various cut conditions
            angClust[3]->Fill(clsmar->clCosTheta);
            prAng[3]->Fill(clsmar->prCosTheta);
            h1Mass[3]->Fill(clsmar->M_prim2);
            h1Eclust[3]->Fill(clsmar->E_prim2);
            h1P[3]->Fill(clsmar->fVertSP/1000.);
	  }
	  std::cout<<" ... Mass2: "<<prim2lv.M2()<<std::endl;
	  std::cout<<" Px, Py, Pz: ["<<clsmar->prim2px<<", "<<clsmar->prim2py<<", "<<clsmar->prim2pz<<"]\n";
	}
	if(clsmar->extraTOF1->size()>=2){
	  Eclust->Fill(clsmar->E_prim2);
	  InvM->Fill(clsmar->M_prim2);
	  h1Pkpi2->Fill(clsmar->fVertSPiplus);
          // fill histos to study various cut conditions
          angClust[2]->Fill(clsmar->clCosTheta);
          prAng[2]->Fill(clsmar->prCosTheta);
          h1Mass[2]->Fill(clsmar->M_prim2);
          h1Eclust[2]->Fill(clsmar->E_prim2);
          h1P[2]->Fill(clsmar->fVertSP/1000.);
	  // Fill various tof1 histograms
          tof1M[2]->Fill(clsmar->vec_extraTOF1->size());
	  for(UInt_t id=0; id<clsmar->vec_extraTOF1->size(); id++){
            // tof1ID
	    tof1ID[2]->Fill((*(clsmar->vec_extraTOF1))[id][0]);
            // tof1 score
	    tof1S[2]->Fill((*(clsmar->vec_extraTOF1))[id][1]);
            tof1IDvM[2]->Fill((*(clsmar->vec_extraTOF1))[id][0],(*(clsmar->vec_extraTOF1))[id][1]);
	  }
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
