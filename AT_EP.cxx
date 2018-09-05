#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TProfile.h>
#include "Analysis.h"
#include "AT_EP.h"

AT_EP::AT_EP() {
  fNpt = 22;//14; //9;
  //float ptbins[15] = {1.0, 1.1, 1.2, 1.3, 1.5,
  //		      1.7, 1.9, 2.2, 2.5, 3.0,
  //		      4.0, 6.0, 10., 15., 20.};
  float ptbins[23] = {1.0, 1.1, 1.2, 1.3, 1.4,
		      1.5, 1.6, 1.7, 1.8, 2.0,
		      2.2, 2.4, 2.6, 2.8, 3.0,
		      3.4, 3.8, 4.4, 5.0, 6.0,
		      8.0, 10., 20.};
  for(int i=0; i!=fNpt+1; ++i) fPtBins[i] = ptbins[i];
  fNma = 15;
  float mabins[16] = {0.040, 0.060, 0.080, 0.100, 0.110,
		      0.120, 0.130, 0.140, 0.150, 0.160,
		      0.170, 0.180, 0.200, 0.220, 0.240, 0.260};
  for(int i=0; i!=fNma+1; ++i) fMassBins[i] = mabins[i];
}

AT_EP::~AT_EP() {
}

void AT_EP::Init() {
  Analysis *ana = Analysis::Instance();
  fCandidates = ana->GetCandidates();
  fCandidates2 = ana->GetCandidates2();
  for(int i=0; i!=4; ++i) fQ[i] = ana->GetQ(i);
  TProfile *fPQ_Cos[4];
  hEta = new TH1F("hEta","hEta",100,-1,+1);
  hEta2 = new TH1F("hEta2","hEta2",100,-1,+1);
  for(int p=0; p!=fNpt; ++p) {
    hMass[p] = new TH1F( Form("hMass_PB%d",p),
			 Form("hMass_PB%d;Mass",p),
			 240,//220, //assuming 0.260-0.040 = 0.220
			 fMassBins[0],fMassBins[fNma] );
    hMass2[p] = new TH1F( Form("hMass2_PB%d",p),
			  Form("hMass2_PB%d;Mass",p),
			  240,//220, //assuming 0.260-0.040 = 0.220
			  fMassBins[0],fMassBins[fNma] );
    for(int n=0; n!=5; ++n) {
      hCos[n][p] = new TProfile( Form("hCos%dDP_PB%d",n,p),
				 Form("hCos%dDP_PB%d;Mass",n,p), 
				 120,//110,
				 fMassBins[0],fMassBins[fNma] );
      hCos2[n][p] = new TProfile( Form("hCos2%dDP_PB%d",n,p),
				  Form("hCos2%dDP_PB%d;Mass",n,p), 
				  120,//110,
				  fMassBins[0],fMassBins[fNma] );
      //				 fNma, fMassBins );
    }
  }

  for(int n=0; n!=4; ++n) {
    hPsi[n] = new TH1F( Form("hPSI%d",n),
			Form("hPSI%d;Mass",n),
			120, -TMath::TwoPi(), +TMath::TwoPi() );
  }

}

void AT_EP::Finish() {
  hEta->Write();
  hEta2->Write();
  for(int p=0; p!=fNpt; ++p) {
    hMass[p]->Write();
    hMass2[p]->Write();
    for(int n=0; n!=5; ++n) {
      hCos[n][p]->Write();
      hCos2[n][p]->Write();
    }
  }
  for(int n=0; n!=4; ++n) {
    hPsi[n]->Write();
  }
}

void AT_EP::Exec() {
  if(fQ[0]->M()<1) return;

  for(int ord=0; ord!=4; ++ord) {
    hPsi[ord]->Fill( fQ[ord]->Psi2Pi() );
  }

  // CANDIDATES 1
  uint npa = fCandidates->size();
  for(int i=0; i!=npa; ++i) {
    TLorentzVector a = fCandidates->at(i);
    double ma = a.M();
    double pt = a.Pt();
    int mb = BinMass( ma );
    int pb = BinPt( pt );
    if(mb<0||pb<0) continue;
    /// recording
    hEta->Fill( a.Eta() );
    hMass[pb]->Fill(ma);
    for(int ord=0; ord!=4; ++ord) {
      int nn = ord+1;
      double dphi = a.Phi() - fQ[ord]->Psi2Pi();
      double cos = TMath::Cos( nn*dphi );
      hCos[ord][pb]->Fill(ma,cos);
      if(nn==4) {
	double dphi = a.Phi() - fQ[1]->Psi2Pi();
	double cos = TMath::Cos( 4*dphi );
	hCos[4][pb]->Fill(ma,cos);
      }
    }
  }

  // CANDIDATES 2
  npa = fCandidates2->size();
  for(int i=0; i!=npa; ++i) {
    TLorentzVector a = fCandidates2->at(i);
    double ma = a.M();
    double pt = a.Pt();
    int mb = BinMass( ma );
    int pb = BinPt( pt );
    if(mb<0||pb<0) continue;
    /// recording
    hEta2->Fill( a.Eta() );
    hMass2[pb]->Fill(ma);
    for(int ord=0; ord!=4; ++ord) {
      int nn = ord+1;
      double dphi = a.Phi() - fQ[ord]->Psi2Pi();
      double cos = TMath::Cos( nn*dphi );
      hCos2[ord][pb]->Fill(ma,cos);
      if(nn==4) {
	double dphi = a.Phi() - fQ[1]->Psi2Pi();
	double cos = TMath::Cos( 4*dphi );
	hCos2[4][pb]->Fill(ma,cos);
      }
    }
  }
}

int AT_EP::BinPt(float pt) {
  int ret = -1;
  if( pt<fPtBins[0] || pt>fPtBins[fNpt] ) return ret;
  for(int p=0; p!=fNpt+1; ++p) 
    if(pt>fPtBins[p]) ret = p;
  return ret;
}

int AT_EP::BinMass(float ma) {
  int ret = -1;
  if( ma<fMassBins[0] || ma>fMassBins[fNma] ) return ret;
  for(int m=0; m!=fNma+1; ++m) 
    if(ma>fMassBins[m]) ret = m;
  return ret;
}
