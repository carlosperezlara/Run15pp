#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include "Analysis.h"
#include "AT_PIDFlow.h"

AT_PIDFlow::AT_PIDFlow() : AT_ReadTree() {
}

AT_PIDFlow::~AT_PIDFlow() {
  if(hPt) delete hPt;
  if(hNTrk) delete hNTrk;
  for(int i=0; i!=3; ++i)
    if(hPtDPhi[i]) delete hPtDPhi[i];
}

void AT_PIDFlow::MyInit() {
  hPt =   new TH1F("hPt",  "hPt", 100,0.0,5.0);
  hNTrk = new TH1F("hNTrk","hNTrk",100,-0.5,99.5);
  hPtDPhi[0] = new TH2F("hPtDPhi1",";pt;phi-Psi1", 100,0.0,5.0,100,-7,7);
  hPtDPhi[1] = new TH2F("hPtDPhi2",";pt;phi-Psi2", 100,0.0,5.0,100,-5,7);
  hPtDPhi[2] = new TH2F("hPtDPhi3",";pt;phi-Psi3", 100,0.0,5.0,100,-3,7);
  hPtDPhiME[0] = new TH2F("hPtDPhi1ME",";pt;phi-Psi1", 100,0.0,5.0,100,-7,7);
  hPtDPhiME[1] = new TH2F("hPtDPhi2ME",";pt;phi-Psi2", 100,0.0,5.0,100,-5,7);
  hPtDPhiME[2] = new TH2F("hPtDPhi3ME",";pt;phi-Psi3", 100,0.0,5.0,100,-3,7);
}

void AT_PIDFlow::MyFinish() {
  hPt->Write();
  hNTrk->Write();
  for(int i=0; i!=3; ++i) {
    hPtDPhi[i]->Write();
    hPtDPhiME[i]->Write();
  }
}

void AT_PIDFlow::MyExec() {
  float vtxz = fGLB.vtxZ;
  float cent = fGLB.cent;
  if(TMath::Abs(vtxz)>20) return;
  if(cent<0||cent>5) return;
  if(ReferenceTracks()<2) return;
  int ntrk=0;
  uint ntrks = pTRKpt->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float pt   = TMath::Abs( pTRKpt->at(itrk) );
    float pz   = pTRKpz->at(itrk);
    float phi  = pTRKphi->at(itrk);
    float ep   = pTRKecore->at(itrk) / TMath::Sqrt(pt*pt+pz*pz);
    float tof  = pTRKetof->at(itrk);
    float chi2 = pTRKchisq->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    float zed  = pTRKzed->at(itrk);
    int qua = pTRKqua->at(itrk);
    if(qua!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    ntrk++;
    hPt->Fill(pt);
    float dphi1 = phi - Psi1_BBC;
    float dphi2 = phi - Psi2_BBC;
    float dphi3 = phi - Psi3_BBC;
    hPtDPhi[0]->Fill(pt,dphi1);
    hPtDPhi[1]->Fill(pt,dphi2);
    hPtDPhi[2]->Fill(pt,dphi3);
    float dphi1me = phi - fPsi1_BBC_PE;
    float dphi2me = phi - fPsi2_BBC_PE;
    float dphi3me = phi - fPsi3_BBC_PE;
    hPtDPhiME[0]->Fill(pt,dphi1me);
    hPtDPhiME[1]->Fill(pt,dphi2me);
    hPtDPhiME[2]->Fill(pt,dphi3me);
  }
  hNTrk->Fill(ntrk);
  fPsi1_BBC_PE = Psi1_BBC;
  fPsi2_BBC_PE = Psi2_BBC;
  fPsi3_BBC_PE = Psi3_BBC;
}
