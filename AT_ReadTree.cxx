#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TH2F.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "Analysis.h"
#include "AT_ReadTree.h"

AT_ReadTree::AT_ReadTree() : AnalysisTask() {
  // -30.0 ==> +30.0 (60+1)
  // 0.5 ==> 60.5 (50+1)
  fNBinsVtx = 40;
  fNBinsCen = 60;
  fMinBinVtx = -20.0;
  fMinBinCen = 0.5;
  for(int bce=0; bce!=fNBinsCen; ++bce) {
    for(int bvt=0; bvt!=fNBinsVtx; ++bvt) {
      for(int ord=0; ord!=6; ++ord) {
	for(int se=0; se!=2; ++se) {
	  bbcm[se][ord][0][bce][bvt] = 0.0;
	  bbcm[se][ord][1][bce][bvt] = 0.0;
	}
      }
      for(int ord=0; ord!=4; ++ord) {
	for(int ior=0; ior!=32; ++ior) {
	  bbcc[ior][ord][bce][bvt] = 0.0;
	  bbcs[ior][ord][bce][bvt] = 0.0;
	}
      }
    }
  }
}

void AT_ReadTree::Init() {
  Analysis *ana = Analysis::Instance();
  TTree *tree = ana->GetTree();
  if(!tree) {
    std::cout << "AT_ReadTree:Init says: Tree not found." << std::endl;
    return;
  }
  //Opening assigning branches
  tree->SetBranchAddress("Event",&fGLB);
  //=
  tree->SetBranchAddress("Q1ex",&pQ1ex);
  tree->SetBranchAddress("Q2ex",&pQ2ex);
  tree->SetBranchAddress("Q3ex",&pQ3ex);
  tree->SetBranchAddress("Q4ex",&pQ4ex);
  tree->SetBranchAddress("Q6ex",&pQ6ex);
  tree->SetBranchAddress("Q8ex",&pQ8ex);
  tree->SetBranchAddress("Q1fv",&pQ1fv);
  tree->SetBranchAddress("Q2fv",&pQ2fv);
  tree->SetBranchAddress("Q3fv",&pQ3fv);
  tree->SetBranchAddress("Q1bb",&pQ1bb);
  tree->SetBranchAddress("Q2bb",&pQ2bb);
  tree->SetBranchAddress("Q3bb",&pQ3bb);
  tree->SetBranchAddress("Q4bb",&pQ4bb);
  tree->SetBranchAddress("Q6bb",&pQ6bb);
  tree->SetBranchAddress("Q8bb",&pQ8bb);
  //=
  tree->SetBranchAddress("EMCid",   &pEMCid);
  tree->SetBranchAddress("EMCtwrid",&pEMCtwrid);
  tree->SetBranchAddress("EMCx",    &pEMCx);
  tree->SetBranchAddress("EMCy",    &pEMCy);
  tree->SetBranchAddress("EMCz",    &pEMCz);
  tree->SetBranchAddress("EMCecore",&pEMCecore);
  tree->SetBranchAddress("EMCecent",&pEMCecent);
  tree->SetBranchAddress("EMCchisq",&pEMCchisq);
  tree->SetBranchAddress("EMCtimef",&pEMCtimef);
  //=
  tree->SetBranchAddress("TRKqua",  &pTRKqua);
  tree->SetBranchAddress("TRKpt",   &pTRKpt);
  tree->SetBranchAddress("TRKphi",  &pTRKphi);
  tree->SetBranchAddress("TRKpz",   &pTRKpz);
  tree->SetBranchAddress("TRKecore",&pTRKecore);
  tree->SetBranchAddress("TRKetof", &pTRKetof);
  tree->SetBranchAddress("TRKplemc",&pTRKplemc);
  tree->SetBranchAddress("TRKtwrid",&pTRKtwrid);
  tree->SetBranchAddress("TRKchisq",&pTRKchisq);
  tree->SetBranchAddress("TRKdphi", &pTRKdphi);
  tree->SetBranchAddress("TRKdz",   &pTRKdz);
  tree->SetBranchAddress("TRKpc3sdphi",&pTRKpc3sdphi);
  tree->SetBranchAddress("TRKpc3sdz",  &pTRKpc3sdz);
  tree->SetBranchAddress("TRKzed",  &pTRKzed);
  tree->SetBranchAddress("TRKdisp", &pTRKdisp);
  tree->SetBranchAddress("TRKprob", &pTRKprob);
  tree->SetBranchAddress("TRKcid",  &pTRKcid);
  //=
  tree->SetBranchAddress("MXSpt",  &pMXSpt);
  tree->SetBranchAddress("MXSpz",  &pMXSpz);
  tree->SetBranchAddress("MXSphi", &pMXSphi);
  tree->SetBranchAddress("MXSflyr",&pMXSflyr);
  tree->SetBranchAddress("MXSsingleD", &pMXSsingleD);
  tree->SetBranchAddress("MXSsingleP", &pMXSsingleP);
  tree->SetBranchAddress("MXSempccent",&pMXSempccent);
  tree->SetBranchAddress("MXSempc3x3", &pMXSempc3x3);

  LoadTableEP();

  hEvents = new TH2F("Events",";vtx;cent",40,-20.0,+20.0,10,0.5,10.5);

  MyInit();
}

void AT_ReadTree::CheckEP1() {
  std::cout << "CheckEP1 called" << std::endl;
  std::cout << "opening runs.dat" << std::endl;
  ifstream fin("runs.dat");
  float bbcqx[2][6][100];
  float bbcqy[2][6][100];
  float runs[100];
  int run;
  int ir = 0;
  for(;ir!=100;++ir) {
    fin >> run;
    LoadTableEP(run);
    runs[ir] = run;
    for(int se=0; se!=2; ++se) {
      for(int ord=0; ord!=6; ++ord) {
	bbcqx[se][ord][ir] = bbcm[se][ord][0][2][20]; // bce=2 bvtx=20
	bbcqy[se][ord][ir] = bbcm[se][ord][1][2][20]; // bce=2 bvtx=20
      }
    }
  }
  TCanvas *main = new TCanvas();
  TGraph *grx[2][6];
  TGraph *gry[2][6];
  main->Divide(6,2);
  int color[2] = {kRed-3,kBlue-3};
  for(int ord=0; ord!=6; ++ord) {
    for(int se=0; se!=2; ++se) {
      grx[se][ord] = new TGraph( ir, runs, bbcqx[se][ord] );
      gry[se][ord] = new TGraph( ir, runs, bbcqy[se][ord] );
      grx[se][ord]->SetMarkerStyle(24);
      gry[se][ord]->SetMarkerStyle(24);
      grx[se][ord]->SetMarkerColor( color[se] );
      gry[se][ord]->SetMarkerColor( color[se] );
      grx[se][ord]->SetLineColor( color[se] );
      gry[se][ord]->SetLineColor( color[se] );
    }
    main->cd(1+ord+0*6);
    grx[0][ord]->Draw("APL");
    grx[1][ord]->Draw("PLSAME");
    grx[0][ord]->GetYaxis()->SetRangeUser(-12,+12);
    main->cd(1+ord+1*6);
    gry[0][ord]->Draw("APL");
    gry[1][ord]->Draw("PLSAME");
    gry[0][ord]->GetYaxis()->SetRangeUser(-2,+2);
  }
  main->SaveAs("CheckEP1.root","root");
}

void AT_ReadTree::CheckEP2() {
  std::cout << "CheckEP2 called" << std::endl;
  std::cout << "opening runs.dat" << std::endl;
  ifstream fin("runs.dat");
  float bbcqc[32][6][100];
  float bbcqs[32][6][100];
  float runs[100];
  int run;
  int ir = 0;
  for(;ir!=100;++ir) {
    fin >> run;
    LoadTableEP(run);
    runs[ir] = run;
    for(int se=0; se!=32; ++se) {
      for(int ord=0; ord!=6; ++ord) {
	bbcqc[se][ord][ir] = bbcc[se][ord][3][20]; // bce=2 bvtx=20
	bbcqs[se][ord][ir] = bbcs[se][ord][3][20]; // bce=2 bvtx=20
      }
    }
  }
  TCanvas *main1 = new TCanvas();
  TCanvas *main2 = new TCanvas();
  TGraph *grx[32][6];
  TGraph *gry[32][6];
  main1->Divide(8,6);
  main2->Divide(8,6);
  int color[4] = { kRed-3, kOrange-3, kCyan-3, kBlue-3};
  for(int ord=0; ord!=6; ++ord) {
    for(int se=0; se!=32; ++se) {
      grx[se][ord] = new TGraph( ir, runs, bbcqc[se][ord] );
      gry[se][ord] = new TGraph( ir, runs, bbcqs[se][ord] );
      grx[se][ord]->SetMarkerStyle(24);
      gry[se][ord]->SetMarkerStyle(24);
      grx[se][ord]->SetMarkerColor( color[se%4] );
      gry[se][ord]->SetMarkerColor( color[se%4] );
      grx[se][ord]->SetLineColor( color[se%4] );
      gry[se][ord]->SetLineColor( color[se%4] );
    }
    for(int se=0; se!=32; ++se) {
      int cvs = se/4;
      main1->cd(1+3*cvs+ord);
      if((se%4)==0) {
	grx[se][ord]->Draw("APL");
	grx[se][ord]->GetYaxis()->SetRangeUser(-0.035,+0.035);
      } else {
	grx[se][ord]->Draw("PLSAME");
      }
      main2->cd(1+3*cvs+ord);
      if((se%4)==0) {
	gry[se][ord]->Draw("APL");
	gry[se][ord]->GetYaxis()->SetRangeUser(-0.035,+0.035);
      } else {
	gry[se][ord]->Draw("PLSAME");
      }
    }
  }
  main1->SaveAs("CheckEP2_1.root","root");
  main2->SaveAs("CheckEP2_2.root","root");
}


void AT_ReadTree::Finish() {
  hEvents->Write();
  MyFinish();
}

AT_ReadTree::~AT_ReadTree() {
  if(hEvents) delete hEvents;
}

void AT_ReadTree::Exec() {
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  unsigned int trigger = fGLB.trig;
  unsigned int kBBCnc = 0x00000008;
  unsigned int kBBCn  = 0x00000010;
  unsigned int  mask = kBBCnc | kBBCn;
  bool trig = false;
  if(trigger & mask) trig = true;
  float frac = fGLB.frac;

  if(cen<0.5||cen>60.5) return;
  if(!trig) return;
  if(frac<0.95) return;
  if(TMath::Abs(vtx)>20) return;
  //std::cout << " " << cen << " " << frac << " " << vtx << std::endl;

  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );
  //std::cout << "  " << bcen << " " << bvtx << std::endl;

  if(bvtx<0 || bcen<0) return;

  hEvents->Fill(vtx,cen);
  //std::cout << hEvents->GetEntries() << std::endl;

  MakeBBCEventPlanes(bcen,bvtx);

  MyExec();
}

//BBC EVENTPLANE
void AT_ReadTree::MakeBBCEventPlanes(int bcen, int bvtx) {
  return;
  qcQ bbq[4][3];
  for(int se=0; se!=2; ++se) {
    bbq[0][se] = pQ1bb->at(se);
    bbq[1][se] = pQ2bb->at(se);
    bbq[2][se] = pQ3bb->at(se);
    bbq[3][se] = pQ4bb->at(se);
    for(int ord=0; ord!=4; ++ord) {
      bbq[ord][se].SetXY( bbq[ord][se].X() - bbcm[se][ord][0][bcen][bvtx],
			  bbq[ord][se].Y() - bbcm[se][ord][1][bcen][bvtx],
			  bbq[ord][se].NP(),
			  bbq[ord][se].M() );
    }
  }
  for(int ord=0; ord!=4; ++ord) {
    bbq[ord][2] = bbq[ord][0] + bbq[ord][1];
  }

  Psi1_BBC = bbq[0][2].Psi();
  Psi2_BBC = bbq[1][2].Psi();
  Psi3_BBC = bbq[2][2].Psi();
  Psi4_BBC = bbq[4][2].Psi();
}

int AT_ReadTree::ReferenceTracks() {
  int ntrk=0;
  uint ntrks = pTRKpt->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float zed  = pTRKzed->at(itrk);
    int qua = pTRKqua->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    if(qua!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    ntrk++;
  }
  return ntrk;
}
int AT_ReadTree::BinVertex(float vtx) {
  int ret=-1;
  for(int i=0; i!=fNBinsVtx+1; ++i) {
    if(vtx<fMinBinVtx+i) {
      ret = i-1;
      break;
    }
  }
  return ret;
}
int AT_ReadTree::BinCentrality(float cen) {
  int ret=-1;
  for(int i=0; i!=fNBinsCen+1; ++i) {
    if(cen<fMinBinCen+i) {
      ret = i-1;
      break;
    }
  }
  return ret;
}

void AT_ReadTree::LoadTableEP( int run ) {
  if(run<0) {
    Analysis *ana = Analysis::Instance();
    run = ana->RunNumber();
    std::cout << " SEGMENT " << ana->SegmentNumber() << std::endl;
  }
  std::cout << " RUN " << run << std::endl;

  ifstream fin;
  int se, ord, xy, bce, bvt;
  float tmp;
  fin.open( Form("BBC_EPC/tables/BBC_%d.dat",run) );
  int nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/9600)%6;
    int xy = (nn/4800)%2;
    int se = (nn/2400)%2;
    int bce = (nn/40)%60;
    int bvt = nn%40;
    bbcm[se][ord][xy][bce][bvt] = tmp*1e-1;
  }
  fin.close();
  std::cout << "   BBC ReCenter coefficients loaded: " << nn << std::endl;
  fin.open( Form("BBC_EPC/tables/BBC_A_%d.dat",run) );
  nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/230400)%4;
    int bce = (nn/3840)%60;
    int bcs = (nn/1920)%2;
    int bor = (nn/40)%32;
    int bvt = nn%40;
    //if(bcs==0) bbcc[bor][ord][bce][bvt] = tmp*1e-3;
    //else bbcs[bor][ord][bce][bvt] = tmp*1e-3;
  }
  std::cout << "   BBC Flattening coefficients loaded: " << nn << std::endl;
}
