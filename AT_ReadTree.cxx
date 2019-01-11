#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "Analysis.h"
#include "AT_ReadTree.h"

AT_ReadTree::AT_ReadTree() : AnalysisTask() {
  // -20.0 ==> +20.0 (40+1)
  // 0.5 ==> 60.5 (60+1)
  fBBCQCal = true;
  Psi_BBC = false;
  fNBinsVtx = 40;
  fNBinsCen = 1;
  fMinBinVtx = -20.0;
  fMinBinCen = 9.5;
  unsigned int kBBC1 = 0x00000001;
  unsigned int kBBC2 = 0x00000002;
  unsigned int kBBCnc = 0x00000008;
  unsigned int kBBCn  = 0x00000010;
  fMask = kBBC1 | kBBC2 | kBBCnc | kBBCn;
  fCentralityMin = 0.0;
  fCentralityMax = 80.0;
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
  pMXSempc3x3 = new std::vector<Float_t>;
  pQ1ex = new std::vector<qcQ>;
  pQ2ex = new std::vector<qcQ>;
  pQ3ex = new std::vector<qcQ>;
  pQ4ex = new std::vector<qcQ>;
  pQ6ex = new std::vector<qcQ>;
  pQ8ex = new std::vector<qcQ>;
  pQ1fv = new std::vector<qcQ>;
  pQ2fv = new std::vector<qcQ>;
  pQ3fv = new std::vector<qcQ>;
  pQ1bb = new std::vector<qcQ>;
  pQ2bb = new std::vector<qcQ>;
  pQ3bb = new std::vector<qcQ>;
  pQ4bb = new std::vector<qcQ>;
  pQ6bb = new std::vector<qcQ>;
  pQ8bb = new std::vector<qcQ>;

  pEMCid = new std::vector<Int_t>;
  pEMCtwrid = new std::vector<Int_t>;
  pEMCx = new std::vector<Float_t>;
  pEMCy = new std::vector<Float_t>;
  pEMCz = new std::vector<Float_t>;
  pEMCecore = new std::vector<Float_t>;
  pEMCecent = new std::vector<Float_t>;
  pEMCchisq = new std::vector<Float_t>;
  pEMCtimef = new std::vector<Float_t>;

  pTRKqua = new std::vector<Int_t>;
  pTRKpt = new std::vector<Float_t>;
  pTRKphi = new std::vector<Float_t>;
  pTRKpz = new std::vector<Float_t>;
  pTRKecore = new std::vector<Float_t>;
  pTRKetof = new std::vector<Float_t>;
  pTRKtwrid = new std::vector<Int_t>;
  pTRKplemc = new std::vector<Float_t>;
  pTRKchisq = new std::vector<Float_t>;
  pTRKdphi = new std::vector<Float_t>;
  pTRKdz = new std::vector<Float_t>;
  pTRKpc3sdphi = new std::vector<Float_t>;
  pTRKpc3sdz = new std::vector<Float_t>;
  pTRKzed = new std::vector<Float_t>;
  pTRKdisp = new std::vector<Float_t>;
  pTRKprob = new std::vector<Float_t>;
  pTRKcid = new std::vector<Int_t>;

  pMXSempccent = new std::vector<Float_t>;
  pMXSempc3x3 = new std::vector<Float_t>;
  pMXSpt = new std::vector<Float_t>;
  pMXSpz = new std::vector<Float_t>;
  pMXSeta = new std::vector<Float_t>;
  pMXSphi = new std::vector<Float_t>;
  pMXSflyr = new std::vector<Int_t>;
  pMXSsingleD = new std::vector<Float_t>;
  pMXSsingleP = new std::vector<Int_t>;

  fSkipPileUpCuts = false;
  fSkipDetails = false;
  fSkipClusters = false;
  fSkipTracks = false;
  fSkipShowers = false;
}

void AT_ReadTree::Init() {
  Analysis *ana = Analysis::Instance();
  fCandidates = ana->GetCandidates();
  fCandidates2= ana->GetCandidates2();
  for(int i=0; i!=4; ++i) fQ[i] = ana->GetQ(i);
  hEvents = new TH1F("hEvents","hEvents",10,-0.5,9.5);
  hEvents->GetXaxis()->SetBinLabel(1,"AllEvents");
  hEvents->GetXaxis()->SetBinLabel(2,"AT_ReadTree");
  hEvents->GetXaxis()->SetBinLabel(3,"AT_ReadTree_PileUp");
  hEvents->GetXaxis()->SetBinLabel(4,"AT_ReadTree");
  hEvents->GetXaxis()->SetBinLabel(5,"AT_ReadTree");
  hEvents->GetXaxis()->SetBinLabel(6,"AT_ReadTree");

  hTriggers0 = new TH1F("hTriggers","hTriggers",10,-0.5,9.5);
  hTriggers0->GetXaxis()->SetBinLabel(1,"0x00000001");
  hTriggers0->GetXaxis()->SetBinLabel(2,"0x00000002");
  hTriggers0->GetXaxis()->SetBinLabel(3,"0x00000008");
  hTriggers0->GetXaxis()->SetBinLabel(4,"0x00000010");

  hCentrality0 = new TH1F("hCentrality","hCentrality",100,-0.5,99.5);
  hVertex0 = new TH1F("hVertex","hVertex",100,-40,+40);
  hPileUpRejectionS = new TH2F("hPileUpRejectionS","PileUpRejection;;BBCs",2,-0.5,1.5,100,0.,200.);
  hPileUpRejectionN = new TH2F("hPileUpRejectionN","PileUpRejection;;BBCn",2,-0.5,1.5,100,0.,200.);
  hCentralitySelection = new TH2F("hCentralitySelection",";BBCs;BBCn",1000,0.,200.,1000,0.,200.);
  
  for (int i = 0; i < 5; i++) {
    hPsi2[i] = new TH1F(Form("hPsi2%d", i), Form("2nd order Psi after calib step %d", i), 200, -6.3, 6.3);
  }

  TTree *tree = ana->GetTree();
  if(!tree) {
    std::cout << "AT_ReadTree:Init says: Tree not found." << std::endl;
    return;
  }
  //Opening assigning branches
  tree->SetBranchAddress("Event",&fGLB);
  tree->SetBranchAddress("EventC",&fGLB2);
  if(!fSkipDetails) {
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
  }
  LoadTableEP();
  LoadTableTime();

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
  float bbcqc[32][4][100];
  float bbcqs[32][4][100];
  float runs[100];
  int run;
  int ir = 0;
  for(;ir!=100;++ir) {
    fin >> run;
    LoadTableEP(run);
    runs[ir] = run;
    for(int se=0; se!=32; ++se) {
      for(int ord=0; ord!=4; ++ord) {
	bbcqc[se][ord][ir] = bbcc[se][ord][3][20]; // bce=2 bvtx=20
	bbcqs[se][ord][ir] = bbcs[se][ord][3][20]; // bce=2 bvtx=20
      }
    }
  }
  TCanvas *main1 = new TCanvas();
  TCanvas *main2 = new TCanvas();
  TGraph *grx[32][4];
  TGraph *gry[32][4];
  main1->Divide(8,4);
  main2->Divide(8,4);
  int color[4] = { kRed-3, kOrange-3, kCyan-3, kBlue-3};
  for(int ord=0; ord!=4; ++ord) {
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
  hTriggers0->Write();
  hCentrality0->Write();
  hVertex0->Write();
  hPileUpRejectionS->Write();
  hPileUpRejectionN->Write();
  hCentralitySelection->Write();
  for (int i = 0; i < 5; i++) {
    hPsi2[i]->Write();
  }
  MyFinish();

}

AT_ReadTree::~AT_ReadTree() {
  if(hEvents) delete hEvents;
  if(hCentrality0) delete hCentrality0;
}

void AT_ReadTree::Exec() {
  //std::cout << "BBB" << std::endl;
  hEvents->Fill(0);
  float vtx = fGLB.vtxZ;
  fGLB.cent = 10; //setting centraliy to 10 for everybody
  float cen = fGLB.cent;
  unsigned int trigger = fGLB.trig;
  bool trig = false;
  if(0x00000001 | trigger) hTriggers0->Fill(0);
  if(0x00000002 | trigger) hTriggers0->Fill(1);
  if(0x00000008 | trigger) hTriggers0->Fill(2);
  if(0x00000010 | trigger) hTriggers0->Fill(3);
  if(trigger & fMask) trig = true;
  if(!trig) return;
    float meanS = fGLB2.bbcsTmean;
    float meanN = fGLB2.bbcnTmean;
  if( TMath::IsNaN(meanS) || TMath::IsNaN(meanN) ) return;

  if(TMath::Abs(vtx)>20) return;
  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );
  if(bvtx<0 || bcen<0) return;

  float sgnS = fGLB.bbcs;
  float sgnN = fGLB2.bbcn;
  hPileUpRejectionS->Fill(0.,sgnS);
  hPileUpRejectionN->Fill(0.,sgnN);
  hEvents->Fill(1);
  if(!fSkipPileUpCuts) {
    double rmsS = fGLB2.bbcsTrms;
    double rmsN = fGLB2.bbcnTrms;
    if( TMath::Abs(meanS-meanN+fTime[0]) > 5*0.60 ) return;
    if( rmsS/fTime[1] + rmsN/fTime[2] > 1.0 ) return;
  }
  hPileUpRejectionS->Fill(1.,sgnS);
  hPileUpRejectionN->Fill(1.,sgnN);
  hCentralitySelection->Fill(sgnS,sgnN);

  hEvents->Fill(2);
  hCentrality0->Fill(cen);
  hVertex0->Fill(vtx);

  if(fBBCQCal&&!fSkipDetails) MakeBBCEventPlanes(0,bvtx);
  MyExec();
}

//BBC EVENTPLANE
void AT_ReadTree::MakeBBCEventPlanes(int bcen, int bvtx) {
  Psi_BBC;
  Psi1_BBC = 0;
  Psi2_BBC = 0;
  Psi3_BBC = 0;
  Psi4_BBC = 0;
  qcQ qvec[4][3];
  for(int se=0; se!=2; ++se) {
    qvec[0][se] = pQ1bb->at(se);
    qvec[1][se] = pQ2bb->at(se);
    qvec[2][se] = pQ3bb->at(se);
    qvec[3][se] = pQ4bb->at(se);
    if(qvec[0][se].M()<1) {
      Psi_BBC = false;
      return;
    }
  }

  //std::cout << "BBC EP Check0 : " << qvec[1][0].Psi2Pi() << std::endl;
  hPsi2[0]->Fill(qvec[1][0].Psi2Pi());

  // ======= STAGE 2: Recentering SubEvents (STEP1)  =======
  for(int k=0; k!=4; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double cn = bbcm[j][k][0][bcen][bvtx];
      double sn = bbcm[j][k][1][bcen][bvtx];
      qvec[k][j].SetXY( x - cn, y - sn, qvec[k][j].NP(), qvec[k][j].M() );
    }
  }
  //std::cout << " LOADED " << bbcm[0][1][0][bcen][bvtx] << " " << bbcm[0][1][1][bcen][bvtx] << std::endl;
  //std::cout << " LOADED " << bbcm[1][1][0][bcen][bvtx] << " " << bbcm[1][1][1][bcen][bvtx] << std::endl;
  //std::cout << "BBC EP Check1 : 1 " << qvec[1][0].Psi2Pi() << std::endl;
  //std::cout << "BBC EP Check1 : 0 " << qvec[1][0].Psi2Pi() << std::endl;
  //std::cout << "BBC EP Check1 : 1 " << qvec[1][1].Psi2Pi() << std::endl;
  hPsi2[1]->Fill(qvec[1][0].Psi2Pi());

  int twon[4] = {1,3,4,5}; // 1,2,3,4,6,8
  // ======= STAGE 4: Twisting SubEvents (STEP2)  =======
  for(int k=0; k!=4; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double c2n = bbcm[j][twon[k]][0][bcen][bvtx] / qvec[k][j].M();
      double s2n = bbcm[j][twon[k]][1][bcen][bvtx] / qvec[k][j].M();
      double ldaSm = s2n/(1.0+c2n);
      double ldaSp = s2n/(1.0-c2n);
      double den = 1.0 - ldaSm*ldaSp;
      qvec[k][j].SetXY( (x-ldaSm*y) / den,
                        (y-ldaSp*x) / den,
                        qvec[k][j].NP(),
                        qvec[k][j].M() );
      if( TMath::IsNaN( qvec[k][j].X() ) || TMath::IsNaN( qvec[k][j].Y() ) ) {
	std::cout << "Error building coefficient ";
	std::cout << " | qvec.M: " << qvec[k][j].M();
	std::cout << " | c2n: " << c2n;
	std::cout << " | s2n: " << s2n;
	std::cout << " | ldaSm: " << ldaSm;
	std::cout << " | ldaSp: " << ldaSp;
	std::cout << " | den: " << den << std::endl;
      }
    }
  }

  //std::cout << "BBC EP Check2 : " << qvec[1][0].Psi2Pi() << std::endl;
  hPsi2[2]->Fill(qvec[1][0].Psi2Pi());
  // ======= STAGE 6: Rescaling SubEvents (STEP3)  =======
  for(int k=0; k!=4; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double c2n = bbcm[j][twon[k]][0][bcen][bvtx] / qvec[k][j].M();
      double a2np = 1.0+c2n;
      double a2nm = 1.0-c2n;
      qvec[k][j].SetXY( x / a2np,
                        y / a2nm,
                        qvec[k][j].NP(),
                        qvec[k][j].M() );
      if( TMath::IsNaN( qvec[k][j].X() ) || TMath::IsNaN( qvec[k][j].Y() ) ) {
	std::cout << "Error building coefficient [2] ";
	std::cout << " | a2np: " << a2np;
	std::cout << " | a2nm: " << a2nm << std::endl;
      }
    }
  }

  //std::cout << "BBC EP Check3 : " << qvec[1][0].Psi2Pi() << std::endl;
  hPsi2[3]->Fill(qvec[1][0].Psi2Pi());
  // ======= STAGE 8: Bulding Full Q and Storing Flattening Coeficients  =======
  double delta[4] = {0,0,0,0};
  for(int k=0; k!=4; ++k) { // order
    qvec[k][2] = qvec[k][0] + qvec[k][1];
    double psi = qvec[k][2].Psi2Pi();
    for(int ik=0; ik!=32; ++ik) { // correction order
      int nn = ik+1;
      delta[k] -= TMath::Cos(nn*psi)*bbcs[ik][k][bcen][bvtx]*2.0/nn;
      delta[k] += TMath::Sin(nn*psi)*bbcc[ik][k][bcen][bvtx]*2.0/nn;
    }
    fQ[k]->CopyFrom( qvec[k][2] );
    double cn = TMath::Cos( (k+1)*delta[k] );
    double sn = TMath::Sin( (k+1)*delta[k] );
    double xprime = fQ[k]->X()*cn - fQ[k]->Y()*sn;
    double yprime = fQ[k]->X()*sn + fQ[k]->Y()*cn;
    fQ[k]->SetXY( xprime, yprime, fQ[k]->NP(),fQ[k]->M() );
  }

  Psi1_BBC = qvec[0][2].Psi2Pi()+delta[0];
  Psi2_BBC = qvec[1][2].Psi2Pi()+delta[1];
  Psi3_BBC = qvec[2][2].Psi2Pi()+delta[2];
  Psi4_BBC = qvec[3][2].Psi2Pi()+delta[3];

  //std::cout << "BBC EP Check4 : " << Psi2_BBC << std::endl;
  hPsi2[4]->Fill(Psi2_BBC);

  /*
  if( (TMath::Abs( fQ[0]->Psi2Pi() - Psi1_BBC ) < 1e-4) ||
      (TMath::Abs( fQ[1]->Psi2Pi() - Psi2_BBC ) < 1e-4) ||
      (TMath::Abs( fQ[2]->Psi2Pi() - Psi3_BBC ) < 1e-4) ||
      (TMath::Abs( fQ[3]->Psi2Pi() - Psi4_BBC ) < 1e-4)
      )
    std::cout << "NOT COMPATIBLE" << std::endl;
  */
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
  double tmp;
  fin.open( Form("BBC_EPC/tables/BBC_%d.dat",run) );
  std::cout <<  Form("LOADING TABLE1 BBC_EPC/tables/BBC_%d.dat",run) << std::endl;
  int nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/160)%6;//9600)%6;
    int xy = (nn/80)%2;//4800)%2;
    int se = (nn/40)%2;//2400)%2;
    int bce = 0; //(nn/40)%60; // centrality
    int bvt = nn%40; // vertex
    bbcm[se][ord][xy][bce][bvt] = tmp;
  }
  fin.close();
  std::cout << "   BBC ReCenter coefficients loaded: " << nn << std::endl;
  fin.open( Form("BBC_EPC/tables/BBC_A_%d.dat",run) );
  std::cout <<  Form("LOADING TABLE1 BBC_EPC/tables/BBC_A_%d.dat",run) << std::endl;
  nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/2560)%4;//153600)%4; //order
    int bce = 0;//(nn/2560)%60; //centrality
    int bcs = (nn/1280)%2; //xy
    int bor = (nn/40)%32; // order of expansion
    int bvt = nn%40; // vertex
    if(bcs==0) bbcc[bor][ord][bce][bvt] = tmp;
    else bbcs[bor][ord][bce][bvt] = tmp;
  }
  std::cout << "   BBC Flattening coefficients loaded: " << nn << std::endl;
}

void AT_ReadTree::LoadTableTime( int run ) {
  if(run<0) {
    Analysis *ana = Analysis::Instance();
    run = ana->RunNumber();
    std::cout << " SEGMENT " << ana->SegmentNumber() << std::endl;
  }
  std::cout << " RUN " << run << std::endl;
  
  ifstream fin;
  int irun;
  float tmp, mean, rmsS, rmsN;
  fin.open( "EventChecker/TimeConstants.dat" );
  int nn=0;
  fTime[0] = -9999;
  fTime[1] = 1;
  fTime[2] = 1;
  for(;;++nn) {
    fin >> irun >> mean >> tmp >> tmp >> tmp >> rmsS >> tmp >> tmp >> rmsN;
    if( irun == run ) {
      fTime[0] = mean;
      fTime[1] = rmsS;
      fTime[2] = rmsN;
      break;
    }
    if(!fin.good()) break;
  }
  fin.close();
  std::cout << " Time Constants " << fTime[0] << " " << fTime[1];
  std::cout << " " << fTime[2] << std::endl;
  return;
}
