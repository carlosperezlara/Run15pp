#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TString.h"
#include "PbScIndexer.h"
#include "PbScIndexer.C"
#include "PbGlIndexer.h"
#include "PbGlIndexer.C"
#include "EmcIndexer.h"
#include "EmcIndexer.C"
#include "TMath.h"

const int NTWRS = 24768;
const double c = 29.979245829979; //[cm/ns]

int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS [NPAS=0]" << std::endl;
    return 1;
  }
  TString srun = argv[1];
  TString snev = argv[2];
  TString spas = argv[3];
  int nev = snev.Atoi();
  int pas = spas.Atoi();

  //===================
  float min[NTWRS]; // minimum and maximum range in rawprof
  float max[NTWRS]; // minimum and maximum range in rawprof
  float p0[NTWRS]; // parameters of walk fit
  float p1[NTWRS]; // parameters of walk fit
  float p2[NTWRS]; // parameters of walk fit
  float MyLC[NTWRS];
  for(int i=0; i!=NTWRS; ++i) {
    min[i] = 0;
    max[i] = 3000;
    p0[i] = 0;
    p1[i] = 0;
    p2[i] = 0;
  }
  std::ifstream rr;

  rr.open("lcs.txt");
  for(int tid;;) {
    rr >> tid;
    if(!rr.good()) break;
    rr >> MyLC[tid];
  }
  rr.close();

  float off[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  if(pas>0) {
    std::cout << "PAS 1: LOADING SECTOR OFFSET PER RUN FROM FILE offs.txt" << std::endl;
    rr.open("offs.txt");
    float tmp;
    for(;;) {
      TString dbrun;
      rr >> dbrun;
      if(!rr.good()) break;
      if(dbrun==srun) {
	rr >> off[0] >> off[1] >> off[2] >> off[3];
	rr >> off[4] >> off[5] >> off[6] >> off[7];
	std::cout << " READ VALUES FOR RUN " << srun.Data() << std::endl; 
      } else {
	rr >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
      }
    }
    rr.close();
  }
  if(pas>1) {
    std::cout << "PAS 2: LOADING RANGES FROM FILE ranges.txt" << std::endl;
    rr.open("ranges.txt");
    for(;;) {
      int tor=0;
      rr >> tor;
      if(!rr.good()) break;
      rr >> min[tor] >> max[tor];
    }
    rr.close();
  }
  if(pas>2) {
    std::cout << "PAS 3: LOADING WALK FIT FROM FILE walks.txt" << std::endl;
    rr.open("walks.txt");
    for(;;) {
      int tor=0;
      rr >> tor;
      if(!rr.good()) break;
      rr >> p0[tor] >> p1[tor] >> p2[tor];
    }
    rr.close();
  }
  float off2[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  if(pas>3) {
    std::cout << "PAS 4: LOADING SECTOR OFFSET PER RUN FROM FILE offset_numbertwo.txt" << std::endl;
    rr.open("offset_numbertwo.txt");
    float tmp;
    for(;;) {
      TString dbrun;
      rr >> dbrun;
      if(!rr.good()) break;
      if(dbrun==srun) {
	rr >> off2[0] >> off2[1] >> off2[2] >> off2[3];
	rr >> off2[4] >> off2[5] >> off2[6] >> off2[7];
	std::cout << " READ VALUES FOR RUN " << srun.Data() << std::endl; 
      } else {
	rr >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
      }
    }
    rr.close();
  }
  float mean2[24768];
  if(pas>4) {
    std::cout << "PAS 5: LOADING TOWERID OFFSET FROM file meansnumbertwo.dat" << std::endl;
    rr.open("meansnumbertwo.dat");
    float tmp;
    for(;;) {
      int tor=0;
      rr >> tor;
      if(!rr.good()) break;
      rr >> mean2[tor] >> tmp;
    }
    rr.close();
  }
  std::cout << "ALL LOADED." << std::endl;

  //===================

  TH2F *hTOF=NULL, *hTOF2=NULL;
  if(pas>2) {
    hTOF = new TH2F("hTOF",";towerID; ToF  [ns]",24768,-0.5,24767.5,900,-200,100);
    hTOF2 = new TH2F("hTOF2",";towerID;ToF  [ns]",24768,-0.5,24767.5,100,-5,+5);
  }
  TH2F* hTvsA[NTWRS];
  TProfile* sTvsA[8];
  for(int i=0;i<8;i++){
    if(pas<2) {
      sTvsA[i]=new TProfile(Form("PROF_TDCvsADC_SEC%d",i),";ADC;<TDC>",
			    200,0.,4096.,0.,4096.,"s");
    } else {
      sTvsA[i] = NULL;
    }
  }
  for(int i=0;i<NTWRS;i++){
    hTvsA[i]=NULL;
    if(pas<4) {
      float rymin = min[i]; //0
      float rymax = max[i]; //3000
      int nbinsx = 100; //default binning
      int nbinsy = 100; //default binning
      hTvsA[i]=new TH2F(Form("H2_TDCvsADC_%d",i),";ADC;TDC",
			nbinsx,0.,4096.,nbinsy,rymin,rymax);
    }
  }
  TString inFile=Form("trees/%s.root",srun.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  TTree* T = (TTree*)f1->Get("nt");
  float TDC,ADC,LC,WK,x,y,z,ecore,ecent,ctof,VTX,BBCT0;
  float towerid;
  T->SetBranchAddress("T0",&BBCT0);
  T->SetBranchAddress("TDC",&TDC);
  T->SetBranchAddress("ADC",&ADC);
  T->SetBranchAddress("LC",&LC);
  T->SetBranchAddress("WK",&WK);
  T->SetBranchAddress("x",&x);
  T->SetBranchAddress("y",&y);
  T->SetBranchAddress("z",&z);
  T->SetBranchAddress("ecore",&ecore);
  T->SetBranchAddress("ecent",&ecent);
  T->SetBranchAddress("ctof",&ctof);
  T->SetBranchAddress("VTX",&VTX);
  T->SetBranchAddress("towerid",&towerid);
  
  Long64_t nevt = T->GetEntries();
  if(nev>0) nevt = nev; // overriding from parameter
  for(Long64_t i1=0; i1<nevt; ++i1) {
    if(i1%2000000 == 0)
      std::cout << "Event:  " << i1 << "/" << nevt << Form(" (%.1f)",i1*100.0/nevt) << std::endl; 
    T->GetEntry(i1);
    int tid = int(towerid+0.1);
    int isc, iz, iy;
    float lc = MyLC[tid];
    //lc=LC;
    EmcIndexer::decodeTowerId(tid,isc,iz,iy);
    float tdcoff = TDC+off[isc]+lc*BBCT0;
    if(hTvsA[tid]) {
      if(pas>2) {
	float mean = p0[tid]+p1[tid]/ADC+p2[tid]/ADC/ADC;
	float mmax = mean + 50;
	float mmin = mean - 50;
	if(tdcoff>mmin&&tdcoff<mmax) // cut around minima
	  hTvsA[tid]->Fill(ADC,tdcoff);
      } else {
	hTvsA[tid]->Fill(ADC,tdcoff);
      }
    }
    if(sTvsA[isc]) {
      sTvsA[isc]->Fill(ADC,tdcoff);
    }
    z-=VTX;
    double d=TMath::Sqrt( x*x + y*y + z*z );
    float light=d/c; // ns
    float tdccor = TDC-p1[tid]/ADC-p2[tid]/ADC/ADC -p0[tid]+off[isc] ;
    float tof=-lc*tdccor-light-BBCT0-off2[isc]-mean2[tid];
    if(hTOF) hTOF->Fill(tid,tof);
    if(hTOF2) hTOF2->Fill(tid,tof);
  }
  TString ooo = Form("histresults_P%d_%s.root",pas,srun.Data());
  std::cout << "SAVING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data(),"RECREATE"); 
  for(int i=0;i<NTWRS;i++) {
    if(hTvsA[i]) hTvsA[i]->Write();
  }
  for(int i=0;i<8;i++) {
    if(sTvsA[i]) sTvsA[i]->Write();
  }
  if(hTOF) hTOF->Write();
  if(hTOF2) hTOF2->Write();
  std::cout << "DONE =)" << std::endl;
  return 0;
}
