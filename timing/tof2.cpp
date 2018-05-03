#include "defs2.h"
#include "configtime.h" //step1_sector[8] isbad[24768] min[24768] max[24768]

int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS" << std::endl;
    return 1;
  }
  LoadLCs();
  LoadBadTowers();
  TString srun = argv[1];
  TString snev = argv[2];
  TString schu = argv[3];
  int nev = snev.Atoi();
  int myrun = srun.Atoi();
  step1(myrun);
  loadranges();
  definebins();
  cout << myrun << endl;

  //===================
  TString inFile=Form("trees/%s.root",srun.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  fTree = (TTree*)f1->Get("TOP");

  fTree->Print();
  #include "branches2.h"
  unsigned int trig = kBBCnc | kBBCn;
  TH2F *hADCTDC[24768];
  for(int i=0;i!=24768;i++){
    float rymin = tdcmin[i];
    float rymax = tdcmax[i];
    hADCTDC[i]=new TH2F(Form("H2_TDCvsADC_%d",i),";ADC;TDC",
                        100,0.,4096.,100,rymin,rymax);
  }

  TH1F *hEvents = new TH1F("Events","",1,-0.5,+0.5);
  Long64_t nevt = fTree->GetEntries();
  if(nev>0) nevt = nev; // overriding from parameter
  int isc,iz,iy;
  int goodevents = 0;
  for(Long64_t i1=0; i1<nevt; ++i1) {
    if(i1%500000 == 0)
      std::cout << "Event:  " << i1 << "/" << nevt << Form(" (%.1f)",i1*100.0/nevt) << std::endl; 
    fTree->GetEntry(i1);
    if(fGLB.frac<0.95) continue;
    if(TMath::Abs(fGLB.vtxZ)>10.0) continue;
    if(fGLB.trig&trig==0) continue;
    if(fGLB.cent>5||fGLB.cent<0) continue;
    float T0 = fGLB.bbct;
    goodevents++;
    unsigned int nnn = pEMCtwrid->size();
    for(unsigned int i=0; i!=nnn ;++i) {
      int twid = pEMCtwrid->at(i);
      if(isbad[twid]!=0) continue;
      int adc = pEMCadc->at(i);
      int tdc = pEMCtdc->at(i);
      float lc = MyLC[twid];
      EmcIndexer::decodeTowerId(twid,isc,iz,iy);
      float ftdc = tdc +lc*T0 -step1_sector[isc] + 2000;
      hADCTDC[twid]->Fill( adc, ftdc );
    }
  }
  hEvents->SetBinContent(1,goodevents);
  TString outfile = Form("out/MyResultsTOF2_%s.root",srun.Data());
  TFile *fout = new TFile( outfile.Data(), "RECREATE" );
  hEvents->Write();
  for(int i=0; i!=24768; ++i) {
    hADCTDC[i]->Write();
  }
  fout->Close();
  std::cout << "OUTPUT IS HERE: " << outfile.Data() << std::endl;
  return 0;
}
