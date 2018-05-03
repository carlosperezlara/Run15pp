#include "defs2.h"
#include "configtime.h" //isbad[24768]

int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS [NPAS=0]" << std::endl;
    return 1;
  }
  LoadBadTowers();
  LoadLCs();
  TString srun = argv[1];
  TString snev = argv[2];
  TString spas = argv[3];
  int nev = snev.Atoi();
  int pas = spas.Atoi();

  //===================
  TString inFile=Form("trees/%s.root",srun.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  fTree = (TTree*)f1->Get("TOP");

  fTree->Print();
  #include "branches2.h"
  unsigned int trig = kBBCnc | kBBCn;
  TProfile *hMeanTDC[8];
  for( int t=0; t!=8; ++t ) {
    hMeanTDC[t] = new TProfile( Form("MeanTDC_SEC%d",t),"",2048,-0.5,4095.5,"S");
  }
  TH1F *hEvents = new TH1F("Events","",1,-0.5,+0.5);
  Long64_t nevt = fTree->GetEntries();
  if(nev>0) nevt = nev; // overriding from parameter
  int isc,iz,iy;
  int goodevents = 0;
  for(Long64_t i1=0; i1<nevt; ++i1) {
    if(i1%200000 == 0)
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
      float tdcsupp = tdc + lc*T0; // con esto calculo offsetsec(run)
      hMeanTDC[isc]->Fill( adc, tdcsupp );
    }
  }
  hEvents->SetBinContent(1,goodevents);
  TString outfile = Form("out/MyResultsTOF1_%s.root",srun.Data());
  TFile *fout = new TFile( outfile.Data(), "RECREATE" );
  hEvents->Write();
  for(int i=0; i!=8; ++i) {
    hMeanTDC[i]->Write();
  }
  fout->Close();
  std::cout << "OUTPUT IS HERE: " << outfile.Data() << std::endl;
  return 0;
}
