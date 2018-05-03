#include "defs2.h"
#include "configtime.h"

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
  int myrun = srun.Atoi();
  step1(myrun);

  //===================
  TString inFile=Form("trees/%s.root",srun.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  fTree = (TTree*)f1->Get("TOP");

  fTree->Print();
  #include "branches2.h"
  unsigned int trig = kBBCnc | kBBCn;
  TProfile2D *hMeanTDC[8];
  definebins();
  //for(int i=0; i!=nbinsx; ++i) {
  //  cout << i << " " << binx[i] << endl;
  //}
  for( int t=0; t!=8; ++t ) {
    //hMeanTDC[t] = new TProfile2D( Form("MeanTDC_SEC%d",t),"",2048,-0.5,4095.5,
    hMeanTDC[t] = new TProfile2D( Form("MeanTDC_SEC%d",t),"",nbinsx-1,binx,
				  secs[t+1]-secs[t],secs[t]-0.5,secs[t+1]-0.5);
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
      //float ene = pEMCene->at(i);
      //if(ene<0.2) continue; // i stored ecent no eclus
      EmcIndexer::decodeTowerId(twid,isc,iz,iy);
      float ftdc = tdc - step1_sector[isc] + 2000 +lc*T0;
      hMeanTDC[isc]->Fill( adc, twid, ftdc );
    }
  }
  hEvents->SetBinContent(1,goodevents);
  TString outfile = Form("out/MyResults2_%s.root",srun.Data());
  TFile *fout = new TFile( outfile.Data(), "RECREATE" );
  hEvents->Write();
  for(int i=0; i!=8; ++i) {
    hMeanTDC[i]->Write();
  }
  fout->Close();
  std::cout << "OUTPUT IS HERE: " << outfile.Data() << std::endl;
  return 0;
}
