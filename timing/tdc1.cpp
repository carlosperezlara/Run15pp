#include "defs2.h"
int isbad[24768];
void LoadBadTowers() {
  for(int i=0; i!=24768; ++i) isbad[i] = 0;
  ifstream badt("badtowers.dat");
  int tid, bd;
  int nnn=0;
  for(;;++nnn) {
    badt >> tid;
    if(!badt.good()) break;
    badt >> bd;
    isbad[tid] = bd;
  }
  badt.close();
  cout << " BAD TOWERS LOADED " << nnn << endl;
}

int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS [NPAS=0]" << std::endl;
    return 1;
  }
  LoadBadTowers();
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
  int secs[9] = { 0, 2592, 5184, 7776, 10368,
		  12960, 15552, 20160, 24768 };
  unsigned int trig = kBBCnc | kBBCn;
  TProfile2D *hMeanTDC[8];
  TProfile *hMeanTDC2[8];
  for( int t=0; t!=8; ++t ) {
    hMeanTDC[t] = new TProfile2D( Form("MeanTDC_SEC%d",t),"",2048,-0.5,4095.5,
				  secs[t+1]-secs[t],secs[t]-0.5,secs[t+1]-0.5);
    hMeanTDC2[t] = new TProfile( Form("MeanTDC2_SEC%d",t),"",
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
    goodevents++;
    unsigned int nnn = pEMCtwrid->size();
    for(unsigned int i=0; i!=nnn ;++i) {
      int twid = pEMCtwrid->at(i);
      if(isbad[twid]!=0) continue;
      int adc = pEMCadc->at(i);
      int tdc = pEMCtdc->at(i);
      EmcIndexer::decodeTowerId(twid,isc,iz,iy);
      hMeanTDC[isc]->Fill( adc, twid, tdc );
      if(adc>2500&&adc<4000) {
	hMeanTDC2[isc]->Fill( twid, tdc );
      }
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
