#include "defs2.h"

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
  TString inFile=Form("trees/%s.root",srun.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  fTree = (TTree*)f1->Get("TOP");

  fTree->Print();
  #include "branches2.h"
  //const int kNTWRS = 24768;
  //unsigned int kBBCnc = 0x00000008;
  //unsigned int kBBCn  = 0x00000010;
  int secs[9] = { 0, 2592, 5184, 7776, 10368,
		  12960, 15552, 20160, 24768 };
  unsigned int trig = kBBCnc | kBBCn;
  TH2F *hADC[8], *hTDC[8], *hTDCADC[8];
  TProfile *hMeanTDC = new TProfile("MeanTDC","",24768,-0.5,24767.5);
  TH1F *hHitTwr = new TH1F("HitTwr","", 24768,-0.5,24767.5);
  TH1F *hEvents = new TH1F("Events","",1,-0.5,+0.5);
  for(int i=0; i!=8; ++i) {
    hADC[i] = new TH2F( Form("ADC_SEC%d",i),"",
			secs[i+1]-secs[i],
			secs[i]-0.5,
			secs[i+1]-0.5,
			100,0,3000);
    hTDC[i] = new TH2F( Form("TDC_SEC%d",i),"",
			secs[i+1]-secs[i],
			secs[i]-0.5,
			secs[i+1]-0.5,
			100,-4000,+4000);
    hTDCADC[i] = new TH2F( Form("TDCADC_SEC%d",i),"",
			   100,0,4096,
			   100,-4000,4000);
  }

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
      int adc = pEMCadc->at(i);
      int tdc = pEMCtdc->at(i);
      EmcIndexer::decodeTowerId(twid,isc,iz,iy);
      hADC[isc]->Fill( twid, adc );
      hTDC[isc]->Fill( twid, tdc );
      hTDCADC[isc]->Fill(adc,tdc);
      if(adc>1000||adc<4096)
	hMeanTDC->Fill( twid, tdc );
      if( pEMCene->at(i)>0.2 && pEMCene->at(i)<6.0 )
	hHitTwr->Fill(twid);
    }
  }
  hEvents->SetBinContent(1,goodevents);
  TH1F *hHitTwrSec[8];
  TString outfile = Form("out/MyResults_%s.root",srun.Data());
  TFile *fout = new TFile( outfile.Data(), "RECREATE" );
  hEvents->Write();
  for(int i=0; i!=8; ++i) {
    hADC[i]->Write();
    hTDC[i]->Write();
    hTDCADC[i]->Write();
  }
  hMeanTDC->Write();
  hHitTwr->Write();
  float ddd[8] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
		   0.2, 0.2 };
  for(int i=0; i!=8; ++i) {
    hHitTwrSec[i] = new TH1F( Form("hHitTwrSec%d",i), "", 100, 0.0, ddd[i] );
    for(int j=secs[i]; j!=secs[i+1]; ++j) {
      hHitTwrSec[i]->Fill( hHitTwr->GetBinContent(j+1)/(goodevents/1000.0) );
    }
    hHitTwrSec[i]->Write();
  }

  fout->Close();
  std::cout << "OUTPUT IS HERE: " << outfile.Data() << std::endl;
  return 0;
}
