#include "../defs.h"
#include "eventplanedefs.h"
//bool fVerbosity = true;
bool fVerbosity = false;

void CreateQCA() {
  q2ca[0].Reset();
  //q2ca[1].Reset();
  //q2ca[2].Reset();
  uint ntrks = pTRKcid->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float pt   = pTRKpt->at(itrk);
    float phi  = pTRKphi->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    float zed  = pTRKzed->at(itrk);
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    q2ca[0].Fill(phi);
    //if(pt>0.5&&pt<3) {
    //  if(phi>1.5 && phi<4.5) {
    //	q2ca[0].Fill(phi);
    //      } else {
    //	q2ca[1].Fill(phi);
    //      }
    //}
  }
  //q2ca[2] = q2ca[0] + q2ca[1];
}

int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS" << std::endl;
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  //======= HISTOS
  TH1F *hEvents = new TH1F("Events","",5,-0.5,4.5);
  hEvents->GetXaxis()->SetBinLabel(1,"All");
  hEvents->GetXaxis()->SetBinLabel(2,"triggered");
  hEvents->GetXaxis()->SetBinLabel(3,"centrality");
  hEvents->GetXaxis()->SetBinLabel(4,"fraction");
  TH1F *hFrac = new TH1F("FRAC","", 100, 0.0, 1.0);
  TH1F *hCent = new TH1F("CENT","", 100, 0.0, 100);
  TH1F *hVtxZ = new TH1F("VTXZ","", 100, -11, +11);

  TH1F *hEp   = new TH1F("Ep","",   100, 0.0, 1.3);
  TH1F *hChi2 = new TH1F("Chi2","", 100, 0.0, 10.);
  TH1F *hDphi = new TH1F("pc3Dphi","", 100,-30, +30);
  TH1F *hphi  = new TH1F("phi","",  100, -TMath::TwoPi(), TMath::TwoPi());
  TH1F *hDz   = new TH1F("pc3Dz","", 100,-30., 30.);
  TH1F *hZED  = new TH1F("ZED","",  100,-100,+100);
  TH2F *hBeta = new TH2F("Beta",";p [GeV/c];[ns]", 300,0,3,1800,-100,+500);

  TProfile2D *hV2[2];
  TH2F *hDP[2];
  TH2F *hDPPE[2];
  hV2[0] = new TProfile2D("hV2_E","", 13,-0.5,12.5, 30,0,3, -1.1,+1.1);
  hV2[1] = new TProfile2D("hV2_W","", 13,-0.5,12.5, 30,0,3, -1.1,+1.1);
  hDP[0] = new TH2F("DPhiE", "", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());
  hDP[1] = new TH2F("DPhiW", "", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());
  hDPPE[0] = new TH2F("DPhiER","", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());
  hDPPE[1] = new TH2F("DPhiWR","", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());

  CreateHistograms();

  LoadCalibs(run);
  //======== OPEN FILE 
  TString inFile=Form("trees/%s.root",run.Data());
  TFile* f1 = new TFile(inFile.Data(),"READ");
  fTree = (TTree*)f1->Get("TOP");
  #include "../branches.h"

  //======= MAIN LOOP
  Long64_t nevt = fTree->GetEntries();
  if(nev>0&&nev<nevt) nevt = nev; // overriding from parameter
  TLorentzVector ii,jj,pp;
  float RndmPSI = 0;
  for(Long64_t i1=0; i1<nevt; ++i1) {
    if(i1%2000000 == 0) {
      std::cout << "Event:  " << i1 << "/" << nevt;
      std::cout << Form(" (%.1f)",i1*100.0/nevt) << std::endl;
    }
    fTree->GetEntry(i1);
    hEvents->Fill(0.0);
    unsigned int trigger = fGLB.trig;
    unsigned int  mask = kBBCnc | kBBCn;
    bool trig = false;
    if(trigger & mask) trig = true;
    if(!trig) continue;
    hEvents->Fill(1.0);
    float cent = fGLB.cent;
    if(cent<0||cent>5) continue;
    hEvents->Fill(2.0);
    float frac = fGLB.frac;
    if(frac<0.95) continue;
    hEvents->Fill(3.0);
    float vtxZ = fGLB.vtxZ;

    if(fVerbosity) std::cout << "===> EVENT ACEPTED" << endl;

    q2bb[0] = pQ2bb->at(0);
    q2bb[1] = pQ2bb->at(1);
    q2bb[2] = q2bb[0] + q2bb[1];
    hBBCcent->Fill( cent, q2bb[2].M() );

    q2fv[0] = pQ2fv->at(0);

    q2ex[0] = pQ2ex->at(0);
    q2ex[1] = pQ2ex->at(1);
    q2ex[2] = pQ2ex->at(2);
    q2ex[3] = pQ2ex->at(3);
    q2ex[4] = pQ2ex->at(4);
    q2ex[5] = pQ2ex->at(5);
    q2ex[6] = pQ2ex->at(6);
    q2ex[7] = pQ2ex->at(7);

    CreateQCA();

    std::cout << " || CA " << q2ca[0].NP() << " " << q2ca[1].NP();
    std::cout << " || BB " << q2bb[0].NP() << " " << q2bb[1].NP();
    std::cout << " || FV " << q2fv[0].NP();
    std::cout << " || EX " << q2ex[0].NP() << " " << q2ex[1].NP();
    std::cout << " " << q2ex[2].NP() << " " << q2ex[3].NP();
    std::cout << " " << q2ex[3].NP() << " " << q2ex[4].NP();
    std::cout << " " << q2ex[5].NP() << " " << q2ex[6].NP();
    std::cout << " " << q2ex[6].NP() << " " << q2ex[7].NP() << std::endl;

    // histograms 
    hVtxZ->Fill( vtxZ );
    hCent->Fill( cent );
    hFrac->Fill( frac );

    //=================================
    //=================================
    // STEP 1:
    StoreQcomponents();
    if(fVerbosity) std::cout << "===> STEP 1 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.0.cpp"
 
    //=================================
    //=================================
    // STEP 2: RECENTER
    Recenterbb();
    Recenterfv();
    Recenterex(); //only EX subdetectors
    if(fVerbosity) std::cout << "===> STEP 2 [DONE]" << endl;
    //=================================
    //=================================


    //=================================
    //=================================
    // STEP 3: CONSTRUCT FULL EX
    q2ex[8] = q2ex[0] + q2ex[1] + q2ex[2] + q2ex[3]; // construct
    q2ex[8] = q2ex[4] + q2ex[5] + q2ex[6] + q2ex[7]; // full EX (after recenter its parts)
    hQx->Fill(12.,q2ex[8].X()); // Qx 12
    hQy->Fill(12.,q2ex[8].Y()); // Qy 12
    //---- ( Qaing, save to remove later)
    hQM->Fill( 12., q2ex[8].M() ); // Qred 12
    hQred0->Fill( 12., q2ex[8].Reduced() ); // Qred 12
    hPSI0->Fill( 12., q2ex[8].Psi2Pi() ); // PSI 12
    hEVC[0][10]->Fill( q2ex[8].Psi2Pi(), q2bb[2].Psi2Pi() ); // EVC 10
    //----
    q2ex[8].SetXY( (q2ex[8].X()-exqx[8][0])/*/exqx[8][1]*/,
		   (q2ex[8].Y()-exqy[8][0])/*/exqy[8][1]*/,
		   q2ex[8].NP(), q2ex[8].M() );
    if(fVerbosity) std::cout << "===> STEP 3 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.1.cpp"

    //=================================
    //=================================
    // STEP 4: FLATTENING
    GetEPbb();
    GetEPfv();
    GetEPex();
    if(fVerbosity) std::cout << "===> STEP 4 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.2.cpp"

    //=== START OF MAIN LOOP
    if(fVerbosity) std::cout << "===> START OF MAIN LOOP" << endl;
    uint ntrks = pTRKcid->size();
    for(uint itrk=0; itrk!=ntrks; ++itrk) {
      float pt   = pTRKpt->at(itrk);
      float pz   = pTRKpz->at(itrk);
      float phi  = pTRKphi->at(itrk);
      float ep   = pTRKecore->at(itrk) / TMath::Sqrt(pt*pt+pz*pz);
      float tof  = pTRKetof->at(itrk);
      float chi2 = pTRKchisq->at(itrk);
      float dphi = pTRKpc3sdphi->at(itrk);
      float dz   = pTRKpc3sdz->at(itrk);
      float zed  = pTRKzed->at(itrk);
      //float beta = tof/TMath::Sqrt(pt*pt+pz*pz);
      if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
      if(TMath::Abs(dphi)>3) continue;
      if(TMath::Abs(dz)>3) continue;

      if(pt>0.5&&pt<3) {
	if(phi>1.5 && phi<4.5) {
	  q2ca[0].Fill(phi);
	} else {
	  q2ca[1].Fill(phi);
	}
      }

      hEp->Fill( ep );
      hChi2->Fill( chi2 );
      hDphi->Fill( dphi );
      hphi->Fill( phi );
      hDz->Fill( dz );
      hZED->Fill( zed );
      hBeta->Fill( TMath::Sqrt(pt*pt+pz*pz), tof/TMath::Sqrt(pt*pt+pz*pz) );
      int ew;
      if(phi>1.5 && phi<4.5) {
	ew = 0; //EAST
      } else {
	ew = 1; //WEST
      }
      hV2[ew]->Fill( 0.0, pt, TMath::Cos( 2.0*(phi-bbEP[0]) ) );
      hV2[ew]->Fill( 1.0, pt, TMath::Cos( 2.0*(phi-bbEP[1]) ) );
      hV2[ew]->Fill( 2.0, pt, TMath::Cos( 2.0*(phi-bbEP[2]) ) );
      hV2[ew]->Fill( 3.0, pt, TMath::Cos( 2.0*(phi-fvEP[0]) ) );
      hV2[ew]->Fill( 4.0, pt, TMath::Cos( 2.0*(phi-exEP[0]) ) );
      hV2[ew]->Fill( 5.0, pt, TMath::Cos( 2.0*(phi-exEP[1]) ) );
      hV2[ew]->Fill( 6.0, pt, TMath::Cos( 2.0*(phi-exEP[2]) ) );
      hV2[ew]->Fill( 7.0, pt, TMath::Cos( 2.0*(phi-exEP[3]) ) );
      hV2[ew]->Fill( 8.0, pt, TMath::Cos( 2.0*(phi-exEP[4]) ) );
      hV2[ew]->Fill( 9.0, pt, TMath::Cos( 2.0*(phi-exEP[5]) ) );
      hV2[ew]->Fill( 10., pt, TMath::Cos( 2.0*(phi-exEP[6]) ) );
      hV2[ew]->Fill( 11., pt, TMath::Cos( 2.0*(phi-exEP[7]) ) );
      hV2[ew]->Fill( 12., pt, TMath::Cos( 2.0*(phi-exEP[8]) ) );
      if(pt>1.0&&pt<1.2) {
	hDP[ew]->Fill( 0.0, phi-bbEP[0] ); hDPPE[ew]->Fill( 0.0, phi-bbEPpe[0] );
	hDP[ew]->Fill( 1.0, phi-bbEP[1] ); hDPPE[ew]->Fill( 1.0, phi-bbEPpe[1] );
	hDP[ew]->Fill( 2.0, phi-bbEP[2] ); hDPPE[ew]->Fill( 2.0, phi-bbEPpe[2] );
	hDP[ew]->Fill( 3.0, phi-fvEP[0] ); hDPPE[ew]->Fill( 3.0, phi-fvEPpe[0] );
	hDP[ew]->Fill( 4.0, phi-exEP[0] ); hDPPE[ew]->Fill( 4.0, phi-exEPpe[0] );
	hDP[ew]->Fill( 5.0, phi-exEP[1] ); hDPPE[ew]->Fill( 5.0, phi-exEPpe[1] );
	hDP[ew]->Fill( 6.0, phi-exEP[2] ); hDPPE[ew]->Fill( 6.0, phi-exEPpe[2] );
	hDP[ew]->Fill( 7.0, phi-exEP[3] ); hDPPE[ew]->Fill( 7.0, phi-exEPpe[3] );
	hDP[ew]->Fill( 8.0, phi-exEP[4] ); hDPPE[ew]->Fill( 8.0, phi-exEPpe[4] );
	hDP[ew]->Fill( 9.0, phi-exEP[5] ); hDPPE[ew]->Fill( 9.0, phi-exEPpe[5] );
	hDP[ew]->Fill( 10., phi-exEP[6] ); hDPPE[ew]->Fill( 10., phi-exEPpe[6] );
	hDP[ew]->Fill( 11., phi-exEP[7] ); hDPPE[ew]->Fill( 11., phi-exEPpe[7] );
	hDP[ew]->Fill( 12., phi-exEP[8] ); hDPPE[ew]->Fill( 12., phi-exEPpe[8] );
      }
    }
    // END OF MAIN LOOP
    bbEPpe[0] = bbEP[0];
    bbEPpe[1] = bbEP[1];
    bbEPpe[2] = bbEP[2];
    fvEPpe[0] = fvEP[0];
    exEPpe[0] = exEP[0];
    exEPpe[1] = exEP[1];
    exEPpe[2] = exEP[2];
    exEPpe[3] = exEP[3];
    exEPpe[4] = exEP[4];
    exEPpe[5] = exEP[5];
    exEPpe[6] = exEP[6];
    exEPpe[7] = exEP[7];
    exEPpe[8] = exEP[8];
    caEPpe[0] = caEP[0];
  }
  //== SAVING EP CALIBS
  StoreEPFlatteningCoeficients();
  SaveCalibFiles(run);

  //======== SAVING OUTPUT
  TString ooo = Form("eventplane/myResults_%s.root",run.Data());
  std::cout << "SAVING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data(),"RECREATE"); 
  hEvents->Write();
  hFrac->Write();
  hVtxZ->Write();
  hCent->Write();

  hEp->Write();
  hChi2->Write();
  hDphi->Write();
  hphi->Write();
  hDz->Write();
  hZED->Write();
  hBeta->Write();

  hV2[0]->Write();
  hV2[1]->Write();
  hDP[0]->Write();
  hDP[1]->Write();
  hDPPE[0]->Write();
  hDPPE[1]->Write();
  
  SaveHistograms();
  
  f2->Close();
  std::cout << "DONE =)" << std::endl;
  return 0;
}
