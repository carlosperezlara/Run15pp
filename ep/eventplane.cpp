#include "../defs.h"
#include "eventplanedefs.h"
//bool fVerbosity = true;
bool fVerbosity = false;

//==================
void CreateQCA() {
  q2ca[0].Reset();
  uint ntrks = pTRKcid->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float pt   = TMath::Abs( pTRKpt->at(itrk) );
    float phi  = pTRKphi->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    float zed  = pTRKzed->at(itrk);
    int fbit = pTRKqua->at(itrk);
    //if(fbit!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    q2ca[0].Fill(phi);
  }
}
//==================
int main(int argc, char *argv[]){
  if(argc<3) {
    std::cout << "launch with args RUN NEVS" << std::endl;
    return 1;
  }
  q2ca[0].SetOrder(2);
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  //======= HISTOS
  TH1F *hEvents = new TH1F("Events","",10,-0.5,9.5);
  hEvents->GetXaxis()->SetBinLabel(1,"All");
  hEvents->GetXaxis()->SetBinLabel(2,"triggered");
  hEvents->GetXaxis()->SetBinLabel(3,"centrality");
  hEvents->GetXaxis()->SetBinLabel(4,"fraction");
  hEvents->GetXaxis()->SetBinLabel(5,"BB");
  hEvents->GetXaxis()->SetBinLabel(6,"FV");
  hEvents->GetXaxis()->SetBinLabel(7,"EX");
  hEvents->GetXaxis()->SetBinLabel(8,"CA");
  TH1F *hFrac = new TH1F("FRAC","", 100, 0.0, 1.0);
  TH1F *hCent = new TH1F("CENT","", 100, 0.0, 100);
  TH1F *hVtxZ = new TH1F("VTXZ","", 100, -11, +11);

  TH1F *hEp   = new TH1F("Ep","",   100, 0.0, 1.3);
  TH1F *hChi2 = new TH1F("Chi2","", 100, 0.0, 10.);
  TH1F *hDphi = new TH1F("pc3Dphi","", 100,-30, +30);
  TH1F *hphi  = new TH1F("phi","",  100, -TMath::TwoPi(), TMath::TwoPi());
  TH1F *hDz   = new TH1F("pc3Dz","", 100,-30., 30.);
  TH1F *hZED  = new TH1F("ZED","",  100,-100,+100);
  TH2F *hBeta = new TH2F("Beta",";p [GeV/c];[cm/ns]", 300,0,3,1800,-10,10);

  TH2F *hQxEX0 = new TH2F("QxEX0","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX1 = new TH2F("QxEX1","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX2 = new TH2F("QxEX2","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX3 = new TH2F("QxEX3","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX4 = new TH2F("QxEX4","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX5 = new TH2F("QxEX5","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX6 = new TH2F("QxEX6","",100,-30,+30,100,-30,+30);
  TH2F *hQxEX7 = new TH2F("QxEX7","",100,-30,+30,100,-30,+30);
  TH2F *hQxFV  = new TH2F("QxFV", "",100,-30,+30,100,-30,+30);

  TProfile2D *hV2[2];
  hV2[0] = new TProfile2D("hV2_E","", 13,-0.5,12.5, 30,0,3, -1.1,+1.1);
  hV2[1] = new TProfile2D("hV2_W","", 13,-0.5,12.5, 30,0,3, -1.1,+1.1);
  TH2F *hDP[2];
  hDP[0] = new TH2F("DPhiE", "", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());
  hDP[1] = new TH2F("DPhiW", "", 13,-0.5,12.5, 100,-TMath::TwoPi(),+TMath::TwoPi());
  TH2F *hDPPE[2];
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
    if(i1%50000 == 0) {
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
    isBBgood = false;
    isFVgood = false;
    isEXgood = false;
    isCAgood = false;
    if(q2bb[2].M()>40) isBBgood = true; 
    if(q2fv[0].M()>10) isFVgood = true; 
    float xm =0;
    for(int i=0; i!=8; ++i) xm += q2ex[i].M();
    if(xm>60) isEXgood = true; 
    if(q2ca[0].M()>2)  isCAgood = true; 

    if(isBBgood) hEvents->Fill(4.0);
    if(isFVgood) hEvents->Fill(5.0);
    if(isEXgood) hEvents->Fill(6.0);
    if(isCAgood) hEvents->Fill(7.0);

    /*
    std::cout << " || CA " << q2ca[0].NP();
    std::cout << " || BB " << q2bb[0].NP() << " " << q2bb[1].NP();
    std::cout << " || FV " << q2fv[0].NP();
    std::cout << " || EX " << q2ex[0].NP() << " " << q2ex[1].NP();
    std::cout << " " << q2ex[2].NP() << " " << q2ex[3].NP();
    std::cout << " " << q2ex[3].NP() << " " << q2ex[4].NP();
    std::cout << " " << q2ex[5].NP() << " " << q2ex[6].NP();
    std::cout << " " << q2ex[6].NP() << " " << q2ex[7].NP() << std::endl;
    */

    // histograms 
    hVtxZ->Fill( vtxZ );
    hCent->Fill( cent );
    hFrac->Fill( frac );

    //=================================
    //=================================
    // STEP 1:
    StoreQcomponents();
    hQxEX0->Fill( q2ex[0].X(), q2bb[2].X() );
    hQxEX1->Fill( q2ex[1].X(), q2bb[2].X() );
    hQxEX2->Fill( q2ex[2].X(), q2bb[2].X() );
    hQxEX3->Fill( q2ex[3].X(), q2bb[2].X() );
    hQxEX4->Fill( q2ex[4].X(), q2bb[2].X() );
    hQxEX5->Fill( q2ex[5].X(), q2bb[2].X() );
    hQxEX6->Fill( q2ex[6].X(), q2bb[2].X() );
    hQxEX7->Fill( q2ex[7].X(), q2bb[2].X() );
    hQxFV->Fill(  q2fv[0].X(), q2bb[2].X() );
    if(fVerbosity) std::cout << "===> STEP 1 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.0.cpp"
 
    //=================================
    //=================================
    // STEP 2: RECENTER
    if(isBBgood) Recenterbb();
    if(isFVgood) Recenterfv();
    if(isEXgood) Recenterex(); //only EX subdetectors
    if(isCAgood) Recenterca();
    if(fVerbosity) std::cout << "===> STEP 2 [DONE]" << endl;
    //=================================
    //=================================


    //=================================
    //=================================
    // STEP 3: CONSTRUCT FULL EX
    if(isEXgood) {
      q2ex[8] = q2ex[0] + q2ex[1] + q2ex[2] + q2ex[3]; // construct
      q2ex[8] = q2ex[4] + q2ex[5] + q2ex[6] + q2ex[7]; // full EX (after recenter its parts)
      //---- ( Qaing, save to remove later)
      hQx->Fill(12.,q2ex[8].X()); // Qx 12
      hQy->Fill(12.,q2ex[8].Y()); // Qy 12
      hQM->Fill( 12., q2ex[8].M() ); // Qred 12
      hQred0->Fill( 12., q2ex[8].Reduced() ); // Qred 12
      hPSI0->Fill( 12., q2ex[8].Psi2Pi() ); // PSI 12
      if(isBBgood)
	hEVC[0][13]->Fill( q2ex[8].Psi2Pi(), q2bb[2].Psi2Pi() ); // EVC 13
      if(isFVgood)
	hEVC[0][22]->Fill( q2ex[8].Psi2Pi(), q2fv[2].Psi2Pi() ); // EVC 22
      if(isCAgood)
	hEVC[0][31]->Fill( q2ex[8].Psi2Pi(), q2ca[2].Psi2Pi() ); // EVC 31
      //----
      q2ex[8].SetXY( (q2ex[8].X()-exqx[8][0])/*/exqx[8][1]*/,
		     (q2ex[8].Y()-exqy[8][0])/*/exqy[8][1]*/,
		     q2ex[8].NP(), q2ex[8].M() );
    }
    if(fVerbosity) std::cout << "===> STEP 3 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.1.cpp"

    //=================================
    //=================================
    // STEP 4: FLATTENING
    //std::cout << bbEP[0] << " " << bbEP[1] << " " << bbEP[2] << " | ";
    if(isBBgood) GetEPbb();
    if(isFVgood) GetEPfv();
    if(isEXgood) GetEPex();
    if(isCAgood) GetEPca();
    //std::cout << bbEP[0] << " " << bbEP[1] << " " << bbEP[2] << std::endl;

    if(fVerbosity) std::cout << "===> STEP 4 [DONE]" << endl;
    //=================================
    //=================================

    #include "eventplane.qa.2.cpp"
    if(isBBgood) StoreEPFlatteningCoeficientsBB();
    if(isFVgood) StoreEPFlatteningCoeficientsFV();
    if(isEXgood) StoreEPFlatteningCoeficientsEX();
    if(isCAgood) StoreEPFlatteningCoeficientsCA();

    //=== START OF MAIN LOOP
    if(fVerbosity) std::cout << "===> START OF MAIN LOOP" << endl;
    uint ntrks = pTRKcid->size();
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
      double dlen = pTRKplemc->at(itrk);
      double time = tof+dlen/kSpeedOfLightinNS;
      double beta = dlen / ( time * kSpeedOfLightinNS );
      if(qua!=63) continue;
      if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
      if(TMath::Abs(dphi)>3) continue;
      if(TMath::Abs(dz)>3) continue;

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
      if(isBBgood) {
	hV2[ew]->Fill( 0.0, pt, TMath::Cos( 2.0*(phi-bbEP[0]) ) );
	hV2[ew]->Fill( 1.0, pt, TMath::Cos( 2.0*(phi-bbEP[1]) ) );
	hV2[ew]->Fill( 2.0, pt, TMath::Cos( 2.0*(phi-bbEP[2]) ) );
      }
      if(isFVgood)
	hV2[ew]->Fill( 3.0, pt, TMath::Cos( 2.0*(phi-fvEP[0]) ) );
      if(isEXgood) {
	hV2[ew]->Fill( 4.0, pt, TMath::Cos( 2.0*(phi-exEP[0]) ) );
	hV2[ew]->Fill( 5.0, pt, TMath::Cos( 2.0*(phi-exEP[1]) ) );
	hV2[ew]->Fill( 6.0, pt, TMath::Cos( 2.0*(phi-exEP[2]) ) );
	hV2[ew]->Fill( 7.0, pt, TMath::Cos( 2.0*(phi-exEP[3]) ) );
	hV2[ew]->Fill( 8.0, pt, TMath::Cos( 2.0*(phi-exEP[4]) ) );
	hV2[ew]->Fill( 9.0, pt, TMath::Cos( 2.0*(phi-exEP[5]) ) );
	hV2[ew]->Fill( 10., pt, TMath::Cos( 2.0*(phi-exEP[6]) ) );
	hV2[ew]->Fill( 11., pt, TMath::Cos( 2.0*(phi-exEP[7]) ) );
	hV2[ew]->Fill( 12., pt, TMath::Cos( 2.0*(phi-exEP[8]) ) );
      }
      if(pt>1.0&&pt<1.2) {
	if(isBBgood) {
	  hDP[ew]->Fill( 0.0, phi-bbEP[0] ); hDPPE[ew]->Fill( 0.0, phi-bbEPpe[0] );
	  hDP[ew]->Fill( 1.0, phi-bbEP[1] ); hDPPE[ew]->Fill( 1.0, phi-bbEPpe[1] );
	  hDP[ew]->Fill( 2.0, phi-bbEP[2] ); hDPPE[ew]->Fill( 2.0, phi-bbEPpe[2] );
	}
	if(isFVgood) {
	  hDP[ew]->Fill( 3.0, phi-fvEP[0] ); hDPPE[ew]->Fill( 3.0, phi-fvEPpe[0] );
	}
	if(isEXgood) {
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
    }
    // END OF MAIN LOOP
    if(isBBgood) {
      bbEPpe[0] = bbEP[0];
      bbEPpe[1] = bbEP[1];
      bbEPpe[2] = bbEP[2];
    }
    if(isFVgood) fvEPpe[0] = fvEP[0];
    if(isEXgood) {
      exEPpe[0] = exEP[0];
      exEPpe[1] = exEP[1];
      exEPpe[2] = exEP[2];
      exEPpe[3] = exEP[3];
      exEPpe[4] = exEP[4];
      exEPpe[5] = exEP[5];
      exEPpe[6] = exEP[6];
      exEPpe[7] = exEP[7];
      exEPpe[8] = exEP[8];
    }
    if(isCAgood) caEPpe[0] = caEP[0];
  }
  //== SAVING EP CALIBS
  SaveCalibFiles(run);

  //======== SAVING OUTPUT
  TString ooo = Form("out/myResults_%s.root",run.Data());
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

  hQxEX0->Write();
  hQxEX1->Write();
  hQxEX2->Write();
  hQxEX3->Write();
  hQxEX4->Write();
  hQxEX5->Write();
  hQxEX6->Write();
  hQxEX7->Write();
  hQxFV->Write();

  SaveHistograms();
  
  f2->Close();
  std::cout << "DONE =)" << std::endl;
  return 0;
}
