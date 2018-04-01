const int kMaxOrder = 36;

//BBC 3 detectors
const int kBBN = 3;
double bba[kBBN][kMaxOrder];
double bbb[kBBN][kMaxOrder];
double bbqx[kBBN][2], bbqy[kBBN][2];
double bbEP[kBBN];
double bbEPpe[kBBN];
qcQ q2bb[kBBN];

// MPCEX 9 detectors
const int kEXN = 9;
double exa[kEXN][kMaxOrder];
double exb[kEXN][kMaxOrder];
double exqx[kEXN][2], exqy[kEXN][2];
double exEP[kEXN];
double exEPpe[kEXN];
qcQ q2ex[kEXN];

// FVTX 1 detector
const int kFVN = 1;
double fva[kFVN][kMaxOrder];
double fvb[kFVN][kMaxOrder];
double fvqx[kFVN][2], fvqy[kFVN][2];
double fvEP[kFVN];
double fvEPpe[kFVN];
qcQ q2fv[kFVN];

// CENTRAL ARM 1 detectors
const int kCAN = 1;
double caa[kCAN][kMaxOrder];
double cab[kCAN][kMaxOrder];
double caqx[kCAN][2], caqy[kCAN][2];
double caEP[kCAN];
double caEPpe[kCAN];
qcQ q2ca[kCAN];

TH2F *hQx;
TH2F *hQy;
TH2F *hQM;
TProfile2D *hPsiC;
TProfile2D *hPsiS;
  
TH2F *hBBCX;
TH2F *hBBCY;
TH2F *hBBCL;
TH2F *hBBCH;
TH2F *hBBCS;
TH2F *hBBCM;
TH2F *hBBCNP;
TH2F *hBBCatan01;
TH2F *hBBCatan02;
TH2F *hBBCatan12;
TH2F *hBBCcent;

TH2F *hQred0;
TH2F *hQred1;
TH2F *hPSI0;
TH2F *hPSI1;
TH2F *hPSI2;

TH2F *hResD;
TProfile *hRes;

TH2F *hEVC[3][32];

//=====================================
void CreateHistograms() {
  // BB3 FV1 EX9 CA1 =14
  TString tit14[14] = {"BBA","BBB","BB", "FV", "EX0","EX1","EX2",
		       "EX3","EX4","EX5","EX6","EX7","EX8","CA"};
  hQx = new TH2F("Qx","", 14,-0.5,13.5, 500,-100,+100);
  hQy = new TH2F("Qy","", 14,-0.5,13.5, 500,-100,+100);
  hQM = new TH2F("QM","", 14,-0.5,13.5, 500,0,+1000);
  hPsiC = new TProfile2D("PsiC","", 14,-0.5,13.5, kMaxOrder,0.5,kMaxOrder+0.5, -1.1,+1.1);
  hPsiS = new TProfile2D("PsiS","", 14,-0.5,13.5, kMaxOrder,0.5,kMaxOrder+0.5, -1.1,+1.1);
  hQred0= new TH2F("Qred0","",14,-0.5,13.5, 300, 0.0, 15.0);
  hQred1= new TH2F("Qred1","",14,-0.5,13.5, 300, 0.0, 15.0);
  hPSI0 = new TH2F("PSI0","", 14,-0.5,13.5, 100, 0.0, TMath::TwoPi());
  hPSI1 = new TH2F("PSI1","", 14,-0.5,13.5, 100, 0.0, TMath::TwoPi());
  hPSI2 = new TH2F("PSI2","", 14,-0.5,13.5, 100, 0.0, TMath::TwoPi());
  for(int i=0; i!=14; ++i) {
    hQx->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hQy->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hQM->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hPsiC->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hPsiS->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hQred0->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hQred1->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hPSI0->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hPSI1->GetXaxis()->SetBinLabel(i+1,tit14[i]);
    hPSI2->GetXaxis()->SetBinLabel(i+1,tit14[i]);
  }
    
  hBBCX = new TH2F("BBCX","", 100,-30,30, 100,-30,30);
  hBBCY = new TH2F("BBCY","", 100,-30,30, 100,-30,30);
  hBBCL = new TH2F("BBCL","", 100,0,30, 100,0,30);
  hBBCH = new TH2F("BBCH","", 100,0,3, 100,0,3);
  hBBCS = new TH2F("BBCS","", 100,-3,3, 100,-3,3);
  hBBCM = new TH2F("BBCM","", 100,0,150, 100,0,150);
  hBBCNP= new TH2F("BBCNP","", 64,-0.5,63.6, 64,-0.5,63.5);
  hBBCatan01 = new TH2F("BBCatan01","", 100,-4,+4, 100,-4,+4);
  hBBCatan02 = new TH2F("BBCatan02","", 100,-4,+4, 100,-4,+4);
  hBBCatan12 = new TH2F("BBCatan12","", 100,-4,+4, 100,-4,+4);
  hBBCcent = new TH2F("BBCcent","", 100,0,100, 100,0,300);
  
  // BBAC BBCBC BBAB 
  // FVBB FVCA
  // EX0BB EX0CA EX0FV (x9)
  // = 32
  TString title32[32] = {";BBSubA;BB", ";BBSubB;BB", ";BBSubA;BBSubB",
			 ";FV;BB",  ";FV;CA",
			 ";EX0;BB", ";EX1;BB", ";EX2;BB", ";EX3;BB",
			 ";EX4;BB", ";EX5;BB", ";EX6;BB", ";EX7;BB", ";EX8;BB"
			 ";EX0;FV", ";EX1;FV", ";EX2;FV", ";EX3;FV",
			 ";EX4;FV", ";EX5;FV", ";EX6;FV", ";EX7;FV", ";EX8;FV"
			 ";EX0;CA", ";EX1;CA", ";EX2;CA", ";EX3;CA",
			 ";EX4;CA", ";EX5;CA", ";EX6;CA", ";EX7;CA", ";EX8;CA"};
  for(int i=0; i!=3; ++i)
    for(int j=0; j!=32; ++j)
      hEVC[i][j] = new TH2F( Form("EVC%d_DET%d",i,j), title32[j].Data(),
			     100,0,TMath::TwoPi(), 100,0,TMath::TwoPi() );

  // 32 +
  // EX04 EX15 EX26 EX37
  // = 36
  TString title36[36] = {"BBA-BB","BBB-BB","BBA-BBB","FB-BB","FV-CA",
			 "EX0-BB","EX1-BB","EX2-BB","EX3-BB","EX4-BB","EX5-BB","EX6-BB","EX7-BB","EX8-BB",
			 "EX0-FV","EX1-FV","EX2-FV","EX3-FV","EX4-FV","EX5-FV","EX6-FV","EX7-FV","EX8-FV",
			 "EX0-CA","EX1-CA","EX2-CA","EX3-CA","EX4-CA","EX5-CA","EX6-CA","EX7-CA","EX8-CA",
			 "EX0-EX4","EX1-EX5","EX2-EX6","EX3-EX7"};
  hResD = new TH2F("ResD","", 36,-0.5,35.5, 100,-1.1,+1.1);
  hRes = new TProfile("Res","", 36,-0.5,35.5, -1.1,+1.1);  
  for(int i=0; i!=36; ++i) {
    hResD->GetXaxis()->SetBinLabel(i+1,title36[i]);
    hRes->GetXaxis()->SetBinLabel(i+1,title36[i]);
  }
}
//=====================================
void SaveHistograms() {
  hBBCX->Write();
  hBBCY->Write();
  hBBCL->Write();
  hBBCH->Write();
  hBBCS->Write();
  hBBCM->Write();
  hBBCNP->Write();
  hBBCcent->Write();
  hBBCatan01->Write();
  hBBCatan02->Write();
  hBBCatan12->Write();

  hQx->Write();
  hQy->Write();
  hQM->Write();
  hPsiC->Write();
  hPsiS->Write();

  hQred0->Write();
  hQred1->Write();
  
  hPSI0->Write();
  hPSI1->Write();
  hPSI2->Write();
  
  hRes->Write();
  hResD->Write();
  for(int i=0; i!=3; ++i) for(int j=0; j!=32; ++j) hEVC[i][j]->Write();
}
//=====================================
void LoadCalibs(TString run) {
  for(int idet=0; idet!=kBBN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      bba[idet][iord] = 0.0;
      bbb[idet][iord] = 0.0;
    }
    bbqx[idet][0] = bbqx[idet][1] = 0.0;
    bbqy[idet][0] = bbqy[idet][1] = 0.0;
  }
  for(int idet=0; idet!=kFVN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      fva[idet][iord] = 0.0;
      fvb[idet][iord] = 0.0;
    }
    fvqx[idet][0] = fvqx[idet][1] = 0.0;
    fvqy[idet][0] = fvqy[idet][1] = 0.0;
  }
  for(int idet=0; idet!=kEXN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      exa[idet][iord] = 0.0;
      exb[idet][iord] = 0.0;
    }
    exqx[idet][0] = exqx[idet][1] = 0.0;
    exqy[idet][0] = exqy[idet][1] = 0.0;
  }
  for(int idet=0; idet!=kCAN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      caa[idet][iord] = 0.0;
      cab[idet][iord] = 0.0;
    }
    caqx[idet][0] = caqx[idet][1] = 0.0;
    caqy[idet][0] = caqy[idet][1] = 0.0;
  }

  int tmp;
  ifstream fin;
  //========
  fin.open( Form("out/calib/bb%s.dat1",run.Data()) );
  fin >> tmp; // ndet
  if(fin.good()) {
    for(int idet=0; idet!=kBBN; ++idet) {
      fin >> bbqx[idet][0] >> bbqy[idet][0];
      fin >> bbqx[idet][1] >> bbqy[idet][1];
    }
  }
  fin.close();
  fin.open( Form("out/calib/fv%s.dat1",run.Data()) );
  fin >> tmp; // ndet
  if(fin.good()) {
    for(int idet=0; idet!=kFVN; ++idet) {
      fin >> fvqx[idet][0] >> fvqy[idet][0];
      fin >> fvqx[idet][1] >> fvqy[idet][1];
    }
  }
  fin.close();
  fin.open( Form("out/calib/ca%s.dat1",run.Data()) );
  fin >> tmp; // ndet
  if(fin.good()) {
    for(int idet=0; idet!=kCAN; ++idet) {
      fin >> caqx[idet][0] >> caqy[idet][0];
      fin >> caqx[idet][1] >> caqy[idet][1];
    }
  }
  fin.close();
  fin.open( Form("out/calib/ex%s.dat1",run.Data()) );
  fin >> tmp; // ndet
  if(fin.good()) {
    for(int idet=0; idet!=kEXN; ++idet) {
      fin >> exqx[idet][0] >> exqy[idet][0];
      fin >> exqx[idet][1] >> exqy[idet][1];
    }
  }
  fin.close();
  //========
  fin.open( Form("out/calib/bb%s.dat2",run.Data()) );
  fin >> tmp; //ord
  if(fin.good()) {
    for(int idet=0; idet!=kBBN; ++idet) {
      for(int iord=0; iord!=kMaxOrder; ++iord) {
	fin >> bba[idet][iord] >> bbb[idet][iord];
      }
    }
  }
  fin.close();
  fin.open( Form("out/calib/fv%s.dat2",run.Data()) );
  fin >> tmp; //ord
  if(fin.good()) {
    for(int idet=0; idet!=kFVN; ++idet) {
      for(int iord=0; iord!=kMaxOrder; ++iord) {
	fin >> fva[idet][iord] >> fvb[idet][iord];
      }
    }
  }
  fin.close();
  fin.open( Form("out/calib/ca%s.dat2",run.Data()) );
  fin >> tmp; //ord
  if(fin.good()) {
    for(int idet=0; idet!=kCAN; ++idet) {
      for(int iord=0; iord!=kMaxOrder; ++iord) {
	fin >> caa[idet][iord] >> cab[idet][iord];
      }
    }
  }
  fin.close();
  fin.open( Form("out/calib/ex%s.dat2",run.Data()) );
  fin >> tmp; //ord
  if(fin.good()) {
    for(int idet=0; idet!=kEXN; ++idet) {
      for(int iord=0; iord!=kMaxOrder; ++iord) {
	fin >> exa[idet][iord] >> exb[idet][iord];
      }
    }
  }
  fin.close();
}
//=====================================
void GetEPbb() {
  for(int idet=0; idet!=kBBN; ++idet) {
    bbEP[idet] = q2bb[idet].Psi2Pi();
    double delta = 0.0;
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      delta += TMath::Cos((iord+1.0)*bbEP[idet])*bbb[idet][iord];
      delta += TMath::Sin((iord+1.0)*bbEP[idet])*bba[idet][iord];
    }
    bbEP[idet] += delta;
  }
}
//=====================================
void GetEPfv() {
  for(int idet=0; idet!=kFVN; ++idet) {
    fvEP[idet] = q2fv[idet].Psi2Pi();
    double delta = 0.0;
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      delta += TMath::Cos((iord+1.0)*fvEP[idet])*fvb[idet][iord];
      delta += TMath::Sin((iord+1.0)*fvEP[idet])*fva[idet][iord];
    }
    fvEP[idet] += delta;
  }
}
//=====================================
void GetEPca() {
  for(int idet=0; idet!=kCAN; ++idet) {
    caEP[idet] = q2ca[idet].Psi2Pi();
    double delta = 0.0;
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      delta += TMath::Cos((iord+1.0)*caEP[idet])*cab[idet][iord];
      delta += TMath::Sin((iord+1.0)*caEP[idet])*caa[idet][iord];
    }
    caEP[idet] += delta;
  }
}
//=====================================
void GetEPex() {
  for(int idet=0; idet!=kEXN; ++idet) {
    exEP[idet] = q2ex[idet].Psi2Pi();
    double delta = 0.0;
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      delta += TMath::Cos((iord+1.0)*exEP[idet])*exb[idet][iord];
      delta += TMath::Sin((iord+1.0)*exEP[idet])*exa[idet][iord];
    }
    exEP[idet] += delta;
  }
}
//=====================================
void Recenterbb() {
  for(int idet=0; idet!=kBBN; ++idet) {
    q2bb[idet].SetXY( (q2bb[idet].X()-bbqx[idet][0])/*/bbqx[idet][1]*/,
		      (q2bb[idet].Y()-bbqy[idet][0])/*/bbqy[idet][1]*/,
		      q2bb[idet].NP(),
		      q2bb[idet].M() );
  }
}
//=====================================
void Recenterfv() {
  for(int idet=0; idet!=kFVN; ++idet) {
    q2fv[idet].SetXY( (q2fv[idet].X()-fvqx[idet][0])/*/fvqx[idet][1]*/,
		      (q2fv[idet].Y()-fvqy[idet][0])/*/fvqy[idet][1]*/,
		      q2fv[idet].NP(),
		      q2fv[idet].M() );
  }
}
//=====================================
void Recenterca() {
  for(int idet=0; idet!=kCAN; ++idet) {
    q2ca[idet].SetXY( (q2ca[idet].X()-caqx[idet][0])/*/caqx[idet][1]*/,
		      (q2ca[idet].Y()-caqy[idet][0])/*/caqy[idet][1]*/,
		      q2ca[idet].NP(),
		      q2ca[idet].M() );
  }
}
//=====================================
void Recenterex() {
  for(int idet=0; idet!=kEXN-1; ++idet) { // on purpose do not touch full EX
    q2ex[idet].SetXY( (q2ex[idet].X()-exqx[idet][0])/*/exqx[idet][1]*/,
		      (q2ex[idet].Y()-exqy[idet][0])/*/exqy[idet][1]*/,
		      q2ex[idet].NP(),
		      q2ex[idet].M() );
  }
}
//=====================================
void StoreQcomponents() {
  hQx->Fill(0.0,q2bb[0].X());
  hQx->Fill(1.0,q2bb[1].X());
  hQx->Fill(2.0,q2bb[2].X());
  hQx->Fill(3.0,q2fv[0].X());
  hQx->Fill(4.0,q2ex[0].X());
  hQx->Fill(5.0,q2ex[1].X());
  hQx->Fill(6.0,q2ex[2].X());
  hQx->Fill(7.0,q2ex[3].X());
  hQx->Fill(8.0,q2ex[4].X());
  hQx->Fill(9.0,q2ex[5].X());
  hQx->Fill(10.,q2ex[6].X());
  hQx->Fill(11.,q2ex[7].X());
  //
  hQx->Fill(13.,q2ca[0].X());

  hQy->Fill(0.0,q2bb[0].Y());
  hQy->Fill(1.0,q2bb[1].Y());
  hQy->Fill(2.0,q2bb[2].Y());
  hQy->Fill(3.0,q2fv[0].Y());
  hQy->Fill(4.0,q2ex[0].Y());
  hQy->Fill(5.0,q2ex[1].Y());
  hQy->Fill(6.0,q2ex[2].Y());
  hQy->Fill(7.0,q2ex[3].Y());
  hQy->Fill(8.0,q2ex[4].Y());
  hQy->Fill(9.0,q2ex[5].Y());
  hQy->Fill(10.,q2ex[6].Y());
  hQy->Fill(11.,q2ex[7].Y());
  //
  hQy->Fill(13.,q2ca[0].Y());

  // Qx 12 not stored yet
  // Qy 12 not stored yet
}
//=====================================
void StoreEPFlatteningCoeficients() {
  for(int iord=0; iord!=kMaxOrder; ++iord) {
    float nn = iord + 1.0;
    hPsiC->Fill(0.0, nn, TMath::Cos(nn*q2bb[0].Psi2Pi()) );
    hPsiC->Fill(1.0, nn, TMath::Cos(nn*q2bb[1].Psi2Pi()) );
    hPsiC->Fill(2.0, nn, TMath::Cos(nn*q2bb[2].Psi2Pi()) );
    hPsiC->Fill(3.0, nn, TMath::Cos(nn*q2fv[0].Psi2Pi()) );
    hPsiC->Fill(4.0, nn, TMath::Cos(nn*q2ex[0].Psi2Pi()) );
    hPsiC->Fill(5.0, nn, TMath::Cos(nn*q2ex[1].Psi2Pi()) );
    hPsiC->Fill(6.0, nn, TMath::Cos(nn*q2ex[2].Psi2Pi()) );
    hPsiC->Fill(7.0, nn, TMath::Cos(nn*q2ex[3].Psi2Pi()) );
    hPsiC->Fill(8.0, nn, TMath::Cos(nn*q2ex[4].Psi2Pi()) );
    hPsiC->Fill(9.0, nn, TMath::Cos(nn*q2ex[5].Psi2Pi()) );
    hPsiC->Fill(10., nn, TMath::Cos(nn*q2ex[6].Psi2Pi()) );
    hPsiC->Fill(11., nn, TMath::Cos(nn*q2ex[7].Psi2Pi()) );
    hPsiC->Fill(12., nn, TMath::Cos(nn*q2ex[8].Psi2Pi()) );
    hPsiC->Fill(13., nn, TMath::Cos(nn*q2ca[0].Psi2Pi()) );
    //=
    hPsiS->Fill(0.0, nn, TMath::Sin(nn*q2bb[0].Psi2Pi()) );
    hPsiS->Fill(1.0, nn, TMath::Sin(nn*q2bb[1].Psi2Pi()) );
    hPsiS->Fill(2.0, nn, TMath::Sin(nn*q2bb[2].Psi2Pi()) );
    hPsiS->Fill(3.0, nn, TMath::Sin(nn*q2fv[0].Psi2Pi()) );
    hPsiS->Fill(4.0, nn, TMath::Sin(nn*q2ex[0].Psi2Pi()) );
    hPsiS->Fill(5.0, nn, TMath::Sin(nn*q2ex[1].Psi2Pi()) );
    hPsiS->Fill(6.0, nn, TMath::Sin(nn*q2ex[2].Psi2Pi()) );
    hPsiS->Fill(7.0, nn, TMath::Sin(nn*q2ex[3].Psi2Pi()) );
    hPsiS->Fill(8.0, nn, TMath::Sin(nn*q2ex[4].Psi2Pi()) );
    hPsiS->Fill(9.0, nn, TMath::Sin(nn*q2ex[5].Psi2Pi()) );
    hPsiS->Fill(10., nn, TMath::Sin(nn*q2ex[6].Psi2Pi()) );
    hPsiS->Fill(11., nn, TMath::Sin(nn*q2ex[7].Psi2Pi()) );
    hPsiS->Fill(12., nn, TMath::Sin(nn*q2ex[8].Psi2Pi()) );
    hPsiS->Fill(13., nn, TMath::Sin(nn*q2ca[0].Psi2Pi()) );
  }
}
//=====================================
void MoveEPtoPreviousEvent() {
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
//=====================================
void SaveCalibFiles(TString run) {
  ofstream fout;
  TH1D *hpx, *hpy;
  //==
  fout.open( Form("out/calib/bb%s.dat1",run.Data()) );
  fout << kBBN << endl;
  for(int idet=0; idet!=kBBN; ++idet) {
    hpx = hQx->ProjectionY( Form("qx_bb%d",idet), 1+idet, 1+idet );
    hpy = hQy->ProjectionY( Form("qy_bb%d",idet), 1+idet, 1+idet );
    fout << hpx->GetMean() << " " << hpy->GetMean() << " ";
    fout << hpx->GetRMS()  << " " << hpy->GetRMS() << endl;
  }
  fout.close();
  //==
  fout.open( Form("out/calib/fv%s.dat1",run.Data()) );
  fout << kFVN << endl;
  hpx = hQx->ProjectionY( "qx_fv", 4, 4 );
  hpy = hQy->ProjectionY( "qy_fv", 4, 4 );
  fout << hpx->GetMean() << " " << hpy->GetMean() << " ";
  fout << hpx->GetRMS()  << " " << hpy->GetRMS() << endl;
  fout.close();
  //==
  fout.open( Form("out/calib/ex%s.dat1",run.Data()) );
  fout << kEXN << endl;
  for(int idet=0; idet!=9; ++idet) {
    hpx = hQx->ProjectionY( Form("qx_ex%d",idet), 5+idet, 5+idet );
    hpy = hQy->ProjectionY( Form("qy_ex%d",idet), 5+idet, 5+idet );
    fout << hpx->GetMean() << " " << hpy->GetMean() << " ";
    fout << hpx->GetRMS()  << " " << hpy->GetRMS() << endl;
  }
  fout.close();
  fout.open( Form("out/calib/ca%s.dat1",run.Data()) );
  fout << kCAN << endl;
  for(int idet=0; idet!=kCAN; ++idet) {
    hpx = hQx->ProjectionY( Form("qx_ca%d",idet), 14+idet, 14+idet );
    hpy = hQy->ProjectionY( Form("qy_ca%d",idet), 14+idet, 14+idet );
    fout << hpx->GetMean() << " " << hpy->GetMean() << " ";
    fout << hpx->GetRMS()  << " " << hpy->GetRMS() << endl;
  }
  fout.close();

  //=======

  fout.open( Form("out/calib/bb%s.dat2",run.Data()) );
  fout << kMaxOrder << endl;
  for(int idet=0; idet!=kBBN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      fout << +2.0/(iord+1.0)*hPsiC->GetBinContent(idet+1,iord+1) << " ";
      fout << -2.0/(iord+1.0)*hPsiS->GetBinContent(idet+1,iord+1) << " ";
    }
    fout << endl;
  }
  fout.close();
  //==
  fout.open( Form("out/calib/fv%s.dat2",run.Data()) );
  fout << kMaxOrder << endl;
  for(int iord=0; iord!=kMaxOrder; ++iord) {
    fout << +2.0/(iord+1.0)*hPsiC->GetBinContent(4,iord+1) << " ";
    fout << -2.0/(iord+1.0)*hPsiS->GetBinContent(4,iord+1) << " ";
  }
  fout << endl;
  fout.close();
  //==
  fout.open( Form("out/calib/ex%s.dat2",run.Data()) );
  fout << kMaxOrder << endl;
  for(int idet=0; idet!=kEXN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      fout << +2.0/(iord+1.0)*hPsiC->GetBinContent(idet+5,iord+1) << " ";
      fout << -2.0/(iord+1.0)*hPsiS->GetBinContent(idet+5,iord+1) << " ";
    }
    fout << endl;
  }
  fout.close();
  //==
  fout.open( Form("out/calib/ca%s.dat2",run.Data()) );
  fout << kMaxOrder << endl;
  for(int idet=0; idet!=kCAN; ++idet) {
    for(int iord=0; iord!=kMaxOrder; ++iord) {
      fout << +2.0/(iord+1.0)*hPsiC->GetBinContent(idet+14,iord+1) << " ";
      fout << -2.0/(iord+1.0)*hPsiS->GetBinContent(idet+14,iord+1) << " ";
    }
    fout << endl;
  }
  fout.close();
  std::cout << "Calibration files have been saved." << std::endl;
  //==
}
