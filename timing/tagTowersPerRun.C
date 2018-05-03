TH2F *hHitDensity[8];
TH2F *hHitDensityFilter[5][8];
TProfile2D *hHitMap[8];
int twr2iz[24768][4];

int secs[9] = { 0, 2592, 5184, 7776, 10368,
		12960, 15552, 20160, 24768 };
int isbad[24768];
float ddd[8] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2 };

void Load() {
  ifstream badt("twid.dat");
  int tid, sec, iy, iz, bd;
  for(int i=0; i!=24768; ++i) {
    twr2iz[i][0] = 0;
    twr2iz[i][1] = 0;
    twr2iz[i][2] = 0;
  }
  for(;;) {
    badt >> tid;
    if(!badt.good()) break;
    badt >> sec >> iy >> iz;
    twr2iz[tid][0] = iy;
    twr2iz[tid][1] = iz;
    twr2iz[tid][2] = sec;
  }
  badt.close();
}

void LoadHitDensity() {
  TFile *map = new TFile("hitdensity.root");
  for(int i=0; i!=8; ++i)
    hHitDensity[i] = (TH2F*) map->Get( Form("HitDensity_SEC%d",i) );
}

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

void CreateHitDensity() {
  LoadBadTowers();
  for(int i=0; i!=8; ++i) {
    hHitDensity[i] = new TH2F( Form("HitDensity_SEC%d",i),"",
			       secs[i+1]-secs[i], secs[i]-0.5, secs[i+1]-0.5,
			       100,0.0,ddd[i]);
    for(int j=0; j!=5; ++j) {
      hHitDensityFilter[j][i] = new TH2F( Form("HitDensityFilter%d_SEC%d",j,i),"",
					  secs[i+1]-secs[i], secs[i]-0.5, secs[i+1]-0.5,
					  100,0.0,ddd[i]);
    }
    hHitMap[i] = new TProfile2D( Form("HitMap_SEC%d",i),"",
				 100,-0.5,99.5, 100,-0.5,99.5);
  }

  ifstream frun( "gruns.dat" );
  TFile *file;
  int run;
  for(;;) {
    frun >> run;
    if(!frun.good()) break;
    cout << "processing " << run << endl;
    file = new TFile( Form("out/MyResults_%d.root",run) );
    TH1F *evt = (TH1F*) file->Get("Events");
    TH1F *hit = (TH1F*) file->Get("HitTwr");
    double evts = evt->GetBinContent(1);
    hit->Scale( 1000.0/evts );
    for(int bin=0; bin!=hit->GetXaxis()->GetNbins(); ++bin) {
      double twr = hit->GetXaxis()->GetBinCenter(bin+1);
      int itwr = int(twr+0.01);
      double hden = hit->GetBinContent(bin+1);
      int isec;
      if(twr<secs[1]) isec = 0;
      else if(twr<secs[2]) isec = 1;
      else if(twr<secs[3]) isec = 2;
      else if(twr<secs[4]) isec = 3;
      else if(twr<secs[5]) isec = 4;
      else if(twr<secs[6]) isec = 5;
      else if(twr<secs[7]) isec = 6;
      else if(twr<secs[8]) isec = 7;
      hHitDensity[isec]->Fill(twr,hden);
      int iy = twr2iz[twr][0];
      int iz = twr2iz[twr][1];
      hHitMap[isec]->Fill(iy,iz+50,hden);
      if( isbad[itwr]==0 ) {
	hHitMap[isec]->Fill(iy,iz,hden);
      }
      hHitDensityFilter[ isbad[itwr] ][isec]->Fill(twr,hden);
    }
  }
  frun.close();
  TFile *map = new TFile("hitdensity.root","RECREATE");
  for(int i=0; i!=8; ++i) {
    hHitDensity[i]->Write();
    for(int j=0; j!=4; ++j) {
      hHitDensityFilter[j][i]->Write();
    }
    hHitMap[i]->Write();
  }
  map->Close();
}

int CreateBadList() {
  LoadHitDensity();

  float histMean[8];
  float meanGaus[8];
  float sigmaGaus[8];

  int tagged[8][5];
  for(int i=0; i!=8; ++i)
    for(int j=0; j!=5; ++j)
      tagged[i][j] = 0;

  ofstream badt("badtowers.dat");

  for(int sss=0; sss!=6; ++sss) {
    cout << "========== S E C T O R   " << sss << endl;
    for(int b=0; b!=hHitDensity[sss]->GetNbinsX(); ++b) {
      TH1D* binh = hHitDensity[sss]->ProjectionY( Form("SEC%d_%d",sss,b), b+1, b+1 );
      double ent = binh->GetEntries();
      double mea = binh->GetMean();
      double rms = binh->GetRMS();
      double inte= binh->Integral();
      int twr =  secs[sss] + b;
      int iz = twr2iz[twr][1];
      int iy = twr2iz[twr][0];
      int clas=0;
      /*if(iz==0||iz==1||iy==0||iy==1||
	 iz==34||iz==35||iy==70||iy==71) clas=4;
	 else*/ if(mea<0.10) clas = 1;
      else if(inte<88||mea>0.4) clas = 2;
      else if(rms>0.035) clas=3;
      if(clas>0) {
	badt << twr << " " << clas << endl;
	if(clas<4) {
	  cout << " tower " << twr << " bin " << b+1 << " ";
	  cout << mea << " " << rms << " " << inte << " ===> " << clas << endl;
	}
      }
      tagged[sss][clas]++;
    }
  }

  for(int sss=6; sss!=8; ++sss) {
    cout << "========== S E C T O R   " << sss << endl;
    for(int b=0; b!=hHitDensity[sss]->GetNbinsX(); ++b) {
      TH1D* binh = hHitDensity[sss]->ProjectionY( Form("SEC%d_%d",sss,b), b+1, b+1 );
      double ent = binh->GetEntries();
      double mea = binh->GetMean();
      double rms = binh->GetRMS();
      double inte= binh->Integral();
      int twr =  secs[sss] + b;
      int iz = twr2iz[twr][1];
      int iy = twr2iz[twr][0];
      int clas=0;
      /*if(iz==0||iz==1||iy==0||iy==1||
	 iz==46||iz==47||iy==94||iy==95) clas=4;
	 else*/ if(mea<0.015) clas = 1;
      else if(inte<88||mea>0.14) clas = 2;
      else if(rms>0.02/*||rms<0.0040*/) clas=3;
      if(clas>0) {
	badt << twr << " " << clas << endl;
	if(clas<4) {
	  cout << " tower " << twr << " bin " << b+1 << " ";
	  cout << mea << " " << rms << " " << inte << " ===> " << clas << endl;
	}
      }
      tagged[sss][clas]++;
    }
  }

  badt.close();

  //==========================
  for(int i=0; i!=8; ++i) {
    cout << "SECTOR " << i << " || ";
    for(int j=0; j!=5; ++j) {
      cout << Form("%.2f ", 100*double(tagged[i][j]) / double(secs[i+1]-secs[i])) << " ";
    }
    cout << endl;
  }

  return 0;
}

void CompareToJaehyeon() {
  CompareTo("Run16dAu200WarnMap_ERT.list");
}
void CompareToVeronica() {
  CompareTo("Run16dAuEmcalDeadMap.txt");
}
void CompareVeronicaJaehyeon() {
  TH2F *hSecMap[8];
  TH2F *hSecMapJ[8];
  for(int i=0; i!=8; ++i) {
    hSecMap[i] = new TH2F( Form("SECMAP%d",i),  ";Y;Z", 100, -0.5, 99.5, 50, -0.5, 49.5 );
    hSecMapJ[i]= new TH2F( Form("SECMAPJ%d",i), ";Y;Z", 100, -0.5, 99.5, 50, -0.5, 49.5 );
  }

  badt.open("Run16dAuEmcalDeadMap.txt");
  for(;;) {
    badt >> sec;
    if(!badt.good()) break;
    badt >> iz >> iy >> bd;
    hSecMap[sec]->Fill(iy,iz,bd+2);
  }
  badt.close();
  badt.open("Run16dAu200WarnMap_ERT.list");
  for(;;) {
    badt >> sec;
    if(!badt.good()) break;
    badt >> iz >> iy >> bd;
    hSecMapJ[sec]->Fill(iy,iz,bd+2);
  }
  badt.close();

  TCanvas *main = new TCanvas();
  main->Divide(4,2,0,0);
  for(int i=0; i!=8; ++i) {
    main->cd(i+1);
    hSecMapJ[i]->Draw("COLZ");
    hSecMap[i]->Draw("TEXT SAME");
  }
}
void CompareTo(TString file) {
  LoadBadTowers();
  TH2F *hSecMap[8];
  TH2F *hSecMapJ[8];
  for(int i=0; i!=8; ++i) {
    hSecMap[i] = new TH2F( Form("SECMAP%d",i),  ";Y;Z", 100, -0.5, 99.5, 50, -0.5, 49.5 );
    hSecMapJ[i]= new TH2F( Form("SECMAPJ%d",i), ";Y;Z", 100, -0.5, 99.5, 50, -0.5, 49.5 );
  }

  for(int i=0; i!=24768; ++i) {
    hSecMap[ twr2iz[i][2] ]->Fill( twr2iz[i][0], twr2iz[i][1], isbad[i] );
  }

  badt.open(file.Data());
  int nnn=0;
  for(;;++nnn) {
    badt >> sec;
    if(!badt.good()) break;
    badt >> iz >> iy >> bd;

    if(sec==6) sec=7;
    else if(sec==7) sec=6;

    if(sec==4) sec=5;
    else if(sec==5) sec=4;

    hSecMapJ[sec]->Fill(iy,iz,bd+2);
  }
  badt.close();

  TCanvas *main = new TCanvas();
  main->Divide(4,2,0,0);
  for(int i=0; i!=8; ++i) {
    main->cd(i+1);
    hSecMap[i]->Draw("COLZ");
    hSecMapJ[i]->Draw("TEXT SAME");
  }
}

void ShowHitDensity() {
  TFile *map = new TFile("hitdensity.root");
  for(int i=0; i!=8; ++i) {
    hHitMap[i] = (TProfile2D*) map->Get( Form("HitMap_SEC%d",i) );
    hHitMap[i]->GetZaxis()->SetRangeUser(0,0.5);
    if(i>5) hHitMap[i]->GetZaxis()->SetRangeUser(0,ddd[i]);
  }

  TCanvas *main = new TCanvas();
  main->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main->cd(i+1);
    hHitMap[i]->Draw("COLZ");
  }
}

int tagTowersPerRun() {
  Load();
  gStyle->SetOptStat(0);
  //CreateHitDensity();
  //CreateBadList();
  //CompareToJaehyeon();
  //CompareToVeronica();
  //CompareVeronicaJaehyeon();
  ShowHitDensity();
}
