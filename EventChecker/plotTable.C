int nrun=0;
double rate[1000];
int run[1000];

void loadrates() {
  ifstream fin("../runs.info.dat");
  for(;;++nrun) {
    fin >> run[nrun];
    if(!fin.good()) break;
    fin >> rate[nrun] >> rate[nrun] >> rate[nrun];
  }
}

int plotTable() {
  loadrates();
  ifstream fin("table.dat");
  TString labx;
  int n = 0;
  Float_t x[10000];
  Float_t xr[10000];
  Float_t ex[10000];
  Float_t bm0[10000];
  Float_t bm1[10000];
  Float_t ebm0[10000];
  Float_t ebm1[10000];
  Float_t tm0[10000];
  Float_t tm1[10000];
  Float_t etm0[10000];
  Float_t etm1[10000];
  for(;;) {
    fin >> labx;
    if(!fin.good()) break;
    fin >> bm0[n] >> ebm0[n] >> bm1[n] >> ebm1[n];
    fin >> tm0[n] >> etm0[n] >> tm1[n] >> etm1[n];
    TObjArray *oarr = labx.Tokenize("_");
    TString srun = ((TObjString*) oarr->At(0))->GetString();
    int trun = srun.Atoi();
    cout << trun << " " << nrun << " ";

    for(int ii=0;ii!=nrun;++ii) {
      if(run[ii]==trun) {
	xr[n] = rate[ii];
	cout << " F ";
	break;
      }
    }
    x[n] = n;
    ex[n] = 0;
    n++;
    cout << labx << endl;
  }
  cout << n << endl;
  TGraphErrors *gBM0 = new TGraphErrors(n,x,bm0,ex,ebm0);
  TGraphErrors *gBM1 = new TGraphErrors(n,x,bm1,ex,ebm1);
  TGraphErrors *gTM0 = new TGraphErrors(n,x,tm0,ex,etm0);
  TGraphErrors *gTM1 = new TGraphErrors(n,x,tm1,ex,etm1);

  gBM0->SetLineColor( kRed-3 );  gBM0->SetMarkerColor( kRed-3 );  gBM0->SetMarkerStyle( 20 );
  gBM1->SetLineColor( kBlue-3 ); gBM1->SetMarkerColor( kBlue-3 ); gBM1->SetMarkerStyle( 20 );
  gTM0->SetLineColor( kRed-3 );  gTM0->SetMarkerColor( kRed-3 );  gTM0->SetMarkerStyle( 20 );
  gTM1->SetLineColor( kBlue-3 ); gTM1->SetMarkerColor( kBlue-3 ); gTM1->SetMarkerStyle( 20 );

  TH2F *axis1 = new TH2F("axis1",";segment  index;BBC  multiplicity",
			 100,-0.5,3999.5,100,0.,126.);
  TH2F *axis2 = new TH2F("axis2",";segment  index;Track  multiplicity",
			 100,-0.5,3999.5,100,0.,26.);

  TCanvas *cV = new TCanvas();
  cV->Divide(2,1);
  cV->cd(1);
  axis1->Draw();
  gBM0->Draw("PSAME");
  gBM1->Draw("PSAME");
  cV->cd(2);
  axis2->Draw();
  gTM0->Draw("PSAME");
  gTM1->Draw("PSAME");

  TCanvas *cv2 = new TCanvas();
  TGraphErrors *gR0 = new TGraphErrors(n,xr,bm0,ex,ebm0);
  TGraphErrors *gR1 = new TGraphErrors(n,xr,bm1,ex,ebm1);
  gR0->SetLineColor( kRed-3 );
  gR0->SetMarkerColor( kRed-3 );
  gR0->SetMarkerStyle( 20 );
  gR1->SetLineColor( kBlue-3 );
  gR1->SetMarkerColor( kBlue-3 );
  gR1->SetMarkerStyle( 20 );
  
  TH2F *axisR = new TH2F("axisR",";rate  [Hz];BBC  multiplicity",
			 100,0.,2.5e6,100,8.,20.);
  axisR->Draw();
  gR0->Draw("PSAME");
  gR1->Draw("PSAME");
}
