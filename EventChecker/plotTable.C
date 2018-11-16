int plotTable() {
  ifstream fin("table.dat");
  TString labx;
  int n = 0;
  Float_t x[1000];
  Float_t ex[1000];
  Float_t bm0[1000];
  Float_t bm1[1000];
  Float_t ebm0[1000];
  Float_t ebm1[1000];
  Float_t tm0[1000];
  Float_t tm1[1000];
  Float_t etm0[1000];
  Float_t etm1[1000];
  for(;;) {
    fin >> labx;
    if(!fin.good()) break;
    fin >> bm0[n] >> ebm0[n] >> bm1[n] >> ebm1[n];
    fin >> tm0[n] >> etm0[n] >> tm1[n] >> etm1[n];
    x[n] = n;
    ex[n] = 0;
    n++;
  }
  TGraphErrors *gBM0 = new TGraphErrors(n,x,bm0,ex,ebm0);
  TGraphErrors *gBM1 = new TGraphErrors(n,x,bm1,ex,ebm1);
  TGraphErrors *gTM0 = new TGraphErrors(n,x,tm0,ex,etm0);
  TGraphErrors *gTM1 = new TGraphErrors(n,x,tm1,ex,etm1);

  gBM0->SetLineColor( kRed-3 );  gBM0->SetMarkerColor( kRed-3 );  gBM0->SetMarkerStyle( 20 );
  gBM1->SetLineColor( kBlue-3 ); gBM1->SetMarkerColor( kBlue-3 ); gBM1->SetMarkerStyle( 20 );
  gTM0->SetLineColor( kRed-3 );  gTM0->SetMarkerColor( kRed-3 );  gTM0->SetMarkerStyle( 20 );
  gTM1->SetLineColor( kBlue-3 ); gTM1->SetMarkerColor( kBlue-3 ); gTM1->SetMarkerStyle( 20 );

  TH2F *axis1 = new TH2F("axis1",";segment  index;BBC  multiplicity",  100,-0.5,999.5,100,0.,16.);
  TH2F *axis2 = new TH2F("axis2",";segment  index;Track  multiplicity",100,-0.5,999.5,100,0.,3.);

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

}
