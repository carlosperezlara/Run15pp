int nptbins[5];
int nptbinsERT[5];

float xxx[5][100];
float exx[5][100];
float yy0[5][100];
float ey0[5][100];
float yy1[5][100];
float ey1[5][100];
float yy2[5][100];
float ey2[5][100];

float xxxERT[5][100];
float exxERT[5][100];
float yy0ERT[5][100];
float ey0ERT[5][100];
float yy1ERT[5][100];
float ey1ERT[5][100];
float yy2ERT[5][100];
float ey2ERT[5][100];

void style(TGraphErrors *gr, int color) {
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetFillColor(kWhite);
}

void loadfiles() {
  int ord, nbins;
  float xmin, xmax;
  ifstream fin("vn.dat");
  for(;;) {
    fin >> ord >> nbins;
    if(!fin.good()) break;
    nptbins[ord] = nbins;
    for(int i=0; i!=nbins; ++i) {
      fin >> xmin >> xmax;
      xxx[ord][i] = 0.5*(xmax+xmin);
      exx[ord][i] = 0.5*(xmax-xmin);
      fin >> yy0[ord][i] >> ey0[ord][i];
      fin >> yy1[ord][i] >> ey1[ord][i];
      fin >> yy2[ord][i] >> ey2[ord][i];
    }
  }
  fin.close();
  fin.open("vnERT.dat");
  for(;;) {
    fin >> ord >> nbins;
    if(!fin.good()) break;
    nptbinsERT[ord] = nbins;
    for(int i=0; i!=nbins; ++i) {
      fin >> xmin >> xmax;
      xxxERT[ord][i] = 0.5*(xmax+xmin);
      exxERT[ord][i] = 0.5*(xmax-xmin);
      fin >> yy0ERT[ord][i] >> ey0ERT[ord][i];
      fin >> yy1ERT[ord][i] >> ey1ERT[ord][i];
      fin >> yy2ERT[ord][i] >> ey2ERT[ord][i];
    }
  }
}

int plot() {
  loadfiles();

  TGraphErrors *gr0 = new TGraphErrors(nptbins[0],xxx[0],yy1[0],exx[0],ey1[0]);
  TGraphErrors *gr1 = new TGraphErrors(nptbins[1],xxx[1],yy1[1],exx[1],ey1[1]);
  TGraphErrors *gr2 = new TGraphErrors(nptbins[2],xxx[2],yy1[2],exx[2],ey1[2]);
  TGraphErrors *gr3 = new TGraphErrors(nptbins[3],xxx[3],yy1[3],exx[3],ey1[3]);
  TGraphErrors *gr4 = new TGraphErrors(nptbins[4],xxx[4],yy1[4],exx[4],ey1[4]);
  style(gr0,kBlue-3);
  style(gr1,kBlue-3);
  style(gr2,kBlue-3);
  style(gr3,kBlue-3);
  style(gr4,kBlue-3);
  TGraphErrors *gr0ERT = new TGraphErrors(nptbinsERT[0],xxxERT[0],yy1ERT[0],exxERT[0],ey1ERT[0]);
  TGraphErrors *gr1ERT = new TGraphErrors(nptbinsERT[1],xxxERT[1],yy1ERT[1],exxERT[1],ey1ERT[1]);
  TGraphErrors *gr2ERT = new TGraphErrors(nptbinsERT[2],xxxERT[2],yy1ERT[2],exxERT[2],ey1ERT[2]);
  TGraphErrors *gr3ERT = new TGraphErrors(nptbinsERT[3],xxxERT[3],yy1ERT[3],exxERT[3],ey1ERT[3]);
  TGraphErrors *gr4ERT = new TGraphErrors(nptbinsERT[4],xxxERT[4],yy1ERT[4],exxERT[4],ey1ERT[4]);
  style(gr0ERT,kRed-3);
  style(gr1ERT,kRed-3);
  style(gr2ERT,kRed-3);
  style(gr3ERT,kRed-3);
  style(gr4ERT,kRed-3);


  TH2F *axis0 = new TH2F("axis0","",100,0,20,100,-0.05,+0.03);
  TH2F *axis1 = new TH2F("axis1","",100,0,20,100,-0.01,+0.04);
  TH2F *axis2 = new TH2F("axis2","",100,0,20,100,-0.02,+0.05);
  TH2F *axis3 = new TH2F("axis3","",100,0,20,100,-0.03,+0.03);
  TH2F *axis4 = new TH2F("axis4","",100,0,20,100,-0.03,+0.03);

  TCanvas *main = new TCanvas("main","main");
  main->Divide(3,2);
  main->cd(1); axis0->Draw(); gr0->Draw("PSAME"); gr0ERT->Draw("PSAME");
  main->cd(2); axis1->Draw(); gr1->Draw("PSAME"); gr1ERT->Draw("PSAME");
  main->cd(3); axis2->Draw(); gr2->Draw("PSAME"); gr2ERT->Draw("PSAME");
  main->cd(4); axis3->Draw(); gr3->Draw("PSAME"); gr3ERT->Draw("PSAME");
  main->cd(5); axis4->Draw(); gr4->Draw("PSAME"); gr4ERT->Draw("PSAME");
  return 0;
}
