#include "configtime.h"

int correctfindmean2(int block=-1, bool justone = true) {
  TFile *file = new TFile("out/allTOF4.root");

  float x[10000];
  float y[10000];
  float ex[10000];
  float ey[10000];
  LoadBadTowers();
  ofstream fout;
  if(!justone) fout.open( Form("offset_B%d.dat",block) );
  int ini = 0;
  int fin = 24768;
  if(justone) {
    ini = block;
    fin = block+1;
  } else {
    ini = block*2500;
    fin = (block+1)*2500;
  }
  if(fin>24768) fin = 24768; 
  int nskip = 1;
  TH2F *pro = (TH2F*) file->Get("hTOF");
  int bmax = pro->GetXaxis()->GetNbins();
  int bymax = pro->GetYaxis()->GetNbins();
  int np = 0;
  TH1F *proj = new TH1F("temp","",
			bymax,
			pro->GetYaxis()->GetBinLowEdge(1),
			pro->GetYaxis()->GetBinLowEdge(bymax+1));
  for(int j=ini;j<fin;j++) {
    if(isbad[j]) {
      if(!justone) {
	fout << j << " " << 0 << " " << 0 << " " << 0 << endl;
	cout << "TRW " << j << " BAD TOWER. SKIPPING... " << endl;
      }
      continue;
    }
    proj->Reset();
    for(int k=0; k!=bymax; ++k) {
      if(pro->GetBinContent(j+1,k+1)>0)
	proj->Fill( pro->GetYaxis()->GetBinCenter(k+1),
		    pro->GetBinContent(j+1,k+1) );
    }
    int nn = proj->GetEntries();
    if(nn<30) continue;
    float maxb;
    float xmin;
    float xmax;
    float scale=1.5;
    if(j<secs[6]) scale=1.0;
    maxb = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() );
    xmin = maxb-1.5;
    xmax = maxb+1.5;
    int rangeswid[4] = {3.5,3.5,3.5,3.0};
    for(int ittr=0; ittr!=0; ++ittr) {
      proj->GetXaxis()->SetRangeUser(xmin,xmax);
      maxb = proj->GetMean();
      xmin = maxb-rangeswid[ittr];
      xmax = maxb+rangeswid[ittr];
    }
    proj->GetXaxis()->SetRangeUser(-10,+10);
    TF1 *fit = new TF1(Form("FIT_%d",j),"[0]*TMath::Gaus(x,[1],[2])");
    fit->SetParameter(2,1.0);
    fit->SetParLimits(2,0.4,3.0);
    float control;
    TString proper = "RMEQ";
    if(justone) {
      proj->Draw();
      proper = "RME";
    }
    cout << "TWR " << j << " " << xmin << " " << xmax << endl;
    for(int tryout=0; ; ++tryout) {
      //cout << xmin << " " << xmax << endl;
      proj->Fit(fit,proper.Data(),"",xmin,xmax);
      proj->Fit(fit,proper.Data(),"",xmin,xmax);
      control = fit->GetChisquare()/fit->GetNDF();
      if(j<secs[6]) {
	if(control<2.5) break;
      } else {
	if(control<2.5) break;
      }
      cout << "COULD NOT FIND ACCEPTABLE CHI2/NDF = " << control;
      cout << " TRYING AGAIN WITH " << tryout+1 << " LESS POINTS.";
      if(tryout>1) cout << " t" << j;
      cout << endl;
      xmin += 0.1*(xmax-xmin);
      xmax -= 0.1*(xmax-xmin);
    }
    if(justone) {
      //fit->Draw("same");
      cout << j << " " << fit->GetParameter(0) << " " << fit->GetParameter(1);
      cout << " " << fit->GetParameter(2) << endl;
    } else {
      fout << j << " " << fit->GetParameter(0) << " " << fit->GetParameter(1);
      fout << " " << fit->GetParameter(2) << endl;
      delete fit;
    }
    cout << " CHI2/NDF " << control;
    cout << endl;
  }
  if(!justone)
    fout.close();
  return 0;
}
