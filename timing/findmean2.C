#include "configtime.h"

int findmean2(int block=-1, bool justone = true, bool hacked=false, float inix=-100) {
  TFile *file = new TFile("out/allTOF4.root");
  TCanvas *main = new TCanvas();
  float x[10000];
  float y[10000];
  float ex[10000];
  float ey[10000];
  LoadBadTowers();
  ofstream fout;
  if(!justone) fout.open( Form("offset_B%05d.dat",block) );
  int ini = 0;
  int fin = 24768;
  if(justone) {
    ini = block;
    fin = block+1;
  } else {
    ini = block*25;
    fin = (block+1)*25;
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
    maxb = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() ); //proj->GetMean();
    if(inix>-99) maxb = inix;
    float peak = maxb;
    xmin = maxb-2.0;
    xmax = maxb+2.0;
    int rangeswid[4] = {2.0,2.0,2.0,2.0};
    float scale=1.5;
    if(j<secs[6]) scale=1.0;
    cout << "TWR " << j << " " << endl;
    for(int ittr=0; ittr!=4; ++ittr) {
      cout << " => " << xmin << " " << xmax << endl;
      proj->GetXaxis()->SetRangeUser(xmin,xmax);
      maxb = proj->GetMean();
      if(!hacked) {
	maxb = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() ); //proj->GetMean();
	xmin = maxb-rangeswid[ittr]*scale;
	xmax = maxb+rangeswid[ittr]*scale*0.8;
      }
    }
    proj->GetXaxis()->SetRangeUser(-10,+10);
    TF1 *fit = new TF1(Form("FIT_%d",j),"[0]*TMath::Gaus(x,[1],[2])");
    float control;
    TString proper = "RMEQ";
    if(justone) {
      proj->Draw();
      proper = "RME";
    }
    for(int tryout=0; ; ++tryout) {
      fit->SetParameter(2,peak);
      if(hacked) {
	fit->SetParLimits(2,0.4,6.0);
      } else {
	fit->SetParLimits(2,0.4,3.0);
      }
      cout << " => " << xmin << " " << xmax << endl;
      proj->Fit(fit,proper.Data(),"",xmin,xmax);
      proj->Fit(fit,proper.Data(),"",xmin,xmax);
      control = fit->GetChisquare()/fit->GetNDF();
      float maxcontrol = 4.5;
      if(hacked) maxcontrol = 3.0;
      if(j<secs[6]) {
	if(control<maxcontrol) break;
      } else {
	if(control<maxcontrol) break;
      }
      cout << "COULD NOT FIND ACCEPTABLE CHI2/NDF = " << control;
      cout << " TRYING AGAIN. TIME " << tryout+1 << ". SHRINKING RANGE BY 10PERCENT";
      if(tryout>1) cout << " t" << j;
      cout << endl;
      xmin += 0.06*(xmax-xmin);
      xmax -= 0.10*(xmax-xmin);
    }
    if(justone) {
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
