#include "configtime.h"

int findmean(int block=-1, bool justone = true) {
  TFile *file = new TFile("out/allTOF2.root");

  float x[10000];
  float y[10000];
  float ex[10000];
  float ey[10000];
  LoadBadTowers();
  ofstream fout;
  if(!justone) fout.open( Form("walks_B%d.dat",block) );
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
  for(int j=ini;j<fin;j++) {
    if(isbad[j]) {
      if(!justone) {
	fout << j << " " << 0 << " " << 0 << " " << 0 << endl;
	cout << "TRW " << j << " BAD TOWER. SKIPPING... " << endl;
      }
      continue;
    }
    TH2F *pro = (TH2F*) file->Get(Form("H2_TDCvsADC_%d",j));
    int bmax = pro->GetXaxis()->GetNbins();
    int bymax = pro->GetYaxis()->GetNbins();
    int np = 0;
    TH1F *proj = new TH1F(Form("temp_%d",j),"",
			  bymax,
			  pro->GetYaxis()->GetBinLowEdge(1),
			  pro->GetYaxis()->GetBinLowEdge(bymax+1));
    for(int b=1; b!=bmax+1; ++b) {
      proj->Reset();
      for(int k=0; k!=bymax+1; ++k) {
	if(pro->GetBinContent(b,k)>0)
	  proj->Fill( pro->GetYaxis()->GetBinCenter(k),
		      pro->GetBinContent(b,k) );
      }
      int nn = proj->GetEntries();
      if(nn<30) continue;
      ++np;
      if(np-1<nskip) {
	continue;
      }
      //cout << "B " << b << " NP " << np << " NSKIP " << nskip << endl;
      //cout << "=== TAKEN " << endl;
      float maxb;
      float xmin;
      float xmax;
      maxb = proj->GetMean();
      xmin = maxb-500;
      xmax = maxb+500;
      int rangeswid[4] = {100,80,50,20};
      for(int ittr=0; ittr!=4; ++ittr) {
	proj->GetXaxis()->SetRangeUser(xmin,xmax);
	maxb = proj->GetMean();
	xmin = maxb-rangeswid[ittr];
	xmax = maxb+rangeswid[ittr];
	//cout << maxb << " " << xmin << " " << xmax << " " << proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() ) << endl;
      }
      //cout << endl;
      int inp = np -nskip-1;
      y[inp] = proj->GetMean();
      //cout << y[inp] << " ";
      //y[inp] = proj->GetMaximum();
      //cout << y[inp] << " ";
      //cout << endl;
      ey[inp] = proj->GetMeanError();
      float maxcenter = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() );
      float widbin = proj->GetXaxis()->GetBinWidth(1);
      if(j<secs[6]) {
	widbin *= 5;
      } else { 
	widbin *= 6;
      }
      if( TMath::Abs(y[inp]-maxcenter)>widbin ) {
	np--;
	continue;
      }

      if( TMath::Abs(maxb-y[inp])>200 ) {
	np--;
	continue;
      }
      if( TMath::Abs(ey[inp])<1e-3 ) {
	np--;
	continue;
      }
      if( TMath::Abs(ey[inp])>10 ) {
	np--;
	continue;
      }
      x[inp] = pro->GetXaxis()->GetBinCenter(b);
      ex[inp] = pro->GetXaxis()->GetBinWidth(1)/2;
      //cout << inp << " " << maxb << " " << xmin << " " << xmax << " || ";
      //cout << nn << " " << y[inp] << " " << ey[inp] << endl;
    }
    cout << "TWR " << j << " NNN===>";
    int nreal = np-nskip;
    cout << nreal-1; // promise that I will remove one point
    if(nreal<10) {
      if(!justone) {
	fout << j << " " << 0 << " " << 0 << " " << 0 << endl;
      }
    } else {
      TGraphErrors *gr = new TGraphErrors(nreal,x,y,ex,ey);
      gr->RemovePoint(0); // always remove the first point
      TF1 *fit = new TF1(Form("FIT_%d",j),"[0]+[1]/x+[2]/x/x",0,5000);
      fit->SetParameter(0,2000); fit->SetParLimits(0,0,4090);
      fit->SetParameter(1,-5e4); fit->SetParLimits(1,-1e8,0);
      fit->SetParameter(2,-0.1); fit->SetParLimits(2,-1e8,+1e8);
      float control;
      TString proper = "RMEQ";
      if(justone) {
	pro->Draw("colz");
	gr->Draw("*Lsame");
	proper = "RME";
      }
      for(int tryout=0; ; ++tryout) {
	gr->Fit(fit,proper.Data(),"",0,3000);
	gr->Fit(fit,proper.Data(),"",0,3000);
	control = fit->GetChisquare()/gr->GetN();
	if(j<secs[6]) {
	  if(control<4.0) break;
	} else {
	  if(control<8.0) break;
	}
	//if(control<1.7) break;
	//if(control<1.0) break;
	cout << "COULD NOT FIND ACCEPTABLE CHI2/NDF = " << control;
	cout << " TRYING AGAIN WITH " << tryout+1 << " LESS POINTS.";
	if(tryout>1) cout << " t" << j;
	cout << endl;
	gr->RemovePoint(0);
	fit->SetParameter(0,2000); fit->SetParLimits(0,0,4090);
	fit->SetParameter(1,-5e4); fit->SetParLimits(1,-1e8,0);
	fit->SetParameter(2,-0.1); fit->SetParLimits(2,-1e8,+1e8);
      }
      if(justone) {
	fit->Draw("same");
      } else {
	fout << j << " " << fit->GetParameter(0) << " " << fit->GetParameter(1);
	fout << " " << fit->GetParameter(2) << endl;
	delete fit;
	delete gr;
      }
      cout << " CHI2/NDF " << control;
    }
    delete proj;
    cout << endl;
  }
  if(!justone)
    fout.close();
  return 0;
}
