#include "configtime.h"

int checkmean(int block=-1, bool justone = true) {
  TFile *file = new TFile("out/allTOF5.root");

  float x[10000];
  float y[10000];
  float ex[10000];
  float ey[10000];
  LoadBadTowers();
  ofstream fout;
  if(!justone) fout.open( Form("check_B%d.dat",block) );
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
	fout << j << " " << 0 << endl;
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
    maxb = 0;
    xmin = maxb-4.5;
    xmax = maxb+4.5;
    proj->GetXaxis()->SetRangeUser(xmin,xmax);
    cout << "TWR " << j << " ";
    float maxat =  proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() );
    if(justone) {
      //proj->Draw();
    } else {
      fout << j << " " << maxat << endl;
    }
    cout << " MAX AT " << maxat;
    cout << endl;
  }
  if(!justone)
    fout.close();
  return 0;
}
