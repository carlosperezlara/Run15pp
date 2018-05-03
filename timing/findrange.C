#include "configtime.h"

int findrange(int thistwr=-1) {
  LoadBadTowers();
  bool drawonly = false;
  if(thistwr>=0) drawonly = true;
  TFile *file = new TFile("out/all2.root");
  ofstream fout;
  if(!drawonly) fout.open(Form("ranges.dat"));
  for(int i=0; i!=8; ++i) {
    TProfile2D *prof = (TProfile2D*) file->Get( Form("MeanTDC_SEC%d",i) );
    for(int j=0; j!=prof->GetYaxis()->GetNbins(); ++j) {
      int twr = secs[i]+j;
      if(drawonly)
	if(twr!=thistwr) continue;
      if(isbad[twr]>0) {
	fout << twr << " " << 0.0 << " " << 4000.0 << endl;
      continue;
      }
      TH1D *pro = (TH1D*) prof->ProjectionX( Form("twr_%d",j), j+1, j+1 );
      TF1 *fit = new TF1("fit","[0]+[1]/x");
      pro->Fit(fit,"QR","",1000,4000);
      int p0 = fit->GetParameter(0);
      fout << twr << " " << p0-400 << " " << p0+100 << endl;
      if(drawonly) {
	pro->Draw();
	break;
      } else {
	delete fit;
	delete pro;
      }
    }
  }
  if(!drawonly) fout.close();
  cout << "ALL DONE." << endl;

  return 0;
}
