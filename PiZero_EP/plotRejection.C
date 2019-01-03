int plotRejection(int run=455053) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TH2F *hist2DS = (TH2F*) file->Get("hPileUpRejectionS");
  TH2F *hist2DN = (TH2F*) file->Get("hPileUpRejectionN");
  TAxis *axis = hist2DS->GetYaxis();

  const int nbinsS = 3;
  double blimitsS[4] = {0,12,24,200};
  const int nbinsN = 3;
  double blimitsN[4] = {0,12,24,200};
  TH1D *hist[20];
  ofstream fout( Form("dat/run%d.dat",run) );
  fout << run << " ";
  for(int i=0; i!=nbinsS; ++i) {
    hist[i] = hist2DS->ProjectionX( Form("histS_%02d%02d",blimitsS[i],blimitsS[i+1]),
				    axis->FindBin(blimitsS[i]+0.1),
				    axis->FindBin(blimitsS[i+1]-0.1));
    fout << hist[i]->GetBinContent(2) << " ";
    fout << hist[i]->GetBinContent(2)/hist[i]->GetBinContent(1) << " ";
  }
  for(int i=0; i!=nbinsN; ++i) {
    hist[i] = hist2DN->ProjectionX( Form("histS_%02d%02d",blimitsN[i],blimitsN[i+1]),
				    axis->FindBin(blimitsN[i]+0.1),
				    axis->FindBin(blimitsN[i+1]-0.1));
    fout << hist[i]->GetBinContent(2) << " ";
    fout << hist[i]->GetBinContent(2)/hist[i]->GetBinContent(1) << " ";
  }
  fout << endl;

  return 0;
}
