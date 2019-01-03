int plotCentrality() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile("out/all.root");
  TH2F *hist2D = (TH2F*) file->Get("hCentralitySelection");
  TAxis *axis = hist2D->GetXaxis();

  const int nbins = 8;
  double blimits[9] = {0,2,3,5,10,20,40,60,80};
  int color[8] = { 
		   kMagenta-3,
		   kRed-3,
		   kOrange-3,
		   kGray,
		   kBlack,
		   kGreen-3,
		   kCyan-3,
		   kBlue-3};
  TCanvas *can = new TCanvas();
  can->SetLogy(1);
  TH2F *axis2D = new TH2F("axis2D","d+Au  200 GeV;BBC  south  signal  (a.u.);counts  (a.u.)",
			  100,0,200,100,0.9e3,1.1e7);
  axis2D->Draw();
  TH1D *hist[nbins];
  for(int i=0; i!=nbins; ++i) {
    hist[i] = hist2D->ProjectionY( Form("hist%02d%02d",blimits[i],blimits[i+1]),
				   axis->FindBin(blimits[i]+0.1),
				   axis->FindBin(blimits[i+1]-0.1));
    hist[i]->SetLineColor( color[i] );
    hist[i]->SetLineWidth(2);
    hist[i]->Draw("same");
  }
  TLegend *leg = new TLegend(0.7,0.49,0.88,0.88,"Centrality Class");
  for(int i=0; i!=8; ++i) {
    leg->AddEntry(hist[i], Form("[%.0f,%.0f]",blimits[i],blimits[i+1]) );
  }
  leg->Draw();
  can->SaveAs("dAu_entrality.pdf","pdf");
  return 0;
}
