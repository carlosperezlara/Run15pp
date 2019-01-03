TF1 *CreateFit(TString name) {
  TF1 *ret = new TF1(name.Data(),"[0]*TMath::Landau(x,[1],[2])",0,0.5);
  ret->SetParameter(1,0.06); ret->SetParLimits(1,0.01,0.1);
  ret->SetParameter(2,0.01); ret->SetParLimits(2,0.001,0.1);
  return ret;
}

int ComputeTimeConstants(int run=428717) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  file->ls();
  TH1F *hM0 = (TH1F*) file->Get("hBBCTimeMean0");
  TH1F *hS0 = (TH1F*) file->Get("hBBCTimeRMS_S0");
  TH1F *hN0 = (TH1F*) file->Get("hBBCTimeRMS_N0");
  TF1 *fS0 = CreateFit("fS0");
  TF1 *fN0 = CreateFit("fN0");
  TF1 *fS0cdf = new TF1("fS0cdf","[0]*ROOT::Math::landau_cdf(x,[2],[1])");
  TF1 *fN0cdf = new TF1("fN0cdf","[0]*ROOT::Math::landau_cdf(x,[2],[1])");
  TF1 *fM0 = new TF1("fM0",
		     "[0]*TMath::Gaus(x,[1],[2])",
		     -3,+3);
  fM0->SetParameter(0,1.0); fM0->SetParLimits(0,0,100000);
  fM0->SetParameter(1,0.0); fM0->SetParLimits(1,-2.,+2.);
  fM0->SetParameter(2,0.2); fM0->SetParLimits(2,0.1,1.0);
  new TCanvas();
  hM0->Draw();
  hM0->Fit("fM0","R","",-3,+3);

  new TCanvas();
  hS0->Draw();
  hS0->Fit("fS0","R","",0.01,0.11);
  fS0->Draw("same");
  fS0cdf->SetParameter(0,1000);//fS0->GetParameter(0));
  fS0cdf->SetParameter(1,fS0->GetParameter(1));
  fS0cdf->SetParameter(2,fS0->GetParameter(2));
  fS0cdf->Draw("same");

  new TCanvas();
  hN0->Draw();
  hN0->Fit("fN0","R","",0.01,0.11);
  fN0->Draw("same");
  fN0cdf->SetParameter(0,1000);//fN0->GetParameter(0));
  fN0cdf->SetParameter(1,fN0->GetParameter(1));
  fN0cdf->SetParameter(2,fN0->GetParameter(2));
  fN0cdf->Draw("same");

  ofstream fout( Form("dat/TimeConstants%d.dat",run) );
  fout << run << " ";
  fout << Form("%0.2f",fM0->GetParameter(1)) << " ";
  fout << Form("%0.3f",fM0->GetParameter(2)) << " ";
  fout << Form("%0.4f",fS0->GetParameter(1)) << " ";
  fout << Form("%0.4f",fS0->GetParameter(2)) << " ";
  fout << Form("%0.3f",fS0cdf->GetX(1000*0.95)) << " ";
  fout << Form("%0.4f",fN0->GetParameter(1)) << " ";
  fout << Form("%0.4f",fN0->GetParameter(2)) << " ";
  fout << Form("%0.3f",fN0cdf->GetX(1000*0.95)) << " " << endl;
  fout.close();

  return 0;
}
