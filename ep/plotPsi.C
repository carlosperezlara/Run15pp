int plotPsi(TString run) {
  gStyle->SetOptStat(0);

  TString ooo = Form("out/%s.root",run.Data());
  std::cout << "READING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data() ); 

  //======= HISTOS
  TH2F *hPSI0 = (TH2F*) f2->Get("PSI0");
  TH2F *hPSI1 = (TH2F*) f2->Get("PSI1");
  TH2F *hPSI2 = (TH2F*) f2->Get("PSI2");
  TH1D *h1[14], *h2[14], *h3[14];
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.2);
  TString tit14[14] = {"BBA","BBB","BB", "FV", "EX0","EX1","EX2",
		       "EX3","EX4","EX5","EX6","EX7","EX8","CA"};
  TCanvas *cpsi = new TCanvas();
  cpsi->Divide( 4, 4 );
  for(int i=0; i!=14-1; ++i) {
    h1[i] = hPSI0->ProjectionY( Form("PSI0_DET_%d",i), i+1, i+1 );
    h2[i] = hPSI1->ProjectionY( Form("PSI1_DET_%d",i), i+1, i+1 );
    h3[i] = hPSI2->ProjectionY( Form("PSI2_DET_%d",i), i+1, i+1 );
    h1[i]->SetLineColor(kBlack);
    h2[i]->SetLineColor(kBlue-3);
    h3[i]->SetLineColor(kRed-3);
    h1[i]->SetLineWidth(2);
    h2[i]->SetLineWidth(2);
    h3[i]->SetLineWidth(2);
    cpsi->cd(i+1);
    h1[i]->GetYaxis()->SetRangeUser(0,h1[i]->GetMaximum());
    h1[i]->GetXaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetLabelSize(0.08);
    h1[i]->DrawCopy();
    h2[i]->DrawCopy("same");
    h3[i]->DrawCopy("same");
    tex->DrawLatex( 0.2, 10e4, tit14[i].Data() );
  }
  cpsi->cd(14);
  //h1[0]->Reset();
  //h1[0]->Draw();
  TLegend *leg = new TLegend(0.1,0.1,0.9,0.9,"#Psi_{2} distribution");
  leg->SetFillColor(kWhite);
  leg->AddEntry(h1[0],"Raw");
  leg->AddEntry(h2[0],"ReCentered");
  leg->AddEntry(h3[0],"Flattened");
  leg->Draw();

  return 0;
}
