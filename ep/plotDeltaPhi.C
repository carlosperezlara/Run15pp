int plotDeltaPhi(TString run) {
  gStyle->SetOptStat(0);

  TString ooo = Form("out/%s.root",run.Data());
  std::cout << "READING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data() ); 

  TH2F *hDP[2];
  TH2F *hDPPE[2];
  hDP[0] = (TH2F*) f2->Get("DPhiE");
  hDP[1] = (TH2F*) f2->Get("DPhiW");
  hDPPE[0] = (TH2F*) f2->Get("DPhiER");
  hDPPE[1] = (TH2F*) f2->Get("DPhiWR");

  TH1D *h1[13], *h2[13], *h3[13];
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.2);

  //======= HISTOS
  TF1 *fit[2][13];
  for(int i=0; i!=2; ++i) {
    for(int j=0; j!=13; ++j) {
      fit[i][j] = new TF1( Form("F_%d_%d",i,j), "1+2*[0]*TMath::Cos(2*x)", -2.5, +2.5 );
      fit[i][j]->SetLineColor(kBlack);
    }
  }
  TString texttt[13] = { "#varphi - #Psi_{2}^{BBC_{A}}",
			 "#varphi - #Psi_{2}^{BBC_{B}}",
			 "#varphi - #Psi_{2}^{BBC}",
			 "#varphi - #Psi_{2}^{FVTX}",
			 "#varphi - #Psi_{2}^{EX0}",
			 "#varphi - #Psi_{2}^{EX1}",
			 "#varphi - #Psi_{2}^{EX2}",
			 "#varphi - #Psi_{2}^{EX3}",
			 "#varphi - #Psi_{2}^{EX4}",
			 "#varphi - #Psi_{2}^{EX5}",
			 "#varphi - #Psi_{2}^{EX6}",
			 "#varphi - #Psi_{2}^{EX7}",
			 "#varphi - #Psi_{2}^{EX}" };
  TCanvas *cqx = new TCanvas();
  cqx->SetLeftMargin(0.15);
  cqx->SetBottomMargin(0.2);
  cqx->Divide( 4, 6, 0, 0 );
  for(int i=0; i!=13; ++i) {
    cqx->cd(i+1);
    if(i>3) cqx->cd(i);
    h1[i] = (TH1D*) hDP[0]->ProjectionY( Form("dp0_%d",i), i+1, i+1 );
    h2[i] = (TH1D*) hDPPE[0]->ProjectionY( Form("dppe0_%d",i), i+1, i+1 );
    h1[i]->SetLineColor( kRed-3 );
    h1[i]->SetLineWidth(2);
    h2[i]->SetLineColor( kBlue-3 );
    h2[i]->SetLineWidth(2);
    h2[i]->GetYaxis()->SetRangeUser(50.4e3,56.6e3);
    h2[i]->GetXaxis()->SetLabelSize(0.12);
    h2[i]->GetYaxis()->SetLabelSize(0.12);
    h2[i]->GetYaxis()->SetNdivisions(505);
    h2[i]->DrawCopy();
    h1[i]->DrawCopy("same");
    tex->DrawLatex( -5.5, 54.8e3, texttt[i].Data() );
  }
  for(int i=0; i!=13; ++i) {
    cqx->cd(i+13);
    if(i>3) cqx->cd(i+12);
    int b1 = h1[i]->GetXaxis()->FindBin(0.0);
    int b2 = h1[i]->GetXaxis()->FindBin(7.0);
    float bgr = h2[i]->Integral(b1,b2);
    float fgr = h1[i]->Integral(b1,b2);
    h2[i]->Scale( fgr/bgr );
    h1[i]->Divide( h2[i] );
    h1[i]->SetLineWidth(2);
    h1[i]->SetLineColor(kOrange-3);
    h1[i]->GetXaxis()->SetLabelSize(0.12);
    h1[i]->GetYaxis()->SetLabelSize(0.12);
    h1[i]->GetYaxis()->SetRangeUser(0.95,1.09);
    h1[i]->GetYaxis()->SetNdivisions(510);
    h1[i]->DrawCopy();
    h1[i]->Fit( fit[0][i], "", "", -2.5, +2.5 );
    tex->DrawLatex( -5.5, 1.05, texttt[i].Data() );
  }

  TCanvas *cqy = new TCanvas();
  cqy->SetLeftMargin(0.15);
  cqy->SetBottomMargin(0.2);
  cqy->Divide( 4, 6, 0, 0 );
  for(int i=0; i!=13; ++i) {
    cqy->cd(i+1);
    if(i>3) cqy->cd(i);
    h1[i] = (TH1D*) hDP[1]->ProjectionY( Form("dp1_%d",i), i+1, i+1 );
    h2[i] = (TH1D*) hDPPE[1]->ProjectionY( Form("dppe1_%d",i), i+1, i+1 );
    h1[i]->SetLineColor( kRed-3 );
    h1[i]->SetLineWidth(2);
    h2[i]->SetLineColor( kBlue-3 );
    h2[i]->SetLineWidth(2);
    h2[i]->GetYaxis()->SetRangeUser(0,70000);
    h2[i]->GetXaxis()->SetLabelSize(0.12);
    h2[i]->GetYaxis()->SetLabelSize(0.12);
    h2[i]->GetYaxis()->SetNdivisions(505);
    h2[i]->DrawCopy();
    h1[i]->DrawCopy("same");
    tex->DrawLatex( -5.5, 50e3, texttt[i].Data() );
  }
  for(int i=0; i!=13; ++i) {
    cqy->cd(i+13);
    if(i>3) cqy->cd(i+12);
    int b1 = h1[i]->GetXaxis()->FindBin(0.0);
    int b2 = h1[i]->GetXaxis()->FindBin(7.0);
    float bgr = h2[i]->Integral(b1,b2);
    float fgr = h1[i]->Integral(b1,b2);
    h2[i]->Scale( fgr/bgr );
    h1[i]->Divide( h2[i] );
    h1[i]->SetLineWidth(2);
    h1[i]->SetLineColor(kOrange-3);
    h1[i]->GetXaxis()->SetLabelSize(0.12);
    h1[i]->GetYaxis()->SetLabelSize(0.12);
    h1[i]->GetYaxis()->SetRangeUser(0.95,1.09);
    h1[i]->GetYaxis()->SetNdivisions(510);
    h1[i]->DrawCopy();
    h1[i]->Fit( fit[1][i], "", "", -4.5, +4.5 );
    fit[1][i]->Draw("SAME");
    tex->DrawLatex( -5.5, 1.05, texttt[i].Data() );
  }

  return 0;
}
