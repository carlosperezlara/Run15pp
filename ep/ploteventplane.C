int ploteventplane(TString run) {
  gStyle->SetOptStat(0);

  ifstream fromSylvia("v2_0_0_FVTX1S.dat");
  //ifstream fromSylvia("v2_0_0_FVTX2p4LS.dat");
  int nnn=0;
  float pfsX[30];
  float pfsY[30];
  float pfsXe[30];
  float pfsYe[30];
  for(;;++nnn) {
    fromSylvia >> pfsX[nnn];
    if(!fromSylvia.good()) break;
    fromSylvia >> pfsY[nnn] >> pfsYe[nnn];
    pfsXe[nnn] = 0;
  }
  fromSylvia.close();
  cout << "READ POINTS: " << nnn << endl;
  TGraphErrors *gFS = new TGraphErrors(nnn,pfsX,pfsY,pfsXe,pfsYe);
  gFS->SetLineColor(kGray);
  gFS->SetMarkerColor(kGray);
  gFS->SetMarkerStyle(25);

  TString ooo = Form("out/myResults_%s.root",run.Data());
  std::cout << "READING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data() ); 

  //======= HISTOS
  TH1F *hEvents = (TH1F*) f2->Get("Events");
  TH1F *hFrac = (TH1F*) f2->Get("FRAC");
  TH1F *hCent = (TH1F*) f2->Get("CENT");
  TH1F *hVtxZ = (TH1F*) f2->Get("VTXZ");

  TH1F *hEp   = (TH1F*) f2->Get("Ep");
  TH1F *hChi2 = (TH1F*) f2->Get("Chi2");
  TH1F *hDphi = (TH1F*) f2->Get("pc3Dphi");
  TH1F *hphi  = (TH1F*) f2->Get("phi");
  TH1F *hDz   = (TH1F*) f2->Get("pc3Dz");
  TH1F *hZED  = (TH1F*) f2->Get("ZED");
  TH2F *hBeta = (TH2F*) f2->Get("Beta");

  TH2F *hQx = (TH2F*) f2->Get("Qx");
  TH2F *hQy = (TH2F*) f2->Get("Qy");
  TProfile2D *hPsiC = (TProfile2D*) f2->Get("PsiC");
  TProfile2D *hPsiS = (TProfile2D*) f2->Get("PsiS");

  TH2F *hPSI  = (TH2F*) f2->Get("PSI0");
  TH2F *hPSI2 = (TH2F*) f2->Get("PSI1");
  TH2F *hPSI3 = (TH2F*) f2->Get("PSI2");
  TH2F *hResD = (TH2F*) f2->Get("ResD");
  TProfile *hRes = (TProfile*) f2->Get("Res");

  TH2F *hEPC[3][11];
  for(int i=0; i!=3; ++i) for(int j=0; j!=11; ++j) {
      hEPC[i][j] = (TH2F*) f2->Get( Form("EVC%d_DET%d",i,j) );
      hEPC[i][j]->GetXaxis()->SetLabelSize(0.12);
      hEPC[i][j]->GetYaxis()->SetLabelSize(0.12);
    }
  TH1F *hreso = new TH1F("hreso","",13,-0.5,12.5);
  hreso->SetBinContent( 1, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(1) ) ) ); // BBC 1
  hreso->SetBinContent( 2, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(1) ) ) ); // BBC 1
  hreso->SetBinContent( 3, TMath::Sqrt( 2*TMath::Abs( hRes->GetBinContent(1) ) ) ); // BBC 1*sqrt2
  hreso->SetBinContent( 4, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(2)*hRes->GetBinContent(20)/
						    hRes->GetBinContent(11) ) ) ); // FVTX 2*20/11
  hreso->SetBinContent( 5, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(3)*hRes->GetBinContent(21)/
						    hRes->GetBinContent(7) ) ) ); // EX0/1B 3*21/7
  hreso->SetBinContent( 6, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(4)*hRes->GetBinContent(22)/
						    hRes->GetBinContent(8) ) ) ); // EX1/0B 4*22/8
  hreso->SetBinContent( 7, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(5)*hRes->GetBinContent(23)/
						    hRes->GetBinContent(9) ) ) ); // EX2/3B 5*23/9
  hreso->SetBinContent( 8, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(6)*hRes->GetBinContent(24)/
						    hRes->GetBinContent(10) ) ) ); // EX3/2B 6*24/10
  hreso->SetBinContent( 9, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(7)*hRes->GetBinContent(21)/
						    hRes->GetBinContent(3) ) ) ); // EX4/5B 7*21/3
  hreso->SetBinContent(10, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(8)*hRes->GetBinContent(22)/
						    hRes->GetBinContent(4) ) ) ); // EX5/4B 8*22/4
  hreso->SetBinContent(11, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(9)*hRes->GetBinContent(23)/
						    hRes->GetBinContent(5)) ) ); // EX6/7B 9*23/5
  hreso->SetBinContent(12, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(10)*hRes->GetBinContent(24)/
						    hRes->GetBinContent(6) ) ) ); // EX7/6B 10*24/6
  hreso->SetBinContent(13, TMath::Sqrt( TMath::Abs( hRes->GetBinContent(11)*hRes->GetBinContent(20)/
						    hRes->GetBinContent(2) ) ) ); // EX/FB 11*20/2
  
  TProfile2D *hV2[2];
  hV2[0] = (TProfile2D*) f2->Get("hV2_E");
  hV2[1] = (TProfile2D*) f2->Get("hV2_W");
  TProfile2D *hV2All = (TProfile2D*) hV2[0]->Clone( "hV2All" );
  hV2All->Add( hV2[1] );

  TH2D *h2V2[2], *h2V2All;
  h2V2[0] = hV2[0]->ProjectionXY("EAST");
  h2V2[1] = hV2[1]->ProjectionXY("WEST");
  h2V2All = hV2All->ProjectionXY("V2RAW");

  TH2F *hDP[2];
  TH2F *hDPPE[2];
  hDP[0] = (TH2F*) f2->Get("DPhiE");
  hDP[1] = (TH2F*) f2->Get("DPhiW");
  hDPPE[0] = (TH2F*) f2->Get("DPhiER");
  hDPPE[1] = (TH2F*) f2->Get("DPhiWR");

  TH1D *h1[13], *h2[13], *h3[13];
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.2);

  TString textt[13] = { "BBC SubA", "BBC SubB", "BBC Full", "FVTX",
			"EX 0i", "EX 1i", "EX 2i", "EX 3i",
			"EX 0o", "EX 1o", "EX 2o", "EX 3o",
			"EX Full"};
  TCanvas *cpsi = new TCanvas();
  cpsi->Divide( 4, 4 );
  for(int i=0; i!=13; ++i) {
    h1[i] = hPSI->ProjectionY( Form("PSI0_DET_%d",i), i+1, i+1 );
    h2[i] = hPSI2->ProjectionY( Form("PSI1_DET_%d",i), i+1, i+1 );
    h3[i] = hPSI3->ProjectionY( Form("PSI2_DET_%d",i), i+1, i+1 );
    h1[i]->SetLineColor(kBlack);
    h2[i]->SetLineColor(kBlue-3);
    h3[i]->SetLineColor(kRed-3);
    cpsi->cd(i+1);
    h1[i]->GetYaxis()->SetRangeUser(0,h1[i]->GetMaximum());
    h1[i]->GetXaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetLabelSize(0.08);
    h1[i]->DrawCopy();
    h2[i]->DrawCopy("same");
    h3[i]->DrawCopy("same");
    tex->DrawLatex( 0.2, 10e4, textt[i].Data() );
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

  TCanvas *cv2raw = new TCanvas();
  cv2raw->Divide( 3, 4, 0, 0 );
  for(int i=0; i!=13; ++i) {
    h1[i] = h2V2[0]->ProjectionY( Form("V2RE_DET_%d",i), i+1, i+1 );
    h2[i] = h2V2[1]->ProjectionY( Form("V2RW_DET_%d",i), i+1, i+1 );
    h1[i]->SetLineColor(kBlue-3);
    h2[i]->SetLineColor(kRed-3);
    h1[i]->GetYaxis()->SetRangeUser(-0.005,0.02);
    if(i<1) continue;
    cv2raw->cd(i);
    h1[i]->DrawCopy();
    h2[i]->DrawCopy("same");
  }

  TCanvas *cv2rawall = new TCanvas();
  cv2rawall->Divide( 3, 4, 0, 0 );
  for(int i=0; i!=13; ++i) {
    h1[i] = h2V2All->ProjectionY( Form("V2R_DET_%d",i), i+1, i+1 );
    h1[i]->SetLineColor(kBlack);
    h1[i]->GetYaxis()->SetRangeUser(-0.005,0.02);
    if(i<1) continue;
    cv2rawall->cd(i);
    h1[i]->DrawCopy();
  }

  int color[13] = { kRed-3, kBlue-3, kMagenta-3, kMagenta-3,
		    kRed-3, kBlue-3, kGreen-3, kOrange-3,
		    kRed-3, kBlue-3, kGreen-3, kOrange-3,
		    kMagenta-3 };
  int style[13] = { 1,1,24,25,
		    1,1,1,1,
		    1,1,1,1,
		    27 };
  TCanvas *cv2 = new TCanvas();
  cv2->Divide( 3, 4, 0, 0 );
  for(int i=0; i!=13; ++i) {
    h1[i]->Scale(1.0/hreso->GetBinContent(i+1)/0.14);
    h1[i]->GetYaxis()->SetRangeUser(-0.005,0.03);
    h1[i]->SetLineColor( color[i] );
    h1[i]->SetMarkerColor( color[i] );
    h1[i]->SetMarkerStyle( style[i] );
    if(i<1) continue;
    cv2->cd(i);
    h1[i]->DrawCopy();
  }

  TCanvas *cres = new TCanvas();
  cres->Divide(2,1);
  cres->cd(1);
  hRes->Draw();
  cres->cd(2);
  hreso->Draw();

  TCanvas *cv2f = new TCanvas();
  cv2f->Divide(2,2);
  cv2f->cd(1);
  gFS->Draw("AP");
  h1[4]->Draw("SAME");
  h1[5]->Draw("SAME");
  h1[6]->Draw("SAME");
  h1[7]->Draw("SAME");
  cv2f->cd(2);
  gFS->Draw("AP");
  h1[8]->Draw("SAME");
  h1[9]->Draw("SAME");
  h1[10]->Draw("SAME");
  h1[11]->Draw("SAME");
  cv2f->cd(3);
  gFS->Draw("AP");
  h1[2]->Draw("SAME");
  h1[3]->Draw("SAME");
  h1[12]->Draw("SAME");
  cv2f->cd(4);
  gFS->Draw("AP");
  h1[0]->Draw("SAME");
  h1[1]->Draw("SAME");
  h1[2]->Draw("SAME");

  TString text[11] = {"#Psi^{BBC_{A}}_{2} :: #Psi^{BBC_{B}}_{2}",
		      "#Psi^{FVTX}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX0i}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX1i}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX2i}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX3i}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX0o}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX1o}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX2o}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EX3o}_{2} :: #Psi^{BBC}_{2}",
		      "#Psi^{EXfull}_{2} :: #Psi^{BBC}_{2}"};

  TCanvas *cepc0 = new TCanvas();
  cepc0->Divide( 3, 4, 0, 0 );
  cepc0->cd(1);
  for(int i=0; i!=11; ++i) {
    hEPC[0][i]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, text[i].Data() );
    cepc0->cd(i+3);
  }

  TCanvas *cepc1 = new TCanvas();
  cepc1->Divide( 3, 4, 0, 0 );
  cepc1->cd(1);
  for(int i=0; i!=11; ++i) {
    hEPC[1][i]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, text[i].Data() );
    cepc1->cd(i+3);
  }

  TCanvas *cepc2 = new TCanvas();
  cepc2->Divide( 3, 4, 0, 0 );
  cepc2->cd(1);
  for(int i=0; i!=11; ++i) {
    hEPC[2][i]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, text[i].Data() );
    cepc2->cd(i+3);
  }

  TString texttt[13] = { "#varphi - #Psi_{2}^{BBC_{A}}",
			 "#varphi - #Psi_{2}^{BBC_{B}}",
			 "#varphi - #Psi_{2}^{BBC}",
			 "#varphi - #Psi_{2}^{FVTX}",
			 "#varphi - #Psi_{2}^{EX0i}",
			 "#varphi - #Psi_{2}^{EX1i}",
			 "#varphi - #Psi_{2}^{EX2i}",
			 "#varphi - #Psi_{2}^{EX3i}",
			 "#varphi - #Psi_{2}^{EX0o}",
			 "#varphi - #Psi_{2}^{EX1o}",
			 "#varphi - #Psi_{2}^{EX2o}",
			 "#varphi - #Psi_{2}^{EX3o}",
			 "#varphi - #Psi_{2}^{EX}" };
  TCanvas *cqx = new TCanvas();
  cqx->Divide( 4, 6, 0, 0 );
  for(int i=0; i!=13; ++i) {
    cqx->cd(i+1);
    if(i>3) cqx->cd(i);
    h1[i] = (TH1D*) hDP[0]->ProjectionY( Form("dp0_%d",i), i+1, i+1 );
    h2[i] = (TH1D*) hDPPE[0]->ProjectionY( Form("dppe0_%d",i), i+1, i+1 );
    h1[i]->SetLineColor( kRed-3 );
    h2[i]->SetLineColor( kBlue-3 );
    h2[i]->GetYaxis()->SetRangeUser(46.5e3,49.5e3);
    h2[i]->GetXaxis()->SetLabelSize(0.08);
    h2[i]->GetYaxis()->SetLabelSize(0.08);
    h2[i]->GetYaxis()->SetNdivisions(505);
    h2[i]->DrawCopy();
    h1[i]->DrawCopy("same");
    tex->DrawLatex( -5.5, 48.5e3, texttt[i].Data() );
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
    h1[i]->GetXaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetRangeUser(0.95,1.05);
    h1[i]->GetYaxis()->SetNdivisions(507);
    h1[i]->DrawCopy();
    tex->DrawLatex( -5.5, 1.02, texttt[i].Data() );
  }

  TCanvas *cqy = new TCanvas();
  cqy->Divide( 4, 6, 0, 0 );
  for(int i=0; i!=13; ++i) {
    cqy->cd(i+1);
    if(i>3) cqy->cd(i);
    h1[i] = (TH1D*) hDP[1]->ProjectionY( Form("dp1_%d",i), i+1, i+1 );
    h2[i] = (TH1D*) hDPPE[1]->ProjectionY( Form("dppe1_%d",i), i+1, i+1 );
    h1[i]->SetLineColor( kRed-3 );
    h2[i]->SetLineColor( kBlue-3 );
    //h2[i]->GetYaxis()->SetRangeUser(46.5e3,49.5e3);
    h2[i]->GetXaxis()->SetLabelSize(0.08);
    h2[i]->GetYaxis()->SetLabelSize(0.08);
    h2[i]->GetYaxis()->SetNdivisions(505);
    h2[i]->DrawCopy();
    h1[i]->DrawCopy("same");
    tex->DrawLatex( -5.5, 32e3, texttt[i].Data() );
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
    h1[i]->GetXaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetLabelSize(0.08);
    h1[i]->GetYaxis()->SetRangeUser(0.95,1.05);
    h1[i]->GetYaxis()->SetNdivisions(507);
    h1[i]->DrawCopy();
    tex->DrawLatex( -5.5, 1.02, texttt[i].Data() );
  }

  return 0;
}
