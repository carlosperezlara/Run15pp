TH1F *tot;
TF1 *sgn;

Double_t SoverN(Double_t *x, Double_t *p) {
  Int_t bin = tot->GetXaxis()->FindBin( x[0] );
  Double_t xmin = tot->GetXaxis()->GetBinLowEdge(bin);
  Double_t wid = tot->GetXaxis()->GetBinWidth(bin);
  Double_t nnn = tot->GetBinContent( bin )*wid;
  Double_t sss = sgn->Integral(xmin,xmin+wid);
  return (sss/nnn);
}

Double_t myFunction(Double_t *x, Double_t *p) {
  // p[0] = v2sgn
  // p[1] = v2bgr0
  // p[2] = v2bgr1 (x)
  // p[3] = v2bgr2 (x^2)
  Double_t ret;
  Double_t v2bgr = p[1] + p[2]*x[0] + p[3]*(2*x[0]*x[0]-1);
  ret = v2bgr + SoverN(x,p) * (p[0]-v2bgr);
  return ret;
}

void style(TH1F *tmp, int color, float factor=1) {
  tmp->SetLineWidth(2);
  tmp->SetLineColor(color);
  tmp->GetYaxis()->SetTitleSize(0.15*factor);
  tmp->GetYaxis()->SetTitleOffset(0.3/factor);
  tmp->GetYaxis()->SetLabelSize(0.10*factor);
  tmp->GetXaxis()->SetTitleSize(0.15*factor);
  tmp->GetXaxis()->SetTitleOffset(1.0);
  tmp->GetXaxis()->SetLabelSize(0.15*factor);
  tmp->GetYaxis()->SetNdivisions(505);
  tmp->GetXaxis()->SetTitle("Mass  [GeV/c^{2}]");
  tmp->SetTitle("");
}

int flow() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile("out/all.root");
  int binsgn = 11;
  int binbgr = 16;
  const int nptbins = 11;
  float ptbins[15] = {1.0, 1.1, 1.2, 1.3, 1.5,
                      1.7, 1.9, 2.2, 2.5, 3.0,
                      4.0, 6.0, 10., 15., 20.};

  //float ptbins[15] = {1.0, 1.1, 1.2, 1.4, 1.6,
  //                    1.8, 2.0, 2.2, 2.5, 3.0,
  //                    4.0, 6.0, 10., 15., 20.};

  //const int nptbins = 6;
  //float ptbins[10] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10., 15., 20.};

  float ss[nptbins];
  float bb[nptbins];
  float vs[nptbins][5];
  float vb[nptbins][5];

  float vny1[5][nptbins], vne1[5][nptbins];
  float vny2[5][nptbins], vne2[5][nptbins];
  float vny3[5][nptbins], vne3[5][nptbins];

  float gaus0n[nptbins], gaus0ne[nptbins];
  float gaus0[nptbins], gaus0e[nptbins];
  float gaus1[nptbins], gaus1e[nptbins];
  float gaus2[nptbins], gaus2e[nptbins];

  TCanvas *main[nptbins];
  TCanvas *main2[nptbins];
  TLatex *tex = new TLatex;
  TLatex *tex2 = new TLatex;
  tex->SetTextSize(0.3);
  tex->SetTextColor(kOrange-3);
  tex2->SetTextSize(0.1);
  tex2->SetTextColor(kOrange-3);
  for(int ord=0; ord!=nptbins; ++ord) {
    main[ord] = new TCanvas( Form("CANVAS_PTB%d",ord), Form("CANVAS_PTB%d",ord) );
    main[ord]->SetLeftMargin(0.11);
    main[ord]->SetBottomMargin(0.3);
    main[ord]->Divide(1,6,0,0);
    main2[ord] = new TCanvas( Form("CANVAS2_PTB%d",ord), Form("CANVAS2_PTB%d",ord) );
    main2[ord]->SetLeftMargin(0.11);
    main2[ord]->SetBottomMargin(0.3);
    main2[ord]->Divide(1,2,0,0);
  } 
  TString str[5] = {"C(#varphi - #Psi_{1})",
		    "C2(#varphi - #Psi_{2})",
		    "C3(#varphi - #Psi_{3})",
		    "C4(#varphi - #Psi_{4})",
		    "C2(#varphi - #Psi_{4})" };
  float v2bin[5][2] = { {-0.05,+0.05},
			{-0.60,+0.60},
			{-0.60,+0.60},
			{-0.60,+0.60},
			{-0.05,+0.05} };
  for(int ptb=0; ptb!=nptbins; ++ptb) {
    TH1F *hm = (TH1F*) file->Get( Form("hMass_PB%d",ptb) );
    style(hm,kBlue-3,1.4);
    hm->GetYaxis()->SetTitle("counts");
    int bin10 = hm->GetXaxis()->FindBin( 0.050 );
    int bin20 = hm->GetXaxis()->FindBin( 0.100 );
    float counts1L = hm->Integral(bin10,bin20);
    int bin11 = hm->GetXaxis()->FindBin( 0.176 );
    int bin21 = hm->GetXaxis()->FindBin( 0.250 );
    float counts1R = hm->Integral(bin11,bin21);
    // BUILDING BGR
    TH1F *hmm = (TH1F*) file->Get( Form("hMass2_PB%d",ptb) );
    TH1F *hmmL = (TH1F*) hmm->Clone( Form("hMassL_PB%d",ptb) );
    TH1F *hmmR = (TH1F*) hmm->Clone( Form("hMassR_PB%d",ptb) );
    float counts2L = hmm->Integral(bin10,bin20);
    float counts2R = hmm->Integral(bin11,bin21);
    hmm->Reset();
    hmmL->Scale(counts1L/counts2L);
    hmmR->Scale(counts1R/counts2R);
    int thL = hm->GetXaxis()->FindBin( 0.110 );
    int thR = hm->GetXaxis()->FindBin( 0.166 );
    for(int iii=0; iii!=hmm->GetXaxis()->GetNbins(); ++iii) {
      float fracL, fracR;
      if(iii<thL) {
	fracL = 1.0;
	fracR = 0.0;
      } else if(iii>thR) {
	fracL = 0.0;
	fracR = 1.0;
      } else {
	fracL = (thR-iii)*1.0/(thR-thL);
	fracR = 1-fracL;
      }
      hmm->SetBinContent( iii,
			  fracL*hmmL->GetBinContent(iii) + 
			  fracR*hmmR->GetBinContent(iii) );
    }
    style(hmmL,kYellow-3);
    style(hmmR,kYellow-3);
    style(hmm,kRed-3);
    //

    main2[ptb]->cd(1);
    hm->DrawCopy();
    hmmL->DrawCopy("same");
    hmmR->DrawCopy("same");
    hmm->DrawCopy("same");
    tex2->DrawTextNDC( 0.6, 0.7, Form("pT  [  %.1f - %.1f ]  GeV/c",ptbins[ptb],ptbins[ptb+1]) );
    main2[ptb]->cd(2);
    TH1F *hmS = (TH1F*) hm->Clone( Form("hMassS_PB%d",ptb) );
    hmS->Add(hmm,-1);
    TF1 *fit = new TF1( Form("FIT_TB%d",ptb), "[0]*TMath::Gaus(x,[1],[2])",
			0.05, 0.25 );
    fit->SetParameter(0,100);
    fit->SetParameter(1,0.139);
    fit->SetParameter(2,0.01);  fit->SetParLimits(2,0.001,0.1);
    hmS->Fit( fit, "NWR", "", 0.115, 0.155 );
    gaus0n[ptb] = fit->GetParameter(0)/(ptbins[ptb+1]-ptbins[ptb]); gaus0ne[ptb] = fit->GetParError(0);
    gaus0[ptb] = fit->GetParameter(0); gaus0e[ptb] = fit->GetParError(0);
    gaus1[ptb] = fit->GetParameter(1); gaus1e[ptb] = fit->GetParError(1);
    gaus2[ptb] = fit->GetParameter(2); gaus2e[ptb] = fit->GetParError(2);

    hmS->DrawCopy();
    fit->Draw("same");
    tex2->DrawTextNDC( 0.6, 0.7, Form("mu = %.1f MeV/c2",1e3*fit->GetParameter(1)) );
    tex2->DrawTextNDC( 0.6, 0.6, Form("sgm = %.1f MeV/c2",1e3*fit->GetParameter(2)) );

    main[ptb]->cd(1);
    hm->DrawCopy();
    hmm->DrawCopy("same");
    tex->DrawTextNDC( 0.6, 0.7,
		      Form("pT  [  %.1f - %.1f ]  GeV/c",ptbins[ptb],ptbins[ptb+1]) );
    ss[ptb] = hm->GetBinContent( binsgn );
    bb[ptb] = hm->GetBinContent( binbgr );
    // setting global sources
    tot = hm;
    sgn = fit;
    TF1 *fitcheck = new TF1( Form("SoN%d",ptb), SoverN, 0,1 );
    main2[ptb]->cd(2);
    fitcheck->Draw("same");
    for(int ord=0; ord!=5; ++ord) {
      TProfile *hv = (TProfile*) file->Get( Form("hCos%dDP_PB%d",ord,ptb) );
      //hv->Rebin(5);
      vs[ptb][ord] = hv->GetBinContent( binsgn );
      vb[ptb][ord] = hv->GetBinContent( binbgr );
      main[ptb]->cd(2+ord);
      hv->GetYaxis()->SetTitle( Form( "%s",str[ord].Data() ) );
      style((TH1F*)hv,kGreen-3,ord==4?1.0:1.4);
      TF1 *fitv = new TF1( Form("fitvn_PT%d_ORD%d",ptb,ord), myFunction, 0.5, 0.25, 4 );
      hv->Fit(fitv);
      vny1[ord][ptb] = fitv->GetParameter(0);
      vne1[ord][ptb] = fitv->GetParError(0);
      fitv->SetParameter(3,0); fitv->SetParLimits(3,+1,-1);
      hv->Fit(fitv);
      vny2[ord][ptb] = fitv->GetParameter(0);
      vne2[ord][ptb] = fitv->GetParError(0);
      fitv->SetParameter(2,0); fitv->SetParLimits(2,+1,-1);
      fitv->SetParameter(3,0); fitv->SetParLimits(3,+1,-1);
      hv->Fit(fitv);
      vny3[ord][ptb] = fitv->GetParameter(0);
      vne3[ord][ptb] = fitv->GetParError(0);
      hv->DrawCopy();
    }
  }

  float zzz[20] = {0,0,0,0,0,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0};
  float xxx[nptbins];
  for(int i=0; i!=nptbins; ++i) xxx[i] = 0.5*(ptbins[i] + ptbins[i+1]);

  TCanvas *main3 = new TCanvas("pi0","pi0");
  main3->Divide(3,1);
  TGraphErrors *gg0n= new TGraphErrors(nptbins,xxx,gaus0n,zzz,gaus0ne);
  TGraphErrors *gg0 = new TGraphErrors(nptbins,xxx,gaus0,zzz,gaus0e);
  TGraphErrors *gg1 = new TGraphErrors(nptbins,xxx,gaus1,zzz,gaus1e);
  TGraphErrors *gg2 = new TGraphErrors(nptbins,xxx,gaus2,zzz,gaus2e);
  main3->cd(1)->SetLogy(1); gg0n->Draw("AP"); gg0n->SetMarkerStyle(20);
  main3->cd(1); gg0->Draw("*SAME");
  main3->cd(2); gg1->Draw("A*"); gg1->GetYaxis()->SetRangeUser(0.130,0.140);
  main3->cd(3); gg2->Draw("A*"); gg2->GetYaxis()->SetRangeUser(0,0.013);

  TCanvas *main4 = new TCanvas("vn","vn");
  main4->Divide(3,2);
  TGraphErrors *gv1[5], *gv2[5], *gv3[5];
  /*
  float ybin[5][2] = { {-0.04, 0},
		       {0, +0.02},
		       {-0.004, +0.004},
		       {-0.004, +0.004},
		       {-0.004, +0.004} };
  */
  float ybin[5][2] = { {-0.1, +0.1},
		       {-0.1, +0.1},
		       {-0.1, +0.1},
		       {-0.1, +0.1},
		       {-0.1, +0.1} };
  TString vss[5] = {"< Cos #varphi-#Psi_{1} >",
		    "< Cos 2(#varphi-#Psi_{2}) >",
		    "< Cos 3(#varphi-#Psi_{3}) >",
		    "< Cos 4(#varphi-#Psi_{4}) >",
		    "< Cos 2(#varphi-#Psi_{4}) >"};
  ofstream fout("vn.dat");
  for(int ord=0; ord!=5; ++ord) {
    fout << ord << " " << nptbins << endl;
    for(int ptb=0; ptb!=nptbins; ++ptb) {
      fout << ptbins[ptb] << " " << ptbins[ptb+1] << " ";
      fout << vny1[ord][ptb] << " " << vne1[ord][ptb] << " ";
      fout << vny2[ord][ptb] << " " << vne2[ord][ptb] << " ";
      fout << vny3[ord][ptb] << " " << vne3[ord][ptb] << endl;
    }
  }
  fout.close();
  for(int ord=0; ord!=5; ++ord) {
    gv1[ord] = new TGraphErrors(nptbins,xxx,vny1[ord],zzz,vne1[ord]);
    gv2[ord] = new TGraphErrors(nptbins,xxx,vny2[ord],zzz,zzz);//vne2[ord]);
    gv3[ord] = new TGraphErrors(nptbins,xxx,vny3[ord],zzz,zzz);///vne3[ord]);
    gv1[ord]->SetMarkerStyle(20);
    gv2[ord]->SetMarkerStyle(20);//24);
    gv3[ord]->SetMarkerStyle(20);//24);
    gv2[ord]->SetMarkerColor(kOrange-3);
    gv3[ord]->SetMarkerColor(kOrange-3);
    main4->cd(ord+1);
    TH2F *axis = new TH2F( Form("tmp%d",ord),
			   Form("%s;p_{T}  [GeV/c]",vss[ord].Data()),
			   100,0,20, 100,ybin[ord][0],ybin[ord][1]);
    axis->Draw();
    gv1[ord]->Draw("PSAME");
    gv2[ord]->Draw("PSAME");
    gv3[ord]->Draw("PSAME");
  }
  return 0;
}
