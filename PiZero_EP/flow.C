//float res[5] = {0.161,0.107,0.064,0.042,0.107};
//float res[5] = {0.161,0.107*1.2,0.064,0.042,0.107};
float res[5] = {1,1,1,1,1};
const int nptbins = 22; //14
float ptbins[23] = {1.0, 1.1, 1.2, 1.3, 1.4,
		    1.5, 1.6, 1.7, 1.8, 2.0,
		    2.2, 2.4, 2.6, 2.8, 3.0,
		    3.4, 3.8, 4.4, 5.0, 6.0,
		    8.0, 10., 20.};
//float ptbins[15] = {1.0, 1.1, 1.2, 1.3, 1.5,
//		    1.7, 1.9, 2.2, 2.5, 3.0,
//		    4.0, 6.0, 10., 15., 20.};
float meanPi0;

TH1F *tot; // they must be loaded
TF1 *sgn; // an actual fit is attempted
TFile *file;

TH1F *hm, *hmm, *hmmL, *hmmR;
TH1F *hmS;
TF1 *fit;
TProfile *hv;
TF1 *fitv;

Double_t SoverN(Double_t *x, Double_t *p) { // requires TOT and SGN
  Int_t bin = tot->GetXaxis()->FindBin( x[0] );
  Double_t xmin = tot->GetXaxis()->GetBinLowEdge(bin);
  Double_t wid = tot->GetXaxis()->GetBinWidth(bin);
  Double_t nnn = tot->GetBinContent( bin )*wid;
  Double_t sss = sgn->Integral(xmin,xmin+wid);
  return (sss/nnn);
}
//=======================================
Double_t Significance(double xmin, double xmax) { // requires TOT and SGN
  Int_t bin0 = tot->GetXaxis()->FindBin( xmin+1e-6 );
  Int_t bin1 = tot->GetXaxis()->FindBin( xmax-1e-6 );
  Double_t nnn = tot->Integral(bin0,bin1);
  Double_t sss = sgn->Integral(xmin,xmax);
  return (sss/TMath::Sqrt(nnn));
}
//=======================================
Double_t myFunction(Double_t *x, Double_t *p) { //requires SoverN
  // p[0] = v2sgn
  // p[1] = v2bgr0
  // p[2] = v2bgr1 (x)
  // p[3] = v2bgr2 (x)^2
  Double_t ret;
  Double_t diff = x[0] - meanPi0;
  Double_t v2bgr = p[1] + p[2]*diff + p[3]*(2*diff*diff-1);
  //Double_t v2bgr = p[1] + p[2]*TMath::Power(diff,p[3]);
  ret = v2bgr + SoverN(x,p) * (p[0]-v2bgr);
  return ret;
}
//=======================================
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
//=======================================
void color(TGraph *g, int c, int t=20) {
  g->SetMarkerColor(c);
  g->SetLineColor(c);
  g->SetFillColor(kWhite);
  g->SetMarkerStyle(t);
}
//=======================================
void LoadVn( TString cut, float v[2][10][50], float er[2][10][50]) { //v[0MB;1ERT][ORD][PTB] //er[0MB;1ERT][ORD][PTB]
  float trash;
  ifstream fin;
  //========
  fin.open( Form("CN%s.dat",cut.Data()) );
  for(;;) {
    int nord, npts;
    fin >> nord >> npts;
    if(!fin.good()) break;
    for(int k=0; k!=npts; ++k) {
      fin >> trash >> trash >> trash >> trash >> v[0][nord][k] >> er[1][nord][k] >> trash >> trash;
    }
  }
  fin.close();
  //========
  fin.open( Form("CNERT%s.dat",cut.Data()) );
  for(;;) {
    int nord, npts;
    fin >> nord >> npts;
    if(!fin.good()) break;
    for(int k=0; k!=npts; ++k) {
      fin >> trash >> trash >> trash >> trash >> v[1][nord][k] >> er[1][nord][k] >> trash >> trash;
    }
  }
  fin.close();
}
//=======================================
void LoadJustSys( TString cut, float v[2][10][50]) { //v[0MB;1ERT][ORD][PTB]
  float er[2][10][50];
  LoadVn( cut, v, er);
}
//=======================================
void LoadFile(TString cut) {
  file = new TFile( Form("all%s.root",cut.Data()) );
}
//=======================================
void CreateBgr(int ptb) {
  hm = (TH1F*) file->Get( Form("hMass_PB%d",ptb) );
  hm->GetYaxis()->SetTitle("counts");
  int bin10 = hm->GetXaxis()->FindBin( 0.050 );
  int bin20 = hm->GetXaxis()->FindBin( 0.100 );
  float counts1L = hm->Integral(bin10,bin20);
  int bin11 = hm->GetXaxis()->FindBin( 0.176 );
  int bin21 = hm->GetXaxis()->FindBin( 0.250 );
  float counts1R = hm->Integral(bin11,bin21);
  // BUILDING BGR
  hmm = (TH1F*) file->Get( Form("hMass2_PB%d",ptb) );
  hmmL = (TH1F*) hmm->Clone( Form("hMassL_PB%d",ptb) );
  hmmR = (TH1F*) hmm->Clone( Form("hMassR_PB%d",ptb) );
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
  style(hm,kBlue-3,0.7);
  style(hmmL,kYellow-3);
  style(hmmR,kYellow-3);
  style(hmm,kRed-3);
}
//=======================================
void FitSignal(int ptb ) {
  hmS = (TH1F*) hm->Clone( Form("hMassS_PB%d",ptb) );
  hmS->Add(hmm,-1);
  fit = new TF1( Form("FIT_TB%d",ptb), "[0]*TMath::Gaus(x,[1],[2])",
		 0.05, 0.25 );
  fit->SetParameter(0,100);
  fit->SetParameter(1,0.139);
  fit->SetParameter(2,0.01);  fit->SetParLimits(2,0.001,0.1);
  hmS->Fit( fit, "NWR", "", 0.115, 0.155 );
  hmS->Fit( fit, "NWR", "", 0.115, 0.155 );
  meanPi0 = fit->GetParameter(1);
}
//=======================================
void LoadVnProfileAndFit(int ptb, int ord, int rebin=1) {
  hv = (TProfile*) file->Get( Form("hCos%dDP_PB%d",ord,ptb) );
  hv->Rebin( rebin );
  hv->Dump();
  fitv = new TF1( Form("fitvn_PT%d_ORD%d",ptb,ord), myFunction, 0.5, 0.25, 4 );
}
//=======================================
void FitVn(bool quiet=true,TString opt="Q",float min=-1,float max=-1) {
  if(min<0) {
    hv->Fit(fitv,opt.Data());
  } else {
    hv->Fit(fitv,opt.Data(),"",min,max);
  }
  if(!quiet) {
    cout << "FIT RESULTS || vsgn evsgn chi2/ndf" << endl;
    cout << fitv->GetParameter(0);
    cout << " ";
    cout << fitv->GetParError(0);
    cout << " ";
    cout << fitv->GetChisquare() / fitv->GetNDF();
    cout << endl;
  }
}
//=======================================
void FitVnConst(bool quiet=true,float min=-1,float max=-1) {
  fitv->SetParameter(2,0); fitv->SetParLimits(2,+1,-1); //fixed
  fitv->SetParameter(3,0); fitv->SetParLimits(3,+1,-1); //fixed
  FitVn(true,"QR",min,max);
  FitVn(quiet,"R",min,max);
}
//=======================================
void FitVnLine(bool quiet=true,float min=-1,float max=-1) {
  FitVnConst(true,min,max);
  fitv->SetParLimits(2,-1,+1); //free
  fitv->SetParLimits(3,+1,-1); //fixed
  FitVn(true,"QR",min,max);
  FitVn(quiet,"R",min,max);
}
//=======================================
void FitVnQuad(bool quiet=true,float min=-1,float max=-1) {
  FitVnLine(true,min,max);
  fitv->SetParLimits(2,-1,+1); // free
  //fitv->SetParLimits(3,-1,+1); // free
  fitv->SetParLimits(3,+1,-1); // still fixed
  FitVn(true,"QR",min,max);
  FitVn(quiet,"R",min,max);
}
//=======================================
int flow(TString cut="",int ptfirst=0, int ptlast=11) {
  gStyle->SetOptStat(0);
  LoadFile(cut);
  float vny1[5][nptbins], vne1[5][nptbins];
  float vny2[5][nptbins], vne2[5][nptbins];
  float vny3[5][nptbins], vne3[5][nptbins];

  float gaus0n[nptbins], gaus0ne[nptbins];
  float gaus0[nptbins], gaus0e[nptbins];
  float gaus1[nptbins], gaus1e[nptbins];
  float gaus2[nptbins], gaus2e[nptbins];

  TCanvas *main[nptbins];
  TCanvas *main2[nptbins];
  TCanvas *main3[nptbins];
  TLatex *tex = new TLatex;
  TLatex *tex2 = new TLatex;
  tex->SetTextSize(0.3);
  tex->SetTextColor(kOrange-3);
  tex2->SetTextSize(0.1);
  tex2->SetTextColor(kOrange-3);
  for(int ptb=ptfirst; ptb!=ptlast; ++ptb) {
    main[ptb] = new TCanvas( Form("pdf/%s_CANVAS_PTB%d",cut.Data(), ptb), Form("CANVAS_PTB%d",ptb) );
    main[ptb]->SetLeftMargin(0.11);
    main[ptb]->SetBottomMargin(0.3);
    main[ptb]->Divide(1,6,0,0);
    main2[ptb] = new TCanvas( Form("pdf/%s_CANVAS2_PTB%d",cut.Data(),ptb), Form("CANVAS2_PTB%d",ptb) );
    main2[ptb]->SetLeftMargin(0.10);
    main2[ptb]->SetBottomMargin(0.2);
    main2[ptb]->Divide(1,2,0,0);
    main3[ptb] = new TCanvas( Form("pdf/%s_CANVAS3_PTB%d",cut.Data(),ptb), Form("CANVAS3_PTB%d",ptb) );
    main3[ptb]->SetLeftMargin(0.11);
    main3[ptb]->SetBottomMargin(0.3);
    main3[ptb]->Divide(1,3,0,0);
  } 
  TString str[5] = {"C(#varphi - #Psi_{1})",
		    "C2(#varphi - #Psi_{2})",
		    "C3(#varphi - #Psi_{3})",
		    "C4(#varphi - #Psi_{4})",
		    "C4(#varphi - #Psi_{2})" };
  float v2bin[5][2] = { {-0.05,+0.05},
			{-0.60,+0.60},
			{-0.60,+0.60},
			{-0.60,+0.60},
			{-0.05,+0.05} };
  ofstream fsgnout( Form("dat/SGN%s.dat",cut.Data()) );
  for(int ptb=ptfirst; ptb!=ptlast; ++ptb) {
    CreateBgr(ptb);//,hm,hmm,hmmL,hmmR);
    FitSignal(ptb);//,hm,hmm,hmS,fit);

    // setting global sources
    tot = hm;
    sgn = fit;

    main2[ptb]->cd(1);
    hm->DrawCopy();
    hmmL->DrawCopy("same");
    hmmR->DrawCopy("same");
    hmm->DrawCopy("same");
    tex2->DrawTextNDC( 0.6, 0.7, Form("pT  [  %.1f - %.1f ]  GeV/c",ptbins[ptb],ptbins[ptb+1]) );
    main2[ptb]->cd(2);
    gaus0n[ptb] = fit->GetParameter(0)/(ptbins[ptb+1]-ptbins[ptb]); gaus0ne[ptb] = fit->GetParError(0);
    gaus0[ptb] = fit->GetParameter(0); gaus0e[ptb] = fit->GetParError(0);
    gaus1[ptb] = fit->GetParameter(1); gaus1e[ptb] = fit->GetParError(1);
    gaus2[ptb] = fit->GetParameter(2); gaus2e[ptb] = fit->GetParError(2);

    hmS->DrawCopy();
    fit->Draw("same");
    tex2->DrawTextNDC( 0.6, 0.7, Form("mu = %.1f MeV/c2",1e3*fit->GetParameter(1)) );
    tex2->DrawTextNDC( 0.6, 0.6, Form("sgm = %.1f MeV/c2",1e3*fit->GetParameter(2)) );

    main3[ptb]->cd(1);
    style(hm,kBlue-3,0.7);
    hm->DrawCopy();
    hmm->DrawCopy("same");
    tex2->SetTextSize(0.2);
    tex2->DrawTextNDC( 0.6, 0.7, Form("pT  [  %.1f - %.1f ]  GeV/c",ptbins[ptb],ptbins[ptb+1]) );
    tex2->SetTextSize(0.1);
  
    main[ptb]->cd(1);
    style(hm,kBlue-3,1.4);
    hm->DrawCopy();
    hmm->DrawCopy("same");
    tex->DrawTextNDC( 0.6, 0.7, Form("pT  [  %.1f - %.1f ]  GeV/c",ptbins[ptb],ptbins[ptb+1]) );
    Double_t signi = Significance( gaus1[ptb] - 3*gaus2[ptb], gaus1[ptb] + 3*gaus2[ptb] );

    fsgnout << ptbins[ptb] << " " << ptbins[ptb+1] << " " << signi << " ";
    fsgnout << sgn->Integral( gaus1[ptb] - 3*gaus2[ptb], gaus1[ptb] + 3*gaus2[ptb] ) << endl;

    TF1 *fitcheck = new TF1( Form("SoN%d",ptb), SoverN, 0,1 );
    main2[ptb]->cd(2);
    fitcheck->Draw("same");
    cout << "******* PT BIN " << ptb << " [" << ptbins[ptb] << ", " << ptbins[ptb+1] << "]" << endl;
    for(int ord=0; ord!=5; ++ord) {
      cout << "    ** Fitting ORD BIN " << ord << endl;
      LoadVnProfileAndFit(ptb,ord); // hv,fitv
      //hv->Rebin(5);
      main[ptb]->cd(2+ord);
      hv->GetYaxis()->SetTitle( Form( "%s",str[ord].Data() ) );
      //==
      FitVnLine();
      vny2[ord][ptb] = fitv->GetParameter(0)/res[ord];
      vne2[ord][ptb] = fitv->GetParError(0)/res[ord];
      //==
      FitVnConst();
      vny3[ord][ptb] = fitv->GetParameter(0)/res[ord];
      vne3[ord][ptb] = fitv->GetParError(0)/res[ord];
      //==
      continue;
      FitVnQuad();
      vny1[ord][ptb] = fitv->GetParameter(0)/res[ord];
      vne1[ord][ptb] = fitv->GetParError(0)/res[ord];
      style((TH1F*)hv,kGreen-3,ord==4?1.0:1.4);
      hv->DrawCopy();
      if(ord>2) continue;
      main3[ptb]->cd(2+ord);
      style((TH1F*)hv,kGreen-3,ord==1?0.60:0.90);
      hv->DrawCopy();
    }
    main[ptb]->SaveAs(  Form("%s.pdf", main[ptb]->GetName()  ), "pdf"  );
    main2[ptb]->SaveAs( Form("%s.pdf", main2[ptb]->GetName() ), "pdf"  );
    main3[ptb]->SaveAs( Form("%s.pdf", main3[ptb]->GetName() ), "pdf"  );
  }
  int nbins = ptlast - ptfirst;
  fsgnout.close();

  //STORING COEFFICIENTS
  ofstream fout( Form("dat/CN%s.dat",cut.Data()) );
  for(int ord=0; ord!=5; ++ord) {
    fout << ord << " " << nbins << endl;
    for(int ptb=ptfirst; ptb!=ptlast; ++ptb) {
      fout << ptbins[ptb] << " " << ptbins[ptb+1] << " ";
      fout << vny1[ord][ptb] << " " << vne1[ord][ptb] << " ";
      fout << vny2[ord][ptb] << " " << vne2[ord][ptb] << " ";
      fout << vny3[ord][ptb] << " " << vne3[ord][ptb] << endl;
    }
  }
  fout.close();

  float zzz[30] = {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0};
  float xxx[30];
  for(int i=ptfirst; i!=ptlast; ++i) xxx[i] = 0.5*(ptbins[i] + ptbins[i+1]);

  TCanvas *main4 = new TCanvas(Form("pdf/%s_pi0",cut.Data()),"pi0");
  main4->Divide(3,1);
  TGraphErrors *gg0n= new TGraphErrors(nbins,xxx,gaus0n,zzz,gaus0ne);
  TGraphErrors *gg0 = new TGraphErrors(nbins,xxx,gaus0,zzz,gaus0e);
  TGraphErrors *gg1 = new TGraphErrors(nbins,xxx,gaus1,zzz,gaus1e);
  TGraphErrors *gg2 = new TGraphErrors(nbins,xxx,gaus2,zzz,gaus2e);
  main4->cd(1)->SetLogy(1); gg0n->Draw("AP"); gg0n->SetMarkerStyle(20);
  main4->cd(1); gg0->Draw("*SAME");
  main4->cd(2); gg1->Draw("A*"); gg1->GetYaxis()->SetRangeUser(0.130,0.140);
  main4->cd(3); gg2->Draw("A*"); gg2->GetYaxis()->SetRangeUser(0,0.013);
  main4->SaveAs(  Form("%s.pdf", main4->GetName()  ), "pdf"  );

  TCanvas *main5 = new TCanvas(Form("pdf/%s_vn",cut.Data()),"vn");
  main5->Divide(3,2);

  TCanvas *main6 = new TCanvas(Form("pdf/%s_v12",cut.Data()),"v12");
  main6->Divide(2,1);

  TGraphErrors *gv1[5], *gv2[5], *gv3[5];
  float ybin[5][2] = { {-0.04, 0},
		       {-0.01, +0.01},
		       {-0.004, +0.004},
		       {-0.004, +0.004},
		       {-0.004, +0.004} };
  TString vss[5] = {"< Cos #varphi-#Psi_{1} >",
		    "< Cos 2(#varphi-#Psi_{2}) >",
		    "< Cos 3(#varphi-#Psi_{3}) >",
		    "< Cos 4(#varphi-#Psi_{4}) >",
		    "< Cos 4(#varphi-#Psi_{2}) >"};
  for(int ord=0; ord!=5; ++ord) {
    gv1[ord] = new TGraphErrors(nbins,xxx,vny1[ord],zzz,vne1[ord]);
    gv2[ord] = new TGraphErrors(nbins,xxx,vny2[ord],zzz,zzz);//vne2[ord]);
    gv3[ord] = new TGraphErrors(nbins,xxx,vny3[ord],zzz,zzz);///vne3[ord]);
    gv1[ord]->SetMarkerStyle(20);
    gv2[ord]->SetMarkerStyle(20);//24);
    gv3[ord]->SetMarkerStyle(20);//24);
    gv2[ord]->SetMarkerColor(kOrange-3);
    gv3[ord]->SetMarkerColor(kOrange-3);
    main5->cd(ord+1);
    TH2F *axis = new TH2F( Form("tmp%d",ord),
			   Form("%s;p_{T}  [GeV/c]",vss[ord].Data()),
			   100,0,20.0, 100,ybin[ord][0],ybin[ord][1]);
    axis->Draw();
    gv1[ord]->Draw("PSAME");
    gv2[ord]->Draw("PSAME");
    gv3[ord]->Draw("PSAME");
    if(ord>1) continue;
    main6->cd(ord+1);
    axis->Draw();
    gv1[ord]->Draw("PSAME");
    gv2[ord]->Draw("PSAME");
    gv3[ord]->Draw("PSAME");
  }
  main5->SaveAs(  Form("%s.pdf", main5->GetName()  ), "pdf"  );

  return 0;
}
//=======================================
void LoadV2Charged() {
  ifstream fin("v2charged.dat");
  for(;;) {
    //fin >> 
  }
}
//=======================================
int fitme(TString cut="NOM",int ptb=0,int ord=1,int vbgr=1,int rebin=1,int skippre=0, int skippos=0) {
  gStyle->SetOptStat(0);
  cout << ptbins[ptb] << " " << ptbins[ptb+1] << endl;
  LoadFile(cut);
  CreateBgr(ptb);//,hm,hmm,hmmL,hmmR);
  TCanvas *main = new TCanvas();
  main->Divide(1,2);
  main->cd(2);
  FitSignal(ptb);//,hm,hmm,hmS,fit);
  hm->Draw();
  // setting global sources
  tot = hm;
  sgn = fit;
  LoadVnProfileAndFit(ptb,ord,rebin); // hv,fitv

  int nbm = hv->GetXaxis()->GetNbins()+1;
  float minr = hv->GetXaxis()->GetBinLowEdge( 1 );
  float maxr = hv->GetXaxis()->GetBinLowEdge( nbm );
  if(skippre>0) minr = hv->GetXaxis()->GetBinLowEdge( 1+skippre );
  if(skippos>0) maxr = hv->GetXaxis()->GetBinLowEdge( nbm-skippos );

  main->cd(1);
  switch(vbgr) {
  case(0):
    FitVnConst(false,minr,maxr);
    break;
  case(1):
    FitVnLine(false,minr,maxr);
    break;
  case(2):
    FitVnQuad(false,minr,maxr);
    break;
  }
  main->SaveAs( Form("pdf/FLOW_%s%s_ptb%d_ord%d.pdf",cut.Data(),vbgr!=1?Form("%d",vbgr):"",ptb,ord) );

  ofstream fout( Form("dat/CN_%s%s_ptb%d_ord%d.dat",cut.Data(),vbgr!=1?Form("%d",vbgr):"",ptb,ord) ); //==================== OUTPUT
  fout << Form("%s%s",cut.Data(),vbgr!=1?Form("%d",vbgr):"") << " "; //===================================================== OUTPUT
  fout << ptbins[ptb] << " " << ptbins[ptb+1] << " "; //==================================================================== OUTPUT
  fout << fitv->GetParameter(0) << " " << fitv->GetParError(0) << " "; //=================================================== OUTPUT
  fout << fitv->GetParameter(1) << " " << fitv->GetParError(1) << " "; //=================================================== OUTPUT
  fout << fitv->GetParameter(2) << " " << fitv->GetParError(2) << " "; //=================================================== OUTPUT
  fout << fitv->GetParameter(3) << " " << fitv->GetParError(3) << " "; //=================================================== OUTPUT
  fout << fitv->GetChisquare() / fitv->GetNDF() << endl; //================================================================= OUTPUT
  fout.close();
  fout.open( Form("dat/SGN_%s%s_ptb%d_ord%d.dat",cut.Data(),vbgr!=1?Form("%d",vbgr):"",ptb,ord) ); //======================= OUTPUT
  fout << ptbins[ptb] << " " << ptbins[ptb+1] << " "; //==================================================================== OUTPUT
  float threem = fit->GetParameter(1) - 3*fit->GetParameter(2);
  float threep = fit->GetParameter(1) + 3*fit->GetParameter(2);
  fout << fit->GetParameter(0) << " " << fit->GetParameter(1) << " " << fit->GetParameter(2) << " "; //===================== OUTPUT
  fout << Significance( threem, threep ) << " " << sgn->Integral( threem, threep ) << endl; //============================== OUTPUT
  fout.close();

  TCanvas *main2 = new TCanvas();
  tot->GetYaxis()->SetRangeUser( 0, tot->GetMaximum()*1.5 );
  tot->GetYaxis()->SetTitleOffset(1.0);
  tot->GetYaxis()->SetLabelSize(0.04);
  tot->GetYaxis()->SetTitleSize(0.04);
  tot->GetXaxis()->SetLabelSize(0.04);
  tot->GetXaxis()->SetTitleSize(0.04);
  tot->Draw();
  sgn->Draw("SAME");
  main2->SaveAs( Form("pdf/TOT_%s%s_ptb%d_ord%d.pdf",cut.Data(),vbgr!=1?Form("%d",vbgr):"",ptb,ord) );

  main->cd(1);
  return 0;
}
//=======================================
void sys() {
  gStyle->SetOptStat(0);
  float pt[2][50], ept[2][50], minpt[2][50], maxpt[2][50], v[2][10][50], ev[2][10][50], vF2[2][10][50], vF0[2][10][50];
  float zzz[50];
  for(int i=0; i!=50; ++i) zzz[i] = 0;
  float trash;

  ifstream fin;
  int n=0, nERT=0;
  int nord, npts;
  //========
  fin.open("CN.dat");
  for(;;) {
    fin >> nord >> npts;
    if(!fin.good()) break;
    n = npts;
    for(int k=0; k!=npts; ++k) {
      fin >> minpt[0][k] >> maxpt[0][k] >> vF2[0][nord][k] >> trash >> v[0][nord][k] >> ev[0][nord][k] >> vF0[0][nord][k] >> trash;
      pt[0][k] = 0.5*(maxpt[0][k] + minpt[0][k]);
      ept[0][k] = 0.5*(maxpt[0][k] - minpt[0][k]);
    }
  }
  fin.close();
  //========
  fin.open("CNERT.dat");
  for(;;) {
    fin >> nord >> npts;
    if(!fin.good()) break;
    nERT = npts;
    for(int k=0; k!=npts; ++k) {
      fin >> minpt[1][k] >> maxpt[1][k] >> vF2[1][nord][k] >> trash >> v[1][nord][k] >> ev[1][nord][k] >> vF0[1][nord][k] >> trash;
      pt[1][k] = 0.5*(maxpt[1][k] + minpt[1][k]);
      ept[1][k] = 0.5*(maxpt[1][k] - minpt[1][k]);
    }
  }
  fin.close();

  float env[50];
  //for(int i=0; i!=50; ++i) env[i] = 0.0002;
  for(int i=0; i!=50; ++i) env[i] = 0.0005;

  TCanvas *main = new TCanvas();
  main->SetLeftMargin(0.12);
  TH2F *av2 = new TH2F("av2",";p_{T}  (GeV/c);#pi^{0}  < Cos 2(#varphi - #Psi_{2}) >",100,0,10,100,0,0.025);
  TGraphErrors *gv2 = new TGraphErrors(n,pt[0],v[0][1],ept[0],ev[0][1]);
  TGraphErrors *gv2ERT = new TGraphErrors(nERT,pt[1],v[1][1],ept[1],ev[1][1]);
  TGraphErrors *gv2env = new TGraphErrors(n,pt[0],v[0][1],ept[0],env);
  gv2env->SetFillColor(kGray);
  gv2env->SetLineWidth(0);

  gv2->SetMarkerStyle(20);
  gv2ERT->SetMarkerStyle(20);
  gv2->SetFillColor(kWhite);
  gv2ERT->SetFillColor(kWhite);
  gv2ERT->SetLineColor(kRed-3);
  gv2ERT->SetMarkerColor(kRed-3);
  av2->Draw();
  gv2env->Draw("E2");
  gv2->Draw("SAMEP");
  gv2ERT->Draw("SAMEP");
  TLegend *leg = new TLegend(0.15,0.15,0.4,0.3);
  leg->AddEntry(gv2,"MBT 0-5%");
  leg->AddEntry(gv2ERT,"ERT 0-5%");
  leg->Draw();

  float vA0[2][10][50], vA1[2][10][50], vD0[2][10][50], vD1[2][10][50], vT0[2][10][50], vT1[2][10][50];
  float vFA0[2][10][50], vFA1[2][10][50], vFD0[2][10][50], vFD1[2][10][50], vFT0[2][10][50], vFT1[2][10][50];
  float vSX[2][10][50], vSXX[2][10][50];

  LoadJustSys("A0",vA0);
  LoadJustSys("A1",vA1);
  LoadJustSys("D0",vD0);
  LoadJustSys("D1",vD1);
  LoadJustSys("T0",vT0);
  LoadJustSys("T1",vT1);

  LoadJustSys("FA0",vFA0);
  LoadJustSys("FA1",vFA1);
  LoadJustSys("FD0",vFD0);
  LoadJustSys("FD1",vFD1);
  LoadJustSys("FT0",vFT0);
  LoadJustSys("FT1",vFT1);

  for(int k=0; k!=n; ++k) {
    vA0[0][1][k] -= v[0][1][k];
    vA1[0][1][k] -= v[0][1][k];
    vD0[0][1][k] -= v[0][1][k];
    vD1[0][1][k] -= v[0][1][k];
    vT0[0][1][k] -= v[0][1][k];
    vT1[0][1][k] -= v[0][1][k];

    vF0[0][1][k] -= v[0][1][k];
    vF2[0][1][k] -= v[0][1][k];

    vFA0[0][1][k] -= v[0][1][k];
    vFA1[0][1][k] -= v[0][1][k];
    vFD0[0][1][k] -= v[0][1][k];
    vFD1[0][1][k] -= v[0][1][k];
    vFT0[0][1][k] -= v[0][1][k];
    vFT1[0][1][k] -= v[0][1][k];

    vSX[0][1][k] = (vA0[0][1][k] + vA1[0][1][k] + vD0[0][1][k] + vD1[0][1][k] +
		    vT0[0][1][k] + vT1[0][1][k] +
		    vFA0[0][1][k]+vFA1[0][1][k] +vFD0[0][1][k] +vFD1[0][1][k] +
		    vFT0[0][1][k]+vFT1[0][1][k] + 
		    vF0[0][1][k] + vF2[0][1][k])/14.0;
    vSXX[0][1][k] = (vA0[0][1][k]*vA0[0][1][k]  +
		     vA1[0][1][k]*vA1[0][1][k] +
		     vD0[0][1][k]*vD0[0][1][k] +
		     vD1[0][1][k]*vD1[0][1][k] + 
		     vT0[0][1][k]*vT0[0][1][k] +
		     vT1[0][1][k]*vT1[0][1][k] +
		     vFA0[0][1][k]*vFA0[0][1][k]  +
		     vFA1[0][1][k]*vFA1[0][1][k] +
		     vFD0[0][1][k]*vFD0[0][1][k] +
		     vFD1[0][1][k]*vFD1[0][1][k] + 
		     vFT0[0][1][k]*vFT0[0][1][k] +
		     vFT1[0][1][k]*vFT1[0][1][k] +
		     vF0[0][1][k]*vF0[0][1][k] +
		     vF2[0][1][k]*vF2[0][1][k])/14.0;
    vSXX[0][1][k] = TMath::Sqrt( vSXX[0][1][k] );
  }

  TGraph *gA0 = new TGraph(n,pt[0],vA0[0][1]);
  TGraph *gA1 = new TGraph(n,pt[0],vA1[0][1]);
  TGraph *gD0 = new TGraph(n,pt[0],vD0[0][1]);
  TGraph *gD1 = new TGraph(n,pt[0],vD1[0][1]);
  TGraph *gT0 = new TGraph(n,pt[0],vT0[0][1]);
  TGraph *gT1 = new TGraph(n,pt[0],vT1[0][1]);

  TGraph *gFA0 = new TGraph(n,pt[0],vFA0[0][1]);
  TGraph *gFA1 = new TGraph(n,pt[0],vFA1[0][1]);
  TGraph *gFD0 = new TGraph(n,pt[0],vFD0[0][1]);
  TGraph *gFD1 = new TGraph(n,pt[0],vFD1[0][1]);
  TGraph *gFT0 = new TGraph(n,pt[0],vFT0[0][1]);
  TGraph *gFT1 = new TGraph(n,pt[0],vFT1[0][1]);

  TGraph *gF0 = new TGraph(n,pt[0],vF0[0][1]);
  TGraph *gF2 = new TGraph(n,pt[0],vF2[0][1]);

  TGraph *gSX = new TGraph(n,pt[0],vSX[0][1]);
  TGraph *gSXX = new TGraph(n,pt[0],vSXX[0][1]);

  //TGraph *gA0ERT = new TGraph(nERT,pt[1],vA0[1][1],ept[1],zzz);
  color(gA0,kGreen-3,24);
  color(gA1,kGreen-3,20);
  color(gD0,kOrange-3,24);
  color(gD1,kOrange-3,20);
  color(gT0,kGray,24);
  color(gT1,kGray,20);

  color(gFA0,kGreen-3,24);
  color(gFA1,kGreen-3,20);
  color(gFD0,kOrange-3,24);
  color(gFD1,kOrange-3,20);
  color(gFT0,kGray,24);
  color(gFT1,kGray,20);

  color(gF0,kCyan-3,24);
  color(gF2,kCyan-3,20);
  color(gSX,kBlue-3,24);
  color(gSXX,kBlue-3,20);

  TCanvas *msys2 = new TCanvas("msys2");
  msys2->Divide(2,1);
  msys2->cd(1);
  TH2F *as2 = new TH2F("as2",";p_{T}  (GeV/c);< Cos 2(#varphi - #Psi_{2}) >_{#pi^{0}}",100,0,10,100,-0.002,+0.002);
  as2->Draw();
  gA0->Draw("SAMEPL");
  gA1->Draw("SAMEPL");
  gD0->Draw("SAMEPL");
  gD1->Draw("SAMEPL");
  gT0->Draw("SAMEPL");
  gT1->Draw("SAMEPL");
  gFA0->Draw("SAMEPL");
  gFA1->Draw("SAMEPL");
  gFD0->Draw("SAMEPL");
  gFD1->Draw("SAMEPL");
  gFT0->Draw("SAMEPL");
  gFT1->Draw("SAMEPL");
  gF0->Draw("SAMEPL");
  gF2->Draw("SAMEPL");
  msys2->cd(2);
  as2->Draw();
  gSX->Draw("SAMEP");
  gSXX->Draw("SAMEP");

  return 0;

  TLegend *leg2 = new TLegend(0.5,0.9,0.5,0.9);
  leg2->AddEntry(gA0, "A0");
  leg2->AddEntry(gA1, "A1");
  leg2->AddEntry(gD0, "D0");
  leg2->AddEntry(gD1, "D1");
  leg2->AddEntry(gT0, "T0");
  leg2->AddEntry(gT1, "T1");
  leg2->AddEntry(gFA0, "FA0");
  leg2->AddEntry(gFA1, "FA1");
  leg2->AddEntry(gFD0, "FD0");
  leg2->AddEntry(gFD1, "FD1");
  leg2->AddEntry(gFT0, "FT0");
  leg2->AddEntry(gFT1, "FT1");
  leg2->SetNColumns(3);
  leg2->Draw();

  return 0;
}
//=======================================
void setupMB( TGraph *gg, int col ) {
  gg->SetLineColor( col );
  gg->SetMarkerColor( col );
  gg->SetMarkerStyle( 20 );
  gg->SetFillColor( kWhite );
  gg->RemovePoint(11); // removing
}
//=======================================
void setupERT( TGraph *gg, int col ) {
  gg->SetLineColor( col );
  gg->SetMarkerColor( col );
  gg->SetMarkerStyle( 24 );
  gg->SetFillColor( kWhite );
  for(int fp=0; fp!=9; ++fp) gg->RemovePoint(0); // removing
}
//=======================================
void setupMB( TGraphErrors *gg, int col ) {
  setupMB( (TGraph*) (gg), col );
}
void setupERT( TGraphErrors *gg, int col ) {
  setupERT( (TGraph*) (gg), col );
}
//=======================================
void final() {
  gStyle->SetOptStat(0);
  const int ncuts = 15;
  const int npts = 12;

  float pt[2][npts], ept[2][npts], minpt[2][npts], maxpt[2][npts];
  float v[2][ncuts][npts], ev[2][ncuts][npts];
  float b0[2][ncuts][npts], eb0[2][ncuts][npts];
  float b1[2][ncuts][npts], eb1[2][ncuts][npts];
  float b2[2][ncuts][npts], eb2[2][ncuts][npts];
  float c[2][ncuts][npts];
  TString cut[ncuts] = {"NOM","A0","A1","D0","D1","T0","T1","FA0","FA1","FD0","FD1","FT0","FT1","NOM0","NOM2"};
  float zeroes[npts];
  for(int i=0; i!=npts; ++i) zeroes[i] = 0;

  for(int ic=0; ic!=ncuts; ++ic) {
    ifstream fin( Form("dat/CN_%s_ord1.dat",cut[ic].Data()) );
    for(;;) {
      TString tcut;
      fin >> tcut;
      if(!fin.good()) break;
      float ptmin, ptmax;
      fin >> ptmin >> ptmax;
      //======= FIND THE PT BIN
      int ptb = 0;
      float ptmed = 0.5*(ptmin+ptmax);
      for(int x=0; x!=npts; ++x) {
	if(ptbins[x+1]>ptmed && ptbins[x]<ptmed) {
	  ptb = x;
	  break;
	}
      }
      minpt[0][ptb] = ptmin;
      maxpt[0][ptb] = ptmax;
      pt[0][ptb] = ptmed;
      //======= FIND THE CUT BIN
      /*
      for(int x=0; x!=ncuts; ++x) {
	TString altcut = Form("ERT%s",cut[x].Data() );
	if(cut[x] == tcut) {
	  fin >> v[0][x][ptb] >> ev[0][x][ptb] >> c[0][x][ptb];
	} else if(altcut == tcut) {
	  fin >> v[1][x][ptb] >> ev[1][x][ptb] >> c[1][x][ptb];
	}
      }
      */
      int ertbin = 0;
      if( tcut.Contains("ERT") ) ertbin = 1;
      fin >> v[ertbin][ic][ptb] >> ev[ertbin][ic][ptb];
      fin >> b0[ertbin][ic][ptb] >> eb0[ertbin][ic][ptb];
      fin >> b1[ertbin][ic][ptb] >> eb1[ertbin][ic][ptb];
      fin >> b2[ertbin][ic][ptb] >> eb2[ertbin][ic][ptb];
      fin >> c[ertbin][ic][ptb];
    }
  }

  TH2F *axis = new TH2F("axis",";pT (GeV);C2",100,0,10, 100,0,0.025);
  TH2F *axisb0 = new TH2F("axisb0",";pT (GeV);b0",100,0,10, 100,-0.045,0.045);
  TH2F *axisb1 = new TH2F("axisb1",";pT (GeV);b1",100,0,10, 100,-0.145,0.145);
  TH2F *axisb2 = new TH2F("axisb2",";pT (GeV);b2",100,0,10, 100,-0.325,0.325);
  int color[ncuts] = {kBlack, kRed-3, kGreen-3, kGreen+3,
		      kRed-3, kRed+3, kCyan-3, kCyan+3,
		      kOrange-3, kOrange+3, kMagenta-3, kMagenta+3,
		      kBlue-3, kBlue+3, kGray};
  TLegend *leg = new TLegend(0.4,0.1,0.9,0.4);
  leg->SetNColumns(3);
  TGraphErrors *gr[2][ncuts];
  TGraph *grb0[2][ncuts];
  TGraph *grb1[2][ncuts];
  TGraph *grb2[2][ncuts];
  TCanvas *cvn = new TCanvas("cvn","cvn");
  axis->Draw();
  TCanvas *cb0 = new TCanvas("cb0","cb0");
  axisb0->Draw();
  TCanvas *cb1 = new TCanvas("cb1","cb1");
  axisb1->Draw();
  TCanvas *cb2 = new TCanvas("cb2","cb2");
  axisb2->Draw();
  for(int x=0; x!=ncuts; ++x) {
    gr[0][x] = new TGraphErrors(npts,pt[0],v[0][x],zeroes,ev[0][x]); setupMB(gr[0][x],color[x]);
    grb0[0][x] = new TGraph(npts,pt[0],b0[0][x]); setupMB(grb0[0][x],color[x]);
    grb1[0][x] = new TGraph(npts,pt[0],b1[0][x]); setupMB(grb1[0][x],color[x]);
    grb2[0][x] = new TGraph(npts,pt[0],b2[0][x]); setupMB(grb2[0][x],color[x]);
    gr[1][x] = new TGraphErrors(npts,pt[0],v[1][x],zeroes,ev[1][x]); setupERT(gr[1][x],color[x]);
    grb0[1][x] = new TGraph(npts,pt[0],b0[1][x]); setupERT(grb0[1][x],color[x]);
    grb1[1][x] = new TGraph(npts,pt[0],b1[1][x]); setupERT(grb1[1][x],color[x]);
    grb2[1][x] = new TGraph(npts,pt[0],b2[1][x]); setupERT(grb2[1][x],color[x]);
    leg->AddEntry( gr[0][x], cut[x].Data() );

    cvn->cd(); gr[0][x]->Draw("PLSAME");   gr[1][x]->Draw("PLSAME");
    cb0->cd(); grb0[0][x]->Draw("PLSAME"); grb0[1][x]->Draw("PLSAME");
    cb1->cd(); grb1[0][x]->Draw("PLSAME"); grb1[1][x]->Draw("PLSAME");
    cb2->cd(); grb2[0][x]->Draw("PLSAME"); grb2[1][x]->Draw("PLSAME");
  }
  cvn->cd(); leg->Draw();
  cb0->cd(); leg->Draw();
  cb1->cd(); leg->Draw();
  cb2->cd(); leg->Draw();

  float sysper[2][npts];
  float ptwidth[npts];
  for(int i=0; i!=npts; ++i) {
    if(i<6) { /// 2percent below 2GeV
      sysper[0][i] = 0.02*v[0][0][i];
    } else { /// 4percent above 2GeV
      sysper[0][i] = 0.04*v[0][0][i];
    }
    sysper[1][i] = 0.18*v[1][0][i];
  }
  for(int i=0; i!=nptbins; ++i) {
    ptwidth[i] = (ptbins[i+1]-ptbins[i])/2;
  }
  TGraphErrors *grSys[2];
  grSys[0] = new TGraphErrors(npts,pt[0],v[0][0],ptwidth,sysper[0]);
  grSys[1] = new TGraphErrors(npts,pt[0],v[1][0],ptwidth,sysper[1]);
  setupMB(grSys[0],color[0]);
  setupERT(grSys[1],color[0]);
  grSys[0]->SetFillColor(kGray);
  grSys[1]->SetFillColor(kGray);
  TCanvas *status = new TCanvas ("results","results");
  axis->Draw();
  grSys[0]->Draw("E2SAME");
  grSys[1]->Draw("E2SAME");
  gr[0][0]->Draw("PSAME");
  gr[1][0]->Draw("PSAME");
  leg->Draw();

  //int whichisuniverse[13] = {0,0,2,0,4,0,6,0,8,0,10,0,12};
  int whichisuniverse[13] = {0,1,0,3,0,5,0,7,0,9,0,11,0};
  float systematics[2][13][npts];
  float systematicsRMS[2][npts];
  TGraph *grSystematics[2][13];
  TH2F *axisSys = new TH2F("axisSys",";pT (GeV);Residuals  /  Nominal",100,0,10,100,-1.3,1.3);
  new TCanvas("sys","sys");
  axisSys->Draw();
  TGraph *grSystematicsRMS[2];
  for(int xme=0; xme!=2; ++xme) {
    for(int xp=0; xp!=npts; ++xp) {
      systematicsRMS[xme][xp] = 0;
      for(int xc=1; xc!=13; ++xc) {
	double vi = v[xme][xc][xp];
	double vo = v[xme][0][xp];
	double evi = ev[xme][xc][xp];
	double evo = ev[xme][0][xp];
	double evb = ev[xme][ whichisuniverse[xc] ][xp];
	systematics[xme][xc][xp] = (vi-vo)/vo;
	systematicsRMS[xme][xp] += systematics[xme][xc][xp]*systematics[xme][xc][xp];
      }
      systematicsRMS[xme][xp] = TMath::Sqrt( systematicsRMS[xme][xp]/12.0 );//13-1
    }
    for(int xc=1; xc!=13; ++xc) {
      grSystematics[xme][xc] = new TGraph(npts,pt[0],systematics[xme][xc]);
      grSystematics[xme][xc]->SetFillColor(kWhite);
      grSystematics[xme][xc]->SetLineColor( color[xc] );
      grSystematics[xme][xc]->SetMarkerColor( color[xc] );
      grSystematics[xme][xc]->SetMarkerStyle(20);
      grSystematics[xme][xc]->Draw("PLSAME");
      if(xme==0)
	setupMB(grSystematics[xme][xc],color[xc]);
      else
	setupERT(grSystematics[xme][xc],color[xc]);
    }
    grSystematicsRMS[xme] = new TGraph(npts,pt[0],systematicsRMS[xme]);
  }
  setupMB(grSystematicsRMS[0],kBlack);
  setupERT(grSystematicsRMS[1],kBlack);
  TH2F *axisSysRMS = new TH2F("axisSysRMS",";pT (GeV);RMS of Residuals/NOM",100,0,10,100,0.0,0.2);
  new TCanvas("rms","rms");
  axisSysRMS->Draw();
  grSystematicsRMS[0]->Draw("PSAME");
  grSystematicsRMS[1]->Draw("PSAME");
  //leg->Draw();

  return 0;

  float sum[2][npts];
  float syst[2][npts];
  float sys2[2][npts];
  float res1[2][ncuts][npts];
  float res2[2][ncuts][npts];
  TGraph *grRes1[2][ncuts];
  TGraph *grRes2[2][ncuts];
  for(int xp=0; xp!=npts-1; ++xp) {
    sum[0][xp] = 0;
    for(int xc=1; xc!=ncuts; ++xc) {
      sum[0][xp] += TMath::Abs(v[0][xc][xp]-v[0][0][xp]);
      res1[0][xc][xp] = (v[0][xc][xp]-v[0][0][xp])/ev[0][0][xp];
      res2[0][xc][xp] = (v[0][xc][xp]-v[0][0][xp])/v[0][0][xp];
    }
    syst[0][xp] = TMath::Abs( sum[0][xp]/(ncuts-1)/ev[0][0][xp] );
    sys2[0][xp] = TMath::Abs( sum[0][xp]/(ncuts-1)/v[0][0][xp] );
    //cout << xp << " " << syst[0][xp] << endl;
  }

  TH2F *axis2 = new TH2F("axis2",";pT (GeV);C2",100,0,10,100,0,0.5);
  new TCanvas();
  axis2->Draw();
  TGraph *grS1 = new TGraph(npts-1,pt[0],syst[0]);
  TGraph *grS2 = new TGraph(npts-1,pt[0],sys2[0]);
  grS1->SetFillColor(kWhite);
  grS1->SetLineColor(kBlack);
  grS1->SetMarkerColor(kBlack);
  grS1->SetMarkerStyle(20);
  grS2->SetFillColor(kWhite);
  grS2->SetLineColor(kBlack);
  grS2->SetMarkerColor(kBlack);
  grS2->SetMarkerStyle(24);
  grS1->Draw("PLSAME");
  grS2->Draw("PLSAME");
  TLegend *leg2 = new TLegend(0.5,0.5,0.9,0.9);
  leg2->AddEntry(grS1,"AvgRes / Sigma");
  leg2->AddEntry(grS2,"AvgRes / V2");
  leg2->Draw();

  TH2F *axisS1 = new TH2F("axisS1",";pT (GeV);Diff  /  StatError",100,0,10,100,-3,+3);
  TH2F *axisS2 = new TH2F("axisS2",";pT (GeV);Diff  /  C2",100,0,10,100,-0.1,+0.1);
  TCanvas *csys1 = new TCanvas();
  axisS1->Draw();
  TCanvas *csys2 = new TCanvas();
  axisS2->Draw();
  for(int xc=1; xc!=ncuts; ++xc) {
    grRes1[0][xc] = new TGraph(npts-1,pt[0],res1[0][xc]);
    grRes2[0][xc] = new TGraph(npts-1,pt[0],res2[0][xc]);
    grRes1[0][xc]->SetFillColor(kWhite);
    grRes1[0][xc]->SetLineColor( color[xc] );
    grRes1[0][xc]->SetMarkerColor( color[xc] );
    grRes1[0][xc]->SetMarkerStyle(20);
    grRes2[0][xc]->SetFillColor(kWhite);
    grRes2[0][xc]->SetLineColor( color[xc] );
    grRes2[0][xc]->SetMarkerColor( color[xc] );
    grRes2[0][xc]->SetMarkerStyle(24);
    csys1->cd();
    grRes1[0][xc]->Draw("PLSAME");
    csys2->cd();
    grRes2[0][xc]->Draw("PLSAME");
  }


  return 0;

}

