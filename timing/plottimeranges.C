int secs[9] = { 0, 2592, 5184, 7776, 10368, 12960, 15552, 20160, 24767 };

void style(TH1F *gr, TH1F *gr2, TH1F *gr3) {
  gr->GetXaxis()->SetNdivisions(504);
  gr->GetXaxis()->SetLabelSize(0.07);
  gr->SetLineColor( kRed-3 );
  gr->SetMarkerStyle( 20 );
  gr2->SetLineColor( kBlue-3 );
  gr2->SetMarkerStyle( 24 );
  gr3->SetLineColor( kGreen-3 );
  gr3->SetMarkerStyle( 24 );
}

void style(TProfile2D *gr) {
  gr->GetXaxis()->SetNdivisions(504);
  gr->GetXaxis()->SetLabelSize(0.07);
  gr->SetMarkerColor( kRed-3 );
  gr->SetMarkerStyle( 20 );
}

void style(TProfile *gr, TProfile *gr2, TProfile *gr3) {
  gr->GetXaxis()->SetNdivisions(504);
  gr->GetXaxis()->SetLabelSize(0.07);
  gr->SetMarkerColor( kRed-3 );
  gr->SetMarkerStyle( 20 );
  gr2->SetMarkerColor( kBlue-3 );
  gr2->SetMarkerStyle( 24 );
  gr3->SetMarkerColor( kGreen-3 );
  gr3->SetMarkerStyle( 24 );
}

void style(TGraph *gr, TGraph *gr2, TGraph *gr3) {
  gr->GetXaxis()->SetNdivisions(504);
  gr->GetXaxis()->SetLabelSize(0.07);
  gr->SetMarkerColor( kRed-3 );
  gr->SetMarkerStyle( 20 );
  gr2->SetMarkerColor( kBlue-3 );
  gr2->SetMarkerStyle( 24 );
  gr3->SetMarkerColor( kGreen-3 );
  gr3->SetMarkerStyle( 24 );
}

int plottimeranges() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile("out/all2.root");
  TProfile2D *pro[8];
  for(int i=0; i!=8; ++i) {
    pro[i] = (TProfile2D*) file->Get(Form("MeanTDC_SEC%d",i));
    style( pro[i] );
  }
  TCanvas *main1 = new TCanvas("main1","",1300,600);
  main1->Divide(4,2);
  int nn;
  float xx[5000];
  float yy[5000];
  float ex[5000];
  float ey[5000];
  for(int i=0; i!=8; ++i) {
    main1->cd(i+1); pro[i]->Draw("colz");
    /*
    for(int b=0; b!=pro[i]->GetNbinsY(); ++b) {
      nn = 0;
      for(int a=0; a!=pro[i]->GetXaxis()->GetNbins(); ++a) {
	if( pro[i]->GetBinError( a+1, b+1 ) < 0.01) continue;
	yy[nn] = pro[i]->GetBinContent( a+1, b+1 );
	ey[nn] = pro[i]->GetBinError( a+1, b+1 );
	ex[nn] = 0.5;
	xx[nn] = pro[i]->GetXaxis()->GetBinCenter( a+1 );
	++nn;
      }
      if(nn<100) continue;
      TGraphErrors *gr = new TGraphErrors(nn,xx,yy,ex,ey);
      TF1 *fitf = new TF1("fitf","[0]+[1]/x");
      gr->Draw("A*");
      gr->Fit(fitf,"MERWWQ","",0,3000);
      cout << "TWR " << b+1 << " A " << fitf->GetParameter(0);
      cout << " B " << fitf->GetParameter(1) << endl;
      //delete fitf;
      //delete gr;
    }
  */
  }
  
  
}
    
