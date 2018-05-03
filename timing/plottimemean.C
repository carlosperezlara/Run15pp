int secs[9] = { 0, 2592, 5184, 7776, 10368, 12960, 15552, 20160, 24768 };

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

int plottimemean() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile("out/all3.root");
  TProfile *pro[4][8];
  for(int i=0; i!=8; ++i) {
    pro[0][i] = (TProfile*) file->Get(Form("MeanTDC1_SEC%d",i));
    pro[1][i] = (TProfile*) file->Get(Form("MeanTDC2_SEC%d",i));
    pro[2][i] = (TProfile*) file->Get(Form("MeanTDC3_SEC%d",i));
    pro[3][i] = (TProfile*) file->Get(Form("MeanTDC4_SEC%d",i));
    style( pro[0][i], pro[1][i], pro[2][i] );
  }
  TCanvas *main1 = new TCanvas("main1","",1300,600);
  main1->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main1->cd(i+1); pro[0][i]->Draw(); pro[1][i]->Draw("SAME"); pro[2][i]->Draw("SAME");
  }
  
  TH1F *sigmas[3][8];
  int bad[8] = {0,0,0,0,0,0,0,0};
  int nnn[8];
  float xx0[8][5000];
  float yy0[8][5000];
  float zz0[8][5000];
  float ww0[8][5000];
  float vv0[8][5000];
  float del[8][5000];
  float min[8] = {1500,1700,1500,1300,1800,1500, 200, 300};
  float max[8] = {3000,6000,6000,6000,6000,6000,1600,6000};
  ofstream fout("badtowerstof.dat");
  for(int i=0; i!=8; ++i) {
    nnn[i] = 0;
    sigmas[0][i] = new TH1F( Form("sigma0_%d",i),"",1000,0,100 );
    sigmas[1][i] = new TH1F( Form("sigma1_%d",i),"",1000,0,100 );
    sigmas[2][i] = new TH1F( Form("sigma2_%d",i),"",1000,0,100 );
    style( sigmas[0][i], sigmas[1][i], sigmas[2][i] );
    cout << "========= S E C T O R     " << i << " ========" << endl;
    for(int b=0; b!=pro[0][i]->GetNbinsX(); ++b) {
      int twr = b + secs[i];
      xx0[i][nnn[i]] = pro[0][i]->GetXaxis()->GetBinCenter(b+1);
      yy0[i][nnn[i]] = pro[0][i]->GetBinContent(b+1);
      zz0[i][nnn[i]] = pro[1][i]->GetBinContent(b+1);
      ww0[i][nnn[i]] = pro[2][i]->GetBinContent(b+1);
      if(yy0[i][nnn[i]]<min[i]||zz0[i][nnn[i]]<min[i]||
      	 yy0[i][nnn[i]]>max[i]||zz0[i][nnn[i]]>max[i]) {
	fout << twr << " " << 1 << endl;
	bad[i]++;
      	continue;
      }
      del[i][nnn[i]] = ( TMath::Abs(ww0[i][nnn[i]]-yy0[i][nnn[i]]) +
			 TMath::Abs(zz0[i][nnn[i]]-ww0[i][nnn[i]]) +
			 TMath::Abs(zz0[i][nnn[i]]-yy0[i][nnn[i]]) );
      if(del[i][nnn[i]]>500) { 
	fout << twr << " " << 2 << endl;
	bad[i]++;
	continue;
      }
      vv0[i][nnn[i]] = pro[3][i]->GetBinContent(b+1);
      if(vv0[i][nnn[i]]<min[i]||vv0[i][nnn[i]]>max[i]) {
	fout << twr << " " << 1 << endl;
	bad[i]++;
      	continue;
      }
      
      sigmas[0][i]->Fill( pro[0][i]->GetBinError(b+1) );
      sigmas[1][i]->Fill( pro[1][i]->GetBinError(b+1) );
      sigmas[2][i]->Fill( pro[2][i]->GetBinError(b+1) );
      nnn[i]++;
    }
    cout << "BAD TOWERS " << bad[i] << endl;
  }
  fout.close();
  TGraph *gr[4][8];
  TGraph *de[8];
  for(int i=0; i!=8; ++i) {
    gr[0][i] = new TGraph(nnn[i],xx0[i],yy0[i]);
    gr[1][i] = new TGraph(nnn[i],xx0[i],zz0[i]);
    gr[2][i] = new TGraph(nnn[i],xx0[i],ww0[i]);
    gr[3][i] = new TGraph(nnn[i],xx0[i],vv0[i]);
    de[i] = new TGraph(nnn[i],xx0[i],del[i]);
  }
  TCanvas *main2 = new TCanvas("main2","",1300,600);
  main2->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main2->cd(i+1); gr[0][i]->Draw("AP"); gr[1][i]->Draw("PSAME"); gr[2][i]->Draw("PSAME");
    style( gr[0][i], gr[1][i], gr[2][i] );
    gr[0][i]->GetYaxis()->SetRangeUser(-100,3000);
  }

  TCanvas *main3 = new TCanvas("main3","",1300,600);
  main3->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main3->cd(i+1)->SetLogy(1); sigmas[0][i]->Draw(); sigmas[1][i]->Draw("SAME"); sigmas[2][i]->Draw("SAME");
  }

  TCanvas *main4 = new TCanvas("main4","",1300,600);
  main4->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main4->cd(i+1); de[i]->Draw("A*L");
    de[i]->GetYaxis()->SetRangeUser(0,1000);
  }

  TCanvas *main5 = new TCanvas("main5","",1300,600);
  main5->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    main5->cd(i+1); gr[3][i]->Draw("A*L");
  }

  return 0;
}
