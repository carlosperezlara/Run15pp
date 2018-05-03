int step1_timemean(int run=454948) {
  TFile *file = new TFile( Form("out/MyResultsTOF1_%d.root",run) );
  TProfile *sec[8];
  sec[0] = (TProfile*) file->Get("MeanTDC_SEC0");
  sec[1] = (TProfile*) file->Get("MeanTDC_SEC1");
  sec[2] = (TProfile*) file->Get("MeanTDC_SEC2");
  sec[3] = (TProfile*) file->Get("MeanTDC_SEC3");
  sec[4] = (TProfile*) file->Get("MeanTDC_SEC4");
  sec[5] = (TProfile*) file->Get("MeanTDC_SEC5");
  sec[6] = (TProfile*) file->Get("MeanTDC_SEC6");
  sec[7] = (TProfile*) file->Get("MeanTDC_SEC7");

  TGraph *gr[8];
  float xx[8][5000];
  float yy[8][5000];
  int nn[8];
  for(int i=0; i!=8; ++i) {
    nn[i] = 0;
    for(int b=0; b!=sec[i]->GetXaxis()->GetNbins(); ++b) {
      double bins = sec[i]->GetBinEntries(b+1);
      if(bins<2) continue;
      yy[i][nn[i]] = sec[i]->GetBinContent(b+1);
      xx[i][nn[i]] = sec[i]->GetXaxis()->GetBinCenter(b+1);
      nn[i]++;
    }
    gr[i] = new TGraph( nn[i], xx[i], yy[i] );
    gr[i]->SetMarkerStyle(4);
    gr[i]->SetMarkerColor(kOrange-3);
  }

  TF1 *fit[8];
  for(int i=0; i!=8; ++i) {
    fit[i] = new TF1( Form("fsec%d",i), "[0]",2000,4000);
    cout << "==========================================" << endl;
    cout << "============ S E C T O R   " << i << " =============" << endl;
    cout << "==========================================" << endl;
    //sec[i]->Fit( fit[i], "LR","",2000,4000);
    //sec[i]->Fit( fit[i], "MEIRWW","",2000,4000);
    //sec[i]->Fit( fit[i], "MEIRWW","",2000,4000);
    gr[i]->Fit( fit[i], "LR","",2000,4000);
    gr[i]->Fit( fit[i], "MEIRWW","",2000,4000);
    gr[i]->Fit( fit[i], "MEIRWW","",2000,4000);
    sec[i]->SetFillColor(kYellow-8);
    sec[i]->SetLineColor(kYellow-8);
    sec[i]->SetMarkerStyle(20);
    sec[i]->SetMarkerColor(kBlue-3);
  }

  ofstream fout( Form("timemean_%d.dat",run) );
  fout << run << " ";
  for(int i=0; i!=8 ; ++i) {
    fout << Form("%.3f ", fit[i]->GetParameter(0));
  }
  fout << endl;
  fout.close();


  TCanvas *main = new TCanvas();
  main->Divide(4,2,0,0);
  main->cd(1); sec[0]->Draw(); gr[0]->Draw("PSAME");
  main->cd(2); sec[1]->Draw(); gr[1]->Draw("PSAME");
  main->cd(3); sec[2]->Draw(); gr[2]->Draw("PSAME");
  main->cd(4); sec[3]->Draw(); gr[3]->Draw("PSAME");
  main->cd(5); sec[4]->Draw(); gr[4]->Draw("PSAME");
  main->cd(6); sec[5]->Draw(); gr[5]->Draw("PSAME");
  main->cd(7); sec[6]->Draw(); gr[6]->Draw("PSAME");
  main->cd(8); sec[7]->Draw(); gr[7]->Draw("PSAME");

  return 0;
}
