int multi(TString sfile="423547_1") {
  TString sin  = "out/out_" + sfile + ".root";
  TString sout = "table.dat";
  TFile *fin = new TFile( sin.Data() );
  //fin->ls();
  TH1F *bm0 = (TH1F*) fin->Get("hBM0");
  TH1F *bm1 = (TH1F*) fin->Get("hBM1");
  TH1F *tm0 = (TH1F*) fin->Get("hTM0");
  TH1F *tm1 = (TH1F*) fin->Get("hTM1");
  ofstream fout( sout.Data(), std::ofstream::out | std::ofstream::app );
  fout << sfile.Data();
  fout << " " << bm0->GetMean() << " " << bm0->GetMeanError();
  fout << " " << bm1->GetMean() << " " << bm1->GetMeanError();
  fout << " " << tm0->GetMean() << " " << tm0->GetMeanError();
  fout << " " << tm1->GetMean() << " " << tm1->GetMeanError();
  fout << endl;
  fout.close();
}
