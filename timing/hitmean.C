int hitmean(int run=454948) {
  TFile *file = new TFile( Form("out/MyResults_%d.root",run) );
  TH1F *sec[8];
  sec[0] = (TH1F*) file->Get("hHitTwrSec0");
  sec[1] = (TH1F*) file->Get("hHitTwrSec1");
  sec[2] = (TH1F*) file->Get("hHitTwrSec2");
  sec[3] = (TH1F*) file->Get("hHitTwrSec3");
  sec[4] = (TH1F*) file->Get("hHitTwrSec4");
  sec[5] = (TH1F*) file->Get("hHitTwrSec5");
  sec[6] = (TH1F*) file->Get("hHitTwrSec6");
  sec[7] = (TH1F*) file->Get("hHitTwrSec7");

  TF1 *fit[8];
  for(int i=0; i!=8; ++i) {
    fit[i] = new TF1( Form("fsec%d",i), "[0]*TMath::Gaus(x,[1],[2],1)" );
    fit[i]->SetParameter(0,sec[i]->GetEntries());
    fit[i]->SetParameter(1,sec[i]->GetMean());
    fit[i]->SetParameter(2,sec[i]->GetRMS());
    fit[i]->SetParLimits(1,sec[i]->GetMean()-sec[i]->GetRMS(),sec[i]->GetMean()+sec[i]->GetRMS());
    fit[i]->SetParLimits(2,0.2*sec[i]->GetRMS(),2*sec[i]->GetRMS());
    cout << "==========================================" << endl;
    cout << "============ S E C T O R   " << i << " =============" << endl;
    cout << "==========================================" << endl;
    sec[i]->Fit( fit[i], "L");
    sec[i]->Fit( fit[i], "MEI");
    sec[i]->Fit( fit[i], "MEI");
    sec[i]->SetFillColor(kYellow-8);
    sec[i]->SetLineColor(kYellow-8);

  }

  ofstream fout( Form("hitmean_%d.dat",run) );
  fout << run << " ";
  for(int i=0; i!=8 ; ++i) {
    fout << Form("%.3f ",sec[i]->GetMean());
  }
  for(int i=0; i!=8 ; ++i) {
    fout << Form("%.3f ", fit[i]->GetParameter(1));
    fout << Form("%.3f ", fit[i]->GetParameter(2));
  }
  fout << endl;
  fout.close();


  TCanvas *main = new TCanvas();
  main->Divide(4,2,0,0);
  main->cd(1); sec[0]->Draw();
  main->cd(2); sec[1]->Draw();
  main->cd(3); sec[2]->Draw();
  main->cd(4); sec[3]->Draw();
  main->cd(5); sec[4]->Draw();
  main->cd(6); sec[5]->Draw();
  main->cd(7); sec[6]->Draw();
  main->cd(8); sec[7]->Draw();

  return 0;
}
