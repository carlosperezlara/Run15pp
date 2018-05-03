#include "configtime.h"

int step4_sectoroffset(int run=454948) {
  LoadBadTowers();
  loadranges();
  LoadLCs();
  TFile *file = new TFile( Form("out/MyResultsTOF3_%d.root",run) );
  TH2F *hTOF = (TH2F*) file->Get("hTOF");
  int bymax = hTOF->GetYaxis()->GetNbins();
  TH1F *proj = new TH1F("temp","",
			bymax,
			hTOF->GetYaxis()->GetBinLowEdge(1),
			hTOF->GetYaxis()->GetBinLowEdge(bymax+1));
  float xx[8][5000];
  float yy[8][5000];
  int nn[8] = {0,0,0,0,0,0,0,0};
  for(int b=0; b!=hTOF->GetXaxis()->GetNbins(); ++b) {
    int twr = b;
    //cout << "TWR " << twr;
    if(isbad[twr]!=0) continue;
    proj->Reset();
    for(int jj=0; jj!=bymax; ++jj) {
      float cont = hTOF->GetBinContent(b+1,jj+1);
      float cent = hTOF->GetYaxis()->GetBinCenter(jj+1);
      proj->Fill(cent,cont);
    }
    int sec = 0;
    for(int i=0; i!=8; ++i)
      if(twr>secs[i]) sec = i;
    //if(sec!=0) break;
    xx[sec][nn[sec]] = twr;
    int rangeswid[4] = {3,3,3,3};
    float maxb;
    float xmin;
    float xmax;
    //float estimate = (tdcmax[twr]-100)*MyLC[twr] - 115;
    float estimate = -23;
    //cout << "TWR " << twr << " "  << estimate << endl;
    xmin = estimate - 3;
    xmax = estimate + 3;
    proj->GetXaxis()->SetRangeUser(xmin,xmax);
    for(int ittr=0; ittr!=0; ++ittr) {
      maxb = proj->GetMean();
      //cout << " " << maxb;
      xmin = maxb-rangeswid[ittr];
      xmax = maxb+rangeswid[ittr];
      proj->GetXaxis()->SetRangeUser(xmin,xmax);
    }
    yy[sec][nn[sec]] = proj->GetMean();
    nn[sec]++;
    //cout << " " <<proj->GetMean();
    //cout << endl;
  }

  TString fileoutname( Form("offset2_%d.dat",run) );
  ofstream fout( fileoutname.Data() );
  TGraph *gr[8];
  TF1 *fit[8];
  fout << run << " ";
  for(int i=0; i!=8; ++i) {
    cout << "sector " << i << " " << nn[i] << endl;
    //if(nn[i]<10) {
    gr[i] = new TGraph(nn[i],xx[i],yy[i]);
    fit[i] = new TF1( Form("FIT%d",i), "[0]");
    gr[i]->Fit(fit[i]);
    fout << fit[i]->GetParameter(0) << " ";
    //} else {
    //fout << 0.0 << " ";
    //}
  }
  fout << endl;
  fout.close();
  cout << "DATA SAVED IN FILE " << fileoutname.Data() << endl;

  TCanvas *main = new TCanvas();
  main->Divide(4,2);
  main->cd(1); gr[0]->Draw("A*L");
  main->cd(2); gr[1]->Draw("A*L");
  main->cd(3); gr[2]->Draw("A*L");
  main->cd(4); gr[3]->Draw("A*L");
  main->cd(5); gr[4]->Draw("A*L");
  main->cd(6); gr[5]->Draw("A*L");
  main->cd(7); gr[6]->Draw("A*L");
  main->cd(8); gr[7]->Draw("A*L");

  

  return 0;
}
