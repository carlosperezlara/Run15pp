void style(TGraph *gG, int col, TString yt, TString xt) {
  gG->SetMarkerStyle(1);
  gG->SetMarkerColor( col );
  gG->SetLineColor( col );
  gG->GetXaxis()->SetTitle( xt.Data() );
  gG->GetYaxis()->SetTitle( yt.Data() );
  gG->GetXaxis()->SetLabelSize(0.08);
  gG->GetYaxis()->SetLabelSize(0.08);
  gG->GetXaxis()->SetTitleSize(0.08);
  gG->GetYaxis()->SetTitleSize(0.08);
  gG->GetYaxis()->SetTitleOffset(0.45);
  gG->SetTitle("");
  gG->GetXaxis()->SetRangeUser(0,3510);
  gG->GetXaxis()->CenterTitle();
  gG->GetYaxis()->CenterTitle();
  gG->GetYaxis()->SetNdivisions(509);
}

int plotFVTXcalibration() {
  gStyle->SetOptStat(0);

  float xx[3500];
  float bbqx[1][3500];
  float bbqy[1][3500];
  float tmp;
  int ns = 0;
  TString segment;
  ifstream fseg("segments.dat");
  ifstream ftmp;
  for(;;++ns) {
    xx[ns] = ns;
    fseg >> segment;
    if(!fseg.good()) break;
    ftmp.open( Form( "out/calib/fv%s.dat1", segment.Data() ) );
    ftmp >> tmp;
    ftmp >> bbqx[0][ns] >> bbqy[0][ns] >> tmp >> tmp;
    ftmp.close();
  }
  fseg.close();

  TGraphErrors *gBBAQX = new TGraphErrors(ns,xx,bbqx[0]);
  style( gBBAQX, kRed-3,  "< QX >", "Segment  index");

  TGraphErrors *gBBAQY = new TGraphErrors(ns,xx,bbqy[0]);
  style( gBBAQY, kRed-3,  "< QY >", "Segment  index");

  TCanvas *main = new TCanvas();
  main->SetLeftMargin( 0.09 );
  main->SetRightMargin( 0.05 );
  main->SetBottomMargin( 0.2 );
  main->Divide(1,2,0,0);

  main->cd(2);
  gBBAQX->Draw("APL");
  gBBAQX->GetYaxis()->SetRangeUser(+1.5,+9.5);

  main->cd(1);
  gBBAQY->Draw("APL");
  gBBAQY->GetYaxis()->SetRangeUser(-2.5,+2.5);
  
  return 0;
}
