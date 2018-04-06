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

int plotBBCcalibration() {
  gStyle->SetOptStat(0);

  float xx[3500];
  float bbqx[3][3500];
  float bbqy[3][3500];
  float tmp;
  int ns = 0;
  TString segment;
  ifstream fseg("segments.dat");
  ifstream ftmp;
  for(;;++ns) {
    xx[ns] = ns;
    fseg >> segment;
    if(!fseg.good()) break;
    ftmp.open( Form( "out/calib/bb%s.dat1", segment.Data() ) );
    ftmp >> tmp;
    ftmp >> bbqx[0][ns] >> bbqy[0][ns] >> tmp >> tmp;
    ftmp >> bbqx[1][ns] >> bbqy[1][ns] >> tmp >> tmp;
    ftmp >> bbqx[2][ns] >> bbqy[2][ns] >> tmp >> tmp;
    ftmp.close();
  }
  fseg.close();

  TGraphErrors *gBBAQX = new TGraphErrors(ns,xx,bbqx[0]);
  TGraphErrors *gBBBQX = new TGraphErrors(ns,xx,bbqx[1]);
  TGraphErrors *gBBCQX = new TGraphErrors(ns,xx,bbqx[2]);
  style( gBBAQX, kRed-3,  "< QX >", "Segment  index");
  style( gBBBQX, kBlue-3, "< QX >", "Segment  index");
  style( gBBCQX, kBlack,  "< QX >", "Segment  index");

  TGraphErrors *gBBAQY = new TGraphErrors(ns,xx,bbqy[0]);
  TGraphErrors *gBBBQY = new TGraphErrors(ns,xx,bbqy[1]);
  TGraphErrors *gBBCQY = new TGraphErrors(ns,xx,bbqy[2]);
  style( gBBAQY, kRed-3,  "< QY >", "Segment  index");
  style( gBBBQY, kBlue-3, "< QY >", "Segment  index");
  style( gBBCQY, kBlack,  "< QY >", "Segment  index");

  TCanvas *main = new TCanvas();
  main->SetLeftMargin( 0.09 );
  main->SetRightMargin( 0.05 );
  main->SetBottomMargin( 0.2 );
  main->Divide(1,2,0,0);

  main->cd(2);
  gBBAQX->Draw("APL");
  gBBBQX->Draw("PLSAME");
  gBBCQX->Draw("PLSAME");
  gBBAQX->GetYaxis()->SetRangeUser(-6.5,+3.5);

  main->cd(1);
  gBBAQY->Draw("APL");
  gBBBQY->Draw("PLSAME");
  gBBCQY->Draw("PLSAME");
  gBBAQY->GetYaxis()->SetRangeUser(-2.5,+2.5);
  
  return 0;
}
