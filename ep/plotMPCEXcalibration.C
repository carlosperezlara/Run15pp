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
  gG->SetFillColor( kWhite );
}

int plotMPCEXcalibration() {
  gStyle->SetOptStat(0);

  float xx[3500];
  float bbqx[9][3500];
  float bbqy[9][3500];
  float tmp;
  int ns = 0;
  TString segment;
  ifstream fseg("segments.dat");
  ifstream ftmp;
  for(;;++ns) {
    xx[ns] = ns;
    fseg >> segment;
    if(!fseg.good()) break;
    ftmp.open( Form( "out/calib/ex%s.dat1", segment.Data() ) );
    ftmp >> tmp;
    ftmp >> bbqx[0][ns] >> bbqy[0][ns] >> tmp >> tmp;
    ftmp >> bbqx[1][ns] >> bbqy[1][ns] >> tmp >> tmp;
    ftmp >> bbqx[2][ns] >> bbqy[2][ns] >> tmp >> tmp;
    ftmp >> bbqx[3][ns] >> bbqy[3][ns] >> tmp >> tmp;
    ftmp >> bbqx[4][ns] >> bbqy[4][ns] >> tmp >> tmp;
    ftmp >> bbqx[5][ns] >> bbqy[5][ns] >> tmp >> tmp;
    ftmp >> bbqx[6][ns] >> bbqy[6][ns] >> tmp >> tmp;
    ftmp >> bbqx[7][ns] >> bbqy[7][ns] >> tmp >> tmp;
    ftmp >> bbqx[8][ns] >> bbqy[8][ns] >> tmp >> tmp;
    ftmp.close();
  }
  fseg.close();

  TGraphErrors *gBBAQX = new TGraphErrors(ns,xx,bbqx[0]);
  TGraphErrors *gBBBQX = new TGraphErrors(ns,xx,bbqx[1]);
  TGraphErrors *gBBCQX = new TGraphErrors(ns,xx,bbqx[2]);
  TGraphErrors *gBBDQX = new TGraphErrors(ns,xx,bbqx[3]);
  TGraphErrors *gBBEQX = new TGraphErrors(ns,xx,bbqx[4]);
  TGraphErrors *gBBFQX = new TGraphErrors(ns,xx,bbqx[5]);
  TGraphErrors *gBBGQX = new TGraphErrors(ns,xx,bbqx[6]);
  TGraphErrors *gBBHQX = new TGraphErrors(ns,xx,bbqx[7]);
  TGraphErrors *gBBIQX = new TGraphErrors(ns,xx,bbqx[8]);
  style( gBBAQX, kRed-3,    "< QX >", "Segment  index");
  style( gBBBQX, kOrange-3, "< QX >", "Segment  index");
  style( gBBCQX, kCyan-3,   "< QX >", "Segment  index");
  style( gBBDQX, kBlue-3,   "< QX >", "Segment  index");
  style( gBBEQX, kRed-3,    "< QX >", "Segment  index");
  style( gBBFQX, kOrange-3, "< QX >", "Segment  index");
  style( gBBGQX, kCyan-3,   "< QX >", "Segment  index");
  style( gBBHQX, kBlue-3,   "< QX >", "Segment  index");
  style( gBBIQX, kBlack,    "< QX >", "Segment  index");

  TGraphErrors *gBBAQY = new TGraphErrors(ns,xx,bbqy[0]);
  TGraphErrors *gBBBQY = new TGraphErrors(ns,xx,bbqy[1]);
  TGraphErrors *gBBCQY = new TGraphErrors(ns,xx,bbqy[2]);
  TGraphErrors *gBBDQY = new TGraphErrors(ns,xx,bbqy[3]);
  TGraphErrors *gBBEQY = new TGraphErrors(ns,xx,bbqy[4]);
  TGraphErrors *gBBFQY = new TGraphErrors(ns,xx,bbqy[5]);
  TGraphErrors *gBBGQY = new TGraphErrors(ns,xx,bbqy[6]);
  TGraphErrors *gBBHQY = new TGraphErrors(ns,xx,bbqy[7]);
  TGraphErrors *gBBIQY = new TGraphErrors(ns,xx,bbqy[8]);
  style( gBBAQY, kRed-3,    "< QY >", "Segment  index");
  style( gBBBQY, kOrange-3, "< QY >", "Segment  index");
  style( gBBCQY, kCyan-3,   "< QY >", "Segment  index");
  style( gBBDQY, kBlue-3,   "< QY >", "Segment  index");
  style( gBBEQY, kRed-3,    "< QY >", "Segment  index");
  style( gBBFQY, kOrange-3, "< QY >", "Segment  index");
  style( gBBGQY, kCyan-3,   "< QY >", "Segment  index");
  style( gBBHQY, kBlue-3,   "< QY >", "Segment  index");
  style( gBBIQY, kBlack,    "< QY >", "Segment  index");


  TCanvas *main = new TCanvas();
  main->SetLeftMargin( 0.09 );
  main->SetRightMargin( 0.05 );
  main->SetBottomMargin( 0.2 );
  main->Divide(2,2,0,0);

  main->cd(1);
  gBBAQY->Draw("APL");
  gBBBQY->Draw("PLSAME");
  gBBCQY->Draw("PLSAME");
  gBBDQY->Draw("PLSAME");
  gBBAQY->GetYaxis()->SetRangeUser(-2.5,+6.5);

  main->cd(2);
  gBBEQY->Draw("APL");
  gBBFQY->Draw("PLSAME");
  gBBGQY->Draw("PLSAME");
  gBBHQY->Draw("PLSAME");
  gBBIQY->Draw("PLSAME");
  gBBEQY->GetYaxis()->SetRangeUser(-2.5,+6.5);

  main->cd(3);
  gBBAQX->Draw("APL");
  gBBBQX->Draw("PLSAME");
  gBBCQX->Draw("PLSAME");
  gBBDQX->Draw("PLSAME");
  gBBAQX->GetYaxis()->SetRangeUser(-8.5,+2.5);

  main->cd(4);
  gBBEQX->Draw("APL");
  gBBFQX->Draw("PLSAME");
  gBBGQX->Draw("PLSAME");
  gBBHQX->Draw("PLSAME");
  gBBIQX->Draw("PLSAME");
  gBBEQX->GetYaxis()->SetRangeUser(-7.5,+1.5);

  TLegend *leg = new TLegend(0.1,0.59,0.9,0.99);
  leg->AddEntry( gBBAQX, "Station 0" );
  leg->AddEntry( gBBBQX, "Station 1" );
  leg->AddEntry( gBBCQX, "Station 2" );
  leg->AddEntry( gBBDQX, "Station 3" );
  leg->AddEntry( gBBIQX, "Full" );
  leg->SetFillColor( kWhite );
  leg->SetNColumns(2);
  main->cd(2);
  leg->Draw();
  
  return 0;
}
