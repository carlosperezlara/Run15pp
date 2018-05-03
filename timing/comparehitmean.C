int comparehitmean(int run1=454774,int run2=455639) {
  TFile *file1 = new TFile( Form("out/MyResults_%d.root",run1) );
  TFile *file2 = new TFile( Form("out/MyResults_%d.root",run2) );

  TH1F *sec0 = (TH1F*) file1->Get("hHitTwrSec0");
  TH1F *sec1 = (TH1F*) file1->Get("hHitTwrSec1");
  TH1F *sec2 = (TH1F*) file1->Get("hHitTwrSec2");
  TH1F *sec3 = (TH1F*) file1->Get("hHitTwrSec3");
  TH1F *sec4 = (TH1F*) file1->Get("hHitTwrSec4");
  TH1F *sec5 = (TH1F*) file1->Get("hHitTwrSec5");
  TH1F *sec6 = (TH1F*) file1->Get("hHitTwrSec6");
  TH1F *sec7 = (TH1F*) file1->Get("hHitTwrSec7");

  TH1F *sec01 = (TH1F*) file2->Get("hHitTwrSec0");
  TH1F *sec11 = (TH1F*) file2->Get("hHitTwrSec1");
  TH1F *sec21 = (TH1F*) file2->Get("hHitTwrSec2");
  TH1F *sec31 = (TH1F*) file2->Get("hHitTwrSec3");
  TH1F *sec41 = (TH1F*) file2->Get("hHitTwrSec4");
  TH1F *sec51 = (TH1F*) file2->Get("hHitTwrSec5");
  TH1F *sec61 = (TH1F*) file2->Get("hHitTwrSec6");
  TH1F *sec71 = (TH1F*) file2->Get("hHitTwrSec7");

  sec0->SetLineColor(kRed-3); sec01->SetLineColor(kBlue-3);
  sec1->SetLineColor(kRed-3); sec11->SetLineColor(kBlue-3);
  sec2->SetLineColor(kRed-3); sec21->SetLineColor(kBlue-3);
  sec3->SetLineColor(kRed-3); sec31->SetLineColor(kBlue-3);
  sec4->SetLineColor(kRed-3); sec41->SetLineColor(kBlue-3);
  sec5->SetLineColor(kRed-3); sec51->SetLineColor(kBlue-3);
  sec6->SetLineColor(kRed-3); sec61->SetLineColor(kBlue-3);
  sec7->SetLineColor(kRed-3); sec71->SetLineColor(kBlue-3);

  TCanvas *main = new TCanvas();
  main->Divide(4,2);
  main->cd(1); sec0->Draw(); sec01->Draw("SAME");
  main->cd(2); sec1->Draw(); sec11->Draw("SAME");
  main->cd(3); sec2->Draw(); sec21->Draw("SAME");
  main->cd(4); sec3->Draw(); sec31->Draw("SAME");
  main->cd(5); sec4->Draw(); sec41->Draw("SAME");
  main->cd(6); sec5->Draw(); sec51->Draw("SAME");
  main->cd(7); sec6->Draw(); sec61->Draw("SAME");
  main->cd(8); sec7->Draw(); sec71->Draw("SAME");

  return 0;
}
