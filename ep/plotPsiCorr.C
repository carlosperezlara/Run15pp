int plotPsiCorr(TString run) {
  gStyle->SetOptStat(0);

  TString ooo = Form("out/%s.root",run.Data());
  std::cout << "READING... " << ooo.Data() << std::endl;
  TFile* f2 = new TFile( ooo.Data() ); 

  //======= HISTOS
  TH2F *hEPC[3][33];
  TString title33[33] = {"#Psi_{2}^{BBA}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{BBB}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{BBA}::#Psi_{2}^{BB}",
                         "#Psi_{2}^{FV}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{FV}::#Psi_{2}^{CA}",
                         "#Psi_{2}^{EX0}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX1}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX2}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX3}::#Psi_{2}^{BB}",
                         "#Psi_{2}^{EX4}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX5}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX6}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX7}::#Psi_{2}^{BB}",
			 "#Psi_{2}^{EX8}::#Psi_{2}^{BB}",
                         "#Psi_{2}^{EX0}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX1}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX2}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX3}::#Psi_{2}^{FV}",
                         "#Psi_{2}^{EX4}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX5}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX6}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX7}::#Psi_{2}^{FV}",
			 "#Psi_{2}^{EX8}::#Psi_{2}^{FV}",
                         "#Psi_{2}^{EX0}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX1}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX2}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX3}::#Psi_{2}^{CA}",
                         "#Psi_{2}^{EX4}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX5}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX6}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX7}::#Psi_{2}^{CA}",
			 "#Psi_{2}^{EX8}::#Psi_{2}^{CA}",
                         "#Psi_{2}^{CA}::#Psi_{2}^{BB}"};
  for(int i=0; i!=3; ++i) for(int j=0; j!=33; ++j) {
      hEPC[i][j] = (TH2F*) f2->Get( Form("EVC%d_DET%d",i,j) );
      hEPC[i][j]->GetXaxis()->SetLabelSize(0.12);
      hEPC[i][j]->GetYaxis()->SetLabelSize(0.12);
    }
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.2);

  TCanvas *cepc0 = new TCanvas();
  cepc0->Divide( 3, 3, 0, 0 );
  cepc0->cd(1);
  for(int i=0; i!=9; ++i) {
    cepc0->cd(i+1);
    hEPC[0][i+5]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[i+5].Data() );
  }

  TCanvas *cepc1 = new TCanvas();
  cepc1->Divide( 3, 3, 0, 0 );
  cepc1->cd(1);
  for(int i=0; i!=9; ++i) {
    cepc1->cd(i+1);
    hEPC[1][i+5]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[i+5].Data() );
  }

  TCanvas *cepc2 = new TCanvas();
  cepc2->Divide( 3, 3, 0, 0 );
  cepc2->cd(1);
  for(int i=0; i!=9; ++i) {
    cepc2->cd(i+1);
    hEPC[2][i+5]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[i+5].Data() );
  }

  int idx[7] = {2,3,32,13,22,31,4};
  TCanvas *cepc4 = new TCanvas("ccc","",1200,600);
  cepc4->Divide( 3, 2, 0, 0 );
  cepc4->cd(1);
  for(int i=0; i!=6; ++i) {
    cepc4->cd(i+1);
    hEPC[0][ idx[i] ]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[ idx[i] ].Data() );
  }

  TCanvas *cepc5 = new TCanvas("ccc","",1200,600);
  cepc5->Divide( 3, 2, 0, 0 );
  cepc5->cd(1);
  for(int i=0; i!=6; ++i) {
    cepc5->cd(i+1);
    hEPC[1][ idx[i] ]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[ idx[i] ].Data() );
  }

  TCanvas *cepc6 = new TCanvas("ccc","",1200,600);
  cepc6->Divide( 3, 2, 0, 0 );
  cepc6->cd(1);
  for(int i=0; i!=6; ++i) {
    cepc6->cd(i+1);
    hEPC[2][ idx[i] ]->Draw("colz");
    tex->DrawLatex( 0.4, 4.5, title33[ idx[i] ].Data() );
  }

  return 0;
}
