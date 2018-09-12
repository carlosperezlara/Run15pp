void res(int run=454777) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TF1 *fit = new TF1( "fit", "[0]+[1]*(x-20)", 10, 30 );
  ofstream fout( Form("tables/BBC_R_%d.dat",run) );
  for(int ord=0; ord!=4; ++ord) {
    for(int i=0; i!=1; ++i) {
      TProfile *QH = (TProfile*) file->Get( Form("BBCRes_Ord%d_Cen%02d",ord,i) );
      QH->Fit( fit, "RL", "", 10.5, 29.5 );
      fout << fit->GetParameter( 0 ) << " " << fit->GetParError( 0 ) << " ";
      fout << fit->GetParameter( 1 ) << " " << fit->GetParError( 1 ) << endl;
    }
    fout << endl;
  }
  fout.close();
  file->Close();
}
