int coef(int run=454777) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TProfile2D *QH;
  float coe;
  ofstream fout( Form("tables/BBC_A_%d.dat",run) );
  for(int ord=0; ord!=4; ++ord) {
    for(int i=0; i!=60; ++i) {
      QH = (TProfile2D*) file->Get( Form("BBCPsiC_Ord%d_Cen%02d",ord,i) );
      //yes, they are reversed!
      int nbins = QH->GetXaxis()->GetNbins(); // vtx
      //cout << " VTX BINS " << nbins << endl;
      // ttable of 32rows x 40columns
      for(int in=0; in!=32; ++in) {
	for(int j=0; j!=nbins; ++j) {
	  //if( QHH->GetBinEntries( j+1 ) < 30 ) coe=0.0;
	  //else
	  coe = QH->GetBinContent( j+1, in+1 );
	  if( TMath::IsNaN( coe ) ) {
	    cout << "ERROR IN " << QH->GetName() << " bins " << j+1 << " " << in+1 << endl;
	  }
	  fout << Form(" %.2f", coe*1e+3);
	}
	fout << endl;
      }
      fout << endl;
      //==
      QH = (TProfile2D*) file->Get( Form("BBCPsiS_Ord%d_Cen%02d",ord,i) );
      //yes, they are reversed!
      int nbins = QH->GetXaxis()->GetNbins(); // vtx
      // ttable of 32rows x 40columns
      for(int in=0; in!=32; ++in) {
	for(int j=0; j!=nbins; ++j) {
	  //if( QHH->GetBinEntries( j+1 ) < 30 ) coe=0.0;
	  //else
	  coe = QH->GetBinContent( j+1, in+1 );
	  if( TMath::IsNaN( coe ) ) {
	    cout << "ERROR IN " << QH->GetName() << " bins " << j+1 << " " << in+1 << endl;
	  }
	  fout << Form(" %.2f", coe*1e+3);
	}
	fout << endl;
      }
      fout << endl;
    }
  }
  return 0;
}
