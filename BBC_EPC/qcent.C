int qcent(int run=424356) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TH2F *QH;
  TH1D *h;
  ofstream fout( Form("tables/BBC_%d.dat",run) );
  char xy[2] = {'x','y'};
  for(int ord=0; ord!=6; ++ord) {
    for(int ix=0; ix!=2; ++ix) {
      for(int se=0; se!=2; ++se) {
	//printing 60rows x 40 columns table
	for(int i=0; i!=1; ++i) {
	  QH = (TH2F*) file->Get( Form("BBCQ%d%c_S%d_CB%02d",ord+1,xy[ix],se,i) );
	  int nbins = QH->GetXaxis()->GetNbins(); // NVtx
	  for(int j=0; j!=nbins; ++j) {
	    h = QH->ProjectionY( Form("%s_P%d",QH->GetName(),j), j+1, j+1 );
	    double mean = h->GetMean();
	    if(h->GetEntries()<100) mean = 0;
	    fout << Form(" %e", mean);
	  }
	  fout << endl;
	}
	fout << endl;
      }
    }
  }
  return 0;
}
