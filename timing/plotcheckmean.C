#include "configtime.h"

int plotcheckmean() {
  ifstream fin( "checkmean.dat" );

  float x[8][10000];
  float y[8][10000];
  int nn[8] = {0,0,0,0,0,0,0,0};
  LoadBadTowers();

  int twr;
  float mea;
  for(;;) {
    fin >> twr;
    if(!fin.good()) break;
    fin >> mea;
    if(isbad[twr]!=0) continue;
    int sec;
    for(int i=0; i!=8; ++i) {
      if(twr>secs[i]) sec = i;
    }
    //cout << twr << " " << sec << endl;
    x[sec][nn[sec]] = twr;
    y[sec][nn[sec]] = mea;
    nn[sec]++;
    if( fabs(mea) > 0.5 ) {
      cout << sec << " || " << twr << " || " << mea << endl;
    }
  }

  TGraph *gr[8];
  TCanvas * main = new TCanvas();
  main->Divide(4,2);
  for(int i=0; i!=8; ++i) {
    gr[i] = new TGraph(nn[i],x[i],y[i]);
    main->cd(i+1);
    gr[i]->Draw("A*L");
  }
  return 0;
}
