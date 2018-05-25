#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TMath.h>

#include "Analysis.h"
#include "AT_BBC_EPC.h"

AT_BBC_EPC::AT_BBC_EPC() : AT_ReadTree() {
}

AT_BBC_EPC::~AT_BBC_EPC() {
  //Should delete histograms, otherwise leak.
  //But this is the end anyway, so it does not matter
}
void AT_BBC_EPC::MyInit() {
  Analysis *ana = Analysis::Instance();
  int run = ana->RunNumber();
  for(int i=0; i!=fNBinsCen; ++i) { //centrality
    for(int k=0; k!=6; ++k) { // order
      for(int j=0; j!=2; ++j) { // subevent
	hQxVtx[k][j][i] = new TH2F( Form("BBCQ%dx_S%d_CB%02d",k+1,j,i), "BBCQx;VtxBin",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 60, -50, +50 );
	hQyVtx[k][j][i] = new TH2F( Form("BBCQ%dy_S%d_CB%02d",k+1,j,i), "BBCQy;VtxBin",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 60, -50, +50 );
      }
    }
    for(int k=0; k!=4; ++k) { // order
      hPsiC[k][i] = new TProfile2D( Form("BBCPsiC_Ord%d_Cen%02d",k,i), "PsiC",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
      hPsiS[k][i] = new TProfile2D( Form("BBCPsiS_Ord%d_Cen%02d",k,i), "PsiS",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
      for(int j=0; j!=2; ++j) { // subevent
	for(int n=0; n!=3; ++n) { //step
	  hQxC[n][k][j][i] = new TH1F( Form("BBCQ%dxC_S%d_CB%02d_ST%d",k+1,j,i,n),
				       "BBCQx", 60, -30, +30 );
	  hQyC[n][k][j][i] = new TH1F( Form("BBCQ%dyC_S%d_CB%02d_ST%d",k+1,j,i,n),
				       "BBCQy", 60, -30, +30 );
	}
      }
    }
  }
}

void AT_BBC_EPC::MyFinish() {
  for(int i=0; i!=fNBinsCen; ++i) { // centrality
    for(int k=0; k!=6; ++k) { // order
      for(int j=0; j!=2; ++j) { // subevent
	hQxVtx[k][j][i]->Write();
	hQyVtx[k][j][i]->Write();
      }
    }
    for(int k=0; k!=4; ++k) { // order
      hPsiC[k][i]->Write();
      hPsiS[k][i]->Write();
      for(int j=0; j!=2; ++j) { // subevent
	for(int n=0; n!=3; ++n) { //step
	  hQxC[n][k][j][i]->Write();
	  hQyC[n][k][j][i]->Write();
	}
      }
    }
  }
}

void AT_BBC_EPC::MyExec() {
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );

  qcQ qvec[6][3];
  for(int se=0; se!=2; ++se) {
    qvec[0][se] = pQ1bb->at(se);
    qvec[1][se] = pQ2bb->at(se);
    qvec[2][se] = pQ3bb->at(se);
    qvec[3][se] = pQ4bb->at(se);
    qvec[4][se] = pQ6bb->at(se);
    qvec[5][se] = pQ8bb->at(se);
    if(qvec[0][se].M()<1) return;
  }

  // ======= STAGE 1: Storing Raw Centroids =======
  for(int j=0; j!=2; ++j) { // subevent
    for(int k=0; k!=6; ++k) { //order
      hQxVtx[k][j][bcen]->Fill( bvtx, qvec[k][j].X() );
      hQyVtx[k][j][bcen]->Fill( bvtx, qvec[k][j].Y() );
    }
  }

  // ======= STAGE 2: Recentering SubEvents (STEP1)  =======
  for(int k=0; k!=6; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double cn = bbcm[j][k][0][bcen][bvtx];
      double sn = bbcm[j][k][1][bcen][bvtx];
      qvec[k][j].SetXY( x - cn, y - sn, qvec[k][j].NP(), qvec[k][j].M() );
    }
  }

  // ======= STAGE 3: Storing Prime Centroids =======
  for(int j=0; j!=2; ++j) { // subevent
    for(int k=0; k!=4; ++k) { //order
      hQxC[0][k][j][bcen]->Fill( qvec[k][j].X() );
      hQyC[0][k][j][bcen]->Fill( qvec[k][j].Y() );
    }
  }

  int twon[4] = {1,3,4,5}; // 1,2,3,4,6,8
  // ======= STAGE 4: Twisting SubEvents (STEP2)  =======
  for(int k=0; k!=4; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double c2n = bbcm[j][twon[k]][0][bcen][bvtx] / qvec[k][j].M();
      double s2n = bbcm[j][twon[k]][1][bcen][bvtx] / qvec[k][j].M();
      double ldaSm = s2n/(1.0+c2n);
      double ldaSp = s2n/(1.0-c2n);
      double den = 1.0 - ldaSm*ldaSp;
      qvec[k][j].SetXY( (x-ldaSm*y) / den,
			(y-ldaSp*x) / den,
			qvec[k][j].NP(),
			qvec[k][j].M() );
      if( TMath::IsNaN( qvec[k][j].X() ) || TMath::IsNaN( qvec[k][j].Y() ) ) {
	std::cout << "Error building coefficient ";
	std::cout << " | qvec.M: " << qvec[k][j].M();
	std::cout << " | c2n: " << c2n;
	std::cout << " | s2n: " << s2n;
	std::cout << " | ldaSm: " << ldaSm;
	std::cout << " | ldaSp: " << ldaSp;
	std::cout << " | den: " << den << std::endl;
      }
    }
  }

  // ======= STAGE 5: Storing Double-Prime Centroids =======
  for(int j=0; j!=2; ++j) { // subevent
    for(int k=0; k!=4; ++k) { //order
      hQxC[1][k][j][bcen]->Fill( qvec[k][j].X() );
      hQyC[1][k][j][bcen]->Fill( qvec[k][j].Y() );
    }
  }

  // ======= STAGE 6: Rescaling SubEvents (STEP3)  =======
  for(int k=0; k!=4; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      double x = qvec[k][j].X();
      double y = qvec[k][j].Y();
      double c2n = bbcm[j][twon[k]][0][bcen][bvtx] / qvec[k][j].M();
      double a2np = 1.0+c2n;
      double a2nm = 1.0-c2n;
      qvec[k][j].SetXY( x / a2np,
			y / a2nm,
			qvec[k][j].NP(),
			qvec[k][j].M() );
      if( TMath::IsNaN( qvec[k][j].X() ) || TMath::IsNaN( qvec[k][j].Y() ) ) {
	std::cout << "Error building coefficient [2] ";
	std::cout << " | a2np: " << a2np;
	std::cout << " | a2nm: " << a2nm << std::endl;
      }
    }
  }

  // ======= STAGE 7: Storing Triple-Prime Centroids =======
  for(int j=0; j!=2; ++j) { // subevent
    for(int k=0; k!=4; ++k) { //order
      hQxC[2][k][j][bcen]->Fill( qvec[k][j].X() );
      hQyC[2][k][j][bcen]->Fill( qvec[k][j].Y() );
    }
  }

  // ======= STAGE 8: Bulding Full Q and Storing Flattening Coeficients  =======
  for(int k=0; k!=4; ++k) { // order
    qvec[k][2] = qvec[k][0] + qvec[k][1];
    for(int ik=0; ik!=32; ++ik) { // correction order
      int nn = 1+ik;
      hPsiC[k][bcen]->Fill(bvtx, nn, TMath::Cos(nn*qvec[k][2].Psi2Pi()) );
      hPsiS[k][bcen]->Fill(bvtx, nn, TMath::Sin(nn*qvec[k][2].Psi2Pi()) );
    }
  }
}
