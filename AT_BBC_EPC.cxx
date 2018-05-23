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
    for(int k=0; k!=3; ++k) { // order
      for(int j=0; j!=2; ++j) { // subevent
	hQxVtx[k][j][i] = new TH2F( Form("BBCQ%dx_S%d_CB%02d",k+1,j,i), "BBCQx;VtxBin",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 60, -30, +30 );
	hQyVtx[k][j][i] = new TH2F( Form("BBCQ%dy_S%d_CB%02d",k+1,j,i), "BBCQy;VtxBin",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 60, -30, +30 );
      }
      hQx_RC[k][i] = new TH1F( Form("BBCQ%dx_RC_CB%02d",k,i), "BBCQx", 60, -30, +30 );
      hQy_RC[k][i] = new TH1F( Form("BBCQ%dy_RC_CB%02d",k,i), "BBCQy", 60, -30, +30 );
      hPsiC[k][i] = new TProfile2D( Form("BBCPsiC_Ord%d_Cen%02d",k,i), "PsiC",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
      hPsiS[k][i] = new TProfile2D( Form("BBCPsiS_Ord%d_Cen%02d",k,i), "PsiS",
				    fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
    }
  }
}

void AT_BBC_EPC::MyFinish() {
  for(int i=0; i!=fNBinsCen; ++i) { // centrality
    for(int k=0; k!=2; ++k) { // order
      for(int j=0; j!=2; ++j) { // subevent
	hQxVtx[k][j][i]->Write();
	hQyVtx[k][j][i]->Write();
      }
      hQx_RC[k][i]->Write();
      hQy_RC[k][i]->Write();
      hPsiC[k][i]->Write();
      hPsiS[k][i]->Write();
    }
  }
}

void AT_BBC_EPC::MyExec() {
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );

  qcQ bbq[3][3];
  for(int se=0; se!=2; ++se) {
    bbq[0][se] = pQ1bb->at(se);
    bbq[1][se] = pQ2bb->at(se);
    bbq[2][se] = pQ3bb->at(se);
  }

  // ======= STAGE 1: Storing Raw Centroids =======
  for(int j=0; j!=2; ++j) { // subevent
    for(int k=0; k!=2; ++k) { //order
      hQxVtx[k][j][bcen]->Fill( bvtx, bbq[k][j].X() );
      hQyVtx[k][j][bcen]->Fill( bvtx, bbq[k][j].Y() );
    }
  }


  // ======= STAGE 2: Recentering SubEvents and Building FullQ  =======
  for(int k=0; k!=3; ++k) { // order
    for(int j=0; j!=2; ++j) { // subevent
      bbq[k][j].SetXY( bbq[k][j].X() - bbcm[j][k][0][bcen][bvtx],
		       bbq[k][j].Y() - bbcm[j][k][1][bcen][bvtx],
		       bbq[k][j].NP(),
		       bbq[k][j].M() );
    }
    bbq[k][2] = bbq[k][0] + bbq[k][1];
    hQx_RC[k][bcen]->Fill( bbq[k][2].X() );
    hQy_RC[k][bcen]->Fill( bbq[k][2].Y() );
  }


  // ======= STAGE 3A: Storing Flattening Coeficients  =======
  for(int k=0; k!=3; ++k) { // order
    for(int ik=0; ik!=32; ++ik) { // correction order
      int nn = 1+ik;
      hPsiC[k][bcen]->Fill(bvtx, nn, TMath::Cos(nn*bbq[k][2].Psi2Pi()) );
      hPsiS[k][bcen]->Fill(bvtx, nn, TMath::Sin(nn*bbq[k][2].Psi2Pi()) );
    }
  }

}
