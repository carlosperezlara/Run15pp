#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TMath.h>

#include "Analysis.h"
#include "AT_EventPlaneCalibrator.h"

AT_EventPlaneCalibrator::AT_EventPlaneCalibrator() : AT_ReadTree() {
}

void AT_EventPlaneCalibrator::MyInit() {
  Analysis *ana = Analysis::Instance();
  int run = ana->RunNumber();
  for(int i=0; i!=fNBinsCen; ++i) {
    for(int j=0; j!=2; ++j) { // subevent
      hBBCQ1xVtx[j][i] = new TH2F( Form("BBCQ1x_%d_S%d_CB%02d",run,j,i),
				   "BBCQ1x;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
      hBBCQ1yVtx[j][i] = new TH2F( Form("BBCQ1y_%d_S%d_CB%02d",run,j,i),
				   "BBCQ1y;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
      hBBCQ2xVtx[j][i] = new TH2F( Form("BBCQ2x_%d_S%d_CB%02d",run,j,i),
				   "BBCQ2x;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
      hBBCQ2yVtx[j][i] = new TH2F( Form("BBCQ2y_%d_S%d_CB%02d",run,j,i),
				   "BBCQ2y;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
      hBBCQ3xVtx[j][i] = new TH2F( Form("BBCQ3x_%d_S%d_CB%02d",run,j,i),
				   "BBCQ3x;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
      hBBCQ3yVtx[j][i] = new TH2F( Form("BBCQ3y_%d_S%d_CB%02d",run,j,i),
				   "BBCQ3y;VtxBin",
				   fNBinsVtx, -0.5, fNBinsVtx-0.5,
				   60, -30, +30 );
    }
    hBBCQ1x_Step2[i] = new TH1F( Form("BBCQ1x_step2_CB%02d",i), "BBCQ1x",
				 60, -30, +30 );
    hBBCQ1y_Step2[i] = new TH1F( Form("BBCQ1y_step2_CB%02d",i), "BBCQ1y",
				 60, -30, +30 );
    hBBCQ2x_Step2[i] = new TH1F( Form("BBCQ2x_step2_CB%02d",i), "BBCQ2x",
				 60, -30, +30 );
    hBBCQ2y_Step2[i] = new TH1F( Form("BBCQ2y_step2_CB%02d",i), "BBCQ2y",
				 60, -30, +30 );
    hBBCQ3x_Step2[i] = new TH1F( Form("BBCQ3x_step2_CB%02d",i), "BBCQ3x",
				 60, -30, +30 );
    hBBCQ3y_Step2[i] = new TH1F( Form("BBCQ3y_step2_CB%02d",i), "BBCQ3y",
				 60, -30, +30 );
    for(int j=0; j!=3; ++j) {
      hBBCPsiC[j][i] = new TProfile2D( Form("BBCPsiC_Ord%d_Cen%02d",i,j), "PsiC",
				       fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
      hBBCPsiS[j][i] = new TProfile2D( Form("BBCPsiS_Ord%d_Cen%02d",i,j), "PsiS",
				       fNBinsVtx, -0.5, fNBinsVtx-0.5, 32, 0.5, 32.5 );
    }
  }
}

void AT_EventPlaneCalibrator::MyFinish() {
  for(int i=0; i!=fNBinsCen; ++i) {
    for(int j=0; j!=2; ++j) {
      hBBCQ1xVtx[j][i]->Write();
      hBBCQ1yVtx[j][i]->Write();
      hBBCQ2xVtx[j][i]->Write();
      hBBCQ2yVtx[j][i]->Write();
      hBBCQ3xVtx[j][i]->Write();
      hBBCQ3yVtx[j][i]->Write();
    }
    hBBCQ1x_Step2[i]->Write();
    hBBCQ1y_Step2[i]->Write();
    hBBCQ2x_Step2[i]->Write();
    hBBCQ2y_Step2[i]->Write();
    hBBCQ3x_Step2[i]->Write();
    hBBCQ3y_Step2[i]->Write();
    for(int j=0; j!=3; ++j) {
      hBBCPsiC[j][i]->Write();
      hBBCPsiS[j][i]->Write();
    }
  }
}

void AT_EventPlaneCalibrator::MyExec() {
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

  hBBCQ1xVtx[0][bcen]->Fill( bvtx, bbq[0][0].X() );
  hBBCQ1yVtx[0][bcen]->Fill( bvtx, bbq[0][0].Y() );
  hBBCQ1xVtx[1][bcen]->Fill( bvtx, bbq[0][1].X() );
  hBBCQ1yVtx[1][bcen]->Fill( bvtx, bbq[0][1].Y() );

  hBBCQ2xVtx[0][bcen]->Fill( bvtx, bbq[1][0].X() );
  hBBCQ2yVtx[0][bcen]->Fill( bvtx, bbq[1][0].Y() );
  hBBCQ2xVtx[1][bcen]->Fill( bvtx, bbq[1][1].X() );
  hBBCQ2yVtx[1][bcen]->Fill( bvtx, bbq[1][1].Y() );

  hBBCQ3xVtx[0][bcen]->Fill( bvtx, bbq[2][0].X() );
  hBBCQ3yVtx[0][bcen]->Fill( bvtx, bbq[2][0].Y() );
  hBBCQ3xVtx[1][bcen]->Fill( bvtx, bbq[2][1].X() );
  hBBCQ3yVtx[1][bcen]->Fill( bvtx, bbq[2][1].Y() );

  for(int ord=0; ord!=3; ++ord) {
    for(int se=0; se!=2; ++se) {
      bbq[ord][se].SetXY( bbq[ord][se].X() - bbcm[se][ord][0][bcen][bvtx],
			  bbq[ord][se].Y() - bbcm[se][ord][1][bcen][bvtx],
			  bbq[ord][se].NP(),
			  bbq[ord][se].M() );
    }
    bbq[ord][2] = bbq[ord][0] + bbq[ord][1];
  }

  hBBCQ1x_Step2[bcen]->Fill( bbq[0][2].X() );
  hBBCQ1y_Step2[bcen]->Fill( bbq[0][2].Y() );
  hBBCQ2x_Step2[bcen]->Fill( bbq[1][2].X() );
  hBBCQ2y_Step2[bcen]->Fill( bbq[1][2].Y() );
  hBBCQ3x_Step2[bcen]->Fill( bbq[2][2].X() );
  hBBCQ3y_Step2[bcen]->Fill( bbq[2][2].Y() );

  Psi1_BBC = bbq[0][2].Psi2Pi();
  Psi2_BBC = bbq[1][2].Psi2Pi();
  Psi3_BBC = bbq[2][2].Psi2Pi();

  for(int ord=0; ord!=3; ++ord) {
    for(int iord=0; iord!=32; ++iord) {
      int nn = 1+iord;
      hBBCPsiC[ord][bcen]->Fill(bvtx, nn, TMath::Cos(nn*bbq[ord][2].Psi2Pi()) );
      hBBCPsiS[ord][bcen]->Fill(bvtx, nn, TMath::Sin(nn*bbq[ord][2].Psi2Pi()) );
    }
  }


}
