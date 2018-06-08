#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include "AT_PiZero.h"
#include "EmcIndexer.h"
#include "EmcIndexer.C"
#include "PbGlIndexer.C"
#include "PbScIndexer.C"

AT_PiZero::AT_PiZero() : AT_ReadTree() {
  for (int i = 0; i < 8; ++i){
    for (int j = 0; j < 48; ++j){
      for (int k = 0; k < 96; ++k){
	EMCMAP[i][j][k] = 0;
      }
    }
  }
  TString fname = "Run16dAu200WarnMap.list";
  std::cout << "AT_PiZero::Ctor === is reading EMCal dead map: ";
  std::cout << fname.Data() << std::endl;
  int armsect = 0, ypos = 0, zpos = 0, status = 0;
  ifstream readmap( fname.Data() );
  while(readmap >> armsect >> ypos >> zpos >> status) {
    EMCMAP[armsect][ypos][zpos] = status;
    //if(status==-1)EMCMAP[armsect][ypos][zpos] = 0; // this is for ERT trigger
  }
  readmap.close();
  fQA = false;
  hVertex = NULL;
  hCentrality = NULL;
  hNClu0 = NULL;
  hNClu1 = NULL;
  for(int j=0; j!=8; ++j) {
    for(int i=0; i!=4; ++i) {
      hPizeroMass[i][j] = NULL;
    }
    hPizeroMixMass[j] = NULL;
  }
}
AT_PiZero::~AT_PiZero() {
}

void AT_PiZero::MyInit() {
  if(fQA) {
    hVertex = new TH1F("Vertex","",100,-30,+30);
    hCentrality = new TH1F("Centrality","",100,0,100);
    hNClu0 = new TProfile("hNClu0","<NClu> exc. bad towers",8,-0.5,7.5);
    hNClu1 = new TProfile("hNClu1","<NClu> exc. bad towers and timing",8,-0.5,7.5);
    for(int j=0; j!=8; ++j) {
      for(int i=0; i!=4; ++i) {
	hPizeroMass[i][j] = new TH2F(Form("PizeroMass_Cut%d_Sector%d",i,j),
				     "Mass;pT", 100,0.0,0.7, 150,0,15);
      }
      hPizeroMixMass[j] = new TH2F(Form("PizeroMixMass_Sector%d",j),
				   "Mass;pT", 100,0.0,0.7, 150,0,15);
    }
  }
}

void AT_PiZero::MyFinish() {
  if(fQA) {
    hVertex->Write();
    hCentrality->Write();
    hNClu0->Write();
    hNClu1->Write();
    for(int j=0; j!=8; ++j) {
      for(int i=0; i!=4; ++i) {
	hPizeroMass[i][j]->Write();
      }
      hPizeroMixMass[j]->Write();
    }
  }
}

void AT_PiZero::MyExec() {
  fCandidates->clear();
  fCandidates2->clear();
  
  //====== EVENT SELECTION ======
  float cent = fGLB.cent;
  float frac = fGLB.frac;
  float vtxZ = fGLB.vtxZ;
  unsigned int trigger = fGLB.trig;
  bool trig = false;
  if(trigger & fMask) trig = true;
  if(cent<fCentralityMin||cent>fCentralityMax) return;
  if(frac<0.95) return;
  if(!trig) return;
  if(fabs(vtxZ)>20) return;
  
  if(fQA) hCentrality->Fill(cent);
  if(fQA) hVertex->Fill(vtxZ);
  //============
  int binvertex = P0_VertexBin(vtxZ);
  if(binvertex<0) return;
  hEvents->Fill(2);
  
  //====== MAIN LOOP ON CLUSTERS ======
  fBuffer.clear();
  int nclu0[8] = {0,0,0,0,0,0,0,0};
  int nclu1[8] = {0,0,0,0,0,0,0,0};
  int isc, jsc;
  int y, z;
  uint nclu = pEMCecore->size();
  for(uint icl=0; icl!=nclu; ++icl) {
    int idx = pEMCtwrid->at(icl);
    float it = pEMCtimef->at(icl);
    EmcIndexer::decodeTowerId(idx,isc,z,y);
    if( IsBad(isc,y,z) ) continue;
    nclu0[isc]++;
    if( fabs(it)<5 ) nclu1[isc]++;
    //=== loading cluster i
    float iecore = pEMCecore->at(icl);
    float ix = pEMCx->at(icl);
    float iy = pEMCy->at(icl);
    float iz = pEMCz->at(icl) - vtxZ;
    double idl = TMath::Sqrt(ix*ix + iy*iy + iz*iz);
    //=== buffering
    FASTCLU clu;
    clu.ecore = iecore;
    clu.idx = idx;
    clu.x = ix;
    clu.y = iy;
    clu.z = iz;
    clu.t = it;
    fBuffer.push_back( clu );
    //=== continue
    TLorentzVector ii;
    ii.SetPx(iecore*ix/idl);
    ii.SetPy(iecore*iy/idl);
    ii.SetPz(iecore*iz/idl);
    ii.SetE(iecore);
    // building foreground
    for(uint jcl=icl+1; jcl<nclu; ++jcl) {
      int jdx = pEMCtwrid->at(jcl);
      float jt = pEMCtimef->at(jcl);
      EmcIndexer::decodeTowerId(jdx,jsc,z,y);
      if(isc!=jsc) continue;
      if( IsBad(jsc,y,z) ) continue;
      //=== loading cluster j
      float jecore = pEMCecore->at(jcl);
      float jidx = pEMCtwrid->at(jcl);
      float jx = pEMCx->at(jcl);
      float jy = pEMCy->at(jcl);
      float jz = pEMCz->at(jcl) - vtxZ;
      double jdl = TMath::Sqrt(jx*jx + jy*jy + jz*jz);
      TLorentzVector jj;
      jj.SetPx(jecore*jx/jdl);
      jj.SetPy(jecore*jy/jdl);
      jj.SetPz(jecore*jz/jdl);
      jj.SetE(jecore);
      
      //=== building pair
      TLorentzVector pp = ii + jj;
      double dist = TMath::Sqrt( +(ix-jx)*(ix-jx)
				 +(iy-jy)*(iy-jy)
				 +(iz-jz)*(iz-jz) );
      float alpha = TMath::Abs(iecore-jecore)/(iecore+jecore);
      if(pp.Pt()<0.8) continue;
      if(pp.Pt()>16.0) continue;
      if(fQA) hPizeroMass[0][isc]->Fill( pp.M(),pp.Pt()); // step0
      if(dist<8) continue;
      if(fQA) hPizeroMass[1][isc]->Fill( pp.M(),pp.Pt()); // step1
      if(alpha>0.8) continue;
      if(fQA) hPizeroMass[2][isc]->Fill( pp.M(),pp.Pt()); // step2
      if( fabs(it)>5 )continue;
      if( fabs(jt)>5 )continue;
      if(fQA) hPizeroMass[3][isc]->Fill( pp.M(),pp.Pt()); // step3
      fCandidates->push_back( pp );
    }
    // building background
    for(uint jcl=0; jcl!=fPrevious[binvertex].size(); ++jcl) {
      float jecore = (fPrevious[binvertex].at(jcl)).ecore;
      float jdx = (fPrevious[binvertex].at(jcl)).idx;
      float jx = (fPrevious[binvertex].at(jcl)).x;
      float jy = (fPrevious[binvertex].at(jcl)).y;
      float jz = (fPrevious[binvertex].at(jcl)).z;
      float jt = (fPrevious[binvertex].at(jcl)).t;
      EmcIndexer::decodeTowerId(jdx,jsc,z,y);
      //std::cout << "GOOD! " << isc << " " << jsc << std::endl;
      if(isc!=jsc) continue;
      if( IsBad(jsc,y,z) ) continue;
      //std::cout << "VGOOD! " << ix << std::endl;
      double jdl = TMath::Sqrt(jx*jx + jy*jy + jz*jz);
      TLorentzVector jj;
      jj.SetPx(jecore*jx/jdl);
      jj.SetPy(jecore*jy/jdl);
      jj.SetPz(jecore*jz/jdl);
      jj.SetE(jecore);
      TLorentzVector pp = ii + jj;
      double dist = TMath::Sqrt( +(ix-jx)*(ix-jx)
                                 +(iy-jy)*(iy-jy)
                                 +(iz-jz)*(iz-jz) );
      float alpha = TMath::Abs(iecore-jecore)/(iecore+jecore);
      if(pp.Pt()<0.8) continue;
      if(pp.Pt()>16.0) continue;
      if(dist<8) continue;
      if(alpha>0.8) continue;
      if( fabs(it)>5 )continue;
      if( fabs(jt)>5 )continue;
      //std::cout << "VVGOOD! " << ix << std::endl;
      if(fQA) hPizeroMixMass[isc]->Fill( pp.M(),pp.Pt());
      fCandidates2->push_back( pp );
    }
  }

  if(fBuffer.size()>0) {
    //d::cout << "Filling vertex " << binvertex << std::endl;
    fPrevious[binvertex].clear();
    for(uint i=0; i!=fBuffer.size(); ++i) {
      fPrevious[binvertex].push_back( fBuffer[i] );
    }
  }

  if(fQA) {
    for(int i=0; i!=8; ++i) {
      hNClu0->Fill( i, nclu0[i] );
      hNClu1->Fill( i, nclu1[i] );
    }
  }
}

bool AT_PiZero::IsBad(int isc, int y, int z) {
  bool ret = false;

  // convert to veronica convention for sectors
  int vSc = isc;
  if(isc==6) vSc=7;
  if(isc==7) vSc=6;
  if(isc==4) vSc=5;
  if(isc==5) vSc=4;
  isc = vSc;
  //

  if( y==0 || z==0 ) ret = true;
  if( isc < 6 && ( y == 35 || z == 71) ) ret = true;
  if( isc > 5 && ( y == 47 || z == 95) ) ret = true;
  if( EMCMAP[isc][y-1][z-1] || EMCMAP[isc][y][z-1] || EMCMAP[isc][y+1][z-1] ||
      EMCMAP[isc][y-1][z]   || EMCMAP[isc][y][z]   || EMCMAP[isc][y+1][z] ||
      EMCMAP[isc][y-1][z+1] || EMCMAP[isc][y][z+1] || EMCMAP[isc][y+1][z+1] )
    ret = true;
  return ret;
}

int AT_PiZero::P0_VertexBin(float vtx) {
  int ret = -1;
  float binning[21] = {-20,-18,-16,-14,-12,-10,-8.,-6.,-4.,-2., 0,
		       +2.,+4.,+6.,+8.,+10,+12,+14,+16,+18,+20};
  for(int i=0; i!=20; ++i) {
    if(vtx>binning[i] && vtx<binning[i+1]) ret = i;
  }
  return ret;
}
