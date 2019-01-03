#ifndef __AT_READTREE_HH__
#define __AT_READTREE_HH__

#include <vector>
#include <TString.h>
#include <TH1F.h>
#include "qcQ.h"
#include "AnalysisTask.h"

class AT_ReadTree : public AnalysisTask {
 public:
  AT_ReadTree();
  virtual ~AT_ReadTree();
  virtual void Init();
  virtual void Exec();
  virtual void Finish();
  virtual void MyInit() {}
  virtual void MyFinish() {}
  virtual void MyExec() {}
  void CheckEP1();
  void CheckEP2();
  int ReferenceTracks();
  void SkipBBCQCal() {fBBCQCal=false;}
  void TriggerMask(unsigned int msk) {fMask=msk;}
  void CentralitySelection(float min, float max)
  {fCentralityMin=min; fCentralityMax=max;}
  void SetSkipDetails() {fSkipDetails = true;}
  void SkipPileUpCuts() {fSkipPileUpCuts = true;}

 private:
  void MakeBBCEventPlanes(int,int);
  void LoadTableEP(int run=-1);
  void LoadTableTime(int run=-1);

 protected:
  int BinVertex(float);
  int BinCentrality(float);

  TH1F *hEvents;
  TH1F *hTriggers0;
  TH1F *hCentrality0;
  TH2F *hPileUpRejectionS;
  TH2F *hPileUpRejectionN;
  TH2F *hCentralitySelection;
  TH1F *hVertex0;
  TH1F *hPsi2[10]; // Check flatness on each EP calib step
  unsigned int fMask;
  float fCentralityMin;
  float fCentralityMax;

  int fNBinsVtx;
  int fNBinsCen;
  float fMinBinVtx;
  float fMinBinCen;

  bool fBBCQCal;
  bool Psi_BBC;
  float Psi1_BBC;
  float Psi2_BBC;
  float Psi3_BBC;
  float Psi4_BBC;
  double bbcm[2][6][2][60][40]; //se ord xy bcen bvtx
  double bbcc[32][4][60][40]; //har ord bcen bvtx
  double bbcs[32][4][60][40]; //har ord bcen bvtx
  double fMXm[8][6][2][60][40]; //se ord xy bcen bvtx
  double fMXc[32][4][60][40]; //har ord bcen bvtx
  double fMXs[32][4][60][40]; //har ord bcen bvtx
  float fTime[3]; //mean rmsS rmsN

  typedef struct MyTreeRegister {
    Float_t vtxZ;
    Float_t cent;
    Float_t bbcs;
    Float_t frac;
    UInt_t  trig;
  } MyTreeRegister_t;
  typedef struct MyTreeRegister2 {
    Float_t bbcn;
    Float_t bbcsTmean;
    Float_t bbcnTmean;
    Double_t bbcsTrms;
    Double_t bbcnTrms;
    UInt_t alltrks;
  } MyTreeRegister2_t;
  MyTreeRegister_t fGLB;
  MyTreeRegister2_t fGLB2;

  std::vector<qcQ> *pQ1ex;
  std::vector<qcQ> *pQ2ex;
  std::vector<qcQ> *pQ3ex;
  std::vector<qcQ> *pQ4ex;
  std::vector<qcQ> *pQ6ex;
  std::vector<qcQ> *pQ8ex;
  std::vector<qcQ> *pQ1fv;
  std::vector<qcQ> *pQ2fv;
  std::vector<qcQ> *pQ3fv;
  std::vector<qcQ> *pQ1bb;
  std::vector<qcQ> *pQ2bb;
  std::vector<qcQ> *pQ3bb;
  std::vector<qcQ> *pQ4bb;
  std::vector<qcQ> *pQ6bb;
  std::vector<qcQ> *pQ8bb;
  
  std::vector<Int_t>   *pEMCid;
  std::vector<Int_t>   *pEMCtwrid;
  std::vector<Float_t> *pEMCx;
  std::vector<Float_t> *pEMCy;
  std::vector<Float_t> *pEMCz;
  std::vector<Float_t> *pEMCecore;
  std::vector<Float_t> *pEMCecent;
  std::vector<Float_t> *pEMCchisq;
  std::vector<Float_t> *pEMCtimef;
  
  //  0 (1)   X1 used
  //  1 (2)   X2 used
  //  2 (4)   UV found
  //  3 (8)   UV unique
  //  4 (16)  PC1 found
  //  5 (48)  PC1 unique
  std::vector<Int_t>   *pTRKqua;
  std::vector<Float_t> *pTRKpt;
  std::vector<Float_t> *pTRKphi;
  std::vector<Float_t> *pTRKpz;
  std::vector<Float_t> *pTRKecore;
  std::vector<Float_t> *pTRKetof;
  std::vector<Int_t>   *pTRKtwrid;
  std::vector<Float_t> *pTRKplemc;
  std::vector<Float_t> *pTRKchisq;// chi2/npe0 -> seme-normalized chi2/dof.
  std::vector<Float_t> *pTRKdphi; // diff in rads btw track model projection and hit in EMC
  std::vector<Float_t> *pTRKdz;   // diff in cm btw track model projection and hit in EMC
  std::vector<Float_t> *pTRKpc3sdphi;
  std::vector<Float_t> *pTRKpc3sdz;
  std::vector<Float_t> *pTRKzed;  // z coord at which track crosses PC1
  std::vector<Float_t> *pTRKdisp; // displacement of ring center wrt Track projection in the RICH PMT array
  std::vector<Float_t> *pTRKprob; // probability that particle shower is electromagnetic
  std::vector<Int_t>   *pTRKcid;
  
  std::vector<Float_t> *pMXSempccent;
  std::vector<Float_t> *pMXSempc3x3;
  std::vector<Float_t> *pMXSpt;
  std::vector<Float_t> *pMXSpz;
  std::vector<Float_t> *pMXSeta;
  std::vector<Float_t> *pMXSphi;
  std::vector<Int_t>   *pMXSflyr;
  std::vector<Float_t> *pMXSsingleD;
  std::vector<Int_t>   *pMXSsingleP;

  bool fSkipDetails;
  bool fSkipPileUpCuts;
  bool fSkipClusters;
  bool fSkipTracks;
  bool fSkipShowers;
};

#endif
