#ifndef __AT_PIDFLOW_HH__
#define __AT_PIDFLOW_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;

class AT_PIDFlow : public AT_ReadTree {
 public:
  AT_PIDFlow();
  virtual ~AT_PIDFlow();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH1F *hPt;
  TH1F *hNTrk;
  TH2F *hPtDPhi[3];
  TH2F *hPtDPhiME[3];
  float fPsi1_BBC_PE;
  float fPsi2_BBC_PE;
  float fPsi3_BBC_PE;

};

#endif
