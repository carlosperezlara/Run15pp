#ifndef __AT_EVENTPLANECALIBRATOR_HH__
#define __AT_EVENTPLANECALIBRATOR_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;
class TProfile2D;

class AT_EventPlaneCalibrator : public AT_ReadTree {
 public:
  AT_EventPlaneCalibrator();
  virtual ~AT_EventPlaneCalibrator() {}
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH2F *hBBCQ1xVtx[2][60];
  TH2F *hBBCQ1yVtx[2][60];
  TH2F *hBBCQ2xVtx[2][60];
  TH2F *hBBCQ2yVtx[2][60];
  TH2F *hBBCQ3xVtx[2][60];
  TH2F *hBBCQ3yVtx[2][60];

  TH1F *hBBCQ1x_Step2[60];
  TH1F *hBBCQ1y_Step2[60];
  TH1F *hBBCQ2x_Step2[60];
  TH1F *hBBCQ2y_Step2[60];
  TH1F *hBBCQ3x_Step2[60];
  TH1F *hBBCQ3y_Step2[60];

  TProfile2D *hBBCPsiC[3][60];
  TProfile2D *hBBCPsiS[3][60];

};

#endif
