#ifndef __AT_EP_HH__
#define __AT_EP_HH__

#include "AnalysisTask.h"

class TH1F;
class TProfile;

class AT_EP : public AnalysisTask {
 public:
  AT_EP();
  virtual ~AT_EP();
  virtual void Init();
  virtual void Exec();
  virtual void Finish();

 private:
  int BinPt(float);
  int BinMass(float);

  int fNpt;
  float fPtBins[100];
  int fNma;
  float fMassBins[100];
  TH1F *hEta;
  TH1F *hEta2;
  TH1F *hMass[100];
  TProfile *hCos[5][100];
  TH1F *hMass2[100];
  TProfile *hCos2[5][100];

};

#endif
