#ifndef __AT_BBC_EPC_HH__
#define __AT_BBC_EPC_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;
class TProfile2D;

class AT_BBC_EPC : public AT_ReadTree {
 public:
  AT_BBC_EPC();
  virtual ~AT_BBC_EPC();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH2F *hQxVtx[3][2][60];
  TH2F *hQyVtx[3][2][60];

  TH1F *hQx_RC[3][60];
  TH1F *hQy_RC[3][60];

  TProfile2D *hPsiC[3][60];
  TProfile2D *hPsiS[3][60];
};

#endif
