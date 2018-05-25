#ifndef __AT_PIZEROMASS_HH__
#define __AT_PIZEROMASS_HH__

#include <vector>
#include "AT_ReadTree.h"
#include "TLorentzVector.h"

class TH1F;
class TH2F;

class AT_PiZeroMass : public AT_ReadTree {
 public:
  AT_PiZeroMass();
  virtual ~AT_PiZeroMass();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  bool IsBad(int sc, int y, int z);
  int  EMCMAP[8][48][96];
  TH1F *hVertex;
  TH1F *hCentrality;
  TH2F *hPizeroMass;
  TH2F *hPizeroM[8];

 protected:
  std::vector<TLorentzVector> fCandidates;
};

#endif
