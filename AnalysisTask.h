#ifndef __ANALYSIS_TASK_HH__
#define __ANALYSIS_TASK_HH__

#include <iostream>
#include <vector>
#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>
#include "qcQ.h"

class AnalysisTask : public TObject { // needed to add to TLists
 public:
  AnalysisTask() {
    fCandidates = NULL;
    fCandidates2 = NULL;
    fQ[0]=fQ[1]=fQ[2]=fQ[3]=NULL;
  }
  virtual ~AnalysisTask() {}
  virtual void Init() {std::cout << "AT::INIT" << std::endl;}
  virtual void Exec() {std::cout << "AT::EXEC" << std::endl;}
  virtual void Finish() {std::cout << "AT::FINISH" << std::endl;}

 protected:
  std::vector<TLorentzVector> *fCandidates;
  std::vector<TLorentzVector> *fCandidates2;
  qcQ *fQ[4];
};

#endif
