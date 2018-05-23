#ifndef __ANALYSIS_TASK_HH__
#define __ANALYSIS_TASK_HH__

#include <iostream>
#include <TObject.h>
#include <TString.h>

class AnalysisTask : public TObject { // needed to add to TLists
 public:
  AnalysisTask() {}
  virtual ~AnalysisTask() {std::cout << "AT::DTOR" << std::endl;}
  virtual void Init() {std::cout << "AT::INIT" << std::endl;}
  virtual void Exec() {std::cout << "AT::EXEC" << std::endl;}
  virtual void Finish() {std::cout << "AT::FINISH" << std::endl;}
};

#endif
