#ifndef __AT_PIZEROFLOW_HH__
#define __AT_PIZEROFLOW_HH__

#include "AT_ReadTree.h"

class AT_PiZeroMass : public AT_ReadTree {
 public:
  AT_PiZeroMass();
  virtual ~AT_PiZeroMass();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

};

#endif
