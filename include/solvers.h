#include <lattice_util.h>

class CG {

 private:
  int maxIter;
  double tol;
  double reliableDelta;
  qkxTM::Gauge &u;

 public:
  CG(LatticeInfo *param , qkxTM::Gauge &gauge);
  ~CG() { ; }

  /////////// functor    
  int MaxIter() { return maxIter; }
  double Tol() { return tol; }
  double ReliableDelta() { return reliableDelta; }
  ///////////////      

  void operator()(qkxTM::Vector &out , qkxTM::Vector &in, LatticeInfo *param);
};

