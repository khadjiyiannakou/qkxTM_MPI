#include <qkxTM.h>
#include <lattice_util.h>
#include <blas_qkxTM.h>
#include <solvers.h>

using namespace qkxTM;

CG::CG(LatticeInfo *param , Gauge &gauge) : u(gauge) {

  maxIter = param->maxIter;
  tol = param->tol;
  reliableDelta = param->reliableDelta;
 
}

void CG::operator()(qkxTM::Vector &out , qkxTM::Vector &in, LatticeInfo *param){
  
  int rank = out.rank;

  Vector *r_P = new Vector(param);
  Vector *p_P = new Vector(param);
  Vector *Ap_P = new Vector(param);
  Vector *tmp_P = new Vector(param);

  // take reference
  Vector &r = *r_P;
  Vector &p = *p_P;
  Vector &Ap = *Ap_P;
  Vector &tmp = *tmp_P;

  r.copy(in);  // r_ 0 = b - A x_0  initial guess zero
  p.copy(r);   // p_0 = r_0  

  double r2;
  double r2_old;
  double stop;
  double src_norm;
  double pAp;
  double alpha;
  double beta;
  int k = 0;           // number of iterations

  r2 = reCdotXdagY(r,r);
  src_norm = reCdotXdagY(in,in);
  stop = src_norm*tol*tol;

  while( (r2 > stop) && (k < maxIter) ){

    printfQKXTM("CG : iteration = %d : r2 = %e\n",k , r2);

    
    Ap.tmDagtm(p,tmp,u,param);
    
    pAp = reCdotXdagY(p,Ap); // p A p
    alpha = r2 / pAp;         // a = r^T r / ( p^T A p)
    r2_old = r2;
    X_eq_X_p_aY(out, alpha, p);  // x = x + a p
    X_eq_X_p_aY(r, -alpha, Ap);   // r = r - a Ap
    r2 = reCdotXdagY(r,r);
    beta = r2/r2_old;       // b = r_(k+1)^T r_(k+1) / ( r_k^T r_k)
    X_eq_aX_p_Y(p, beta, r);  // p = r + beta p
    
    k++;
    

  } //close while loop

  if(k == maxIter){
    printfQKXTM("Warning solver reached number of max iterations \n");
  }
  else
    {
      printfQKXTM("CG converged after iteration %d\n",k);
    }

  delete r_P;
  delete p_P;
  delete Ap_P;
  delete tmp_P;

}
