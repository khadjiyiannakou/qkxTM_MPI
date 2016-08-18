#include <blas_qkxTM.h>
#include <qkxTM.h>

double reCdotXdagY(qkxTM::Vector &x , qkxTM::Vector &y){

  Complex res;
  Complex globalRes;

  res.real() = 0.;
  res.imag() = 0.;

  globalRes.real() = 0.;
  globalRes.imag() = 0.;

  int lV4d;
  lV4d = x.lV4d;

  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        res = res + conj(x.M[iv][mu][ic]) * y.M[iv][mu][ic];
      }

  MPI_Allreduce(&(res.real()) , &(globalRes.real()) , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  return globalRes.real();
}


void X_eq_X_p_aY(qkxTM::Vector &x , double a , qkxTM::Vector &y){

  int lV4d;
  lV4d = x.lV4d;

  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        x.M[iv][mu][ic] = x.M[iv][mu][ic] + a * y.M[iv][mu][ic];
      }

}


void X_eq_aX_p_Y(qkxTM::Vector &x , double a , qkxTM::Vector &y){

  int lV4d;
  lV4d = x.lV4d;

  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        x.M[iv][mu][ic] = a * x.M[iv][mu][ic] + y.M[iv][mu][ic];
      }

}
