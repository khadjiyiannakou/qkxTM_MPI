#ifndef _QKXTM_H
#define _QKXTM_H

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <blas.h>

typedef std::complex<double> Complex;

/*
inline Complex operator*(const double a , const Complex b)
{
  Complex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}

inline Complex operator*(const Complex b , const double a)
{
  Complex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}

inline Complex operator*(const int a , const Complex b)
{
  Complex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}

inline Complex operator*(const Complex b , const int a)
{
  Complex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}
*/
///////////////////// Classes for matrix vector operators for 3x3 ////////////
class Vector_3 {
public:
  Complex M[3];
  Vector_3(){;}
  Vector_3(Complex a[3]){
    M[0]=a[0];
    M[1]=a[1];
    M[2]=a[2];
  }
  ~Vector_3(){;}
  void print(){
    for(int c = 0 ; c<3 ; c++)
      printf("%f %f\n",M[c].real(),M[c].imag());
  }
  Vector_3 operator+(Vector_3 v){
    Vector_3 res;
    res.M[0]=M[0]+v.M[0];
    res.M[1]=M[1]+v.M[1];
    res.M[2]=M[2]+v.M[2];
    return res;
  }
};
  
class Matrix_3x3 {
public:
  Complex M[3][3];
  Matrix_3x3(){;}
  Matrix_3x3(Complex a[3][3]){
    for(int c1 = 0 ; c1 < 3 ; c1++)
      for(int c2 = 0 ; c2 < 3 ; c2++)
	M[c1][c2] = a[c1][c2];
  }
  ~Matrix_3x3(){;}
  void print(){
    for(int c1 = 0; c1 < 3 ; c1++){
      printf("(%f,%f) \t (%f,%f) \t (%f,%f)\n",M[c1][0].real(),M[c1][0].imag(),M[c1][1].real(),M[c1][1].imag(),M[c1][2].real(),M[c1][2].imag());
    }
  }
  Vector_3 operator*(Vector_3 v){
    Vector_3 res;
    res.M[0]=M[0][0]*v.M[0]+M[0][1]*v.M[1]+M[0][2]*v.M[2];
    res.M[1]=M[1][0]*v.M[0]+M[1][1]*v.M[1]+M[1][2]*v.M[2];
    res.M[2]=M[2][0]*v.M[0]+M[2][1]*v.M[1]+M[2][2]*v.M[2];
    return res;
  }
};

template<typename T>
Vector_3 operator*(T a, Vector_3 x){
  Vector_3 res;
  res.M[0]=a*x.M[0];
  res.M[1]=a*x.M[1];
  res.M[2]=a*x.M[2];
  return res;
}

template<typename T>
Matrix_3x3 operator*(T a, Matrix_3x3 X){
  Matrix_3x3 res;
  for(int c1 = 0 ; c1 < 3 ; c1++)
    for(int c2 = 0 ; c2 < 3 ; c2++)
      res.M[c1][c2] = a*X.M[c1][c2];
  return res;
}

/////////////////////////////////////


#define NDIM 4  // the timeSpace has always 4 dimensions
#define NSPINS 4
#define NCOLORS 3
#define NDF 2 // for real and imag part
#define MOM_MAX 1000
#define NSPINOR (NSPINS*NCOLORS*NDF)
#define NLINKS (NDIM*NCOLORS*NCOLORS*NDF)
#define NPROPS (NSPINS*NSPINS*NCOLORS*NCOLORS*NDF)

#define MS(spin,color) ( (spin)*NCOLORS + color )                          // macro for spinor
#define MG(v,color1,color2) ( (v)*NCOLORS*NCOLORS + (color1)*NCOLORS + color2 )                           // macro for gauge field
#define MU(color1,color2) ( (color1)*NCOLORS + color2 )

#define LEXIC(it,iz,iy,ix,L) ( (it)*L[0]*L[1]*L[2] + (iz)*L[0]*L[1] + (iy)*L[0] + (ix) )
#define LEXIC_TZY(it,iz,iy,L) ( (it)*L[1]*L[2] + (iz)*L[1] + (iy) )
#define LEXIC_TZX(it,iz,ix,L) ( (it)*L[0]*L[2] + (iz)*L[0] + (ix) )
#define LEXIC_TYX(it,iy,ix,L) ( (it)*L[0]*L[1] + (iy)*L[0] + (ix) )
#define LEXIC_ZYX(iz,iy,ix,L) ( (iz)*L[0]*L[1] + (iy)*L[0] + (ix) )
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define errorQKXTM(...) do {			\
  fprintf(stderr,__VA_ARGS__);			\
  MPI_Abort(MPI_COMM_WORLD,-1);			\
  exit(-1);					\
  }while(0)

#define printfQKXTM(...) do {			\
  if(rank == 0)					\
    fprintf(stdout,__VA_ARGS__);		\
  } while(0)





///////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////
struct LatticeInfo{

  // index 0 -> x direction
  // index 1 -> y direction
  // index 2 -> z direction
  // index 3 -> t direction

  int L[4]; // stores the lattice length in each direction
  int P[4]; // process in each direction
  int x_src[4];
  int Nprocs;
  int rank;
  double kappa;
  double csw; // clover coefficient
  double mu;
  int twistSign; // +1 or -1 only
  double tol;        // need it when try inverter
  double maxIter;
  double reliableDelta;
  int NsmearAPE;
  int NsmearGauss;
  double alphaAPE;
  double alphaGauss;
  int boundaryT;   // +1 for periodic boundary conditions -1 for anti-periodic boundary conditions
  int Q_sq;

  LatticeInfo()
  {
    L[0] = 0 ; L[1] = 0 ; L[2] = 0 ; L[3] = 0;
    P[0] = 0 ; P[1] = 0 ; P[2] = 0 ; P[3] = 0;
    x_src[0] = 0 ; x_src[1] = 0 ; x_src[2] = 0 ; x_src[3] = 0 ;
    Nprocs = 0;
    rank = -1;
    kappa=0.16;
    mu=0.008;
    twistSign=+1;
    tol = 1e-1;
    maxIter = 1000;
    reliableDelta=0.1;
    NsmearAPE = 20;
    NsmearGauss = 50;
    alphaAPE = 0.5;
    alphaGauss = 4.0;
    boundaryT = -1;
    Q_sq = 0;
    csw=1.57551;
  }

};

struct Pointers{
  Complex *val;
  int **ind;
  Complex *val2[3];
  int *ind2[3][4];

  Pointers(){
    val = NULL;
    ind = NULL;
    for(int i = 0 ; i < 3 ; i++){
      val2[i] = NULL;
      for(int mu = 0 ; mu < 4 ; mu++)
	ind2[i][mu] = NULL;
    }
  }
  
};
//////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////

void setGrid(int x[], int iv , const int gridDim[]);
void* safe_malloc(size_t size);
void getCoord(const int iv, int x[], const int L[]);
int getNeighborPlus(int it, int iz, int iy, int ix, const int L[], int mu);
int getNeighborMinus(int it, int iz, int iy, int ix, const int L[], int mu);
//int qkxTM_getLimeMessage(char *fname);
//int qkxTM_getGaugeLime(char *fname, Complex **u, LatticeInfo *latInfo);



#endif
