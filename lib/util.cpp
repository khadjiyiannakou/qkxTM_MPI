#include <stdio.h>
#include <stdlib.h>
#include <qkxTM.h>

void setGrid(int x[], int iv , const int gridDim[]){
  for(int i = 0 ; i < 3 ; i++){
    x[i] = iv % gridDim[i];
    iv -= x[i];
    iv /= gridDim[i];
  }
  x[3] = iv;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void* safe_malloc(size_t size){						
  void *ptr;
  ptr = malloc(size);							
  if(ptr == NULL){							
    errorQKXTM("Error out of memory");
  }
  return ptr;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void getCoord(const int iv, int x[], const int L[]){
  int r,r2;
  r = iv/L[0];
  r2 = r*L[0];
  x[0] = iv - r2;
  r = iv/(L[0]*L[1]);
  r2 = r*(L[0]*L[1]);
  x[1] = (iv - r2 - x[0])/L[0];
  x[3] = iv/(L[0]*L[1]*L[2]);
  r = iv/(L[0]*L[1]);
  x[2] = r - x[3]*L[2];
} 

int getNeighborPlus(int it, int iz, int iy, int ix, const int L[], int mu){
  int res;
  switch(mu){
  case 0:
    res = LEXIC(it,iz,iy,(ix+1)%L[0],L);
    break;
  case 1:
    res = LEXIC(it,iz,(iy+1)%L[1],ix,L);
    break;
  case 2:
    res = LEXIC(it,(iz+1)%L[2],iy,ix,L);
    break;
  case 3:
    res = LEXIC((it+1)%L[3],iz,iy,ix,L);
  }
  return res;
}

int getNeighborMinus(int it, int iz, int iy, int ix, const int L[], int mu){
  int res;
  switch(mu){
  case 0:
    res = LEXIC(it,iz,iy,(ix-1+L[0])%L[0],L);
    break;
  case 1:
    res = LEXIC(it,iz,(iy-1+L[1])%L[1],ix,L);
    break;
  case 2:
    res = LEXIC(it,(iz-1+L[2])%L[2],iy,ix,L);
    break;
  case 3:
    res = LEXIC((it-1+L[3])%L[3],iz,iy,ix,L);
  }
  return res;
}
