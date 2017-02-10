#include <qkxTM.h>
#include <lattice_util.h>
#define PI 3.14159265359

//#include <blas_qkxTM.h>

#include <hopping_term.h>
#include <hopping_term_chiral.h>
#define DEBUG_MODE

using namespace qkxTM;


static const int eps[6][3]=
  {
    {0,1,2},
    {2,0,1},
    {1,2,0},
    {2,1,0},
    {0,2,1},
    {1,0,2}
  };

static const int sgn_eps[6]=
  {
    +1,+1,+1,-1,-1,-1
  };







MPI_Group fullGroup , spaceGroup , timeGroup;
MPI_Comm spaceComm , timeComm;
int localRank;
int localSize;
int timeRank;
int timeSize;
bool initializeGroupsFlag = false;


LatticeGeometry::LatticeGeometry(LatticeInfo *latInfo){

  rank = latInfo->rank;

  for(int i = 0 ; i < NDIM ; i++){
    L[i] = latInfo->L[i];
    P[i] = latInfo->P[i];
  }


  if( (L[0]%P[0] !=0 ) || (L[1]%P[1] != 0) || (L[2]%P[2] != 0) || (L[3]%P[3] !=0) ) errorQKXTM("Error lattice cant subdivide to this numbers of processors\n");
  
  if( P[0]*P[1]*P[2]*P[3] != latInfo->Nprocs) errorQKXTM("Error MPI processors dont agree with process you pass\n");
  


  V4d = L[0]*L[1]*L[2]*L[3];
  V3d = L[0]*L[1]*L[2];
  for(int i = 0 ; i < NDIM ; i++) lL[i] = L[i] / P[i];
  lV4d = lL[0]*lL[1]*lL[2]*lL[3];
  lV3d = lL[0]*lL[1]*lL[2];


  setGrid(procPos,rank,P);            // align procs in a 4d grid and give them coordinates
  Nprocs = latInfo->Nprocs;

  procPlus[0] = LEXIC(procPos[3], procPos[2], procPos[1], (procPos[0] + 1)%P[0], P );           // find the next proc in lexicographic order
  procPlus[1] = LEXIC(procPos[3], procPos[2], (procPos[1] + 1)%P[1] , procPos[0], P );
  procPlus[2] = LEXIC(procPos[3], (procPos[2] + 1)%P[2] ,procPos[1] , procPos[0], P );
  procPlus[3] = LEXIC( (procPos[3] + 1)%P[3], procPos[2] ,procPos[1] , procPos[0], P );
  procMinus[0] = LEXIC(procPos[3], procPos[2], procPos[1], (procPos[0] - 1 + P[0])%P[0], P );    // same for previous
  procMinus[1] = LEXIC(procPos[3], procPos[2], (procPos[1] - 1 + P[1])%P[1] , procPos[0], P );
  procMinus[2] = LEXIC(procPos[3], (procPos[2] - 1 + P[2])%P[2] ,procPos[1] , procPos[0], P );
  procMinus[3] = LEXIC( (procPos[3] - 1 + P[3])%P[3], procPos[2] ,procPos[1] , procPos[0], P );


  for(int i = 0 ; i < NDIM ; i++){
    pointPlus[i] = (int*) safe_malloc(lV4d*sizeof(int));  // perform safe malloc for local lattice points indexing
    pointMinus[i] = (int*) safe_malloc(lV4d*sizeof(int));
  }

  
  for(int it = 0 ; it < lL[3] ; it++)
    for(int iz = 0 ; iz < lL[2] ; iz++)
      for(int iy = 0 ; iy < lL[1] ; iy++)
	for(int ix = 0 ; ix < lL[0] ; ix++){
	  int iv = LEXIC(it,iz,iy,ix,lL);
	  pointPlus[0][iv] = LEXIC(it,iz,iy,(ix+1)%lL[0],lL);                               // each lattice point memorize the index of plus and minus point
	  pointPlus[1][iv] = LEXIC(it,iz,(iy+1)%lL[1],ix,lL);
	  pointPlus[2][iv] = LEXIC(it,(iz+1)%lL[2],iy,ix,lL);
	  pointPlus[3][iv] = LEXIC((it+1)%lL[3],iz,iy,ix,lL);
	  pointMinus[0][iv] = LEXIC(it,iz,iy,(ix-1+lL[0])%lL[0],lL);
	  pointMinus[1][iv] = LEXIC(it,iz,(iy-1+lL[1])%lL[1],ix,lL);
	  pointMinus[2][iv] = LEXIC(it,(iz-1+lL[2])%lL[2],iy,ix,lL);
	  pointMinus[3][iv] = LEXIC((it-1+lL[3])%lL[3],iz,iy,ix,lL);
	}

  
  // check if we are using more than one proc we need change index for edge points

  for(int i = 0 ; i < NDIM ; i++){
    placeBoundaryStartPlus[i] = 0;
    placeBoundaryStartMinus[i] = 0;
  }

  long int lastIndex = lV4d;

  if(lL[0] < L[0]){  // this is true always except for npx = 1
    // x is the edge points
    for(int it = 0 ; it < lL[3] ; it++)
      for(int iz = 0 ; iz < lL[2] ; iz++)
	for(int iy = 0 ; iy < lL[1] ; iy++){
	  int iv = LEXIC(it,iz,iy,lL[0]-1,lL);
	  pointPlus[0][iv] = LEXIC_TZY(it,iz,iy,lL) + lastIndex;
	}

    placeBoundaryStartPlus[0] = lastIndex;
    lastIndex += lL[3]*lL[2]*lL[1];
    for(int it = 0 ; it < lL[3] ; it++)
      for(int iz = 0 ; iz < lL[2] ; iz++)
        for(int iy = 0 ; iy < lL[1] ; iy++){
          int iv = LEXIC(it,iz,iy,0,lL);
          pointMinus[0][iv] = LEXIC_TZY(it,iz,iy,lL) + lastIndex;
        }
    placeBoundaryStartMinus[0] = lastIndex;
    lastIndex += lL[3]*lL[2]*lL[1];

  } // close if 

  if(lL[1] < L[1]){  // this is true always except for npy = 1
    // x is the edge points
    for(int it = 0 ; it < lL[3] ; it++)
      for(int iz = 0 ; iz < lL[2] ; iz++)
	for(int ix = 0 ; ix < lL[0] ; ix++){
	  int iv = LEXIC(it,iz,lL[1]-1,ix,lL);
	  pointPlus[1][iv] = LEXIC_TZX(it,iz,ix,lL) + lastIndex;
	}

    placeBoundaryStartPlus[1] = lastIndex ;
    lastIndex += lL[3]*lL[2]*lL[0];

    for(int it = 0 ; it < lL[3] ; it++)
      for(int iz = 0 ; iz < lL[2] ; iz++)
        for(int ix = 0 ; ix < lL[0] ; ix++){
          int iv = LEXIC(it,iz,0,ix,lL);
          pointMinus[1][iv] = LEXIC_TZX(it,iz,ix,lL) + lastIndex ;
        }
    placeBoundaryStartMinus[1] = lastIndex;
    lastIndex += lL[3]*lL[2]*lL[0];

  } // close if 


  if(lL[2] < L[2]){  // this is true always except for npz = 1
    // x is the edge points
    for(int it = 0 ; it < lL[3] ; it++)
      for(int iy = 0 ; iy < lL[1] ; iy++)
	for(int ix = 0 ; ix < lL[0] ; ix++){
	  int iv = LEXIC(it,lL[2]-1,iy,ix,lL);
	  pointPlus[2][iv] = LEXIC_TYX(it,iy,ix,lL) + lastIndex;
	}

    placeBoundaryStartPlus[2] = lastIndex ;
    lastIndex += lL[3]*lL[1]*lL[0];

    for(int it = 0 ; it < lL[3] ; it++)
      for(int iy = 0 ; iy < lL[1] ; iy++)
        for(int ix = 0 ; ix < lL[0] ; ix++){
          int iv = LEXIC(it,0,iy,ix,lL);
          pointMinus[2][iv] = LEXIC_TYX(it,iy,ix,lL) + lastIndex ;
        }
    placeBoundaryStartMinus[2] = lastIndex;
    lastIndex += lL[3]*lL[1]*lL[0];

  } // close if 


  if(lL[3] < L[3]){  // this is true always except for npt = 1
    // x is the edge points
    for(int iz = 0 ; iz < lL[2] ; iz++)
      for(int iy = 0 ; iy < lL[1] ; iy++)
	for(int ix = 0 ; ix < lL[0] ; ix++){
	  int iv = LEXIC(lL[3]-1,iz,iy,ix,lL);
	  pointPlus[3][iv] = LEXIC_ZYX(iz,iy,ix,lL) + lastIndex;
	}

    placeBoundaryStartPlus[3] = lastIndex ;
    lastIndex += lL[2]*lL[1]*lL[0];

    for(int iz = 0 ; iz < lL[2] ; iz++)
      for(int iy = 0 ; iy < lL[1] ; iy++)
        for(int ix = 0 ; ix < lL[0] ; ix++){
          int iv = LEXIC(0,iz,iy,ix,lL);
          pointMinus[3][iv] = LEXIC_ZYX(iz,iy,ix,lL) + lastIndex ;
        }
    placeBoundaryStartMinus[3] = lastIndex;
    lastIndex += lL[2]*lL[1]*lL[0];

  } // close if 


  lSurfaceBoundary[0] = lL[1]*lL[2]*lL[3];
  lSurfaceBoundary[1] = lL[0]*lL[2]*lL[3];
  lSurfaceBoundary[2] = lL[0]*lL[1]*lL[3];
  lSurfaceBoundary[3] = lL[0]*lL[1]*lL[2];

  NtotalBoundary = lastIndex - lV4d;

  // because data are contiguous in memory we need packing for the directions except the slower running (in out case t direction)

  // for vector style
  MPI_Type_vector(lL[1]*lL[2]*lL[3], NSPINOR, lL[0]*NSPINOR,MPI_DOUBLE, &(stypeV[0]));                           
  MPI_Type_vector(lL[2]*lL[3], NSPINOR*lL[0], lL[0]*lL[1]*NSPINOR,MPI_DOUBLE, &(stypeV[1]));
  MPI_Type_vector(lL[3], NSPINOR*lL[0]*lL[1], lL[0]*lL[1]*lL[2]*NSPINOR,MPI_DOUBLE, &(stypeV[2]));

  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NSPINOR, MPI_DOUBLE, &(stypeV[3]));
  MPI_Type_contiguous(lL[1]*lL[2]*lL[3]*NSPINOR, MPI_DOUBLE, &(rtypeV[0]));
  MPI_Type_contiguous(lL[0]*lL[2]*lL[3]*NSPINOR, MPI_DOUBLE, &(rtypeV[1]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[3]*NSPINOR, MPI_DOUBLE, &(rtypeV[2]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NSPINOR, MPI_DOUBLE, &(rtypeV[3]));


  // for gauge style
  MPI_Type_vector(lL[1]*lL[2]*lL[3], NLINKS, lL[0]*NLINKS,MPI_DOUBLE, &(stypeU[0]));
  MPI_Type_vector(lL[2]*lL[3], NLINKS*lL[0], lL[0]*lL[1]*NLINKS,MPI_DOUBLE, &(stypeU[1]));
  MPI_Type_vector(lL[3], NLINKS*lL[0]*lL[1], lL[0]*lL[1]*lL[2]*NLINKS,MPI_DOUBLE, &(stypeU[2]));

  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NLINKS, MPI_DOUBLE, &(stypeU[3]));
  MPI_Type_contiguous(lL[1]*lL[2]*lL[3]*NLINKS, MPI_DOUBLE, &(rtypeU[0]));
  MPI_Type_contiguous(lL[0]*lL[2]*lL[3]*NLINKS, MPI_DOUBLE, &(rtypeU[1]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[3]*NLINKS, MPI_DOUBLE, &(rtypeU[2]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NLINKS, MPI_DOUBLE, &(rtypeU[3]));


  // for propagator style
  MPI_Type_vector(lL[1]*lL[2]*lL[3], NPROPS, lL[0]*NPROPS,MPI_DOUBLE, &(stypeP[0]));
  MPI_Type_vector(lL[2]*lL[3], NPROPS*lL[0], lL[0]*lL[1]*NPROPS,MPI_DOUBLE, &(stypeP[1]));
  MPI_Type_vector(lL[3], NPROPS*lL[0]*lL[1], lL[0]*lL[1]*lL[2]*NPROPS,MPI_DOUBLE, &(stypeP[2]));

  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NPROPS, MPI_DOUBLE, &(stypeP[3]));
  MPI_Type_contiguous(lL[1]*lL[2]*lL[3]*NPROPS, MPI_DOUBLE, &(rtypeP[0]));
  MPI_Type_contiguous(lL[0]*lL[2]*lL[3]*NPROPS, MPI_DOUBLE, &(rtypeP[1]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[3]*NPROPS, MPI_DOUBLE, &(rtypeP[2]));
  MPI_Type_contiguous(lL[0]*lL[1]*lL[2]*NPROPS, MPI_DOUBLE, &(rtypeP[3]));


  //mpi_commit inside vector,propagator,gauge



  // to create mpi groups


    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++){
	Ccg[mu][nu].real()=0.;
	Ccg[mu][nu].imag()=0.;
	Ccg_bar[mu][nu].real()=0.;
	Ccg_bar[mu][nu].imag()=0.;
	for(int i = 0 ; i < 5 ; i++){
	gammaM[i][mu][nu].real()=0.;
	gammaM[i][mu][nu].imag()=0.;
	CgammaM[i][mu][nu].real()=0.;
	CgammaM[i][mu][nu].imag()=0.;
	CgammaM_bar[i][mu][nu].real()=0.;
	CgammaM_bar[i][mu][nu].imag()=0.;
	}
      }
       
  // 0 -> g1
  gammaM[0][0][3].imag() = 1.;gammaM[0][1][2].imag() = 1.; gammaM[0][2][1].imag() = -1.; gammaM[0][3][0].imag() = -1.;
  // 1 -> g2
  gammaM[1][0][3].real() = 1.;gammaM[1][1][2].real() = -1.; gammaM[1][2][1].real() = -1.; gammaM[1][3][0].real() = 1.;
  // 2 -> g3
  gammaM[2][0][2].imag() = 1.;gammaM[2][1][3].imag() = -1.; gammaM[2][2][0].imag() = -1.; gammaM[2][3][1].imag() = 1.;
  // 3 -> g4
  gammaM[3][0][0].real() = 1.;gammaM[3][1][1].real() = 1.; gammaM[3][2][2].real() = -1.; gammaM[3][3][3].real() = -1.;
  // 4 -> g5
  gammaM[4][0][2].real() = 1.;gammaM[4][1][3].real() = 1.; gammaM[4][2][0].real() = 1.; gammaM[4][3][1].real() = 1.;

  // C = g2g4
  

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++)
      for(int ku = 0 ; ku < 4 ; ku++)
	Ccg[mu][nu] = Ccg[mu][nu] + gammaM[1][mu][ku]*gammaM[3][ku][nu];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++)
      for(int ku = 0 ; ku < 4 ; ku++)
	for(int lu = 0 ; lu < 4 ; lu++)
	  Ccg_bar[mu][nu] = Ccg_bar[mu][nu] + gammaM[3][mu][lu]*conj(Ccg[ku][lu])*gammaM[3][ku][nu];

  for(int i = 0 ; i < 5 ; i++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int ku = 0 ; ku < 4 ; ku++)
	  CgammaM[i][mu][nu] = CgammaM[i][mu][nu] + Ccg[mu][ku]*gammaM[i][ku][nu];

  for(int i = 0 ; i < 5 ; i++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int ku = 0 ; ku < 4 ; ku++)
	  for(int lu = 0 ; lu < 4 ; lu++)
	    CgammaM_bar[i][mu][nu] = CgammaM_bar[i][mu][nu] + gammaM[3][mu][lu]*conj(CgammaM[i][ku][lu])*gammaM[3][ku][nu]; //g4 Cg^\dagger g4 
  /*
#ifdef DEBUG_MODE
  printfQKXTM("Using UKQCD the gamma matrices are the following\n");
  for(int i = 0 ; i < 5 ; i++){
    printfQKXTM("\n");
    printfQKXTM("g%d=\n",i+1);
    for(int mu = 0 ; mu < 4 ; mu++)
      printfQKXTM("(%+f,%+f) \t (%+f,%+f) \t (%+f,%+f) \t (%+f,%+f)\n",gammaM[i][mu][0].real(),gammaM[i][mu][0].imag(),gammaM[i][mu][1].real(),gammaM[i][mu][1].imag(),gammaM[i][mu][2].real(),gammaM[i][mu][2].imag(),gammaM[i][mu][3].real(),gammaM[i][mu][3].imag());
  }
  printfQKXTM("\n");
  printfQKXTM("Th C matrix charge conjucate is\n");
    for(int mu = 0 ; mu < 4 ; mu++)
      printfQKXTM("(%+f,%+f) \t (%+f,%+f) \t (%+f,%+f) \t (%+f,%+f)\n",Ccg[mu][0].real(),Ccg[mu][0].imag(),Ccg[mu][1].real(),Ccg[mu][1].imag(),Ccg[mu][2].real(),Ccg[mu][2].imag(),Ccg[mu][3].real(),Ccg[mu][3].imag());

  printfQKXTM("\n");
  printfQKXTM("Th C matrix charge conjucate bar is\n");
    for(int mu = 0 ; mu < 4 ; mu++)
      printfQKXTM("(%+f,%+f) \t (%+f,%+f) \t (%+f,%+f) \t (%+f,%+f)\n",Ccg_bar[mu][0].real(),Ccg_bar[mu][0].imag(),Ccg_bar[mu][1].real(),Ccg_bar[mu][1].imag(),Ccg_bar[mu][2].real(),Ccg_bar[mu][2].imag(),Ccg_bar[mu][3].real(),Ccg_bar[mu][3].imag());

  printfQKXTM("The C gamma matrix is \n");
  for(int i = 0 ; i < 5 ; i++){
    printfQKXTM("\n");
    printfQKXTM("Cg%d = \n",i+1);
    for(int mu = 0 ; mu < 4 ; mu++)
      printfQKXTM("(%+f,%+f) \t (%+f,%+f) \t (%+f,%+f) \t (%+f,%+f)\n",CgammaM[i][mu][0].real(),CgammaM[i][mu][0].imag(),CgammaM[i][mu][1].real(),CgammaM[i][mu][1].imag(),CgammaM[i][mu][2].real(),CgammaM[i][mu][2].imag(),CgammaM[i][mu][3].real(),CgammaM[i][mu][3].imag());
  }
  printfQKXTM("\n");
  printfQKXTM("The C gamma matrix bar is \n");
  for(int i = 0 ; i < 5 ; i++){
    printfQKXTM("\n");
    printfQKXTM("Cg%d = \n",i+1);
    for(int mu = 0 ; mu < 4 ; mu++)
      printfQKXTM("(%+f,%+f) \t (%+f,%+f) \t (%+f,%+f) \t (%+f,%+f)\n",CgammaM_bar[i][mu][0].real(),CgammaM_bar[i][mu][0].imag(),CgammaM_bar[i][mu][1].real(),CgammaM_bar[i][mu][1].imag(),CgammaM_bar[i][mu][2].real(),CgammaM_bar[i][mu][2].imag(),CgammaM_bar[i][mu][3].real(),CgammaM_bar[i][mu][3].imag());
  }

#endif
  */

  if(initializeGroupsFlag != true){

  MPI_Comm_group(MPI_COMM_WORLD, &fullGroup);
  int space3D_proc;
  space3D_proc = latInfo->P[0] * latInfo->P[1] * latInfo->P[2];
  int *ranks = (int*) malloc(space3D_proc*sizeof(int));

  for(int i= 0 ; i < space3D_proc ; i++)
    ranks[i] = procPos[3] * space3D_proc + i;

  MPI_Group_incl(fullGroup,space3D_proc,ranks,&spaceGroup);
  MPI_Group_rank(spaceGroup,&localRank);
  MPI_Group_size(spaceGroup,&localSize);
  MPI_Comm_create(MPI_COMM_WORLD, spaceGroup , &spaceComm);

  free(ranks);

  int *ranksTime = (int*) malloc(latInfo->P[3]*sizeof(int));
  for(int i=0 ; i < latInfo->P[3] ; i++)
    ranksTime[i] = i*space3D_proc;

  MPI_Group_incl(fullGroup,latInfo->P[3], ranksTime, &timeGroup);
  MPI_Group_rank(timeGroup, &timeRank);
  MPI_Group_size(timeGroup, &timeSize);
  MPI_Comm_create(MPI_COMM_WORLD, timeGroup, &timeComm);

  free(ranksTime);
  }

  initializeGroupsFlag = true;
} // close constructor



Vector::Vector(LatticeInfo *latInfo) : LatticeGeometry(latInfo) {

  size = lV4d ;
  size_total = lV4d;

  for(int i = 0 ; i < NDIM ; i++)
    if(lL[i] < L[i]) size_total += 2*lSurfaceBoundary[i];
  
  M = (Complex(*)[NSPINS][NCOLORS]) safe_malloc(NSPINOR * size_total * sizeof(double));


  for(int dir = 0 ; dir < NDIM ; dir++){
    if(lL[dir] < L[dir]){
      ptr_boundaryPlusStart[dir] = &(M[placeBoundaryStartPlus[dir]][0][0].real());
      ptr_boundaryMinusStart[dir] = &(M[placeBoundaryStartMinus[dir]][0][0].real());
    }
  }

  for(int dir = 0; dir < NDIM; dir++)
    {
      if( lL[dir] < L[dir] )
	{

	  MPI_Type_commit(&(stypeV[dir]));
	  MPI_Type_commit(&(rtypeV[dir]));
	}
    }

    

  initialized = true;
}

Gauge::~Gauge(){
  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }
  free(M);
}                

Propagator::~Propagator(){
  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }
  free(M);
}                


Vector::~Vector(){
  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }
  free(M);
}                

void Vector::zero(){

  memset(&(M[0][0][0].real()) , 0 , NSPINOR*size_total*sizeof(double) );

}

void Vector::copy(Vector &src){

  memcpy(&(this->M[0][0][0].real()) , &(src.M[0][0][0].real()) , lV4d*NSPINOR*sizeof(double) );

}

void Vector::copyPropagator(Propagator &prop,int nu , int c2){
  for(long int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int c1 = 0 ; c1 < 3 ; c1++){
	//      memcpy(&(this->M[0][mu][c1].real()) , &(prop.M[0][mu][nu][c1][c2].real()) , lV4d*2*sizeof(double));
	M[iv][mu][c1] = prop.M[iv][mu][nu][c1][c2];
    }
}

double Vector::norm2(){
  double localResult = 0;
  double globalResult = 0;

  for(long int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
	localResult += norm(M[iv][mu][ic]);
      }

  MPI_Reduce(&localResult,&globalResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return globalResult;

}

void Vector::communicateToPlus(){

  //  MPI_Request request_sent;
  // MPI_Request request_recv;
  int startPos;

  for(int dir = 0 ; dir < NDIM ; dir++){

    startPos = lL[dir] - 1;
    for(int i = 0 ; i < dir ; i++)
      startPos *= lL[i];

    if( lL[dir] < L[dir]){
      
      MPI_Sendrecv(&(M[startPos][0][0].real()), 1 , stypeV[dir] , procPlus[dir] , dir , ptr_boundaryMinusStart[dir] , 1 ,rtypeV[dir] , procMinus[dir] , dir, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
      //      MPI_Isend(&(M[startPos][0][0].real()), 1, stypeV[dir], procPlus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      //      MPI_Irecv(ptr_boundaryMinusStart[dir], 1, rtypeV[dir], procMinus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&(request_recv),MPI_STATUSES_IGNORE);

    }

  } 


}


void Vector::communicateToMinus(){

  //  MPI_Request request_sent;
  // MPI_Request request_recv;


  for(int dir = 0 ; dir < NDIM ; dir++){

    if( lL[dir] < L[dir]){
      MPI_Sendrecv(&(M[0][0][0].real()) , 1 , stypeV[dir] , procMinus[dir] , dir , ptr_boundaryPlusStart[dir] , 1 , rtypeV[dir] , procPlus[dir] , dir , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
      //MPI_Isend(&(M[0][0][0].real()), 1, stypeV[dir], procMinus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      // MPI_Irecv(ptr_boundaryPlusStart[dir], 1, rtypeV[dir], procPlus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&request_recv,MPI_STATUSES_IGNORE);
    }

  } 

 


}

Propagator::Propagator(LatticeInfo *latInfo) : LatticeGeometry(latInfo) {


  size = lV4d ;
  size_total = lV4d;

  for(int i = 0 ; i < NDIM ; i++)
    if(lL[i] < L[i]) size_total += 2*lSurfaceBoundary[i];

  M = (Complex(*)[NSPINS][NSPINS][NCOLORS][NCOLORS]) safe_malloc(NPROPS * size_total * sizeof(double));

  for(int dir = 0 ; dir < NDIM ; dir++){
    if(lL[dir] < L[dir]){
      ptr_boundaryPlusStart[dir] = &(M[placeBoundaryStartPlus[dir]][0][0][0][0].real());
      ptr_boundaryMinusStart[dir] = &(M[placeBoundaryStartMinus[dir]][0][0][0][0].real());
    }
  }

  for(int dir = 0; dir < NDIM; dir++)
    {
      if( lL[dir] < L[dir] )
	{

	  MPI_Type_commit(&(stypeP[dir]));
	  MPI_Type_commit(&(rtypeP[dir]));
	}
    }


  initialized = true;

}

void Propagator::zero(){

  memset(&(M[0][0][0][0][0].real()) , 0 , NPROPS*size_total*sizeof(double) );

}

void Propagator::copy(Propagator &src){

  memcpy(&(this->M[0][0][0][0][0].real()) , &(src.M[0][0][0][0][0].real()) , lV4d*NPROPS*sizeof(double) );

}


void Propagator::copyVector(Vector &v,int nu , int c2){
  for(long int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int c1 = 0 ; c1 < 3 ; c1++){
	//      memcpy(&(this->M[0][mu][c1].real()) , &(prop.M[0][mu][nu][c1][c2].real()) , lV4d*2*sizeof(double));
	M[iv][mu][nu][c1][c2] = v.M[iv][mu][c1];
    }
}

double Propagator::norm2(){
  double localResult = 0;
  double globalResult = 0;

  for(long int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int nu = 0 ; nu < NSPINS ; nu++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	  for(int c2 = 0 ; c2 < NCOLORS ; c2++){
	    localResult += norm(M[iv][mu][nu][c1][c2]);
	  }

  MPI_Reduce(&localResult,&globalResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return globalResult;

}



void Propagator::communicateToPlus(){

  //  MPI_Request request_sent;
  // MPI_Request request_recv;
  int startPos;

  for(int dir = 0 ; dir < NDIM ; dir++){

    startPos = lL[dir] - 1;
    for(int i = 0 ; i < dir ; i++)
      startPos *= lL[i];

    if( lL[dir] < L[dir]){
      MPI_Sendrecv(&(M[startPos][0][0][0][0].real()) , 1 , stypeP[dir] , procPlus[dir] , dir , ptr_boundaryMinusStart[dir] , 1 , rtypeP[dir] , procMinus[dir] , dir, MPI_COMM_WORLD , MPI_STATUSES_IGNORE);
      //      MPI_Isend(&(M[startPos][0][0][0][0].real()), 1, stypeP[dir], procPlus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      // MPI_Irecv(ptr_boundaryMinusStart[dir], 1, rtypeP[dir], procMinus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&request_recv,MPI_STATUSES_IGNORE);
    }

  } 


}


void Propagator::communicateToMinus(){

  //  MPI_Request request_sent;
  // MPI_Request request_recv;

  for(int dir = 0 ; dir < NDIM ; dir++){

    if( lL[dir] < L[dir]){
      MPI_Sendrecv(&(M[0][0][0][0][0].real()) , 1 , stypeP[dir] , procMinus[dir] , dir , ptr_boundaryPlusStart[dir], 1, rtypeP[dir], procPlus[dir], dir, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      //      MPI_Isend(&(M[0][0][0][0][0].real()), 1, stypeP[dir], procMinus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      //MPI_Irecv(ptr_boundaryPlusStart[dir], 1, rtypeP[dir], procPlus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&request_recv,MPI_STATUSES_IGNORE);

    }

  } 


}



Gauge::Gauge(LatticeInfo *latInfo) : LatticeGeometry(latInfo) , commToPlusOnce(false) , commToMinusOnce(false) , constGauge(false) {

  size = lV4d ;
  size_total = lV4d;

  for(int i = 0 ; i < NDIM ; i++)
    if(lL[i] < L[i]) size_total += 2*lSurfaceBoundary[i];

  M = (Complex(*)[NDIM][NCOLORS][NCOLORS]) safe_malloc(NLINKS * size_total * sizeof(double));


  for(int dir = 0 ; dir < NDIM ; dir++){
    if(lL[dir] < L[dir]){
      ptr_boundaryPlusStart[dir] = &(M[placeBoundaryStartPlus[dir]][0][0][0].real());
      ptr_boundaryMinusStart[dir] = &(M[placeBoundaryStartMinus[dir]][0][0][0].real());
    }
  }

  for(int dir = 0; dir < NDIM; dir++)
    {
      if( lL[dir] < L[dir] )
	{
	  MPI_Type_commit(&(stypeU[dir]));
	  MPI_Type_commit(&(rtypeU[dir]));
	}
    }


  initialized = true;
}

void Gauge::zero(){

  memset(&(M[0][0][0][0].real()) , 0 , NLINKS*size_total*sizeof(double) );

}

void Gauge::copy(Gauge &src){

  memcpy(&(this->M[0][0][0][0].real()) , &(src.M[0][0][0][0].real()) , lV4d*NLINKS*sizeof(double) );

}

double Gauge::norm2(){
  double localResult = 0;
  double globalResult = 0;

  for(long int iv = 0 ; iv < lV4d ; iv++)
    for(int dir = 0 ; dir < NDIM ; dir++)
      for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	for(int c2 = 0 ; c2 < NCOLORS ; c2++){
	localResult += norm(M[iv][dir][c1][c2]);
      }

  MPI_Reduce(&localResult,&globalResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return globalResult;

}

void Gauge::applyBoundaryCondition(LatticeInfo *latInfo){
  int sign = latInfo->boundaryT;
  bool last_node_in_t = (procPos[3] == P[3]-1) ? true : false;

  if( (sign == -1) && (last_node_in_t == true) ){
    printfQKXTM("Flip sign in gauge link\n");
    for(long int iv = lL[0]*lL[1]*lL[2]*(lL[3]-1) ; iv < lV4d ; iv++)
      for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	for(int c2 = 0 ; c2 < NCOLORS ; c2++){
	  M[iv][3][c1][c2].real() *= -1;
	  M[iv][3][c1][c2].imag() *= -1;
	}
  }
  else
    return;

}

void Gauge::communicateToPlus(){

  //  MPI_Request request_sent;
  // MPI_Request request_recv;
  int startPos;

  if( (constGauge == true) && (commToPlusOnce == true) ) return;

  for(int dir = 0 ; dir < NDIM ; dir++){

    startPos = lL[dir] - 1;
    for(int i = 0 ; i < dir ; i++)
      startPos *= lL[i];

    if( lL[dir] < L[dir]){
      MPI_Sendrecv(&(M[startPos][0][0][0].real()), 1, stypeU[dir], procPlus[dir] , dir,ptr_boundaryMinusStart[dir], 1, rtypeU[dir], procMinus[dir], dir, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      //  MPI_Isend(&(M[startPos][0][0][0].real()), 1, stypeU[dir], procPlus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      //MPI_Irecv(ptr_boundaryMinusStart[dir], 1, rtypeU[dir], procMinus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&request_recv,MPI_STATUSES_IGNORE);
    }

  } 

  commToPlusOnce = true;

}


void Gauge::communicateToMinus(){

  // MPI_Request request_sent;
  //  MPI_Request request_recv;

  if( (constGauge == true) && (commToMinusOnce == true) ) return;

  for(int dir = 0 ; dir < NDIM ; dir++){

    if( lL[dir] < L[dir]){
      MPI_Sendrecv(&(M[0][0][0][0].real()), 1, stypeU[dir], procMinus[dir] , dir, ptr_boundaryPlusStart[dir], 1, rtypeU[dir], procPlus[dir], dir, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      //  MPI_Isend(&(M[0][0][0][0].real()), 1, stypeU[dir], procMinus[dir] , dir, MPI_COMM_WORLD, &(request_sent) );
      //MPI_Irecv(ptr_boundaryPlusStart[dir], 1, rtypeU[dir], procPlus[dir], dir, MPI_COMM_WORLD, &(request_recv));
      // MPI_Wait(&request_recv,MPI_STATUSES_IGNORE);
    }

  } 

  commToMinusOnce = true;
}

/*
double Gauge::calculatePlaquette(){

  // first we need the point infront of us so we need communication
    communicateToMinus(); // 
  //    communicateToPlus(); // 

  Complex *a = 0;
  Complex *b = 0;
  Complex *c = 0;
  Complex *d = 0;
  
  register Complex prod1[NCOLORS][NCOLORS][NCOLORS][NCOLORS];
  register Complex prod2[NCOLORS][NCOLORS][NCOLORS][NCOLORS];
  register double plaq = 0.;

  double fullPlaq = 0.;

  //  for(int it = 0 ; it < lL[3] ; it++)
  // for(int iz =0 ; iz < lL[2] ; iz++)
  //   for(int iy = 0 ; iy < lL[1] ; iy++)
  //	for(int ix = 0 ; ix < lL[0] ; ix++){
  
  //	  int iv = LEXIC(it,iz,iy,ix,lL);
  for(int iv = 0 ; iv < lV4d ; iv++){
    for(int mu = 0 ; mu < NSPINS -1 ; mu++)
	    for(int nu = mu+1 ; nu < NSPINS ; nu++){
	      a = &(M[iv][mu][0][0]);
	      d = &(M[iv][nu][0][0]);
	      b = &(M[pointPlus[mu][iv]][nu][0][0]);
	      c = &(M[pointPlus[nu][iv]][mu][0][0]);


	      for(int c1 = 0 ; c1 < NCOLORS ; c1++)
		for(int c2 = 0 ; c2 < NCOLORS ; c2++)
		  for(int c3 = 0 ; c3 < NCOLORS ; c3++)
		    for(int c4 = 0 ; c4 < NCOLORS ; c4++){
		      prod1[c1][c2][c3][c4] = a[c1*NCOLORS+c2] * conj(d[c3*NCOLORS+c4]);
		    }



	      for(int c1 = 0 ; c1 < NCOLORS ; c1++)
		for(int c2 = 0 ; c2 < NCOLORS ; c2++)
		  for(int c3 = 0 ; c3 < NCOLORS ; c3++)
		    for(int c4 = 0 ; c4 < NCOLORS ; c4++){
		      prod2[c1][c2][c3][c4] = b[c1*NCOLORS+c2] * conj(c[c3*NCOLORS+c4]);
		    }

	      plaq += (prod1[0][0][0][0] * prod2[0][0][0][0]).real();
	      plaq += (prod1[0][1][0][0] * prod2[1][0][0][0]).real();
	      plaq += (prod1[0][2][0][0] * prod2[2][0][0][0]).real();
	      plaq += (prod1[0][0][0][0] * prod2[0][1][0][1]).real();
	      plaq += (prod1[0][1][0][0] * prod2[1][1][0][1]).real();
	      plaq += (prod1[0][2][0][0] * prod2[2][1][0][1]).real();
	      plaq += (prod1[0][0][0][0] * prod2[0][2][0][2]).real();
	      plaq += (prod1[0][1][0][0] * prod2[1][2][0][2]).real();
	      plaq += (prod1[0][2][0][0] * prod2[2][2][0][2]).real();
	      plaq += (prod1[0][0][0][1] * prod2[0][0][1][0]).real();
	      plaq += (prod1[0][1][0][1] * prod2[1][0][1][0]).real();
	      plaq += (prod1[0][2][0][1] * prod2[2][0][1][0]).real();

	      plaq += (prod1[0][0][0][1] * prod2[0][1][1][1]).real();
	      plaq += (prod1[0][1][0][1] * prod2[1][1][1][1]).real();
	      plaq += (prod1[0][2][0][1] * prod2[2][1][1][1]).real();
	      plaq += (prod1[0][0][0][1] * prod2[0][2][1][2]).real();
	      plaq += (prod1[0][1][0][1] * prod2[1][2][1][2]).real();
	      plaq += (prod1[0][2][0][1] * prod2[2][2][1][2]).real();
	      plaq += (prod1[0][0][0][2] * prod2[0][0][2][0]).real();
	      plaq += (prod1[0][1][0][2] * prod2[1][0][2][0]).real();
	      plaq += (prod1[0][2][0][2] * prod2[2][0][2][0]).real();
	      plaq += (prod1[0][0][0][2] * prod2[0][1][2][1]).real();
	      plaq += (prod1[0][1][0][2] * prod2[1][1][2][1]).real();
	      plaq += (prod1[0][2][0][2] * prod2[2][1][2][1]).real();

	      plaq += (prod1[0][0][0][2] * prod2[0][2][2][2]).real();
	      plaq += (prod1[0][1][0][2] * prod2[1][2][2][2]).real();
	      plaq += (prod1[0][2][0][2] * prod2[2][2][2][2]).real();
	      plaq += (prod1[1][0][1][0] * prod2[0][0][0][0]).real();
	      plaq += (prod1[1][1][1][0] * prod2[1][0][0][0]).real();
	      plaq += (prod1[1][2][1][0] * prod2[2][0][0][0]).real();
	      plaq += (prod1[1][0][1][0] * prod2[0][1][0][1]).real();
	      plaq += (prod1[1][1][1][0] * prod2[1][1][0][1]).real();
	      plaq += (prod1[1][2][1][0] * prod2[2][1][0][1]).real();
	      plaq += (prod1[1][0][1][0] * prod2[0][2][0][2]).real();
	      plaq += (prod1[1][1][1][0] * prod2[1][2][0][2]).real();
	      plaq += (prod1[1][2][1][0] * prod2[2][2][0][2]).real();

	      plaq += (prod1[1][0][1][1] * prod2[0][0][1][0]).real();
	      plaq += (prod1[1][1][1][1] * prod2[1][0][1][0]).real();
	      plaq += (prod1[1][2][1][1] * prod2[2][0][1][0]).real();
	      plaq += (prod1[1][0][1][1] * prod2[0][1][1][1]).real();
	      plaq += (prod1[1][1][1][1] * prod2[1][1][1][1]).real();
	      plaq += (prod1[1][2][1][1] * prod2[2][1][1][1]).real();
	      plaq += (prod1[1][0][1][1] * prod2[0][2][1][2]).real();
	      plaq += (prod1[1][1][1][1] * prod2[1][2][1][2]).real();
	      plaq += (prod1[1][2][1][1] * prod2[2][2][1][2]).real();
	      plaq += (prod1[1][0][1][2] * prod2[0][0][2][0]).real();
	      plaq += (prod1[1][1][1][2] * prod2[1][0][2][0]).real();
	      
	      plaq += (prod1[1][2][1][2] * prod2[2][0][2][0]).real();
	      plaq += (prod1[1][0][1][2] * prod2[0][1][2][1]).real();
	      plaq += (prod1[1][1][1][2] * prod2[1][1][2][1]).real();
	      plaq += (prod1[1][2][1][2] * prod2[2][1][2][1]).real();
	      plaq += (prod1[1][0][1][2] * prod2[0][2][2][2]).real();
	      plaq += (prod1[1][1][1][2] * prod2[1][2][2][2]).real();
	      plaq += (prod1[1][2][1][2] * prod2[2][2][2][2]).real();
	      plaq += (prod1[2][0][2][0] * prod2[0][0][0][0]).real();
	      plaq += (prod1[2][1][2][0] * prod2[1][0][0][0]).real();
	      plaq += (prod1[2][2][2][0] * prod2[2][0][0][0]).real();
	      plaq += (prod1[2][0][2][0] * prod2[0][1][0][1]).real();

	      plaq += (prod1[2][1][2][0] * prod2[1][1][0][1]).real();
	      plaq += (prod1[2][2][2][0] * prod2[2][1][0][1]).real();
	      plaq += (prod1[2][0][2][0] * prod2[0][2][0][2]).real();
	      plaq += (prod1[2][1][2][0] * prod2[1][2][0][2]).real();
	      plaq += (prod1[2][2][2][0] * prod2[2][2][0][2]).real();
	      plaq += (prod1[2][0][2][1] * prod2[0][0][1][0]).real();
	      plaq += (prod1[2][1][2][1] * prod2[1][0][1][0]).real();
	      plaq += (prod1[2][2][2][1] * prod2[2][0][1][0]).real();
	      plaq += (prod1[2][0][2][1] * prod2[0][1][1][1]).real();
	      plaq += (prod1[2][1][2][1] * prod2[1][1][1][1]).real();
	      plaq += (prod1[2][2][2][1] * prod2[2][1][1][1]).real();

	      plaq += (prod1[2][0][2][1] * prod2[0][2][1][2]).real();
	      plaq += (prod1[2][1][2][1] * prod2[1][2][1][2]).real();
	      plaq += (prod1[2][2][2][1] * prod2[2][2][1][2]).real();
	      plaq += (prod1[2][0][2][2] * prod2[0][0][2][0]).real();
	      plaq += (prod1[2][1][2][2] * prod2[1][0][2][0]).real();
	      plaq += (prod1[2][2][2][2] * prod2[2][0][2][0]).real();
	      plaq += (prod1[2][0][2][2] * prod2[0][1][2][1]).real();
	      plaq += (prod1[2][1][2][2] * prod2[1][1][2][1]).real();
	      plaq += (prod1[2][2][2][2] * prod2[2][1][2][1]).real();
	      plaq += (prod1[2][0][2][2] * prod2[0][2][2][2]).real();
	      plaq += (prod1[2][1][2][2] * prod2[1][2][2][2]).real();
	      plaq += (prod1[2][2][2][2] * prod2[2][2][2][2]).real();
	    } // close mu nu

	} // close space


  MPI_Reduce(&plaq,&fullPlaq,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return fullPlaq / (6*NCOLORS*V4d);

}
*/


double Gauge::calculatePlaquette(){ 

  communicateToMinus();
  Complex plaq[3][3],tmp[3][3];
  double meanplaq=0;
  double result=0;

  for(long int iv = 0 ; iv < lV4d ; iv++){                                                                                                                                     
    for(int mu = 0 ; mu < NSPINS -1 ; mu++)                                                                                                                                 
      for(int nu = mu+1 ; nu < NSPINS ; nu++){       

	qcd_MUL3x3(plaq,M[iv][mu],M[pointPlus[mu][iv]][nu]);
	qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[nu][iv]][mu]);
	qcd_MULADJOINT3x3(plaq,tmp,M[iv][nu]);
	meanplaq += qcd_SU3TRACER(plaq);
      }

  }
    MPI_Allreduce(&meanplaq, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return(result/(V4d * 3.0 * 6.0)); // volume * N_color * number of mu-nu combinations            



}


void Gauge::APE_smearing(Gauge &U_tmp, LatticeInfo *latInfo){

  int nsmear = latInfo->NsmearAPE;
  double alpha = latInfo->alphaAPE;
  Complex stapleForward[3][3];
  Complex stapleTmp[3][3];
  if(nsmear == 0) return;
  
  //  Gauge *null_object = NULL;

  Propagator *prop = new Propagator(latInfo);
  Propagator &prp = *prop;
  
  Gauge *in = NULL;
  Gauge *out = NULL;
  Gauge *tmp = NULL;
  
  in = &(U_tmp);
  out = &(*this);

  printfQKXTM("Perform APE smearing\n");  
  for(int iter = 0 ; iter < nsmear ; iter++){

    in->communicateToMinus(); // so has the plus points


    for(int iv = 0 ; iv < lV4d ; iv++)
      for(int mu = 0 ; mu < NSPINS -1 ; mu++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	  for(int c2 = 0 ; c2 < NCOLORS ; c2++){
	    out->M[iv][mu][c1][c2].real() = 0.; out->M[iv][mu][c1][c2].imag() = 0.;
	  }

    for(int mu = 0 ; mu < NDIM -1 ; mu++)
      for(int nu = 0 ; nu < NDIM -1 ; nu++)
	if( nu != mu)
	  for(int iv = 0 ; iv < lV4d ; iv++){

	    qcd_MUL3x3(stapleTmp, in->M[iv][mu], in->M[pointPlus[mu][iv]][nu]);
	    qcd_ADJOINTMUL3x3(prp.M[iv][mu][nu],in->M[iv][nu],stapleTmp);

	    /*
	    if( (iv == 0) && (mu == 0) && (nu==1) ){
	      for(int c1 = 0 ; c1 < 3 ; c1++)
		for(int c2 = 0 ; c2 < 3 ; c2++)
		  printf("%e %e\n",prp.M[iv][mu][nu][c1][c2].real(), prp.M[iv][mu][nu][c1][c2].imag());
	    }
	    */
	    qcd_MUL3x3(stapleTmp, in->M[iv][nu],in->M[pointPlus[nu][iv]][mu]);
	    qcd_MULADJOINT3x3(stapleForward, stapleTmp, in->M[pointPlus[mu][iv]][nu]);

	    for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	      for(int c2 = 0 ; c2 < NCOLORS ; c2++)
		out->M[iv][mu][c1][c2] = out->M[iv][mu][c1][c2] + stapleForward[c1][c2];
	    /*
	    if( (iv == 0) && (mu == 1) && (nu==2) ){
	      for(int c1 = 0 ; c1 < 3 ; c1++)
		for(int c2 = 0 ; c2 < 3 ; c2++)
		  printf("%e %e\n",out->M[iv][mu][c1][c2].real(), out->M[iv][mu][c1][c2].imag());                                                                        
	    }                
	    */

	  }
    
    /*
    for(int mu = 0 ; mu < 3 ; mu++)
    for(int c1 = 0 ; c1 < 3 ; c1++)
      for(int c2 = 0 ; c2 < 3 ; c2++)
	printf("%e %e\n",out->M[0][mu][c1][c2].real(), out->M[0][mu][c1][c2].imag());                                                                        
    */
    /*
    for(int mu = 0 ; mu < 3 ; mu++)
      for(int nu = 0 ; nu < 3 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3 ; c2++)
	    printf("%e %e\n",prp.M[0][mu][nu][c1][c2].real(), prp.M[0][mu][nu][c1][c2].imag());                                                                        
    */
    
    prp.communicateToPlus();
    
    for(int mu = 0 ; mu < NDIM -1 ; mu++)
      for(int nu = 0 ; nu < NDIM -1 ; nu++)
        if( nu != mu)
          for(int iv = 0 ; iv < lV4d ; iv++)
	    for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	      for(int c2 = 0 ; c2 < NCOLORS ; c2++){
		out->M[iv][mu][c1][c2] = out->M[iv][mu][c1][c2] + prp.M[pointMinus[nu][iv]][mu][nu][c1][c2];
		
		//		if(iv == 0 && mu == 0 && nu == 2){
		  //printf("%e %e\n",out->M[0][mu][c1][c2].real(), out->M[0][mu][c1][c2].imag());                                                        
		// printf("%e %e\n",prp.M[pointMinus[nu][iv]][mu][nu][c1][c2].real(), prp.M[pointMinus[nu][iv]][mu][nu][c1][c2].imag());                                                        
		//	}

	      }

    /*
    for(int dir = 0 ; dir < 1 ; dir++)
      for(int c1 = 0 ; c1 < 3 ; c1++)
	for(int c2 = 0 ; c2 < 3 ; c2++)
	  printf("%e %e\n",out->M[0][dir][c1][c2].real(), out->M[0][dir][c1][c2].imag());                                                        
    */

    for(int iv = 0 ; iv < lV4d ; iv++)
      for(int mu = 0 ; mu < NSPINS -1 ; mu++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++)
	  for(int c2 = 0 ; c2 < NCOLORS ; c2++){
	    out->M[iv][mu][c1][c2] = alpha*out->M[iv][mu][c1][c2] + in->M[iv][mu][c1][c2];
	  }

    qcd_projectSU33d(*out);

    tmp=in;
    in=out;
    out=tmp;

    
  } // close iteration

  if( (nsmear%2) == 0){
    memcpy( &(this->M[0][0][0][0].real()) , &(U_tmp.M[0][0][0][0].real()) , lV4d  * NLINKS * sizeof(double) );
  }
  else
    {
      for(int iv = 0 ; iv < lV4d ; iv++)
	memcpy(&(this->M[iv][3][0][0].real()) , &(U_tmp.M[iv][3][0][0].real()) , NCOLORS*NCOLORS*NDF * sizeof(double) );
    }

  delete prop;
}

void Vector::Gaussian_smearing(Vector &vec_tmp, Gauge &gauge, LatticeInfo *latInfo){

  int nsmear = latInfo->NsmearGauss;
  double alpha = latInfo->alphaGauss;
  double normalize = 1./(1 + 6*alpha);

  double *uu;
  double *psi;
  double upsi[24];
  double udaggerpsi[24];
  double *total;

  if(nsmear == 0)return;

  //communicate gauge field 
  gauge.communicateToPlus();

  Vector *in = NULL;
  Vector *out = NULL;
  Vector *tmp = NULL;

  in = &(vec_tmp);
  out = &(*this);

  printfQKXTM("Perform Gaussian smearing\n");
  for(int iter = 0 ; iter < nsmear ; iter++){


    in->communicateToPlus();      // for every iteration we need communication
    in->communicateToMinus();
    out->zero();

    for(int iv = 0 ; iv < lV4d ; iv++){
      
      for(int mu = 0 ; mu < NDIM -1 ; mu++){
	
	uu = &(gauge.M[iv][mu][0][0].real());
	psi = &(in->M[pointPlus[mu][iv]][0][0].real());
	qcd_APPLY_U(uu,upsi,psi);

	uu = &(gauge.M[pointMinus[mu][iv]][mu][0][0].real());
        psi = &(in->M[pointMinus[mu][iv]][0][0].real());
	qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

	total = &(out->M[iv][0][0].real());
	qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
      }
      
    } // close local volume 

    for(int iv = 0 ; iv < lV4d ; iv++)
      for(int mu = 0 ; mu < NSPINS ; mu++)
	for(int ic = 0 ; ic < NCOLORS ; ic++){
	  out->M[iv][mu][ic] = normalize*(in->M[iv][mu][ic] + alpha*out->M[iv][mu][ic]);
	}

    tmp = in;
    in = out;
    out = tmp;

  } // smearing iteration

  if( (nsmear%2) == 0){
    memcpy( &(this->M[0][0][0].real()),&(in->M[0][0][0]), lV4d*NSPINOR*sizeof(double) );
  }   

}

static int getOddBit(int latt_cord,int nx, int ny, int nz){ // return 1 if is odd and 0 if is even
  int x,y,z,t;
  t = latt_cord/(nz*ny*nx);
  z = latt_cord/(ny*nx) - t*nz;
  y = latt_cord/nx - t*nz*ny - z*ny;
  x = latt_cord - t*nz*ny*nx - z*ny*nx - y*nx;
  return (x+y+z+t) & 1;
}


void Vector::applyDslash(Vector &psi, Gauge &u){

  Complex phi[NSPINS*NCOLORS];
  Complex R[NSPINS*NCOLORS];
  Complex xi[NSPINS*NCOLORS];

  // ATTENTION !!!!! remember to communicate Gauge to Plus outside

  psi.communicateToPlus();
  psi.communicateToMinus();
  u.communicateToPlus();

  for(int iv = 0 ; iv < lV4d ; iv++){

    ZERO(R);
    //plus X direction                                                                                                                                                      
    PROJ_MINUS_X(phi,psi);
    APPLY_LINK(xi,u,phi,0);
    COLLECT_MINUS_X(R,xi);
    //minus X direction                                                                                                                                                      
    PROJ_PLUS_X(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,0);
    COLLECT_PLUS_X(R,xi);
    //plus Y direction                                                                                                                                                       
    PROJ_MINUS_Y(phi,psi);
    APPLY_LINK(xi,u,phi,1);
    COLLECT_MINUS_Y(R,xi);
    //minus Y direction                                                                                                                                                     
    PROJ_PLUS_Y(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,1);
    COLLECT_PLUS_Y(R,xi);
    //plus Z direction                                                                                                                                                        
    PROJ_MINUS_Z(phi,psi);
    APPLY_LINK(xi,u,phi,2);
    COLLECT_MINUS_Z(R,xi);
    //minus Z direction                                                                                                                                         
    PROJ_PLUS_Z(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,2);
    COLLECT_PLUS_Z(R,xi);
    //plus T direction                                                                                                                                 
    PROJ_MINUS_T(phi,psi);
    APPLY_LINK(xi,u,phi,3);
    COLLECT_MINUS_T(R,xi);
    //minus T direction                                                                                                                           
    PROJ_PLUS_T(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,3);
    COLLECT_PLUS_T(R,xi);
    //////////////////////////////////////
    // save data to memory
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
        this->M[iv][mu][ic] = R[MS(mu,ic)];

  } // local volume loop

}

void Vector::applyDslash_chiral(Vector &psi, Gauge &u){

  Complex phi[NSPINS*NCOLORS];
  Complex R[NSPINS*NCOLORS];
  Complex xi[NSPINS*NCOLORS];

  // ATTENTION !!!!! remember to communicate Gauge to Plus outside

  psi.communicateToPlus();
  psi.communicateToMinus();
  u.communicateToPlus();

  for(int iv = 0 ; iv < lV4d ; iv++){

    ZERO(R);
    //plus X direction                                                                                                                                                      
    PROJ_MINUS_X_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,0);
    COLLECT_MINUS_X_CHIRAL(R,xi);
    //minus X direction                                                                                                                                                      
    PROJ_PLUS_X_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,0);
    COLLECT_PLUS_X_CHIRAL(R,xi);
    //plus Y direction                                                                                                                                                       
    PROJ_MINUS_Y_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,1);
    COLLECT_MINUS_Y_CHIRAL(R,xi);
    //minus Y direction                                                                                                                                                     
    PROJ_PLUS_Y_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,1);
    COLLECT_PLUS_Y_CHIRAL(R,xi);
    //plus Z direction                                                                                                                                                        
    PROJ_MINUS_Z_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,2);
    COLLECT_MINUS_Z_CHIRAL(R,xi);
    //minus Z direction                                                                                                                                         
    PROJ_PLUS_Z_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,2);
    COLLECT_PLUS_Z_CHIRAL(R,xi);
    //plus T direction                                                                                                                                 
    PROJ_MINUS_T_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,3);
    COLLECT_MINUS_T_CHIRAL(R,xi);
    //minus T direction                                                                                                                           
    PROJ_PLUS_T_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,3);
    COLLECT_PLUS_T_CHIRAL(R,xi);
    //////////////////////////////////////
    // save data to memory
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
        this->M[iv][mu][ic] = R[MS(mu,ic)];

  } // local volume loop

}

void Vector::applyDslash_xx(Vector &psi, Gauge &u,PARITY parity){

  Complex phi[NSPINS*NCOLORS];
  Complex R[NSPINS*NCOLORS];
  Complex xi[NSPINS*NCOLORS];

  // ATTENTION !!!!! remember to communicate Gauge to Plus outside

  psi.communicateToPlus();
  psi.communicateToMinus();
  u.communicateToPlus();

  for(int iv = 0 ; iv < lV4d ; iv++){

    int oddBit = getOddBit(iv,lL[0],lL[1],lL[2]);
    if(parity == OE){
      if(!oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else if(parity == EO){
      if(oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else{
      errorQKXTM("Error wrong parity for hopping term\n");
    }

    ZERO(R);
    //plus X direction                                                                                                                                                      
    PROJ_MINUS_X(phi,psi);
    APPLY_LINK(xi,u,phi,0);
    COLLECT_MINUS_X(R,xi);
    //minus X direction                                                                                                                                                      
    PROJ_PLUS_X(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,0);
    COLLECT_PLUS_X(R,xi);
    //plus Y direction                                                                                                                                                       
    PROJ_MINUS_Y(phi,psi);
    APPLY_LINK(xi,u,phi,1);
    COLLECT_MINUS_Y(R,xi);
    //minus Y direction                                                                                                                                                     
    PROJ_PLUS_Y(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,1);
    COLLECT_PLUS_Y(R,xi);
    //plus Z direction                                                                                                                                                        
    PROJ_MINUS_Z(phi,psi);
    APPLY_LINK(xi,u,phi,2);
    COLLECT_MINUS_Z(R,xi);
    //minus Z direction                                                                                                                                         
    PROJ_PLUS_Z(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,2);
    COLLECT_PLUS_Z(R,xi);
    //plus T direction                                                                                                                                 
    PROJ_MINUS_T(phi,psi);
    APPLY_LINK(xi,u,phi,3);
    COLLECT_MINUS_T(R,xi);
    //minus T direction                                                                                                                           
    PROJ_PLUS_T(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,3);
    COLLECT_PLUS_T(R,xi);
    //////////////////////////////////////
    // save data to memory
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
        this->M[iv][mu][ic] = R[MS(mu,ic)];

  } // local volume loop

}

void Vector::applyDslash_xx_chiral(Vector &psi, Gauge &u,PARITY parity){

  Complex phi[NSPINS*NCOLORS];
  Complex R[NSPINS*NCOLORS];
  Complex xi[NSPINS*NCOLORS];

  // ATTENTION !!!!! remember to communicate Gauge to Plus outside

  psi.communicateToPlus();
  psi.communicateToMinus();
  u.communicateToPlus();

  for(int iv = 0 ; iv < lV4d ; iv++){

    int oddBit = getOddBit(iv,lL[0],lL[1],lL[2]);
    if(parity == OE){
      if(!oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else if(parity == EO){
      if(oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else{
      errorQKXTM("Error wrong parity for hopping term\n");
    }

    ZERO(R);
    //plus X direction                                                                                                                                                      
    PROJ_MINUS_X_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,0);
    COLLECT_MINUS_X_CHIRAL(R,xi);
    //minus X direction                                                                                                                                                      
    PROJ_PLUS_X_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,0);
    COLLECT_PLUS_X_CHIRAL(R,xi);
    //plus Y direction                                                                                                                                                       
    PROJ_MINUS_Y_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,1);
    COLLECT_MINUS_Y_CHIRAL(R,xi);
    //minus Y direction                                                                                                                                                     
    PROJ_PLUS_Y_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,1);
    COLLECT_PLUS_Y_CHIRAL(R,xi);
    //plus Z direction                                                                                                                                                        
    PROJ_MINUS_Z_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,2);
    COLLECT_MINUS_Z_CHIRAL(R,xi);
    //minus Z direction                                                                                                                                         
    PROJ_PLUS_Z_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,2);
    COLLECT_PLUS_Z_CHIRAL(R,xi);
    //plus T direction                                                                                                                                 
    PROJ_MINUS_T_CHIRAL(phi,psi);
    APPLY_LINK(xi,u,phi,3);
    COLLECT_MINUS_T_CHIRAL(R,xi);
    //minus T direction                                                                                                                           
    PROJ_PLUS_T_CHIRAL(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,3);
    COLLECT_PLUS_T_CHIRAL(R,xi);
    //////////////////////////////////////
    // save data to memory
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
        this->M[iv][mu][ic] = R[MS(mu,ic)];

  } // local volume loop

}

void Vector::applyTwistAddDslash(Vector &dslash, Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = (-latInfo->kappa)*dslash.M[iv][mu][ic] + psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}

void Vector::applyTwistAddDslash_chiral(Vector &dslash, Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][0][ic];
      phi[MS(1,ic)] = psi.M[iv][1][ic];
      phi[MS(2,ic)] = -psi.M[iv][2][ic];
      phi[MS(3,ic)] = -psi.M[iv][3][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = (-latInfo->kappa)*dslash.M[iv][mu][ic] + psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}

void Vector::applyTwistAddDslash_PC(Vector &dslash_PC, Vector &psi, LatticeInfo *latInfo, PARITY parity){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){

    int oddBit = getOddBit(iv,lL[0],lL[1],lL[2]);
    if(parity == OO){
      if(!oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else if(parity == EE){
      if(oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else{
      errorQKXTM("Error wrong parity for hopping term\n");
    }

    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = (-latInfo->kappa*latInfo->kappa)*dslash_PC.M[iv][mu][ic] + (psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)]);
      }

  } // close volume loop                                     

}

void Vector::applyTwistAddDslash_PC_chiral(Vector &dslash_PC, Vector &psi, LatticeInfo *latInfo, PARITY parity){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){

    int oddBit = getOddBit(iv,lL[0],lL[1],lL[2]);
    if(parity == OO){
      if(!oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else if(parity == EE){
      if(oddBit){
	memset(&(this->M[iv][0][0]),0,4*3*sizeof(Complex));
	continue;
      }
    }
    else{
      errorQKXTM("Error wrong parity for hopping term\n");
    }

    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][0][ic];
      phi[MS(1,ic)] = psi.M[iv][1][ic];
      phi[MS(2,ic)] = -psi.M[iv][2][ic];
      phi[MS(3,ic)] = -psi.M[iv][3][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = (-latInfo->kappa*latInfo->kappa)*dslash_PC.M[iv][mu][ic] + (psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)]);
      }

  } // close volume loop                                     

}

void Vector::applyDeltaAddDslash(Vector &dslash, Vector &psi, LatticeInfo *latInfo){

  for(int iv = 0 ; iv < lV4d; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic =0 ; ic < 3 ; ic++)
	this->M[iv][mu][ic] = psi.M[iv][mu][ic] - (latInfo->kappa)*dslash.M[iv][mu][ic];

}



void Vector::applyTwist(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] =  psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}


void Vector::applyTwistInv(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);
  double normalize = (1 + (2*latInfo->kappa*latInfo->mu)*(2*latInfo->kappa*latInfo->mu) );

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] =  psi.M[iv][mu][ic] - twistCoeff*phi[MS(mu,ic)];
	this->M[iv][mu][ic].real() /= normalize; 
	this->M[iv][mu][ic].imag() /= normalize; 
      }

  } // close volume loop                                     

}

void Vector::applyTwistInv_chiral(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = 2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);
  double normalize = (1 + (2*latInfo->kappa*latInfo->mu)*(2*latInfo->kappa*latInfo->mu) );

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][0][ic];
      phi[MS(1,ic)] = psi.M[iv][1][ic];
      phi[MS(2,ic)] = -psi.M[iv][2][ic];
      phi[MS(3,ic)] = -psi.M[iv][3][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] =  psi.M[iv][mu][ic] - twistCoeff*phi[MS(mu,ic)];
	this->M[iv][mu][ic].real() /= normalize; 
	this->M[iv][mu][ic].imag() /= normalize; 
      }

  } // close volume loop                                     

}


void Vector::applyDslashDag(Vector &psi, Gauge &u){

  Complex phi[NSPINS*NCOLORS];
  Complex R[NSPINS*NCOLORS];
  Complex xi[NSPINS*NCOLORS];

  // ATTENTION !!!!! remember to communicate Gauge to Plus outside

  psi.communicateToPlus();
  psi.communicateToMinus();
  u.communicateToPlus();

  for(int iv = 0 ; iv < lV4d ; iv++){

    ZERO(R);
    //plus X direction                                                                                                                                                      
    DAG_PROJ_PLUS_X(phi,psi);
    APPLY_LINK(xi,u,phi,0);
    COLLECT_PLUS_X(R,xi);
    //minus X direction                                                                                                                                                      
    DAG_PROJ_MINUS_X(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,0);
    COLLECT_MINUS_X(R,xi);
    //plus Y direction                                                                                                                                                       
    DAG_PROJ_PLUS_Y(phi,psi);
    APPLY_LINK(xi,u,phi,1);
    COLLECT_PLUS_Y(R,xi);
    //minus Y direction                                                                                                                                                     
    DAG_PROJ_MINUS_Y(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,1);
    COLLECT_MINUS_Y(R,xi);
    //plus Z direction                                                                                                                                                        
    DAG_PROJ_PLUS_Z(phi,psi);
    APPLY_LINK(xi,u,phi,2);
    COLLECT_PLUS_Z(R,xi);
    //minus Z direction                                                                                                                                         
    DAG_PROJ_MINUS_Z(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,2);
    COLLECT_MINUS_Z(R,xi);
    //plus T direction                                                                                                                                 
    DAG_PROJ_PLUS_T(phi,psi);
    APPLY_LINK(xi,u,phi,3);
    COLLECT_PLUS_T(R,xi);
    //minus T direction                                                                                                                           
    DAG_PROJ_MINUS_T(phi,psi);
    APPLY_LINK_DAG(xi,u,phi,3);
    COLLECT_MINUS_T(R,xi);
    //////////////////////////////////////
    // save data to memory
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
        this->M[iv][mu][ic] = R[MS(mu,ic)];

  } // local volume loop

}


void Vector::applyTwistDagAddDslashDag(Vector &dslashDag, Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = -2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = (-latInfo->kappa)*dslashDag.M[iv][mu][ic] + psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}


void Vector::applyTwistDag(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = -2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = psi.M[iv][mu][ic] + twistCoeff*phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}

void Vector::applyTwistDagInv(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];
  Complex twistCoeff;
  twistCoeff.real() = 0.;
  twistCoeff.imag() = -2*(latInfo->kappa)*(latInfo->mu)*(latInfo->twistSign);
  double normalize = (1 + (2*latInfo->kappa*latInfo->mu)*(2*latInfo->kappa*latInfo->mu) );

  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] =  psi.M[iv][mu][ic] - twistCoeff*phi[MS(mu,ic)];
	this->M[iv][mu][ic].real() /= normalize; 
	this->M[iv][mu][ic].imag() /= normalize; 
      }

  } // close volume loop                                     

}


void Vector::applyGamma5(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];


  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][2][ic];
      phi[MS(1,ic)] = psi.M[iv][3][ic];
      phi[MS(2,ic)] = psi.M[iv][0][ic];
      phi[MS(3,ic)] = psi.M[iv][1][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}

void Vector::applyGamma5_chiral(Vector &psi, LatticeInfo *latInfo){
  
  Complex phi[NSPINS*NCOLORS];


  for(int iv = 0 ; iv < lV4d; iv++){
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      phi[MS(0,ic)] = psi.M[iv][0][ic];
      phi[MS(1,ic)] = psi.M[iv][1][ic];
      phi[MS(2,ic)] = -psi.M[iv][2][ic];
      phi[MS(3,ic)] = -psi.M[iv][3][ic];
    } // gamma_5 rotation

    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
        this->M[iv][mu][ic] = phi[MS(mu,ic)];
      }

  } // close volume loop                                     

}

void Vector::daggerVectorGamma5(Vector &psi, LatticeInfo *latInfo){

  for(int iv = 0 ; iv < lV4d; iv++)
    for(int ic = 0 ; ic < NCOLORS ; ic++){
      M[iv][0][ic] = conj(psi.M[iv][2][ic]);
      M[iv][1][ic] = conj(psi.M[iv][3][ic]);
      M[iv][2][ic] = conj(psi.M[iv][0][ic]);
      M[iv][3][ic] = conj(psi.M[iv][1][ic]);      
    } 
  
}

void Vector::rotateFromChiralToUKQCD(){
  Complex transMatrix[4][4];
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      transMatrix[mu][nu].real();
      transMatrix[mu][nu].imag();
    }
  double value = 1./sqrt(2.);

  transMatrix[0][0].real() = -value;
  transMatrix[1][1].real() = -value;
  transMatrix[2][2].real() = +value;
  transMatrix[3][3].real() = +value;

  transMatrix[0][2].real() = +value;
  transMatrix[1][3].real() = +value;
  transMatrix[2][0].real() = +value;
  transMatrix[3][1].real() = +value;

  Complex tmp[4];

  for(int iv = 0 ; iv < this->lV4d ; iv++)
      for(int ic = 0 ; ic < 3 ; ic++){
        memset(tmp,0,4*2*sizeof(double));
        for(int mu = 0 ; mu < 4 ; mu++)
          for(int nu = 0 ; nu < 4 ; nu++)
            tmp[mu] = tmp[mu] + transMatrix[mu][nu] * this->M[iv][nu][ic];
        for(int mu = 0 ; mu < 4 ; mu++)
          this->M[iv][mu][ic] = tmp[mu];
      }
}

void Vector::multiply_by_phase(){
  Complex phase;

    for(int t=0; t < lL[3];t++)
      for(int z=0; z < lL[2];z++)
        for(int y=0; y < lL[1];y++)
          for(int x=0; x < lL[0];x++)
            for(int mu=0; mu<4; mu++)
              for(int c1=0; c1<3; c1++)
                {
		  
		  phase.real() = cos(PI*(t+procPos[3]*lL[3])/((double) L[3]));
		  phase.imag() = sin(PI*(t+procPos[3]*lL[3])/((double) L[3]));
		  M[LEXIC(t,z,y,x,lL)][mu][c1] = phase*M[LEXIC(t,z,y,x,lL)][mu][c1];
                }
	      }
	      
	      
void Vector::sendEvenToOdd(){
  int rv = -1;
  for(int t=0; t < lL[3];t++)
    for(int z=0; z < lL[2];z++)
      for(int y=0; y < lL[1];y++)
	for(int x=0; x < lL[0];x++){
	  int latt_cord = LEXIC(t,z,y,x,lL);
	  int oddBit_x = (x+y+z+t)&1;
	  int oddBit_xp1 = getOddBit(latt_cord+1,lL[0],lL[1],lL[2]);
	  for(int mu=0; mu<4; mu++)
	    for(int c1=0; c1<3; c1++)
	      {
		if(oddBit_x)
		  M[LEXIC(t,z,y,x,lL)][mu][c1] = M[LEXIC(t,z,y,x,lL)+rv][mu][c1];
	      }
	  if(oddBit_x == oddBit_xp1)
	    rv *= -1;
	}

  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu=0; mu<4; mu++)
      for(int c1=0; c1<3; c1++){
	int oddBit = getOddBit(iv,lL[0],lL[1],lL[2]);
	if(!oddBit){
	  M[iv][mu][c1].real() = 0.;
	  M[iv][mu][c1].imag() = 0.;
	}
      }

}

void Vector::normalize(){
  double norm_square = norm2();
  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu=0; mu<4; mu++)
      for(int c1=0; c1<3; c1++){
	M[iv][mu][c1].real() /= sqrt(norm_square);
	M[iv][mu][c1].imag() /= sqrt(norm_square);
      }
}

void Vector::tm(Vector &psi, Gauge &u, LatticeInfo *latInfo){
  this->applyDslash(psi,u);
  this->applyTwistAddDslash(*this,psi,latInfo);
}

void Vector::tm_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo){
  this->applyDslash_chiral(psi,u);
  this->applyTwistAddDslash_chiral(*this,psi,latInfo);
}

void Vector::tmDag(Vector &psi, Gauge &u, LatticeInfo *latInfo){
  this->applyDslashDag(psi,u);
  this->applyTwistDagAddDslashDag(*this,psi,latInfo);
}

void Vector::tmDagtm(Vector &psi,Vector &vec_tmp ,Gauge &u, LatticeInfo *latInfo){

  vec_tmp.applyDslash(psi,u);
  vec_tmp.applyTwistAddDslash(vec_tmp,psi,latInfo);
  this->applyDslashDag(vec_tmp,u);
  this->applyTwistDagAddDslashDag(*this,vec_tmp,latInfo);

}

void Vector::Wilson(Vector &psi, Gauge &u, LatticeInfo *latInfo){ // 1-k*Dslash
  this->applyDslash(psi,u);
  this->applyDeltaAddDslash(*this,psi,latInfo);
}

void Vector::applyCloverTerm(Vector &psi, FieldStrength &fs,LatticeInfo *latInfo){
  Complex I;
  I.real()=0.;
  I.imag()=1.;
  double  coeff=-latInfo->kappa*latInfo->csw;
  // the chiral basis is faster but I do not want to rotate now
  for(int iv = 0 ; iv < lV4d ; iv++){
    Matrix_3x3 F0(fs.M[iv][0]);
    Matrix_3x3 F1(fs.M[iv][1]);
    Matrix_3x3 F2(fs.M[iv][2]);
    Matrix_3x3 F3(fs.M[iv][3]);
    Matrix_3x3 F4(fs.M[iv][4]);
    Matrix_3x3 F5(fs.M[iv][5]);
    
    Vector_3 s0(psi.M[iv][0]);
    Vector_3 s1(psi.M[iv][1]);
    Vector_3 s2(psi.M[iv][2]);
    Vector_3 s3(psi.M[iv][3]);

    Vector_3 r0 = coeff*( (-I)*(F0*s0) + (F1*s1+(-I)*F2*s1) + I*F5*s2 + (I*F3*s3+F4*s3)  );
    Vector_3 r1 = coeff*( ((-1.)*F1*s0+(-I)*F2*s0) + (I)*F0*s1 + ((I)*F3*s2+(-1.)*F4*s2) + (-I)*F5*s3 );
    Vector_3 r2 = coeff*( (I)*F5*s0 + ((I)*F3*s1+F4*s1) + (-I)*F0*s2 + (F1*s3 + (-I)*F2*s3));
    Vector_3 r3 = coeff*( ((I)*F3*s0+(-1.)*F4*s0) + (-I)*F5*s1 + ((-1.)*F1*s2+(-I)*F2*s2) + I*F0*s3);

    for(int c = 0 ; c < 3 ; c++) M[iv][0][c] = r0.M[c];
    for(int c = 0 ; c < 3 ; c++) M[iv][1][c] = r1.M[c];
    for(int c = 0 ; c < 3 ; c++) M[iv][2][c] = r2.M[c];
    for(int c = 0 ; c < 3 ; c++) M[iv][3][c] = r3.M[c];
  }
}

void Vector::Clover(Vector &psi, Gauge &u, FieldStrength &fs, LatticeInfo *latInfo){
  Vector *tmpVec0 = new Vector(latInfo);
  Vector *tmpVec1 = new Vector(latInfo);

  tmpVec0->Wilson(psi,u,latInfo);
  tmpVec1->applyCloverTerm(psi,fs,latInfo);
  
  for(int iv = 0 ; iv < lV4d; iv++)
    for(int s = 0 ; s < 4 ; s++)
      for(int c = 0 ; c < 3 ; c++)
	M[iv][s][c] = tmpVec0->M[iv][s][c] + tmpVec1->M[iv][s][c];

  delete tmpVec0;
  delete tmpVec1;
}

void Vector::tm_xx(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity){
  Vector *tmp = new Vector(latInfo);
  tmp->copy(psi);
  this->zero();

  if( parity == EE){
    tmp->applyDslash_xx(psi,u,OE);
    this->applyTwistInv(*tmp,latInfo);
    tmp->applyDslash_xx(*this,u,EO);
    this->applyTwistAddDslash_PC(*tmp,psi,latInfo,EE);
  }
  else if(parity == OO){
    tmp->applyDslash_xx(psi,u,EO);
    this->applyTwistInv(*tmp,latInfo);
    tmp->applyDslash_xx(*this,u,OE);
    this->applyTwistAddDslash_PC(*tmp,psi,latInfo,OO);
  }
  else
    errorQKXTM("Error wrong parity choise");

  delete tmp;
}

void Vector::tm_xxDagTm_xx(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity){
  Vector *tmp = new Vector(latInfo);
  LatticeInfo latInfo_flip_sign = *latInfo;
  latInfo_flip_sign.twistSign = -latInfo_flip_sign.twistSign;

  tmp->tm_xx(psi,u,latInfo,parity);
  this->applyGamma5(*tmp,latInfo);
  tmp->tm_xx(*this,u,&latInfo_flip_sign,parity);
  this->applyGamma5(*tmp,latInfo);

  delete tmp;
}

void Vector::tm_xx_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity){
  Vector *tmp = new Vector(latInfo);
  tmp->copy(psi);
  this->zero();

  if( parity == EE){
    tmp->applyDslash_xx_chiral(psi,u,OE);
    this->applyTwistInv_chiral(*tmp,latInfo);
    tmp->applyDslash_xx_chiral(*this,u,EO);
    this->applyTwistAddDslash_PC_chiral(*tmp,psi,latInfo,EE);
  }
  else if(parity == OO){
    tmp->applyDslash_xx_chiral(psi,u,EO);
    this->applyTwistInv_chiral(*tmp,latInfo);
    tmp->applyDslash_xx_chiral(*this,u,OE);
    this->applyTwistAddDslash_PC_chiral(*tmp,psi,latInfo,OO);
  }
  else
    errorQKXTM("Error wrong parity choise");

  delete tmp;
}

void Vector::tm_xxDagTm_xx_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity){
  Vector *tmp = new Vector(latInfo);
  LatticeInfo latInfo_flip_sign = *latInfo;
  latInfo_flip_sign.twistSign = -latInfo_flip_sign.twistSign;

  tmp->tm_xx_chiral(psi,u,latInfo,parity);
  this->applyGamma5_chiral(*tmp,latInfo);
  tmp->tm_xx_chiral(*this,u,&latInfo_flip_sign,parity);
  this->applyGamma5_chiral(*tmp,latInfo);

  delete tmp;
}

void Propagator::absorbVector(Vector &phi ,int nu , int c2 ){

  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int c1 = 0 ; c1 < NCOLORS ; c1++){
	M[iv][mu][nu][c1][c2] = phi.M[iv][mu][c1];
      }

}


void Vector::seqSourceFixSinkPart1(Propagator &prop1,Propagator &prop2,int timeslice,int nu, int c2 , int flagProj){

  int my_timeslice = timeslice - procPos[3]*lL[3];

  if(my_timeslice >= 0 && my_timeslice < lL[3]){
  if(flagProj != 1 && flagProj != 2)errorQKXTM("Error flag must be 1 or 2");

  zero();

  Complex projector[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      projector[mu][nu].real() = 0 ; projector[mu][nu].imag() = 0;
    }

  if(flagProj == 1){
    projector[0][0].real() = 0.5;
    projector[1][1].real() = 0.5;
  }
  else{
    projector[0][0].real() = 0.5;
    projector[0][1].real() = 0.5;  projector[0][1].imag() = -0.5;
    projector[1][0].real() = 0.5;  projector[1][0].imag() = 0.5;
    projector[1][1].real() = -0.5;
  }

  Complex Cg5[4][4];
  Complex Cg5_bar[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      Cg5[mu][nu].real() = 0.; Cg5[mu][nu].imag()=0.;
      Cg5_bar[mu][nu].real() = 0.; Cg5_bar[mu][nu].imag()=0.;
    }

  Cg5[0][1].real() = 1.;
  Cg5[1][0].real() = -1.;
  Cg5[2][3].real() = 1.;
  Cg5[3][2].real() = -1.;
                         
  Cg5_bar[0][1].real() = -1.;
  Cg5_bar[1][0].real() = 1.;
  Cg5_bar[2][3].real() = -1.;
  Cg5_bar[3][2].real() = 1.;


  Complex C_temp;
  Complex Cg5Cg5bar_val[16*16];
  unsigned short int Cg5Cg5bar_ind[16*16][4];
  int counter = 0;
  // adopt Thomas convention
  for(unsigned short int ksi = 0 ; ksi < 4 ; ksi++)
    for(unsigned short int pi = 0 ; pi < 4 ; pi++)
      for(unsigned short int theta = 0 ; theta < 4 ; theta++)
	for(unsigned short int rho = 0 ; rho < 4 ; rho++){
	  C_temp = Cg5[ksi][pi] * Cg5_bar[theta][rho];
	  if( norm(C_temp) > 1e-3 ){
	    Cg5Cg5bar_val[counter] = C_temp;
	    Cg5Cg5bar_ind[counter][0] = ksi;
	    Cg5Cg5bar_ind[counter][1] = pi;
	    Cg5Cg5bar_ind[counter][2] = theta;
	    Cg5Cg5bar_ind[counter][3] = rho;
	    counter++;
	  }

	}
  Complex factor;

  for(int cc1 = 0 ; cc1 < 6 ; cc1++){
    int a = eps[cc1][0];
    int b = eps[cc1][1];
    int c = eps[cc1][2];
    for(int cc2 = 0 ; cc2 < 6 ; cc2++){
      int ap = eps[cc2][0];
      int bp = eps[cc2][1];
      int cp = eps[cc2][2];
      if(ap == c2){
	for(int idx = 0 ; idx < counter ; idx++){
	  int ksi = Cg5Cg5bar_ind[idx][0];
	  int pi = Cg5Cg5bar_ind[idx][1];
	  int theta = Cg5Cg5bar_ind[idx][2];
	  int rho = Cg5Cg5bar_ind[idx][3];

	  for(int delta = 0 ; delta < 4 ; delta++)
	    for(int epsilon = 0 ; epsilon < 4 ; epsilon++){
	      if( norm(projector[delta][epsilon]) > 1e-3){
		
		for(int mu = 0 ; mu < 4 ; mu++){
		  factor = sgn_eps[cc1]*sgn_eps[cc2]*Cg5Cg5bar_val[idx] * projector[delta][epsilon];
		  
		  for(int iz = 0 ; iz < lL[2] ; iz++)
		    for(int iy = 0 ; iy < lL[1] ; iy++)
		      for(int ix = 0 ; ix < lL[0] ; ix++){
			int iv = LEXIC(timeslice,iz,iy,ix,lL);
			
			if(delta == nu && epsilon == mu) M[iv][mu][a] = M[iv][mu][a] + factor * prop2.M[iv][pi][theta][b][bp] * prop1.M[iv][ksi][rho][c][cp];
			if(rho == nu && epsilon == mu) M[iv][mu][a] = M[iv][mu][a] + factor * prop2.M[iv][pi][theta][b][bp] * prop1.M[iv][ksi][delta][c][cp];
			if(delta == nu && ksi == mu) M[iv][mu][a] = M[iv][mu][a] + factor * prop2.M[iv][pi][theta][b][bp] * prop1.M[iv][epsilon][rho][c][cp];
			if(rho == nu && ksi == mu) M[iv][mu][a] = M[iv][mu][a] + factor * prop2.M[iv][pi][theta][b][bp] * prop1.M[iv][epsilon][delta][c][cp];
			
		      } // volume
		    
		}// mu statement

	      } // if statement

	    } // epsilon loop

	} // spin
      } // if statement
    } // color
  } // color

  }
}


void Vector::seqSourceFixSinkPart2(Propagator &prop1,int timeslice,int nu, int c2 , int flagProj){

  if(flagProj != 1 && flagProj != 2)errorQKXTM("Error flag must be 1 or 2");

  zero();

  Complex projector[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      projector[mu][nu].real() = 0 ; projector[mu][nu].imag() = 0;
    }

  if(flagProj == 1){
    projector[0][0].real() = 0.5;
    projector[1][1].real() = 0.5;
  }
  else{
    projector[0][0].real() = 0.5;
    projector[0][1].real() = 0.5;  projector[0][1].imag() = -0.5;
    projector[1][0].real() = 0.5;  projector[1][0].imag() = 0.5;
    projector[1][1].real() = -0.5;
  }

  Complex Cg5[4][4];
  Complex Cg5_bar[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      Cg5[mu][nu].real() = 0.; Cg5[mu][nu].imag()=0.;
      Cg5_bar[mu][nu].real() = 0.; Cg5_bar[mu][nu].imag()=0.;
    }

  Cg5[0][1].real() = 1.;
  Cg5[1][0].real() = -1.;
  Cg5[2][3].real() = 1.;
  Cg5[3][2].real() = -1.;
                         
  Cg5_bar[0][1].real() = -1.;
  Cg5_bar[1][0].real() = 1.;
  Cg5_bar[2][3].real() = -1.;
  Cg5_bar[3][2].real() = 1.;


  Complex C_temp;
  Complex Cg5Cg5bar_val[16*16];
  unsigned short int Cg5Cg5bar_ind[16*16][4];
  int counter = 0;

  for(unsigned short int ksi = 0 ; ksi < 4 ; ksi++)
    for(unsigned short int pi = 0 ; pi < 4 ; pi++)
      for(unsigned short int theta = 0 ; theta < 4 ; theta++)
	for(unsigned short int rho = 0 ; rho < 4 ; rho++){
	  C_temp = Cg5[ksi][pi] * Cg5_bar[theta][rho];
	  if( norm(C_temp) > 1e-3 ){
	    Cg5Cg5bar_val[counter] = C_temp;
	    Cg5Cg5bar_ind[counter][0] = ksi;
	    Cg5Cg5bar_ind[counter][1] = pi;
	    Cg5Cg5bar_ind[counter][2] = theta;
	    Cg5Cg5bar_ind[counter][3] = rho;
	    counter++;
	  }

	}
  Complex factor;

  for(int cc1 = 0 ; cc1 < 6 ; cc1++){
    int a = eps[cc1][0];
    int b = eps[cc1][1];
    int c = eps[cc1][2];
    for(int cc2 = 0 ; cc2 < 6 ; cc2++){
      int ap = eps[cc2][0];
      int bp = eps[cc2][1];
      int cp = eps[cc2][2];
      if(ap == c2){
	for(int idx = 0 ; idx < counter ; idx++){
	  int alpha = Cg5Cg5bar_ind[idx][0];
	  int mu = Cg5Cg5bar_ind[idx][1];
	  int mup = Cg5Cg5bar_ind[idx][2];
	  int delta = Cg5Cg5bar_ind[idx][3];

	  if(mup == nu){
	    for(int beta = 0 ; beta < 4 ; beta++)
	      for(int gamma = 0 ; gamma < 4 ; gamma++){
		if( norm(projector[beta][gamma]) > 1e-3){
		 
		  factor = sgn_eps[cc1]*sgn_eps[cc2]*Cg5Cg5bar_val[idx] * projector[beta][gamma];
		  
		  for(int iz = 0 ; iz < lL[2] ; iz++)
		    for(int iy = 0 ; iy < lL[1] ; iy++)
		      for(int ix = 0 ; ix < lL[0] ; ix++){
			int iv = LEXIC(timeslice,iz,iy,ix,lL);
			
		        M[iv][mu][a] = M[iv][mu][a] + factor * prop1.M[iv][alpha][beta][b][bp] * prop1.M[iv][gamma][delta][c][cp];
			M[iv][mu][a] = M[iv][mu][a] + factor * prop1.M[iv][gamma][beta][b][bp] * prop1.M[iv][alpha][delta][c][cp];
			
		      } // volume
		    
		}// if statement
		
	      } // gamma loop

	  } // if statement
	  
	} // spin
      } // if statement
    } // color
  } // color
  
}


void Vector::rescale(double a){
  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
	M[iv][mu][ic].real() = M[iv][mu][ic].real() * a;
	M[iv][mu][ic].imag() = M[iv][mu][ic].imag() * a;
      }
}


void Propagator::rotateToPhysicalBasePlus(){

  Complex temp[4][4][3][3];
  Complex temp2[4][4][3][3];

  Complex Rotate[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      Rotate[mu][nu].real() =0.;
      Rotate[mu][nu].imag() =0.;
    }

  Rotate[0][0].real()= 1.;
  Rotate[1][1].real()= 1.;
  Rotate[2][2].real()= 1.;
  Rotate[3][3].real()= 1.;

  Rotate[0][2].imag()= 1.;
  Rotate[1][3].imag()= 1.;
  Rotate[2][0].imag()= 1.;
  Rotate[3][1].imag()= 1.;

  for(int iv = 0 ; iv < lV4d ; iv++){

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3; c2++){
	    temp[mu][nu][c1][c2] = M[iv][mu][nu][c1][c2];
	    temp2[mu][nu][c1][c2].real()=0.;	    temp2[mu][nu][c1][c2].imag()=0.;
	  }
    
    for(int c1 = 0 ; c1 < 3 ; c1++)
      for(int c2 = 0 ; c2 < 3 ; c2++){

	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int mup = 0 ; mup < 4 ; mup++)
	    for(int nup = 0 ; nup < 4 ; nup++)
	      for(int nu = 0 ; nu < 4 ; nu++){
		temp2[mu][nu][c1][c2] = temp2[mu][nu][c1][c2] + Rotate[mu][mup] * temp[mup][nup][c1][c2] * Rotate[nup][nu];
	      }
	
      }

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3; c2++){
	    M[iv][mu][nu][c1][c2].real() = 0.5 * temp2[mu][nu][c1][c2].real(); 
	    M[iv][mu][nu][c1][c2].imag() = 0.5 * temp2[mu][nu][c1][c2].imag(); 
	  }
    

  }



}


void Propagator::rotateToPhysicalBaseMinus(){


  Complex temp[4][4][3][3];
  Complex temp2[4][4][3][3];

  Complex Rotate[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      Rotate[mu][nu].real() =0.;
      Rotate[mu][nu].imag() =0.;
    }

  Rotate[0][0].real()= 1.;
  Rotate[1][1].real()= 1.;
  Rotate[2][2].real()= 1.;
  Rotate[3][3].real()= 1.;

  Rotate[0][2].imag()= -1.;
  Rotate[1][3].imag()= -1.;
  Rotate[2][0].imag()= -1.;
  Rotate[3][1].imag()= -1.;

  for(int iv = 0 ; iv < lV4d ; iv++){

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3; c2++){
	    temp[mu][nu][c1][c2] = M[iv][mu][nu][c1][c2];
	    temp2[mu][nu][c1][c2].real()=0.;	    temp2[mu][nu][c1][c2].imag()=0.;
	  }
    
    for(int c1 = 0 ; c1 < 3 ; c1++)
      for(int c2 = 0 ; c2 < 3 ; c2++){

	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int mup = 0 ; mup < 4 ; mup++)
	    for(int nup = 0 ; nup < 4 ; nup++)
	      for(int nu = 0 ; nu < 4 ; nu++){
		temp2[mu][nu][c1][c2] = temp2[mu][nu][c1][c2] + Rotate[mu][mup] * temp[mup][nup][c1][c2] * Rotate[nup][nu];
	      }
	
      }

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3; c2++){
	    M[iv][mu][nu][c1][c2].real() = 0.5 * temp2[mu][nu][c1][c2].real(); 
	    M[iv][mu][nu][c1][c2].imag() = 0.5 * temp2[mu][nu][c1][c2].imag(); 
	  }
    

  }




}


////////////////////////////////////////////////////////////////

// the following are tests for disconnected loops

// for mpi use global static variables                                                                                                                                         


void Vector::standardEndTrick_sigmaTerm(Complex *loops, LatticeInfo *latInfo ){
  
  Complex *local_loops = (Complex*) calloc(lL[3]*7,sizeof(Complex));
  Complex *total_loops = (Complex*) calloc(lL[3]*7,sizeof(Complex));

  Complex phase;
  double tmp;

  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;

  for(int lt = 0 ; lt < lL[3] ; lt++){
    for(int imom = 0 ; imom < 7 ; imom++){

      for(int lx = 0 ; lx < lL[0] ; lx++)
	for(int ly = 0 ; ly < lL[1] ; ly++)
	  for(int lz = 0 ; lz < lL[2] ; lz++){
	    int iv = LEXIC(lt,lz,ly,lx,lL); 

	    x = lx + procPos[0]*lL[0];
	    y = ly + procPos[1]*lL[1];
	    z = lz + procPos[2]*lL[2];
	    tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	    phase = (Complex) {cos(tmp) , -sin(tmp)};


	    for(int alpha = 0 ; alpha < 4 ; alpha++)	    
	      for(int a = 0 ; a < 3 ; a++)
		local_loops[lt*7+imom] = local_loops[lt*7+imom] + conj(this->M[iv][alpha][a]) * this->M[iv][alpha][a] * phase;
		
	      
	  } // spatial
    } // momentum

  } // time

  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                                                  
    MPI_Gather(total_loops,lL[3]*7*2,MPI_DOUBLE,loops,lL[3]*7*2,MPI_DOUBLE,0,timeComm);
  }

  if(rank == 0){

    for(int it = 0 ; it < L[3] ; it++)
      for(int imom = 0 ; imom < 7 ; imom++){
	loops[it*7+imom].real() = -8*(latInfo->mu)*(latInfo->kappa)*(latInfo->kappa)*loops[it*7+imom].real();
	loops[it*7+imom].imag() = -8*(latInfo->mu)*(latInfo->kappa)*(latInfo->kappa)*loops[it*7+imom].imag();
      }

  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);

}

void Vector::applyConjugation(){
  for(int iv = 0 ; iv < lV4d ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int c1 = 0 ; c1 < NCOLORS ; c1++)
        M[iv][mu][c1] = conj(M[iv][mu][c1]);
}

void Vector::applyg5gm(int mu){
  Complex cmplx_unit;
  Complex temp[4];

  cmplx_unit = (Complex) {0.,1};
  switch(mu)
    {
    case 0:  // g5g1
      for(int iv = 0 ; iv < lV4d ; iv++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++){

	  for(int mu = 0 ; mu < NSPINS ; mu++)
	    temp[mu] = M[iv][mu][c1];
	  
	  M[iv][0][c1] = -cmplx_unit*temp[1];
	  M[iv][1][c1] = -cmplx_unit*temp[0];
	  M[iv][2][c1] = cmplx_unit*temp[3];
	  M[iv][3][c1] = cmplx_unit*temp[2];
	}
      break;
    case 1: // g5g2
      for(int iv = 0 ; iv < lV4d ; iv++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++){

	  for(int mu = 0 ; mu < NSPINS ; mu++)
	    temp[mu] = M[iv][mu][c1];
	  
	  M[iv][0][c1] = -temp[1];
	  M[iv][1][c1] = temp[0];
	  M[iv][2][c1] = temp[3];
	  M[iv][3][c1] = -temp[2];
	}

      break;
    case 2: // g5g3
      for(int iv = 0 ; iv < lV4d ; iv++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++){

	  for(int mu = 0 ; mu < NSPINS ; mu++)
	    temp[mu] = M[iv][mu][c1];
	  
	  M[iv][0][c1] = -cmplx_unit*temp[0];
	  M[iv][1][c1] = cmplx_unit*temp[1];
	  M[iv][2][c1] = cmplx_unit*temp[2];
	  M[iv][3][c1] = -cmplx_unit*temp[3];
	}

      break;
    case 3: // g5g4
      for(int iv = 0 ; iv < lV4d ; iv++)
	for(int c1 = 0 ; c1 < NCOLORS ; c1++){

	  for(int mu = 0 ; mu < NSPINS ; mu++)
	    temp[mu] = M[iv][mu][c1];
	  
	  M[iv][0][c1] = -temp[2];
	  M[iv][1][c1] = -temp[3];
	  M[iv][2][c1] = temp[0];
	  M[iv][3][c1] = temp[1];
	}

      break;
    }

}



void Vector::generalizedEndTrick_axial(Gauge &u,Complex *loops, LatticeInfo *latInfo ){
  


  // prepare vector
  Vector *tmpVec0 = new Vector(latInfo);
  Vector *tmpVec1 = new Vector(latInfo);
  Vector *tmpVec2 = new Vector(latInfo);
  Vector *tmpVec3 = new Vector(latInfo);
  Vector *tmpVec = new Vector(latInfo);

  tmpVec->copy(*this);
  tmpVec0->Wilson(*tmpVec,u,latInfo);
  tmpVec1->Wilson(*tmpVec,u,latInfo);
  tmpVec2->Wilson(*tmpVec,u,latInfo);
  tmpVec3->Wilson(*tmpVec,u,latInfo);

  tmpVec0->applyg5gm(0);
  tmpVec1->applyg5gm(1);
  tmpVec2->applyg5gm(2);
  tmpVec3->applyg5gm(3);


  Complex *local_loops = (Complex*) calloc(lL[3]*4*7,sizeof(Complex));
  Complex *total_loops = (Complex*) calloc(lL[3]*4*7,sizeof(Complex));
  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
    for(int mu = 0 ; mu < 4 ; mu++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      

	      for(int alpha = 0 ; alpha < 4 ; alpha++)	    
		for(int a = 0 ; a < 3 ; a++){
		  if(mu == 0) local_loops[lt*4*7+mu*7+imom] = local_loops[lt*4*7+mu*7+imom] + conj(M[iv][alpha][a]) * tmpVec0->M[iv][alpha][a] * phase;
		  if(mu == 1) local_loops[lt*4*7+mu*7+imom] = local_loops[lt*4*7+mu*7+imom] + conj(M[iv][alpha][a]) * tmpVec1->M[iv][alpha][a] * phase;
		  if(mu == 2) local_loops[lt*4*7+mu*7+imom] = local_loops[lt*4*7+mu*7+imom] + conj(M[iv][alpha][a]) * tmpVec2->M[iv][alpha][a] * phase;
		  if(mu == 3) local_loops[lt*4*7+mu*7+imom] = local_loops[lt*4*7+mu*7+imom] + conj(M[iv][alpha][a]) * tmpVec3->M[iv][alpha][a] * phase;
		}

	    } // spatial
      } // momentum
    } // mu
  } // time

  
  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*7*2,MPI_DOUBLE,loops,lL[3]*4*7*2,MPI_DOUBLE,0,timeComm);
  }

  if(rank == 0){


    for(int it = 0 ; it < L[3] ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < 7 ; imom++){
	  loops[it*4*7+mu*7+imom].real() = 4*(latInfo->kappa)*loops[it*4*7+mu*7+imom].real();
	  loops[it*4*7+mu*7+imom].imag() = 4*(latInfo->kappa)*loops[it*4*7+mu*7+imom].imag();
	}

  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);
  

  delete tmpVec0;
  delete tmpVec1;
  delete tmpVec2;
  delete tmpVec3;
  delete tmpVec;

}


void Vector::generalizedEndTrick_axial_openIndices(Gauge &u,Complex *loops, LatticeInfo *latInfo ){
  


  // prepare vector
  Vector *tmpVecRight = new Vector(latInfo);
  Vector *tmpVecLeft = new Vector(latInfo);
  
  tmpVecRight->Wilson(*this,u,latInfo);
  tmpVecRight->applyGamma5(*tmpVecRight,latInfo);
  tmpVecLeft->daggerVectorGamma5(*this,latInfo);

  Complex *local_loops = (Complex*) calloc(lL[3]*4*4*7,sizeof(Complex));
  Complex *total_loops = (Complex*) calloc(lL[3]*4*4*7,sizeof(Complex));
  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      

	      for(int alpha = 0 ; alpha < 4 ; alpha++)
		for(int beta = 0 ; beta < 4 ; beta++)
		  for(int a = 0 ; a < 3 ; a++){
		    local_loops[lt*4*4*7+alpha*4*7+beta*7+imom] = local_loops[lt*4*4*7+alpha*4*7+beta*7+imom] + tmpVecLeft->M[iv][alpha][a] * tmpVecRight->M[iv][beta][a] * phase;
		  }

	    } // spatial
      } // momentum
  } // time

  
  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*4*7*2,MPI_DOUBLE,loops,lL[3]*4*4*7*2,MPI_DOUBLE,0,timeComm);
  }

  /*
  if(rank == 0){


    for(int it = 0 ; it < L[3] ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < 7 ; imom++){
	  loops[it*4*7+mu*7+imom].real() = 4*(latInfo->kappa)*loops[it*4*7+mu*7+imom].real();
	  loops[it*4*7+mu*7+imom].imag() = 4*(latInfo->kappa)*loops[it*4*7+mu*7+imom].imag();
	}

  }
  */

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);
  

  delete tmpVecLeft;
  delete tmpVecRight;

}

void Vector::generalizedEndTrick_covDer(Gauge &u,Complex *loops, Complex *loops_cv , LatticeInfo *latInfo ){

  //  this->communicateToPlus();
  //  this->communicateToMinus();
  u.communicateToPlus();
  
  Vector *tmpVecRight = new Vector(latInfo);
  Vector *tmpVecLeft = new Vector(latInfo);
  FieldStrength *fs1 = new FieldStrength(latInfo,u);
  tmpVecRight->Clover(*this,u,*fs1,latInfo);
  //  tmpVecRight->Wilson(*this,u,latInfo);
  tmpVecRight->applyGamma5(*tmpVecRight,latInfo);
  tmpVecLeft->daggerVectorGamma5(*this,latInfo);

  tmpVecRight->communicateToPlus();
  tmpVecRight->communicateToMinus();
  tmpVecLeft->communicateToPlus();
  tmpVecLeft->communicateToMinus();
  
  Complex *local_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex *local_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      
	      for(int nu = 0 ; nu < 4 ; nu++)
		for(int alpha = 0 ; alpha < 4 ; alpha++)
		  for(int beta = 0 ; beta < 4 ; beta++)
		    for(int a = 0 ; a < 3 ; a++)
		      for(int b = 0 ; b < 3 ; b++){
			local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( tmpVecLeft->M[iv][alpha][a] * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   - tmpVecLeft->M[iv][alpha][a] * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   - tmpVecLeft->M[pointPlus[nu][iv]][alpha][a] * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + tmpVecLeft->M[pointMinus[nu][iv]][alpha][a] * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

			local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( tmpVecLeft->M[iv][alpha][a] * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   + tmpVecLeft->M[iv][alpha][a] * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   + tmpVecLeft->M[pointPlus[nu][iv]][alpha][a] * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + tmpVecLeft->M[pointMinus[nu][iv]][alpha][a] * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

		  }


	    } // spatial
      } // momentum
  } // time

  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  MPI_Reduce(&(local_loops_cv[0]) , &(total_loops_cv[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
    MPI_Gather(total_loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);

  free(local_loops_cv);
  free(total_loops_cv);
  

  delete tmpVecLeft;
  delete tmpVecRight;
  delete fs1;

}

void Vector::volumeSource_ultralocal(Vector &xi, Gauge &u,Complex *loops, LatticeInfo *latInfo ){
  


  // prepare vector
  Vector *tmpVecRight = new Vector(latInfo);
  Vector *tmpVecLeft = new Vector(latInfo);
  
  tmpVecRight->copy(*this);
  tmpVecLeft->copy(xi);

  Complex *local_loops = (Complex*) calloc(lL[3]*4*4*7,sizeof(Complex));
  Complex *total_loops = (Complex*) calloc(lL[3]*4*4*7,sizeof(Complex));
  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      

	      for(int alpha = 0 ; alpha < 4 ; alpha++)
		for(int beta = 0 ; beta < 4 ; beta++)
		  for(int a = 0 ; a < 3 ; a++){
		    local_loops[lt*4*4*7+alpha*4*7+beta*7+imom] = local_loops[lt*4*4*7+alpha*4*7+beta*7+imom] + conj(tmpVecLeft->M[iv][alpha][a]) * tmpVecRight->M[iv][beta][a] * phase;
		  }

	    } // spatial
      } // momentum
  } // time

  
  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*4*7*2,MPI_DOUBLE,loops,lL[3]*4*4*7*2,MPI_DOUBLE,0,timeComm);
  }


  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);
  

  delete tmpVecLeft;
  delete tmpVecRight;

}


void Vector::volumeSource_covDer(Vector &xi, Gauge &u, Complex *loops, Complex *loops_cv, LatticeInfo *latInfo ){

  u.communicateToPlus();
  
  Vector *tmpVecRight = new Vector(latInfo);
  Vector *tmpVecLeft = new Vector(latInfo);
  
  tmpVecRight->copy(*this);
  tmpVecLeft->copy(xi);

  tmpVecRight->communicateToPlus();
  tmpVecRight->communicateToMinus();
  tmpVecLeft->communicateToPlus();
  tmpVecLeft->communicateToMinus();
  
  Complex *local_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex *local_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      
	      for(int nu = 0 ; nu < 4 ; nu++)
		for(int alpha = 0 ; alpha < 4 ; alpha++)
		  for(int beta = 0 ; beta < 4 ; beta++)
		    for(int a = 0 ; a < 3 ; a++)
		      for(int b = 0 ; b < 3 ; b++){
			local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( conj(tmpVecLeft->M[iv][alpha][a]) * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   - conj(tmpVecLeft->M[iv][alpha][a]) * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   - conj(tmpVecLeft->M[pointPlus[nu][iv]][alpha][a]) * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + conj(tmpVecLeft->M[pointMinus[nu][iv]][alpha][a]) * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

			local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( conj(tmpVecLeft->M[iv][alpha][a]) * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   + conj(tmpVecLeft->M[iv][alpha][a]) * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   + conj(tmpVecLeft->M[pointPlus[nu][iv]][alpha][a]) * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + conj(tmpVecLeft->M[pointMinus[nu][iv]][alpha][a]) * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

		  }


	    } // spatial
      } // momentum
  } // time

  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  MPI_Reduce(&(local_loops_cv[0]) , &(total_loops_cv[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
    MPI_Gather(total_loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);

  free(local_loops_cv);
  free(total_loops_cv);
  

  delete tmpVecLeft;
  delete tmpVecRight;


}

void Vector::standardEndTrick_covDer(Gauge &u,Complex *loops, Complex *loops_cv , LatticeInfo *latInfo ){

  //  this->communicateToPlus();
  //  this->communicateToMinus();
  u.communicateToPlus();
  
  Vector *tmpVecRight = new Vector(latInfo);
  Vector *tmpVecLeft = new Vector(latInfo);
  
  tmpVecRight->copy(*this);
  tmpVecLeft->daggerVectorGamma5(*this,latInfo);

  tmpVecRight->communicateToPlus();
  tmpVecRight->communicateToMinus();
  tmpVecLeft->communicateToPlus();
  tmpVecLeft->communicateToMinus();
  
  Complex *local_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex *local_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex)); // time,covDer index,dirac,dirac,momentum
  Complex *total_loops_cv = (Complex*) calloc(lL[3]*4*4*4*7,sizeof(Complex));

  Complex phase;
  double tmp;
  int x,y,z;

  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;
  

  for(int lt = 0 ; lt < lL[3] ; lt++){
      for(int imom = 0 ; imom < 7 ; imom++){

	for(int lx = 0 ; lx < lL[0] ; lx++)
	  for(int ly = 0 ; ly < lL[1] ; ly++)
	    for(int lz = 0 ; lz < lL[2] ; lz++){
	      int iv = LEXIC(lt,lz,ly,lx,lL); 
	      
	      x = lx + procPos[0]*lL[0];
	      y = ly + procPos[1]*lL[1];
	      z = lz + procPos[2]*lL[2];
	      tmp = ( ((double)momList[0][imom]*x)/L[0] + ((double)momList[1][imom]*y)/L[1] + ((double)momList[2][imom]*z)/L[2]  )*2*PI;
	      phase = (Complex) {cos(tmp) , -sin(tmp)};
	      
	      for(int nu = 0 ; nu < 4 ; nu++)
		for(int alpha = 0 ; alpha < 4 ; alpha++)
		  for(int beta = 0 ; beta < 4 ; beta++)
		    for(int a = 0 ; a < 3 ; a++)
		      for(int b = 0 ; b < 3 ; b++){
			local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( tmpVecLeft->M[iv][alpha][a] * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   - tmpVecLeft->M[iv][alpha][a] * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   - tmpVecLeft->M[pointPlus[nu][iv]][alpha][a] * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + tmpVecLeft->M[pointMinus[nu][iv]][alpha][a] * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

			local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom] = local_loops_cv[lt*4*4*4*7+nu*4*4*7+alpha*4*7+beta*7+imom]
			  + 0.25*( tmpVecLeft->M[iv][alpha][a] * u.M[iv][nu][a][b] * tmpVecRight->M[pointPlus[nu][iv]][beta][b] 
				   + tmpVecLeft->M[iv][alpha][a] * conj(u.M[pointMinus[nu][iv]][nu][b][a]) * tmpVecRight->M[pointMinus[nu][iv]][beta][b]
				   + tmpVecLeft->M[pointPlus[nu][iv]][alpha][a] * conj(u.M[iv][nu][b][a]) * tmpVecRight->M[iv][beta][b]
				   + tmpVecLeft->M[pointMinus[nu][iv]][alpha][a] * u.M[pointMinus[nu][iv]][nu][a][b] * tmpVecRight->M[iv][beta][b]
				   )* phase;

		  }


	    } // spatial
      } // momentum
  } // time

  MPI_Reduce(&(local_loops[0]) , &(total_loops[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  MPI_Reduce(&(local_loops_cv[0]) , &(total_loops_cv[0]) , lL[3]*4*4*4*7*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);
  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(total_loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
    MPI_Gather(total_loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,loops_cv,lL[3]*4*4*4*7*2,MPI_DOUBLE,0,timeComm);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(local_loops);
  free(total_loops);

  free(local_loops_cv);
  free(total_loops_cv);
  

  delete tmpVecLeft;
  delete tmpVecRight;


}



void Gauge::calculateGluonLoop(Complex *gluonLoop){ 

  communicateToMinus();
  Complex *localGluon = (Complex*)calloc(lL[3],sizeof(Complex));
  Complex *totalGluon = (Complex*)calloc(lL[3],sizeof(Complex));

  Complex plaq[3][3],tmp[3][3];
  Complex tracePlaq;
  tracePlaq.real() = 0.;
  tracePlaq.imag() = 0.;

  //+++++++++++++++++++++++++++++++++++++   XY  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][0],M[pointPlus[0][iv]][1]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[1][iv]][0]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][1]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+0*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++   XZ  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][0],M[pointPlus[0][iv]][2]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[2][iv]][0]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][2]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+1*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++   XT  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][0],M[pointPlus[0][iv]][3]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[3][iv]][0]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][3]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+2*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


  //+++++++++++++++++++++++++++++++++++++   YZ  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][1],M[pointPlus[1][iv]][2]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[2][iv]][1]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][2]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+3*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++   YT  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][1],M[pointPlus[1][iv]][3]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[3][iv]][1]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][3]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+4*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++   ZT  ++++++++++++++++++++++++++++++++++++++++++++++//
  for(int lt = 0 ; lt < lL[3] ; lt++)  
    for(int lx = 0 ; lx < lL[0] ; lx++)
      for(int ly = 0 ; ly < lL[1] ; ly++)
	for(int lz = 0 ; lz < lL[2] ; lz++){

	  int iv = LEXIC(lt,lz,ly,lx,lL); 
	  qcd_MUL3x3(plaq,M[iv][2],M[pointPlus[2][iv]][3]);
	  qcd_MULADJOINT3x3(tmp,plaq,M[pointPlus[3][iv]][2]);
	  qcd_MULADJOINT3x3(plaq,tmp,M[iv][3]);

	  localGluon[lt] = localGluon[lt] + plaq[0][0] + plaq[1][1] + plaq[2][2];

	}

  MPI_Reduce(&(localGluon[0]) , &(totalGluon[0]) , lL[3]*2 , MPI_DOUBLE , MPI_SUM , 0 , spaceComm);  
  if(timeRank >= 0 && timeRank < P[3] ){                                                                                                
    MPI_Gather(totalGluon,lL[3]*2,MPI_DOUBLE,gluonLoop+5*L[3],lL[3]*2,MPI_DOUBLE,0,timeComm);
  }

  for(int lt = 0 ; lt < lL[3] ; lt++){
    localGluon[lt].real() = 0.;
    localGluon[lt].imag() = 0.;
    totalGluon[lt].real() = 0.;
    totalGluon[lt].imag() = 0.;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  free(localGluon);
  free(totalGluon);

  MPI_Barrier(MPI_COMM_WORLD);


}


static void createMom(int *Nmom, int momElem[][3],int Q_sq){
  int counter = 0;

  for(int iQ = 0 ; iQ <= Q_sq ; iQ++){
    for(int nx = iQ ; nx >= -iQ ; nx--)
      for(int ny = iQ ; ny >= -iQ ; ny--)
        for(int nz = iQ ; nz >= -iQ ; nz--){
          if( nx*nx + ny*ny + nz*nz == iQ ){
            momElem[counter][0] = nx;
            momElem[counter][1] = ny;
            momElem[counter][2] = nz;
            counter++;
          }
        }
  }


  *Nmom = counter;

  if(*Nmom > MOM_MAX){
    errorQKXTM("Error Nmom can't be greater than the maximum number of momenta\n");
  }

}



void contractNucleon(Propagator &prop1, Propagator &prop2, char *output, int Q_sq, int x_src[4]){

  int Nmom;
  int momElem[MOM_MAX][3];

  double t1;
  double t2;
  
  /*
  Complex *block[4][4];
  for(int gamma = 0 ; gamma < 4 ; gamma++)
    for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++){
      block[gamma][gamma1] = (Complex*)malloc(prop1.lV3d*sizeof(Complex));
      if(block[gamma][gamma1] == NULL) errorQKXTM("Error out of memory\n");
    }
  */
  Complex *block = (Complex*)malloc(prop1.lV3d*4*4*sizeof(Complex));
  if(block == NULL)errorQKXTM("Error out of memory\n");

  createMom(&Nmom,momElem,Q_sq);


  if(prop1.P[3] != 1){
    errorQKXTM("Error number of process in time direction must be 1\n");
  }
  FILE *ptr_out = NULL;
  if(prop1.rank == 0){
    ptr_out = fopen(output,"w");
    if(ptr_out == NULL){
      errorQKXTM("Error open file for writting \n");
    }
  }

  Complex Cg5[4][4];
  Complex Cg5_bar[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      Cg5[mu][nu].real() = 0.; Cg5[mu][nu].imag()=0.;
      Cg5_bar[mu][nu].real() = 0.; Cg5_bar[mu][nu].imag()=0.;
    }

  Cg5[0][1].real() = 1.;
  Cg5[1][0].real() = -1.;
  Cg5[2][3].real() = 1.;
  Cg5[3][2].real() = -1.;
                         
  Cg5_bar[0][1].real() = -1.;
  Cg5_bar[1][0].real() = 1.;
  Cg5_bar[2][3].real() = -1.;
  Cg5_bar[3][2].real() = 1.;


  Complex C_temp;
  Complex Cg5Cg5bar_val[16*16];
  unsigned short int Cg5Cg5bar_ind[16*16][4];
  int counter = 0;

  for(unsigned short int alpha = 0 ; alpha < 4 ; alpha++)
    for(unsigned short int beta = 0 ; beta < 4 ; beta++)
      for(unsigned short int beta1 = 0 ; beta1 < 4 ; beta1++)
	for(unsigned short int alpha1 = 0 ; alpha1 < 4 ; alpha1++){
	  C_temp = Cg5[alpha][beta] * Cg5_bar[beta1][alpha1];
	  if( norm(C_temp) > 1e-3 ){
	    Cg5Cg5bar_val[counter] = C_temp;
	    Cg5Cg5bar_ind[counter][0] = alpha;
	    Cg5Cg5bar_ind[counter][1] = beta;
	    Cg5Cg5bar_ind[counter][2] = beta1;
	    Cg5Cg5bar_ind[counter][3] = alpha1;
	    counter++;
	  }
	}

  unsigned short int a,b,c,a1,b1,c1;
  Complex factor;
  Complex test;
  Complex test2;
  Complex *corr = (Complex*) calloc(Nmom*4*4,sizeof(Complex));;
  Complex *corr2 = (Complex*) calloc(Nmom*4*4,sizeof(Complex));

  for(int it = 0 ; it < prop1.L[3] ; it++){
    if(prop1.rank == 0 ) printf("it = %d ",it);
    t1 = MPI_Wtime();
    memset(&(block[0].real()),0,prop1.lV3d*4*4*sizeof(Complex));
    test.real() = 0.; test.imag() = 0.;
    for(int gamma = 0 ; gamma < 4 ; gamma++)       
      for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++)
	for(int idx = 0 ; idx < counter ; idx++){

	  unsigned short int alpha = Cg5Cg5bar_ind[idx][0];
	  unsigned short int beta = Cg5Cg5bar_ind[idx][1];
	  unsigned short int beta1 = Cg5Cg5bar_ind[idx][2];
	  unsigned short int alpha1 = Cg5Cg5bar_ind[idx][3];

	  for(int cc1 = 0 ; cc1 < 6 ; cc1++){
	    a = eps[cc1][0];
	    b = eps[cc1][1];
	    c = eps[cc1][2];
	    for(int cc2 = 0 ; cc2 < 6 ; cc2++){
	      a1 = eps[cc2][0];
	      b1 = eps[cc2][1];
	      c1 = eps[cc2][2];

	      factor = sgn_eps[cc1] * sgn_eps[cc2] * Cg5Cg5bar_val[idx];

	      for(int lz = 0 ; lz < prop1.lL[2] ; lz++)
		for(int ly = 0 ; ly < prop1.lL[1] ; ly++)
		  for(int lx = 0 ; lx < prop1.lL[0] ; lx++){
		    int v3 = LEXIC_ZYX(lz,ly,lx,prop1.lL);
		    int v = LEXIC(it,lz,ly,lx,prop1.lL);
		    
		    block[v3*4*4 + gamma*4 + gamma1] = block[v3*4*4 + gamma*4 + gamma1] + factor*prop2.M[v][beta][beta1][b][b1] * prop1.M[v][alpha][alpha1][a][a1] * prop1.M[v][gamma][gamma1][c][c1];
		    block[v3*4*4 + gamma*4 + gamma1] = block[v3*4*4 + gamma*4 + gamma1] - factor*prop2.M[v][beta][beta1][b][b1] * prop1.M[v][alpha][gamma1][a][c1] * prop1.M[v][gamma][alpha1][c][a1];
		    /*
		    if(v3 == 0 && gamma ==0 && gamma1 ==0){
		      test =test + factor*prop2.M[v][beta][beta1][b][b1] * prop1.M[v][alpha][alpha1][a][a1] * prop1.M[v][gamma][gamma1][c][c1];
		      test2 = factor*prop2.M[v][beta][beta1][b][b1] * prop1.M[v][alpha][alpha1][a][a1] * prop1.M[v][gamma][gamma1][c][c1];
		    }
		    if(prop1.rank == 0 && v3 == 0 && gamma ==0 && gamma1 ==0)printf("%d %+e %+e %+e %+e\n",it,test.real(),test.imag(),test2.real(),test2.imag());
		    */
		  } // volume

	    } // cc2
	  } // cc1
	} // gamma indices

    //    if(prop1.rank == 0 )printf("%+e %+e %+e %+e\n",test.real(),test.imag(),block[0].real(),block[0].imag());

    double tmp;
    memset(corr,0,Nmom*4*4*sizeof(Complex));
    for(int imom = 0 ; imom < Nmom ; imom++){

      for(int gamma = 0 ; gamma < 4 ; gamma++) 
	for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++){

	  for(int lz = 0 ; lz < prop1.lL[2] ; lz++)
	    for(int ly = 0 ; ly < prop1.lL[1] ; ly++)
	      for(int lx = 0 ; lx < prop1.lL[0] ; lx++){	
		int v3 = LEXIC_ZYX(lz,ly,lx,prop1.lL);
		int x = lx + prop1.procPos[0]*prop1.lL[0] - x_src[0];
		int y = ly + prop1.procPos[1]*prop1.lL[1] - x_src[1];
		int z = lz + prop1.procPos[2]*prop1.lL[2] - x_src[2];
		tmp = (((double) momElem[imom][0]*x)/prop1.L[0] + ((double) momElem[imom][1]*y)/prop1.L[1] + ((double) momElem[imom][2]*z)/prop1.L[2])*2*PI;
		Complex C2 (cos(tmp),-sin(tmp));
		corr[imom*4*4 + gamma*4 + gamma1] = corr[imom*4*4 + gamma*4 + gamma1] + block[v3*4*4 + gamma*4 + gamma1] * C2;
	      }
	  //	  MPI_Reduce(&(corr.real()), &(corr2.real()), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	} // close dirac indices
    } // close fourier

    MPI_Reduce(corr, corr2, Nmom*4*4*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(prop1.rank == 0){
      for(int imom = 0 ; imom < Nmom ; imom++)
	for(int gamma = 0 ; gamma < 4 ; gamma++)
	  for(int gamma1 = 0 ; gamma1 < 4 ; gamma1++)
	    fprintf(ptr_out,"%d   %+d %+d %+d   %d %d  \t %+e %+e\n",it,momElem[imom][0],momElem[imom][1],momElem[imom][2],gamma,gamma1,corr2[imom*4*4 + gamma*4 + gamma1].real(),corr2[imom*4*4 + gamma*4 + gamma1].imag());
    }

    t2 = MPI_Wtime();
    if(prop1.rank == 0) printf(" took %f sec \n",t2-t1);
  } // time

  if(prop1.rank == 0)fclose(ptr_out);
  free(block);
  free(corr);
  free(corr2);

}





FieldStrength::FieldStrength(LatticeInfo *latInfo, Gauge &u) : LatticeGeometry(latInfo) {
  size = lV4d;
  M = (Complex(*)[6][NCOLORS][NCOLORS]) safe_malloc(6*NCOLORS*NCOLORS*2*size*sizeof(double));
  zero();
  //============== Serial version of calculation of field tensor =========//
  if (lV4d != V4d)
    errorQKXTM("Error field tensor supports only serial version for now\n");

  Complex tmp1[3][3];
  Complex tmp2[3][3];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < mu; nu++){

      for(int it = 0 ; it < L[3]; it++)
	for(int iz = 0 ; iz < L[2]; iz++)
	  for(int iy = 0 ; iy < L[1]; iy++)
	    for(int ix = 0 ; ix < L[0]; ix++){
	      int iv = LEXIC(it,iz,iy,ix,L);
	      // (mu,nu) -> index,  (1,0) -> 0, (2,0) -> 1, (2,1) -> 2, (3,0) -> 3, (3,1) -> 4, (3,2) -> 5
	      int munu = (mu*(mu-1))/2+nu; 
	      int itmp;
	      int xtmp[4];
	      int PnuMmu;
	      int MnuMmu;
	      int MnuPmu;
	      //== U_\mu(x) * U_\nu(x+\mu) * U^dag_\mu(x+\nu) * U^dag_\nu(x) ==//
	      qcd_MUL3x3(tmp1,u.M[iv][mu],u.M[ getNeighborPlus(it,iz,iy,ix,L,mu) ][nu]);
	      qcd_MULADJOINT3x3(tmp2,tmp1,u.M[ getNeighborPlus(it,iz,iy,ix,L,nu) ][mu]);
	      qcd_MULADJOINT3x3(tmp1,tmp2,u.M[iv][nu]);
	      qcd_ACCUM_3x3(this->M[iv][munu],tmp1);
	      //== U_\nu(x) * U^dag_\mu(x+\nu-\mu) * U^\dag_\nu(x-\mu) * U_\mu(x-\mu) ===//
	      itmp=getNeighborPlus(it,iz,iy,ix,L,nu);
	      getCoord(itmp,xtmp,L);
	      PnuMmu = getNeighborMinus(xtmp[3],xtmp[2],xtmp[1],xtmp[0],L,mu);
	      qcd_MULADJOINT3x3(tmp1,u.M[iv][nu],u.M[PnuMmu][mu]);
	      qcd_MULADJOINT3x3(tmp2,tmp1,u.M[ getNeighborMinus(it,iz,iy,ix,L,mu) ][nu]);
	      qcd_MUL3x3(tmp1,tmp2,u.M[getNeighborMinus(it,iz,iy,ix,L,mu)][mu]);
	      qcd_ACCUM_3x3(this->M[iv][munu],tmp1);
	      //== U^\dag_\mu(x-\mu) * U^dag_\nu(x-\nu-\mu) * U_\mu(x-nu-mu) * U_\nu(x-\nu) ==//
	      // here we use that (U_\nu(x-\nu-\mu) * U_\mu(x-\mu))^dag = U^\dag_\mu(x-\mu) * U^dag_\nu(x-\nu-\mu)
	      itmp=getNeighborMinus(it,iz,iy,ix,L,nu);
	      getCoord(itmp,xtmp,L);
	      MnuMmu=getNeighborMinus(xtmp[3],xtmp[2],xtmp[1],xtmp[0],L,mu);
	      qcd_MUL3x3(tmp1,u.M[MnuMmu][nu],u.M[getNeighborMinus(it,iz,iy,ix,L,mu)][mu]);
	      qcd_ADJOINTMUL3x3(tmp2,tmp1,u.M[MnuMmu][mu]);
	      qcd_MUL3x3(tmp1,tmp2,u.M[getNeighborMinus(it,iz,iy,ix,L,nu)][nu]);
	      qcd_ACCUM_3x3(this->M[iv][munu],tmp1);
	      //=== U^dag_\nu(x-\nu) * U_\mu(x-\nu) * U_\nu(x-\nu+\mu) * U^dag_\mu(x)
	      itmp=getNeighborMinus(it,iz,iy,ix,L,nu);
              getCoord(itmp,xtmp,L);
	      MnuPmu=getNeighborPlus(xtmp[3],xtmp[2],xtmp[1],xtmp[0],L,mu);
	      qcd_ADJOINTMUL3x3(tmp1,u.M[getNeighborMinus(it,iz,iy,ix,L,nu)][nu], u.M[getNeighborMinus(it,iz,iy,ix,L,nu)][mu]);
	      qcd_MUL3x3(tmp2,tmp1,u.M[MnuPmu][nu]);
	      qcd_MULADJOINT3x3(tmp1,tmp2,u.M[iv][mu]);
	      qcd_ACCUM_3x3(this->M[iv][munu],tmp1);

	      // calculate the h.c
	      for(int c1 = 0 ; c1 < 3 ; c1++)
		for(int c2 = 0 ; c2 < 3 ; c2++)
		  tmp1[c1][c2] = conj(this->M[iv][munu][c2][c1]);
	      
	      // complete the field strength tensor
	      for(int c1 = 0 ; c1 < 3 ; c1++)
		for(int c2 = 0 ; c2 < 3 ; c2++)
		  this->M[iv][munu][c1][c2] = (1./8.) * (this->M[iv][munu][c1][c2] - tmp1[c1][c2]);
	      
	    }
    }
  
}

FieldStrength::~FieldStrength(){
  free(M);
}

void FieldStrength::zero(){
  memset(&(M[0][0][0][0].real()) , 0 , 6*NCOLORS*NCOLORS*2*size*sizeof(double) );
}

