#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  latInfo.L[0] = 8;
  latInfo.L[1] = 8;
  latInfo.L[2] = 8;
  latInfo.L[3] = 8;

  latInfo.P[0] = 2;
  latInfo.P[1] = 1;
  latInfo.P[2] = 2;
  latInfo.P[3] = 2;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                


  LatticeGeometry *geo = new LatticeGeometry(&latInfo);

  printfQKXTM("rank is %d\n",geo->rank);
  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",geo->L[0],geo->L[1],geo->L[2] , geo->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",geo->lL[0],geo->lL[1],geo->lL[2] , geo->lL[3]);
  printfQKXTM("global volume = %d\n",geo->V4d);
  printfQKXTM("local volume = %d\n",geo->lV4d);

  //  for(int i = 0 ; i < NDIM ; i++)
  // printf("%d %d %d\n",geo->rank, geo->procPlus[i] , geo->procMinus[i]);
  int x[4];
  int xp[4];

  for(int dir = 0 ; dir < 4 ; dir++)
  for(int it = 0 ; it < geo->lL[3] ; it++)
    for(int iz = 0 ; iz < geo->lL[2] ; iz++)
      for(int iy = 0 ; iy < geo->lL[1] ; iy++)
	for(int ix = 0 ; ix < geo->lL[0] ; ix++){
	  int iv = LEXIC(it,iz,iy,ix,geo->lL);
	  getCoord(iv,x,geo->lL);
	  getCoord(geo->pointMinus[dir][iv],xp,geo->lL);

	  if(geo->pointMinus[dir][iv] < geo->lV4d){
	    if(rank == 0)printf("%d (%d,%d,%d,%d) , (%d,%d,%d,%d)\n",iv,x[0],x[1],x[2],x[3],xp[0],xp[1],xp[2],xp[3]);
	  }

	  else{
	    if(rank == 0)printf("%d (%d,%d,%d,%d) %d out \n",iv,x[0],x[1],x[2],x[3],geo->pointMinus[dir][iv]);
	  }
	}

  /*
  for(int i = 0 ; i < geo->lV4d ; i++){
    getCoord(i,x,geo->L);
    getCoord(geo->pointPlus[1][i],xp,geo->L);

    if(geo->pointPlus[1][i] < geo->lV4d){
      if(rank == 0)printf("%d (%d,%d,%d,%d) , (%d,%d,%d,%d)\n",i,x[0],x[1],x[2],x[3],xp[0],xp[1],xp[2],xp[3]);
    }

    else{
      if(rank == 0)printf("%d (%d,%d,%d,%d) %d out \n",i,x[0],x[1],x[2],x[3],geo->pointPlus[1][i]);
    }

  }
  */

  MPI_Finalize();

}
