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
  latInfo.P[1] = 2;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                


  LatticeGeometry *geo = new LatticeGeometry(&latInfo);
  Vector *vec = new Vector(&latInfo);
  Propagator *prop = new Propagator(&latInfo);
  Gauge *gauge = new Gauge(&latInfo);


  //printfQKXTM("rank is %d\n",geo->rank);
  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",geo->L[0],geo->L[1],geo->L[2] , geo->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",geo->lL[0],geo->lL[1],geo->lL[2] , geo->lL[3]);
  // printfQKXTM("global volume = %d\n",geo->V4d);
  //printfQKXTM("local volume = %d\n",geo->lV4d);
  //  printf("rank is %d (%d,%d,%d,%d) local procPlus is (%d,%d,%d,%d) \n",rank,geo->procPos[0],geo->procPos[1],geo->procPos[2],geo->procPos[3], geo->procPlus[0],geo->procPlus[1],geo->procPlus[2],geo->procPlus[3]);
  //  printf("rank is %d (%d,%d,%d,%d) local procPlus is (%d,%d,%d,%d) \n",rank,geo->procPlus[0],geo->procPlus[1],geo->procPlus[2],geo->procPlus[3]);
  //  printfQKXTM("size is %d and total size is %d\n",vec->size,vec->size_total);

   vec->zero();
   prop->zero();
   gauge->zero();

   vec->communicateToPlus();
   vec->communicateToMinus();
    
   prop->communicateToPlus();
   prop->communicateToMinus();

   gauge->communicateToPlus();
   gauge->communicateToMinus();
  

   delete geo;
   delete vec;
   delete prop;
   delete gauge;

  MPI_Finalize();

}
