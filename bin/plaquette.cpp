#include <qkxTM.h>
#include <lattice_util.h>
#include <string.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  if(argc != 3)errorQKXTM("wrong number of inputs\n");

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 48;
  latInfo.L[1] = 48;
  latInfo.L[2] = 48;
  latInfo.L[3] = 96;

  latInfo.P[0] = 2;
  latInfo.P[1] = 4;
  latInfo.P[2] = 6;
  latInfo.P[3] = 1;

  latInfo.NsmearAPE = 50;
  latInfo.alphaAPE = 0.5;


  char file_conf[257];
  char file_conf_smeared[257]; 
  strcpy(file_conf,argv[1]);
  strcpy(file_conf_smeared,argv[2]);
  //  char message[] = "Smeared configuration using APE smearing 3D, alpa = 0.4 , nsmear = 50 , kappa = 0.156361 , mu = 0.001500";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                

  if(rank == 0){
    //printf("read from %s\n",file_conf);
    // printf("write on %s\n",file_conf_smeared);
  }


  Gauge *gauge = new Gauge(&latInfo);
  Gauge *gauge_tmp = new Gauge(&latInfo);
  Gauge *gauge_ape = new Gauge(&latInfo);

  //  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  //  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

   gauge->zero();

   // start_time = MPI_Wtime();
   qkxTM_MPI_getGaugeLime(file_conf,*gauge);
   //end_time = MPI_Wtime();
   // if(rank == 0)printf("read conf in time is %e in sec\n",end_time-start_time);

   gauge_tmp->copy(*gauge);

   // start_time = MPI_Wtime();
   double plaq = gauge->calculatePlaquette();
   //   end_time = MPI_Wtime();

   if(rank == 0)printf("%s , plaquette is %8.7lf \n",file_conf,plaq);

   
   delete gauge;

  MPI_Finalize();

}
