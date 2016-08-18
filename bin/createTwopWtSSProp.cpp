#include <qkxTM.h>
#include <lattice_util.h>
#include <string.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  if(argc != 9)errorQKXTM("wrong number on inputs\n");

  LatticeInfo latInfo;

  double start_time;
  double end_time;


  latInfo.L[0] = 48;
  latInfo.L[1] = 48;
  latInfo.L[2] = 48;
  latInfo.L[3] = 96;

  latInfo.P[0] = 3;
  latInfo.P[1] = 4;
  latInfo.P[2] = 8;
  latInfo.P[3] = 1;  

  /*        
  latInfo.L[0] = 16;
  latInfo.L[1] = 16;
  latInfo.L[2] = 16;
  latInfo.L[3] = 32;

  latInfo.P[0] = 2;
  latInfo.P[1] = 2;
  latInfo.P[2] = 2;
  latInfo.P[3] = 1;  
  */

  //  char file_prop_up[] = "/home/khadjiyiannakou/test_write_props/prop_SS.up";
  //  char file_prop_down[] = "/home/khadjiyiannakou/test_write_props/prop_SS.down";

  char file_prop_up[257];
  char file_prop_down[257];

  strcpy(file_prop_up,argv[1]);
  strcpy(file_prop_down,argv[2]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                

  Propagator *prop_up = new Propagator(&latInfo);
  Propagator *prop_down = new Propagator(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",prop_up->L[0],prop_up->L[1],prop_up->L[2] , prop_up->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",prop_up->lL[0],prop_up->lL[1],prop_up->lL[2] , prop_up->lL[3]);
  printfQKXTM("Got %d process\n",latInfo.Nprocs);

   start_time = MPI_Wtime();
   qkxTM_MPI_getPropLime(file_prop_up,*prop_up,&latInfo);
   qkxTM_MPI_getPropLime(file_prop_down,*prop_down,&latInfo);
   end_time = MPI_Wtime();
   if(rank == 0)printf("read props in time is %e in sec\n",end_time-start_time);
   

   prop_up->rotateToPhysicalBasePlus();
   prop_down->rotateToPhysicalBaseMinus();

   
   char twop_proton[257] ;//= "/home/khadjiyiannakou/qkxTM_MPI/bin/twop_proton.dat";
   char twop_neutron[257] ;// = "/home/khadjiyiannakou/qkxTM_MPI/bin/twop_neutron.dat";

   strcpy(twop_proton,argv[3]);
   strcpy(twop_neutron,argv[4]);

   int x_src[4];
   x_src[0] = atoi(argv[5]);
   x_src[1] = atoi(argv[6]);;
   x_src[2] = atoi(argv[7]);;
   x_src[3] = atoi(argv[8]);

   contractNucleon(*prop_up,*prop_down,twop_proton,16,x_src);
   contractNucleon(*prop_down,*prop_up,twop_neutron,16,x_src);
   

   delete prop_up;
   delete prop_down;

  MPI_Finalize();

}
