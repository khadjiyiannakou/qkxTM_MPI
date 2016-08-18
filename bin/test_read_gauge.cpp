#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 24;
  latInfo.L[1] = 24;
  latInfo.L[2] = 24;
  latInfo.L[3] = 48;

  latInfo.P[0] = 2;
  latInfo.P[1] = 2;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  char file_conf[] = "/home/krikitos3/Desktop/QCD/qkxTM/confs/conf.1000";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);


  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

   gauge->zero();

   start_time = MPI_Wtime();
   qkxTM_MPI_getGaugeLime(file_conf,*gauge);
   end_time = MPI_Wtime();
   if(rank == 0)printf("read conf in time is %e in sec\n",end_time-start_time);

   start_time = MPI_Wtime();
   double plaq = gauge->calculatePlaquette();
   end_time = MPI_Wtime();


   if(rank == 0)printf("plaquette is %8.7lf and time is %e in sec\n",plaq,end_time-start_time);

   /*      
   if(rank == 1){
     FILE *ptr;
     ptr = fopen("kale.dat","w");
     for(int iv = 0 ; iv < gauge->lV4d ; iv++)
       for(int mu = 0 ; mu < 4 ; mu++)
	 for(int c1 = 0 ; c1 < 3 ; c1++)
	   for(int c2 = 0 ; c2 < 3 ; c2++)
	     fprintf(ptr,"%+e %+e\n",gauge->M[iv][mu][c1][c2].real(),gauge->M[iv][mu][c1][c2].imag());
   }
   */
   //   printf("from rank = %d , first gauge element is %+e,%+e\n",rank,gauge->M[0][0][0][0].real(),gauge->M[0][0][0][0].imag());

   delete gauge;

  MPI_Finalize();

}
