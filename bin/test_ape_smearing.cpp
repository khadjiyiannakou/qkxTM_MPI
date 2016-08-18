#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 48;
  latInfo.L[1] = 48;
  latInfo.L[2] = 48;
  latInfo.L[3] = 96;

  latInfo.P[0] = 4;
  latInfo.P[1] = 2;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-6;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = +1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 50;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 10;
  latInfo.alphaGauss = 4.0;


  char file_conf[] = "/users/krikitos/scratch/cA2.09.48/confs/conf.2042";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);
  Gauge *gauge_tmp = new Gauge(&latInfo);
  Gauge *gauge_ape = new Gauge(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

   gauge->zero();

   // start_time = MPI_Wtime();
   qkxTM_MPI_getGaugeLime(file_conf,*gauge);
   //end_time = MPI_Wtime();
   // if(rank == 0)printf("read conf in time is %e in sec\n",end_time-start_time);

   gauge_tmp->copy(*gauge);

   // start_time = MPI_Wtime();
   double plaq = gauge->calculatePlaquette();
   //   end_time = MPI_Wtime();

   if(rank == 0)printf("plaquette is %8.7lf \n",plaq);

   start_time = MPI_Wtime();
   gauge_ape->APE_smearing(*gauge_tmp,&latInfo);
   end_time = MPI_Wtime();
   if(rank == 0)printf("finish smearing in %e secs\n",end_time-start_time);
   plaq = gauge_ape->calculatePlaquette();

   if(rank == 0)printf("plaquette smeared is %8.7lf \n",plaq);

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
