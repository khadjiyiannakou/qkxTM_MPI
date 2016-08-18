#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 32;
  latInfo.L[1] = 32;
  latInfo.L[2] = 32;
  latInfo.L[3] = 64;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 8;

  latInfo.kappa = 0.161236;
  latInfo.mu = 0.005500;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-6;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = +1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 0;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 0;
  latInfo.alphaGauss = 4.0;


  char file_conf[] = "/home/khadjiyiannakou/confs/L32T64/conf.1000";
  char file_out[] = "/home/khadjiyiannakou/qkxTM_MPI/bin/gluonLoop.1000.dat";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;

  Complex *gluonLoop = (Complex*)calloc(6*latInfo.L[3],sizeof(Complex));

  Gauge *gauge = new Gauge(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

   gauge->zero();

   qkxTM_MPI_getGaugeLime(file_conf,*gauge);
   double plaq = gauge->calculatePlaquette();


   if(rank == 0)printf("plaquette is %8.7lf \n",plaq);
   gauge->calculateGluonLoop(gluonLoop);


   FILE *ptr_out = NULL;
   if(rank == 0){
     ptr_out = fopen(file_out,"w");
     for(int i = 0 ; i < 6 ; i++)
       for(int it = 0 ; it < latInfo.L[3] ; it++)
	 fprintf(ptr_out,"%d %d %+.10e %+.10e\n",it,i,gluonLoop[i*latInfo.L[3]+it].real(),gluonLoop[i*latInfo.L[3]+it].imag());
     fclose(ptr_out);
   }

   free(gluonLoop);

   delete gauge;

  MPI_Finalize();

}
