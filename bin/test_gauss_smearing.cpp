#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 16;
  latInfo.L[1] = 16;
  latInfo.L[2] = 16;
  latInfo.L[3] = 32;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-6;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = +1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 20;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 50;
  latInfo.alphaGauss = 4.0;


  char file_conf[] = "/home/khadjiyiannakou/confs/L16T32/conf.1505";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);
  Gauge *gauge_tmp = new Gauge(&latInfo);
  Gauge *gauge_ape = new Gauge(&latInfo);

  Vector *vec = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);

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
  delete gauge_tmp; // dont need it anymore
  end_time = MPI_Wtime();
  if(rank == 0)printf("finish smearing in %e secs\n",end_time-start_time);
  plaq = gauge_ape->calculatePlaquette();
  
  if(rank == 0)printf("plaquette smeared is %8.7lf \n",plaq);

  vec_tmp->zero();

  if(rank == 0)  vec_tmp->M[0][0][0].real() = 1.;

  vec->Gaussian_smearing(*vec_tmp,*gauge_ape,&latInfo);
  delete vec_tmp;

  /*
  if(rank == 0){
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	printf("%e %e\n",vec->M[0][mu][ic].real(),vec->M[0][mu][ic].imag());

  }
  */
  double vec_norm;
  vec_norm = vec->norm2();
  if(rank == 0) printf("the norm of vector is %e\n",vec_norm);

   if(rank == 0){
     FILE *ptr;
     ptr = fopen("kale.dat","w");
     for(int z = 0 ; z < latInfo.L[2] ; z++)
       for(int y = 0 ; y < latInfo.L[1] ; y++)
	 for(int x = 0 ; x < latInfo.L[0] ; x++)
	   for(int mu = 0 ; mu < 4 ; mu++)
	     for(int c1 = 0 ; c1 < 3 ; c1++)
	       fprintf(ptr,"%d %d %d %+e %+e\n",x,y,z,vec->M[LEXIC_ZYX(z,y,x,latInfo.L)][mu][c1].real(),vec->M[LEXIC_ZYX(z,y,x,latInfo.L)][mu][c1].imag());
   }

   //   printf("from rank = %d , first gauge element is %+e,%+e\n",rank,gauge->M[0][0][0][0].real(),gauge->M[0][0][0][0].imag());

   delete gauge;
   delete gauge_ape;
   delete vec;
  MPI_Finalize();

}
