#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 24;
  latInfo.L[1] = 24;
  latInfo.L[2] = 24;
  latInfo.L[3] = 48;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 8;

  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = -1;
  latInfo.tol = 1e-12;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = -1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 20;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 50;
  latInfo.alphaGauss = 4.0;


  char file_conf[] = "/home/khadjiyiannakou/qudaQKXTM_package/confs/L24T48/conf.1000";


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);
  Gauge *gauge_tmp = new Gauge(&latInfo);
  Gauge *gauge_ape = new Gauge(&latInfo);

  Vector *vec_out = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);
  Vector *vec_in = new Vector(&latInfo);

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
  vec_out->zero();
  for(int mu = 0 ; mu < NSPINS ; mu++)
    for(int ic =0 ; ic < NCOLORS ; ic++)
      if(rank == 0)  vec_tmp->M[0][mu][ic].real() = 1.;

  vec_out->Gaussian_smearing(*vec_tmp,*gauge_ape,&latInfo);
  //vec_out->copy(*vec_tmp);
  
  double vec_norm;
  vec_in->tmDag(*vec_out,*gauge,&latInfo);
  vec_out->zero();
  //  vec_norm = vec_in->norm2();
  // if(rank == 0)printf("the norm of vector is %e\n",vec_norm);
  //  MPI_Abort();
  //  vec_in->applyDslashDag(*vec_tmp, *gauge);
  // vec_in->applyTwistDagAddDslashDag(*vec_in,*vec_tmp,&latInfo);
  
  //delete vec_tmp;
  


  // apply boundary conditions on gauge

  gauge->applyBoundaryCondition(&latInfo);
  gauge->constGauge = true;



  CG *cg_P = new CG(&latInfo,*gauge);
  CG &cg = *cg_P;
  //  vec_tmp->zero();
  cg(*vec_out,*vec_in,&latInfo);

  vec_norm = vec_out->norm2();
  if(rank == 0)printf("the norm of vector is %e\n",vec_norm);

  //  vec_tmp->tmDag(*vec,*gauge,&latInfo);
  //vec_norm = vec_tmp->norm2();
  //if(rank == 0) printf("the norm of vector is %e\n",vec_norm);
  
//  vec_norm = vec->norm2();
  // if(rank == 0) printf("the norm of vector is %e\n",vec_norm);
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

  delete vec_tmp;
  delete gauge;
  delete gauge_ape;
  delete vec_out;
  MPI_Finalize();

}
