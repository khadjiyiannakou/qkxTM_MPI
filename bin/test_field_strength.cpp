#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  latInfo.L[0] = 16;
  latInfo.L[1] = 16;
  latInfo.L[2] = 16;
  latInfo.L[3] = 32;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.160856;
  latInfo.mu = 0.00400;
  latInfo.twistSign = -1;
  latInfo.tol = 1e-12;
  latInfo.maxIter = 50000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = -1;

  latInfo.NsmearAPE = 0;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 0;
  latInfo.alphaGauss = 4.0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;

  char file_conf[257] = "/users/krikitos/scratch/test_cov_der/16x32/conf.1505";
  char file_out[257] = "/users/krikitos/scratch/test_cov_der/16x32/field_strength.dat";
  Gauge *gauge = new Gauge(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();

   // start_time = MPI_Wtime();
  qkxTM_MPI_getGaugeLime(file_conf,*gauge);
  double plaq = gauge->calculatePlaquette();
  if(rank == 0)printf("plaquette is %8.7lf \n",plaq);

  FILE *out_ptr = fopen(file_out,"w");
  if(out_ptr == NULL){
    fprintf(stderr,"Error open output file\n");
    exit(-1);
  }

  FieldStrength *fs1 = new FieldStrength(&latInfo,*gauge);
  gauge->applyBoundaryCondition(&latInfo);
  FieldStrength *fs2 = new FieldStrength(&latInfo,*gauge);

  for(int iv = 0 ; iv < gauge->V4d; iv++)
    for(int dir = 0 ; dir < 6 ; dir++)
      for(int c1 = 0 ; c1 < 3 ; c1++)
	for(int c2 = 0 ; c2 < 3 ; c2++)
	  fprintf(out_ptr,"%+e %+e \t %+e %+e\n",fs1->M[iv][dir][c1][c2].real(), fs1->M[iv][dir][c1][c2].imag(),
		  fs2->M[iv][dir][c1][c2].real(), fs2->M[iv][dir][c1][c2].imag());
  gauge->constGauge = true;

  delete gauge;
  delete fs1;
  delete fs2;
  MPI_Finalize();

}
