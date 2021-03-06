#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 6;
  latInfo.L[1] = 6;
  latInfo.L[2] = 6;
  latInfo.L[3] = 6;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.156;
  latInfo.mu = 0.1;
  latInfo.twistSign = -1;
  latInfo.tol = 1e-6;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = +1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 0;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 0;
  latInfo.alphaGauss = 4.0;


  char file_conf[] = "/users/krikitos/scratch/16x32/eigenVectors/conf.0000";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);
  Vector *vec = new Vector(&latInfo);
  Vector *eigenVec = new Vector(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();
  qkxTM_MPI_getGaugeLime(file_conf,*gauge);
  double plaq = gauge->calculatePlaquette();
  printfQKXTM("Plaquette is %e\n",plaq);
  gauge->applyBoundaryCondition(&latInfo);

  char file_eigen[] = "/users/krikitos/scratch/16x32/eigenVectors/theta0/ev.0000.00000";
  //char file_eigen[] = "/users/krikitos/scratch/16x32/eigenVectors/eigenVectors_tmLQCD/twisted_mass/ev.0000.00000";
  qkxTM_MPI_readEigenVectors(file_eigen,*eigenVec);
  eigenVec->sendEvenToOdd();
  eigenVec->rotateFromChiralToUKQCD();
  //  eigenVec->multiply_by_phase();
  //  eigenVec->normalize();
  
  //  vec->tm_xxDagTm_xx(*eigenVec,*gauge,&latInfo,OO);


  for(int iv = 0 ; iv < eigenVec->lV4d; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int c1 = 0 ; c1 < 3 ; c1++)
	printf("%+e %+e\n",eigenVec->M[iv][mu][c1].real(),eigenVec->M[iv][mu][c1].imag());

  delete gauge;
  delete eigenVec;
  delete vec;
  MPI_Finalize();

}
