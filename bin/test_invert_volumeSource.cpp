#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  //  double start_time;
  //double end_time;

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
  latInfo.boundaryT = -1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 0;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 0;
  latInfo.alphaGauss = 4.0;




  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;

  char file_conf[257] = "/users/krikitos/run/test_loop/conf.1505";
  char file_xi[257] = "/users/krikitos/run/test_loop/volumeSource_kx.In";
  char file_phi[257] = "/users/krikitos/run/test_loop/volumeSource_kx.Out";
  char file_ultra_local[257] ="/users/krikitos/run/test_loop/ultralocal.dat";
  char file_oneD[257] ="/users/krikitos/run/test_loop/oneD.dat";
  char file_noe[257] ="/users/krikitos/run/test_loop/noe.dat";

  Gauge *gauge = new Gauge(&latInfo);

  Vector *vec_phi = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);
  Vector *vec_xi = new Vector(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();

   // start_time = MPI_Wtime();
  qkxTM_MPI_getGaugeLime(file_conf,*gauge);
  double plaq = gauge->calculatePlaquette();
  if(rank == 0)printf("plaquette is %8.7lf \n",plaq);

  
  FILE *ptr_xi = NULL;
  FILE *ptr_phi = NULL;
  FILE *ptr_ultra_local = NULL;
  FILE *ptr_oneD = NULL;
  FILE *ptr_noe = NULL;

  ptr_xi = fopen(file_xi,"r");
  ptr_phi = fopen(file_phi,"r");
  ptr_ultra_local = fopen(file_ultra_local,"w");
  ptr_oneD = fopen(file_oneD,"w");
  ptr_noe = fopen(file_noe,"w");

  //  int ierr;
  if(ptr_xi == NULL || ptr_phi == NULL){
    fprintf(stderr,"Error open input files\n");
    exit(-1);
  }

  if(rank == 0){
    if(ptr_ultra_local == NULL || ptr_oneD == NULL || ptr_noe == NULL){
      fprintf(stderr,"Error open output file\n");
      exit(-1);
    }
  }
//MPI_Abort(MPI_COMM_WORLD,ierr);
 
  int local_volume = gauge->lL[0] * gauge->lL[1] * gauge->lL[2] * gauge->lL[3];

  MPI_Barrier(MPI_COMM_WORLD);

  gauge->applyBoundaryCondition(&latInfo);
  gauge->constGauge = true;
  // apply boundary conditions on gauge


  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(ptr_xi,"%lf %lf",&(vec_xi->M[iv][mu][ic].real()) , &(vec_xi->M[iv][mu][ic].imag()) );

  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(ptr_phi,"%lf %lf",&(vec_phi->M[iv][mu][ic].real()) , &(vec_phi->M[iv][mu][ic].imag()) );


  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;

  Complex *loops_ultra_local = (Complex*)malloc(latInfo.L[3]*4*4*7*sizeof(Complex));
  Complex *loops_oneD = (Complex*)malloc(latInfo.L[3]*4*4*4*7*sizeof(Complex));
  Complex *loops_noe = (Complex*)malloc(latInfo.L[3]*4*4*4*7*sizeof(Complex));

  vec_phi->volumeSource_ultralocal(*vec_xi,*gauge,loops_ultra_local, &latInfo );
  vec_phi->volumeSource_covDer(*vec_xi,*gauge,loops_oneD,loops_noe,&latInfo);

  if(rank == 0){

    for(int imom = 0 ; imom < 7 ; imom++)
      for(int it = 0 ; it < latInfo.L[3] ; it++)
	for(int ioper = 0 ; ioper < 16 ; ioper++)
	  fprintf(ptr_ultra_local,"%d %d %d %+d %+d %+d %+e %+e\n",1,it,ioper,momList[0][imom],momList[1][imom],momList[2][imom],loops_ultra_local[it*16*7+ioper*7+imom].real() , loops_ultra_local[it*16*7+ioper*7+imom].imag() );

    for(int imom = 0 ; imom < 7 ; imom++)
      for(int it = 0 ; it < latInfo.L[3] ; it++)
	for(int ioper = 0 ; ioper < 16 ; ioper++)
	  for(int nu = 0 ; nu < 4 ; nu++)
	    fprintf(ptr_oneD,"%d %d %d %d %+d %+d %+d %+e %+e\n",1,it,nu,ioper,momList[0][imom],momList[1][imom],momList[2][imom],loops_oneD[it*4*16*7+nu*16*7+ioper*7+imom].real() , loops_oneD[it*4*16*7+nu*16*7+ioper*7+imom].imag() );

    for(int imom = 0 ; imom < 7 ; imom++)
      for(int it = 0 ; it < latInfo.L[3] ; it++)
	for(int ioper = 0 ; ioper < 16 ; ioper++)
	  for(int nu = 0 ; nu < 4 ; nu++)
	    fprintf(ptr_noe,"%d %d %d %d %+d %+d %+d %+e %+e\n",1,it,nu,ioper,momList[0][imom],momList[1][imom],momList[2][imom],loops_noe[it*4*16*7+nu*16*7+ioper*7+imom].real() , loops_noe[it*4*16*7+nu*16*7+ioper*7+imom].imag() );
    
  }


  free(loops_ultra_local);
  free(loops_oneD);
  free(loops_noe);

  delete vec_xi;
  delete vec_tmp;
  delete gauge;
  delete vec_phi;
  MPI_Finalize();

}
