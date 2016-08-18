#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  latInfo.L[0] = 48;
  latInfo.L[1] = 48;
  latInfo.L[2] = 48;
  latInfo.L[3] = 96;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.137290;
  latInfo.mu = 0.000900;
  latInfo.twistSign = -1;
  latInfo.tol = 1e-12;
  latInfo.maxIter = 50000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = -1;                   // remember write program that applies boundaries
  latInfo.csw = 1.57551;

  latInfo.NsmearAPE = 0;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 0;
  latInfo.alphaGauss = 4.0;




  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;

  char file_conf[257] = "/homec/ecy00/ecy011/work/test_cov_der/conf.1202";
  char file_in[257] = "/homec/ecy00/ecy011/work/test_cov_der/inverted_source.1";
  char file_out_st[257] = "/homec/ecy00/ecy011/work/test_cov_der/sigmaTerm.dat"; // give sigma term
  char file_out_ax[257] = "/homec/ecy00/ecy011/work/test_cov_der/generalizedLoops.dat"; 
  char file_out_gD[257] = "/homec/ecy00/ecy011/work/test_cov_der/generalizedDerLoops.dat"; 
  char file_out_cD[257] = "/homec/ecy00/ecy011/work/test_cov_der/generalizedConservedLoops.dat"; 

  sprintf(file_in,"%s_%d",file_in,rank);

  Gauge *gauge = new Gauge(&latInfo);

  Vector *vec_out = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);
  Vector *vec_in = new Vector(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();


   // start_time = MPI_Wtime();
  qkxTM_MPI_getGaugeLime(file_conf,*gauge);
  double plaq = gauge->calculatePlaquette();
  if(rank == 0)printf("plaquette is %8.7lf \n",plaq);
  fflush(stdout);
  
  FILE *in_file = NULL;
  FILE *out_file_st = NULL;
  FILE *out_file_ax = NULL;
  FILE *out_file_gD = NULL;
  FILE *out_file_cD = NULL;

  in_file = fopen(file_in,"r");
  printf("%d %s\n",rank,file_in);

  if(rank == 0){
    //    out_file_st = fopen(file_out_st,"w");
    //    out_file_ax = fopen(file_out_ax,"w");
    out_file_gD = fopen(file_out_gD,"w");
    out_file_cD = fopen(file_out_cD,"w");
  }

  //  int ierr;
  if(in_file == NULL){
    fprintf(stderr,"Error open input file\n");
    exit(-1);
  }

  if(rank == 0){
    if(out_file_gD == NULL || out_file_cD == NULL){
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

  /*
  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(in_file,"%lf %lf",&(vec_tmp->M[iv][mu][ic].real()) , &(vec_tmp->M[iv][mu][ic].imag()) );


  vec_in->tmDag(*vec_tmp,*gauge,&latInfo);
  vec_out->zero();



  CG *cg_P = new CG(&latInfo,*gauge);
  CG &cg = *cg_P;
  cg(*vec_out,*vec_in,&latInfo);
  */

  int dummy;
  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(in_file,"%d %d %d %d %d %d %lf %lf",&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&(vec_out->M[iv][mu][ic].real()) , &(vec_out->M[iv][mu][ic].imag()) );


  int momList[3][7];
  momList[0][0] = 0; momList[1][0] = 0 ; momList[2][0] = 0;
  momList[0][1] = 1; momList[1][1] = 0 ; momList[2][1] = 0;
  momList[0][2] = -1; momList[1][2] = 0 ; momList[2][2] = 0;
  momList[0][3] = 0; momList[1][3] = 1 ; momList[2][3] = 0;
  momList[0][4] = 0; momList[1][4] = -1 ; momList[2][4] = 0;
  momList[0][5] = 0; momList[1][5] = 0 ; momList[2][5] = 1;
  momList[0][6] = 0; momList[1][6] = 0 ; momList[2][6] = -1;

  Complex *loops_st = (Complex*)malloc(latInfo.L[3]*7*sizeof(Complex));
  Complex *loops_ax = (Complex*)malloc(latInfo.L[3]*4*4*7*sizeof(Complex));
  Complex *loops_gD = (Complex*)malloc(latInfo.L[3]*4*4*4*7*sizeof(Complex));
  Complex *loops_cD = (Complex*)malloc(latInfo.L[3]*4*4*4*7*sizeof(Complex));

  //  vec_out->standardEndTrick_sigmaTerm(loops_st, &latInfo );
  //vec_out->generalizedEndTrick_axial_openIndices(*gauge,loops_ax, &latInfo );
  vec_out->generalizedEndTrick_covDer(*gauge,loops_gD,loops_cD,&latInfo);

  if(rank == 0){

//     for(int it = 0 ; it < latInfo.L[3] ; it++)
//       for(int imom = 0 ; imom < 7 ; imom++)
// 	fprintf(out_file_st,"%d %+d %+d %+d %+e %+e\n",it,momList[0][imom],momList[1][imom],momList[2][imom],loops_st[it*7+imom].real() , loops_st[it*7+imom].imag() );

    

//     for(int imom = 0 ; imom < 7 ; imom++)
//       for(int it = 0 ; it < latInfo.L[3] ; it++)
// 	for(int ioper = 0 ; ioper < 16 ; ioper++)
// 	  fprintf(out_file_ax,"%d %d %d %+d %+d %+d %+e %+e\n",1,it,ioper,momList[0][imom],momList[1][imom],momList[2][imom],loops_ax[it*16*7+ioper*7+imom].real() , loops_ax[it*16*7+ioper*7+imom].imag() );

    for(int nu = 0 ; nu < 4 ; nu++)
      for(int imom = 0 ; imom < 7 ; imom++)
	for(int it = 0 ; it < latInfo.L[3] ; it++)
	  for(int ioper = 0 ; ioper < 16 ; ioper++)
	    fprintf(out_file_gD,"%02d %02d %02d %+d %+d %+d %+16.15e %+16.15e\n",it,ioper,nu,momList[0][imom],momList[1][imom],momList[2][imom],loops_gD[it*4*16*7+nu*16*7+ioper*7+imom].real() , loops_gD[it*4*16*7+nu*16*7+ioper*7+imom].imag() );

    for(int nu = 0 ; nu < 4 ; nu++)
      for(int imom = 0 ; imom < 7 ; imom++)
	for(int it = 0 ; it < latInfo.L[3] ; it++)
	  for(int ioper = 0 ; ioper < 16 ; ioper++)
	    fprintf(out_file_cD,"%02d %02d %02d %+d %+d %+d %+16.15e %+16.15e\n",it,ioper,nu,momList[0][imom],momList[1][imom],momList[2][imom],loops_cD[it*4*16*7+nu*16*7+ioper*7+imom].real() , loops_cD[it*4*16*7+nu*16*7+ioper*7+imom].imag() );
    
  }

  /*
  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fprintf(out_file,"%+e %+e\n",vec_out->M[iv][mu][ic].real() , vec_out->M[iv][mu][ic].imag() );
  */

  free(loops_st);
  free(loops_ax);
  free(loops_gD);
  free(loops_cD);

  delete vec_in;
  delete vec_tmp;
  delete gauge;
  delete vec_out;
  MPI_Finalize();

}
