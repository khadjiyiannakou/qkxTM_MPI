#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  //  double start_time;
  //double end_time;

  latInfo.L[0] = 24;
  latInfo.L[1] = 24;
  latInfo.L[2] = 24;
  latInfo.L[3] = 48;

  latInfo.P[0] = 1;
  latInfo.P[1] = 1;
  latInfo.P[2] = 1;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
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

  char file_conf[257] = "/home/khadjiyiannakou/qudaQKXTM_package/confs/L24T48/conf.1000";

  char file_in[257] = "/home/khadjiyiannakou/qkxTM_MPI/test/Spinor.In";
  char file_out[257] = "/home/khadjiyiannakou/qkxTM_MPI/test/outSpinor";

  //  sprintf(file_in,"%s.%d.dat",file_in,rank);
  sprintf(file_out,"%s.%d.dat",file_out,rank);

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
  
  FILE *in_file;
  FILE *out_file;
  
  in_file = fopen(file_in,"r");
  printf("%d %s\n",rank,file_in);


  int ierr;
  if(in_file == NULL)MPI_Abort(MPI_COMM_WORLD,ierr);
  out_file = fopen(file_out,"w");
  int local_volume = gauge->lL[0] * gauge->lL[1] * gauge->lL[2] * gauge->lL[3];

  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(in_file,"%lf %lf",&(vec_tmp->M[iv][mu][ic].real()) , &(vec_tmp->M[iv][mu][ic].imag()) );

  MPI_Barrier(MPI_COMM_WORLD);

  gauge->applyBoundaryCondition(&latInfo);
  gauge->constGauge = true;

  vec_in->tmDag(*vec_tmp,*gauge,&latInfo);
  vec_out->zero();

  // apply boundary conditions on gauge



  CG *cg_P = new CG(&latInfo,*gauge);
  CG &cg = *cg_P;
  cg(*vec_out,*vec_in,&latInfo);

  for(int iv = 0 ; iv < local_volume ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int ic = 0 ; ic < 3 ; ic++)
	fscanf(out_file,"%lf %lf",&(vec_out->M[iv][mu][ic].real()) , &(vec_out->M[iv][mu][ic].imag()) );


  delete vec_in;
  delete vec_tmp;
  delete gauge;
  delete vec_out;
  MPI_Finalize();

}
