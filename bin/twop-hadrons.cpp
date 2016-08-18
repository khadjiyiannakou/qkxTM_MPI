#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  if(argc != 20)errorQKXTM("wrong number on inputs\n");

  LatticeInfo latInfo;

  double start_time;
  double end_time;


  latInfo.L[0] = atoi(argv[1]);
  latInfo.L[1] = atoi(argv[2]);
  latInfo.L[2] = atoi(argv[3]);
  latInfo.L[3] = atoi(argv[4]);

  latInfo.P[0] = atoi(argv[5]);
  latInfo.P[1] = atoi(argv[6]);
  latInfo.P[2] = atoi(argv[7]);
  latInfo.P[3] = atoi(argv[8]);


  char file_prop_up[257];
  char file_prop_down[257];
  char file_gauge[257];
  
  strcpy(file_gauge,argv[9]);
  strcpy(file_prop_up,argv[10]);
  strcpy(file_prop_down,argv[11]);

  char filename_out[257];
  strcpy(filename_out,argv[12]);

  latInfo.x_src[0] = atoi(argv[13]);
  latInfo.x_src[1] = atoi(argv[14]);;
  latInfo.x_src[2] = atoi(argv[15]);;
  latInfo.x_src[3] = atoi(argv[16]);

  latInfo.NsmearGauss = atoi(argv[17]);
  latInfo.alphaGauss = atof(argv[18]);
  
  latInfo.Q_sq = atoi(argv[19]);  


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                


  Gauge *gauge_ape = new Gauge(&latInfo);

  Vector *vec = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);

  Propagator *prop_up = new Propagator(&latInfo);
  Propagator *prop_down = new Propagator(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",vec->L[0],vec->L[1],vec->L[2] , vec->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",vec->lL[0],vec->lL[1],vec->lL[2] , vec->lL[3]);


  gauge_ape->zero();
  qkxTM_MPI_getGaugeLime(file_gauge,*gauge_ape);
  double plaq = gauge_ape->calculatePlaquette();
  if(rank == 0)printf("plaquette smeared is %8.7lf \n",plaq);

   start_time = MPI_Wtime();
   qkxTM_MPI_getPropLime(file_prop_up,*prop_up,&latInfo);
   qkxTM_MPI_getPropLime(file_prop_down,*prop_down,&latInfo);
   end_time = MPI_Wtime();
   if(rank == 0)printf("read props in time is %f in sec\n",end_time-start_time);


  prop_up->rotateToPhysicalBasePlus();
  prop_down->rotateToPhysicalBaseMinus();

  if(rank == 0)printf("Perform contractions using SL propagator\n");
  // contractions for SL propagators

  char file_mesons_SL[257];
  char file_baryons_SL[257];
  char file_baryons_dec_SL[257];

  char file_mesons_SS[257];
  char file_baryons_SS[257];
  char file_baryons_dec_SS[257];

  sprintf(file_mesons_SL,"%s_%s.dat",filename_out,"mesons_SL");
  sprintf(file_baryons_SL,"%s_%s.dat",filename_out,"baryons_SL");
  sprintf(file_baryons_dec_SL,"%s_%s.dat",filename_out,"baryons_dec_SL");

  sprintf(file_mesons_SS,"%s_%s.dat",filename_out,"mesons_SS");
  sprintf(file_baryons_SS,"%s_%s.dat",filename_out,"baryons_SS");
  sprintf(file_baryons_dec_SS,"%s_%s.dat",filename_out,"baryons_dec_SS");


  // mesons
  start_time = MPI_Wtime();
  //  contract_Mesons_omp(prop_up,prop_down,file_mesons_SL,&latInfo);
  contractMesons(prop_up,prop_down,file_mesons_SL,&latInfo);
  contractBaryons(prop_up,prop_down,file_baryons_SL,&latInfo);
  contractBaryonsDec(prop_up,prop_down,file_baryons_dec_SL,&latInfo);
  end_time = MPI_Wtime();
  if(rank == 0)printf("perform meson contractions in %f in sec\n",end_time-start_time);

  if(rank == 0)printf("Apply Gaussian smearing\n");

  /*

  for(int nu = 0 ; nu < 4 ; nu++)
    for(int c2 = 0 ; c2 < 3 ; c2++){
      vec_tmp->copyPropagator(*prop_up,nu,c2);
      vec->Gaussian_smearing(*vec_tmp,*gauge_ape,&latInfo);
      prop_up->copyVector(*vec,nu,c2);
      
      vec_tmp->copyPropagator(*prop_down,nu,c2);
      vec->Gaussian_smearing(*vec_tmp,*gauge_ape,&latInfo);
      prop_down->copyVector(*vec,nu,c2);   
    }
  
  if(rank == 0)printf("Perform contractions using SS propagator\n");


  start_time = MPI_Wtime();
  contractMesons(prop_up,prop_down,file_mesons_SS,&latInfo);
  end_time = MPI_Wtime();
  if(rank == 0)printf("perform meson contractions in %f in sec\n",end_time-start_time);

  start_time = MPI_Wtime();
  contractBaryons(prop_up,prop_down,file_baryons_SS,&latInfo);
  end_time = MPI_Wtime();
  if(rank == 0)printf("perform baryon contractions in %f in sec\n",end_time-start_time);

  start_time = MPI_Wtime();
  contractBaryonsDec(prop_up,prop_down,file_baryons_dec_SS,&latInfo);
  end_time = MPI_Wtime();
  if(rank == 0)printf("perform baryon dec contractions in %f in sec\n",end_time-start_time);
  */
   delete prop_up;
   delete prop_down;
   delete gauge_ape;
   // delete gauge;
   delete vec;
   delete vec_tmp;
  MPI_Finalize();

}
