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

  latInfo.P[0] = 2;
  latInfo.P[1] = 2;
  latInfo.P[2] = 2;
  latInfo.P[3] = 1;

  latInfo.kappa = 0.160856;
  latInfo.mu = 0.004000;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-6;
  latInfo.maxIter = 5000;
  latInfo.reliableDelta = 0.01;
  latInfo.boundaryT = -1;                   // remember write program that applies boundaries

  latInfo.NsmearAPE = 20;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 50;
  latInfo.alphaGauss = 4.0;
  
  

  char file_prop_up[] = "/home/khadjiyiannakou/test_write_props/prop_SS.up";
  char file_prop_down[] = "/home/khadjiyiannakou/test_write_props/prop_SS.down";
  char file_gauge[] = "/home/khadjiyiannakou/confs/L16T32/conf.1505";
  

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

  Propagator *prop_up = new Propagator(&latInfo);
  Propagator *prop_down = new Propagator(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",vec->L[0],vec->L[1],vec->L[2] , vec->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",vec->lL[0],vec->lL[1],vec->lL[2] , vec->lL[3]);

  gauge->zero();
  qkxTM_MPI_getGaugeLime(file_gauge,*gauge);
  gauge_tmp->copy(*gauge);
  double plaq = gauge->calculatePlaquette();
  if(rank == 0)printf("plaquette is %8.7lf \n",plaq);


  start_time = MPI_Wtime();
  gauge_ape->APE_smearing(*gauge_tmp,&latInfo);
  delete gauge_tmp; // dont need it anymore
  end_time = MPI_Wtime();
  if(rank == 0)printf("finish smearing in %e secs\n",end_time-start_time);
  plaq = gauge_ape->calculatePlaquette();
  
  if(rank == 0)printf("plaquette smeared is %8.7lf \n",plaq);


   start_time = MPI_Wtime();
   qkxTM_MPI_getPropLime(file_prop_up,*prop_up,&latInfo);
   qkxTM_MPI_getPropLime(file_prop_down,*prop_down,&latInfo);
   end_time = MPI_Wtime();
   if(rank == 0)printf("read props in time is %e in sec\n",end_time-start_time);



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
   */
   

   prop_up->rotateToPhysicalBasePlus();
   prop_down->rotateToPhysicalBaseMinus();

   
   char twop_proton[257] = "/home/khadjiyiannakou/qkxTM_MPI/bin/twop_proton.dat";
   char twop_neutron[257] = "/home/khadjiyiannakou/qkxTM_MPI/bin/twop_neutron.dat";
   int x_src[4];
   x_src[0] = 0;
   x_src[1] = 0;
   x_src[2] = 0;
   x_src[3] = 0;
   contractNucleon(*prop_up,*prop_down,twop_proton,16,x_src);
   contractNucleon(*prop_down,*prop_up,twop_neutron,16,x_src);
   

   delete prop_up;
   delete prop_down;
   delete gauge_ape;
   delete gauge;
   delete vec;
   delete vec_tmp;
  MPI_Finalize();

}
