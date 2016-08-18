#include <qkxTM.h>
#include <lattice_util.h>
#include <solvers.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;
  
  //double start_time;
  //double end_time;

  latInfo.L[0] = 24;
  latInfo.L[1] = 24;
  latInfo.L[2] = 24;
  latInfo.L[3] = 48;

  latInfo.P[0] = 1;
  latInfo.P[1] = 2;
  latInfo.P[2] = 4;
  latInfo.P[3] = 8;

  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = 0; // change later
  latInfo.tol = 1e-6;
  latInfo.maxIter = 50000;
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
  Gauge *gauge_APE = new Gauge(&latInfo);
  Vector *vec_out = new Vector(&latInfo);
  Vector *vec_tmp = new Vector(&latInfo);
  Vector *vec_in = new Vector(&latInfo);
  Propagator *uprop = new Propagator(&latInfo);
  Propagator *dprop = new Propagator(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();

  qkxTM_MPI_getGaugeLime(file_conf,*gauge);
  gauge_tmp->copy(*gauge);
  double plaq = gauge->calculatePlaquette();

  if(rank == 0)printf("plaquette is %8.7lf \n",plaq);

  gauge_APE->APE_smearing(*gauge_tmp,&latInfo); 
  plaq = gauge_APE->calculatePlaquette();
  if(rank == 0)printf("plaquette smeared is %8.7lf \n",plaq);

  gauge->applyBoundaryCondition(&latInfo);
  gauge->constGauge = true;


  CG *cg_P = new CG(&latInfo,*gauge);
  CG &cg = *cg_P;
  


  for(int mu = 0 ; mu < 4 ; mu++)
    for(int ic = 0 ; ic < 3 ; ic++){
      vec_tmp->zero();
      vec_out->zero();
      latInfo.twistSign = +1;
      if(rank == 0)  vec_tmp->M[0][mu][ic].real() = 1.;
      vec_out->Gaussian_smearing(*vec_tmp,*gauge_APE,&latInfo);

      vec_in->tmDag(*vec_out,*gauge,&latInfo);
      vec_out->zero();

      cg(*vec_out,*vec_in,&latInfo);
      if(rank==0)printf("mu=%d,c=%d up\n",mu,ic);
      vec_out->rescale(2*latInfo.kappa);
      uprop->absorbVector(*vec_out,mu,ic);
      //      vec_in->Gaussian_smearing(*vec_out,*gauge_APE,&latInfo);
      //uprop->absorbVector(*vec_in,mu,ic);
    }


  
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int ic = 0 ; ic < 3 ; ic++){
      vec_tmp->zero();
      vec_out->zero();
      latInfo.twistSign = -1;
      if(rank == 0)  vec_tmp->M[0][mu][ic].real() = 1.;
      vec_out->Gaussian_smearing(*vec_tmp,*gauge_APE,&latInfo);

      vec_in->tmDag(*vec_out,*gauge,&latInfo);
      vec_out->zero();
      
      cg(*vec_out,*vec_in,&latInfo);
      if(rank==0)printf("mu=%d,c=%d down\n",mu,ic);
      vec_out->rescale(2*latInfo.kappa);
      dprop->absorbVector(*vec_out,mu,ic);
      //vec_in->Gaussian_smearing(*vec_out,*gauge_APE,&latInfo);
      // dprop->absorbVector(*vec_in,mu,ic);
    }
  
  /*  
  char file_up[257] = "/home/khadjiyiannakou/qkxTM_MPI/bin/test/propu";
  char file_down[257] = "/home/khadjiyiannakou/qkxTM_MPI/bin/test/propd";

  sprintf(file_up,"%s.%d.dat",file_up,rank);
  sprintf(file_down,"%s.%d.dat",file_down,rank);

  FILE *p_u,*p_d;
  p_u = fopen(file_up,"r");
  p_d = fopen(file_down,"r");
  int i_int;
  for(int iv = 0 ;iv < 24*24*24*6 ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3 ; c2++){
	    fscanf(p_u,"%d %d %d %d %d %lf %lf",&i_int,&i_int,&i_int,&i_int,&i_int,&(uprop->M[iv][mu][nu][c1][c2].real()) , &(uprop->M[iv][mu][nu][c1][c2].imag()) );
	    fscanf(p_d,"%d %d %d %d %d %lf %lf",&i_int,&i_int,&i_int,&i_int,&i_int,&(dprop->M[iv][mu][nu][c1][c2].real()) , &(dprop->M[iv][mu][nu][c1][c2].imag()) );
	  }
  */

  uprop->rotateToPhysicalBasePlus();
  dprop->rotateToPhysicalBaseMinus();


  for(int nu = 0 ; nu < 4 ; nu++)
    for(int c2 = 0 ; c2 < 3 ; c2++){
      vec_tmp->copyPropagator(*uprop,nu,c2);
      vec_out->Gaussian_smearing(*vec_tmp ,*gauge_APE, &latInfo );
      uprop->absorbVector(*vec_out,nu,c2);

      vec_tmp->copyPropagator(*dprop,nu,c2);
      vec_out->Gaussian_smearing(*vec_tmp ,*gauge_APE, &latInfo );
      dprop->absorbVector(*vec_out,nu,c2);
                                                                                                                                                                             
    }


    
    
 
    /*
  if(rank == 0){
    printf("NOW PRINT THE ELEMENTS OF PROP FOR FIRST POINT ZERO COLORS\n");
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int nu = 0 ; nu < NSPINS ; nu++)
	printf("%+e %+e %+e %+e\n",uprop->M[0][mu][nu][0][0].real() , uprop->M[0][mu][nu][0][0].imag(), dprop->M[0][mu][nu][0][0].real() , dprop->M[0][mu][nu][0][0].imag());
  }
  */
  /*  
  if(rank == 0){
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3 ; c2++)
	    printf("%+e %+e\n",dprop->M[0][mu][nu][c1][c2].real(),dprop->M[0][mu][nu][c1][c2].imag() );
    
  }
  */
      
  vec_tmp->zero();
  vec_out->zero();

  vec_tmp->seqSourceFixSinkPart2(*uprop,0,0,0,1); 
  MPI_Barrier(MPI_COMM_WORLD);
  vec_out->Gaussian_smearing(*vec_tmp,*gauge_APE,&latInfo);


  if(rank == 0){
    printf("NOW PRINT THE ELEMENTS OF PROP FOR FIRST POINT ZERO COLORS\n");
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++)
	printf("mu=%d ic=%d %+e %+e\n",mu,ic,vec_out->M[0][mu][ic].real() , vec_out->M[0][mu][ic].imag());
  }
  

  // create seq source

  /*
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int ic = 0 ; ic < 3 ; ic++){
      vec_tmp->zero();
      vec_out->zero();
      latInfo.twistSign = -1;
      if(rank == 0)  vec_tmp->M[0][mu][ic].real() = 1.;
      vec_out->Gaussian_smearing(*vec_tmp,*gauge_APE,&latInfo);

      vec_in->tmDag(*vec_out,*gauge,&latInfo);
      vec_out->zero();

      cg(*vec_out,*vec_in,&latInfo);
      vec_in->Gaussian_smearing(*vec_out,*gauge_APE,&latInfo);
      dprop->absorbVector(*vec_in,mu,ic);
    }

  */

  

  //prop->absorbVector(*vec_out,mu,ic);
  /*
  char filename_out[257] = "/home/khadjiyiannakou/qkxTM_MPI/bin/test/prop";
  sprintf(filename_out,"%s.%d.dat",filename_out,rank);
  FILE *ptr_out;
  ptr_out = fopen(filename_out,"w");

  for(int iv = 0 ; iv < prop->lV4d ; iv++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int c1 = 0 ; c1 < 3 ; c1++)
	  for(int c2 = 0 ; c2 < 3 ; c2++)
	    fprintf(ptr_out,"%d %d %d %d %d %+e %+e\n",iv,mu,nu,c1,c2,prop->M[iv][mu][nu][c1][c2].real() , prop->M[iv][mu][nu][c1][c2].imag());
  */
      
  /*
  if(rank == 0){
    printf("NOW PRINT THE ELEMENTS OF PROP FOR FIRST POINT ZERO COLORS\n");
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int nu = 0 ; nu < NSPINS ; nu++)
	printf("mu=%d nu=%d %+e %+e\n",mu,nu,prop->M[0][mu][nu][0][0].real() , prop->M[0][mu][nu][0][0].imag());
  }
  */
  delete vec_in;
  delete vec_tmp;
  delete gauge;
  delete vec_out;
  delete uprop;
  delete dprop;

  MPI_Finalize();

}
