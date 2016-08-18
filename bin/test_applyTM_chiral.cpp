#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

int main(int argc,char* argv[]){

  LatticeInfo latInfo;

  double start_time;
  double end_time;

  latInfo.L[0] = 8;
  latInfo.L[1] = 8;
  latInfo.L[2] = 8;
  latInfo.L[3] = 8;

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


  char file_conf[] = "/users/krikitos/scratch/8to4/conf/conf_8c8.0000";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&(latInfo.Nprocs));         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&(latInfo.rank));             // each process gets its ID
  int rank = latInfo.rank;
  //  MPI_Get_processor_name(processor_name,&namelen); //                



  Gauge *gauge = new Gauge(&latInfo);
  Vector *vec = new Vector(&latInfo);
  Vector *tmp = new Vector(&latInfo);
  Vector *vecIn = new Vector(&latInfo);

  printfQKXTM("global lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->L[0],gauge->L[1],gauge->L[2] , gauge->L[3]);
  printfQKXTM("local lattice (x,y,z,t) =  %d x %d x %d x %d\n",gauge->lL[0],gauge->lL[1],gauge->lL[2] , gauge->lL[3]);

  gauge->zero();
  qkxTM_MPI_getGaugeLime(file_conf,*gauge);

  //  for(int c1 = 0 ; c1 < 3 ; c1++)                                                                                                                                                            
  //  for(int c2 = 0 ; c2 < 3 ; c2++){                                                                                                                                                         
  //    printf("%+9.8f %+9.8f\n",gauge->M[LEXIC(1,0,0,1,gauge->lL)][3][c1][c2].real(),gauge->M[LEXIC(1,0,0,1,gauge->lL)][3][c1][c2].imag());
  //  }

  /*
  for(int t=0; t < gauge->lL[3];t++)
    for(int z=0; z < gauge->lL[2];z++)
      for(int y=0; y < gauge->lL[1];y++)
	for(int x=0; x < gauge->lL[0];x++)
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int c1 = 0 ; c1 < 3 ; c1++)
	      for(int c2 = 0 ; c2 < 3 ; c2++){
		gauge->M[LEXIC(t,z,y,x,gauge->lL)][mu][c1][c2].real() = c1==c2 ? 1. : 0.;
	      }
  */
  double plaq = gauge->calculatePlaquette();
  printfQKXTM("Plaquette is %e\n",plaq);
  gauge->applyBoundaryCondition(&latInfo);

  /*
  char filename[]="/users/krikitos/scratch/8to4/vec";
  char filename_out[257];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int c = 0 ; c < 3 ; c++){
      int i = mu*3+c;
      vecIn->zero();
      vec->zero();
      vecIn->M[LEXIC(2,0,0,1,vecIn->lL)][mu][c].real() = 1.;
      vec->tm_xxDagTm_xx_chiral(*vecIn,*gauge,&latInfo,OO);
      //  tmp->applyDslash_xx_chiral(*vecIn,*gauge,EO);
      // vec->applyTwistInv_chiral(*tmp,&latInfo);

      FILE *ptr_out;
      sprintf(filename_out,"%s_%02d.dat",filename,i);
      ptr_out = fopen(filename_out,"w");

      for(int t=0; t < vecIn->lL[3];t++)
	for(int z=0; z < vecIn->lL[2];z++)
	  for(int y=0; y < vecIn->lL[1];y++)
	    for(int x=0; x < vecIn->lL[0];x++)
	      for(int nu = 0 ; nu < 4 ; nu++)
		for(int c2 = 0 ; c2 < 3 ; c2++)
		  fprintf(ptr_out,"%d %d %d %d %+e %+e\n",t,z,y,x,vec->M[LEXIC(t,z,y,x,vecIn->lL)][nu][c2].real(),vec->M[LEXIC(t,z,y,x,vecIn->lL)][nu][c2].imag());
      fclose(ptr_out);
    }
  */
  
  /*
  vecIn->zero();
  for(int iv = 0 ; iv < vecIn->lV4d ; iv++)
    vecIn->M[iv][0][0].real() = 1.;

  vecIn->rotateFromChiralToUKQCD();
  vec->tm_xxDagTm_xx_chiral(*vecIn,*gauge,&latInfo,OO);
  vec->rotateFromChiralToUKQCD();
  */

  char file_eigen[] = "/users/krikitos/scratch/8to4/eigenVec_new/ev.0000.00000";
  qkxTM_MPI_readEigenVectors(file_eigen,*vecIn);
  //  vecIn->sendEvenToOdd();
  vec->tm_xxDagTm_xx_chiral(*vecIn,*gauge,&latInfo,OO);
  
  for(int t=0; t < vecIn->lL[3];t++)                                                                                                                                                               
    for(int z=0; z < vecIn->lL[2];z++)                                                                                                                                                             
      for(int y=0; y < vecIn->lL[1];y++)                                                                                                                                                           
	for(int x=0; x < vecIn->lL[0];x++)       
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int c1 = 0 ; c1 < 3 ; c1++)
	      printf("%d %d %d %d %+e %+e \t %+e %+e\n",t,z,y,x,vecIn->M[LEXIC(t,z,y,x,vecIn->lL)][mu][c1].real(),vecIn->M[LEXIC(t,z,y,x,vecIn->lL)][mu][c1].imag(), vec->M[LEXIC(t,z,y,x,vecIn->lL)][mu][c1].real(),vec->M[LEXIC(t,z,y,x,vecIn->lL)][mu][c1].imag());
  
  delete gauge;
  delete vecIn;
  delete vec;
  delete tmp;
  MPI_Finalize();

}
