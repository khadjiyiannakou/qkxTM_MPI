#include <qkxTM.h>
#include <lattice_util.h>
#include <omp.h>

#define PI 3.14159265359
#define DEBUG_MODE


extern MPI_Group fullGroup , spaceGroup , timeGroup;
extern MPI_Comm spaceComm , timeComm;
extern int localRank;
extern int localSize;
extern int timeRank;
extern int timeSize;
extern bool initializeGroupsFlag;


using namespace qkxTM;

static const int eps[6][3]=
  {
    {0,1,2},
    {2,0,1},
    {1,2,0},
    {2,1,0},
    {0,2,1},
    {1,0,2}
  };

static const int sgn_eps[6]=
  {
    +1,+1,+1,-1,-1,-1
  };

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// CLASS MESONS /////////////////
static int createMomenta(int momElem [][3], int Q_sq){
  int counter = 0;



  for(int iQ = 0 ; iQ <= Q_sq ; iQ++){
    for(int nx = iQ ; nx >= -iQ ; nx--)
      for(int ny = iQ ; ny >= -iQ ; ny--)
        for(int nz = iQ ; nz >= -iQ ; nz--){
          if( nx*nx + ny*ny + nz*nz == iQ ){
            momElem[counter][0] = nx;
            momElem[counter][1] = ny;
            momElem[counter][2] = nz;
            counter++;
          }
        }
  }
  if(counter > MOM_MAX)errorQKXTM("Error exceeded max number of momenta\n");
  /*
#ifdef DEBUG_MODE
  printfQKXTM("With Q^2 = %d I create %d momenta combinations\n",Q_sq,Nmom);
  for(int imom = 0 ; imom < Nmom ; imom++)
    printfQKXTM("%d \t %+d %+d %+d\n",imom,momElem[imom][0],momElem[imom][1],momElem[imom][2]);
#endif
  */
  return counter;
}

Meson::Meson(LatticeInfo *latInfo,char filename_twop[]) : LatticeGeometry(latInfo) {

  if(P[3] != 1){
    errorQKXTM("Error number of process in time direction must be 1 for this class I can change it later if necessary\n");
  }

  size = lV4d;
  Q_sq = latInfo->Q_sq;
  for(int i = 0 ; i < 4 ; i ++) x_src[i]=latInfo->x_src[i];
  M = (Complex*) safe_malloc(size*2* sizeof(double));
  memset(M,0,size*2*sizeof(double));
  Nmom = createMomenta(momElem,Q_sq);
 
  F = (Complex*) safe_malloc(L[3]*Nmom*2*sizeof(double));

  if(rank == 0){
    ptr_file = fopen(filename_twop,"w");
    if(ptr_file == NULL){
      errorQKXTM("Error cannot open file for writting meson correlators\n");
    }
  }
  else{
    ptr_file=NULL;
  }
  
}

Meson::Meson(LatticeInfo *latInfo): LatticeGeometry(latInfo){
}

Meson::~Meson(){

  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }
  free(M); M = NULL;
  free(F); F = NULL;

  if(ptr_file != NULL){
    fclose(ptr_file);
    ptr_file = NULL;
  }

}


void Meson::zero_M(){
  memset(M,0,size*2*sizeof(double));
}

void Meson::zero_F(){
  memset(F,0,L[3]*Nmom*2*sizeof(double));
}

void Meson::fourier(FOURIER_DIRECTION Direction){

  int sign = ( Direction == FORWARD ) ? -1 : +1;
  Complex *corr = (Complex*)safe_malloc(L[3]*Nmom*2*sizeof(double));
  memset(corr,0,L[3]*Nmom*2*sizeof(double));
  double tmp;

  for(int it = 0 ; it < L[3] ; it++){
    for(int imom = 0 ; imom < Nmom ; imom++){
      for(int lz = 0 ; lz < lL[2] ; lz++)
	for(int ly = 0 ; ly < lL[1] ; ly++)
	  for(int lx = 0 ; lx < lL[0] ; lx++){
	    int v = LEXIC(it,lz,ly,lx,lL);
	    int x = lx + procPos[0]*lL[0] - x_src[0];
	    int y = ly + procPos[1]*lL[1] - x_src[1];
	    int z = lz + procPos[2]*lL[2] - x_src[2];
	    tmp = (((double) momElem[imom][0]*x)/L[0] + ((double) momElem[imom][1]*y)/L[1] + ((double) momElem[imom][2]*z)/L[2])*2*PI;
	    Complex C2 (cos(tmp),sign*sin(tmp));
	    corr[it*Nmom+imom] = corr[it*Nmom+imom] + M[v] * C2;
	  }

    } // momenta
  } // time
  MPI_Reduce(corr, F, L[3]*Nmom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  free(corr);
}

void Meson::dumpData(MESONS_OPERATOR OPER){

  if( rank == 0 ){
    for(int it = 0 ; it < L[3] ; it++)
      for(int imom = 0 ; imom < Nmom ; imom++)
	fprintf(ptr_file,"%d \t %d \t %+d %+d %+d \t %+e %+e\n",OPER,it,momElem[imom][0],momElem[imom][1],momElem[imom][2],F[it*Nmom+imom].real(),F[it*Nmom+imom].imag());
  }

}

int Meson::tabulateIndices(Complex val[],int ind[][4],MESONS_OPERATOR OPER){
  int counter = 0;
  Complex g[4][4];
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      g[mu][nu].real() = 0.; g[mu][nu].imag() = 0.;
    }

    switch(OPER){
    case ONE :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++){
	  if(mu == nu){
	    g[mu][nu].real() = 1.; g[mu][nu].imag() = 0.;
	  }
	  else{
	    g[mu][nu].real() = 0.; g[mu][nu].imag() = 0.;
	  }
	}
      break;
    case G5 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  g[mu][nu] = gammaM[4][mu][nu];
      break;
    case G1 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  g[mu][nu] = gammaM[0][mu][nu];
      break;
    case G2 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  g[mu][nu] = gammaM[1][mu][nu];
      break;
    case G3 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  g[mu][nu] = gammaM[2][mu][nu];
      break;
    case G4 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  g[mu][nu] = gammaM[3][mu][nu];
      break;
    case G5G1 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  for(int ku = 0 ; ku < 4 ; ku++)
	    g[mu][nu] = g[mu][nu] + gammaM[4][mu][ku] * gammaM[0][ku][nu];
      break;
    case G5G2 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  for(int ku = 0 ; ku < 4 ; ku++)
	    g[mu][nu] = g[mu][nu] + gammaM[4][mu][ku] * gammaM[1][ku][nu];
      break;
    case G5G3 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  for(int ku = 0 ; ku < 4 ; ku++)
	    g[mu][nu] = g[mu][nu] + gammaM[4][mu][ku] * gammaM[2][ku][nu];
      break;
    case G5G4 :
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int nu = 0 ; nu < 4 ; nu++)
	  for(int ku = 0 ; ku < 4 ; ku++)
	    g[mu][nu] = g[mu][nu] + gammaM[4][mu][ku] * gammaM[3][ku][nu];
      break;
    default:
      errorQKXTM("Something gone wrong in the switch\n");
    }
    
    Complex C_temp;
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int gamma = 0 ; gamma < 4 ; gamma++)
	for(int delta = 0 ; delta < 4 ; delta++)
	  for(int alpha = 0 ; alpha < 4 ; alpha++){
	    C_temp = g[beta][gamma] * g[delta][alpha];
	    if( norm(C_temp) > 1e-3 ){
	      val[counter] = C_temp;
	      ind[counter][0] = beta;
	      ind[counter][1] = gamma;
	      ind[counter][2] = delta;
	      ind[counter][3] = alpha;
	      counter++;
	    }
	  }
  
  return counter;
}

void Meson::contract(Propagator &prop, MESONS_OPERATOR OPER){

  Complex pp[4][4][3][3];
  Complex g5_pp_dag_g5[4][4][3][3];
  Complex val[16*16];
  int ind[16*16][4];

  int count = tabulateIndices(val,ind,OPER);
  zero_M();

  /*
#pragma omp parallel
  {
    for(int it = 0 ; it < L[3] ; it++)
      for(int iz = 0 ; iz < lL[2] ; iz++)
#pragma omp for collapse(2)
	for(int iy = 0 ; iy < lL[1] ; iy++)
	  for(int ix = 0 ; ix < lL[0] ; ix++){
	    int iv = LEXIC(it,iz,iy,ix,lL);
  */

#pragma omp parallel for
  for(int iv = 0 ; iv < lV4d ; iv++){

      memcpy(pp,&(prop.M[iv][0][0][0][0].real()),4*4*3*3*2*sizeof(double));

      
      for(int c1 = 0 ; c1 < 3 ; c1++)
	for(int c2 = 0 ; c2 < 3 ; c2++)
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int nu = 0 ; nu < 4 ; nu++){
	      g5_pp_dag_g5[mu][nu][c1][c2].real() = 0; g5_pp_dag_g5[mu][nu][c1][c2].imag() = 0;
	      for(int ku = 0 ; ku < 4 ; ku++)
		for(int lu = 0 ; lu < 4 ; lu++)
		  g5_pp_dag_g5[mu][nu][c1][c2] = g5_pp_dag_g5[mu][nu][c1][c2] + gammaM[4][mu][ku] * conj(pp[lu][ku][c2][c1]) * gammaM[4][lu][nu];
	    }
      
      for(int a = 0 ; a < 3 ; a++)
	for(int b = 0 ; b < 3 ; b++)
	  for(int is = 0 ; is < count ; is++){
	    Complex gammaValue = val[is];
	    int beta = ind[is][0];
	    int gamma = ind[is][1];
	    int delta = ind[is][2];
	    int alpha = ind[is][3];
	    M[iv] = M[iv] + gammaValue * pp[alpha][beta][a][b] * g5_pp_dag_g5[gamma][delta][b][a] ;
	  }
  } // close volume


}

Baryon::Baryon(LatticeInfo *latInfo,char filename_twop[]) : LatticeGeometry(latInfo) {

  if(P[3] != 1){
    errorQKXTM("Error number of process in time direction must be 1 for this class I can change it later if necessary\n");
  }

  size = lV4d;
  Q_sq = latInfo->Q_sq;
  for(int i = 0 ; i < 4 ; i ++) x_src[i]=latInfo->x_src[i];

  M = (Complex(*)[4][4]) safe_malloc(size*4*4*2*sizeof(double));
  memset(M,0,size*4*4*2*sizeof(double));
  Nmom = createMomenta(momElem,Q_sq);
  F = (Complex(*)[4][4]) safe_malloc(L[3]*Nmom*4*4*2*sizeof(double));

  if(rank == 0){
    ptr_file = fopen(filename_twop,"w");
    if(ptr_file == NULL){
      errorQKXTM("Error cannot open file for writting Baryon correlators\n");
    }
  }
  else{
    ptr_file=NULL;
  }
  
}


Baryon::~Baryon(){
  
  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }
  
  free(M); M = NULL;
  free(F); F = NULL;

  if(ptr_file != NULL){
    fclose(ptr_file);
    ptr_file = NULL;
  }
}

void Baryon::zero_M(){
  memset(M,0,size*4*4*2*sizeof(double));
}

void Baryon::zero_F(){
  memset(F,0,L[3]*Nmom*4*4*2*sizeof(double));
}

void Baryon::fourier(FOURIER_DIRECTION Direction){

  int sign = ( Direction == FORWARD ) ? -1 : +1;
  Complex *corr = (Complex*)safe_malloc(L[3]*Nmom*4*4*2*sizeof(double));
  memset(corr,0,L[3]*Nmom*4*4*2*sizeof(double));
  double tmp;

  for(int it = 0 ; it < L[3] ; it++){
    for(int imom = 0 ; imom < Nmom ; imom++){
      for(int gamma = 0 ; gamma < 4 ; gamma++)
	for(int gammap = 0 ; gammap < 4 ; gammap++)
	  for(int lz = 0 ; lz < lL[2] ; lz++)
	    for(int ly = 0 ; ly < lL[1] ; ly++)
	      for(int lx = 0 ; lx < lL[0] ; lx++){
		int v = LEXIC(it,lz,ly,lx,lL);
		int x = lx + procPos[0]*lL[0] - x_src[0];
		int y = ly + procPos[1]*lL[1] - x_src[1];
		int z = lz + procPos[2]*lL[2] - x_src[2];
		tmp = (((double) momElem[imom][0]*x)/L[0] + ((double) momElem[imom][1]*y)/L[1] + ((double) momElem[imom][2]*z)/L[2])*2*PI;
		Complex C2 (cos(tmp),sign*sin(tmp));
		corr[it*Nmom*4*4+imom*4*4+gamma*4+gammap] = corr[it*Nmom*4*4+imom*4*4+gamma*4+gammap] + M[v][gamma][gammap] * C2;
	      }
      
    } // momenta
  } // time
  MPI_Reduce(corr, F, L[3]*Nmom*4*4*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  free(corr);

}

int Baryon::tabulateIndices_NTN(Complex val[4*4*4*4],int ind[8][4*4*4*4]){
  int counter = 0;
  Complex C_temp;

  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int betap = 0 ; betap < 4 ; betap++)
	for(int alphap = 0 ; alphap < 4 ; alphap++){
	  C_temp = CgammaM[4][alpha][beta] * CgammaM_bar[4][betap][alphap];
	  if( norm(C_temp) > 1e-3 ){
	    val[counter] = C_temp;
	    ind[0][counter] = alpha;
	    ind[1][counter] = beta;
	    ind[2][counter] = betap;
	    ind[3][counter] = alphap;
	    counter++;
	  }
	}

  return counter;
}

int Baryon::tabulateIndices_NTR(Complex val[4*4*4*4],int ind[8][4*4*4*4]){
  int counter = 0;
  Complex C_temp;

  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int betap = 0 ; betap < 4 ; betap++)
	for(int alphap = 0 ; alphap < 4 ; alphap++)
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    for(int deltap = 0 ; deltap < 4 ; deltap++){
	      C_temp = CgammaM[4][alpha][beta] * Ccg_bar[betap][alphap] * gammaM[4][gammap][deltap];
	      if( norm(C_temp) > 1e-3 ){
		val[counter] = C_temp;
		ind[0][counter] = alpha;
		ind[1][counter] = beta;
		ind[2][counter] = betap;
		ind[3][counter] = alphap;
		ind[4][counter] = gammap;
		ind[5][counter] = deltap;
		counter++;
	      }
	    }

  
  return counter;
}

int Baryon::tabulateIndices_RTN(Complex val[4*4*4*4],int ind[8][4*4*4*4]){
  int counter = 0;
  Complex C_temp;

  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int betap = 0 ; betap < 4 ; betap++)
	for(int alphap = 0 ; alphap < 4 ; alphap++)
	  for(int gamma = 0 ; gamma < 4 ; gamma++)
	    for(int delta = 0 ; delta < 4 ; delta++){
	      C_temp = Ccg[alpha][beta] * CgammaM_bar[4][betap][alphap] * gammaM[4][gamma][delta];
	      if( norm(C_temp) > 1e-3 ){
		val[counter] = C_temp;
		ind[0][counter] = alpha;
		ind[1][counter] = beta;
		ind[2][counter] = betap;
		ind[3][counter] = alphap;
		ind[4][counter] = gamma;
		ind[5][counter] = delta;
		counter++;
	      }
	    }
  
  return counter;
}

int Baryon::tabulateIndices_RTR(Complex val[4*4*4*4],int ind[8][4*4*4*4]){
  int counter = 0;
  Complex C_temp;

  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int betap = 0 ; betap < 4 ; betap++)
	for(int alphap = 0 ; alphap < 4 ; alphap++)
	  for(int gamma = 0 ; gamma < 4 ; gamma++)
	    for(int delta = 0 ; delta < 4 ; delta++)
	      for(int gammap = 0 ; gammap < 4 ; gammap++)
		for(int deltap = 0 ; deltap < 4 ; deltap++){
		  C_temp = Ccg[alpha][beta] * Ccg_bar[betap][alphap] * gammaM[4][gamma][delta] * gammaM[4][gammap][deltap];
		  if( norm(C_temp) > 1e-3 ){
		    val[counter] = C_temp;
		    ind[0][counter] = alpha;
		    ind[1][counter] = beta;
		    ind[2][counter] = betap;
		    ind[3][counter] = alphap;
		    ind[4][counter] = gamma;
		    ind[5][counter] = delta;
		    ind[6][counter] = gammap;
		    ind[7][counter] = deltap;
		counter++;
	      }
	    }
  
  return counter;
}



void Baryon::contract(Propagator &prop1,Propagator &prop2,BARYON_OPERATOR OPER){
  Complex val[4*4*4*4];
  int ind[8][4*4*4*4];;

  zero_M();
  Pointers pointer;
  
  int counter;
  switch (OPER){
  case NTN:
    counter = tabulateIndices_NTN(val,ind);
    break;
  case NTR:
    counter = tabulateIndices_NTR(val,ind);
    break;
  case RTN:
    counter = tabulateIndices_RTN(val,ind);
    break;
  case RTR:
    counter = tabulateIndices_RTR(val,ind);
    break;
  default:
    errorQKXTM("Error with the switch");
  }//close switch
  


  int a,b,c,a1,b1,c1;
  switch (OPER){
  case NTN:

    for(int it = 0 ; it < L[3] ; it++)
      for(int gamma = 0 ; gamma < 4 ; gamma++)       
	for(int gammap = 0 ; gammap < 4 ; gammap++)
	  for(int idx = 0 ; idx < counter ; idx++){
	    int alpha = ind[0][idx];
	    int beta = ind[1][idx];
	    int betap = ind[2][idx];
	    int alphap = ind[3][idx];
	    for(int cc1 = 0 ; cc1 < 6 ; cc1++){
	      a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
	      for(int cc2 = 0 ; cc2 < 6 ; cc2++){
		a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
		Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[idx];
		for(int lz = 0 ; lz < lL[2] ; lz++)
		  for(int ly = 0 ; ly < lL[1] ; ly++)
		    for(int lx = 0 ; lx < lL[0] ; lx++){
		      int iv = LEXIC(it,lz,ly,lx,lL);
		      M[iv][gamma][gammap] = M[iv][gamma][gammap] + factor * prop2.M[iv][beta][betap][b][b1] * 
		      	( prop1.M[iv][alpha][alphap][a][a1] * prop1.M[iv][gamma][gammap][c][c1] - prop1.M[iv][alpha][gammap][a][c1] * prop1.M[iv][gamma][alphap][c][a1] );

		    }
	      }
	    }
	  }

    break;
  case NTR:
    for(int it = 0 ; it < L[3] ; it++)
      for(int gamma = 0 ; gamma < 4 ; gamma++)       
	for(int idx = 0 ; idx < counter ; idx++){
	  int alpha = ind[0][idx];
	  int beta = ind[1][idx];
	  int betap = ind[2][idx];
	  int alphap = ind[3][idx];
	  int gammap = ind[4][idx];
	  int deltap = ind[5][idx];
	  for(int cc1 = 0 ; cc1 < 6 ; cc1++){
	    a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
	    for(int cc2 = 0 ; cc2 < 6 ; cc2++){
	      a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
	      Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[idx];
	      for(int lz = 0 ; lz < lL[2] ; lz++)
		for(int ly = 0 ; ly < lL[1] ; ly++)
		  for(int lx = 0 ; lx < lL[0] ; lx++){
		    int iv = LEXIC(it,lz,ly,lx,lL);
		    M[iv][gamma][gammap] = M[iv][gamma][gammap] - factor * prop2.M[iv][beta][betap][b][b1] * 
		      ( prop1.M[iv][alpha][alphap][a][a1] * prop1.M[iv][gamma][deltap][c][c1] - prop1.M[iv][alpha][deltap][a][c1] * prop1.M[iv][gamma][alphap][c][a1] );
		    if(it == 0 && gamma == 0 && iv == 0) printf("%d %d %d %d %d %d %d %d %d %d %d %+e %+e\n",gamma,gammap,idx,alpha,beta,betap,alphap,gammap,deltap,cc1,cc2,M[iv][gamma][gammap].real(),M[iv][gamma][gammap].imag());
		  }
	      }
	    }
	  }

    break;
  case RTN:
    
    for(int it = 0 ; it < L[3] ; it++)
      for(int gammap = 0 ; gammap < 4 ; gammap++)       
	for(int idx = 0 ; idx < counter ; idx++){
	  int alpha = ind[0][idx];
	  int beta = ind[1][idx];
	  int betap = ind[2][idx];
	  int alphap = ind[3][idx];
	  int gamma = ind[4][idx];
	  int delta = ind[5][idx];
	  for(int cc1 = 0 ; cc1 < 6 ; cc1++){
	    a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
	    for(int cc2 = 0 ; cc2 < 6 ; cc2++){
	      a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
	      Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[idx];
	      for(int lz = 0 ; lz < lL[2] ; lz++)
		for(int ly = 0 ; ly < lL[1] ; ly++)
		  for(int lx = 0 ; lx < lL[0] ; lx++){
		    int iv = LEXIC(it,lz,ly,lx,lL);
		    M[iv][gamma][gammap] = M[iv][gamma][gammap] + factor * prop2.M[iv][beta][betap][b][b1] * 
		      ( prop1.M[iv][alpha][alphap][a][a1] * prop1.M[iv][delta][gammap][c][c1] - prop1.M[iv][alpha][gammap][a][c1] * prop1.M[iv][delta][alphap][c][a1] );
		  }
	      }
	    }
	  }
    break;
  case RTR:
    for(int it = 0 ; it < L[3] ; it++)
      for(int idx = 0 ; idx < counter ; idx++){
	int alpha = ind[0][idx];
	int beta = ind[1][idx];
	int betap = ind[2][idx];
	int alphap = ind[3][idx];
	int gamma = ind[4][idx];
	int delta = ind[5][idx];
	int gammap = ind[6][idx];
	int deltap = ind[7][idx];
	for(int cc1 = 0 ; cc1 < 6 ; cc1++){
	  a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
	  for(int cc2 = 0 ; cc2 < 6 ; cc2++){
	    a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
	    Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[idx];
	    for(int lz = 0 ; lz < lL[2] ; lz++)
	      for(int ly = 0 ; ly < lL[1] ; ly++)
		for(int lx = 0 ; lx < lL[0] ; lx++){
		  int iv = LEXIC(it,lz,ly,lx,lL);
		  M[iv][gamma][gammap] = M[iv][gamma][gammap] - factor * prop2.M[iv][beta][betap][b][b1] * 
		    ( prop1.M[iv][alpha][alphap][a][a1] * prop1.M[iv][delta][deltap][c][c1] - prop1.M[iv][alpha][deltap][a][c1] * prop1.M[iv][delta][alphap][c][a1] );
		}
	  }
	}
      }
    
    break;
  default:
    errorQKXTM("Problem with the switch\n");
  }
  

}

void Baryon::dumpData(BARYON_OPERATOR OPER){

  if( rank == 0 ){
    for(int it = 0 ; it < L[3] ; it++)
      for(int imom = 0 ; imom < Nmom ; imom++)
	for(int gamma = 0 ; gamma < 4 ; gamma++)
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    fprintf(ptr_file,"%d \t %d \t %+d %+d %+d \t %d %d \t %+e %+e\n",OPER,it,momElem[imom][0],momElem[imom][1],momElem[imom][2],gamma,gammap,F[it*Nmom+imom][gamma][gammap].real(),F[it*Nmom+imom][gamma][gammap].imag());
  }

}

//////// Baryon_Dec
Baryon_Dec::Baryon_Dec(LatticeInfo *latInfo,char filename_twop[]) : LatticeGeometry(latInfo) {

  if(P[3] != 1){
    errorQKXTM("Error number of process in time direction must be 1 for this class I can change it later if necessary\n");
  }

  size = lV4d;
  Q_sq = latInfo->Q_sq;
  for(int i = 0 ; i < 4 ; i ++) x_src[i]=latInfo->x_src[i];
  M = (Complex(*)[3][4][4]) safe_malloc(size*3*4*4*2*sizeof(double));
  memset(M,0,size*3*4*4*2*sizeof(double));
  Nmom = createMomenta(momElem,Q_sq);
  
  F = (Complex(*)[3][4][4]) safe_malloc(L[3]*Nmom*3*4*4*2*sizeof(double));

  if(rank == 0){
    ptr_file = fopen(filename_twop,"w");
    if(ptr_file == NULL){
      errorQKXTM("Error cannot open file for writting Baryon decuplet correlators\n");
    }
  }
  else{
    ptr_file=NULL;
  }
  
}


Baryon_Dec::~Baryon_Dec(){
  for(int i = 0 ; i < NDIM ; i++)
    if( lL[i] < L[i] ){
      MPI_Type_free( (MPI_Datatype*) &(stypeV[i]));
      MPI_Type_free( (MPI_Datatype*) &(rtypeV[i]));
    }

  free(M);
  free(F);

  if(ptr_file != NULL)  fclose(ptr_file);  
}

void Baryon_Dec::zero_M(){
  memset(M,0,size*3*4*4*2*sizeof(double));
}

void Baryon_Dec::zero_F(){
  memset(F,0,L[3]*Nmom*3*4*4*2*sizeof(double));
}

void Baryon_Dec::fourier(FOURIER_DIRECTION Direction){

  int sign = ( Direction == FORWARD ) ? -1 : +1;
  Complex *corr = (Complex*)safe_malloc(L[3]*Nmom*4*4*2*sizeof(double));
  Complex *F_tmp = (Complex*)safe_malloc(L[3]*Nmom*4*4*2*sizeof(double));
  double tmp;

  for(int i = 0 ; i < 3 ; i++){
    memset(corr,0,L[3]*Nmom*4*4*2*sizeof(double));
    memset(F_tmp,0,L[3]*Nmom*4*4*2*sizeof(double));
    for(int it = 0 ; it < L[3] ; it++){
      for(int imom = 0 ; imom < Nmom ; imom++){
	for(int gamma = 0 ; gamma < 4 ; gamma++)
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    for(int lz = 0 ; lz < lL[2] ; lz++)
	      for(int ly = 0 ; ly < lL[1] ; ly++)
		for(int lx = 0 ; lx < lL[0] ; lx++){
		  int v = LEXIC(it,lz,ly,lx,lL);
		  int x = lx + procPos[0]*lL[0] - x_src[0];
		  int y = ly + procPos[1]*lL[1] - x_src[1];
		  int z = lz + procPos[2]*lL[2] - x_src[2];
		  tmp = (((double) momElem[imom][0]*x)/L[0] + ((double) momElem[imom][1]*y)/L[1] + ((double) momElem[imom][2]*z)/L[2])*2*PI;
		  Complex C2 (cos(tmp),sign*sin(tmp));
		  corr[it*Nmom*4*4+imom*4*4+gamma*4+gammap] = corr[it*Nmom*4*4+imom*4*4+gamma*4+gammap] + M[v][i][gamma][gammap] * C2;
		}
      
      } // momenta
    } // time
    MPI_Reduce(corr, F_tmp, L[3]*Nmom*4*4*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for(int it = 0 ; it < L[3] ; it++)
      for(int imom = 0 ; imom < Nmom ; imom++)
	for(int gamma = 0 ; gamma < 4 ; gamma++)
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    F[it*Nmom+imom][i][gamma][gammap] = F_tmp[it*Nmom*4*4+imom*4*4+gamma*4+gammap];
  }
  free(F_tmp);
  free(corr);

}

void Baryon_Dec::tabulateIndices_Delta(Complex val[3][4*4],int ind[3][4][4*4],int *counter){

  counter[0] = 0;  counter[1] = 0;   counter[2] = 0;

  Complex C_temp;

  for(int i = 0 ; i < 3 ; i++)
    for(int alpha = 0 ; alpha < 4 ; alpha++)
      for(int beta = 0 ; beta < 4 ; beta++)
	for(int betap = 0 ; betap < 4 ; betap++)
	  for(int alphap = 0 ; alphap < 4 ; alphap++){
	    C_temp = CgammaM[i][alpha][beta] * CgammaM_bar[i][betap][alphap];
	    if( norm(C_temp) > 1e-3 ){
	      val[i][counter[i]] = C_temp;
	      ind[i][0][counter[i]] = alpha;
	      ind[i][1][counter[i]] = beta;
	      ind[i][2][counter[i]] = betap;
	      ind[i][3][counter[i]] = alphap;
	      counter[i]++;
	    }
	  }
  
}


void Baryon_Dec::contract(Propagator &prop1,Propagator &prop2,BARYON_DEC_OPERATOR OPER){
  int counter[3];
  Complex val[3][4*4];
  int ind[3][4][4*4];

  zero_M();
  tabulateIndices_Delta(val,ind,counter);


  int a,b,c,a1,b1,c1;
  switch(OPER){
  case DELTA_ISO1:
    for(int i = 0 ; i < 3 ; i++)
      for(int it = 0 ; it < L[3] ; it++)
	for(int gamma = 0 ; gamma < 4 ; gamma++)       
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    for(int idx = 0 ; idx < counter[i] ; idx++){
	      int alpha = ind[i][0][idx];
	      int beta = ind[i][1][idx];
	      int betap = ind[i][2][idx];
	      int alphap = ind[i][3][idx];
	      for(int cc1 = 0 ; cc1 < 6 ; cc1++){
		a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
		for(int cc2 = 0 ; cc2 < 6 ; cc2++){
		  a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
		  Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[i][idx];
		  for(int lz = 0 ; lz < lL[2] ; lz++)
		    for(int ly = 0 ; ly < lL[1] ; ly++)
		      for(int lx = 0 ; lx < lL[0] ; lx++){
			int iv = LEXIC(it,lz,ly,lx,lL);
			M[iv][i][gamma][gammap] = M[iv][i][gamma][gammap] + factor * 
			  (
			   -prop1.M[iv][alpha][gammap][a][c1]*prop1.M[iv][beta][betap][b][b1]*prop1.M[iv][gamma][alphap][c][a1]
			   +prop1.M[iv][alpha][betap][a][b1]*prop1.M[iv][beta][gammap][b][c1]*prop1.M[iv][gamma][alphap][c][a1]
			   +prop1.M[iv][alpha][gammap][a][c1]*prop1.M[iv][beta][alphap][b][a1]*prop1.M[iv][gamma][betap][c][b1]
			   -prop1.M[iv][alpha][alphap][a][a1]*prop1.M[iv][beta][gammap][b][c1]*prop1.M[iv][gamma][betap][c][b1]
			   -prop1.M[iv][alpha][betap][a][b1]*prop1.M[iv][beta][alphap][b][a1]*prop1.M[iv][gamma][gammap][c][c1]
			   +prop1.M[iv][alpha][alphap][a][a1]*prop1.M[iv][beta][betap][b][b1]*prop1.M[iv][gamma][gammap][c][c1]
			   );
		      }
		}
	      }
	    }
    break;
  case DELTA_ISO1O2:

    for(int i = 0 ; i < 3 ; i++)
      for(int it = 0 ; it < L[3] ; it++)
	for(int gamma = 0 ; gamma < 4 ; gamma++)       
	  for(int gammap = 0 ; gammap < 4 ; gammap++)
	    for(int idx = 0 ; idx < counter[i] ; idx++){
	      int alpha = ind[i][0][idx];
	      int beta = ind[i][1][idx];
	      int betap = ind[i][2][idx];
	      int alphap = ind[i][3][idx];
	      for(int cc1 = 0 ; cc1 < 6 ; cc1++){
		a = eps[cc1][0]; b = eps[cc1][1]; c = eps[cc1][2];
		for(int cc2 = 0 ; cc2 < 6 ; cc2++){
		  a1 = eps[cc2][0]; b1 = eps[cc2][1]; c1 = eps[cc2][2];
		  Complex factor = sgn_eps[cc1] * sgn_eps[cc2] * val[i][idx];
		  for(int lz = 0 ; lz < lL[2] ; lz++)
		    for(int ly = 0 ; ly < lL[1] ; ly++)
		      for(int lx = 0 ; lx < lL[0] ; lx++){
			int iv = LEXIC(it,lz,ly,lx,lL);
			M[iv][i][gamma][gammap] = M[iv][i][gamma][gammap] + factor * (1./3.) * 
			  (
			   -4*prop1.M[iv][alpha][gammap][a][c1]*prop2.M[iv][beta][betap][b][b1]*prop1.M[iv][gamma][alphap][c][a1]
			   +2*prop1.M[iv][alpha][betap][a][b1]*prop2.M[iv][beta][gammap][b][c1]*prop1.M[iv][gamma][alphap][c][a1]
			   +2*prop1.M[iv][alpha][gammap][a][c1]*prop1.M[iv][beta][alphap][b][a1]*prop2.M[iv][gamma][betap][c][b1]
			   -2*prop1.M[iv][alpha][alphap][a][a1]*prop1.M[iv][beta][gammap][b][c1]*prop2.M[iv][gamma][betap][c][b1]
			   -2*prop1.M[iv][alpha][alphap][a][a1]*prop2.M[iv][beta][gammap][b][c1]*prop1.M[iv][gamma][betap][c][b1]
			   -prop1.M[iv][alpha][betap][a][b1]*prop1.M[iv][beta][alphap][b][a1]*prop2.M[iv][gamma][gammap][c][c1]
			   +prop1.M[iv][alpha][alphap][a][a1]*prop1.M[iv][beta][betap][b][b1]*prop2.M[iv][gamma][gammap][c][c1]
			   +4*prop1.M[iv][alpha][alphap][a][a1]*prop2.M[iv][beta][betap][b][b1]*prop1.M[iv][gamma][gammap][c][c1]
			   );
		      }
		}
	      }
	    }

    break;
  default:
    errorQKXTM("Error with the switch\n");
  }



}

void Baryon_Dec::dumpData(BARYON_DEC_OPERATOR OPER){

  if( rank == 0 ){
    for(int it = 0 ; it < L[3] ; it++)
      for(int imom = 0 ; imom < Nmom ; imom++)
	for(int i = 0 ; i < 3 ; i++)
	  for(int gamma = 0 ; gamma < 4 ; gamma++)
	    for(int gammap = 0 ; gammap < 4 ; gammap++)
	      fprintf(ptr_file,"%d \t %d \t %+d %+d %+d \t %d \t %d %d \t %+e %+e\n",OPER,it,momElem[imom][0],momElem[imom][1],momElem[imom][2],i,gamma,gammap,F[it*Nmom+imom][i][gamma][gammap].real(),F[it*Nmom+imom][i][gamma][gammap].imag());
  }

}


//////////////////////////////////// FUNCTIONS 

void contractMesons(Propagator *prop1,Propagator *prop2, char *twop_filename, LatticeInfo *latInfo){

  Meson *meson = new Meson(latInfo,twop_filename);
  enum MESONS_OPERATOR Oper[10]={ONE,G5,G1,G2,G3,G4,G5G1,G5G2,G5G3,G5G4};

  for(int i = 0 ; i < 10 ; i++){
    meson->contract(*prop1,Oper[i]);
    meson->fourier(FORWARD);
    meson->dumpData(Oper[i]);
  }

  for(int i = 0 ; i < 10 ; i++){
    meson->contract(*prop2,Oper[i]);
    meson->fourier(FORWARD);
    meson->dumpData(Oper[i]);
  }

  delete meson;
}

void contractBaryons(Propagator *prop1,Propagator *prop2, char *twop_filename, LatticeInfo *latInfo){


  Baryon *baryon = new Baryon(latInfo,twop_filename);
  //  baryon->contract(*prop1,*prop2,NTN);

  enum BARYON_OPERATOR Oper[4]={NTN,NTR,RTN,RTR};
  
  for(int i = 0 ; i < 4 ; i++){
    baryon->contract(*prop1,*prop2,Oper[i]);
    baryon->fourier(FORWARD);
    baryon->dumpData(Oper[i]);
  }

  for(int i = 0 ; i < 4 ; i++){
    baryon->contract(*prop2,*prop1,Oper[i]);
    baryon->fourier(FORWARD);
    baryon->dumpData(Oper[i]);
  }

  delete baryon;
}

void contractBaryonsDec(Propagator *prop1,Propagator *prop2, char *twop_filename, LatticeInfo *latInfo){


  Baryon_Dec *baryon_dec = new Baryon_Dec(latInfo,twop_filename);
  //  baryon->contract(*prop1,*prop2,NTN);

  enum BARYON_DEC_OPERATOR Oper[2]={DELTA_ISO1,DELTA_ISO1O2};
  
  for(int i = 0 ; i < 2 ; i++){
    baryon_dec->contract(*prop1,*prop2,Oper[i]);
    baryon_dec->fourier(FORWARD);
    baryon_dec->dumpData(Oper[i]);
  }

  for(int i = 0 ; i < 2 ; i++){
    baryon_dec->contract(*prop2,*prop1,Oper[i]);
    baryon_dec->fourier(FORWARD);
    baryon_dec->dumpData(Oper[i]);
  }

  delete baryon_dec;
}

void contract_Mesons_omp(Propagator *prop1,Propagator *prop2, char *twop_filename, LatticeInfo *latInfo){

  int rank = prop1->rank;
  printfQKXTM("You use a function with omp support\n");
  int n_cores = omp_get_num_procs();
  printfQKXTM("The system you are running has %d cores\n",n_cores);
  int max_n_threads = omp_get_max_threads();
  printfQKXTM("You asked for %d threads\n",max_n_threads);

}
