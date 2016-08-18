/*
void contractBaryonsDec(Propagator *prop1,Propagator *prop2, char *twop_filename, LatticeInfo *latInfo){

  Baryon_Dec *baryon = new Baryon_Dec(latInfo,twop_filename);
  enum BARYON_DEC_OPERATOR Oper[4]={DELTA_ISO1,DELTA_ISO1O2};
  
  for(int i = 0 ; i < 2 ; i++){
    baryon->contract(*prop1,*prop2,Oper[i]);
    baryon->fourier(FORWARD);
    baryon->dumpData(Oper[i]);
  }

  for(int i = 0 ; i < 2 ; i++){
    baryon->contract(*prop2,*prop1,Oper[i]);
    baryon->fourier(FORWARD);
    baryon->dumpData(Oper[i]);
  }

  delete baryon;
}

#define N_MESONS 10

static __inline__ Complex prop_tr(Complex prop[4][4][3][3]){
  Complex trace;
  trace.real() = 0.; trace.imag() = 0.;

#pragma unroll_and_jam (4)
  for(int mu = 0 ; mu < 4 ; mu++)
#pragma unroll_and_jam (3)
    for(int ic = 0 ; ic < 3 ; ic++){
      trace = trace + prop[mu][mu][ic][ic];
    }
  return trace;
}

void contract_mesons_cache(Propagator &prop1, Propagator &prop2, char *filename_twop_mesons, LatticeInfo latInfo){

  if(initializeGroupsFlag == false)errorQKXTM("You must initialize the group Communicator first\n");
    
  
  register Complex pp1[4][4][3][3];
  register Complex pp2[4][4][3][3];

  register Complex pp1_dag[4][4][3][3];
  register Complex pp2_dag[4][4][3][3];

  register Complex pp_temp1[4][4][3][3];
  register Complex pp_temp2[4][4][3][3];

  register Complex pp[4][4][3][3];
  register Complex aux[4][4][3][3];

  register Complex trace;
  register Complex mesons[2*N_MESONS];

  int momElem[MOM_MAX][3];
  int Nmom = createMomenta(momElem,latInfo.Q_sq);

  Complex *corr_mom = (Complex*)calloc(latInfo.L[3]*2*N_MESONS*Nmom,sizeof(Complex));
  Complex *corr_mom_local = (Complex*)calloc(prop1.lL[3]*2*N_MESONS*Nmom,sizeof(Complex));
  Complex *corr_mom_local_reduced = (Complex*)calloc(prop1.lL[3]*2*N_MESONS*Nmom,sizeof(Complex));

  for(int it = 0 ; it < prop1.lL[3] ; it++){

    for(int iz = 0 ; iz < prop1.lL[2] ; iz++)
      for(int iy = 0 ; iy < prop1.lL[1] ; iy++)
	for(int ix = 0 ; ix < prop1.lL[0] ; ix++){
	  int iv = LEXIC(it,iz,iy,ix,prop1.lL);
	  memcpy(pp1,&(prop1.M[iv][0][0][0][0]),144*2*sizeof(double)); // store to registers
	  memcpy(pp2,&(prop2.M[iv][0][0][0][0]),144*2*sizeof(double));
	  prop_G_dag(pp1_dag,pp1);
	  prop_G_dag(pp2_dag,pp2);

	  ///////////////// uprop //////////////////////////////
	  // ONE 
	  memcpy(pp_temp1,pp1,144*2*sizeof(double));
	  memcpy(pp_temp2,pp1_dag,144*2*sizeof(double));
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[0] = trace;
	  // G5 
	  prop_gamma_5_G(pp_temp1,pp1);
	  prop_gamma_5_G(pp_temp2,pp1_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[1] = trace;
	  // G1
	  prop_gamma_x_G(pp_temp1,pp1);
	  prop_gamma_x_G(pp_temp2,pp1_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[2] = trace;
	  // G2
	  prop_gamma_y_G(pp_temp1,pp1);
	  prop_gamma_y_G(pp_temp2,pp1_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[3] = trace;
	  // G3
	  prop_gamma_z_G(pp_temp1,pp1);
	  prop_gamma_z_G(pp_temp2,pp1_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[4] = trace;
	  // G4
	  prop_gamma_t_G(pp_temp1,pp1);
	  prop_gamma_t_G(pp_temp2,pp1_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[5] = trace;
	  // G5G1
	  prop_gamma_x_G(aux,pp1);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_x_G(aux,pp1_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[6] = trace;
	  // G5G2
	  prop_gamma_y_G(aux,pp1);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_y_G(aux,pp1_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[7] = trace;
	  // G5G3
	  prop_gamma_z_G(aux,pp1);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_z_G(aux,pp1_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[8] = trace;
	  // G5G4
	  prop_gamma_t_G(aux,pp1);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_t_G(aux,pp1_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[9] = trace;
	  ///////////////////////////////////////////////////

	  ///////////////// uprop //////////////////////////////
	  // ONE 
	  memcpy(pp_temp1,pp2,144*2*sizeof(double));
	  memcpy(pp_temp2,pp2_dag,144*2*sizeof(double));
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[10] = trace;
	  // G5 
	  prop_gamma_5_G(pp_temp1,pp2);
	  prop_gamma_5_G(pp_temp2,pp2_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[11] = trace;
	  // G1
	  prop_gamma_x_G(pp_temp1,pp2);
	  prop_gamma_x_G(pp_temp2,pp2_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[12] = trace;
	  // G2
	  prop_gamma_y_G(pp_temp1,pp2);
	  prop_gamma_y_G(pp_temp2,pp2_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[13] = trace;
	  // G3
	  prop_gamma_z_G(pp_temp1,pp2);
	  prop_gamma_z_G(pp_temp2,pp2_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[14] = trace;
	  // G4
	  prop_gamma_t_G(pp_temp1,pp2);
	  prop_gamma_t_G(pp_temp2,pp2_dag);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[15] = trace;
	  // G5G1
	  prop_gamma_x_G(aux,pp2);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_x_G(aux,pp2_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[16] = trace;
	  // G5G2
	  prop_gamma_y_G(aux,pp2);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_y_G(aux,pp2_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[17] = trace;
	  // G5G3
	  prop_gamma_z_G(aux,pp2);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_z_G(aux,pp2_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[18] = trace;
	  // G5G4
	  prop_gamma_t_G(aux,pp2);
	  prop_gamma_5_G(pp_temp1,aux);
	  prop_gamma_t_G(aux,pp2_dag);
	  prop_gamma_5_G(pp_temp2,aux);
	  prop_G_G(pp, pp_temp1 ,pp_temp2);
	  trace=prop_tr(pp);
	  mesons[19] = trace;
	  ///////////////////////////////////////////////////
	  // fourier
	  int x = ix + prop1.procPos[0]*prop1.lL[0] - latInfo.x_src[0];
	  int y = iy + prop1.procPos[1]*prop1.lL[1] - latInfo.x_src[1];
	  int z = iz + prop1.procPos[2]*prop1.lL[2] - latInfo.x_src[2];
	  
#pragma unroll_and_jam(2*N_MESONS)
	  for(int ip = 0 ; ip < 2*N_MESONS ; ip++){
	    for(int imom = 0 ; imom < Nmom ; imom++){
	      double tmp = (((double) momElem[imom][0]*x)/latInfo.L[0] + ((double) momElem[imom][1]*y)/latInfo.L[1] + ((double) momElem[imom][2]*z)/latInfo.L[2])*2*PI;
	      Complex C2 (cos(tmp),-sin(tmp));
	      corr_mom_local[it*2*N_MESONS*Nmom+ip*Nmom+imom] = corr_mom_local[it*2*N_MESONS*Nmom+ip*Nmom+imom] + mesons[ip]*C2;
	    }
	  }

	} // close spatial loop
  } // close time loop

  // collect from all process
  MPI_Reduce(corr_mom_local,corr_mom_local_reduced,prop1.lL[3]*2*N_MESONS*Nmom*2,MPI_DOUBLE,MPI_SUM,0, spaceComm);
  if( timeRank >= 0 && timeRank < latInfo.P[3] ){
    MPI_Gather(corr_mom_local_reduced,prop1.lL[3]*2*N_MESONS*Nmom*2,MPI_DOUBLE,corr_mom,prop1.lL[3]*2*N_MESONS*Nmom*2,MPI_DOUBLE,0,timeComm);
  }

  FILE *ptr_out;
  if(prop1.rank == 0){
    ptr_out = fopen(filename_twop_mesons,"w");
    for(int ip = 0 ; ip < 2*N_MESONS ; ip++)
      for(int it = 0 ; it < latInfo.L[3] ; it++)
	for(int imom = 0 ; imom < Nmom ; imom++)
	  fprintf(ptr_out,"%d \t %d \t %+d %+d %+d \t %+e %+e\n",ip,it,momElem[imom][0],momElem[imom][1],momElem[imom][2],corr_mom[it*2*N_MESONS*Nmom+ip*Nmom+imom].real(),corr_mom[it*2*N_MESONS*Nmom+ip*Nmom+imom].imag());
    fclose(ptr_out);
  }

  free(corr_mom);
  free(corr_mom_local);
  free(corr_mom_local_reduced);
}

*/
