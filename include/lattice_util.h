#ifndef _LATTICE_UTIL_H
#define _LATTICE_UTIL_H

#include <qkxTM.h>

// need to create a class which include all information about the lattice


namespace qkxTM {



  // forward definition
  class LatticeGeometry;
  class Gauge;
  class Vector;
  class Propagator;
  class Meson;
  class Baryon;
  class Baryon_Dec;
  class FieldStrength;

  enum PARITY{EO,OE,OO,EE};

  class LatticeGeometry {

  public:
    long int V4d;
    long int V3d;
    long int lV4d;
    long int lV3d;
    int L[NDIM];
    int P[NDIM];
    int lL[NDIM];
    int procPos[NDIM];        // position of processor in processors grid
    int rank;                    // id of my processor
    int Nprocs;                   // total number of processors
    
    int procPlus[NDIM];          // id of the next processor in each direction
    int procMinus[NDIM];          // id of the previous processor in each direction
    
    int *pointPlus[NDIM];   // index of the forward point 
    int *pointMinus[NDIM];   // index of the backward point

    int placeBoundaryStartPlus[NDIM];   // to find the pointer where the boundary starts
    int placeBoundaryStartMinus[NDIM];   // to find the pointer where the boundary starts
    
    int lSurfaceBoundary[NDIM];

    int NtotalBoundary;

    MPI_Datatype	stypeV[NDIM];	// data-types that allow sends of strided data       
    MPI_Datatype	rtypeV[NDIM];	//      vector-fields 
    MPI_Datatype        stypeU[NDIM];
    MPI_Datatype	rtypeU[NDIM];	//      gauge-fields
    MPI_Datatype	stypeP[NDIM];	//
    MPI_Datatype	rtypeP[NDIM];	//      propagator-fields                                 

    bool initialized;                   // check initialization
    
    // constructor
    LatticeGeometry(LatticeInfo *latInfo);
    // destructor

    virtual ~LatticeGeometry(){
      for(int i  = 0 ; i < NDIM ; i++){
	free(pointPlus[i]);
	free(pointMinus[i]);
      }
    }

    Complex gammaM[5][4][4];
    Complex Ccg[4][4];
    Complex Ccg_bar[4][4];
    Complex CgammaM[5][4][4];
    Complex CgammaM_bar[5][4][4];
  };

  ///////////////////////////////////////
  class Vector : public LatticeGeometry {

  public:
    long int size;
    long int size_total;
  

    Vector(LatticeInfo *latInfo);
    virtual ~Vector();
    
    Complex (*M)[NSPINS][NCOLORS];          // elements of the vector with ghost zone
    void *ptr_boundaryPlusStart[4];
    void *ptr_boundaryMinusStart[4];

    void zero();
    void copy(Vector &src);
    double norm2();
    void normalize();
    void communicateToPlus();
    void communicateToMinus();

    void Gaussian_smearing(Vector &vec_tmp,Gauge &gauge, LatticeInfo *latInfo);
    void applyDslash(Vector &psi, Gauge &u);
    void applyDslash_chiral(Vector &psi, Gauge &u);
    void applyDslash_xx(Vector &psi, Gauge &u, PARITY parity);
    void applyDslash_xx_chiral(Vector &psi, Gauge &u, PARITY parity);
    void applyDeltaAddDslash(Vector &dslash, Vector &psi, LatticeInfo *latInfo);
    void applyTwistAddDslash(Vector &dslash, Vector &psi, LatticeInfo *latInfo);
    void applyTwistAddDslash_chiral(Vector &dslash, Vector &psi, LatticeInfo *latInfo);
    void applyTwistAddDslash_PC(Vector &dslash, Vector &psi, LatticeInfo *latInfo, PARITY parity);
    void applyTwistAddDslash_PC_chiral(Vector &dslash, Vector &psi, LatticeInfo *latInfo, PARITY parity);

    void applyDslashDag(Vector &psi, Gauge &u);
    void applyTwistDagAddDslashDag(Vector &dslashDag, Vector &psi, LatticeInfo *latInfo);
    void applyTwist(Vector &psi, LatticeInfo *latInfo);
    void applyTwistDag(Vector &psi, LatticeInfo *latInfo);
    void applyTwistInv(Vector &psi, LatticeInfo *latInfo);
    void applyTwistInv_chiral(Vector &psi, LatticeInfo *latInfo);
    void applyTwistDagInv(Vector &psi, LatticeInfo *latInfo);
    void applyGamma5(Vector &psi, LatticeInfo *latInfo);
    void applyGamma5_chiral(Vector &psi, LatticeInfo *latInfo);
    void daggerVectorGamma5(Vector &psi, LatticeInfo *latInfo);
    
    void rotateFromChiralToUKQCD();
    void multiply_by_phase();
    void sendEvenToOdd();

    void tm(Vector &psi, Gauge &u, LatticeInfo *latInfo);
    void tm_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo);
    void tmDag(Vector &psi, Gauge &u, LatticeInfo *latInfo);
    void tmDagtm(Vector &psi,Vector &vec_tmp ,Gauge &u, LatticeInfo *latInfo);
    void Wilson(Vector &psi, Gauge &u, LatticeInfo *latInfo);

    void tm_xx(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity);
    void tm_xxDagTm_xx(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity);

    void tm_xx_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity);
    void tm_xxDagTm_xx_chiral(Vector &psi, Gauge &u, LatticeInfo *latInfo, PARITY parity);

    void seqSourceFixSinkPart1(Propagator &prop1,Propagator &prop2,int timeslice,int nu, int c2, int flagProj);
    void seqSourceFixSinkPart2(Propagator &prop1,int timeslice,int nu, int c2, int flagProj);
    void rescale(double a);
    void applyConjugation();
    void applyg5gm(int mu);
    void copyPropagator(Propagator &prop,int nu , int c2);

    void standardEndTrick_sigmaTerm(Complex *loops, LatticeInfo *latInfo );
    void generalizedEndTrick_axial(Gauge &u,Complex *loops, LatticeInfo *latInfo );
    void generalizedEndTrick_axial_openIndices(Gauge &u,Complex *loops, LatticeInfo *latInfo );
    void generalizedEndTrick_covDer(Gauge &u,Complex *loops, Complex *loops_cv, LatticeInfo *latInfo );
    void standardEndTrick_covDer(Gauge &u,Complex *loops, Complex *loops_cv , LatticeInfo *latInfo );
    void volumeSource_ultralocal(Vector &xi, Gauge &u,Complex *loops, LatticeInfo *latInfo );
    void volumeSource_covDer(Vector &xi, Gauge &u, Complex *loops, Complex *loops_cv, LatticeInfo *latInfo );

    void applyCloverTerm(Vector &psi, FieldStrength &fs, LatticeInfo *latInfo);
    void Clover(Vector &psi, Gauge &u, FieldStrength &fs, LatticeInfo *latInfo);
  };
  //////////////////////////////////
  class Propagator : public LatticeGeometry {

  public:
    long int size;
    long int size_total;


    Propagator(LatticeInfo *latInfo);
    virtual ~Propagator();

    Complex (*M)[NSPINS][NSPINS][NCOLORS][NCOLORS];
    void *ptr_boundaryPlusStart[4];
    void *ptr_boundaryMinusStart[4];

    void zero();
    void copy(Propagator &src);
    void copyVector(Vector &v,int nu , int c2);
    void communicateToPlus();
    double norm2();
    void communicateToMinus();

    void absorbVector(Vector &phi ,int nu , int c2 );

    void rotateToPhysicalBasePlus();
    void rotateToPhysicalBaseMinus();

  };
  ////////////////////////////////////////////////////
  class Gauge : public LatticeGeometry {

  public:
    long int size;
    long int size_total;
    bool commToPlusOnce;  // flag for constGauge
    bool commToMinusOnce; // flag for constGauge
    bool constGauge; //         flag if its true means already we have done comms

    Gauge(LatticeInfo *latInfo);
    virtual ~Gauge();

    Complex (*M)[NDIM][NCOLORS][NCOLORS];
    void *ptr_boundaryPlusStart[4];
    void *ptr_boundaryMinusStart[4];

    void zero();
    void copy(Gauge &src);
    double norm2();
    void applyBoundaryCondition(LatticeInfo *latInfo);
    void communicateToPlus();
    void communicateToMinus();

    double calculatePlaquette();
    void APE_smearing( Gauge &tmp , LatticeInfo *latInfo);
    void calculateGluonLoop(Complex *gluonLoop);
  };

  class FieldStrength : public LatticeGeometry {
  public:
    int size;
    FieldStrength(LatticeInfo *latInfo, Gauge &u);
    virtual ~FieldStrength();
    Complex (*M)[6][NCOLORS][NCOLORS];
    void zero();
  };

  enum FOURIER_DIRECTION{FORWARD,BACKWARD};
  enum MESONS_OPERATOR{ONE=0,G5,G1,G2,G3,G4,G5G1,G5G2,G5G3,G5G4};
  // for example ONE is the scalar meson and G5 is the pion
  // Giannis for ONE is pion !! attention
  ///////////////////////////////////////////////////
  class Meson : public LatticeGeometry {

  public:
    Complex *M;
    Complex *F;

    void zero_M();
    void zero_F();
    int tabulateIndices(Complex val[],int ind[][4], MESONS_OPERATOR OPER);

    int size;
    int Q_sq;
    int x_src[4];
    int Nmom;
    int momElem[MOM_MAX][3];
    FILE *ptr_file;

    Meson(LatticeInfo *latInfo,char filename_twop[]);
    Meson(LatticeInfo *latInfo);
    virtual ~Meson();

    void contract(Propagator &prop,MESONS_OPERATOR);
    void fourier(FOURIER_DIRECTION);
    void dumpData(MESONS_OPERATOR);


  };

  ///////////////////////////////////////////////////
  enum BARYON_OPERATOR{NTN=0,NTR,RTN,RTR};


  class Baryon : public LatticeGeometry {

  public:
    int size;
    int Q_sq;
    int x_src[4];
    int Nmom;
    int momElem[MOM_MAX][3];
    FILE *ptr_file;

    void contract(Propagator &prop1,Propagator &prop2,BARYON_OPERATOR);
    int tabulateIndices_NTN(Complex val[4*4*4*4],int ind[8][4*4*4*4]);
    int tabulateIndices_NTR(Complex val[4*4*4*4],int ind[8][4*4*4*4]);
    int tabulateIndices_RTN(Complex val[4*4*4*4],int ind[8][4*4*4*4]);
    int tabulateIndices_RTR(Complex val[4*4*4*4],int ind[8][4*4*4*4]);
    Complex (*M)[4][4];
    Complex (*F)[4][4];

    Baryon(LatticeInfo *latInfo, char filename_twop[]);
    virtual ~Baryon();
    
    void zero_M();
    void zero_F();
    void fourier(FOURIER_DIRECTION);
    
    void dumpData(BARYON_OPERATOR);
  };

  enum BARYON_DEC_OPERATOR{DELTA_ISO1=0,DELTA_ISO1O2};

  class Baryon_Dec : public LatticeGeometry {

  public:
    int size;
    int Q_sq;
    int x_src[4];
    int Nmom;
    int momElem[MOM_MAX][3];
    FILE *ptr_file;

    Baryon_Dec(LatticeInfo *latInfo, char filename_twop[]);
    virtual ~Baryon_Dec();

    Complex (*M)[3][4][4];
    Complex (*F)[3][4][4];
    
    void zero_M();
    void zero_F();
    void fourier(FOURIER_DIRECTION);
    
    void dumpData(BARYON_DEC_OPERATOR);

    void contract(Propagator &prop1,Propagator &prop2,BARYON_DEC_OPERATOR);
    void tabulateIndices_Delta(Complex val[3][4*4],int ind[3][4][4*4],int *);
  };

} // close namespace

  // functions
int qkxTM_MPI_getGaugeLime(char *fname, qkxTM::Gauge &u);
int qkxTM_writeGaugeLime(char *fname, qkxTM::Gauge &u, char* message);
int qkxTM_MPI_getVectorLime(char *fname, qkxTM::Vector &v);
int qkxTM_MPI_getPropLime(char *fname, qkxTM::Propagator &prop, LatticeInfo *latInfo);
void qkxTM_MPI_readEigenVectors(char *fname, qkxTM::Vector &v);

void qcd_projectSU33d(qkxTM::Gauge &gf);
void contractNucleon(qkxTM::Propagator &prop1, qkxTM::Propagator &prop2, char *output, int Q_sq, int x_src[4]);
void contractMesons(qkxTM::Propagator *prop1,qkxTM::Propagator *prop2, char *twop_filename, LatticeInfo *latInfo);
void contractBaryons(qkxTM::Propagator *prop1,qkxTM::Propagator *prop2, char *twop_filename, LatticeInfo *latInfo);
void contractBaryonsDec(qkxTM::Propagator *prop1,qkxTM::Propagator *prop2, char *twop_filename, LatticeInfo *latInfo);
void contract_mesons_cache(qkxTM::Propagator &prop1, qkxTM::Propagator &prop2, char *filename_twop_mesons, LatticeInfo latInfo);
void contract_Mesons_omp(qkxTM::Propagator *prop1,qkxTM::Propagator *prop2, char *twop_filename, LatticeInfo *latInfo);
#endif
