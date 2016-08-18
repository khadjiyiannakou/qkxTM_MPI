#define PROJ_MINUS_X_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointPlus[0][iv]][0][0].real() - psi.M[pointPlus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() =psi.M[pointPlus[0][iv]][0][0].imag() + psi.M[pointPlus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() =psi.M[pointPlus[0][iv]][0][1].real() - psi.M[pointPlus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() =psi.M[pointPlus[0][iv]][0][1].imag() + psi.M[pointPlus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() =psi.M[pointPlus[0][iv]][0][2].real() - psi.M[pointPlus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() =psi.M[pointPlus[0][iv]][0][2].imag() + psi.M[pointPlus[0][iv]][3][2].real(); \
  phi[MS(1,0)].real() =psi.M[pointPlus[0][iv]][1][0].real() - psi.M[pointPlus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() =psi.M[pointPlus[0][iv]][1][0].imag() + psi.M[pointPlus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() =psi.M[pointPlus[0][iv]][1][1].real() - psi.M[pointPlus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() =psi.M[pointPlus[0][iv]][1][1].imag() + psi.M[pointPlus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() =psi.M[pointPlus[0][iv]][1][2].real() - psi.M[pointPlus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() =psi.M[pointPlus[0][iv]][1][2].imag() + psi.M[pointPlus[0][iv]][2][2].real(); 

#define PROJ_PLUS_X_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointMinus[0][iv]][0][0].real() + psi.M[pointMinus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() =psi.M[pointMinus[0][iv]][0][0].imag() - psi.M[pointMinus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() =psi.M[pointMinus[0][iv]][0][1].real() + psi.M[pointMinus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() =psi.M[pointMinus[0][iv]][0][1].imag() - psi.M[pointMinus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() =psi.M[pointMinus[0][iv]][0][2].real() + psi.M[pointMinus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() =psi.M[pointMinus[0][iv]][0][2].imag() - psi.M[pointMinus[0][iv]][3][2].real(); \
  phi[MS(1,0)].real() =psi.M[pointMinus[0][iv]][1][0].real() + psi.M[pointMinus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() =psi.M[pointMinus[0][iv]][1][0].imag() - psi.M[pointMinus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() =psi.M[pointMinus[0][iv]][1][1].real() + psi.M[pointMinus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() =psi.M[pointMinus[0][iv]][1][1].imag() - psi.M[pointMinus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() =psi.M[pointMinus[0][iv]][1][2].real() + psi.M[pointMinus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() =psi.M[pointMinus[0][iv]][1][2].imag() - psi.M[pointMinus[0][iv]][2][2].real(); 

#define PROJ_MINUS_Y_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointPlus[1][iv]][0][0].real() + psi.M[pointPlus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() =psi.M[pointPlus[1][iv]][0][0].imag() + psi.M[pointPlus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() =psi.M[pointPlus[1][iv]][0][1].real() + psi.M[pointPlus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() =psi.M[pointPlus[1][iv]][0][1].imag() + psi.M[pointPlus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() =psi.M[pointPlus[1][iv]][0][2].real() + psi.M[pointPlus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() =psi.M[pointPlus[1][iv]][0][2].imag() + psi.M[pointPlus[1][iv]][3][2].imag(); \
  phi[MS(1,0)].real() =psi.M[pointPlus[1][iv]][1][0].real() - psi.M[pointPlus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() =psi.M[pointPlus[1][iv]][1][0].imag() - psi.M[pointPlus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() =psi.M[pointPlus[1][iv]][1][1].real() - psi.M[pointPlus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() =psi.M[pointPlus[1][iv]][1][1].imag() - psi.M[pointPlus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() =psi.M[pointPlus[1][iv]][1][2].real() - psi.M[pointPlus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() =psi.M[pointPlus[1][iv]][1][2].imag() - psi.M[pointPlus[1][iv]][2][2].imag(); 
 
#define PROJ_PLUS_Y_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointMinus[1][iv]][0][0].real() - psi.M[pointMinus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() =psi.M[pointMinus[1][iv]][0][0].imag() - psi.M[pointMinus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() =psi.M[pointMinus[1][iv]][0][1].real() - psi.M[pointMinus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() =psi.M[pointMinus[1][iv]][0][1].imag() - psi.M[pointMinus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() =psi.M[pointMinus[1][iv]][0][2].real() - psi.M[pointMinus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() =psi.M[pointMinus[1][iv]][0][2].imag() - psi.M[pointMinus[1][iv]][3][2].imag(); \
  phi[MS(1,0)].real() =psi.M[pointMinus[1][iv]][1][0].real() + psi.M[pointMinus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() =psi.M[pointMinus[1][iv]][1][0].imag() + psi.M[pointMinus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() =psi.M[pointMinus[1][iv]][1][1].real() + psi.M[pointMinus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() =psi.M[pointMinus[1][iv]][1][1].imag() + psi.M[pointMinus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() =psi.M[pointMinus[1][iv]][1][2].real() + psi.M[pointMinus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() =psi.M[pointMinus[1][iv]][1][2].imag() + psi.M[pointMinus[1][iv]][2][2].imag(); 

#define PROJ_MINUS_Z_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointPlus[2][iv]][0][0].real() - psi.M[pointPlus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() =psi.M[pointPlus[2][iv]][0][0].imag() + psi.M[pointPlus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() =psi.M[pointPlus[2][iv]][0][1].real() - psi.M[pointPlus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() =psi.M[pointPlus[2][iv]][0][1].imag() + psi.M[pointPlus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() =psi.M[pointPlus[2][iv]][0][2].real() - psi.M[pointPlus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() =psi.M[pointPlus[2][iv]][0][2].imag() + psi.M[pointPlus[2][iv]][2][2].real(); \
  phi[MS(1,0)].real() =psi.M[pointPlus[2][iv]][1][0].real() + psi.M[pointPlus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() =psi.M[pointPlus[2][iv]][1][0].imag() - psi.M[pointPlus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() =psi.M[pointPlus[2][iv]][1][1].real() + psi.M[pointPlus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() =psi.M[pointPlus[2][iv]][1][1].imag() - psi.M[pointPlus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() =psi.M[pointPlus[2][iv]][1][2].real() + psi.M[pointPlus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() =psi.M[pointPlus[2][iv]][1][2].imag() - psi.M[pointPlus[2][iv]][3][2].real(); 
  
#define PROJ_PLUS_Z_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointMinus[2][iv]][0][0].real() + psi.M[pointMinus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() =psi.M[pointMinus[2][iv]][0][0].imag() - psi.M[pointMinus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() =psi.M[pointMinus[2][iv]][0][1].real() + psi.M[pointMinus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() =psi.M[pointMinus[2][iv]][0][1].imag() - psi.M[pointMinus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() =psi.M[pointMinus[2][iv]][0][2].real() + psi.M[pointMinus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() =psi.M[pointMinus[2][iv]][0][2].imag() - psi.M[pointMinus[2][iv]][2][2].real(); \
  phi[MS(1,0)].real() =psi.M[pointMinus[2][iv]][1][0].real() - psi.M[pointMinus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() =psi.M[pointMinus[2][iv]][1][0].imag() + psi.M[pointMinus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() =psi.M[pointMinus[2][iv]][1][1].real() - psi.M[pointMinus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() =psi.M[pointMinus[2][iv]][1][1].imag() + psi.M[pointMinus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() =psi.M[pointMinus[2][iv]][1][2].real() - psi.M[pointMinus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() =psi.M[pointMinus[2][iv]][1][2].imag() + psi.M[pointMinus[2][iv]][3][2].real(); 
 
#define PROJ_MINUS_T_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointPlus[3][iv]][0][0].real() + psi.M[pointPlus[3][iv]][2][0].real(); \
  phi[MS(0,0)].imag() =psi.M[pointPlus[3][iv]][0][0].imag() + psi.M[pointPlus[3][iv]][2][0].imag(); \
  phi[MS(0,1)].real() =psi.M[pointPlus[3][iv]][0][1].real() + psi.M[pointPlus[3][iv]][2][1].real(); \
  phi[MS(0,1)].imag() =psi.M[pointPlus[3][iv]][0][1].imag() + psi.M[pointPlus[3][iv]][2][1].imag(); \
  phi[MS(0,2)].real() =psi.M[pointPlus[3][iv]][0][2].real() + psi.M[pointPlus[3][iv]][2][2].real(); \
  phi[MS(0,2)].imag() =psi.M[pointPlus[3][iv]][0][2].imag() + psi.M[pointPlus[3][iv]][2][2].imag(); \
  phi[MS(1,0)].real() =psi.M[pointPlus[3][iv]][1][0].real() + psi.M[pointPlus[3][iv]][3][0].real(); \
  phi[MS(1,0)].imag() =psi.M[pointPlus[3][iv]][1][0].imag() + psi.M[pointPlus[3][iv]][3][0].imag(); \
  phi[MS(1,1)].real() =psi.M[pointPlus[3][iv]][1][1].real() + psi.M[pointPlus[3][iv]][3][1].real(); \
  phi[MS(1,1)].imag() =psi.M[pointPlus[3][iv]][1][1].imag() + psi.M[pointPlus[3][iv]][3][1].imag(); \
  phi[MS(1,2)].real() =psi.M[pointPlus[3][iv]][1][2].real() + psi.M[pointPlus[3][iv]][3][2].real(); \
  phi[MS(1,2)].imag() =psi.M[pointPlus[3][iv]][1][2].imag() + psi.M[pointPlus[3][iv]][3][2].imag(); 

#define PROJ_PLUS_T_CHIRAL(phi,psi)					\
  phi[MS(0,0)].real() =psi.M[pointMinus[3][iv]][0][0].real() - psi.M[pointMinus[3][iv]][2][0].real(); \
  phi[MS(0,0)].imag() =psi.M[pointMinus[3][iv]][0][0].imag() - psi.M[pointMinus[3][iv]][2][0].imag(); \
  phi[MS(0,1)].real() =psi.M[pointMinus[3][iv]][0][1].real() - psi.M[pointMinus[3][iv]][2][1].real(); \
  phi[MS(0,1)].imag() =psi.M[pointMinus[3][iv]][0][1].imag() - psi.M[pointMinus[3][iv]][2][1].imag(); \
  phi[MS(0,2)].real() =psi.M[pointMinus[3][iv]][0][2].real() - psi.M[pointMinus[3][iv]][2][2].real(); \
  phi[MS(0,2)].imag() =psi.M[pointMinus[3][iv]][0][2].imag() - psi.M[pointMinus[3][iv]][2][2].imag(); \
  phi[MS(1,0)].real() =psi.M[pointMinus[3][iv]][1][0].real() - psi.M[pointMinus[3][iv]][3][0].real(); \
  phi[MS(1,0)].imag() =psi.M[pointMinus[3][iv]][1][0].imag() - psi.M[pointMinus[3][iv]][3][0].imag(); \
  phi[MS(1,1)].real() =psi.M[pointMinus[3][iv]][1][1].real() - psi.M[pointMinus[3][iv]][3][1].real(); \
  phi[MS(1,1)].imag() =psi.M[pointMinus[3][iv]][1][1].imag() - psi.M[pointMinus[3][iv]][3][1].imag(); \
  phi[MS(1,2)].real() =psi.M[pointMinus[3][iv]][1][2].real() - psi.M[pointMinus[3][iv]][3][2].real(); \
  phi[MS(1,2)].imag() =psi.M[pointMinus[3][iv]][1][2].imag() - psi.M[pointMinus[3][iv]][3][2].imag(); 

#define COLLECT_PLUS_X_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
							\
  R[MS(2,0)].real() = R[MS(2,0)].real() - xi[MS(1,0)].imag();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() + xi[MS(1,0)].real();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() - xi[MS(1,1)].imag();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() + xi[MS(1,1)].real();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() - xi[MS(1,2)].imag();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() + xi[MS(1,2)].real();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() - xi[MS(0,0)].imag();  \
  R[MS(3,0)].imag() = R[MS(3,0)].imag() + xi[MS(0,0)].real();  \
  R[MS(3,1)].real() = R[MS(3,1)].real() - xi[MS(0,1)].imag();  \
  R[MS(3,1)].imag() = R[MS(3,1)].imag() + xi[MS(0,1)].real();  \
  R[MS(3,2)].real() = R[MS(3,2)].real() - xi[MS(0,2)].imag();  \
  R[MS(3,2)].imag() = R[MS(3,2)].imag() + xi[MS(0,2)].real();  

#define COLLECT_MINUS_X_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
							\
  R[MS(2,0)].real() = R[MS(2,0)].real() + xi[MS(1,0)].imag();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() - xi[MS(1,0)].real();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() + xi[MS(1,1)].imag();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() - xi[MS(1,1)].real();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() + xi[MS(1,2)].imag();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() - xi[MS(1,2)].real();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() + xi[MS(0,0)].imag();  \
  R[MS(3,0)].imag() = R[MS(3,0)].imag() - xi[MS(0,0)].real();  \
  R[MS(3,1)].real() = R[MS(3,1)].real() + xi[MS(0,1)].imag();  \
  R[MS(3,1)].imag() = R[MS(3,1)].imag() - xi[MS(0,1)].real();  \
  R[MS(3,2)].real() = R[MS(3,2)].real() + xi[MS(0,2)].imag();  \
  R[MS(3,2)].imag() = R[MS(3,2)].imag() - xi[MS(0,2)].real();  


#define COLLECT_PLUS_Y_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
							\
  R[MS(2,0)].real() = R[MS(2,0)].real() + xi[MS(1,0)].real();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() + xi[MS(1,0)].imag();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() + xi[MS(1,1)].real();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() + xi[MS(1,1)].imag();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() + xi[MS(1,2)].real();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() + xi[MS(1,2)].imag();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() - xi[MS(0,0)].real();	\
  R[MS(3,0)].imag() = R[MS(3,0)].imag() - xi[MS(0,0)].imag();	\
  R[MS(3,1)].real() = R[MS(3,1)].real() - xi[MS(0,1)].real();	\
  R[MS(3,1)].imag() = R[MS(3,1)].imag() - xi[MS(0,1)].imag();	\
  R[MS(3,2)].real() = R[MS(3,2)].real() - xi[MS(0,2)].real();	\
  R[MS(3,2)].imag() = R[MS(3,2)].imag() - xi[MS(0,2)].imag();	


#define COLLECT_MINUS_Y_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
							\
  R[MS(2,0)].real() = R[MS(2,0)].real() - xi[MS(1,0)].real();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() - xi[MS(1,0)].imag();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() - xi[MS(1,1)].real();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() - xi[MS(1,1)].imag();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() - xi[MS(1,2)].real();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() - xi[MS(1,2)].imag();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() + xi[MS(0,0)].real();	\
  R[MS(3,0)].imag() = R[MS(3,0)].imag() + xi[MS(0,0)].imag();	\
  R[MS(3,1)].real() = R[MS(3,1)].real() + xi[MS(0,1)].real();	\
  R[MS(3,1)].imag() = R[MS(3,1)].imag() + xi[MS(0,1)].imag();	\
  R[MS(3,2)].real() = R[MS(3,2)].real() + xi[MS(0,2)].real();	\
  R[MS(3,2)].imag() = R[MS(3,2)].imag() + xi[MS(0,2)].imag();	

#define COLLECT_PLUS_Z_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
									\
  R[MS(2,0)].real() = R[MS(2,0)].real() - xi[MS(0,0)].imag();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() + xi[MS(0,0)].real();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() - xi[MS(0,1)].imag();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() + xi[MS(0,1)].real();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() - xi[MS(0,2)].imag();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() + xi[MS(0,2)].real();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() + xi[MS(1,0)].imag();  \
  R[MS(3,0)].imag() = R[MS(3,0)].imag() - xi[MS(1,0)].real();  \
  R[MS(3,1)].real() = R[MS(3,1)].real() + xi[MS(1,1)].imag();  \
  R[MS(3,1)].imag() = R[MS(3,1)].imag() - xi[MS(1,1)].real();  \
  R[MS(3,2)].real() = R[MS(3,2)].real() + xi[MS(1,2)].imag();  \
  R[MS(3,2)].imag() = R[MS(3,2)].imag() - xi[MS(1,2)].real();  


#define COLLECT_MINUS_Z_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
									\
  R[MS(2,0)].real() = R[MS(2,0)].real() + xi[MS(0,0)].imag();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() - xi[MS(0,0)].real();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() + xi[MS(0,1)].imag();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() - xi[MS(0,1)].real();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() + xi[MS(0,2)].imag();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() - xi[MS(0,2)].real();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() - xi[MS(1,0)].imag();  \
  R[MS(3,0)].imag() = R[MS(3,0)].imag() + xi[MS(1,0)].real();  \
  R[MS(3,1)].real() = R[MS(3,1)].real() - xi[MS(1,1)].imag();  \
  R[MS(3,1)].imag() = R[MS(3,1)].imag() + xi[MS(1,1)].real();  \
  R[MS(3,2)].real() = R[MS(3,2)].real() - xi[MS(1,2)].imag();  \
  R[MS(3,2)].imag() = R[MS(3,2)].imag() + xi[MS(1,2)].real();  
  

#define COLLECT_PLUS_T_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	\
							\
  R[MS(2,0)].real() = R[MS(2,0)].real() - xi[MS(0,0)].real();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() - xi[MS(0,0)].imag();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() - xi[MS(0,1)].real();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() - xi[MS(0,1)].imag();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() - xi[MS(0,2)].real();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() - xi[MS(0,2)].imag();	\
									\
  R[MS(3,0)].real() = R[MS(3,0)].real() - xi[MS(1,0)].real();	\
  R[MS(3,0)].imag() = R[MS(3,0)].imag() - xi[MS(1,0)].imag();	\
  R[MS(3,1)].real() = R[MS(3,1)].real() - xi[MS(1,1)].real();	\
  R[MS(3,1)].imag() = R[MS(3,1)].imag() - xi[MS(1,1)].imag();	\
  R[MS(3,2)].real() = R[MS(3,2)].real() - xi[MS(1,2)].real();	\
  R[MS(3,2)].imag() = R[MS(3,2)].imag() - xi[MS(1,2)].imag();	

#define COLLECT_MINUS_T_CHIRAL(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
  						\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];		\
  								\
  R[MS(2,0)].real() = R[MS(2,0)].real() + xi[MS(0,0)].real();	\
  R[MS(2,0)].imag() = R[MS(2,0)].imag() + xi[MS(0,0)].imag();	\
  R[MS(2,1)].real() = R[MS(2,1)].real() + xi[MS(0,1)].real();	\
  R[MS(2,1)].imag() = R[MS(2,1)].imag() + xi[MS(0,1)].imag();	\
  R[MS(2,2)].real() = R[MS(2,2)].real() + xi[MS(0,2)].real();	\
  R[MS(2,2)].imag() = R[MS(2,2)].imag() + xi[MS(0,2)].imag();		\
  									\
  R[MS(3,0)].real() = R[MS(3,0)].real() + xi[MS(1,0)].real();	\
  R[MS(3,0)].imag() = R[MS(3,0)].imag() + xi[MS(1,0)].imag();	\
  R[MS(3,1)].real() = R[MS(3,1)].real() + xi[MS(1,1)].real();	\
  R[MS(3,1)].imag() = R[MS(3,1)].imag() + xi[MS(1,1)].imag();	\
  R[MS(3,2)].real() = R[MS(3,2)].real() + xi[MS(1,2)].real();	\
  R[MS(3,2)].imag() = R[MS(3,2)].imag() + xi[MS(1,2)].imag();	

