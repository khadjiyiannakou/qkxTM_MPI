#define ZERO(R)              \
  R[MS(0,0)].real() = 0.; \
  R[MS(0,0)].imag() = 0.; \
  R[MS(0,1)].real() = 0.; \
  R[MS(0,1)].imag() = 0.; \
  R[MS(0,2)].real() = 0.; \
  R[MS(0,2)].imag() = 0.; \
                             \
  R[MS(1,0)].real() = 0.; \
  R[MS(1,0)].imag() = 0.; \
  R[MS(1,1)].real() = 0.; \
  R[MS(1,1)].imag() = 0.; \
  R[MS(1,2)].real() = 0.; \
  R[MS(1,2)].imag() = 0.; \
                             \
  R[MS(2,0)].real() = 0.; \
  R[MS(2,0)].imag() = 0.; \
  R[MS(2,1)].real() = 0.; \
  R[MS(2,1)].imag() = 0.; \
  R[MS(2,2)].real() = 0.; \
  R[MS(2,2)].imag() = 0.; \
                             \
  R[MS(3,0)].real() = 0.; \
  R[MS(3,0)].imag() = 0.; \
  R[MS(3,1)].real() = 0.; \
  R[MS(3,1)].imag() = 0.; \
  R[MS(3,2)].real() = 0.; \
  R[MS(3,2)].imag() = 0.; 

#define PROJ_MINUS_X(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointPlus[0][iv]][0][0].real() + psi.M[pointPlus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[0][iv]][0][0].imag() - psi.M[pointPlus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[0][iv]][0][1].real() + psi.M[pointPlus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[0][iv]][0][1].imag() - psi.M[pointPlus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[0][iv]][0][2].real() + psi.M[pointPlus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[0][iv]][0][2].imag() - psi.M[pointPlus[0][iv]][3][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[0][iv]][1][0].real() + psi.M[pointPlus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[0][iv]][1][0].imag() - psi.M[pointPlus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[0][iv]][1][1].real() + psi.M[pointPlus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[0][iv]][1][1].imag() - psi.M[pointPlus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[0][iv]][1][2].real() + psi.M[pointPlus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[0][iv]][1][2].imag() - psi.M[pointPlus[0][iv]][2][2].real(); 


#define PROJ_PLUS_X(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointMinus[0][iv]][0][0].real() - psi.M[pointMinus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[0][iv]][0][0].imag() + psi.M[pointMinus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[0][iv]][0][1].real() - psi.M[pointMinus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[0][iv]][0][1].imag() + psi.M[pointMinus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[0][iv]][0][2].real() - psi.M[pointMinus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[0][iv]][0][2].imag() + psi.M[pointMinus[0][iv]][3][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[0][iv]][1][0].real() - psi.M[pointMinus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[0][iv]][1][0].imag() + psi.M[pointMinus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[0][iv]][1][1].real() - psi.M[pointMinus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[0][iv]][1][1].imag() + psi.M[pointMinus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[0][iv]][1][2].real() - psi.M[pointMinus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[0][iv]][1][2].imag() + psi.M[pointMinus[0][iv]][2][2].real(); 



#define PROJ_MINUS_Y(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointPlus[1][iv]][0][0].real() - psi.M[pointPlus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[1][iv]][0][0].imag() - psi.M[pointPlus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[1][iv]][0][1].real() - psi.M[pointPlus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[1][iv]][0][1].imag() - psi.M[pointPlus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[1][iv]][0][2].real() - psi.M[pointPlus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[1][iv]][0][2].imag() - psi.M[pointPlus[1][iv]][3][2].imag(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[1][iv]][1][0].real() + psi.M[pointPlus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[1][iv]][1][0].imag() + psi.M[pointPlus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[1][iv]][1][1].real() + psi.M[pointPlus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[1][iv]][1][1].imag() + psi.M[pointPlus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[1][iv]][1][2].real() + psi.M[pointPlus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[1][iv]][1][2].imag() + psi.M[pointPlus[1][iv]][2][2].imag(); 


#define PROJ_PLUS_Y(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointMinus[1][iv]][0][0].real() + psi.M[pointMinus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[1][iv]][0][0].imag() + psi.M[pointMinus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[1][iv]][0][1].real() + psi.M[pointMinus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[1][iv]][0][1].imag() + psi.M[pointMinus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[1][iv]][0][2].real() + psi.M[pointMinus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[1][iv]][0][2].imag() + psi.M[pointMinus[1][iv]][3][2].imag(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[1][iv]][1][0].real() - psi.M[pointMinus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[1][iv]][1][0].imag() - psi.M[pointMinus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[1][iv]][1][1].real() - psi.M[pointMinus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[1][iv]][1][1].imag() - psi.M[pointMinus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[1][iv]][1][2].real() - psi.M[pointMinus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[1][iv]][1][2].imag() - psi.M[pointMinus[1][iv]][2][2].imag(); 


#define PROJ_MINUS_Z(phi,psi)		      \
  phi[MS(0,0)].real() = psi.M[pointPlus[2][iv]][0][0].real() + psi.M[pointPlus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[2][iv]][0][0].imag() - psi.M[pointPlus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[2][iv]][0][1].real() + psi.M[pointPlus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[2][iv]][0][1].imag() - psi.M[pointPlus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[2][iv]][0][2].real() + psi.M[pointPlus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[2][iv]][0][2].imag() - psi.M[pointPlus[2][iv]][2][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[2][iv]][1][0].real() - psi.M[pointPlus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[2][iv]][1][0].imag() + psi.M[pointPlus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[2][iv]][1][1].real() - psi.M[pointPlus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[2][iv]][1][1].imag() + psi.M[pointPlus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[2][iv]][1][2].real() - psi.M[pointPlus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[2][iv]][1][2].imag() + psi.M[pointPlus[2][iv]][3][2].real(); 


#define PROJ_PLUS_Z(phi,psi)		      \
  phi[MS(0,0)].real() = psi.M[pointMinus[2][iv]][0][0].real() - psi.M[pointMinus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[2][iv]][0][0].imag() + psi.M[pointMinus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[2][iv]][0][1].real() - psi.M[pointMinus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[2][iv]][0][1].imag() + psi.M[pointMinus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[2][iv]][0][2].real() - psi.M[pointMinus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[2][iv]][0][2].imag() + psi.M[pointMinus[2][iv]][2][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[2][iv]][1][0].real() + psi.M[pointMinus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[2][iv]][1][0].imag() - psi.M[pointMinus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[2][iv]][1][1].real() + psi.M[pointMinus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[2][iv]][1][1].imag() - psi.M[pointMinus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[2][iv]][1][2].real() + psi.M[pointMinus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[2][iv]][1][2].imag() - psi.M[pointMinus[2][iv]][3][2].real(); 


#define PROJ_MINUS_T(phi,psi)			\
  phi[MS(0,0)].real() = 2*psi.M[pointPlus[3][iv]][2][0].real();	\
  phi[MS(0,0)].imag() = 2*psi.M[pointPlus[3][iv]][2][0].imag();	\
  phi[MS(0,1)].real() = 2*psi.M[pointPlus[3][iv]][2][1].real();	\
  phi[MS(0,1)].imag() = 2*psi.M[pointPlus[3][iv]][2][1].imag();	\
  phi[MS(0,2)].real() = 2*psi.M[pointPlus[3][iv]][2][2].real();	\
  phi[MS(0,2)].imag() = 2*psi.M[pointPlus[3][iv]][2][2].imag();	\
								\
  phi[MS(1,0)].real() = 2*psi.M[pointPlus[3][iv]][3][0].real();	\
  phi[MS(1,0)].imag() = 2*psi.M[pointPlus[3][iv]][3][0].imag();	\
  phi[MS(1,1)].real() = 2*psi.M[pointPlus[3][iv]][3][1].real();	\
  phi[MS(1,1)].imag() = 2*psi.M[pointPlus[3][iv]][3][1].imag();	\
  phi[MS(1,2)].real() = 2*psi.M[pointPlus[3][iv]][3][2].real();	\
  phi[MS(1,2)].imag() = 2*psi.M[pointPlus[3][iv]][3][2].imag();	


#define PROJ_PLUS_T(phi,psi)			\
  phi[MS(0,0)].real() = 2*psi.M[pointMinus[3][iv]][0][0].real();	\
  phi[MS(0,0)].imag() = 2*psi.M[pointMinus[3][iv]][0][0].imag();	\
  phi[MS(0,1)].real() = 2*psi.M[pointMinus[3][iv]][0][1].real();	\
  phi[MS(0,1)].imag() = 2*psi.M[pointMinus[3][iv]][0][1].imag();	\
  phi[MS(0,2)].real() = 2*psi.M[pointMinus[3][iv]][0][2].real();	\
  phi[MS(0,2)].imag() = 2*psi.M[pointMinus[3][iv]][0][2].imag();	\
								\
  phi[MS(1,0)].real() = 2*psi.M[pointMinus[3][iv]][1][0].real();	\
  phi[MS(1,0)].imag() = 2*psi.M[pointMinus[3][iv]][1][0].imag();	\
  phi[MS(1,1)].real() = 2*psi.M[pointMinus[3][iv]][1][1].real();	\
  phi[MS(1,1)].imag() = 2*psi.M[pointMinus[3][iv]][1][1].imag();	\
  phi[MS(1,2)].real() = 2*psi.M[pointMinus[3][iv]][1][2].real();	\
  phi[MS(1,2)].imag() = 2*psi.M[pointMinus[3][iv]][1][2].imag();	



//////////////////////////////////////////////////////////////////////

#define APPLY_LINK(xi,u,phi,mu)			\
  xi[MS(0,0)] = u.M[iv][mu][0][0] * phi[MS(0,0)] + u.M[iv][mu][0][1] * phi[MS(0,1)] + u.M[iv][mu][0][2] * phi[MS(0,2)]; \
  xi[MS(0,1)] = u.M[iv][mu][1][0] * phi[MS(0,0)] + u.M[iv][mu][1][1] * phi[MS(0,1)] + u.M[iv][mu][1][2] * phi[MS(0,2)]; \
  xi[MS(0,2)] = u.M[iv][mu][2][0] * phi[MS(0,0)] + u.M[iv][mu][2][1] * phi[MS(0,1)] + u.M[iv][mu][2][2] * phi[MS(0,2)]; \
									\
  xi[MS(1,0)] = u.M[iv][mu][0][0] * phi[MS(1,0)] + u.M[iv][mu][0][1] * phi[MS(1,1)] + u.M[iv][mu][0][2] * phi[MS(1,2)]; \
  xi[MS(1,1)] = u.M[iv][mu][1][0] * phi[MS(1,0)] + u.M[iv][mu][1][1] * phi[MS(1,1)] + u.M[iv][mu][1][2] * phi[MS(1,2)]; \
  xi[MS(1,2)] = u.M[iv][mu][2][0] * phi[MS(1,0)] + u.M[iv][mu][2][1] * phi[MS(1,1)] + u.M[iv][mu][2][2] * phi[MS(1,2)]; 

#define APPLY_LINK_DAG(xi,u,phi,mu)			\
  xi[MS(0,0)] = conj(u.M[pointMinus[mu][iv]][mu][0][0]) * phi[MS(0,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][0]) * phi[MS(0,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][0]) * phi[MS(0,2)]; \
  xi[MS(0,1)] = conj(u.M[pointMinus[mu][iv]][mu][0][1]) * phi[MS(0,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][1]) * phi[MS(0,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][1]) * phi[MS(0,2)]; \
  xi[MS(0,2)] = conj(u.M[pointMinus[mu][iv]][mu][0][2]) * phi[MS(0,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][2]) * phi[MS(0,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][2]) * phi[MS(0,2)]; \
									\
  xi[MS(1,0)] = conj(u.M[pointMinus[mu][iv]][mu][0][0]) * phi[MS(1,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][0]) * phi[MS(1,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][0]) * phi[MS(1,2)]; \
  xi[MS(1,1)] = conj(u.M[pointMinus[mu][iv]][mu][0][1]) * phi[MS(1,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][1]) * phi[MS(1,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][1]) * phi[MS(1,2)]; \
  xi[MS(1,2)] = conj(u.M[pointMinus[mu][iv]][mu][0][2]) * phi[MS(1,0)] + conj(u.M[pointMinus[mu][iv]][mu][1][2]) * phi[MS(1,1)] + conj(u.M[pointMinus[mu][iv]][mu][2][2]) * phi[MS(1,2)]; 



////////////////////////////////////////////////////////////////////////////////

#define COLLECT_MINUS_X(R,xi)			\
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

#define COLLECT_PLUS_X(R,xi)			\
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


#define COLLECT_MINUS_Y(R,xi)			\
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


#define COLLECT_PLUS_Y(R,xi)			\
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

#define COLLECT_MINUS_Z(R,xi)			\
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


#define COLLECT_PLUS_Z(R,xi)			\
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
  

#define COLLECT_MINUS_T(R,xi)			\
  R[MS(2,0)] = R[MS(2,0)] + xi[MS(0,0)];	\
  R[MS(2,1)] = R[MS(2,1)] + xi[MS(0,1)];	\
  R[MS(2,2)] = R[MS(2,2)] + xi[MS(0,2)];	\
							\
  R[MS(3,0)] = R[MS(3,0)] + xi[MS(1,0)];	\
  R[MS(3,1)] = R[MS(3,1)] + xi[MS(1,1)];	\
  R[MS(3,2)] = R[MS(3,2)] + xi[MS(1,2)];	

#define COLLECT_PLUS_T(R,xi)			\
  R[MS(0,0)] = R[MS(0,0)] + xi[MS(0,0)];	\
  R[MS(0,1)] = R[MS(0,1)] + xi[MS(0,1)];	\
  R[MS(0,2)] = R[MS(0,2)] + xi[MS(0,2)];	\
							\
  R[MS(1,0)] = R[MS(1,0)] + xi[MS(1,0)];	\
  R[MS(1,1)] = R[MS(1,1)] + xi[MS(1,1)];	\
  R[MS(1,2)] = R[MS(1,2)] + xi[MS(1,2)];	



//////////////////////////////////////////// for dagger /////////////////////////////////////////

#define DAG_PROJ_MINUS_X(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointMinus[0][iv]][0][0].real() + psi.M[pointMinus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[0][iv]][0][0].imag() - psi.M[pointMinus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[0][iv]][0][1].real() + psi.M[pointMinus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[0][iv]][0][1].imag() - psi.M[pointMinus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[0][iv]][0][2].real() + psi.M[pointMinus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[0][iv]][0][2].imag() - psi.M[pointMinus[0][iv]][3][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[0][iv]][1][0].real() + psi.M[pointMinus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[0][iv]][1][0].imag() - psi.M[pointMinus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[0][iv]][1][1].real() + psi.M[pointMinus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[0][iv]][1][1].imag() - psi.M[pointMinus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[0][iv]][1][2].real() + psi.M[pointMinus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[0][iv]][1][2].imag() - psi.M[pointMinus[0][iv]][2][2].real(); 


#define DAG_PROJ_PLUS_X(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointPlus[0][iv]][0][0].real() - psi.M[pointPlus[0][iv]][3][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[0][iv]][0][0].imag() + psi.M[pointPlus[0][iv]][3][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[0][iv]][0][1].real() - psi.M[pointPlus[0][iv]][3][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[0][iv]][0][1].imag() + psi.M[pointPlus[0][iv]][3][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[0][iv]][0][2].real() - psi.M[pointPlus[0][iv]][3][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[0][iv]][0][2].imag() + psi.M[pointPlus[0][iv]][3][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[0][iv]][1][0].real() - psi.M[pointPlus[0][iv]][2][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[0][iv]][1][0].imag() + psi.M[pointPlus[0][iv]][2][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[0][iv]][1][1].real() - psi.M[pointPlus[0][iv]][2][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[0][iv]][1][1].imag() + psi.M[pointPlus[0][iv]][2][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[0][iv]][1][2].real() - psi.M[pointPlus[0][iv]][2][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[0][iv]][1][2].imag() + psi.M[pointPlus[0][iv]][2][2].real(); 


#define DAG_PROJ_MINUS_Y(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointMinus[1][iv]][0][0].real() - psi.M[pointMinus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[1][iv]][0][0].imag() - psi.M[pointMinus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[1][iv]][0][1].real() - psi.M[pointMinus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[1][iv]][0][1].imag() - psi.M[pointMinus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[1][iv]][0][2].real() - psi.M[pointMinus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[1][iv]][0][2].imag() - psi.M[pointMinus[1][iv]][3][2].imag(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[1][iv]][1][0].real() + psi.M[pointMinus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[1][iv]][1][0].imag() + psi.M[pointMinus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[1][iv]][1][1].real() + psi.M[pointMinus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[1][iv]][1][1].imag() + psi.M[pointMinus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[1][iv]][1][2].real() + psi.M[pointMinus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[1][iv]][1][2].imag() + psi.M[pointMinus[1][iv]][2][2].imag(); 


#define DAG_PROJ_PLUS_Y(phi,psi)						\
  phi[MS(0,0)].real() = psi.M[pointPlus[1][iv]][0][0].real() + psi.M[pointPlus[1][iv]][3][0].real(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[1][iv]][0][0].imag() + psi.M[pointPlus[1][iv]][3][0].imag(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[1][iv]][0][1].real() + psi.M[pointPlus[1][iv]][3][1].real(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[1][iv]][0][1].imag() + psi.M[pointPlus[1][iv]][3][1].imag(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[1][iv]][0][2].real() + psi.M[pointPlus[1][iv]][3][2].real(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[1][iv]][0][2].imag() + psi.M[pointPlus[1][iv]][3][2].imag(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[1][iv]][1][0].real() - psi.M[pointPlus[1][iv]][2][0].real(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[1][iv]][1][0].imag() - psi.M[pointPlus[1][iv]][2][0].imag(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[1][iv]][1][1].real() - psi.M[pointPlus[1][iv]][2][1].real(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[1][iv]][1][1].imag() - psi.M[pointPlus[1][iv]][2][1].imag(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[1][iv]][1][2].real() - psi.M[pointPlus[1][iv]][2][2].real(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[1][iv]][1][2].imag() - psi.M[pointPlus[1][iv]][2][2].imag(); 


#define DAG_PROJ_MINUS_Z(phi,psi)		      \
  phi[MS(0,0)].real() = psi.M[pointMinus[2][iv]][0][0].real() + psi.M[pointMinus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointMinus[2][iv]][0][0].imag() - psi.M[pointMinus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointMinus[2][iv]][0][1].real() + psi.M[pointMinus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointMinus[2][iv]][0][1].imag() - psi.M[pointMinus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointMinus[2][iv]][0][2].real() + psi.M[pointMinus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointMinus[2][iv]][0][2].imag() - psi.M[pointMinus[2][iv]][2][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointMinus[2][iv]][1][0].real() - psi.M[pointMinus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointMinus[2][iv]][1][0].imag() + psi.M[pointMinus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointMinus[2][iv]][1][1].real() - psi.M[pointMinus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointMinus[2][iv]][1][1].imag() + psi.M[pointMinus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointMinus[2][iv]][1][2].real() - psi.M[pointMinus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointMinus[2][iv]][1][2].imag() + psi.M[pointMinus[2][iv]][3][2].real(); 


#define DAG_PROJ_PLUS_Z(phi,psi)		      \
  phi[MS(0,0)].real() = psi.M[pointPlus[2][iv]][0][0].real() - psi.M[pointPlus[2][iv]][2][0].imag(); \
  phi[MS(0,0)].imag() = psi.M[pointPlus[2][iv]][0][0].imag() + psi.M[pointPlus[2][iv]][2][0].real(); \
  phi[MS(0,1)].real() = psi.M[pointPlus[2][iv]][0][1].real() - psi.M[pointPlus[2][iv]][2][1].imag(); \
  phi[MS(0,1)].imag() = psi.M[pointPlus[2][iv]][0][1].imag() + psi.M[pointPlus[2][iv]][2][1].real(); \
  phi[MS(0,2)].real() = psi.M[pointPlus[2][iv]][0][2].real() - psi.M[pointPlus[2][iv]][2][2].imag(); \
  phi[MS(0,2)].imag() = psi.M[pointPlus[2][iv]][0][2].imag() + psi.M[pointPlus[2][iv]][2][2].real(); \
									\
  phi[MS(1,0)].real() = psi.M[pointPlus[2][iv]][1][0].real() + psi.M[pointPlus[2][iv]][3][0].imag(); \
  phi[MS(1,0)].imag() = psi.M[pointPlus[2][iv]][1][0].imag() - psi.M[pointPlus[2][iv]][3][0].real(); \
  phi[MS(1,1)].real() = psi.M[pointPlus[2][iv]][1][1].real() + psi.M[pointPlus[2][iv]][3][1].imag(); \
  phi[MS(1,1)].imag() = psi.M[pointPlus[2][iv]][1][1].imag() - psi.M[pointPlus[2][iv]][3][1].real(); \
  phi[MS(1,2)].real() = psi.M[pointPlus[2][iv]][1][2].real() + psi.M[pointPlus[2][iv]][3][2].imag(); \
  phi[MS(1,2)].imag() = psi.M[pointPlus[2][iv]][1][2].imag() - psi.M[pointPlus[2][iv]][3][2].real(); 


#define DAG_PROJ_MINUS_T(phi,psi)			\
  phi[MS(0,0)].real() = 2*psi.M[pointMinus[3][iv]][2][0].real();	\
  phi[MS(0,0)].imag() = 2*psi.M[pointMinus[3][iv]][2][0].imag();	\
  phi[MS(0,1)].real() = 2*psi.M[pointMinus[3][iv]][2][1].real();	\
  phi[MS(0,1)].imag() = 2*psi.M[pointMinus[3][iv]][2][1].imag();	\
  phi[MS(0,2)].real() = 2*psi.M[pointMinus[3][iv]][2][2].real();	\
  phi[MS(0,2)].imag() = 2*psi.M[pointMinus[3][iv]][2][2].imag();	\
								\
  phi[MS(1,0)].real() = 2*psi.M[pointMinus[3][iv]][3][0].real();	\
  phi[MS(1,0)].imag() = 2*psi.M[pointMinus[3][iv]][3][0].imag();	\
  phi[MS(1,1)].real() = 2*psi.M[pointMinus[3][iv]][3][1].real();	\
  phi[MS(1,1)].imag() = 2*psi.M[pointMinus[3][iv]][3][1].imag();	\
  phi[MS(1,2)].real() = 2*psi.M[pointMinus[3][iv]][3][2].real();	\
  phi[MS(1,2)].imag() = 2*psi.M[pointMinus[3][iv]][3][2].imag();	


#define DAG_PROJ_PLUS_T(phi,psi)			\
  phi[MS(0,0)].real() = 2*psi.M[pointPlus[3][iv]][0][0].real();	\
  phi[MS(0,0)].imag() = 2*psi.M[pointPlus[3][iv]][0][0].imag();	\
  phi[MS(0,1)].real() = 2*psi.M[pointPlus[3][iv]][0][1].real();	\
  phi[MS(0,1)].imag() = 2*psi.M[pointPlus[3][iv]][0][1].imag();	\
  phi[MS(0,2)].real() = 2*psi.M[pointPlus[3][iv]][0][2].real();	\
  phi[MS(0,2)].imag() = 2*psi.M[pointPlus[3][iv]][0][2].imag();	\
								\
  phi[MS(1,0)].real() = 2*psi.M[pointPlus[3][iv]][1][0].real();	\
  phi[MS(1,0)].imag() = 2*psi.M[pointPlus[3][iv]][1][0].imag();	\
  phi[MS(1,1)].real() = 2*psi.M[pointPlus[3][iv]][1][1].real();	\
  phi[MS(1,1)].imag() = 2*psi.M[pointPlus[3][iv]][1][1].imag();	\
  phi[MS(1,2)].real() = 2*psi.M[pointPlus[3][iv]][1][2].real();	\
  phi[MS(1,2)].imag() = 2*psi.M[pointPlus[3][iv]][1][2].imag();	
