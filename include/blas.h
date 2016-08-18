// macros from qcd package be carefull


#define qcd_MUL3x3(a,b,c) ({ a[0][0].real() =  b[0][0].real()*c[0][0].real() - b[0][0].imag()*c[0][0].imag() \
                                          +b[0][1].real()*c[1][0].real() - b[0][1].imag()*c[1][0].imag() \
                                          +b[0][2].real()*c[2][0].real() - b[0][2].imag()*c[2][0].imag();\
                             a[0][0].imag() =  b[0][0].real()*c[0][0].imag() + b[0][0].imag()*c[0][0].real() \
                                          +b[0][1].real()*c[1][0].imag() + b[0][1].imag()*c[1][0].real() \
                                          +b[0][2].real()*c[2][0].imag() + b[0][2].imag()*c[2][0].real();\
                             a[0][1].real() =  b[0][0].real()*c[0][1].real() - b[0][0].imag()*c[0][1].imag() \
                                          +b[0][1].real()*c[1][1].real() - b[0][1].imag()*c[1][1].imag() \
                                          +b[0][2].real()*c[2][1].real() - b[0][2].imag()*c[2][1].imag();\
                             a[0][1].imag() =  b[0][0].real()*c[0][1].imag() + b[0][0].imag()*c[0][1].real() \
                                          +b[0][1].real()*c[1][1].imag() + b[0][1].imag()*c[1][1].real() \
                                          +b[0][2].real()*c[2][1].imag() + b[0][2].imag()*c[2][1].real();\
                             a[0][2].real() =  b[0][0].real()*c[0][2].real() - b[0][0].imag()*c[0][2].imag() \
                                          +b[0][1].real()*c[1][2].real() - b[0][1].imag()*c[1][2].imag() \
                                          +b[0][2].real()*c[2][2].real() - b[0][2].imag()*c[2][2].imag();\
                             a[0][2].imag() =  b[0][0].real()*c[0][2].imag() + b[0][0].imag()*c[0][2].real() \
                                          +b[0][1].real()*c[1][2].imag() + b[0][1].imag()*c[1][2].real() \
                                          +b[0][2].real()*c[2][2].imag() + b[0][2].imag()*c[2][2].real();\
                                                                                         \
                             a[1][0].real() =  b[1][0].real()*c[0][0].real() - b[1][0].imag()*c[0][0].imag() \
                                          +b[1][1].real()*c[1][0].real() - b[1][1].imag()*c[1][0].imag() \
                                          +b[1][2].real()*c[2][0].real() - b[1][2].imag()*c[2][0].imag();\
                             a[1][0].imag() =  b[1][0].real()*c[0][0].imag() + b[1][0].imag()*c[0][0].real() \
                                          +b[1][1].real()*c[1][0].imag() + b[1][1].imag()*c[1][0].real() \
                                          +b[1][2].real()*c[2][0].imag() + b[1][2].imag()*c[2][0].real();\
                             a[1][1].real() =  b[1][0].real()*c[0][1].real() - b[1][0].imag()*c[0][1].imag() \
                                          +b[1][1].real()*c[1][1].real() - b[1][1].imag()*c[1][1].imag() \
                                          +b[1][2].real()*c[2][1].real() - b[1][2].imag()*c[2][1].imag();\
                             a[1][1].imag() =  b[1][0].real()*c[0][1].imag() + b[1][0].imag()*c[0][1].real() \
                                          +b[1][1].real()*c[1][1].imag() + b[1][1].imag()*c[1][1].real() \
                                          +b[1][2].real()*c[2][1].imag() + b[1][2].imag()*c[2][1].real();\
                             a[1][2].real() =  b[1][0].real()*c[0][2].real() - b[1][0].imag()*c[0][2].imag() \
                                          +b[1][1].real()*c[1][2].real() - b[1][1].imag()*c[1][2].imag() \
                                          +b[1][2].real()*c[2][2].real() - b[1][2].imag()*c[2][2].imag();\
                             a[1][2].imag() =  b[1][0].real()*c[0][2].imag() + b[1][0].imag()*c[0][2].real() \
                                          +b[1][1].real()*c[1][2].imag() + b[1][1].imag()*c[1][2].real() \
                                          +b[1][2].real()*c[2][2].imag() + b[1][2].imag()*c[2][2].real();\
                                                                                         \
                             a[2][0].real() =  b[2][0].real()*c[0][0].real() - b[2][0].imag()*c[0][0].imag() \
                                          +b[2][1].real()*c[1][0].real() - b[2][1].imag()*c[1][0].imag() \
                                          +b[2][2].real()*c[2][0].real() - b[2][2].imag()*c[2][0].imag();\
                             a[2][0].imag() =  b[2][0].real()*c[0][0].imag() + b[2][0].imag()*c[0][0].real() \
                                          +b[2][1].real()*c[1][0].imag() + b[2][1].imag()*c[1][0].real() \
                                          +b[2][2].real()*c[2][0].imag() + b[2][2].imag()*c[2][0].real();\
                             a[2][1].real() =  b[2][0].real()*c[0][1].real() - b[2][0].imag()*c[0][1].imag() \
                                          +b[2][1].real()*c[1][1].real() - b[2][1].imag()*c[1][1].imag() \
                                          +b[2][2].real()*c[2][1].real() - b[2][2].imag()*c[2][1].imag();\
                             a[2][1].imag() =  b[2][0].real()*c[0][1].imag() + b[2][0].imag()*c[0][1].real() \
                                          +b[2][1].real()*c[1][1].imag() + b[2][1].imag()*c[1][1].real() \
                                          +b[2][2].real()*c[2][1].imag() + b[2][2].imag()*c[2][1].real();\
                             a[2][2].real() =  b[2][0].real()*c[0][2].real() - b[2][0].imag()*c[0][2].imag() \
                                          +b[2][1].real()*c[1][2].real() - b[2][1].imag()*c[1][2].imag() \
                                          +b[2][2].real()*c[2][2].real() - b[2][2].imag()*c[2][2].imag();\
                             a[2][2].imag() =  b[2][0].real()*c[0][2].imag() + b[2][0].imag()*c[0][2].real() \
                                          +b[2][1].real()*c[1][2].imag() + b[2][1].imag()*c[1][2].real() \
                                          +b[2][2].real()*c[2][2].imag() + b[2][2].imag()*c[2][2].real();\
})
 

/* a = b* adjoint(c), whereal() a,b,c areal() 3x3 matrices */
#define qcd_MULADJOINT3x3(a,b,c) ({\
                             a[0][0].real() =  b[0][0].real()*c[0][0].real() + b[0][0].imag()*c[0][0].imag() \
                                          +b[0][1].real()*c[0][1].real() + b[0][1].imag()*c[0][1].imag() \
                                          +b[0][2].real()*c[0][2].real() + b[0][2].imag()*c[0][2].imag();\
                             a[0][0].imag() = -b[0][0].real()*c[0][0].imag() + b[0][0].imag()*c[0][0].real() \
                                          -b[0][1].real()*c[0][1].imag() + b[0][1].imag()*c[0][1].real() \
                                          -b[0][2].real()*c[0][2].imag() + b[0][2].imag()*c[0][2].real();\
                             a[0][1].real() =  b[0][0].real()*c[1][0].real() + b[0][0].imag()*c[1][0].imag() \
                                          +b[0][1].real()*c[1][1].real() + b[0][1].imag()*c[1][1].imag() \
                                          +b[0][2].real()*c[1][2].real() + b[0][2].imag()*c[1][2].imag();\
                             a[0][1].imag() = -b[0][0].real()*c[1][0].imag() + b[0][0].imag()*c[1][0].real() \
                                          -b[0][1].real()*c[1][1].imag() + b[0][1].imag()*c[1][1].real() \
                                          -b[0][2].real()*c[1][2].imag() + b[0][2].imag()*c[1][2].real();\
                             a[0][2].real() =  b[0][0].real()*c[2][0].real() + b[0][0].imag()*c[2][0].imag() \
                                          +b[0][1].real()*c[2][1].real() + b[0][1].imag()*c[2][1].imag() \
                                          +b[0][2].real()*c[2][2].real() + b[0][2].imag()*c[2][2].imag();\
                             a[0][2].imag() = -b[0][0].real()*c[2][0].imag() + b[0][0].imag()*c[2][0].real() \
                                          -b[0][1].real()*c[2][1].imag() + b[0][1].imag()*c[2][1].real() \
                                          -b[0][2].real()*c[2][2].imag() + b[0][2].imag()*c[2][2].real();\
                                                                                         \
                             a[1][0].real() =  b[1][0].real()*c[0][0].real() + b[1][0].imag()*c[0][0].imag() \
                                          +b[1][1].real()*c[0][1].real() + b[1][1].imag()*c[0][1].imag() \
                                          +b[1][2].real()*c[0][2].real() + b[1][2].imag()*c[0][2].imag();\
                             a[1][0].imag() = -b[1][0].real()*c[0][0].imag() + b[1][0].imag()*c[0][0].real() \
                                          -b[1][1].real()*c[0][1].imag() + b[1][1].imag()*c[0][1].real() \
                                          -b[1][2].real()*c[0][2].imag() + b[1][2].imag()*c[0][2].real();\
                             a[1][1].real() =  b[1][0].real()*c[1][0].real() + b[1][0].imag()*c[1][0].imag() \
                                          +b[1][1].real()*c[1][1].real() + b[1][1].imag()*c[1][1].imag() \
                                          +b[1][2].real()*c[1][2].real() + b[1][2].imag()*c[1][2].imag();\
                             a[1][1].imag() = -b[1][0].real()*c[1][0].imag() + b[1][0].imag()*c[1][0].real() \
                                          -b[1][1].real()*c[1][1].imag() + b[1][1].imag()*c[1][1].real() \
                                          -b[1][2].real()*c[1][2].imag() + b[1][2].imag()*c[1][2].real();\
                             a[1][2].real() =  b[1][0].real()*c[2][0].real() + b[1][0].imag()*c[2][0].imag() \
                                          +b[1][1].real()*c[2][1].real() + b[1][1].imag()*c[2][1].imag() \
                                          +b[1][2].real()*c[2][2].real() + b[1][2].imag()*c[2][2].imag();\
                             a[1][2].imag() = -b[1][0].real()*c[2][0].imag() + b[1][0].imag()*c[2][0].real() \
                                          -b[1][1].real()*c[2][1].imag() + b[1][1].imag()*c[2][1].real() \
                                          -b[1][2].real()*c[2][2].imag() + b[1][2].imag()*c[2][2].real();\
                                                                                         \
                             a[2][0].real() =  b[2][0].real()*c[0][0].real() + b[2][0].imag()*c[0][0].imag() \
                                          +b[2][1].real()*c[0][1].real() + b[2][1].imag()*c[0][1].imag() \
                                          +b[2][2].real()*c[0][2].real() + b[2][2].imag()*c[0][2].imag();\
                             a[2][0].imag() = -b[2][0].real()*c[0][0].imag() + b[2][0].imag()*c[0][0].real() \
                                          -b[2][1].real()*c[0][1].imag() + b[2][1].imag()*c[0][1].real() \
                                          -b[2][2].real()*c[0][2].imag() + b[2][2].imag()*c[0][2].real();\
                             a[2][1].real() =  b[2][0].real()*c[1][0].real() + b[2][0].imag()*c[1][0].imag() \
                                          +b[2][1].real()*c[1][1].real() + b[2][1].imag()*c[1][1].imag() \
                                          +b[2][2].real()*c[1][2].real() + b[2][2].imag()*c[1][2].imag();\
                             a[2][1].imag() = -b[2][0].real()*c[1][0].imag() + b[2][0].imag()*c[1][0].real() \
                                          -b[2][1].real()*c[1][1].imag() + b[2][1].imag()*c[1][1].real() \
                                          -b[2][2].real()*c[1][2].imag() + b[2][2].imag()*c[1][2].real();\
                             a[2][2].real() =  b[2][0].real()*c[2][0].real() + b[2][0].imag()*c[2][0].imag() \
                                          +b[2][1].real()*c[2][1].real() + b[2][1].imag()*c[2][1].imag() \
                                          +b[2][2].real()*c[2][2].real() + b[2][2].imag()*c[2][2].imag();\
                             a[2][2].imag() = -b[2][0].real()*c[2][0].imag() + b[2][0].imag()*c[2][0].real() \
                                          -b[2][1].real()*c[2][1].imag() + b[2][1].imag()*c[2][1].real() \
                                          -b[2][2].real()*c[2][2].imag() + b[2][2].imag()*c[2][2].real();\
})


/* a = adjoint(b)*c, whereal() a,b,c areal() 3x3 matrices */ 
#define qcd_ADJOINTMUL3x3(a,b,c) ({\
                             a[0][0].real() =  b[0][0].real()*c[0][0].real() + b[0][0].imag()*c[0][0].imag() \
                                          +b[1][0].real()*c[1][0].real() + b[1][0].imag()*c[1][0].imag() \
                                          +b[2][0].real()*c[2][0].real() + b[2][0].imag()*c[2][0].imag();\
                             a[0][0].imag() =  b[0][0].real()*c[0][0].imag() - b[0][0].imag()*c[0][0].real() \
                                          +b[1][0].real()*c[1][0].imag() - b[1][0].imag()*c[1][0].real() \
                                          +b[2][0].real()*c[2][0].imag() - b[2][0].imag()*c[2][0].real();\
                             a[0][1].real() =  b[0][0].real()*c[0][1].real() + b[0][0].imag()*c[0][1].imag() \
                                          +b[1][0].real()*c[1][1].real() + b[1][0].imag()*c[1][1].imag() \
                                          +b[2][0].real()*c[2][1].real() + b[2][0].imag()*c[2][1].imag();\
                             a[0][1].imag() =  b[0][0].real()*c[0][1].imag() - b[0][0].imag()*c[0][1].real() \
                                          +b[1][0].real()*c[1][1].imag() - b[1][0].imag()*c[1][1].real() \
                                          +b[2][0].real()*c[2][1].imag() - b[2][0].imag()*c[2][1].real();\
                             a[0][2].real() =  b[0][0].real()*c[0][2].real() + b[0][0].imag()*c[0][2].imag() \
                                          +b[1][0].real()*c[1][2].real() + b[1][0].imag()*c[1][2].imag() \
                                          +b[2][0].real()*c[2][2].real() + b[2][0].imag()*c[2][2].imag();\
                             a[0][2].imag() =  b[0][0].real()*c[0][2].imag() - b[0][0].imag()*c[0][2].real() \
                                          +b[1][0].real()*c[1][2].imag() - b[1][0].imag()*c[1][2].real() \
                                          +b[2][0].real()*c[2][2].imag() - b[2][0].imag()*c[2][2].real();\
                                                                                         \
                             a[1][0].real() =  b[0][1].real()*c[0][0].real() + b[0][1].imag()*c[0][0].imag() \
                                          +b[1][1].real()*c[1][0].real() + b[1][1].imag()*c[1][0].imag() \
                                          +b[2][1].real()*c[2][0].real() + b[2][1].imag()*c[2][0].imag();\
                             a[1][0].imag() =  b[0][1].real()*c[0][0].imag() - b[0][1].imag()*c[0][0].real() \
                                          +b[1][1].real()*c[1][0].imag() - b[1][1].imag()*c[1][0].real() \
                                          +b[2][1].real()*c[2][0].imag() - b[2][1].imag()*c[2][0].real();\
                             a[1][1].real() =  b[0][1].real()*c[0][1].real() + b[0][1].imag()*c[0][1].imag() \
                                          +b[1][1].real()*c[1][1].real() + b[1][1].imag()*c[1][1].imag() \
                                          +b[2][1].real()*c[2][1].real() + b[2][1].imag()*c[2][1].imag();\
                             a[1][1].imag() =  b[0][1].real()*c[0][1].imag() - b[0][1].imag()*c[0][1].real() \
                                          +b[1][1].real()*c[1][1].imag() - b[1][1].imag()*c[1][1].real() \
                                          +b[2][1].real()*c[2][1].imag() - b[2][1].imag()*c[2][1].real();\
                             a[1][2].real() =  b[0][1].real()*c[0][2].real() + b[0][1].imag()*c[0][2].imag() \
                                          +b[1][1].real()*c[1][2].real() + b[1][1].imag()*c[1][2].imag() \
                                          +b[2][1].real()*c[2][2].real() + b[2][1].imag()*c[2][2].imag();\
                             a[1][2].imag() =  b[0][1].real()*c[0][2].imag() - b[0][1].imag()*c[0][2].real() \
                                          +b[1][1].real()*c[1][2].imag() - b[1][1].imag()*c[1][2].real() \
                                          +b[2][1].real()*c[2][2].imag() - b[2][1].imag()*c[2][2].real();\
                                                                                         \
                             a[2][0].real() =  b[0][2].real()*c[0][0].real() + b[0][2].imag()*c[0][0].imag() \
                                          +b[1][2].real()*c[1][0].real() + b[1][2].imag()*c[1][0].imag() \
                                          +b[2][2].real()*c[2][0].real() + b[2][2].imag()*c[2][0].imag();\
                             a[2][0].imag() =  b[0][2].real()*c[0][0].imag() - b[0][2].imag()*c[0][0].real() \
                                          +b[1][2].real()*c[1][0].imag() - b[1][2].imag()*c[1][0].real() \
                                          +b[2][2].real()*c[2][0].imag() - b[2][2].imag()*c[2][0].real();\
                             a[2][1].real() =  b[0][2].real()*c[0][1].real() + b[0][2].imag()*c[0][1].imag() \
                                          +b[1][2].real()*c[1][1].real() + b[1][2].imag()*c[1][1].imag() \
                                          +b[2][2].real()*c[2][1].real() + b[2][2].imag()*c[2][1].imag();\
                             a[2][1].imag() =  b[0][2].real()*c[0][1].imag() - b[0][2].imag()*c[0][1].real() \
                                          +b[1][2].real()*c[1][1].imag() - b[1][2].imag()*c[1][1].real() \
                                          +b[2][2].real()*c[2][1].imag() - b[2][2].imag()*c[2][1].real();\
                             a[2][2].real() =  b[0][2].real()*c[0][2].real() + b[0][2].imag()*c[0][2].imag() \
                                          +b[1][2].real()*c[1][2].real() + b[1][2].imag()*c[1][2].imag() \
                                          +b[2][2].real()*c[2][2].real() + b[2][2].imag()*c[2][2].imag();\
                             a[2][2].imag() =  b[0][2].real()*c[0][2].imag() - b[0][2].imag()*c[0][2].real() \
                                          +b[1][2].real()*c[1][2].imag() - b[1][2].imag()*c[1][2].real() \
                                          +b[2][2].real()*c[2][2].imag() - b[2][2].imag()*c[2][2].real();\
})


#define qcd_ACCUM_3x3(x,y) \
  ({						\
    x[0][0] += y[0][0];				\
    x[0][1] += y[0][1];				\
    x[0][2] += y[0][2];				\
    x[1][0] += y[1][0];				\
    x[1][1] += y[1][1];				\
    x[1][2] += y[1][2];				\
    x[2][0] += y[2][0];				\
    x[2][1] += y[2][1];				\
    x[2][2] += y[2][2];				\
  })


#define qcd_SU3TRACER(u) ( (double) (u[0][0].real() + u[1][1].real() + u[2][2].real()))



#define qcd_APPLY_U(u,x,y)  \
   ({ x[0] = u[0] *y[0]  - u[1] *y[1]  + u[2] *y[2]  - u[3] *y[3]  + u[4] *y[4]  - u[5] *y[5];     \
      x[1] = u[0] *y[1]  + u[1] *y[0]  + u[2] *y[3]  + u[3] *y[2]  + u[4] *y[5]  + u[5] *y[4];     \
      x[2] = u[6] *y[0]  - u[7] *y[1]  + u[8] *y[2]  - u[9] *y[3]  + u[10]*y[4]  - u[11]*y[5];     \
      x[3] = u[6] *y[1]  + u[7] *y[0]  + u[8] *y[3]  + u[9] *y[2]  + u[10]*y[5]  + u[11]*y[4];     \
      x[4] = u[12]*y[0]  - u[13]*y[1]  + u[14]*y[2]  - u[15]*y[3]  + u[16]*y[4]  - u[17]*y[5];     \
      x[5] = u[12]*y[1]  + u[13]*y[0]  + u[14]*y[3]  + u[15]*y[2]  + u[16]*y[5]  + u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  - u[1] *y[7]  + u[2] *y[8]  - u[3] *y[9]  + u[4] *y[10] - u[5] *y[11];    \
      x[7] = u[0] *y[7]  + u[1] *y[6]  + u[2] *y[9]  + u[3] *y[8]  + u[4] *y[11] + u[5] *y[10];    \
      x[8] = u[6] *y[6]  - u[7] *y[7]  + u[8] *y[8]  - u[9] *y[9]  + u[10]*y[10] - u[11]*y[11];    \
      x[9] = u[6] *y[7]  + u[7] *y[6]  + u[8] *y[9]  + u[9] *y[8]  + u[10]*y[11] + u[11]*y[10];    \
      x[10]= u[12]*y[6]  - u[13]*y[7]  + u[14]*y[8]  - u[15]*y[9]  + u[16]*y[10] - u[17]*y[11];    \
      x[11]= u[12]*y[7]  + u[13]*y[6]  + u[14]*y[9]  + u[15]*y[8]  + u[16]*y[11] + u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] - u[1] *y[13] + u[2] *y[14] - u[3] *y[15] + u[4] *y[16] - u[5] *y[17];    \
      x[13]= u[0] *y[13] + u[1] *y[12] + u[2] *y[15] + u[3] *y[14] + u[4] *y[17] + u[5] *y[16];    \
      x[14]= u[6] *y[12] - u[7] *y[13] + u[8] *y[14] - u[9] *y[15] + u[10]*y[16] - u[11]*y[17];    \
      x[15]= u[6] *y[13] + u[7] *y[12] + u[8] *y[15] + u[9] *y[14] + u[10]*y[17] + u[11]*y[16];    \
      x[16]= u[12]*y[12] - u[13]*y[13] + u[14]*y[14] - u[15]*y[15] + u[16]*y[16] - u[17]*y[17];    \
      x[17]= u[12]*y[13] + u[13]*y[12] + u[14]*y[15] + u[15]*y[14] + u[16]*y[17] + u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] - u[1] *y[19] + u[2] *y[20] - u[3] *y[21] + u[4] *y[22] - u[5] *y[23];    \
      x[19]= u[0] *y[19] + u[1] *y[18] + u[2] *y[21] + u[3] *y[20] + u[4] *y[23] + u[5] *y[22];    \
      x[20]= u[6] *y[18] - u[7] *y[19] + u[8] *y[20] - u[9] *y[21] + u[10]*y[22] - u[11]*y[23];    \
      x[21]= u[6] *y[19] + u[7] *y[18] + u[8] *y[21] + u[9] *y[20] + u[10]*y[23] + u[11]*y[22];    \
      x[22]= u[12]*y[18] - u[13]*y[19] + u[14]*y[20] - u[15]*y[21] + u[16]*y[22] - u[17]*y[23];    \
      x[23]= u[12]*y[19] + u[13]*y[18] + u[14]*y[21] + u[15]*y[20] + u[16]*y[23] + u[17]*y[22];    \
   })

#define qcd_APPLY_U_DAGGER(u,x,y)  \
   ({ x[0] = u[0] *y[0]  + u[1] *y[1]  + u[6] *y[2]  + u[7] *y[3]  + u[12]*y[4]  + u[13]*y[5];     \
      x[1] = u[0] *y[1]  - u[1] *y[0]  + u[6] *y[3]  - u[7] *y[2]  + u[12]*y[5]  - u[13]*y[4];     \
      x[2] = u[2] *y[0]  + u[3] *y[1]  + u[8] *y[2]  + u[9] *y[3]  + u[14]*y[4]  + u[15]*y[5];     \
      x[3] = u[2] *y[1]  - u[3] *y[0]  + u[8] *y[3]  - u[9] *y[2]  + u[14]*y[5]  - u[15]*y[4];     \
      x[4] = u[4] *y[0]  + u[5] *y[1]  + u[10]*y[2]  + u[11]*y[3]  + u[16]*y[4]  + u[17]*y[5];     \
      x[5] = u[4] *y[1]  - u[5] *y[0]  + u[10]*y[3]  - u[11]*y[2]  + u[16]*y[5]  - u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  + u[1] *y[7]  + u[6] *y[8]  + u[7] *y[9]  + u[12]*y[10] + u[13]*y[11];    \
      x[7] = u[0] *y[7]  - u[1] *y[6]  + u[6] *y[9]  - u[7] *y[8]  + u[12]*y[11] - u[13]*y[10];    \
      x[8] = u[2] *y[6]  + u[3] *y[7]  + u[8] *y[8]  + u[9] *y[9]  + u[14]*y[10] + u[15]*y[11];    \
      x[9] = u[2] *y[7]  - u[3] *y[6]  + u[8] *y[9]  - u[9] *y[8]  + u[14]*y[11] - u[15]*y[10];    \
      x[10]= u[4] *y[6]  + u[5] *y[7]  + u[10]*y[8]  + u[11]*y[9]  + u[16]*y[10] + u[17]*y[11];    \
      x[11]= u[4] *y[7]  - u[5] *y[6]  + u[10]*y[9]  - u[11]*y[8]  + u[16]*y[11] - u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] + u[1] *y[13] + u[6] *y[14] + u[7] *y[15] + u[12]*y[16] + u[13]*y[17];    \
      x[13]= u[0] *y[13] - u[1] *y[12] + u[6] *y[15] - u[7] *y[14] + u[12]*y[17] - u[13]*y[16];    \
      x[14]= u[2] *y[12] + u[3] *y[13] + u[8] *y[14] + u[9] *y[15] + u[14]*y[16] + u[15]*y[17];    \
      x[15]= u[2] *y[13] - u[3] *y[12] + u[8] *y[15] - u[9] *y[14] + u[14]*y[17] - u[15]*y[16];    \
      x[16]= u[4] *y[12] + u[5] *y[13] + u[10]*y[14] + u[11]*y[15] + u[16]*y[16] + u[17]*y[17];    \
      x[17]= u[4] *y[13] - u[5] *y[12] + u[10]*y[15] - u[11]*y[14] + u[16]*y[17] - u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] + u[1] *y[19] + u[6] *y[20] + u[7] *y[21] + u[12]*y[22] + u[13]*y[23];    \
      x[19]= u[0] *y[19] - u[1] *y[18] + u[6] *y[21] - u[7] *y[20] + u[12]*y[23] - u[13]*y[22];    \
      x[20]= u[2] *y[18] + u[3] *y[19] + u[8] *y[20] + u[9] *y[21] + u[14]*y[22] + u[15]*y[23];    \
      x[21]= u[2] *y[19] - u[3] *y[18] + u[8] *y[21] - u[9] *y[20] + u[14]*y[23] - u[15]*y[22];    \
      x[22]= u[4] *y[18] + u[5] *y[19] + u[10]*y[20] + u[11]*y[21] + u[16]*y[22] + u[17]*y[23];    \
      x[23]= u[4] *y[19] - u[5] *y[18] + u[10]*y[21] - u[11]*y[20] + u[16]*y[23] - u[17]*y[22];    \
   })



#define qcd_SUM_UP_HOPP(t,a,b)  \
  ({ t[0]  += a[0] + b[0];  \
  t[1]  += a[1] + b[1];  \
  t[2]  += a[2] + b[2];  \
  t[3]  += a[3] + b[3];  \
  t[4]  += a[4] + b[4];  \
  t[5]  += a[5] + b[5];  \
  t[6]  += a[6] + b[6];  \
  t[7]  += a[7] + b[7];  \
  t[8]  += a[8] + b[8];  \
  t[9]  += a[9] + b[9];  \
  t[10] += a[10]+ b[10]; \
  t[11] += a[11]+ b[11]; \
  t[12] += a[12]+ b[12]; \
  t[13] += a[13]+ b[13]; \
  t[14] += a[14]+ b[14]; \
  t[15] += a[15]+ b[15]; \
  t[16] += a[16]+ b[16]; \
  t[17] += a[17]+ b[17]; \
  t[18] += a[18]+ b[18]; \
  t[19] += a[19]+ b[19]; \
  t[20] += a[20]+ b[20]; \
  t[21] += a[21]+ b[21]; \
  t[22] += a[22]+ b[22]; \
  t[23] += a[23]+ b[23]; \
  })




#define qcd_CONJ(x)    ( (Complex) {x.real(),-x.imag()} )
#define qcd_CMUL(x,y)  ( (Complex) {x.real() * y.real() - x.imag() * y.imag(), x.real() * y.imag() + x.imag() * y.real() } )
#define qcd_CMULR(x,y) ( (double) (x.real() * y.real() - x.imag() * y.imag()) )
#define qcd_CMULI(x,y) ( (double) (x.real() * y.imag() + x.imag() * y.real()) )
#define qcd_CADJOINTMUL(x,y)  ( (Complex) {x.real() * y.real() + x.imag() * y.imag(), x.real() * y.imag() - x.imag() * y.real() } )

#define qcd_CADD(x,y)  ( (Complex) {x.real()+y.real(), x.imag()+y.imag()})
#define qcd_CADDR(x,y) ( (double) (x.real()+y.real()))
#define qcd_CADDI(x,y) ( (double) (x.imag()+y.imag()))

#define qcd_CSUB(x,y)  ( (Complex) {x.real()-y.real(), x.imag()-y.imag()})
#define qcd_CSUBR(x,y) ( (double) (x.real()-y.real()))
#define qcd_CSUBI(x,y) ( (double) (x.imag()-y.imag()))

#define qcd_CSCALE(x,a)  ( (Complex) {x.real()*(a), x.imag()*(a)})
#define qcd_CSCALER(x,a) ( (double) x.real()*(a))
#define qcd_CSCALEI(x,a) ( (double) x.imag()*(a))

#define qcd_CDEV(x,y)  ( (Complex) {(x.real() * y.real() + x.imag() * y.imag())/(y.real()*y.real() + y.imag()*y.imag()), (x.imag() * y.real() - x.real() * y.imag())/(y.real()*y.real() + y.imag()*y.imag()) } )
#define qcd_CDEVR(x,y) ( (double)     ((x.real() * y.real() + x.imag() * y.imag())/(y.real()*y.real() + y.imag()*y.imag()) ))
#define qcd_CDEVI(x,y) ( (double)     ((x.imag() * y.real() - x.real() * y.imag())/(y.real()*y.real() + y.imag()*y.imag()) ))

#define qcd_NORM(x)    ( (double) sqrt(x.real() * x.real() + x.imag() * x.imag()))
#define qcd_NORMSQUAREAL()D(x) ( (double) (x.real() * x.real() + x.imag() * x.imag()))
#define qcd_ARG(x)     ( (double) atan2(x.imag(),x.real()))
#define qcd_CPOW(x,a)  ( (Complex) {pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)), pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))})
#define qcd_CPOWR(x,a) ( (double) pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)))
#define qcd_CPOWI(x,a) ( (double) pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a)))

