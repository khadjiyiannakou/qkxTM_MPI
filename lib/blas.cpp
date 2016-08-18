#include <qkxTM.h>
#include <lattice_util.h>

using namespace qkxTM;

void qcd_projectSU33d(Gauge &gf)
{

   int iv,mu,c1,c2;
   Complex H[3][3],detM,U[3][3],M[3][3];
   Complex b,D,v[3][3],vr[3][3];
   double a,ThirdRoot_18,ThirdRoot_12,ThirdRoot_2_3;
   double trace,e[3],de,cor[3];
   Complex w,ThirdRootOne[2];
   double  sum;
   double  norm,phase;

   /*
     Constants Used:
   */

   ThirdRootOne[0].real()= 1;
   ThirdRootOne[0].imag()= sqrt(3.);

   ThirdRootOne[1].real()= 1;
   ThirdRootOne[1].imag()=-sqrt(3.);

   ThirdRoot_12 =pow(12.,1./3.);
   ThirdRoot_18 =pow(18.,1./3.);
   ThirdRoot_2_3=pow((2./3.),1./3.);

   for(iv=0; iv < gf.lV4d ; iv++)
   for(mu=0; mu<3; mu++) 
   {
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
         M[c1][c2] = gf.M[iv][mu][c1][c2];

      
      detM = qcd_CADD(qcd_CADD( qcd_CMUL(M[0][0],qcd_CMUL(M[1][1],M[2][2])),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][2],M[2][0]))),
                                qcd_CMUL(M[0][2],qcd_CMUL(M[1][0],M[2][1])));
                                
      detM = qcd_CSUB(detM,
             qcd_CADD(qcd_CADD( qcd_CMUL(M[0][2],qcd_CMUL(M[1][1],M[2][0])),
                                qcd_CMUL(M[0][0],qcd_CMUL(M[1][2],M[2][1]))),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][0],M[2][2]))));


      phase = qcd_ARG(detM)/3.;



      H[0][0].real()= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].real()= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].real()= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][2]);
      H[0][0].imag()= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].imag()= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].imag()= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][2]);

      H[1][0].real()= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].real()= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].real()= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][2]);
      H[1][0].imag()= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].imag()= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].imag()= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][2]);

      H[2][0].real()= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].real()= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].real()= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][2]);
      H[2][0].imag()= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].imag()= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].imag()= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][2]);



      /*
        Assureal() Hermiticity:
      */
      H[0][1].real() = (H[0][1].real() + H[1][0].real())/2.;
      H[0][1].imag() = (H[0][1].imag() - H[1][0].imag())/2.;

      H[1][0] =  qcd_CONJ(H[0][1]);

      H[0][2].real() = (H[0][2].real() + H[2][0].real())/2.;
      H[0][2].imag() = (H[0][2].imag() - H[2][0].imag())/2.;

      H[2][0] =  qcd_CONJ(H[0][2]);

      H[1][2].real() = (H[1][2].real() + H[2][1].real())/2.;
      H[1][2].imag() = (H[1][2].imag() - H[2][1].imag())/2.;

      H[2][1] =  qcd_CONJ(H[1][2]);



      /*
        If H^2 is alreal()ad diagonal skip diagonalization and
        calculate U direal()ctly
      */
      sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);

      if(sum<=1e-08)
      {
         e[0]=1./sqrt(H[0][0].real());
         e[1]=1./sqrt(H[1][1].real());
         e[2]=1./sqrt(H[2][2].real());

         U[0][0] = (Complex) { M[0][0].real()*e[0], M[0][0].imag()*e[0] };
         U[0][1] = (Complex) { M[0][1].real()*e[0], M[0][1].imag()*e[0] };
         U[0][2] = (Complex) { M[0][2].real()*e[0], M[0][2].imag()*e[0] };

         U[1][0] = (Complex) { M[1][0].real()*e[1], M[1][0].imag()*e[1] };
         U[1][1] = (Complex) { M[1][1].real()*e[1], M[1][1].imag()*e[1] };
         U[1][2] = (Complex) { M[1][2].real()*e[1], M[1][2].imag()*e[1] };

         U[2][0] = (Complex) { M[2][0].real()*e[2], M[2][0].imag()*e[2] };
         U[2][1] = (Complex) { M[2][1].real()*e[2], M[2][1].imag()*e[2] };
         U[2][2] = (Complex) { M[2][2].real()*e[2], M[2][2].imag()*e[2] };

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf.M[iv][mu][c1][c2] = U[c1][c2]; 
      }
      else
      {
         /*
           Make traceless to elimag()inate second order term in eigenvalue equation,
           i.e. eig^3 + a eig + b = 0, when H is traceless.
         */
         trace=(H[0][0].real()+H[1][1].real()+H[2][2].real())/3.;

	 
	 


         H[0][0].real()-=trace;
         H[1][1].real()-=trace;
         H[2][2].real()-=trace;


         /*
           Solve for eigenvalues:
           e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
           e^3 + a*e + b = 0,

           a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
           b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

           D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

           e = D/(18^(1/3)) - ((2/3)^(1/3))/D
           e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
           e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
         */
         a = -(H[2][2].real()*H[2][2].real() - H[0][0].real()*H[1][1].real() + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));


         b.real()  = - H[0][0].real()*H[1][1].real()*H[2][2].real() + H[2][2].real()*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].real()*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
         b.imag()  =                                      H[2][2].real()*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].real()*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

         b.real() +=   H[0][0].real()*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
         b.imag() +=   H[0][0].real()*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

         w.real()=qcd_CPOWR(((Complex){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
         w.imag()=qcd_CPOWI(((Complex){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

         D=qcd_CPOW(((Complex){-9.*b.real() + w.real(),-9.*b.imag() + w.imag()}),1./3.); 

         e[0] = D.real()/(ThirdRoot_18) - qcd_CDEVR(((Complex){a*ThirdRoot_2_3,0}),D);
         e[1] = a*qcd_CDEVR(ThirdRootOne[0],((Complex){D.real()*ThirdRoot_12,D.imag()*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
         e[2] = -e[0]-e[1];

         e[0]+= trace;
         e[1]+= trace;
         e[2]+= trace;

         H[0][0].real()+=trace;
         H[1][1].real()+=trace;
         H[2][2].real()+=trace;

         /*
           Eigenvectors:
           v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
           v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
           v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
         */

         v[0][0].real() = -(e[0]*H[2][0].real() - H[2][0].real()*H[1][1].real() + qcd_CMULR(H[1][0],H[2][1]));
         v[0][0].imag() = -(e[0]*H[2][0].imag() - H[2][0].imag()*H[1][1].real() + qcd_CMULI(H[1][0],H[2][1]));

         v[0][1].real() = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].real() - H[0][0].real()*H[2][1].real());
         v[0][1].imag() = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].imag() - H[0][0].real()*H[2][1].imag());

         v[0][2].real() =-e[0]*e[0] + e[0]*H[0][0].real() + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].real() - H[0][0].real()*H[1][1].real();
         v[0][2].imag() = 0.;

         v[1][0].real() = -(e[1]*H[2][0].real() - H[2][0].real()*H[1][1].real() + qcd_CMULR(H[1][0],H[2][1]));
         v[1][0].imag() = -(e[1]*H[2][0].imag() - H[2][0].imag()*H[1][1].real() + qcd_CMULI(H[1][0],H[2][1]));

         v[1][1].real() = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].real() - H[0][0].real()*H[2][1].real());
         v[1][1].imag() = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].imag() - H[0][0].real()*H[2][1].imag());

         v[1][2].real() =-e[1]*e[1] + e[1]*H[0][0].real() + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].real() - H[0][0].real()*H[1][1].real();;
         v[1][2].imag() = 0.;

         /*
           Assureal() eigenvectors orthonormality:
           norm =  inner product v1.v1
           w    = (inner product v1.v2)/norm
           v2   = w*v1
         */

         norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
         w.real()  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
         w.imag()  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
         w.real() /= norm;
         w.imag() /= norm;

         v[1][0].real()-= qcd_CMULR(w,v[0][0]);
         v[1][0].imag()-= qcd_CMULI(w,v[0][0]);

         v[1][1].real()-= qcd_CMULR(w,v[0][1]);
         v[1][1].imag()-= qcd_CMULI(w,v[0][1]);

         v[1][2].real()-= qcd_CMULR(w,v[0][2]);
         v[1][2].imag()-= qcd_CMULI(w,v[0][2]);

         norm=1./sqrt(norm);

         /*
           Normalize first and second eigenvector:
         */

         v[0][0].real()*= norm;
         v[0][0].imag()*= norm;

         v[0][1].real()*= norm;
         v[0][1].imag()*= norm;

         v[0][2].real()*= norm;
         v[0][2].imag()*= norm;


         norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

         norm=1./sqrt(norm);

         v[1][0].real()*= norm;
         v[1][0].imag()*= norm;

         v[1][1].real()*= norm;
         v[1][1].imag()*= norm;

         v[1][2].real()*= norm;
         v[1][2].imag()*= norm;

         /*
           v3 = v1 x v2
         */


         v[2][0].real() =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
         v[2][0].imag() = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

         v[2][1].real() = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
         v[2][1].imag() = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

         v[2][2].real() =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
         v[2][2].imag() = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

         de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
         /*
         cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
         cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
         cor[2] = - cor[0] - cor[1];
         */
         //to be compatible with Greal()noble & Paris, don't apply correal()ctions
         cor[0]=0;
         cor[1]=0;
         cor[2]=0;

         de = 1./sqrt(e[0]);
         b.real() = de*cos(phase-cor[0]);
         b.imag() =-de*sin(phase-cor[0]);
         vr[0][0] = qcd_CMUL(b,v[0][0]);
         vr[0][1] = qcd_CMUL(b,v[0][1]);
         vr[0][2] = qcd_CMUL(b,v[0][2]);

         de = 1./sqrt(e[1]);
         b.real() = de*cos(phase-cor[1]);
         b.imag() =-de*sin(phase-cor[1]);

         vr[1][0] = qcd_CMUL(b,v[1][0]);
         vr[1][1] = qcd_CMUL(b,v[1][1]);
         vr[1][2] = qcd_CMUL(b,v[1][2]);

         de = 1./sqrt(e[2]);
         b.real() = de*cos(phase-cor[2]);
         b.imag() =-de*sin(phase-cor[2]);

         vr[2][0] = qcd_CMUL(b,v[2][0]);
         vr[2][1] = qcd_CMUL(b,v[2][1]);
         vr[2][2] = qcd_CMUL(b,v[2][2]);


         H[0][0].real()= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].real()= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].real()= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

         H[0][0].imag()= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].imag()= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].imag()= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


         H[1][0].real()= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].real()= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].real()= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

         H[1][0].imag()= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].imag()= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].imag()= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


         H[2][0].real()= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].real()= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].real()= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

         H[2][0].imag()= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].imag()= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].imag()= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

         U[0][0].real()= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
         U[0][1].real()= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
         U[0][2].real()= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

         U[0][0].imag()= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
         U[0][1].imag()= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
         U[0][2].imag()= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


         U[1][0].real()= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
         U[1][1].real()= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
         U[1][2].real()= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

         U[1][0].imag()= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
         U[1][1].imag()= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
         U[1][2].imag()= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


         U[2][0].real()= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
         U[2][1].real()= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
         U[2][2].real()= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

         U[2][0].imag()= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
         U[2][1].imag()= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
         U[2][2].imag()= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


         /*
           w    = inner product: col1.col2
           norm = inner product: col1.col1
         */

         norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
         w.real()  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
         w.imag()  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
         w.real() /= norm;
         w.imag() /= norm;


         U[0][1].real()-=qcd_CMULR(w,U[0][0]);
         U[0][1].imag()-=qcd_CMULI(w,U[0][0]);

         U[1][1].real()-=qcd_CMULR(w,U[1][0]);
         U[1][1].imag()-=qcd_CMULI(w,U[1][0]);

         U[2][1].real()-=qcd_CMULR(w,U[2][0]);
         U[2][1].imag()-=qcd_CMULI(w,U[2][0]);

         norm = 1./sqrt(norm);

         U[0][0].real()*= norm;
         U[0][0].imag()*= norm;
         U[1][0].real()*= norm;
         U[1][0].imag()*= norm;
         U[2][0].real()*= norm;
         U[2][0].imag()*= norm;

         norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
         norm = 1./sqrt(norm);

         U[0][1].real()*= norm;
         U[0][1].imag()*= norm;
         U[1][1].real()*= norm;
         U[1][1].imag()*= norm;
         U[2][1].real()*= norm;
         U[2][1].imag()*= norm;

         /*
           col3 = col1 x col2
         */
         U[0][2].real() =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
         U[0][2].imag() = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

         U[1][2].real() = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
         U[1][2].imag() =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

         U[2][2].real() =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
         U[2][2].imag() = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf.M[iv][mu][c1][c2] = U[c1][c2]; 

	 //	 if(iv==0)printf("%18.16e\n",gf->D[iv][mu][0][0].real()); 
      //exit(-1);

      }
   }//end volume-mu-loop
   return;
}//end qcd_projectSU33d

