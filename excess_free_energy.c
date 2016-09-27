#include "headdefs.h"

real F_MFA(real rho){
     real t;
	 // the attractive part of the correlation between repulsion and attraction 
	 t = -16.0*M_PI*Beta*Atom[0].epslion*Sqr(rho)*Cube(Atom[0].sigma)/9.0;
     return Beta*t;
}


real F_HS(real rho, real sigma){
     real eta = M_PI*Cube(sigma)*rho/6.0;
	 real t;
	 t = -1.0*log(1-eta)+6*eta/(1.0-eta)+Sqr(3.0*eta)/(2.0*Sqr(1.0-eta)) + (eta-1.0)*(1.0+eta+Sqr(eta))/Cube(1.0-eta)+1.0;
	 return t;
}

/*The excess chemical potential for HS calculated from SPT*/
real P1(real sigma, real rho_0) {
	real P0;
	real t0;
	real xi_m1, xi_m2, xi_m3;
	
	xi_m1 = M_PI*rho_0*sigma/6.0;
	xi_m2 = M_PI*rho_0*sigma*sigma/6.0;
	xi_m3 = M_PI*rho_0*sigma*sigma*sigma/6.0;
	
	
	P0 =rho_0 * (1.0+xi_m3+Sqr(xi_m3))/Cube(1.0-xi_m3);

	t0 = -log(1.0-xi_m3)+3.0*xi_m2/(1.0-xi_m3)*sigma+(3.0*xi_m1/(1.0-xi_m3)+9.0/2.0*Sqr(xi_m2)/Sqr(1.0-xi_m3))*sigma*sigma+M_PI*sigma*sigma*sigma*P0/6.0;

	return t0;	
}

/*The excess chemical potential as calculated from the EOS of Johnson
 *  Molecular Physics 78(3):591-618 · February 1993
 * All quantities are in reduced units.
 * */
real P2(real sigma, real epslion, real rho, char name[50]){
	int i;
	real x[32] = {
	0.8623085097507421,\
	2.976218765822098,\
	-8.402230115796038,\
	0.1054136629203555,\
	-0.8564583828174598,\
	1.582759470107601,\
	0.7639421948305453,\
	1.753173414312048,\
	2.798291772190376E3,\
	-4.8394220260857657E-2,\
	0.9963265197721935,\
	-3.698000291272493E1,\
	2.084012299434647E1,\
	8.305402124717285E1,\
	-9.574799715203068E2,\
	-1.477746229234994E2,\
	6.398607852471505E1,\
	1.603993673294834E1,\
	6.805916615864377E1,\
	-2.791293578795945E3,\
	-6.245128304568454,\
	-8.116836104958410E3,\
	1.488735559561229E1,\
	-1.059346754655084E4,\
	-1.131607632802822E2,\
	-8.867771540418822E3,\
	-3.986982844450543E1,\
	-4.689270299917261E3,\
	2.593535277438717E2,\
	-2.694523589434903E3,\
	-7.218487631550215E2,\
	1.721802063863269E2};
	
	real gama = 3.0;
	real F = exp(-gama*Sqr(rho));
	real a[8], b[6], G[6];
	real P, A, miu;
	
	a[0] = x[0]*Temperature + x[1]*sqrt(Temperature) + x[2] + x[3]/Temperature + x[4]/Sqr(Temperature);
	a[1] = x[5]*Temperature + x[6] + x[7]/Temperature + x[8]/Sqr(Temperature);
	a[2] = x[9]*Temperature + x[10] + x[11]/Temperature;
	a[3] = x[12];
	a[4] = x[13]/Temperature + x[14]/Sqr(Temperature);
	a[5] = x[15]/Temperature;
	a[6] = x[16]/Temperature + x[17]/Sqr(Temperature);
	a[7] = x[18]/Sqr(Temperature);
	
    b[0] = x[19]/Sqr(Temperature) + x[20]/Cube(Temperature);
    b[1] = x[21]/Sqr(Temperature) + x[22]/Sqr(Sqr(Temperature));
    b[2] = x[23]/Sqr(Temperature) + x[24]/Cube(Temperature);
    b[3] = x[25]/Sqr(Temperature) + x[26]/Sqr(Sqr(Temperature));
    b[4] = x[27]/Sqr(Temperature) + x[28]/Cube(Temperature);
    b[5] = x[29]/Sqr(Temperature) + x[30]/Cube(Temperature) + x[31]/Sqr(Sqr(Temperature));	
	
	G[0] = (1.0 - F)/(2.0*gama);
	for(i=1;i<6;i++){
		G[i] = -1.0*(F*pow(rho,2*i) - 2*i*G[i-1])/(2.0*gama);
		}
	P=0.0;
	A=0.0;	
	for(i=0;i<8;i++){
		if(i<6){
	    P += a[i]*pow(rho,i+2) + F*b[i]*pow(rho,2*i+3);
	    A += a[i]*pow(rho,i+1)/(i+1) + b[i]*G[i];  	
	    }
	    else{
			 P += a[i]*pow(rho,i+2);
			 A += a[i]*pow(rho,i+1)/(i+1); 
		 }
	}
	
	P += rho*Temperature;	
	miu = A + P/rho - Temperature;
	
	if(strcmp(name,"Chemical potential")==0)
		return Beta*miu*Atom[0].epslion;
	else if(strcmp(name,"Excess free energy")==0)
		return Beta*A*Atom[0].epslion;	
}
