#include "headdefs.h"

fftw_plan p;
fftw_plan p_rho_bar;
fftw_plan p_n0;
fftw_plan p_n1;
fftw_plan p_n2;
fftw_plan p_n3;
fftw_plan p_n1V_x;
fftw_plan p_n1V_y;
fftw_plan p_n1V_z;
fftw_plan p_n2V_x;
fftw_plan p_n2V_y;
fftw_plan p_n2V_z;
fftw_plan p_Phi0;
fftw_plan p_Phi1;
fftw_plan p_Phi2;
fftw_plan p_Phi3;
fftw_plan p_Phi1V_x;
fftw_plan p_Phi1V_y;
fftw_plan p_Phi1V_z;
fftw_plan p_Phi2V_x;
fftw_plan p_Phi2V_y;
fftw_plan p_Phi2V_z;
fftw_plan p_Miu_ex;

/*The Integrand of Eq.(17) in JCP 124,144709(2006)*/
real Integrand(real r, real k){
	real t;
	t = Cube(Atom[0].sigma/r);
	if(k!=0)
	return 4.0*M_PI*sin(k*r)*r*(4.0*Atom[0].epslion*(Sqr(Sqr(t))-Sqr(t)))/k;
	else
	return 4.0*M_PI*Sqr(r)*(4.0*Atom[0].epslion*(Sqr(Sqr(t))-Sqr(t)));	
	}

/*Simpson method for calculating the Integrand defined above*/
real simp(real a,real b,real c,real eps)
  { int n,k;
    real h,t1,t2,s1,s2,ep,p,x;
    n=1; h=b-a;
    t1=h*(Integrand(a,c)+Integrand(b,c))/2.0;
    s1=t1;
    ep=eps+1.0;
    while (ep>=eps)
      { p=0.0;
        for (k=0;k<=n-1;k++)
          { x=a+(k+0.5)*h;
            p=p+Integrand(x,c);
          }
        t2=(t1+h*p)/2.0;
        s2=(4.0*t2-t1)/3.0;
        ep=fabs(s2-s1);
        t1=t2; s1=s2; n=n+n; h=h/2.0;
      }
    return(s2);
  }

void Set_FFT(){
	
	int i,j,u,loop;
	struct Vector K;
	int Length = Pts.x*Pts.y*(Pts.z/2+1);
	
	F_Density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_rho_bar = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
    AllocMem(F_Ulj, Pts.x * Pts.y * (Pts.z/2+1), real);
	
	F_n0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
	F_n1V_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n1V_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n1V_z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
    F_n2V_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n2V_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_n2V_z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
	n0 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n1 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n2 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
    n3 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
    
    n1V_x = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n1V_y = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n1V_z = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	
	n2V_x = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n2V_y = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	n2V_z = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	
	
	F_Phi0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
	F_Phi1V_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi1V_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi1V_z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
    F_Phi2V_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi2V_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	F_Phi2V_z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
	Phi0 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi1 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi2 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
    Phi3 = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
    
    Phi1V_x = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi1V_y = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi1V_z = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	
	Phi2V_x = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi2V_y = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	Phi2V_z = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	
	Miu_ex = (real*)malloc(sizeof(real)*Pts.x * Pts.y * Pts.z);
	F_Miu_ex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
	p = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Density,F_Density,FFTW_ESTIMATE);
	

    // make a table for the fourie transform of the LJ 12-6 potential.  
    // There is no analytical form for the fourie transform of the LJ 12-6 potential, I could only prepare it at the very beginning of the simulation.
    if(AtomType==LJ){
	for(loop=0;loop<Length;loop++){
		u = loop % (Pts.z/2+1);
		j = (loop / (Pts.z/2+1)) % Pts.y;
		i = (loop / (Pts.z/2+1)) / Pts.y;
		
		K.x = 2.0 * M_PI * i / Size.x;
		K.y = 2.0 * M_PI * j / Size.y;
		K.z = 2.0 * M_PI * u / Size.z;
		
		F_Ulj[loop] = simp(Atom[0].sigma,Rcut,VLen(K),Stop);
		}
	}
		
	p_n0 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n0,n0,FFTW_ESTIMATE);
	p_n1 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1,n1,FFTW_ESTIMATE);
	p_n2 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n2,n2,FFTW_ESTIMATE);
	p_n3 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n3,n3,FFTW_ESTIMATE);
	
	p_rho_bar = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_rho_bar,rho_bar,FFTW_ESTIMATE);
	
	p_n1V_x = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_x,n1V_x,FFTW_ESTIMATE);
	p_n1V_y = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_y,n1V_y,FFTW_ESTIMATE);
	p_n1V_z = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_z,n1V_z,FFTW_ESTIMATE);
	
	p_n2V_x = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n2V_x,n2V_x,FFTW_ESTIMATE);
	p_n2V_y = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n2V_y,n2V_y,FFTW_ESTIMATE);
	p_n2V_z = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n2V_z,n2V_z,FFTW_ESTIMATE);
	
	p_Phi0 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi0,F_Phi0,FFTW_ESTIMATE);
	p_Phi1 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1,F_Phi1,FFTW_ESTIMATE);
	p_Phi2 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi2,F_Phi2,FFTW_ESTIMATE);
	p_Phi3 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi3,F_Phi3,FFTW_ESTIMATE);
	
	p_Phi1V_x = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_x,F_Phi1V_x,FFTW_ESTIMATE);
	p_Phi1V_y = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_y,F_Phi1V_y,FFTW_ESTIMATE);
	p_Phi1V_z = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_z,F_Phi1V_z,FFTW_ESTIMATE);
	
	p_Phi2V_x = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi2V_x,F_Phi2V_x,FFTW_ESTIMATE);
	p_Phi2V_y = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi2V_y,F_Phi2V_y,FFTW_ESTIMATE);
	p_Phi2V_z = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi2V_z,F_Phi2V_z,FFTW_ESTIMATE);
	
	p_Miu_ex = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_Miu_ex,Miu_ex,FFTW_ESTIMATE);
	printf("Setting FFT is done\n");
	return;
}

void Cal_Miu_HS_ex(){
	int i,j,u,loop;
	struct Vector K;
	int Length = Pts.x*Pts.y*(Pts.z/2+1);
	fftw_execute(p); /* repeat as needed */

	for(loop=0;loop<Length;loop++){
		u = loop % (Pts.z/2+1);
		j = (loop / (Pts.z/2+1)) % Pts.y;
		i = (loop / (Pts.z/2+1)) / Pts.y;
		
		K.x = 2.0 * M_PI * i / Size.x;
		K.y = 2.0 * M_PI * j / Size.y;
		K.z = 2.0 * M_PI * u / Size.z;
		
		
		F_n0[loop] = F_Density[loop] * Omega0(VLen(K), Radius) / VProd(Pts);
		F_n1[loop] = F_Density[loop] * Omega1(VLen(K), Radius) / VProd(Pts);
		F_n2[loop] = F_Density[loop] * Omega2(VLen(K), Radius) / VProd(Pts);
		F_n3[loop] = F_Density[loop] * Omega3(VLen(K), Radius) / VProd(Pts);
		
		F_rho_bar[loop] = (3.0/(4.0*M_PI*Cube(2.0*Radius)))*F_Density[loop] * Omega3(VLen(K), 2.0*Radius) / VProd(Pts);
		
		F_n1V_x[loop] = F_Density[loop] * Omega1V(K, Radius).x / VProd(Pts);
		F_n1V_y[loop] = F_Density[loop] * Omega1V(K, Radius).y / VProd(Pts);
		F_n1V_z[loop] = F_Density[loop] * Omega1V(K, Radius).z / VProd(Pts);
		
		F_n2V_x[loop] = F_Density[loop] * Omega2V(K, Radius).x / VProd(Pts);
		F_n2V_y[loop] = F_Density[loop] * Omega2V(K, Radius).y / VProd(Pts);
		F_n2V_z[loop] = F_Density[loop] * Omega2V(K, Radius).z / VProd(Pts);	
		}

	fftw_execute(p_n0); /* repeat as needed */
	fftw_execute(p_n1);
	fftw_execute(p_n2);
	fftw_execute(p_n3);
	
	fftw_execute(p_rho_bar);
	
	fftw_execute(p_n1V_x);
	fftw_execute(p_n1V_y);
	fftw_execute(p_n1V_z);
	
	fftw_execute(p_n2V_x);
	fftw_execute(p_n2V_y);
	fftw_execute(p_n2V_z);

		
	for(loop=0;loop<VProd(Pts);loop++){
		Phi0[loop] = -1.0 * log(1.0 - n3[loop]);
		Phi1[loop] = n2[loop] / (1.0 - n3[loop]);

/*		// Shuangliang Zhao
		Phi2[loop] = n1[loop] / (1.0 - n3[loop]) + (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop])) * (Sqr(n2[loop]) - Sqr(n2V_x[loop]) - Sqr(n2V_y[loop]) - Sqr(n2V_z[loop])) / (12.0 * M_PI * n3[loop]);
		Phi3[loop] = -1.0*(log(1.0-n3[loop])/(18.0*M_PI*Cube(n3[loop])) + (1.0 - 3.0*n3[loop] + Sqr(1.0 - n3[loop]))/(36.0*M_PI*Sqr(n3[loop])*Cube(1.0 - n3[loop]))) * (Cube(n2[loop]) - 3.0*n2[loop]*(Sqr(n2V_x[loop])+Sqr(n2V_y[loop])+Sqr(n2V_z[loop]))) \
		+ n0[loop] / (1.0 - n2[loop]) + (n1[loop]*n2[loop] - n1V_x[loop]*n2V_x[loop] - n1V_y[loop]*n2V_y[loop] - n1V_z[loop]*n2V_z[loop])/Sqr(1.0 - n3[loop]);
*/		
        // M. P. Sears
        Phi2[loop] = n1[loop] / (1.0 - n3[loop]) + Sqr(n2[loop])/(8.0*M_PI*Sqr(1.0-n3[loop]))-(Sqr(n2V_x[loop]) + Sqr(n2V_y[loop]) + Sqr(n2V_z[loop]))/(8.0*M_PI*Sqr(1.0-n3[loop]));
        Phi3[loop] = n0[loop] / (1.0 - n3[loop]) + n1[loop]*n2[loop]/Sqr(1.0-n3[loop]) + Cube(n2[loop])/(12.0*M_PI*Cube(1.0-n3[loop])) - (n1V_x[loop]*n2V_x[loop] + n1V_y[loop]*n2V_y[loop] + n1V_z[loop]*n2V_z[loop]) / Sqr(1.0-n3[loop])\
         - n2[loop]* (Sqr(n2V_x[loop]) + Sqr(n2V_y[loop]) + Sqr(n2V_z[loop]))/(4.0*M_PI*Cube(1.0-n3[loop]));
        
         
		Phi1V_x[loop] = -1.0 * n2V_x[loop] / (1.0 - n3[loop]);
		Phi1V_y[loop] = -1.0 * n2V_y[loop] / (1.0 - n3[loop]);
		Phi1V_z[loop] = -1.0 * n2V_z[loop] / (1.0 - n3[loop]);
/*		
        // Shuabngliang Zhao
		Phi2V_x[loop] = -1.0 * n1V_x[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_x[loop]/(6.0*M_PI*n3[loop]);
		Phi2V_y[loop] = -1.0 * n1V_y[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_y[loop]/(6.0*M_PI*n3[loop]);
		Phi2V_z[loop] = -1.0 * n1V_z[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_z[loop]/(6.0*M_PI*n3[loop]);
*/		
        // M.P. Sears
		Phi2V_x[loop] = -1.0 * n1V_x[loop] / (1.0 - n3[loop]) - n2[loop]*n2V_x[loop]/(4.0*M_PI*Sqr(1.0-n3[loop]));
		Phi2V_y[loop] = -1.0 * n1V_y[loop] / (1.0 - n3[loop]) - n2[loop]*n2V_y[loop]/(4.0*M_PI*Sqr(1.0-n3[loop]));
		Phi2V_z[loop] = -1.0 * n1V_z[loop] / (1.0 - n3[loop]) - n2[loop]*n2V_z[loop]/(4.0*M_PI*Sqr(1.0-n3[loop]));
		}
	
	fftw_execute(p_Phi0); /* repeat as needed */
	fftw_execute(p_Phi1);
	fftw_execute(p_Phi2);
	fftw_execute(p_Phi3);
	
	fftw_execute(p_Phi1V_x);
	fftw_execute(p_Phi1V_y);
	fftw_execute(p_Phi1V_z);
	
	fftw_execute(p_Phi2V_x);
	fftw_execute(p_Phi2V_y);
	fftw_execute(p_Phi2V_z);
	
	for(loop=0;loop<Length;loop++){
		u = loop % (Pts.z/2+1);
		j = (loop / (Pts.z/2+1)) % Pts.y;
		i = (loop / (Pts.z/2+1)) / Pts.y;
		
		K.x = 2 * M_PI * i / Size.x;
		K.y = 2 * M_PI * j / Size.y;
		K.z = 2 * M_PI * u / Size.z;
		
		if(AtomType == LJ){
		F_Miu_ex[loop] = F_Phi0[loop]*Omega0(VLen(K), Radius) +F_Phi1[loop]*Omega1(VLen(K), Radius) + F_Phi2[loop]*Omega2(VLen(K), Radius) + F_Phi3[loop]*Omega3(VLen(K), Radius) \
		+ F_Phi1V_x[loop]*Omega1V(K, Radius).x + F_Phi1V_y[loop]*Omega1V(K, Radius).y + F_Phi1V_z[loop]*Omega1V(K, Radius).z \
		+ F_Phi2V_x[loop]*Omega2V(K, Radius).x + F_Phi2V_y[loop]*Omega2V(K, Radius).y + F_Phi2V_z[loop]*Omega2V(K, Radius).z \
		+ Beta*F_Density[loop]*F_Ulj[loop];
	    } 
	    else if(AtomType == HS){
		F_Miu_ex[loop] = F_Phi0[loop]*Omega0(VLen(K), Radius) +F_Phi1[loop]*Omega1(VLen(K), Radius) + F_Phi2[loop]*Omega2(VLen(K), Radius) + F_Phi3[loop]*Omega3(VLen(K), Radius) \
		+ F_Phi1V_x[loop]*Omega1V(K, Radius).x + F_Phi1V_y[loop]*Omega1V(K, Radius).y + F_Phi1V_z[loop]*Omega1V(K, Radius).z \
		+ F_Phi2V_x[loop]*Omega2V(K, Radius).x + F_Phi2V_y[loop]*Omega2V(K, Radius).y + F_Phi2V_z[loop]*Omega2V(K, Radius).z;			
		}
		F_Miu_ex[loop] /= VProd(Pts);
		}
		
	fftw_execute(p_Miu_ex);
	return;
}

void Clean_FFT(){
	fftw_destroy_plan (p);
	fftw_destroy_plan (p_rho_bar);
	fftw_destroy_plan (p_n0);
	fftw_destroy_plan (p_n1);
	fftw_destroy_plan (p_n2);
	fftw_destroy_plan (p_n3);
	fftw_destroy_plan (p_n1V_x);
	fftw_destroy_plan (p_n1V_y);
	fftw_destroy_plan (p_n1V_z);
	fftw_destroy_plan (p_n2V_x);
	fftw_destroy_plan (p_n2V_y);
	fftw_destroy_plan (p_n2V_z);
	fftw_destroy_plan (p_Phi0);
	fftw_destroy_plan (p_Phi1);
	fftw_destroy_plan (p_Phi2);
	fftw_destroy_plan (p_Phi3);
	fftw_destroy_plan (p_Phi1V_x);
	fftw_destroy_plan (p_Phi1V_y);
	fftw_destroy_plan (p_Phi1V_z);
	fftw_destroy_plan (p_Phi2V_x);
	fftw_destroy_plan (p_Phi2V_y);
	fftw_destroy_plan (p_Phi2V_z);
	fftw_destroy_plan (p_Miu_ex);

    fftw_cleanup();
    
    free(n0);
    free(n1);
    free(n2);
    free(n3);
    free(n1V_x);
    free(n1V_y);
    free(n1V_z);
    free(n2V_x);
    free(n2V_y);
    free(n2V_z);
    
    free(Phi0);
    free(Phi1);
    free(Phi2);
    free(Phi3);
    free(Phi1V_x);
    free(Phi1V_y);
    free(Phi1V_z);
    free(Phi2V_x);
    free(Phi2V_y);
    free(Phi2V_z);
    
    fftw_free(F_Density);
    fftw_free(F_rho_bar);
    free(F_Ulj);
    fftw_free(F_n0);
    fftw_free(F_n1);
    fftw_free(F_n2);
    fftw_free(F_n3);
    fftw_free(F_n1V_x);
    fftw_free(F_n1V_y);
    fftw_free(F_n1V_z);
    fftw_free(F_n2V_x);
    fftw_free(F_n2V_y);
    fftw_free(F_n2V_z);

    fftw_free(F_Phi0);
    fftw_free(F_Phi1);
    fftw_free(F_Phi2);
    fftw_free(F_Phi3);
    fftw_free(F_Phi1V_x);
    fftw_free(F_Phi1V_y);
    fftw_free(F_Phi1V_z);
    fftw_free(F_Phi2V_x);
    fftw_free(F_Phi2V_y);
    fftw_free(F_Phi2V_z);
    
    free(Miu_ex);
    fftw_free(F_Miu_ex);
    
    return;
	}
