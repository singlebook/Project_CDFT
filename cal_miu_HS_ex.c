#include "headdefs.h"

fftw_plan p;
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

void Set_FFT(){
	F_Density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Pts.x * Pts.y * (Pts.z/2+1));
	
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
	
	p_n0 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n0,n0,FFTW_ESTIMATE);
	p_n1 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1,n1,FFTW_ESTIMATE);
	p_n2 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n2,n2,FFTW_ESTIMATE);
	p_n3 = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n3,n3,FFTW_ESTIMATE);
	
	p_n1V_x = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_x,n1V_x,FFTW_ESTIMATE);
	p_n1V_y = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_y,n1V_y,FFTW_ESTIMATE);
	p_n1V_z = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_z,n1V_z,FFTW_ESTIMATE);
	
	p_n2V_x = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_x,n1V_x,FFTW_ESTIMATE);
	p_n2V_y = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_y,n1V_y,FFTW_ESTIMATE);
	p_n2V_z = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_n1V_z,n1V_z,FFTW_ESTIMATE);
	
	p_Phi0 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi0,F_Phi0,FFTW_ESTIMATE);
	p_Phi1 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1,F_Phi1,FFTW_ESTIMATE);
	p_Phi2 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi2,F_Phi2,FFTW_ESTIMATE);
	p_Phi3 = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi3,F_Phi3,FFTW_ESTIMATE);
	
	p_Phi1V_x = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_x,F_Phi1V_x,FFTW_ESTIMATE);
	p_Phi1V_y = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_y,F_Phi1V_y,FFTW_ESTIMATE);
	p_Phi1V_z = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_z,F_Phi1V_z,FFTW_ESTIMATE);
	
	p_Phi2V_x = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_x,F_Phi1V_x,FFTW_ESTIMATE);
	p_Phi2V_y = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_y,F_Phi1V_y,FFTW_ESTIMATE);
	p_Phi2V_z = fftw_plan_dft_r2c_3d(Pts.x,Pts.y,Pts.z,Phi1V_z,F_Phi1V_z,FFTW_ESTIMATE);
	
	p_Miu_ex = fftw_plan_dft_c2r_3d(Pts.x,Pts.y,Pts.z,F_Miu_ex,Miu_ex,FFTW_ESTIMATE);
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
		
		K.x = 2 * M_PI * i / Pts.x;
		K.y = 2 * M_PI * j / Pts.y;
		K.z = 2 * M_PI * u / Pts.z;
		
		
		F_n0[loop] = F_Density[loop] * Omega0(VLen(K), Radius) / VProd(Pts);
		F_n1[loop] = F_Density[loop] * Omega1(VLen(K), Radius) / VProd(Pts);
		F_n2[loop] = F_Density[loop] * Omega2(VLen(K), Radius) / VProd(Pts);
		F_n3[loop] = F_Density[loop] * Omega3(VLen(K), Radius) / VProd(Pts);
		
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
	
	fftw_execute(p_n1V_x);
	fftw_execute(p_n1V_y);
	fftw_execute(p_n1V_z);
	
	fftw_execute(p_n2V_x);
	fftw_execute(p_n2V_y);
	fftw_execute(p_n2V_z);
	
	for(loop=0;loop<VProd(Pts);loop++){
		Phi0[loop] = -1.0 * log(1.0 - n3[loop]);
		Phi1[loop] = n2[loop] / (1.0 - n3[loop]);
		Phi2[loop] = n1[loop] / (1.0 - n3[loop]) + (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop])) * (Sqr(n2[loop]) - Sqr(n2V_x[loop]) - Sqr(n2V_y[loop]) - Sqr(n2V_z[loop])) / (12.0 * M_PI * n3[loop]);
		Phi3[loop] = -1.0*(log(1.0-n3[loop])/(18.0*M_PI*Cube(n3[loop])) + (1.0 - 3.0*n3[loop] + Sqr(1.0 - n3[loop]))/(36.0*M_PI*Sqr(n3[loop])*Cube(1.0 - n3[loop]))) * (Cube(n2[loop]) - 3.0*n2[loop]*(Sqr(n2V_x[loop])+Sqr(n2V_y[loop])+Sqr(n2V_z[loop]))) \
		+ n0[loop] / (1.0 - n2[loop]) + (n1[loop]*n2[loop] - n1V_x[loop]*n2V_x[loop] - n1V_y[loop]*n2V_y[loop] - n1V_z[loop]*n2V_z[loop])/Sqr(1.0 - n3[loop]);
		
		Phi1V_x[loop] = -1.0 * n2V_x[loop] / (1.0 - n3[loop]);
		Phi1V_y[loop] = -1.0 * n2V_y[loop] / (1.0 - n3[loop]);
		Phi1V_z[loop] = -1.0 * n2V_z[loop] / (1.0 - n3[loop]);
		
		Phi2V_x[loop] = -1.0 * n1V_x[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_x[loop]/(6.0*M_PI*n3[loop]);
		Phi2V_y[loop] = -1.0 * n1V_y[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_y[loop]/(6.0*M_PI*n3[loop]);
		Phi2V_z[loop] = -1.0 * n1V_z[loop] / (1.0 - n3[loop]) - (log(1.0 - n3[loop])/n3[loop] + 1.0 / Sqr(1.0 - n3[loop]))*n2[loop]*n2V_z[loop]/(6.0*M_PI*n3[loop]);
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
		
		K.x = 2 * M_PI * i / Pts.x;
		K.y = 2 * M_PI * j / Pts.y;
		K.z = 2 * M_PI * u / Pts.z;
		
		F_Miu_ex[loop] = F_Phi0[loop]*Omega0(VLen(K), Radius) + F_Phi1[loop]*Omega1(VLen(K), Radius) + F_Phi2[loop]*Omega2(VLen(K), Radius) + F_Phi3[loop]*Omega3(VLen(K), Radius) \
		+ F_Phi1V_x[loop]*Omega1V(K, Radius).x + F_Phi1V_y[loop]*Omega1V(K, Radius).y + F_Phi1V_z[loop]*Omega1V(K, Radius).z \
		+ F_Phi2V_x[loop]*Omega2V(K, Radius).x + F_Phi2V_y[loop]*Omega2V(K, Radius).y + F_Phi2V_z[loop]*Omega2V(K, Radius).z;
		
		F_Miu_ex[loop] /= VProd(Pts);
		}
		
	fftw_execute(p_Miu_ex);
	
	return;
}

void Clean_FFT(){
	fftw_destroy_plan (p);
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
