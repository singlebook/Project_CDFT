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
	return;
}

void FFT_R2C(){
	int i,j,u,loop;
	struct Vector K
	int Lenght = Pts.x*Pts.y*(Pts.z/2+1);
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
	return;
}

void FFT_C2R(){
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
	}
