#include "headdefs.h" 

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

void Setup(){
	int i;
	int LCount = 0;
	FILE * fp;
    
    fp = fopen("LJParas.dat", "r");

    if (fp == NULL) 
    {
    fprintf(stderr, "Can't open input file!\n");
    exit(1);
    }
    
	Delta.x = Size.x / (Pts.x - 1.0);
	Delta.y = Size.y / (Pts.y - 1.0);
	Delta.z = Size.z / (Pts.z - 1.0);
	
	
	

//	Beta = 1.0 / (Kb * Na * Temperature);

    Beta = 1.0 / (Kb * Temperature);

	
	AllocMem(ThermWaveLength, NumAtomType, real);
	
	AllocMem(Atom, NumAtomType, struct PropertyAtom);
	
	while (1){
	if (feof (fp)) break;
    fscanf(fp, "%lf %lf %lf %lf %lf\n", &(Atom[LCount].mass), &(Atom[LCount].charge), &(Atom[LCount].sigma), &(Atom[LCount].epslion), &(Atom[LCount].density));
    Atom[LCount].miu = P1(Atom[LCount].sigma, Atom[LCount].density);
    ++ LCount;
	}

	for(i=0;i<NumAtomType;i++){
//		ThermWaveLength[i] = sqrt(Beta * Sqr(h) / (2.0 * M_PI * Atom[i].mass));
        ThermWaveLength[i] = 1.0;
		}
	if(AtomType==LJ)
		Radius = 0.5 *  Atom[0].sigma * (1.0+0.2977*Temperature) / (1.0+0.33163*Temperature+0.0010477*Sqr(Temperature));		
	else if(AtomType==HS)
	    Radius = 0.5 *  Atom[0].sigma;
	    
	AllocMem(Vext, VProd(Pts), real);	
	AllocMem(Density, VProd(Pts), real);
	AllocMem(Ulj, VProd(Pts), real);
	
	SetVext();
	SetInitialDensity();
	SetLJ();
	Set_FFT();
	
	printf ("Setting the calculation is done\n ");	
	return;
	}

void Free_memory(){
	free(ThermWaveLength);
	free(Atom);
	free(Vext);
	free(Density);
	free(Ulj);
	Clean_FFT();
	}
