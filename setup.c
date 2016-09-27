#include "headdefs.h" 


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
    if(AtomType==LJ)
		Atom[LCount].miu = P2(Atom[LCount].sigma, Atom[LCount].epslion, Atom[LCount].density,"Chemical potential");		
	else if(AtomType==HS)
	    Atom[LCount].miu = P1(Atom[LCount].sigma, Atom[LCount].density);
    ++ LCount;
	}

	for(i=0;i<NumAtomType;i++){
//		ThermWaveLength[i] = sqrt(Beta * Sqr(h) / (2.0 * M_PI * Atom[i].mass));
        ThermWaveLength[i] = 1.0;
		}
	if(AtomType==LJ){
		Radius = 0.5 *  Atom[0].sigma * (1.0+0.2977*Temperature) / (1.0+0.33163*Temperature+0.0010477*Sqr(Temperature));
	}		
	else if(AtomType==HS){
	    Radius = 0.5 *  Atom[0].sigma;
	}
	       
	AllocMem(Vext, VProd(Pts), real);	
	AllocMem(Density, VProd(Pts), real);
	
	SetVext();
	SetInitialDensity();
	Set_FFT();
	
	printf ("Setting the calculation is done\n ");	
	return;
	}

void Free_memory(){
	free(ThermWaveLength);
	free(Atom);
	free(Vext);
	free(Density);
	Clean_FFT();
	}
