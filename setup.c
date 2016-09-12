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
	
	Beta = 1.0 / (Kb * Na * Temperature);
	
	AllocMem(ThermWaveLength, NumAtomType, real);
	
	AllocMem(Atom, NumAtomType, struct PropertyAtom);
	
	while (1){
	if (feof (fp)) break;
    fscanf(fp, "%lf %lf %lf %lf %lf\n", &(Atom[LCount].mass), &(Atom[LCount].charge), &(Atom[LCount].sigma), &(Atom[LCount].epslion), &(Atom[LCount].miu));
    ++ LCount;
	}
	
	for(i=0;i<NumAtomType;i++){
		ThermWaveLength[i] = sqrt(Beta * Sqr(h) / (2.0 * M_PI * Atom[i].mass));
		}
		
	AllocMem(Vext, VProd(Pts), real);	
	AllocMem(Density, VProd(Pts), real);	
	return;
	}

void Free_memory(){
	free(ThermWaveLength);
	free(Atom);
	free(Vext);
	free(Density);
	}
