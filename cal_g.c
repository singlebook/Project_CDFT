#include "headdefs.h"

void Cal_g(){
	real r;
	real count;
	int i, j, k;
	FILE * fp;
	
    fp = fopen("g_r.dat", "w");

    if (fp == NULL) 
    {
    fprintf(stderr, "Can't open input file!\n");
    exit(1);
    }	
    

	for(k=0;k<Pts.z;k++){
		r = k * Delta.z;
		count = 0.0;
		if(r < Rcut){
			for(j=0;j<Pts.y;j++){
				for(i=0;i<Pts.x;i++){
					count += Density[k+Pts.z*(j+Pts.y*i)];
					}
				}
			fprintf(fp, "%lf  %lf\n", r, count/(Pts.x*Pts.y));	
			}
		else break;
		}
	printf("The radial distribution function is calculated.\n");
	return;
	}
