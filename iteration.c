/*
 * The fixed point problem of the density is solved using a Picard iteration.
 * rho_new = (1.0-alpha) * rho_iteration + alpha * rho_old,
 * where alpha is less than 0.9.
 * */

#include "headdefs.h"

void Iteration(){
	int loop;
	real Density_temp;
	real Sum;
	
	while(1){
		Sum = 0.0;
		Cal_Miu_HS_ex();
		for(loop=0;loop<VProd(Pts);loop++){
			Density_temp = Atom[0].density * exp(Atom[0].miu - Miu_ex[loop] - Beta * Vext[loop]);
			Sum +=Sqr((1.0-Alpha) * Density_temp + Alpha * Density[loop] - Density[loop]);
			Density[loop] = (1.0-Alpha) * Density_temp + Alpha * Density[loop];
			}
		printf("Sum = %lf\n", Sum);	
		if (Sum < Stop){
			printf("The iteration is successfully finished!\n");
			break;
			}
		}
	return;	
	}
