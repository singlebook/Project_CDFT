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
			if(AtomType == HS) Density_temp = Atom[0].density * exp(Atom[0].miu - Miu_ex[loop] - Beta * Vext[loop]);
			else if (AtomType == LJ) Density_temp = Atom[0].density * exp(Atom[0].miu - (Miu_ex[loop]+P2(Atom[0].sigma, Atom[0].epslion, rho_bar[loop], "Excess free energy")-F_HS(rho_bar[loop], 2.0*Radius) - F_MFA(rho_bar[loop])) \
			- Beta * Vext[loop]);
			printf("%lf %lf\n",Atom[0].miu , +P2(Atom[0].sigma, Atom[0].epslion, rho_bar[loop], "Excess free energy")-F_HS(rho_bar[loop], 2.0*Radius) - F_MFA(rho_bar[loop]));exit(0);
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
