#include "headdefs.h"

// Set up the external potential
void SetVext(){
	
	int loop;
	/*I set here the hard wall potential as a test.
	 * The hard wall is put at the bottom of the simulation box.*/
	
	for(loop=0;loop<VProd(Pts);loop++){
		if(loop%Pts.z == 0) Vext[loop] = 20.0;
		else Vext[loop] = 0.0; 
		}
	printf("Setting the external potential is done\n");
	return;
	}
	

void SetInitialDensity(){
	int loop;
	
	for(loop=0;loop<VProd(Pts);loop++){
	Density[loop] = exp(Beta*(Atom[0].miu - Vext[loop]));
		}
	printf("Setting the initial density is done\n");
	return;
	}
