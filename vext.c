#include "headdefs.h"

// Set up the external potential
void SetVext(){
	
	int loop;
	int i,j,k;
	real rr;
	/*I set here the hard wall potential as a test.
	 * The hard wall is put at the bottom of the simulation box.*/
	 
/*  //Hard Wall	
	for(loop=0;loop<VProd(Pts);loop++){
		if(loop%Pts.z == 0) Vext[loop] = 100.0;
		else Vext[loop] = 0.0; 
   }
*/
    // LJ 12-6 Wall
	for(k=0;k<Pts.z;k++){
		rr = Sqr(k*Delta.z);
		for(j=0;j<Pts.y;j++){
			for(i=0;i<Pts.x;i++){
				if( rr > Sqr(Rcut)){
					Vext[k+Pts.z*(j+Pts.y*i)] = 0.0;
                  }
		       else if(rr <= Sqr(Rcut) && rr>0.0){
					Vext[k+Pts.z*(j+Pts.y*i)] = 4.0 * Atom[0].epslion * (pow(Atom[0].sigma, 12)/pow(rr, 6) - pow(Atom[0].sigma, 6)/pow(rr, 3));
					} 
			   else if(rr == 0.0){
					Vext[k+Pts.z*(j+Pts.y*i)] = 1000.0;
					}						
				}
			}
	}
	
	printf("Setting the external potential is done\n");
	return;
	}
	

void SetInitialDensity(){
	int loop;
	for(loop=0;loop<VProd(Pts);loop++){
        Density[loop] = 1.1*Atom[0].density;
		}
	printf("Setting the initial density is done\n");
	return;
	}

void SetLJ(){
	int i,j,k;
	real rr;
	
	for(k=0;k<Pts.z;k++){
		for(j=0;j<Pts.y;j++){
			for(i=0;i<Pts.x;i++){
			    rr = Sqr(i*Delta.x) + Sqr(j*Delta.y) + Sqr(k*Delta.z);
			    if(rr< Sqr(Atom[0].sigma) || rr > Sqr(Rcut)){
					Ulj[k+Pts.z*(j+Pts.y*i)] = 0.0;
					}
				else{
					Ulj[k+Pts.z*(j+Pts.y*i)] = 4.0 * Atom[0].epslion * (pow(Atom[0].sigma, 12)/pow(rr, 6) - pow(Atom[0].sigma, 6)/pow(rr, 3));
					} 	
				}
			}
		}
	return;
	}	
