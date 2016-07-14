#include "common.h"


// defin the grand potential
template<typename T>
class GrandPotential : public cppoptlib::Problem<T> {
  public:
    using typename cppoptlib::Problem<T>::TVector;
    // this is just the objective (NOT optional)
    T value(const TVector &density) {
        long index;
        T term = 0;
        T term1;
        T term2;
        int i, j, k;
        int up, down;
        int left, right;
        int front, back;
        long length = Nx * Ny * Nz;
        for(index = 0; index < length; index ++){
			term1 = BOLTZMANN_K * Temperature * (density[index] * log(density[index]) + (1.0 - density[index]) * log(1.0 - density[index])) + (Vext[index] - miu) * density[index];
            // index = i + Nx * (j + Ny * k)
            i = index%Nx;
            j = (index/Nx)%Ny;
            k = index/Nx/Ny;
            
            front = i + 1;
            if(front > Nx -1 && BoundaryCondition_X_Front == 0){
				front -= Nx;
				}
			back = i - 1;
            if(back < 0 && BoundaryCondition_X_Back == 0){
				back += Nx;
				}
			
			right = j + 1;
			if(right > Ny -1 && BoundaryCondition_Y_Left == 0){
				right -= Ny;
				}
			left = j -1;
			if(left < 0 && BoundaryCondition_Y_Right == 0){
				left += Ny;
				}
			
			up = k + 1;
			if(up > Nz -1 && BoundaryCondition_Z_Up == 0){
				up -= Nz;
				}
			down = k -1;
			if(down < 0 && BoundaryCondition_Z_Down == 0) {
				down += Nz;
				}
			
			if(BoundaryCondition_Z_Up == 0 && BoundaryCondition_Z_Down == 0){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
				}
				
			else if(BoundaryCondition_Z_Up == 1 && BoundaryCondition_Z_Down == 1){
				if(k != 0 && k != Nz - 1){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == 0){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] ) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] \\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == Nz -1){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)]  + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				}
			
			else if	(BoundaryCondition_Z_Up == 2 && BoundaryCondition_Z_Down == 1){
				if(k != 0 && k != Nz - 1){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == 0){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] \\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] \\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == Nz -1){
				term2 = -0.5 * density[index] * ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + Density_Bulk + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (4.0 * Density_Bulk + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )				
					}				
				}
			term += term1 + term2;	
			}
        return term;
    };

    // if you calculated the derivative by hand
    // you can implement it here (OPTIONAL)
    // otherwise it will fall back to (bad) numerical finite differences
    void gradient(const TVector &density, TVector &grad) {
        T term1;
        T term2;
        int i, j, k;
        int up, down;
        int left, right;
        int front, back;
        long length = Nx * Ny * Nz;
		for(index = 0; index < length; index ++){ 
			term1 = BOLTZMANN_K * Temperature * log(density[index] / (1.0 - density[index])) + 	Vext[index] - miu;
            // index = i + Nx * (j + Ny * k)
            i = index%Nx;
            j = (index/Nx)%Ny;
            k = index/Nx/Ny;
            
            front = i + 1;
            if(front > Nx -1 && BoundaryCondition_X_Front == 0){
				front -= Nx;
				}
			back = i - 1;
            if(back < 0 && BoundaryCondition_X_Back == 0){
				back += Nx;
				}
			
			right = j + 1;
			if(right > Ny -1 && BoundaryCondition_Y_Left == 0){
				right -= Ny;
				}
			left = j -1;
			if(left < 0 && BoundaryCondition_Y_Right == 0){
				left += Ny;
				}
			
			up = k + 1;
			if(up > Nz -1 && BoundaryCondition_Z_Up == 0){
				up -= Nz;
				}
			down = k -1;
			if(down < 0 && BoundaryCondition_Z_Down == 0) {
				down += Nz;
				}
			
			if(BoundaryCondition_Z_Up == 0 && BoundaryCondition_Z_Down == 0){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
				}
				
			else if(BoundaryCondition_Z_Up == 1 && BoundaryCondition_Z_Down == 1){
				if(k != 0 && k != Nz - 1){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == 0){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] ) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] \\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == Nz -1){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)]  + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				}
			
			else if	(BoundaryCondition_Z_Up == 2 && BoundaryCondition_Z_Down == 1){
				if(k != 0 && k != Nz - 1){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)] + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == 0){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + density[i + Nx * (j + Ny * up)]) \\
				                     + 0.25 * Esplion * (density[front + Nx * (j + Ny * up)] + density[back + Nx * (j + Ny * up)] \\
				                     + density[i + Nx * (left + Ny * up)] + density[i + Nx * (right + Ny * up)] \\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )
					}
				else if(k == Nz -1){
				term2 = -0.5 *  ( \\
				           Esplion * (density[front + Nx * (j + Ny * k)] + density[back + Nx * (j + Ny * k)]  + density[i + Nx * (left + Ny * k)] + density[i + Nx * (right + Ny * k)] + Density_Bulk + density[i + Nx * (j + Ny * down)]) \\
				                     + 0.25 * Esplion * (4.0 * Density_Bulk + density[front + Nx * (j + Ny * down)] + density[back + Nx * (j + Ny * down)]\\
				                     + density[i + Nx * (left + Ny * down)] + density[i + Nx * (right + Ny * down)]\\				                    
				                     + density[front + Nx * (right + Ny * k)] + density[front + Nx * (left + Ny * k)]  + density[back + Nx * (left + Ny * k)] + density[back + Nx * (right + Ny * k)] \\				                       
				                     ) 
				       )				
					}		               
			}
		grad[index] = term1 + term2;
		}
	};
};
