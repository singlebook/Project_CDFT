#include "common.h"

const double BOLTZMANN_K = 0.008314463;
const double ELECTROSTATIC_K = 1389.354325379097;

int Nx;
int Ny;
int Nz;

int BoundaryCondition_X_Front; 
int BoundaryCondition_X_Back; 

int BoundaryCondition_Y_Left; 
int BoundaryCondition_Y_Right; 

int BoundaryCondition_Z_Up; 
int BoundaryCondition_Z_Down;

double Temperature;
double Miu;
double Density_Bulk;
double Esplion; 

double * Vext;
