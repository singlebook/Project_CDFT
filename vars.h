 /*
 * Created on Jul 12, 2016
 * Author: Wei Chen
 */

#pragma once
#ifndef VARS_H
#define	VARS_H

typedef std::complex<double> dcomplex;


// physical constants
extern const double BOLTZMANN_K; // units of kJ  / (mol* K) = (1.38064E-26) * 6.02214179E23
extern const double ELECTROSTATIC_K; // 1/(4*PI*vacuum_permittivity) units of (kJ*Angstrom)/(mol * elementary_charge^2) (see below for notes on obtaining this constant)

// system parameters
extern int Nx; // The number of grids on X direction 
extern int Ny; // The number of grids on Y direction
extern int Nz; // The number of grids on Z direction

extern int BoundaryCondition_X_Front; // 0 is Periodic Boundary, 1 is wall, 2 is bulk
extern int BoundaryCondition_X_Back; // 0 is Periodic Boundary, 1 is wall, 2 is bulk

extern int BoundaryCondition_Y_Left; // 0 is Periodic Boundary, 1 is wall, 2 is bulk
extern int BoundaryCondition_Y_Right; // 0 is Periodic Boundary, 1 is wall, 2 is bulk

extern int BoundaryCondition_Z_Up; // 0 is Periodic Boundary, 1 is wall, 2 is bulk
extern int BoundaryCondition_Z_Down; // 0 is Periodic Boundary, 1 is wall, 2 is bulk

extern double Temperature;
extern double Miu; // chemical potential
extern double Density_Bulk;
extern double Esplion; // the strength of the interactions between the nearest-neighbor sites

extern double * Vext;
