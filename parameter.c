#include "headdefs.h" 


 struct PropertyAtom * Atom;
 
 int NumAtomType;
 
 real Temperature;
 
 real Beta;
 
 real Radius; // The radius of the hard sphere.

 const real Kb = 1.0 ;// 1.38064852E-26; //Boltzmann constant, unit in kJ / K

 const real Na = 6.022140857E23; //Avogradro constant, unit in mol^-1

 const real Eps = 8.854187817E-12; //The vacuum permittivity, unit in F/m (farads per metre)

 const real Q = 1.6021766208E-19; //The elementary charge, usually denoted as e or sometimes q, is the electric charge carried by a single proton, or equivalently, the magnitude of the electric charge carried by a single electron, which has charge -e. Unit in coulombs.

 const real h = 6.626070040E-37; //Planck constant  KJ * S, energy multiplied by time
 
 real Rcut; // cutoff for LJ attraction

 struct Vector Size;

 real * ThermWaveLength; //The thermal wavelenght "Lambda" = (beta * h * h / (2.0 * Pi * mass))

 struct VectorInt Pts; // Discretize the space into a given number of grid points, pts[0] is X dirextion, pts[1] is Y dirextion, pts[2] is Z dirextion

 struct Vector Delta; //The distance between two grid nodes
 
 int AtomType; // 0 is HS, 1 is LJ
 
 real * Density;
 
 real * F_Ulj;
 
 real Alpha; // Alpha is the line search parameter in the Picard iteration.
 
 real Stop;
 
 real_complex * F_Density; // The output arrays of the FFT transform of the density.
 
 real * n0, * n1, * n2, * n3, * n1V_x, * n1V_y, * n1V_z, * n2V_x, * n2V_y, * n2V_z;
 
 real_complex * F_n0, * F_n1, * F_n2, * F_n3, * F_n1V_x, * F_n1V_y, * F_n1V_z, * F_n2V_x, * F_n2V_y, * F_n2V_z; // The output arrays of the FFT transform of the n0, n1, n2, n3, n1v, n2v.
 
 real * Phi0, * Phi1, * Phi2, * Phi3, * Phi1V_x, * Phi1V_y, * Phi1V_z, * Phi2V_x, * Phi2V_y, * Phi2V_z; // The excess free energy density
 
 real_complex * F_Phi0, * F_Phi1, * F_Phi2, * F_Phi3, * F_Phi1V_x, * F_Phi1V_y, * F_Phi1V_z, * F_Phi2V_x, * F_Phi2V_y, * F_Phi2V_z; // The excess free energy density
 
 real * Miu_ex;
 
 real_complex * F_Miu_ex;
 
 real * Vext;

 int HS_size_approximation; 
 /*
	approximation for estimating effective HS size
   0 BH :   Baker-Handerson method
   1 WCA:   Weitheim-Chandler-Anderson method
   2 OCB:   Oxtoby
	*/

 int Approx_excess_free_energy;

/*
   	approximation for attractive excess helmholtz free energy
   0 MSA:  Mean Spherical Approximation
   1 FMSA: First order Mean Spherical Approximation
   2 MFA: Mean field approximation

    reasonable combination:
    BH  + MFA
    BH  + FMSA
 */
 


