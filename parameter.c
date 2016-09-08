#include "headdefs.h" 


 struct PropertyMol * Atom;
 
 int NumAtomType;
 
 real Temperature;
 
 real Beta;
 
 real Radius; // The radius of the hard sphere.

 const real Kb = 1.38064852E-26; //Boltzmann constant, unit in kJ / K

 const real Na = 6.022140857E23; //Avogradro constant, unit in mol^-1

 const real Eps = 8.854187817E-12; //The vacuum permittivity, unit in F/m (farads per metre)

 const real Q = 1.6021766208E-19; //The elementary charge, usually denoted as e or sometimes q, is the electric charge carried by a single proton, or equivalently, the magnitude of the electric charge carried by a single electron, which has charge -e. Unit in coulombs.

 const real h = 6.626070040E-37; //Planck constant  KJ * S, energy multiplied by time
 
 real Rcut; // cutoff for LJ attraction

 struct Vector Size;

 real * ThermWaveLength; //The thermal wavelenght "Lambda" = (beta * h * h / (2.0 * Pi * mass))

 struct VectorInt Pts; // Discretize the space into a given number of grid points, pts[0] is X dirextion, pts[1] is Y dirextion, pts[2] is Z dirextion

 struct Vextor Delta; //The distance between two grid nodes
 
 int AtomType; // 0 is HS, 1 is LJ
 
 real * Density;
 
 real_complex * F_Density; // The output arrays of the FFT transform of the density.
 
 real_complex * F_n0, * F_n1, * F_n2, * F_n3, * F_n1V_x, * F_n1V_y, * F_n1V_z, * F_n2V_x, * F_n2V_y, * F_n2V_z; // The output arrays of the FFT transform of the n0, n1, n2, n3, n1v, n2v.
 
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
 
 NameList nameList[NUM_PARAMETER] = {
  NameI (NumMol),
  NameR (deltaT),
  NameI (HeadType),
  NameI (EndType),
  NameI (StepDump),
  NameI (num_his_bars),
  NameR (MolLength),
  NameR (LengthX),
  NameR (LengthY),
};

