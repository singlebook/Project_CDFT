#ifndef EXCESS_FREE_ENERGY_H
#define EXCESS_FREE_ENERGY_H
real F_MFA(real rho);
real F_HS(real rho, real sigma);
real P1(real sigma, real rho_0);
real P2(real sigma, real epslion, real rho, char name[20]);
#endif
