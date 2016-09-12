#include "headdefs.h" 

real Omega0(real K, real R){
	return sin(K*R) / (K*R);
}

real Omega1(real K, real R){
	return sin(K*R) / K;
}

real Omega2(real K, real R){
	return 4.0 * M_PI * R * sin(K*R)/K;
}

real Omega3(real K, real R){
	return (4.0 * M_PI / Cube(K)) * (sin(K*R) - K*R*cos(K*R));
}
	
struct VectorComplex Omega1V(struct Vector K, real R){
 real kk;
 struct VectorComplex t; 

 kk = VLen(K);
 VSCopy(t, -1.0 * I *  Omega3(kk, R) / (4.0 * M_PI * R), K);
 
 return t;
}


struct VectorComplex Omega2V(struct Vector K, real R){
  struct VectorComplex t; 
  VSCopy(t, 4.0 * M_PI * R, Omega1V(K, R));
  return t;
}
