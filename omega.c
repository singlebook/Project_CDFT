#include "headdefs.h" 

real Omega0(real K, real R){
	if(K != 0)
	return sin(K*R) / (K*R);
	else 
	return 1.0;
}

real Omega1(real K, real R){
	if(K != 0)
	return sin(K*R) / K;
	else
	return R;
}

real Omega2(real K, real R){
	if(K != 0)
	return 4.0 * M_PI * R * sin(K*R)/K;
	else
	return 4.0 * M_PI * Sqr(R);
}

real Omega3(real K, real R){
	if(K != 0)
	return (4.0 * M_PI / Cube(K)) * (sin(K*R) - K*R*cos(K*R));
	else
	return 4.0 * M_PI * Cube(R) / 3.0 ;
}
	
struct VectorComplex Omega1V(struct Vector K, real R){
 real kk;
 struct VectorComplex t; 

 kk = VLen(K);
 if (kk != 0)
 VSCopy(t, -1.0 * I *  Omega3(kk, R) / (4.0 * M_PI * R), K);
 else
 VSCopy(t, 0.0, K);
 
 return t;
}


struct VectorComplex Omega2V(struct Vector K, real R){
  struct VectorComplex t; 
  VSCopy(t, 4.0 * M_PI * R, Omega1V(K, R));
  return t;
}
