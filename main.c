#include "headdefs.h"

int main (int argc, char **argv){
	GetNameList (argc, argv);
	PrintNameList (stdout);
	Setup();
	SetVext();
	SetInitialDensity();
	Set_FFT();
	Free_memory();
	Clean_FFT();
	return 0;
}
