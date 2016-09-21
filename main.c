#include "headdefs.h"

int main (int argc, char **argv){
	GetNameList (argc, argv);
	PrintNameList (stdout);
	Setup();
	Iteration();
	Cal_g();
	Free_memory();
	return 0;
}
