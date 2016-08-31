#ifndef PARALLEL_H
#define PARALLEL_H

extern int myid;
extern int numprocs;

/* this subroutine initializes parallel processing*/
void InitMPI(int argc, char **argv);

/* this subroutine terminates parallel processing*/
void ExitMPI();

#endif /*PARALLEL_H*/

