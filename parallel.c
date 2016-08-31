#include "headdefs.h"

int myid;
int numprocs;


void InitMPI(int argc, char **argv) {
   int ierror, flag;
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (!flag) {
/* initialize the MPI execution environment */
      ierror = MPI_Init(&argc,&argv);
      if (ierror) exit(1);
   }
/* determine the rank of the calling process in the communicator */
   ierror = MPI_Comm_rank(MPI_COMM_WORLD,&myid);
/* determine the size of the group associated with a communicator */
   ierror = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   return;
}

void ExitMPI() {
   int ierror, flag;
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (flag) {
/* synchronize processes */
      ierror = MPI_Barrier(MPI_COMM_WORLD);
/* terminate MPI execution environment */
      ierror = MPI_Finalize();
   }
   return;
}
