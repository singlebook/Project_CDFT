#ifndef NAMELIST_H
#define	NAMELIST_H

typedef enum {N_I, N_R} VType;

#define NameI(x)  {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x)  {#x, &x, N_R, sizeof (x) / sizeof (double)}

typedef struct {
  char *vName;
  void *vPtr;
  VType vType;
  int vLen, vStatus;
} NameList;

#define ValI(x)  {&x, N_I, sizeof (x) / sizeof (int)}
#define ValR(x)  {&x, N_R, sizeof (x) / sizeof (double)}

typedef struct {
  void *vPtr;
  VType vType;
  int vLen;
} ValList;

int GetNameList (int argc, char **argv);
void PrintNameList (FILE *fp);



#endif	/* NAMELIST_H */


