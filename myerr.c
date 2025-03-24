#include <stdio.h>
#include <stdlib.h>
//#include "common.h"
#include "myerr.h"

/* die: print msg to stderr then exit */
void die(char *msg)
{
  fprintf(stderr, "mckaol: %s\n", msg);
  exit(1);
}



