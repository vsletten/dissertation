#include <stdio.h>
#include <stdlib.h>
#include "myerr.hpp"

/* die: print msg to stderr then exit */
void Myerr::die(const char *msg)
{
  fprintf(stderr, "mckaol: %s\n", msg);
  exit(1);
}



