#include <iostream>
#include <stdlib.h>
#include "myerr.hpp"

/* die: print msg to stderr then exit */
void Myerr::die(const char *msg)
{
  std::cerr << "mckaol: " << msg << std::endl;
  exit(1);
}



