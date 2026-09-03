/* envrn.h: functions for checking the environment of a site */

#ifndef envrn_h
#define envrn_h

#include "lattice.h"

int checkEnv(Lattice lattice, int site);
int isActive(int site, Lattice lattice, int rxn);

#endif
