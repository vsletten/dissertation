/* output.h: functions for writing out various data files */

#ifndef output_h
#define output_h

#include "lattice.h"
#include "ucell.h"

void initDatafile(void);
void writeData(Lattice lattice, int n, float t);
void writeMSI(Lattice l,  unitCell cell, char *name, int bonds);
void writeSurf(Lattice l, unitCell cell);
void writeXYZ(Lattice l, unitCell cell);


#endif
