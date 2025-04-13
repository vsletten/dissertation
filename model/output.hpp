/* output.h: functions for writing out various data files */

#ifndef output_h
#define output_h

#include "lattice.hpp"
#include "ucell.hpp"

class output {
public:
  static void initDatafile(void);
  static void writeData(Lattice *lattice, int n, float t);
  static void writeMSI(Lattice *lattice,  ucell::unitCell cell, const char *name, int bonds);
  static void writeSurf(Lattice *lattice, ucell::unitCell cell);
  static void writeXYZ(Lattice *lattice, ucell::unitCell cell);
};


#endif
