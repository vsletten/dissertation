/* output.h: functions for writing out various data files */

#ifndef output_h
#define output_h

#include "lattice.hpp"
#include "ucell.hpp"

class output {
public:
  static void initDatafile(void);
  static void writeData(lattice::Lattice lattice, int n, float t);
  static void writeMSI(lattice::Lattice l,  ucell::unitCell cell, const char *name, int bonds);
  static void writeSurf(lattice::Lattice l, ucell::unitCell cell);
  static void writeXYZ(lattice::Lattice l, ucell::unitCell cell);
};


#endif
