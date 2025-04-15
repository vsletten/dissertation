/* output.h: functions for writing out various data files */

#ifndef output_h
#define output_h

#include "lattice.hpp"

class output {
public:
  static void initDatafile(void);
  static void writeData(Lattice *lattice, int n, float t);
  static void writeMSI(Lattice *lattice, const char *name, int bonds);
  static void writeSurf(Lattice *lattice);
  static void writeXYZ(Lattice *lattice);
};


#endif
