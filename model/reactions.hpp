/* reactions.h: functions for updating as result of reaction */

#ifndef reactions_h
#define reactions_h

#include "lattice.hpp"

class Reactions
{
public:
  static void adsorbAl(lattice::Lattice lattice, int site);
  static void adsorbSi(lattice::Lattice lattice, int site);
  static void desorbAl(lattice::Lattice lattice, int site);
  static void desorbSi(lattice::Lattice lattice, int site);
  static void diffuse(lattice::Lattice lattice, int source, int target);
  static void r0(lattice::Lattice lattice, int site);
  static void r1(lattice::Lattice lattice, int site);
  static void r2(lattice::Lattice lattice, int site);
  static void r3(lattice::Lattice lattice, int site);
  static void r4(lattice::Lattice lattice, int site);
  static void r5(lattice::Lattice lattice, int site);
  static void r6(lattice::Lattice lattice, int site);
  static void r7(lattice::Lattice lattice, int site);
  static void r8(lattice::Lattice lattice, int site);
  static void r9(lattice::Lattice lattice, int site);
  static void r10(lattice::Lattice lattice, int site);
  static void r11(lattice::Lattice lattice, int site);
  static void r12(lattice::Lattice lattice, int site);
  static void r13(lattice::Lattice lattice, int site);
  static void r14(lattice::Lattice lattice, int site);
  static void r15(lattice::Lattice lattice, int site);
};


#endif
