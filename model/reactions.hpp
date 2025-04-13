/* reactions.h: functions for updating as result of reaction */

#ifndef reactions_h
#define reactions_h

#include "lattice.hpp"

class Reactions
{
public:
  static void adsorbAl(Lattice *lattice, int site);
  static void adsorbSi(Lattice *lattice, int site);
  static void desorbAl(Lattice *lattice, int site);
  static void desorbSi(Lattice *lattice, int site);
  static void diffuse(Lattice *lattice, int source, int target);
  static void r0(Lattice *lattice, int site);
  static void r1(Lattice *lattice, int site);
  static void r2(Lattice *lattice, int site);
  static void r3(Lattice *lattice, int site);
  static void r4(Lattice *lattice, int site);
  static void r5(Lattice *lattice, int site);
  static void r6(Lattice *lattice, int site);
  static void r7(Lattice *lattice, int site);
  static void r8(Lattice *lattice, int site);
  static void r9(Lattice *lattice, int site);
  static void r10(Lattice *lattice, int site);
  static void r11(Lattice *lattice, int site);
  static void r12(Lattice *lattice, int site);
  static void r13(Lattice *lattice, int site);
  static void r14(Lattice *lattice, int site);
  static void r15(Lattice *lattice, int site);
};


#endif
