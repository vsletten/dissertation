/* reactions.h: functions for updating as result of reaction */

#ifndef reactions_h
#define reactions_h

#include "lattice.h"

void adsorbAl(Lattice lattice, int site);
void adsorbSi(Lattice lattice, int site);
void desorbAl(Lattice lattice, int site);
void desorbSi(Lattice lattice, int site);
void diffuse(Lattice lattice, int source, int target);
void r0(Lattice lattice, int site);
void r1(Lattice lattice, int site);
void r2(Lattice lattice, int site);
void r3(Lattice lattice, int site);
void r4(Lattice lattice, int site);
void r5(Lattice lattice, int site);
void r6(Lattice lattice, int site);
void r7(Lattice lattice, int site);
void r8(Lattice lattice, int site);
void r9(Lattice lattice, int site);
void r10(Lattice lattice, int site);
void r11(Lattice lattice, int site);
void r12(Lattice lattice, int site);
void r13(Lattice lattice, int site);
void r14(Lattice lattice, int site);
void r15(Lattice lattice, int site);


#endif
