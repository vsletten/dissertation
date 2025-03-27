/* actions.h: exports functions for doing something */

#ifndef actions_h
#define actions_h

#include "lattice.h"
#include "evtlist.h"


float doEvent(Lattice lattice, eventList el);
void doReaction(Lattice lattice, int site, int rxn);


#endif
