/* actions.h: exports functions for doing state transitions */

#ifndef actions_h
#define actions_h

#include "lattice.h"
#include "evtlist.h"


/* doEvent: randomly pick event and update state accordingly
 *          return time increment for the event */
float doEvent(Lattice lattice, eventList el);


/* doReaction: update state of site and its nbrs based on reaction rxn */
void doReaction(Lattice lattice, int site, int rxn);


#endif
