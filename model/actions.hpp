/* actions.h: exports functions for doing something */

#ifndef actions_h
#define actions_h

#include "lattice.hpp"
#include "evtlist.hpp"

class Actions {
public:
    static float doEvent(Lattice *lattice, EventList *eventList);
    static void doReaction(Lattice *lattice, int site, int rxn);
};


#endif
