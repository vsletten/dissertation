/* actions.h: exports functions for doing something */

#ifndef actions_h
#define actions_h

#include "lattice.hpp"
#include "evtlist.hpp"


class Actions {
public:
    static float doEvent(lattice::Lattice lattice, evtlist::eventList el);
    static void doReaction(lattice::Lattice lattice, int site, int rxn);
};


#endif
