/* evtlist.h: exports structs and functions for the event list */

#ifndef evtlist_h
#define evtlist_h

#include "lattice.hpp"
#include "rxnlist.hpp"
#include "sim.hpp"

class evtlist {
public:
  /* event: Uniquely identifies an event which may occur. */
  struct event {
    int site,   /* site where reaction happens */
        rxn;    /* reaction which happens
                 * for diffusion, this will contain OFFSET + index
                 * of target for jump */
    float rate; /* specific rate for this event (e.g., adsorption at
                 * a specific site) */
    struct event *next;
  };

  typedef struct event *eventList;

  static evtlist::eventList new_evtList(lattice::Lattice l, rxnlist::reactionList rl, int nsites);
  static void free_evtList(evtlist::eventList el);
};

#endif
