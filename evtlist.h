/* evtlist.h: exports structs and functions for the event list */

#ifndef evtlist_h
#define evtlist_h

#include "lattice.h"
#include "sim.h"
#include "rxnlist.h"

/* event: Uniquely identifies an event which may occur. */
struct event {
  int   site,     /* site where reaction happens */
        rxn;      /* reaction which happens
                   * for diffusion, this will contain OFFSET + index 
                   * of target for jump */
  float rate;     /* specific rate for this event (e.g., adsorption at
                   * a specific site) */
  struct event *next;
};

typedef struct event *eventList;

eventList new_evtList(Lattice l, reactionList rl, int nsites);
void free_evtList(eventList el);



#endif
