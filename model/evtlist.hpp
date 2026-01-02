/* evtlist.h: exports structs and functions for the event list */

#ifndef evtlist_h
#define evtlist_h

#include "envrn.hpp"
#include "lattice.hpp"
#include "rxnlist.hpp"

  /* event: Uniquely identifies an event which may occur. */
class Event {
public:
  int site,   /* site where reaction happens */
      rxn;    /* reaction which happens
               * for diffusion, this will contain OFFSET + index
               * of target for jump */
  float rate; /* specific rate for this event (e.g., adsorption at
               * a specific site) */
};

class EventList : public Event {
public:
  EventList *next;

  static EventList *CreateEventList(Lattice *lattice, Environment *environment, Reaction *rxList);

  static void DisposeEventList(EventList *&eventList);
};

#endif
