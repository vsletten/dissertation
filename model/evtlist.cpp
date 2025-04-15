#include "evtlist.hpp"
#include "envrn.hpp"
#include "myerr.hpp"
#include "rxnlist.hpp"
#include <stdio.h>
#include <stdlib.h>

EventList *EventList::CreateEventList(Lattice *lattice, Reaction *rxList,
                           int nSites) {
  EventList *eventList = nullptr;
  EventList *currEvent = nullptr;
  int s, i, env, lo, hi;
  char msg[100];

  for (s = 0; s < nSites; s++) { /* build event list */
    if ((lattice->sites[s].state % 100 == 0) && lattice->sites[s].state > 200) {
      continue;
    } else if (lattice->sites[s].state > 500) {
      lo = N400;
      hi = NHYD;
    } else if (lattice->sites[s].state > 400) {
      lo = N300;
      hi = N400;
    } else if (lattice->sites[s].state > 300) {
      lo = 0;
      hi = N300;
    } else {
      lo = NHYD;
      hi = NRXN;
    }
    for (i = lo; i < hi; i++) {
      if (lattice->sites[s].state == rxList[i].reactant /* matches reactant and is active? */
          && Environment::isActive(s, lattice, i)) {
        currEvent = new EventList();
        if (currEvent == nullptr)
          Myerr::die("out of memory in CreateEventList");
        currEvent->next = eventList;
        currEvent->site = s;
        currEvent->rxn = i;
        env = Environment::checkEnv(lattice, s);
        if (env < 0 || env >= rxList[i].nrates) {
          sprintf(msg, "invalid environment: site %d state %d env %d nrates %d",
                  s, lattice->sites[s].state, env, rxList[i].nrates);
          Myerr::die(msg);
        }
        currEvent->rate = rxList[i].rate[env];
        eventList = currEvent;
      }
    }
  }
  return eventList;
}

void EventList::DisposeEventList(EventList *eventList) {
  EventList *current = eventList;
  EventList *next = nullptr;  

  while (current != nullptr) {
    next = current->next;
    delete current;
    current = next;
  }
}
