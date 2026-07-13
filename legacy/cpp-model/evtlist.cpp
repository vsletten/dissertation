#include "evtlist.hpp"
#include "envrn.hpp"
#include "myerr.hpp"
#include "rxnlist.hpp"
#include <iostream>
#include <stdlib.h>

EventList *EventList::CreateEventList(Lattice *lattice, Environment *environment, Reaction *rxList) {
  EventList *eventList = nullptr;
  EventList *currEvent = nullptr;
  int s, i, env, lo, hi;

  /* build event list */
  for (s = 0; s < lattice->GetNsites(); s++) {
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
       /* matches reactant and is active? */
      if (lattice->sites[s].state ==
              rxList[i].reactant
          && environment->IsActive(s, i)) 
      {
        currEvent = new EventList();
        if (currEvent == nullptr)
        {
          Myerr::die("out of memory in CreateEventList");
        }
        currEvent->next = eventList;
        currEvent->site = s;
        currEvent->rxn = i;
        env = environment->CheckEnv(s);
        if (env < 0 || env >= rxList[i].nrates) 
        {
          std::cerr << "invalid environment: site " << s 
            << " state " << lattice->sites[s].state 
            << " env " << env 
            << " nrates " << rxList[i].nrates << std::endl;
          return nullptr;
        }
        currEvent->rate = rxList[i].rate[env];
        eventList = currEvent;
      }
    }
  }
  return eventList;
}

void EventList::DisposeEventList(EventList *&eventList) {
  EventList *current = eventList;
  EventList *next = nullptr;

  while (current != nullptr) {
    next = current->next;
    delete current;
    current = next;
  }
  eventList = nullptr;
}
