#include <stdlib.h>
#include <stdio.h>
#include "evtlist.hpp"
#include "envrn.hpp"
#include "myerr.hpp"
#include "rxnlist.hpp"

evtlist::eventList evtlist::new_evtList(lattice::Lattice l, rxnlist::reactionList rl, int nsites)
{
  evtlist::eventList el = NULL;
  struct event *e;
  int s, i, env, lo, hi;
  char msg[100];

  for (s = 0; s < nsites; s++) {   /* build event list */
    if ((l[s].state % 100 == 0) && l[s].state > 200) {
      continue;
    } else if (l[s].state > 500) {
      lo = N400; hi = NHYD;
    } else if (l[s].state > 400) {
      lo = N300; hi = N400;
    } else if (l[s].state > 300) {
      lo = 0; hi = N300;
    } else {
      lo = NHYD; hi = NRXN;
    }
    for (i = lo; i < hi; i++) {
      if (l[s].state == rl[i].reactant   /* matches reactant and is active? */
	  && Environment::isActive(s, l, i)) {
	if (!(e = (struct event *) malloc (sizeof(struct event))))
	  Myerr::die("malloc failed");
	e->next = el;
	e->site = s;
	e->rxn = i;
	env = Environment::checkEnv(l, s);
	if (env < 0 || env >= rl[i].nrates) {
	  sprintf(msg, "invalid environment: site %d state %d env %d nrates %d"
                  , s, l[s].state, env, rl[i].nrates);
	  Myerr::die(msg);
	}
	e->rate = rl[i].rate[env];
	el = e;
      }
    }
  }
  return el;
}



void evtlist::free_evtList(evtlist::eventList el)
{
  struct event *e;

  while(el) {
    e = el->next;
    free(el);
    el = e;
  }
}




