//#include <stdio.h>
//#include "common.hpp"
#include "actions.hpp"
#include "reactions.hpp"
#include "lattice.hpp"
#include "ran2.hpp"
#include "myerr.hpp"
#include <cmath>


/* doEvent: randomly pick event and update state accordingly
 *          return time increment for the event */
float Actions::doEvent(Lattice *lattice, EventList *eventList)
{
  float ratesum, partsum, eps, dt;
  EventList *e;

  e = eventList;
  ratesum = 0;

  while (e->next) {
    ratesum += e->rate;
    e = e->next;
  }
  ratesum += e->rate;

  if (ratesum == 0.0)
    Myerr::die("ratesum is zero!");
  dt = - log(ran2()) / ratesum;

  eps = ran2();
  partsum = 0;
  e = eventList;
  while (e->next) {
    partsum += e->rate / ratesum;
    if (eps <= partsum)
      break;
    else
      e = e->next;
  }

  doReaction(lattice, e->site, e->rxn);
  return dt;
}



/* doReaction: update state of site and its nbrs based on reaction rxn */
void Actions::doReaction(Lattice *lattice, int site, int rxn)
{
  switch (rxn) {
    case 0: Reactions::r0(lattice, site); break;
    case 1: Reactions::r1(lattice, site); break;
    case 2: Reactions::r2(lattice, site); break;
    case 3: Reactions::r3(lattice, site); break;
    case 4: Reactions::r4(lattice, site); break;
    case 5: Reactions::r5(lattice, site); break;
    case 6: Reactions::r6(lattice, site);break;
    case 7: Reactions::r7(lattice, site); break;
    case 8: Reactions::r8(lattice, site); break;
    case 9: Reactions::r9(lattice, site); break;
    case 10: Reactions::r10(lattice, site); break;
    case 11: Reactions::r11(lattice, site); break;
    case 12: Reactions::r12(lattice, site); break;
    case 13: Reactions::r13(lattice, site); break;
    case 14: Reactions::r14(lattice, site); break;
    case 15: Reactions::r15(lattice, site); break;
    case 16: Reactions::adsorbAl(lattice, site); break;
    case 17: Reactions::adsorbAl(lattice, site); break;
    case 18: Reactions::adsorbSi(lattice, site); break;
    case 19: Reactions::adsorbSi(lattice, site); break;
    case 20: Reactions::desorbAl(lattice, site); break;
    case 21: Reactions::desorbSi(lattice, site); break;
    case 22: Reactions::desorbSi(lattice, site); break;
    case 23: Reactions::desorbAl(lattice, site); break;
             lattice->RemoveUnattachedClusters();
             break;
    default: Reactions::diffuse(lattice, site, rxn); 
             lattice->RemoveUnattachedClusters();
             break;
  }
}
