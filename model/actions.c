//#include <stdio.h>
//#include "common.h"
#include "actions.h"
#include "reactions.h"
#include "ran2.h"
#include "myerr.h"
#include "math.h"


/* doEvent: randomly pick event and update state accordingly
 *          return time increment for the event */
float doEvent(Lattice lattice, eventList el)
{
  float ratesum, partsum, eps, dt;
  struct event *e;

  e = el;
  ratesum = 0;

  while (e->next) {
    ratesum += e->rate;
    e = e->next;
  }
  ratesum += e->rate;

  if (ratesum == 0.0)
    die("ratesum is zero!");
  dt = - log(ran2()) / ratesum;

  eps = ran2();
  partsum = 0;
  e = el;
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
void doReaction(Lattice lattice, int site, int rxn)
{
  switch (rxn) {
    case 0: r0(lattice, site); break;
    case 1: r1(lattice, site); break;
    case 2: r2(lattice, site); break;
    case 3: r3(lattice, site); break;
    case 4: r4(lattice, site); break;
    case 5: r5(lattice, site); break;
    case 6: r6(lattice, site);break;
    case 7: r7(lattice, site); break;
    case 8: r8(lattice, site); break;
    case 9: r9(lattice, site); break;
    case 10: r10(lattice, site); break;
    case 11: r11(lattice, site); break;
    case 12: r12(lattice, site); break;
    case 13: r13(lattice, site); break;
    case 14: r14(lattice, site); break;
    case 15: r15(lattice, site); break;
    case 16: adsorbAl(lattice, site); break;
    case 17: adsorbAl(lattice, site); break;
    case 18: adsorbSi(lattice, site); break;
    case 19: adsorbSi(lattice, site); break;
    case 20: desorbAl(lattice, site);
             clusters(lattice);
             break;
    case 21: desorbSi(lattice, site);
             clusters(lattice);
             break;
    case 22: desorbSi(lattice, site);
             clusters(lattice);
             break;
    case 23: desorbAl(lattice, site);
             clusters(lattice);
             break;
    default: diffuse(lattice, site, rxn); 
             clusters(lattice);
             break;
  }
}
