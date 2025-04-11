/* Libraries */
#include <stdio.h>
#include "rxnlist.h"
#include "evtlist.h"
#include "lattice.h"
#include "ucell.h"
#include "myerr.h"
#include "actions.h"
#include "sim.h"
#include "output.h"
#include "ran2.h"

int main(int argc, char *argv[])
{
  Simulation sc;
  reactionList rl;
  eventList el;
  unitCell uc;
  Lattice l;
  int i, nsites;
  float time = 0;
  char fname[50];

  sc = readCond();
  initran2(&sc->ranseed);
  initDatafile();
  rl = readRxns();
  uc = readCell();

  l = makeLattice(uc);
  nsites = getNsites();
  findPairs(l);
  populateSolid(l);
  terminateSurface(l);
  terminateLattice(l);

  writeMSI(l, uc, "start", sc->drawbonds);
  for (i = 0; i < sc->nsteps; i++) {
    el = new_evtList(l, rl, nsites);
    if (!el) {                  /* no events possible */
      writeMSI(l, uc, "end", sc->drawbonds);
      if (sc->wsteps)
	writeData(l, nsites, time);
      die("out of events");
    }
    time += doEvent(l, el);
    free_evtList(el);            /* check if time to write stuff */
    if (sc->wsteps && (i % sc->wsteps == 0 || i == 0))
      writeData(l, nsites, time);
    if (sc->msteps && i && i % sc->msteps == 0) {
      sprintf(fname, "step%d", i);
      writeMSI(l, uc, fname, sc->drawbonds);
    }
  }
  writeSurf(l, uc);
  writeXYZ(l, uc);
  writeMSI(l, uc, "end", sc->drawbonds);
  free_rxnList(rl);
  free_lattice(l);
  free_cell(uc);
  free_Simulation(sc);
  return 0;
}
