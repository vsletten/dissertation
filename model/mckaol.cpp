/* Libraries */
#include <stdio.h>
#include "rxnlist.hpp"
#include "evtlist.hpp"
#include "lattice.hpp"
#include "ucell.hpp"
#include "myerr.hpp"
#include "actions.hpp"
#include "sim.hpp"
#include "output.hpp"
#include "ran2.hpp"

int main(int argc, char *argv[])
{
  sim::Simulation sc;
  rxnlist::reactionList rl;
  evtlist::eventList el;
  ucell::unitCell uc;
  lattice::Lattice l;
  int i, nsites;
  float time = 0;
  char fname[50];

  sc = sim::readCond();
  initran2(&sc->ranseed);
  output::initDatafile();
  rl = rxnlist::readRxns();
  uc = ucell::readCell();

  l = lattice::makeLattice(uc);
  nsites = lattice::getNsites();
  lattice::findPairs(l);
  lattice::populateSolid(l);
  lattice::terminateSurface(l);
  lattice::terminateLattice(l);

  output::writeMSI(l, uc, "start", sc->drawbonds);
  for (i = 0; i < sc->nsteps; i++) {
    el = evtlist::new_evtList(l, rl, nsites);
    if (!el) {                  /* no events possible */
      output::writeMSI(l, uc, "end", sc->drawbonds);
      if (sc->wsteps)
	output::writeData(l, nsites, time);
      Myerr::die("out of events");
    }
    time += Actions::doEvent(l, el);
    evtlist::free_evtList(el);            /* check if time to write stuff */
    if (sc->wsteps && (i % sc->wsteps == 0 || i == 0))
      output::writeData(l, nsites, time);
    if (sc->msteps && i && i % sc->msteps == 0) {
      sprintf(fname, "step%d", i);
      output::writeMSI(l, uc, fname, sc->drawbonds);
    }
  }
  output::writeSurf(l, uc);
  output::writeXYZ(l, uc);
  output::writeMSI(l, uc, "end", sc->drawbonds);
  rxnlist::free_rxnList(rl);
  lattice::free_lattice(l);
  ucell::free_cell(uc);
  sim::free_Simulation(sc);
  return 0;
}
