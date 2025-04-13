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

int main()
{
  sim::Simulation sc;
  rxnlist::reactionList rl;
  EventList *el;
  ucell::unitCell uc;
  Lattice *lattice;
  int i, nsites;
  float time = 0;
  char fname[50];

  sc = sim::readCond();
  initran2(&sc->ranseed);
  output::initDatafile();
  rl = rxnlist::readRxns();
  uc = ucell::readCell();

  lattice = Lattice::CreateLattice(uc);
  nsites = lattice->GetNsites();
  lattice->FindPairs();
  lattice->PopulateSolid();
  lattice->TerminateSurface();
  lattice->TerminateLattice();

  output::writeMSI(lattice, uc, "start", sc->drawbonds);
  for (i = 0; i < sc->nsteps; i++) {
    el = EventList::CreateEventList(lattice, rl, nsites);
    if (!el) {                  /* no events possible */
      output::writeMSI(lattice, uc, "end", sc->drawbonds);
      if (sc->wsteps)
	  output::writeData(lattice, nsites, time);
      Myerr::die("out of events");
    }
    time += Actions::doEvent(lattice, el);
    EventList::DisposeEventList(el);            /* check if time to write stuff */
    if (sc->wsteps && (i % sc->wsteps == 0 || i == 0))
      output::writeData(lattice, nsites, time);
    if (sc->msteps && i && i % sc->msteps == 0) {
      sprintf(fname, "step%d", i);
      output::writeMSI(lattice, uc, fname, sc->drawbonds);
    }
  }
  output::writeSurf(lattice, uc);
  output::writeXYZ(lattice, uc);
  output::writeMSI(lattice, uc, "end", sc->drawbonds);
  rxnlist::free_rxnList(rl);
  lattice->TerminateLattice();
  ucell::free_cell(uc);
  sim::free_Simulation(sc);
  return 0;
}
