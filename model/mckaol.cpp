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
  Simulation *sc;
  ReactionList *reactionList;
  EventList *eventList;
  UnitCell *unitCell;
  Lattice *lattice;
  int i, nsites;
  float time = 0;
  char fname[50];

  sc = Simulation::CreateSimulation();
  initran2(&sc->ranseed);
  output::initDatafile();
  reactionList = ReactionList::CreateReactionList();
  unitCell = UnitCell::CreateUnitCell();

  lattice = Lattice::CreateLattice(*unitCell);
  nsites = lattice->GetNsites();
  lattice->FindPairs();
  lattice->PopulateSolid(reactionList->GetSiPotential(), reactionList->GetAlPotential());
  lattice->TerminateSurface();
  lattice->TerminateLattice();
  output::writeMSI(lattice, "start", sc->drawbonds);
  for (i = 0; i < sc->nsteps; i++) {
    eventList = EventList::CreateEventList(lattice, reactionList->GetReactions(), nsites);
    Actions actions(lattice);
    if (!eventList) {                  /* no events possible */
      output::writeMSI(lattice, "end", sc->drawbonds);
      if (sc->wsteps)
	  output::writeData(lattice, nsites, time);
      Myerr::die("out of events");
    }
    time += actions.DoEvent(eventList);
    EventList::DisposeEventList(eventList);            /* check if time to write stuff */
    if (sc->wsteps && (i % sc->wsteps == 0 || i == 0))
      output::writeData(lattice, nsites, time);
    if (sc->msteps && i && i % sc->msteps == 0) {
      sprintf(fname, "step%d", i);
      output::writeMSI(lattice, fname, sc->drawbonds);
    }
  }
  output::writeSurf(lattice);
  output::writeXYZ(lattice);
  output::writeMSI(lattice, "end", sc->drawbonds);
  ReactionList::DisposeReactionList(reactionList);
  lattice->TerminateLattice();
  UnitCell::DisposeUnitCell(unitCell);
  Simulation::DisposeSimulation(sc);
  return 0;
}
