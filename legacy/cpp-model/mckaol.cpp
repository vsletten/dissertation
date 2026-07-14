/* Libraries */
#include "actions.hpp"
#include "evtlist.hpp"
#include "lattice.hpp"
#include "output.hpp"
#include "ran2.hpp"
#include "rxnlist.hpp"
#include "sim.hpp"
#include "ucell.hpp"
#include <sstream>
#include <iostream>

int main() {
  Simulation *sc;
  ReactionList *reactionList;
  EventList *eventList;
  UnitCell *unitCell;
  Lattice *lattice;
  int i;
  float time = 0;

  /* Initialize simulation */
  sc = Simulation::CreateSimulation();
  initran2(&sc->ranseed);
  output::initDatafile();

  reactionList = ReactionList::CreateReactionList();
  if (!reactionList) {
    std::cerr << "failed to create reaction list. Exiting." << std::endl;
    return 1;
  }
  
  unitCell = UnitCell::CreateUnitCell();
  if (!unitCell) {
    std::cerr << "failed to create unit cell. Exiting." << std::endl;
    return 1;
  }

  /* Create lattice */
  lattice = Lattice::CreateLattice(*unitCell);
  if (!lattice) {
    std::cerr << "failed to create lattice. Exiting." << std::endl;
    return 1;
  }
  lattice->FindPairs();
  lattice->PopulateSolid(reactionList->GetSiPotential(),
                         reactionList->GetAlPotential());
  lattice->TerminateSurface();
  lattice->TerminateLattice();  

  /* Initialize actions */
  Actions actions(lattice);

  /* Initialize environment */
  Environment *environment = new Environment(lattice);

  /* Write initial state */
  output::writeMSI(lattice, "start", sc->drawbonds);

  /* Run simulation */
  for (i = 0; i < sc->nsteps; i++) 
  {
    /* Create list of possible events for each simulation step */
    eventList = EventList::CreateEventList(
        lattice, environment, reactionList->GetReactions());

    /* Exit if no events possible or termination condition met */
    if (!eventList) 
    { 
      std::cerr << "failed to create event list. Exiting after " << i << " steps." << std::endl;
      break;
    }

    /* Perform event */
    float dt = actions.DoEvent(eventList);
    if (dt < 0.0)
    {
      std::cerr << "failed to do event. Exiting after " << i << " steps." << std::endl;
      break;
    }
    time += dt;

    /* Write data if necessary */
    if (sc->wsteps && (i % sc->wsteps == 0 || i == 0))
    {
      std::stringstream fname;
      fname << "step" << i << ".dat";
      output::writeData(fname.str().c_str(), lattice, lattice->GetNsites(), time);
    }
    
    if (sc->msteps && i && i % sc->msteps == 0) 
    {
      std::stringstream fname;
      fname << "step" << i;
      output::writeMSI(lattice, fname.str().c_str(), sc->drawbonds);
    }

    EventList::DisposeEventList(eventList);
  }
  output::writeData("end.dat", lattice, lattice->GetNsites(), time);
  output::writeSurf(lattice);
  output::writeXYZ(lattice);
  output::writeMSI(lattice, "end", sc->drawbonds);

  /* Clean up */
  delete environment;
  if (eventList != nullptr) {
    EventList::DisposeEventList(eventList);
  }
  ReactionList::DisposeReactionList(reactionList);
  UnitCell::DisposeUnitCell(unitCell);
  Lattice::DisposeLattice(lattice);
  Simulation::DisposeSimulation(sc);
  return 0;
}
