#include <stdlib.h>
#include "futil.hpp"
#include "myerr.hpp"
#include "sim.hpp"

/* readCond: read in some more simulation conditions */
Simulation *Simulation::CreateSimulation(void)
{
  std::ifstream f;
  Simulation *s;

  f = Futil::OpenInputFile("data.sim");
  s = new Simulation();
  if (!(f >> s->nsteps))
    Myerr::die("invalid number of steps in input file");
  Futil::EatComment(f, '#');
  if (!(f >> s->wsteps))
    Myerr::die("invalid number of data write steps in input file");
  Futil::EatComment(f, '#');
  if (!(f >> s->msteps))
    Myerr::die("invalid number of movie write steps in input file");
  if (s->msteps == -1)
    s->msteps = s->nsteps + 1;
  Futil::EatComment(f, '#');
  if (!(f >> s->drawbonds))
    Myerr::die("invalid number of steps in input file");
  Futil::EatComment(f, '#');
  if (!(f >> s->drawbonds))
    Myerr::die("invalid draw bonds? parameter in input file");

  Futil::CloseFile(f);
  return s;
}


void Simulation::DisposeSimulation(Simulation *&s)  
{
  if (s == nullptr) {
    return;
  }
  
  delete s;
  s = nullptr;
}

