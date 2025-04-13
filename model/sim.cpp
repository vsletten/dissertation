#include <stdio.h>
#include <stdlib.h>
#include "futil.hpp"
#include "myerr.hpp"
#include "sim.hpp"

/* readCond: read in some more simulation conditions */
sim::Simulation sim::readCond(void)
{
  FILE *f;
  Simulation s;

  f = Futil::openFile("data.sim","r");
  s = (Simulation) malloc (sizeof(*s));
  if (fscanf(f, " %d ", &s->nsteps) != 1)
    Myerr::die("invalid number of steps in input file");
  Futil::eatComment(f, '#');
  if (fscanf(f, " %d ", &s->wsteps) != 1)
    Myerr::die("invalid number of data write steps in input file");
  Futil::eatComment(f, '#');
  if (fscanf(f, " %d ", &s->msteps) != 1)
    Myerr::die("invalid number of movie write steps in input file");
  if (s->msteps == -1)
    s->msteps = s->nsteps + 1;
  Futil::eatComment(f, '#');
  if (fscanf(f, " %ld ", &s->ranseed) != 1)
    Myerr::die("invalid number of steps in input file");
  Futil::eatComment(f, '#');
  if (fscanf(f, " %d ", &s->drawbonds) != 1)
    Myerr::die("invalid draw bonds? parameter in input file");

  return s;
  Futil::closeFile(f);
}


void sim::free_Simulation(Simulation s)
{
  free(s);
}

