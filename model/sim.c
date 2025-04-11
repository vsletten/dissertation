#include <stdio.h>
#include <stdlib.h>
#include "futil.h"
#include "myerr.h"
#include "sim.h"

/* readCond: read in some more simulation conditions */
Simulation readCond(void)
{
  FILE *f;
  Simulation s;

  f = openFile("data.sim","r");
  s = (Simulation) malloc (sizeof(*s));
  if (fscanf(f, " %d ", &s->nsteps) != 1)
    die("invalid number of steps in input file");
  eatComment(f, '#');
  if (fscanf(f, " %d ", &s->wsteps) != 1)
    die("invalid number of data write steps in input file");
  eatComment(f, '#');
  if (fscanf(f, " %d ", &s->msteps) != 1)
    die("invalid number of movie write steps in input file");
  if (s->msteps == -1)
    s->msteps = s->nsteps + 1;
  eatComment(f, '#');
  if (fscanf(f, " %ld ", &s->ranseed) != 1)
    die("invalid number of steps in input file");
  eatComment(f, '#');
  if (fscanf(f, " %d ", &s->drawbonds) != 1)
    die("invalid draw bonds? parameter in input file");

  return s;
  closeFile(f);
}


void free_Simulation(Simulation s)
{
  free(s);
}

