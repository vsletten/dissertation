/* ucell.c: functions for unit cell */
#include <stdlib.h>
#include <stdio.h>
//#include "common.hpp"
#include "ucell.hpp"
#include "futil.hpp"


static float Alpha, Beta, Gamma, A, B, C;  /* unit cell params */
static int Npos;


void ucell::free_cell(unitCell c)
{
  free(c);
}


ucell::cellDim *ucell::getCelldim(void)
{
  ucell::cellDim *cd;

  cd = (ucell::cellDim *) malloc (sizeof(ucell::cellDim));
  cd->a = A; cd->b = B; cd->c = C;
  cd->alpha = Alpha; cd->beta = Beta; cd->gamma = Gamma;
  return cd;
}


int ucell::getNpos(void)
{
  return Npos;
}

/* getNumNbrs: returns the number of neighbors that a site should
 *             have based on what kind it is */
int ucell::getNumNbrs(int state)
{
  int nbrs;

  switch (state / 100) {
    case 1: nbrs = 6; break;          /* Al */
    case 2: nbrs = 4; break;          /* Si */
    case 3: case 5: nbrs = 2; break;  /* Si-O-Si or Al-OH-Al type O */
    case 4: nbrs = 3; break;          /* Si-O-Al2 type O */
    default: nbrs = -1; break;
  }
  return nbrs;
}




/* readCell: read in unit cell parameters, and some simulation conditions */
ucell::unitCell ucell::readCell(void)
{
  int i, j;
  ucell::unitCell c;
  FILE *f;

  f = Futil::openFile("data.cell", "r");
  fscanf(f, " %f %f %f ", &A, &B, &C);
  Futil::eatComment(f, '#');
  fscanf(f, " %f %f %f ", &Alpha, &Beta, &Gamma);
  Futil::eatComment(f, '#');
  fscanf(f, " %d ", &Npos);
  Futil::eatComment(f, '#');
  c = (ucell::unitCell) malloc (sizeof(struct ucell::cellSite) * (Npos + 1));
  for (i = 0; i < Npos; i++) {
    fscanf(f, " %d %d %f %f %f ", &c[i].n, &c[i].state, &c[i].x, &c[i].y, 
	   &c[i].z);
    for (j = 0; j < 6; j++) {
      fscanf(f, "%d %d %d %d ", &c[i].nbr[j].n, &c[i].nbr[j].a, &c[i].nbr[j].b,
	     &c[i].nbr[j].c);
    }
  }
  c[Npos].n = -1;
  return c;
  Futil::closeFile(f);
}





