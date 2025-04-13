#include <stdlib.h>
#include <stdio.h>
#include "futil.hpp"
#include "myerr.hpp"
#include "bfsearch.hpp"
#include "common.hpp"
#include "lattice.hpp"
#include "rxnlist.hpp"

static int Acells;    /* number of cells in a direction */
static int Bcells;    /* number of cells in b direction */
static int Surfplane; /* 0 for ac surface, 1 for bc surface */
static int Nsites;    /* total number of sites in simulation */



/* clusters: check for and remove any unattached clusters */
void lattice::clusters(lattice::Lattice lattice)
{
  int i;

  BFSearch::BFS(lattice, 0, Nsites);
  for (i = 0; i < Nsites; i++) 
  {
    if (lattice[i].color == WHITE &&
	      lattice[i].state != 9)
      lattice[i].state = (int) (lattice[i].state / 100) * 100;
    else if (lattice[i].color == GRAY) 
    {
      printf("bad color at %d  %d\n", i, lattice[i].state);
      Myerr::die("error in clusters");
    } 
  }
}



/* countNbrs: returns actual number of neighbors 
 *            (for lattice edge detection) */
int lattice::countNbrs(lattice::Lattice lattice, int site)
{
  int i;

  for (i = 0; i < 6 && lattice[site].nbr[i] >= 0; i++)
    ;
  return i;
}




/* findPairs: sets up the double bridge information */
void lattice::findPairs(lattice::Lattice lattice)
{
  int j, al1, al2, al3, al4, o1, o2, found;

  for (o1 = 0; o1 < Nsites; o1++) 
  {
    if (lattice[o1].state < 400 || lattice[o1].pair >= 0)
      continue;
    al1 = lattice[o1].nbr[0];
    al2 = lattice[o1].nbr[1];
    if (al1 == -1 || al2 == -1)
      continue;
    found = FALSE;
    for (j = 0; j < 6 && !found; j++) 
    {
      if ((o2 = lattice[al1].nbr[j]) < 0) /* no neighbor, must be edge */
	      found = TRUE;
      else 
      {
        al3 = lattice[o2].nbr[0];
        al4 = lattice[o2].nbr[1];
        if ((al1 == al3 || al1 == al4) &&
            (al2 == al3 || al2 == al4)) 
        {
          found = TRUE;
          lattice[o1].pair = o2;
          lattice[o2].pair = o1;
        }
      }
    }
  }
}



void lattice::free_lattice(lattice::Lattice l)
{
  free(l);
}


void lattice::getDim(int *a, int *b)
{
  *a = Acells;
  *b = Bcells;
}

/* getNeighbor: determine lattice index of nbr j for site i 
 *              depends upon index generation scheme in makeLattice */
int lattice::getNeighbor(lattice::Lattice l, ucell::unitCell c, int i, int j)
{
  int a, b, nbr, npos;

  npos = ucell::getNpos();
  a = l[i].a + c[l[i].n].nbr[j].a;
  b = l[i].b + c[l[i].n].nbr[j].b;
  if (c[l[i].n].nbr[j].n < 0) 
  {      /* no jth neighbor */
    nbr = c[l[i].n].nbr[j].n;
  } 
  else if (((a == Acells || a < 0) && Surfplane) ||
	        ((b == Bcells || b < 0) && !Surfplane)) 
  {
    nbr = -1;
  } 
  else 
  {
    if (a >= Acells )                   /* boundary conditions - a */
      a = 0;
    else if (a < 0)
      a = Acells - 1;
    if (b >= Bcells)                    /* boundary conditions - b */
      b = 0;
    else if (b < 0)
      b = Bcells - 1;
    nbr = a * Bcells * npos + b * npos + c[l[i].n].nbr[j].n;
  }
  return nbr;
}




int lattice::getNsites(void)
{
  return Nsites;
}



/* makeLattice: repeat unit cell to make the lattice and
 *              initialize state of sites 
 *              index generation scheme important for getNeighbor */
lattice::Lattice lattice::makeLattice(ucell::unitCell c)
{
  int a, b, n, i, j, npos;
  FILE *f;
  lattice::Lattice lattice;

  f = Futil::openFile("data.lattice","r");
  fscanf(f, " %d %d %d", &Acells, &Bcells, &Surfplane);
  Futil::closeFile(f);

  npos = ucell::getNpos();
  Nsites = Acells * Bcells * npos;
  lattice = (lattice::Lattice) malloc (sizeof(*lattice) * Nsites);
  for (a = 0; a < Acells; a++) 
  {
    for (b = 0; b < Bcells; b++) 
    {
      i = a * Bcells * npos + b * npos;
      for (n = 0; n < npos; n++) 
      {
        lattice[i].a = a;
        lattice[i].b = b;
        lattice[i].n = n;
        lattice[i].state = c[n].state;
        lattice[i].pair = -1;
        lattice[i].lostal = -1;
        for (j = 0; j < 6; j++) 
        {
          lattice[i].nbr[j] = getNeighbor(lattice, c, i, j);
        }
        i++;
      }
    }
  }
  return lattice;
}



/* populateSolid: update state of sites which start as part of the solid
 *                as determined by input parameters */
void lattice::populateSolid(lattice::Lattice lattice)
{
  int z, x, n, i, top, len, npos;
  float dms, dma, frac;
  char s[20];

  npos = ucell::getNpos();
  rxnlist::getChem(&dms, &dma);
  if (dms + dma > 0.5) 
  {
    frac = 0.3;
    sprintf(s, "supersaturated");
  } 
  else if (dms + dma < -0.5) 
  {
    frac = 0.7;
    sprintf(s, "undersaturated");
  } 
  else 
  {
    frac = 0.5;
    sprintf(s, "near equilibrium");
  }
  top = (Surfplane ? (Acells * frac) : (Bcells * frac));
  printf("Detected %s conditions -- filling %d unit cells\n", s, top);
  len = (Surfplane ? Bcells : Acells);
  for (z = 0; z < top; z++) 
  {
    for (x = 0; x < len; x++) 
    {  /* if bc, then x = b else x = a */
      i = (Surfplane ? (z * Bcells * npos + x * npos) : 
	   (x * Bcells * npos + z * npos));
      for (n = 0; n < npos; n++) 
      {
        lattice[i].state++;
        i++;
      }
    }
  }
}




/* terminateLattice: set to EDGE state of any atoms at very top or 
 *                   bottom of the lattice */
void lattice::terminateLattice(lattice::Lattice lattice)
{
  int top, len, x, i, n, npos;

  npos = ucell::getNpos();
  top = (Surfplane ? (Acells - 1) : (Bcells - 1));
  len = (Surfplane ? Bcells : Acells);
  for (x = 0; x < len; x++) 
  {  /* if bc, then x = b else x = a */
    i = (Surfplane ? 
        (top * Bcells * npos + x * npos) :
        (x * Bcells * npos + top * npos));
    for (n = 0; n < npos; n++, i++) 
      if (ucell::getNumNbrs(lattice[i].state) != countNbrs(lattice, i))
      	lattice[i].state = EDGE;
    i = (Surfplane ? (x * npos) :
	      (x * Bcells * npos));
    for (n = 0; n < npos; n++, i++)
      if (ucell::getNumNbrs(lattice[i].state) != countNbrs(lattice, i))
	      lattice[i].state = EDGE;
  }
}




/* terminateSurface: set state of O sites bonded to cations at surface 
 *                   so as to terminate dangling bonds with OH's */
void lattice::terminateSurface(lattice::Lattice lattice)
{
  int i, i2, j, type;

  /* "dangle" the O's which are missing cations */
  for (i = 0; i < Nsites; i++) 
  {
    if (lattice[i].state % 100 == 0 || lattice[i].state / 100 < 3 ||
      	lattice[i].state == EDGE)
      continue;
    for (j = 0; j < 6; j++) 
    {
      i2 = lattice[i].nbr[j];
      if (i2 >= 0 && (type = lattice[i2].state % 100) == 0) 
      {
        switch (lattice[i].state) 
        {
          case 401: lattice[i].state = (type == 2) ? 404 : 406; break;
          case 404: lattice[i].state = 409; break;
          case 406: lattice[i].state = (type == 2) ? 408 : 409; break;
          case 409: case 408: lattice[i].state = 400; break;
          case 301: lattice[i].state = 303; break;
          case 303: lattice[i].state = 300; break;
          case 501: lattice[i].state = 503; break;
          case 503: lattice[i].state = 500; break;
          default: break;
        }
      }
    }
  }

  /* turn on dangling oxygens to terminate dangling bonds */
  for (i = 0; i < Nsites; i++) 
  {
    if ((type = lattice[i].state / 100) > 2 || 
	      lattice[i].state % 100 != 1) 
      continue;
    for (j = 0; j < 6; j++) 
    {
      i2 = lattice[i].nbr[j];
      if (i2 >= 0 &&lattice[i2].state % 100 != 1)
      switch (lattice[i2].state) 
      {
        case 300: lattice[i2].state = 303; break;
        case 400: lattice[i2].state = ((type == 2) ? 408 : 409); break;
        case 408: lattice[i2].state = 406; break;
        case 409: lattice[i2].state = ((type == 2) ? 406 : 404); break;
        case 404: lattice[i2].state = 401; break;
        case 406: lattice[i2].state = 401; break;
        case 500: lattice[i2].state = 503; break;
        case 503: lattice[i2].state = 501; break;
        default: break;
      }
    }
  }

  /* update state of Si and Al to reflect terminal OH's */
  for (i = 0; i < Nsites; i++) 
  {  
    if (lattice[i].state % 100 == 0 ||  /* unoccupied */
	      lattice[i].state / 100 > 2 ||   /* O or OH */
	      lattice[i].state == EDGE)       /* edge */
      continue;
    if (lattice[i].state / 100 == 1) 
    {  /* Al */
      for (j = 0; j < 6; j++) 
      {
        i2 = lattice[i].nbr[j];
        if (i2 < 0)
          continue;
        switch (lattice[i2].state) 
        {
          case 503: case 409: lattice[i].state++; break;
          default: break;
        }
      }
    } 
    else 
    {
      for (j = 0; j < 6; j++) 
      {
        i2 = lattice[i].nbr[j];
        if (i2 < 0)
          continue;
        switch (lattice[i2].state) 
        {
          case 303: case 408: lattice[i].state++; break;
          default: break;
        }
      }
    }
  }
}





