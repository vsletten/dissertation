#include "lattice.hpp"
#include "bfsearch.hpp"
#include "common.hpp"
#include "futil.hpp"
#include "myerr.hpp"
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <string>

using namespace std;
/* clusters: check for and remove any unattached clusters */
void Lattice::RemoveUnattachedClusters() {
  int i;

  BreadthFirstSearch::ColorNodes(this, 0);
  for (i = 0; i < Num_Sites; i++) {
    if (this->sites[i].color == UNREACHABLE && this->sites[i].state != 9)
      this->sites[i].state = (int)(this->sites[i].state / 100) * 100;
    else if (this->sites[i].color == ENQUEUED) {
      std::cout << "bad color at " << i << "; state is " << this->sites[i].state << std::endl;
      Myerr::die("error in clusters");
    }
  }
}

/* countNbrs: returns actual number of neighbors
 *            (for lattice edge detection) */
int Lattice::CountNbrs(int site) {
  int i;

  for (i = 0; i < 6 && this->sites[site].nbr[i] >= 0; i++)
    ;
  return i;
}

/* findPairs: sets up the double bridge information */
void Lattice::FindPairs() {
  int j, al1, al2, al3, al4, o1, o2, found;

  for (o1 = 0; o1 < Num_Sites; o1++) {
    if (this->sites[o1].state < 400 || this->sites[o1].pair >= 0)
      continue;
    al1 = this->sites[o1].nbr[0];
    al2 = this->sites[o1].nbr[1];
    if (al1 == -1 || al2 == -1)
      continue;
    found = FALSE;
    for (j = 0; j < 6 && !found; j++) {
      if ((o2 = this->sites[al1].nbr[j]) < 0) /* no neighbor, must be edge */
        found = TRUE;
      else {
        al3 = this->sites[o2].nbr[0];
        al4 = this->sites[o2].nbr[1];
        if ((al1 == al3 || al1 == al4) && (al2 == al3 || al2 == al4)) {
          found = TRUE;
          this->sites[o1].pair = o2;
          this->sites[o2].pair = o1;
        }
      }
    }
  }
}

void Lattice::GetDim(int *a, int *b) {
  *a = Num_aCells;
  *b = Num_bCells;
}

/* getNeighbor: determine lattice index of nbr j for site i
 *              depends upon index generation scheme in makeLattice */
int Lattice::GetNeighbor(int i, int j) {
  int a, b, nbr, npos;

  npos = this->unitCell->GetNumPositions();
  a = this->sites[i].a + this->unitCell->GetCellSites()[this->sites[i].n].nbr[j].a;
  b = this->sites[i].b + this->unitCell->GetCellSites()[this->sites[i].n].nbr[j].b;
  if (this->unitCell->GetCellSites()[this->sites[i].n].nbr[j].n < 0) { /* no jth neighbor */
    nbr = this->unitCell->GetCellSites()[this->sites[i].n].nbr[j].n;
  } else if (((a == Num_aCells || a < 0) && SurfacePlane) ||
             ((b == Num_bCells || b < 0) && !SurfacePlane)) {
    nbr = -1;
  } else {
    if (a >= Num_aCells) /* boundary conditions - a */
      a = 0;
    else if (a < 0)
      a = Num_aCells - 1;
    if (b >= Num_bCells) /* boundary conditions - b */
      b = 0;
    else if (b < 0)
      b = Num_bCells - 1;
    nbr = a * Num_bCells * npos + b * npos + this->unitCell->GetCellSites()[this->sites[i].n].nbr[j].n;
  }
  return nbr;
}

int Lattice::GetNsites(void) { return this->Num_Sites; }

/* makeLattice: repeat unit cell to make the lattice and
 *              initialize state of sites
 *              index generation scheme important for getNeighbor */
Lattice *Lattice::CreateLattice(UnitCell &unitCell) {
  int a, b, n, i, j, npos;
  ifstream f;
  Lattice *lattice = new Lattice();
  lattice->unitCell = &unitCell;

  f = Futil::OpenInputFile("data.lattice");
  f >> lattice->Num_aCells >> lattice->Num_bCells >> lattice->SurfacePlane;
  Futil::CloseFile(f);

  npos = lattice->unitCell->GetNumPositions();
  lattice->Num_Sites = lattice->Num_aCells * lattice->Num_bCells * npos;
  lattice->sites = new LatticeSite[lattice->Num_Sites];
  for (a = 0; a < lattice->Num_aCells; a++) {
    for (b = 0; b < lattice->Num_bCells; b++) {
      i = a * lattice->Num_bCells * npos + b * npos;
      for (n = 0; n < npos; n++) {
        lattice->sites[i].a = a;
        lattice->sites[i].b = b;
        lattice->sites[i].n = n;
        lattice->sites[i].state = lattice->unitCell->GetCellSites()[n].state;
        lattice->sites[i].pair = -1;
        lattice->sites[i].lostal = -1;
        for (j = 0; j < 6; j++) {
          lattice->sites[i].nbr[j] = lattice->GetNeighbor(i, j);
        }
        i++;
      }
    }
  }
  return lattice;
}

/* populateSolid: update state of sites which start as part of the solid
 *                as determined by input parameters */
void Lattice::PopulateSolid(float dmSi, float dmAl) {
  int z, x, n, i, top, len, npos;
  float frac;
  char s[20];

  npos = this->unitCell->GetNumPositions();
  if (dmSi + dmAl > 0.5) {
    frac = 0.3;
    sprintf(s, "supersaturated");
  } else if (dmSi + dmAl < -0.5) {
    frac = 0.7;
    sprintf(s, "undersaturated");
  } else {
    frac = 0.5;
    sprintf(s, "near equilibrium");
  }
  top = (SurfacePlane ? (Num_aCells * frac) : (Num_bCells * frac));
  printf("Detected %s conditions -- filling %d unit cells\n", s, top);
  len = (SurfacePlane ? Num_bCells : Num_aCells);
  for (z = 0; z < top; z++) {
    for (x = 0; x < len; x++) { /* if bc, then x = b else x = a */
      i = (SurfacePlane ? (z * Num_bCells * npos + x * npos)
                        : (x * Num_bCells * npos + z * npos));
      for (n = 0; n < npos; n++) {
        this->sites[i].state++;
        i++;
      }
    }
  }
}

/* terminateLattice: set to EDGE state of any atoms at very top or
 *                   bottom of the lattice */
void Lattice::TerminateLattice() {
  int top, len, x, i, n, npos;

  npos = this->unitCell->GetNumPositions();
  top = (SurfacePlane ? (Num_aCells - 1) : (Num_bCells - 1));
  len = (SurfacePlane ? Num_bCells : Num_aCells);
  for (x = 0; x < len; x++) { /* if bc, then x = b else x = a */
    i = (SurfacePlane ? (top * Num_bCells * npos + x * npos)
                      : (x * Num_bCells * npos + top * npos));
    for (n = 0; n < npos; n++, i++)
      if (this->unitCell->GetNumNeighbors(this->sites[i].state) != this->CountNbrs(i))
        this->sites[i].state = EDGE;
    i = (SurfacePlane ? (x * npos) : (x * Num_bCells * npos));
    for (n = 0; n < npos; n++, i++)
      if (this->unitCell->GetNumNeighbors(this->sites[i].state) != this->CountNbrs(i))
        this->sites[i].state = EDGE;
  }
}

/* terminateSurface: set state of O sites bonded to cations at surface
 *                   so as to terminate dangling bonds with OH's */
void Lattice::TerminateSurface() {
  int i, i2, j, type;

  /* "dangle" the O's which are missing cations */
  for (i = 0; i < this->Num_Sites; i++) {
    if (this->sites[i].state % 100 == 0 || this->sites[i].state / 100 < 3 ||
        this->sites[i].state == EDGE)
      continue;
    for (j = 0; j < 6; j++) {
      i2 = this->sites[i].nbr[j];
      if (i2 >= 0 && (type = this->sites[i2].state % 100) == 0) {
        switch (this->sites[i].state) {
        case 401:
          this->sites[i].state = (type == 2) ? 404 : 406;
          break;
        case 404:
          this->sites[i].state = 409;
          break;
        case 406:
          this->sites[i].state = (type == 2) ? 408 : 409;
          break;
        case 409:
        case 408:
          this->sites[i].state = 400;
          break;
        case 301:
          this->sites[i].state = 303;
          break;
        case 303:
          this->sites[i].state = 300;
          break;
        case 501:
          this->sites[i].state = 503;
          break;
        case 503:
          this->sites[i].state = 500;
          break;
        default:
          break;
        }
      }
    }
  }

  /* turn on dangling oxygens to terminate dangling bonds */
  for (i = 0; i < this->Num_Sites; i++) {
    if ((type = this->sites[i].state / 100) > 2 ||
        this->sites[i].state % 100 != 1)
      continue;
    for (j = 0; j < 6; j++) {
      i2 = this->sites[i].nbr[j];
      if (i2 >= 0 && this->sites[i2].state % 100 != 1)
        switch (this->sites[i2].state) {
        case 300:
          this->sites[i2].state = 303;
          break;
        case 400:
          this->sites[i2].state = ((type == 2) ? 408 : 409);
          break;
        case 408:
          this->sites[i2].state = 406;
          break;
        case 409:
          this->sites[i2].state = ((type == 2) ? 406 : 404);
          break;
        case 404:
          this->sites[i2].state = 401;
          break;
        case 406:
          this->sites[i2].state = 401;
          break;
        case 500:
          this->sites[i2].state = 503;
          break;
        case 503:
          this->sites[i2].state = 501;
          break;
        default:
          break;
        }
    }
  }

  /* update state of Si and Al to reflect terminal OH's */
  for (i = 0; i < this->Num_Sites; i++) {
    if (this->sites[i].state % 100 == 0 || /* unoccupied */
        this->sites[i].state / 100 > 2 ||  /* O or OH */
        this->sites[i].state == EDGE)      /* edge */
      continue;
    if (this->sites[i].state / 100 == 1) { /* Al */
      for (j = 0; j < 6; j++) {
        i2 = this->sites[i].nbr[j];
        if (i2 < 0)
          continue;
        switch (this->sites[i2].state) {
        case 503:
        case 409:
          this->sites[i].state++;
          break;
        default:
          break;
        }
      }
    } else {
      for (j = 0; j < 6; j++) {
        i2 = this->sites[i].nbr[j];
        if (i2 < 0)
          continue;
        switch (this->sites[i2].state) {
        case 303:
        case 408:
          this->sites[i].state++;
          break;
        default:
          break;
        }
      }
    }
  }
}
