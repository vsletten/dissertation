/* lattice.h: exports struct for a lattice site and functions 
 *            for changing or obtaining info about a lattice */

#ifndef lattice_h
#define lattice_h

#include "ucell.h"
#include "sim.h"

#define EDGE 9
#define WRONG 99
#define AL     1
#define SI     2
#define SIOSI  3
#define SIOAL  4
#define ALOHAL 5
#define OFFSET 100


/* macros */
#define ISOCC(A)       ((A).state % 100 > 0)
#define EXISTS(A)      ((A) >= 0)
#define ISEDGE(A)      ((A).state == EDGE)
#define TYPE(A)        ((A).state / 100)
#define ISSI(A)        (TYPE(A) == SI)
#define ISAL(A)        (TYPE(A) == AL)
#define ISOX(A)        ((A).state > 299)
#define WRONGCATION(A) ((A).state % 100 == WRONG)


/* latticeSite: Uniquely identifies a lattice site and describes
 *              its current state. */
struct latticeSite {
  int a, b,       /* unit cell coordinates (for single sheet) */
      n,          /* site within unit cell */
      state,      /* current state */
      color,      /* for BFS for clusters, WHITE, GRAY, BLACK */
      pair,       /* for 400s and 500s, other part of double bridge */
      lostal,     /* for 400s, when lose one al which it is */
      nbr[6];     /* neighbor list */
};
typedef struct latticeSite *Lattice;


void clusters(Lattice lattice);
int countNbrs(Lattice l, int s);
void findPairs(Lattice l);
void free_lattice(Lattice l);
int getNeighbor(Lattice l, unitCell c, int i, int j);
int getNsites(void);
Lattice makeLattice(unitCell uc);
void populateSolid(Lattice lattice);
void terminateLattice(Lattice lattice);
void terminateSurface(Lattice lattice);
void getDim(int *a, int *b);

#endif
