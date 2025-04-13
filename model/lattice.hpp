/* lattice.h: exports struct for a lattice site and functions
 *            for changing or obtaining info about a lattice */

#ifndef lattice_h
#define lattice_h

#include "sim.hpp"
#include "ucell.hpp"

#define EDGE 9
#define WRONG 99
#define AL 1
#define SI 2
#define SIOSI 3
#define SIOAL 4
#define ALOHAL 5
#define OFFSET 100

/* macros */
#define ISOCC(A) ((A).state % 100 > 0)
#define EXISTS(A) ((A) >= 0)
#define ISEDGE(A) ((A).state == EDGE)
#define TYPE(A) ((A).state / 100)
#define ISSI(A) (TYPE(A) == SI)
#define ISAL(A) (TYPE(A) == AL)
#define ISOX(A) ((A).state > 299)
#define WRONGCATION(A) ((A).state % 100 == WRONG)

class lattice {
public:
  /* latticeSite: Uniquely identifies a lattice site and describes
   *              its current state. */
  struct latticeSite {
    int a, b,   /* unit cell coordinates (for single sheet) */
        n,      /* site within unit cell */
        state,  /* current state */
        color,  /* for BFS for clusters, WHITE, GRAY, BLACK */
        pair,   /* for 400s and 500s, other part of double bridge */
        lostal, /* for 400s, when lose one al which it is */
        nbr[6]; /* neighbor list */
  };
  typedef struct latticeSite *Lattice;

  static void clusters(Lattice lattice);
  static int countNbrs(Lattice l, int s);
  static void findPairs(Lattice l);
  static void free_lattice(Lattice l);
  static int getNeighbor(Lattice l, ucell::unitCell c, int i, int j);
  static int getNsites(void);
  static Lattice makeLattice(ucell::unitCell uc);
  static void populateSolid(Lattice lattice);
  static void terminateLattice(Lattice lattice);
  static void terminateSurface(Lattice lattice);
  static void getDim(int *a, int *b);
};

#endif
