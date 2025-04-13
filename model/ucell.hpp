/* ucell.h: exports structs and routines for dealing with the unit cell */

#ifndef ucell_h
#define ucell_h

#include "sim.hpp"

class ucell {
public:
  /* neighborSite: Only for unit cell template; specifies a neighbor of
   *               a unit cell position. */
  struct neighborSite {
    int a, b, c, /* unit cell offsets */
        n;       /* unit cell site number */
  };

  /* cellSite: Uniquely identifies a unit cell position and specifies
   *           what type of site and how to find neighbors. */
  struct cellSite {
    float x, y, z;              /* coordinate of site in unit cell */
    int n,                      /* site number */
        state;                  /* initial state (which kind of site?) */
    struct neighborSite nbr[6]; /* neighbor list */
  };
  typedef struct cellSite *unitCell;

  struct celldim {
    float a, b, c, alpha, beta, gamma;
  };
  typedef struct celldim cellDim;

  static void free_cell(unitCell c);
  static cellDim *getCelldim(void);
  static int getNpos(void);
  static int getNumNbrs(int state);
  static unitCell readCell(void);
};

#endif
