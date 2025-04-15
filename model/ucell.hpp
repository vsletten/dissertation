/* ucell.h: exports structs and routines for dealing with the unit cell */

#ifndef ucell_h
#define ucell_h

/* neighborSite: Only for unit cell template; specifies a neighbor of
 *               a unit cell position. */
struct NeighborSite {
  int a, b, c, /* unit cell offsets */
      n;       /* unit cell site number */
};

/* cellSite: Uniquely identifies a unit cell position and specifies
 *           what type of site and how to find neighbors. */
struct CellSite {
  float x, y, z;              /* coordinate of site in unit cell */
  int n,                      /* site number */
      state;                  /* initial state (which kind of site?) */
  struct NeighborSite nbr[6]; /* neighbor list */
};

struct CellDimensions {
  float a, b, c, alpha, beta, gamma;
};


class UnitCell {
private:
  CellSite *sites = nullptr;
  float Alpha, Beta, Gamma, A, B, C;  /* unit cell params */
  int Npos;

public:
  static UnitCell *CreateUnitCell(void);
  static void DisposeUnitCell(UnitCell *uc);

  CellDimensions *GetCellDimensions(void);
  int GetNumPositions(void);
  int GetNumNeighbors(int state);
  CellSite *GetCellSites(void) { return this->sites; }
};

#endif
