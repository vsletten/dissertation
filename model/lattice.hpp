/* lattice.h: exports struct for a lattice site and functions
 *            for changing or obtaining info about a lattice */

#ifndef lattice_h
#define lattice_h

#include "ucell.hpp"
#include "common.hpp"

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

class LatticeSite {
public:
  /* latticeSite: Uniquely identifies a lattice site and describes
   *              its current state. */
  int a, b,   /* unit cell coordinates (for single sheet) */
      n,      /* site within unit cell */
      state,  /* current state */
      pair,   /* for 400s and 500s, other part of double bridge */
      lostal, /* for 400s, when lose one al which it is */
      nbr[6]; /* neighbor list */
  Color color;/* for BFS for clusters, WHITE, GRAY, BLACK */

};

class Lattice {
private:
  int Num_aCells;   /* number of cells in a direction */
  int Num_bCells;   /* number of cells in b direction */
  int Num_Sites;    /* total number of sites in simulation */
  int SurfacePlane; /* 0 for ac surface, 1 for bc surface */
  UnitCell *unitCell;
  Lattice()
      : Num_aCells(0), Num_bCells(0), Num_Sites(0),
        SurfacePlane(0), unitCell(nullptr) {}

public:
  static Lattice *CreateLattice(UnitCell &uc);
  static void DisposeLoattice();

  LatticeSite *sites = nullptr;

  ~Lattice()
  {
    if (this->sites != nullptr) { 
      delete[] this->sites;
      this->sites = nullptr;
    }
  }
  
  void RemoveUnattachedClusters();
  UnitCell *GetUnitCell(void) { return this->unitCell; }  
  int CountNbrs(int s);
  void FindPairs();
  int GetNeighbor(int i, int j);
  int GetNsites(void);
  void PopulateSolid(float dmSi, float dmAl);
  void TerminateLattice();
  void TerminateSurface();
  void GetDim(int *a, int *b);
};

#endif
