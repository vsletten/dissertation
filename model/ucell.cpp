/* ucell.c: functions for unit cell */
#include <stdlib.h>
#include "ucell.hpp"
#include "futil.hpp"


void UnitCell::DisposeUnitCell(UnitCell *unitCell)
{
  if (unitCell != nullptr) {
    delete[] unitCell->sites;
    delete unitCell;
  }
}


CellDimensions *UnitCell::GetCellDimensions(void)
{
  CellDimensions *cd;

  cd = new CellDimensions();
  cd->a = this->A; cd->b = this->B; cd->c = this->C;
  cd->alpha = this->Alpha; cd->beta = this->Beta; cd->gamma = this->Gamma;
  return cd;
}


int UnitCell::GetNumPositions(void)
{
  return this->Npos;
}

/* getNumNbrs: returns the number of neighbors that a site should
 *             have based on what kind it is */
int UnitCell::GetNumNeighbors(int state)
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
UnitCell *UnitCell::CreateUnitCell(void)
{
  int i, j;
  UnitCell *unitCell = new UnitCell();
  std::ifstream f;

  f = Futil::OpenInputFile("data.cell");
  f >> unitCell->A >> unitCell->B >> unitCell->C;
  Futil::EatComment(f, '#');
  f >> unitCell->Alpha >> unitCell->Beta >> unitCell->Gamma;
  Futil::EatComment(f, '#');
  f >> unitCell->Npos;
  unitCell->sites = new CellSite[unitCell->Npos + 1];
  Futil::EatComment(f, '#');
  for (i = 0; i < unitCell->Npos; i++) {
    f >> unitCell->sites[i].n >> unitCell->sites[i].state >> unitCell->sites[i].x >> unitCell->sites[i].y >> unitCell->sites[i].z;
    for (j = 0; j < 6; j++) {
      f >> unitCell->sites[i].nbr[j].n >> unitCell->sites[i].nbr[j].a >> unitCell->sites[i].nbr[j].b >> unitCell->sites[i].nbr[j].c;
    }
  }
  unitCell->sites[unitCell->Npos].n = -1;
  Futil::CloseFile(f);
  return unitCell;
}





