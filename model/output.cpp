#include "output.hpp"
#include "futil.hpp"
#include "myerr.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

void output::initDatafile(void) { std::remove("results.dat"); }

void output::writeData(const char *filename, Lattice *lattice, int n, float t) {
  std::ofstream f;
  int i, nsi[5], nal[7], oh, type, sitot, altot;

  f = Futil::OpenOutputFile(filename);
  sitot = altot = 0;
  for (i = 0; i < 5; i++)
    nsi[i] = 0;
  for (i = 0; i < 7; i++)
    nal[i] = 0;

  for (i = 0; i < n; i++) {
    if (lattice->sites[i].state % 100 == 0) /*skip empty sites */
      continue;
    if ((type = lattice->sites[i].state / 100) == 1) {
      oh = lattice->sites[i].state - 101;
      nal[oh]++;
      altot++;
    } else if (type == 2) {
      oh = lattice->sites[i].state - 201;
      nsi[oh]++;
      sitot++;
    }
  }

  f << t;
  f << "," << sitot;
  for (i = 0; i < 5; i++)
    f << "," << nsi[i];
  f << "," << altot;
  for (i = 0; i < 7; i++)
    f << "," << nal[i];
  f << "\n";

  Futil::CloseFile(f);
}

void output::writeXYZ(Lattice *lattice) {
  int idxSite, nAtoms, numAcells, numBcells;
  float x, y, z, al, bl, cl;
  std::ofstream xyzFile;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&numAcells, &numBcells);
  /* open file and print header */
  xyzFile = Futil::OpenOutputFile("start.xyz");

  nAtoms = lattice->GetNsites();

  for (idxSite = 0; idxSite < nAtoms; idxSite++) {
    if ((lattice->sites[idxSite].state > 200 && lattice->sites[idxSite].state < 206) ||
        (lattice->sites[idxSite].state > 100 && lattice->sites[idxSite].state < 108) ||
        (lattice->sites[idxSite].state > 300 && lattice->sites[idxSite].state < 400) ||
        (lattice->sites[idxSite].state > 400 && lattice->sites[idxSite].state < 500) ||
        (lattice->sites[idxSite].state > 500 && lattice->sites[idxSite].state < 600)) {
      x = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].x / al + lattice->sites[idxSite].a) * cd[0][0] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].y / bl + lattice->sites[idxSite].b) * cd[0][1] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].z / cl) * cd[0][2];
      y = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].y / bl + lattice->sites[idxSite].b) * cd[1][1] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].z / cl) * cd[1][2];
      z = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[idxSite].n].z / cl) * cd[2][2];
      if (lattice->sites[idxSite].state / 100 == 1)
        xyzFile << 13 << "," << x << "," << y << "," << z << "\n";
      else if (lattice->sites[idxSite].state / 100 == 2)
        xyzFile << 14 << "," << x << "," << y << "," << z << "\n";
      else
        xyzFile << 8 << "," << x << "," << y << "," << z << "\n";
    }
  }
  Futil::CloseFile(xyzFile);
}

void output::writeSurf(Lattice *lattice) {
  int i, j, nbr, flag, natom, acells, bcells;
  float x, y, al, bl, cl;
  // float z; // for 3D
  std::ofstream siSurfaceFile, alSurfaceFile;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&acells, &bcells);
  /* open file and print header */
  siSurfaceFile = Futil::OpenOutputFile("surfSi.out");
  alSurfaceFile = Futil::OpenOutputFile("surfAl.out");

  natom = lattice->GetNsites();

  for (i = 0; i < natom; i++) {
    if ((lattice->sites[i].state > 201 && lattice->sites[i].state < 206) ||
        (lattice->sites[i].state > 101 && lattice->sites[i].state < 108)) {
      flag = 0;
      for (j = 0; j < 6 && lattice->sites[i].nbr[j] >= 0; j++) {
        nbr = lattice->sites[i].nbr[j];
        if (lattice->sites[nbr].state == 303 || (lattice->sites[nbr].state > 403 && lattice->sites[nbr].state < 410) ||
            (lattice->sites[nbr].state == 503))
          flag = 1;
      }
      if (flag == 1) {
        x = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].x / al + lattice->sites[i].a) * cd[0][0] +
            (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[0][1] +
            (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].z / cl) * cd[0][2];
        y = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[1][1] +
            (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].z / cl) * cd[1][2];
        // z = (cell[lattice->sites[i].n].z / cl) * cd[2][2]; // for 3D
        if (lattice->sites[i].state > 200) {
          siSurfaceFile << x << "," << y << "\n";
        } else {
          alSurfaceFile << x << "," << y << "\n";
        }
      }
    }
  }

  Futil::CloseFile(alSurfaceFile);
  Futil::CloseFile(siSurfaceFile);
}

void output::writeMSI(Lattice *lattice, const char *name, int bonds) {
  int i, j, natom, nthing = 1, *id, i2, acells, bcells;
  std::stringstream fname;
  float x, y, z, al, bl, cl;
  std::ofstream fxyz;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&acells, &bcells);
  /* open file and print header */
  fname << name << ".msi";
  fxyz = Futil::OpenOutputFile(fname.str().c_str());
  fxyz << "# MSI CERIUS2 DataModel File Version 3 5\n";
  fxyz << "(1 Model\n";
  fxyz << "(A I Id 1)\n";
  fxyz << " (A C Label \"" << name << "\")\n";

  natom = lattice->GetNsites();
  if (!(id = (int *)malloc(sizeof(int) * natom)))
    Myerr::die("malloc failed");

  for (i = 0; i < natom; i++) {
    if (lattice->sites[i].state % 100 > 0 && lattice->sites[i].state != EDGE) {
      fxyz << " (" << ++nthing << " Atom\n";
      id[i] = nthing;
      fxyz << "  (A C ACL ";
      switch (lattice->sites[i].state / 100) {
      case 0:
        fxyz << "\"10 Ne\")\n";
        break;
      case 1:
        fxyz << "\"13 Al\")\n";
        break;
      case 2:
        fxyz << "\"14 Si\")\n";
        break;
      case 3:
      case 4:
        if (lattice->sites[i].state % 100 > 1)
          fxyz << "\"1 H\")\n";
        else
          fxyz << "\"8 O\")\n";
        break;
      case 5:
        if (lattice->sites[i].state % 100 > 1)
          fxyz << "\"1 H\")\n";
        else
          fxyz << "\"9 F\")\n";
        break;
      default:
        break;
      }
      x = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].x / al + lattice->sites[i].a) * cd[0][0] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[0][1] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].z / cl) * cd[0][2];
      y = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[1][1] +
          (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].z / cl) * cd[1][2];
      z = (lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].z / cl) * cd[2][2];
      fxyz << "  (A D XYZ (" << x << " " << y << " " << z << "))\n";
      fxyz << "  (A I Id " << nthing - 1 << ")\n";
      fxyz << "  (A C Label \"";
      switch (lattice->sites[i].state / 100) {
      case 0:
        fxyz << "Edge";
        break;
      case 1:
        fxyz << "Al";
        break;
      case 2:
        fxyz << "Si";
        break;
      case 3:
      case 4:
        fxyz << "O";
        break;
      case 5:
        fxyz << "OH";
        break;
      default:
        break;
      }
      fxyz << i << "\")\n";
      fxyz << "  (A I LabelType 0)\n";
      fxyz << " )\n";
    }
  }

  for (i = 0; i < natom; i++) {
    if (bonds && (lattice->sites[i].state % 100 > 0 && lattice->sites[i].state != EDGE)) {
      for (j = 0; j < 6; j++) {
        if ((i2 = lattice->sites[i].nbr[j]) >= 0 && i2 < i && lattice->sites[i2].state % 100 != 0 &&
            lattice->sites[i2].state != EDGE &&
            !(lattice->sites[i].a == 0 && lattice->sites[i2].a == acells - 1 &&
              lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].nbr[j].a == -1) &&
            !(lattice->sites[i].a == acells - 1 && lattice->sites[i2].a == 0 &&
              lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].nbr[j].a == 1) &&
            !(lattice->sites[i].b == 0 && lattice->sites[i2].b == bcells - 1 &&
              lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].nbr[j].b == -1) &&
            !(lattice->sites[i].b == bcells - 1 && lattice->sites[i2].b == 0 &&
              lattice->GetUnitCell()->GetCellSites()[lattice->sites[i].n].nbr[j].b == 1)) {
          fxyz << " (" << ++nthing << " Bond\n";
          fxyz << "  (A O Atom1 " << id[i] << ")\n";
          fxyz << "  (A O Atom2 " << id[i2] << ")\n";
          fxyz << " )\n";
        }
      }
    }
  }

  fxyz << ")\n";
  Futil::CloseFile(fxyz);

  free(id);
}
