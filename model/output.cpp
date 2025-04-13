#include "output.hpp"
#include "futil.hpp"
#include "myerr.hpp"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void output::initDatafile(void) { unlink("results.dat"); }

void output::writeData(Lattice *lattice, int n, float t) {
  FILE *f;
  int i, nsi[5], nal[7], oh, type, sitot, altot;

  f = Futil::openFile("results.dat", "a");
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

  fprintf(f, "%f", t);
  fprintf(f, ",%d", sitot);
  for (i = 0; i < 5; i++)
    fprintf(f, ",%d", nsi[i]);
  fprintf(f, ",%d", altot);
  for (i = 0; i < 7; i++)
    fprintf(f, ",%d", nal[i]);
  fprintf(f, "\n");

  Futil::closeFile(f);
}

void output::writeXYZ(Lattice *lattice, ucell::unitCell cell) {
  int i, j, natom, acells, bcells;
  float x, y, z, al, bl, cl;
  FILE *xyz;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&acells, &bcells);
  /* open file and print header */
  xyz = Futil::openFile("start.xyz", "w");

  natom = lattice->GetNsites();

  for (i = 0; i < natom; i++) {
    if ((lattice->sites[i].state > 200 && lattice->sites[i].state < 206) ||
        (lattice->sites[i].state > 100 && lattice->sites[i].state < 108) ||
        (lattice->sites[i].state > 300 && lattice->sites[i].state < 400) ||
        (lattice->sites[i].state > 400 && lattice->sites[i].state < 500) ||
        (lattice->sites[i].state > 500 && lattice->sites[i].state < 600)) {
      x = (cell[lattice->sites[i].n].x / al + lattice->sites[i].a) * cd[0][0] +
          (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[0][1] +
          (cell[lattice->sites[i].n].z / cl) * cd[0][2];
      y = (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[1][1] +
          (cell[lattice->sites[i].n].z / cl) * cd[1][2];
      z = (cell[lattice->sites[i].n].z / cl) * cd[2][2];
      if (lattice->sites[i].state / 100 == 1)
        fprintf(xyz, "%d,%f,%f,%f\n", 13, x, y, z);
      else if (lattice->sites[i].state / 100 == 2)
        fprintf(xyz, "%d,%f,%f,%f\n", 14, x, y, z);
      else
        fprintf(xyz, "%d,%f,%f,%f\n", 8, x, y, z);
    }
  }
  Futil::closeFile(xyz);
}

void output::writeSurf(Lattice *lattice, ucell::unitCell cell) {
  int i, j, nbr, flag, natom, acells, bcells;
  float x, y, z, al, bl, cl;
  FILE *fsi, *fal;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&acells, &bcells);
  /* open file and print header */
  fsi = Futil::openFile("surfSi.out", "w");
  fal = Futil::openFile("surfAl.out", "w");

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
        x = (cell[lattice->sites[i].n].x / al + lattice->sites[i].a) * cd[0][0] +
            (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[0][1] +
            (cell[lattice->sites[i].n].z / cl) * cd[0][2];
        y = (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[1][1] +
            (cell[lattice->sites[i].n].z / cl) * cd[1][2];
        z = (cell[lattice->sites[i].n].z / cl) * cd[2][2];
        if (lattice->sites[i].state > 200) {
          fprintf(fsi, "%f,%f\n", x, y);
        } else {
          fprintf(fal, "%f,%f\n", x, y);
        }
      }
    }
  }

  Futil::closeFile(fal);
  Futil::closeFile(fsi);
}

void output::writeMSI(Lattice *lattice, ucell::unitCell cell, const char *name,
                      int bonds) {
  int i, j, natom, nthing = 1, *id, i2, acells, bcells;
  char fname[100];
  float x, y, z, al, bl, cl;
  FILE *fxyz;
  float cd[][3] = {{4.9725, -0.0262374, -1.3362},
                   {0.0, 8.92893, -0.30084},
                   {0.0, 0.0, 7.384}};
  al = sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1] + cd[0][2] * cd[0][2]);
  bl = sqrt(cd[1][1] * cd[1][1] + cd[1][2] * cd[1][2]);
  cl = cd[2][2];

  lattice->GetDim(&acells, &bcells);
  /* open file and print header */
  sprintf(fname, "%s.msi", name);
  fxyz = Futil::openFile(fname, "w");
  fprintf(fxyz, "# MSI CERIUS2 DataModel File Version 3 5\n");
  fprintf(fxyz, "(1 Model\n");
  fprintf(fxyz, "(A I Id 1)\n");
  fprintf(fxyz, " (A C Label \"%s\")\n", name);

  natom = lattice->GetNsites();
  if (!(id = (int *)malloc(sizeof(int) * natom)))
    Myerr::die("malloc failed");

  for (i = 0; i < natom; i++) {
    if (lattice->sites[i].state % 100 > 0 && lattice->sites[i].state != EDGE) {
      fprintf(fxyz, " (%d Atom\n", ++nthing);
      id[i] = nthing;
      fprintf(fxyz, "  (A C ACL ");
      switch (lattice->sites[i].state / 100) {
      case 0:
        fprintf(fxyz, "\"10 Ne\")\n");
        break;
      case 1:
        fprintf(fxyz, "\"13 Al\")\n");
        break;
      case 2:
        fprintf(fxyz, "\"14 Si\")\n");
        break;
      case 3:
      case 4:
        if (lattice->sites[i].state % 100 > 1)
          fprintf(fxyz, "\"1 H\")\n");
        else
          fprintf(fxyz, "\"8 O\")\n");
        break;
      case 5:
        if (lattice->sites[i].state % 100 > 1)
          fprintf(fxyz, "\"1 H\")\n");
        else
          fprintf(fxyz, "\"9 F\")\n");
        break;
      default:
        break;
      }
      x = (cell[lattice->sites[i].n].x / al + lattice->sites[i].a) * cd[0][0] +
          (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[0][1] +
          (cell[lattice->sites[i].n].z / cl) * cd[0][2];
      y = (cell[lattice->sites[i].n].y / bl + lattice->sites[i].b) * cd[1][1] +
          (cell[lattice->sites[i].n].z / cl) * cd[1][2];
      z = (cell[lattice->sites[i].n].z / cl) * cd[2][2];
      fprintf(fxyz, "  (A D XYZ (%f %f %f))\n", x, y, z);
      fprintf(fxyz, "  (A I Id %d)\n", nthing - 1);
      fprintf(fxyz, "  (A C Label \"");
      switch (lattice->sites[i].state / 100) {
      case 0:
        fprintf(fxyz, "Edge");
        break;
      case 1:
        fprintf(fxyz, "Al");
        break;
      case 2:
        fprintf(fxyz, "Si");
        break;
      case 3:
      case 4:
        fprintf(fxyz, "O");
        break;
      case 5:
        fprintf(fxyz, "OH");
        break;
      default:
        break;
      }
      fprintf(fxyz, "%d\")\n", i);
      fprintf(fxyz, "  (A I LabelType 0)\n");
      fprintf(fxyz, " )\n");
    }
  }

  for (i = 0; i < natom; i++) {
    if (bonds && (lattice->sites[i].state % 100 > 0 && lattice->sites[i].state != EDGE)) {
      for (j = 0; j < 6; j++) {
        if ((i2 = lattice->sites[i].nbr[j]) >= 0 && i2 < i && lattice->sites[i2].state % 100 != 0 &&
            lattice->sites[i2].state != EDGE &&
            !(lattice->sites[i].a == 0 && lattice->sites[i2].a == acells - 1 &&
              cell[lattice->sites[i].n].nbr[j].a == -1) &&
            !(lattice->sites[i].a == acells - 1 && lattice->sites[i2].a == 0 &&
              cell[lattice->sites[i].n].nbr[j].a == 1) &&
            !(lattice->sites[i].b == 0 && lattice->sites[i2].b == bcells - 1 &&
              cell[lattice->sites[i].n].nbr[j].b == -1) &&
            !(lattice->sites[i].b == bcells - 1 && lattice->sites[i2].b == 0 &&
              cell[lattice->sites[i].n].nbr[j].b == 1)) {
          fprintf(fxyz, " (%d Bond\n", ++nthing);
          fprintf(fxyz, "  (A O Atom1 %d)\n", id[i]);
          fprintf(fxyz, "  (A O Atom2 %d)\n", id[i2]);
          fprintf(fxyz, " )\n");
        }
      }
    }
  }

  fprintf(fxyz, ")\n");
  Futil::closeFile(fxyz);

  free(id);
}
