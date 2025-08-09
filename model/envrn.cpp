#include "envrn.hpp"
#include "common.hpp"
#include "rxnlist.hpp"
#include <iostream>

/* checkEnv: returns index which tells environment for determination
 *           of appropriate rate constant */
int Environment::CheckEnv(int site) {
  int status;

  switch (this->lattice->sites[site].state / 100) {
  case 1:
    if (this->lattice->sites[site].state == 100)
      status = 0;
    else
      status = Check100(site);
    break;
  case 2:
    if (this->lattice->sites[site].state == 200)
      status = 0;
    else
      status = Check200(site);
    break;
  case 3:
    status = Check300(site);
    break;
  case 4:
    if (this->lattice->sites[site].state == 404 ||
        this->lattice->sites[site].state == 405)
      status = Check500(site);
    else
      status = Check400(site);
    break;
  case 5:
    status = Check500(site);
    break;
  default:
    status = -1;
    std::cerr << "invalid environment: site " << site << ", state "
              << this->lattice->sites[site].state << std::endl;
    break;
  }
  return status;
}

/* check100: returns index for environment of an Al */
int Environment::Check100(int site) {
  int i, nbr, x, y;

  /* x = double 500 occupied ?
   * y = number of occupied 400's (0 - 2) */

  x = y = 0;
  for (i = 0; i < 6; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    if (this->lattice->sites[nbr].state == EDGE) {
      std::cerr << "ran into lattice edge in check100" << std::endl;
      return -1;
    }
    switch (this->lattice->sites[nbr].state) {
    case 502:
      x++;
      break;
    case 403:
    case 405:
    case 407:
    case 409:
    case 410:
      y++;
      break;
    default:
      break;
    }
  }
  return ((x + y) / 2);
}

/* check200: returns index for environment of an Si */
int Environment::Check200(int site) {
  int i, nbr, x, y;

  /* x = 400 occupied ?
   * y = number of occupied 300's (0 - 3) */

  y = 0;
  x = 1;
  for (i = 0; i < 4; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    if (this->lattice->sites[nbr].state == EDGE) {
      std::cerr << "ran into lattice edge in check200" << std::endl;
      return -1;
    }
    switch (this->lattice->sites[nbr].state) {
    case 408:
      x = 0;
      break;
    case 302:
      y++;
      break;
    default:
      break;
    }
  }
  return (x + y);
}

/* check300: returns index for environment of Si-O-Si type O */
int Environment::Check300(int site) {
  int si1, si2, sio1, sio2, x, y;

  /* x = number of broken 400s (0 - 2)
   * y = number of broken 300s (0 - 4) */

  x = y = 0;
  si1 = this->lattice->sites[site].nbr[0];
  si2 = this->lattice->sites[site].nbr[1];
  
  if (si1 < 0 || si2 < 0) {
    std::cerr << "ran into lattice edge in check300" << std::endl;
    std::cerr << "DEBUG: site " << site << " has missing Si neighbors: si1=" << si1 << " si2=" << si2 << std::endl;
    return -1;
  }
  
  sio1 = this->lattice->sites[si1].nbr[0]; /* Si-O-Al is Si nbr[0] */
  sio2 = this->lattice->sites[si2].nbr[0];
  
  if (sio1 < 0 || sio2 < 0 ||
      this->lattice->sites[si1].state == EDGE ||
      this->lattice->sites[si2].state == EDGE ||
      this->lattice->sites[sio1].state == EDGE ||
      this->lattice->sites[sio2].state == EDGE) {
    std::cerr << "ran into lattice edge in check300" << std::endl;
    return -1;
  }

  if (this->lattice->sites[sio1].state == 402 || /* 400s */
      this->lattice->sites[sio1].state == 403 ||
      this->lattice->sites[sio1].state == 407 ||
      this->lattice->sites[sio1].state == 408) {
    x++;
  }
  if (this->lattice->sites[sio2].state == 402 || /* 400s */
      this->lattice->sites[sio2].state == 403 ||
      this->lattice->sites[sio2].state == 407 ||
      this->lattice->sites[sio2].state == 408) {
    x++;
  }

  y = (this->lattice->sites[si1].state - 201) +
      (this->lattice->sites[si2].state - 201) - x; /* 300s */
  if (this->lattice->sites[site].state > 301) {
    y -= 2;
  }

  return (5 * x + y);
}

/* check400: returns index for environment of Si-O-Al2 type O */
int Environment::Check400(int site) {
  int i, x, y, al1, al2, si, o1, p;

  /* x = number of broken 300s (0 - 3)
   * y = number of broken 400s and 500s (0 - 9) */

  x = y = 0;
  al1 = this->lattice->sites[site].nbr[0];
  al2 = this->lattice->sites[site].nbr[1];
  si = this->lattice->sites[site].nbr[2];
  p = this->lattice->sites[site].pair;
  if (al1 < 0 || al2 < 0 || si < 0 ||
      (p = this->lattice->sites[site].pair) < 0 ||
      this->lattice->sites[al1].state == EDGE ||
      this->lattice->sites[al2].state == EDGE ||
      this->lattice->sites[si].state == EDGE ||
      this->lattice->sites[p].state == EDGE) {
    std::cerr << "ran into lattice edge in check400" << std::endl;
    return -1;
  }

  for (i = 1; i < 3; i++) {                         /* 300s */
    if ((o1 = this->lattice->sites[si].nbr[i]) < 0) /* 400 site is Si nbr[0] */
    {
      std::cerr << "ran into lattice edge in check400" << std::endl;
      return -1;
    }
    if (this->lattice->sites[o1].state > 301)
      x++;
  }

  y = (this->lattice->sites[al1].state - 101) +
      (this->lattice->sites[al2].state - 101);
  if (this->lattice->sites[site].state == 403 || /* subtract out the site */
      this->lattice->sites[site].state == 405) {
    y -= 2;
  }
  if (this->lattice->sites[site].state == 407 ||
      this->lattice->sites[site].state == 410 ||
      this->lattice->sites[p].state == 502) {
    y--;
  }

  if (this->lattice->sites[p].state == 502) {
    y--;
  }

  if (this->lattice->sites[site].state == 406 ||
      this->lattice->sites[site].state == 407) {
    return (x * 5 + y);
  } else {
    return (x * 10 + y);
  }
}

/* check500: returns index of environment for Al-OH-Al type O */
int Environment::Check500(int site) {
  int al1, al2, x, y, p;

  /* x = 0 if normal double bridge else 1
   * y = number of broken 500s (0 - 9) */

  y = 0;
  al1 = this->lattice->sites[site].nbr[0];
  al2 = this->lattice->sites[site].nbr[1];
  if ((p = this->lattice->sites[site].pair) < 0 || al1 < 0 || al2 < 0 ||
      this->lattice->sites[al1].state == EDGE ||
      this->lattice->sites[al2].state == EDGE ||
      this->lattice->sites[p].state == EDGE) {
    std::cerr << "ran into lattice edge in check500" << std::endl;
    return -1;
  }

  if (this->lattice->sites[p].state == 502 ||
      this->lattice->sites[p].state == 403 ||
      this->lattice->sites[p].state == 405 ||
      this->lattice->sites[p].state == 410) {
    x = 1;
  } else {
    x = 0;
  }

  y = (this->lattice->sites[al1].state - 101) +
      (this->lattice->sites[al2].state - 101);
  if (this->lattice->sites[p].state ==
          502 || /* subtract out doubly counted bridges */
      this->lattice->sites[p].state == 403 ||
      this->lattice->sites[p].state == 405) {
    y--;
  }
  if (this->lattice->sites[site].state == 502 || /* subtract out site */
      this->lattice->sites[site].state == 405) {
    y -= 2;
  }

  return (x * 9 + y);
}

/* isActive: determine whether site is active to do reaction rxn
 *           return TRUE or FALSE */
int Environment::IsActive(int site, int rxn) {
  int result, nbr, nbr2, i, j;

  result = FALSE;
  if (rxn < NHYD && rxn % 2 == 0) { /* hydrolysis
                                     * only allowed at "surface" */
    for (i = 0; (nbr = this->lattice->sites[site].nbr[i]) >= 0 && !result;
         i++) 
    {
      for (j = 0;
           (nbr2 = this->lattice->sites[nbr].nbr[j]) >= 0 && j < 6 && !result;
           j++) 
      {
        switch (this->lattice->sites[nbr2].state) 
        {
          case 303:
          case 404:
          case 405:
          case 406:
          case 408:
          case 409:
          {
            result = TRUE;
            break;
          }
          default:
            break;
        }
      }
    }
  }
  else if (rxn < NHYD) 
  { /* reverse hydrolysis
     * always allowed*/
    result = TRUE;
  }
  else if (rxn == 16 || rxn == 19) 
  { /* adsorption to correct site */
        result = FALSE;                  /* must be at least one occupied nbr */
        for (i = 0; (nbr = this->lattice->sites[site].nbr[i]) >= 0 && i < 6;
             i++)
        {
          if (ISOCC(this->lattice->sites[nbr]) &&
              !ISEDGE(this->lattice->sites[nbr]))
          {
            result = TRUE;
          }
        }
  }
  else if (rxn == 20 || rxn == 22) 
  { /* desorption */
    result = TRUE;
  }
  else 
  { /* diffusion */
    result = FALSE;
  }
  return result;
  return result;
}
