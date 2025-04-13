#include <stdio.h>
#include "common.hpp"
#include "envrn.hpp"
#include "myerr.hpp"
#include "rxnlist.hpp"


/* checkEnv: returns index which tells environment for determination
 *           of appropriate rate constant */
int Environment::checkEnv(Lattice *lattice, int site)
{
  int status;
  char msg[100];

  switch (lattice->sites[site].state / 100) {
    case 1: 
      if (lattice->sites[site].state  == 100)
	status = 0;
      else
	status = check100(lattice, site);
      break;
    case 2:
      if (lattice->sites[site].state == 200)
	status = 0;
      else
	status = check200(lattice, site);
      break;
    case 3: status = check300(lattice, site); break;
    case 4: 
      if (lattice->sites[site].state == 404 || lattice->sites[site].state == 405)
	status = check500(lattice, site);
      else
	status = check400(lattice, site);
      break;
    case 5: status = check500(lattice, site); break;
    default: 
      status = -1;
      sprintf(msg, "invalid environment: site %d, state %d", site,
	      lattice->sites[site].state);
      Myerr::die(msg);
      break;
  }
  return status;
}



/* check100: returns index for environment of an Al */
int Environment::check100(Lattice *lattice, int site)
{
  int i, nbr, x, y;

  /* x = double 500 occupied ?
   * y = number of occupied 400's (0 - 2) */

  x = y = 0;
  for (i = 0; i < 6; i++) {
    nbr = lattice->sites[site].nbr[i];
    if (lattice->sites[nbr].state == EDGE)
      Myerr::die("ran into lattice edge in check100");
    switch (lattice->sites[nbr].state) {
      case 502: x++; break;
      case 403: case 405: case 407: case 409: case 410: y++; break;
      default: break;
    }
  }
  return ((x + y) / 2);
}



/* check200: returns index for environment of an Si */
int Environment::check200(Lattice *lattice, int site)
{
  int i, nbr, x, y;

  /* x = 400 occupied ?
   * y = number of occupied 300's (0 - 3) */

  y = 0;
  x = 1;
  for (i = 0; i < 4; i++) {
    nbr = lattice->sites[site].nbr[i];
    if (lattice->sites[nbr].state == EDGE)
      Myerr::die("ran into lattice edge in check200");
    switch (lattice->sites[nbr].state) {
      case 408: x = 0; break;
      case 302: y++; break;
      default: break;
    }
  }
  return (x + y);
}



/* check300: returns index for environment of Si-O-Si type O */
int Environment::check300(Lattice *lattice, int site)
{
  int si1, si2, sio1, sio2, x, y;

  /* x = number of broken 400s (0 - 2)
   * y = number of broken 300s (0 - 4) */

  x = y = 0;
  si1 = lattice->sites[site].nbr[0];
  sio1 = lattice->sites[si1].nbr[0];      /* Si-O-Al is Si nbr[0] */
  si2 = lattice->sites[site].nbr[1];
  sio2 = lattice->sites[si2].nbr[0];
  if (si1 < 0 || si2 < 0 || sio1 < 0 || sio2 < 0 ||
      lattice->sites[si1].state == EDGE || lattice->sites[si2].state ==EDGE ||
      lattice->sites[sio1].state == EDGE || lattice->sites[sio2].state == EDGE)
    Myerr::die("ran into lattice edge in check300");

  if (lattice->sites[sio1].state == 402 ||    /* 400s */
      lattice->sites[sio1].state == 403 ||
      lattice->sites[sio1].state == 407 ||
      lattice->sites[sio1].state == 408)
    x++;
  if (lattice->sites[sio2].state == 402 ||    /* 400s */
      lattice->sites[sio2].state == 403 ||
      lattice->sites[sio2].state == 407 ||
      lattice->sites[sio2].state == 408)
    x++;

  y = (lattice->sites[si1].state - 201) + (lattice->sites[si2].state - 201) - x;  /* 300s */
  if (lattice->sites[site].state > 301)
    y -= 2;

  return (5 * x + y);
}



/* check400: returns index for environment of Si-O-Al2 type O */
int Environment::check400(Lattice *lattice, int site)
{
  int i, x, y, al1, al2, si, o1, p;

  /* x = number of broken 300s (0 - 3)
   * y = number of broken 400s and 500s (0 - 9) */

  x = y = 0;
  al1 = lattice->sites[site].nbr[0];
  al2 = lattice->sites[site].nbr[1];
  si = lattice->sites[site].nbr[2];
  p = lattice->sites[site].pair;
  if (al1 < 0 || al2 < 0 || si < 0 || (p = lattice->sites[site].pair) < 0 ||
      lattice->sites[al1].state == EDGE || lattice->sites[al2].state == EDGE ||
      lattice->sites[si].state == EDGE || lattice->sites[p].state == EDGE)
    Myerr::die("ran into lattice edge in check400");

  for (i = 1; i < 3; i++) {             /* 300s */
    if ((o1 = lattice->sites[si].nbr[i]) < 0)  /* 400 site is Si nbr[0] */
      Myerr::die("ran into lattice edge");
    if (lattice->sites[o1].state > 301)
      x++;
  }

  y = (lattice->sites[al1].state - 101) + (lattice->sites[al2].state - 101);
  if (lattice->sites[site].state == 403 ||     /* subtract out the site */
      lattice->sites[site].state == 405)
    y -= 2;
  if (lattice->sites[site].state == 407 ||
      lattice->sites[site].state == 410 ||
      lattice->sites[p].state == 502)
    y--;

  if (lattice->sites[p].state == 502)
    y--;

  if (lattice->sites[site].state == 406 || lattice->sites[site].state == 407)
    return (x * 5 + y);
  else 
    return (x * 10 + y);
}



/* check500: returns index of environment for Al-OH-Al type O */
int Environment::check500(Lattice *lattice, int site)
{
  int al1, al2, x, y, p;

  /* x = 0 if normal double bridge else 1 
   * y = number of broken 500s (0 - 9) */

  y = 0;
  al1 = lattice->sites[site].nbr[0];
  al2 = lattice->sites[site].nbr[1];
  if ((p = lattice->sites[site].pair) < 0 || al1 < 0 || al2 < 0 ||
      lattice->sites[al1].state == EDGE || lattice->sites[al2].state == EDGE ||
      lattice->sites[p].state == EDGE)
    Myerr::die("ran into lattice edge in check500");

  if (lattice->sites[p].state == 502 ||
      lattice->sites[p].state == 403 ||
      lattice->sites[p].state == 405 ||
      lattice->sites[p].state == 410)
    x = 1;
  else 
    x = 0;
  
  y = (lattice->sites[al1].state - 101) + (lattice->sites[al2].state - 101);
  if (lattice->sites[p].state == 502 ||     /* subtract out doubly counted bridges */
      lattice->sites[p].state == 403 ||
      lattice->sites[p].state == 405)
    y--;
  if (lattice->sites[site].state == 502 ||  /* subtract out site */
      lattice->sites[site].state == 405)
    y -= 2;

  return (x * 9 + y);
}      
  


/* isActive: determine whether site is active to do reaction rxn 
 *           return TRUE or FALSE */
int Environment::isActive(int site, Lattice *lattice, int rxn)
{
  int result, nbr, nbr2, i, j;

  result = FALSE;                        
  if (rxn < NHYD && rxn % 2 == 0) {         /* hydrolysis
                                             * only allowed at "surface" */
    for (i = 0; (nbr = lattice->sites[site].nbr[i]) >= 0 && !result; i++)
      for (j = 0; (nbr2 = lattice->sites[nbr].nbr[j]) >= 0 && j < 6 && !result; j++)
	switch (lattice->sites[nbr2].state) {
	  case 303: case 404: case 405: case 406: case 408: case 409:
	            result = TRUE; break;
	  default: break;
	}
  } else if (rxn < NHYD) {             /* reverse hydrolysis
					* always allowed*/
    result = TRUE;
    
  } else if (rxn == 16 || rxn == 19) {  /* adsorption to correct site */
    result = FALSE;                    /* must be at least one occupied nbr */
    for (i = 0; (nbr = lattice->sites[site].nbr[i]) >= 0 && i < 6; i++)
      if (ISOCC(lattice->sites[nbr]) && !ISEDGE(lattice->sites[nbr]))
	result = TRUE;   
  } else if (rxn == 20 || rxn == 22) {  /* desorption */
    result = TRUE;
  } else {                              /* diffusion */
    result = FALSE;
  }
  return result;
}




