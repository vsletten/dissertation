#include <stdio.h>
#include "common.h"
#include "envrn.h"
#include "myerr.h"
#include "rxnlist.h"


static int check100(Lattice lattice, int site);
static int check200(Lattice lattice, int site);
static int check300(Lattice lattice, int site);
static int check400(Lattice lattice, int site);
static int check500(Lattice lattice, int site);


/* checkEnv: returns index which tells environment for determination
 *           of appropriate rate constant */
int checkEnv(Lattice lattice, int site)
{
  int status;
  char msg[100];

  switch (lattice[site].state / 100) {
    case 1: 
      if (lattice[site].state  == 100)
	status = 0;
      else
	status = check100(lattice, site);
      break;
    case 2:
      if (lattice[site].state == 200)
	status = 0;
      else
	status = check200(lattice, site);
      break;
    case 3: status = check300(lattice, site); break;
    case 4: 
      if (lattice[site].state == 404 || lattice[site].state == 405)
	status = check500(lattice, site);
      else
	status = check400(lattice, site);
      break;
    case 5: status = check500(lattice, site); break;
    default: 
      status = -1;
      sprintf(msg, "invalid environment: site %d, state %d", site,
	      lattice[site].state);
      die(msg);
      break;
  }
  return status;
}



/* check100: returns index for environment of an Al */
int check100(Lattice lattice, int site)
{
  int i, nbr, x, y;

  /* x = double 500 occupied ?
   * y = number of occupied 400's (0 - 2) */

  x = y = 0;
  for (i = 0; i < 6; i++) {
    nbr = lattice[site].nbr[i];
    if (lattice[nbr].state == EDGE)
      die("ran into lattice edge in check100");
    switch (lattice[nbr].state) {
      case 502: x++; break;
      case 403: case 405: case 407: case 409: case 410: y++; break;
      default: break;
    }
  }
  return ((x + y) / 2);
}



/* check200: returns index for environment of an Si */
int check200(Lattice lattice, int site)
{
  int i, nbr, x, y;

  /* x = 400 occupied ?
   * y = number of occupied 300's (0 - 3) */

  y = 0;
  x = 1;
  for (i = 0; i < 4; i++) {
    nbr = lattice[site].nbr[i];
    if (lattice[nbr].state == EDGE)
      die("ran into lattice edge in check200");
    switch (lattice[nbr].state) {
      case 408: x = 0; break;
      case 302: y++; break;
      default: break;
    }
  }
  return (x + y);
}



/* check300: returns index for environment of Si-O-Si type O */
int check300(Lattice lattice, int site)
{
  int si1, si2, sio1, sio2, x, y;

  /* x = number of broken 400s (0 - 2)
   * y = number of broken 300s (0 - 4) */

  x = y = 0;
  si1 = lattice[site].nbr[0];
  sio1 = lattice[si1].nbr[0];      /* Si-O-Al is Si nbr[0] */
  si2 = lattice[site].nbr[1];
  sio2 = lattice[si2].nbr[0];
  if (si1 < 0 || si2 < 0 || sio1 < 0 || sio2 < 0 ||
      lattice[si1].state == EDGE || lattice[si2].state ==EDGE ||
      lattice[sio1].state == EDGE || lattice[sio2].state == EDGE)
    die("ran into lattice edge in check300");

  if (lattice[sio1].state == 402 ||    /* 400s */
      lattice[sio1].state == 403 ||
      lattice[sio1].state == 407 ||
      lattice[sio1].state == 408)
    x++;
  if (lattice[sio2].state == 402 ||    /* 400s */
      lattice[sio2].state == 403 ||
      lattice[sio2].state == 407 ||
      lattice[sio2].state == 408)
    x++;

  y = (lattice[si1].state - 201) + (lattice[si2].state - 201) - x;  /* 300s */
  if (lattice[site].state > 301)
    y -= 2;

  return (5 * x + y);
}



/* check400: returns index for environment of Si-O-Al2 type O */
int check400(Lattice lattice, int site)
{
  int i, x, y, al1, al2, si, o1, p;

  /* x = number of broken 300s (0 - 3)
   * y = number of broken 400s and 500s (0 - 9) */

  x = y = 0;
  al1 = lattice[site].nbr[0];
  al2 = lattice[site].nbr[1];
  si = lattice[site].nbr[2];
  p = lattice[site].pair;
  if (al1 < 0 || al2 < 0 || si < 0 || (p = lattice[site].pair) < 0 ||
      lattice[al1].state == EDGE || lattice[al2].state == EDGE ||
      lattice[si].state == EDGE || lattice[p].state == EDGE)
    die("ran into lattice edge in check400");

  for (i = 1; i < 3; i++) {             /* 300s */
    if ((o1 = lattice[si].nbr[i]) < 0)  /* 400 site is Si nbr[0] */
      die("ran into lattice edge");
    if (lattice[o1].state > 301)
      x++;
  }

  y = (lattice[al1].state - 101) + (lattice[al2].state - 101);
  if (lattice[site].state == 403 ||     /* subtract out the site */
      lattice[site].state == 405)
    y -= 2;
  if (lattice[site].state == 407 ||
      lattice[site].state == 410 ||
      lattice[p].state == 502)
    y--;

  if (lattice[p].state == 502)
    y--;

  if (lattice[site].state == 406 || lattice[site].state == 407)
    return (x * 5 + y);
  else 
    return (x * 10 + y);
}



/* check500: returns index of environment for Al-OH-Al type O */
int check500(Lattice lattice, int site)
{
  int al1, al2, x, y, p;

  /* x = 0 if normal double bridge else 1 
   * y = number of broken 500s (0 - 9) */

  y = 0;
  al1 = lattice[site].nbr[0];
  al2 = lattice[site].nbr[1];
  if ((p = lattice[site].pair) < 0 || al1 < 0 || al2 < 0 ||
      lattice[al1].state == EDGE || lattice[al2].state == EDGE ||
      lattice[p].state == EDGE)
    die("ran into lattice edge in check500");

  if (lattice[p].state == 502 ||
      lattice[p].state == 403 ||
      lattice[p].state == 405 ||
      lattice[p].state == 410)
    x = 1;
  else 
    x = 0;
  
  y = (lattice[al1].state - 101) + (lattice[al2].state - 101);
  if (lattice[p].state == 502 ||     /* subtract out doubly counted bridges */
      lattice[p].state == 403 ||
      lattice[p].state == 405)
    y--;
  if (lattice[site].state == 502 ||  /* subtract out site */
      lattice[site].state == 405)
    y -= 2;

  return (x * 9 + y);
}      
  


/* isActive: determine whether site is active to do reaction rxn 
 *           return TRUE or FALSE */
int isActive(int site, Lattice lattice, int rxn)
{
  int result, nbr, nbr2, i, j;

  result = FALSE;                        
  if (rxn < NHYD && rxn % 2 == 0) {         /* hydrolysis
                                             * only allowed at "surface" */
    for (i = 0; (nbr = lattice[site].nbr[i]) >= 0 && !result; i++)
      for (j = 0; (nbr2 = lattice[nbr].nbr[j]) >= 0 && j < 6 && !result; j++)
	switch (lattice[nbr2].state) {
	  case 303: case 404: case 405: case 406: case 408: case 409:
	            result = TRUE; break;
	  default: break;
	}
  } else if (rxn < NHYD) {             /* reverse hydrolysis
					* always allowed*/
    result = TRUE;
    
  } else if (rxn == 16 || rxn == 19) {  /* adsorption to correct site */
    result = FALSE;                    /* must be at least one occupied nbr */
    for (i = 0; (nbr = lattice[site].nbr[i]) >= 0 && i < 6; i++)
      if (ISOCC(lattice[nbr]) && !ISEDGE(lattice[nbr]))
	result = TRUE;   
  } else if (rxn == 20 || rxn == 22) {  /* desorption */
    result = TRUE;
  } else {                              /* diffusion */
    result = FALSE;
  }
  return result;
}




