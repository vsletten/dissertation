#include <stdio.h>
//#include "common.hpp"
#include "reactions.hpp"
#include "myerr.hpp"
#include "ran2.hpp"

/* adsorbAl: update state of O site upon adsorbing Al(OH,H2O)6 */
void Reactions::adsorbAl(Lattice *lattice, int site)
{
  int i, nbr;

  if (ISAL(lattice->sites[site])) {
    lattice->sites[site].state = 107;
    for (i = 0; i < 6; i++) {
      nbr = lattice->sites[site].nbr[i];
      switch (lattice->sites[nbr].state) {
        case 503: lattice->sites[nbr].state = 502; break;
	case 500: lattice->sites[nbr].state = 503; break;
	case 407: lattice->sites[nbr].state = 403; break;
	case 409: lattice->sites[nbr].state = 405; break;
	case 408: lattice->sites[nbr].state = 407; break;
	case 400: lattice->sites[nbr].state = 409; break;
	case 406: 
	  lattice->sites[nbr].state = 410; 
	  lattice->sites[nbr].lostal = site;
	  break;
	default: 
	  printf("state %d\n", lattice->sites[nbr].state);
	  Myerr::die("invalid state in adsorbAl");
	  break;
      }
    }
  } else {
    lattice->sites[site].state = 200 + WRONG;
  }
}



/* adsorbSi: update state of O site upon adsorbing Si(OH)4 */
void Reactions::adsorbSi(Lattice *lattice, int site)
{
  int i, nbr;

  if (ISSI(lattice->sites[site])) {
    lattice->sites[site].state = 205;
    for (i = 0; i < 4; i++) {
      nbr = lattice->sites[site].nbr[i];
      switch (lattice->sites[nbr].state) {
        case 300: lattice->sites[nbr].state = 303; break;
	case 303: lattice->sites[nbr].state = 302; break;
	case 404: lattice->sites[nbr].state = 402; break;
	case 405: lattice->sites[nbr].state = 403; break;
	case 409: lattice->sites[nbr].state = 407; break;
	case 400: lattice->sites[nbr].state = 408; break;
	default:
	  printf("state %d\n", lattice->sites[nbr].state);
	  Myerr::die("invalid state in adsorbSi");
	  break;
      }
    }
  } else {
    lattice->sites[site].state = 100 + WRONG;
  }
}


/* desorbAl: update state of O site upon desorbing Al(OH,H2O)6 */
void Reactions::desorbAl(Lattice *lattice, int site)
{
  int i, nbr;

  if (ISAL(lattice->sites[site])) {
    lattice->sites[site].state = 100;
    for (i = 0; i < 6; i++) {
      nbr = lattice->sites[site].nbr[i];             /* nbr is an o site */
      switch (lattice->sites[nbr].state) {
        case 502: lattice->sites[nbr].state = 503; break;
	      case 503: lattice->sites[nbr].state = 500; break;
	      case 403: lattice->sites[nbr].state = 407; break;
	      case 405: lattice->sites[nbr].state = 409; break;
	      case 407: lattice->sites[nbr].state = 408; break;
	      case 409: lattice->sites[nbr].state = 400; break;
	      case 410: 
	        lattice->sites[nbr].state = 406; 
	        lattice->sites[nbr].lostal = -1;
	        break;
	      default: break;
      }
    }
  } else {
    lattice->sites[site].state = 200;
  }
}



/* desorbSi: update state of O site upon desorbing Si(OH)4 */
void Reactions::desorbSi(Lattice *lattice, int site)
{
  int i, nbr;

  if (ISSI(lattice->sites[site])) {
    lattice->sites[site].state = 200;
    for (i = 0; i < 4; i++) {
      nbr = lattice->sites[site].nbr[i];
      switch (lattice->sites[nbr].state) {
        case 303: lattice->sites[nbr].state = 300; break;
	case 302: lattice->sites[nbr].state = 303; break;
	case 402: lattice->sites[nbr].state = 404; break;
	case 403: lattice->sites[nbr].state = 405; break;
	case 407: lattice->sites[nbr].state = 409; break;
	case 408: lattice->sites[nbr].state = 400; break;
	default: break;
      }
    }
  } else {
    lattice->sites[site].state = 100;
  }
}




void Reactions::diffuse(Lattice *lattice, int source, int target)
{
  if (lattice->sites[source].state == 107 || lattice->sites[source].state == 299) {
    desorbAl(lattice, source);
    adsorbAl(lattice, target);
  } else {
    desorbSi(lattice, source);
    adsorbSi(lattice, target);
  }
}



/* r0: =Si-O-Si= + H2O => 2(=Si-OH) */
void Reactions::r0(Lattice *lattice, int site)
{
  int si, i;

  lattice->sites[site].state = 302;
  for (i = 0;  i < 2; i++) {
    si = lattice->sites[site].nbr[i];
    lattice->sites[si].state++;
  }
}



/* r1: 2(=Si-OH) => =Si-O-Si= + H2O */
void Reactions::r1(Lattice *lattice, int site)
{
  int si, i;

  lattice->sites[site].state = 301;
  for (i = 0;  i < 2; i++) {
    si = lattice->sites[site].nbr[i];
    lattice->sites[si].state--;
  }
}



/* r2: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
void Reactions::r2(Lattice *lattice, int site)
{
  int nbr;

  lattice->sites[site].state = 402;
  nbr = lattice->sites[site].nbr[2];
  lattice->sites[nbr].state++;  /* increment for si only */
}



/* r3: =Si-OH + =Al-OH-Al= => =Si-O-2(Al=) + H2O */
void Reactions::r3(Lattice *lattice, int site)
{
  int nbr;

  lattice->sites[site].state = 401;
  nbr = lattice->sites[site].nbr[2];
  lattice->sites[nbr].state--;  /* decrement for si only */
}



/* r4: =Si-O-2(Al=) + H2O => =Si-OH-Al= + =Al-OH */
void Reactions::r4(Lattice *lattice, int site)
{
  int nbr;
  float r;

  lattice->sites[site].state = 410;
  r = ran2();
  if (r < 0.5)
    nbr = lattice->sites[site].nbr[0];
  else 
    nbr = lattice->sites[site].nbr[1];
  lattice->sites[nbr].state++;
  lattice->sites[site].lostal = nbr;
}



/* r5: =Si-OH-Al= + =Al-OH => =Si-O-2(Al=) + H2O */
void Reactions::r5(Lattice *lattice, int site)
{
  int nbr;

  lattice->sites[site].state = 401;
  nbr = lattice->sites[site].lostal;
  lattice->sites[nbr].state--;
  lattice->sites[site].lostal = -1;
}



/* r6: =Si-OH + =Al-OH-Al= + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Reactions::r6(Lattice *lattice, int site)
{
  int i, nbr;

  lattice->sites[site].state = 403;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state++;
  }
}



/* r7: =Si-OH + =Al-OH + =Al-H2O => =Si-OH + =Al-OH-Al= + H2O */
void Reactions::r7(Lattice *lattice, int site)
{
  int i, nbr;

  lattice->sites[site].state = 402;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state--;
  }
}



/* r8: =Si-OH-Al= + =Al-OH + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Reactions::r8(Lattice *lattice, int site)
{
  int nbr;

  lattice->sites[site].state = 403;
  nbr = lattice->sites[site].nbr[2];   /* si */
  lattice->sites[nbr].state++;
  nbr = (lattice->sites[site].lostal == lattice->sites[site].nbr[0]) ? 
         lattice->sites[site].nbr[1] : lattice->sites[site].nbr[0];
  lattice->sites[nbr].state++;
  lattice->sites[site].lostal = -1;
}




/* r9: =Si-OH + =Al-OH + =Al-H2O => =Si-OH-Al= + =Al-OH + H2O */
void Reactions::r9(Lattice *lattice, int site)
{
  int nbr, lost;
  float r;

  lattice->sites[site].state = 410;
  r = ran2();
  nbr = lattice->sites[site].nbr[2];   /* si */
  lattice->sites[nbr].state--;
  if (r < 0.5) {
    nbr = lattice->sites[site].nbr[0];   /* al */
    lost = lattice->sites[site].nbr[1];
  } else {
    nbr = lattice->sites[site].nbr[1];
    lost = lattice->sites[site].nbr[0];
  }
  lattice->sites[nbr].state--;
  lattice->sites[site].lostal = lost;
}



/* r10: =Si-OH-Al= + H2O => =Si-OH + =Al=H2O */
void Reactions::r10(Lattice *lattice, int site)
{
  int nbr, al1, al2;

  lattice->sites[site].state = 407;
  nbr = lattice->sites[site].nbr[2];   /* si */
  lattice->sites[nbr].state++;
  al1 = lattice->sites[site].nbr[0];
  al2 = lattice->sites[site].nbr[1];
  nbr = (lattice->sites[al1].state == 100) ? al2 : al1;
  lattice->sites[nbr].state++;         /* al */
}



/* r11: =Si-OH + =Al=H2O => =Si-OH-Al= + H2O */
void Reactions::r11(Lattice *lattice, int site)
{
  int nbr, al1, al2;

  lattice->sites[site].state = 406;
  nbr = lattice->sites[site].nbr[2];
  lattice->sites[nbr].state--;
  al1 = lattice->sites[site].nbr[0];
  al2 = lattice->sites[site].nbr[1];
  nbr = (lattice->sites[al1].state == 100) ? al2 : al1;
  lattice->sites[nbr].state--;         /* al */
}



/* r12: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (400s) */
void Reactions::r12(Lattice *lattice, int site)
{
  int nbr, i;

  lattice->sites[site].state = 405;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state++;
  }
}



/* r13:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (400s) */
void Reactions::r13(Lattice *lattice, int site)
{
  int nbr, i;

  lattice->sites[site].state = 404;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state--;
  }
}



/* r14: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (500s) */
void Reactions::r14(Lattice *lattice, int site)
{
  int nbr, i;

  lattice->sites[site].state = 502;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state++;
  }
}



/* r15:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (500s) */
void Reactions::r15(Lattice *lattice, int site)
{
  int nbr, i;

  lattice->sites[site].state = 501;
  for (i = 0; i < 2; i++) {
    nbr = lattice->sites[site].nbr[i];
    lattice->sites[nbr].state--;
  }
}
