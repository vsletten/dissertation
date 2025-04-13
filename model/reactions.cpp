#include <stdio.h>
//#include "common.hpp"
#include "reactions.hpp"
#include "myerr.hpp"
#include "ran2.hpp"

/* adsorbAl: update state of O site upon adsorbing Al(OH,H2O)6 */
void Reactions::adsorbAl(lattice::Lattice lattice, int site)
{
  int i, nbr;

  if (ISAL(lattice[site])) {
    lattice[site].state = 107;
    for (i = 0; i < 6; i++) {
      nbr = lattice[site].nbr[i];
      switch (lattice[nbr].state) {
        case 503: lattice[nbr].state = 502; break;
	case 500: lattice[nbr].state = 503; break;
	case 407: lattice[nbr].state = 403; break;
	case 409: lattice[nbr].state = 405; break;
	case 408: lattice[nbr].state = 407; break;
	case 400: lattice[nbr].state = 409; break;
	case 406: 
	  lattice[nbr].state = 410; 
	  lattice[nbr].lostal = site;
	  break;
	default: 
	  printf("state %d\n", lattice[nbr].state);
	  Myerr::die("invalid state in adsorbAl");
	  break;
      }
    }
  } else {
    lattice[site].state = 200 + WRONG;
  }
}



/* adsorbSi: update state of O site upon adsorbing Si(OH)4 */
void Reactions::adsorbSi(lattice::Lattice lattice, int site)
{
  int i, nbr;

  if (ISSI(lattice[site])) {
    lattice[site].state = 205;
    for (i = 0; i < 4; i++) {
      nbr = lattice[site].nbr[i];
      switch (lattice[nbr].state) {
        case 300: lattice[nbr].state = 303; break;
	case 303: lattice[nbr].state = 302; break;
	case 404: lattice[nbr].state = 402; break;
	case 405: lattice[nbr].state = 403; break;
	case 409: lattice[nbr].state = 407; break;
	case 400: lattice[nbr].state = 408; break;
	default:
	  printf("state %d\n", lattice[nbr].state);
	  Myerr::die("invalid state in adsorbSi");
	  break;
      }
    }
  } else {
    lattice[site].state = 100 + WRONG;
  }
}


/* desorbAl: update state of O site upon desorbing Al(OH,H2O)6 */
void Reactions::desorbAl(lattice::Lattice lattice, int site)
{
  int i, nbr;

  if (ISAL(lattice[site])) {
    lattice[site].state = 100;
    for (i = 0; i < 6; i++) {
      nbr = lattice[site].nbr[i];             /* nbr is an o site */
      switch (lattice[nbr].state) {
        case 502: lattice[nbr].state = 503; break;
	      case 503: lattice[nbr].state = 500; break;
	      case 403: lattice[nbr].state = 407; break;
	      case 405: lattice[nbr].state = 409; break;
	      case 407: lattice[nbr].state = 408; break;
	      case 409: lattice[nbr].state = 400; break;
	      case 410: 
	        lattice[nbr].state = 406; 
	        lattice[nbr].lostal = -1;
	        break;
	      default: break;
      }
    }
  } else {
    lattice[site].state = 200;
  }
}



/* desorbSi: update state of O site upon desorbing Si(OH)4 */
void Reactions::desorbSi(lattice::Lattice lattice, int site)
{
  int i, nbr;

  if (ISSI(lattice[site])) {
    lattice[site].state = 200;
    for (i = 0; i < 4; i++) {
      nbr = lattice[site].nbr[i];
      switch (lattice[nbr].state) {
        case 303: lattice[nbr].state = 300; break;
	case 302: lattice[nbr].state = 303; break;
	case 402: lattice[nbr].state = 404; break;
	case 403: lattice[nbr].state = 405; break;
	case 407: lattice[nbr].state = 409; break;
	case 408: lattice[nbr].state = 400; break;
	default: break;
      }
    }
  } else {
    lattice[site].state = 100;
  }
}




void Reactions::diffuse(lattice::Lattice lattice, int source, int target)
{
  if (lattice[source].state == 107 || lattice[source].state == 299) {
    desorbAl(lattice, source);
    adsorbAl(lattice, target);
  } else {
    desorbSi(lattice, source);
    adsorbSi(lattice, target);
  }
}



/* r0: =Si-O-Si= + H2O => 2(=Si-OH) */
void Reactions::r0(lattice::Lattice lattice, int site)
{
  int si, i;

  lattice[site].state = 302;
  for (i = 0;  i < 2; i++) {
    si = lattice[site].nbr[i];
    lattice[si].state++;
  }
}



/* r1: 2(=Si-OH) => =Si-O-Si= + H2O */
void Reactions::r1(lattice::Lattice lattice, int site)
{
  int si, i;

  lattice[site].state = 301;
  for (i = 0;  i < 2; i++) {
    si = lattice[site].nbr[i];
    lattice[si].state--;
  }
}



/* r2: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
void Reactions::r2(lattice::Lattice lattice, int site)
{
  int nbr;

  lattice[site].state = 402;
  nbr = lattice[site].nbr[2];
  lattice[nbr].state++;  /* increment for si only */
}



/* r3: =Si-OH + =Al-OH-Al= => =Si-O-2(Al=) + H2O */
void Reactions::r3(lattice::Lattice lattice, int site)
{
  int nbr;

  lattice[site].state = 401;
  nbr = lattice[site].nbr[2];
  lattice[nbr].state--;  /* decrement for si only */
}



/* r4: =Si-O-2(Al=) + H2O => =Si-OH-Al= + =Al-OH */
void Reactions::r4(lattice::Lattice lattice, int site)
{
  int nbr;
  float r;

  lattice[site].state = 410;
  r = ran2();
  if (r < 0.5)
    nbr = lattice[site].nbr[0];
  else 
    nbr = lattice[site].nbr[1];
  lattice[nbr].state++;
  lattice[site].lostal = nbr;
}



/* r5: =Si-OH-Al= + =Al-OH => =Si-O-2(Al=) + H2O */
void Reactions::r5(lattice::Lattice lattice, int site)
{
  int nbr;

  lattice[site].state = 401;
  nbr = lattice[site].lostal;
  lattice[nbr].state--;
  lattice[site].lostal = -1;
}



/* r6: =Si-OH + =Al-OH-Al= + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Reactions::r6(lattice::Lattice lattice, int site)
{
  int i, nbr;

  lattice[site].state = 403;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state++;
  }
}



/* r7: =Si-OH + =Al-OH + =Al-H2O => =Si-OH + =Al-OH-Al= + H2O */
void Reactions::r7(lattice::Lattice lattice, int site)
{
  int i, nbr;

  lattice[site].state = 402;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state--;
  }
}



/* r8: =Si-OH-Al= + =Al-OH + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Reactions::r8(lattice::Lattice lattice, int site)
{
  int nbr;

  lattice[site].state = 403;
  nbr = lattice[site].nbr[2];   /* si */
  lattice[nbr].state++;
  nbr = (lattice[site].lostal == lattice[site].nbr[0]) ? 
         lattice[site].nbr[1] : lattice[site].nbr[0];
  lattice[nbr].state++;
  lattice[site].lostal = -1;
}




/* r9: =Si-OH + =Al-OH + =Al-H2O => =Si-OH-Al= + =Al-OH + H2O */
void Reactions::r9(lattice::Lattice lattice, int site)
{
  int nbr, lost;
  float r;

  lattice[site].state = 410;
  r = ran2();
  nbr = lattice[site].nbr[2];   /* si */
  lattice[nbr].state--;
  if (r < 0.5) {
    nbr = lattice[site].nbr[0];   /* al */
    lost = lattice[site].nbr[1];
  } else {
    nbr = lattice[site].nbr[1];
    lost = lattice[site].nbr[0];
  }
  lattice[nbr].state--;
  lattice[site].lostal = lost;
}



/* r10: =Si-OH-Al= + H2O => =Si-OH + =Al=H2O */
void Reactions::r10(lattice::Lattice lattice, int site)
{
  int nbr, al1, al2;

  lattice[site].state = 407;
  nbr = lattice[site].nbr[2];   /* si */
  lattice[nbr].state++;
  al1 = lattice[site].nbr[0];
  al2 = lattice[site].nbr[1];
  nbr = (lattice[al1].state == 100) ? al2 : al1;
  lattice[nbr].state++;         /* al */
}



/* r11: =Si-OH + =Al=H2O => =Si-OH-Al= + H2O */
void Reactions::r11(lattice::Lattice lattice, int site)
{
  int nbr, al1, al2;

  lattice[site].state = 406;
  nbr = lattice[site].nbr[2];
  lattice[nbr].state--;
  al1 = lattice[site].nbr[0];
  al2 = lattice[site].nbr[1];
  nbr = (lattice[al1].state == 100) ? al2 : al1;
  lattice[nbr].state--;         /* al */
}



/* r12: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (400s) */
void Reactions::r12(lattice::Lattice lattice, int site)
{
  int nbr, i;

  lattice[site].state = 405;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state++;
  }
}



/* r13:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (400s) */
void Reactions::r13(lattice::Lattice lattice, int site)
{
  int nbr, i;

  lattice[site].state = 404;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state--;
  }
}



/* r14: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (500s) */
void Reactions::r14(lattice::Lattice lattice, int site)
{
  int nbr, i;

  lattice[site].state = 502;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state++;
  }
}



/* r15:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (500s) */
void Reactions::r15(lattice::Lattice lattice, int site)
{
  int nbr, i;

  lattice[site].state = 501;
  for (i = 0; i < 2; i++) {
    nbr = lattice[site].nbr[i];
    lattice[nbr].state--;
  }
}
