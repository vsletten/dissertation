//#include <stdio.h>
//#include "common.hpp"
#include "actions.hpp"
#include "ran2.hpp"
#include "myerr.hpp"
#include <cmath>
#include <iostream>

/* doEvent: randomly pick event and update state accordingly
 *          return time increment for the event */
float Actions::DoEvent(EventList *eventList)
{
  float ratesum, partsum, eps, dt;
  EventList *currentEvent;

  currentEvent = eventList;
  ratesum = 0;

  while (currentEvent->next) {
    ratesum += currentEvent->rate;
    currentEvent = currentEvent->next;
  }
  ratesum += currentEvent->rate;

  if (ratesum == 0.0)
    Myerr::die("ratesum is zero!");
  dt = - log(ran2()) / ratesum;

  eps = ran2();
  partsum = 0;
  currentEvent = eventList;
  while (currentEvent->next) {
    partsum += currentEvent->rate / ratesum;
    if (eps <= partsum)
      break;
    else
      currentEvent = currentEvent->next;
  }

  DoReaction(currentEvent->site, currentEvent->rxn);
  return dt;
}



/* doReaction: update state of site and its nbrs based on reaction rxn */
void Actions::DoReaction(int site, int rxn)
{
  switch (rxn) {
    case 0: this->DoReaction0(site); break;
    case 1: this->DoReaction1(site); break;
    case 2: this->DoReaction2(site); break;
    case 3: this->DoReaction3(site); break;
    case 4: this->DoReaction4(site); break;
    case 5: this->DoReaction5(site); break;
    case 6: this->DoReaction6(site); break;
    case 7: this->DoReaction7(site); break;
    case 8: this->DoReaction8(site); break;
    case 9: this->DoReaction9(site); break;
    case 10: this->DoReaction10(site); break;
    case 11: this->DoReaction11(site); break;
    case 12: this->DoReaction12(site); break;
    case 13: this->DoReaction13(site); break;
    case 14: this->DoReaction14(site); break;
    case 15: this->DoReaction15(site); break;
    case 16: this->AdsorbAl(site); break;
    case 17: this->AdsorbAl(site); break;
    case 18: this->AdsorbSi(site); break;
    case 19: this->AdsorbSi(site); break;
    case 20: this->DesorbAl(site); break;
    case 21: this->DesorbSi(site); break;
    case 22: this->DesorbSi(site); break;
    case 23: this->DesorbAl(site); break;
    default: this->Diffuse(site, rxn); 
             this->lattice->RemoveUnattachedClusters();
             break;
  }
}

/* adsorbAl: update state of O site upon adsorbing Al(OH,H2O)6 */
void Actions::AdsorbAl(int site)
{
  int i, nbr;

  if (ISAL(this->lattice->sites[site])) {
    this->lattice->sites[site].state = 107;
    for (i = 0; i < 6; i++) {
      nbr = this->lattice->sites[site].nbr[i];
      switch (this->lattice->sites[nbr].state) {
        case 503: this->lattice->sites[nbr].state = 502; break;
	case 500: this->lattice->sites[nbr].state = 503; break;
	case 407: this->lattice->sites[nbr].state = 403; break;
	case 409: this->lattice->sites[nbr].state = 405; break;
	case 408: this->lattice->sites[nbr].state = 407; break;
	case 400: this->lattice->sites[nbr].state = 409; break;
	case 406: 
	  this->lattice->sites[nbr].state = 410; 
	  this->lattice->sites[nbr].lostal = site;
	  break;
	default: 
	  std::cout << "state " << this->lattice->sites[nbr].state << std::endl;
	  Myerr::die("invalid state in adsorbAl");
	  break;
      }
    }
  } else {
    this->lattice->sites[site].state = 200 + WRONG;
  }
}



/* adsorbSi: update state of O site upon adsorbing Si(OH)4 */
void Actions::AdsorbSi(int site)
{
  int i, nbr;

  if (ISSI(this->lattice->sites[site])) {
    this-> lattice->sites[site].state = 205;
    for (i = 0; i < 4; i++) {
      nbr = this->lattice->sites[site].nbr[i];
      switch (this->lattice->sites[nbr].state) {
        case 300: this->lattice->sites[nbr].state = 303; break;
	case 303: this->lattice->sites[nbr].state = 302; break;
	case 404: this->lattice->sites[nbr].state = 402; break;
	case 405: this->lattice->sites[nbr].state = 403; break;
	case 409: this->lattice->sites[nbr].state = 407; break;
	case 400: this->lattice->sites[nbr].state = 408; break;
	default:
	  std::cout << "state " << this->lattice->sites[nbr].state << std::endl;
	  Myerr::die("invalid state in adsorbSi");
	  break;
      }
    }
  } else {
    this->lattice->sites[site].state = 100 + WRONG;
  }
}


/* desorbAl: update state of O site upon desorbing Al(OH,H2O)6 */
void Actions::DesorbAl(int site)
{
  int i, nbr;

  if (ISAL(this->lattice->sites[site])) {
    this->lattice->sites[site].state = 100;
    for (i = 0; i < 6; i++) {
      nbr = this->lattice->sites[site].nbr[i];             /* nbr is an o site */
      switch (this->lattice->sites[nbr].state) {
        case 502: this->lattice->sites[nbr].state = 503; break;
	      case 503: this->lattice->sites[nbr].state = 500; break;
	      case 403: this->lattice->sites[nbr].state = 407; break;
	      case 405: this->lattice->sites[nbr].state = 409; break;
	      case 407: this->lattice->sites[nbr].state = 408; break;
	      case 409: this->lattice->sites[nbr].state = 400; break;
	      case 410: 
	        this->lattice->sites[nbr].state = 406; 
	        this->lattice->sites[nbr].lostal = -1;
	        break;
	      default: break;
      }
    }
  } else {
    this->lattice->sites[site].state = 200;
  }
}



/* desorbSi: update state of O site upon desorbing Si(OH)4 */
void Actions::DesorbSi(int site)
{
  int i, nbr;

  if (ISSI(this->lattice->sites[site])) {
    this->lattice->sites[site].state = 200;
    for (i = 0; i < 4; i++) {
      nbr = this->lattice->sites[site].nbr[i];
      switch (this->lattice->sites[nbr].state) {
        case 303: this->lattice->sites[nbr].state = 300; break;
	case 302: this->lattice->sites[nbr].state = 303; break;
	case 402: this->lattice->sites[nbr].state = 404; break;
	case 403: this->lattice->sites[nbr].state = 405; break;
	case 407: this->lattice->sites[nbr].state = 409; break;
	case 408: this->lattice->sites[nbr].state = 400; break;
	default: break;
      }
    }
  } else {
    this->lattice->sites[site].state = 100;
  }
}




void Actions::Diffuse(int source, int target)
{
  if (this->lattice->sites[source].state == 107 || this->lattice->sites[source].state == 299) {
    DesorbAl(source);
    AdsorbAl(target);
  } else {
    DesorbSi(source);
    AdsorbSi(target);
  }
}



/* r0: =Si-O-Si= + H2O => 2(=Si-OH) */
void Actions::DoReaction0(int site)
{
  int si, i;

  this->lattice->sites[site].state = 302;
  for (i = 0;  i < 2; i++) {
    si = this->lattice->sites[site].nbr[i];
    this->lattice->sites[si].state++;
  }
}



/* r1: 2(=Si-OH) => =Si-O-Si= + H2O */
void Actions::DoReaction1(int site)
{
  int si, i;

  this->lattice->sites[site].state = 301;
  for (i = 0;  i < 2; i++) {
    si = this->lattice->sites[site].nbr[i];
    this->lattice->sites[si].state--;
  }
}



/* r2: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
void Actions::DoReaction2(int site)
{
  int nbr;

  this->lattice->sites[site].state = 402;
  nbr = this->lattice->sites[site].nbr[2];
  this->lattice->sites[nbr].state++;  /* increment for si only */
}



/* r3: =Si-OH + =Al-OH-Al= => =Si-O-2(Al=) + H2O */
void Actions::DoReaction3(int site)
{
  int nbr;

  this->lattice->sites[site].state = 401;
  nbr = this->lattice->sites[site].nbr[2];
  this->lattice->sites[nbr].state--;  /* decrement for si only */
}



/* r4: =Si-O-2(Al=) + H2O => =Si-OH-Al= + =Al-OH */
void Actions::DoReaction4(int site)
{
  int nbr;
  float r;

  this->lattice->sites[site].state = 410;
  r = ran2();
  if (r < 0.5)
    nbr = this->lattice->sites[site].nbr[0];
  else 
    nbr = this->lattice->sites[site].nbr[1];
  this->lattice->sites[nbr].state++;
  this->lattice->sites[site].lostal = nbr;
}



/* r5: =Si-OH-Al= + =Al-OH => =Si-O-2(Al=) + H2O */
void Actions::DoReaction5(int site)
{
  int nbr;

  this->lattice->sites[site].state = 401;
  nbr = this->lattice->sites[site].lostal;
  this->lattice->sites[nbr].state--;
  this->lattice->sites[site].lostal = -1;
}



/* r6: =Si-OH + =Al-OH-Al= + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Actions::DoReaction6(int site)
{
  int i, nbr;

  this->lattice->sites[site].state = 403;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state++;
  }
}



/* r7: =Si-OH + =Al-OH + =Al-H2O => =Si-OH + =Al-OH-Al= + H2O */
void Actions::DoReaction7(int site)
{
  int i, nbr;

  this->lattice->sites[site].state = 402;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state--;
  }
}



/* r8: =Si-OH-Al= + =Al-OH + H2O => =Si-OH + =Al-OH + =Al-H2O */
void Actions::DoReaction8(int site)
{
  int nbr;

  this->lattice->sites[site].state = 403;
  nbr = this->lattice->sites[site].nbr[2];   /* si */
  this->lattice->sites[nbr].state++;
  nbr = (this->lattice->sites[site].lostal == this->lattice->sites[site].nbr[0]) ? 
         this->lattice->sites[site].nbr[1] : this->lattice->sites[site].nbr[0];
  this->lattice->sites[nbr].state++;
  this->lattice->sites[site].lostal = -1;
}




/* r9: =Si-OH + =Al-OH + =Al-H2O => =Si-OH-Al= + =Al-OH + H2O */
void Actions::DoReaction9(int site)
{
  int nbr, lost;
  float r;

  this->lattice->sites[site].state = 410;
  r = ran2();
  nbr = this->lattice->sites[site].nbr[2];   /* si */
  this->lattice->sites[nbr].state--;
  if (r < 0.5) {
    nbr = this->lattice->sites[site].nbr[0];   /* al */
    lost = this->lattice->sites[site].nbr[1];
  } else {
    nbr = this->lattice->sites[site].nbr[1];
    lost = this->lattice->sites[site].nbr[0];
  }
  this->lattice->sites[nbr].state--;
  this->lattice->sites[site].lostal = lost;
}



/* r10: =Si-OH-Al= + H2O => =Si-OH + =Al=H2O */
void Actions::DoReaction10(int site)
{
  int nbr, al1, al2;

  this->lattice->sites[site].state = 407;
  nbr = this->lattice->sites[site].nbr[2];   /* si */
  this->lattice->sites[nbr].state++;
  al1 = this->lattice->sites[site].nbr[0];
  al2 = this->lattice->sites[site].nbr[1];
  nbr = (this->lattice->sites[al1].state == 100) ? al2 : al1;
  this->lattice->sites[nbr].state++;         /* al */
}



/* r11: =Si-OH + =Al=H2O => =Si-OH-Al= + H2O */
void Actions::DoReaction11(int site)
{
  int nbr, al1, al2;

  this->lattice->sites[site].state = 406;
  nbr = this->lattice->sites[site].nbr[2];
  this->lattice->sites[nbr].state--;
  al1 = this->lattice->sites[site].nbr[0];
  al2 = this->lattice->sites[site].nbr[1];
  nbr = (this->lattice->sites[al1].state == 100) ? al2 : al1;
  this->lattice->sites[nbr].state--;         /* al */
}



/* r12: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (400s) */
void Actions::DoReaction12(int site)
{
  int nbr, i;

  this->lattice->sites[site].state = 405;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state++;
  }
}



/* r13:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (400s) */
void Actions::DoReaction13(int site)
{
  int nbr, i;

  this->lattice->sites[site].state = 404;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state--;
  }
}



/* r14: =Al-OH-Al= + H2O => =Al-OH + =Al-H2O  (500s) */
void Actions::DoReaction14(int site)
{
  int nbr, i;

  this->lattice->sites[site].state = 502;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state++;
  }
}



/* r15:  =Al-OH + =Al-H2O => =Al-OH-Al= + H2O  (500s) */
void Actions::DoReaction15(int site)
{
  int nbr, i;

  this->lattice->sites[site].state = 501;
  for (i = 0; i < 2; i++) {
    nbr = this->lattice->sites[site].nbr[i];
    this->lattice->sites[nbr].state--;
  }
}
