#include <stdio.h>
#include <stdlib.h>
//#include "common.h"
#include "futil.h"
#include "rxnlist.h"
#include "math.h"

static float dmsi, dmal;

void free_rxnList(reactionList rl)
{
  int i;

  for (i = 0; i < NRXN; i++)
    free(rl[i].rate);
  free(rl);
}



/* readRxns: read in rates, dE's for forward rxns
 *           calculate rates for reverse
 *           complete reactionList */
reactionList readRxns(void)
{
  FILE *f;
  reactionList l;
  int rct, prd, i, j, nrts;
  float kp, de, t, a, rt, kt;


  f = openFile("data.rxn","r");
  fscanf(f, "%f", &t); eatComment(f, '#');
  fscanf(f, "%f", &dmsi); eatComment(f, '#');
  fscanf(f, "%f", &dmal); eatComment(f, '#');
  l = (reactionList) malloc (sizeof(*l) * NRXN);
  rt = R * t;
  kt = Kb * t;

  /* hydrolysis reactions 
   * de is in kcal/mol/K */
  for (i = 0; i < NHYD; i += 2) {
    fscanf(f, " %d %d", &rct, &prd); eatComment(f, '#');
    l[i].reactant = rct;
    l[i+1].reactant = prd;
    fscanf(f, "%d", &nrts); eatComment(f, '#');
    l[i].nrates = l[i+1].nrates = nrts;
    l[i].rate = (float *) malloc (sizeof(float) * nrts);
    l[i+1].rate = (float *) malloc (sizeof(float) * nrts);
    for (j = 0; j < nrts; j++) {
      fscanf(f, " %f %f", &kp, &de);
      l[i].rate[j] = kp;
      l[i+1].rate[j] = kp * exp(de / rt);
    }
  }

  /* adsorption : a is \nu * exp(\mu_solid / RT) 
   * dm's are in kcal/mol/K                      */
  eatComment(f, '#');
  for (i = NHYD; i < NADS; i++) {
    fscanf(f, " %d %d", &rct, &prd); eatComment(f, '#');
    fscanf(f, " %f %f", &a, &de); eatComment(f, '#');
    l[i].reactant = rct;
    l[i].info = prd;
    l[i].nrates = 1;
    l[i].rate = (float *) malloc (sizeof(float));
    if (i < 18)  /* si adsorption or al adsorption */
      l[i].rate[0] = a * exp(de / rt) * exp(dmsi / rt);
    else
      l[i].rate[0] = a * exp(de / rt) * exp(dmal / rt);
  }

  /* desorption */
  eatComment(f, '#');
  for (i = NADS; i < NDES; i++) {
    fscanf(f, " %d", &rct); eatComment(f, '#');
    l[i].reactant = rct;
    fscanf(f, " %d", &nrts);
    l[i].nrates = nrts;
    l[i].rate = (float *) malloc (sizeof(float) * nrts);
    for (j = 0; j < nrts; j++) {
      fscanf(f," %f %f", &a, &de);
      l[i].rate[j] = a * exp(- de / rt);
    }
  }
  
  /* diffusion */
  for (i = NDES; i < NRXN; i++) {
    l[i].reactant = l[i-4].reactant;
    l[i].nrates = l[i-4].nrates;
    l[i].rate = (float *) malloc (sizeof(float) * nrts);
    for (j = 0; j < nrts; j++)
      l[i].rate[j] = l[i-4].rate[j];
  }
  

  fclose(f);
  return l;
}



void getChem(float *si, float *al)
{
  *si = dmsi;
  *al = dmal;
}
