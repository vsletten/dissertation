#include <cmath>
#include <stdlib.h>
#include "futil.hpp"
#include "myerr.hpp"
#include "rxnlist.hpp"


void ReactionList::DisposeReactionList(ReactionList  *rl)
{
  int i;

  for (i = 0; i < NRXN; i++)
    delete[] rl->reactions[i].rate;
  delete[] rl->reactions;
  delete rl;
}


/* readRxns: read in rates, dE's for forward rxns
 *           calculate rates for reverse
 *           complete reactionList */
ReactionList *ReactionList::CreateReactionList(void)
{
  std::ifstream f;
  Reaction *reactions;
  int reactant, product, i, j, numRates;
  float kp, de, t, a, rt, dmsi, dmal;
  // float kt; // TODO: Not sure what this is for


  f = Futil::OpenInputFile("data.rxn");
  if (!f.is_open()) {
    Myerr::die("Error: Could not open data.rxn");
  }
  f >> t; 
  Futil::EatComment(f, '#');
  f >> dmsi; 
  Futil::EatComment(f, '#');
  f >> dmal; 
  Futil::EatComment(f, '#');
  reactions = new Reaction[NRXN];
  rt = R * t;
  // kt = Kb * t;

  /* hydrolysis reactions 
   * de is in kcal/mol/K */
  for (i = 0; i < NHYD; i += 2) {
    f >> reactant >> product; 
    Futil::EatComment(f, '#');
    reactions[i].reactant = reactant;
    reactions[i+1].reactant = product;
    f >> numRates; 
    Futil::EatComment(f, '#');
    reactions[i].nrates = reactions[i+1].nrates = numRates;
    reactions[i].rate = new float[numRates]; 
    reactions[i+1].rate = new float[numRates]; 
    for (j = 0; j < numRates; j++) {
      f >> kp >> de;
      reactions[i].rate[j] = kp;
      reactions[i+1].rate[j] = kp * exp(de / rt);
    }
  }

  /* adsorption : a is \nu * exp(\mu_solid / RT) 
   * dm's are in kcal/mol/K                      */
  Futil::EatComment(f, '#');
  for (i = NHYD; i < NADS; i++) {
    f >> reactant >> product; 
    Futil::EatComment(f, '#');
    f >> a >> de; 
    Futil::EatComment(f, '#');
    reactions[i].reactant = reactant;
    reactions[i].info = product;
    reactions[i].nrates = 1;
    reactions[i].rate = new float[1];
    if (i < 18)  /* si adsorption or al adsorption */
      reactions[i].rate[0] = a * exp(de / rt) * exp(dmsi / rt);
    else
      reactions[i].rate[0] = a * exp(de / rt) * exp(dmal / rt);
  }

  /* desorption */
  Futil::EatComment(f, '#');
  for (i = NADS; i < NDES; i++) {
    f >> reactant; Futil::EatComment(f, '#');
    reactions[i].reactant = reactant;
    f >> numRates;
    reactions[i].nrates = numRates;
    reactions[i].rate = new float[numRates];
    for (j = 0; j < numRates; j++) {
      f >> a >> de;
      reactions[i].rate[j] = a * exp(- de / rt);
    }
  }
  
  /* diffusion */
  for (i = NDES; i < NRXN; i++) {
    reactions[i].reactant = reactions[i-4].reactant;
    reactions[i].nrates = reactions[i-4].nrates;
    reactions[i].rate = new float[numRates];
    for (j = 0; j < numRates; j++)
      reactions[i].rate[j] = reactions[i-4].rate[j];
  }

  f.close();
  ReactionList *rl = new ReactionList(reactions);
  rl->dmsi = dmsi;
  rl->dmal = dmal;
  return rl;
}
