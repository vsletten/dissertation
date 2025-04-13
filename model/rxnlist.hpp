/* rxnlist.h: exports types and functions for reaction list */

#ifndef rxnlist_h
#define rxnlist_h

#define N300 2
#define N400 14
#define NHYD 16
#define NADS 20
#define NDES 24
#define NRXN 28
#define Kb 1.38066e-23 /* boltzman in J/mol/K */
#define R 1.987e-3     /* R in kcal/mol/K */
#define JtoCAL 4.1842  /* multiply cal by JtoCAL to get J */

class rxnlist {
public:
  /* reaction: Gives reactant and product states for given reaction,
   *           and reaction rate. */
  struct reaction {
    int reactant, /* reactant state */
        info,     /* multi-purpose info */
        nrates;   /* number of actual rates */
    float *rate;  /* reaction rate list since actual rate
                   * depends on next neighbors' states */
  };
  typedef struct reaction *reactionList;

  static reactionList readRxns(void);
  static void free_rxnList(reactionList rl);
  static void getChem(float *si, float *al);
};

#endif
