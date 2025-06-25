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

class Reaction {
public:
  int reactant, /* reactant state */
      info,     /* multi-purpose info */
      nrates;   /* number of actual rates */
  float *rate;  /* reaction rate list since actual rate
                 * depends on next neighbors' states */
};

class ReactionList {
private:
  float dmsi, dmal;
  Reaction *reactions = nullptr;

  ReactionList(Reaction *reactions) : reactions(reactions) {}
  
public:
  static ReactionList *CreateReactionList(void);
  static void DisposeReactionList(ReactionList *&rl);

  float GetSiPotential(void) { return this->dmsi; }
  float GetAlPotential(void) { return this->dmal; }
  Reaction *GetReactions(void) { return this->reactions; }
};

#endif
