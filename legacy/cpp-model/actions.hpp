/* actions.h: exports functions for doing something */

#ifndef actions_h
#define actions_h

#include "lattice.hpp"
#include "evtlist.hpp"

class Actions {
private:
  Lattice *lattice = nullptr;

  bool DoReaction(int site, int rxn);

/* adsorbAl: update state of O site upon adsorbing Al(OH,H2O)6 */
  bool AdsorbAl(int site);

/* adsorbSi: update state of O site upon adsorbing Si(OH)4 */
  bool AdsorbSi(int site);

/* desorbAl: update state of O site upon desorbing Al(OH,H2O)6 */
  void DesorbAl(int site);

/* desorbSi: update state of O site upon desorbing Si(OH)4 */
  void DesorbSi(int site);

/* diffuse: update state of O site upon diffusing */
  void Diffuse(int source, int target);

/* r0: =Si-O-Si= + H2O => 2(=Si-OH) */
  void DoReaction0(int site);

/* r1: 2(=Si-OH) => =Si-O-Si= + H2O */
  void DoReaction1(int site);

/* r2: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction2(int site);

/* r3: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction3(int site);

/* r4: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction4(int site);

/* r5: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction5(int site);

/* r6: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction6(int site);

/* r7: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction7(int site);

/* r8: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction8(int site);

/* r9: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction9(int site);

/* r10: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction10(int site);

/* r11: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction11(int site);

/* r12: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction12(int site);

/* r13: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction13(int site);

/* r14: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction14(int site);

/* r15: =Si-O-2(Al=) + H2O => =Si-OH + =Al-OH-Al= */
  void DoReaction15(int site);

public:
    Actions(Lattice *lattice) : lattice(lattice) { }

/* doEvent: randomly pick event and update lattice state accordingly
 *          return time increment for the event */
    float DoEvent(EventList *eventList);
};


#endif
