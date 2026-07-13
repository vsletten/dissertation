/* bfsearch.h: exports routine to do breadth-first search
 *             on a lattice                               */

#ifndef bfsearch_h
#define bfsearch_h

#include "lattice.hpp"

class BreadthFirstSearch {
public:
  static void ColorNodes(Lattice *lattice, int startSiteIdx);

private:
  /* Qnode: node of a queue used in BFS for clusters */
  struct Qnode {
    int val;
    struct Qnode *next;
  };

  static int deQ(struct Qnode **q);
  static void enQ(struct Qnode **q, int val);
  static int headQ(struct Qnode *q);
};

#endif
