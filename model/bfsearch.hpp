/* bfsearch.h: exports routine to do breadth-first search
 *             on a lattice                               */

#ifndef bfsearch_h
#define bfsearch_h

#include "lattice.hpp"

#define WHITE 0
#define GRAY 1
#define BLACK 2

class BreadthFirstSearch {
public:
  static void ColorNodes(Lattice *lattice, int s, int n);

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
