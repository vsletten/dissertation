/* bfsearch.h: exports routine to do breadth-first search 
 *             on a lattice                               */

#ifndef bfsearch_h
#define bfsearch_h

#include "lattice.h"

#define WHITE 0
#define GRAY 1
#define BLACK 2


void BFS(Lattice lattice, int s, int n);

#endif
