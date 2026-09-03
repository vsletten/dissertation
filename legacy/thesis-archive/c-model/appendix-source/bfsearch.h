/* bfsearch.h: exports routine to do breadth-first search 
 *             on a lattice                               */

#ifndef bfsearch_h
#define bfsearch_h

#include "lattice.h"

#define WHITE 0
#define GRAY 1
#define BLACK 2


/* BFS: perform breadth-first search to mark all "solid" nodes
        start in a known occupied node, s, and discover all 
        occupied neighbors. (discovered => BLACK)
        At the end, undiscovered, occupied nodes are unbonded
        clusters. */
void BFS(Lattice lattice, int s, int n);

#endif
