#include <stdlib.h>
#include "myerr.hpp"
#include "bfsearch.hpp"
#include "lattice.hpp"

/* BFS: perform breadth-first search to mark all "solid" nodes
        start in a known occupied node and discover all 
        occupied neighbors. (discovered => BLACK)
        At the end, undiscovered, occupied nodes are unbonded
        clusters. */
void BFSearch::BFS(lattice::Lattice lattice, int s, int n)
{
  int u, v, i;
  struct Qnode *Q = NULL;

  for (u = 0; u < n; u++)
    lattice[u].color = WHITE;
  lattice[s].color = GRAY;
  enQ(&Q, s);
  while (Q) {
    u = headQ(Q);
    for (i = 0; i < 6; i++) {
      v = lattice[u].nbr[i];
      if (!EXISTS(v) ||          /* no ith neighbor */
	  !ISOCC(lattice[v]))    /* unoccupied */
	continue;
      else if (lattice[v].color == WHITE) {
	lattice[v].color = GRAY;
	enQ(&Q, v);
      }
    }
    deQ(&Q);
    lattice[u].color = BLACK;
  }
}
  
/* enQ: add val to tail of queue */
void BFSearch::enQ(struct Qnode **q, int val)
{
  struct BFSearch::Qnode *newNode; 

  if (!(newNode = (struct Qnode *) malloc (sizeof(struct Qnode))))
    Myerr::die("malloc failed");
  newNode->val = val;
  if (NULL == *q) {
    newNode->next = newNode;                         /* single el't; tail = head  */
  } else { 
    newNode->next = (*q)->next;                  /* tail points to head */
    (*q)->next = newNode;                        /* old tail points to new tail */
  }
  *q = newNode;                                  /* queue points to new tail */ 
}


/* headQ: return val for head of queue */
int BFSearch::headQ(struct BFSearch::Qnode *q)
{
  if (NULL == q) {                           /* empty queue */
    return -1;
  } else {
    return q->next->val;
  }
}


/* deQ: remove head of queue */
int BFSearch::deQ(struct BFSearch::Qnode **q)
{
  int val;
  struct BFSearch::Qnode *head;

  if (NULL == *q)                          /* empty queue */
    return -1;
  head = (*q)->next;
  val = head->val;
  if (*q == head) {              /* was only el't in queue, so */
    *q = NULL;                             /*  make queue empty */
  } else {
    (*q)->next = head->next;
  }
  free(head);
  return val;
}
  

