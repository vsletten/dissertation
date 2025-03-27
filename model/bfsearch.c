#include <stdlib.h>
#include "myerr.h"
//#include "common.h"
#include "bfsearch.h"
#include "lattice.h"

/* Qnode: node of a queue used in BFS for clusters */
struct Qnode {
  int val;
  struct Qnode *next;
};

static int deQ(struct Qnode **q);
static void enQ(struct Qnode **q, int val);
static int headQ(struct Qnode *q);

/* BFS: perform breadth-first search to mark all "solid" nodes
        start in a known occupied node and discover all 
        occupied neighbors. (discovered => BLACK)
        At the end, undiscovered, occupied nodes are unbonded
        clusters. */
void BFS(Lattice lattice, int s, int n)
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
void enQ(struct Qnode **q, int val)
{
  struct Qnode *new;

  if (!(new = (struct Qnode *) malloc (sizeof(struct Qnode))))
    die("malloc failed");
  new->val = val;
  if (NULL == *q) {
    new->next = new;                         /* single el't; tail = head  */
  } else { 
    new->next = (*q)->next;                  /* tail points to head */
    (*q)->next = new;                        /* old tail points to new tail */
  }
  *q = new;                                  /* queue points to new tail */ 
}


/* headQ: return val for head of queue */
int headQ(struct Qnode *q)
{
  if (NULL == q) {                           /* empty queue */
    return -1;
  } else {
    return q->next->val;
  }
}


/* deQ: remove head of queue */
int deQ(struct Qnode **q)
{
  int val;
  struct Qnode *head;

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
  

