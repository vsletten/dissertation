#include "bfsearch.hpp"
#include "lattice.hpp"
#include "myerr.hpp"
#include <cmath>
#include <cstddef>
#include <stdlib.h>

/* BFS: perform breadth-first search to mark all "solid" nodes
        start in a known occupied node and discover all
        occupied neighbors. (discovered => BLACK)
        At the end, undiscovered, occupied nodes are unbonded
        clusters. */
void BreadthFirstSearch::ColorNodes(Lattice *lattice, int startSiteIdx) {
  int u, v, i;
  struct Qnode *currSite = NULL;

  for (u = 0; u < lattice->GetNsites(); u++) {
    lattice->sites[u].color = UNREACHABLE;
  }

  lattice->sites[startSiteIdx].color = ENQUEUED;
  enQ(&currSite, startSiteIdx);
  while (currSite) {
    u = headQ(currSite);
    for (i = 0; i < 6; i++) {
      v = lattice->sites[u].nbr[i];
      if (!EXISTS(v) ||              /* no ith neighbor */
          !ISOCC(lattice->sites[v])) /* unoccupied */
        continue;
      else if (lattice->sites[v].color == UNREACHABLE) {
        lattice->sites[v].color = ENQUEUED;
        enQ(&currSite, v);
      }
    }
    deQ(&currSite);
    lattice->sites[u].color = DISCOVERED;
  }
}

/* enQ: add val to tail of queue */
void BreadthFirstSearch::enQ(struct Qnode **q, int val) {
  struct BreadthFirstSearch::Qnode *newNode = new Qnode();

  if (newNode == nullptr)
    Myerr::die("malloc failed");
  newNode->val = val;
  if (nullptr == *q) {
    newNode->next = newNode; /* single el't; tail = head  */
  } else {
    newNode->next = (*q)->next; /* tail points to head */
    (*q)->next = newNode;       /* old tail points to new tail */
  }
  *q = newNode; /* queue points to new tail */
}

/* headQ: return val for head of queue */
int BreadthFirstSearch::headQ(struct BreadthFirstSearch::Qnode *q) {
  if (NULL == q) { /* empty queue */
    return -1;
  } else {
    return q->next->val;
  }
}

/* deQ: remove head of queue */
int BreadthFirstSearch::deQ(struct BreadthFirstSearch::Qnode **q) {
  int val;
  struct BreadthFirstSearch::Qnode *head;

  if (q == nullptr || *q == nullptr) /* empty queue */
    return -1;
  head = (*q)->next;
  val = head->val;
  if (*q == head) { /* was only el't in queue, so */
    *q = nullptr;   /*  make queue empty */
  } else {
    (*q)->next = head->next;
  }
  delete head;
  return val;
}
