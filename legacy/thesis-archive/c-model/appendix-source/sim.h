/* sim.h: exports struct and function for simulation parameters */

#ifndef sim_h
#define sim_h

struct simulation {
  int nsteps,      /* total number of steps in simulation */
      wsteps,      /* number of steps between data writes */
      msteps,      /* number of steps between movie frames */
      drawbonds;   /* TRUE = draw bonds, FALSE = don't */
  long ranseed;    /* seed for initran2 */
};
typedef struct simulation *Simulation;

Simulation readCond(void);
void free_Simulation(Simulation s);

#endif
