/* sim.h: exports struct and function for simulation parameters */

#ifndef sim_h
#define sim_h

class sim {
public:
  struct simulation {
    int nsteps,      /* total number of steps in simulation */
      wsteps,      /* number of steps between data writes */
      msteps,      /* number of steps between movie frames */
      drawbonds;   /* TRUE = draw bonds, FALSE = don't */
  long ranseed;    /* seed for initran2 */
  };
  typedef struct simulation *Simulation;

  static Simulation readCond(void);
  static void free_Simulation(Simulation s);
};

#endif
