/* sim.h: exports struct and function for simulation parameters */

#ifndef sim_h
#define sim_h

class Simulation {
public:
  int nsteps,      /* total number of steps in simulation */
      wsteps,      /* number of steps between data writes */
      msteps,      /* number of steps between movie frames */
      drawbonds;   /* TRUE = draw bonds, FALSE = don't */
  long ranseed;    /* seed for initran2 */

  static Simulation *CreateSimulation(void);
  static void DisposeSimulation(Simulation *s);
};

#endif
