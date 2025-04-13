/* envrn.h: functions for checking the environment of a site */

#ifndef envrn_h
#define envrn_h

#include "lattice.hpp"

class Environment {
public:
    static int checkEnv(Lattice *lattice, int site);
    static int isActive(int site, Lattice *lattice, int rxn);

private:
    static int check100(Lattice *lattice, int site);
    static int check200(Lattice *lattice, int site);
    static int check300(Lattice *lattice, int site);
    static int check400(Lattice *lattice, int site);
    static int check500(Lattice *lattice, int site);
};

#endif
