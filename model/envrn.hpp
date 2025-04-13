/* envrn.h: functions for checking the environment of a site */

#ifndef envrn_h
#define envrn_h

#include "lattice.hpp"

class Environment {
public:
    static int checkEnv(lattice::Lattice lattice, int site);
    static int isActive(int site, lattice::Lattice lattice, int rxn);

private:
    static int check100(lattice::Lattice lattice, int site);
    static int check200(lattice::Lattice lattice, int site);
    static int check300(lattice::Lattice lattice, int site);
    static int check400(lattice::Lattice lattice, int site);
    static int check500(lattice::Lattice lattice, int site);
};

#endif
