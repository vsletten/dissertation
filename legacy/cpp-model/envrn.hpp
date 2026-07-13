/* envrn.h: functions for checking the environment of a site */

#ifndef envrn_h
#define envrn_h

#include "lattice.hpp"

class Environment {
public:
    Environment(Lattice *lattice) : lattice(lattice) { }

    int CheckEnv(int site);
    int IsActive(int site, int rxn);

private:
    Lattice *lattice;

    int Check100(int site);
    int Check200(int site);
    int Check300(int site);
    int Check400(int site);
    int Check500(int site);
};

#endif
