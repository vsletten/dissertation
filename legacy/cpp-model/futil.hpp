/* futil.h: exports functions for dealing with file I/O */

#ifndef futil_h
#define futil_h

#include <fstream>
#include <string>

class Futil {
public:
    static void EatComment(std::ifstream &f, char flag);
    static std::ifstream OpenInputFile(const std::string &name);
    static std::ofstream OpenOutputFile(const std::string &name);
    static int CloseFile(std::ifstream &f);
    static int CloseFile(std::ofstream &f);
};

#endif
