/* futil.h: exports functions for dealing with file I/O */

#ifndef futil_h
#define futil_h

#include <stdio.h>

class Futil {
public:
    static void eatComment(FILE *f, char flag);
    static FILE *openFile(const char *name, const char *mode);
    static int closeFile(FILE *f);
};

#endif
