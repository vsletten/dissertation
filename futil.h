/* futil.h: exports functions for dealing with file I/O */

#ifndef futil_h
#define futil_h

#include <stdio.h>

void eatComment(FILE *f, char flag);
FILE *openFile(const char *name, const char *mode);
int closeFile(FILE *f);

#endif
