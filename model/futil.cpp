#include <ctype.h>
#include "common.hpp"
#include "futil.hpp"
#include "myerr.hpp"

/* eatComment: advance file pointer to beginning of next line 
 *             flag is first character of all comments       */
void Futil::eatComment(FILE *f, char flag)
{
  int i;
  char c;
  
  i = TRUE;
  while (i) {              /* do this for multiline comments */
    while (isspace(c = fgetc(f)))  /* eat leading whitespace */
      ;
    if (c == flag) {
      while ((c = fgetc(f)) != '\n')  /* comment, read to end of line */
	;
    } else {
      ungetc(c, f);
      i = FALSE;
    }
  }
}


/* openFile: try to open a file die if unsuccessful */
FILE *Futil::openFile(const char *name, const char *mode)
{
  FILE *f;
  char msg[100];

  if (!(f = fopen(name, mode))) {
    sprintf(msg, "could not open %s with mode %s", name, mode);
    Myerr::die(msg);
  }
  return f;
}


int Futil::closeFile(FILE *f)
{
  return (fclose(f));
}


