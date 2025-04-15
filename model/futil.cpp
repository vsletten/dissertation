#include <cctype>
#include <iostream>
#include "common.hpp"
#include "futil.hpp"
#include "myerr.hpp"

using namespace std;

/* eatComment: advance file pointer to beginning of next line 
 *             flag is first character of all comments       */
void Futil::EatComment(ifstream &f, char flag)
{
  int i;
  char c;
  
  i = TRUE;
  while (i) {              /* do this for multiline comments */
    while (f.get(c) && isspace(c))  /* eat leading whitespace */
      ;
    if (c == flag) {
      while (f.get(c) && c != '\n')  /* comment, read to end of line */
	;
    } else {
      f.putback(c);
      i = FALSE;
    }
  }
}


/* openFile: try to open a file die if unsuccessful */
ifstream Futil::OpenInputFile(const string &name)
{
  ifstream f(name);

  if (!f) {
    string msg = "could not open " + name + " for reading";
    Myerr::die(msg.c_str());
  }
  return f;
}

ofstream Futil::OpenOutputFile(const string &name)
{
  ofstream f(name);

  if (!f) {
    string msg = "could not open " + name + " for writing";
    Myerr::die(msg.c_str());
  }
  return f;
}

int Futil::CloseFile(ifstream &f)
{
  f.close();
  return 0;
}

int Futil::CloseFile(ofstream &f)
{
  f.close();
  return 0;
}
