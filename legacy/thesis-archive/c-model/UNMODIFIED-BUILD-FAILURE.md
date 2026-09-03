# Frozen unmodified-build failure

Observed on 2026-09-02 with Apple clang 21.0.0 on arm64 macOS.

```text
$ make clean all
rm -f *.o *~ *.Addr *.Counts *.pixie mon.out logfile
gcc -O3 -ffast-math -o mckaol mckaol.c actions.c bfsearch.c envrn.c evtlist.c futil.c lattice.c myerr.c output.c ran2.c reactions.c rxnlist.c sim.c ucell.c -lm
myerr.c:9:3: error: call to undeclared library function 'exit'
output.c:13:3: error: call to undeclared function 'unlink'
output.c:168:22: error: call to undeclared library function 'malloc'
output.c:244:3: error: call to undeclared library function 'free'
make: *** [mckaol] Error 1
```

The curated files are intentionally not patched. The conformance harness supplies
`-std=gnu89 -include stdlib.h -include unistd.h`; this is the smallest compatibility
layer and leaves the A8 byte-copy invariant intact. The historical optimization flags
`-O3 -ffast-math` remain unchanged.
