# Legacy portability limitations

The curated 1999 source is preserved byte-for-byte. Modern execution exposed two
initialization-time out-of-bounds reads that depended on historical allocator
layout. They are evidence about the model as implemented, not changes to make
silently in the archive.

## Reproduction

The digest-pinned `linux/amd64` container built the unmodified source with GCC
12.2.0 and the forced standard-header compatibility flags. All five fixtures
terminated with `SIGSEGV` (`returncode: -11`) before writing outputs.

An AddressSanitizer diagnostic on `hotrox/935077498` first reported:

```text
ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 4
#0 readRxns /archive/c-model/source/rxnlist.c:93
#1 main /archive/c-model/source/mckaol.c:28
0 bytes to the right of 16-byte region
allocated by readRxns /archive/c-model/source/rxnlist.c:80
```

`rxnlist.c:91–93` uses the stale final desorption `nrts` value for both the
diffusion allocation and copy loop. In these fixtures that value is five, so the
first two diffusion copies read a fifth float from four-float source arrays.

After adding suffix-only allocation slack externally, AddressSanitizer exposed a
second independent pre-buffer read:

```text
ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 4
#0 terminateSurface /archive/c-model/source/lattice.c:271
#1 main /archive/c-model/source/mckaol.c:35
40 bytes to the left of the lattice allocation
allocated by makeLattice /archive/c-model/source/lattice.c:147
```

## Compatibility decision

`compat/historical_malloc_slack.c` is an external linker wrapper around `malloc`
and `free`. It supplies 64 aligned bytes before and after each allocation, making
both historical undefined accesses non-fatal without changing curated source
bytes. The harness verifies and records the shim SHA-256 and links it explicitly
with `--wrap=malloc --wrap=free`.

This is not a claim that the out-of-bounds reads are correct. The report keeps the
shim visible as compatibility provenance, and the defect is carried as a follow-on
issue for the Rust-port semantics decision. The compatibility shim exists only to
replay the historical executable oracle; it must not be copied into the Rust port.
