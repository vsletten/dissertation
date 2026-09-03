# 1999 C KMC conformance report

Generated: `2026-09-03T06:04:58.109620+00:00`

## Provenance

- Git commit: `unavailable; content manifests govern`
- Source manifest SHA-256: `73ee92b8cda415f490112bba398c1b7e5e8e010dad733fcbab8dd5ecc5e85482`
- Fixture-input manifest SHA-256: `3612d4136f72aedf32e5fbace72726c720b07fea1da1565b07ea87cec89f4c48`
- Platform: `Linux-6.8.0-117-generic-aarch64-with-glibc2.36`
- Architecture: `aarch64`
- libc: `glibc 2.36`
- Compiler: `gcc (Debian 12.2.0-14+deb12u1) 12.2.0`
- Compiler target: `aarch64-linux-gnu`
- Canonical container: `python:3.11.15-bookworm@sha256:a8f8fbe1a0edc9e4dddafa64ba73f7e04be7be5ebc23f332362e779e0a2e4e52`
- Compiler lock matched: `true`
- Compatibility flags: `-O3 -ffast-math -std=gnu89 -include stdlib.h -include unistd.h`
- Allocator compatibility shim SHA-256: `e224b3c9e5d5f731135a6216b63b3171041682da5091399666f08449521b30d5`
- Diffusion IDs 24–27: `pinned_disabled`

## Fixture outcomes

| fixture | outcome | results parity | first divergence | surfAl rows | surfSi rows | seconds |
|---|---|---|---:|---:|---:|---:|
| hotrox/935077498 | behavioral_mismatch | diverged | 97 | 317 | 194 | 1384.411 |
| hotrox/936930575 | behavioral_mismatch | diverged | 1 | 234 | 164 | 1528.028 |
| hotrox/937172019 | behavioral_mismatch | diverged | 1 | 235 | 185 | 1476.897 |
| jasper/933892971 | behavioral_mismatch | diverged | 1 | 234 | 164 | 1523.819 |
| jasper/935835145 | behavioral_mismatch | diverged | 38 | 345 | 218 | 1215.047 |

`byte_parity` is exact SHA-256 equality. A
`compiler_prng_drift_candidate` is explicitly non-canonical evidence: content
manifests match, the compiler/architecture differs, every schema and row count
matches, and at least ten initial trajectory rows are exact. It remains a
candidate rather than a normalized fact and never passes the canonical gate.
Anything else is a
`behavioral_mismatch` and fails the gate.

## Independent controls

- Historical identical-input Hotrox/Jasper outputs byte-equal: `true`.
- Late-row sabotage detected without being labeled drift: `true`.
