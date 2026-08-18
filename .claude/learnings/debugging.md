# Debugging

### Saddle searches slide to trivial saddles from uncrossed-ridge guesses (2026-08-18)
First live hydrolysis TS hunt: the r(Si-Ow) scan ended (1.9 A) with energy
still rising, its endpoint became the Sella guess, and the search relaxed
back to a 208i cm^-1 water-wag saddle — mechanically valid (1 imaginary
mode, IRC converges) but chemically wrong (dG‡ 2.5 vs ~120 kJ/mol lit.).
Diagnosis signature: the saddle's driven coordinate far past the guess and
both quick-IRC ends identical to the reactant complex. Fixes in quarry:
scan_to_maximum (extend until interior max) + escaped-channel guard in
the phase-1 driver. A verified saddle is not necessarily the right saddle.
