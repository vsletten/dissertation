import { describe, expect, it } from 'vitest';

import { loadMsi } from '../src/loaders/msi';

/** A faithful slice of the KMC golden start.msi grammar (kmc-rs data/golden). */
const MSI_FIXTURE = `# MSI CERIUS2 DataModel File Version 3 5
(1 Model
(A I Id 1)
 (A C Label "start")
 (2 Atom
  (A C ACL "13 Al")
  (A D XYZ (1.17181 4.38134 2.84548))
  (A I Id 2)
  (A C Label "Al0")
  (A I LabelType 0)
 )
 (3 Atom
  (A C ACL "14 Si")
  (A D XYZ (1.30217 7.28290 2.58231))
  (A I Id 3)
  (A C Label "Si1")
 )
 (4 Atom
  (A C ACL "8 O")
  (A D XYZ (0.24775 3.10361 0.225645))
  (A I Id 4)
  (A C Label "O2")
  (A F Charge -0.5)
 )
 (1002 Bond
  (A O Atom1 2)
  (A O Atom2 3)
 )
 (1003 Bond
  (A O Atom1 3)
  (A O Atom2 4)
 )
 (1004 Bond
  (A O Atom1 2)
  (A O Atom2 999)
 )
)
`;

describe('msi loader', () => {
  it('parses atoms with element, position, label, charge', () => {
    const g = loadMsi(MSI_FIXTURE);
    expect(g.count).toBe(3);
    expect(Array.from(g.ids)).toEqual([2, 3, 4]);
    expect(g.typeDict).toEqual(['Al', 'Si', 'O']);
    expect(g.positions![0]).toBeCloseTo(1.17181, 4);
    expect(g.positions![5]).toBeCloseTo(2.58231, 4);
    expect((g.nodeColumns['label'].data as string[])[0]).toBe('Al0');
    expect((g.nodeColumns['charge'].data as Float32Array)[2]).toBeCloseTo(-0.5, 5);
  });

  it('remaps bond object-ids to node indices and drops dangling bonds', () => {
    const g = loadMsi(MSI_FIXTURE);
    // bond 1004 references missing atom 999 → dropped
    expect(g.edgeCount).toBe(2);
    expect(Array.from(g.edgeSrc)).toEqual([0, 1]);
    expect(Array.from(g.edgeDst)).toEqual([1, 2]);
  });

  it('returns an empty graph for non-MSI text', () => {
    const g = loadMsi('hello\nworld\n');
    expect(g.count).toBe(0);
    expect(g.edgeCount).toBe(0);
  });
});
