// Test file for the Cerius2 parser
import { parseModelFile, extractAtoms } from './cerius2Parser';

// Simple methane example
const methaneContent = `# MSI CERIUS2 DataModel File Version 2 0
(1 Model
 (A C Label methane)
 (2 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A C Label H1)
 )
 (3 Atom
  (A I ACL "6 C")
  (A F Charge -0.11)
  (A D XYZ (1.087 0 0))
  (A C Label C2)
 )
 (4 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.4493 1.02483 0))
  (A C Label H3)
 )
 (5 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.44936 -0.51245 -0.88749))
  (A C Label H4)
 )
 (6 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.44934 -0.51236 0.88756))
  (A C Label H5)
 )
 (7 Bond
  (A O Atom1 2)
  (A O Atom2 3)
 )
 (8 Bond
  (A O Atom1 3)
  (A O Atom2 4)
 )
 (9 Bond
  (A O Atom1 3)
  (A O Atom2 5)
 )
 (10 Bond
  (A O Atom1 3)
  (A O Atom2 6)
 )
)`;

// Run this via node with ts-node
export function runTest() {
  // Test the parser
  const model = parseModelFile(methaneContent);
  console.log('Parsed model:', model ? 'Success' : 'Failed');

  if (model) {
    // Extract atoms and bonds
    const { atoms, bonds } = extractAtoms(model);
    console.log('Extracted atoms:', atoms.length);
    console.log('Extracted bonds:', bonds.length);
    
    // Output the first atom and bond
    console.log('First atom:', atoms[0]);
    console.log('First bond:', bonds[0]);
  }
}

// For TypeScript safety
declare const module: any;

// Run the test if this file is executed directly
if (typeof window === 'undefined') {
  runTest();
}