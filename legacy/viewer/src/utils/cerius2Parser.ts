// Super simple parser that only cares about Atoms and Bonds

// Output types
type Atom = {
  id: number;
  element: string;
  position: [number, number, number];
  label: string;
  charge: number;
};

type Bond = {
  from: number;
  to: number;
};

/**
 * Parse an MSI file to extract atoms and bonds
 */
export function parseModelFile(content: string): { atoms: Atom[], bonds: Bond[] } | null {
  try {
    // Remove comments and split by lines
    const lines = content.split('\n')
      .map(line => line.trim())
      .filter(line => line && !line.startsWith('#'));
    
    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    
    let currentObject: 'atom' | 'bond' | null = null;
    let currentId = 0;
    let currentElement = '';
    let currentPosition: [number, number, number] = [0, 0, 0];
    let currentLabel = '';
    let currentCharge = 0;
    let currentAtom1 = 0;
    let currentAtom2 = 0;
    
    // Process line by line
    for (const line of lines) {
      // Check for start of atom
      if (line.match(/^\(\d+\s+Atom/)) {
        // Save previous atom if there was one
        if (currentObject === 'atom') {
          atoms.push({
            id: currentId,
            element: currentElement,
            position: currentPosition,
            label: currentLabel,
            charge: currentCharge
          });
        }
        
        // Start a new atom
        currentObject = 'atom';
        currentId = parseInt(line.match(/^\((\d+)/)?.[1] || '0', 10);
        currentElement = '';
        currentPosition = [0, 0, 0];
        currentLabel = '';
        currentCharge = 0;
      }
      // Check for start of bond
      else if (line.match(/^\(\d+\s+Bond/)) {
        // Save previous atom or bond if there was one
        if (currentObject === 'atom') {
          atoms.push({
            id: currentId,
            element: currentElement,
            position: currentPosition,
            label: currentLabel,
            charge: currentCharge
          });
        } else if (currentObject === 'bond' && currentAtom1 > 0 && currentAtom2 > 0) {
          bonds.push({
            from: currentAtom1,
            to: currentAtom2
          });
        }
        
        // Start a new bond
        currentObject = 'bond';
        currentId = parseInt(line.match(/^\((\d+)/)?.[1] || '0', 10);
        currentAtom1 = 0;
        currentAtom2 = 0;
      }
      // Check for end of object
      else if (line === ')') {
        // Save previous atom or bond if there was one
        if (currentObject === 'atom') {
          atoms.push({
            id: currentId,
            element: currentElement,
            position: currentPosition,
            label: currentLabel,
            charge: currentCharge
          });
          currentObject = null;
        } else if (currentObject === 'bond' && currentAtom1 > 0 && currentAtom2 > 0) {
          bonds.push({
            from: currentAtom1,
            to: currentAtom2
          });
          currentObject = null;
        }
      }
      // Parse atom attributes
      else if (currentObject === 'atom') {
        // Parse ACL for element (can be "A C ACL" or "A I ACL")
        if (line.includes('ACL')) {
          const match = line.match(/"([^"]+)"/);
          if (match) {
            const aclValue = match[1];
            const elemMatch = aclValue.match(/\d+\s+([A-Za-z]+)/);
            if (elemMatch) {
              currentElement = elemMatch[1];
            }
          }
        }
        // Parse XYZ for position
        else if (line.includes('(A D XYZ')) {
          const match = line.match(/\(([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\)/);
          if (match) {
            currentPosition = [
              parseFloat(match[1]),
              parseFloat(match[2]),
              parseFloat(match[3])
            ];
          }
        }
        // Parse Label
        else if (line.includes('(A C Label')) {
          const match = line.match(/"([^"]+)"/);
          if (match) {
            currentLabel = match[1];
          }
        }
        // Parse Charge
        else if (line.includes('(A F Charge')) {
          const match = line.match(/Charge\s+([-\d.]+)/);
          if (match) {
            currentCharge = parseFloat(match[1]);
          }
        }
      }
      // Parse bond attributes
      else if (currentObject === 'bond') {
        // Parse Atom1
        if (line.includes('(A O Atom1')) {
          const match = line.match(/Atom1\s+(\d+)/);
          if (match) {
            currentAtom1 = parseInt(match[1], 10);
          }
        }
        // Parse Atom2
        else if (line.includes('(A O Atom2')) {
          const match = line.match(/Atom2\s+(\d+)/);
          if (match) {
            currentAtom2 = parseInt(match[1], 10);
          }
        }
      }
    }
    
    // Add the last atom or bond if there is one
    if (currentObject === 'atom') {
      atoms.push({
        id: currentId,
        element: currentElement,
        position: currentPosition,
        label: currentLabel,
        charge: currentCharge
      });
    } else if (currentObject === 'bond' && currentAtom1 > 0 && currentAtom2 > 0) {
      bonds.push({
        from: currentAtom1,
        to: currentAtom2
      });
    }
    
    return { atoms, bonds };
  } catch (error) {
    console.error('Error parsing MSI file:', error);
    return null;
  }
}

// For backward compatibility with the existing UI
export function extractAtoms(model: any) {
  // Check if model is already in the simple format
  if (model && Array.isArray(model.atoms) && Array.isArray(model.bonds)) {
    return model;
  }
  
  // Return empty arrays as fallback
  return { atoms: [], bonds: [] };
}