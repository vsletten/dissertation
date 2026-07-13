# CLAUDE.md - Guidelines for Agentic Coding Assistants

## Project Overview

This is a React/TypeScript web application for visualizing molecular structures from Cerius2 MSI files. It's part of a dissertation project studying kaolin mineral surface reactions.

## Build, Lint and Test Commands
- Build: `npm run build`
- Lint: `npm run lint`
- Dev server: `npm run dev` (runs on http://localhost:5173)
- Preview prod build: `npm run preview`

## Architecture

### Core Data Flow
1. User uploads .msi file via `FileUploader.tsx`
2. `App.tsx` calls `cerius2Parser.parseModelFile()` to extract atoms and bonds
3. Parsed data passed to `MoleculeViewer.tsx` for 3D rendering
4. `ModelInfo.tsx` displays metadata (composition, dimensions, charges)

### Key Components

**MoleculeViewer.tsx** - Three.js/R3F scene with:
- `AtomSphere` - Renders atoms as colored spheres with element-based sizing
- `Bond` - Renders bonds as cylinders using quaternion rotation for proper orientation
- `TextLabel` - Billboard labels using drei's `Text` component
- Display modes: ball-and-stick, space-filling, wireframe
- Color schemes: element-based (CPK) or charge-based

**cerius2Parser.ts** - Simple line-by-line MSI parser that extracts:
- Atoms: id, element (from ACL field), position (XYZ), label, charge
- Bonds: from/to atom IDs (Atom1, Atom2 fields)

**elementColors.ts** - CPK coloring scheme and van der Waals radii for elements

### MSI File Format Notes
- Atoms defined as `(N Atom` blocks with attributes like `(A I ACL "14 Si")`
- The ACL attribute type can be 'I' (Integer) or 'C' (Character)
- Bonds defined as `(N Bond` blocks with `(A O Atom1 X)` and `(A O Atom2 Y)`
- See `msi-schema.md` for detailed format documentation

## Code Style Guidelines
- **Formatting**: Use Prettier with default settings
- **Imports**: Group imports (1. React/framework 2. External libraries 3. Internal modules)
- **Types**: Prefer TypeScript interfaces over types, explicit return types on functions
- **Naming**: camelCase for variables/functions, PascalCase for classes/components
- **Components**: Functional components with React hooks
- **3D Rendering**: Use @react-three/fiber and @react-three/drei for Three.js integration

## Common Issues and Solutions

### Bond cylinders not connecting atoms
The cylinder must be rotated using quaternions, not Euler angles:
```typescript
const yAxis = new THREE.Vector3(0, 1, 0);
const quat = new THREE.Quaternion();
quat.setFromUnitVectors(yAxis, direction.clone().normalize());
```

### Element not being parsed from MSI files
The ACL attribute can have different type codes. Check for `line.includes('ACL')` not just `(A C ACL`.

### Labels not rendering
Use drei's `Billboard` + `Text` components instead of canvas textures for reliable text rendering.

## Dependencies
- React 19, TypeScript, Vite
- three, @react-three/fiber, @react-three/drei

Update this file as the project evolves.
