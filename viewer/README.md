# Cerius2 Molecular Viewer

A web-based 3D viewer for Cerius2 molecular model files (.msi), designed for visualizing kaolinite mineral surfaces from Monte Carlo simulations.

## Quick Start

```bash
npm install
npm run dev
```

Then open http://localhost:5173 and upload an MSI file.

## Features

- **3D Visualization**: Interactive molecular structure rendering using Three.js
- **File Support**: Cerius2 Model Structure Interface (.msi) files
- **Display Modes**:
  - Ball-and-stick (default)
  - Space-filling
  - Wireframe
- **Coloring**: Element-based (CPK) or charge-based
- **Controls**:
  - Left-click + drag: Rotate
  - Right-click + drag: Pan
  - Scroll: Zoom
- **Atom Labels**: Toggle element/label display
- **Model Info**: Composition, dimensions, and charge statistics

## Project Structure

```
viewer/
├── src/
│   ├── components/
│   │   ├── FileUploader.tsx   # Drag-and-drop file upload
│   │   ├── MoleculeViewer.tsx # Three.js 3D scene
│   │   └── ModelInfo.tsx      # Structure statistics panel
│   ├── utils/
│   │   ├── cerius2Parser.ts   # MSI file parser
│   │   └── elementColors.ts   # CPK colors and atomic radii
│   ├── types/
│   │   └── cerius2.ts         # TypeScript interfaces
│   ├── App.tsx                # Main application
│   └── main.tsx               # Entry point
├── public/
│   └── example_kaolinite.msi  # Sample model file
├── sample_data/               # Additional test files
└── msi-schema.md              # MSI file format documentation
```

## MSI File Format

The viewer parses Cerius2 MSI files with the following structure:

```
# MSI CERIUS2 DataModel File Version 2 0
(1 Model
 (2 Atom
  (A I ACL "14 Si")           # Element (atomic number + symbol)
  (A F Charge 1.1)            # Partial charge
  (A D XYZ (0.0 0.0 0.0))     # 3D coordinates
  (A C Label Si1)             # Atom label
 )
 (34 Bond
  (A O Atom1 2)               # First atom ID
  (A O Atom2 6)               # Second atom ID
 )
)
```

See `msi-schema.md` for detailed format documentation.

## Development

```bash
npm run dev      # Start dev server with hot reload
npm run build    # Production build
npm run lint     # Run ESLint
npm run preview  # Preview production build
```

## Technologies

- **React 19** - UI framework
- **TypeScript** - Type safety
- **Vite** - Build tool
- **Three.js** - 3D rendering
- **@react-three/fiber** - React renderer for Three.js
- **@react-three/drei** - Useful helpers (OrbitControls, Text, Billboard)

## Kaolinite Context

This viewer is designed for the dissertation project studying kaolinite (Al2Si2O5(OH)4) mineral surfaces. Kaolinite has a layered structure with:
- **Silica tetrahedral sheet**: Si atoms (tan) bonded to O atoms (red)
- **Alumina octahedral sheet**: Al atoms (grey) coordinated with O and OH groups
- **Interlayer**: Hydroxyl groups

The viewer supports visualizing output from the Monte Carlo simulation in the `model/` directory.

## License

ISC
