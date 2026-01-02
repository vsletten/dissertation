# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a dissertation project studying kaolin mineral surface reactions through Monte Carlo simulation and visualization. The repository contains two main components:

1. **model/** - C++ kinetic Monte Carlo (KMC) simulation for modeling surface reactions on kaolin minerals
2. **viewer/** - React/TypeScript web application for visualizing molecular structures from simulation output

## Build, Test, and Development Commands

### Model (C++ Simulation)
```bash
cd model
make                # Build the mckaol executable
make clean          # Remove build artifacts
./mckaol            # Run simulation (requires input files)
```

### Viewer (Web Application)
```bash
cd viewer
npm install         # Install dependencies
npm run dev         # Start development server at http://localhost:5173
npm run build       # Build for production (TypeScript compilation + Vite build)
npm run lint        # Run ESLint on TypeScript/TSX files
npm run preview     # Preview production build
```

## Model Architecture

The Monte Carlo simulation models chemical reactions at the atomic level for alumina and silica sites on mineral surfaces using a lattice-based approach.

### Core Data Flow
1. **Initialization** (mckaol.cpp:20-98): Load simulation parameters → create unit cell → build lattice → populate solid structure → terminate surface
2. **Simulation Loop** (mckaol.cpp:101-139): For each step, generate event list → select and execute random event → update time → write output
3. **Output** (mckaol.cpp:140-143): Write final state files (MSI, XYZ, surface analysis)

### Key Components

**Lattice System** (lattice.hpp/cpp):
- LatticeSite stores state (100-599 encoding cation type and coordination), unit cell coordinates, neighbor list, and bond pair information
- State encoding: first digit = site type (1=Al, 2=Si, 3=Si-O-Si, 4=Si-O-Al2, 5=Al-OH-Al), last two digits = coordination/occupancy
- Lattice class manages 2D array of sites with periodic boundary conditions

**Reaction System**:
- ReactionList (rxnlist.hpp/cpp): Loads 28 reaction types with temperature-dependent rates from data.rxn
- Environment (envrn.hpp/cpp): Checks local environment (next-neighbor states) for each site to determine valid reactions
- EventList (evtlist.hpp/cpp): Linked list of possible events with specific rates based on local environment

**Execution** (actions.hpp/cpp):
- DoEvent() selects random event weighted by rate, executes state changes, returns time increment
- 16 hydrolysis reactions (R0-R15) plus adsorption/desorption for Si and Al species
- Reactions update site states and propagate changes to neighbors

### Site State Classification
See common.hpp:1-160 for complete scheme. Key categories:
- 100s: Al sites (100=empty, 101-107=Al with varying OH/H2O coordination)
- 200s: Si sites (200=empty, 201-205=Si with varying OH coordination)
- 300s: Si-O-Si bridging oxygen sites
- 400s: Si-O-Al2 bridging sites (most complex, 10 states)
- 500s: Al-OH-Al bridging sites

### Input Files (all in model/)
- **data.sim**: MC steps, output frequency, random seed
- **data.rxn**: Temperature, chemical potentials, rate constants for 28 reactions
- **data.cell**: Unit cell atomic positions and connectivity
- **data.lattice**: Lattice dimensions (a_cells, b_cells), surface plane orientation

### Output Files
- **start.msi, end.msi**: Cerius2 Model Structure Interface format (initial/final structures)
- **start.xyz**: XYZ Cartesian coordinates
- **results.dat**: Time series of site type populations
- **surfSi.out, surfAl.out**: Surface composition analysis for Si and Al sites
- **stepN.msi**: Optional intermediate snapshots at intervals

### Known Issues
- mckaol.cpp contains merge conflict markers (lines 5-70) that need resolution
- Code has both old and new initialization paths visible in conflict

## Viewer Architecture

React application using Three.js (via @react-three/fiber) to render molecular structures from simulation output files.

### Core Components

**FileUploader** (components/FileUploader.tsx): Accepts .msi, .car, .pdb, .xyz files via drag-and-drop or file picker

**MoleculeViewer** (components/MoleculeViewer.tsx): Three.js scene with interactive camera controls (orbit, zoom, pan), renders atoms as spheres and bonds as cylinders

**ModelInfo** (components/ModelInfo.tsx): Displays structure statistics, element counts, atom/bond information

### Data Pipeline
1. User uploads file → FileUploader reads content
2. App.tsx calls cerius2Parser.parseModelFile()
3. Parser extracts atoms (id, element, position, label, charge) and bonds (from, to)
4. Data passed to MoleculeViewer for 3D rendering
5. ModelInfo displays structural metadata

### File Format Support
Parser handles Cerius2 .msi files (see msi-schema.md for format details). The MSI format uses keyword-value pairs with sections like (1 Atom), (1 Bond) for defining molecular structure.

### Styling
Element colors defined in utils/elementColors.ts, supports charge-based and element-based coloring schemes

## Development Notes

### Model
- Uses C++11 standard with aggressive optimization (-O3 -ffast-math) in production builds
- Comment out optimization and add -ggdb in Makefile for debugging
- Breadth-first search (bfsearch.hpp/cpp) used to identify and remove unattached atomic clusters
- Random number generation via ran2.cpp (Numerical Recipes algorithm)

### Viewer
- Uses Vite for fast development and optimized builds
- TypeScript strict mode enabled
- React 19 with functional components and hooks
- Three.js r174+ with @react-three/fiber for declarative 3D
- No tests currently configured (npm test will fail)

### Code Style (from viewer/CLAUDE.md)
- Imports: Group as (1) React/framework (2) External libraries (3) Internal modules
- TypeScript: Prefer interfaces over types, explicit return types on functions
- Naming: camelCase for variables/functions, PascalCase for classes/components
- Use Prettier defaults for formatting
