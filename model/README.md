# Monte Carlo Simulation for Kaolin Mineral Surface Reactions

A kinetic Monte Carlo simulation tool for modeling surface reactions and transformations on kaolin minerals.

## Overview

This software simulates chemical reactions at the atomic level for alumina and silica sites on mineral surfaces. It focuses on processes such as hydrolysis, adsorption, and desorption of Si and Al atoms in aluminum-silicate systems using a stochastic Monte Carlo approach.

## Features

- Lattice-based representation of mineral surfaces
- Realistic reaction kinetics based on temperature and chemical potentials
- Tracking of site transformations over time
- Visualization output in MSI and XYZ formats
- Surface composition analysis

## Building

```bash
make
```

The executable `mckaol` will be created.

## Usage

1. Prepare input files (see Input section below)
2. Run the simulation:
   ```bash
   ./mckaol
   ```

## Input Files

The simulation requires four input files:

- **data.sim**: Simulation parameters (MC steps, output frequency, random seed)
- **data.rxn**: Reaction parameters (temperature, chemical potentials, rate constants)
- **data.cell**: Unit cell configuration (atomic positions, connectivity)
- **data.lattice**: Lattice parameters (dimensions, surface plane)

Template input files are provided in the `input/` directory.

## Output Files

- **results.dat**: Time series data of site type counts
- **start.msi, end.msi**: Initial and final structure files
- **start.xyz**: XYZ format file of atomic positions
- **surfSi.out, surfAl.out**: Surface analysis data

## Code Structure

- **mckaol.c**: Main program
- **lattice.c/h**: Lattice data structure
- **reactions.c/h**: Chemical reaction definitions
- **rxnlist.c/h**: Reaction list management
- **evtlist.c/h**: Monte Carlo event list
- **output.c/h**: Output file generation
- **sim.c/h**: Simulation parameter handling
- **ucell.c/h**: Unit cell configuration
- **actions.c/h**: Lattice modification actions
- **bfsearch.c/h**: Breadth-first search for connected components

## License

Copyright information for this code should be added here.

## Citation

If using this code in published work, please cite: [Citation information should be added here]