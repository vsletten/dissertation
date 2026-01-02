# Kinetic Monte Carlo Simulation for Kaolin Mineral Surface Reactions

A stochastic simulation tool for modeling surface reactions and transformations on kaolin minerals at the atomic level.

## Overview

This software implements a kinetic Monte Carlo (KMC) simulation using the Gillespie algorithm to model chemical reactions on aluminum-silicate mineral surfaces. It tracks hydrolysis, adsorption, and desorption processes for alumina (Al) and silica (Si) sites using a lattice-based representation with environment-dependent reaction rates.

### Key Features

- Lattice-based representation of mineral surfaces with periodic boundary conditions
- 28 reaction types including hydrolysis, adsorption, desorption, and diffusion
- Temperature-dependent kinetics via Boltzmann rate calculations
- Environment-dependent rate selection based on local neighbor configurations
- Multiple output formats: MSI (Cerius2), XYZ, and surface composition analysis

## Building

### Requirements
- C++ compiler with C++11 support (g++ recommended)
- GNU Make

### Compilation
```bash
make
```

This creates the `mckaol` executable with aggressive optimization (`-O3 -ffast-math`).

### Debug Build
Edit the Makefile to enable debug symbols:
```makefile
# Comment out production flags:
# CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -ffast-math

# Uncomment debug flags:
CXXFLAGS = -std=c++11 -Wall -Wextra -ggdb
```

Then rebuild:
```bash
make clean && make
```

## Usage

1. Prepare the four required input files (see [Input Files](#input-files))
2. Run the simulation:
   ```bash
   ./mckaol
   ```
3. Results are written to output files in the current directory

## Input Files

The simulation requires four data files in the working directory:

### data.sim - Simulation Parameters
```
5000000     # Number of Monte Carlo steps
1000        # Data output frequency (0 = no intermediate output)
1000000     # MSI snapshot frequency (0 = none, -1 = final only)
-2          # Random seed (negative value initializes the generator)
1           # Draw bonds: 0 = only occupied atoms, 1 = full lattice structure
```

### data.lattice - Lattice Configuration
```
20 3 0      # Number of unit cells: a-direction, b-direction, surface plane
```

The surface plane parameter controls which direction has open boundaries:
- `0`: ac surface (b-direction has open boundary)
- `1`: bc surface (a-direction has open boundary)

### data.rxn - Reaction Parameters and Rate Constants
```
8000.0      # Temperature (Kelvin)
-1.0        # Chemical potential delta for Si (kcal/mol)
-1.0        # Chemical potential delta for Al (kcal/mol)

# Hydrolysis reaction pairs (forward and reverse)
301 302     # Reactant state, Product state
15          # Number of environment-dependent rate variants
1 2.6 1 2.6 ...  # Pre-exponential (k+) and activation energy (dE) pairs

# Additional reactions follow same pattern...
```

The chemical potentials control initial lattice filling:
- Positive values (supersaturated): ~30% filled
- Near zero (equilibrium): ~50% filled
- Negative values (undersaturated): ~70% filled

### data.cell - Unit Cell Definition
```
5.140 8.930 7.370           # Cell dimensions a, b, c (angstroms)
-0.03141 -0.25038 0.0       # Cell angles alpha, beta, gamma (radians)
26                          # Number of atomic positions

# For each position:
# position_id  site_type  x  y  z
0 100 1.76027 4.49982 2.84548
# Followed by 6 neighbor lines:
# neighbor_id  delta_a  delta_b  delta_c
  8   0  0  0
  9   0  0  0
  ...
```

Site types:
- `100`: Aluminum site
- `200`: Silicon site
- `300`: Si-O-Si bridging oxygen
- `400`: Si-O-Al2 bridging oxygen
- `500`: Al-OH-Al bridging oxygen

## Output Files

### results.dat
Comma-separated time series of site populations:
```
time, Si_total, Si(OH)0, Si(OH)1, Si(OH)2, Si(OH)3, Si(OH)4, Al_total, Al(OH,H2O)0, ..., Al(OH,H2O)6
```

### Structure Files
- **start.msi**: Initial structure in Cerius2 MSI format
- **end.msi**: Final structure in MSI format
- **start.xyz**: Initial structure in XYZ format
- **stepN.msi**: Intermediate snapshots (if msteps > 0)

### Surface Analysis
- **surfSi.out**: 2D surface coordinates of Si atoms (x, y)
- **surfAl.out**: 2D surface coordinates of Al atoms (x, y)

## Site State Encoding

The simulation uses a three-digit integer encoding for site states:

| Range | Type | Description |
|-------|------|-------------|
| 100-107 | Al | Al site with 0-6 OH/H2O ligands |
| 199 | Al | Si adsorbed on Al site |
| 200-205 | Si | Si site with 0-4 OH groups |
| 299 | Si | Al adsorbed on Si site |
| 300-303 | Si-O-Si | Bridging oxygen between Si atoms |
| 400-410 | Si-O-Al2 | Bridging oxygen between Si and 2 Al |
| 500-503 | Al-OH-Al | Bridging oxygen between Al atoms |

See `common.hpp` for the complete state encoding scheme.

## Reaction Types

### Hydrolysis Reactions (0-15)
Eight forward/reverse reaction pairs for bond breaking and formation:
- **R0/R1**: Si-O-Si bond hydrolysis/formation
- **R2-R15**: Si-O-Al and Al-O-Al bond transformations

### Adsorption (16-19)
- Al species adsorption onto Al sites
- Si species adsorption onto Si sites

### Desorption (20-23)
- Al species desorption from Al sites
- Si species desorption from Si sites

## Algorithm

The simulation uses the Gillespie algorithm (first-reaction method):

1. **Event Generation**: Scan all lattice sites and enumerate possible reactions based on current state and local environment
2. **Event Selection**: Randomly select an event weighted by reaction rates
3. **Execution**: Apply the selected reaction, update site states and neighbors
4. **Time Advancement**: Calculate time increment as `-ln(rand) / total_rate`
5. **Cleanup**: Remove unattached clusters using breadth-first search
6. **Repeat**: Continue for specified number of MC steps

## Code Structure

```
model/
├── mckaol.cpp      # Main program
├── lattice.cpp/hpp # Lattice data structure and initialization
├── ucell.cpp/hpp   # Unit cell configuration
├── sim.cpp/hpp     # Simulation parameters
├── rxnlist.cpp/hpp # Reaction database and rate calculation
├── evtlist.cpp/hpp # Event list management
├── envrn.cpp/hpp   # Environment checking for rate selection
├── actions.cpp/hpp # Reaction execution
├── output.cpp/hpp  # Output file generation
├── bfsearch.cpp/hpp# Graph connectivity algorithms
├── ran2.cpp/hpp    # Random number generation
├── futil.cpp/hpp   # File I/O utilities
├── myerr.cpp/hpp   # Error handling
├── common.hpp      # Constants and state encoding documentation
└── Makefile        # Build configuration
```

## Sample Input Files

Sample input files for testing are provided in the `input/` directory (if available) or can be generated based on the format specifications above.

## Performance Notes

- Memory usage scales with lattice size: O(aCells × bCells × positions_per_cell)
- Larger lattices provide better statistics but require more computation
- Event list regeneration is the main computational cost per step
- Use the `-O3 -ffast-math` optimization flags for production runs

## Visualization

Output MSI files can be visualized using:
- The companion **viewer** application in this repository
- Cerius2 molecular modeling software
- Any molecular visualization tool that supports MSI format

XYZ files are compatible with most molecular viewers including VMD, OVITO, and ASE.

## Citation

If you use this software in published research, please cite appropriately.

## License

Copyright information to be added.
