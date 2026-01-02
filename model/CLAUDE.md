# CLAUDE.md - Model Component

This file provides Claude Code guidance for the kinetic Monte Carlo simulation component.

## Quick Reference

### Build Commands
```bash
make           # Build mckaol executable (production: -O3 -ffast-math)
make clean     # Remove object files and executable
./mckaol       # Run simulation (requires data.* input files in current directory)
```

### Debug Build
Edit Makefile to comment out optimization flags and add `-ggdb`:
```makefile
# CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -ffast-math
CXXFLAGS = -std=c++11 -Wall -Wextra -ggdb
```

## Source File Map

| File | Purpose | Key Functions/Classes |
|------|---------|----------------------|
| mckaol.cpp | Main entry point | `main()` - initialization and simulation loop |
| lattice.cpp/hpp | Lattice structure | `Lattice`, `LatticeSite` classes |
| ucell.cpp/hpp | Unit cell template | `UnitCell`, `CellSite` classes |
| sim.cpp/hpp | Simulation parameters | `Simulation` class |
| rxnlist.cpp/hpp | Reaction database | `ReactionList`, `Reaction` classes |
| evtlist.cpp/hpp | Event list | `EventList`, `Event` classes |
| envrn.cpp/hpp | Environment checking | `Environment` class with Check100-500 methods |
| actions.cpp/hpp | Reaction execution | `Actions` class with DoEvent, DoReaction methods |
| output.cpp/hpp | File I/O | `output` static class |
| bfsearch.cpp/hpp | Graph algorithms | `BreadthFirstSearch` class |
| ran2.cpp/hpp | Random numbers | `initran2()`, `ran2()` |
| futil.cpp/hpp | File utilities | `Futil` static class |
| myerr.cpp/hpp | Error handling | `Myerr::die()` |
| common.hpp | Constants and enums | Site state encoding documentation |

## Core Data Structures

### LatticeSite (lattice.hpp)
```cpp
struct LatticeSite {
    int a, b;        // Unit cell coordinates
    int n;           // Position within unit cell (0 to Npos-1)
    int state;       // Site state (100-599 encoding)
    int pair;        // For 400/500 sites, paired bridging oxygen index
    int lostal;      // Al coordination tracking
    int nbr[6];      // Six neighbors (index or -1)
    Color color;     // BFS marking: UNREACHABLE/ENQUEUED/DISCOVERED
};
```

### Site State Encoding (common.hpp)
Three-digit integer: `[Type][Coordination]`
- **1xx (Al sites)**: 100=empty, 101-107=Al with 0-6 OH/H2O ligands
- **2xx (Si sites)**: 200=empty, 201-205=Si with 0-4 OH groups
- **3xx (Si-O-Si)**: 300=empty, 301-303=bridging oxygen states
- **4xx (Si-O-Al2)**: 400=empty, 401-410=complex bridging states
- **5xx (Al-OH-Al)**: 500=empty, 501-503=Al bridging states
- **9 (EDGE)**: Boundary site
- **99 (WRONG)**: Error state

### Constants (lattice.hpp)
```cpp
const int AL = 1;      // Type 1: Al sites
const int SI = 2;      // Type 2: Si sites
const int SIOSI = 3;   // Type 3: Si-O-Si bridging oxygen
const int SIOAL = 4;   // Type 4: Si-O-Al2 bridging oxygen
const int ALOHAL = 5;  // Type 5: Al-OH-Al bridging oxygen
const int EDGE = 9;    // Boundary marker
const int WRONG = 99;  // Error marker
```

## Simulation Algorithm

### Main Loop (mckaol.cpp:61-99)
```
For each MC step:
  1. CreateEventList() - Build linked list of all possible reactions
  2. DoEvent() - Select random event weighted by rate, execute it
  3. RemoveUnattachedClusters() - BFS cleanup
  4. Update time: dt = -ln(rand) / sum(rates)
  5. Write output at specified intervals
```

### Event Selection (actions.cpp:DoEvent)
Gillespie algorithm (first-reaction method):
1. Sum all event rates
2. Generate random number eps in [0,1]
3. Iterate events, accumulating normalized rates
4. Select event where cumulative sum >= eps

### Environment-Dependent Rates
Each reaction has multiple rate variants selected by local neighbor configuration:
- `Check100()`: 3 variants based on 400/500 neighbor states
- `Check200()`: 4 variants based on 300/400 neighbor states
- `Check300()`: 10 variants based on neighboring 300/400 states
- `Check400()`: Up to 40 variants (most complex)
- `Check500()`: 18 variants based on Al coordination

## Input File Formats

### data.sim
```
5000000     # Number of MC steps
1000        # Data write frequency (0=none)
1000000     # MSI snapshot frequency (0=none, -1=final only)
-2          # Random seed (negative = initialize)
1           # Draw bonds (0=atoms only, 1=full lattice)
```

### data.lattice
```
20 3 0      # aCells bCells SurfacePlane
            # SurfacePlane: 0=ac plane (b open), 1=bc plane (a open)
```

### data.rxn
```
8000.0      # Temperature (K)
-1.0        # Delta mu for Si
-1.0        # Delta mu for Al
# Hydrolysis reactions (pairs of forward/reverse):
301 302     # Reactant Product
15          # Number of rate variants
1 2.6 1 2.6 ...  # (k+, dE) pairs for each environment
# Adsorption reactions:
100 107     # Reactant Product
100.0 -14.5 # Pre-exponential, energy offset
# Desorption reactions:
107         # Reactant (product is 100)
4           # Number of variants
1998.0 0.0 1998.0 12.0 ...  # (a, dE) pairs
```

### data.cell
```
5.140 8.930 7.370    # Cell parameters a, b, c (angstroms)
-0.03141 -0.25038 0.0  # Cell angles alpha, beta, gamma (radians)
26                   # Number of positions
# For each position:
0 100 1.76027 4.49982 2.84548   # id type x y z
  8   0  0  0    # neighbor: nbr_id da db dc
  9   0  0  0
  ...
```

## Reaction Types

### Hydrolysis (Reactions 0-15)
Forward/reverse pairs for bond breaking/forming:
- R0/R1: Si-O-Si hydrolysis (301 ↔ 302)
- R2/R3: Si-O-Al2 first step (401 ↔ 402)
- R4/R5: Si-O-Al2 alternate path (401 ↔ 410)
- R6-R15: Complex multi-step Al coordination changes

### Adsorption (Reactions 16-19)
- Al adsorption: 100 → 107
- Si adsorption: 200 → 205
- Rate = k * exp(dE/RT) * exp(dm/RT)

### Desorption (Reactions 20-23)
- Al desorption: 107 → 100
- Si desorption: 205 → 200
- Rate = k * exp(-dE/RT)

## Output Files

| File | Content |
|------|---------|
| results.dat | Time series: time, Si counts, Al counts by coordination |
| start.msi | Initial structure (Cerius2 MSI format) |
| end.msi | Final structure |
| start.xyz | Initial coordinates (XYZ format) |
| surfSi.out | Surface Si positions (x,y pairs) |
| surfAl.out | Surface Al positions (x,y pairs) |
| stepN.msi | Intermediate snapshots at step N |

## Key Algorithms

### Lattice Initialization
1. `CreateLattice()`: Allocate sites, copy unit cell template with periodic boundaries
2. `FindPairs()`: Link bridging oxygen pairs (400/500 sites)
3. `PopulateSolid()`: Fill sites based on chemical potential
4. `TerminateSurface()`: Create surface hydroxyl groups
5. `TerminateLattice()`: Mark boundary sites as EDGE

### Breadth-First Search (bfsearch.cpp)
Used to identify and remove unattached atomic clusters after reactions. Marks all connected atoms as DISCOVERED, leaving floating atoms as UNREACHABLE for removal.

### Boltzmann Rate Calculation (rxnlist.cpp)
```cpp
k_forward = k+
k_reverse = k+ * exp(dE / RT)  // RT = 1.987e-3 * T
k_ads = a * exp(dE/RT) * exp(dm/RT)
k_des = a * exp(-dE/RT)
```

## Common Modifications

### Adding a New Reaction
1. Add entry to data.rxn with reactant, product, rate data
2. Increment reaction count in rxnlist.cpp
3. Add DoReactionN() method in actions.cpp
4. Add environment check in envrn.cpp if needed
5. Update event creation in evtlist.cpp

### Changing Lattice Size
Edit data.lattice: `aCells bCells SurfacePlane`
Memory scales as O(aCells * bCells * Npos)

### Adjusting Output Frequency
Edit data.sim:
- Line 2: data write frequency (wsteps)
- Line 3: MSI snapshot frequency (msteps)

## Code Style

- C++11 standard
- Classes use PascalCase, methods use PascalCase
- Member variables without prefix
- Static classes for utilities (output, Futil, Myerr)
- Constants use UPPER_CASE
- Error handling via Myerr::die()
