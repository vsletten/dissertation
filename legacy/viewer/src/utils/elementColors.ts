/**
 * CPK coloring scheme for chemical elements
 * Common atoms in molecular visualizations
 */
export const ELEMENT_COLORS: Record<string, string> = {
  H: '#FFFFFF', // White
  C: '#909090', // Light grey
  N: '#3050F8', // Blue
  O: '#FF0D0D', // Red
  F: '#90E050', // Light green
  P: '#FF8000', // Orange
  S: '#FFFF30', // Yellow
  Cl: '#1FF01F', // Green
  Si: '#F0C8A0', // Tan/beige
  Al: '#BFA6A6', // Light purple/grey
  Na: '#AB5CF2', // Purple
  Mg: '#8AFF00', // Bright green
  K: '#8F40D4',  // Purple
  Ca: '#3DFF00', // Green
  Fe: '#E06633', // Orange/brown
  // Additional elements common in kaolinite and clay minerals
  Mn: '#9C7AC7', // Purple
  Ti: '#BFC2C7', // Grey
  Zn: '#7D80B0', // Blue-grey
};

// Default color for elements not in the map
export const DEFAULT_ELEMENT_COLOR = '#A0A0A0';

/**
 * Get color for a chemical element
 */
export function getElementColor(element: string): string {
  return ELEMENT_COLORS[element] || DEFAULT_ELEMENT_COLOR;
}

/**
 * Covalent radii for chemical elements (in Angstroms)
 * Used for ball-and-stick visualization to represent bonding distances
 * Source: Cordero et al. (2008) "Covalent radii revisited"
 */
export const COVALENT_RADII: Record<string, number> = {
  H: 0.31,
  C: 0.76,
  N: 0.71,
  O: 0.66,
  F: 0.57,
  P: 1.07,
  S: 1.05,
  Cl: 1.02,
  Si: 1.11,
  Al: 1.21,
  Na: 1.66,
  Mg: 1.41,
  K: 2.03,
  Ca: 1.76,
  Fe: 1.32,
  Mn: 1.39,
  Ti: 1.60,
  Zn: 1.22,
};

/**
 * Van der Waals radii for chemical elements (in Angstroms)
 * Used for space-filling visualization to represent atomic volume
 * Source: Bondi (1964), Mantina et al. (2009)
 */
export const VDW_RADII: Record<string, number> = {
  H: 1.20,
  C: 1.70,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  P: 1.80,
  S: 1.80,
  Cl: 1.75,
  Si: 2.10,
  Al: 1.84,
  Na: 2.27,
  Mg: 1.73,
  K: 2.75,
  Ca: 2.31,
  Fe: 2.04,
  Mn: 2.05,
  Ti: 2.11,
  Zn: 2.01,
};

// Default radii for elements not in the map
export const DEFAULT_COVALENT_RADIUS = 0.75;
export const DEFAULT_VDW_RADIUS = 1.50;

/**
 * Get covalent radius for a chemical element (for ball-and-stick mode)
 */
export function getCovalentRadius(element: string): number {
  return COVALENT_RADII[element] || DEFAULT_COVALENT_RADIUS;
}

/**
 * Get van der Waals radius for a chemical element (for space-filling mode)
 */
export function getVdwRadius(element: string): number {
  return VDW_RADII[element] || DEFAULT_VDW_RADIUS;
}

/**
 * Get display radius for a chemical element
 * @deprecated Use getCovalentRadius or getVdwRadius instead
 */
export function getElementRadius(element: string): number {
  return getCovalentRadius(element);
}