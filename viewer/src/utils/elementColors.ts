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
 * Get radius for a chemical element (in Angstroms)
 * These are approximations of van der Waals radii
 */
export const ELEMENT_RADII: Record<string, number> = {
  H: 0.32,
  C: 0.75,
  N: 0.71,
  O: 0.60,
  F: 0.50,
  P: 1.10,
  S: 1.02,
  Cl: 0.99,
  Si: 1.15,
  Al: 1.18,
  Na: 1.16,
  Mg: 1.10,
  K: 1.52,
  Ca: 1.26,
  Fe: 1.16,
  Mn: 1.17,
  Ti: 1.32,
  Zn: 1.09,
};

// Default radius for elements not in the map
export const DEFAULT_ELEMENT_RADIUS = 0.75;

/**
 * Get display radius for a chemical element
 */
export function getElementRadius(element: string): number {
  return ELEMENT_RADII[element] || DEFAULT_ELEMENT_RADIUS;
}