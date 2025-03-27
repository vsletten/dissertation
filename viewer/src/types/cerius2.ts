// Simple types for atoms and bonds
export type AttributeType = 'B' | 'C' | 'D' | 'F' | 'I' | 'O' | 'S' | 'T';

export interface Attribute {
  type: AttributeType;
  name: string;
  value: any;
}

export interface Atom {
  id: number;
  element: string;
  position: [number, number, number];
  label: string;
  charge: number;
}

export interface Bond {
  from: number;
  to: number;
}

// Dummy type for backward compatibility, not used anymore
export interface Cerius2Model {
  id: number;
  type: string;
}

export interface Cerius2Object {
  id: number;
  type: string;
}