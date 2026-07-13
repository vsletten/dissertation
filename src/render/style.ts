/**
 * Type → color/size styling. Known chemistry gets the classic palette; any
 * other categorical type gets a stable color hashed from its name, so an
 * arbitrary property graph is legible with zero configuration.
 */

const ELEMENT_STYLES: Record<string, { color: string; radius: number }> = {
  H: { color: '#e8e8e8', radius: 0.55 },
  C: { color: '#7a7a7a', radius: 0.75 },
  N: { color: '#4d7dff', radius: 0.72 },
  O: { color: '#ff5c47', radius: 0.7 },
  OH: { color: '#7fd98b', radius: 0.7 },
  F: { color: '#8be08b', radius: 0.65 },
  Al: { color: '#e3c236', radius: 0.95 },
  Si: { color: '#57aef2', radius: 0.9 },
  P: { color: '#ffa04d', radius: 0.85 },
  S: { color: '#ffe34d', radius: 0.85 },
  Fe: { color: '#c96a3d', radius: 0.95 },
  Mg: { color: '#77e0c8', radius: 0.9 },
  Na: { color: '#b56ee0', radius: 0.95 },
  K: { color: '#9a5bd1', radius: 1.0 },
  Ca: { color: '#66d17a', radius: 0.95 },
};

/** Golden-angle hue walk seeded by the name — distinct, stable, theme-friendly. */
function hashedStyle(name: string): { color: string; radius: number } {
  let h = 2166136261;
  for (let i = 0; i < name.length; i++) {
    h ^= name.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  const hue = ((h >>> 0) * 137.508) % 360;
  const sat = 62 + ((h >>> 8) % 20);
  const light = 56 + ((h >>> 16) % 14);
  return { color: `hsl(${hue.toFixed(1)}, ${sat}%, ${light}%)`, radius: 0.8 };
}

export interface TypeStyle {
  name: string;
  color: string;
  radius: number;
}

export function stylesFor(typeDict: string[]): TypeStyle[] {
  return typeDict.map(name => {
    const s = ELEMENT_STYLES[name] ?? hashedStyle(name);
    return { name, color: s.color, radius: s.radius };
  });
}
