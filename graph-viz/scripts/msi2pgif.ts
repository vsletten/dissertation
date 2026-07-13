/**
 * One-way door out of CERIUS2: convert legacy .msi files to PGIF.
 *
 *   npm run msi2pgif -- input.msi [output.pgif.bin]
 *   npm run msi2pgif -- input.msi --json [output.pgif.json]
 *
 * Binary by default (typed-array columns, zero-parse load); --json for the
 * human-readable form. Uses the exact loader/encoder the viewer ships.
 */

import { readFileSync, writeFileSync } from 'node:fs';

import { loadMsi } from '../src/loaders/msi';
import { encodeBinary, encodeJson } from '../src/pgif/encode';

const args = process.argv.slice(2).filter(a => a !== '--');
const json = args.includes('--json');
const positional = args.filter(a => !a.startsWith('--'));

if (positional.length < 1) {
  console.error('usage: npm run msi2pgif -- <input.msi> [--json] [output]');
  process.exit(2);
}

const input = positional[0];
const output = positional[1] ?? input.replace(/\.msi$/i, '') + (json ? '.pgif.json' : '.pgif.bin');

const graph = loadMsi(readFileSync(input, 'utf-8'));
if (graph.count === 0) {
  console.error(`no atoms found in ${input} — not an MSI file?`);
  process.exit(1);
}

if (json) {
  writeFileSync(output, encodeJson(graph));
} else {
  writeFileSync(output, encodeBinary(graph));
}
console.log(
  `${input} → ${output}: ${graph.count.toLocaleString()} nodes, ` +
    `${graph.edgeCount.toLocaleString()} edges, types [${graph.typeDict.join(', ')}]`,
);
