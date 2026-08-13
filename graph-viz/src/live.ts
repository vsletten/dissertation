/**
 * Live simulation: petra compiled to WebAssembly, running in this page.
 *
 * The sim produces exactly the trajectory-sidecar rows of docs/PGIF.md §6a,
 * generated on demand into a *growing* Trajectory that the ordinary
 * TrajectoryPlayer consumes — so a live run gets pause, backward scrubbing
 * through its full history, and rate control for free, and "watching live"
 * is just keeping the playhead at the growing edge.
 *
 * The wasm package is served from public/petra-wasm/ (rebuild with
 * `npm run build:wasm`); determinism carries over — a (deck, seed) pair
 * reproduces the native CLI trajectory exactly (IEEE f64 + PCG64).
 */

import type { Trajectory, TrajectoryEvent, TrajectoryHeader } from './traj';

interface WasmSimInstance {
  snapshot_json(): string;
  header_json(): string;
  step_batch(n: number): string;
  stop_reason(): string | undefined;
  step_count(): number;
  time(): number;
  n_sites(): number;
}

interface WasmModule {
  default(input?: string | URL): Promise<unknown>;
  WasmSim: new (deckToml: string, seed: number) => WasmSimInstance;
}

let wasmModule: WasmModule | null = null;

async function loadWasm(): Promise<WasmModule> {
  if (wasmModule) return wasmModule;
  // Served statically from public/ — no bundler involvement, so the
  // viewer's build stays independent of the Rust toolchain. The import is
  // constructed at runtime (Function) so neither TS nor Vite's import
  // analysis touches it — Vite otherwise rewrites even @vite-ignore'd
  // dynamic imports of public assets and 500s on them.
  const dynImport = new Function('u', 'return import(u)') as (u: string) => Promise<unknown>;
  const mod = (await dynImport(new URL('/petra-wasm/petra_wasm.js', location.origin).href)) as WasmModule;
  await mod.default('/petra-wasm/petra_wasm_bg.wasm');
  wasmModule = mod;
  return mod;
}

/** Compact §6a row as emitted by the sim. */
type Row = [number, number, number, [number, number, number][]];

export class LiveSim {
  private constructor(
    private readonly sim: WasmSimInstance,
    readonly header: TrajectoryHeader,
    readonly traj: Trajectory,
  ) {}

  /** Events generated per ensure() slice — a frame's worth at top speed. */
  static readonly CHUNK = 2048;

  static async create(deckToml: string, seed: number): Promise<{ live: LiveSim; snapshotText: string }> {
    const mod = await loadWasm();
    const sim = new mod.WasmSim(deckToml, seed);
    const snapshotText = sim.snapshot_json();
    const header = JSON.parse(sim.header_json()) as TrajectoryHeader;
    const traj: Trajectory = { header, events: [] };
    return { live: new LiveSim(sim, header, traj), snapshotText };
  }

  /** null while the sim can still advance, else the stop reason. */
  get stopped(): string | null {
    return this.sim.stop_reason() ?? null;
  }

  /**
   * Generate events until the trajectory holds at least `target` (or the
   * sim stops). Returns the new event count.
   */
  ensure(target: number): number {
    const events = this.traj.events;
    while (events.length < target && !this.stopped) {
      const want = Math.min(LiveSim.CHUNK, target - events.length);
      const rows = JSON.parse(this.sim.step_batch(want)) as Row[];
      for (const [step, time, rxn, changes] of rows) {
        events.push({ step, time, rxn, changes } satisfies TrajectoryEvent);
      }
      if (rows.length < want) break; // stopped mid-batch
    }
    return events.length;
  }
}
