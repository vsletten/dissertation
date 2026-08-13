/**
 * Trajectory playback: the petra event-log sidecar to a PGIF snapshot
 * (docs/PGIF.md, trajectory section). The log carries (site, old, new)
 * state deltas per KMC event, so the player can seek both directions —
 * forward applies `new`, backward restores `old`.
 *
 * The player mutates the loaded Graph's `state` column and `typeIndex`
 * in place (the columnar graph stays the single source of truth) and
 * reports which nodes it touched; the scene applies those to the GPU
 * instance buffers.
 */

import type { Graph } from './graph';

export interface TrajectoryHeader {
  petra_traj: 1;
  deck: string;
  seed: number;
  n_sites: number;
  /** State names, indexed by the dense state ids used in events. */
  states: string[];
  /** Display type per state (element name or 'vacant'). */
  state_types: string[];
  /** Reaction names, indexed by the event rows' reaction ids. */
  reactions: string[];
}

export interface TrajectoryEvent {
  step: number;
  time: number;
  rxn: number;
  /** [site, oldState, newState] — actual changes only. */
  changes: [number, number, number][];
}

export interface Trajectory {
  header: TrajectoryHeader;
  events: TrajectoryEvent[];
}

export function parseTrajectory(text: string): Trajectory {
  const lines = text.split('\n');
  let i = 0;
  while (i < lines.length && lines[i].trim() === '') i++;
  if (i >= lines.length) throw new Error('trajectory: empty file');
  const header = JSON.parse(lines[i]) as TrajectoryHeader;
  if (header.petra_traj !== 1) {
    throw new Error(`trajectory: unsupported version ${String((header as { petra_traj?: unknown }).petra_traj)}`);
  }
  if (!Array.isArray(header.states) || !Array.isArray(header.reactions)) {
    throw new Error('trajectory: header missing states/reactions');
  }
  const events: TrajectoryEvent[] = [];
  for (i += 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (line === '') continue;
    const row = JSON.parse(line) as [number, number, number, [number, number, number][]];
    events.push({ step: row[0], time: row[1], rxn: row[2], changes: row[3] });
  }
  return { header, events };
}

/**
 * Bind a trajectory to a loaded graph. Throws if they don't match
 * (site count, state dict). Call `prepareGraphForTrajectory` BEFORE
 * `scene.setGraph` so any display types that never occur at t=0 (e.g. a
 * fully occupied crystal that dissolves) still get typeDict entries.
 */
export function prepareGraphForTrajectory(graph: Graph, header: TrajectoryHeader): Int32Array {
  if (header.n_sites !== graph.count) {
    throw new Error(`trajectory: ${header.n_sites} sites vs graph's ${graph.count} nodes`);
  }
  const stateCol = graph.nodeColumns['state'];
  if (!stateCol || stateCol.type !== 'categorical' || !stateCol.dict) {
    throw new Error("trajectory: graph has no categorical 'state' column (not a petra snapshot?)");
  }
  if (stateCol.dict.length !== header.states.length) {
    throw new Error(
      `trajectory: state dict mismatch (${stateCol.dict.length} in snapshot, ${header.states.length} in log)`,
    );
  }
  // type per state id, extending the graph's typeDict for types unseen at t=0
  const typeOfState = new Int32Array(header.states.length);
  header.state_types.forEach((name, s) => {
    let t = graph.typeDict.indexOf(name);
    if (t < 0) {
      graph.typeDict.push(name);
      t = graph.typeDict.length - 1;
    }
    typeOfState[s] = t;
  });
  return typeOfState;
}

export class TrajectoryPlayer {
  /** Number of events currently applied (0 = snapshot state). */
  private cursor = 0;
  private warned = false;

  constructor(
    private readonly traj: Trajectory,
    private readonly graph: Graph,
    private readonly typeOfState: Int32Array,
  ) {}

  get length(): number {
    return this.traj.events.length;
  }

  get position(): number {
    return this.cursor;
  }

  /** Simulated time and step of the last applied event (0 at snapshot). */
  get time(): number {
    return this.cursor === 0 ? 0 : this.traj.events[this.cursor - 1].time;
  }

  get step(): number {
    return this.cursor === 0 ? 0 : this.traj.events[this.cursor - 1].step;
  }

  get lastReaction(): string | null {
    if (this.cursor === 0) return null;
    return this.traj.header.reactions[this.traj.events[this.cursor - 1].rxn] ?? null;
  }

  /**
   * Move to `target` events applied, mutating the graph's state column and
   * typeIndex; touched node indices are added to `touched`.
   */
  seek(target: number, touched: Set<number>): void {
    const t = Math.max(0, Math.min(this.length, Math.floor(target)));
    const states = this.graph.nodeColumns['state'].data as Int32Array;
    while (this.cursor < t) {
      for (const [site, old, next] of this.traj.events[this.cursor].changes) {
        if (states[site] !== old && !this.warned) {
          console.warn(`trajectory: site ${site} expected state ${old}, found ${states[site]}`);
          this.warned = true;
        }
        states[site] = next;
        this.graph.typeIndex[site] = this.typeOfState[next];
        touched.add(site);
      }
      this.cursor++;
    }
    while (this.cursor > t) {
      this.cursor--;
      const changes = this.traj.events[this.cursor].changes;
      for (let i = changes.length - 1; i >= 0; i--) {
        const [site, old, next] = changes[i];
        if (states[site] !== next && !this.warned) {
          console.warn(`trajectory: site ${site} expected state ${next}, found ${states[site]}`);
          this.warned = true;
        }
        states[site] = old;
        this.graph.typeIndex[site] = this.typeOfState[old];
        touched.add(site);
      }
    }
  }
}
