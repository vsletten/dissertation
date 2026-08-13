/**
 * Trajectory sidecar: parse, bind to a snapshot graph, seek both
 * directions (docs/PGIF.md trajectory section).
 */
import { describe, expect, it } from 'vitest';

import { loadPgif } from '../src/pgif/decode';
import { parseTrajectory, prepareGraphForTrajectory, TrajectoryPlayer } from '../src/traj';

/** A 3-site petra-style snapshot: all occupied 'M', one bond chain. */
const SNAPSHOT = JSON.stringify({
  pgif: 1,
  meta: { producer: 'petra', kind: 'kmc-lattice' },
  nodes: {
    count: 3,
    columns: {
      x: { type: 'f32', data: [0, 2, 4] },
      y: { type: 'f32', data: [0, 0, 0] },
      z: { type: 'f32', data: [0, 0, 0] },
      type: { type: 'categorical', dict: ['M'], data: [0, 0, 0] },
      state: { type: 'categorical', dict: ['X.occ', 'X.gone'], data: [0, 0, 0] },
    },
  },
  edges: { count: 2, src: [0, 1], dst: [1, 2] },
});

/** Two events: site 2 dissolves, then site 1. */
const EVENTS = [
  JSON.stringify({
    petra_traj: 1,
    deck: 'mini',
    seed: 1,
    n_sites: 3,
    states: ['X.occ', 'X.gone'],
    state_types: ['M', 'vacant'],
    reactions: ['leave'],
  }),
  '[1,1.5e-1,0,[[2,0,1]]]',
  '[2,4.0e-1,0,[[1,0,1]]]',
].join('\n');

describe('trajectory', () => {
  it('parses header and compact event rows', () => {
    const traj = parseTrajectory(EVENTS);
    expect(traj.header.deck).toBe('mini');
    expect(traj.events).toHaveLength(2);
    expect(traj.events[0]).toEqual({ step: 1, time: 0.15, rxn: 0, changes: [[2, 0, 1]] });
  });

  it('extends typeDict with types unseen at t=0', () => {
    const g = loadPgif(SNAPSHOT);
    expect(g.typeDict).toEqual(['M']);
    const typeOfState = prepareGraphForTrajectory(g, parseTrajectory(EVENTS).header);
    expect(g.typeDict).toEqual(['M', 'vacant']);
    expect(Array.from(typeOfState)).toEqual([0, 1]);
  });

  it('seeks forward and backward, mutating state column and typeIndex', () => {
    const g = loadPgif(SNAPSHOT);
    const traj = parseTrajectory(EVENTS);
    const typeOfState = prepareGraphForTrajectory(g, traj.header);
    const player = new TrajectoryPlayer(traj, g, typeOfState);
    const states = g.nodeColumns['state'].data as Int32Array;
    const touched = new Set<number>();

    player.seek(2, touched);
    expect(Array.from(states)).toEqual([0, 1, 1]);
    expect(Array.from(g.typeIndex)).toEqual([0, 1, 1]);
    expect([...touched].sort()).toEqual([1, 2]);
    expect(player.time).toBeCloseTo(0.4);
    expect(player.lastReaction).toBe('leave');

    touched.clear();
    player.seek(0, touched);
    expect(Array.from(states)).toEqual([0, 0, 0]);
    expect(Array.from(g.typeIndex)).toEqual([0, 0, 0]);
    expect(player.position).toBe(0);
    expect(player.time).toBe(0);

    // clamped seeks
    player.seek(99, touched);
    expect(player.position).toBe(2);
    player.seek(-5, touched);
    expect(player.position).toBe(0);
  });

  it('rejects mismatched site counts and missing state column', () => {
    const g = loadPgif(SNAPSHOT);
    const badHeader = { ...parseTrajectory(EVENTS).header, n_sites: 7 };
    expect(() => prepareGraphForTrajectory(g, badHeader)).toThrow(/7 sites/);
  });
});
