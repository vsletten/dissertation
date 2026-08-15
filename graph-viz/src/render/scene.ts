/**
 * The render core: one InstancedMesh for nodes, one for edges, adaptive
 * geometry detail by graph size. This is the architecture the bench proved —
 * a fixed handful of draw calls no matter the node count.
 *
 * Type visibility filtering rebuilds the instance buffers from the columnar
 * source arrays (O(N), a few ms at 1M) rather than tracking per-instance
 * state — the columnar Graph is the single source of truth.
 */

import {
  AmbientLight,
  Color,
  CylinderGeometry,
  DirectionalLight,
  Group,
  HemisphereLight,
  InstancedMesh,
  Matrix4,
  MeshStandardMaterial,
  PerspectiveCamera,
  Quaternion,
  Scene,
  SphereGeometry,
  Vector3,
  WebGLRenderer,
} from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

import { computeBounds, medianEdgeLength, type Bounds, type Graph } from '../graph';
import { stylesFor, type TypeStyle } from './style';

export interface SceneStats {
  fps: number;
  frameMs: number;
  drawCalls: number;
  triangles: number;
  visibleNodes: number;
  visibleEdges: number;
}

export class GraphScene {
  readonly renderer: WebGLRenderer;
  readonly scene: Scene;
  readonly camera: PerspectiveCamera;
  readonly controls: OrbitControls;

  private group = new Group();
  private nodeMesh: InstancedMesh | null = null;
  private edgeMesh: InstancedMesh | null = null;
  private graph: Graph | null = null;
  private styles: TypeStyle[] = [];
  private hiddenTypes = new Set<number>();
  /** instance slot → node index (visibility filtering compacts the buffer) */
  private slotToNode = new Uint32Array(0);
  /** edge index → edge instance slot (-1 when hidden by type filtering) */
  private edgeSlotOf = new Int32Array(0);
  /** CSR incidence: edges touching each node (for in-place updates) */
  private incOff = new Uint32Array(0);
  private incEdge = new Uint32Array(0);
  /** `seam` edge column (periodic wrap bonds — never drawn), if present */
  private seam: Uint8Array | null = null;
  private visibleNodes = 0;
  private visibleEdges = 0;
  private bounds: Bounds | null = null;
  private radiusScale = 1;

  private frameMsWindow: number[] = [];
  private lastTs = 0;
  private stats: SceneStats = { fps: 0, frameMs: 0, drawCalls: 0, triangles: 0, visibleNodes: 0, visibleEdges: 0 };
  onFrame?: (stats: SceneStats) => void;

  constructor(container: HTMLElement) {
    this.renderer = new WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.setSize(container.clientWidth, container.clientHeight);
    this.renderer.setClearColor(new Color('#0b1018'), 1);
    container.appendChild(this.renderer.domElement);

    this.scene = new Scene();
    this.camera = new PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 50_000);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;

    this.scene.add(new HemisphereLight(0xcfe0ff, 0x1a2030, 0.9));
    const key = new DirectionalLight(0xffffff, 1.6);
    key.position.set(1, 1.4, 0.8);
    this.scene.add(key);
    const fill = new DirectionalLight(0x88aaff, 0.35);
    fill.position.set(-1, -0.4, -0.7);
    this.scene.add(fill);
    this.scene.add(new AmbientLight(0xffffff, 0.25));
    this.scene.add(this.group);

    new ResizeObserver(() => {
      const w = container.clientWidth, h = container.clientHeight;
      if (!w || !h) return;
      this.camera.aspect = w / h;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(w, h);
    }).observe(container);

    this.renderer.setAnimationLoop(ts => this.frame(ts));
  }

  setGraph(graph: Graph): void {
    this.graph = graph;
    this.styles = stylesFor(graph.typeDict);
    this.hiddenTypes.clear();
    if (!graph.positions) throw new Error('GraphScene.setGraph requires positions (run layout first)');
    this.bounds = computeBounds(graph.positions, graph.count);
    // radius that reads well at this graph's density
    this.radiusScale = Math.max(0.05, Math.min(2.5, medianEdgeLength(graph) * 0.16));
    const seamCol = graph.edgeColumns['seam'];
    this.seam = seamCol && seamCol.type === 'bool' ? (seamCol.data as Uint8Array) : null;
    this.buildIncidence(graph);
    this.rebuild();
    this.frameCamera();
  }

  /** CSR edge lists per node — the in-place update path's lookup. */
  private buildIncidence(g: Graph): void {
    const deg = new Uint32Array(g.count);
    for (let e = 0; e < g.edgeCount; e++) {
      deg[g.edgeSrc[e]]++;
      deg[g.edgeDst[e]]++;
    }
    this.incOff = new Uint32Array(g.count + 1);
    for (let i = 0; i < g.count; i++) this.incOff[i + 1] = this.incOff[i] + deg[i];
    this.incEdge = new Uint32Array(this.incOff[g.count]);
    const cursor = Uint32Array.from(this.incOff.subarray(0, g.count));
    for (let e = 0; e < g.edgeCount; e++) {
      this.incEdge[cursor[g.edgeSrc[e]]++] = e;
      this.incEdge[cursor[g.edgeDst[e]]++] = e;
    }
  }

  getStyles(): TypeStyle[] {
    return this.styles;
  }

  isTypeHidden(t: number): boolean {
    return this.hiddenTypes.has(t);
  }

  toggleType(t: number): void {
    if (this.hiddenTypes.has(t)) this.hiddenTypes.delete(t);
    else this.hiddenTypes.add(t);
    this.rebuild();
  }

  /** Node index for an instance slot (picking). */
  nodeForSlot(slot: number): number {
    return this.slotToNode[slot];
  }

  getGraph(): Graph | null {
    return this.graph;
  }

  getBounds(): Bounds | null {
    return this.bounds;
  }

  getRadiusScale(): number {
    return this.radiusScale;
  }

  frameCamera(): void {
    if (!this.bounds) return;
    const [cx, cy, cz] = this.bounds.center;
    const d = this.bounds.extent * 1.1 + 4;
    this.controls.target.set(cx, cy, cz);
    this.camera.position.set(cx + d * 0.85, cy + d * 0.65, cz + d * 0.55);
    this.camera.near = Math.max(0.01, d / 5000);
    this.camera.far = Math.max(10_000, d * 40);
    this.camera.updateProjectionMatrix();
    this.controls.update();
  }

  /** Rebuild instance buffers from the columnar graph, honoring hidden types. */
  private rebuild(): void {
    const g = this.graph;
    if (!g || !g.positions) return;

    this.disposeMeshes();

    // node detail budget: heavy spheres for small graphs, light for huge ones
    const detail = g.count <= 20_000 ? [24, 16] : g.count <= 200_000 ? [14, 10] : [8, 6];
    const bondDetail = g.count <= 200_000 ? 8 : 5;

    // visible nodes
    let visible = 0;
    for (let i = 0; i < g.count; i++) if (!this.hiddenTypes.has(g.typeIndex[i])) visible++;
    this.visibleNodes = visible;
    this.slotToNode = new Uint32Array(visible);

    const nodeGeo = new SphereGeometry(1, detail[0], detail[1]);
    const nodeMat = new MeshStandardMaterial({ roughness: 0.42, metalness: 0.08 });
    const nodes = new InstancedMesh(nodeGeo, nodeMat, visible);
    const m = new Matrix4();
    const color = new Color();
    const colorCache = this.styles.map(s => new Color(s.color));
    let slot = 0;
    for (let i = 0; i < g.count; i++) {
      const t = g.typeIndex[i];
      if (this.hiddenTypes.has(t)) continue;
      const r = this.styles[t].radius * this.radiusScale;
      m.makeScale(r, r, r);
      m.setPosition(g.positions[i * 3], g.positions[i * 3 + 1], g.positions[i * 3 + 2]);
      nodes.setMatrixAt(slot, m);
      nodes.setColorAt(slot, colorCache[t] ?? color.set('#ffffff'));
      this.slotToNode[slot] = i;
      slot++;
    }
    nodes.instanceMatrix.needsUpdate = true;
    if (nodes.instanceColor) nodes.instanceColor.needsUpdate = true;
    this.nodeMesh = nodes;
    this.group.add(nodes);

    // visible edges: both endpoints' TYPES visible. Slots are allocated for
    // all of these; seam edges and edges with a radius-0 (vacant) endpoint
    // keep their slot but get a zero matrix, so trajectory playback can
    // flip them in place without recompacting.
    let evis = 0;
    this.edgeSlotOf = new Int32Array(g.edgeCount).fill(-1);
    for (let e = 0; e < g.edgeCount; e++) {
      if (!this.hiddenTypes.has(g.typeIndex[g.edgeSrc[e]]) && !this.hiddenTypes.has(g.typeIndex[g.edgeDst[e]])) {
        this.edgeSlotOf[e] = evis++;
      }
    }
    if (evis > 0) {
      const bondRadius = this.radiusScale * 0.18;
      const edgeGeo = new CylinderGeometry(bondRadius, bondRadius, 1, bondDetail, 1, true);
      const edgeMat = new MeshStandardMaterial({ color: '#7d8ba1', roughness: 0.6, metalness: 0.05 });
      const edges = new InstancedMesh(edgeGeo, edgeMat, evis);
      this.edgeMesh = edges;
      let drawn = 0;
      for (let e = 0; e < g.edgeCount; e++) {
        if (this.edgeSlotOf[e] >= 0 && this.applyEdgeMatrix(e)) drawn++;
      }
      edges.instanceMatrix.needsUpdate = true;
      this.visibleEdges = drawn;
      this.group.add(edges);
    } else {
      this.visibleEdges = 0;
    }
  }

  private static readonly edgeUp = new Vector3(0, 1, 0);
  private edgeTmp = {
    m: new Matrix4(),
    q: new Quaternion(),
    a: new Vector3(),
    b: new Vector3(),
    mid: new Vector3(),
    axis: new Vector3(),
    scale: new Vector3(),
  };

  /**
   * Write edge `e`'s instance matrix into its slot: oriented cylinder, or a
   * zero matrix when the edge is a seam or an endpoint is invisible
   * (radius 0). Returns whether the edge is drawn.
   */
  private applyEdgeMatrix(e: number): boolean {
    const g = this.graph;
    const slot = this.edgeSlotOf[e];
    if (!g || !g.positions || !this.edgeMesh || slot < 0) return false;
    const { m, q, a, b, mid, axis, scale } = this.edgeTmp;
    const si = g.edgeSrc[e], di = g.edgeDst[e];
    const hidden =
      (this.seam && this.seam[e]) ||
      this.styles[g.typeIndex[si]].radius <= 0 ||
      this.styles[g.typeIndex[di]].radius <= 0;
    if (hidden) {
      m.makeScale(0, 0, 0);
      this.edgeMesh.setMatrixAt(slot, m);
      return false;
    }
    a.set(g.positions[si * 3], g.positions[si * 3 + 1], g.positions[si * 3 + 2]);
    b.set(g.positions[di * 3], g.positions[di * 3 + 1], g.positions[di * 3 + 2]);
    mid.addVectors(a, b).multiplyScalar(0.5);
    axis.subVectors(b, a);
    const len = Math.max(axis.length(), 1e-6);
    q.setFromUnitVectors(GraphScene.edgeUp, axis.normalize());
    m.makeRotationFromQuaternion(q);
    m.scale(scale.set(1, len, 1));
    m.setPosition(mid);
    this.edgeMesh.setMatrixAt(slot, m);
    return true;
  }

  /**
   * In-place refresh of specific nodes after their type/state changed
   * (trajectory playback). Fast path requires no type filtering (slots are
   * then identity); with hidden types it falls back to a full rebuild.
   */
  updateNodeStates(touched: Iterable<number>): void {
    const g = this.graph;
    if (!g || !g.positions || !this.nodeMesh) return;
    if (this.hiddenTypes.size > 0) {
      this.rebuild();
      return;
    }
    const m = this.edgeTmp.m;
    const color = new Color();
    for (const i of touched) {
      const t = g.typeIndex[i];
      const style = this.styles[t];
      const r = (style?.radius ?? 0.8) * this.radiusScale;
      m.makeScale(r, r, r);
      m.setPosition(g.positions[i * 3], g.positions[i * 3 + 1], g.positions[i * 3 + 2]);
      this.nodeMesh.setMatrixAt(i, m);
      this.nodeMesh.setColorAt(i, color.set(style?.color ?? '#ffffff'));
      for (let k = this.incOff[i]; k < this.incOff[i + 1]; k++) {
        this.applyEdgeMatrix(this.incEdge[k]);
      }
    }
    this.nodeMesh.instanceMatrix.needsUpdate = true;
    if (this.nodeMesh.instanceColor) this.nodeMesh.instanceColor.needsUpdate = true;
    if (this.edgeMesh) this.edgeMesh.instanceMatrix.needsUpdate = true;
  }

  private disposeMeshes(): void {
    for (const mesh of [this.nodeMesh, this.edgeMesh]) {
      if (!mesh) continue;
      this.group.remove(mesh);
      mesh.geometry.dispose();
      (mesh.material as MeshStandardMaterial).dispose();
      mesh.dispose();
    }
    this.nodeMesh = null;
    this.edgeMesh = null;
  }

  private frame(ts: number): void {
    if (this.lastTs) {
      const dt = ts - this.lastTs;
      this.frameMsWindow.push(dt);
      if (this.frameMsWindow.length > 60) this.frameMsWindow.shift();
      const avg = this.frameMsWindow.reduce((s, v) => s + v, 0) / this.frameMsWindow.length;
      this.stats = {
        fps: avg > 0 ? 1000 / avg : 0,
        frameMs: avg,
        drawCalls: this.renderer.info.render.calls,
        triangles: this.renderer.info.render.triangles,
        visibleNodes: this.visibleNodes,
        visibleEdges: this.visibleEdges,
      };
      this.onFrame?.(this.stats);
    }
    this.lastTs = ts;
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  getStats(): SceneStats {
    return this.stats;
  }
}
