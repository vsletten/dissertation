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
    this.rebuild();
    this.frameCamera();
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

    // visible edges: both endpoints visible
    let evis = 0;
    for (let e = 0; e < g.edgeCount; e++) {
      if (!this.hiddenTypes.has(g.typeIndex[g.edgeSrc[e]]) && !this.hiddenTypes.has(g.typeIndex[g.edgeDst[e]])) evis++;
    }
    this.visibleEdges = evis;
    if (evis > 0) {
      const bondRadius = this.radiusScale * 0.18;
      const edgeGeo = new CylinderGeometry(bondRadius, bondRadius, 1, bondDetail, 1, true);
      const edgeMat = new MeshStandardMaterial({ color: '#7d8ba1', roughness: 0.6, metalness: 0.05 });
      const edges = new InstancedMesh(edgeGeo, edgeMat, evis);
      const up = new Vector3(0, 1, 0);
      const q = new Quaternion();
      const a = new Vector3(), b = new Vector3(), mid = new Vector3(), axis = new Vector3();
      let eslot = 0;
      for (let e = 0; e < g.edgeCount; e++) {
        const si = g.edgeSrc[e], di = g.edgeDst[e];
        if (this.hiddenTypes.has(g.typeIndex[si]) || this.hiddenTypes.has(g.typeIndex[di])) continue;
        a.set(g.positions[si * 3], g.positions[si * 3 + 1], g.positions[si * 3 + 2]);
        b.set(g.positions[di * 3], g.positions[di * 3 + 1], g.positions[di * 3 + 2]);
        mid.addVectors(a, b).multiplyScalar(0.5);
        axis.subVectors(b, a);
        const len = Math.max(axis.length(), 1e-6);
        q.setFromUnitVectors(up, axis.normalize());
        m.makeRotationFromQuaternion(q);
        m.scale(new Vector3(1, len, 1));
        m.setPosition(mid);
        edges.setMatrixAt(eslot, m);
        eslot++;
      }
      edges.instanceMatrix.needsUpdate = true;
      this.edgeMesh = edges;
      this.group.add(edges);
    }
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
