// geo.js
import * as THREE from 'three';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

export const LINE_WIDTH = .1;

export function createLabelRenderer(container){
  const r = new CSS2DRenderer();
  r.setSize(container.clientWidth, container.clientHeight);
  r.domElement.style.position = 'absolute';
  r.domElement.style.top = '0';
  r.domElement.style.left = '0';
  r.domElement.style.pointerEvents = 'none';
  container.appendChild(r.domElement);
  function resize(){
    const rect = container.getBoundingClientRect();
    r.setSize(rect.width, rect.height);
  }
  return { renderer: r, resize };
}

export function removeCubeAxes(scene){
  const oldLines = scene.getObjectByName('axes-lines');
  const oldLabels = scene.getObjectByName('axes-labels');
  if (oldLines) oldLines.removeFromParent();
  if (oldLabels) {
    // Remove CSS2D label DOM elements before removing from scene
    oldLabels.traverse((obj) => {
      if (obj instanceof CSS2DObject && obj.element && obj.element.parentNode) {
        obj.element.parentNode.removeChild(obj.element);
      }
    });
    oldLabels.removeFromParent();
  }
}

export function addCubeEdgeAxes(scene, mapper, cfg){
  // cfg = {
  //   x: {domain:[min,max], ticks:number[], label:(v)=>string}
  //   y: {domain:[min,max], ticks:number[], label:(v)=>string}
  //   z: {domain:[min,max], ticks:number[], label:(v)=>string}
  //   tickLen?: number, color?: number
  // }
  const color = cfg.color ?? 0x8b91a6;
  const tickColor = cfg.tickColor ?? color;
  const showEdges = cfg.showEdges ?? true;
  const tickOutwardOnly = cfg.tickOutwardOnly ?? false;
  const lineWidth = cfg.lineWidth ?? LINE_WIDTH * 0.25;
  const tickRadius = 1.0 * lineWidth;
  const tickLen = cfg.tickLen ?? 0.5;
  const labelOffset = cfg.labelOffset ?? 2.6;
  const labelColor = cfg.labelColor ?? '#dfe3ee';

  const linesGroup = new THREE.Group(); linesGroup.name = 'axes-lines';
  const labelsGroup = new THREE.Group(); labelsGroup.name = 'axes-labels';
  scene.add(linesGroup); scene.add(labelsGroup);

  const lineMat = new THREE.LineBasicMaterial({ color, transparent:true, opacity:0.9, linewidth: lineWidth });
  const tickMeshMat = new THREE.MeshBasicMaterial({ color: tickColor, transparent:true, opacity:0.9, depthWrite:false });

  const dom = {
    x0: cfg.x.domain[0], x1: cfg.x.domain[1],
    y0: cfg.y.domain[0], y1: cfg.y.domain[1],
    z0: cfg.z.domain[0], z1: cfg.z.domain[1],
  };

  // 12 cube edges as parametric lines: varying one axis, fixing the other two at min/max
  const edges = [
    // X-parallel edges (y,z fixed)
    { vary:'x', fix:{y:'y0', z:'z0'} }, { vary:'x', fix:{y:'y0', z:'z1'} },
    { vary:'x', fix:{y:'y1', z:'z0'} }, { vary:'x', fix:{y:'y1', z:'z1'} },
    // Y-parallel edges
    { vary:'y', fix:{x:'x0', z:'z0'} }, { vary:'y', fix:{x:'x0', z:'z1'} },
    { vary:'y', fix:{x:'x1', z:'z0'} }, { vary:'y', fix:{x:'x1', z:'z1'} },
    // Z-parallel edges
    { vary:'z', fix:{x:'x0', y:'y0'} }, { vary:'z', fix:{x:'x0', y:'y1'} },
    { vary:'z', fix:{x:'x1', y:'y0'} }, { vary:'z', fix:{x:'x1', y:'y1'} },
  ];

  function worldOf(x,y,z){ return mapper(x,y,z); }

  // Draw each edge line (optional)
  if (showEdges) {
    for (const e of edges){
      const p0 = { x: dom.x0, y: dom.y0, z: dom.z0 };
      const p1 = { x: dom.x0, y: dom.y0, z: dom.z0 };
      p0[e.vary] = dom[ e.vary + '0' ];
      p1[e.vary] = dom[ e.vary + '1' ];
      if (e.fix.x) { p0.x = dom[e.fix.x]; p1.x = dom[e.fix.x]; }
      if (e.fix.y) { p0.y = dom[e.fix.y]; p1.y = dom[e.fix.y]; }
      if (e.fix.z) { p0.z = dom[e.fix.z]; p1.z = dom[e.fix.z]; }
      const g = new THREE.BufferGeometry().setFromPoints([worldOf(p0.x,p0.y,p0.z), worldOf(p1.x,p1.y,p1.z)]);
      linesGroup.add(new THREE.Line(g, lineMat));
    }
  }

  // For each axis, put ticks on all 4 edges parallel to it
  function addTicksForAxis(axis){
    const vals = cfg[axis].ticks;
    const label = cfg[axis].label ?? ((v)=>String(v));
    // four edges for this axis:
    const combos = axis === 'x'
      ? [{y:'y0',z:'z0'},{y:'y0',z:'z1'},{y:'y1',z:'z0'},{y:'y1',z:'z1'}]
      : axis === 'y'
      ? [{x:'x0',z:'z0'},{x:'x0',z:'z1'},{x:'x1',z:'z0'},{x:'x1',z:'z1'}]
      : [{x:'x0',y:'y0'},{x:'x0',y:'y1'},{x:'x1',y:'y0'},{x:'x1',y:'y1'}];

    const cam = scene.userData?.camera;
    const camPos = new THREE.Vector3();
    if (cam?.getWorldPosition) cam.getWorldPosition(camPos);

    for (const v of vals){
      // choose only the edge combo closest to camera for this tick value
      let chosen = combos[0];
      if (cam) {
        let bestDist = Infinity;
        for (const c of combos) {
          const p = { x: dom.x0, y: dom.y0, z: dom.z0 };
          p[axis] = v;
          if (c.x) p.x = dom[c.x];
          if (c.y) p.y = dom[c.y];
          if (c.z) p.z = dom[c.z];
          const P = worldOf(p.x, p.y, p.z);
          const d = P.distanceToSquared(camPos);
          if (d < bestDist) {
            bestDist = d;
            chosen = c;
          }
        }
      }

      // base point (on chosen edge)
      const p = { x: dom.x0, y: dom.y0, z: dom.z0 };
      p[axis] = v;
      if (chosen.x) p.x = dom[chosen.x];
      if (chosen.y) p.y = dom[chosen.y];
      if (chosen.z) p.z = dom[chosen.z];
      const P = worldOf(p.x,p.y,p.z);

      // outward tick direction: from cube center to edge, normalized & scaled
      const center = worldOf( (dom.x0+dom.x1)/2, (dom.y0+dom.y1)/2, (dom.z0+dom.z1)/2 );
      const outward = new THREE.Vector3().subVectors(P, center).normalize().multiplyScalar(tickLen);

      // tick segment
      const tickP0 = tickOutwardOnly ? P.clone() : P.clone().sub(outward);
      const tickP1 = P.clone().add(outward);
      const tickDir = new THREE.Vector3().subVectors(tickP1, tickP0);
      const tickLenWorld = tickDir.length();
      if (tickLenWorld > 1e-9) {
        const tick = new THREE.Mesh(
          new THREE.CylinderGeometry(tickRadius, tickRadius, tickLenWorld, 8),
          tickMeshMat
        );
        tick.position.copy(tickP0.clone().add(tickP1).multiplyScalar(0.5));
        tick.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), tickDir.normalize());
        linesGroup.add(tick);
      }

      // crisp label using CSS2D
      const el = document.createElement('div');
      const labelText = label(v);
      el.innerHTML = labelText;
      el.style.color = labelColor;
      el.style.font = '30px ui-monospace, Menlo, Consolas, monospace';
      el.style.padding = '0';
      el.style.background = 'transparent';
      el.style.borderRadius = '0';
      const obj = new CSS2DObject(el);
      obj.position.copy(P.clone().add(outward.multiplyScalar(labelOffset))); // farther from axis
      labelsGroup.add(obj);
    }
  }

  addTicksForAxis('x');
  addTicksForAxis('y');
  addTicksForAxis('z');

  return { linesGroup, labelsGroup };
}



// ==================== GEOMETRY HELPERS ====================
export function buildBezierPath(points, tension = 0.6) {
  const path = new THREE.CurvePath();
  if (points.length === 2) {
    const [p0, p1] = points;
    const d = new THREE.Vector3().subVectors(p1, p0).multiplyScalar(1 / 3);
    path.add(new THREE.CubicBezierCurve3(p0, p0.clone().add(d), p1.clone().sub(d), p1));
    return path;
  }
  const tangents = points.map(() => new THREE.Vector3());
  for (let i = 0; i < points.length; i++) {
    const a = points[Math.max(0, i - 1)], b = points[Math.min(points.length - 1, i + 1)];
    tangents[i].copy(new THREE.Vector3().subVectors(b, a).multiplyScalar(0.5 * tension));
  }
  for (let i = 0; i < points.length - 1; i++) {
    const p0 = points[i], p1 = points[i + 1];
    const c0 = p0.clone().add(tangents[i].clone().multiplyScalar(1 / 3));
    const c1 = p1.clone().sub(tangents[i + 1].clone().multiplyScalar(1 / 3));
    path.add(new THREE.CubicBezierCurve3(p0.clone(), c0, c1, p1.clone()));
  }
  return path;
}

export function tubeFromPathZGradient(path, {
  radius = 0.10,
  tubularSegments = 200,
  radialSegments = 8,
  openEnded = false,
  colorStart = 0xffffff,
  colorEnd = 0x3da3ff
} = {}) {
  const geom = new THREE.TubeGeometry(path, tubularSegments, radius, radialSegments, openEnded);
  const pos = geom.getAttribute('position');
  const count = pos.count;
  let zMin = +Infinity, zMax = -Infinity;

  for (let i = 0; i < count; i++) {
    const z = pos.getZ(i);
    if (z < zMin) zMin = z;
    if (z > zMax) zMax = z;
  }
  const span = Math.max(1e-9, zMax - zMin);

  const c0 = new THREE.Color(colorStart);
  const c1 = new THREE.Color(colorEnd);
  const colors = new Float32Array(count * 3);

  for (let i = 0; i < count; i++) {
    const z = pos.getZ(i);
    const t = 1 - (z - zMin) / span;
    const r = c0.r + (c1.r - c0.r) * t;
    const g = c0.g + (c1.g - c0.g) * t;
    const b = c0.b + (c1.b - c0.b) * t;
    colors[i * 3 + 0] = r;
    colors[i * 3 + 1] = g;
    colors[i * 3 + 2] = b;
  }
  geom.setAttribute('color', new THREE.BufferAttribute(colors, 3));

  const mat = new THREE.MeshPhongMaterial({
    vertexColors: true,
    shininess: 30,
    specular: 0x666666,
    transparent: true
  });

  return new THREE.Mesh(geom, mat);
}

export function addBoundsBox(
  scene,
  size = 36,
  colors = { x: 0x555a6b, y: 0x555a6b, z: 0x555a6b },
  opts = {}
) {
  const yOnly = opts.yOnly ?? false;
  const topBottomOnly = opts.topBottomOnly ?? false;
  const lineWidth = opts.lineWidth ?? LINE_WIDTH;
  const opacity = opts.opacity ?? 1.0;
  const old = scene.getObjectByName('bounds-box');
  if (old) old.removeFromParent();

  const g = new THREE.Group();
  g.name = 'bounds-box';
  scene.add(g);

  const width = Number.isFinite(Number(opts.width)) ? Number(opts.width) : size;
  const height = Number.isFinite(Number(opts.height)) ? Number(opts.height) : size;
  const depth = Number.isFinite(Number(opts.depth)) ? Number(opts.depth) : size;
  const hx = width / 2;
  const hy = height / 2;
  const hz = depth / 2;
  const X = [-hx, hx], Y = [-hy, hy], Z = [-hz, hz];

  const matX = new THREE.LineBasicMaterial({ color: colors.x, transparent: true, opacity: 0.9 * opacity, linewidth: lineWidth });
  const matY = new THREE.LineBasicMaterial({
    color: colors.y,
    transparent: true,
    opacity: 1.0 * opacity,
    linewidth: lineWidth,
    depthTest: false,
    depthWrite: false
  });
  const matZ = new THREE.LineBasicMaterial({ color: colors.z, transparent: true, opacity: 0.95 * opacity, linewidth: lineWidth });

  function addEdge(p0, p1, mat) {
    const geom = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(...p0), new THREE.Vector3(...p1)]);
    const line = new THREE.Line(geom, mat);
    line.renderOrder = 20;
    g.add(line);
  }

  if (!topBottomOnly) {
    // Y-parallel edges
    addEdge([X[0], Y[0], Z[0]], [X[0], Y[1], Z[0]], matY);
    addEdge([X[1], Y[0], Z[0]], [X[1], Y[1], Z[0]], matY);
    addEdge([X[0], Y[0], Z[1]], [X[0], Y[1], Z[1]], matY);
    addEdge([X[1], Y[0], Z[1]], [X[1], Y[1], Z[1]], matY);
  }

  if (yOnly) {
    return g;
  }

  // X-parallel edges
  addEdge([X[0], Y[0], Z[0]], [X[1], Y[0], Z[0]], matX);
  addEdge([X[0], Y[0], Z[1]], [X[1], Y[0], Z[1]], matX);
  addEdge([X[0], Y[1], Z[0]], [X[1], Y[1], Z[0]], matX);
  addEdge([X[0], Y[1], Z[1]], [X[1], Y[1], Z[1]], matX);

  // Z-parallel edges
  addEdge([X[0], Y[0], Z[0]], [X[0], Y[0], Z[1]], matZ);
  addEdge([X[1], Y[0], Z[0]], [X[1], Y[0], Z[1]], matZ);
  addEdge([X[0], Y[1], Z[0]], [X[0], Y[1], Z[1]], matZ);
  addEdge([X[1], Y[1], Z[0]], [X[1], Y[1], Z[1]], matZ);

  return g;
}