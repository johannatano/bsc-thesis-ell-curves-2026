import * as THREE from 'three';
import { buildBezierPath, tubeFromPathZGradient } from '/static/js/geo.js';

export class EllipticCurveRenderer {
  constructor(P, N, q) {
    this.P = P;
    this.N = N;
    this.q = q;
    this.curves = [];
  }

  setCurves(curves) {
    this.curves = curves;
  }

  renderCurves(div, mapCurveXYZ, opts = {}) {
    const R_MAIN = opts.radius ?? 0.1;
    const items = [];
    const meshes = [];

    this.curves.forEach(([A, B, Ncur, j, torsion], idx) => {
      let trank = 0;
      const t = this.q + 1 - Ncur;
      if (div > 1) {
        const entry = torsion?.[String(div)];
        const list = entry?.gens ?? [];
        if (!Array.isArray(list) || !list.length) return;
        trank = 1;
        if (list.length >= 2) {
          trank = 2;
        }
      }
      items.push({ A, B, Ncur, j, torsion, t, idx, trank });
    });

    const byT = new Map();
    for (const it of items) {
      if (!byT.has(it.t)) byT.set(it.t, []);
      byT.get(it.t).push(it);
    }

    const rankMap = new Map();
    for (const [tVal, arr] of byT) {
      arr.sort((a, b) => {
        const A = a.j, B = b.j;
        const L = Math.max(A.length, B.length);
        for (let i = 0; i < L; i++) {
          const av = A[i] ?? 0, bv = B[i] ?? 0;
          if (av !== bv) return av - bv;
        }
        return 0;
      });
      const count = arr.length;
      arr.forEach((it, i) => {
        rankMap.set(`${tVal}|${it.j.join(',')}`, { rank: i, count });
      });
    }

    for (const it of items) {
      const { A, B, Ncur, j, t, torsion, idx, trank } = it;

      const xs = [...j, ...new Array(this.N).fill(0)].slice(0, this.N);
      const ys = new Array(this.N).fill(t);
      const zs = [...Array(this.N).keys()];
      const pts = xs.map((x, i) => mapCurveXYZ(x, ys[i], zs[i]));

      const key = `${t}|${j.join(',')}`;
      const { rank, count } = rankMap.get(key) || { rank: 0, count: 1 };

      let col1 = trank > 1 ? new THREE.Color(0, 1, 0) : new THREE.Color(1, 1, 1);
      let col2 = trank > 1 ? new THREE.Color(0, 1, 0) : new THREE.Color(1, 1, 1);

      if (t % this.P == 0) {
        col1 = trank > 1 ? new THREE.Color(0, 1, 0) : new THREE.Color(1, 0, 0);
        col2 = new THREE.Color(1, 0, 0);
      }

      if (this.N === 1) {
        const pos = mapCurveXYZ(j[0] ?? 0, t, 0);
        const sphere = new THREE.Mesh(
          new THREE.SphereGeometry(R_MAIN * 1.8, 24, 24),
          new THREE.MeshPhongMaterial({
            color: col1,
            transparent: true,
            opacity: 0.95,
            shininess: 50,
            specular: 0x999999
          })
        );
        sphere.userData = { idx, j, N: Ncur, t, torsion, A, B };
        sphere.position.copy(pos);
        meshes.push(sphere);
      } else {
        const path = buildBezierPath(pts, 0.6);
        const mesh = tubeFromPathZGradient(path, {
          radius: R_MAIN,
          colorStart: col1,
          colorEnd: col2
        });
        mesh.material.transparent = true;
        mesh.material.depthWrite = false;
        mesh.material.opacity = trank > 1 ? 1 : 1;
        mesh.userData = { idx, j, N: Ncur, t, torsion, A, B };
        meshes.push(mesh);
      }
    }

    return meshes;
  }

  renderTorsionPoints(mesh, div, mapTorsXYZ, opts = {}) {
    if (!mesh?.userData?.torsion) return null;

    const entry = mesh.userData.torsion[String(div)];
    const ptsJson = entry?.points ?? [];
    if (!ptsJson.length) return null;

    const group = new THREE.Group();
    const radius = opts.radius ?? 0.1;
    const tension = opts.tension ?? 0.6;
    const opacity = opts.opacity ?? 0.95;
    const baseCol = (opts.color instanceof THREE.Color) ? opts.color : new THREE.Color(opts.color ?? 0x8bd0ff);

    ptsJson.forEach((pt, i) => {
      if (this.N === 1) {
        const pos = mapTorsXYZ(pt[0] ?? 0, pt[1] ?? 0, 0);
        const sphere = new THREE.Mesh(
          new THREE.SphereGeometry(radius * 1.5, 16, 16),
          new THREE.MeshPhongMaterial({
            color: baseCol,
            transparent: true,
            opacity
          })
        );
        sphere.position.copy(pos);
        group.add(sphere);
      } else {
        const pts = Array.from({ length: this.N }, (_, k) =>
          mapTorsXYZ(pt[0][k] ?? 0, pt[1][k] ?? 0, k)
        );
        const path = buildBezierPath(pts, tension);
        const tubeMesh = tubeFromPathZGradient(path, {
          radius,
          colorStart: 0xffffff,
          colorEnd: 0x3da3ff
        });
        tubeMesh.material.transparent = true;
        tubeMesh.material.depthWrite = false;
        tubeMesh.material.opacity = opacity;
        group.add(tubeMesh);
      }
    });

    return group;
  }

  computeStats(div) {
    let totalOrbits = 0;
    let traces = [];

    for (const [A, B, Ncur, j, torsion] of this.curves) {
      const entry = torsion?.[String(div)];
      if (!entry) continue;
      const n_orbits = entry?.n_orbits ?? 0;
      if (n_orbits > 0) {
        const t = this.q + 1 - Ncur;
        if (!traces.includes(t)) {
          traces.push(t);
        }
        totalOrbits += n_orbits;
      }
    }

    return { totalOrbits, traces };
  }

  getSelectedOrbits(selectedMesh, div) {
    if (!selectedMesh?.userData) return 0;

    const { j } = selectedMesh.userData;
    for (const [A, B, Ncur, jCurve, torsion] of this.curves) {
      if (j.join(', ') === jCurve.join(', ')) {
        const entry = torsion?.[String(div)];
        return entry?.n_orbits ?? 0;
      }
    }
    return 0;
  }
}
