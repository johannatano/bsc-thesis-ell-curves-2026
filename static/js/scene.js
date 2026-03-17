import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { addCubeEdgeAxes, removeCubeAxes } from '/static/js/geo.js';

export class SceneManager {
  constructor(leftContainer, rightContainer) {
    this.left = this.setupScene(leftContainer);
    this.right = this.setupScene(rightContainer);
    this.scenes = [this.left, this.right];
    
    this.mapCurveXYZ = null;
    this.mapTorsXYZ = null;
    
    this.P = null;
    this.N = null;
    this.q = null;

    this.boundsBoxColors = { x: 0x555a6b, y: 0x555a6b, z: 0x555a6b };
    this.boundsBoxSize = 36 * 0.9;
    this.rightBoundsOpacity = 0.90;
  }

  setBoundsBoxColors(colors = {}) {
    this.boundsBoxColors = {
      ...this.boundsBoxColors,
      ...colors
    };
    const oldRightBounds = this.right.scene.getObjectByName('bounds-box');
    if (oldRightBounds) oldRightBounds.removeFromParent();
  }

  setupScene(container) {
    const scene = new THREE.Scene();
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    container.appendChild(renderer.domElement);
    const camera = new THREE.PerspectiveCamera(60, 1, 0.1, 5000);
    camera.position.set(0, 0, 60);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    scene.userData.camera = camera;

    // Stronger base lighting for better overall visibility.
    scene.add(new THREE.AmbientLight(0xffffff, 1.2));

    // Main key light.
    const key = new THREE.DirectionalLight(0xffffff, 1.3);
    key.position.set(-1, -1, -1);
    scene.add(key);

    // Fill lights to reduce harsh dark sides.
    const fillA = new THREE.DirectionalLight(0xffffff, 1.3);
    fillA.position.set(1, -0.5, -0.8);
    scene.add(fillA);

    const fillB = new THREE.DirectionalLight(0xffffff, 1.3);
    fillB.position.set(-0.3, 1, 0.6);
    scene.add(fillB);

    function resize() {
      const r = container.getBoundingClientRect();
      renderer.setSize(r.width, r.height, false);
      camera.aspect = r.width / r.height;
      camera.updateProjectionMatrix();
    }
    resize();

    return { scene, camera, renderer, controls, resize };
  }

  init(P, N, q, box = 36, pad = 0.9, preserveCamera = false) {
    this.P = P;
    this.N = N;
    this.q = q;

    this.boundsBoxSize = box * pad;

    // Setup mappers
    this.setupMappers(P, N, q, box, pad, preserveCamera);

    // Ensure no bounds box is shown in the right scene.
    const oldRightBounds = this.right.scene.getObjectByName('bounds-box');
    if (oldRightBounds) oldRightBounds.removeFromParent();

    // Setup axes
    removeCubeAxes(this.left.scene);
    removeCubeAxes(this.right.scene);

    const tMin = -2 * Math.sqrt(q), tMax = 2 * Math.sqrt(q);
    addCubeEdgeAxes(this.left.scene, this.mapCurveXYZ, {
      x: { domain: [0, P - 1], ticks: [...Array(P).keys()], label: v => v },
      y: { domain: [tMin, tMax], ticks: [-2, -1, 0, 1, 2], label: v => v },
      z: { domain: [0, N - 1], ticks: [...Array(N).keys()], label: v => v },
      tickLen: 1.0, color: 0x6e7387
    });

    const X0 = -Math.floor(P / 2), X1 = Math.ceil(P / 2) - 1;
    const rangeInts = (a, b) => { const out = []; for (let v = a; v <= b; v++) out.push(v); return out; };
    addCubeEdgeAxes(this.right.scene, this.mapTorsXYZ, {
      x: { domain: [X0, X1], ticks: rangeInts(X0, X1), label: v => v },
      y: { domain: [X0, X1], ticks: rangeInts(X0, X1), label: v => v },
      z: { domain: [0, N - 1], ticks: [...Array(N).keys()], label: v => v },
      tickLen: 1.0, color: 0x6e7387, showEdges: false
    });
  }

  setupMappers(P, N, q, box = 36, pad = 0.9, preserveCamera = false) {
    // CURVES cube
    {
      const X0 = 0, X1 = P - 1;
      const Y0 = -2 * Math.sqrt(q), Y1 = 2 * Math.sqrt(q);
      const Z0 = 0, Z1 = N - 1;
      const cx = (X0 + X1) / 2, cy = (Y0 + Y1) / 2, cz = (Z0 + Z1) / 2;
      const sx = (box * pad) / Math.max(1e-9, X1 - X0);
      const sy = (box * pad) / Math.max(1e-9, Y1 - Y0);
      const sz = (box * pad) / Math.max(1e-9, Z1 - Z0);
      this.mapCurveXYZ = (x, y, z) => new THREE.Vector3((x - cx) * sx, (y - cy) * sy, -(z - cz) * sz);
    }

    // TORSION cube
    {
      const X0 = -Math.floor(P / 2), X1 = Math.ceil(P / 2) - 1;
      const Y0 = X0, Y1 = X1;
      const Z0 = 0, Z1 = N - 1;
      const cx = (X0 + X1) / 2;
      const cy = (Y0 + Y1) / 2;
      const cz = (Z0 + Z1) / 2;
      const sx = (box * pad) / Math.max(1e-9, X1 - X0);
      const sy = (box * pad) / Math.max(1e-9, Y1 - Y0);
      const sz = (box * pad) / Math.max(1e-9, Z1 - Z0);

      this.mapTorsXYZ = (xCoeff, yCoeff, layerIndex) => {
        const X = this.symModP(xCoeff, P);
        const Y = this.symModP(yCoeff, P);
        const Z = layerIndex;
        return new THREE.Vector3(
          (X - cx) * sx,
          (Y - cy) * sy,
          -((Z - Z0) * sz - (sz * (Z1 - Z0)) / 2)
        );
      };
    }

    if (!preserveCamera) {
      const diag = Math.sqrt((box * pad) ** 2 * 3);
      this.left.camera.position.set(0, 0, diag * 1.2);
      this.left.controls.target.set(0, 0, 0);
      this.right.camera.position.set(0, 0, diag * 1.2);
      this.right.controls.target.set(0, 0, 0);
      this.left.camera.updateProjectionMatrix();
      this.right.camera.updateProjectionMatrix();
    }
  }

  symModP(c, p) {
    c = ((c % p) + p) % p;
    const half = Math.floor(p / 2);
    return (c <= half) ? c : (c - p);
  }

  resize() {
    this.scenes.forEach(s => s.resize());
  }

  animate(labelRenderer = null) {
    this.scenes.forEach(s => {
      s.controls.update();
      s.renderer.render(s.scene, s.camera);
    });
    
    // Render CSS2D labels for left scene
    if (labelRenderer) {
      labelRenderer.renderer.render(this.left.scene, this.left.camera);
    }
  }

  startAnimation(labelRenderer = null, rightLabelRenderer = null) {
    const loop = () => {
      requestAnimationFrame(loop);
      this.animate(labelRenderer);
      
      // Render right scene labels separately
      if (rightLabelRenderer) {
        rightLabelRenderer.renderer.render(this.right.scene, this.right.camera);
      }
    };
    loop();
  }
}