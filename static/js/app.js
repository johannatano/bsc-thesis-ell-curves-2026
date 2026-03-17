import * as THREE from 'three';
import { SceneManager } from '/static/js/scene.js';
import { EllipticCurveRenderer } from '/static/js/elliptic.js';
import { UI } from '/static/js/ui.js';
import { parseNumberFieldData, mergeCurvesIntoFields } from '/static/js/data.js';
import { buildBezierPath, tubeFromPathZGradient, createLabelRenderer, addCubeEdgeAxes, removeCubeAxes, LINE_WIDTH } from '/static/js/geo.js';
import { CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

// ==================== CONFIG ====================
const urlParams = new URLSearchParams(window.location.search);
const P = Number(urlParams.get('p')) || 13;
let N = Number(urlParams.get('n')) || 2;
const _q = P ** N;

// URL parameters for auto-selection
const URL_D = urlParams.get('D') ? parseInt(urlParams.get('D')) : null;
const URL_J = urlParams.get('j') || null;

// ==================== THEMES ====================
const THEMES = {
  dark: {
    rendererBg: 0x0e0e10,
    bodyBg: '#0e0e10',
    themeColor1: '#000000',
    themeColor2: '#ff00ff',
    sphere: {
      supersingular: new THREE.Color(1, 0, 0),
      eisenstein:    new THREE.Color(0.7, 0.7, 0.7),
      gaussian:      new THREE.Color(0.7, 0.7, 0.7),
      ordinary:      new THREE.Color(0.7, 0.7, 0.7),
      piConductor:   new THREE.Color(0x3da3ff),
      aboveFloor:    new THREE.Color(0x3da3ff),
      floorLevel:    new THREE.Color(0.7, 0.7, 0.7),
      hover:         new THREE.Color(1, 0, 1),
      selected:      new THREE.Color(0xffff00),
    },
    tube:   { colorStart: 0xffffff, colorEnd: 0xffffff },
    axis:   0x6e7387,
    uiBg:   'rgba(0,0,0,0.45)',
    uiColor:'#ddd',
    labelBg:'rgba(0,0,0,0.4)',
    labelColor: '#aaa',
    inputBg: '#2a2a30',
    inputColor: '#ddd',
    inputBorder: '#444',
    pillBg: '#1f1f26',
  },
  light: {
    rendererBg: 0xffffff,
    bodyBg: '#ffffff',
    themeColor1: '#000000',
    themeColor2: '#ff00ff',
    sphere: {
      supersingular: new THREE.Color(0.8, 0.1, 0.1),
      eisenstein:    new THREE.Color(.7, .7, .7),
      gaussian:      new THREE.Color(.7, .7, .7),
      ordinary:      new THREE.Color(1.0, 1.0, 1.0),
      piConductor:   new THREE.Color(0x3da3ff),
      aboveFloor:    new THREE.Color(0x1a5bbf),
      floorLevel:    new THREE.Color(0.9, 0.9, 0.9),
      hover:         new THREE.Color(0.7, 0, 0.7),
      selected:      new THREE.Color(0.8, 0.6, 0),
    },
    tube:   { colorStart: 0x222222, colorEnd: 0x222222 },
    axis:   0x444444,
    uiBg:   'rgba(255,255,255,0.92)',
    uiColor:'#111',
    labelBg:'rgba(240,240,240,0.9)',
    labelColor: '#333',
    inputBg: '#f0f0f0',
    inputColor: '#111',
    inputBorder: '#aaa',
    pillBg: '#e0e0e0',
  },
};
let currentTheme = 'light';

// ==================== STATE ====================
const SELECTED = new Set();

const OPQ_ACTIVE = 0.1;   // 10% opacity for all planes when something selected
const OPQ_BASE = 0;       // Invisible by default
const OPQ_DIM = 0;        // Not used anymore
const SHOW_ORDER_CONNECTOR_TUBES = false;
const INACTIVE_DOT_OPACITY = 1.0;
const GUIDE_RADIUS = 0.5 * LINE_WIDTH;
const GUIDE_TICK_RADIUS = .25 * LINE_WIDTH;

let CURVES;
let FIELDS;

let fieldData = null;
let q = null;

let sceneManager = null;
let ellipticRenderer = null;
let ui = null;

// Layout mode: 'circle' or 'scatter'
let layoutMode = 'circle';

const nrFieldsGroup = new THREE.Group();
let curveMeshes = [];
let numberFieldGroups = [];  // Array of {group, discriminant}
let hiddenObjects = [];  // Track objects removed during selection
let hoveredDiscriminant = null;  // Track which field is being hovered
let hoveredConductor = null;  // Track conductor of hovered dot
let selectedDiscriminant = null;  // Track the selected discriminant across N changes
let selectedVolcanoTrace = null;  // Trace for the volcano explicitly clicked in status panel
let originalColors = new Map();  // Store original colors for hover restoration
let fieldInfoEl = null;  // Field info display element
let lastActiveSceneKey = 'left';

const selectedFieldGroup = new THREE.Group();
const discriminantLabelsGroup = new THREE.Group();
discriminantLabelsGroup.name = 'discriminant-labels';

// Crosshair for left box
const leftCrosshairGroup = new THREE.Group();
leftCrosshairGroup.name = 'left-crosshair';

// Store global conductor normalization (used for positioning and axis labels)
let globalMinLogCond = 0;
let globalMaxLogCond = 1;
let globalLogCondRange = 1;

let labelRenderer = null;

// Toggle connected-component spatial splitting in volcano view.
// false: keep edge-aware ordering, but do not separate components in XZ.
const DISTRIBUTE_COMPONENTS = false;

// ==================== SCENE SETUP ====================
sceneManager = new SceneManager(
  document.getElementById('left'),
  document.getElementById('right')
);
sceneManager.setBoundsBoxColors({
  x: 0x555a6b,
  y: 0x555a6b,
  z: 0x555a6b
});
// Apply initial theme renderer background
sceneManager.scenes.forEach(s => s.renderer.setClearColor(THEMES[currentTheme].rendererBg, 1));

sceneManager.left.scene.add(nrFieldsGroup);
sceneManager.left.scene.add(discriminantLabelsGroup);
sceneManager.left.scene.add(leftCrosshairGroup);
sceneManager.right.scene.add(selectedFieldGroup);


// ==================== DATA FETCHING ====================


async function loadFields(p) {
  const url = `/fields/${p}`;
  const r = await fetch(url, { cache: "no-store" });
  if (!r.ok) throw new Error(`${url} → ${r.status} ${r.statusText}`);
  const json = await r.json();
  return parseNumberFieldData(json);
}

async function loadCurves(p, n) {
  const url = `/curves/${p}/${n}`;
  const r = await fetch(url, { cache: "no-store" });
  if (!r.ok) throw new Error(`${url} → ${r.status} ${r.statusText}`);
  return r.json();
}

// ==================== CURVE RENDERING ====================
// ==================== NUMBER FIELD RENDERING ====================

// ==================== THEME ====================
function applyTheme(name) {
  currentTheme = name;
  const t = THEMES[name];

  // Renderer backgrounds
  sceneManager.scenes.forEach(s => {
    s.renderer.setClearColor(t.rendererBg, 1);
  });
  sceneManager.setBoundsBoxColors({
    x: 0x555a6b,
    y: 0x555a6b,
    z: 0x555a6b
  });

  // CSS body + UI panel
  document.body.style.background = t.bodyBg;
  const uiEl = document.getElementById('ui');
  if (uiEl) {
    uiEl.style.background = t.uiBg;
    uiEl.style.color = t.uiColor;
  }
  // Update all inputs/selects/buttons inside UI
  document.querySelectorAll('#ui input, #ui select, #ui button, #ui label').forEach(el => {
    el.style.background = t.inputBg;
    el.style.color = t.inputColor;
    el.style.borderColor = t.inputBorder;
  });
  document.querySelectorAll('#ui .pill').forEach(el => {
    el.style.background = t.pillBg;
    el.style.color = t.uiColor;
  });
  document.querySelectorAll('.label').forEach(el => {
    el.style.background = t.labelBg;
    el.style.color = t.labelColor;
  });
  if (fieldInfoEl) {
    fieldInfoEl.style.color = t.labelColor;
  }

  // Update theme button label
  const btn = document.getElementById('theme_toggle');
  if (btn) btn.textContent = name === 'dark' ? '☀ Light' : '☾ Dark';

  // Rebuild the scene so all meshes/tubes use the new palette
  if (typeof rebuildScene === 'function' && typeof ui !== 'undefined' && ui) {
    ui.setPiConductorColor?.(THEMES[currentTheme].sphere.piConductor.getStyle());
    rebuildScene(ui.activeDiv);
    if (selectedDiscriminant !== null) drawCurvesForField(selectedDiscriminant);
  }
}

function rebuildScene(div) {
  nrFieldsGroup.clear();
  discriminantLabelsGroup.clear();
  curveMeshes.length = 0;
  numberFieldGroups.length = 0;

  SELECTED.clear();
  selectedFieldGroup.clear();
  
  // Clear CSS2D labels from DOM
  if (labelRenderer && labelRenderer.domElement) {
    while (labelRenderer.domElement.firstChild) {
      labelRenderer.domElement.removeChild(labelRenderer.domElement.firstChild);
    }
  }

  // Render planes for each number field
  const numberFields = fieldData;
  const getAbsTraceClasses = (classes = []) => {
    const seen = new Set();
    return classes.filter(ic => {
      const key = Math.abs(Number(ic.trace));
      if (!Number.isFinite(key)) return true;
      if (seen.has(key)) return false;
      seen.add(key);
      return true;
    });
  };
  
  // Get discriminant range for normalization
  const discriminants = numberFields.map(nf => Number(nf.discriminant));
  const minDisc = Math.min(...discriminants);
  const maxDisc = Math.max(...discriminants);
  
  // Fixed maximum discriminant for consistent scaling across N values
  const FIXED_MAX_DISC = 10000000;
  const Z_LOG_COMPRESS = 1.0; // >1 compresses z-spread
  
  // Apply log scale to discriminants (use absolute values since discriminants are negative)
  // log(|D| + 1) to handle D=0 case
  const logDiscriminants = discriminants.map(d => Math.log(Math.abs(d) + 1));
  const maxLogDisc = Math.log(FIXED_MAX_DISC + 1) * Z_LOG_COMPRESS;
  
  // Get trace t range for x-axis normalization
  const xValues = [];
  const qVal = Number.isFinite(q) ? q : (P ** N);
  numberFields.forEach(nf => {
    getAbsTraceClasses(nf.isogeny_classes).forEach(ic => {
      xValues.push(Number(ic.trace));
    });
  });
  const minX = Math.min(...xValues);
  const maxX = Math.max(...xValues);
  const xRange = maxX - minX || 1;
  
  // Get conductor range for normalization
  const conductors = [];
  numberFields.forEach(nf => {
    getAbsTraceClasses(nf.isogeny_classes).forEach(ic => {
      if (ic.orders) {
        ic.orders.forEach(order => {
          if (order.conductor !== undefined) conductors.push(order.conductor);
        });
      }
    });
  });
  const minCond = Math.min(...conductors);
  const maxCond = Math.max(...conductors);
  
  // Apply log scale to conductors for better spacing
  const logConductors = conductors.map(c => Math.log(c + 1));
  const minLogCond = Math.log(minCond + 1);
  const maxLogCond = Math.log(maxCond + 1);
  const logCondRange = maxLogCond - minLogCond || 1;
  
  // Store global conductor normalization for axis labels
  globalMinLogCond = minLogCond;
  globalMaxLogCond = maxLogCond;
  globalLogCondRange = logCondRange;

  numberFields.forEach(nf => {
    const discriminant = Number(nf.discriminant);
    // Create a group for this number field
    const fieldGroup = new THREE.Group();
    fieldGroup.userData = { discriminant };
    numberFieldGroups.push({ group: fieldGroup, discriminant });
    
    // Use fixed box dimensions (the scene box is 36 units with 0.9 padding)
    const boxSize = 36 * 0.9;
    const width = boxSize;
    let height = boxSize; // Make height equal to width for proper rectangles
    
    // Apply log scale to discriminant and remap to [-0.5, 0.5] range
    // D=0 -> back (-0.5), D=-FIXED_MAX_DISC -> front (0.5)
    const logDisc = Math.log(Math.abs(discriminant) + 1);
    const zValue = -0.5 + (logDisc / maxLogDisc);

    // Scale directly by height (don't use mapCurveXYZ which may apply additional scaling)
    const zPos = zValue * height;

    let baseRad = 0.15;
    // Draw conductor spheres and tubes for each isogeny class
    getAbsTraceClasses(nf.isogeny_classes).forEach(ic => {
      const trace = ic.trace;
      const dPi = Number.isFinite(Number(ic.D_pi))
        ? Number(ic.D_pi)
        : (Number(trace) ** 2 - 4 * qVal);
      const fPiRaw = Number(ic.f_pi);
      const fPiFromOrders = (ic.orders && Array.isArray(ic.orders) && ic.orders.length > 0)
        ? Math.max(...ic.orders.map(o => Number(o.conductor)).filter(Number.isFinite))
        : NaN;
      const fPiTarget = Number.isFinite(fPiRaw) ? fPiRaw : fPiFromOrders;
      
      // Normalize trace t to fit within plane width
      const xPos = (Number(trace) - minX) / xRange * width - width/2;
      
      // Mark conductors of all orders
      if (ic.orders && Array.isArray(ic.orders)) {
        ic.orders.forEach(order => {
        const conductor = order.conductor;
        // Apply log scale to conductor for Y position (inverted: higher conductor = lower Y)
        const logCond = Math.log(conductor + 1);
        const normalizedY = -((logCond - minLogCond) / logCondRange * boxSize - boxSize/2);  // Inverted to match right box
        
        // q-normalized step scaling (non-logarithmic):
        // low q -> small step threshold (e.g. +1 curve per step),
        // large q -> larger threshold (e.g. +10 curves per step).
        const curveCount = Math.max(1, Number(order.class_number) || 1);
        const qScale = Math.max(1, Number(qVal) || 1);
        const stepThreshold = Math.max(1, Math.round(Math.sqrt(qScale) / 10));
        const scaleIncr = 1.0;
        const scaleLevel = Math.floor((curveCount - 1) / stepThreshold);
        const sphereRadius = Math.min(baseRad * (1 + scaleLevel * scaleIncr), 2.0);
        
        const isPiConductor = Number.isFinite(fPiTarget) && Number(conductor) === fPiTarget;
        let sphereColor;
        const sp = THEMES[currentTheme].sphere;
        if (!ic.ordinary) {
          sphereColor = sp.supersingular.clone();
        } else if (discriminant === -3) {
          sphereColor = sp.eisenstein.clone();
        } else if (discriminant === -4) {
          sphereColor = sp.gaussian.clone();
        } else {
          sphereColor = sp.ordinary.clone();
        }
        const orderGeometry = new THREE.SphereGeometry(sphereRadius, 16, 16);
        const orderMaterial = new THREE.MeshPhongMaterial({
          color: sphereColor,
          transparent: true,
          opacity: INACTIVE_DOT_OPACITY,
          shininess: 30
        });
        
        const orderSphere = new THREE.Mesh(orderGeometry, orderMaterial);
        orderSphere.position.set(xPos, normalizedY, zPos);
        orderSphere.userData = {
          discriminant: nf.discriminant,
          trace: ic.trace,
          D_pi: dPi,
          conductor: conductor,
          piConductor: fPiTarget,
          isPiConductor,
          type: 'order',
          class_number: order.class_number
        };
        
        fieldGroup.add(orderSphere);
        curveMeshes.push(orderSphere);
        });
      }
      
      // Draw bezier tubes connecting all conductor spheres in this isogeny class
      const conductorPositions = [];
      // Add all order conductor positions
      if (ic.orders && Array.isArray(ic.orders)) {
        ic.orders.forEach(order => {
          const logCond = Math.log(order.conductor + 1);
          const normalizedY = -((logCond - minLogCond) / logCondRange * boxSize - boxSize/2);
          conductorPositions.push(new THREE.Vector3(xPos, normalizedY, zPos));
        });
      }
      // Sort by y position for cleaner connections
      conductorPositions.sort((a, b) => a.y - b.y);
      // Create bezier tube if we have at least 2 points
      if (conductorPositions.length >= 2 && SHOW_ORDER_CONNECTOR_TUBES) {
        const path = buildBezierPath(conductorPositions, 0.3);
        const tubeMesh = tubeFromPathZGradient(path, {
          radius: 0.02,
          tubularSegments: 100,
          radialSegments: 6,
          colorStart: THEMES[currentTheme].tube.colorStart,
          colorEnd: THEMES[currentTheme].tube.colorEnd
        });
        tubeMesh.material.transparent = false;
        tubeMesh.material.opacity = 1.0;
        fieldGroup.add(tubeMesh);
      }
    });
    
    // Add the complete number field group to the scene
    nrFieldsGroup.add(fieldGroup);
  });

  // Add axes to left scene without conductor ticks (will be added on hover/select)
  const boxSize = 36 * 0.9;
  const halfBox = boxSize / 2;
  
  removeCubeAxes(sceneManager.left.scene);
  addCubeEdgeAxes(sceneManager.left.scene, (x, y, z) => new THREE.Vector3(x, y, z), {
    x: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    y: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    z: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    tickLen: 1.0,
    showEdges: false,
    color: THEMES[currentTheme].axis
  });

  // No default conductor planes; grid is drawn only for highlighted conductor.

  // Re-select the previously selected discriminant if it exists in current data  if (selectedDiscriminant !== null) {
    const meshToSelect = curveMeshes.find(m => Number(m.userData.discriminant) === selectedDiscriminant);
    if (meshToSelect) {
      SELECTED.add(meshToSelect);
      highlightField(
        selectedDiscriminant,
        Number(meshToSelect.userData.piConductor ?? meshToSelect.userData.conductor),
        true,
        Number(meshToSelect.userData.D_pi)
      );
      showFieldInfo(selectedDiscriminant);
      //drawTorsionForMesh(meshToSelect, div);
    }
    // Keep selectedDiscriminant stored even if not found, so it can be re-selected later
  }

// ==================== SELECTION ====================
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

// Track selected curves in right box (can be multiple with same j)
let selectedCurves = [];
let hoveredCurves = [];

// Helper function to compare j-invariants (arrays or primitives)
function jInvariantsEqual(j1, j2) {
  if (Array.isArray(j1) && Array.isArray(j2)) {
    // Compare arrays ignoring trailing zeros
    // Trim trailing zeros from both arrays
    const trim = (arr) => {
      let lastNonZero = -1;
      for (let i = 0; i < arr.length; i++) {
        if (arr[i] !== 0) lastNonZero = i;
      }
      return arr.slice(0, lastNonZero + 1);
    };
    
    const trimmed1 = trim(j1);
    const trimmed2 = trim(j2);
    
    if (trimmed1.length !== trimmed2.length) return false;
    return trimmed1.every((val, idx) => val === trimmed2[idx]);
  }
  return j1 === j2;
}

// Find all curves with the same j-invariant
function findCurvesWithJ(targetJ) {
  const matching = [];
  selectedFieldGroup.traverse(obj => {
    if (obj.isMesh && obj.userData.curveID && obj.userData.j) {
      if (jInvariantsEqual(obj.userData.j, targetJ)) {
        matching.push(obj);
      }
    }
  });
  return matching;
}

function getThemeColor1() {
  return new THREE.Color(THEMES[currentTheme].themeColor1);
}

function getCurveBaseColorByHeight(height) {
  return height > 0
    ? THEMES[currentTheme].sphere.aboveFloor.clone()
    : THEMES[currentTheme].sphere.floorLevel.clone();
}

function setCurveVisualState(curve, state = 'base') {
  if (!curve?.material?.color) return;
  if (state === 'selected' || state === 'hover') {
    curve.material.color.copy(getThemeColor1());
    return;
  }
  const height = curve.userData?.height || 0;
  curve.material.color.copy(getCurveBaseColorByHeight(height));
}

function updateVolcanoEdgeHighlights() {
  const selectedIDs = new Set((selectedCurves || []).map(c => c?.userData?.curveID).filter(Boolean));
  const hoveredIDs = new Set((hoveredCurves || []).map(c => c?.userData?.curveID).filter(Boolean));
  const activeIDs = selectedIDs.size > 0 ? selectedIDs : hoveredIDs;

  const activeColor = getThemeColor1();

  selectedFieldGroup.traverse(obj => {
    if (!(obj instanceof THREE.Line)) return;
    if (!obj.userData?.edgeFrom || !obj.userData?.edgeTo || !obj.material?.color) return;

    const baseColor = obj.userData.baseEdgeColor ?? 0x444444;
    const baseOpacity = obj.userData.baseEdgeOpacity ?? 0.4;
    const isActive = activeIDs.size > 0 && (activeIDs.has(obj.userData.edgeFrom) || activeIDs.has(obj.userData.edgeTo));

    if (isActive) {
      obj.material.color.copy(activeColor);
      obj.material.opacity = 0.95;
    } else {
      obj.material.color.setHex(baseColor);
      obj.material.opacity = baseOpacity;
    }
  });
}

function onClick(ev) {
  lastActiveSceneKey = 'left';
  const rect = sceneManager.left.renderer.domElement.getBoundingClientRect();
  mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera(mouse, sceneManager.left.camera);
  const hits = raycaster.intersectObjects(curveMeshes, false);
  if (!hits.length) return;

  const mesh = hits[0].object;
  toggleSelection(mesh);
}

function onMouseMove(ev) {
  lastActiveSceneKey = 'left';
  const rect = sceneManager.left.renderer.domElement.getBoundingClientRect();
  
  // Check if mouse is within the left canvas
  if (ev.clientX >= rect.left && ev.clientX <= rect.right && 
      ev.clientY >= rect.top && ev.clientY <= rect.bottom) {
    
    mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
    
    // Update crosshair
    raycaster.setFromCamera(mouse, sceneManager.left.camera);
    
    // Create an invisible plane at the center of the scene to raycast against
    const boxSize = 36 * 0.9;
    const planeGeometry = new THREE.PlaneGeometry(boxSize, boxSize);
    const planeMaterial = new THREE.MeshBasicMaterial({ visible: false });
    const plane = new THREE.Mesh(planeGeometry, planeMaterial);
    plane.rotation.y = Math.PI / 2; // Orient plane perpendicular to camera view
    
    const planeHits = raycaster.intersectObject(plane);
    
    if (planeHits.length > 0) {
      const hitPoint = planeHits[0].point;
      updateLeftCrosshair(hitPoint);
    } else {
      clearLeftCrosshair();
    }
    
    // Only do mesh hover detection if nothing is selected
    if (SELECTED.size === 0) {
      const hits = raycaster.intersectObjects(curveMeshes, false);
      const newHoveredDiscriminant = hits.length > 0 ? Number(hits[0].object.userData.discriminant) : null;
      const newHoveredConductor = hits.length > 0
        ? Number(hits[0].object.userData.piConductor ?? hits[0].object.userData.conductor)
        : null;
      
      if (newHoveredDiscriminant !== hoveredDiscriminant || newHoveredConductor !== hoveredConductor) {
        // Restore previous hover
        if (hoveredDiscriminant !== null && newHoveredDiscriminant !== hoveredDiscriminant) {
          restoreFieldColors(hoveredDiscriminant);
        }
        
        // Apply new hover
        if (newHoveredDiscriminant !== null) {
          if (newHoveredDiscriminant !== hoveredDiscriminant) {
            highlightField(newHoveredDiscriminant, newHoveredConductor, false);
          } else {
            updateLeftAxisForField(newHoveredDiscriminant, newHoveredConductor);
          }
          showFieldInfo(newHoveredDiscriminant);
        } else {
          hideFieldInfo();
        }
        
        hoveredDiscriminant = newHoveredDiscriminant;
        hoveredConductor = newHoveredConductor;
      }
    }
  } else {
    // Mouse outside left canvas - clear crosshair
    clearLeftCrosshair();
  }
}

// ==================== CROSSHAIR ====================
function updateLeftCrosshair(point) {
  clearLeftCrosshair();
  
  const boxSize = 36 * 0.9;
  const halfBox = boxSize / 2;
  
  // Main horizontal line at z=0 along X axis at Y height
  const horizontalLine = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(-halfBox, point.y, 0),
    new THREE.Vector3(halfBox, point.y, 0)
  ]);
  
  // Line along Z axis at left edge (x = -halfBox)
  const leftZLine = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(-halfBox, point.y, -halfBox),
    new THREE.Vector3(-halfBox, point.y, halfBox)
  ]);
  
  // Line along Z axis at right edge (x = +halfBox)
  const rightZLine = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(halfBox, point.y, -halfBox),
    new THREE.Vector3(halfBox, point.y, halfBox)
  ]);
  
  const lineMaterial = new THREE.LineBasicMaterial({ 
    color: 0xff00ff, 
    opacity: 0.5, 
    transparent: true,
    linewidth: 3
  });
  
  const hLine = new THREE.Line(horizontalLine, lineMaterial);
  const lLine = new THREE.Line(leftZLine, lineMaterial);
  const rLine = new THREE.Line(rightZLine, lineMaterial);
  
  leftCrosshairGroup.add(hLine);
  leftCrosshairGroup.add(lLine);
  leftCrosshairGroup.add(rLine);
}

function clearLeftCrosshair() {
  while (leftCrosshairGroup.children.length > 0) {
    const child = leftCrosshairGroup.children[0];
    if (child.geometry) child.geometry.dispose();
    if (child.material) child.material.dispose();
    leftCrosshairGroup.remove(child);
  }
}

// ==================== CURVE SELECTION (RIGHT BOX) ====================
function onCurveClick(ev) {
  lastActiveSceneKey = 'right';
  const rect = sceneManager.right.renderer.domElement.getBoundingClientRect();
  mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera(mouse, sceneManager.right.camera);
  const curveMeshes = [];
  selectedFieldGroup.traverse(obj => {
    if (obj.isMesh && obj.userData.curveID) {
      curveMeshes.push(obj);
    }
  });
  
  const hits = raycaster.intersectObjects(curveMeshes, false);
  if (!hits.length) {
    // Clicked on empty space - do nothing (don't clear selection)
    return;
  }

  const mesh = hits[0].object;
  const targetJ = mesh.userData.j;
  
  // Find all curves with the same j-invariant
  const matchingCurves = findCurvesWithJ(targetJ);
  
  // Check if already selected
  const isAlreadySelected = selectedCurves.some(c => jInvariantsEqual(c.userData.j, targetJ));
  
  if (isAlreadySelected) {
    // Deselect all with this j
    selectedCurves.forEach(curve => {
      setCurveVisualState(curve, 'base');
    });
    selectedCurves = [];
    updateVolcanoEdgeHighlights();
  } else {
    // Clear previous selection
    selectedCurves.forEach(curve => {
      setCurveVisualState(curve, 'base');
    });
    // Select all with this j
    selectedCurves = matchingCurves;
    matchingCurves.forEach(curve => {
      setCurveVisualState(curve, 'selected');
    });
    updateVolcanoEdgeHighlights();
    console.log(`Selected ${matchingCurves.length} curve(s) with j = ${mesh.userData.jFormatted}:`, 
                matchingCurves.map(c => c.userData.curveID));
  }
}

function onCurveMouseMove(ev) {
  lastActiveSceneKey = 'right';
  if (selectedCurves.length > 0) return; // Don't hover if something is selected
  
  const rect = sceneManager.right.renderer.domElement.getBoundingClientRect();
  mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera(mouse, sceneManager.right.camera);
  const curveMeshes = [];
  selectedFieldGroup.traverse(obj => {
    if (obj.isMesh && obj.userData.curveID) {
      curveMeshes.push(obj);
    }
  });
  
  const hits = raycaster.intersectObjects(curveMeshes, false);
  
  const newHoveredMesh = hits.length > 0 ? hits[0].object : null;
  
  // Check if we're hovering over a different j-invariant
  const hoveredJChanged = !newHoveredMesh || hoveredCurves.length === 0 || 
                          !jInvariantsEqual(newHoveredMesh.userData.j, hoveredCurves[0].userData.j);
  
  if (hoveredJChanged) {
    // Restore previous hover
    hoveredCurves.forEach(curve => {
      setCurveVisualState(curve, 'base');
    });
    hoveredCurves = [];
    
    // Apply new hover
    if (newHoveredMesh) {
      const matchingCurves = findCurvesWithJ(newHoveredMesh.userData.j);
      hoveredCurves = matchingCurves;
      matchingCurves.forEach(curve => {
        setCurveVisualState(curve, 'hover');
      });
      const aFormatted = formatJInvariant(newHoveredMesh.userData.a, P);
      const bFormatted = formatJInvariant(newHoveredMesh.userData.b, P);
      console.log(`Hovering ${matchingCurves.length} curve(s) with j=${newHoveredMesh.userData.jFormatted}, A=${aFormatted}, B=${bFormatted}:`,
                  matchingCurves.map(c => c.userData.curveID));
    }

    updateVolcanoEdgeHighlights();
  }
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

let _html2canvasPromise = null;
function ensureHtml2Canvas() {
  if (window.html2canvas) return Promise.resolve(window.html2canvas);
  if (_html2canvasPromise) return _html2canvasPromise;

  _html2canvasPromise = new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js';
    script.async = true;
    script.onload = () => {
      if (window.html2canvas) resolve(window.html2canvas);
      else reject(new Error('html2canvas failed to initialize'));
    };
    script.onerror = () => reject(new Error('Failed to load html2canvas'));
    document.head.appendChild(script);
  });

  return _html2canvasPromise;
}

async function exportHighDpiScene(side = 'active', scale = 2, mimeType = 'image/png') {
  const which = side === 'active' ? lastActiveSceneKey : side;
  const slot = which === 'right' ? sceneManager.right : sceneManager.left;
  if (!slot?.renderer || !slot?.scene || !slot?.camera) return;

  const panel = document.getElementById(which === 'right' ? 'right' : 'left');
  const canvas = slot.renderer.domElement;
  const cssW = Math.max(1, Math.floor(canvas.clientWidth));
  const cssH = Math.max(1, Math.floor(canvas.clientHeight));
  const outW = Math.max(1, Math.floor(cssW * scale));
  const outH = Math.max(1, Math.floor(cssH * scale));
  const ext = mimeType === 'image/png' ? 'png' : 'webp';

  const outCanvas = document.createElement('canvas');
  outCanvas.width = outW;
  outCanvas.height = outH;
  const ctx = outCanvas.getContext('2d');
  if (!ctx) return;

  // 1) Draw current WebGL frame (scene geometry)
  slot.renderer.render(slot.scene, slot.camera);
  ctx.drawImage(canvas, 0, 0, outW, outH);

  try {
    // 2) Draw HTML/CSS overlays (labels/tooltips), but ignore canvases
    const h2c = await ensureHtml2Canvas();
    const overlayCanvas = await h2c(panel || canvas, {
      scale,
      useCORS: true,
      backgroundColor: null,
      logging: false,
      ignoreElements: (el) => (el?.tagName || '').toLowerCase() === 'canvas'
    });

    ctx.drawImage(overlayCanvas, 0, 0);

    const ts = new Date().toISOString().replace(/[:.]/g, '-');
    const filename = `scene-${which}-${outCanvas.width}x${outCanvas.height}-${ts}.${ext}`;
    outCanvas.toBlob((blob) => {
      if (!blob) return;
      downloadBlob(blob, filename);
      console.log(`Exported ${which} scene (with labels): ${filename}`);
    }, mimeType, 0.98);
    return;
  } catch (err) {
    console.warn('Overlay capture failed; exporting scene only:', err);
  }

  const ts = new Date().toISOString().replace(/[:.]/g, '-');
  const filename = `scene-${which}-${outW}x${outH}-${ts}.${ext}`;
  outCanvas.toBlob((blob) => {
    if (!blob) return;
    downloadBlob(blob, filename);
    console.log(`Exported ${which} scene (fallback): ${filename}`);
  }, mimeType, 0.98);
}

function drawConductorGrid(scene, {
  boxSize,
  boxWidth,
  boxHeight,
  minLogCond,
  logCondRange,
  conductor,
  groupName,
  color,
  outlineOnly = false
}) {
  const old = scene.getObjectByName(groupName);
  if (old) old.removeFromParent();

  const conductorsToDraw = Array.isArray(conductor)
    ? conductor.map(x => Number(x)).filter(Number.isFinite)
    : (Number.isFinite(Number(conductor)) ? [Number(conductor)] : []);
  if (conductorsToDraw.length === 0) return;

  const g = new THREE.Group();
  g.name = groupName;
  scene.add(g);

  const width = Number.isFinite(Number(boxWidth)) ? Number(boxWidth) : boxSize;
  const height = Number.isFinite(Number(boxHeight)) ? Number(boxHeight) : boxSize;
  const halfW = width / 2;
  const halfH = height / 2;

  if (outlineOnly) {
    const outlineMat = new THREE.LineBasicMaterial({
      color,
      transparent: true,
      opacity: 0.75,
      depthWrite: false
    });

    for (const c of conductorsToDraw) {
      const logCond = Math.log(c + 1);
      const y = -((logCond - minLogCond) / logCondRange * height - halfH);
      if (!Number.isFinite(y) || y < -halfH || y > halfH) continue;

      const outlineGeom = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-halfW, y, -halfW),
        new THREE.Vector3( halfW, y, -halfW),
        new THREE.Vector3( halfW, y,  halfW),
        new THREE.Vector3(-halfW, y,  halfW)
      ]);
      const outline = new THREE.LineLoop(outlineGeom, outlineMat);
      outline.renderOrder = 6;
      g.add(outline);
    }
    return;
  }

  // Build a visibly thick grid (LineBasicMaterial linewidth is ignored on most WebGL drivers).
  const divisions = 12;
  const step = width / divisions;
  const gridRadius = .25 * LINE_WIDTH;
  const gridMat = new THREE.MeshBasicMaterial({
    color,
    //color: new THREE.Color(0.0),
    transparent: true,
    opacity: 0.55,
    depthWrite: false
  });

  for (const c of conductorsToDraw) {
    const logCond = Math.log(c + 1);
    const y = -((logCond - minLogCond) / logCondRange * height - halfH);
    if (!Number.isFinite(y) || y < -halfH || y > halfH) continue;

    for (let i = 0; i <= divisions; i++) {
      const coord = -halfW + i * step;

      // Line parallel to X at z = coord
      const lineX = new THREE.Mesh(
        new THREE.CylinderGeometry(gridRadius, gridRadius, width, 8),
        gridMat
      );
      lineX.position.set(0, y, coord);
      lineX.rotation.z = Math.PI / 2;
      lineX.renderOrder = 6;
      g.add(lineX);

      // Line parallel to Z at x = coord
      const lineZ = new THREE.Mesh(
        new THREE.CylinderGeometry(gridRadius, gridRadius, width, 8),
        gridMat
      );
      lineZ.position.set(coord, y, 0);
      lineZ.rotation.x = Math.PI / 2;
      lineZ.renderOrder = 6;
      g.add(lineZ);
    }
  }
}

function clearLeftGuidePlanes() {
  const roof = sceneManager.left.scene.getObjectByName('left-roof-ok-outline');
  if (roof) {
    roof.traverse((node) => {
      if (node instanceof CSS2DObject && node.element?.parentNode) {
        node.element.parentNode.removeChild(node.element);
      }
    });
    roof.removeFromParent();
  }

  const cond = sceneManager.left.scene.getObjectByName('selected-conductor-grid-left');
  if (cond) cond.removeFromParent();
}

function updateLeftAxisForField(discriminant, highlightedConductor = null) {
  if (!fieldData) return;
  
  const nf = fieldData.find(f => Number(f.discriminant) === discriminant);
  if (!nf) return;
  
  // Get conductors for this specific field
  const conductors = [];
  nf.isogeny_classes.forEach(ic => {
    if (ic.orders) {
      ic.orders.forEach(order => {
        if (order.conductor !== undefined) conductors.push(order.conductor);
      });
    }
  });
  
  if (conductors.length === 0) return;
  
  const boxSize = 36 * 0.9;
  const halfBox = boxSize / 2;
  
  // Show only highlighted conductor tick if provided, else all
  let uniqueConductors = [...new Set(conductors)].sort((a, b) => a - b);
  if (highlightedConductor !== null && Number.isFinite(Number(highlightedConductor))) {
    const hc = Number(highlightedConductor);
    uniqueConductors = uniqueConductors.includes(hc) ? [hc] : [];
  }
  if (uniqueConductors.length === 0) {
    clearLeftAxisTicks();
    return;
  }
  
  // Use GLOBAL normalization (not field-specific)
  const conductorTicks = uniqueConductors.map(c => {
    const logCond = Math.log(c + 1);
    return -((logCond - globalMinLogCond) / globalLogCondRange * boxSize - boxSize/2);
  });
  
  removeCubeAxes(sceneManager.left.scene);
  addCubeEdgeAxes(sceneManager.left.scene, (x, y, z) => new THREE.Vector3(x, y, z), {
    x: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    y: { 
      domain: [-halfBox, halfBox], 
      ticks: [],
      label: () => ''
    },
    z: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    tickLen: 4.8, 
    labelOffset: 1.5,
    labelColor: THEMES[currentTheme].themeColor1,
    tickColor: THEMES[currentTheme].themeColor1,
    tickOutwardOnly: true,
    showEdges: false,
    color: THEMES[currentTheme].axis
  });

  clearLeftGuidePlanes();
}

function clearLeftAxisTicks() {
  const boxSize = 36 * 0.9;
  const halfBox = boxSize / 2;
  
  removeCubeAxes(sceneManager.left.scene);
  addCubeEdgeAxes(sceneManager.left.scene, (x, y, z) => new THREE.Vector3(x, y, z), {
    x: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    y: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    z: { domain: [-halfBox, halfBox], ticks: [], label: v => '' },
    tickLen: 1.0,
    showEdges: false,
    color: 0x6e7387
  });

  clearLeftGuidePlanes();
}

function drawSelectedDKVerticalTubes(discriminant, highlightedConductor = null, highlightedDpi = null) {
  const old = sceneManager.left.scene.getObjectByName('selected-dk-vertical-tubes');
  if (old) {
    old.traverse((node) => {
      if (node instanceof CSS2DObject && node.element?.parentNode) {
        node.element.parentNode.removeChild(node.element);
      }
    });
    old.removeFromParent();
  }

  const entry = numberFieldGroups.find(({ discriminant: d }) => d === discriminant);
  if (!entry?.group) return;

  const g = new THREE.Group();
  g.name = 'selected-dk-vertical-tubes';
  sceneManager.left.scene.add(g);

  const boxSize = 36 * 0.9;
  const halfBox = boxSize / 2;
  const yMin = -halfBox;
  const yMax = halfBox;
  const tubeHeight = yMax - yMin;

  const mat = new THREE.MeshBasicMaterial({
    color: THEMES[currentTheme].sphere.hover.getHex(),
    transparent: true,
    opacity: 0.9,
    depthWrite: false
  });

  const seen = new Set();
  const selectedPoints = [];
  const allOrderPoints = [];
  entry.group.traverse((obj) => {
    if (!obj.isMesh || obj.userData?.type !== 'order') return;
    const conductor = Number(obj.userData?.conductor);
    const dPi = Number(obj.userData?.D_pi);
    allOrderPoints.push({
      x: obj.position.x,
      y: obj.position.y,
      z: obj.position.z,
      conductor,
      dPi
    });

    const x = obj.position.x;
    const z = obj.position.z;
    const key = `${x.toFixed(4)}|${z.toFixed(4)}`;
    if (seen.has(key)) return;
    seen.add(key);
    selectedPoints.push({
      x,
      y: obj.position.y,
      z,
      dPi,
      conductor
    });
  });

  selectedPoints.forEach(({ x, z }) => {
    const dPiTubeExtraHeight = 4.8;
    const dPiTubeHeight = tubeHeight + dPiTubeExtraHeight;
    const dPiTopY = yMax + dPiTubeExtraHeight;

    const tube = new THREE.Mesh(
      new THREE.CylinderGeometry(GUIDE_RADIUS, GUIDE_RADIUS, dPiTubeHeight, 16),
      mat
    );
    tube.position.set(x, yMin + dPiTubeHeight / 2, z);
    tube.renderOrder = 8;
    g.add(tube);

    // Keep D_K value on top (as before).
    const labelDiv = document.createElement('div');
    labelDiv.style.color = '#ff00ff';
    labelDiv.style.fontSize = '30px';
    labelDiv.style.fontFamily = 'ui-monospace, monospace';
    labelDiv.style.background = 'transparent';
    labelDiv.style.padding = '0';
    labelDiv.style.borderRadius = '0';
    labelDiv.style.pointerEvents = 'none';
    labelDiv.innerHTML = `D<sub>K</sub>: ${discriminant}`;

    const label = new CSS2DObject(labelDiv);
    label.position.set(x, dPiTopY + 4.5, z);
    g.add(label);
  });

  const anchorUseCounts = new Map();
  const getNearestFaceDirection = (anchor) => {
    const faceCandidates = [
      { dist: halfBox - anchor.x, dir: new THREE.Vector3(1, 0, 0) },
      { dist: anchor.x + halfBox, dir: new THREE.Vector3(-1, 0, 0) },
      { dist: halfBox - anchor.z, dir: new THREE.Vector3(0, 0, 1) },
      { dist: anchor.z + halfBox, dir: new THREE.Vector3(0, 0, -1) }
    ];
    let best = faceCandidates[0];
    for (let i = 1; i < faceCandidates.length; i++) {
      if (faceCandidates[i].dist < best.dist) best = faceCandidates[i];
    }
    return {
      direction: best.dir,
      distanceToFace: Math.max(0, best.dist)
    };
  };

  const addDotTooltip = (point, html, color, keyTag = '', opts = {}) => {
    if (!point) return;
    const anchor = new THREE.Vector3(point.x, point.y, point.z);

    const anchorKey = `${anchor.x.toFixed(3)}|${anchor.y.toFixed(3)}|${anchor.z.toFixed(3)}`;
    const usage = anchorUseCounts.get(anchorKey) ?? 0;
    anchorUseCounts.set(anchorKey, usage + 1);

    const seed = `${keyTag}|${anchorKey}`;
    let hash = 0;
    for (let i = 0; i < seed.length; i++) hash = ((hash * 31) + seed.charCodeAt(i)) % 104729;
    const phase = (hash / 104729) * Math.PI * 2;

    const nearestFace = getNearestFaceDirection(anchor);
    const tickDir = nearestFace.direction.clone().normalize();
    const minOutsideDistance = 1.2;
    const requestedTickLen = Number.isFinite(Number(opts.tickLen)) ? Number(opts.tickLen) : 4.8;
    const tickLen = Math.max(requestedTickLen, nearestFace.distanceToFace + minOutsideDistance);

    // Stable side direction for label separation.
    let side = new THREE.Vector3().crossVectors(tickDir, new THREE.Vector3(0, 1, 0));
    if (side.lengthSq() < 1e-9) side = new THREE.Vector3(1, 0, 0);
    side.normalize();
    const sidePhase = 0.9 * Math.sin(phase + usage * 1.1);

    const tickMat = new THREE.MeshBasicMaterial({
      color: new THREE.Color(color),
      transparent: true,
      opacity: 0.95,
      depthWrite: false
    });

    const tick = new THREE.Mesh(
      new THREE.CylinderGeometry(GUIDE_TICK_RADIUS, GUIDE_TICK_RADIUS, tickLen, 12),
      tickMat
    );
    tick.position.copy(anchor.clone().add(tickDir.clone().multiplyScalar(tickLen / 2)));
    tick.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), tickDir);
    tick.renderOrder = 8;
    g.add(tick);

    const labelDiv = document.createElement('div');
    labelDiv.style.color = color;
    labelDiv.style.fontSize = '30px';
    labelDiv.style.fontFamily = 'ui-monospace, monospace';
    labelDiv.style.background = 'transparent';
    labelDiv.style.padding = '0';
    labelDiv.style.borderRadius = '0';
    labelDiv.style.pointerEvents = 'none';
    labelDiv.innerHTML = html;

    const label = new CSS2DObject(labelDiv);
    const labelOffset = tickDir.clone().multiplyScalar(tickLen + 4.5)
      .add(side.clone().multiplyScalar((usage * 1.4) + sidePhase));
    label.position.copy(anchor.clone().add(labelOffset));
    g.add(label);
  };

  const f1Candidates = allOrderPoints.filter(p => Number(p.conductor) === 1);
  const f1Point = f1Candidates.length > 0
    ? f1Candidates.reduce((best, p) => (p.y > best.y ? p : best), f1Candidates[0])
    : null;

  const fPi = Number(highlightedConductor);
  const fPiCandidates = Number.isFinite(fPi)
    ? allOrderPoints.filter(p => Number(p.conductor) === fPi)
    : [];
  const targetDpi = Number(highlightedDpi);
  const fPiPoint = (Number.isFinite(targetDpi)
    ? fPiCandidates.find(p => Number(p.dPi) === targetDpi)
    : null)
    || (fPiCandidates.length > 0
      ? fPiCandidates.reduce((best, p) => (p.y > best.y ? p : best), fPiCandidates[0])
      : null);

  addDotTooltip(f1Point, '𝒪<sub>K</sub>', THEMES[currentTheme].themeColor2, 'ok');
  addDotTooltip(fPiPoint, 'Z[π]', THEMES[currentTheme].themeColor2, 'zpi', { baseLift: -0.35, tickLen: 7.2 });

  // D_K guide plane/tick label temporarily disabled.
}

function highlightField(discriminant, highlightedConductor = null, selectedOnlyTube = false, highlightedDpi = null) {
  const magentaColor = new THREE.Color(1, 0, 1);
  
  // Update left axis to show conductors for this field
  updateLeftAxisForField(discriminant, highlightedConductor);

  const oldTubes = sceneManager.left.scene.getObjectByName('selected-dk-vertical-tubes');
  if (oldTubes) {
    oldTubes.traverse((node) => {
      if (node instanceof CSS2DObject && node.element?.parentNode) {
        node.element.parentNode.removeChild(node.element);
      }
    });
    oldTubes.removeFromParent();
  }
  if (selectedOnlyTube) {
    drawSelectedDKVerticalTubes(discriminant, highlightedConductor, highlightedDpi);
  }
  
  numberFieldGroups.forEach(({ group, discriminant: groupDisc }) => {
    if (groupDisc === discriminant && group) {
      group.traverse(obj => {
        if (obj.material && obj.material.color) {
          if (obj.userData?.type === 'order' && obj.scale?.setScalar) {
            obj.scale.setScalar(selectedOnlyTube ? 1.2 : 1.0);
          }

          // Store original color
          if (!originalColors.has(obj.uuid)) {
            originalColors.set(obj.uuid, obj.material.color.clone());
          }
          const isGuideDot = selectedOnlyTube
            && obj.userData?.type === 'order'
            && (
              Number(obj.userData.conductor) === 1
              || (
                Number.isFinite(Number(highlightedConductor))
                && Number(obj.userData.conductor) === Number(highlightedConductor)
              )
            );

          if (obj.userData?.type === 'order') {
            obj.material.transparent = true;
            obj.material.opacity = isGuideDot ? 1.0 : INACTIVE_DOT_OPACITY;
          }

          if (isGuideDot) {
            obj.material.color.set(THEMES[currentTheme].themeColor2);
          } else {
            obj.material.color.copy(magentaColor);
          }
        }
      });
    }
  });
}

function restoreFieldColors(discriminant) {
  numberFieldGroups.forEach(({ group, discriminant: groupDisc }) => {
    if (groupDisc === discriminant && group) {
      group.traverse(obj => {
        if (obj.userData?.type === 'order' && obj.scale?.setScalar) {
          obj.scale.setScalar(1.0);
        }
        if (obj.userData?.type === 'order' && obj.material) {
          obj.material.transparent = true;
          obj.material.opacity = INACTIVE_DOT_OPACITY;
        }
        if (obj.material && obj.material.color) {
          const originalColor = originalColors.get(obj.uuid);
          if (originalColor) {
            obj.material.color.copy(originalColor);
            originalColors.delete(obj.uuid);
          }
        }
      });
    }
  });
  
  // Clear conductor ticks when unhighlighting
  clearLeftAxisTicks();

  const oldTubes = sceneManager.left.scene.getObjectByName('selected-dk-vertical-tubes');
  if (oldTubes) {
    oldTubes.traverse((node) => {
      if (node instanceof CSS2DObject && node.element?.parentNode) {
        node.element.parentNode.removeChild(node.element);
      }
    });
    oldTubes.removeFromParent();
  }
}

function showFieldInfo(discriminant) {

  if (!fieldInfoEl || !fieldData) return;
  const nf = fieldData.find(f => Number(f.discriminant) === discriminant);
  if (!nf) return;
  
  const qDisplay = Number.isFinite(Number(q)) ? Number(q) : (P ** N);
  let info = `<strong><span style="font-size: 32px;">𝔽<sub>${qDisplay}</sub></span></strong><br><strong>D = ${discriminant}</strong>`;
  if (nf.name) info += ` (${nf.name})`;
  
  // Get trace values and D_pi values
  const traces = nf.isogeny_classes.map(ic => ic.trace).sort((a, b) => a - b);
  const uniqueTraces = [...new Set(traces)];
  
  info += `<br>${nf.isogeny_classes.length} isogeny class${nf.isogeny_classes.length !== 1 ? 'es' : ''}`;
  info += ` (t = ${uniqueTraces.join(', ')})`;
  
  // Show D_pi values for each isogeny class
  const dpiValues = nf.isogeny_classes
    .filter(ic => ic.D_pi !== undefined)
    .map(ic => `t=${ic.trace}: D_π=${ic.D_pi}`)
    .join(', ');
  if (dpiValues) {
    info += `<br>${dpiValues}`;
  }
  
  // Group by conductor and collect class numbers
  const conductorMap = new Map();
  nf.isogeny_classes.forEach(ic => {
    ic.orders.forEach(order => {
      const conductor = order.conductor;
      if (!conductorMap.has(conductor)) {
        conductorMap.set(conductor, []);
      }
      conductorMap.get(conductor).push(order.class_number);
    });
  });
  
  if (conductorMap.size > 0) {
    info += `<br><br><strong>Conductors:</strong>`;
    const sortedConductors = Array.from(conductorMap.keys()).sort((a, b) => a - b);
    sortedConductors.forEach(conductor => {
      const classNumbers = conductorMap.get(conductor);
      const uniqueClassNumbers = [...new Set(classNumbers)].sort((a, b) => a - b);
      info += `<br>f = ${conductor}: num curves = ${uniqueClassNumbers.join(', ')}`;
    });
  }
  
  fieldInfoEl.innerHTML = info;
  fieldInfoEl.style.display = 'block';

  ui.updateStatus(nf);
}

function hideFieldInfo() {
  if (fieldInfoEl) {
    fieldInfoEl.style.display = 'none';
  }
  ui.updateStatus(null);
}

async function ensureCurvesLoaded(discriminant) {
  // Find the number field
  const nf = fieldData.find(f => Number(f.discriminant) === discriminant);
  if (!nf) return;
  
  // Check if any isogeny class has curves loaded
  const hasCurves = nf.isogeny_classes.some(ic => 
    ic.orders && ic.orders.some(o => o.curves && o.curves.length > 0)
  );
  
  if (hasCurves) {
    return; // Already loaded
  }
  
  // Fetch curves via HCP API
  console.log(`Fetching curves for D=${discriminant} via HCP API...`);
  
  try {
    console.log(`Requesting /get_curves/${P}/${N}?D=${discriminant}`);
    const response = await fetch(`/get_curves/${P}/${N}?D=${discriminant}`);
    if (!response.ok) {
      throw new Error(`API returned ${response.status}`);
    }
    
    const data = await response.json();
    
    if (!data.success) {
      throw new Error(data.error || 'Unknown error');
    }
    // Use existing mergeCurvesIntoFields function - same format as precomputed JSON!
    mergeCurvesIntoFields(FIELDS, data.catalogue, N);
    
    console.log(`Successfully loaded curves for D=${discriminant}`);
  } catch (err) {
    console.error(`Failed to load curves for D=${discriminant}:`, err);
  }
}

async function toggleSelection(mesh) {
  // Clear any hover state
  if (hoveredDiscriminant !== null) {
    restoreFieldColors(hoveredDiscriminant);
    hoveredDiscriminant = null;
  }
  
  if (SELECTED.has(mesh)) {
    clearAllSelections();
    return;
  }

  clearAllSelections();
  SELECTED.add(mesh);
  selectedDiscriminant = Number(mesh.userData.discriminant);

  // Just highlight the selected field with magenta, don't hide others
  highlightField(
    selectedDiscriminant,
    Number(mesh.userData.piConductor ?? mesh.userData.conductor),
    true,
    Number(mesh.userData.D_pi)
  );

  updateEllButtons(selectedDiscriminant);
  
  // Check if curves are loaded for this field, if not fetch via API
  await ensureCurvesLoaded(selectedDiscriminant);
  
  drawCurvesForField(selectedDiscriminant);
  showFieldInfo(selectedDiscriminant);
}

function clearAllSelections() {
  SELECTED.forEach(m => {
    const discriminant = Number(m.userData.discriminant);
    restoreFieldColors(discriminant);
  });
  SELECTED.clear();
  selectedDiscriminant = null;
  hoveredConductor = null;

  // Clear curve visualizations
  selectedFieldGroup.clear();

  // Reset ell buttons to show all ells
  updateEllButtons();

  hideFieldInfo();
}


// ==================== RESIZE ====================
addEventListener('resize', () => {
  sceneManager.resize();
  if (labelRenderer) {
    labelRenderer.resize();
  }
});

// ==================== INITIALIZATION ====================
(async function init() {
  // Create label renderers for both scenes
  labelRenderer = createLabelRenderer(document.getElementById('left'));
  const rightLabelRenderer = createLabelRenderer(document.getElementById('right'));
  
  // Create field info display element
  fieldInfoEl = document.createElement('div');
  fieldInfoEl.id = 'field-info';
  fieldInfoEl.style.cssText = `
    position: absolute;
    left: 12px;
    bottom: 12px;
    background: transparent;
    color: ${THEMES[currentTheme].labelColor};
    padding: 0;
    border-radius: 0;
    font-family: ui-monospace, Menlo, Consolas, monospace;
    font-size: 22px;
    line-height: 1.4;
    max-width: 400px;
    display: none;
    z-index: 30;
    pointer-events: none;
  `;
  const leftPanel = document.getElementById('left');
  if (leftPanel) {
    leftPanel.appendChild(fieldInfoEl);
  } else {
    document.body.appendChild(fieldInfoEl);
  }
  FIELDS = await loadFields(P);

  console.log(`Loaded fields data for p=${P}, n=${N}:`, FIELDS);
  
  // Try to load curves data and merge with fields
  try {
    const curvesJson = await loadCurves(P, N);
    mergeCurvesIntoFields(FIELDS, curvesJson, N);
  } catch (err) {
    console.warn(`Curves data not available for p=${P}, n=${N}:`, err.message);
  }
  
  // Extract only the data for current N from all fields
  q = P**N;
  
  // Get filtered data for current N
  fieldData = FIELDS.catalogue.getFieldsByN(N);
  
  for (const nf of fieldData) {
    nf.isogeny_classes.forEach(ic => {
      /*for (const v of ic.volcanoes?.values() || []) {
        // Check volcano with h > 0
        const h = v.getMaxHeight();
        if (h > 0) {
          console.log(`  Volcano found with h=${h} at ell=${v.ell} for discriminant D=${nf.discriminant}`);
        }
      }*/
    });
  }

  

  ui = new UI(P, N);
    ui.setPiConductorColor?.(THEMES[currentTheme].sphere.piConductor.getStyle());
  ui.onDivChange = (div) => rebuildScene(div);
  ui.onNChange = async (newN) => {
    try {
      // Update global N variable
      N = newN;
      
      // Save current input values before clearing
      const currentD = ui.dInput ? ui.dInput.value : null;
      const currentJ = ui.jInput ? ui.jInput.value : null;
      
      // Clear all selections when switching N
      clearAllSelections();
      selectedCurves = [];
      hoveredCurves = [];
      
      // Load curves data for new N
      try {
        const curvesJson = await loadCurves(P, newN);
        mergeCurvesIntoFields(FIELDS, curvesJson, newN);
      } catch (err) {
        console.warn(`Curves data not available for p=${P}, n=${newN}:`, err.message);
      }
      
      // Get filtered data for new N
      fieldData = FIELDS.catalogue.getFieldsByN(newN);
      q = P ** newN;
      sceneManager.init(P, newN, q, 36, 0.9, true);
      updateEllButtons();
      rebuildScene(ui.activeDiv);
      
      // Re-apply selections after new data is loaded
      setTimeout(() => {
        if (currentJ) {
          // Prioritize j-invariant search
          ui.onJSelect(currentJ);
        } else if (currentD) {
          // Fall back to D if no j
          const d = parseInt(currentD);
          if (!isNaN(d)) {
            const nf = fieldData.find(field => Number(field.discriminant) === d);
            if (nf) {
              ui.onDSelect(d);
            } else {
              console.warn(`Discriminant D=${d} not found in new N=${newN}, clearing D input`);
              if (ui.dInput) ui.dInput.value = '';
            }
          }
        }
      }, 150);

    } catch (err) {
      console.error('Failed to update view for n=' + newN, err);
      ui.showError?.('Failed to update view for n=' + newN + '\n' + (err?.message || err));
    }
  };
  ui.onEllChange = (ell) => {
    selectedVolcanoTrace = null;
    // Redraw curves for the currently selected field with new ell filter
    if (selectedDiscriminant !== null) {
      drawCurvesForField(selectedDiscriminant);
    }
  };
  
  // Listen for volcano row clicks
  window.addEventListener('volcanoRowClick', (e) => {
    const { ell, trace } = e.detail;
    ui.selectedEll = ell;
    ui.ellSelect.value = ell;
    selectedVolcanoTrace = trace;
    if (selectedDiscriminant !== null) {
      drawCurvesForField(selectedDiscriminant, trace);
    }
  });
  
  ui.onDSelect = (d) => {
    // Find the number field with this discriminant
    const nf = fieldData.find(field => Number(field.discriminant) === d);
    if (nf) {
      // Find a mesh belonging to this field to trigger selection
      const fieldMesh = curveMeshes.find(mesh => Number(mesh.userData.discriminant) === d);
      if (fieldMesh) {
        // Clear previous selection and select this field
        clearAllSelections();
        toggleSelection(fieldMesh);
      } else {
        console.warn('No mesh found for discriminant D =', d);
      }
    } else {
      console.warn('No number field found with discriminant D =', d);
    }
  };
  
  ui.onJSelect = (jStr) => {
    // Parse j-invariant - could be an integer or array notation [a,b,c]
    let targetJ = null;
    
    // Try parsing as array notation like [1,2,3]
    if (jStr.startsWith('[') && jStr.endsWith(']')) {
      try {
        const arr = jStr.slice(1, -1).split(',').map(s => parseInt(s.trim()));
        if (arr.every(n => !isNaN(n))) {
          targetJ = arr;
        }
      } catch (e) {
        console.warn('Failed to parse j-invariant as array:', e);
      }
    } else {
      // Try parsing as integer - convert to base-p representation
      const jInt = parseInt(jStr);
      if (!isNaN(jInt)) {
        targetJ = intToBaseP(jInt, P);
      }
    }
    
    if (targetJ) {
      // Search through ALL fields to find curves with this j-invariant
      let foundField = null;
      let foundCurves = [];
      let checkedCount = 0;
      
      for (const nf of fieldData) {
        for (const ic of nf.isogeny_classes) {
          if (ic.orders) {
            for (const order of ic.orders) {
              if (order.curves) {
                for (const curve of order.curves) {
                  checkedCount++;
                  if (jInvariantsEqual(curve.j, targetJ)) {
                    foundCurves.push({
                      discriminant: nf.discriminant,
                      curveID: curve.ID,  // Use uppercase ID to match userData.curveID
                      j: curve.j
                    });
                    if (!foundField) foundField = nf.discriminant;
                  }
                }
              }
            }
          }
        }
      }
      if (foundCurves.length > 0) {
        // Check if we need to switch to a different field
        if (selectedDiscriminant !== foundField) {
          console.log(`Switching from field D=${selectedDiscriminant} to D=${foundField}`);
        }
        // Select the field containing these curves (will switch if different)
        const fieldMesh = curveMeshes.find(mesh => Number(mesh.userData.discriminant) === foundField);
        if (fieldMesh) {
          clearAllSelections();
          toggleSelection(fieldMesh);
          
          // Wait a bit for drawCurvesForField to complete, then highlight the matching curves
          setTimeout(() => {
            const curvesToHighlight = [];
            selectedFieldGroup.traverse(obj => {
              if (obj.isMesh && obj.userData.curveID) {
                const matchingCurve = foundCurves.find(fc => fc.curveID === obj.userData.curveID);
                if (matchingCurve) {
                  curvesToHighlight.push(obj);
                }
              }
            });
            
            // Clear previous curve selection
            selectedCurves.forEach(curve => {
              setCurveVisualState(curve, 'base');
            });
            
            // Highlight matching curves
            curvesToHighlight.forEach(curve => {
              setCurveVisualState(curve, 'selected');
            });
            selectedCurves = curvesToHighlight;
            updateVolcanoEdgeHighlights();
          }, 100);
        } else {
          console.warn(`Could not find mesh for discriminant D=${foundField}`);
        }
      } else {
        console.warn('No curves found with j-invariant:', targetJ, 'in any number field');
      }
    } else {
      console.warn('Invalid j-invariant format. Use integer or [a,b,c] notation');
    }
  };
  
  ui.onASelect = (aStr) => {
    let targetA = null;
    
    if (aStr.startsWith('[') && aStr.endsWith(']')) {
      const arr = aStr.slice(1, -1).split(',').map(s => parseInt(s.trim()));
      if (arr.every(n => !isNaN(n))) targetA = arr;
    } else {
      const aInt = parseInt(aStr);
      if (!isNaN(aInt)) targetA = intToBaseP(aInt, P);
    }
    
    if (targetA) {
      const curvesToSelect = [];
      selectedFieldGroup.traverse(obj => {
        if (obj.isMesh && obj.userData.curveID && jInvariantsEqual(obj.userData.a, targetA)) {
          curvesToSelect.push(obj);
        }
      });
      
      if (curvesToSelect.length > 0) {
        selectedCurves.forEach(curve => {
          setCurveVisualState(curve, 'base');
        });
        curvesToSelect.forEach(curve => setCurveVisualState(curve, 'selected'));
        selectedCurves = curvesToSelect;
        updateVolcanoEdgeHighlights();
      } else {
        console.warn('No curves found with A:', targetA);
      }
    }
  };
  
  ui.onBSelect = (bStr) => {
    let targetB = null;
    
    if (bStr.startsWith('[') && bStr.endsWith(']')) {
      const arr = bStr.slice(1, -1).split(',').map(s => parseInt(s.trim()));
      if (arr.every(n => !isNaN(n))) targetB = arr;
    } else {
      const bInt = parseInt(bStr);
      if (!isNaN(bInt)) targetB = intToBaseP(bInt, P);
    }
    
    if (targetB) {
      const curvesToSelect = [];
      selectedFieldGroup.traverse(obj => {
        if (obj.isMesh && obj.userData.curveID && jInvariantsEqual(obj.userData.b, targetB)) {
          curvesToSelect.push(obj);
        }
      });
      
      if (curvesToSelect.length > 0) {
        selectedCurves.forEach(curve => {
          setCurveVisualState(curve, 'base');
        });
        curvesToSelect.forEach(curve => setCurveVisualState(curve, 'selected'));
        selectedCurves = curvesToSelect;
        updateVolcanoEdgeHighlights();
      } else {
        console.warn('No curves found with B:', targetB);
      }
    }
  };
  
  ui.init();
  updateEllButtons();
  
  // Add layout toggle handler
  ui.onLayoutToggle = () => {
    layoutMode = layoutMode === 'circle' ? 'scatter' : 'circle';
    if (ui.layoutToggle) {
      ui.layoutToggle.textContent = `Layout: ${layoutMode === 'circle' ? 'Circle' : 'Scatter'}`;
    }
    console.log(`Layout mode changed to: ${layoutMode}`);
    
    // Redraw the right scene with the current selection
    if (selectedDiscriminant !== null) {
      drawNumberField(selectedDiscriminant);
    }
  };

  // Theme toggle
  const themeBtn = document.getElementById('theme_toggle');
  if (themeBtn) {
    themeBtn.addEventListener('click', () => {
      applyTheme(currentTheme === 'dark' ? 'light' : 'dark');
    });
  }

  const exportBtn = document.getElementById('export_png');
  if (exportBtn) {
    exportBtn.addEventListener('click', () => {
      exportHighDpiScene('active', 2, 'image/png');
    });
  }
  
  // Apply URL parameters for auto-selection
  if (URL_D !== null || URL_J !== null) {
    // Populate input fields
    if (URL_D !== null && ui.dInput) {
      ui.dInput.value = URL_D;
    }
    if (URL_J !== null && ui.jInput) {
      ui.jInput.value = URL_J;
    }
    
    // Auto-trigger selections after a short delay to ensure everything is loaded
    setTimeout(() => {
      if (URL_D !== null) {
        // First select by D
        ui.onDSelect(URL_D);
        
        // Check if j is valid for this D
        if (URL_J !== null) {
          setTimeout(() => {
            // Check if the selected field contains curves with this j
            const nf = fieldData.find(field => Number(field.discriminant) === URL_D);
            let jFoundInField = false;
            
            if (nf) {
              // Parse URL_J to target format
              let targetJ = null;
              if (URL_J.startsWith('[') && URL_J.endsWith(']')) {
                const arr = URL_J.slice(1, -1).split(',').map(s => parseInt(s.trim()));
                if (arr.every(n => !isNaN(n))) targetJ = arr;
              } else {
                const jInt = parseInt(URL_J);
                if (!isNaN(jInt)) targetJ = intToBaseP(jInt, P);
              }
              
              if (targetJ) {
                // Check if this j exists in the selected field
                for (const ic of nf.isogeny_classes) {
                  if (ic.orders) {
                    for (const order of ic.orders) {
                      if (order.curves) {
                        for (const curve of order.curves) {
                          if (jInvariantsEqual(curve.j, targetJ)) {
                            jFoundInField = true;
                            break;
                          }
                        }
                      }
                      if (jFoundInField) break;
                    }
                  }
                  if (jFoundInField) break;
                }
              }
              
              if (jFoundInField) {
                // j exists in this field, select it
                ui.onJSelect(URL_J);
              } else {
                // j doesn't exist in this field, clear j input
                console.warn(`j=${URL_J} not found in field D=${URL_D}, clearing j input`);
                if (ui.jInput) ui.jInput.value = '';
              }
            }
          }, 200);
        }
      } else if (URL_J !== null) {
        // Only j specified, no D
        ui.onJSelect(URL_J);
      }
    }, 300);
  }

  //ellipticRenderer = new EllipticCurveRenderer(P, N, q);
  //ellipticRenderer.setCurves(CURVES);

  sceneManager.init(P, N, q);
  sceneManager.left.renderer.domElement.addEventListener('click', onClick);
  sceneManager.left.renderer.domElement.addEventListener('mousemove', onMouseMove);
  sceneManager.left.renderer.domElement.addEventListener('mouseenter', () => { lastActiveSceneKey = 'left'; });
  sceneManager.right.renderer.domElement.addEventListener('click', onCurveClick);
  sceneManager.right.renderer.domElement.addEventListener('mousemove', onCurveMouseMove);
  sceneManager.right.renderer.domElement.addEventListener('mouseenter', () => { lastActiveSceneKey = 'right'; });

  // Export high-DPI frame from active pane: Ctrl/Cmd+Shift+E or Shift+P.
  window.addEventListener('keydown', (ev) => {
    const active = document.activeElement;
    const tag = (active?.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || active?.isContentEditable) return;

    const hotkeyA = (ev.ctrlKey || ev.metaKey) && ev.shiftKey && ev.code === 'KeyE';
    const hotkeyB = ev.shiftKey && !ev.ctrlKey && !ev.metaKey && ev.code === 'KeyP';
    if (!hotkeyA && !hotkeyB) return;

    ev.preventDefault();
    exportHighDpiScene('active', 2, 'image/png');
  });

  rebuildScene(ui.activeDiv);
  sceneManager.startAnimation(labelRenderer, rightLabelRenderer);

})().catch(err => {
  if (ui) {
    ui.showError('Failed to load curves.json.\n' + (err?.message || err));
  }
  console.error(err);
});

// Collect all ell values from loaded data and populate UI
// If discriminant is provided, only show ells from that specific field
function updateEllButtons(discriminant = null) {
  const ellSet = new Set();
  
  if (fieldData) {
    const fieldsToScan = discriminant !== null 
      ? fieldData.filter(nf => Number(nf.discriminant) === discriminant)
      : fieldData;
    
    fieldsToScan.forEach(nf => {
      nf.isogeny_classes.forEach(ic => {
        if (ic.volcanoes) {
          for (const volcano of ic.volcanoes) {
            ellSet.add(volcano.ell);
          }
        }
      });
    });
  }
  
  ui.populateEllSelect([...ellSet]);
}

// ==================== CURVE VISUALIZATION ====================

// Convert j-invariant array representation to integer in base p
// j = [a₀, a₁, a₂, ...] represents a₀ + a₁*p + a₂*p² + ...
function jInvariantToInteger(jArray, p) {
  if (!Array.isArray(jArray)) return jArray;
  
  let result = 0;
  let power = 1;
  for (const coeff of jArray) {
    result += coeff * power;
    power *= p;
  }
  return result;
}

// Convert integer to base-p array representation
function intToBaseP(n, p) {
  if (n === 0) return [0];
  const result = [];
  let remaining = n;
  while (remaining > 0) {
    result.push(remaining % p);
    remaining = Math.floor(remaining / p);
  }
  return result;
}

// Convert base-p array representation to integer
function baseP_to_int(arr, p) {
  if (!Array.isArray(arr)) return arr;
  let result = 0;
  for (let i = 0; i < arr.length; i++) {
    result += arr[i] * Math.pow(p, i);
  }
  return result;
}

// Convert j-invariant to a compact string representation
function formatJInvariant(jArray, p) {
  if (!Array.isArray(jArray)) return String(jArray);
  
  // Remove trailing zeros
  let trimmed = [...jArray];
  while (trimmed.length > 1 && trimmed[trimmed.length - 1] === 0) {
    trimmed.pop();
  }
  
  // Convert to integer in base p
  const intValue = jInvariantToInteger(trimmed, p);
  
  // Show both representations: integer value and coefficients
  if (trimmed.length === 1) {
    return `${intValue}`;
  } else {
    return `${intValue} [${trimmed.join(',')}]`;
  }
}

function _spreadDesiredXs(desiredXs, minGap, maxOffset) {
  if (!Array.isArray(desiredXs) || desiredXs.length === 0) return [];

  const indexed = desiredXs.map((x, i) => ({ x, i })).sort((a, b) => a.x - b.x);
  const placed = [];
  let prev = -Infinity;

  for (const item of indexed) {
    let x = item.x;
    if (x < prev + minGap) x = prev + minGap;
    placed.push({ i: item.i, x });
    prev = x;
  }

  const meanDesired = desiredXs.reduce((a, b) => a + b, 0) / desiredXs.length;
  const meanPlaced = placed.reduce((s, p) => s + p.x, 0) / placed.length;
  const shift = meanDesired - meanPlaced;

  const out = new Array(desiredXs.length).fill(0);
  for (const p of placed) {
    const shifted = p.x + shift;
    out[p.i] = Math.max(-maxOffset, Math.min(maxOffset, shifted));
  }
  return out;
}

function buildVolcanoGuidedCircleRanks(volcano) {
  const rankByCurve = new Map();
  if (!volcano || !Array.isArray(volcano.levels) || volcano.levels.length === 0) return rankByCurve;

  const normID = (v) => String(v);

  const levels = [...volcano.levels]
    .filter(l => Array.isArray(l.vertices) && l.vertices.length > 0)
    .sort((a, b) => b.h - a.h);
  if (levels.length === 0) return rankByCurve;

  const levelOf = new Map();
  const adjacency = new Map();

  const addAdj = (a, b) => {
    if (!adjacency.has(a)) adjacency.set(a, new Set());
    adjacency.get(a).add(b);
  };

  for (const level of levels) {
    for (const v of [...new Set(level.vertices)]) {
      levelOf.set(normID(v), level.h);
    }
    for (const [a, b] of (level.edges ?? [])) {
      addAdj(normID(a), normID(b));
      addAdj(normID(b), normID(a));
    }
  }

  const siblingStep = 0.55;   // compact spacing within a sibling cluster
  const clusterGap = 1.8;     // stronger separation between parent clusters

  const orderByLevel = new Map();
  const topVertices = [...new Set(levels[0].vertices.map(v => normID(v)))].sort((a, b) => String(a).localeCompare(String(b)));
  orderByLevel.set(levels[0].h, topVertices);
  topVertices.forEach((v, i) => rankByCurve.set(v, i * siblingStep));

  for (let li = 1; li < levels.length; li++) {
    const parentLevel = levels[li - 1];
    const currLevel = levels[li];
    const currVertices = [...new Set(currLevel.vertices.map(v => normID(v)))];
    if (currVertices.length === 0) continue;

    const orderedParents = (orderByLevel.get(parentLevel.h) ?? [...new Set(parentLevel.vertices.map(v => normID(v)))]);

    // Assign each child to exactly one primary parent (first in ordered parent sequence)
    const primaryParent = new Map();
    orderedParents.forEach((parent, pIdx) => {
      const children = [...(adjacency.get(parent) ?? [])]
        .filter(v => levelOf.get(v) === currLevel.h)
        .sort((a, b) => String(a).localeCompare(String(b)));
      for (const child of children) {
        if (!primaryParent.has(child)) {
          primaryParent.set(child, pIdx);
        }
      }
    });

    const blocks = new Map();
    orderedParents.forEach((_, i) => blocks.set(i, []));
    for (const child of currVertices) {
      if (!primaryParent.has(child)) continue;
      const pIdx = primaryParent.get(child);
      blocks.get(pIdx).push(child);
    }
    for (const [pIdx, arr] of blocks) {
      arr.sort((a, b) => String(a).localeCompare(String(b)));
      blocks.set(pIdx, arr);
    }

    const levelOrder = [];

    // Center each sibling cluster around its parent rank so fans are symmetric.
    // Keep sibling spacing compact and separate only between parent clusters.
    let prevPlaced = null;
    orderedParents.forEach((parent, pIdx) => {
      const kids = blocks.get(pIdx) ?? [];
      if (kids.length === 0) return;
      const desiredParentRank = Number.isFinite(rankByCurve.get(parent))
        ? rankByCurve.get(parent)
        : (pIdx * siblingStep);
      const halfWidth = ((kids.length - 1) * siblingStep) / 2;
      let parentRank = desiredParentRank;
      if (prevPlaced !== null) {
        const minCenter = prevPlaced.center + prevPlaced.halfWidth + clusterGap + halfWidth;
        if (parentRank < minCenter) parentRank = minCenter;
      }
      prevPlaced = { center: parentRank, halfWidth };

      const center = (kids.length - 1) / 2;
      kids.forEach((child, i) => {
        const childRank = parentRank + (i - center) * siblingStep;
        levelOrder.push(child);
        rankByCurve.set(child, childRank);
      });
    });

    // Any unattached vertices at this level go last.
    const inOrder = new Set(levelOrder);
    const leftovers = currVertices.filter(v => !inOrder.has(v)).sort((a, b) => String(a).localeCompare(String(b)));
    let cursor = (levelOrder.length > 0)
      ? ((rankByCurve.get(levelOrder[levelOrder.length - 1]) ?? 0) + clusterGap)
      : 0;
    for (const child of leftovers) {
      levelOrder.push(child);
      rankByCurve.set(child, cursor);
      cursor += siblingStep;
    }

    orderByLevel.set(currLevel.h, levelOrder);
  }

  return rankByCurve;
}

function buildVolcanoPrimaryParentMap(volcano, rankByCurve = null) {
  const primaryParentByCurve = new Map();
  if (!volcano || !Array.isArray(volcano.levels) || volcano.levels.length === 0) return primaryParentByCurve;

  const normID = (v) => String(v);

  const levels = [...volcano.levels]
    .filter(l => Array.isArray(l.vertices) && l.vertices.length > 0)
    .sort((a, b) => b.h - a.h);
  if (levels.length === 0) return primaryParentByCurve;

  const levelOf = new Map();
  const adjacency = new Map();
  const addAdj = (a, b) => {
    if (!adjacency.has(a)) adjacency.set(a, new Set());
    adjacency.get(a).add(b);
  };

  for (const level of levels) {
    for (const v of [...new Set(level.vertices)]) {
      levelOf.set(normID(v), level.h);
    }
    for (const [a, b] of (level.edges ?? [])) {
      addAdj(normID(a), normID(b));
      addAdj(normID(b), normID(a));
    }
  }

  const compareByRankThenID = (a, b) => {
    const ra = rankByCurve?.get(a);
    const rb = rankByCurve?.get(b);
    const aHas = Number.isFinite(ra);
    const bHas = Number.isFinite(rb);
    if (aHas && bHas && ra !== rb) return ra - rb;
    if (aHas && !bHas) return -1;
    if (!aHas && bHas) return 1;
    return String(a).localeCompare(String(b));
  };

  for (let li = 1; li < levels.length; li++) {
    const parentLevel = levels[li - 1];
    const currLevel = levels[li];
    const orderedParents = [...new Set(parentLevel.vertices.map(v => normID(v)))].sort(compareByRankThenID);

    for (const parent of orderedParents) {
      const children = [...(adjacency.get(parent) ?? [])]
        .filter(v => levelOf.get(v) === currLevel.h)
        .sort(compareByRankThenID);
      for (const child of children) {
        if (!primaryParentByCurve.has(child)) {
          primaryParentByCurve.set(child, parent);
        }
      }
    }
  }

  return primaryParentByCurve;
}

function getVolcanoConnectedComponents(volcano) {
  const componentByCurve = new Map();
  const components = [];
  if (!volcano || !Array.isArray(volcano.levels)) {
    return { componentByCurve, components };
  }

  const adjacency = new Map();
  const allVertices = new Set();

  const touch = (v) => {
    if (!adjacency.has(v)) adjacency.set(v, new Set());
    allVertices.add(v);
  };
  const link = (a, b) => {
    touch(a);
    touch(b);
    adjacency.get(a).add(b);
    adjacency.get(b).add(a);
  };

  volcano.levels.forEach(level => {
    (level.vertices ?? []).forEach(v => touch(v));
    (level.edges ?? []).forEach(([a, b]) => link(a, b));
  });

  const visited = new Set();
  for (const start of allVertices) {
    if (visited.has(start)) continue;
    const comp = [];
    const stack = [start];
    visited.add(start);
    while (stack.length > 0) {
      const v = stack.pop();
      comp.push(v);
      for (const u of (adjacency.get(v) ?? [])) {
        if (!visited.has(u)) {
          visited.add(u);
          stack.push(u);
        }
      }
    }
    components.push(comp.sort((a, b) => String(a).localeCompare(String(b))));
  }

  // Deterministic order: larger components first, then lexicographic.
  components.sort((a, b) => {
    if (b.length !== a.length) return b.length - a.length;
    return String(a[0] ?? '').localeCompare(String(b[0] ?? ''));
  });

  components.forEach((comp, idx) => {
    comp.forEach(v => componentByCurve.set(v, idx));
  });

  return { componentByCurve, components };
}

function computeComponentOffsets(componentRadii, boxSize) {
  const componentCount = Array.isArray(componentRadii) ? componentRadii.length : 0;
  if (componentCount <= 1) return [{ x: 0, z: 0 }];

  // User-requested behavior: split space into equal cells and place each
  // component at the center of a distinct cell.
  // (Do not squeeze by component footprint.)
  const cols = Math.ceil(Math.sqrt(componentCount));
  const rows = Math.ceil(componentCount / cols);
  const span = boxSize * 0.78;
  const stepX = cols > 1 ? span / (cols - 1) : 0;
  const stepZ = rows > 1 ? span / (rows - 1) : 0;
  const x0 = -((cols - 1) * stepX) / 2;
  const z0 = -((rows - 1) * stepZ) / 2;

  const out = new Array(componentCount).fill(null);
  for (let i = 0; i < componentCount; i++) {
    const row = Math.floor(i / cols);
    const col = i % cols;
    out[i] = { x: x0 + col * stepX, z: z0 + row * stepZ };
  }
  return out;
}

function drawCurvesForField(discriminant, traceFilter = null) {
  // Clear previous curve visualizations
  selectedFieldGroup.clear();

  const oldConductorLabels = sceneManager.right.scene.getObjectByName('right-conductor-rim-labels');
  if (oldConductorLabels) {
    oldConductorLabels.traverse((node) => {
      if (node instanceof CSS2DObject && node.element?.parentNode) {
        node.element.parentNode.removeChild(node.element);
      }
    });
    oldConductorLabels.removeFromParent();
  }
  
  // Clear and rebuild right scene axes
  removeCubeAxes(sceneManager.right.scene);
  
  // Find the number field with this discriminant
  const nf = fieldData.find(f => Number(f.discriminant) === discriminant);
  if (!nf) return;
  
  const boxBase = 36 * 0.9;
  const rightCurveWidth = boxBase * 1.5;
  const rightBoxWidth = boxBase * 2.0;
  const rightBoxHeight = boxBase;
  const rightBoxDepth = boxBase;
  const boxSize = rightCurveWidth; // keep curve radius formulas tied to 150% width

  const oldRightBounds = sceneManager.right.scene.getObjectByName('bounds-box');
  if (oldRightBounds) oldRightBounds.removeFromParent();
  const sphereRadius = 0.4;
  const circleRadius = 5.0; // Small radius for the circle around each trace position
  const traceCircleRadius = rightCurveWidth * 0.32; // Larger base radius for default order/twist view
  const maxCircleRadius = rightCurveWidth * 0.34; // Larger maximum default circle radius
  const minCircleRadius = 0.9; // Larger minimum default circle radius
  
  // Get selected ell filter
  const selectedEll = ui.selectedEll;
  const activeVolcanoTrace = traceFilter !== null ? Number(traceFilter) : selectedVolcanoTrace;
  const volcanoFocusMode = layoutMode === 'circle' && selectedEll !== 'ALL' && activeVolcanoTrace !== null;
  // Calculate conductor range for log scaling (same as left box)
  const conductors = [];
  nf.isogeny_classes.forEach(ic => {
    ic.orders.forEach(order => {
      conductors.push(order.conductor);
    });
  });
  const minCond = Math.min(...conductors);
  const maxCond = Math.max(...conductors);
  const minLogCond = Math.log(minCond + 1);
  const maxLogCond = Math.log(maxCond + 1);
  const logCondRange = maxLogCond - minLogCond || 1;
  
  // Get unique conductor values for axis ticks
  const uniqueConductors = [...new Set(conductors)].sort((a, b) => a - b);
  // Add conductor-based axes to right scene
  // Y axis maps conductor values (log-scaled, inverted)
  const conductorTicks = uniqueConductors.map(c => {
    const logCond = Math.log(c + 1);
    return -((logCond - minLogCond) / logCondRange * rightBoxHeight - rightBoxHeight/2);
  });
  
  // Create a map for exact lookup
  const yPosToCondMap = new Map();
  uniqueConductors.forEach((c, i) => {
    yPosToCondMap.set(conductorTicks[i], c);
  });
  
  const halfW = rightBoxWidth / 2;
  const halfH = rightBoxHeight / 2;
  const conductorLabelForY = (yVal) => {
    if (!Number.isFinite(Number(yVal)) || conductorTicks.length === 0) return '';
    let bestIdx = 0;
    let bestDist = Math.abs(conductorTicks[0] - yVal);
    for (let i = 1; i < conductorTicks.length; i++) {
      const d = Math.abs(conductorTicks[i] - yVal);
      if (d < bestDist) {
        bestDist = d;
        bestIdx = i;
      }
    }
    return `f = ${uniqueConductors[bestIdx]}`;
  };
  addCubeEdgeAxes(sceneManager.right.scene, (x, y, z) => new THREE.Vector3(x, y, z), {
    x: { domain: [-halfW, halfW], ticks: [], label: v => '' },
    y: { 
      domain: [-halfH, halfH], 
      ticks: conductorTicks,
      label: conductorLabelForY
    },
    z: { domain: [-halfW, halfW], ticks: [], label: v => '' },
    tickLen: 4.8, 
    labelOffset: 1.5,
    labelColor: THEMES[currentTheme].themeColor2,
    tickColor: THEMES[currentTheme].themeColor2,
    tickOutwardOnly: true,
    showEdges: false,
    color: THEMES[currentTheme].axis
  });

  const axisLines = sceneManager.right.scene.getObjectByName('axes-lines');
  const axisLabels = sceneManager.right.scene.getObjectByName('axes-labels');

  // Clear any previous right-side conductor grid group.
  drawConductorGrid(sceneManager.right.scene, {
    boxSize: rightBoxWidth,
    boxWidth: rightBoxWidth,
    boxHeight: rightBoxHeight,
    minLogCond,
    logCondRange,
    conductor: null,
    groupName: 'selected-conductor-grid-right',
    color: THEMES[currentTheme].themeColor1
  });
  
  // Calculate ranges for normalization
  const maxJIndex = Math.pow(P, N) - 1; // j-invariants range from 0 to p^n - 1
  
  // Use Hasse bound for trace: |t| ≤ 2√q, so t ∈ [-2√q, 2√q]
  const q = Math.pow(P, N);
  const hasseBound = 2 * Math.sqrt(q);
  const minTraceHasse = -hasseBound;
  const maxTraceHasse = hasseBound;
  const traceRange = maxTraceHasse - minTraceHasse; // This is 4√q
  
  // Circle layout parameters
  const totalIsogenyClasses = nf.isogeny_classes.length;
  const traceAngleStep = (2 * Math.PI) / Math.max(totalIsogenyClasses, 1);

  // Precompute trace centers used for radial placement of conductor ticks/labels.
  // Final placement is applied after curve dots are built, so we can use
  // per-conductor actual radii from the rendered levels.
  const traceCentersXZ = nf.isogeny_classes.map((_, traceIndex) => {
    if (volcanoFocusMode) {
      return new THREE.Vector2(0, 0);
    }
    if (layoutMode === 'scatter') {
      const trace = nf.isogeny_classes[traceIndex].trace;
      const z = traceRange > 0
        ? ((trace - minTraceHasse) / traceRange) * rightCurveWidth - rightCurveWidth / 2
        : 0;
      return new THREE.Vector2(0, z);
    }
    const traceAngle = traceIndex * traceAngleStep;
    return new THREE.Vector2(
      traceCircleRadius * Math.cos(traceAngle),
      traceCircleRadius * Math.sin(traceAngle)
    );
  });

  const rightAxisTickLen = 4.8;
  const rightAxisLabelExtraOutward = rightAxisTickLen * 1.0;
  const conductorOuterRadiusByValue = new Map();

  const placeAlongNearestPoleDirection = (obj, distanceFromPoleCenter = 0) => {
    if (!obj?.position || traceCentersXZ.length === 0) return;
    const current = new THREE.Vector2(Number(obj.position.x) || 0, Number(obj.position.z) || 0);

    let nearest = traceCentersXZ[0];
    let bestD2 = current.distanceToSquared(nearest);
    for (let i = 1; i < traceCentersXZ.length; i++) {
      const d2 = current.distanceToSquared(traceCentersXZ[i]);
      if (d2 < bestD2) {
        bestD2 = d2;
        nearest = traceCentersXZ[i];
      }
    }

    const fromCenter = current.clone().sub(nearest);
    if (fromCenter.lengthSq() < 1e-9) {
      fromCenter.set(-1, 0);
    }
    fromCenter.normalize().multiplyScalar(Math.max(0, Number(distanceFromPoleCenter) || 0));
    const target = nearest.clone().add(fromCenter);

    obj.position.x = target.x;
    obj.position.z = target.y;
  };

  // Make right-axis conductor label text smaller.
  axisLabels?.traverse?.((obj) => {
    if (!(obj instanceof CSS2DObject) || !obj.element?.style) return;
    const currentSize = parseFloat(obj.element.style.fontSize || '30');
    if (!Number.isFinite(currentSize) || currentSize <= 0) return;
    obj.element.style.fontSize = `${currentSize * 0.60}px`;
    obj.element.style.lineHeight = '1';
    obj.element.style.whiteSpace = 'nowrap';
    obj.center.set(0.5, 0.5);
    obj.element.querySelector('.axis-label-center-dot')?.remove();
  });
  
  console.log(`Normalization ranges: P=${P}, N=${N}, q=${q}, maxJIndex=${maxJIndex}, Hasse bounds: t ∈ [${minTraceHasse.toFixed(2)}, ${maxTraceHasse.toFixed(2)}], boxSize=${boxSize}, layoutMode=${layoutMode}`);
  
  // Draw trace markers for both modes
  const traces = nf.isogeny_classes.map(ic => ic.trace);
  const uniqueTraces = [...new Set(traces)].sort((a, b) => a - b);
  
  nf.isogeny_classes.forEach((ic, traceIndex) => {
    const trace = ic.trace;

    const isActiveVolcanoTrace = volcanoFocusMode && Number(trace) === Number(activeVolcanoTrace);
    if (volcanoFocusMode && !isActiveVolcanoTrace) {
      return;
    }
    
    let xPos, zPos;
    
    if (volcanoFocusMode) {
      // In volcano focus mode, pin the active trace indicator at center.
      xPos = 0;
      zPos = 0;
    } else if (layoutMode === 'scatter') {
      // In scatter mode, place markers along the Z axis
      xPos = 0;
      zPos = traceRange > 0 ? ((trace - minTraceHasse) / traceRange) * rightCurveWidth - rightCurveWidth / 2 : 0;
    } else {
      // In circle mode, place markers at the radial trace positions
      const traceAngle = traceIndex * traceAngleStep;
      xPos = traceCircleRadius * Math.cos(traceAngle);
      zPos = traceCircleRadius * Math.sin(traceAngle);
    }
    
    // Match D_K pole vertical span exactly:
    // bottom at box floor, top at (box top + 4.8).
    const dkExtraTop = 4.8;
    const poleBottomY = -rightBoxHeight / 2;
    const poleTopY = rightBoxHeight / 2 + dkExtraTop;
    const poleLen = Math.max(0.001, poleTopY - poleBottomY);
    const poleMaterial = new THREE.MeshBasicMaterial({
      color: THEMES[currentTheme].sphere.hover.getHex(),
      transparent: true,
      opacity: 0.9,
      depthWrite: false
    });
    const pole = new THREE.Mesh(
      new THREE.CylinderGeometry(GUIDE_RADIUS, GUIDE_RADIUS, poleLen, 16),
      poleMaterial
    );
    pole.position.set(xPos, (poleBottomY + poleTopY) / 2, zPos);
    pole.renderOrder = 8;
    selectedFieldGroup.add(pole);
    
    // Add label at the top
    const labelDiv = document.createElement('div');
    labelDiv.style.color = THEMES[currentTheme].themeColor2;
    labelDiv.style.fontSize = '30px';
    labelDiv.style.fontFamily = 'ui-monospace, monospace';
    labelDiv.style.pointerEvents = 'none';
    labelDiv.textContent = `t=${trace}`;
    
    const label = new CSS2DObject(labelDiv);
    // Match D_K label offset above pole top.
    label.position.set(xPos, poleTopY + 4.5, zPos);
    selectedFieldGroup.add(label);
  });
  
  // Map to store sphere meshes by curve ID for edge drawing
  const curveIdToMesh = new Map();
  
  // Loop over all isogeny classes (traces) in this field
  nf.isogeny_classes.forEach((ic, traceIndex) => {
    const trace = ic.trace;
    
    // Skip if trace filter is set and doesn't match
    if (traceFilter !== null && trace !== traceFilter) {
      return;
    }
    
    // Build a set of curve IDs that belong to selected ell volcano(es)
    // Also build a map of curveID -> height for highlighting
    let allowedCurveIDs = null;
    let selectedVolcano = null;
    const curveHeights = new Map(); // curveID -> h value
    
    if (selectedEll !== 'ALL' && ic.volcanoes) {
      allowedCurveIDs = new Set();
      for (const volcano of ic.volcanoes) {
        if (volcano.ell === Number(selectedEll)) {
          if (activeVolcanoTrace !== null && Number(ic.trace) === Number(activeVolcanoTrace)) {
            selectedVolcano = volcano;
          }
          volcano.levels.forEach(level => {
            level.vertices.forEach(curveID => {
              allowedCurveIDs.add(curveID);
              curveHeights.set(curveID, level.h);
            });
          });
        }
      }
    }

    // Keep default layout unless a volcano row was explicitly clicked.
    // In volcano mode, align lower-level vertices under connected upper-level vertices.
    const volcanoCircleRanks = (selectedVolcano && selectedEll !== 'ALL')
      ? buildVolcanoGuidedCircleRanks(selectedVolcano)
      : null;
    const volcanoPrimaryParentByCurve = (selectedVolcano && selectedEll !== 'ALL')
      ? buildVolcanoPrimaryParentMap(selectedVolcano, volcanoCircleRanks)
      : null;
    const componentData = (selectedVolcano && selectedEll !== 'ALL')
      ? getVolcanoConnectedComponents(selectedVolcano)
      : null;

    // Count curves per (component, conductor level).
    const componentConductorCounts = new Map();
    if (componentData) {
      ic.orders.forEach(order => {
        const c = Number(order.conductor);
        if (!Number.isFinite(c)) return;
        (order.curves ?? []).forEach(curve => {
          const curveID = curve.ID;
          if (allowedCurveIDs !== null && !allowedCurveIDs.has(curveID)) return;
          const compIdx = DISTRIBUTE_COMPONENTS
            ? (componentData.componentByCurve.get(curveID) ?? 0)
            : 0;
          const key = `${compIdx}|${c}`;
          componentConductorCounts.set(key, (componentConductorCounts.get(key) || 0) + 1);
        });
      });
    }

    let componentOffsets = [{ x: 0, z: 0 }];

    const compCountGlobal = Math.max(1, componentData?.components?.length ?? 1);
    const desiredSpacingPre = volcanoFocusMode
      ? (DISTRIBUTE_COMPONENTS ? (1.35 * sphereRadius) : (1.2 * sphereRadius))
      : (3.6 * sphereRadius);
    const focusRadiusScalePre = (0.26 / Math.sqrt(compCountGlobal)) + 0.05;
    const maxRadiusHerePre = volcanoFocusMode
      ? (boxSize * (DISTRIBUTE_COMPONENTS
        ? Math.max(0.16, Math.min(0.49, focusRadiusScalePre))
        : 0.40))
      : maxCircleRadius;
    const minRadiusHerePre = volcanoFocusMode
      ? (DISTRIBUTE_COMPONENTS
        ? Math.max(0.25, 0.55 / Math.sqrt(compCountGlobal))
        : Math.max(0.6, sphereRadius * 3.0))
      : minCircleRadius;

    if (componentData && DISTRIBUTE_COMPONENTS) {
      const componentMaxRadii = new Array(compCountGlobal).fill(minRadiusHerePre);
      const levelGap = Math.max(0.12, sphereRadius * 0.6);
      const growthFactor = 1.10;

      for (let compIdx = 0; compIdx < compCountGlobal; compIdx++) {
        const conductors = [...new Set(
          [...componentConductorCounts.keys()]
            .map(k => String(k).split('|'))
            .filter(([c]) => Number(c) === compIdx)
            .map(([, c]) => Number(c))
            .filter(Number.isFinite)
        )].sort((a, b) => a - b);

        let prevR = minRadiusHerePre / growthFactor;
        for (const c of conductors) {
          const count = Math.max(1, Number(componentConductorCounts.get(`${compIdx}|${c}`)) || 1);
          const idealR = (count * desiredSpacingPre) / (2 * Math.PI);
          const required = Math.max(prevR + levelGap, prevR * growthFactor);
          const r = Math.min(maxRadiusHerePre, Math.max(minRadiusHerePre, idealR, required));
          componentMaxRadii[compIdx] = Math.max(componentMaxRadii[compIdx], r);
          prevR = r;
        }
      }

      componentOffsets = computeComponentOffsets(componentMaxRadii, boxSize);
    }

    // Also enforce monotone widening by displayed conductor row (top -> bottom).
    const conductorBaseRadius = new Map(); // key: conductor -> enforced radius
    const conductorDesiredSpacing = new Map(); // key: conductor -> target spacing on that row
    const conductorGapByLevel = new Map(); // key: conductor -> additive radial gap vs previous row
    if (volcanoFocusMode && selectedVolcano) {
      const minDotPadding = DISTRIBUTE_COMPONENTS
        ? Math.max(0.14, 1.25 * sphereRadius)
        : Math.max(0.22, 1.80 * sphereRadius);
      const desiredSpacingByConductor = DISTRIBUTE_COMPONENTS ? (1.15 * sphereRadius) : (0.9 * sphereRadius);
      const focusRadiusScaleByConductor = (0.26 / Math.sqrt(compCountGlobal)) + 0.05;
      const maxRadiusByConductor = boxSize * (DISTRIBUTE_COMPONENTS
        ? Math.max(0.16, Math.min(0.49, focusRadiusScaleByConductor))
        : 0.52);
      const minRadiusByConductor = DISTRIBUTE_COMPONENTS
        ? Math.max(0.25, 0.55 / Math.sqrt(compCountGlobal))
        : Math.max(0.6, sphereRadius * 3.0);
      const conductorGap = DISTRIBUTE_COMPONENTS ? Math.max(0.12, sphereRadius * 0.6) : Math.max(0.18, sphereRadius * 0.9);
      const conductorGrowthFactor = 1.10;

      const curvesCountByConductor = new Map();
      ic.orders.forEach(order => {
        const c = Number(order.conductor);
        if (!Number.isFinite(c)) return;
        const count = (allowedCurveIDs !== null)
          ? (order.curves ?? []).filter(curve => allowedCurveIDs.has(curve.ID)).length
          : (order.curves ?? []).length;
        curvesCountByConductor.set(c, count);
      });

      const sortedConductors = [...curvesCountByConductor.keys()].sort((a, b) => a - b);
      if (!DISTRIBUTE_COMPONENTS && sortedConductors.length > 0) {
        const topTargetRadius = boxSize * 0.25;
        const bottomTargetRadius = boxSize * 0.50;
        const topPadding = Math.max(minDotPadding, 2.20 * sphereRadius);
        const bottomPadding = Math.max(minDotPadding * 0.85, 1.55 * sphereRadius);
        for (let idx = 0; idx < sortedConductors.length; idx++) {
          const c = sortedConductors[idx];
          const t = (sortedConductors.length <= 1)
            ? 0
            : (idx / (sortedConductors.length - 1));
          const cnt = Math.max(1, curvesCountByConductor.get(c) || 1);
          const targetRadius = topTargetRadius + (bottomTargetRadius - topTargetRadius) * t;
          const spacingFromTargetRadius = (2 * Math.PI * targetRadius) / cnt;
          const rowMinPadding = topPadding + (bottomPadding - topPadding) * t;
          const rowSpacing = Math.max(rowMinPadding, spacingFromTargetRadius);
          conductorDesiredSpacing.set(c, rowSpacing);
          conductorGapByLevel.set(c, Math.max(0.08, rowSpacing * 0.22));
        }
      }
      let prevR = minRadiusByConductor / conductorGrowthFactor;
      for (let idx = 0; idx < sortedConductors.length; idx++) {
        const c = sortedConductors[idx];
        const cnt = Math.max(1, curvesCountByConductor.get(c) || 1);
        const rowSpacing = conductorDesiredSpacing.get(c) || desiredSpacingByConductor;
        const rowGap = conductorGapByLevel.get(c) || conductorGap;
        const idealR = (cnt * rowSpacing) / (2 * Math.PI);
        const stepsLeft = sortedConductors.length - idx;
        const maxPossibleGrowth = (prevR > 0 && stepsLeft > 0)
          ? Math.pow(Math.max(1.0, maxRadiusByConductor / prevR), 1 / stepsLeft)
          : conductorGrowthFactor;
        const adaptiveGrowth = Math.max(1.001, Math.min(conductorGrowthFactor, maxPossibleGrowth));
        const requiredByConductor = Math.max(prevR + rowGap, prevR * adaptiveGrowth);
        const r = Math.min(maxRadiusByConductor, Math.max(minRadiusByConductor, idealR, requiredByConductor));
        conductorBaseRadius.set(c, r);
        prevR = r;
      }
    }
    
    // Loop over all orders (conductors) in this isogeny class
    ic.orders.forEach(order => {
      const conductor = order.conductor;
      // Y position based on conductor (inverted: higher conductor = lower Y)
      // Apply log scale for consistent spacing with left box
      const logCond = Math.log(conductor + 1);
      const yPos = -((logCond - minLogCond) / logCondRange * rightBoxHeight - rightBoxHeight/2);

      const curvesToDraw = (allowedCurveIDs !== null)
        ? order.curves.filter(c => allowedCurveIDs.has(c.ID))
        : order.curves;

      let orderedCurves = curvesToDraw;
      if (layoutMode === 'circle' && volcanoCircleRanks && curvesToDraw.length > 1) {
        const withMeta = curvesToDraw.map(curve => {
          const curveID = curve.ID;
          const curveKey = String(curveID);
          const match = curveID.match(/_(\d+)$/);
          const localIndex = match ? parseInt(match[1]) : 0;
          const compIdx = componentData?.componentByCurve?.get(curveID) ?? 0;
          const guidedRank = volcanoCircleRanks.has(curveKey)
            ? volcanoCircleRanks.get(curveKey)
            : Number.POSITIVE_INFINITY;
          return { curve, localIndex, guidedRank, compIdx };
        });
        withMeta.sort((a, b) => {
          if (DISTRIBUTE_COMPONENTS && a.compIdx !== b.compIdx) return a.compIdx - b.compIdx;
          if (a.guidedRank !== b.guidedRank) return a.guidedRank - b.guidedRank;
          return a.localIndex - b.localIndex;
        });
        orderedCurves = withMeta.map(x => x.curve);
      }

      // Per-component circular slots/counts at this conductor for split display.
      const compToCurves = new Map();
      if (DISTRIBUTE_COMPONENTS) {
        for (const c of orderedCurves) {
          const compIdx = componentData?.componentByCurve?.get(c.ID) ?? 0;
          if (!compToCurves.has(compIdx)) compToCurves.set(compIdx, []);
          compToCurves.get(compIdx).push(c);
        }
      }
      const compSlotByID = new Map();
      const compCountByID = new Map();
      for (const [compIdx, arr] of compToCurves) {
        arr.forEach((c, i) => {
          compSlotByID.set(c.ID, i);
          compCountByID.set(c.ID, arr.length);
        });
      }

      // Keep radial placement simple: order only, no extra angle remapping.
      const guidedAngleByID = new Map();
      
      // Draw a dot for each curve in this order
      orderedCurves.forEach((curve, circleSlot) => {
        const curveID = curve.ID;
        const height = curveHeights.get(curveID) || 0;
        const compIdx = componentData?.componentByCurve?.get(curveID) ?? 0;
        const compOffset = componentOffsets[compIdx] ?? { x: 0, z: 0 };
        const compSlot = DISTRIBUTE_COMPONENTS ? (compSlotByID.get(curveID) ?? circleSlot) : circleSlot;
        const compCount = DISTRIBUTE_COMPONENTS ? (compCountByID.get(curveID) ?? orderedCurves.length) : orderedCurves.length;
        
        // Parse curve ID to extract local index (e.g., "f1_0" -> 0, "f2_1" -> 1)
        const match = curveID.match(/_(\d+)$/);
        const localIndex = match ? parseInt(match[1]) : 0;
        
        let xPos, zPos;
        
        if (layoutMode === 'scatter') {
          // SCATTER MODE: map j-invariant to X, trace to Z
          const jIndex = baseP_to_int(curve.j, P);
          
          // Map j-index to X position
          if (maxJIndex > 0) {
            xPos = (jIndex / maxJIndex) * rightCurveWidth - rightCurveWidth / 2;
          } else {
            xPos = 0;
          }
          
          // Map trace to Z position using Hasse bound
          if (traceRange > 0) {
            zPos = ((trace - minTraceHasse) / traceRange) * rightCurveWidth - rightCurveWidth / 2;
          } else {
            zPos = 0;
          }

          if (componentData && DISTRIBUTE_COMPONENTS) {
            xPos += compOffset.x;
            zPos += compOffset.z;
          }
        } else {
          // CIRCLE MODE: radial arrangement by trace, circles around each trace
          const traceAngle = traceIndex * traceAngleStep;
          const baseXPos = (volcanoFocusMode ? 0 : traceCircleRadius * Math.cos(traceAngle))
            + (DISTRIBUTE_COMPONENTS ? compOffset.x : 0);
          const baseZPos = (volcanoFocusMode ? 0 : traceCircleRadius * Math.sin(traceAngle))
            + (DISTRIBUTE_COMPONENTS ? compOffset.z : 0);
          
          if (order.class_number === 1) {
            xPos = baseXPos;
            zPos = baseZPos;
          } else {
            const conductorKey = `${DISTRIBUTE_COMPONENTS ? compIdx : 0}|${Number(conductor)}`;
            const conductorCount = componentConductorCounts.get(conductorKey) || 0;
            const ringPopulationForRadius = Math.max(1, conductorCount || compCount);
            const totalCurvesAtConductor = Math.max(1, compCount);
            // Calculate circle radius so spheres are evenly spaced
            // For n spheres on a circle, spacing = R * 2π/n
            // In non-separated mode we need a much wider ring footprint.
            const conductorNumber = Number(conductor);
            const defaultDesiredSpacing = volcanoFocusMode
              ? (DISTRIBUTE_COMPONENTS
                ? Math.max(0.14, 1.25 * sphereRadius)
                : Math.max(0.20, 1.65 * sphereRadius))
              : (3.6 * sphereRadius);
            const desiredSpacing = (volcanoFocusMode && conductorDesiredSpacing.has(conductorNumber))
              ? (conductorDesiredSpacing.get(conductorNumber) || defaultDesiredSpacing)
              : defaultDesiredSpacing;
            const idealRadius = (ringPopulationForRadius * desiredSpacing) / (2 * Math.PI);
            const compCountGlobal = componentData?.components?.length ?? 1;
            const focusRadiusScale = (0.26 / Math.sqrt(compCountGlobal)) + 0.05;
            const maxRadiusHere = volcanoFocusMode
              ? (boxSize * (DISTRIBUTE_COMPONENTS
                ? Math.max(0.16, Math.min(0.49, focusRadiusScale))
                : 0.52))
              : maxCircleRadius;
            const minRadiusHere = volcanoFocusMode
              ? (DISTRIBUTE_COMPONENTS
                ? Math.max(0.25, 0.55 / Math.sqrt(compCountGlobal))
                : Math.max(0.6, sphereRadius * 3.0))
              : minCircleRadius;
            let dynamicCircleRadius = Math.max(minRadiusHere, Math.min(maxRadiusHere, idealRadius));

            // Enforce monotone widening by conductor row as well.
            const enforcedConductorRadius = conductorBaseRadius.get(Number(conductor));
            if (Number.isFinite(enforcedConductorRadius)) {
              dynamicCircleRadius = Math.max(dynamicCircleRadius, enforcedConductorRadius);
            }
            
            const angleStep = (2 * Math.PI) / Math.max(totalCurvesAtConductor, 1);
            const angle = Number.isFinite(guidedAngleByID.get(String(curveID)))
              ? guidedAngleByID.get(String(curveID))
              : (compSlot * angleStep);
            
            xPos = baseXPos + dynamicCircleRadius * Math.cos(angle);
            zPos = baseZPos + dynamicCircleRadius * Math.sin(angle);
          }
        }
        
        const geometry = new THREE.SphereGeometry(sphereRadius, 16, 16);
        
        // Check if this curve is at height h > 0
        const isAboveFloor = height > 0;
        
        // Use different color/material for curves above floor
        let sphereColor, sphereOpacity;
        const sp2 = THEMES[currentTheme].sphere;
        if (isAboveFloor) {
          sphereColor = sp2.aboveFloor.clone();
          sphereOpacity = 1.0;
        } else {
          sphereColor = sp2.floorLevel.clone();
          sphereOpacity = 0.8;
        }
        
        const material = new THREE.MeshPhongMaterial({
          color: sphereColor,
          transparent: true,
          opacity: sphereOpacity
        });
        
        const sphere = new THREE.Mesh(geometry, material);
        sphere.position.set(xPos, yPos, zPos);
        sphere.userData = {
          curveID: curve.ID,
          j: curve.j,
          jFormatted: formatJInvariant(curve.j, P),
          a: curve.A,
          b: curve.B,
          conductor: conductor,
          discriminant: discriminant,
          trace: ic.trace,
          localIndex: localIndex,
          height: height  // Store height for reference
        };
        
        selectedFieldGroup.add(sphere);

        // Track the outer radius reached by this conductor level,
        // measured from its trace-center in the XZ plane.
        const centerXZ = traceCentersXZ[traceIndex] || new THREE.Vector2(0, 0);
        const radial = Math.hypot(xPos - centerXZ.x, zPos - centerXZ.y);
        const prevR = conductorOuterRadiusByValue.get(Number(conductor)) ?? 0;
        if (radial > prevR) {
          conductorOuterRadiusByValue.set(Number(conductor), radial);
        }

        // Store sphere by curve ID for edge drawing
        curveIdToMesh.set(curve.ID, sphere);
      });
    });
    
    // Draw edges for selected volcano
    if (selectedVolcano && selectedEll !== 'ALL') {
      // First pass: count duplicate edges
      const edgeCounts = new Map();
      const edgeIndices = new Map();
      
      selectedVolcano.levels.forEach(level => {
        level.edges.forEach(edge => {
          const [fromID, toID] = edge;
          const edgeKey = `${fromID}->${toID}`;
          edgeCounts.set(edgeKey, (edgeCounts.get(edgeKey) || 0) + 1);
        });
      });
      
      // Second pass: draw edges with offsets for duplicates
      selectedVolcano.levels.forEach(level => {
        level.edges.forEach(edge => {
          const [fromID, toID] = edge;
          const edgeKey = `${fromID}->${toID}`;
          const totalCount = edgeCounts.get(edgeKey);
          const currentIndex = edgeIndices.get(edgeKey) || 0;
          edgeIndices.set(edgeKey, currentIndex + 1);
          
          const fromMesh = curveIdToMesh.get(fromID);
          const toMesh = curveIdToMesh.get(toID);
          
          if (fromMesh && toMesh) {
            const baseEdgeColor = 0x444444;
            const baseEdgeOpacity = 0.4;
            const lineMaterial = new THREE.LineBasicMaterial({
              color: baseEdgeColor,
              linewidth: 1,
              opacity: baseEdgeOpacity,
              transparent: true
            });
            
            let lineGeometry;
            
            // Check if it's a self-loop
            if (fromID === toID) {
              // Create a circular loop coming out from the vertex
              const center = fromMesh.position.clone();
              const loopRadius = 0.5; // Much larger loop
              const segments = 64;
              const points = [];
              
              // Calculate rotation angle for duplicate self-loops
              const rotationAngle = totalCount > 1 ? (currentIndex / totalCount) * Math.PI * 2 : 0;
              
              // Create a circle extending outward from the vertex
              // Offset the center of the circle away from the vertex
              const offsetCenter = center.clone();
              const offsetDir = new THREE.Vector3(
                Math.cos(rotationAngle),
                Math.sin(rotationAngle),
                0.5
              ).normalize();
              offsetCenter.add(offsetDir.multiplyScalar(loopRadius));
              
              for (let i = 0; i <= segments; i++) {
                const angle = (i / segments) * Math.PI * 2;
                // Rotate the loop plane based on rotationAngle
                const localX = loopRadius * Math.cos(angle);
                const localY = loopRadius * Math.sin(angle);
                const x = offsetCenter.x + localX * Math.cos(rotationAngle) - localY * Math.sin(rotationAngle);
                const y = offsetCenter.y + localX * Math.sin(rotationAngle) + localY * Math.cos(rotationAngle);
                const z = offsetCenter.z;
                points.push(new THREE.Vector3(x, y, z));
              }
              
              lineGeometry = new THREE.BufferGeometry().setFromPoints(points);
              
              // Keep same base styling as regular edges; highlight is handled dynamically.
              lineMaterial.color.setHex(baseEdgeColor);
              lineMaterial.opacity = baseEdgeOpacity;
              
              // Add arrow for self-loop (at top of circle)
              // DISABLED: Arrow cones
              // const arrowPos = new THREE.Vector3(offsetCenter.x, offsetCenter.y, offsetCenter.z + loopRadius);
              // const arrowDir = new THREE.Vector3(1, 0, 0); // Tangent direction
              // const arrowLength = 0.5;
              // const arrowColor = 0xff00ff;
              // const arrowHelper = new THREE.ArrowHelper(arrowDir, arrowPos, arrowLength, arrowColor, 0.4, 0.3);
              // selectedFieldGroup.add(arrowHelper);
            } else {
              // Regular edge between two different vertices
              const direction = new THREE.Vector3().subVectors(toMesh.position, fromMesh.position);
              const edgeLength = direction.length();
              direction.normalize();
              
              // Calculate perpendicular offset for duplicate edges
              let offsetVector = new THREE.Vector3(0, 0, 0);
              if (totalCount > 1) {
                // Create perpendicular vector for offset
                const perpendicular = new THREE.Vector3(-direction.y, direction.x, 0).normalize();
                if (perpendicular.length() < 0.01) {
                  // If edge is vertical, use a different perpendicular
                  perpendicular.set(-direction.z, 0, direction.x).normalize();
                }
                // Spread edges symmetrically around the center
                const offsetAmount = ((currentIndex - (totalCount - 1) / 2) * 0.15);
                offsetVector = perpendicular.multiplyScalar(offsetAmount);
              }
              
              const points = [
                fromMesh.position.clone().add(offsetVector),
                toMesh.position.clone().add(offsetVector)
              ];
              lineGeometry = new THREE.BufferGeometry().setFromPoints(points);
              
              // Add arrow head at destination
              // DISABLED: Arrow cones
              // Position arrow so it stops at the sphere surface (radius = 0.5)
              // const arrowPos = toMesh.position.clone().add(offsetVector).sub(direction.clone().multiplyScalar(sphereRadius + 0.3));
              // const arrowLength = 0.5;
              // const arrowColor = 0xffffff;
              // const arrowHelper = new THREE.ArrowHelper(direction, arrowPos, arrowLength, arrowColor, 0.4, 0.3);
              // selectedFieldGroup.add(arrowHelper);
            }
            
            const line = new THREE.Line(lineGeometry, lineMaterial);
            line.userData = {
              edgeFrom: fromID,
              edgeTo: toID,
              ell: selectedVolcano.ell,
              baseEdgeColor,
              baseEdgeOpacity
            };
            selectedFieldGroup.add(line);
          }
        });
      });
    }
  });

  // Position ticks/labels per conductor level at 10% outside that level radius.
  const nearestConductorForY = (yVal) => {
    if (!Number.isFinite(Number(yVal)) || conductorTicks.length === 0) return null;
    let bestIdx = 0;
    let bestDist = Math.abs(conductorTicks[0] - yVal);
    for (let i = 1; i < conductorTicks.length; i++) {
      const d = Math.abs(conductorTicks[i] - yVal);
      if (d < bestDist) {
        bestDist = d;
        bestIdx = i;
      }
    }
    return uniqueConductors[bestIdx];
  };

  const nearestCenterAndDirXZ = (obj) => {
    const current = new THREE.Vector2(Number(obj?.position?.x) || 0, Number(obj?.position?.z) || 0);
    let nearest = traceCentersXZ[0] || new THREE.Vector2(0, 0);
    let bestD2 = current.distanceToSquared(nearest);
    for (let i = 1; i < traceCentersXZ.length; i++) {
      const d2 = current.distanceToSquared(traceCentersXZ[i]);
      if (d2 < bestD2) {
        bestD2 = d2;
        nearest = traceCentersXZ[i];
      }
    }
    const dir2 = current.clone().sub(nearest);
    if (dir2.lengthSq() < 1e-9) dir2.set(-1, 0);
    dir2.normalize();
    return { nearest, dir2 };
  };

  const tickPlacements = [];
  axisLines?.children?.forEach((obj) => {
    const conductorForObj = nearestConductorForY(obj?.position?.y);
    const levelRadius = conductorOuterRadiusByValue.get(Number(conductorForObj)) ?? 0;
    const baseOutsideRadius = Math.max(sphereRadius, levelRadius) * 1.10; // 10% outside level radius
    const lineExtraTowardLabel = rightAxisLabelExtraOutward;
    const outsideRadius = baseOutsideRadius + lineExtraTowardLabel;
    const { nearest, dir2 } = nearestCenterAndDirXZ(obj);

    // Tick must start at pole center and end at outsideRadius.
    const baseLen = (obj?.geometry?.parameters?.height ?? rightAxisTickLen);
    if (Number.isFinite(baseLen) && baseLen > 1e-9) {
      const currentScaleX = Number(obj?.scale?.x) || 1;
      const currentScaleZ = Number(obj?.scale?.z) || 1;
      obj.scale.set(currentScaleX, outsideRadius / baseLen, currentScaleZ);
    }

    const dir3 = new THREE.Vector3(dir2.x, 0, dir2.y);
    obj.position.set(
      nearest.x + dir2.x * (outsideRadius * 0.5),
      Number(obj?.position?.y) || 0,
      nearest.y + dir2.y * (outsideRadius * 0.5)
    );
    obj.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir3.normalize());

    tickPlacements.push({
      conductor: Number(conductorForObj),
      nearest: nearest.clone(),
      dir2: dir2.clone(),
      outsideRadius,
      y: Number(obj?.position?.y) || 0
    });
  });

  const rightAxisLabelMargin = 2.4;
  axisLabels?.children?.forEach((obj) => {
    if (obj instanceof CSS2DObject) {
      // Use CSS2D anchor at the exact text-box center: translate(-50%, -50%).
      obj.center.set(0.5, 0.5);
    }
    const labelY = Number(obj?.position?.y) || 0;
    let tickPlacement = null;
    let bestDist = Infinity;
    for (const tp of tickPlacements) {
      const d = Math.abs((Number(tp?.y) || 0) - labelY);
      if (d < bestDist) {
        bestDist = d;
        tickPlacement = tp;
      }
    }

    const conductorForObj = nearestConductorForY(labelY);
    const levelRadius = conductorOuterRadiusByValue.get(Number(conductorForObj)) ?? 0;
    const baseOutsideRadius = Math.max(sphereRadius, levelRadius) * 1.10;
    const lineExtraTowardLabel = rightAxisLabelExtraOutward;
    const outsideRadius = baseOutsideRadius + lineExtraTowardLabel;
    const nearest = tickPlacement?.nearest ?? nearestCenterAndDirXZ(obj).nearest;
    const dir2 = tickPlacement?.dir2 ?? nearestCenterAndDirXZ(obj).dir2;
    const labelDistance = (tickPlacement?.outsideRadius ?? outsideRadius) + rightAxisLabelMargin;
    const yVal = Number.isFinite(tickPlacement?.y) ? tickPlacement.y : (Number(obj?.position?.y) || 0);
    obj.position.set(
      nearest.x + dir2.x * labelDistance,
      yVal,
      nearest.y + dir2.y * labelDistance
    );
  });

  updateVolcanoEdgeHighlights();
}


