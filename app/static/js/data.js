// Data structures for elliptic curve catalogue



export class Curve {
  constructor(data) {
    this.ID = data.ID;
    this.j = data.j;  // j-invariant coefficients
    this.j_minpoly = data.j_minpoly ?? data.jMinpoly ?? null;
    this.A = data.A;  // Weierstrass A coefficients
    this.B = data.B;  // Weierstrass B coefficients
    this.inv = data.inv;  // Group structure invariants
    this.torsion_subgroups = data.torsion_subgroups ?? {};
  }
}

export class EndomorphismOrder {
  constructor(data) {
    this.conductor = data.f ?? data.conductor;
    this.class_number = data.cn ?? data.class_number;
    this.curves = (data.curves ?? []).map(c => new Curve(c));
  }

  getCurve(id) {
    return this.curves.find(c => c.ID === id);
  }

  addCurves(curvesArray) {
    this.curves = curvesArray.map(c => new Curve(c));
  }
}

export class VolcanoLevel {
  constructor(data) {
    this.h = data.h;  // Height/level in volcano
    this.vertices = data.vertrices ?? data.vertices ?? [];  // Curve IDs (note: JSON has typo "vertrices")
    this.edges = data.edges ?? [];  // [from_id, to_id] pairs
  }
}

export class Volcano {
  constructor(data) {
    this.ell = data.ell;  // Prime ℓ
    // Backward/forward compatible roots:
    // - old shape: fx_roots
    // - new shape: +/- roots with active sign chosen by trace variant
    this.fx_roots = data.fx_roots
      ?? (data._activeSign === '+' ? data['+fx_roots'] : data['-fx_roots'])
      ?? [];
    this.levels = (data.levels ?? []).map(l => new VolcanoLevel(l));
  }

  getLevel(h) {
    return this.levels.find(l => l.h === h);
  }

  getMaxHeight() {
    return Math.max(...this.levels.map(l => l.h));
  }
}

// Merge curves data into existing FIELDS structure
export function mergeCurvesIntoFields(fieldsData, curvesJson, n) {
  
  const curvesNF = curvesJson?.catalogue?.number_fields
    || curvesJson?.number_fields
    || [];

  for (const curveField of curvesNF) {
    const D = curveField.D;
    const targetField = fieldsData.catalogue.getNumberField(String(D));
    
    if (!targetField) {
      console.warn(`No matching field found for discriminant ${D}`);
      continue;
    }
    
    // Get the tree for this n
    const targetTree = targetField.tree.find(t => t.n === n);
    if (!targetTree) {
      console.warn(`No tree found for n=${n} in field D=${D}`);
      continue;
    }
    
    const sourceTree = (curveField.tree ?? []).find(t => Number(t.n) === Number(n));
    const sourceClasses = sourceTree
      ? (sourceTree.I_t ?? sourceTree.ic ?? [])
      : (curveField.I_t ?? curveField.ic ?? []);
    if (!sourceClasses.length) continue;

    // Merge each isogeny class
    for (const curveIC of sourceClasses) {
      const sourceTrace = Number(curveIC.t ?? curveIC.trace);

      // Generator may store one class per |t|. Always merge curves into both ±t
      // so j-invariants are present in orders on D-select.
      const candidateTraces = Number.isFinite(sourceTrace)
        ? (sourceTrace === 0 ? [0] : [sourceTrace, -sourceTrace])
        : [curveIC.t ?? curveIC.trace];

      const uniqueTargetTraces = [...new Set(candidateTraces)];

      for (const targetTrace of uniqueTargetTraces) {
        const targetIC = targetTree.isogeny_classes.find(ic => Number(ic.trace) === Number(targetTrace));
        if (!targetIC) {
          continue;
        }

        // Keep volcano data sign-aware per target trace.
        if (curveIC.volcanoes) {
          targetIC.addVolcanoes(curveIC.volcanoes, targetTrace);
        }

        // Add curves to orders using EndomorphismOrder objects (old shape)
        if (curveIC.orders) {
          for (const curveOrder of curveIC.orders) {
            const targetOrder = targetIC.getOrder(curveOrder.f);
            if (targetOrder && curveOrder.curves) {
              targetOrder.addCurves(curveOrder.curves);
            }
          }
        }

        // Add curves to orders (new shape: curves map keyed by conductor)
        if (curveIC.curves) {
          for (const targetOrder of targetIC.orders) {
            const byF = curveIC.curves[String(targetOrder.conductor)]
              ?? curveIC.curves[targetOrder.conductor];
            if (Array.isArray(byF)) {
              targetOrder.addCurves(byF);
            }
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////// NR FIELD DATA STRUCTURES AND LOADING

export class IsogenyClass {
  constructor(data) {
    this.trace = data.t;  // Frobenius trace t
    this.D_pi = data.D_pi;  // Discriminant of Frobenius endomorphism (optional)
    this.f_pi = data.f_pi;  // Conductor of Frobenius (optional)
    // Backward/forward compatible:
    // - old shape: data.orders
    // - new shape: data.O
    this.orders = (data.orders ?? data.O ?? []).map(o => new EndomorphismOrder(o));

    // New shape may store curves grouped by conductor at class-level:
    // data.curves = { "1": [...], "3": [...] }
    const curvesByOrder = data.curves ?? data.curves_by_order ?? null;
    if (curvesByOrder) {
      this.orders.forEach(order => {
        const c = curvesByOrder[String(order.conductor)] ?? curvesByOrder[order.conductor];
        if (Array.isArray(c)) order.addCurves(c);
      });
    }

    this.volcanoes = this._buildVolcanoes(data.volcanoes ?? [], this.trace);
  }

  _buildVolcanoes(rawVolcanoes, traceValue) {
    const activeSign = Number(traceValue) < 0 ? '+' : '-';
    return (rawVolcanoes ?? [])
      .filter(v => (v?.[activeSign] === undefined ? true : !!v?.[activeSign]))
      .map(v => new Volcano({ ...v, _activeSign: activeSign }));
  }

  get conductor() {
    return this.f_pi;
  }

  get ordinary() {
    // If f_pi is present (old format), use it. Otherwise infer ordinary from having orders.
    if (this.f_pi !== undefined && this.f_pi !== null) return true;
    return this.orders.length > 0;
  }

  getOrder(conductor) {
    return this.orders.find(o => o.conductor === conductor);
  }

  addVolcanoes(volcanoesArray, sourceTrace = this.trace) {
    this.volcanoes = this._buildVolcanoes(volcanoesArray, sourceTrace);
  }

  getVolcano(ell) {
    return this.volcanoes.find(v => v.ell === ell);
  }
}

export class PrimePowerTree {
  constructor(data) {
    this.n = data.n;  // Power n in p^n
    // Backward/forward compatible:
    // - old shape: data.ic
    // - new shape: data.I_t
    const rawIsogenyClasses = (data.ic ?? data.I_t ?? []);
    const expanded = [];
    const seenTraces = new Set();

    for (const ic of rawIsogenyClasses) {
      const t = Number(ic?.t);

      // If trace is missing/non-numeric, keep entry as-is.
      if (!Number.isFinite(t)) {
        expanded.push(ic);
        continue;
      }

      if (!seenTraces.has(t)) {
        expanded.push(ic);
        seenTraces.add(t);
      }

      // Your generator stores one class per |t|. Mirror it to -t on load.
      if (t !== 0 && !seenTraces.has(-t)) {
        const mirrored = JSON.parse(JSON.stringify(ic));
        mirrored.t = -t;
        expanded.push(mirrored);
        seenTraces.add(-t);
      }
    }

    this.isogeny_classes = expanded.map(ic => new IsogenyClass(ic));
  }

  getIsogenyClass(trace) {
    return this.isogeny_classes.find(ic => ic.trace === trace);
  }

  getTraces() {
    return this.isogeny_classes.map(ic => ic.trace).sort((a, b) => a - b);
  }

  getUniqueTraces() {
    return [...new Set(this.getTraces())];
  }
}

export class NumberField {
  constructor(data) {
    this.discriminant = data.D;
    this.tree = (data.tree ?? []).map(t => new PrimePowerTree(t));
  }

  getTreeByN(n) {
    const branch = this.tree.find(t => t.n === n);
    // Return just the isogeny classes, not the full tree structure
    return branch ? branch.isogeny_classes : null;
  }

  getAllIsogenyClasses() {
    return this.tree.flatMap(t => t.isogeny_classes);
  }

  getTracesForN(n) {
    const isogenyClasses = this.getTreeByN(n);
    return isogenyClasses ? [...new Set(isogenyClasses.map(ic => ic.trace))].sort((a, b) => a - b) : [];
  }
}

export class NumberFieldCatalogue {
  constructor(data) {
    this.number_fields = [];
    
    if (data?.nf) {
      for (const nf_data of data.nf) {
        const field = new NumberField(nf_data);
        this.number_fields.push(field);
      }
    }
  }

  getNumberField(discriminant) {
    return this.number_fields.find(nf => String(nf.discriminant) === String(discriminant));
  }

  getAllNumberFields() {
    return this.number_fields;
  }

  getFieldsByN(n) {
    const fields = [];
    for (const field of this.number_fields) {
      const isogenyClasses = field.getTreeByN(n);
      if (isogenyClasses && isogenyClasses.length > 0) {
        // Return a new object with only the data for this n
        fields.push({
          discriminant: field.discriminant,
          isogeny_classes: isogenyClasses
        });
      }
    }
    return fields;
  }

  getDiscriminants() {
    return this.number_fields.map(nf => nf.discriminant).sort((a, b) => Math.abs(a) - Math.abs(b));
  }
}

export class NumberFieldData {
  constructor(data) {
    this.char = data.char;  // Characteristic p
    this.catalogue = new NumberFieldCatalogue(data.nr_fields);
  }

  getNumberField(discriminant) {
    return this.catalogue.getNumberField(discriminant);
  }

  getAllNumberFields() {
    return this.catalogue.getAllNumberFields();
  }

  getFieldsByN(n) {
    return this.catalogue.getFieldsByN(n);
  }

  getDiscriminants() {
    return this.catalogue.getDiscriminants();
  }

  // Get all available n values across all fields
  getAvailableNValues() {
    const nValues = new Set();
    for (const field of this.getAllNumberFields()) {
      for (const tree of field.tree) {
        nValues.add(tree.n);
      }
    }
    return Array.from(nValues).sort((a, b) => a - b);
  }

  // Statistics
  getTotalIsogenyClasses(n = null) {
    if (n === null) {
      return this.getAllNumberFields().reduce(
        (sum, field) => sum + field.getAllIsogenyClasses().length,
        0
      );
    } else {
      return this.getFieldsByN(n).reduce(
        (sum, field) => sum + field.isogeny_classes.length,
        0
      );
    }
  }

  getTracesForN(n) {
    const traces = new Set();
    for (const field of this.getFieldsByN(n)) {
      for (const ic of field.isogeny_classes) {
        traces.add(ic.trace);
      }
    }
    return Array.from(traces).sort((a, b) => a - b);
  }
}

// Factory function to parse JSON into structured data
export function parseNumberFieldData(json) {
  return new NumberFieldData(json);
}

