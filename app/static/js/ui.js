import { sievePrimes, setStatus } from '/static/js/utils.js';

export class UI {
  constructor(P, N) {
    this.P = P;
    this.N = N;
    this.q = P ** N;
    this.selectedEll = 'ALL';
    
    this.statusEl = document.getElementById('status');
    this.nSelect = document.getElementById('n_select');
    this.ellSelect = document.getElementById('ell_select');
    this.selectedFieldBarEl = document.getElementById('selected_field_bar');
    this.selectedFieldLabelEl = document.getElementById('selected_field_label');
    this.selectedFieldValueEl = document.getElementById('selected_field_value');
    this.selectedCurveBarEl = document.getElementById('selected_curve_bar');
    this.selectedEllLabelEl = document.getElementById('selected_ell_label');
    this.selectedEllValueEl = document.getElementById('selected_ell_value');
    this.selectedCurveLabelEl = document.getElementById('selected_curve_label');
    this.selectedCurveValueEl = document.getElementById('selected_curve_value');
    this.divVal = document.getElementById('divVal');
    this.dInput = document.getElementById('d_input');
    this.jInput = document.getElementById('j_input');
    this.aInput = document.getElementById('a_input');
    this.bInput = document.getElementById('b_input');
    this.layoutToggle = document.getElementById('layout_toggle');
    
    const pPill = document.getElementById('pPill');
    if (pPill) pPill.textContent = ` ${P} `;
    
    const nSup = document.getElementById('nSup');
    if (nSup) nSup.textContent = ` ${N} `;
    
    this.activeDiv = undefined;
    this.onDivChange = null; // callback
    this.onNChange = null; // callback for n value change
    this.onEllChange = null; // callback for ell value change
    this.onVolcanoRowClick = null; // callback for volcano row click (ell, trace)
    this.onDSelect = null; // callback for D input
    this.onJSelect = null; // callback for j input
    this.onASelect = null; // callback for A input
    this.onBSelect = null; // callback for B input
    this.onLayoutToggle = null; // callback for layout mode toggle
    this.piConductorColor = '#3da3ff';
    this.activeVolcanoEll = null;
    this.activeVolcanoTrace = null;
  }

  setActiveVolcano(ell = null, trace = null) {
    this.activeVolcanoEll = ell === null || ell === undefined || ell === 'ALL' ? null : Number(ell);
    this.activeVolcanoTrace = trace === null || trace === undefined ? null : Number(trace);
  }

  setAvailableNValues(nValues) {
    if (!this.nSelect) return;
    this.nSelect.innerHTML = '';
    const rawValues = [...new Set((nValues || []).map(Number).filter(Number.isFinite))].sort((a, b) => a - b);
    const maxN = rawValues.length ? Math.max(...rawValues) : this.N;
    const values = Array.from({ length: Math.max(1, maxN) }, (_, i) => i + 1);
    values.forEach((n) => {
      const option = document.createElement('option');
      option.value = String(n);
      option.textContent = `𝔽_${this.P ** n} (p=${this.P}, n=${n})`;
      if (n === this.N) option.selected = true;
      this.nSelect.appendChild(option);
    });
  }

  setPiConductorColor(color) {
    if (color) this.piConductorColor = color;
  }

  init() {
    if (this.nSelect) {
      this.nSelect.addEventListener('change', (e) => {
        this.handleNSelectChange(Number(e.target.value));
      });
    }
    this.updateDisplay();
    
    // Add event listener for ℓ select dropdown
    if (this.ellSelect) {
      this.ellSelect.addEventListener('change', (e) => this.handleEllChange(e.target.value));
    }
    
    // Add event listeners for D and j inputs
    if (this.dInput) {
      this.dInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          const d = parseInt(this.dInput.value);
          if (!isNaN(d) && this.onDSelect) {
            this.onDSelect(d);
          }
        }
      });
    }
    
    if (this.jInput) {
      this.jInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          const jStr = this.jInput.value.trim();
          if (jStr && this.onJSelect) {
            this.onJSelect(jStr);
          }
        }
      });
    }
    
    if (this.aInput) {
      this.aInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          const aStr = this.aInput.value.trim();
          if (aStr && this.onASelect) {
            this.onASelect(aStr);
          }
        }
      });
    }
    
    if (this.bInput) {
      this.bInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          const bStr = this.bInput.value.trim();
          if (bStr && this.onBSelect) {
            this.onBSelect(bStr);
          }
        }
      });
    }
    
    // Add event listener for layout toggle button
    if (this.layoutToggle) {
      this.layoutToggle.addEventListener('click', () => {
        if (this.onLayoutToggle) {
          this.onLayoutToggle();
        }
      });
    }
  }

  handleNSelectChange(newN) {
    this.updateDisplay();
    
    // Only trigger reload if N actually changed
    if (newN !== this.N && this.onNChange) {
      this.N = newN;
      this.q = this.P ** newN;  // Update q when N changes
      this.onNChange(newN);
    }
  }

  updateDisplay() {
    const n = this.nSelect ? Number(this.nSelect.value || this.N) : this.N;
    this.q = this.P ** n;
    if (this.divVal) this.divVal.textContent = ` ${this.q} `;
    const nSup = document.getElementById('nSup');
    if (nSup) nSup.textContent = ` ${n} `;
  }

  updateStatus(nf) {

    if (!nf) {
      setStatus(this.statusEl, ``);
      return;
    }
    const { isogeny_classes } = nf;

    if (!isogeny_classes || isogeny_classes.length === 0) {
      setStatus(this.statusEl, 'No isogeny classes');
      return;
    }

    // Helper function to compute prime factorization
    const primeFactorization = (n) => {
      const factors = [];
      let d = 2;
      while (d * d <= n) {
        while (n % d === 0) {
          factors.push(d);
          n /= d;
        }
        d++;
      }
      if (n > 1) factors.push(n);
      
      // Group factors by power: e.g., [2,2,3] -> "2²×3"
      const grouped = {};
      factors.forEach(f => grouped[f] = (grouped[f] || 0) + 1);
      
      // Map for superscript digits
      const superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'};
      
      return Object.entries(grouped)
        .map(([p, exp]) => {
          if (exp === 1) return p;
          const expStr = String(exp).split('').map(d => superscripts[d]).join('');
          return `${p}${expStr}`;
        })
        .join('·');
    };

    const classCount = Math.max(1, isogeny_classes.length);
    const availableWidth = Math.max(260, Math.floor(window.innerWidth * 0.5) - 40);
    const gapPx = Math.max(2, Math.min(8, Math.floor(availableWidth / Math.max(16 * classCount, 1))));
    const columnWidthPx = Math.max(92, Math.floor((availableWidth - gapPx * (classCount - 1)) / classCount));

    let infoStr = '';
    infoStr += `<div style="display: flex; gap: ${gapPx}px; flex-wrap: nowrap; align-items: flex-start; margin-top: 4px; width: 100%; white-space: nowrap;">`;
    
    [...isogeny_classes].reverse().forEach(ic => {
      const t = ic.trace;
      const N_pts = this.q + 1 - t;
      const factorization = primeFactorization(N_pts);
      // Create a column for each isogeny class
      infoStr += `<div style="background: transparent; padding: 2px 4px; border-radius: 6px; width: ${columnWidthPx}px; min-width: ${columnWidthPx}px; white-space: nowrap; overflow: hidden; text-align: center;">`;
      infoStr += `<strong>t: ${t}</strong>`;
      infoStr += ` | #E: ${factorization}`;
      
      
      // Show volcanoes for this isogeny class
      if (ic.volcanoes && ic.volcanoes.length > 0) {
        infoStr += `<div style="margin-top: 4px;">`;
        ic.volcanoes.forEach(volcano => {
          const fx_roots = volcano.fx_roots ? `[${volcano.fx_roots.join(', ')}]` : '[]';
          // Calculate max height from levels
          const maxHeight = volcano.levels && volcano.levels.length > 0 
            ? Math.max(...volcano.levels.map(level => level.h))
            : 0;
          // Use theme-aware pi-conductor color if height > 0
          const isActiveVolcano = Number(this.activeVolcanoEll) === Number(volcano.ell)
            && Number(this.activeVolcanoTrace) === Number(t);
          const rowColor = isActiveVolcano
            ? '#ff00ff'
            : (maxHeight > 0 ? this.piConductorColor : 'inherit');
          const rowStyle = `color: ${rowColor}; cursor: pointer; text-decoration: underline; font-size: 11px;`;
          infoStr += `<span class="volcano-row" data-ell="${volcano.ell}" data-trace="${t}" style="${rowStyle}">ℓ=${volcano.ell}: ${fx_roots}, h=${maxHeight}</span><br>`;
        });
        infoStr += `</div>`;
      }
      infoStr += '</div>';
    });
    
    infoStr += '</div></div>';
    
    setStatus(this.statusEl, infoStr);
  }

  showError(message) {
    setStatus(this.statusEl, message);
  }

  populateEllSelect(ellValues) {
    if (!this.ellSelect) return;
    
    // Remove all options except ALL
    while (this.ellSelect.options.length > 1) {
      this.ellSelect.remove(1);
    }

    if (this.ellSelect.options.length > 0) {
      this.ellSelect.options[0].textContent = 'ℓ';
    }
    
    // Add options for each ell value
    ellValues.sort((a, b) => a - b).forEach(ell => {
      const option = document.createElement('option');
      option.value = ell;
      option.textContent = `ℓ=${ell}`;
      
      this.ellSelect.appendChild(option);
    });
    
    // Reset to ALL
    this.ellSelect.value = 'ALL';
    this.selectedEll = 'ALL';
  }

  handleEllChange(ell) {
    this.selectedEll = ell;
    if (this.onEllChange) {
      this.onEllChange(ell);
    }
  }

  updateSelectedBarEmptyState() {
    if (!this.selectedCurveBarEl) return;
    const hasCurve = Boolean(this.selectedCurveValueEl && this.selectedCurveValueEl.textContent && this.selectedCurveValueEl.textContent !== '—');
    const hasEll = Boolean(this.selectedEllValueEl && this.selectedEllValueEl.textContent && this.selectedEllValueEl.textContent !== '—');
    this.selectedCurveBarEl.classList.toggle('is-empty', !hasCurve && !hasEll);
  }

  setSelectedEll(value, label = 'Volcano', { isHtml = false } = {}) {
    if (this.selectedEllLabelEl) {
      this.selectedEllLabelEl.textContent = label;
    }
    if (!this.selectedEllValueEl) return;

    const hasValue = typeof value === 'string' && value.trim().length > 0;
    if (isHtml && hasValue) {
      this.selectedEllValueEl.innerHTML = value;
    } else {
      this.selectedEllValueEl.textContent = hasValue ? value : '—';
    }
    this.updateSelectedBarEmptyState();
  }

  clearSelectedEll() {
    this.setSelectedEll('');
  }

  setSelectedField(value, label = 'Selected field') {
    if (this.selectedFieldLabelEl) {
      this.selectedFieldLabelEl.textContent = label;
    }
    if (!this.selectedFieldBarEl || !this.selectedFieldValueEl) return;

    const hasValue = typeof value === 'string' && value.trim().length > 0;
    this.selectedFieldValueEl.textContent = hasValue ? value : '—';
    this.selectedFieldBarEl.classList.toggle('is-empty', !hasValue);
  }

  clearSelectedField() {
    this.setSelectedField('');
  }

  setSelectedCurveJInvariant(value, label = 'Selected curve j') {
    if (this.selectedCurveLabelEl) {
      this.selectedCurveLabelEl.textContent = label;
    }
    if (!this.selectedCurveBarEl || !this.selectedCurveValueEl) return;

    const hasValue = typeof value === 'string' && value.trim().length > 0;
    this.selectedCurveValueEl.textContent = hasValue ? value : '—';
    this.updateSelectedBarEmptyState();
  }

  clearSelectedCurveJInvariant() {
    this.setSelectedCurveJInvariant('');
  }
}
