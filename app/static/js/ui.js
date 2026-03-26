import { sievePrimes, setStatus } from '/static/js/utils.js';

export class UI {
  constructor(P, N) {
    this.P = P;
    this.N = N;
    this.q = P ** N;
    this.selectedEll = 'ALL';
    
    this.statusEl = document.getElementById('status');
    this.nButtons = document.querySelectorAll('input[name="n_select"]');
    this.ellSelect = document.getElementById('ell_select');
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
  }

  setPiConductorColor(color) {
    if (color) this.piConductorColor = color;
  }

  init() {
    // Set initial selected button
    this.nButtons.forEach(btn => {
      if (Number(btn.value) === this.N) {
        btn.checked = true;
      }
      btn.addEventListener('change', () => this.handleButtonChange(btn));
    });
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

  handleButtonChange(btn) {
    const newN = Number(btn.value);
    this.updateDisplay();
    
    // Only trigger reload if N actually changed
    if (newN !== this.N && this.onNChange) {
      this.N = newN;
      this.q = this.P ** newN;  // Update q when N changes
      this.onNChange(newN);
    }
  }

  updateDisplay() {
    const selectedBtn = Array.from(this.nButtons).find(btn => btn.checked);
    const n = selectedBtn ? Number(selectedBtn.value) : this.N;
    if (this.divVal) this.divVal.textContent = ` ${this.P}^${n} `;
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

    let infoStr = '<strong>Isogeny Classes:</strong><br>';
    infoStr += '<div style="display: flex; gap: 16px; flex-wrap: wrap; margin-top: 6px;">';
    
    isogeny_classes.forEach(ic => {
      const t = ic.trace;
      const N_pts = this.q + 1 - t;
      const f_pi = ic.f_pi || 'N/A';
      const factorization = primeFactorization(N_pts);
      const factorization_f_pi = primeFactorization(f_pi);
      // Create a column for each isogeny class
      infoStr += '<div style="background: rgba(255,255,255,0.05); padding: 8px; border-radius: 6px; min-width: 240px;">';
      infoStr += `<strong>t: ${t}</strong>`;
      infoStr += ` | f_π: ${f_pi} = ${factorization_f_pi}<br>`;
      infoStr += `#E: ${N_pts} = `;
      infoStr += `${factorization}`;
      
      
      // Show volcanoes for this isogeny class
      if (ic.volcanoes && ic.volcanoes.length > 0) {
        infoStr += `<div style="margin-top: 4px;">`;
        infoStr += `<strong>Volcanoes:</strong><br>`;
        ic.volcanoes.forEach(volcano => {
          const fx_roots = volcano.fx_roots ? `[${volcano.fx_roots.join(', ')}]` : '[]';
          // Calculate max height from levels
          const maxHeight = volcano.levels && volcano.levels.length > 0 
            ? Math.max(...volcano.levels.map(level => level.h))
            : 0;
          // Calculate boolean checks: ℓ | (q-1) and ℓ | (t-2)
          const dividesQMinus1 = (this.q - 1) % volcano.ell === 0;
          const dividesTMinus2 = (t - 2) % volcano.ell === 0;
          const check1 = dividesQMinus1 ? '✓' : '✗';
          const check2 = dividesTMinus2 ? '✓' : '✗';
          // Use theme-aware pi-conductor color if height > 0
          const rowColor = maxHeight > 0 ? this.piConductorColor : 'inherit';
          infoStr += `<span class="volcano-row" data-ell="${volcano.ell}" data-trace="${t}" style="color: ${rowColor}; cursor: pointer; text-decoration: underline;">ℓ=${volcano.ell}: ${fx_roots}, h=${maxHeight}, ℓ|(q-1): ${check1}, ℓ|(t-2): ${check2}</span><br>`;
        });
        infoStr += `</div>`;
      }
      infoStr += '</div>';
    });
    
    infoStr += '</div>';
    
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
    
    // Add options for each ell value
    ellValues.sort((a, b) => a - b).forEach(ell => {
      const option = document.createElement('option');
      option.value = ell;
      
      // Check if ell divides q-1 and mark with indicator
      const dividesQMinus1 = (this.q - 1) % ell === 0;
      option.textContent = dividesQMinus1 ? `ℓ=${ell} ✓` : `ℓ=${ell}`;
      
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
}
