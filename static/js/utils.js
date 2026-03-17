
// ==================== UTILITIES ====================
export function sievePrimes(max) {
  const sieve = new Uint8Array(max + 1), out = [1];
  for (let i = 2; i <= max; i++) if (!sieve[i]) { out.push(i); for (let j = i * i; j <= max; j += i) sieve[j] = 1; }
  return out;
}

export function setStatus(statusEl, txt) { 
  statusEl.innerHTML = txt;
  
  // Add click handlers for volcano rows after content is set
  setTimeout(() => {
    const volcanoRows = statusEl.querySelectorAll('.volcano-row');
    volcanoRows.forEach(row => {
      row.addEventListener('click', () => {
        const ell = row.getAttribute('data-ell');
        const trace = row.getAttribute('data-trace');
        // Dispatch custom event that app.js will listen to
        const event = new CustomEvent('volcanoRowClick', {
          detail: { ell: Number(ell), trace: Number(trace) }
        });
        window.dispatchEvent(event);
      });
    });
  }, 0);
}

export function padToN(arr, n) {
  const out = arr.slice(0, n);
  while (out.length < n) out.push(0);
  return out;
}