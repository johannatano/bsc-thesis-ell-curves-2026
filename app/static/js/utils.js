
// ==================== UTILITIES ====================
export function sievePrimes(max) {
  const sieve = new Uint8Array(max + 1), out = [1];
  for (let i = 2; i <= max; i++) if (!sieve[i]) { out.push(i); for (let j = i * i; j <= max; j += i) sieve[j] = 1; }
  return out;
}

export function setStatus(statusEl, txt) { 
  statusEl.innerHTML = txt;

  if (!statusEl.dataset.volcanoClickBound) {
    statusEl.addEventListener('click', (ev) => {
      const row = ev.target?.closest?.('.volcano-row');
      if (!row || !statusEl.contains(row)) return;
      const ell = row.getAttribute('data-ell');
      const trace = row.getAttribute('data-trace');
      const event = new CustomEvent('volcanoRowClick', {
        detail: { ell: Number(ell), trace: Number(trace) }
      });
      window.dispatchEvent(event);
    });
    statusEl.dataset.volcanoClickBound = 'true';
  }
}

export function padToN(arr, n) {
  const out = arr.slice(0, n);
  while (out.length < n) out.push(0);
  return out;
}