/* ============================================================
   PSYCHOMETRIC CALCULATORS · Clinical Suite
   ============================================================ */

/* ---------- STATISTICS ---------- */
// Standard normal CDF using Abramowitz & Stegun erf approximation
function erf(x){
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.exp(-x*x);
  return sign * y;
}
function normCDF(z){ return 0.5 * (1 + erf(z / Math.SQRT2)); }
function normPDF(z){ return Math.exp(-0.5*z*z) / Math.sqrt(2*Math.PI); }

/* ---- Student's t-distribution helpers (used by Crawford RCI) ---- */
// log-gamma via Lanczos approximation
function logGamma(z){
  const g = 7;
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
             -176.61502916214059, 12.507343278686905, -0.13857109526572012,
             9.9843695780195716e-6, 1.5056327351493116e-7];
  if (z < 0.5){
    return Math.log(Math.PI / Math.sin(Math.PI*z)) - logGamma(1 - z);
  }
  z -= 1;
  let x = c[0];
  for (let i = 1; i < g + 2; i++) x += c[i] / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
// regularised incomplete beta I_x(a,b) - continued-fraction form (Numerical Recipes-style)
function _betacf(x, a, b){
  const MAXIT = 200, EPS = 3e-7, FPMIN = 1e-30;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - qab*x/qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1/d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++){
    const m2 = 2*m;
    let aa = m*(b-m)*x / ((qam+m2)*(a+m2));
    d = 1 + aa*d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa/c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1/d; h *= d*c;
    aa = -(a+m)*(qab+m)*x / ((a+m2)*(qap+m2));
    d = 1 + aa*d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa/c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1/d;
    const del = d*c; h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}
function ibeta(x, a, b){
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lnBT = logGamma(a+b) - logGamma(a) - logGamma(b) + a*Math.log(x) + b*Math.log(1-x);
  if (x < (a+1)/(a+b+2)) return Math.exp(lnBT) * _betacf(x, a, b) / a;
  return 1 - Math.exp(lnBT) * _betacf(1-x, b, a) / b;
}
// CDF of Student's t with df degrees of freedom
function tCDF(t, df){
  if (df <= 0 || !isFinite(df)) return NaN;
  const x = df / (df + t*t);
  const ib = ibeta(x, df/2, 0.5);
  return t >= 0 ? 1 - 0.5*ib : 0.5*ib;
}
// inverse CDF (quantile) of Student's t - bisection
function tInv(p, df){
  if (p <= 0 || p >= 1 || df <= 0) return NaN;
  if (p === 0.5) return 0;
  // bracket: |t| <= 50 covers any practical df
  let lo = -50, hi = 50;
  for (let i = 0; i < 80; i++){
    const mid = 0.5*(lo+hi);
    if (tCDF(mid, df) < p) lo = mid; else hi = mid;
  }
  return 0.5*(lo+hi);
}
// Acklam's inverse normal CDF
function normInv(p){
  if (p <= 0 || p >= 1) return p <= 0 ? -Infinity : Infinity;
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const plow = 0.02425, phigh = 1 - plow;
  let q, r;
  if (p < plow){
    q = Math.sqrt(-2*Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  if (p <= phigh){
    q = p - 0.5; r = q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
  q = Math.sqrt(-2*Math.log(1-p));
  return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
}

/* ---------- SCORE CONVERSIONS ---------- */
// Convert any input score type to z
function toZ(value, type){
  if (value === '' || value == null || isNaN(value)) return null;
  const v = parseFloat(value);
  switch(type){
    case 'z': return v;
    case 't': return (v - 50) / 10;
    case 'standard': return (v - 100) / 15;
    case 'scaled': return (v - 10) / 3;
    case 'percentile':
      if (v <= 0 || v >= 100) return null;
      return normInv(v / 100);
    default: return null;
  }
}
function fromZ(z, type){
  if (z == null || isNaN(z)) return null;
  switch(type){
    case 'z': return z;
    case 't': return 50 + 10*z;
    case 'standard': return 100 + 15*z;
    case 'scaled': return 10 + 3*z;
    case 'percentile': return normCDF(z) * 100;
    default: return null;
  }
}
function fmt(n, dp = 2){
  if (n == null || isNaN(n)) return '-';
  if (Math.abs(n) < 0.0001 && n !== 0) return n.toExponential(2);
  return n.toFixed(dp);
}
function fmtPct(p){
  if (p == null || isNaN(p)) return '-';
  if (p < 0.1) return p.toFixed(2);
  if (p > 99.9) return p.toFixed(2);
  if (p < 1 || p > 99) return p.toFixed(1);
  return Math.round(p).toString();
}
function fmtP(p){
  if (p == null || isNaN(p)) return '-';
  if (p < 0.001) return '< .001';
  return p.toFixed(3).replace(/^0\./, '.');
}

/* ---------- DESCRIPTORS ---------- */
// Wechsler classification (based on Standard Score)
function wechslerDesc(ss){
  if (ss == null) return '-';
  if (ss >= 130) return 'Very Superior';
  if (ss >= 120) return 'Superior';
  if (ss >= 110) return 'High Average';
  if (ss >= 90)  return 'Average';
  if (ss >= 80)  return 'Low Average';
  if (ss >= 70)  return 'Borderline';
  return 'Extremely Low';
}
// AACN classification (Guilmette et al., 2020)
function aanDesc(ss){
  if (ss == null) return '-';
  if (ss >= 130) return 'Exceptionally High';
  if (ss >= 120) return 'Above Average';
  if (ss >= 110) return 'High Average';
  if (ss >= 90)  return 'Average';
  if (ss >= 80)  return 'Low Average';
  if (ss >= 70)  return 'Below Average';
  return 'Exceptionally Low';
}

// Descriptor carousel ───────────────────────────────────────────────────────
const WECHSLER_LEVELS = [
  { label:'Extremely Low', range:'< 70'    },
  { label:'Borderline',    range:'70–79'   },
  { label:'Low Average',   range:'80–89'   },
  { label:'Average',       range:'90–109'  },
  { label:'High Average',  range:'110–119' },
  { label:'Superior',      range:'120–129' },
  { label:'Very Superior', range:'≥ 130'   }
];
const AACN_LEVELS = [
  { label:'Exceptionally Low',  range:'< 70'    },
  { label:'Below Average',      range:'70–79'   },
  { label:'Low Average',        range:'80–89'   },
  { label:'Average',            range:'90–109'  },
  { label:'High Average',       range:'110–119' },
  { label:'Above Average',      range:'120–129' },
  { label:'Exceptionally High', range:'≥ 130'   }
];
const DESC_PILL_W = 116; // must match CSS flex-basis on .conv-desc-pill
const DESC_MID    = 3;   // index of the middle pill (0-based, 7 items → index 3)
// Red → neutral → green scale across the 7 descriptor bands
const DESC_COLOURS = [
  '#9C3D2A', // Extremely Low   - deep red
  '#B5631C', // Borderline      - burnt orange
  '#A88818', // Low Average     - amber
  '#6B7A5C', // Average         - olive/neutral
  '#3D7550', // High Average    - muted green
  '#2A6640', // Superior        - medium green
  '#1A5430', // Very Superior   - deep green
];
function ssToDescIndex(ss){
  if (ss < 70)  return 0;
  if (ss < 80)  return 1;
  if (ss < 90)  return 2;
  if (ss < 110) return 3;
  if (ss < 120) return 4;
  if (ss < 130) return 5;
  return 6;
}
function buildDescCarousels(){
  buildDescCarousel('conv-wechsler-block', WECHSLER_LEVELS);
  buildDescCarousel('conv-aan-block', AACN_LEVELS);
}
function buildDescCarousel(id, levels){
  const block = document.getElementById(id);
  if (!block) return;
  const track = document.createElement('div');
  track.className = 'desc-carousel-track';
  levels.forEach(l => {
    const pill = document.createElement('div');
    pill.className = 'conv-desc-pill';
    pill.innerHTML = `<span class="pill-label">${l.label}</span><span class="pill-range">${l.range}</span>`;
    track.appendChild(pill);
  });
  block.innerHTML = '';
  block.appendChild(track);
}
function updateDescCarousel(id, activeIdx){
  const block = document.getElementById(id);
  if (!block) return;
  const track = block.querySelector('.desc-carousel-track');
  if (!track) return;
  track.querySelectorAll('.conv-desc-pill').forEach((p, i) => {
    const d = Math.abs(i - activeIdx);
    p.className = 'conv-desc-pill' + (d === 0 ? ' active' : d === 1 ? ' adj-1' : d === 2 ? ' adj-2' : '');
    p.style.color = d === 0 ? (DESC_COLOURS[activeIdx] || '') : '';
  });
  // Track starts centred (CSS justify-content:center); shift so activeIdx pill lands at centre
  track.style.transform = `translateX(${(DESC_MID - activeIdx) * DESC_PILL_W}px)`;
}

/* ---------- NAVIGATION ---------- */
function setNavGroupOpen(group, isOpen){
  if (!group) return;
  group.classList.toggle('is-collapsed', !isOpen);
  const label = group.querySelector('.nav-label');
  if (label) label.setAttribute('aria-expanded', String(isOpen));
}

function openOnlyNavGroup(group){
  document.querySelectorAll('.nav-group').forEach(g => setNavGroupOpen(g, g === group));
}

const activeNavGroup = document.querySelector('.nav-item.active')?.closest('.nav-group');
if (activeNavGroup) openOnlyNavGroup(activeNavGroup);

const CUSTOM_TESTS_PASSWORD = 'unlock';
const CUSTOM_TESTS_UNLOCK_KEY = 'customTestsUnlocked';
function isCustomTestsUnlocked(){
  return sessionStorage.getItem(CUSTOM_TESTS_UNLOCK_KEY) === 'true';
}
function requestCustomTestsUnlock(){
  const entered = window.prompt('Custom Tests is password-protected. Enter password:');
  if (entered == null) return false;
  if (entered === CUSTOM_TESTS_PASSWORD){
    sessionStorage.setItem(CUSTOM_TESTS_UNLOCK_KEY, 'true');
    showToast('✓ Custom Tests unlocked for this session');
    return true;
  }
  showToast('Incorrect password', true);
  return false;
}

document.querySelectorAll('.nav-label').forEach(label => {
  label.addEventListener('click', () => {
    const group = label.closest('.nav-group');
    const willOpen = group.classList.contains('is-collapsed');
    if (willOpen) {
      openOnlyNavGroup(group);
    } else {
      setNavGroupOpen(group, false);
    }
  });
});

document.querySelectorAll('.nav-item').forEach(el => {
  el.addEventListener('click', () => {
    const target = el.dataset.target;
    if (target === 'custom-tests' && !isCustomTestsUnlocked()){
      const ok = requestCustomTestsUnlock();
      if (!ok) return;
    }
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));
    document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
    el.classList.add('active');
    const targetSection = document.getElementById(target);
    if (targetSection) targetSection.classList.add('active');
    openOnlyNavGroup(el.closest('.nav-group'));
    document.querySelector('.main').scrollTop = 0;
  });
});

/* ---------- TOAST ---------- */
let toastTimer;
function showToast(msg, isError){
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = 'toast show' + (isError ? ' error' : '');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => t.className = 'toast', 2200);
}

/* ---------- COPY TO CLIPBOARD (rich HTML) ---------- */
async function copyApaTable(containerId){
  const container = document.getElementById(containerId);
  if (!container || !container.innerHTML.trim()){
    showToast('No table to copy yet - fill in some data first', true);
    return;
  }
  // Build standalone HTML with inline styles for Word/Docs paste
  const html = buildStandaloneHtml(container);
  const plain = htmlToPlain(container);
  try {
    if (navigator.clipboard && window.ClipboardItem){
      const item = new ClipboardItem({
        'text/html': new Blob([html], { type: 'text/html' }),
        'text/plain': new Blob([plain], { type: 'text/plain' })
      });
      await navigator.clipboard.write([item]);
    } else {
      await navigator.clipboard.writeText(plain);
    }
    showToast('✓ Table copied - ready to paste into your report');
    if (typeof ReportBundle !== 'undefined' && ReportBundle.showKofiPrompt) ReportBundle.showKofiPrompt();
  } catch(e){
    console.error(e);
    showToast('Copy failed - try selecting and copying manually', true);
  }
}
function buildStandaloneHtml(container){
  // Inline styles so Word/Docs render APA formatting
  const clone = container.cloneNode(true);
  // Title block
  clone.querySelectorAll('.apa-table-num').forEach(el => {
    el.setAttribute('style', "font-family:'Times New Roman',serif;font-size:11pt;font-style:italic;font-weight:normal;color:#000;margin:0 0 2pt 0;");
  });
  clone.querySelectorAll('.apa-table-title').forEach(el => {
    el.setAttribute('style', "font-family:'Times New Roman',serif;font-size:11pt;font-weight:normal;font-style:italic;color:#000;margin:0 0 8pt 0;line-height:1.4;");
  });
  // Table
  clone.querySelectorAll('.apa-table').forEach(t => {
    t.setAttribute('style', "border-collapse:collapse;font-family:'Times New Roman',serif;font-size:11pt;color:#000;width:auto;");
    t.setAttribute('cellpadding', '3');
    t.setAttribute('cellspacing', '0');
  });
  // Header rows: top double rule via top border, single bottom rule under header
  clone.querySelectorAll('.apa-table thead tr').forEach((tr, i, arr) => {
    if (i === 0){
      tr.querySelectorAll('th').forEach(th => {
        th.setAttribute('style', "border-top:1.5pt solid #000;padding:3pt 7pt;font-weight:normal;text-align:left;font-family:'Times New Roman',serif;line-height:1.05;");
      });
    }
    if (i === arr.length - 1){
      tr.querySelectorAll('th').forEach(th => {
        const existing = th.getAttribute('style') || '';
        th.setAttribute('style', existing + 'border-bottom:0.5pt solid #000;padding:3pt 7pt;font-weight:normal;text-align:left;font-family:\'Times New Roman\',serif;line-height:1.05;');
      });
    }
  });
  // Numeric headers center
  clone.querySelectorAll('.apa-table th.num').forEach(th => {
    const existing = th.getAttribute('style') || '';
    th.setAttribute('style', existing + 'text-align:center;');
  });
  // Body cells
  clone.querySelectorAll('.apa-table tbody td').forEach(td => {
    td.setAttribute('style', "padding:2.5pt 7pt;border:none;font-family:'Times New Roman',serif;color:#000;line-height:1.05;");
  });
  clone.querySelectorAll('.apa-table tbody td.num').forEach(td => {
    const existing = td.getAttribute('style') || '';
    td.setAttribute('style', existing + 'text-align:center;');
  });
  // Group separator rows (bold italic, no border)
  clone.querySelectorAll('.apa-table tbody tr.apa-group td').forEach(td => {
    td.setAttribute('style', "padding:5pt 7pt 2.5pt;border:none;font-family:'Times New Roman',serif;color:#000;font-style:italic;font-weight:bold;line-height:1.05;");
  });
  // Indent the first cell of grouped subtest rows
  clone.querySelectorAll('.apa-table tbody tr.apa-grouped-row td:first-child').forEach(td => {
    const existing = td.getAttribute('style') || '';
    td.setAttribute('style', existing + 'padding-left:24pt;');
  });
  // Bottom border on last body row
  const lastRows = clone.querySelectorAll('.apa-table tbody tr:last-child td');
  lastRows.forEach(td => {
    const existing = td.getAttribute('style') || '';
    td.setAttribute('style', existing + 'border-bottom:1.5pt solid #000;padding-bottom:3pt;');
  });
  // Preserve emphasis on the expected-range star markers, but DO NOT bold
  // .bat-class-extreme - extremely-low scores should render in normal weight
  // alongside every other classification cell.
  clone.querySelectorAll('.bat-expected-stars').forEach(el => {
    const existing = el.getAttribute('style') || '';
    el.setAttribute('style', existing + 'font-weight:bold;background:transparent;');
  });
  clone.querySelectorAll('.bat-class-extreme').forEach(el => {
    const existing = el.getAttribute('style') || '';
    el.setAttribute('style', existing + 'font-weight:normal;background:transparent;');
  });
  // Note
  clone.querySelectorAll('.apa-note').forEach(n => {
    n.setAttribute('style', "font-family:'Times New Roman',serif;font-size:10pt;font-style:italic;color:#000;margin-top:8pt;line-height:1.4;");
  });
  // Wrap
  return `<html><head><meta charset="utf-8"></head><body>${clone.innerHTML}</body></html>`;
}
function htmlToPlain(container){
  // Tab-separated, with title on first lines
  const out = [];
  const num = container.querySelector('.apa-table-num');
  const title = container.querySelector('.apa-table-title');
  if (num) out.push(num.textContent);
  if (title) out.push(title.textContent);
  out.push('');
  const table = container.querySelector('.apa-table');
  if (table){
    table.querySelectorAll('tr').forEach(tr => {
      const cells = [...tr.querySelectorAll('th,td')].map(c => c.textContent.trim().replace(/\s+/g,' '));
      out.push(cells.join('\t'));
    });
  }
  const note = container.querySelector('.apa-note');
  if (note){ out.push(''); out.push(note.textContent); }
  return out.join('\n');
}
function csvEscape(value){
  const s = String(value ?? '').replace(/\s+/g, ' ').trim();
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
function slugifyFilename(value, fallback){
  const slug = String(value || fallback || 'apa-table')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '')
    .slice(0, 72);
  return slug || fallback || 'apa-table';
}
function buildApaCsv(container){
  const rows = [];
  const num = container.querySelector('.apa-table-num');
  const title = container.querySelector('.apa-table-title');
  if (num) rows.push([num.textContent]);
  if (title) rows.push([title.textContent]);
  if (rows.length) rows.push([]);

  const table = container.querySelector('.apa-table');
  if (table){
    table.querySelectorAll('tr').forEach(tr => {
      rows.push([...tr.querySelectorAll('th,td')].map(c => c.textContent));
    });
  }

  const note = container.querySelector('.apa-note');
  if (note){
    rows.push([]);
    rows.push([note.textContent]);
  }
  return '\ufeff' + rows.map(row => row.map(csvEscape).join(',')).join('\r\n');
}
function downloadApaTableCsv(containerId){
  const container = document.getElementById(containerId);
  if (!container || !container.querySelector('.apa-table')){
    showToast('No table to download yet - fill in some data first', true);
    return;
  }
  const title = container.querySelector('.apa-table-title')?.textContent || containerId;
  const blob = new Blob([buildApaCsv(container)], { type:'text/csv;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `${slugifyFilename(title, containerId)}.csv`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
  showToast('✓ CSV downloaded - opens in Excel');
  if (typeof ReportBundle !== 'undefined' && ReportBundle.showKofiPrompt) ReportBundle.showKofiPrompt();
}
function enhanceApaToolbars(){
  document.querySelectorAll('[data-copy]').forEach(btn => {
    const outId = btn.dataset.copy;
    btn.classList.add('apa-action-btn');
    btn.childNodes.forEach(node => {
      if (node.nodeType === Node.TEXT_NODE && node.textContent.includes('Copy table')){
        node.textContent = ' Copy table';
      }
    });
    if (!document.querySelector(`[data-download="${outId}"]`)){
      const dl = document.createElement('button');
      dl.type = 'button';
      dl.className = 'btn apa-action-btn';
      dl.dataset.download = outId;
      dl.title = 'Download CSV for Excel';
      dl.innerHTML = `
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><path d="M7 10l5 5 5-5"></path><path d="M12 15V3"></path></svg>
        Download CSV
      `;
      btn.insertAdjacentElement('afterend', dl);
    }
  });
}
document.addEventListener('click', e => {
  const btn = e.target.closest('[data-copy]');
  if (btn) copyApaTable(btn.dataset.copy);
  const dl = e.target.closest('[data-download]');
  if (dl) downloadApaTableCsv(dl.dataset.download);
});

const examples = {
  'rci-basic': {name:'Example index score',sd:'15',r:'0.90',t1:'100',t2:'89'},
  'rci-practice': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89'},
  'rci-srb': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89'},
  'rci-crawford': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',n:'100',t1:'100',t2:'89'}
};
function renderConverter(){
  const type = document.getElementById('conv-type').value;
  const val = document.getElementById('conv-value').value;
  const out = document.getElementById('conv-output');
  if (val === '' || isNaN(val)){ out.style.display = 'none'; return; }
  out.style.display = 'block';
  const z = toZ(val, type);
  if (z == null){ out.style.display = 'none'; return; }
  // Sync slider and readout
  const slider = document.getElementById('conv-slider');
  if (slider) slider.value = Math.max(-3, Math.min(3, z)).toFixed(2);
  const pct = normCDF(z) * 100;
  const pctStr = pct < 0.1 ? '<0.1' : pct > 99.9 ? '>99.9' : pct.toFixed(1);
  const zEl = document.getElementById('conv-z-display');
  const pEl = document.getElementById('conv-pct-display');
  if (zEl) zEl.textContent = z.toFixed(2);
  if (pEl) pEl.textContent = pctStr;
  const convTypes = [
    { key:'z',          label:'Z-Score',       v: fmt(fromZ(z,'z'),2) },
    { key:'t',          label:'T-Score',        v: fmt(fromZ(z,'t'),1) },
    { key:'scaled',     label:'Scaled Score',   v: fmt(fromZ(z,'scaled'),1) },
    { key:'standard',   label:'Standard Score', v: fmt(fromZ(z,'standard'),1) },
    { key:'percentile', label:'Percentile',     v: fmtPct(fromZ(z,'percentile')) }
  ];
  document.getElementById('conv-grid').innerHTML = convTypes.map(c => {
    const active = c.key === type;
    return `<div class="conv-score-item${active ? ' conv-score-active' : ''}">
      <span class="conv-score-label">${c.label}</span>
      <span class="conv-score-value">${c.v}</span>
    </div>`;
  }).join('');
  const ss = fromZ(z, 'standard');
  const activeIdx = ssToDescIndex(ss);
  updateDescCarousel('conv-wechsler-block', activeIdx);
  updateDescCarousel('conv-aan-block', activeIdx);

  drawCurve(z, type);
}
function drawCurve(z, scoreType){
  const svg = document.getElementById('conv-curve');
  const W = 640, H = 408, padL = 94, padR = 48;
  const bottomPad = 124;
  const topPad = 76;
  const curveH = H - topPad - bottomPad;
  const xMin = -3.5, xMax = 3.5;
  const xScale = x => padL + (x - xMin) / (xMax - xMin) * (W - padL - padR);
  const yMax = 0.42;
  const base = H - bottomPad;
  const yScale = y => base - y / yMax * curveH;
  const zClamp = Math.max(xMin, Math.min(xMax, z));

  // Subtle band fills aligned to ±1 / ±2 SD
  function bandPath(z1, z2) {
    const za = Math.max(xMin, z1), zb = Math.min(xMax, z2);
    if (za >= zb) return '';
    let d = `M${xScale(za).toFixed(1)},${base} `;
    for (let i = 0; i <= 80; i++) {
      const x = za + (zb - za) * i / 80;
      d += `L${xScale(x).toFixed(1)},${yScale(normPDF(x)).toFixed(1)} `;
    }
    return d + `L${xScale(zb).toFixed(1)},${base} Z`;
  }
  const bands = [
    { z1:-4, z2:-2, fill:'rgba(156,61,42,0.09)'  },
    { z1:-2, z2:-1, fill:'rgba(195,95,40,0.07)'  },
    { z1:-1, z2: 0, fill:'rgba(190,155,45,0.06)' },
    { z1: 0, z2: 1, fill:'rgba(65,115,70,0.05)'  },
    { z1: 1, z2: 2, fill:'rgba(50,95,175,0.07)'  },
    { z1: 2, z2: 4, fill:'rgba(75,45,158,0.09)'  },
  ];
  const bandsSvg = bands.map(b => {
    const d = bandPath(b.z1, b.z2);
    return d ? `<path d="${d}" fill="${b.fill}" stroke="none"/>` : '';
  }).join('');

  // Left-tail shading
  let area = `M${xScale(xMin)},${base} `;
  for (let i = 0; i <= 200; i++){
    const x = xMin + (zClamp - xMin) * (i / 200);
    area += `L${xScale(x).toFixed(1)},${yScale(normPDF(x)).toFixed(1)} `;
  }
  area += `L${xScale(zClamp)},${base} Z`;

  // Vertical SD lines with labels above
  const sdConfig = [
    { z:-2, label:'−2 SD', bold:false },
    { z:-1, label:'−1 SD', bold:false },
    { z: 0, label:'Mean',  bold:true  },
    { z: 1, label:'+1 SD', bold:false },
    { z: 2, label:'+2 SD', bold:false },
  ];
  let sdLineSegs = '', sdLabelTexts = '';
  sdConfig.forEach(({ z: zv, label, bold }) => {
    const lx = xScale(zv);
    const curveY = yScale(normPDF(zv));
    sdLineSegs += `<line x1="${lx}" y1="${curveY.toFixed(1)}" x2="${lx}" y2="${base}"
      stroke="${bold ? '#555' : '#999'}" stroke-width="${bold ? 1.2 : 0.8}" opacity="0.55"/>`;
    sdLabelTexts += `<text x="${lx}" y="${(curveY - 9).toFixed(1)}"
      font-family="IBM Plex Sans" font-size="${bold ? 10 : 9}" fill="${bold ? '#444' : '#777'}"
      text-anchor="middle" font-weight="${bold ? '600' : '400'}">${label}</text>`;
  });
  // Short ticks at ±3
  [-3, 3].forEach(zv => {
    const lx = xScale(zv);
    sdLineSegs += `<line x1="${lx}" y1="${base}" x2="${lx}" y2="${base+4}" stroke="#B0A89E" stroke-width="0.6"/>`;
  });

  // Bell curve
  let path = '';
  for (let i = 0; i <= 260; i++){
    const x = xMin + (xMax - xMin) * (i / 260);
    path += (i===0?'M':'L') + xScale(x).toFixed(1) + ',' + yScale(normPDF(x)).toFixed(1) + ' ';
  }

  // Band percentage labels
  const bandPcts = [
    { zMid:-2.5, pct:'2.14%'  },
    { zMid:-1.5, pct:'13.59%' },
    { zMid:-0.5, pct:'34.13%' },
    { zMid: 0.5, pct:'34.13%' },
    { zMid: 1.5, pct:'13.59%' },
    { zMid: 2.5, pct:'2.14%'  },
  ];
  let pctLabels = '';
  bandPcts.forEach(({ zMid, pct: p }) => {
    const bx = xScale(zMid);
    const curveY = yScale(normPDF(zMid));
    const labelY = curveY + (base - curveY) * 0.45;
    pctLabels += `<text x="${bx.toFixed(1)}" y="${labelY.toFixed(1)}"
      font-family="IBM Plex Mono" font-size="8.5" fill="#888" text-anchor="middle" opacity="0.85">${p}</text>`;
  });

  // Score scale rows
  const scaleRows = [
    { key:'standard',     label:'Standard Score', vals:[55,70,85,100,115,130,145], fmt:v=>String(v) },
    { key:'scaled',       label:'Scaled Score',   vals:[1,4,7,10,13,16,19],        fmt:v=>String(v) },
    { key:'t',            label:'T-Score',        vals:[20,30,40,50,60,70,80],     fmt:v=>String(v) },
    { key:'z',            label:'z-Score',        vals:[-3,-2,-1,0,1,2,3],         fmt:v=>v===0?'0':(v>0?'+':'')+v },
    { key:'percentile',   label:'Percentile',     vals:[-3,-2,-1,0,1,2,3],
      fmt:v=>{const p=normCDF(v)*100; return p<0.5?'<1':p>99.5?'>99':String(Math.round(p));} },
    { key:'descriptor',   label:'Classification', vals:[-3,-2,-1,0,1,2,3], isDesc:true,
      fmt:v=>{const ss=v*15+100; if(ss>=130)return 'V.Superior'; if(ss>=120)return 'Superior'; if(ss>=110)return 'High Avg'; if(ss>=90)return 'Average'; if(ss>=80)return 'Low Avg'; if(ss>=70)return 'Borderline'; return 'Ext.Low';} },
  ];
  let rowsSvg = '';
  scaleRows.forEach((row, i) => {
    const rowY = base + 22 + i * 15;
    const active = !row.isDesc && row.key === scoreType;
    const col = row.isDesc ? '#909090' : (active ? '#9C3D2A' : '#A8A29E');
    const wt  = active ? '600' : '400';
    if (i > 0) rowsSvg += `<line x1="${padL}" y1="${base+22+i*15-9}" x2="${W-padR}" y2="${base+22+i*15-9}" stroke="#EAE4D8" stroke-width="0.5"/>`;
    rowsSvg += `<text x="${padL-6}" y="${rowY}" font-family="IBM Plex Mono" font-size="7.5" fill="${col}" text-anchor="end" font-weight="${wt}">${row.label}</text>`;
    row.vals.forEach((v, j) => {
      rowsSvg += `<text x="${xScale(j-3).toFixed(1)}" y="${rowY}" font-family="IBM Plex Mono" font-size="${row.isDesc?'8':'9'}" fill="${col}" text-anchor="middle" font-weight="${active?'500':'400'}">${row.fmt(v)}</text>`;
    });
  });

  // Position marker
  const yAtZ = yScale(normPDF(zClamp));
  const cx   = xScale(zClamp);
  const zLine  = `<line x1="${cx}" y1="${yAtZ}" x2="${cx}" y2="${base}" stroke="#9C3D2A" stroke-width="1.8" stroke-dasharray="4,3" opacity="0.85"/>`;
  const zDot   = `<circle cx="${cx}" cy="${yAtZ}" r="5.5" fill="#9C3D2A"/>`;
  const zInner = `<circle cx="${cx}" cy="${yAtZ}" r="2.2" fill="#fff"/>`;

  // Classification & percentile
  const pct = normCDF(z) * 100;
  const ss  = z * 15 + 100;
  let classification = 'Average';
  if      (ss >= 130) classification = 'Very Superior';
  else if (ss >= 120) classification = 'Superior';
  else if (ss >= 110) classification = 'High Average';
  else if (ss >= 90)  classification = 'Average';
  else if (ss >= 80)  classification = 'Low Average';
  else if (ss >= 70)  classification = 'Borderline';
  else                classification = 'Extremely Low';

  const pctRound = Math.round(pct);
  const pctOrd = (() => {
    if (pct < 0.5) return '<1st';
    if (pct > 99.5) return '>99th';
    const sfx = (pctRound>=11&&pctRound<=13)?'th':pctRound%10===1?'st':pctRound%10===2?'nd':pctRound%10===3?'rd':'th';
    return pctRound + sfx;
  })();

  // Callout box
  const boxW = 152, boxH = 52;
  let bx = cx - boxW / 2;
  bx = Math.max(padL, Math.min(W - padR - boxW, bx));
  const by = Math.max(4, yAtZ - boxH - 18);
  const stemX = Math.max(bx + 10, Math.min(bx + boxW - 10, cx));
  const callout = `
    <rect x="${bx}" y="${by}" width="${boxW}" height="${boxH}" rx="5"
      fill="#FAF7F1" stroke="#DDD6CC" stroke-width="1" filter="url(#cshadow)"/>
    <line x1="${stemX}" y1="${by+boxH}" x2="${cx}" y2="${yAtZ-7}" stroke="#CEC8BE" stroke-width="1"/>
    <text x="${bx+boxW/2}" y="${by+20}" font-family="Source Serif 4,serif"
      font-style="italic" font-size="14.5" font-weight="500" fill="#9C3D2A" text-anchor="middle">${classification}</text>
    <text x="${bx+boxW/2}" y="${by+38}" font-family="IBM Plex Mono"
      font-size="9" fill="#6B6B6B" text-anchor="middle">${pctOrd} percentile</text>`;

  // Two-tail annotation
  const belowN = pct < 1 ? '<1' : pct > 99 ? '>99' : String(pctRound);
  const aboveN = (100-pct) < 1 ? '<1' : (100-pct) > 99 ? '>99' : String(100-pctRound);
  let twoTail = '';
  if (cx > padL + 70)
    twoTail += `<text x="${cx-10}" y="${base-8}" font-family="IBM Plex Mono" font-size="8" fill="#BFB2A5" text-anchor="end">${belowN}% scored below ◂</text>`;
  if (cx < W - padR - 70)
    twoTail += `<text x="${cx+10}" y="${base-8}" font-family="IBM Plex Mono" font-size="8" fill="#BFB2A5" text-anchor="start">▸ ${aboveN}% scored above</text>`;

  const baseLine = `<line x1="${padL}" y1="${base}" x2="${W-padR}" y2="${base}" stroke="#2A2A2A" stroke-width="0.8"/>`;

  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = `
    <defs>
      <linearGradient id="tail-grad" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%"  stop-color="#9C3D2A" stop-opacity="0.28"/>
        <stop offset="100%" stop-color="#9C3D2A" stop-opacity="0.04"/>
      </linearGradient>
      <filter id="cshadow" x="-20%" y="-30%" width="140%" height="160%">
        <feDropShadow dx="0" dy="1" stdDeviation="3" flood-color="rgba(60,40,20,0.12)"/>
      </filter>
    </defs>
    ${bandsSvg}
    <path d="${area}" fill="url(#tail-grad)" stroke="none"/>
    ${sdLineSegs}
    <path d="${path}" fill="none" stroke="#2A2A2A" stroke-width="1.8"/>
    ${sdLabelTexts}
    ${baseLine}
    ${pctLabels}
    ${twoTail}
    ${rowsSvg}
    ${zLine}${zDot}${zInner}${callout}
  `;
}
function updateSliderTicks(type) {
  const scaleParams = {
    standard: { mean:100, sd:15, fmt: v => Math.round(v) },
    t:        { mean:50,  sd:10, fmt: v => Math.round(v) },
    scaled:   { mean:10,  sd:3,  fmt: v => Math.round(v) },
    z:        { mean:0,   sd:1,  fmt: v => (v > 0 ? '+' : '') + v.toFixed(1) }
  };
  const zVals = [-3, -2, -1, 0, 1, 2, 3];
  zVals.forEach((z, i) => {
    const el = document.getElementById('conv-tick-' + i);
    if (!el) return;
    let label;
    if (type === 'percentile') {
      const p = normCDF(z) * 100;
      label = p < 0.5 ? '<1%' : p > 99.5 ? '>99%' : Math.round(p) + '%';
    } else {
      const s = scaleParams[type] || scaleParams.standard;
      label = s.fmt(s.mean + z * s.sd);
    }
    el.textContent = label;
    el.className = z === 0 ? 'conv-slider-center' : '';
  });
}

document.getElementById('conv-type').addEventListener('change', function() {
  updateSliderTicks(this.value);
  renderConverter();
});
document.getElementById('conv-value').addEventListener('input', renderConverter);

// Slider - syncs back to the value input and re-renders
document.getElementById('conv-slider').addEventListener('input', function(){
  const z = parseFloat(this.value);
  const type = document.getElementById('conv-type').value;
  let displayVal;
  if (type === 'percentile') {
    displayVal = (normCDF(z) * 100).toFixed(1);
  } else {
    const scaleParams = {standard:{mean:100,sd:15},t:{mean:50,sd:10},scaled:{mean:10,sd:3},z:{mean:0,sd:1}};
    const s = scaleParams[type] || scaleParams.standard;
    const raw = s.mean + z * s.sd;
    displayVal = (type === 'z') ? raw.toFixed(2) : Math.round(raw * 10) / 10;
  }
  document.getElementById('conv-value').value = displayVal;
  renderConverter();
});

// Initialise slider ticks for default score type
updateSliderTicks(document.getElementById('conv-type').value);

// Click-and-drag on the bell curve to set z-score directly
(function setupCurveDrag(){
  const svg = document.getElementById('conv-curve');
  if (!svg) return;
  const W = 640, padL = 94, padR = 48;
  const xMin = -3.5, xMax = 3.5;
  let dragging = false;

  function clientToZ(clientX){
    const rect = svg.getBoundingClientRect();
    const localX = (clientX - rect.left) * (W / rect.width);
    const z = xMin + (localX - padL) / (W - padL - padR) * (xMax - xMin);
    return Math.max(-3, Math.min(3, z));
  }
  function applyZ(z){
    const type = document.getElementById('conv-type').value;
    const scaleParams = {standard:{mean:100,sd:15},t:{mean:50,sd:10},scaled:{mean:10,sd:3},z:{mean:0,sd:1}};
    let displayVal;
    if (type === 'percentile'){
      displayVal = (normCDF(z) * 100).toFixed(1);
    } else {
      const s = scaleParams[type] || scaleParams.standard;
      const raw = s.mean + z * s.sd;
      displayVal = (type === 'z') ? raw.toFixed(2) : Math.round(raw * 10) / 10;
    }
    document.getElementById('conv-value').value = displayVal;
    renderConverter();
  }

  svg.addEventListener('pointerdown', e => {
    dragging = true;
    try { svg.setPointerCapture(e.pointerId); } catch(_) {}
    applyZ(clientToZ(e.clientX));
    e.preventDefault();
  });
  svg.addEventListener('pointermove', e => {
    if (!dragging) return;
    applyZ(clientToZ(e.clientX));
  });
  ['pointerup','pointercancel'].forEach(ev =>
    svg.addEventListener(ev, e => {
      dragging = false;
      try { svg.releasePointerCapture(e.pointerId); } catch(_) {}
    })
  );
})();

// Slider tick marks
(function initSliderMarks(){
  const wrap = document.getElementById('conv-slider-marks');
  if (!wrap) return;
  let html = '';
  for (let v = -3; v <= 3; v += 0.5){
    const pct = ((v + 3) / 6) * 100;
    const major = Number.isInteger(v);
    html += `<div class="mark${major ? ' major' : ''}" style="left:${pct}%"></div>`;
  }
  wrap.innerHTML = html;
})();

/* ============================================================
   02 · BATTERY TABLE
   ============================================================ */
let batteryRows = [];
function batteryAddRow(initial){
  batteryRows.push(initial || { name:'', raw:'', score:'' });
  renderBattery();
}
function batteryRemove(i){ batteryRows.splice(i, 1); renderBattery(); }
function batteryRemoveGroup(group){
  batteryRows = batteryRows.filter(r => r.group !== group);
  renderBattery();
}
window.batteryRemove = batteryRemove;
window.batteryRemoveGroup = batteryRemoveGroup;

// Infer score type from a family name in the database.
// "Indices" / "Index Scores" / "IQ" → standard. Otherwise → scaled.
function inferScoreType(familyName){
  const n = (familyName || '').toLowerCase();
  if (n.includes('indic') || /\bindex\b/.test(n) || /\biq\b/.test(n)) return 'standard';
  return 'scaled';
}
// Strip age-band suffixes for compact group labels in editable and APA tables.
function stripAgeRange(name){
  if (!name) return name;
  return name
    .replace(/\s*·\s*Ages?\s+[\d–\-–]+\s*$/i, '')
    .replace(/\s*·\s*All\s+Ages\s*$/i, '')
    .replace(/\s*\(all ages\)\s*$/i, '');
}
function scoreTypeLabel(type){
  return {scaled:'Scaled Score', standard:'Standard Score', t:'T-Score', z:'Z-Score'}[type] || 'Score';
}
function updateBatteryScoreHeader(){
  const type = document.getElementById('bat-type')?.value;
  const label = scoreTypeLabel(type);
  document.querySelectorAll('#bat-table .bat-score-head').forEach(th => {
    th.innerHTML = `<span class="bat-score-head-label">${escapeHtml(label)}</span>`;
    th.title = `Manual rows use the selected default score type: ${label}. Auto-filled test families show their inferred score type in the group row.`;
  });
}
function rowScoreType(r){
  return r.scoreType || document.getElementById('bat-type').value;
}

function batteryPremorbidThresholdLabel(v){
  const n = parseFloat(v);
  return Number.isInteger(n) ? String(n) : String(n).replace(/0+$/,'').replace(/\.$/,'');
}
function syncBatteryPremorbidControls(){
  const enabled = document.getElementById('bat-prem-enable')?.checked;
  ['bat-prem-score','bat-prem-threshold'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.disabled = !enabled;
  });
}
function getBatteryPremorbidComparison(){
  const enabled = document.getElementById('bat-prem-enable')?.checked;
  const scoreEl = document.getElementById('bat-prem-score');
  const estimate = parseFloat(scoreEl?.value);
  if (!enabled || isNaN(estimate)) return null;
  return { estimate };
}
function batteryPremorbidStars(ss, prem){
  if (!prem || !Number.isFinite(ss)) return '';
  const diffSd = (prem.estimate - ss) / 15;
  if (diffSd >= 2) return '***';
  if (diffSd >= 1.5) return '**';
  if (diffSd >= 1) return '*';
  return '';
}
function batteryClassificationDetails(r, cls){
  const z = toZ(r.score, rowScoreType(r));
  if (z == null) return { text:'', html:'', className:'' };
  const ss = fromZ(z, 'standard');
  const desc = cls === 'wechsler' ? wechslerDesc(ss) : aanDesc(ss);
  const prem = getBatteryPremorbidComparison();
  const extreme = ss < 70;
  const stars = batteryPremorbidStars(ss, prem);
  const className = `${extreme ? ' bat-class-extreme' : ''}`.trim();
  const descriptor = `${escapeHtml(desc)}${stars ? `<span class="bat-expected-stars">${stars}</span>` : ''}`;
  return {
    text: `${desc}${stars}`,
    html: className ? `<span class="${className}">${descriptor}</span>` : descriptor,
    className,
    ss,
    prem,
    stars
  };
}
function setupBatteryContextualTabbing(tbody){
  if (!tbody) return;

  tbody.querySelectorAll('input[data-f]').forEach(el => {
    el.addEventListener('keydown', e => {
      if (e.key !== 'Tab' || e.altKey || e.ctrlKey || e.metaKey) return;

      const currentField = el.dataset.f;
      if (!currentField) return;

      const focusables = Array.from(tbody.querySelectorAll(`input[data-f="${currentField}"]`))
        .filter(input =>
          input.offsetParent !== null &&
          !input.disabled
        );

      const idx = focusables.indexOf(el);
      if (idx === -1 || focusables.length < 2) return;

      e.preventDefault();

      const nextIdx = e.shiftKey
        ? (idx - 1 + focusables.length) % focusables.length
        : (idx + 1) % focusables.length;

      focusables[nextIdx].focus();

      if (typeof focusables[nextIdx].select === 'function') {
        focusables[nextIdx].select();
      }
    });
  });
}
function renderBattery(){
  syncBatteryPremorbidControls();
  updateBatteryScoreHeader();
  const cls = document.getElementById('bat-class').value;
  const tbody = document.querySelector('#bat-table tbody');
  tbody.innerHTML = '';
  let lastGroup = null;
  batteryRows.forEach((r, i) => {
    // Inject a group header when group changes
    if (r.group && r.group !== lastGroup){
      const ghr = document.createElement('tr');
      ghr.className = 'group-header';
      const stLabel = scoreTypeLabel(r.scoreType || inferScoreType(r.group));
      ghr.innerHTML = `<td colspan="7">${escapeHtml(stripAgeRange(r.group))} <span class="type-badge">· ${stLabel}</span><button class="group-remove" data-rm-group="${escapeAttr(r.group)}" title="Remove group">×</button></td>`;
      tbody.appendChild(ghr);
      lastGroup = r.group;
    } else if (!r.group){
      lastGroup = null;
    }
    const rowType = rowScoreType(r);
    const z = toZ(r.score, rowType);
    const pct = z == null ? '' : fmtPct(normCDF(z) * 100);
    const details = batteryClassificationDetails(r, cls);
    const tr = document.createElement('tr');
    if (r.group) tr.className = 'in-group';
    tr.innerHTML = `
      <td class="row-num">${i+1}</td>
      <td><input type="text" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Subtest name"></td>
      <td><input type="number" step="any" data-r="${i}" data-f="raw" value="${escapeAttr(r.raw)}"></td>
      <td><input type="number" step="any" data-r="${i}" data-f="score" value="${escapeAttr(r.score)}"></td>
      <td class="computed">${pct}</td>
      <td class="computed ${details.className}">${details.html}</td>
      <td class="row-actions"><button onclick="batteryRemove(${i})" title="Remove">×</button></td>
    `;
    tbody.appendChild(tr);
  });
  // Wire group-remove buttons
  tbody.querySelectorAll('[data-rm-group]').forEach(b => {
    b.addEventListener('click', () => batteryRemoveGroup(b.dataset.rmGroup));
  });
  // In-place updates while typing
  tbody.querySelectorAll('input').forEach(inp => {
    inp.addEventListener('input', e => {
      const i = +e.target.dataset.r, f = e.target.dataset.f;
      batteryRows[i][f] = e.target.value;
      const tr = e.target.closest('tr');
      const rowType = rowScoreType(batteryRows[i]);
      const z = toZ(batteryRows[i].score, rowType);
      const cells = tr.querySelectorAll('.computed');
      cells[0].textContent = z == null ? '' : fmtPct(normCDF(z) * 100);
      const details = batteryClassificationDetails(batteryRows[i], cls);
      cells[1].className = `computed ${details.className}`.trim();
      cells[1].innerHTML = details.html;
      renderBatteryApa();
    });
  });
    setupBatteryContextualTabbing(tbody);
  renderBatteryApa();
}


const apaColumnState = {};
function getApaVisibleColumns(outId, columns){
  const allKeys = columns.map(c => c.key);
  const defaultKeys = columns.filter(c => c.defaultVisible !== false).map(c => c.key);
  if (!apaColumnState[outId]) apaColumnState[outId] = new Set(defaultKeys.length ? defaultKeys : allKeys);
  apaColumnState[outId] = new Set([...apaColumnState[outId]].filter(k => allKeys.includes(k)));
  if (apaColumnState[outId].size === 0) apaColumnState[outId] = new Set(allKeys);
  return columns.filter(c => apaColumnState[outId].has(c.key));
}
function updateApaColumnControls(outId, columns, renderFn){
  const btn = document.querySelector(`[data-copy="${outId}"]`);
  if (!btn) return;
  let controls = document.querySelector(`.apa-column-controls[data-for="${outId}"]`);
  if (!controls){
    const disclosure = document.createElement('details');
    disclosure.className = 'apa-column-disclosure';
    disclosure.innerHTML = '<summary>Toggle columns</summary>';
    controls = document.createElement('div');
    controls.className = 'apa-column-controls';
    controls.dataset.for = outId;
    disclosure.appendChild(controls);
    btn.parentNode.insertBefore(disclosure, btn);
  }
  const visible = new Set(getApaVisibleColumns(outId, columns).map(c => c.key));
  controls.innerHTML = columns.map(c => `
    <label><input type="checkbox" value="${escapeAttr(c.key)}" ${visible.has(c.key) ? 'checked' : ''}>${escapeHtml(c.label.replace(/<[^>]*>/g,''))}</label>
  `).join('');
  controls.querySelectorAll('input[type="checkbox"]').forEach(cb => {
    cb.addEventListener('change', () => {
      let next = new Set([...controls.querySelectorAll('input[type="checkbox"]:checked')].map(x => x.value));
      if (next.size === 0){ cb.checked = true; next.add(cb.value); }
      apaColumnState[outId] = next;
      renderFn();
    });
  });
}
function buildApaTableFromColumns(outId, columns, rows, groupLabelFn){
  const visible = getApaVisibleColumns(outId, columns);
  const header = `<tr>${visible.map(c => {
    const cls = `${c.num ? 'num ' : ''}col-${c.key}`.trim();
    return `<th class="${cls}">${c.label}</th>`;
  }).join('')}</tr>`;
  let body = '';
  let lastGroup = null;
  let inGroup = false;
  rows.forEach(r => {
    const group = groupLabelFn ? groupLabelFn(r) : '';
    if (group && group !== lastGroup){
      body += `<tr class="apa-group"><td colspan="${visible.length}">${escapeHtml(stripAgeRange(group))}</td></tr>`;
      lastGroup = group;
      inGroup = true;
    } else if (!group){
      lastGroup = null;
      inGroup = false;
    }
    const cls = inGroup ? ' class="apa-grouped-row"' : '';
    body += `<tr${cls}>${visible.map(c => {
      const tdCls = `${c.num ? 'num ' : ''}col-${c.key}`.trim();
      return `<td class="${tdCls}">${c.render(r)}</td>`;
    }).join('')}</tr>`;
  });
  return `<table class="apa-table"><thead>${header}</thead><tbody>${body}</tbody></table>`;
}

function renderBatteryApa(){
  const cls = document.getElementById('bat-class').value;
  const title = document.getElementById('bat-title').value || 'Test scores';
  const out = document.getElementById('bat-apa');
  const valid = batteryRows.filter(r => r.name);
  const completed = valid.filter(r => r.score !== '' && !isNaN(r.score));
  const types = new Set((completed.length ? completed : valid).map(r => rowScoreType(r)));
  const headerLabel = types.size === 1 ? scoreTypeLabel([...types][0]) : 'Score';
  const columns = [
    { key:'subtest', label:'Subtest', num:false, render:r => escapeHtml(r.name) },
    { key:'raw', label:'Raw Score', group:'Scores', num:true, render:r => escapeHtml(r.raw || '-') },
    { key:'score', label:headerLabel, group:'Scores', num:true, render:r => escapeHtml(r.score || '') },
    { key:'percentile', label:'Percentile', group:'Scores', num:true, render:r => { const z = toZ(r.score, rowScoreType(r)); return z == null ? '' : fmtPct(normCDF(z) * 100); }},
    { key:'classification', label:'Classification', group:'Interpretation', num:false, render:r => batteryClassificationDetails(r, cls).html }
  ];
  updateApaColumnControls('bat-apa', columns, renderBatteryApa);
  if (valid.length === 0){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Add or select at least one subtest to preview the APA table.</div>';
    return;
  }
  const prem = getBatteryPremorbidComparison();
  const premNote = prem ? ` Premorbid estimate markers: * = at least 1 SD below the premorbid estimate, ** = at least 1.5 SD below, *** = at least 2 SD below (premorbid estimate Standard Score = ${fmt(prem.estimate, 1)}).` : '';
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(title)}</div>
    ${buildApaTableFromColumns('bat-apa', columns, valid, r => r.group)}
    <div class="apa-note"><strong>Note.</strong> Classification follows ${cls === 'wechsler' ? 'Wechsler conventions' : 'Guilmette et al. (2020)'}.${types.size > 1 ? ' Subtest scores are reported in their native standardised metric (scaled-score subtests vs. standard-score indices).' : ''}${premNote}</div>
  `;
}

// Auto-fill: append subtest names from a family as a new group
function loadFamilyIntoBattery(family){
  const db = getMergedDB();
  if (!db[family]) return;
  // Avoid duplicating an already-loaded group
  if (batteryRows.some(r => r.group === family)){
    showToast(`${family} is already loaded`, true);
    return;
  }
  const inferredType = inferScoreType(family);
  const names = Object.keys(db[family]);
  names.forEach(name => {
    batteryRows.push({ name, raw:'', score:'', group:family, scoreType:inferredType });
  });
  renderBattery();
  // Toast suppressed - the working-report pill is the single feedback channel
  // for "things added to the report". (Old toast was a duplicate.)
}
function clearBattery(){
  batteryRows.length = 0;
  document.getElementById('bat-family-input').value = '';
  renderBattery();
}

// Wire up the battery family combobox (after DOM nodes exist)
function wireBatteryAutofill(){
  const inp = document.getElementById('bat-family-input');
  const list = document.getElementById('bat-family-list');
  if (!inp || !list) return;
  inp.addEventListener('focus', () => { list.classList.add('show'); filterFamilyListEl(list, ''); });
  inp.addEventListener('input', () => { list.classList.add('show'); filterFamilyListEl(list, inp.value); });
  inp.addEventListener('keydown', e => {
    if (e.key === 'Escape') list.classList.remove('show');
    if (e.key === 'Enter') {
      const add = list.querySelector('.combo-add:not([disabled])');
      if (add){ e.preventDefault(); add.click(); }
    }
  });
  inp.addEventListener('blur', () => setTimeout(() => {
    if (!list.matches(':hover')) list.classList.remove('show');
  }, 180));
}

function comboCustomTag(isCustom){
  return '';
}
function comboFooterHtml(){
  return '<div class="combo-footer"><span class="combo-count">0 selected</span><button class="btn btn-ghost combo-clear" type="button">Clear</button><button class="btn btn-primary combo-add" type="button" disabled>Add selected tests</button></div>';
}
function comboAgeBandNoteHtml(){
  return '<div class="combo-ageband-note"><span class="combo-ageband-note-icon">ℹ</span><span><strong>Specific age bands</strong> offer greater normative precision but rest on smaller samples, which reduces the stability of <em>r</em>. <strong>All Ages</strong> norms draw on larger <em>N</em>, yielding a more robust <em>r</em>, at the cost of age specificity.</span></div>';
}
function comboCheckboxItemHtml(f, isCustom, indented, groupKey, displayLabel){
  const cls = 'combo-item combo-check' + (indented ? ' combo-indented' : '');
  const label = displayLabel
    ? escapeHtml(displayLabel)
    : (indented ? ageBandLabel(f) : escapeHtml(f));
  const groupAttr = groupKey ? ` data-group="${escapeAttr(groupKey)}"` : '';
  return `<label class="${cls}" data-family="${escapeAttr(f)}"${groupAttr}><input type="checkbox" value="${escapeAttr(f)}"><span class="combo-check-text">${label}${comboCustomTag(isCustom)}</span></label>`;
}
function comboOptionsHtml(itemsHtml){
  return `<div class="combo-options">${itemsHtml}</div>`;
}
function updateComboSelectionState(list){
  const selected = list.querySelectorAll('.combo-check input:checked').length;
  const count = list.querySelector('.combo-count');
  const add = list.querySelector('.combo-add');
  if (count) count.textContent = `${selected} selected`;
  if (add) add.disabled = selected === 0;
}
function selectedComboFamilies(list){
  return Array.from(list.querySelectorAll('.combo-check input:checked')).map(cb => cb.value);
}
function clearComboSelections(list){
  list.querySelectorAll('.combo-check input:checked').forEach(cb => { cb.checked = false; });
  updateComboSelectionState(list);
}
function wireMultiSelectFamilyList(list, onAdd){
  list.classList.add('is-multiselect');
  if (list.dataset.multiselectReady === 'true'){
    list._comboOnAdd = onAdd;
    updateComboSelectionState(list);
    return;
  }
  list.dataset.multiselectReady = 'true';
  list._comboOnAdd = onAdd;

  // Keep footer clicks from blurring/closing the dropdown before the action runs.
  // Checkbox rows are deliberately left alone so the browser can handle native
  // checkbox and label-click behaviour reliably.
  list.addEventListener('mousedown', e => {
    if (e.target.closest('.combo-footer')) e.preventDefault();
  });

  // Update the selected count/Add button after native checkbox toggling.
  list.addEventListener('change', e => {
    if (e.target.matches('.combo-check input[type="checkbox"]')) {
      updateComboSelectionState(list);
    }
  });

  list.addEventListener('click', e => {
    const clear = e.target.closest('.combo-clear');
    if (clear){
      e.preventDefault();
      clearComboSelections(list);
      return;
    }

    const add = e.target.closest('.combo-add');
    if (add){
      e.preventDefault();
      const families = selectedComboFamilies(list);
      if (!families.length) return;
      if (typeof list._comboOnAdd === 'function') list._comboOnAdd(families);
      clearComboSelections(list);
      // Keep the dropdown open so the user can keep adding tests without
      // re-clicking the input. They can dismiss with Esc, blur, or by
      // clicking outside (the existing global outside-click handler at
      // line ~2132 still closes it when appropriate).
      list.classList.add('show');
      return;
    }
  });

  updateComboSelectionState(list);
}
function rebuildBatteryFamilyList(){
  const list = document.getElementById('bat-family-list');
  if (!list) return;
  const db = getMergedDB();
  const families = Object.keys(db).sort();
  // Battery page: collapse age bands to a single entry per family, no
  // age-band note - norms don't affect the resulting table here.
  list.innerHTML = comboFooterHtml() + buildFamilyListHtml(families, { flat: true });
  wireMultiSelectFamilyList(list, families => {
    families.forEach(loadFamilyIntoBattery);
    const inp = document.getElementById('bat-family-input');
    if (inp){ inp.value = ''; inp.focus(); }
  });
}

document.getElementById('bat-add').addEventListener('click', () => batteryAddRow());
document.getElementById('bat-type').addEventListener('change', renderBattery);
document.getElementById('bat-class').addEventListener('change', renderBattery);
document.getElementById('bat-title').addEventListener('input', renderBatteryApa);
document.getElementById('bat-prem-enable').addEventListener('change', renderBattery);
document.getElementById('bat-prem-score').addEventListener('input', renderBattery);
document.getElementById('bat-prem-threshold').addEventListener('change', renderBattery);
document.getElementById('bat-clear').addEventListener('click', clearBattery);

/* ============================================================
   03 · SDI
   ============================================================ */
let sdiRows = [];
const sdiLabelState = { d1:'Test', d2:'Retest' };
function sdiAddRow(initial){ sdiRows.push(initial || { name:'', t1:'', t2:'', sd:'' }); renderSdi(); }
function sdiRemove(i){ sdiRows.splice(i, 1); renderSdi(); }
function sdiRemoveGroup(group){
  sdiRows = sdiRows.filter(r => r.group !== group);
  renderSdi();
}
window.sdiRemoveGroup = sdiRemoveGroup;
function sdiSdUnit(type){ return {scaled:3, standard:15, t:10, z:1}[type]; }
function sdiMode(){ return document.getElementById('sdi-mode').value; }
function sdiCvHit(change, cv){
  if (cv === 0.90) return Math.abs(change) >= 1.645;
  if (cv === 0.95) return Math.abs(change) >= 1.96;
  return Math.abs(change) >= cv;
}
function sdiDateLabel(which){
  const fallback = which === 'd1' ? 'Test' : 'Retest';
  return (sdiLabelState[which] || '').trim() || fallback;
}
function sdiComputeChange(r){
  if (r.t1 === '' || r.t2 === '' || isNaN(r.t1) || isNaN(r.t2)) return null;
  if (sdiMode() === 'raw'){
    if (r.sd === '' || isNaN(r.sd) || parseFloat(r.sd) <= 0) return null;
    return (parseFloat(r.t2) - parseFloat(r.t1)) / parseFloat(r.sd);
  }
  return (parseFloat(r.t2) - parseFloat(r.t1)) / sdiSdUnit(document.getElementById('sdi-type').value);
}
function renderSdiHead(){
  const raw = sdiMode() === 'raw';
  const d1 = sdiDateLabel('d1');
  const d2 = sdiDateLabel('d2');
  document.getElementById('sdi-table-head').innerHTML = `
    <tr class="table-group-row">
      <th colspan="2"></th>
      <th colspan="${raw ? 3 : 2}">Scores</th>
      <th colspan="3">Results</th>
      <th></th>
    </tr>
    <tr>
      <th class="row-num">#</th>
      <th style="min-width:200px">Subtest</th>
      <th style="width:100px"><input class="change-date-head-input" id="sdi-d1-head" type="text" value="${escapeAttr(d1)}" aria-label="Date 1 column label"></th>
      <th style="width:100px"><input class="change-date-head-input" id="sdi-d2-head" type="text" value="${escapeAttr(d2)}" aria-label="Date 2 column label"></th>
      ${raw ? '<th style="width:100px;text-align:right">SD</th>' : ''}
      <th style="width:90px;text-align:right">SD Δ</th>
      <th style="width:90px;text-align:right"><i>p</i></th>
      <th style="width:160px">Significance</th>
      <th class="row-actions"></th>
    </tr>`;
  document.getElementById('sdi-d1-head')?.addEventListener('input', e => {
    sdiLabelState.d1 = e.target.value;
    renderSdiApa();
  });
  document.getElementById('sdi-d2-head')?.addEventListener('input', e => {
    sdiLabelState.d2 = e.target.value;
    renderSdiApa();
  });
}function setupSdiContextualTabbing(tbody){
  if (!tbody) return;

  function fieldsFor(field){
    const raw = sdiMode() === 'raw';

    if (field === 'name') {
      return ['name'];
    }

    if (field === 't1' || field === 't2' || field === 'sd') {
      return raw ? ['t1', 't2', 'sd'] : ['t1', 't2'];
    }

    return [];
  }

  tbody.querySelectorAll('input[data-f]').forEach(el => {
    el.addEventListener('keydown', e => {
      if (e.key !== 'Tab' || e.altKey || e.ctrlKey || e.metaKey) return;

      const allowed = fieldsFor(el.dataset.f);
      if (!allowed.length) return;

      const focusables = Array.from(tbody.querySelectorAll('input[data-f]'))
        .filter(input =>
          allowed.includes(input.dataset.f) &&
          input.offsetParent !== null &&
          !input.disabled
        );

      const idx = focusables.indexOf(el);
      if (idx === -1 || focusables.length < 2) return;

      e.preventDefault();

      const nextIdx = e.shiftKey
        ? (idx - 1 + focusables.length) % focusables.length
        : (idx + 1) % focusables.length;

      focusables[nextIdx].focus();

      if (typeof focusables[nextIdx].select === 'function') {
        focusables[nextIdx].select();
      }
    });
  });
}
function renderSdi(){
  const raw = sdiMode() === 'raw';
  const cv = parseFloat(document.getElementById('sdi-cv').value);
  document.getElementById('sdi-type-field').classList.toggle('is-hidden', raw);
  document.getElementById('sdi-raw-help').style.display = raw ? 'block' : 'none';
  renderSdiHead();
  const tbody = document.querySelector('#sdi-table tbody');
  tbody.innerHTML = '';
  let lastGroup = null;
  sdiRows.forEach((r, i) => {
    if (r.group && r.group !== lastGroup){
      const ghr = document.createElement('tr');
      ghr.className = 'group-header';
      const colspan = raw ? 9 : 8;
      ghr.innerHTML = `<td colspan="${colspan}">${escapeHtml(stripAgeRange(r.group))}<button class="group-remove" data-rm-sdi-group="${escapeAttr(r.group)}" title="Remove group">×</button></td>`;
      tbody.appendChild(ghr);
      lastGroup = r.group;
    } else if (!r.group){
      lastGroup = null;
    }
    let sdChange = '', pStr = '', sigStr = '', sigCls = '';
    const change = sdiComputeChange(r);
    if (change === null){
      sigStr = sdiProblem(r);
      sigCls = sigStr === 'Awaiting values' ? 'status-awaiting' : 'status-check';
    }
    if (change !== null){
      sdChange = fmt(change, 2);
      const p = 2 * (1 - normCDF(Math.abs(change)));
      pStr = fmtP(p);
      const sig = sdiCvHit(change, cv);
      sigStr = sig ? 'Significant change' : 'No significant change';
      sigCls = sig ? 'sig-yes' : 'sig-no';
    }
    const tr = document.createElement('tr');
    if (r.group) tr.className = 'in-group';
    if (change === null && hasAnyRowValue(r)) tr.classList.add('row-check');
    else if (change === null) tr.classList.add('row-awaiting');
    tr.innerHTML = `
      <td class="row-num">${i+1}</td>
      <td><input type="text" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Subtest name"></td>
      <td><input type="number" step="any" data-r="${i}" data-f="t1" value="${escapeAttr(r.t1)}"></td>
      <td><input type="number" step="any" data-r="${i}" data-f="t2" value="${escapeAttr(r.t2)}"></td>
      ${raw ? `<td><input type="number" min="0" step="any" data-r="${i}" data-f="sd" value="${escapeAttr(r.sd || '')}" placeholder="SD"></td>` : ''}
      <td class="computed">${sdChange}</td>
      <td class="computed">${pStr}</td>
      <td class="computed ${sigCls}">${sigStr}</td>
      <td class="row-actions"><button onclick="sdiRemove(${i})">×</button></td>
    `;
    tbody.appendChild(tr);
  });
  tbody.querySelectorAll('[data-rm-sdi-group]').forEach(b => {
    b.addEventListener('click', () => sdiRemoveGroup(b.dataset.rmSdiGroup));
  });
  tbody.querySelectorAll('input').forEach(inp => {
    inp.addEventListener('input', e => {
      const i = +e.target.dataset.r, f = e.target.dataset.f;
      sdiRows[i][f] = e.target.value;
      updateSdiRow(i, e.target.closest('tr'));
      renderSdiApa();
    });
  });
  setupSdiContextualTabbing(tbody);
  renderSdiApa();
  addColumnTitles();
}

function updateSdiRow(i, tr){
  const r = sdiRows[i];
  const cv = parseFloat(document.getElementById('sdi-cv').value);
  const cells = tr.querySelectorAll('.computed');
  cells.forEach(c => { c.className = 'computed'; });
  const change = sdiComputeChange(r);
  tr.classList.remove('row-awaiting','row-check');
  if (change === null){
    const status = sdiProblem(r);
    cells.forEach(c => c.textContent = '');
    setOutcomeStatus(cells[cells.length - 1], status, status === 'Awaiting values' ? 'awaiting' : 'check');
    tr.classList.add(status === 'Awaiting values' && !hasAnyRowValue(r) ? 'row-awaiting' : 'row-check');
    return;
  }
  if (change !== null){
    const p = 2 * (1 - normCDF(Math.abs(change)));
    const sig = sdiCvHit(change, cv);
    cells[0].textContent = fmt(change, 2);
    cells[1].textContent = fmtP(p);
    cells[2].textContent = sig ? 'Significant change' : 'No significant change';
    cells[2].classList.add(sig ? 'sig-yes' : 'sig-no');
  } else {
    cells.forEach(c => c.textContent = '');
  }
}
function renderSdiApa(){
  const raw = sdiMode() === 'raw';
  const cv = parseFloat(document.getElementById('sdi-cv').value);
  const title = 'Test-retest comparison';
  const d1 = sdiDateLabel('d1');
  const d2 = sdiDateLabel('d2');
  const out = document.getElementById('sdi-apa');
  const named = sdiRows.filter(r => r.name);
  if (named.length === 0){ out.innerHTML = '<div class="apa-empty"><strong>APA-formatted output</strong>Add or select at least one test to preview the report-ready table.</div>'; return; }
  const cvDesc = cv === 0.90 ? '90% critical value (1.645)' : cv === 0.95 ? '95% critical value (1.96)' : `${cv} SD threshold`;
  const colCount = raw ? 7 : 6;
  // Insert apa-group separator rows when the test family changes - gives the
  // table visible grouping AND lets the working-report pill detect the most
  // recently added family (it reads the last `tr.apa-group`).
  let lastGroup = null;
  const rows = named.map(r => {
    const group = r.group || null;
    let prefix = '';
    if (group && group !== lastGroup){
      prefix = `<tr class="apa-group"><td colspan="${colCount}">${escapeHtml(stripAgeRange(group))}</td></tr>`;
      lastGroup = group;
    } else if (!group){
      lastGroup = null;
    }
    const inGroup = !!group;
    const trCls = inGroup ? ' class="apa-grouped-row"' : '';
    const change = sdiComputeChange(r);
    if (change === null){
      return prefix + `<tr${trCls}><td>${escapeHtml(r.name)}</td><td class="num">${escapeHtml(r.t1 || '')}</td><td class="num">${escapeHtml(r.t2 || '')}</td>${raw ? `<td class="num">${escapeHtml(r.sd || '')}</td>` : ''}<td class="num"></td><td class="num"></td><td></td></tr>`;
    }
    const p = 2 * (1 - normCDF(Math.abs(change)));
    const sig = sdiCvHit(change, cv);
    return prefix + `<tr${trCls}><td>${escapeHtml(r.name)}</td><td class="num">${escapeHtml(r.t1)}</td><td class="num">${escapeHtml(r.t2)}</td>${raw ? `<td class="num">${escapeHtml(r.sd)}</td>` : ''}<td class="num">${fmt(change, 2)}</td><td class="num">${fmtP(p)}</td><td>${sig ? 'Significant' : 'Not Significant'}</td></tr>`;
  }).join('');
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(title)}</div>
    <table class="apa-table">
      <thead>
        <tr><th>Subtest</th><th class="num">${escapeHtml(raw ? `${d1} raw` : d1)}</th><th class="num">${escapeHtml(raw ? `${d2} raw` : d2)}</th>${raw ? '<th class="num">SD</th>' : ''}<th class="num">SD Δ</th><th class="num"><i>p</i></th><th>Significance</th></tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>
    <div class="apa-note"><strong>Note.</strong> SD Δ = standard-deviation change between assessments${raw ? ', calculated as (retest raw score − test raw score) ÷ SD' : ''}. <i>p</i>-values are two-tailed assuming normality. Significance threshold = ${cvDesc}.</div>
  `;
}
function sdiNormSd(p){
  const sd = p && p.sd1 != null ? Number(p.sd1) : NaN;
  return Number.isFinite(sd) && sd > 0 ? String(sd) : '';
}
function loadFamilyIntoSdi(family){
  const db = getMergedDB();
  if (!db[family]) return;
  if (sdiRows.some(r => r.group === family)){
    showToast(`${family} is already loaded`, true);
    return;
  }
  const raw = sdiMode() === 'raw';
  const subtests = Object.entries(db[family]);
  subtests.forEach(([name, p]) => {
    sdiRows.push({ name, t1:'', t2:'', sd: raw ? sdiNormSd(p) : '', group:family });
  });
  renderSdi();
  // Toast suppressed - the working-report pill is the single feedback channel
}
function clearSdi(){
  sdiRows = [];
  const inp = document.getElementById('sdi-family-input');
  if (inp) inp.value = '';
  renderSdi();
}
function rebuildSdiFamilyList(){
  const list = document.getElementById('sdi-family-list');
  if (!list) return;
  const db = getMergedDB();
  const families = Object.keys(db).sort();
  list.innerHTML = comboFooterHtml() + buildFamilyListHtml(families);
  wireMultiSelectFamilyList(list, families => {
    families.forEach(loadFamilyIntoSdi);
    const inp = document.getElementById('sdi-family-input');
    if (inp){ inp.value = ''; inp.focus(); }
  });
}
function wireSdiAutofill(){
  const inp = document.getElementById('sdi-family-input');
  const list = document.getElementById('sdi-family-list');
  if (!inp || !list) return;
  inp.addEventListener('focus', () => { list.classList.add('show'); filterFamilyListEl(list, ''); });
  inp.addEventListener('input', () => { list.classList.add('show'); filterFamilyListEl(list, inp.value); });
  inp.addEventListener('keydown', e => {
    if (e.key === 'Escape') list.classList.remove('show');
    if (e.key === 'Enter') {
      const add = list.querySelector('.combo-add:not([disabled])');
      if (add){ e.preventDefault(); add.click(); }
    }
  });
  inp.addEventListener('blur', () => setTimeout(() => {
    if (!list.matches(':hover')) list.classList.remove('show');
  }, 180));
}
document.getElementById('sdi-add').addEventListener('click', () => sdiAddRow());
document.getElementById('sdi-mode').addEventListener('change', renderSdi);
document.getElementById('sdi-cv').addEventListener('change', renderSdi);
document.getElementById('sdi-type').addEventListener('change', renderSdi);
const sdiTitleEl = document.getElementById('sdi-title');
if (sdiTitleEl) sdiTitleEl.addEventListener('input', renderSdiApa);
document.getElementById('sdi-clear').addEventListener('click', clearSdi);

/* ============================================================
   04-07 · RCI CALCULATORS (Basic / Practice / SRB / Crawford)
   ============================================================ */
// Defaults for the corrected-r toggle:
//   - Simple RCI / Practice-Adjusted: ON  (corrected r is the better reliability estimate)
//   - SRB / Crawford:                  OFF (raw r matches the regression equations as published)
const rciState = {
  'rci-basic':    { rows:[], cv:0.95, useCorrectedR:true,  d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-practice': { rows:[], cv:0.95, useCorrectedR:true,  d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-srb':      { rows:[], cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-crawford': { rows:[], cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' }
};
function newRciRow(method){
  if (method === 'rci-basic') return { name:'', sd:'', r:'', t1:'', t2:'' };
  if (method === 'rci-crawford') return { name:'', m1:'', sd1:'', m2:'', sd2:'', r:'', n:'', t1:'', t2:'' };
  return { name:'', m1:'', sd1:'', m2:'', sd2:'', r:'', t1:'', t2:'' };
}
function rciAddRow(method){ rciState[method].rows.push(newRciRow(method)); renderRci(method); }
function rciRemove(method, i){ rciState[method].rows.splice(i, 1); renderRci(method); }
function rciRemoveGroup(method, group){
  rciState[method].rows = rciState[method].rows.filter(r => r.group !== group);
  renderRci(method);
}
window.rciRemove = rciRemove;
window.rciRemoveGroup = rciRemoveGroup;

function calcBasicRow(r, method){
  method = method || 'rci-basic';
  const sd = parseFloat(r.sd);
  const rEff = rciEffectiveR(method, r);
  const rel = rEff.value;
  const t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);
  if (isNaN(sd) || isNaN(rel) || isNaN(t1) || isNaN(t2) || rel < 0 || rel >= 1 || sd <= 0) return null;
  const sem = sd * Math.sqrt(1 - rel);
  const sed = Math.sqrt(2 * sem * sem);
  const rci = (t2 - t1) / sed;
  const p = 2 * (1 - normCDF(Math.abs(rci)));
  return { sem, sed, rci, p, usedR: rel, usedCorrected: rEff.fromCorrected, rFallback: !!rEff.fallbackBecauseMissing };
}
function calcPracticeRow(r, method){
  method = method || 'rci-practice';
  const m1 = parseFloat(r.m1), sd1 = parseFloat(r.sd1), m2 = parseFloat(r.m2), sd2 = parseFloat(r.sd2);
  const rEff = rciEffectiveR(method, r);
  const rel = rEff.value;
  const t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);
  if ([m1,sd1,m2,sd2,rel,t1,t2].some(isNaN) || rel < 0 || rel >= 1 || sd1 <= 0 || sd2 <= 0) return null;
  const sem1 = sd1 * Math.sqrt(1 - rel);
  const sem2 = sd2 * Math.sqrt(1 - rel);
  const sdiff = Math.sqrt(sem1*sem1 + sem2*sem2);
  const rci = ((t2 - t1) - (m2 - m1)) / sdiff;
  const p = 2 * (1 - normCDF(Math.abs(rci)));
  return { sem1, sem2, sdiff, rci, p, usedR: rel, usedCorrected: rEff.fromCorrected, rFallback: !!rEff.fallbackBecauseMissing };
}
/* Session-level "use corrected r" toggle. Defaults differ by method:
   ON for Simple RCI + Practice-Adjusted, OFF for SRB + Crawford. The toggle
   is stored on rciState[method].useCorrectedR, falling back to true. */
function rciEffectiveR(method, row){
  const wantCorrected = !rciState[method] || rciState[method].useCorrectedR !== false;
  const rCorr = parseFloat(row.rCorrected);
  const r     = parseFloat(row.r);
  if (wantCorrected && Number.isFinite(rCorr)) return { value: rCorr, fromCorrected: true };
  // Fallback to plain r (covers: toggle off, OR row has no corrected r available)
  return { value: r, fromCorrected: false, fallbackBecauseMissing: wantCorrected && !Number.isFinite(rCorr) };
}
function calcSrbRow(r, method){
  method = method || 'rci-srb';
  const m1 = parseFloat(r.m1), sd1 = parseFloat(r.sd1), m2 = parseFloat(r.m2), sd2 = parseFloat(r.sd2);
  const rEff = rciEffectiveR(method, r);
  const rel = rEff.value;
  const t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);
  if ([m1,sd1,m2,sd2,rel,t1,t2].some(isNaN) || rel < 0 || rel >= 1 || sd1 <= 0 || sd2 <= 0) return null;
  const slope = rel * (sd2 / sd1);
  const intercept = m2 - slope * m1;
  const predicted = intercept + slope * t1;
  const see = sd2 * Math.sqrt(1 - rel*rel);
  const rci = (t2 - predicted) / see;
  const p = 2 * (1 - normCDF(Math.abs(rci)));
  return { slope, intercept, predicted, see, rci, p, usedR: rel, usedCorrected: rEff.fromCorrected, rFallback: !!rEff.fallbackBecauseMissing };
}
function calcCrawfordRow(r, method){
  method = method || 'rci-crawford';
  const m1 = parseFloat(r.m1), sd1 = parseFloat(r.sd1), m2 = parseFloat(r.m2), sd2 = parseFloat(r.sd2);
  const rEff = rciEffectiveR(method, r);
  const rel = rEff.value;
  const n = parseFloat(r.n), t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);
  if ([m1,sd1,m2,sd2,rel,n,t1,t2].some(isNaN) || rel < 0 || rel >= 1 || sd1 <= 0 || sd2 <= 0 || n < 3) return null;
  const slope = rel * (sd2 / sd1);
  const intercept = m2 - slope * m1;
  const predicted = intercept + slope * t1;
  const see = sd2 * Math.sqrt(1 - rel*rel);
  // Crawford & Garthwaite (2007) sample-size-adjusted standard error of prediction
  const sePred = see * Math.sqrt(1 + 1/n + Math.pow(t1 - m1, 2) / ((n - 1) * sd1 * sd1));
  const df = n - 2;
  const tStat = (t2 - predicted) / sePred;
  const p = 2 * (1 - tCDF(Math.abs(tStat), df));
  return { slope, intercept, predicted, see, sePred, df, rci: tStat, p, usedR: rel, usedCorrected: rEff.fromCorrected, rFallback: !!rEff.fallbackBecauseMissing };
}
function rciOutcome(rci, cv, df){
  // For t-distributed RCIs (Crawford method), use t-quantile with df; otherwise z
  let crit;
  if (df != null && isFinite(df) && df > 0){
    crit = tInv(cv === 0.95 ? 0.975 : 0.95, df);
  } else {
    crit = cv === 0.95 ? 1.96 : 1.645;
  }
  if (Math.abs(rci) < crit) return { label:'No reliable change', cls:'sig-no' };
  return rci > 0
    ? { label:'Reliable improvement', cls:'sig-improve' }
    : { label:'Reliable decline', cls:'sig-decline' };
}



// Make keyboard entry match the user's current task.
// If the cursor is in Date 1/Date 2, Tab stays within the patient-score cells.
// If the cursor is in the test/normative fields, Tab stays within those fields.
function setupRciContextualTabbing(tbody){
  if (!tbody) return;
  const fieldsFor = current => (current === 't1' || current === 't2')
    ? ['t1','t2']
    : ['name','sd','r','m1','sd1','m2','sd2','n'];

  tbody.querySelectorAll('input[data-f]').forEach(el => {
    const field = el.dataset.f;
    el.tabIndex = 0;
    if (field === 't1' || field === 't2'){
      el.classList.add('patient-score-tab-target');
      el.title = 'Fast entry: Tab moves between Date 1 and Date 2 score cells.';
      el.setAttribute('aria-label', field === 't1' ? 'Date 1 score' : 'Date 2 score');
    } else {
      el.classList.add('norm-tab-target');
      el.title = 'Norm entry: Tab moves through the test and normative fields.';
    }

    el.addEventListener('keydown', e => {
      if (e.key !== 'Tab' || e.altKey || e.ctrlKey || e.metaKey) return;
      const allowed = fieldsFor(field);
      const focusables = Array.from(tbody.querySelectorAll('input[data-f]'))
        .filter(input => allowed.includes(input.dataset.f) && input.offsetParent !== null && !input.disabled);
      const idx = focusables.indexOf(el);
      if (idx === -1 || focusables.length < 2) return;
      e.preventDefault();
      const nextIdx = e.shiftKey
        ? (idx - 1 + focusables.length) % focusables.length
        : (idx + 1) % focusables.length;
      focusables[nextIdx].focus();
      if (typeof focusables[nextIdx].select === 'function') focusables[nextIdx].select();
    });
  });
}

function renderRciDateHeader(method){
  const st = rciState[method];
  const table = document.getElementById(`${method}-table`);
  if (!st || !table) return;

  const headerHtml = (which, value) => `
    <input
      class="change-date-head-input"
      type="text"
      data-date-head="${which}"
      data-target="${method}"
      value="${escapeAttr(value)}"
      aria-label="${which === 'd1' ? 'Date 1' : 'Date 2'} column label"
    >`;

  table.querySelectorAll('.d1-head').forEach(th => { th.innerHTML = headerHtml('d1', st.d1); });
  table.querySelectorAll('.d2-head').forEach(th => { th.innerHTML = headerHtml('d2', st.d2); });
  table.querySelectorAll('.change-date-head-input').forEach(input => {
    input.addEventListener('input', e => {
      const target = e.target.dataset.target;
      const which = e.target.dataset.dateHead;
      if (!rciState[target] || !which) return;
      rciState[target][which] = e.target.value;
      renderRciApa(target);
    });
  });
}

function renderRci(method){
  const st = rciState[method];
  const cv = st.cv;
  const tbody = document.querySelector(`#${method}-table tbody`);
  tbody.innerHTML = '';
  // Update header labels
  renderRciDateHeader(method);

  const colCount = document.querySelector(`#${method}-table thead tr:last-child`)?.children.length || 10;
  let lastGroup = null;
  st.rows.forEach((r, i) => {
    if (r.group && r.group !== lastGroup){
      const ghr = document.createElement('tr');
      ghr.className = 'group-header';
      ghr.innerHTML = `<td colspan="${colCount}">${escapeHtml(stripAgeRange(r.group))}<button class="group-remove" data-rm-rci-group="${escapeAttr(r.group)}" data-method="${method}" title="Remove group">×</button></td>`;
      tbody.appendChild(ghr);
      lastGroup = r.group;
    } else if (!r.group){
      lastGroup = null;
    }
    let calc = null;
    if (method === 'rci-basic') calc = calcBasicRow(r, method);
    else if (method === 'rci-practice') calc = calcPracticeRow(r, method);
    else if (method === 'rci-srb') calc = calcSrbRow(r, method);
    else if (method === 'rci-crawford') calc = calcCrawfordRow(r, method);

    const rciStr = calc ? fmt(calc.rci, 2) : '';
    const pStr = calc ? fmtP(calc.p) : '';
    const status = !calc ? numericProblem(r, method) : '';
    const outcome = calc ? rciOutcome(calc.rci, cv, calc.df) : { label: status, cls: (status === 'Awaiting values' ? 'status-awaiting' : 'status-check') };
    const predStr = calc ? fmt(calc.predicted, 2) : '';

    const tr = document.createElement('tr');
    if (r.group) tr.className = 'in-group';
    if (!calc && hasAnyRowValue(r)) tr.classList.add('row-check');
    else if (!calc) tr.classList.add('row-awaiting');
    if (method === 'rci-basic'){
      tr.innerHTML = `
        <td class="row-num">${i+1}</td>
        <td><input type="text" data-m="${method}" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Test name"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd" value="${escapeAttr(r.sd)}"></td>
        <td><input type="number" step="0.01" min="0" max="1" data-m="${method}" data-r="${i}" data-f="r" value="${escapeAttr(r.r)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t1" value="${escapeAttr(r.t1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t2" value="${escapeAttr(r.t2)}"></td>
        <td class="computed">${rciStr}</td>
        <td class="computed">${pStr}</td>
        <td class="computed ${outcome.cls}">${outcome.label}</td>
        <td class="row-actions"><button onclick="rciRemove('${method}',${i})">×</button></td>
      `;
    } else if (method === 'rci-practice'){
      tr.innerHTML = `
        <td class="row-num">${i+1}</td>
        <td><input type="text" data-m="${method}" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Test name"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m1" value="${escapeAttr(r.m1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd1" value="${escapeAttr(r.sd1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m2" value="${escapeAttr(r.m2)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd2" value="${escapeAttr(r.sd2)}"></td>
        <td><input type="number" step="0.01" min="0" max="1" data-m="${method}" data-r="${i}" data-f="r" value="${escapeAttr(r.r)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t1" value="${escapeAttr(r.t1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t2" value="${escapeAttr(r.t2)}"></td>
        <td class="computed">${rciStr}</td>
        <td class="computed">${pStr}</td>
        <td class="computed ${outcome.cls}">${outcome.label}</td>
        <td class="row-actions"><button onclick="rciRemove('${method}',${i})">×</button></td>
      `;
    } else if (method === 'rci-srb'){
      tr.innerHTML = `
        <td class="row-num">${i+1}</td>
        <td><input type="text" data-m="${method}" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Test name"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m1" value="${escapeAttr(r.m1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd1" value="${escapeAttr(r.sd1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m2" value="${escapeAttr(r.m2)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd2" value="${escapeAttr(r.sd2)}"></td>
        <td><input type="number" step="0.01" min="0" max="1" data-m="${method}" data-r="${i}" data-f="r" value="${escapeAttr(r.r)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t1" value="${escapeAttr(r.t1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t2" value="${escapeAttr(r.t2)}"></td>
        <td class="computed">${predStr}</td>
        <td class="computed">${rciStr}</td>
        <td class="computed">${pStr}</td>
        <td class="computed ${outcome.cls}">${outcome.label}</td>
        <td class="row-actions"><button onclick="rciRemove('${method}',${i})">×</button></td>
      `;
    } else if (method === 'rci-crawford'){
      tr.innerHTML = `
        <td class="row-num">${i+1}</td>
        <td><input type="text" data-m="${method}" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Test name"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m1" value="${escapeAttr(r.m1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd1" value="${escapeAttr(r.sd1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="m2" value="${escapeAttr(r.m2)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="sd2" value="${escapeAttr(r.sd2)}"></td>
        <td><input type="number" step="0.01" min="0" max="1" data-m="${method}" data-r="${i}" data-f="r" value="${escapeAttr(r.r)}"></td>
        <td><input type="number" step="1" min="3" data-m="${method}" data-r="${i}" data-f="n" value="${escapeAttr(r.n)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t1" value="${escapeAttr(r.t1)}"></td>
        <td><input type="number" step="any" data-m="${method}" data-r="${i}" data-f="t2" value="${escapeAttr(r.t2)}"></td>
        <td class="computed">${predStr}</td>
        <td class="computed">${rciStr}</td>
        <td class="computed">${pStr}</td>
        <td class="computed ${outcome.cls}">${outcome.label}</td>
        <td class="row-actions"><button onclick="rciRemove('${method}',${i})">×</button></td>
      `;
    }
    tbody.appendChild(tr);
  });
  tbody.querySelectorAll('[data-rm-rci-group]').forEach(b => {
    b.addEventListener('click', () => rciRemoveGroup(b.dataset.method, b.dataset.rmRciGroup));
  });
  tbody.querySelectorAll('input').forEach(inp => {
    inp.addEventListener('input', e => {
      const m = e.target.dataset.m, i = +e.target.dataset.r, f = e.target.dataset.f;
      rciState[m].rows[i][f] = e.target.value;
      // Update only this row's computed cells (don't destroy focus)
      updateRciRow(m, i, e.target.closest('tr'));
      renderRciApa(m);
    });
  });
  setupRciContextualTabbing(tbody);
  applyRciGroupedHeaders();
  addColumnTitles();
  renderRciApa(method);
}

function updateRciRow(method, i, tr){
  const r = rciState[method].rows[i];
  const cv = rciState[method].cv;
  let calc = null;
  if (method === 'rci-basic') calc = calcBasicRow(r, method);
  else if (method === 'rci-practice') calc = calcPracticeRow(r, method);
  else if (method === 'rci-srb') calc = calcSrbRow(r, method);
  else if (method === 'rci-crawford') calc = calcCrawfordRow(r, method);
  const cells = tr.querySelectorAll('.computed');
  cells.forEach(c => { c.className = 'computed'; });
  tr.classList.remove('row-awaiting','row-check');
  if (!calc){
    const status = numericProblem(r, method);
    const outcomeCell = cells[cells.length - 1];
    cells.forEach(c => c.textContent = '');
    setOutcomeStatus(outcomeCell, status, status === 'Awaiting values' ? 'awaiting' : 'check');
    tr.classList.add(status === 'Awaiting values' && !hasAnyRowValue(r) ? 'row-awaiting' : 'row-check');
    return;
  }
  if (method === 'rci-srb' || method === 'rci-crawford'){
    // Cells: predicted, RCI/t, p, outcome
    if (calc){
      cells[0].textContent = fmt(calc.predicted, 2);
      cells[1].textContent = fmt(calc.rci, 2);
      cells[2].textContent = fmtP(calc.p);
      const oc = rciOutcome(calc.rci, cv, calc.df);
      cells[3].textContent = oc.label;
      if (oc.cls) cells[3].classList.add(oc.cls);
    } else {
      cells.forEach(c => c.textContent = '');
    }
  } else {
    // basic / practice - cells: RCI, p, outcome
    if (calc){
      cells[0].textContent = fmt(calc.rci, 2);
      cells[1].textContent = fmtP(calc.p);
      const oc = rciOutcome(calc.rci, cv);
      cells[2].textContent = oc.label;
      if (oc.cls) cells[2].classList.add(oc.cls);
    } else {
      cells.forEach(c => c.textContent = '');
    }
  }
}

function renderRciApa(method){
  const st = rciState[method];
  const outId = `${method}-apa`;
  const out = document.getElementById(outId);
  const valid = st.rows.filter(r => r.name);
  const cvLabel = st.cv === 0.95 ? '95%' : '90%';
  const cvLabelZ = st.cv === 0.95 ? '95% (z = 1.96)' : '90% (z = 1.645)';
  const methodNote = {
    'rci-basic':    `RCI (z) is the reliable-change statistic expressed as a standard-normal z value, computed per Jacobson and Truax (1991). Reliable change threshold = ${cvLabelZ}.`,
    'rci-practice': `RCI (z) is the reliable-change statistic expressed as a standard-normal z value, computed per Iverson (2001) and adjusted for practice effects. Reliable change threshold = ${cvLabelZ}.`,
    'rci-srb':      `McSweeney Regression-Based RCI (z) per McSweeney et al. (1993). Ŷ₂ = predicted retest score from the regression model; RCI (z) tests whether Date 2 differs reliably from Ŷ₂. Reliable change threshold = ${cvLabelZ}.`,
    'rci-crawford': `<i>t</i>(RB) is the Crawford regression-based reliable-change statistic. Reliable change threshold = ${cvLabel} CI.`
  }[method];
  const safe = (calc, prop, digits=2) => calc ? fmt(calc[prop], digits) : '';
  const safeP = calc => calc ? fmtP(calc.p) : '';
  const safeOutcome = calc => calc ? escapeHtml(rciOutcome(calc.rci, st.cv, calc.df).label) : '';
  const baseColumns = [
    { key:'subtest', label:'Subtest', num:false, render:r => escapeHtml(r.name) },
    { key:'d1', label:escapeHtml(st.d1), group:'Scores', num:true, render:r => escapeHtml(r.t1 || '') },
    { key:'d2', label:escapeHtml(st.d2), group:'Scores', num:true, render:r => escapeHtml(r.t2 || '') }
  ];
  let columns;
  if (method === 'rci-basic'){
    columns = baseColumns.concat([
      { key:'rci', label:'RCI (z)', group:'Results', num:true, render:r => safe(calcBasicRow(r, method), 'rci') },
      { key:'p', label:'<i>p</i>', group:'Results', num:true, render:r => safeP(calcBasicRow(r, method)) },
      { key:'outcome', label:'Outcome', group:'Outcome', num:false, render:r => safeOutcome(calcBasicRow(r, method)) }
    ]);
  } else if (method === 'rci-practice'){
    columns = baseColumns.concat([
      { key:'delta', label:'Δ', group:'Results', num:true, defaultVisible:false, render:r => (r.t1 !== '' && r.t2 !== '' && !isNaN(r.t1) && !isNaN(r.t2)) ? fmt(parseFloat(r.t2) - parseFloat(r.t1), 1) : '' },
      { key:'rci', label:'RCI (z)', group:'Results', num:true, render:r => safe(calcPracticeRow(r, method), 'rci') },
      { key:'p', label:'<i>p</i>', group:'Results', num:true, render:r => safeP(calcPracticeRow(r, method)) },
      { key:'outcome', label:'Outcome', group:'Outcome', num:false, render:r => safeOutcome(calcPracticeRow(r, method)) }
    ]);
  } else if (method === 'rci-srb'){
    columns = baseColumns.concat([
      { key:'predicted', label:'Ŷ₂', group:'Results', num:true, defaultVisible:false, render:r => safe(calcSrbRow(r, method), 'predicted') },
      { key:'rci', label:'RCI (z)', group:'Results', num:true, render:r => safe(calcSrbRow(r, method), 'rci') },
      { key:'p', label:'<i>p</i>', group:'Results', num:true, render:r => safeP(calcSrbRow(r, method)) },
      { key:'outcome', label:'Outcome', group:'Outcome', num:false, render:r => safeOutcome(calcSrbRow(r, method)) }
    ]);
  } else {
    columns = baseColumns.concat([
      { key:'predicted', label:'Ŷ₂', group:'Results', num:true, defaultVisible:false, render:r => safe(calcCrawfordRow(r, method), 'predicted') },
      { key:'trb', label:'<i>t</i>(RB)', group:'Results', num:true, render:r => safe(calcCrawfordRow(r, method), 'rci') },
      { key:'p', label:'<i>p</i>', group:'Results', num:true, render:r => safeP(calcCrawfordRow(r, method)) },
      { key:'outcome', label:'Outcome', group:'Outcome', num:false, render:r => safeOutcome(calcCrawfordRow(r, method)) }
    ]);
  }
  updateApaColumnControls(outId, columns, () => renderRciApa(method));
  if (valid.length === 0){
    out.innerHTML = '<div class="apa-empty"><strong>APA-formatted output</strong>Add or select at least one test to preview the report-ready table.</div>';
    return;
  }
  // Footnote on which r flavour was used (all 4 RCI methods). When the user
  // has the toggle ON but the active rows include any without a corrected r
  // value, append a "fell back to raw r" qualifier.
  let rNote = '';
  if (method === 'rci-basic' || method === 'rci-practice' || method === 'rci-srb' || method === 'rci-crawford'){
    const defaultOn = method === 'rci-basic' || method === 'rci-practice';
    const wantCorrected = st.useCorrectedR === undefined ? defaultOn : st.useCorrectedR !== false;
    if (wantCorrected){
      const fallbackTests = valid
        .filter(r => !Number.isFinite(parseFloat(r.rCorrected)))
        .map(r => r.name);
      if (fallbackTests.length === 0){
        rNote = ' Calculations used the corrected (attenuation-adjusted) test-retest correlation.';
      } else {
        const list = fallbackTests.join(', ');
        rNote = ' Calculations used the corrected (attenuation-adjusted) test-retest correlation where available; raw <i>r</i> was used for the following tests because no corrected <i>r</i> is published: ' + escapeHtml(list) + '.';
      }
    } else {
      rNote = ' Calculations used the raw test-retest correlation.';
    }
  }
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(st.title)}</div>
    ${buildApaTableFromColumns(outId, columns, valid, r => r.group)}
    <div class="apa-note"><strong>Note.</strong> ${methodNote}${rNote} <i>p</i>-values are two-tailed.</div>
  `;
}

// Wire up RCI controls
document.querySelectorAll('.rci-cv').forEach(sel => {
  sel.addEventListener('change', e => { rciState[e.target.dataset.target].cv = parseFloat(e.target.value); renderRci(e.target.dataset.target); });
});
document.querySelectorAll('.rci-d1').forEach(inp => {
  inp.addEventListener('input', e => { rciState[e.target.dataset.target].d1 = e.target.value; renderRci(e.target.dataset.target); });
});
document.querySelectorAll('.rci-d2').forEach(inp => {
  inp.addEventListener('input', e => { rciState[e.target.dataset.target].d2 = e.target.value; renderRci(e.target.dataset.target); });
});
document.querySelectorAll('.rci-title').forEach(inp => {
  inp.addEventListener('input', e => { rciState[e.target.dataset.target].title = e.target.value; renderRciApa(e.target.dataset.target); });
});
// Use corrected r toggle (SRB + Crawford only)
document.querySelectorAll('.rci-use-corrected-r').forEach(cb => {
  cb.addEventListener('change', e => {
    const m = e.target.dataset.target;
    if (rciState[m]){
      rciState[m].useCorrectedR = !!e.target.checked;
      renderRci(m);
    }
  });
});
document.querySelectorAll('[data-add]').forEach(btn => {
  btn.addEventListener('click', () => rciAddRow(btn.dataset.add));
});

/* ============================================================
   PER-METHOD AUTO-FILL FROM NORMATIVE DATABASE
   ============================================================ */

function getMergedDB(){
  const custom = JSON.parse(localStorage.getItem('customTests') || '{}');
  // Merge: built-in first, then custom (custom overrides on conflict)
  return { ...normDB, ...custom };
}
// Age-band helpers ─────────────────────────────────────────────────────────
// Strip "· Ages X-Y", "· Age N", or "· All Ages" suffix to get the base test name
function familyBaseName(f){
  return f.replace(/\s+·\s+(Ages\s+[\d–\-]+|Age\s+\d+|All\s+Ages)\s*$/i, '').trim();
}
// True for entries that carry any age-band suffix (All Ages, Ages X-Y, or Age N)
function hasAgeBandSuffix(f){
  return /·\s*(All\s+Ages|Ages\s+[\d–\-]+|Age\s+\d+)\s*$/i.test(f);
}
// Return just the age-band portion for display in indented items
function ageBandLabel(f){
  const m = f.match(/·\s*((?:All\s+Ages|Ages\s+[\d–\-]+|Age\s+\d+))\s*$/i);
  return m ? escapeHtml(m[1]) : escapeHtml(f);
}
// Build grouped HTML: group header + indented items for age-banded families,
// plain items for everything else.
//
// `opts.flat = true` collapses age-banded families into a single plain entry
// labelled by the base name. The chosen value is the "All Ages" variant when
// available, otherwise the first member. Used on the Neuropsych Tables page
// where age bands don't change the resulting table and just clutter the UI.
function buildFamilyListHtml(families, opts){
  const flat = !!(opts && opts.flat);
  const isCustom = f => !normDB[f];
  const order = [];
  const groups = {};
  families.forEach(f => {
    const base = hasAgeBandSuffix(f) ? familyBaseName(f) : f;
    if (!groups[base]){ groups[base] = []; order.push(base); }
    groups[base].push(f);
  });

  let html = '';
  order.forEach(base => {
    const members = groups[base];
    const groupKey = `grp:${base}`;
    if (members.length === 1 && !hasAgeBandSuffix(members[0])){
      html += comboCheckboxItemHtml(members[0], isCustom(members[0]), false, groupKey);
    } else if (flat){
      // Pick the "All Ages" canonical entry if present, else the first.
      const canon = members.find(m => /·\s*All\s+Ages\s*$/i.test(m)) || members[0];
      html += comboCheckboxItemHtml(canon, isCustom(canon), false, groupKey, base);
    } else {
      html += `<div class="combo-group-heading" data-group="${escapeAttr(groupKey)}">${escapeHtml(base)}</div>`;
      // Wrap age-banded variants in a flex row so they render side-by-side
      // as compact pills instead of stacking - cuts vertical scroll roughly
      // in half on long family lists (e.g. CVLT-3 INDICES + TRIALS, etc.).
      html += `<div class="combo-indented-row" data-group="${escapeAttr(groupKey)}">`;
      members.forEach(f => {
        html += comboCheckboxItemHtml(f, isCustom(f), true, groupKey);
      });
      html += `</div>`;
    }
  });
  return comboOptionsHtml(html);
}
function rebuildAllFamilyLists(){
  document.querySelectorAll('.rci-family-list').forEach(list => populateFamilyList(list));
}
// Returns true iff every subtest in the family carries a valid numeric N (≥3)
function familyHasN(familyEntries){
  if (!familyEntries || Object.keys(familyEntries).length === 0) return false;
  return Object.values(familyEntries).every(p => {
    const n = (p && p.n != null) ? Number(p.n) : NaN;
    return Number.isFinite(n) && n >= 3;
  });
}
function populateFamilyList(list){
  const db = getMergedDB();
  const method = list.dataset.method;
  let families = Object.keys(db).sort();
  // Crawford method requires N for the t-distributed test statistic - hide
  // families where no subtest carries a usable N. CVLT-3 now ships an N=100
  // holding value so it qualifies; WAIS-IV age-band data still doesn't.
  if (method === 'rci-crawford'){
    families = families.filter(f => familyHasN(db[f]));
  }
  list.innerHTML = comboFooterHtml() + comboAgeBandNoteHtml() + buildFamilyListHtml(families);
  wireMultiSelectFamilyList(list, families => {
    families.forEach(family => loadFamilyIntoMethod(method, family));
    const inp = document.querySelector(`.rci-family-input[data-method="${method}"]`);
    if (inp){ inp.value = ''; inp.focus(); }
  });
}
function loadFamilyIntoMethod(method, family){
  const db = getMergedDB();
  if (!db[family]) return;
  const subtests = Object.entries(db[family]);
  let newRows;
  if (method === 'rci-basic'){
    newRows = subtests.map(([name, p]) => ({
      name, group:family,
      sd: String(p.sd1),
      r: String(p.r),
      // Carry corrected r when present so the user-level toggle can pick it
      rCorrected: (p.rCorrected != null ? String(p.rCorrected) : ''),
      t1:'', t2:''
    }));
  } else if (method === 'rci-crawford'){
    newRows = subtests.map(([name, p]) => ({
      name,
      group:family,
      m1: String(p.m1), sd1: String(p.sd1),
      m2: String(p.m2), sd2: String(p.sd2),
      r: String(p.r),
      // Carry corrected r when present so the user-level toggle can pick it
      rCorrected: (p.rCorrected != null ? String(p.rCorrected) : ''),
      n: (p.n != null ? String(p.n) : ''),
      t1:'', t2:''
    }));
  } else { // rci-practice or rci-srb
    newRows = subtests.map(([name, p]) => ({
      name,
      group:family,
      m1: String(p.m1), sd1: String(p.sd1),
      m2: String(p.m2), sd2: String(p.sd2),
      r: String(p.r),
      rCorrected: (p.rCorrected != null ? String(p.rCorrected) : ''),
      t1:'', t2:''
    }));
  }
  // Append new auto-filled tests rather than replacing any tests already entered.
  rciState[method].rows = rciState[method].rows.concat(newRows);
  renderRci(method);
  // Toast suppressed - the working-report pill is the single feedback channel
}
function clearMethodRows(method){
  rciState[method].rows = [];
  renderRci(method);
  const inp = document.querySelector(`.rci-family-input[data-method="${method}"]`);
  if (inp) inp.value = '';
}
// Wire up the comboboxes
document.querySelectorAll('.rci-family-input').forEach(inp => {
  const method = inp.dataset.method;
  const list = document.querySelector(`.rci-family-list[data-method="${method}"]`);
  inp.addEventListener('focus', () => { list.classList.add('show'); filterFamilyListEl(list, ''); });
  inp.addEventListener('input', e => { list.classList.add('show'); filterFamilyListEl(list, e.target.value); });
  inp.addEventListener('keydown', e => {
    if (e.key === 'Escape') list.classList.remove('show');
    if (e.key === 'Enter') {
      const add = list.querySelector('.combo-add:not([disabled])');
      if (add){ e.preventDefault(); add.click(); }
    }
  });
  inp.addEventListener('blur', () => setTimeout(() => {
    if (!list.matches(':hover')) list.classList.remove('show');
  }, 180));
});
function filterFamilyListEl(list, q){
  const ql = q.toLowerCase().trim();
  // Show/hide individual items; also match against the full family value (base name + age band)
  list.querySelectorAll('.combo-item').forEach(it => {
    const familyVal = (it.dataset.family || '').toLowerCase();
    it.style.display = !ql || familyVal.includes(ql) ? '' : 'none';
  });
  // Show a group heading only when at least one of its sibling items is visible
  list.querySelectorAll('.combo-group-heading').forEach(heading => {
    const groupKey = heading.dataset.group;
    const anyVisible = Array.from(list.querySelectorAll('.combo-item'))
      .filter(it => it.dataset.group === groupKey)
      .some(it => it.style.display !== 'none');
    heading.style.display = anyVisible ? '' : 'none';
  });
}
document.addEventListener('mousedown', e => {
  if (!e.target.closest('.combo')) {
    document.querySelectorAll('.combo-list.show').forEach(list => list.classList.remove('show'));
  }
});
// Wire up clear buttons
document.querySelectorAll('[data-clear]').forEach(btn => {
  btn.addEventListener('click', () => clearMethodRows(btn.dataset.clear));
});

/* ============================================================
   08 · CUSTOM TESTS DATABASE MANAGEMENT
   ============================================================ */
function getCustom(){ return JSON.parse(localStorage.getItem('customTests') || '{}'); }
function saveCustom(c){ localStorage.setItem('customTests', JSON.stringify(c)); }

function refreshFamilySelect(selectedFamily){
  const sel = document.getElementById('ct-family-select');
  const current = selectedFamily || sel.value;
  const merged = getMergedDB();
  const custom = getCustom();
  // Allow adding subtests to any family (built-in or custom). Built-in additions become custom overrides.
  const families = Object.keys(merged).sort();
  sel.innerHTML = families.map(f =>
    `<option value="${escapeAttr(f)}">${escapeHtml(f)}${custom[f] ? ' (custom)' : ''}</option>`
  ).join('');
  if (current && families.includes(current)){
    sel.value = current;
  }
}
function makeFamilyName(baseName, ageMode, ageMinRaw, ageMaxRaw){
  const base = String(baseName || '').trim();
  if (!base) return '';
  const mode = String(ageMode || 'none');
  if (mode === 'all'){
    return `${base} · All Ages`;
  }
  if (mode === 'range'){
    const min = parseInt(ageMinRaw, 10);
    const max = parseInt(ageMaxRaw, 10);
    if (!Number.isFinite(min) || !Number.isFinite(max)) return '';
    if (min < 0 || max < 0 || min > max) return '';
    return `${base} · Ages ${min}-${max}`;
  }
  return base;
}
function renderDbList(){
  const search = document.getElementById('ct-search').value.toLowerCase().trim();
  const merged = getMergedDB();
  const custom = getCustom();
  const list = document.getElementById('db-list');
  const families = Object.keys(merged).sort();
  list.innerHTML = '';
  families.forEach(family => {
    const isCustom = !!custom[family];
    const subtests = merged[family];
    const familyMatch = !search || family.toLowerCase().includes(search);
    const matchedSubs = Object.entries(subtests).filter(([n]) => familyMatch || n.toLowerCase().includes(search));
    if (!familyMatch && matchedSubs.length === 0) return;
    const subsToShow = familyMatch ? Object.entries(subtests) : matchedSubs;
    const fam = document.createElement('div');
    fam.className = 'db-test-family';
    fam.innerHTML = `
      <div class="db-test-family-header">
        <span>${escapeHtml(family)} ${isCustom ? '<span class="custom-tag">Custom</span>' : ''}</span>
        ${isCustom ? `<button class="btn btn-ghost" data-del-family="${escapeAttr(family)}">Delete family</button>` : ''}
      </div>
      <div class="db-subtests">
        <div class="db-subtest-row" style="font-size:10px;color:var(--faint);text-transform:uppercase;letter-spacing:0.1em;border-bottom:1px solid var(--border-soft);padding-bottom:6px">
          <span>Subtest</span><span class="num-label">M₁</span><span class="num-label">SD₁</span><span class="num-label">M₂</span><span class="num-label">SD₂</span><span class="num-label">r</span><span class="num-label" title="Attenuation-corrected test–retest correlation">corr. r</span><span class="num-label">N</span><span></span>
        </div>
        ${subsToShow.map(([name, p]) => `
          <div class="db-subtest-row">
            <span class="name">${escapeHtml(name)}</span>
            <span class="num">${fmt(p.m1, 2)}</span>
            <span class="num">${fmt(p.sd1, 2)}</span>
            <span class="num">${fmt(p.m2, 2)}</span>
            <span class="num">${fmt(p.sd2, 2)}</span>
            <span class="num">${fmt(p.r, 2)}</span>
            <span class="num">${p.rCorrected == null ? '-' : fmt(p.rCorrected, 2)}</span>
            <span class="num">${p.n == null ? '-' : escapeHtml(String(p.n))}</span>
            ${isCustom ? `<button class="btn btn-ghost btn-icon" data-del-sub="${escapeAttr(family)}::${escapeAttr(name)}" title="Remove subtest">×</button>` : '<span></span>'}
          </div>
        `).join('')}
      </div>
    `;
    list.appendChild(fam);
  });
  list.querySelectorAll('[data-del-family]').forEach(b => b.addEventListener('click', () => {
    if (!confirm(`Delete custom family "${b.dataset.delFamily}" and all its subtests?`)) return;
    const c = getCustom(); delete c[b.dataset.delFamily]; saveCustom(c);
    refreshAll();
  }));
  list.querySelectorAll('[data-del-sub]').forEach(b => b.addEventListener('click', () => {
    const [family, sub] = b.dataset.delSub.split('::');
    const c = getCustom();
    if (c[family]){ delete c[family][sub]; if (Object.keys(c[family]).length === 0) delete c[family]; saveCustom(c); }
    refreshAll();
  }));
}
function refreshAll(selectedFamily){
  refreshFamilySelect(selectedFamily);
  renderDbList();
  rebuildAllFamilyLists();
  rebuildBatteryFamilyList();
  rebuildSdiFamilyList();
  refreshReportWriterOptions();
}
function ctEntryRowHtml(i, seed){
  const row = seed || {};
  return `<tr>
    <td class="ct-row-num">${i + 1}</td>
    <td><input class="ct-subtest-name" type="text" data-k="name" placeholder="Test name" value="${escapeAttr(row.name || '')}"></td>
    <td><input type="number" step="any" data-k="m1" value="${row.m1 != null ? escapeAttr(row.m1) : ''}"></td>
    <td><input type="number" step="any" data-k="sd1" value="${row.sd1 != null ? escapeAttr(row.sd1) : ''}"></td>
    <td><input type="number" step="any" data-k="m2" value="${row.m2 != null ? escapeAttr(row.m2) : ''}"></td>
    <td><input type="number" step="any" data-k="sd2" value="${row.sd2 != null ? escapeAttr(row.sd2) : ''}"></td>
    <td><input type="number" step="0.01" min="0" max="0.99" data-k="r" value="${row.r != null ? escapeAttr(row.r) : ''}"></td>
    <td><input type="number" step="1" min="3" data-k="n" value="${row.n != null ? escapeAttr(row.n) : ''}"></td>
    <td><button type="button" class="btn btn-ghost btn-icon" data-ct-del-row title="Delete row">×</button></td>
  </tr>`;
}
function ctRenumberRows(){
  const body = document.getElementById('ct-entry-body');
  if (!body) return;
  body.querySelectorAll('tr').forEach((tr, i) => {
    const num = tr.querySelector('.ct-row-num');
    if (num) num.textContent = String(i + 1);
  });
}
function ctRenderRows(rows){
  const body = document.getElementById('ct-entry-body');
  if (!body) return;
  body.innerHTML = rows.map((r, i) => ctEntryRowHtml(i, r)).join('');
}
function ctReadRows(){
  const body = document.getElementById('ct-entry-body');
  if (!body) return [];
  const out = [];
  body.querySelectorAll('tr').forEach(tr => {
    const get = k => tr.querySelector(`[data-k="${k}"]`)?.value?.trim() || '';
    const name = get('name');
    const m1s = get('m1'), sd1s = get('sd1'), m2s = get('m2'), sd2s = get('sd2'), rs = get('r'), ns = get('n');
    const any = [name,m1s,sd1s,m2s,sd2s,rs,ns].some(v => v !== '');
    if (!any) return;
    out.push({ name, m1s, sd1s, m2s, sd2s, rs, ns });
  });
  return out;
}
function ctInitEntryRows(){
  ctRenderRows([{ name:'Example index score', m1:100, sd1:15, m2:103, sd2:15, r:0.90, n:100 }, {}]);
  const body = document.getElementById('ct-entry-body');
  const addBlankRow = () => {
    const rows = ctReadRows().map(r => ({ name:r.name, m1:r.m1s, sd1:r.sd1s, m2:r.m2s, sd2:r.sd2s, r:r.rs, n:r.ns }));
    rows.push({});
    ctRenderRows(rows);
  };
  document.getElementById('ct-add-row')?.addEventListener('click', () => {
    addBlankRow();
  });
  document.getElementById('ct-load-example')?.addEventListener('click', () => {
    const rows = ctReadRows().map(r => ({ name:r.name, m1:r.m1s, sd1:r.sd1s, m2:r.m2s, sd2:r.sd2s, r:r.rs, n:r.ns }));
    rows.push({ name:'Example trial', m1:10, sd1:3, m2:11, sd2:3, r:0.80, n:100 });
    ctRenderRows(rows);
  });
  document.getElementById('ct-clear-rows')?.addEventListener('click', () => ctRenderRows([{}]));
  body?.addEventListener('click', e => {
    const del = e.target.closest('[data-ct-del-row]');
    if (!del) return;
    const tr = del.closest('tr');
    if (!tr) return;
    tr.remove();
    if (!body.querySelector('tr')){
      ctRenderRows([{}]);
    } else {
      ctRenumberRows();
    }
  });
  body?.addEventListener('keydown', e => {
    if (e.key !== 'Tab' || e.shiftKey) return;
    const target = e.target;
    if (!(target instanceof HTMLInputElement)) return;
    const tr = target.closest('tr');
    if (!tr) return;
    const rows = Array.from(body.querySelectorAll('tr'));
    if (rows[rows.length - 1] !== tr) return;
    const inputs = Array.from(tr.querySelectorAll('input'));
    if (inputs[inputs.length - 1] !== target) return;
    e.preventDefault();
    addBlankRow();
    const nextRows = body.querySelectorAll('tr');
    const newLast = nextRows[nextRows.length - 1];
    const firstInput = newLast?.querySelector('input[data-k="name"]');
    if (firstInput) firstInput.focus();
  });
}
document.getElementById('ct-add-family').addEventListener('click', () => {
  const base = document.getElementById('ct-family').value.trim();
  const ageMode = document.getElementById('ct-family-age-mode')?.value || 'none';
  const ageMin = document.getElementById('ct-family-age-min')?.value || '';
  const ageMax = document.getElementById('ct-family-age-max')?.value || '';
  const name = makeFamilyName(base, ageMode, ageMin, ageMax);
  if (!base) return showToast('Enter a family name', true);
  if (!name) return showToast('Check age-band fields (use valid min/max ages)', true);
  const c = getCustom();
  if (c[name] || normDB[name]) return showToast('A family with that name already exists', true);
  c[name] = {}; saveCustom(c);
  document.getElementById('ct-family').value = '';
  if (document.getElementById('ct-family-age-mode')) document.getElementById('ct-family-age-mode').value = 'none';
  if (document.getElementById('ct-family-age-min')) document.getElementById('ct-family-age-min').value = '';
  if (document.getElementById('ct-family-age-max')) document.getElementById('ct-family-age-max').value = '';
  showToast(`✓ Added family "${name}"`);
  refreshAll(name);
});
document.getElementById('ct-add-subtest').addEventListener('click', () => {
  const family = document.getElementById('ct-family-select').value;
  const rows = ctReadRows();
  if (!family) return showToast('Select a family first', true);
  if (rows.length === 0) return showToast('Add at least one row', true);
  const c = getCustom();
  // If user is adding to a built-in family, clone it first into custom (so we don't lose built-ins on conflict)
  if (!c[family]){
    if (normDB[family]) c[family] = { ...normDB[family] };
    else c[family] = {};
  }
  let added = 0;
  for (const r of rows){
    const subtest = r.name;
    const m1 = parseFloat(r.m1s), sd1 = parseFloat(r.sd1s), m2 = parseFloat(r.m2s), sd2 = parseFloat(r.sd2s), rel = parseFloat(r.rs);
    const n = r.ns === '' ? null : parseInt(r.ns, 10);
    if (!subtest) return showToast('Each used row needs a subtest name', true);
    if ([m1, sd1, m2, sd2, rel].some(isNaN)) return showToast(`Fill all numeric fields for "${subtest}"`, true);
    if (rel < 0 || rel >= 1) return showToast(`Reliability out of range for "${subtest}"`, true);
    if (sd1 <= 0 || sd2 <= 0) return showToast(`SD must be positive for "${subtest}"`, true);
    if (n !== null && (!Number.isFinite(n) || n < 3)) return showToast(`N must be ≥ 3 for "${subtest}"`, true);
    const entry = { m1, sd1, m2, sd2, r: rel };
    if (n !== null) entry.n = n;
    c[family][subtest] = entry;
    added += 1;
  }
  saveCustom(c);
  ctRenderRows([{}]);
  showToast(`✓ Added ${added} row${added === 1 ? '' : 's'} to ${family}`);
  refreshAll();
});
document.getElementById('ct-search').addEventListener('input', renderDbList);
const ctAgeMode = document.getElementById('ct-family-age-mode');
const ctAgeRange = document.getElementById('ct-family-age-range');
if (ctAgeMode && ctAgeRange){
  const syncCtAgeRange = () => { ctAgeRange.style.display = ctAgeMode.value === 'range' ? 'grid' : 'none'; };
  ctAgeMode.addEventListener('change', syncCtAgeRange);
  syncCtAgeRange();
}
ctInitEntryRows();
document.getElementById('ct-export').addEventListener('click', () => {
  const blob = new Blob([JSON.stringify(getCustom(), null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'custom-tests.json';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
});
document.getElementById('ct-import').addEventListener('change', e => {
  const file = e.target.files[0]; if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const data = JSON.parse(reader.result);
      if (typeof data !== 'object') throw 0;
      // Merge with existing custom
      const c = getCustom();
      Object.entries(data).forEach(([fam, subs]) => {
        c[fam] = { ...(c[fam] || {}), ...subs };
      });
      saveCustom(c);
      showToast('✓ Imported custom tests');
      refreshAll();
    } catch (err) {
      showToast('Invalid JSON file', true);
    }
  };
  reader.readAsText(file);
  e.target.value = '';
});

/* ============================================================
   08 · PREMORBID ESTIMATE
   Ports the user's reference implementation. All formulas/data
   from WAIS-IV/WMS-IV manuals + Crawford & Allan (2001) + OPIE-4.
   ============================================================ */

function preModelCell(label, tipKey){
  const tip = PRE_MODEL_TOOLTIPS[tipKey] || 'Hover information is not available for this model.';
  return `<td class="model model-has-tip">${escapeHtml(label)}<span class="model-info-dot" data-tooltip="${escapeAttr(tip)}" aria-label="${escapeAttr(label + '. ' + tip)}">?</span></td>`;
}

// Premorbid state - separate from the input fields so we can re-render APA reliably
const preState = { achieved: {}, opieAchieved: {} }; // achieved[idx] = string; opieAchieved[modelKey] = string

function preGet(id){ return document.getElementById(id); }
function preNum(id){
  const el = preGet(id);
  if (!el || el.value === '' || el.value == null) return null;
  const n = parseFloat(el.value);
  return isNaN(n) ? null : n;
}
function preStr(id){
  const el = preGet(id);
  return el && el.value ? el.value : null;
}
function preCiMult(){
  const ci = preStr('pre-ci');
  return ci === '0.95' ? 1.96 : 1.645;
}
function fmtIntOrDash(v){
  if (v == null || isNaN(v)) return '-';
  return Math.round(v).toString();
}
function fmtPctBr(v){
  if (v == null || isNaN(v)) return '-';
  return (v*100).toFixed(2) + '%';
}

// Tab switching for premorbid section
const PRE_TAB_ORDER = ['inputs', 'estimates', 'predict', 'opiepredict'];
const PRE_TAB_LABELS = {
  inputs: 'Inputs',
  estimates: 'Premorbid Estimates',
  predict: 'ToPF vs WAIS-IV',
  opiepredict: 'OPIE-4 vs WAIS-IV'
};

function switchPreTab(tabId){
  document.querySelectorAll('[data-pre-tab]').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.pre-tab-content').forEach(c => c.classList.remove('active'));
  const tabBtn = document.querySelector(`[data-pre-tab="${tabId}"]`);
  if (tabBtn) tabBtn.classList.add('active');
  const content = document.getElementById('pre-' + tabId);
  if (content) content.classList.add('active');
  updatePreNav(tabId);

  // Smart scroll: only scroll the tab strip into view if it's currently NOT visible.
  // This prevents the page from jumping when the user is already looking at the tabs,
  // but ensures they see the new content if they were scrolled far below.
  const tabsStrip = document.querySelector('#premorbid .pre-tabs');
  if (tabsStrip){
    const rect = tabsStrip.getBoundingClientRect();
    const topbarH = 108; // fixed topbar height
    // If tabs are above viewport (negative top) or below the fold, bring them just under the topbar
    if (rect.top < topbarH + 8 || rect.top > window.innerHeight){
      const offset = window.scrollY + rect.top - topbarH - 16;
      window.scrollTo({ top: offset, behavior: 'smooth' });
    }
    // Otherwise: do nothing - user stays exactly where they were, content swaps in place
  }
}

function updatePreNav(activeTabId){
  const idx = PRE_TAB_ORDER.indexOf(activeTabId);
  if (idx === -1) return;
  const back = document.querySelector('.pre-nav-back');
  const next = document.querySelector('.pre-nav-next');
  const step = document.getElementById('pre-nav-step');
  if (!back || !next) return;

  // Step counter ("Step 2 of 4")
  if (step) step.textContent = `Step ${idx + 1} of ${PRE_TAB_ORDER.length}`;

  // Back button
  if (idx > 0){
    const prevId = PRE_TAB_ORDER[idx - 1];
    back.style.visibility = 'visible';
    back.dataset.goTab = prevId;
    const lbl = back.querySelector('.pre-nav-back-label');
    if (lbl) lbl.textContent = 'Back: ' + PRE_TAB_LABELS[prevId];
  } else {
    back.style.visibility = 'hidden';
  }

  // Next button
  if (idx < PRE_TAB_ORDER.length - 1){
    const nextId = PRE_TAB_ORDER[idx + 1];
    next.style.visibility = 'visible';
    next.dataset.goTab = nextId;
    const lbl = next.querySelector('.pre-nav-next-label');
    if (lbl) lbl.textContent = 'Next: ' + PRE_TAB_LABELS[nextId];
  } else {
    next.style.visibility = 'hidden';
  }
}

function setupPreTabs(){
  document.querySelectorAll('[data-pre-tab]').forEach(tab => {
    tab.addEventListener('click', () => switchPreTab(tab.dataset.preTab));
  });
  const back = document.querySelector('.pre-nav-back');
  const next = document.querySelector('.pre-nav-next');
  if (back) back.addEventListener('click', () => { if (back.dataset.goTab) switchPreTab(back.dataset.goTab); });
  if (next) next.addEventListener('click', () => { if (next.dataset.goTab) switchPreTab(next.dataset.goTab); });
  // Initialise nav state for the default-active tab
  updatePreNav('inputs');
}

// Build the achieved-input table for ToPF vs WAIS-IV
function buildPredictTable(){
  const tbody = document.querySelector('#pre-predict-table tbody');
  tbody.innerHTML = '';

  // WAIS-IV header
  const waisH = document.createElement('tr');
  waisH.className = 'group-row';
  waisH.innerHTML = '<td colspan="7">WAIS-IV</td>';
  tbody.appendChild(waisH);
  WAIS_COEF.forEach(c => {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      ${preModelCell(c.label, 'predictWais')}
      <td class="num" id="pred-${c.idx}">-</td>
      <td class="num" id="pred-${c.idx}-lo">-</td>
      <td class="num" id="pred-${c.idx}-hi">-</td>
      <td class="achieved-cell"><input type="number" min="40" max="160" step="1" data-pre-ach="${c.idx}" value="${escapeAttr(preState.achieved[c.idx] || '')}"></td>
      <td class="num diff" id="diff-${c.idx}">-</td>
      <td class="num" id="br-${c.idx}">-</td>
    `;
    tbody.appendChild(tr);
  });

  // WMS-IV header
  const wmsH = document.createElement('tr');
  wmsH.className = 'group-row';
  wmsH.innerHTML = '<td colspan="7">WMS-IV</td>';
  tbody.appendChild(wmsH);
  WMS_COEF.forEach(c => {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      ${preModelCell(c.label, 'predictWms')}
      <td class="num" id="pred-${c.idx}">-</td>
      <td class="num" id="pred-${c.idx}-lo">-</td>
      <td class="num" id="pred-${c.idx}-hi">-</td>
      <td class="achieved-cell"><input type="number" min="40" max="160" step="1" data-pre-ach="${c.idx}" value="${escapeAttr(preState.achieved[c.idx] || '')}"></td>
      <td class="num diff" id="diff-${c.idx}">-</td>
      <td class="num" id="br-${c.idx}">-</td>
    `;
    tbody.appendChild(tr);
  });

  // Wire achieved inputs
  document.querySelectorAll('[data-pre-ach]').forEach(inp => {
    inp.addEventListener('input', e => {
      preState.achieved[e.target.dataset.preAch] = e.target.value;
      calcPredict();
    });
  });
}

// === Premorbid Estimates calculation ===
function calcPremorbid(){
  const topf = preNum('pre-topf');
  const vc   = preNum('pre-vc');
  const mr   = preNum('pre-mr');
  const sex  = preStr('pre-sex');
  const edu  = preNum('pre-edu');
  const age  = preNum('pre-age');
  const occ  = preStr('pre-occ');
  const ci   = preStr('pre-ci') || '0.90';
  // Sex coding differs between models:
  //   TOPF / WAIS-IV manual demographic equations:  Female = 1, Male = 2
  //   OPIE-4 (Schoenberg et al., 2011):             Female = 0, Male = 1
  const sexC_topf = sex === 'Female' ? 1 : sex === 'Male' ? 2 : null;
  const sexC_opie = sex === 'Female' ? 0 : sex === 'Male' ? 1 : null;
  const occC = occ ? OCC_CODE[occ] : null;
  const mult = preCiMult();
  const ciPct = ci === '0.95' ? '95%' : '90%';

  preGet('pre-ci-lo-hdr').textContent = `Lower ${ciPct}`;
  preGet('pre-ci-hi-hdr').textContent = `Upper ${ciPct}`;

  // Build four model rows
  const rows = [];

  // 1. ToPF Raw Only (lookup)
  let v1 = null;
  if (topf != null && topf >= 0 && topf <= 70 && Number.isFinite(topf)){
    v1 = TOPF_TO_FSIQ[Math.round(topf)];
  }
  rows.push({ name:'ToPF Raw Score', val:v1, see:9.867, r:0.72, tipKey:'topfRaw' });

  // 2. ToPF + Demographics  (uses TOPF sex coding F=1, M=2)
  let v2 = null;
  if (topf != null && edu != null && sexC_topf != null){
    v2 = 29.991 + 2.09426*topf + (-0.0404559)*topf*topf
       + 0.000340705*Math.pow(topf,3) + 1.4617126*edu + 4.925*sexC_topf;
  }
  rows.push({ name:'Demographic Adjusted ToPF', val:v2, see:8.441, r:0.81, tipKey:'topfDemo' });

  // 3. Crawford & Allan (2001)
  let v3 = null;
  if (occC != null && edu != null && age != null){
    v3 = 87.14 - 5.21*occC + 1.78*edu + 0.18*age;
  }
  rows.push({ name:'Crawford & Allan (2001) Demographic', val:v3, see:9.11, r:0.73, tipKey:'crawfordAllan' });

  // 4. OPIE-4 - prorated FSIQ, uses OPIE sex coding F=0, M=1
  // Label, R and SEE update as soon as subtest inputs are present (branch alone).
  // FSIQ computation also requires age; sex term contributes 0 if not entered.
  const sexEffect = sexC_opie != null ? sexC_opie : 0;
  const hasVC = vc != null, hasMR = mr != null;
  let branch = null;
  if (hasVC && hasMR) branch = 'VC_MR';
  else if (hasVC) branch = 'VC';
  else if (hasMR) branch = 'MR';

  let v4 = null, s4 = null, name4 = 'OPIE-4: Vocab and/or Matrix', tipKey4 = 'opieDefault';
  if (branch != null){
    const c = OPIE_PRORATED_FSIQ[branch];
    // Update label, tooltip and R/SEE from the branch alone (no age required)
    s4 = { see: c.see, r: c.r };
    if (branch === 'VC_MR'){ name4 = 'OPIE-4: Vocab + Matrix'; tipKey4 = 'opieVCMR'; }
    else if (branch === 'VC') { name4 = 'OPIE-4: Vocab only';  tipKey4 = 'opieVC'; }
    else if (branch === 'MR') { name4 = 'OPIE-4: Matrix only'; tipKey4 = 'opieMR'; }
    // Compute prediction only once age is also available
    if (age != null){
      let pred = c.intercept + (c.age != null ? c.age * age : 0) + c.age3 * Math.pow(age, 3) + c.sex * sexEffect;
      if (c.vc != null && hasVC) pred += c.vc * vc;
      if (c.mr != null && hasMR) pred += c.mr * mr;
      v4 = pred;
    }
  }
  rows.push({ name:name4, val:v4, see:s4 ? s4.see : null, r:s4 ? s4.r : null, tipKey:tipKey4 });

  // Render results table
  const tbody = document.querySelector('#pre-results-table tbody');
  tbody.innerHTML = rows.map(row => {
    const fsiq = fmtIntOrDash(row.val);
    const lo   = (row.val == null || row.see == null) ? '-' : fmtIntOrDash(row.val - mult*row.see);
    const hi   = (row.val == null || row.see == null) ? '-' : fmtIntOrDash(row.val + mult*row.see);
    const rStr = row.r != null ? row.r.toFixed(2) : '-';
    const seeStr = row.see != null ? row.see.toFixed(2) : '-';
    return `<tr>
      ${preModelCell(row.name, row.tipKey)}
      <td class="fsiq">${fsiq}</td>
      <td class="num">${lo}</td>
      <td class="num">${hi}</td>
      <td class="psy">${rStr}</td>
      <td class="psy">${seeStr}</td>
    </tr>`;
  }).join('');

  preState.estimateRows = rows; // cache for APA
  preState.ciMult = mult;
  preState.ciPct = ciPct;
  renderPreEstimatesApa();
}

// === Predicted vs Actual calculation (ToPF-based) ===
function calcPredict(){
  const topf = preNum('pre-topf');
  const sex  = preStr('pre-sex');
  const edu  = preNum('pre-edu');
  const age  = preNum('pre-age');
  const ci   = preStr('pre-ci') || '0.90';
  // ToPF-predicted WAIS-IV indices use the WAIS-IV manual sex coding: F = 1, M = 2.
  const sexC_topf = sex === 'Female' ? 1 : sex === 'Male' ? 2 : null;
  const mult = preCiMult();
  const ciPct = ci === '0.95' ? '95%' : '90%';

  document.querySelectorAll('.pre-ci-lo-hdr-2').forEach(e => e.textContent = `Lower ${ciPct}`);
  document.querySelectorAll('.pre-ci-hi-hdr-2').forEach(e => e.textContent = `Upper ${ciPct}`);

  // WAIS-IV predictions
  WAIS_COEF.forEach(c => {
    let pred = null;
    if (topf != null && edu != null && sexC_topf != null){
      pred = c.intercept + c.b1*topf + c.b2*topf*topf + c.b3*Math.pow(topf,3) + c.edu*edu + c.sex*sexC_topf;
    }
    updatePredictRow(c.idx, pred, mult, c.see);
  });
  // WMS-IV predictions
  WMS_COEF.forEach(c => {
    let pred = null;
    if (topf != null && age != null){
      pred = c.intercept + c.b1*topf + c.age*age;
    }
    updatePredictRow(c.idx, pred, mult, c.see);
  });

  preState.ciPct = ciPct;
  renderPrePredictApa();
}

function updatePredictRow(idx, pred, mult, see){
  const predEl = preGet('pred-'+idx);
  const loEl = preGet('pred-'+idx+'-lo');
  const hiEl = preGet('pred-'+idx+'-hi');
  if (!predEl) return;
  predEl.textContent = pred == null ? '-' : fmtIntOrDash(pred);
  loEl.textContent = pred == null ? '-' : fmtIntOrDash(pred - mult*see);
  hiEl.textContent = pred == null ? '-' : fmtIntOrDash(pred + mult*see);

  const ach = preState.achieved[idx];
  const achNum = (ach != null && ach !== '') ? parseFloat(ach) : null;
  const diffEl = preGet('diff-'+idx);
  const brEl = preGet('br-'+idx);
  diffEl.className = 'num diff';
  if (pred == null || achNum == null || isNaN(achNum)){
    diffEl.textContent = '-';
    brEl.textContent = '-';
    return;
  }
  const diff = Math.round(achNum - pred);
  diffEl.textContent = (diff > 0 ? '+' : '') + diff;
  if (diff < 0) diffEl.classList.add('neg');
  else if (diff > 0) diffEl.classList.add('pos');
  // Base rate only meaningful for negative discrepancies
  const row = BASE_RATES[String(diff)];
  if (row && row[idx] != null) brEl.textContent = fmtPctBr(row[idx]);
  else brEl.textContent = diff < 0 ? '< 0.01%' : '-';
}

// === APA Output: Premorbid Estimates ===
function renderPreEstimatesApa(){
  const out = preGet('pre-estimates-apa');
  if (!preState.estimateRows){ out.innerHTML = ''; return; }
  const valid = preState.estimateRows.filter(r => r.val != null);
  if (valid.length === 0){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Enter at least the ToPF raw score to generate estimates.</div>';
    return;
  }
  const title = preStr('pre-title') || 'Premorbid cognitive estimate';
  const ciPct = preState.ciPct || '90%';
  const mult = preState.ciMult || 1.645;
  const rows = valid.map(r => {
    const lo = (r.see != null) ? Math.round(r.val - mult*r.see) : '-';
    const hi = (r.see != null) ? Math.round(r.val + mult*r.see) : '-';
    const rStr = r.r != null ? r.r.toFixed(2) : '-';
    const seeStr = r.see != null ? r.see.toFixed(2) : '-';
    return `<tr><td>${escapeHtml(r.name)}</td><td class="num">${Math.round(r.val)}</td><td class="num">${lo}</td><td class="num">${hi}</td><td class="num">${rStr}</td><td class="num">${seeStr}</td></tr>`;
  }).join('');
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(title)}</div>
    <table class="apa-table">
      <thead>
        <tr><th>Model</th><th class="num">FSIQ</th><th class="num">Lower ${ciPct}</th><th class="num">Upper ${ciPct}</th><th class="num"><i>r</i></th><th class="num">SEE</th></tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>
    <div class="apa-note"><strong>Note.</strong> FSIQ = Full Scale IQ estimate. CI = confidence interval based on ${ciPct === '95%' ? '1.96' : '1.645'} × SEE. <i>r</i> = correlation between predictor(s) and criterion. SEE = standard error of estimate.</div>
  `;
}

// === APA Output: ToPF vs WAIS-IV ===
function renderPrePredictApa(){
  const out = preGet('pre-predict-apa');
  if (!out) return;
  const title = preStr('pre-title') || 'Premorbid cognitive estimate';
  const ciPct = preState.ciPct || '90%';

  // Filter to indices where the patient has provided an Achieved score
  const allCoef = [...WAIS_COEF.map(c => ({...c, family:'WAIS-IV'})), ...WMS_COEF.map(c => ({...c, family:'WMS-IV'}))];
  const withAch = allCoef.filter(c => {
    const a = preState.achieved[c.idx];
    return a != null && a !== '' && !isNaN(parseFloat(a));
  });
  if (withAch.length === 0){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Enter at least one Achieved score to generate the discrepancy table.</div>';
    return;
  }

  // Group by family
  const byFamily = {};
  withAch.forEach(c => { (byFamily[c.family] = byFamily[c.family] || []).push(c); });

  let body = '';
  Object.entries(byFamily).forEach(([fam, items]) => {
    body += `<tr><td colspan="7" style="font-style:italic;padding-top:6px;padding-bottom:2px">${fam}</td></tr>`;
    items.forEach(c => {
      const pred = preGet('pred-'+c.idx).textContent;
      const lo = preGet('pred-'+c.idx+'-lo').textContent;
      const hi = preGet('pred-'+c.idx+'-hi').textContent;
      const ach = preState.achieved[c.idx];
      const diff = preGet('diff-'+c.idx).textContent;
      const br = preGet('br-'+c.idx).textContent;
      body += `<tr><td>&nbsp;&nbsp;${escapeHtml(c.label)}</td><td class="num">${pred}</td><td class="num">${lo}</td><td class="num">${hi}</td><td class="num">${ach}</td><td class="num">${diff}</td><td class="num">${br}</td></tr>`;
    });
  });

  out.innerHTML = `
    <div class="apa-table-num">Table 2</div>
    <div class="apa-table-title">ToPF-predicted versus achieved index scores</div>
    <table class="apa-table">
      <thead>
        <tr><th>Index</th><th class="num">Predicted score</th><th class="num">Lower ${ciPct}</th><th class="num">Upper ${ciPct}</th><th class="num">Achieved score</th><th class="num">Difference</th><th class="num">Base Rate</th></tr>
      </thead>
      <tbody>${body}</tbody>
    </table>
    <div class="apa-note"><strong>Note.</strong> WAIS-IV indices predicted from ToPF, education, and sex; WMS-IV indices predicted from ToPF and age. Difference = Achieved − Predicted. Base rate = % of standardisation sample with discrepancy ≤ this value (negative discrepancies only).</div>
  `;
}

// === OPIE-4 vs WAIS-IV calculation ===
// Builds a list of OPIE-4 prediction rows for whichever combinations of inputs are available,
// then renders them into the third premorbid tab. Achieved values are kept in preState.opieAchieved
// keyed by a stable model key (e.g. "FSIQ_VC_MR").
function getOpiePredictions(){
  const vc  = preNum('pre-vc');
  const mr  = preNum('pre-mr');
  const age = preNum('pre-age');
  const sex = preStr('pre-sex');
  const sexC_opie = sex === 'Female' ? 0 : sex === 'Male' ? 1 : null;

  const sexEffect = sexC_opie != null ? sexC_opie : 0;
  const hasVC = vc != null, hasMR = mr != null;

  function predict(c){
    let pred = c.intercept + (c.age != null ? c.age * age : 0) + c.age3 * Math.pow(age, 3) + c.sex * sexEffect;
    if (c.vc != null && hasVC) pred += c.vc * vc;
    if (c.mr != null && hasMR) pred += c.mr * mr;
    return pred;
  }

  const list = [];
  const canUseAge = age != null;
  function pushModel(key, label, c, needVC, needMR, tipKey){
    const canPredict = canUseAge && (!needVC || hasVC) && (!needMR || hasMR);
    list.push({ key, label, val: canPredict ? predict(c) : null, see: c.see, r: c.r, tipKey });
  }

  // FSIQ models
  pushModel('FSIQ_VC_MR', 'FSIQ - Vocab + Matrix', OPIE_PRORATED_FSIQ.VC_MR, true, true, 'opiePredFSIQ_VCMR');
  pushModel('FSIQ_VC', 'FSIQ - Vocab only', OPIE_PRORATED_FSIQ.VC, true, false, 'opiePredFSIQ_VC');
  pushModel('FSIQ_MR', 'FSIQ - Matrix only', OPIE_PRORATED_FSIQ.MR, false, true, 'opiePredFSIQ_MR');

  // GAI models
  pushModel('GAI_VC_MR', 'GAI - Vocab + Matrix', OPIE_PRORATED_GAI.VC_MR, true, true, 'opiePredGAI_VCMR');
  pushModel('GAI_VC', 'GAI - Vocab only', OPIE_PRORATED_GAI.VC, true, false, 'opiePredGAI_VC');

  return list;
}

function opieBaseRateFor(rowKey, diff){
  if (!diff) return '-';
  const key = diff > 0 ? `+${diff}` : String(diff);
  const row = OPIE_BASE_RATES[key];
  if (row && row[rowKey] != null) return fmtPctBr(row[rowKey]);
  // Beyond the range of the published table - for negative discrepancies this means
  // the base rate is smaller than the minimum tabulated value (< 0.1%)
  if (diff < 0) return '< 0.1%';
  return '-';
}

function calcOpiePredict(){
  const tbody = preGet('pre-opiepredict-tbody');
  if (!tbody) return;

  const ci = preStr('pre-ci') || '0.90';
  const mult = preCiMult();
  const ciPct = ci === '0.95' ? '95%' : '90%';
  document.querySelectorAll('.pre-ci-lo-hdr-3').forEach(e => e.textContent = `Lower ${ciPct}`);
  document.querySelectorAll('.pre-ci-hi-hdr-3').forEach(e => e.textContent = `Upper ${ciPct}`);

  const rows = getOpiePredictions();
  preState.opieRows = rows;
  preState.ciPct = ciPct;

  // Group: FSIQ rows, then GAI rows
  const fsiq = rows.filter(r => r.key.startsWith('FSIQ_'));
  const gai  = rows.filter(r => r.key.startsWith('GAI_'));

  function modelCell(row){
    const tip = PRE_MODEL_TOOLTIPS[row.tipKey] || '';
    return `<td class="model model-has-tip">${escapeHtml(row.label)}<span class="model-info-dot" data-tooltip="${escapeAttr(tip)}" aria-label="${escapeAttr(row.label + '. ' + tip)}">?</span></td>`;
  }

  function rowHtml(row){
    const ach = preState.opieAchieved[row.key];
    const achVal = (ach != null && ach !== '') ? parseFloat(ach) : null;
    const hasPred = row.val != null && Number.isFinite(row.val);
    const lo = hasPred ? fmtIntOrDash(row.val - mult * row.see) : '-';
    const hi = hasPred ? fmtIntOrDash(row.val + mult * row.see) : '-';
    let diffHtml = '<td class="num diff">-</td>';
    let brHtml = '<td class="num opie-base-rate">-</td>';
    if (hasPred && achVal != null && !isNaN(achVal)){
      const diff = Math.round(achVal - row.val);
      const sign = diff > 0 ? '+' : '';
      const cls = diff < 0 ? 'neg' : (diff > 0 ? 'pos' : '');
      diffHtml = `<td class="num diff ${cls}">${sign}${diff}</td>`;
      brHtml = `<td class="num opie-base-rate">${opieBaseRateFor(row.key, diff)}</td>`;
    }
    const achInputVal = ach != null ? escapeAttr(ach) : '';
    return `<tr>
      ${modelCell(row)}
      <td class="num">${hasPred ? fmtIntOrDash(row.val) : '-'}</td>
      <td class="num">${lo}</td>
      <td class="num">${hi}</td>
      <td class="achieved-cell"><input type="number" min="40" max="160" step="1" data-pre-opie-ach="${row.key}" value="${achInputVal}"></td>
      ${diffHtml}
      ${brHtml}
    </tr>`;
  }

  let html = '';
  if (fsiq.length){
    html += '<tr class="group-row"><td colspan="7">FSIQ predictions</td></tr>';
    html += fsiq.map(rowHtml).join('');
  }
  if (gai.length){
    html += '<tr class="group-row"><td colspan="7">GAI predictions</td></tr>';
    html += gai.map(rowHtml).join('');
  }
  tbody.innerHTML = html;

  // Wire the achieved inputs (re-wired each render because rows are rebuilt)
  tbody.querySelectorAll('[data-pre-opie-ach]').forEach(inp => {
    inp.addEventListener('input', e => {
      preState.opieAchieved[e.target.dataset.preOpieAch] = e.target.value;
      const tr = e.target.closest('tr');
      if (!tr) return;
      const row = (preState.opieRows || []).find(r => r.key === e.target.dataset.preOpieAch);
      if (!row) return;
      const diffCell = tr.querySelector('.diff');
      const brCell = tr.querySelector('.opie-base-rate');
      const achVal = e.target.value !== '' ? parseFloat(e.target.value) : null;
      diffCell.className = 'num diff';
      if (achVal == null || isNaN(achVal)){
        diffCell.textContent = '-';
        if (brCell) brCell.textContent = '-';
      } else {
        if (row.val == null || !Number.isFinite(row.val)){
          diffCell.textContent = '-';
          if (brCell) brCell.textContent = '-';
        } else {
          const diff = Math.round(achVal - row.val);
          diffCell.textContent = (diff > 0 ? '+' : '') + diff;
          if (diff < 0) diffCell.classList.add('neg');
          else if (diff > 0) diffCell.classList.add('pos');
          if (brCell) brCell.textContent = opieBaseRateFor(row.key, diff);
        }
      }
      renderOpiePredictApa();
    });
  });

  renderOpiePredictApa();
}

// === APA Output: OPIE-4 vs WAIS-IV ===
function renderOpiePredictApa(){
  const out = preGet('pre-opiepredict-apa');
  if (!out) return;
  const ciPct = preState.ciPct || '90%';
  const mult = preCiMult();
  const rows = preState.opieRows || [];

  const byGroup = { FSIQ:[], GAI:[] };
  rows.forEach(r => {
    if (r.key.startsWith('FSIQ_')) byGroup.FSIQ.push(r);
    else if (r.key.startsWith('GAI_')) byGroup.GAI.push(r);
  });

  let body = '';
  ['FSIQ','GAI'].forEach(group => {
    if (!byGroup[group].length) return;
    body += `<tr><td colspan="7" style="font-style:italic;padding-top:6px;padding-bottom:2px">Prorated ${group}</td></tr>`;
    byGroup[group].forEach(r => {
      const achRaw = preState.opieAchieved[r.key];
      const ach = achRaw != null && achRaw !== '' && !isNaN(parseFloat(achRaw)) ? parseFloat(achRaw) : null;
      const hasPred = r.val != null && Number.isFinite(r.val);
      const lo = hasPred ? Math.round(r.val - mult * r.see) : '-';
      const hi = hasPred ? Math.round(r.val + mult * r.see) : '-';
      const diff = (!hasPred || ach == null) ? null : Math.round(ach - r.val);
      const diffText = diff == null ? '-' : `${diff > 0 ? '+' : ''}${diff}`;
      const br = diff == null ? '-' : opieBaseRateFor(r.key, diff);
      const achText = ach == null ? '-' : String(ach);
      const predText = hasPred ? String(Math.round(r.val)) : '-';
      body += `<tr><td>&nbsp;&nbsp;${escapeHtml(r.label)}</td><td class="num">${predText}</td><td class="num">${lo}</td><td class="num">${hi}</td><td class="num">${achText}</td><td class="num">${diffText}</td><td class="num">${br}</td></tr>`;
    });
  });

  out.innerHTML = `
    <div class="apa-table-num">Table 3</div>
    <div class="apa-table-title">OPIE-4-predicted versus achieved prorated index scores</div>
    <table class="apa-table">
      <thead>
        <tr><th>Model</th><th class="num">Predicted score</th><th class="num">Lower ${ciPct}</th><th class="num">Upper ${ciPct}</th><th class="num">Achieved score</th><th class="num">Difference</th><th class="num">Base Rate</th></tr>
      </thead>
      <tbody>${body}</tbody>
    </table>
    <div class="apa-note"><strong>Note.</strong> OPIE-4 prorated FSIQ and GAI predicted from age plus Vocabulary and/or Matrix Reasoning. Difference = Achieved − Predicted. Base rate = % of standardisation sample with discrepancy ≤ this value. Based on US norms.</div>
  `;
}

// Wire all premorbid input listeners
function setupPremorbidListeners(){
  ['pre-topf','pre-vc','pre-mr','pre-sex','pre-edu','pre-age','pre-occ','pre-ci'].forEach(id => {
    const el = preGet(id);
    if (!el) return;
    el.addEventListener('input', () => { calcPremorbid(); calcPredict(); calcOpiePredict(); });
    el.addEventListener('change', () => { calcPremorbid(); calcPredict(); calcOpiePredict(); });
  });
  const titleEl = preGet('pre-title');
  if (titleEl){
    titleEl.addEventListener('input', () => { renderPreEstimatesApa(); renderPrePredictApa(); renderOpiePredictApa(); });
  }
}


function escapeHtml(s){ return String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'})[c]); }
function escapeAttr(s){ return escapeHtml(s); }


function applyCalculatorPolish(){
  ['sdi','rci-basic','rci-practice','rci-srb','rci-crawford'].forEach(id => {
    const sec = document.getElementById(id);
    if (!sec) return;
    const firstPanel = sec.querySelector('.panel:not(.autofill-panel)');
    if (firstPanel && !firstPanel.classList.contains('settings-panel')){
      firstPanel.classList.add('settings-panel');
      if (!firstPanel.querySelector('.panel-kicker')){
        const kicker = document.createElement('div');
        kicker.className = 'panel-kicker';
        kicker.textContent = 'Analysis settings';
        firstPanel.insertBefore(kicker, firstPanel.firstChild);
      }
    }
    const table = sec.querySelector('.input-table');
    if (table) table.classList.add('calculator-table');
  });
  applyRciGroupedHeaders();
}
function applyRciGroupedHeaders(){
  const specs = {
    'rci-basic-table':    [{t:'',c:2},{t:'Norms',c:2},{t:'Scores',c:2},{t:'Results',c:3},{t:'',c:1}],
    'rci-practice-table': [{t:'',c:2},{t:'Norms',c:5},{t:'Scores',c:2},{t:'Results',c:3},{t:'',c:1}],
    'rci-srb-table':      [{t:'',c:2},{t:'Norms',c:5},{t:'Scores',c:2},{t:'Results',c:4},{t:'',c:1}],
    'rci-crawford-table': [{t:'',c:2},{t:'Norms',c:6},{t:'Scores',c:2},{t:'Results',c:4},{t:'',c:1}]
  };
  Object.entries(specs).forEach(([id, groups]) => {
    const table = document.getElementById(id);
    if (!table || !table.tHead || table.tHead.querySelector('.table-group-row')) return;
    const row = document.createElement('tr');
    row.className = 'table-group-row';
    row.innerHTML = groups.map(g => `<th colspan="${g.c}">${g.t}</th>`).join('');
    table.tHead.insertBefore(row, table.tHead.firstChild);
  });
  // After group rows exist, set up the visual lock toggle for each table
  // (skip rci-basic - only 2 norm cells, the toggle adds clutter without much benefit)
  Object.keys(specs)
    .filter(id => id !== 'rci-basic-table')
    .forEach(id => setupNormsLockToggle(id));
}

/* ---- Norms visual-lock toggle ----
   Aesthetic only: when "locked", norm columns visually fade and their
   inputs are disabled. Layout is untouched - no display:none, no width
   changes. The user can toggle via a small pill button next to the
   "Test data & patient scores" heading. */
function setupNormsLockToggle(tableId){
  const table = document.getElementById(tableId);
  if (!table || !table.tHead) return;
  const groupRow = table.tHead.querySelector('.table-group-row');
  if (!groupRow) return;
  // Compute the column range covered by the "Norms" group cell
  let startCol = 0, endCol = 0, normsCell = null;
  for (const c of groupRow.children){
    const span = parseInt(c.colSpan || 1, 10);
    if ((c.textContent || '').trim() === 'Norms'){
      normsCell = c;
      endCol = startCol + span;
      break;
    }
    startCol += span;
  }
  if (!normsCell) return;
  // Tag norm cells (header + body) so CSS can target them
  function tagAll(){
    Array.from(table.querySelectorAll('thead tr:not(.table-group-row), tbody tr')).forEach(row => {
      Array.from(row.children).forEach((cell, i) => {
        if (i >= startCol && i < endCol) cell.dataset.normCell = 'true';
        else delete cell.dataset.normCell;
      });
    });
    normsCell.dataset.normCell = 'true';
    // If currently locked, ensure newly-added inputs are also disabled
    if (table.classList.contains('norms-locked')){
      table.querySelectorAll('[data-norm-cell="true"] input').forEach(inp => { inp.disabled = true; });
    }
  }
  tagAll();
  // Re-tag when tbody is re-rendered (add row, autofill, etc.)
  if (!table.dataset.normsLockObserver){
    new MutationObserver(() => tagAll()).observe(table.tBodies[0] || table, { childList: true, subtree: true });
    table.dataset.normsLockObserver = '1';
  }
  // Inject a two-segment "Lock | Unlock" pill INSIDE the Norms group cell
  if (!normsCell.querySelector('.norms-toggle-pill')){
    const labelText = (normsCell.textContent || 'Norms').trim();
    normsCell.innerHTML = `
      <span class="norms-group-label">${labelText}</span>
      <div class="norms-toggle-pill" role="group" aria-label="Norm column lock state" data-table-id="${tableId}">
        <button type="button" class="norms-pill-segment is-active" data-state="locked" aria-pressed="true">
          <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
            <rect x="3" y="7.5" width="10" height="6.5" rx="1"/>
            <path d="M5 7.5V5.5a3 3 0 0 1 6 0v2"/>
          </svg>
          <span>Lock</span>
        </button>
        <button type="button" class="norms-pill-segment" data-state="unlocked" aria-pressed="false">
          <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
            <rect x="3" y="7.5" width="10" height="6.5" rx="1"/>
            <path d="M5 7.5V5.5a3 3 0 0 1 6 0"/>
          </svg>
          <span>Unlock</span>
        </button>
      </div>
    `;
    normsCell.querySelectorAll('.norms-pill-segment').forEach(btn => {
      btn.addEventListener('click', e => {
        e.stopPropagation();
        toggleNormsLock(tableId, btn.dataset.state === 'locked');
      });
    });
  }
  // Default state: locked (norms faded, inputs disabled)
  if (!table.dataset.normsLockInitialized){
    toggleNormsLock(tableId, true);
    table.dataset.normsLockInitialized = '1';
  }
}

function toggleNormsLock(tableId, force){
  const table = document.getElementById(tableId);
  if (!table) return;
  const willLock = (typeof force === 'boolean') ? force : !table.classList.contains('norms-locked');
  table.classList.toggle('norms-locked', willLock);
  // Disable/enable inputs in norm columns
  table.querySelectorAll('[data-norm-cell="true"] input').forEach(inp => { inp.disabled = willLock; });
  // Sync segmented pill state - highlight the segment matching current state
  const pill = table.querySelector('.norms-toggle-pill[data-table-id="' + tableId + '"]');
  if (pill){
    pill.querySelectorAll('.norms-pill-segment').forEach(seg => {
      const isActive = (seg.dataset.state === 'locked') === willLock;
      seg.classList.toggle('is-active', isActive);
      seg.setAttribute('aria-pressed', String(isActive));
    });
  }
}


function hasAnyRowValue(row){
  return Object.keys(row || {}).some(k => String(row[k] ?? '').trim() !== '');
}
function numericProblem(row, method){
  const req = method === 'rci-basic' ? ['sd','r','t1','t2'] : method === 'rci-crawford' ? ['m1','sd1','m2','sd2','r','n','t1','t2'] : ['m1','sd1','m2','sd2','r','t1','t2'];
  for (const k of req){
    if (String(row[k] ?? '').trim() === '') return 'Awaiting values';
    const v = Number(row[k]);
    if (!Number.isFinite(v)) return 'Check values';
  }
  const r = Number(row.r);
  if (r < 0 || r >= 1) return 'Reliability must be 0–.99';
  if (method === 'rci-basic' && Number(row.sd) <= 0) return 'SD must be > 0';
  if (method !== 'rci-basic' && (Number(row.sd1) <= 0 || Number(row.sd2) <= 0)) return 'SD must be > 0';
  if (method === 'rci-crawford' && Number(row.n) < 3) return 'N must be ≥ 3';
  return 'Check values';
}
function sdiProblem(row){
  const raw = sdiMode && sdiMode() === 'raw';
  const req = raw ? ['t1','t2','sd'] : ['t1','t2'];
  for (const k of req){
    if (String(row[k] ?? '').trim() === '') return 'Awaiting values';
    const v = Number(row[k]);
    if (!Number.isFinite(v)) return 'Check values';
  }
  if (raw && Number(row.sd) <= 0) return 'SD must be > 0';
  return 'Check values';
}
function clearOutcomeStatus(td){
  td.classList.remove('status-awaiting','status-check');
}
function setOutcomeStatus(td, label, mode){
  td.textContent = label;
  td.classList.add(mode === 'check' ? 'status-check' : 'status-awaiting');
}
function insertStepBefore(el, num, title, copy){
  // Workflow step labels removed to keep calculator pages visually lighter.
  return;
}

function buildColumnGuide(method){
  return '';
}
function enhanceCalculatorWorkflow(){
  ['sdi','rci-basic','rci-practice','rci-srb','rci-crawford'].forEach((id) => {
    const sec = document.getElementById(id);
    if (!sec) return;
    const firstPanel = sec.querySelector('.panel');
    insertStepBefore(firstPanel, 1, 'Configure analysis', 'Set the confidence threshold and report labels before entering patient scores.');
    if (firstPanel) firstPanel.classList.add('settings-panel');
    const tableTitle = Array.from(sec.querySelectorAll('h2.block-title')).find(h => /Test data|patient scores/i.test(h.textContent));
    if (tableTitle){
      insertStepBefore(tableTitle, 2, 'Enter scores', 'Use auto-fill when normative values are available, or enter values manually.');
    }
    const apa = sec.querySelector('.apa-wrap');
    insertStepBefore(apa, 3, 'Review and copy output', 'The report-ready APA table updates once at least one valid row is complete.');
    const add = sec.querySelector('.add-row-btn');
    if (add){
      let wrap = add.parentElement;
      if (!wrap.classList.contains('add-row-actions')){
        wrap = document.createElement('div');
        wrap.className = 'add-row-actions';
        add.parentNode.insertBefore(wrap, add);
        wrap.appendChild(add);
      }
      if (!wrap.querySelector(`[data-example="${id}"]`)){
        const ex = document.createElement('button');
        ex.type = 'button'; ex.className = 'add-row-btn btn-example';
        ex.dataset.example = id;
        ex.textContent = '+ Load example row';
        const clear = wrap.querySelector('.btn-clear');
        wrap.insertBefore(ex, clear || null);
      }
    }
  });
  document.querySelectorAll('[data-example]').forEach(btn => {
    if (btn.dataset.wired) return;
    btn.dataset.wired = '1';
    btn.addEventListener('click', () => loadExampleRow(btn.dataset.example));
  });
  addColumnTitles();
}
function addColumnTitles(){
  const titleMap = {
    'M₁':'Normative mean at baseline / Date 1', 'SD₁':'Normative standard deviation at baseline / Date 1',
    'M₂':'Normative mean at retest / Date 2', 'SD₂':'Normative standard deviation at retest / Date 2',
    'r':'Test–retest reliability coefficient', 'R':'Test–retest reliability coefficient', 'N':'Normative sample size',
    'Ŷ₂':'Predicted Date 2 score', 't(RB)':'Regression-based t statistic', 'RCI':'Reliable Change Index', 'RCI (z)':'Reliable Change Index as a z statistic',
    'p':'Two-tailed p value', 'Outcome':'Clinical interpretation of reliable change', 'SD':'Standard deviation', 'SD Δ':'Standard-deviation change'
  };
  document.querySelectorAll('.input-table th').forEach(th => {
    const txt = th.textContent.trim().replace(/\s+/g,' ');
    const key = txt.replace(/[()]/g, m=>m);
    const plain = txt.replace(/[^A-Za-z0-9₁₂ŶΔ]/g,'');
    if (titleMap[txt]) th.title = titleMap[txt];
    else if (txt === 'RCI (z)') th.title = titleMap['RCI (z)'];
    else if (txt.toLowerCase() === 'p') th.title = titleMap.p;
    else if (txt.includes('t') && txt.includes('RB')) th.title = titleMap['t(RB)'];
  });
}
function loadExampleRow(method){
  if (method === 'sdi'){
    const example = sdiMode() === 'raw' ? {name:'Example memory score',t1:'42',t2:'36',sd:'8'} : {name:'Example memory score',t1:'9',t2:'6'};
    sdiRows.push(example); renderSdi(); showToast('Example row added'); return;
  }
  if (method === 'battery'){
    const example = {name:'Example subtest', raw:'25', score:'10'};
    batteryRows.push(example); renderBattery(); showToast('Example row added'); return;
  }
  const examples = {
    'rci-basic': {name:'Example index score',sd:'15',r:'0.90',t1:'100',t2:'89'},
    'rci-practice': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89'},
    'rci-srb': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89'},
    'rci-crawford': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',n:'100',t1:'100',t2:'89'}
  };
  if (rciState[method]){ rciState[method].rows.push(examples[method]); renderRci(method); showToast('Example row added'); }
}

/* ============================================================
   REPORT WRITER · descriptive narrative builder
   ============================================================ */
const REPORT_DOMAINS = [
  { id:'attention_working_memory', label:'Attention and Working Memory', ability:'attention and working memory', proseLabel:'Attention and Working Memory' },
  { id:'processing_speed', label:'Processing Speed', ability:'processing speed', proseLabel:'Processing Speed' },
  { id:'verbal_learning_memory', label:'Verbal Learning and Memory', ability:'verbal learning and memory', proseLabel:'Verbal Learning and Memory' },
  { id:'visual_learning_memory', label:'Visual Learning and Memory', ability:'visual learning and memory', proseLabel:'Visual Learning and Memory' },
  { id:'language', label:'Language', ability:'language', proseLabel:'Language ability' },
  { id:'visuospatial', label:'Visuospatial / Constructional Skills', ability:'visuospatial and constructional ability', proseLabel:'Visuospatial and Constructional ability' },
  { id:'executive', label:'Executive Functioning', ability:'executive functioning', proseLabel:'Executive Functioning' },
  { id:'intellectual', label:'Intellectual Functioning', ability:'intellectual functioning', proseLabel:'Intellectual Functioning' },
  { id:'mood', label:'Mood / Symptoms', ability:'mood and symptom reporting', proseLabel:'mood and symptom reporting' },
  { id:'other', label:'Other / Custom', ability:'the selected domain', proseLabel:'this domain' }
];
const REPORT_SCORE_TYPES = [
  { id:'standard', label:'Standard Score', short:'SS', sentence:'standard score' },
  { id:'scaled', label:'Scaled Score', short:'ScS', sentence:'scaled score' },
  { id:'t', label:'T-Score', short:'T', sentence:'T-score' },
  { id:'z', label:'Z-Score', short:'z', sentence:'z-score' },
  { id:'percentile', label:'Percentile', short:'%ile', sentence:'percentile' }
];
const REPORT_CONSTRUCT_LABELS = {
  common:'Commonly used',
  auditory_attention:'Auditory attention',
  working_memory:'Working memory / mental manipulation',
  sequencing_attention:'Sequencing / attention',
  graphomotor_speed:'Visual scanning and speed',
  naming_reading_speed:'Rapid naming / reading speed',
  verbal_learning:'Verbal learning',
  verbal_recall:'Verbal recall',
  verbal_recognition:'Verbal recognition',
  visual_learning:'Visual learning',
  visual_recall:'Visual recall',
  visual_recognition:'Visual recognition',
  word_knowledge:'Word knowledge',
  verbal_reasoning:'Verbal reasoning',
  verbal_fluency:'Verbal fluency / generativity',
  visuoconstruction:'Visuoconstruction',
  visuospatial_reasoning:'Visuospatial reasoning',
  inhibition:'Inhibition',
  switching:'Cognitive flexibility / switching',
  planning_problem_solving:'Planning / problem-solving',
  abstract_reasoning:'Abstract reasoning / concept formation',
  global_indices:'Global indices',
  mood_symptoms:'Mood / symptom measures',
  other:'Other mapped measures'
};
const REPORT_DOMAIN_GROUP_ORDER = {
  attention_working_memory:['auditory_attention','working_memory','sequencing_attention','other'],
  processing_speed:['graphomotor_speed','naming_reading_speed','sequencing_attention','other'],
  verbal_learning_memory:['verbal_learning','verbal_recall','verbal_recognition','other'],
  visual_learning_memory:['visual_learning','visual_recall','visual_recognition','other'],
  language:['word_knowledge','verbal_reasoning','verbal_fluency','other'],
  visuospatial:['visuoconstruction','visuospatial_reasoning','graphomotor_speed','other'],
  executive:['inhibition','switching','planning_problem_solving','verbal_fluency','abstract_reasoning','other'],
  intellectual:['global_indices','word_knowledge','visuospatial_reasoning','working_memory','graphomotor_speed','other'],
  mood:['mood_symptoms','other'],
  other:['other']
};
const REPORT_FAMILY_ORDER = [
  'WAIS-IV Indices','WAIS-IV Core Subtests','WISC-V Indices','WISC-V Subtests','WMS-IV Indices','WMS-IV Subtests',
  'CVLT-3 Indices','CVLT-3 Trials','RBANS Indices','RBANS Subtests',
  'D-KEFS Trail Making Test','D-KEFS Colour-Word Interference','D-KEFS Verbal Fluency','D-KEFS Design Fluency',
  'D-KEFS Sorting Test','D-KEFS Tower Test','D-KEFS Word Context Test','D-KEFS Word Proverb Test'
];
const REPORT_PICKER_LIMIT = 160;
const REPORT_STORAGE_KEY = 'reportWriterState_v2';

// Reporting-role bands shown within each test in the picker. Order matters - indices first.
const REPORT_ROLE_ORDER = ['index','subtest','learning-trial','delayed-recall','recognition','intrusion','condition','measure'];
const REPORT_ROLE_LABELS = {
  'index':'Indices',
  'subtest':'Subtests',
  'learning-trial':'Encoding / Learning',
  'delayed-recall':'Delayed Recall',
  'recognition':'Recognition',
  'intrusion':'Intrusions / Monitoring',
  'condition':'Conditions',
  'measure':'Measures'
};

// Test catalog - clinical administration units. Each entry maps families → measures with
// reporting role + domain mapping + curated "core" measures per cognitive domain.
// Measures not in the catalog fall back to construct-cluster inference (REPORT_CLUSTER_RULES).
const REPORT_TEST_CATALOG = [
  {
    id:'wais-iv', name:'WAIS-IV', longName:'Wechsler Adult Intelligence Scale-Fourth Edition',
    families:['WAIS-IV Indices','WAIS-IV Core Subtests'],
    measures:{
      'Full Scale IQ':              { role:'index', cluster:'global_indices', construct:'overall intellectual functioning', domains:['intellectual'] },
      'Verbal Comprehension Index': { role:'index', cluster:'global_indices', construct:'verbal comprehension', domains:['intellectual','language'] },
      'Perceptual Reasoning Index': { role:'index', cluster:'global_indices', construct:'perceptual reasoning', domains:['intellectual','visuospatial'] },
      'Working Memory Index':       { role:'index', cluster:'global_indices', construct:'working memory and mental manipulation', domains:['intellectual','attention_working_memory'] },
      'Processing Speed Index':     { role:'index', cluster:'global_indices', construct:'processing speed', domains:['intellectual','processing_speed'] },
      'Block Design':         { role:'subtest', cluster:'visuoconstruction', construct:'visuospatial reasoning and constructional ability', domains:['visuospatial','intellectual'] },
      'Similarities':         { role:'subtest', cluster:'verbal_reasoning', construct:'verbal reasoning and concept formation', domains:['language','executive','intellectual'] },
      'Digit Span':           { role:'subtest', cluster:'working_memory', construct:'working memory and mental manipulation', domains:['attention_working_memory','intellectual'] },
      'Matrix Reasoning':     { role:'subtest', cluster:'visuospatial_reasoning', construct:'visuospatial reasoning', domains:['visuospatial','intellectual'] },
      'Vocabulary':           { role:'subtest', cluster:'word_knowledge', construct:'word knowledge and crystallised verbal ability', domains:['language','intellectual'] },
      'Arithmetic':           { role:'subtest', cluster:'working_memory', construct:'working memory and mental manipulation', domains:['attention_working_memory','intellectual'] },
      'Symbol Search':        { role:'subtest', cluster:'graphomotor_speed', construct:'graphomotor and visual scanning speed', domains:['processing_speed','intellectual'] },
      'Visual Puzzles':       { role:'subtest', cluster:'visuospatial_reasoning', construct:'visuospatial reasoning', domains:['visuospatial','intellectual'] },
      'Information':          { role:'subtest', cluster:'word_knowledge', construct:'word knowledge and crystallised verbal ability', domains:['language','intellectual'] },
      'Coding':               { role:'subtest', cluster:'graphomotor_speed', construct:'graphomotor and visual scanning speed', domains:['processing_speed','intellectual'] }
    },
    coreByDomain:{
      attention_working_memory:['Working Memory Index','Digit Span','Arithmetic'],
      processing_speed:['Processing Speed Index','Coding','Symbol Search'],
      visuospatial:['Perceptual Reasoning Index','Block Design','Matrix Reasoning','Visual Puzzles'],
      language:['Verbal Comprehension Index','Vocabulary','Similarities','Information'],
      intellectual:['Full Scale IQ','Verbal Comprehension Index','Perceptual Reasoning Index','Working Memory Index','Processing Speed Index']
    },
    indexComposition:{
      'Full Scale IQ':['Verbal Comprehension Index','Perceptual Reasoning Index','Working Memory Index','Processing Speed Index'],
      'Verbal Comprehension Index':['Vocabulary','Similarities','Information'],
      'Perceptual Reasoning Index':['Block Design','Matrix Reasoning','Visual Puzzles'],
      'Working Memory Index':['Digit Span','Arithmetic'],
      'Processing Speed Index':['Coding','Symbol Search']
    }
  },
  {
    id:'wisc-v', name:'WISC-V', longName:'Wechsler Intelligence Scale for Children-Fifth Edition',
    families:['WISC-V Indices','WISC-V Subtests'],
    measures:{
      'Full Scale IQ':              { role:'index', cluster:'global_indices', construct:'overall intellectual functioning', domains:['intellectual'] },
      'Verbal Comprehension Index': { role:'index', cluster:'global_indices', construct:'verbal comprehension', domains:['intellectual','language'] },
      'Visuospatial Index':         { role:'index', cluster:'global_indices', construct:'visuospatial reasoning', domains:['intellectual','visuospatial'] },
      'Fluid Reasoning Index':      { role:'index', cluster:'global_indices', construct:'fluid reasoning', domains:['intellectual','visuospatial'] },
      'Working Memory Index':       { role:'index', cluster:'global_indices', construct:'working memory and mental manipulation', domains:['intellectual','attention_working_memory'] },
      'Processing Speed Index':     { role:'index', cluster:'global_indices', construct:'processing speed', domains:['intellectual','processing_speed'] },
      'Similarities':         { role:'subtest', cluster:'verbal_reasoning', construct:'verbal reasoning and concept formation', domains:['language','executive','intellectual'] },
      'Vocabulary':           { role:'subtest', cluster:'word_knowledge', construct:'word knowledge and crystallised verbal ability', domains:['language','intellectual'] },
      'Information':          { role:'subtest', cluster:'word_knowledge', construct:'word knowledge and crystallised verbal ability', domains:['language','intellectual'] },
      'Comprehension':        { role:'subtest', cluster:'verbal_reasoning', construct:'verbal reasoning and social/practical knowledge', domains:['language','intellectual'] },
      'Block Design':         { role:'subtest', cluster:'visuoconstruction', construct:'visuospatial reasoning and constructional ability', domains:['visuospatial','intellectual'] },
      'Visual Puzzles':       { role:'subtest', cluster:'visuospatial_reasoning', construct:'visuospatial reasoning', domains:['visuospatial','intellectual'] },
      'Matrix Reasoning':     { role:'subtest', cluster:'visuospatial_reasoning', construct:'visuospatial reasoning', domains:['visuospatial','intellectual'] },
      'Figure Weights':       { role:'subtest', cluster:'visuospatial_reasoning', construct:'quantitative and analogical reasoning', domains:['visuospatial','intellectual'] },
      'Picture Completion':   { role:'subtest', cluster:'visuospatial_reasoning', construct:'visual perception and detail recognition', domains:['visuospatial','intellectual'] },
      'Arithmetic':           { role:'subtest', cluster:'working_memory', construct:'working memory and mental manipulation', domains:['attention_working_memory','intellectual'] },
      'Digit Span':           { role:'subtest', cluster:'working_memory', construct:'working memory and mental manipulation', domains:['attention_working_memory','intellectual'] },
      'Picture Span':         { role:'subtest', cluster:'working_memory', construct:'visual working memory', domains:['attention_working_memory','intellectual'] },
      'Letter-Number Sequencing':{ role:'subtest', cluster:'working_memory', construct:'working memory and mental manipulation', domains:['attention_working_memory','intellectual'] },
      'Coding':               { role:'subtest', cluster:'graphomotor_speed', construct:'graphomotor and visual scanning speed', domains:['processing_speed','intellectual'] },
      'Symbol Search':        { role:'subtest', cluster:'graphomotor_speed', construct:'graphomotor and visual scanning speed', domains:['processing_speed','intellectual'] },
      'Cancellation':         { role:'subtest', cluster:'graphomotor_speed', construct:'visual scanning and cancellation speed', domains:['processing_speed','intellectual'] }
    },
    coreByDomain:{
      attention_working_memory:['Working Memory Index','Digit Span','Picture Span'],
      processing_speed:['Processing Speed Index','Coding','Symbol Search'],
      visuospatial:['Visuospatial Index','Block Design','Visual Puzzles'],
      language:['Verbal Comprehension Index','Vocabulary','Similarities'],
      intellectual:['Full Scale IQ','Verbal Comprehension Index','Visuospatial Index','Fluid Reasoning Index','Working Memory Index','Processing Speed Index']
    },
    indexComposition:{
      'Full Scale IQ':['Verbal Comprehension Index','Visuospatial Index','Fluid Reasoning Index','Working Memory Index','Processing Speed Index'],
      'Verbal Comprehension Index':['Similarities','Vocabulary'],
      'Visuospatial Index':['Block Design','Visual Puzzles'],
      'Fluid Reasoning Index':['Matrix Reasoning','Figure Weights'],
      'Working Memory Index':['Digit Span','Picture Span'],
      'Processing Speed Index':['Coding','Symbol Search']
    }
  },
  {
    id:'wms-iv', name:'WMS-IV', longName:'Wechsler Memory Scale-Fourth Edition',
    families:['WMS-IV Indices','WMS-IV Subtests'],
    measures:{
      'Auditory Memory Index':      { role:'index', cluster:'global_indices', construct:'auditory memory', domains:['verbal_learning_memory'] },
      'Visual Memory Index':        { role:'index', cluster:'global_indices', construct:'visual memory', domains:['visual_learning_memory'] },
      'Visual Working Memory Index':{ role:'index', cluster:'global_indices', construct:'visual working memory', domains:['attention_working_memory'] },
      'Immediate Memory Index':     { role:'index', cluster:'global_indices', construct:'immediate memory', domains:['verbal_learning_memory','visual_learning_memory'] },
      'Delayed Memory Index':       { role:'index', cluster:'global_indices', construct:'delayed memory', domains:['verbal_learning_memory','visual_learning_memory'] },
      'Logical Memory I':           { role:'subtest', cluster:'verbal_learning', construct:'contextual verbal memory (immediate)', domains:['verbal_learning_memory'] },
      'Logical Memory II':          { role:'subtest', cluster:'verbal_recall', construct:'contextual verbal memory (delayed)', domains:['verbal_learning_memory'] },
      'Verbal Paired Associates I': { role:'subtest', cluster:'verbal_learning', construct:'verbal paired-associate learning', domains:['verbal_learning_memory'] },
      'Verbal Paired Associates II':{ role:'subtest', cluster:'verbal_recall', construct:'delayed verbal paired-associate recall', domains:['verbal_learning_memory'] },
      'Verbal Paired Associates II - Word Recall':{ role:'subtest', cluster:'verbal_recall', construct:'delayed verbal paired-associate word recall', domains:['verbal_learning_memory'] },
      'Designs I':                  { role:'subtest', cluster:'visual_learning', construct:'visual design memory (immediate)', domains:['visual_learning_memory'] },
      'Designs I - Content':        { role:'subtest', cluster:'visual_learning', construct:'visual design content memory (immediate)', domains:['visual_learning_memory'] },
      'Designs I - Spatial':        { role:'subtest', cluster:'visual_learning', construct:'visual design spatial memory (immediate)', domains:['visual_learning_memory'] },
      'Designs II':                 { role:'subtest', cluster:'visual_recall', construct:'delayed visual design memory', domains:['visual_learning_memory'] },
      'Designs II - Content':       { role:'subtest', cluster:'visual_recall', construct:'delayed visual design content memory', domains:['visual_learning_memory'] },
      'Designs II - Spatial':       { role:'subtest', cluster:'visual_recall', construct:'delayed visual design spatial memory', domains:['visual_learning_memory'] },
      'Visual Reproduction I':      { role:'subtest', cluster:'visual_learning', construct:'visual reproduction (immediate)', domains:['visual_learning_memory'] },
      'Visual Reproduction II':     { role:'subtest', cluster:'visual_recall', construct:'delayed visual reproduction', domains:['visual_learning_memory'] },
      'Symbol Span':                { role:'subtest', cluster:'working_memory', construct:'visual working memory', domains:['attention_working_memory'] }
    },
    coreByDomain:{
      verbal_learning_memory:['Auditory Memory Index','Logical Memory I','Logical Memory II','Verbal Paired Associates I','Verbal Paired Associates II'],
      visual_learning_memory:['Visual Memory Index','Visual Reproduction I','Visual Reproduction II','Designs I','Designs II'],
      attention_working_memory:['Visual Working Memory Index','Symbol Span']
    },
    indexComposition:{
      'Auditory Memory Index':['Logical Memory I','Logical Memory II','Verbal Paired Associates I','Verbal Paired Associates II'],
      'Visual Memory Index':['Designs I','Designs II','Visual Reproduction I','Visual Reproduction II'],
      'Visual Working Memory Index':['Symbol Span'],
      'Immediate Memory Index':['Logical Memory I','Verbal Paired Associates I','Designs I','Visual Reproduction I'],
      'Delayed Memory Index':['Logical Memory II','Verbal Paired Associates II','Designs II','Visual Reproduction II']
    }
  },
  {
    id:'cvlt-3', name:'CVLT-3', longName:'California Verbal Learning Test-Third Edition',
    families:['CVLT-3 Indices','CVLT-3 Trials'],
    measures:{
      'T1-5 Correct':           { role:'index', cluster:'global_indices', construct:'total verbal learning across trials 1–5', domains:['verbal_learning_memory'] },
      'Delayed Recall Correct': { role:'index', cluster:'global_indices', construct:'delayed verbal recall', domains:['verbal_learning_memory'] },
      'Total Recall Correct':   { role:'index', cluster:'global_indices', construct:'overall verbal recall', domains:['verbal_learning_memory'] },
      'Trial 1':                { role:'learning-trial', cluster:'verbal_learning', construct:'initial verbal encoding', domains:['verbal_learning_memory'] },
      'Trial 2':                { role:'learning-trial', cluster:'verbal_learning', construct:'verbal learning across trials', domains:['verbal_learning_memory'] },
      'Trial 3':                { role:'learning-trial', cluster:'verbal_learning', construct:'verbal learning across trials', domains:['verbal_learning_memory'] },
      'Trial 4':                { role:'learning-trial', cluster:'verbal_learning', construct:'verbal learning across trials', domains:['verbal_learning_memory'] },
      'Trial 5':                { role:'learning-trial', cluster:'verbal_learning', construct:'final verbal learning trial', domains:['verbal_learning_memory'] },
      'List B Correct':         { role:'learning-trial', cluster:'verbal_learning', construct:'interference list recall', domains:['verbal_learning_memory'] },
      'Short Delay Free Recall':{ role:'delayed-recall', cluster:'verbal_recall', construct:'short-delay free recall', domains:['verbal_learning_memory'] },
      'Short Delay Cued Recall':{ role:'delayed-recall', cluster:'verbal_recall', construct:'short-delay cued recall', domains:['verbal_learning_memory'] },
      'Long Delay Free Recall': { role:'delayed-recall', cluster:'verbal_recall', construct:'long-delay free recall', domains:['verbal_learning_memory'] },
      'Long Delay Cued Recall': { role:'delayed-recall', cluster:'verbal_recall', construct:'long-delay cued recall', domains:['verbal_learning_memory'] },
      'Recognition':            { role:'recognition', cluster:'verbal_recognition', construct:'verbal recognition memory', domains:['verbal_learning_memory'] },
      'Recognition False Positive':{ role:'recognition', cluster:'verbal_recognition', construct:'recognition false-positive errors', domains:['verbal_learning_memory'] },
      'Recognition Discrimination':{ role:'recognition', cluster:'verbal_recognition', construct:'recognition discriminability', domains:['verbal_learning_memory'] },
      'Discrimination Nonparametric':{ role:'recognition', cluster:'verbal_recognition', construct:'recognition discriminability (nonparametric)', domains:['verbal_learning_memory'] },
      'Total Intrusions':       { role:'intrusion', cluster:'verbal_recall', construct:'verbal recall intrusions', domains:['verbal_learning_memory'] },
      'Total Repetitions':      { role:'intrusion', cluster:'verbal_recall', construct:'verbal recall repetitions', domains:['verbal_learning_memory'] }
    },
    coreByDomain:{
      verbal_learning_memory:['T1-5 Correct','Long Delay Free Recall','Recognition Discrimination']
    },
    indexComposition:{
      'T1-5 Correct':['Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'],
      'Delayed Recall Correct':['Short Delay Free Recall','Short Delay Cued Recall','Long Delay Free Recall','Long Delay Cued Recall'],
      'Total Recall Correct':['T1-5 Correct','Delayed Recall Correct']
    }
  },
  {
    id:'rbans', name:'RBANS', longName:'Repeatable Battery for the Assessment of Neuropsychological Status',
    families:['RBANS Indices','RBANS Subtests'],
    measures:{
      'Total Scale':                  { role:'index', cluster:'global_indices', construct:'overall cognitive functioning', domains:['intellectual'] },
      'Immediate Memory':             { role:'index', cluster:'global_indices', construct:'immediate verbal memory', domains:['verbal_learning_memory'] },
      'Visuospatial/Constructional':  { role:'index', cluster:'global_indices', construct:'visuospatial and constructional ability', domains:['visuospatial'] },
      'Attention':                    { role:'index', cluster:'global_indices', construct:'attention', domains:['attention_working_memory','processing_speed'] },
      'Language':                     { role:'index', cluster:'global_indices', construct:'language', domains:['language'] },
      'Delayed Memory':               { role:'index', cluster:'global_indices', construct:'delayed verbal and visual memory', domains:['verbal_learning_memory','visual_learning_memory'] },
      'List Learning':       { role:'subtest', cluster:'verbal_learning', construct:'verbal list-learning', domains:['verbal_learning_memory'] },
      'Story Memory':        { role:'subtest', cluster:'verbal_learning', construct:'contextual verbal memory (immediate)', domains:['verbal_learning_memory'] },
      'Figure Copy':         { role:'subtest', cluster:'visuoconstruction', construct:'visuoconstructional ability', domains:['visuospatial'] },
      'Line Orientation':    { role:'subtest', cluster:'visuospatial_reasoning', construct:'visual orientation judgement', domains:['visuospatial'] },
      'Picture Naming':      { role:'subtest', cluster:'word_knowledge', construct:'confrontation naming', domains:['language'] },
      'Semantic Fluency':    { role:'subtest', cluster:'verbal_fluency', construct:'semantic fluency', domains:['language','executive'] },
      'Digit Span':          { role:'subtest', cluster:'working_memory', construct:'auditory attention and working memory', domains:['attention_working_memory'] },
      'Coding':              { role:'subtest', cluster:'graphomotor_speed', construct:'graphomotor and visual scanning speed', domains:['processing_speed'] },
      'List Recall':         { role:'subtest', cluster:'verbal_recall', construct:'delayed verbal list recall', domains:['verbal_learning_memory'] },
      'List Recognition':    { role:'subtest', cluster:'verbal_recognition', construct:'verbal recognition memory', domains:['verbal_learning_memory'] },
      'Story Recall':        { role:'subtest', cluster:'verbal_recall', construct:'delayed contextual verbal memory', domains:['verbal_learning_memory'] },
      'Figure Recall':       { role:'subtest', cluster:'visual_recall', construct:'delayed visual reproduction', domains:['visual_learning_memory'] }
    },
    coreByDomain:{
      verbal_learning_memory:['Immediate Memory','Delayed Memory','List Learning','Story Memory','List Recall','Story Recall'],
      visual_learning_memory:['Figure Recall'],
      visuospatial:['Visuospatial/Constructional','Figure Copy','Line Orientation'],
      language:['Language','Picture Naming','Semantic Fluency'],
      attention_working_memory:['Attention','Digit Span'],
      processing_speed:['Attention','Coding'],
      intellectual:['Total Scale']
    },
    indexComposition:{
      'Total Scale':['Immediate Memory','Visuospatial/Constructional','Attention','Language','Delayed Memory'],
      'Immediate Memory':['List Learning','Story Memory'],
      'Visuospatial/Constructional':['Figure Copy','Line Orientation'],
      'Attention':['Digit Span','Coding'],
      'Language':['Picture Naming','Semantic Fluency'],
      'Delayed Memory':['List Recall','List Recognition','Story Recall','Figure Recall']
    }
  },
  {
    id:'dkefs-tmt', name:'D-KEFS Trail Making', longName:'D-KEFS Trail Making Test',
    families:['D-KEFS Trail Making Test'],
    measures:{
      'Visual Scanning':           { role:'condition', cluster:'graphomotor_speed', construct:'visual scanning speed', domains:['processing_speed'] },
      'Number Sequencing':         { role:'condition', cluster:'sequencing_attention', construct:'number sequencing speed', domains:['processing_speed','attention_working_memory'] },
      'Letter Sequencing':         { role:'condition', cluster:'sequencing_attention', construct:'letter sequencing speed', domains:['processing_speed','attention_working_memory'] },
      'Switching':                 { role:'condition', cluster:'switching', construct:'set-shifting and cognitive flexibility', domains:['executive'] },
      'Motor Speed':               { role:'condition', cluster:'graphomotor_speed', construct:'motor speed', domains:['processing_speed'] },
      'Combined Number + Letter':  { role:'condition', cluster:'sequencing_attention', construct:'combined number and letter sequencing', domains:['processing_speed','attention_working_memory'] }
    },
    coreByDomain:{
      processing_speed:['Visual Scanning','Number Sequencing','Letter Sequencing'],
      attention_working_memory:['Number Sequencing','Letter Sequencing'],
      executive:['Switching']
    }
  },
  {
    id:'dkefs-vf', name:'D-KEFS Verbal Fluency', longName:'D-KEFS Verbal Fluency Test',
    families:['D-KEFS Verbal Fluency'],
    measures:{
      'Letter Fluency':       { role:'condition', cluster:'verbal_fluency', construct:'phonemic verbal fluency', domains:['executive','language'] },
      'Category Fluency':     { role:'condition', cluster:'verbal_fluency', construct:'semantic verbal fluency', domains:['executive','language'] },
      'Category Switching':   { role:'condition', cluster:'switching', construct:'category switching fluency', domains:['executive'] },
      'Switching Accuracy':   { role:'condition', cluster:'switching', construct:'category switching accuracy', domains:['executive'] }
    },
    coreByDomain:{
      executive:['Letter Fluency','Category Fluency','Category Switching'],
      language:['Letter Fluency','Category Fluency']
    }
  },
  {
    id:'dkefs-cwi', name:'D-KEFS Colour-Word', longName:'D-KEFS Colour-Word Interference Test',
    families:['D-KEFS Colour-Word Interference'],
    measures:{
      'Colour Naming':         { role:'condition', cluster:'naming_reading_speed', construct:'rapid colour naming', domains:['processing_speed'] },
      'Word Reading':          { role:'condition', cluster:'naming_reading_speed', construct:'rapid word reading', domains:['processing_speed'] },
      'Inhibition':            { role:'condition', cluster:'inhibition', construct:'inhibitory control', domains:['executive'] },
      'Inhibition/Switching':  { role:'condition', cluster:'switching', construct:'inhibition with set-shifting', domains:['executive'] }
    },
    coreByDomain:{
      processing_speed:['Colour Naming','Word Reading'],
      executive:['Inhibition','Inhibition/Switching']
    }
  },
  {
    id:'dkefs-df', name:'D-KEFS Design Fluency', longName:'D-KEFS Design Fluency Test',
    families:['D-KEFS Design Fluency'],
    measures:{
      'Filled Dots':  { role:'condition', cluster:'verbal_fluency', construct:'nonverbal generative fluency (filled)', domains:['executive'] },
      'Empty Dots':   { role:'condition', cluster:'verbal_fluency', construct:'nonverbal generative fluency (empty)', domains:['executive'] },
      'Switching':    { role:'condition', cluster:'switching', construct:'design fluency with set-shifting', domains:['executive'] }
    },
    coreByDomain:{ executive:['Filled Dots','Empty Dots','Switching'] }
  },
  {
    id:'dkefs-sort', name:'D-KEFS Sorting', longName:'D-KEFS Sorting Test',
    families:['D-KEFS Sorting Test'],
    measures:{
      'Free Sorting Confirmed Sorts':{ role:'condition', cluster:'planning_problem_solving', construct:'concept formation (sorting)', domains:['executive'] },
      'Free Sorting Total Score':{ role:'condition', cluster:'planning_problem_solving', construct:'sorting description quality', domains:['executive','language'] },
      'Free Sorting Description Total Score':{ role:'condition', cluster:'planning_problem_solving', construct:'sorting description quality', domains:['executive','language'] },
      'Sort Recognition Score':{ role:'condition', cluster:'planning_problem_solving', construct:'recognition of sorting concepts', domains:['executive','language'] },
      'Sort Recognition Total Description Score':{ role:'condition', cluster:'planning_problem_solving', construct:'recognition of sorting concepts', domains:['executive','language'] },
      'Total Weighted Achievement':{ role:'condition', cluster:'planning_problem_solving', construct:'overall sorting achievement', domains:['executive'] },
      'Initial Abstraction Score':{ role:'condition', cluster:'verbal_reasoning', construct:'initial abstraction', domains:['executive','language'] }
    },
    coreByDomain:{ executive:['Free Sorting Confirmed Sorts','Sort Recognition Score','Total Weighted Achievement'] }
  },
  {
    id:'dkefs-tower', name:'D-KEFS Tower', longName:'D-KEFS Tower Test',
    families:['D-KEFS Tower Test'],
    measures:{
      'Total Achievement Score':{ role:'condition', cluster:'planning_problem_solving', construct:'planning and rule-following', domains:['executive'] }
    },
    coreByDomain:{ executive:['Total Achievement Score'] }
  },
  {
    id:'dkefs-wc', name:'D-KEFS Word Context', longName:'D-KEFS Word Context Test',
    families:['D-KEFS Word Context Test'],
    measures:{
      'Total First Trial Consistently Correct':{ role:'condition', cluster:'verbal_reasoning', construct:'verbal deductive reasoning', domains:['executive','language'] }
    },
    coreByDomain:{ executive:['Total First Trial Consistently Correct'], language:['Total First Trial Consistently Correct'] }
  },
  {
    id:'dkefs-prov', name:'D-KEFS Proverb', longName:'D-KEFS Proverb Test',
    families:['D-KEFS Word Proverb Test'],
    measures:{
      'Total Achievement Score: Free Inquiry':{ role:'condition', cluster:'verbal_reasoning', construct:'abstract verbal reasoning', domains:['language','executive'] },
      'Total Achievement Score':{ role:'condition', cluster:'verbal_reasoning', construct:'abstract verbal reasoning', domains:['language','executive'] }
    },
    coreByDomain:{ language:['Total Achievement Score'], executive:['Total Achievement Score'] }
  }
];
const REPORT_FAMILY_TO_TEST = (() => {
  const map = {};
  REPORT_TEST_CATALOG.forEach(t => t.families.forEach(f => { map[f] = t; }));
  return map;
})();
function reportTestById(id){ return REPORT_TEST_CATALOG.find(t => t.id === id) || null; }
function reportTestForFamily(family){ return REPORT_FAMILY_TO_TEST[family] || null; }
function reportTestMeasureMeta(test, measureName){ return test?.measures?.[measureName] || null; }
// Construct cluster lookup - first match wins. Replaces the previous regex chain.
const REPORT_CLUSTER_RULES = [
  { m:/\b(hads|phq|gad|bdi|bai|dass|pcl|ies|stai|geriatric depression|brief symptom|insomnia)\b|mood|anxiety|depression/i,
    g:'mood_symptoms', c:'mood and symptom reporting', d:['mood'] },
  { m:/digit span forward|auditory attention/i,
    g:'auditory_attention', c:'basic auditory attention', d:['attention_working_memory'] },
  { m:/digit span backward|digit span sequencing|letter[- ]number|arithmetic|spatial addition|symbol span|working memory/i,
    g:'working_memory', c:'working memory and mental manipulation', d:['attention_working_memory'] },
  { m:/(combined )?number sequencing|letter sequencing|sequencing/i,
    g:'sequencing_attention', c:'sequencing and divided attention', d:['attention_working_memory','processing_speed'] },
  { m:/coding|symbol search|cancellation|visual scanning|motor speed|processing speed/i,
    g:'graphomotor_speed', c:'graphomotor and visual scanning speed', d:['processing_speed'] },
  { m:/colou?r naming|word reading/i,
    g:'naming_reading_speed', c:'rapid naming and reading speed', d:['processing_speed'] },
  { m:/inhibition\/switching|inhibition.switching|combined inhibition|number-letter switching|category switching/i,
    g:'switching', c:'set-shifting and cognitive flexibility', d:['executive'] },
  { m:/inhibition/i,
    g:'inhibition', c:'inhibitory control', d:['executive'] },
  { m:/tower|sort recognition|free sorting|confirmed sorts|total achievement|planning|problem.?solving/i,
    g:'planning_problem_solving', c:'planning and executive problem-solving', d:['executive'] },
  { m:/letter fluency|category fluency|verbal fluency|naming.*fluency/i,
    g:'verbal_fluency', c:'verbal fluency and lexical retrieval', d:['executive','language'] },
  { m:/word context|proverb|abstract|similarities/i,
    g:'verbal_reasoning', c:'verbal reasoning and concept formation', d:['language','executive'] },
  { m:/vocabulary|information|comprehension|naming|semantic/i,
    g:'word_knowledge', c:'word knowledge and crystallised verbal ability', d:['language','intellectual'] },
  { m:/cvlt|word list|verbal paired|list (a|learning|trials)/i,
    g:'verbal_learning', c:'verbal list-learning', d:['verbal_learning_memory'] },
  { m:/(short|long)[- ]?delay|delayed|free recall|cued recall|story|logical memory|narrative/i,
    g:'verbal_recall', c:'delayed verbal recall', d:['verbal_learning_memory'] },
  { m:/recognition|discriminability/i,
    g:'verbal_recognition', c:'verbal recognition memory', d:['verbal_learning_memory'] },
  { m:/visual reproduction|design memory|figure(?! copy)/i,
    g:'visual_learning', c:'visual learning and memory', d:['visual_learning_memory'] },
  { m:/figure copy|copy(?!ing)|construction/i,
    g:'visuoconstruction', c:'visuoconstructional ability', d:['visuospatial'] },
  { m:/block design|matrix|visual puzzles|line orientation|perceptual reasoning|visuospatial/i,
    g:'visuospatial_reasoning', c:'visuospatial reasoning', d:['visuospatial','intellectual'] },
  { m:/fsiq|gai|\biq\b|vci|pri|wmi|psi|imi|dmi|index|indices/i,
    g:'global_indices', c:'overall index-level functioning', d:['intellectual'] }
];

const REPORT_DEFAULT_PREFS = {
  reference:'mrs', descriptorSystem:'wechsler', locked:false
};
const REPORT_REFERENCE_MAP = {
  mrs: { name:'Mrs Doe', subject:'She', object:'her', possessive:'Her', firstPerson:false },
  mr:  { name:'Mr Doe',  subject:'He', object:'him', possessive:'His', firstPerson:false },
  you: { name:'You',     subject:'You', object:'you', possessive:'Your', firstPerson:true }
};
function reportReference(){
  return REPORT_REFERENCE_MAP[reportState.prefs.reference] || REPORT_REFERENCE_MAP.mrs;
}

const reportState = {
  sections:[], activeId:null, activeTab:'setup', nextSectionId:1, nextRowId:1,
  options:[], activeFilter:'common', selectedOptionKeys:new Set(), visibleOptionKeys:[],
  prefs:{ ...REPORT_DEFAULT_PREFS },
  outputHtml:'', stale:false, hydrated:false, suppressSave:false
};

/* ---------- Persistence ---------- */
function reportSaveState(){
  if (reportState.suppressSave) return;
  try {
    localStorage.setItem(REPORT_STORAGE_KEY, JSON.stringify({
      sections: reportState.sections,
      activeId: reportState.activeId,
      activeTab: reportState.activeTab,
      nextSectionId: reportState.nextSectionId,
      nextRowId: reportState.nextRowId,
      prefs: reportState.prefs,
      outputHtml: reportState.outputHtml
    }));
  } catch(e){ /* ignore quota / privacy */ }
}
function reportLoadState(){
  try {
    const raw = localStorage.getItem(REPORT_STORAGE_KEY);
    if (!raw) return;
    const payload = JSON.parse(raw);
    if (!payload || !Array.isArray(payload.sections)) return;
    reportState.sections = payload.sections.map(s => ({ ...s, rows: Array.isArray(s.rows) ? s.rows : [] }));
    reportState.activeId = payload.activeId || (reportState.sections[0]?.id ?? null);
    reportState.activeTab = ['setup','build','scores'].includes(payload.activeTab) ? payload.activeTab : 'setup';
    reportState.nextSectionId = payload.nextSectionId || (reportState.sections.length + 1);
    reportState.nextRowId = payload.nextRowId || 1;
    reportState.prefs = { ...REPORT_DEFAULT_PREFS, ...(payload.prefs || {}) };
    // Migrate legacy `pronouns` field → `reference`
    if (!reportState.prefs.reference && payload.prefs?.pronouns){
      const legacyMap = { she:'mrs', he:'mr', they:'mrs' };
      reportState.prefs.reference = legacyMap[payload.prefs.pronouns] || 'mrs';
    }
    delete reportState.prefs.pronouns;
    reportState.outputHtml = payload.outputHtml || '';
    reportState.hydrated = true;
  } catch(e){ /* ignore corrupt */ }
}

/* ---------- Helpers ---------- */
function reportCleanText(value){
  return String(value || '')
    .replace(/Â·/g, '·')
    .replace(/â€“/g, '–')
    .replace(/â€”/g, '-')
    .replace(/â€™/g, "'")
    .replace(/\s+/g, ' ')
    .trim();
}
function reportDomain(id){ return REPORT_DOMAINS.find(d => d.id === id) || REPORT_DOMAINS[REPORT_DOMAINS.length - 1]; }
function reportScoreType(id){ return REPORT_SCORE_TYPES.find(t => t.id === id) || REPORT_SCORE_TYPES[0]; }
function reportPossessiveName(name){
  const n = reportCleanText(name) || 'the patient';
  return /s$/i.test(n) ? `${n}'` : `${n}'s`;
}
function reportSubject(form, useName){
  const r = reportReference();
  // First-person ignores useName entirely (always "you/your")
  if (r.firstPerson || !useName){
    return form === 'possessive' ? r.possessive : r.subject;
  }
  return form === 'possessive' ? reportPossessiveName(r.name) : r.name;
}
function reportLowerFirst(value){
  const text = String(value || '');
  return text ? text.charAt(0).toLowerCase() + text.slice(1) : text;
}
function reportDisplayName(){
  return reportReference().name;
}
function reportOrdinal(value){
  const raw = Number(value);
  if (raw > 0 && raw < 1) return '&lt;1<sup>st</sup>';
  if (raw > 99 && raw < 100) return '&gt;99<sup>th</sup>';
  const n = Math.round(raw);
  if (!Number.isFinite(n)) return '';
  const mod100 = n % 100;
  let suffix = 'th';
  if (mod100 < 11 || mod100 > 13){
    if (n % 10 === 1) suffix = 'st';
    else if (n % 10 === 2) suffix = 'nd';
    else if (n % 10 === 3) suffix = 'rd';
  }
  return `${n}<sup>${suffix}</sup>`;
}
function reportFormatList(items){
  const list = [...new Set(items.map(reportCleanText).filter(Boolean))];
  if (list.length === 0) return '';
  if (list.length === 1) return list[0];
  if (list.length === 2) return `${list[0]} and ${list[1]}`;
  return `${list.slice(0, -1).join(', ')}, and ${list[list.length - 1]}`;
}
function reportDescriptorFor(ss){
  return reportState.prefs.descriptorSystem === 'aan' ? aanDesc(ss) : wechslerDesc(ss);
}
/* ---------- Cluster lookup (replaces regex chain) ---------- */
function reportLookupCluster(family, subtest, domainHint){
  const text = `${family || ''} ${subtest || ''}`;
  for (const rule of REPORT_CLUSTER_RULES){
    if (rule.m.test(text)) return rule;
  }
  return { g:'other', c: reportDomain(domainHint || 'other').ability, d:['other'] };
}
function reportInferDomains(family, subtest){ return reportLookupCluster(family, subtest).d; }
function reportInferConstruct(family, subtest, domainHint){ return reportLookupCluster(family, subtest, domainHint).c; }
function reportInferConstructGroup(family, subtest, domainHint){ return reportLookupCluster(family, subtest, domainHint).g; }
function reportIsCommonOption(family, subtest, domainId, groupId){
  const text = `${family} ${subtest}`.toLowerCase();
  if (domainId === 'executive') return /inhibition|switch|category switching|tower|sorting|total achievement|confirmed sorts/.test(text);
  if (domainId === 'processing_speed') return /coding|symbol search|visual scanning|motor speed|colou?r naming|word reading/.test(text);
  if (domainId === 'attention_working_memory') return /digit span|letter[- ]number|arithmetic|working memory|symbol span/.test(text);
  if (domainId === 'verbal_learning_memory') return /cvlt|list a|short delay|long delay|recognition|logical memory|story/.test(text);
  if (domainId === 'visual_learning_memory') return /visual reproduction|figure|design memory/.test(text);
  if (domainId === 'language') return /vocabulary|similarities|naming|letter fluency|category fluency|word context/.test(text);
  if (domainId === 'visuospatial') return /block design|matrix|visual puzzles|figure copy|line orientation/.test(text);
  if (domainId === 'intellectual') return /fsiq|gai|vci|pri|wmi|psi|index|indices|vocabulary|matrix|block design|coding/.test(text);
  return groupId !== 'other';
}
/* ---------- Score conversion ---------- */
function reportScoreToStandard(row){
  const scoreType = row.scoreType || 'standard';
  const score = parseFloat(row.score);
  if (scoreType === 'percentile'){
    const pct = parseFloat(row.percentile || row.score);
    const z = toZ(pct, 'percentile');
    return z == null ? null : fromZ(z, 'standard');
  }
  if (!Number.isFinite(score)){
    const pct = parseFloat(row.percentile);
    const z = toZ(pct, 'percentile');
    return z == null ? null : fromZ(z, 'standard');
  }
  const z = toZ(score, scoreType);
  return z == null ? null : fromZ(z, 'standard');
}
function reportUpdateComputed(row){
  const ss = reportScoreToStandard(row);
  if (ss == null || !Number.isFinite(ss)){
    row.descriptor = '';
    row.standardScore = null;
    return row;
  }
  const z = toZ(ss, 'standard');
  if (row.scoreType !== 'percentile') row.percentile = fmtPct(normCDF(z) * 100);
  row.descriptor = reportDescriptorFor(ss);
  row.standardScore = ss;
  return row;
}

/* ---------- Parenthetical (long form, supports CI) ---------- */
function reportScoreParenthetical(row){
  const parts = [];
  const type = reportScoreType(row.scoreType);
  if (row.scoreType === 'percentile'){
    const pct = String(row.percentile || row.score || '').trim();
    if (pct) parts.push(`${reportOrdinal(pct)} percentile`);
  } else {
    if (String(row.score || '').trim() !== '') parts.push(`${type.sentence} = ${escapeHtml(row.score)}`);
    if (String(row.ci || '').trim() !== '') parts.push(`95% CI [${escapeHtml(row.ci)}]`);
    if (String(row.percentile || '').trim() !== '') parts.push(`${reportOrdinal(row.percentile)} percentile`);
  }
  return parts.length ? ` (${parts.join(', ')})` : '';
}
function reportMeasureLabel(row){
  const test = reportCleanText(row.displayName || row.test);
  const subtest = reportCleanText(row.subtest);
  if (test && subtest && !test.toLowerCase().includes(subtest.toLowerCase())) return `${test}: ${subtest}`;
  return test || subtest || 'the selected measure';
}

/* ---------- Descriptor rank (for spread analysis) ---------- */
const REPORT_RANKS = {
  'extremely low':0, 'exceptionally low':0,
  'borderline':1,    'below average':1,
  'low average':2,
  'average':3,
  'high average':4,
  'superior':5,      'above average':5,
  'very superior':6, 'exceptionally high':6
};
function reportDescriptorRank(desc){
  const d = reportCleanText(desc).toLowerCase();
  return Object.prototype.hasOwnProperty.call(REPORT_RANKS, d) ? REPORT_RANKS[d] : null;
}
function reportRankToDescriptor(rank){
  const sys = reportState.prefs.descriptorSystem;
  const wechsler = ['extremely low','borderline','low average','average','high average','superior','very superior'];
  const aan = ['exceptionally low','below average','low average','average','high average','above average','exceptionally high'];
  return (sys === 'aan' ? aan : wechsler)[rank] || '';
}

/* ---------- Grouping ---------- */
function reportRowGroupKey(row){
  // Cluster by descriptor + score type only - measures with the same descriptor
  // (e.g. four CVLT trials all in the Average range) cluster into one sentence,
  // even if their exact scores differ. The parenthetical shows a range when needed.
  return [
    reportCleanText(row.descriptor).toLowerCase(),
    reportCleanText(row.scoreType)
  ].join('|');
}
/* ---------- Filtering: only rows with usable scores enter the prose ---------- */
function reportFilterScoredRows(rows){
  return rows.filter(row => reportCleanText(row.test) && reportScoreToStandard(row) != null);
}

/* ---------- Test-grouped list helper ---------- */
// Format a list of rows so that 2+ measures from the same test read as
// "A and B from the [Test]". Single measures stay bare unless multiple
// test groups appear in the same list (in which case all groups label).
// Compact a sequential numbered series (e.g. ["Trial 1","Trial 2","Trial 3"]
// → ["Trials 1–3"]). Returns the original list if not a clean run.
function reportCompactSequentialNames(names){
  if (names.length < 3) return names;
  const matches = names.map(n => String(n).match(/^(.*?)\s+(\d+)$/));
  if (matches.some(m => !m)) return names;
  const prefix = matches[0][1].trim();
  if (!prefix) return names;
  if (matches.some(m => m[1].trim() !== prefix)) return names;
  const numbers = matches.map(m => parseInt(m[2], 10));
  const sorted = [...numbers].sort((a, b) => a - b);
  for (let i = 1; i < sorted.length; i++){
    if (sorted[i] !== sorted[i-1] + 1) return names;
  }
  // De-duped sort already ensured consecutive
  const plural = /s$/i.test(prefix) ? prefix : `${prefix}s`;
  return [`${plural} ${sorted[0]}–${sorted[sorted.length - 1]}`];
}

function reportFormatTestGroupedList(rows){
  const order = [];
  const byTest = new Map();
  rows.forEach(r => {
    const key = r.testId || `__row_${r.id}__`;
    if (!byTest.has(key)){
      byTest.set(key, { testName: r.testName || r.testLongName || '', rows:[] });
      order.push(key);
    }
    byTest.get(key).rows.push(r);
  });
  const onlyOneGroup = order.length === 1;
  const parts = order.map(key => {
    const g = byTest.get(key);
    const rawLabels = g.rows.map(r => reportMeasureLabel(r));
    const labels = reportCompactSequentialNames(rawLabels);
    if (labels.length === 1 && g.rows.length === 1){
      if (onlyOneGroup || !g.testName) return labels[0];
      return `${labels[0]} from the ${g.testName}`;
    }
    if (g.testName) return `${reportFormatList(labels)} from the ${g.testName}`;
    return reportFormatList(labels);
  });
  return reportFormatList(parts);
}

// Parenthetical for a measure-group cluster.
// - 1 row: existing per-row parenthetical (with CI if present)
// - 2+ rows, all scores identical: same as per-row (CI dropped - per-row only)
// - 2+ rows, scores vary: range syntax, e.g. (scaled scores 9–11, 37th–63rd percentile)
function reportClusterParenthetical(rows){
  if (rows.length === 1) return reportScoreParenthetical(rows[0]);
  const type = reportScoreType(rows[0].scoreType);
  const isPercentileType = rows[0].scoreType === 'percentile';
  const scores = rows.map(r => parseFloat(r.score)).filter(Number.isFinite);
  const pcts   = rows.map(r => parseFloat(r.percentile)).filter(Number.isFinite);
  const fmtS = n => Number.isInteger(n) ? String(n) : (Math.round(n * 10) / 10).toString();

  if (isPercentileType){
    if (!pcts.length) return '';
    const min = Math.min(...pcts), max = Math.max(...pcts);
    if (Math.round(min) === Math.round(max)) return ` (${reportOrdinal(min)} percentile)`;
    return ` (${reportOrdinal(min)} to ${reportOrdinal(max)} percentile)`;
  }

  if (!scores.length){
    if (!pcts.length) return '';
    const min = Math.min(...pcts), max = Math.max(...pcts);
    if (Math.round(min) === Math.round(max)) return ` (${reportOrdinal(min)} percentile)`;
    return ` (${reportOrdinal(min)} to ${reportOrdinal(max)} percentile)`;
  }

  const minS = Math.min(...scores), maxS = Math.max(...scores);
  const sameScore = minS === maxS;
  const parts = [];
  parts.push(sameScore
    ? `${type.sentence} = ${fmtS(minS)}`
    : `${type.sentence}s ${fmtS(minS)}–${fmtS(maxS)}`);

  if (pcts.length){
    const minP = Math.min(...pcts), maxP = Math.max(...pcts);
    if (Math.round(minP) === Math.round(maxP)){
      parts.push(`${reportOrdinal(minP)} percentile`);
    } else {
      parts.push(`${reportOrdinal(Math.round(minP))} to ${reportOrdinal(Math.round(maxP))} percentile`);
    }
  }
  return ` (${parts.join(', ')})`;
}

/* ---------- Connector ---------- */
function reportComputeConnector(descriptor, rank, previous, counts, isFinal){
  if (isFinal) return 'Finally,';
  if (!previous) return '';
  const cur = reportCleanText(descriptor || '').toLowerCase();
  const prev = reportCleanText(previous.descriptor || '').toLowerCase();
  if (cur && prev && cur === prev){
    const opts = ['Similarly,', 'Also,'];
    const text = opts[counts.same % opts.length];
    counts.same += 1;
    return text;
  }
  if (rank != null && previous.rank != null && Math.abs(rank - previous.rank) >= 2){
    const opts = ['However,', 'Alternatively,', 'Conversely,'];
    const text = opts[counts.contrast % opts.length];
    counts.contrast += 1;
    return text;
  }
  return '';
}
function reportSentencePrefix(connector, useName){
  const possessive = reportSubject('possessive', useName);
  const ref = reportReference();
  if (!connector) return possessive;
  const after = (useName || ref.firstPerson) ? possessive : reportLowerFirst(possessive);
  return `${connector} ${after}`;
}

/* ---------- Sentence renderers ---------- */
// Detect if a row-list compacts via trial naming ("Trial 1, 2, 3" → "Trials 1–3").
// Only when ALL rows compact into a single label do we use the group-level range
// parenthetical; otherwise each measure carries its own inline parenthetical.
function reportListingShape(rows){
  if (rows.length <= 1) return { compacted: false, list: reportFormatTestGroupedList(rows) };
  const labels = rows.map(r => reportMeasureLabel(r));
  const compacted = reportCompactSequentialNames(labels);
  return {
    compacted: compacted.length === 1 && labels.length > 1,
    list: reportFormatTestGroupedList(rows)
  };
}

function reportRenderMeasureGroupSentence(group, previous, counts, useName, isFinal){
  const first = group.rows[0];
  const descriptor = (first.descriptor || '').toLowerCase();
  const desc = escapeHtml(descriptor || 'unclassified');
  const connector = reportComputeConnector(group.descriptor, group.rank, previous, counts, isFinal);
  const ref = reportReference();
  const size = group.rows.length;

  if (size === 1){
    const tests = escapeHtml(reportFormatTestGroupedList(group.rows));
    const paren = reportScoreParenthetical(first);
    const prefix = reportSentencePrefix(connector, useName);
    return `${prefix} performance on ${tests} fell within the <em>${desc}</em> range${paren}.`;
  }

  const shape = reportListingShape(group.rows);
  if (shape.compacted){
    const tests = escapeHtml(shape.list);
    const paren = reportClusterParenthetical(group.rows);
    if (size >= 3){
      if (useName || ref.firstPerson){
        const prefix = reportSentencePrefix(connector, useName);
        return `${prefix} performance on ${tests} was consistently within the <em>${desc}</em> range${paren}.`;
      }
      const lead = connector ? `${connector} ` : '';
      return `${lead}${tests} all fell within the <em>${desc}</em> range${paren}.`;
    }
    const prefix = reportSentencePrefix(connector, useName);
    return `${prefix} performance on ${tests} fell within the <em>${desc}</em> range${paren}.`;
  }

  // Multi-row, not compacted - inline parens per measure
  const inlineItems = group.rows.map(r =>
    `${escapeHtml(reportMeasureLabel(r))}${reportScoreParenthetical(r)}`
  );
  const tests = reportFormatList(inlineItems);
  if (size >= 3){
    if (useName || ref.firstPerson){
      const prefix = reportSentencePrefix(connector, useName);
      return `${prefix} performance on ${tests} was consistently within the <em>${desc}</em> range.`;
    }
    const lead = connector ? `${connector} ` : '';
    return `${lead}${tests} all fell within the <em>${desc}</em> range.`;
  }
  const prefix = reportSentencePrefix(connector, useName);
  return `${prefix} performance on ${tests} fell within the <em>${desc}</em> range.`;
}
/* ---------- Plan + render whole section paragraph ---------- */
function reportRowGroupKey(row){
  // Cluster by descriptor + score type only - measures with the same descriptor
  // (e.g. four CVLT trials all in the Average range) cluster into one sentence,
  // even if their exact scores differ. The parenthetical shows a range when needed.
  return [
    reportCleanText(row.descriptor).toLowerCase(),
    reportCleanText(row.scoreType)
  ].join('|');
}
// Plan sentences: sort indices first, then group adjacent rows by descriptor + score type.
function reportPlanSentences(rows){
  const filtered = reportFilterScoredRows(rows);
  if (!filtered.length) return [];
  filtered.forEach(reportUpdateComputed);
  const sorted = [
    ...filtered.filter(r => r.role === 'index'),
    ...filtered.filter(r => r.role !== 'index')
  ];
  const groups = [];
  sorted.forEach(row => {
    const key = reportRowGroupKey(row);
    const previous = groups[groups.length - 1];
    if (previous && previous.key === key){
      previous.rows.push(row);
    } else {
      groups.push({ key, rows:[row], descriptor:row.descriptor || '', rank:reportDescriptorRank(row.descriptor) });
    }
  });
  return groups.map(g => ({ type:'group', group:g }));
}

function reportRenderSentences(sentences, opts = {}){
  if (!sentences.length) return '';
  const counts = { same:0, contrast:0 };
  const total = sentences.length;
  const startWithName = !!opts.startWithName;
  const ref = reportReference();
  return sentences.map((s, i) => {
    const useName = ref.firstPerson
      ? false
      : ((i === 0 && startWithName) || (i > 0 && i % 3 === 0));
    const isFinal = (i === total - 1) && total >= 3;
    const previous = i > 0 ? { descriptor: sentences[i - 1].group.descriptor, rank: sentences[i - 1].group.rank } : null;
    return reportRenderMeasureGroupSentence(s.group, previous, counts, useName, isFinal);
  }).join(' ');
}

function reportTextForRows(rows, opts = {}){
  return reportRenderSentences(reportPlanSentences(rows), opts);
}
function reportTextForRow(row){ return reportTextForRows([row]); }

/* ---------- Report-level Summary (Overall + top + bottom scores) ---------- */
function reportBuildSummary(){
  const allRows = reportState.sections.flatMap(s => s.rows.map(r => ({ ...r, section: s })))
    .filter(r => reportCleanText(r.test) && Number.isFinite(reportScoreToStandard(r)));
  if (allRows.length < 3) return ''; // not enough data for a meaningful summary

  // Compute standard-score equivalents and rank
  const enriched = allRows.map(r => ({
    row: r,
    ss: reportScoreToStandard(r),
    descriptor: (function(){ reportUpdateComputed(r); return r.descriptor || ''; })()
  })).filter(e => Number.isFinite(e.ss));
  if (enriched.length < 3) return '';

  const sortedDesc = [...enriched].sort((a, b) => b.ss - a.ss);
  const top = sortedDesc.slice(0, 3);
  const bottom = sortedDesc.slice(-3).reverse();

  const subjectName = reportSubject('possessive', true); // first mention - full name
  const ref = reportReference();

  const renderItem = e => `${escapeHtml(reportMeasureLabel(e.row))}${reportScoreParenthetical(e.row)}`;
  const topList = reportFormatList(top.map(renderItem));
  const bottomList = reportFormatList(bottom.map(renderItem));

  // Overall sentence - pick a phrasing based on overall central tendency (mean SS).
  // reportSubject already returns the correctly-cased token: lowercase 'your' for
  // first-person, proper-cased "Mrs Doe's" for named reports.
  const meanSs = enriched.reduce((a, e) => a + e.ss, 0) / enriched.length;
  const overallDescriptor = reportDescriptorFor(meanSs).toLowerCase();
  const overallLine = `Overall, ${reportSubject('possessive', true)} performance across the assessed domains was broadly within the <em>${escapeHtml(overallDescriptor)}</em> range.`;

  // Top / bottom sentences use pronouns rather than restating the name:
  //   "Her highest performance was observed on X."
  //   "However, she performed poorest on Y."
  // For first-person: "Your highest..." / "However, you performed poorest on Y."
  const subjPronoun = reportSubject('subject', false);
  const possPronoun = reportSubject('possessive', false);
  const possCap = possPronoun.charAt(0).toUpperCase() + possPronoun.slice(1);
  const topNoun = top.length === 1 ? 'highest performance' : 'highest performances';
  const topVerb = top.length === 1 ? 'was' : 'were';
  const topLine = `${possCap} ${topNoun} ${topVerb} observed on ${topList}.`;
  const bottomLine = `However, ${subjPronoun} performed poorest on ${bottomList}.`;

  return `<h3>Summary</h3><p>${overallLine} ${topLine} ${bottomLine}</p>`;
}

/* ---------- Section rendering - single paragraph, indices first then subtests ---------- */
function reportRenderSectionHtml(section, rows, scored){
  const out = [];
  out.push(`<h3>${escapeHtml(reportSectionHeading(section))}</h3>`);
  if (!scored.length) return out.join('');
  const opener = reportSectionOpener(section, rows);
  const detail = reportTextForRows(scored, { startWithName: false });
  out.push(`<p>${opener}${opener && detail ? ' ' : ''}${detail}</p>`);
  return out.join('');
}

/* ---------- Section opener - "[Subject's] [domain] was assessed using [list]." ---------- */
function reportSectionOpener(section, rows){
  const ability = reportDomain(section.domain).ability;
  const subj = reportSubject('possessive', true); // first mention - name
  const valid = reportFilterScoredRows(rows);
  if (!valid.length) return '';
  const list = escapeHtml(reportFormatTestGroupedList(valid));
  if (section.domain === 'other'){
    return `${subj} performance was assessed using ${list}.`;
  }
  return `${subj} ${escapeHtml(ability)} was assessed using ${list}.`;
}

function reportFamilyOrder(family){
  const clean = reportCleanText(family).toLowerCase();
  const idx = REPORT_FAMILY_ORDER.findIndex(item => clean.includes(item.toLowerCase()));
  return idx === -1 ? 999 : idx;
}
function reportGroupOrder(domainId, groupId){
  const order = REPORT_DOMAIN_GROUP_ORDER[domainId] || REPORT_DOMAIN_GROUP_ORDER.other;
  const idx = order.indexOf(groupId);
  return idx === -1 ? 999 : idx;
}
function reportOptionOrder(option, domainId){
  const text = `${option.family} ${option.displayName}`.toLowerCase();
  let subOrder = 50;
  [
    /visual scanning/, /number sequencing/, /letter sequencing/, /coding/, /symbol search/,
    /colour naming|color naming/, /word reading/, /inhibition$/, /inhibition\/switching|inhibition switching/,
    /letter fluency/, /category fluency/, /category switching/, /switching accuracy/,
    /tower|total achievement/, /sorting|confirmed sorts|sort recognition/,
    /vocabulary/, /similarities/, /matrix/, /block design/
  ].some((pattern, index) => {
    if (pattern.test(text)){ subOrder = index; return true; }
    return false;
  });
  return (reportGroupOrder(domainId, option.groupId) * 1000) + (option.common ? 0 : 100) + (reportFamilyOrder(option.family) * 10) + subOrder;
}
function reportBuildMeasureOptions(){
  const db = getMergedDB();
  const seen = new Set();
  const options = [];
  Object.entries(db).forEach(([family, tests]) => {
    const cleanFamily = reportCleanText(family);
    const baseFamily = reportCleanText(stripAgeRange(cleanFamily));
    const catalogTest = reportTestForFamily(baseFamily);
    Object.keys(tests || {}).forEach(subtest => {
      const cleanSub = reportCleanText(subtest);
      const key = `${baseFamily}::${cleanSub}`;
      if (seen.has(key)) return;
      seen.add(key);
      // Catalog metadata wins; fall back to cluster-rule inference for measures not in the catalog.
      const meta = reportTestMeasureMeta(catalogTest, cleanSub);
      let role, cluster, construct, domains;
      if (meta){
        role = meta.role;
        cluster = meta.cluster;
        construct = meta.construct;
        domains = meta.domains && meta.domains.length ? meta.domains : ['other'];
      } else {
        const fallback = reportLookupCluster(baseFamily, cleanSub);
        role = 'measure';
        cluster = fallback.g;
        construct = fallback.c;
        domains = fallback.d.length ? fallback.d : ['other'];
      }
      const testId = catalogTest?.id || `family:${baseFamily}`;
      const testName = catalogTest?.name || baseFamily;
      const testLongName = catalogTest?.longName || baseFamily;
      options.push({
        key,
        family:baseFamily,
        subtest:cleanSub,
        label:`${baseFamily}: ${cleanSub}`,
        displayName:cleanSub,
        testId, testName, testLongName,
        role,
        domains,
        domain:domains[0],
        groupId:cluster,
        groupLabel:REPORT_CONSTRUCT_LABELS[cluster] || REPORT_CONSTRUCT_LABELS.other,
        construct,
        common:domains.some(domainId => reportIsCommonOption(baseFamily, cleanSub, domainId, cluster)),
        aliases:[baseFamily, cleanSub, construct, testName, REPORT_CONSTRUCT_LABELS[cluster] || ''].join(' '),
        scoreType:inferScoreType(baseFamily)
      });
    });
  });
  options.sort((a,b) => a.label.localeCompare(b.label));
  reportState.options = options;
}
function reportTestCoreKeysForDomain(testId, domainId){
  const test = reportTestById(testId);
  if (!test || !test.coreByDomain) return new Set();
  const names = test.coreByDomain[domainId] || [];
  // Build keys by matching against options (a measure may live in multiple families/age bands;
  // resolve to the option keys present in the picker).
  const keys = new Set();
  reportState.options.forEach(opt => {
    if (opt.testId !== testId) return;
    if (names.includes(opt.subtest)) keys.add(opt.key);
  });
  return keys;
}
function reportActiveSection(){ return reportState.sections.find(s => s.id === reportState.activeId) || null; }
function renderReportDomainChecklist(){
  const list = document.getElementById('rw-domain-checklist');
  if (!list) return;
  const used = new Set(reportState.sections.map(s => s.domain));
  const available = REPORT_DOMAINS.filter(d => !used.has(d.id));
  if (!available.length){
    list.innerHTML = '<div class="rw-domain-checklist-empty">All domains already added.</div>';
    return;
  }
  list.innerHTML = available.map(d => `
    <label class="rw-domain-check">
      <input type="checkbox" data-rw-domain="${d.id}">
      <span>${escapeHtml(d.label)}</span>
    </label>
  `).join('');
}
function renderReportSectionList(){
  const list = document.getElementById('rw-section-list');
  if (!list) return;
  if (reportState.sections.length === 0){
    list.innerHTML = '<div class="report-empty-state">No report sections yet. Add a heading/subheading above.</div>';
    return;
  }
  const last = reportState.sections.length - 1;
  list.innerHTML = reportState.sections.map((section, i) => {
    const num = String(i + 1).padStart(2, '0');
    const count = section.rows.length;
    return `
    <div class="report-section-item ${section.id === reportState.activeId ? 'active' : ''}" data-rw-section="${section.id}" draggable="true">
      <span class="report-section-grip" aria-hidden="true">≡</span>
      <div class="report-section-content">
        <span class="report-section-num">${num}</span>
        <span class="report-section-title">${escapeHtml(reportSectionHeading(section))}</span>
        <span class="report-section-meta">${count} measure${count === 1 ? '' : 's'}</span>
      </div>
      <div class="report-section-actions">
        <button type="button" data-rw-move="up" data-rw-section="${section.id}" title="Move up" ${i === 0 ? 'disabled' : ''}>▲</button>
        <button type="button" data-rw-move="down" data-rw-section="${section.id}" title="Move down" ${i === last ? 'disabled' : ''}>▼</button>
        <button class="report-section-remove" type="button" data-rw-remove-section="${section.id}" title="Remove section">×</button>
      </div>
    </div>
    `;
  }).join('');
}
function reportOptionAlreadyAdded(section, key){
  return !!section && !!key && section.rows.some(row => row.sourceKey === key);
}
function reportResetPickerState(clearSearch = false){
  reportState.activeFilter = 'common';
  reportState.selectedOptionKeys.clear();
  reportState.visibleOptionKeys = [];
  if (clearSearch){
    const search = document.getElementById('rw-test-search');
    if (search) search.value = '';
  }
}
function reportUpdateFloatingAddBar(){
  const bar = document.getElementById('rw-picker-floating-add');
  const count = document.getElementById('rw-picker-floating-count');
  const addBtn = document.getElementById('rw-add-mapped');
  if (!bar) return;
  const n = reportState.selectedOptionKeys.size;
  if (n === 0){ bar.hidden = true; return; }
  bar.hidden = false;
  if (count) count.textContent = `${n} selected`;
  if (addBtn) addBtn.textContent = `Add ${n} →`;
}
function renderReportTestOptions(){
  const browser = document.getElementById('rw-test-list');
  const chips = document.getElementById('rw-filter-chips');
  if (!browser || !chips) return;
  const section = reportActiveSection();
  const query = reportCleanText(document.getElementById('rw-test-search')?.value).toLowerCase();
  if (!section){
    chips.innerHTML = '';
    reportState.selectedOptionKeys.clear();
    reportState.visibleOptionKeys = [];
    reportUpdateFloatingAddBar();
    browser.innerHTML = '<div class="report-empty-state">Add a section first, then mapped tests will appear here.</div>';
    return;
  }
  const domainId = section.domain;
  [...reportState.selectedOptionKeys].forEach(key => {
    if (reportOptionAlreadyAdded(section, key)) reportState.selectedOptionKeys.delete(key);
  });

  // All options that map to this domain
  const allForDomain = reportState.options.filter(option =>
    option.domains.includes(domainId) || domainId === 'other'
  );

  // Pre-compute which options are "core" for the current domain (per-test curated set)
  const coreKeySet = new Set();
  REPORT_TEST_CATALOG.forEach(t => reportTestCoreKeysForDomain(t.id, domainId).forEach(k => coreKeySet.add(k)));

  // Tests present in this domain (preserve catalog order, then non-catalog families alphabetically)
  const testsPresent = [];
  const seenTests = new Set();
  REPORT_TEST_CATALOG.forEach(t => {
    if (allForDomain.some(o => o.testId === t.id)){
      testsPresent.push({ id:t.id, name:t.name, longName:t.longName, isCatalog:true });
      seenTests.add(t.id);
    }
  });
  const nonCatalogIds = [...new Set(allForDomain.filter(o => !seenTests.has(o.testId)).map(o => o.testId))]
    .sort((a,b) => (allForDomain.find(o => o.testId === a)?.testName || '').localeCompare(allForDomain.find(o => o.testId === b)?.testName || ''));
  nonCatalogIds.forEach(id => {
    const sample = allForDomain.find(o => o.testId === id);
    if (sample) testsPresent.push({ id, name:sample.testName, longName:sample.testLongName, isCatalog:false });
  });

  // Resolve active filter (Core / All / per-test). Default to Core if any core measures exist for this domain.
  const hasCore = coreKeySet.size > 0;
  const validFilters = new Set(['all', ...(hasCore ? ['common'] : []), ...testsPresent.map(t => `test:${t.id}`)]);
  if (!validFilters.has(reportState.activeFilter)) reportState.activeFilter = hasCore ? 'common' : 'all';

  // Build chip row
  const chipDefs = [
    ...(hasCore ? [{ id:'common', label:'Core' }] : []),
    { id:'all', label:'All' },
    ...testsPresent.map(t => ({ id:`test:${t.id}`, label:t.name, title:t.longName }))
  ];
  chips.innerHTML = chipDefs.map(c =>
    `<button class="report-filter-chip ${reportState.activeFilter === c.id ? 'active' : ''}" type="button" data-rw-filter="${escapeAttr(c.id)}"${c.title ? ` title="${escapeAttr(c.title)}"` : ''}>${escapeHtml(c.label)}</button>`
  ).join('');

  // Filter according to active filter + search
  const filterId = reportState.activeFilter;
  const filtered = allForDomain.filter(option => {
    const searchText = `${option.displayName} ${option.family} ${option.construct} ${option.testName} ${option.aliases}`.toLowerCase();
    if (query && !searchText.includes(query)) return false;
    if (filterId === 'all') return true;
    if (filterId === 'common') return coreKeySet.has(option.key);
    if (filterId.startsWith('test:')) return option.testId === filterId.slice(5);
    return true;
  });

  // Group by test then sort by role within each test
  const byTest = new Map();
  filtered.forEach(option => {
    if (!byTest.has(option.testId)){
      byTest.set(option.testId, { name:option.testName, longName:option.testLongName, options:[] });
    }
    byTest.get(option.testId).options.push(option);
  });
  // Test display order matches testsPresent
  const orderedTestIds = testsPresent.map(t => t.id).filter(id => byTest.has(id));

  // Visible cap
  const visible = filtered.slice(0, REPORT_PICKER_LIMIT);
  const visibleSet = new Set(visible.map(o => o.key));
  reportState.visibleOptionKeys = visible.map(o => o.key);

  reportUpdateFloatingAddBar();

  if (!filtered.length){
    browser.innerHTML = '<div class="report-empty-state">No mapped measures match this filter. Try a different test or search term, or add a custom row below.</div>';
    return;
  }

  const limitNote = filtered.length > visible.length
    ? `<div class="report-picker-limit-note">Showing the first ${visible.length} of ${filtered.length} measures. Narrow with a test chip or search.</div>`
    : '';

  // Render - one section per test, role bands inside
  const sectionsHtml = orderedTestIds.map(testId => {
    const group = byTest.get(testId);
    const optionsInTest = group.options.filter(o => visibleSet.has(o.key));
    if (!optionsInTest.length) return '';
    // Pre-compute Add core set / Add all set for this test (within current domain)
    const coreKeys = reportTestCoreKeysForDomain(testId, domainId);
    const allInDomainForTest = allForDomain.filter(o => o.testId === testId);
    const coreAvailable = [...coreKeys].filter(k => !reportOptionAlreadyAdded(section, k)).length;
    const allAvailable = allInDomainForTest.filter(o => !reportOptionAlreadyAdded(section, o.key)).length;
    // Sort by role then display name
    optionsInTest.sort((a,b) => {
      const ra = REPORT_ROLE_ORDER.indexOf(a.role); const rb = REPORT_ROLE_ORDER.indexOf(b.role);
      if (ra !== rb) return (ra === -1 ? 999 : ra) - (rb === -1 ? 999 : rb);
      return a.displayName.localeCompare(b.displayName);
    });
    // Group within test by role
    const byRole = new Map();
    optionsInTest.forEach(o => {
      if (!byRole.has(o.role)) byRole.set(o.role, []);
      byRole.get(o.role).push(o);
    });
    const roleHtml = [...byRole.entries()].map(([role, opts]) => `
      <div class="report-role-band">
        <div class="report-role-band-head">${escapeHtml(REPORT_ROLE_LABELS[role] || 'Measures')}</div>
        ${opts.map(option => {
          const isAdded = reportOptionAlreadyAdded(section, option.key);
          // Checkbox now reflects "in this section?" - tick to add, untick to remove
          const checked = isAdded || reportState.selectedOptionKeys.has(option.key);
          const badge = isAdded ? 'Added' : coreKeySet.has(option.key) ? 'Core' : '';
          return `
            <label class="report-test-option ${isAdded ? 'is-added' : ''} ${checked && !isAdded ? 'is-selected' : ''}">
              <input type="checkbox" data-rw-option="${escapeAttr(option.key)}" ${checked ? 'checked' : ''}>
              <span>
                <span class="report-test-name">${escapeHtml(option.displayName)}</span>
                <span class="report-test-family">${escapeHtml(option.construct)}</span>
              </span>
              ${badge ? `<span class="report-test-badge">${badge}</span>` : '<span></span>'}
            </label>
          `;
        }).join('')}
      </div>
    `).join('');
    return `
      <section class="report-test-group" data-rw-test-group="${escapeAttr(testId)}">
        <div class="report-test-group-head">
          <span title="${escapeAttr(group.longName || group.name)}">${escapeHtml(group.name)}</span>
          <span class="report-test-group-actions">
            ${coreAvailable ? `<button type="button" class="report-test-add" data-rw-add-test-core="${escapeAttr(testId)}" title="Add the curated core measures for ${escapeHtml(reportDomain(domainId).label)}">+ Core (${coreAvailable})</button>` : ''}
            ${allAvailable ? `<button type="button" class="report-test-add" data-rw-add-test-all="${escapeAttr(testId)}" title="Add every mapped measure from this test for ${escapeHtml(reportDomain(domainId).label)}">+ All (${allAvailable})</button>` : ''}
          </span>
        </div>
        ${roleHtml}
      </section>
    `;
  }).join('');

  browser.innerHTML = `${limitNote}${sectionsHtml}`;
}
function renderReportActivePanel(){
  const section = reportActiveSection();
  const title = document.getElementById('rw-active-title');
  const copy = document.getElementById('rw-active-copy');
  const kicker = document.getElementById('rw-active-kicker');
  if (!title || !copy || !kicker) return;
  if (!section){
    kicker.textContent = 'Selected Section';
    title.textContent = 'Add a section to begin';
    copy.textContent = 'Choose a mapped ability above, then add tests/trials that belong in that section.';
    return;
  }
  kicker.textContent = 'Selected Section';
  title.textContent = reportSectionHeading(section);
  copy.textContent = `Mapped ability wording: ${reportDomain(section.domain).ability}.`;
}
/* ---------- Subtab strip (Build / Scores tabs) ---------- */
function renderReportSubtabs(){
  const strips = [
    { id:'rw-build-subtabs', prefix:'2' },
    { id:'rw-scores-subtabs', prefix:'3' }
  ];
  strips.forEach(({ id, prefix }) => {
    const el = document.getElementById(id); if (!el) return;
    if (!reportState.sections.length){
      el.innerHTML = '<span class="rw-subtabs-empty">No sections yet - add one in Setup.</span>';
      return;
    }
    el.innerHTML = reportState.sections.map((section, i) => {
      const count = section.rows.length;
      const isEmpty = count === 0;
      const active = section.id === reportState.activeId;
      return `<button type="button" class="rw-subtab ${active ? 'active' : ''} ${isEmpty ? 'is-empty' : ''}" data-rw-subtab="${section.id}" title="${escapeAttr(reportSectionHeading(section))}">
        <span class="rw-subtab-num">${prefix}.${i + 1}</span>
        <span class="rw-subtab-label">${escapeHtml(reportSectionHeading(section))}</span>
        <span class="rw-subtab-count">${count}</span>
        <span class="rw-subtab-remove" data-rw-remove-subtab="${section.id}" title="Delete domain" aria-label="Delete ${escapeAttr(reportSectionHeading(section))}">×</span>
      </button>`;
    }).join('');
  });
  const progress = document.getElementById('rw-build-progress');
  if (progress){
    if (!reportState.sections.length){
      progress.textContent = '';
    } else {
      const idx = Math.max(0, reportState.sections.findIndex(s => s.id === reportState.activeId));
      progress.textContent = `Section ${idx + 1} of ${reportState.sections.length}`;
    }
  }
  const setupProg = document.getElementById('rw-setup-progress');
  if (setupProg){
    const n = reportState.sections.length;
    setupProg.textContent = n ? `${n} section${n === 1 ? '' : 's'} ready` : 'Add a section to continue';
  }
}

/* ---------- Tab 2: Selected measures (compact list, no scores) ---------- */
function renderReportSelectedList(){
  const list = document.getElementById('rw-selected-list');
  if (!list) return;
  const section = reportActiveSection();
  if (!section || section.rows.length === 0){
    list.innerHTML = `<div class="report-selected-empty">${section ? 'No measures yet - pick from the list above, or add a custom row.' : 'Add a section first in Setup.'}</div>`;
    return;
  }
  list.innerHTML = section.rows.map(row => {
    const tooltipParts = [row.source, row.construct].filter(Boolean).map(s => reportCleanText(s)).filter(Boolean);
    const title = tooltipParts.join(' · ');
    return `
      <div class="report-selected-row" data-rw-row="${row.id}" draggable="true" title="${escapeAttr(title)}">
        <span class="report-selected-grip" aria-hidden="true">≡</span>
        <span class="report-selected-name">${escapeHtml(row.test)}</span>
        <span class="report-selected-tools">
          <button type="button" class="report-selected-remove" data-rw-remove-row="${row.id}" title="Remove">×</button>
        </span>
      </div>
    `;
  }).join('');
}

/* ---------- Tab 3: Score-entry table ---------- */
function renderReportScoresTable(){
  const tbody = document.querySelector('#rw-scores-table tbody');
  const titleEl = document.getElementById('rw-scores-title');
  const section = reportActiveSection();
  if (!tbody) return;
  if (titleEl) titleEl.textContent = section ? reportSectionHeading(section) : 'Add measures first';
  if (!section){
    tbody.innerHTML = '<tr><td colspan="8" class="computed muted">Add a section in Setup.</td></tr>';
    return;
  }
  if (section.rows.length === 0){
    tbody.innerHTML = '<tr><td colspan="8" class="computed muted">No measures in this section yet - add some in Build.</td></tr>';
    return;
  }
  tbody.innerHTML = section.rows.map(row => {
    reportUpdateComputed(row);
    return `
      <tr data-rw-row="${row.id}" draggable="true">
        <td class="report-scores-grip" aria-hidden="true">≡</td>
        <td class="report-scores-test" title="${escapeAttr(row.source || row.test)}">${escapeHtml(row.test)}</td>
        <td><select data-rw-k="scoreType">${REPORT_SCORE_TYPES.map(t => `<option value="${t.id}" ${row.scoreType === t.id ? 'selected' : ''}>${escapeHtml(t.label)}</option>`).join('')}</select></td>
        <td><input type="number" step="any" data-rw-k="score" data-rw-col="score" value="${escapeAttr(row.score || '')}"></td>
        <td><input type="text" data-rw-k="percentile" data-rw-col="percentile" value="${escapeAttr(row.percentile || '')}"></td>
        <td><input type="text" data-rw-k="ci" data-rw-col="ci" value="${escapeAttr(row.ci || '')}" placeholder="e.g. 95–107"></td>
        <td class="report-scores-descriptor"><em>${escapeHtml((row.descriptor || '').toLowerCase())}</em></td>
        <td class="row-actions"><button type="button" data-rw-remove-row="${row.id}" title="Remove">×</button></td>
      </tr>
    `;
  }).join('');
}
function renderReportOutput(force){
  const out = document.getElementById('rw-output');
  if (!out) return;
  if (reportState.prefs.locked && !force){
    if (!out.innerHTML.trim() && reportState.outputHtml) out.innerHTML = reportState.outputHtml;
    reportState.stale = true;
    out.classList.add('is-stale');
    return;
  }
  reportState.stale = false;
  out.classList.remove('is-stale');
  const parts = [];
  reportState.sections.forEach(section => {
    const rows = section.rows.filter(row => reportCleanText(row.test));
    if (!rows.length) return;
    const scored = reportFilterScoredRows(rows);
    if (!scored.length) return;
    parts.push(reportRenderSectionHtml(section, rows, scored));
  });
  // Summary at the end of the report (Overall + top + bottom)
  const summary = reportBuildSummary();
  if (summary) parts.push(summary);
  out.innerHTML = parts.length ? parts.join('') : `<p>Add a heading, select tests/trials, and enter scores to generate descriptive report text.</p>`;
  reportState.outputHtml = out.innerHTML;
}
function renderReportStructure(){
  renderReportDomainChecklist();
  renderReportSectionList();
  renderReportActivePanel();
  renderReportSubtabs();
  renderReportTestOptions();
  renderReportSelectedList();
  renderReportScoresTable();
  reportUpdateBuildInstruction();
  reportUpdateWizardNav();
}
function renderReportWriter(){
  renderReportStructure();
  renderReportOutput();
  reportSaveState();
}
// Add one section for the given domain id and re-render. Triggered by ticking
// a checkbox on the Setup tab - the just-added domain falls out of the
// checklist (renderReportDomainChecklist filters out already-used domains).
function reportAddDomainSection(domain){
  if (!domain) return;
  const section = { id:`rw-sec-${reportState.nextSectionId++}`, domain, rows:[] };
  reportState.sections.push(section);
  reportState.activeId = section.id;
  reportResetPickerState(true);
  renderReportWriter();
}
function reportSectionHeading(section){
  return reportDomain(section.domain).label;
}
/* ---------- Sort sections by clinical-convention domain order ---------- */
const REPORT_DOMAIN_CONVENTION_ORDER = [
  'intellectual','attention_working_memory','processing_speed',
  'language','visuospatial','verbal_learning_memory','visual_learning_memory',
  'executive','mood','other'
];
function reportSortByConvention(){
  reportState.sections.sort((a, b) => {
    const ai = REPORT_DOMAIN_CONVENTION_ORDER.indexOf(a.domain);
    const bi = REPORT_DOMAIN_CONVENTION_ORDER.indexOf(b.domain);
    return (ai === -1 ? 999 : ai) - (bi === -1 ? 999 : bi);
  });
  renderReportWriter();
  showToast('✓ Sorted by convention');
}

function reportSwitchTab(name){
  const root = document.getElementById('report-writer');
  if (!root) return;
  reportState.activeTab = name;
  const tabs = root.querySelectorAll('.rw-tabs .tab[data-rw-tab]');
  const panels = root.querySelectorAll('.rw-tab-content');
  tabs.forEach(t => t.classList.toggle('active', t.dataset.rwTab === name));
  panels.forEach(p => p.classList.toggle('active', p.id === `rw-tab-${name}`));
  reportUpdateWizardNav();
  reportUpdateBuildInstruction();
}

/* ---------- Strictly sequential wizard walker ---------- */
// Step layout: 0 = Setup, 1..N = Build / section i, N+1..2N = Scores / section i.
function reportTotalSteps(){ return 1 + 2 * reportState.sections.length; }
function reportCurrentStepIdx(){
  const N = reportState.sections.length;
  if (reportState.activeTab === 'setup') return 0;
  const sIdx = reportState.sections.findIndex(s => s.id === reportState.activeId);
  const idx = sIdx === -1 ? 0 : sIdx;
  if (reportState.activeTab === 'build') return 1 + idx;
  if (reportState.activeTab === 'scores') return 1 + N + idx;
  return 0;
}
function reportApplyStepIdx(stepIdx){
  const N = reportState.sections.length;
  if (stepIdx === 0){
    reportSwitchTab('setup');
  } else if (stepIdx >= 1 && stepIdx <= N){
    const section = reportState.sections[stepIdx - 1];
    if (section) reportState.activeId = section.id;
    reportResetPickerState(false);
    reportSwitchTab('build');
    renderReportStructure();
  } else if (stepIdx >= N + 1 && stepIdx <= 2 * N){
    const section = reportState.sections[stepIdx - N - 1];
    if (section) reportState.activeId = section.id;
    reportSwitchTab('scores');
    renderReportStructure();
  }
  reportSaveState();
}
function reportWalkNext(){
  const total = reportTotalSteps();
  const current = reportCurrentStepIdx();
  if (current + 1 < total){
    reportApplyStepIdx(current + 1);
  } else if (current === 0 && reportState.sections.length === 0){
    showToast('Add a section first', true);
  }
}
function reportWalkPrev(){
  const current = reportCurrentStepIdx();
  if (current > 0) reportApplyStepIdx(current - 1);
}
// Friendly label for any step idx ("Setup", "Build: X", "Scores: X")
function reportStepLabel(idx){
  const N = reportState.sections.length;
  if (idx <= 0) return 'Setup';
  if (idx <= N){
    const section = reportState.sections[idx - 1];
    return section ? `Build: ${reportSectionHeading(section)}` : 'Build';
  }
  const sIdx = idx - N - 1;
  const section = reportState.sections[sIdx];
  return section ? `Scores: ${reportSectionHeading(section)}` : 'Scores';
}

function reportUpdateWizardNav(){
  const nav = document.getElementById('rw-wizard-nav');
  if (!nav) return;
  const N = reportState.sections.length;
  const total = reportTotalSteps();
  const current = reportCurrentStepIdx();
  const prevBtn = nav.querySelector('[data-rw-wizard="prev"]');
  const nextBtn = nav.querySelector('[data-rw-wizard="next"]');
  const progress = nav.querySelector('#rw-wizard-progress');
  if (prevBtn) prevBtn.disabled = current === 0;
  if (nextBtn) nextBtn.disabled = (current === total - 1) || (current === 0 && N === 0);

  // Update Back/Next button destination labels
  if (prevBtn){
    const lbl = prevBtn.querySelector('.rw-wizard-prev-label');
    if (lbl){
      lbl.textContent = current > 0 ? 'Back: ' + reportStepLabel(current - 1) : 'Back';
    }
  }
  if (nextBtn){
    const lbl = nextBtn.querySelector('.rw-wizard-next-label');
    if (lbl){
      lbl.textContent = (current < total - 1) ? 'Next: ' + reportStepLabel(current + 1) : 'Next';
    }
  }

  // Progress label - concise step-of-total summary
  if (progress){
    let label;
    if (current === 0){
      label = N ? `Step 1 of ${total} · Setup (${N} section${N === 1 ? '' : 's'} ready)` : `Step 1 of ${total} · Setup (add a section to continue)`;
    } else {
      label = `Step ${current + 1} of ${total} · ${reportStepLabel(current)}`;
    }
    progress.textContent = label;
  }
}

/* ---------- Build tab instruction line ---------- */
function reportUpdateBuildInstruction(){
  const el = document.getElementById('rw-build-instruction');
  if (!el) return;
  const section = reportActiveSection();
  if (!section){
    el.innerHTML = 'Add a section in <strong>Setup</strong> first.';
    return;
  }
  el.innerHTML = `Add tests to the <strong>${escapeHtml(reportSectionHeading(section))}</strong> domain.`;
}

/* ---------- DEV: autofill plausible random scores across every row ---------- */
function reportRandomZ(){
  // Box–Muller, clamped to ±2.5 SD for sane scores
  const u1 = Math.random(), u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1 || 1e-9)) * Math.cos(2 * Math.PI * u2);
  return Math.max(-2.5, Math.min(2.5, z));
}
function reportRandomScore(scoreType){
  const z = reportRandomZ();
  switch(scoreType){
    case 'standard':   return String(Math.max(40, Math.round(100 + 15 * z)));
    case 'scaled':     return String(Math.max(1, Math.round(10 + 3 * z)));
    case 't':          return String(Math.max(10, Math.round(50 + 10 * z)));
    case 'z':          return z.toFixed(1);
    case 'percentile': return String(Math.max(1, Math.min(99, Math.round(normCDF(z) * 100))));
    default:           return String(Math.max(1, Math.round(10 + 3 * z)));
  }
}
function reportAutofillAllScores(){
  let filled = 0;
  reportState.sections.forEach(section => {
    section.rows.forEach(row => {
      const type = row.scoreType || 'scaled';
      if (type === 'percentile') row.percentile = reportRandomScore(type);
      else { row.score = reportRandomScore(type); row.percentile = ''; }
      row.ci = '';
      reportUpdateComputed(row);
      filled += 1;
    });
  });
  renderReportWriter();
  showToast(filled ? `✓ Autofilled ${filled} row${filled === 1 ? '' : 's'}` : 'No rows to fill', !filled);
}

/* ---------- Reset all (clear sections, prefs, localStorage) ---------- */
function reportResetAll(){
  if (typeof confirm === 'function' && !confirm('Clear all sections, measures, and saved state? This cannot be undone.')) return;
  reportState.sections = [];
  reportState.activeId = null;
  reportState.nextSectionId = 1;
  reportState.nextRowId = 1;
  reportState.outputHtml = '';
  reportState.selectedOptionKeys = new Set();
  reportState.visibleOptionKeys = [];
  reportState.activeFilter = 'common';
  reportState.prefs = { ...REPORT_DEFAULT_PREFS };
  try { localStorage.removeItem(REPORT_STORAGE_KEY); } catch(_){}
  reportApplyPrefsToDom();
  reportSwitchTab('setup');
  renderReportWriter();
  showToast('✓ Reset');
}
function reportSectionMoveTo(id, targetIndex){
  const idx = reportState.sections.findIndex(s => s.id === id);
  if (idx === -1 || idx === targetIndex || targetIndex < 0 || targetIndex > reportState.sections.length) return;
  const arr = reportState.sections;
  const [moved] = arr.splice(idx, 1);
  // Adjust target index if removing from before it
  const insertAt = targetIndex > idx ? targetIndex - 1 : targetIndex;
  arr.splice(insertAt, 0, moved);
  renderReportWriter();
}
function reportSectionMove(id, direction){
  const idx = reportState.sections.findIndex(s => s.id === id);
  if (idx === -1) return;
  const swap = direction === 'up' ? idx - 1 : idx + 1;
  if (swap < 0 || swap >= reportState.sections.length) return;
  const arr = reportState.sections;
  [arr[idx], arr[swap]] = [arr[swap], arr[idx]];
  renderReportWriter();
}
function reportAddOption(option, settings = {}){
  const section = reportActiveSection();
  if (!section || !option) return false;
  if (reportOptionAlreadyAdded(section, option.key)) return false;
  section.rows.push({
    id:`rw-row-${reportState.nextRowId++}`,
    sourceKey:option.key || '',
    test:option.displayName || option.subtest || option.label,
    source:option.label,
    subtest:option.subtest,
    construct:option.construct,
    scoreType:option.scoreType || 'scaled',
    score:'', percentile:'', ci:'',
    descriptor:'',
    testId: option.testId || '',
    testName: option.testName || '',
    testLongName: option.testLongName || '',
    role: option.role || 'measure',
    cluster: option.groupId || 'other'
  });
  if (option.key) reportState.selectedOptionKeys.delete(option.key);
  if (settings.render !== false) renderReportWriter();
  return true;
}
function reportOptionsForKeys(keys){
  const byKey = new Map(reportState.options.map(option => [option.key, option]));
  return (keys || []).map(key => byKey.get(key)).filter(Boolean);
}
function reportAddOptions(options){
  const section = reportActiveSection();
  if (!section) return showToast('Add a report section first', true);
  const list = (options || []).filter(Boolean);
  if (!list.length) return showToast('Select one or more mapped tests first', true);
  let added = 0;
  list.forEach(option => { if (reportAddOption(option, { render:false })) added += 1; });
  reportState.selectedOptionKeys.clear();
  renderReportWriter();
  if (!added) return showToast('Those tests are already in this section', true);
  showToast(`✓ Added ${added} mapped test${added === 1 ? '' : 's'}`);
}
function reportAddManual(){
  const section = reportActiveSection();
  if (!section) return showToast('Add a report section first', true);
  const testEl = document.getElementById('rw-manual-test');
  const constructEl = document.getElementById('rw-manual-construct');
  const test = reportCleanText(testEl?.value);
  if (!test) return showToast('Enter a custom test or trial name first', true);
  section.rows.push({
    id:`rw-row-${reportState.nextRowId++}`,
    sourceKey:'',
    test,
    subtest:'',
    construct:reportCleanText(constructEl?.value) || reportDomain(section.domain).ability,
    scoreType:'standard',
    score:'', percentile:'', ci:'',
    descriptor:'',
    testId:'', testName:'', testLongName:'',
    role:'measure',
    cluster:'other'
  });
  if (testEl) testEl.value = '';
  if (constructEl) constructEl.value = '';
  renderReportWriter();
}
function refreshReportWriterOptions(){
  if (!document.getElementById('report-writer')) return;
  reportBuildMeasureOptions();
  renderReportTestOptions();
}

/* ---------- Prefs sync (DOM <-> state) ---------- */
function reportSyncPrefsFromDom(){
  const $ = id => document.getElementById(id);
  const get = (id, key, type = 'value') => {
    const el = $(id); if (!el) return;
    reportState.prefs[key] = type === 'checked' ? !!el.checked : el.value;
  };
  get('rw-reference', 'reference');
  get('rw-descriptor-system', 'descriptorSystem');
  get('rw-lock', 'locked', 'checked');
}
function reportApplyPrefsToDom(){
  const $ = id => document.getElementById(id);
  const set = (id, key, type = 'value') => {
    const el = $(id); if (!el) return;
    if (type === 'checked') el.checked = !!reportState.prefs[key];
    else el.value = reportState.prefs[key];
  };
  set('rw-reference', 'reference');
  set('rw-descriptor-system', 'descriptorSystem');
  set('rw-lock', 'locked', 'checked');
  const out = $('rw-output');
  if (out) out.classList.toggle('is-locked', !!reportState.prefs.locked);
}

function setupReportWriter(){
  if (!document.getElementById('report-writer')) return;
  renderReportDomainChecklist();
  reportBuildMeasureOptions();

  // Hydrate state from localStorage before wiring listeners (avoid save during apply)
  reportState.suppressSave = true;
  reportLoadState();
  reportApplyPrefsToDom();
  reportState.suppressSave = false;
  // Restore the active tab in the DOM (default is 'setup')
  reportSwitchTab(reportState.activeTab || 'setup');

  // Domain checklist - tick to add a section instantly (no button)
  document.getElementById('rw-domain-checklist')?.addEventListener('change', e => {
    const cb = e.target.closest('input[type="checkbox"][data-rw-domain]');
    if (!cb || !cb.checked) return;
    reportAddDomainSection(cb.dataset.rwDomain);
  });

  const sectionList = document.getElementById('rw-section-list');
  sectionList?.addEventListener('click', e => {
    const move = e.target.closest('[data-rw-move]');
    if (move){ reportSectionMove(move.dataset.rwSection, move.dataset.rwMove); return; }
    const remove = e.target.closest('[data-rw-remove-section]');
    if (remove){
      reportState.sections = reportState.sections.filter(section => section.id !== remove.dataset.rwRemoveSection);
      if (reportState.activeId === remove.dataset.rwRemoveSection){
        reportState.activeId = reportState.sections[0]?.id || null;
      }
      reportResetPickerState(true);
      renderReportWriter();
      return;
    }
    // Don't navigate if user clicked on the drag grip
    if (e.target.closest('.report-section-grip')) return;
    const item = e.target.closest('[data-rw-section]');
    if (!item) return;
    reportState.activeId = item.dataset.rwSection;
    reportResetPickerState(true);
    renderReportStructure();
    reportSaveState();
    reportSwitchTab('build'); // jump to Build tab so the user can pick tests
  });
  // Drag-to-reorder for section cards
  let dragSourceId = null;
  sectionList?.addEventListener('dragstart', e => {
    const item = e.target.closest('[data-rw-section]');
    if (!item) return;
    dragSourceId = item.dataset.rwSection;
    item.classList.add('is-dragging');
    e.dataTransfer.effectAllowed = 'move';
    try { e.dataTransfer.setData('text/plain', dragSourceId); } catch(_){}
  });
  sectionList?.addEventListener('dragover', e => {
    if (!dragSourceId) return;
    const item = e.target.closest('[data-rw-section]');
    if (!item || item.dataset.rwSection === dragSourceId) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = 'move';
    sectionList.querySelectorAll('.report-section-item.is-drop-target').forEach(el => {
      if (el !== item) el.classList.remove('is-drop-target');
    });
    item.classList.add('is-drop-target');
  });
  sectionList?.addEventListener('dragleave', e => {
    const item = e.target.closest('[data-rw-section]');
    if (item && !item.contains(e.relatedTarget)) item.classList.remove('is-drop-target');
  });
  sectionList?.addEventListener('drop', e => {
    if (!dragSourceId) return;
    const item = e.target.closest('[data-rw-section]');
    if (!item) return;
    e.preventDefault();
    const targetId = item.dataset.rwSection;
    if (targetId !== dragSourceId){
      const targetIdx = reportState.sections.findIndex(s => s.id === targetId);
      reportSectionMoveTo(dragSourceId, targetIdx);
    }
    sectionList.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    sectionList.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
    dragSourceId = null;
  });
  sectionList?.addEventListener('dragend', () => {
    dragSourceId = null;
    sectionList.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    sectionList.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
  });

  // Tab strip
  document.querySelector('#report-writer .rw-tabs')?.addEventListener('click', e => {
    const tab = e.target.closest('.tab[data-rw-tab]');
    if (!tab) return;
    reportSwitchTab(tab.dataset.rwTab);
  });
  document.querySelector('#report-writer .rw-tabs')?.addEventListener('keydown', e => {
    if (e.key !== 'Enter' && e.key !== ' ') return;
    const tab = e.target.closest('.tab[data-rw-tab]');
    if (!tab) return;
    e.preventDefault();
    reportSwitchTab(tab.dataset.rwTab);
  });

  // Search input + clear-button visibility
  const searchInput = document.getElementById('rw-test-search');
  const searchClear = document.getElementById('rw-test-search-clear');
  function syncSearchClear(){
    if (searchClear) searchClear.hidden = !(searchInput && searchInput.value.length);
  }
  searchInput?.addEventListener('input', () => {
    syncSearchClear();
    renderReportTestOptions();
  });
  searchClear?.addEventListener('click', () => {
    if (!searchInput) return;
    searchInput.value = '';
    syncSearchClear();
    renderReportTestOptions();
    searchInput.focus();
  });
  syncSearchClear();
  document.getElementById('rw-filter-chips')?.addEventListener('click', e => {
    const chip = e.target.closest('[data-rw-filter]');
    if (!chip) return;
    reportState.activeFilter = chip.dataset.rwFilter || 'all';
    renderReportTestOptions();
  });
  document.getElementById('rw-test-list')?.addEventListener('change', e => {
    const checkbox = e.target.closest('[data-rw-option]');
    if (!checkbox) return;
    const section = reportActiveSection();
    if (!section){
      checkbox.checked = false;
      return showToast('Add a report section first', true);
    }
    const key = checkbox.dataset.rwOption;
    if (checkbox.checked){
      // Tick → add immediately to the active section
      const option = reportOptionsForKeys([key])[0];
      reportState.selectedOptionKeys.delete(key);
      if (option) reportAddOption(option); // calls renderReportWriter
    } else {
      // Untick → remove all rows in this section that came from this option
      section.rows = section.rows.filter(row => row.sourceKey !== key);
      reportState.selectedOptionKeys.delete(key);
      renderReportWriter();
    }
  });
  document.getElementById('rw-test-list')?.addEventListener('click', e => {
    const section = reportActiveSection();
    if (!section) return;
    const addCoreBtn = e.target.closest('[data-rw-add-test-core]');
    if (addCoreBtn){
      const testId = addCoreBtn.dataset.rwAddTestCore;
      const keys = [...reportTestCoreKeysForDomain(testId, section.domain)]
        .filter(k => !reportOptionAlreadyAdded(section, k));
      if (!keys.length) return showToast('Core measures already in this section', true);
      reportAddOptions(reportOptionsForKeys(keys));
      return;
    }
    const addAllBtn = e.target.closest('[data-rw-add-test-all]');
    if (addAllBtn){
      const testId = addAllBtn.dataset.rwAddTestAll;
      const keys = reportState.options
        .filter(o => o.testId === testId && (o.domains.includes(section.domain) || section.domain === 'other'))
        .map(o => o.key)
        .filter(k => !reportOptionAlreadyAdded(section, k));
      if (!keys.length) return showToast('All mapped measures from this test are already added', true);
      reportAddOptions(reportOptionsForKeys(keys));
      return;
    }
  });
  document.getElementById('rw-add-mapped')?.addEventListener('click', () => {
    reportAddOptions(reportOptionsForKeys([...reportState.selectedOptionKeys]));
  });
  document.getElementById('rw-add-visible')?.addEventListener('click', () => {
    const section = reportActiveSection();
    if (!section) return showToast('Add a report section first', true);
    const keys = reportState.visibleOptionKeys.filter(key => !reportOptionAlreadyAdded(section, key));
    if (!keys.length) return showToast('All visible tests are already in this section', true);
    reportAddOptions(reportOptionsForKeys(keys));
  });
  document.getElementById('rw-clear-selection')?.addEventListener('click', () => {
    reportState.selectedOptionKeys.clear();
    renderReportTestOptions();
  });
  document.getElementById('rw-add-manual')?.addEventListener('click', reportAddManual);

  // Tab 2 (Build) - Selected measures list: remove
  const selectedList = document.getElementById('rw-selected-list');
  selectedList?.addEventListener('click', e => {
    const section = reportActiveSection();
    if (!section) return;
    const btn = e.target.closest('[data-rw-remove-row]');
    if (!btn) return;
    section.rows = section.rows.filter(row => row.id !== btn.dataset.rwRemoveRow);
    renderReportWriter();
  });
  // Drag-to-reorder rows within the active section's selected list
  let dragRowId = null;
  selectedList?.addEventListener('dragstart', e => {
    const row = e.target.closest('[data-rw-row]');
    if (!row) return;
    dragRowId = row.dataset.rwRow;
    row.classList.add('is-dragging');
    e.dataTransfer.effectAllowed = 'move';
    try { e.dataTransfer.setData('text/plain', dragRowId); } catch(_){}
  });
  selectedList?.addEventListener('dragover', e => {
    if (!dragRowId) return;
    const row = e.target.closest('[data-rw-row]');
    if (!row || row.dataset.rwRow === dragRowId) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = 'move';
    selectedList.querySelectorAll('.is-drop-target').forEach(el => { if (el !== row) el.classList.remove('is-drop-target'); });
    row.classList.add('is-drop-target');
  });
  selectedList?.addEventListener('dragleave', e => {
    const row = e.target.closest('[data-rw-row]');
    if (row && !row.contains(e.relatedTarget)) row.classList.remove('is-drop-target');
  });
  selectedList?.addEventListener('drop', e => {
    if (!dragRowId) return;
    const section = reportActiveSection();
    const targetEl = e.target.closest('[data-rw-row]');
    if (!section || !targetEl) return;
    e.preventDefault();
    const targetId = targetEl.dataset.rwRow;
    if (targetId !== dragRowId){
      const arr = section.rows;
      const fromIdx = arr.findIndex(r => r.id === dragRowId);
      const toIdx = arr.findIndex(r => r.id === targetId);
      if (fromIdx !== -1 && toIdx !== -1){
        const [moved] = arr.splice(fromIdx, 1);
        const insertAt = toIdx > fromIdx ? toIdx - 1 : toIdx;
        arr.splice(insertAt, 0, moved);
        renderReportWriter();
      }
    }
    selectedList.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    selectedList.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
    dragRowId = null;
  });
  selectedList?.addEventListener('dragend', () => {
    dragRowId = null;
    selectedList.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    selectedList.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
  });

  // Tab 3 (Scores) - score table inputs / select / remove + column-major Tab key
  const scoresTable = document.getElementById('rw-scores-table');
  scoresTable?.addEventListener('input', e => {
    const rowEl = e.target.closest('[data-rw-row]');
    const key = e.target.dataset.rwK;
    const section = reportActiveSection();
    if (!rowEl || !key || !section) return;
    const row = section.rows.find(item => item.id === rowEl.dataset.rwRow);
    if (!row) return;
    row[key] = e.target.value;
    reportUpdateComputed(row);
    const pct = rowEl.querySelector('[data-rw-k="percentile"]');
    if ((key === 'score' || key === 'scoreType') && pct) pct.value = row.percentile || '';
    const desc = rowEl.querySelector('td.report-scores-descriptor');
    if (desc) desc.innerHTML = `<em>${escapeHtml((row.descriptor || '').toLowerCase())}</em>`;
    renderReportOutput();
    reportSaveState();
  });
  scoresTable?.addEventListener('change', e => {
    if (e.target.matches('select[data-rw-k]')) e.target.dispatchEvent(new Event('input', { bubbles:true }));
  });
  scoresTable?.addEventListener('click', e => {
    const section = reportActiveSection();
    if (!section) return;
    const btn = e.target.closest('[data-rw-remove-row]');
    if (!btn) return;
    section.rows = section.rows.filter(row => row.id !== btn.dataset.rwRemoveRow);
    renderReportWriter();
  });
  // Drag-to-reorder rows in the scores table
  let dragScoreRowId = null;
  scoresTable?.addEventListener('dragstart', e => {
    const tr = e.target.closest('tr[data-rw-row]');
    if (!tr) return;
    // Don't initiate drag from inputs/selects - only from the row itself or the grip cell
    if (e.target.matches('input,select,button')) { e.preventDefault(); return; }
    dragScoreRowId = tr.dataset.rwRow;
    tr.classList.add('is-dragging');
    e.dataTransfer.effectAllowed = 'move';
    try { e.dataTransfer.setData('text/plain', dragScoreRowId); } catch(_){}
  });
  scoresTable?.addEventListener('dragover', e => {
    if (!dragScoreRowId) return;
    const tr = e.target.closest('tr[data-rw-row]');
    if (!tr || tr.dataset.rwRow === dragScoreRowId) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = 'move';
    scoresTable.querySelectorAll('tr.is-drop-target').forEach(el => { if (el !== tr) el.classList.remove('is-drop-target'); });
    tr.classList.add('is-drop-target');
  });
  scoresTable?.addEventListener('dragleave', e => {
    const tr = e.target.closest('tr[data-rw-row]');
    if (tr && !tr.contains(e.relatedTarget)) tr.classList.remove('is-drop-target');
  });
  scoresTable?.addEventListener('drop', e => {
    if (!dragScoreRowId) return;
    const section = reportActiveSection();
    const tr = e.target.closest('tr[data-rw-row]');
    if (!section || !tr) return;
    e.preventDefault();
    const targetId = tr.dataset.rwRow;
    if (targetId !== dragScoreRowId){
      const arr = section.rows;
      const fromIdx = arr.findIndex(r => r.id === dragScoreRowId);
      const toIdx = arr.findIndex(r => r.id === targetId);
      if (fromIdx !== -1 && toIdx !== -1){
        const [moved] = arr.splice(fromIdx, 1);
        const insertAt = toIdx > fromIdx ? toIdx - 1 : toIdx;
        arr.splice(insertAt, 0, moved);
        renderReportWriter();
      }
    }
    scoresTable.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    scoresTable.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
    dragScoreRowId = null;
  });
  scoresTable?.addEventListener('dragend', () => {
    dragScoreRowId = null;
    scoresTable.querySelectorAll('.is-drop-target').forEach(el => el.classList.remove('is-drop-target'));
    scoresTable.querySelectorAll('.is-dragging').forEach(el => el.classList.remove('is-dragging'));
  });
  // Column-major Tab navigation: Score₁ → Score₂ → … → Pct₁ → Pct₂ → … → CI₁ → CI₂ → …
  scoresTable?.addEventListener('keydown', e => {
    if (e.key !== 'Tab') return;
    const cell = e.target.closest('[data-rw-col]');
    if (!cell) return;
    const cols = ['score','percentile','ci'];
    const currentCol = cell.dataset.rwCol;
    const colIdx = cols.indexOf(currentCol);
    if (colIdx === -1) return;
    const allInCol = (colName) => Array.from(scoresTable.querySelectorAll(`[data-rw-col="${colName}"]`));
    const colNodes = allInCol(currentCol);
    const rowIdx = colNodes.indexOf(cell);
    let nextNode = null;
    if (e.shiftKey){
      if (rowIdx > 0) nextNode = colNodes[rowIdx - 1];
      else if (colIdx > 0){ const prevCol = allInCol(cols[colIdx - 1]); nextNode = prevCol[prevCol.length - 1] || null; }
    } else {
      if (rowIdx < colNodes.length - 1) nextNode = colNodes[rowIdx + 1];
      else if (colIdx < cols.length - 1){ const nextCol = allInCol(cols[colIdx + 1]); nextNode = nextCol[0] || null; }
    }
    if (nextNode){ e.preventDefault(); nextNode.focus(); nextNode.select?.(); }
  });

  // Subtabs (Build / Scores) - clicking switches active section, × deletes the domain
  document.querySelectorAll('#rw-build-subtabs, #rw-scores-subtabs').forEach(strip => {
    strip.addEventListener('click', e => {
      const remove = e.target.closest('[data-rw-remove-subtab]');
      if (remove){
        e.stopPropagation();
        const sectionId = remove.dataset.rwRemoveSubtab;
        reportState.sections = reportState.sections.filter(s => s.id !== sectionId);
        if (reportState.activeId === sectionId){
          reportState.activeId = reportState.sections[0]?.id || null;
        }
        reportResetPickerState(false);
        renderReportWriter();
        return;
      }
      const btn = e.target.closest('[data-rw-subtab]');
      if (!btn) return;
      reportState.activeId = btn.dataset.rwSubtab;
      reportResetPickerState(false);
      renderReportStructure();
      reportSaveState();
    });
  });

  // Persistent wizard nav - strictly sequential prev/next
  document.getElementById('rw-wizard-nav')?.addEventListener('click', e => {
    const btn = e.target.closest('[data-rw-wizard]');
    if (!btn || btn.disabled) return;
    if (btn.dataset.rwWizard === 'prev') reportWalkPrev();
    else reportWalkNext();
  });

  // Reset all
  document.getElementById('rw-reset')?.addEventListener('click', reportResetAll);
  // Sort sections by clinical convention
  document.getElementById('rw-sort-convention')?.addEventListener('click', reportSortByConvention);
  // DEV - autofill scores
  document.getElementById('rw-autofill')?.addEventListener('click', reportAutofillAllScores);

  // Reference (Mrs Doe / Mr Doe / You): regenerates output
  document.getElementById('rw-reference')?.addEventListener('change', () => {
    reportSyncPrefsFromDom(); renderReportOutput(); reportSaveState();
  });
  // Descriptor system: re-render scores table (descriptor column) + output
  document.getElementById('rw-descriptor-system')?.addEventListener('change', () => {
    reportSyncPrefsFromDom(); renderReportScoresTable(); renderReportOutput(); reportSaveState();
  });
  document.getElementById('rw-lock')?.addEventListener('change', e => {
    reportState.prefs.locked = !!e.target.checked;
    const out = document.getElementById('rw-output');
    if (out) out.classList.toggle('is-locked', reportState.prefs.locked);
    if (!reportState.prefs.locked) renderReportOutput();
    reportSaveState();
  });
  document.getElementById('rw-regenerate')?.addEventListener('click', () => {
    renderReportOutput(true);
    reportSaveState();
    showToast('✓ Output regenerated');
  });
  document.getElementById('rw-copy-report')?.addEventListener('click', async () => {
    const out = document.getElementById('rw-output');
    if (!out || !out.textContent.trim()) return showToast('No report text to copy yet', true);
    try {
      if (navigator.clipboard && window.ClipboardItem){
        await navigator.clipboard.write([
          new ClipboardItem({
            'text/html': new Blob([out.innerHTML], { type:'text/html' }),
            'text/plain': new Blob([out.innerText], { type:'text/plain' })
          })
        ]);
      } else {
        await navigator.clipboard.writeText(out.innerText);
      }
      showToast('✓ Report text copied');
    } catch(e){
      showToast('Copy failed. Select the draft text and copy manually.', true);
    }
  });
  // Save edits made directly inside the contenteditable output
  document.getElementById('rw-output')?.addEventListener('input', () => {
    const out = document.getElementById('rw-output');
    reportState.outputHtml = out ? out.innerHTML : '';
    reportSaveState();
  });

  renderReportWriter();
}

/* ---------- INITIAL POPULATION ---------- */
// Add a starter example row and one blank row to each editable table
batteryRows = [{name:'Example subtest', raw:'25', score:'10'}, {name:'', raw:'', score:''}];
sdiRows = [sdiMode() === 'raw' ? {name:'Example memory score',t1:'42',t2:'36',sd:'8'} : {name:'Example memory score',t1:'9',t2:'6'}, {}];
['rci-basic', 'rci-practice', 'rci-srb', 'rci-crawford'].forEach(m => { rciState[m].rows = [examples[m], {}]; renderRci(m); });
renderBattery();
renderSdi();
applyCalculatorPolish();
enhanceCalculatorWorkflow();
enhanceApaToolbars();
setupReportWriter();

buildDescCarousels();
renderConverter();

/* Premorbid panel - fade fields when comparison is off (kept always visible) */
(function setupBatteryPremorbidDisable(){
  const block = document.getElementById('bat-premorbid-block');
  const checkbox = document.getElementById('bat-prem-enable');
  if (!block || !checkbox) return;
  function sync(){ block.dataset.enabled = String(checkbox.checked); }
  checkbox.addEventListener('change', sync);
  sync();
})();

// Premorbid setup
setupPreTabs();
buildPredictTable();
setupPremorbidListeners();
calcPremorbid();
calcPredict();
calcOpiePredict();

// Battery autofill input wiring (combobox listeners)
wireBatteryAutofill();
wireSdiAutofill();

// Final initialization
refreshAll();

/* ---------- GLOBAL CLEAR — Topbar "Clear all tables" button ----------
   Wipes every tool's session data in one action: battery rows, SDI rows,
   the four RCI methods, premorbid inputs, working report items.

   The auto-add MutationObserver on each APA container is suppressed during
   the clear cascade so the residual table HTML doesn't sneak items back
   into the bundle after we wipe it. */
(function wireGlobalClear(){
  const btn = document.getElementById('topbar-clear-all');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const ok = confirm('Are you sure?\n\nThis clears every table you\'ve been working on across all tools, including the Working Report. It cannot be undone.');
    if (!ok) return;

    // Suppress the observers BEFORE any clearing happens. The 350ms debounce
    // means handlers from changes triggered during the cascade can fire well
    // after we've finished, so we keep suppression on for ~700ms to be safe.
    if (typeof ReportBundle !== 'undefined' && ReportBundle.setSuppressed){
      ReportBundle.setSuppressed(true);
    }

    // Battery / Neuropsych Tables
    try { batteryRows.length = 0; renderBattery(); } catch(e){}
    // SDI
    try { if (typeof clearSdi === 'function') clearSdi(); } catch(e){}
    // All four RCI methods
    try {
      ['rci-basic','rci-practice','rci-srb','rci-crawford'].forEach(m => {
        if (typeof clearMethodRows === 'function') clearMethodRows(m);
      });
    } catch(e){}
    // Premorbid inputs (predictors + demographics)
    try {
      ['pre-topf','pre-vc','pre-mr','pre-sex','pre-occ','pre-edu','pre-age'].forEach(id => {
        const el = document.getElementById(id);
        if (!el) return;
        el.value = '';
        el.dispatchEvent(new Event('input',  { bubbles:true }));
        el.dispatchEvent(new Event('change', { bubbles:true }));
      });
      // Premorbid achieved-score cells (per-index)
      document.querySelectorAll('[data-pre-ach]').forEach(inp => {
        inp.value = '';
        inp.dispatchEvent(new Event('input',  { bubbles:true }));
        inp.dispatchEvent(new Event('change', { bubbles:true }));
      });
      document.querySelectorAll('[data-pre-opie-ach]').forEach(inp => {
        inp.value = '';
        inp.dispatchEvent(new Event('input',  { bubbles:true }));
        inp.dispatchEvent(new Event('change', { bubbles:true }));
      });
    } catch(e){}

    // Working Report bundle - clear after the tools so any pending observer
    // tasks fire against an empty state.
    try {
      if (typeof ReportBundle !== 'undefined' && ReportBundle.clearSilent){
        ReportBundle.clearSilent();
      }
    } catch(e){}

    // Re-clear after the debounce window in case any observer queued an
    // add during the cascade. Then turn observers back on.
    setTimeout(() => {
      try {
        if (typeof ReportBundle !== 'undefined' && ReportBundle.clearSilent){
          ReportBundle.clearSilent();
        }
        if (typeof ReportBundle !== 'undefined' && ReportBundle.setSuppressed){
          ReportBundle.setSuppressed(false);
        }
      } catch(e){}
    }, 700);

    if (typeof showToast === 'function') showToast('All tables cleared');
  });
})();


/* ============================================================
   AUTH OVERLAY · prototype-ready login/register behaviour
   ============================================================ */
(function(){
  const overlay = document.getElementById('auth-overlay');
  if (!overlay) return;

  const SESSION_KEY = 'paAuthPrototypeSession';
  const USER_KEY = 'paAuthPrototypeUser';

  const loginForm = document.getElementById('auth-login-form');
  const registerForm = document.getElementById('auth-register-form');
  const loginFeedback = document.getElementById('auth-login-feedback');
  const registerFeedback = document.getElementById('auth-register-feedback');

  function isValidEmail(value){
    return /^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(String(value || '').trim());
  }

  function setFeedback(el, message, type){
    if (!el) return;
    el.textContent = message || '';
    el.className = 'auth-feedback' + (type ? ' ' + type : '');
  }

  function hideOverlay(){
    overlay.classList.add('is-hidden');
    overlay.setAttribute('aria-hidden', 'true');
  }

  function showOverlay(){
    overlay.classList.remove('is-hidden');
    overlay.removeAttribute('aria-hidden');
  }

  function getUser(){
    try { return JSON.parse(localStorage.getItem(USER_KEY) || 'null'); }
    catch(e){ return null; }
  }

  function activateHomePage(){
    const homeNav = document.querySelector('.nav-item[data-target="home"]');
    const homeSection = document.getElementById('home');
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));
    document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
    if (homeNav) {
      homeNav.classList.add('active');
      if (typeof openOnlyNavGroup === 'function') openOnlyNavGroup(homeNav.closest('.nav-group'));
    }
    if (homeSection) homeSection.classList.add('active');
    const main = document.querySelector('.main');
    if (main) main.scrollTop = 0;
  }

  function setSession(user){
    localStorage.setItem(SESSION_KEY, 'true');
    if (user) localStorage.setItem(USER_KEY, JSON.stringify(user));
    activateHomePage();
    hideOverlay();
    renderAuthChip();
    if (typeof showToast === 'function') showToast('✓ Prototype access granted');
  }

  function clearSession(){
    localStorage.removeItem(SESSION_KEY);
    activateHomePage();
    showOverlay();
    renderAuthChip();
  }

  function renderAuthChip(){
    let chip = document.getElementById('auth-user-chip');
    const sidebar = document.querySelector('.sidebar');
    if (!sidebar) return;

    if (!chip){
      chip = document.createElement('div');
      chip.id = 'auth-user-chip';
      chip.className = 'auth-user-chip';
      sidebar.appendChild(chip);
    }

    const active = localStorage.getItem(SESSION_KEY) === 'true';
    const user = getUser();
    if (!active){
      chip.classList.remove('show');
      chip.innerHTML = '';
      return;
    }

    const label = user && user.email ? user.email : 'Demo mode';
    chip.classList.add('show');
    chip.innerHTML = `
      <span>Signed in as <strong>${escapeHtml(label)}</strong></span>
      <button class="btn btn-ghost" type="button" id="auth-logout">Log out</button>
    `;
    const logout = document.getElementById('auth-logout');
    if (logout) logout.addEventListener('click', clearSession);
  }

  document.querySelectorAll('[data-auth-tab]').forEach(tab => {
    tab.addEventListener('click', () => {
      const target = tab.dataset.authTab;
      document.querySelectorAll('[data-auth-tab]').forEach(t => t.classList.toggle('active', t === tab));
      loginForm.classList.toggle('active', target === 'login');
      registerForm.classList.toggle('active', target === 'register');
      setFeedback(loginFeedback, '');
      setFeedback(registerFeedback, '');
    });
  });

  loginForm.addEventListener('submit', e => {
    e.preventDefault();
    const email = document.getElementById('auth-login-email').value.trim();
    const password = document.getElementById('auth-login-password').value;
    if (!isValidEmail(email)) return setFeedback(loginFeedback, 'Enter a valid email address.', 'error');
    if (!password) return setFeedback(loginFeedback, 'Enter your password.', 'error');

    // TODO later:
    // Replace this local prototype action with your provider call, for example:
    // await supabase.auth.signInWithPassword({ email, password });
    setSession({ email, name: email.split('@')[0], mode: 'prototype-login' });
  });

  registerForm.addEventListener('submit', e => {
    e.preventDefault();
    const name = document.getElementById('auth-register-name').value.trim();
    const email = document.getElementById('auth-register-email').value.trim();
    const password = document.getElementById('auth-register-password').value;
    if (!name) return setFeedback(registerFeedback, 'Enter your name.', 'error');
    if (!isValidEmail(email)) return setFeedback(registerFeedback, 'Enter a valid email address.', 'error');
    if (password.length < 8) return setFeedback(registerFeedback, 'Use at least 8 characters for the password.', 'error');

    // TODO later:
    // Replace this local prototype action with your provider call, for example:
    // await supabase.auth.signUp({ email, password, options:{ data:{ full_name:name } } });
    setSession({ email, name, mode: 'prototype-register' });
  });

  document.getElementById('auth-demo').addEventListener('click', () => {
    setSession({ email: 'Demo mode', name: 'Demo user', mode: 'demo' });
  });

  document.getElementById('auth-forgot').addEventListener('click', () => {
    setFeedback(loginFeedback, 'Password reset is not connected yet. This link is ready to wire to a real provider later.', 'success');
  });

  renderAuthChip();
  activateHomePage();
  if (localStorage.getItem(SESSION_KEY) === 'true') hideOverlay();
  else showOverlay();
})();

/* =====================================================================
   REDESIGN - Top-bar navigation bucket sync
   Each section ID maps to the topnav bucket (data-bucket) that should
   be highlighted when that section is active.
   ===================================================================== */
(function(){
  const TOPNAV_BUCKETS = {
    home: 'home',
    converter: 'converter',
    battery: 'battery',
    'report-writer': 'report',
    effectsize: 'effectsize',
    // All change-analysis methods map to the "change" bucket
    sdi: 'change',
    'rci-basic': 'change',
    'rci-practice': 'change',
    'rci-srb': 'change',
    'rci-crawford': 'change',
    'change-analysis': 'change',
    premorbid: 'premorbid'
    // Norms (custom-tests) and Reference (about) live in the footer, not the top nav
  };

  // Page-title labels for the brand row (mirrors Fisherman's session-title slot)
  const PAGE_TITLES = {
    home: 'Home',
    converter: 'Score Converter',
    battery: 'Neuropsych Report Tables',
    'report-writer': 'Report Writer',
    effectsize: 'Effect Size Tools',
    sdi: 'Standard Deviation Index',
    'rci-basic': 'Simple Reliable Change',
    'rci-practice': 'Practice Effect-Adjusted',
    'rci-srb': 'McSweeney Regression-Based',
    'rci-crawford': 'Crawford Regression-Based',
    'change-analysis': 'Change Analysis',
    premorbid: 'Premorbid Estimation',
    'custom-tests': 'Norms Database',
    about: 'Methods & References'
  };

  function syncTopnav(){
    const active = document.querySelector('.section.active');
    const id = active ? active.id : 'home';
    const bucket = TOPNAV_BUCKETS[id] || 'home';
    // Highlight the matching topnav button (and its parent group, if dropdown)
    document.querySelectorAll('.topnav-item').forEach(b => {
      b.classList.toggle('active', b.dataset.bucket === bucket);
    });
    // Also mark the active dropdown item so it shows highlighted on hover
    document.querySelectorAll('.topnav-drop-item').forEach(b => {
      b.classList.toggle('active', b.dataset.target === id);
    });
    // Mark active footer link (Norms / Methods & References)
    document.querySelectorAll('.site-footer-link').forEach(b => {
      b.classList.toggle('active', b.dataset.target === id);
    });
    // Update the brand-row page title (right of the wordmark divider)
    const titleEl = document.getElementById('topbar-page-title');
    if (titleEl) titleEl.textContent = PAGE_TITLES[id] || '';
  }

  // Watch for any section becoming active → re-sync top-nav
  document.querySelectorAll('.section').forEach(s => {
    new MutationObserver(syncTopnav).observe(s, { attributes:true, attributeFilter:['class'] });
  });

  // Initial sync, plus a small deferred re-sync to catch the JS-built change-analysis section
  syncTopnav();
  setTimeout(syncTopnav, 200);
})();

/* =====================================================================
   Score Converter - view-mode tabs (Equivalents / Distribution)
   ===================================================================== */
(function(){
  document.querySelectorAll('#converter .conv-view-tab').forEach(btn => {
    btn.addEventListener('click', () => {
      const view = btn.dataset.convView;
      document.querySelectorAll('#converter .conv-view-tab').forEach(t => {
        const active = t.dataset.convView === view;
        t.classList.toggle('is-active', active);
        t.setAttribute('aria-selected', String(active));
      });
      document.querySelectorAll('#converter .conv-view-pane').forEach(p => {
        p.classList.toggle('is-active', p.dataset.pane === view);
      });
    });
  });
})();

/* ================================================================
   WORKING REPORT BUNDLE v2
   - Renders full APA tables in a wide drawer (replaces inline panels)
   - Auto-updates via MutationObserver on each tool's APA container
   - Dedupes by sourceId - re-Add refreshes the existing entry
   - Drag-to-reorder items
   - Toggle layout: floating popover ↔ docked side panel
   - Persists across reloads via localStorage
   ================================================================ */
const ReportBundle = (function(){
  const STORAGE_KEY = 'workingReport_v1';
  const SOURCE_LABELS = {
    'bat-apa':           'Neuropsych Tables',
    'sdi-apa':           'Standard Deviation Index',
    'rci-basic-apa':     'Simple Reliable Change',
    'rci-practice-apa':  'Practice-Adjusted RCI',
    'rci-srb-apa':       'McSweeney Regression-Based',
    'rci-crawford-apa':  'Crawford Regression-Based',
    'pre-estimates-apa':    'Premorbid · Estimates',
    'pre-predict-apa':      'Premorbid · ToPF Predicted',
    'pre-opiepredict-apa':  'Premorbid · OPIE-4 Predicted'
  };
  /* Method / tool names - combined with the detected test family to produce
     intelligent table titles like "Crawford Regression-Based Change: WAIS-IV". */
  const SOURCE_METHOD_NAMES = {
    'bat-apa':              'Cognitive Outcomes',
    'sdi-apa':              'Standard-Deviation Discrepancy',
    'rci-basic-apa':        'Reliable Change (Jacobson & Truax)',
    'rci-practice-apa':     'Practice-Adjusted Reliable Change',
    'rci-srb-apa':          'McSweeney Regression-Based Change',
    'rci-crawford-apa':     'Crawford Regression-Based Change',
    'pre-estimates-apa':    'Premorbid Cognitive Estimate',
    'pre-predict-apa':      'ToPF-Predicted vs Achieved',
    'pre-opiepredict-apa':  'OPIE-4-Predicted vs Achieved'
  };
  /* Backwards alias - SOURCE_TITLES still referenced in a couple of places */
  const SOURCE_TITLES = SOURCE_METHOD_NAMES;
  /* Build an intelligent title:  "Method: Family"  if family detected, else just method name.
     For split items, pass the family name explicitly via the second arg. */
  function buildIntelligentTitle(sourceId, html, explicitFamily){
    // Strip the "::groupname" suffix to get the parent's method
    const parentId = (sourceId || '').split('::')[0];
    const method = SOURCE_METHOD_NAMES[parentId] || null;

    // Premorbid sources: the table CONTENT mentions WAIS-IV / WMS-IV (those are
    // the predicted indices), but those aren't the "test family" - the test is
    // ToPF / OPIE, already named in the method. Skip family detection here so
    // we don't end up with "ToPF-Predicted vs Achieved: WAIS-IV".
    if (parentId && parentId.startsWith('pre-')){
      return method || 'APA Table';
    }

    const family = explicitFamily || detectTestFamily(html);
    if (method && family) return `${method}: ${family}`;
    if (method) return method;
    return family || 'APA Table';
  }
  const SOURCE_IDS = Object.keys(SOURCE_LABELS);

  let state = { items: [], minimized: true, drawerWidth: null, onboardingSeen: false, maximised: false };
  /* Last seen body-row count per source. Used to decide whether the pill
     should fire - we only want it on actual row additions, not score edits. */
  const lastRowCount = {};
  function countTableRows(container){
    return container ? container.querySelectorAll('table tbody tr:not(.apa-group)').length : 0;
  }
  let rootEl = null;
  const observers = new Map(); // sourceId -> MutationObserver
  let dragId = null;
  let lastChangedItemId = null; // for auto-scroll on render
  const KOFI_SEEN_KEY = 'kofiPromptSeen_v1'; // localStorage flag - persists across tabs/sessions
  let kofiToastShown = false;   // in-memory guard for current page load
  function hasSeenKofi(){
    try { return !!localStorage.getItem(KOFI_SEEN_KEY); } catch(e){ return false; }
  }
  function markKofiSeen(){
    try { localStorage.setItem(KOFI_SEEN_KEY, String(Date.now())); } catch(e){}
  }

  /* Post-export thank-you modal. Auto-fires the FIRST time EVER (across tabs
     and sessions) that the user exports/copies a table or spends 5 minutes
     on the site. Can also be opened manually (e.g. from the header "Buy me a
     coffee" button) by passing { force: true }, which bypasses the seen flag.
     Persisted via localStorage. Center-screen with a dimming backdrop.
     Dismissable via backdrop click, the X, or Escape. */
  function maybeShowKofiToast(opts){
    const force = !!(opts && opts.force);
    // Don't double-render while a modal is already on screen
    if (kofiToastShown) return;
    // Auto-triggers (export, copy, timer) respect the localStorage flag.
    // Force-opens (header button) always proceed regardless.
    if (!force && hasSeenKofi()) return;
    kofiToastShown = true;
    // Only mark the persistent "seen" flag for AUTO-triggers - that's what the
    // flag is for (don't repeatedly nudge passive users). Manual header-button
    // clicks are user-initiated and shouldn't burn the budget.
    if (!force) markKofiSeen();

    const card = document.createElement('aside');
    card.className = 'rb-kofi-toast';
    card.setAttribute('role', 'status');
    card.setAttribute('aria-live', 'polite');
    card.setAttribute('aria-label', 'Support note');
    card.innerHTML = `
      <button class="rb-kofi-toast-close" type="button" aria-label="Dismiss">
        <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" aria-hidden="true"><path d="M3 3l6 6M9 3l-6 6"/></svg>
      </button>
      <div class="rb-kofi-toast-body">
        <span class="rb-kofi-toast-icon" aria-hidden="true">
          <svg viewBox="0 0 24 24" fill="none" stroke="#FFFFFF" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">
            <path d="M11.5 4.5c-1 0-2 .3-2.7 1c-1 .2-1.8 1-2 2c-1.5.4-2.6 1.7-2.6 3.2c0 .5.1 1 .3 1.4c-.6.5-1 1.3-1 2.2c0 1 .5 1.9 1.3 2.4c-.2.4-.3.8-.3 1.3c0 1.7 1.4 3.1 3.1 3.1c.5 0 .9-.1 1.3-.3c.6.5 1.4.7 2.2.7c.6 0 1.2-.2 1.7-.5"/>
            <path d="M12.5 4.5c1 0 2 .3 2.7 1c1 .2 1.8 1 2 2c1.5.4 2.6 1.7 2.6 3.2c0 .5-.1 1-.3 1.4c.6.5 1 1.3 1 2.2c0 1-.5 1.9-1.3 2.4c.2.4.3.8.3 1.3c0 1.7-1.4 3.1-3.1 3.1c-.5 0-.9-.1-1.3-.3c-.6.5-1.4.7-2.2.7c-.6 0-1.2-.2-1.7-.5"/>
            <line x1="12" y1="5" x2="12" y2="20"/>
            <path d="M8 11c1-.6 2-.6 3 0" stroke-opacity="0.7"/>
            <path d="M16 11c-1-.6-2-.6-3 0" stroke-opacity="0.7"/>
          </svg>
        </span>
        <div class="rb-kofi-toast-msg">
          <strong class="rb-kofi-toast-title">Saved you time today?</strong>
          <span class="rb-kofi-toast-sub">Free and ad&#8209;free for clinicians. Contributions help keep it that way.</span>
        </div>
      </div>
      <div class="rb-kofi-toast-foot">
        <a class="rb-kofi-toast-cta" href="https://ko-fi.com/clinpsyry" target="_blank" rel="noopener noreferrer">
          <svg class="rb-kofi-toast-cta-cup" viewBox="0 0 18 18" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
            <path d="M3 7.5h10V13a3 3 0 0 1-3 3H6a3 3 0 0 1-3-3V7.5z"/>
            <path d="M13 9h1.5a2 2 0 0 1 0 4H13"/>
            <path d="M5.5 3.5v1.4M8 3v1.4M10.5 3.5v1.4"/>
          </svg>
          Support this work
        </a>
        <span class="rb-kofi-toast-credit">via ko&#8209;fi.com</span>
      </div>
    `;
    document.body.appendChild(card);
    requestAnimationFrame(() => card.classList.add('is-visible'));

    let autoDismiss = setTimeout(dismiss, 14000);
    function dismiss(){
      clearTimeout(autoDismiss);
      card.classList.remove('is-visible');
      setTimeout(() => {
        card.remove();
        kofiToastShown = false;
      }, 280);
    }

    /* Pause auto-dismiss while the user's pointer is over the toast,
       so they have time to read or click without it disappearing. */
    card.addEventListener('mouseenter', () => clearTimeout(autoDismiss));
    card.addEventListener('mouseleave', () => { autoDismiss = setTimeout(dismiss, 6000); });

    card.querySelector('.rb-kofi-toast-close').addEventListener('click', e => {
      e.preventDefault(); e.stopPropagation(); dismiss();
    });
    card.querySelector('.rb-kofi-toast-cta').addEventListener('click', () => {
      setTimeout(dismiss, 200);
    });
  }

  /* ---------- persistence ---------- */
  function load(){
    try {
      const raw = localStorage.getItem(STORAGE_KEY);
      if (raw){
        const parsed = JSON.parse(raw);
        const w = Number(parsed.drawerWidth);
        state = {
          items: Array.isArray(parsed.items) ? parsed.items : [],
          minimized: parsed.minimized !== false,
          drawerWidth: (Number.isFinite(w) && w >= 320) ? w : null,
          onboardingSeen: parsed.onboardingSeen === true,
          maximised: parsed.maximised === true
        };
      }
    } catch(e){ state = { items:[], minimized:true, drawerWidth:null, onboardingSeen:false, maximised:false }; }
  }
  function save(){
    try { localStorage.setItem(STORAGE_KEY, JSON.stringify(state)); } catch(e){}
  }

  /* ---------- helpers ---------- */
  function newId(){ return 'rb-' + Date.now().toString(36) + '-' + Math.random().toString(36).slice(2,7); }
  function escapeHtmlLocal(s){
    return String(s == null ? '' : s)
      .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')
      .replace(/"/g,'&quot;').replace(/'/g,'&#39;');
  }
  function debounce(fn, wait){
    let t = null;
    return function(){
      const args = arguments;
      clearTimeout(t);
      t = setTimeout(() => fn.apply(this, args), wait);
    };
  }
  function formatRelative(iso){
    if (!iso) return '';
    const t = new Date(iso).getTime();
    const now = Date.now();
    const sec = Math.max(1, Math.round((now - t) / 1000));
    if (sec < 60) return 'just now';
    const min = Math.round(sec / 60);
    if (min < 60) return `${min} min ago`;
    const hr = Math.round(min / 60);
    if (hr < 24) return `${hr} hr ago`;
    return `${Math.round(hr / 24)}d ago`;
  }
  function getTitleFromContainer(container, fallback){
    const titleEl = container.querySelector('.apa-table-title');
    const numEl = container.querySelector('.apa-table-num');
    return (titleEl?.textContent || numEl?.textContent || fallback || 'APA Table').trim();
  }
  /* Replace the captured table's title text with an intelligent one of the form
     "Method: Test family" (e.g. "Crawford Regression-Based Change: WAIS-IV"). */
  function applyMeaningfulTitle(html, sourceId){
    const title = buildIntelligentTitle(sourceId, html);
    if (!title) return html;
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    const titleEl = tmp.querySelector('.apa-table-title');
    if (titleEl) titleEl.textContent = title;
    return tmp.innerHTML;
  }
  /* Read the column labels from a captured APA table for the column-toggle UI.
     Uses the LAST header row (skipping spanner rows above it). */
  function getItemColumns(item){
    const tmp = document.createElement('div');
    tmp.innerHTML = item.html;
    const headerRows = tmp.querySelectorAll('table thead tr');
    if (!headerRows.length) return [];
    const lastRow = headerRows[headerRows.length - 1];
    return [...lastRow.children].map((cell, idx) => ({
      idx,
      label: cell.textContent.replace(/\s+/g, ' ').trim() || `Column ${idx + 1}`
    }));
  }
  /* Apply the per-item hiddenColumns filter to a captured-HTML string by
     setting display:none on cells at hidden indices. Inline style so this
     also carries through to clipboard / Word export. Skips spanner rows
     in <thead> (we only target the last header row, where actual column
     labels live, plus all body rows). */
  function applyHiddenColumns(html, hiddenColumns){
    if (!hiddenColumns || !hiddenColumns.length) return html;
    const hidden = new Set(hiddenColumns);
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    const hideCell = cell => {
      const existing = cell.getAttribute('style') || '';
      if (!/display\s*:\s*none/i.test(existing)){
        cell.setAttribute('style', existing + ';display:none;');
      }
    };
    tmp.querySelectorAll('table').forEach(table => {
      const headRows = table.querySelectorAll('thead tr');
      if (headRows.length){
        const lastHead = headRows[headRows.length - 1];
        [...lastHead.children].forEach((cell, idx) => {
          if (hidden.has(idx)) hideCell(cell);
        });
      }
      table.querySelectorAll('tbody tr').forEach(tr => {
        [...tr.children].forEach((cell, idx) => {
          if (hidden.has(idx)) hideCell(cell);
        });
      });
    });
    return tmp.innerHTML;
  }
  /* Apply per-item header label overrides to the LAST <thead> row of the
     captured HTML. Skips entries that are null/empty so the original captured
     text shows through. */
  function applyHeaderOverrides(html, overrides){
    if (!overrides || !overrides.length) return html;
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    const headerRows = tmp.querySelectorAll('table thead tr');
    if (!headerRows.length) return html;
    const lastRow = headerRows[headerRows.length - 1];
    [...lastRow.children].forEach((cell, idx) => {
      const v = overrides[idx];
      if (v != null && String(v).trim() !== '') cell.textContent = v;
    });
    return tmp.innerHTML;
  }

  /* Replace the captured "Table 1" number text with the item's position in
     the report (Table 1 / Table 2 / Table 3 …). */
  function renumberTable(html, num){
    if (!Number.isFinite(num)) return html;
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    const numEl = tmp.querySelector('.apa-table-num');
    if (numEl) numEl.textContent = `Table ${num}`;
    return tmp.innerHTML;
  }
  function effectiveItemHtml(item, indexInReport){
    let html = applyHiddenColumns(item.html, item.hiddenColumns || []);
    html = applyHeaderOverrides(html, item.headerOverrides || []);
    if (typeof indexInReport === 'number'){
      html = renumberTable(html, indexInReport + 1);
    }
    return html;
  }

  /* Detect the test family name (CVLT-3, WAIS-IV, etc.) from a captured APA
     table's content. Used to label the "Added to report" pill with the
     actual test rather than the analysis type. */
  const TEST_FAMILY_PATTERNS = [
    'CVLT-3', 'CVLT-II', 'CVLT',
    'D-KEFS',
    'RBANS',
    'WAIS-IV', 'WAIS-V', 'WAIS-III', 'WAIS',
    'WMS-IV', 'WMS-V', 'WMS-III', 'WMS',
    'WISC-V', 'WISC-IV', 'WISC',
    'TOMM',
    'ToPF',
    'OPIE-4', 'OPIE-3', 'OPIE'
  ];
  function detectTestFamily(html){
    const tmp = document.createElement('div');
    tmp.innerHTML = html || '';
    // 1. Use the LAST group separator row - that's the most recently added
    //    family. (Iterating from the first would echo whichever family was
    //    added earliest forever.)
    const groupCells = tmp.querySelectorAll('table tbody tr.apa-group td');
    if (groupCells.length){
      for (let i = groupCells.length - 1; i >= 0; i--){
        const text = (groupCells[i].textContent || '').trim();
        if (!text) continue;
        const cleaned = text
          .replace(/\s*[·•]\s*Ages?\s+[\w\d-]+(\s+\w+)?\s*$/i, '')
          .replace(/\s*[·•]\s*All\s+Ages\s*$/i, '')
          .trim();
        if (cleaned) return cleaned;
      }
    }
    // 2. Scan the full text for any known abbreviation as a fallback
    const allText = tmp.textContent || '';
    for (const fam of TEST_FAMILY_PATTERNS){
      if (allText.includes(fam)) return fam;
    }
    return null;
  }
  function pillLabelFor(html, sourceId){
    const parentId = (sourceId || '').split('::')[0];
    // Premorbid sources: don't pattern-match WAIS-IV/WMS-IV from the table -
    // those are predicted outcomes, not the test family.
    if (parentId && parentId.startsWith('pre-')){
      return SOURCE_LABELS[parentId] || null;
    }
    const family = detectTestFamily(html);
    if (family) return family;
    return SOURCE_LABELS[parentId] || null;
  }

  /* Walk the captured table for group separators (apa-group rows). If 2+
     groups are present, return a per-group array of { name, html } so each
     test family becomes its own table in the report. Returns [] when 0 or 1
     groups (in which case the caller treats it as a single merged table). */
  function extractGroupsFromHtml(html, sourceId){
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    const tbody = tmp.querySelector('table tbody');
    if (!tbody) return [];

    const sections = [];
    let current = null;
    [...tbody.children].forEach(row => {
      if (row.classList.contains('apa-group')){
        if (current) sections.push(current);
        current = {
          name: ((row.querySelector('td')?.textContent) || '').trim(),
          rows: []
        };
      } else if (current){
        current.rows.push(row);
      }
    });
    if (current) sections.push(current);

    if (sections.length < 2) return [];

    return sections.map(section => {
      const cloneTmp = document.createElement('div');
      cloneTmp.innerHTML = html;
      const cloneTbody = cloneTmp.querySelector('table tbody');
      if (!cloneTbody) return null;
      cloneTbody.innerHTML = '';
      section.rows.forEach(r => cloneTbody.appendChild(r.cloneNode(true)));
      // Build an intelligent title combining the parent method name with this
      // group's family name - e.g. "Crawford Regression-Based Change: WAIS-IV".
      const titleText = buildIntelligentTitle(sourceId, cloneTmp.innerHTML, section.name);
      const titleEl = cloneTmp.querySelector('.apa-table-title');
      if (titleEl) titleEl.textContent = titleText;
      return { name: section.name, title: titleText, html: cloneTmp.innerHTML };
    }).filter(Boolean);
  }

  /* When a parent source produces 2+ groups, replace any merged parent item
     with per-group items. Match existing per-group items by sourceId so that
     user customisations (column toggles, header overrides) survive.

     Pill firing logic:
     - Splits already present in state get silent updates.
     - When transitioning merged→split, the group that the merged item
       represented is treated as a continuation (no pill).
     - Only TRULY new groups fire pills.
  */
  function splitAndUpsert(parentSourceId, splits){
    // Identify the previously-merged item's group name (if any) so we can
    // treat that group as a continuation rather than a fresh addition.
    const merged = state.items.find(i => i.sourceId === parentSourceId);
    let prevGroupName = null;
    if (merged){
      const tmp = document.createElement('div');
      tmp.innerHTML = merged.html;
      const firstGroup = tmp.querySelector('table tbody tr.apa-group td');
      if (firstGroup) prevGroupName = (firstGroup.textContent || '').trim();
    }

    // Drop the merged parent item
    state.items = state.items.filter(i => i.sourceId !== parentSourceId);

    const validIds = new Set();
    splits.forEach(split => {
      const splitId = `${parentSourceId}::${split.name}`;
      validIds.add(splitId);
      const existing = state.items.find(i => i.sourceId === splitId);
      if (existing){
        if (existing.html !== split.html || existing.title !== split.title){
          existing.html = split.html;
          existing.title = split.title || split.name;
          existing.sourceTool = split.name;
          existing.updatedAt = new Date().toISOString();
          lastChangedItemId = existing.id;
        }
      } else {
        const newItem = {
          id: newId(),
          title: split.title || split.name,
          sourceTool: split.name,
          sourceId: splitId,
          html: split.html,
          addedAt: new Date().toISOString(),
          updatedAt: new Date().toISOString(),
          hiddenColumns: [],
          headerOverrides: []
        };
        state.items.push(newItem);
        lastChangedItemId = newItem.id;
        // Only pill for genuinely new content - not for the group that was
        // already represented by the merged item.
        const isContinuation = prevGroupName && split.name === prevGroupName;
        if (!isContinuation){
          showAddPrompt(split.name, splitId);
        }
      }
    });

    // Remove orphaned per-group items (groups no longer present)
    state.items = state.items.filter(i =>
      !i.sourceId.startsWith(`${parentSourceId}::`) || validIds.has(i.sourceId)
    );

    save();
    render();
  }

  /* ---------- core API ---------- */
  function addOrReplace({title, sourceTool, sourceId, html}){
    const now = new Date().toISOString();
    const existingIdx = state.items.findIndex(i => i.sourceId === sourceId);
    if (existingIdx >= 0){
      const existing = state.items[existingIdx];
      state.items[existingIdx] = {
        ...existing,
        title: title || existing.title,
        html,
        updatedAt: now
      };
    } else {
      state.items.push({
        id: newId(),
        title: (title || 'APA Table').trim().slice(0, 200),
        sourceTool: sourceTool || 'Unknown',
        sourceId,
        html,
        addedAt: now,
        updatedAt: now
      });
      ensureObserver(sourceId);
    }
    save();
    render();
    flashChip();
    if (typeof showToast === 'function'){
      const verb = existingIdx >= 0 ? '↻ Updated' : '✓ Added';
      showToast(`${verb} in working report`);
    }
  }
  function remove(id){
    state.items = state.items.filter(i => i.id !== id);
    save();
    render();
  }
  function clear(){
    if (!state.items.length) return;
    if (!window.confirm('Start a new report?\n\nThis clears all collected tables - the change cannot be undone.')) return;
    state.items = [];
    save();
    render();
    if (typeof showToast === 'function') showToast('New report started');
  }
  function clearSilent(){
    state.items = [];
    save();
    render();
  }
  /* Observer suppression - lets a global "clear all" wipe state without
     the auto-add MutationObserver immediately re-populating from the
     residual APA table HTML during the clear cascade. */
  let _suppressed = false;
  function setSuppressed(v){ _suppressed = !!v; }
  function isSuppressed(){ return _suppressed; }
  function moveItem(fromId, toId, dropAfter){
    const fromIdx = state.items.findIndex(i => i.id === fromId);
    let toIdx = state.items.findIndex(i => i.id === toId);
    if (fromIdx < 0 || toIdx < 0 || fromId === toId) return;
    const [moved] = state.items.splice(fromIdx, 1);
    if (toIdx > fromIdx) toIdx--;
    state.items.splice(toIdx + (dropAfter ? 1 : 0), 0, moved);
    save();
    render();
  }

  /* ---------- auto-add + auto-update via MutationObserver ----------
     Tables now flow into the working report automatically as you enter data.
     - First time a tool's APA container has a real table → auto-create item.
     - Subsequent edits → silently update the existing item.
     - No toast / chip-pulse on updates (would be spam during typing).
     - Subtle chip pulse only on first auto-add for a source.
     - When a container becomes empty, the bundle entry is PRESERVED - clearing
       inputs by accident shouldn't wipe your collected report.
  */
  function ensureObserver(sourceId){
    if (observers.has(sourceId)) return;
    const container = document.getElementById(sourceId);
    if (!container) return;
    // Seed the row count from the current DOM state so existing example rows
    // don't trigger a phantom "added" pill on first edit.
    if (lastRowCount[sourceId] == null){
      lastRowCount[sourceId] = countTableRows(container);
    }
    const handler = debounce(() => {
      // Bail out during a global clear so we don't immediately re-populate
      // the bundle from residual APA HTML during the clear cascade.
      if (_suppressed) return;
      // Always update the row-count tracker first, even if we're about to
      // early-return, so cleared tables don't leave stale counts behind.
      const currentRowCount = countTableRows(container);
      const rowsAdded = currentRowCount > (lastRowCount[sourceId] || 0);
      lastRowCount[sourceId] = currentRowCount;

      if (!container.querySelector('.apa-table')) return;
      if (typeof buildStandaloneHtml !== 'function') return;
      const rawHtml = buildStandaloneHtml(container);
      const html = applyMeaningfulTitle(rawHtml, sourceId);

      // AUTO-SPLIT: if the captured table has 2+ test families (group separators),
      // each becomes its own item in the working report.
      const splits = extractGroupsFromHtml(html, sourceId);
      if (splits.length >= 2){
        splitAndUpsert(sourceId, splits);
        return;
      }
      // Source is back to a single (or no) group - drop any per-group items
      // for this parent so the merged version takes over cleanly.
      const orphans = state.items.filter(i => i.sourceId.startsWith(`${sourceId}::`));
      if (orphans.length){
        state.items = state.items.filter(i => !i.sourceId.startsWith(`${sourceId}::`));
      }

      const item = state.items.find(i => i.sourceId === sourceId);
      const title = buildIntelligentTitle(sourceId, html);

      if (item){
        if (html === item.html && title === item.title) return;
        item.html = html;
        item.title = title;
        item.updatedAt = new Date().toISOString();
        lastChangedItemId = item.id;
        save();
        render();
        if (rowsAdded){
          showAddPrompt(pillLabelFor(html, sourceId), sourceId);
        }
      } else {
        // Auto-add (first time this tool produces data)
        const now = new Date().toISOString();
        const newItem = {
          id: newId(),
          title: (title || 'APA Table').trim().slice(0, 200),
          sourceTool: SOURCE_LABELS[sourceId] || sourceId,
          sourceId,
          html,
          addedAt: now,
          updatedAt: now,
          hiddenColumns: []
        };
        state.items.push(newItem);
        lastChangedItemId = newItem.id;
        save();
        render();
        // Show the confirmation pill - it'll fly into the chip after a hold
        // and pulse the chip on impact, so we don't need a separate flashChip().
        showAddPrompt(pillLabelFor(html, sourceId), sourceId);
      }
    }, 350);
    const obs = new MutationObserver(handler);
    obs.observe(container, { childList: true, subtree: true, characterData: true });
    observers.set(sourceId, obs);
  }
  function setupAllObservers(){
    SOURCE_IDS.forEach(id => ensureObserver(id));
  }

  /* ---------- exports ---------- */
  function buildReportHtmlBody(){
    if (!state.items.length) return '<p style="color:#888;font-style:italic;">No items in the working report yet.</p>';
    return state.items.map((it, i) => `<div style="page-break-inside:avoid;margin:0 0 28pt;">${effectiveItemHtml(it, i)}</div>`).join('');
  }
  async function copyAll(){
    if (!state.items.length){
      if (typeof showToast === 'function') showToast('Working report is empty', true);
      return;
    }
    const html = `<div style="font-family:'Times New Roman',serif;font-size:11pt;color:#000;">${buildReportHtmlBody()}</div>`;
    const plain = state.items.map((it, i) => {
      const tmp = document.createElement('div');
      tmp.innerHTML = effectiveItemHtml(it, i);
      return tmp.textContent.replace(/\s+/g, ' ').trim();
    }).join('\n\n---\n\n');
    try {
      if (navigator.clipboard && window.ClipboardItem){
        await navigator.clipboard.write([new ClipboardItem({
          'text/html':  new Blob([html],  { type:'text/html' }),
          'text/plain': new Blob([plain], { type:'text/plain' })
        })]);
      } else {
        await navigator.clipboard.writeText(plain);
      }
      if (typeof showToast === 'function') showToast(`✓ ${state.items.length} table${state.items.length===1?'':'s'} copied`);
      maybeShowKofiToast();
    } catch(e){
      console.error(e);
      if (typeof showToast === 'function') showToast('Copy failed - try selecting manually', true);
    }
  }
  function exportExcel(){
    if (!state.items.length){
      if (typeof showToast === 'function') showToast('Working report is empty', true);
      return;
    }
    // Generate CSV - Excel opens it natively without any "format mismatch"
    // warning, and the data round-trips cleanly. Multiple tables are stacked
    // with blank lines between them.
    const csvCell = s => {
      const v = String(s == null ? '' : s);
      return /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;
    };
    const sections = state.items.map((it, i) => {
      const tmp = document.createElement('div');
      tmp.innerHTML = effectiveItemHtml(it, i);
      const table = tmp.querySelector('table');
      if (!table) return '';
      const numText = (tmp.querySelector('.apa-table-num')?.textContent || `Table ${i + 1}`).trim();
      const titleText = (tmp.querySelector('.apa-table-title')?.textContent || '').trim();
      const heading = titleText ? `${numText}: ${titleText}` : numText;
      const lines = [csvCell(heading), ''];
      table.querySelectorAll('thead tr, tbody tr').forEach(tr => {
        const cells = [...tr.children].map(td =>
          csvCell((td.textContent || '').replace(/\s+/g, ' ').trim())
        );
        lines.push(cells.join(','));
      });
      const noteEl = tmp.querySelector('.apa-note, .apa-table-note, .apa-footer');
      if (noteEl){
        const note = (noteEl.textContent || '').replace(/\s+/g, ' ').trim();
        if (note){ lines.push(''); lines.push(csvCell(note)); }
      }
      return lines.join('\n');
    }).filter(Boolean).join('\n\n');

    // BOM + CRLF so Excel recognises UTF-8 and respects line endings on Windows.
    const csvBlob = new Blob(['﻿', sections.replace(/\n/g, '\r\n')], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(csvBlob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `assessment-report-${new Date().toISOString().slice(0,10)}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 100);
    if (typeof showToast === 'function') showToast('✓ CSV downloaded - opens in Excel');
    maybeShowKofiToast();
  }
  function exportWord(){
    if (!state.items.length){
      if (typeof showToast === 'function') showToast('Working report is empty', true);
      return;
    }
    const dateStr = new Date().toLocaleDateString(undefined, { year:'numeric', month:'long', day:'numeric' });
    const doc = `<html xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:w="urn:schemas-microsoft-com:office:word" xmlns="http://www.w3.org/TR/REC-html40">
<head>
<meta charset="utf-8">
<title>Assessment Report</title>
<!--[if gte mso 9]><xml><w:WordDocument><w:View>Print</w:View><w:Zoom>100</w:Zoom><w:DoNotOptimizeForBrowser/></w:WordDocument></xml><![endif]-->
<style>
  @page { margin: 1in; }
  body { font-family:'Times New Roman',serif; font-size:11pt; color:#000; line-height:1.4; }
  h1 { font-family:'Times New Roman',serif; font-size:14pt; font-weight:bold; margin:0 0 12pt; }
  .rb-doc-meta { color:#666; font-size:10pt; font-style:italic; margin:0 0 24pt; }
</style>
</head>
<body>
<h1>Assessment Report</h1>
<p class="rb-doc-meta">Compiled ${escapeHtmlLocal(dateStr)} · ${state.items.length} table${state.items.length===1?'':'s'}</p>
${buildReportHtmlBody()}
</body>
</html>`;
    const blob = new Blob(['﻿', doc], { type: 'application/msword' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `assessment-report-${new Date().toISOString().slice(0,10)}.doc`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 100);
    if (typeof showToast === 'function') showToast('✓ Word document downloaded');
    maybeShowKofiToast();
  }

  /* ---------- UI rendering ---------- */
  function injectUI(){
    if (document.getElementById('report-bundle-root')) return;
    const root = document.createElement('div');
    root.className = 'rb-root';
    root.id = 'report-bundle-root';
    root.innerHTML = `
      <button class="rb-chip" data-rb-action="toggle" type="button" aria-label="Toggle working report">
        <span class="rb-chip-icon" aria-hidden="true">
          <svg viewBox="0 0 18 18" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">
            <rect x="3.25" y="2.25" width="9" height="11" rx="1.25"/>
            <rect x="5.75" y="4.75" width="9" height="11" rx="1.25" fill="currentColor" fill-opacity="0.18"/>
            <line x1="8" y1="8" x2="13" y2="8"/>
            <line x1="8" y1="11" x2="11.5" y2="11"/>
          </svg>
        </span>
        <span class="rb-chip-label">Working Report <span class="rb-chip-sub">APA Tables</span></span>
        <span class="rb-chip-count" data-rb-count>0</span>
      </button>
      <div class="rb-onboarding" data-rb-onboarding hidden aria-hidden="true">
        <div class="rb-onboarding-bubble">
          <div class="rb-onboarding-title">Your Working Report</div>
          <div class="rb-onboarding-body">
            <p>Every APA-formatted table you generate across the suite is <strong>auto-saved into one place</strong> as you work.</p>
            <p>Open the panel to <strong>reorder, hide columns, edit titles, and export</strong> the whole bundle to Word, Excel, or your clipboard - paste straight into your report.</p>
            <p>Everything lives in your browser. <strong>Nothing is uploaded, nothing leaves your device.</strong></p>
          </div>
          <div class="rb-onboarding-cta">Click the button below to open it ↓</div>
          <button class="rb-onboarding-dismiss" data-rb-action="dismiss-onboarding" type="button" aria-label="Dismiss">
            <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" aria-hidden="true"><path d="M3 3l6 6M9 3l-6 6"/></svg>
          </button>
          <div class="rb-onboarding-tail" aria-hidden="true"></div>
        </div>
      </div>
      <div class="rb-drawer" hidden>
        <div class="rb-drawer-head">
          <div class="rb-drawer-title">
            <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
              <path d="M3.5 1.5h7l3 3v10h-10z"/>
              <path d="M10.5 1.5v3h3"/>
              <line x1="5.5" y1="8" x2="11.5" y2="8"/>
              <line x1="5.5" y1="11" x2="9.5" y2="11"/>
            </svg>
            <span>Working Report</span>
            <span class="rb-drawer-count" data-rb-count-text>0 items</span>
          </div>
          <div class="rb-drawer-head-actions">
            <div class="rb-mode-toggle" role="tablist" aria-label="View mode">
              <button class="rb-mode-btn is-active" data-rb-action="mode-edit" type="button" role="tab" aria-selected="true">
                <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                  <path d="M2 11.5h10"/><path d="M2 7.5h10"/><path d="M2 3.5h10"/>
                </svg>
                Edit
              </button>
              <button class="rb-mode-btn" data-rb-action="mode-preview" type="button" role="tab" aria-selected="false">
                <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                  <rect x="3" y="2" width="8" height="10" rx="0.5"/>
                  <line x1="5" y1="5" x2="9" y2="5"/>
                  <line x1="5" y1="7.5" x2="9" y2="7.5"/>
                  <line x1="5" y1="10" x2="7.5" y2="10"/>
                </svg>
                Preview
              </button>
            </div>
            <button class="rb-head-btn" data-rb-action="maximise" type="button" aria-label="Maximise" title="Maximise">
              <svg class="rb-icon-maximise" viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                <rect x="2" y="2" width="10" height="10" rx="1"/>
              </svg>
              <svg class="rb-icon-restore" viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                <rect x="4" y="4" width="8" height="8" rx="1"/>
                <path d="M2 9V3a1 1 0 0 1 1-1h6"/>
              </svg>
            </button>
            <button class="rb-head-btn" data-rb-action="minimise" type="button" aria-label="Minimise" title="Minimise">
              <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" aria-hidden="true">
                <path d="M3 11h8"/>
              </svg>
            </button>
          </div>
        </div>
        <div class="rb-drawer-body" data-rb-body></div>
        <div class="rb-drawer-actions">
          <button class="btn btn-clear rb-action-clear" data-rb-action="clear" type="button">
            <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M2 4h10"/><path d="M5.5 4V2.5h3V4"/><path d="M3 4l1 8h6l1-8"/></svg>
            New report
          </button>
          <div class="rb-actions-right">
            <button class="btn rb-action-copy" data-rb-action="copy" type="button">
              <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="4" y="4" width="8" height="8" rx="1"/><path d="M2 9V3a1 1 0 0 1 1-1h6"/></svg>
              Copy all tables
            </button>
            <button class="btn rb-action-export-excel" data-rb-action="export-excel" type="button">
              <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="2" y="2.5" width="10" height="9" rx="1"/><path d="M2 6h10"/><path d="M5.5 2.5v9"/></svg>
              Export to Excel
            </button>
            <button class="btn rb-action-export" data-rb-action="export-word" type="button">
              <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M12 9v2.5a1 1 0 0 1-1 1H3a1 1 0 0 1-1-1V9"/><path d="M4 6l3 3 3-3"/><path d="M7 9V1.5"/></svg>
              Export to Word
            </button>
          </div>
        </div>
      </div>
    `;
    document.body.appendChild(root);
    rootEl = root;
  }

  function bindEvents(){
    if (!rootEl) return;
    // Delegated click handler at the DOCUMENT level
    document.addEventListener('click', e => {
      // Per-item Options menu close-on-outside-click
      const insideMenu = e.target.closest('.rb-item-menu');
      const onMenuTrigger = e.target.closest('[data-rb-item-options]');
      if (!insideMenu && !onMenuTrigger){
        rootEl?.querySelectorAll('.rb-item-menu.is-open').forEach(m => m.classList.remove('is-open'));
      }

      if (!e.target.closest('#report-bundle-root')) return;
      const actionBtn = e.target.closest('[data-rb-action]');
      if (actionBtn){
        e.preventDefault();
        e.stopPropagation();
        const action = actionBtn.dataset.rbAction;
        if (action === 'toggle')        toggle();
        else if (action === 'minimise') close();
        else if (action === 'maximise') toggleMaximised();
        else if (action === 'clear')    {
          // Forward to the topbar's "Clear all tables" button so it does
          // the full session wipe (every tool + the bundle), with a single
          // confirm. Falls back to local-only clear if the topbar button
          // isn't available for any reason.
          const topbarBtn = document.getElementById('topbar-clear-all');
          if (topbarBtn) topbarBtn.click();
          else clear();
        }
        else if (action === 'copy')     copyAll();
        else if (action === 'export-word') exportWord();
        else if (action === 'export-excel') exportExcel();
        else if (action === 'dismiss-onboarding') dismissOnboarding();
        else if (action === 'mode-edit')    setPreviewMode(false);
        else if (action === 'mode-preview') setPreviewMode(true);
        return;
      }
      // Per-item: Copy single table
      const copyOne = e.target.closest('[data-rb-item-copy]');
      if (copyOne){
        e.preventDefault();
        e.stopPropagation();
        copyItem(copyOne.dataset.rbItemCopy);
        return;
      }
      // Per-item: Options menu trigger
      const optsTrigger = e.target.closest('[data-rb-item-options]');
      if (optsTrigger){
        e.preventDefault();
        e.stopPropagation();
        const id = optsTrigger.dataset.rbItemOptions;
        const menu = rootEl.querySelector(`.rb-item-menu[data-rb-item-menu="${id}"]`);
        rootEl.querySelectorAll('.rb-item-menu.is-open').forEach(m => {
          if (m !== menu) m.classList.remove('is-open');
        });
        menu?.classList.toggle('is-open');
        return;
      }
      // Per-item: column toggles inside the Options menu
      const colToggle = e.target.closest('[data-rb-col-toggle]');
      if (colToggle){
        e.preventDefault();
        e.stopPropagation();
        const id = colToggle.dataset.rbItemId;
        const colIdx = colToggle.dataset.rbColToggle;
        toggleColumn(id, colIdx);
        // Re-open the menu after render so the user can toggle multiple columns
        requestAnimationFrame(() => {
          const menu = rootEl?.querySelector(`.rb-item-menu[data-rb-item-menu="${id}"]`);
          if (menu) menu.classList.add('is-open');
        });
        return;
      }
      // Per-item: Options menu items
      const menuItem = e.target.closest('[data-rb-item-action]');
      if (menuItem){
        e.preventDefault();
        e.stopPropagation();
        const action = menuItem.dataset.rbItemAction;
        const id = menuItem.dataset.rbItemId;
        if (action === 'up')          moveUp(id);
        else if (action === 'down')   moveDown(id);
        else if (action === 'top')    moveToTop(id);
        else if (action === 'bottom') moveToBottom(id);
        else if (action === 'refresh') refreshItem(id);
        else if (action === 'remove') remove(id);
        rootEl.querySelectorAll('.rb-item-menu.is-open').forEach(m => m.classList.remove('is-open'));
        return;
      }
      const removeBtn = e.target.closest('[data-rb-remove]');
      if (removeBtn){
        e.preventDefault();
        e.stopPropagation();
        remove(removeBtn.dataset.rbRemove);
        return;
      }
    });
    // Drag-and-drop reorder
    const body = rootEl.querySelector('[data-rb-body]');
    // Helper: lock the drawer's current height so it doesn't shrink as items
    // collapse to header strips during reorder.
    function lockDrawerHeight(){
      const drawerEl = rootEl.querySelector('.rb-drawer');
      if (!drawerEl) return;
      drawerEl.style.height = drawerEl.offsetHeight + 'px';
    }
    function unlockDrawerHeight(){
      const drawerEl = rootEl.querySelector('.rb-drawer');
      if (!drawerEl) return;
      drawerEl.style.height = '';
    }

    body.addEventListener('dragstart', e => {
      // `draggable` lives on the grip handle only, so dragstart only ever
      // fires when the user grabbed the grip. No guard needed.
      const grip = e.target.closest('.rb-item-grip');
      if (!grip) return;
      const item = grip.closest('.rb-item');
      if (!item) return;
      dragId = item.dataset.rbId;
      item.classList.add('is-dragging');
      lockDrawerHeight();
      body.classList.add('is-reorder-mode');
      e.dataTransfer.effectAllowed = 'move';
      try { e.dataTransfer.setData('text/plain', dragId); } catch(_){}
      // Show the entire item as the drag preview, not the small grip handle.
      try {
        const rect = item.getBoundingClientRect();
        e.dataTransfer.setDragImage(item, e.clientX - rect.left, e.clientY - rect.top);
      } catch(_){}
    });
    body.addEventListener('dragend', e => {
      body.querySelectorAll('.rb-item.is-dragging').forEach(el => el.classList.remove('is-dragging'));
      body.classList.remove('is-reorder-mode');
      unlockDrawerHeight();
      body.querySelectorAll('.is-drop-before, .is-drop-after').forEach(el => el.classList.remove('is-drop-before','is-drop-after'));
      dragId = null;
    });
    body.addEventListener('dragover', e => {
      const target = e.target.closest('.rb-item');
      if (!target || !dragId || target.dataset.rbId === dragId) return;
      e.preventDefault();
      const rect = target.getBoundingClientRect();
      const after = (e.clientY - rect.top) > rect.height / 2;
      body.querySelectorAll('.is-drop-before, .is-drop-after').forEach(el => el.classList.remove('is-drop-before','is-drop-after'));
      target.classList.add(after ? 'is-drop-after' : 'is-drop-before');
      e.dataTransfer.dropEffect = 'move';
    });
    body.addEventListener('drop', e => {
      body.classList.remove('is-reorder-mode');
      unlockDrawerHeight();
      body.querySelectorAll('.is-drop-before, .is-drop-after').forEach(el =>
        el.classList.remove('is-drop-before', 'is-drop-after')
      );
      const target = e.target.closest('.rb-item');
      if (!target || !dragId || target.dataset.rbId === dragId){ dragId = null; return; }
      e.preventDefault();
      const rect = target.getBoundingClientRect();
      const after = (e.clientY - rect.top) > rect.height / 2;
      const fromId = dragId;
      dragId = null;
      moveItem(fromId, target.dataset.rbId, after);
    });

    // Keyboard reorder - focus a grip handle, then Alt+ArrowUp/Down to
    // shift that item one position. Standard accessibility pattern.
    body.addEventListener('keydown', e => {
      if (!(e.altKey && (e.key === 'ArrowUp' || e.key === 'ArrowDown'))) return;
      const grip = e.target.closest('.rb-item-grip');
      if (!grip) return;
      const item = grip.closest('.rb-item');
      if (!item) return;
      e.preventDefault();
      const id = item.dataset.rbId;
      const idx = state.items.findIndex(x => x.id === id);
      if (idx < 0) return;
      const dir = e.key === 'ArrowUp' ? -1 : 1;
      const target = idx + dir;
      if (target < 0 || target >= state.items.length) return;
      const [moved] = state.items.splice(idx, 1);
      state.items.splice(target, 0, moved);
      lastChangedItemId = id;
      save();
      render();
      // Restore focus to the grip on the moved item so chained Alt+Arrow works
      requestAnimationFrame(() => {
        const newGrip = rootEl?.querySelector(`.rb-item[data-rb-id="${CSS.escape(id)}"] .rb-item-grip`);
        newGrip?.focus();
      });
    });
    // ESC closes (and blurs any in-progress header edit first)
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && !state.minimized){
        // If editing a header, just blur - don't close the drawer
        const editing = document.activeElement?.closest('.rb-editable-header');
        if (editing){ editing.blur(); return; }
        close();
      }
      // Enter on a header cell → commit (blur to save)
      if (e.key === 'Enter'){
        const cell = e.target.closest?.('.rb-editable-header');
        if (cell){ e.preventDefault(); cell.blur(); }
      }
    });
    // Save header override on blur (focusout bubbles, blur doesn't)
    document.addEventListener('focusout', e => {
      const cell = e.target.closest?.('.rb-editable-header');
      if (!cell) return;
      const itemId = cell.dataset.rbItemId;
      const colIdx = parseInt(cell.dataset.rbColIdx, 10);
      saveHeaderOverride(itemId, colIdx, cell.textContent);
    });
    // Click anywhere outside the bundle root → minimize. Uses composedPath()
    // which is captured at dispatch time, so it sees the original ancestor
    // chain even if the inner click handler removed the target from DOM
    // (e.g. clicking ✕ to remove an item triggers a full re-render).
    document.addEventListener('click', e => {
      if (state.minimized) return;
      const path = (typeof e.composedPath === 'function') ? e.composedPath() : [];
      for (const el of path){
        if (el && el.id === 'report-bundle-root') return;
      }
      // Fallback for browsers without composedPath, when target is still in DOM
      if (e.target.closest && e.target.closest('#report-bundle-root')) return;
      close();
    });

    // Resize handle removed - drawer is now fixed-size.
  }

  /* ---------- state transitions ---------- */
  function open(){
    state.minimized = false;
    state.onboardingSeen = true; // first open dismisses the hint
    hideAddPrompt();             // dismiss the transient "View live report" prompt
    save();
    render();
  }
  function close(){
    // Play the bubble-out animation (drawer shrinks back into the chip), then
    // commit the state change once the animation finishes.
    if (rootEl && !state.minimized){
      rootEl.classList.add('is-closing');
      setTimeout(() => {
        state.minimized = true;
        save();
        render(); // removes is-open + hides drawer
        rootEl?.classList.remove('is-closing'); // cleanup last so the drawer never re-animates
      }, 240); // matches rb-drawer-bubble-out duration
    } else {
      state.minimized = true;
      save();
      render();
    }
  }
  function toggle(){ state.minimized ? open() : close(); }
  function toggleMaximised(){
    state.maximised = !state.maximised;
    save();
    render();
  }
  /* Preview mode - toggles the drawer body into a Word-style print
     preview where APA tables sit on a white "page" with proper margins,
     editing chrome (drag, remove, options) hidden. Session-only; not
     persisted, since it's a transient view, not a setting. */
  let previewMode = false;
  function setPreviewMode(on){
    previewMode = !!on;
    if (!rootEl) return;
    const body = rootEl.querySelector('[data-rb-body]');
    if (body) body.classList.toggle('is-preview-mode', previewMode);
    rootEl.querySelectorAll('.rb-mode-btn').forEach(btn => {
      const isPreviewBtn = btn.dataset.rbAction === 'mode-preview';
      const shouldBeActive = isPreviewBtn ? previewMode : !previewMode;
      btn.classList.toggle('is-active', shouldBeActive);
      btn.setAttribute('aria-selected', String(shouldBeActive));
    });
  }
  function flashChip(){
    if (!rootEl) return;
    rootEl.classList.add('rb-flash');
    setTimeout(() => rootEl?.classList.remove('rb-flash'), 600);
  }

  /* "Added to report" pills - one per source, stackable. Each pill pops in,
     bobs gently, then flies into the chip. Multiple sources => stacked above
     each other with the newest closest to the chip. Same source updating
     => existing pill stays and its hold timer resets (no duplicate pile-up
     during continuous typing). */
  const PILL_HEIGHT = 40;       // approx pill height
  const PILL_GAP    = 8;        // gap between stacked pills
  const PILL_BASE_BOTTOM = 82;  // distance from viewport bottom for pill #1
  const CHIP_CENTER_BOTTOM = 24; // chip's vertical centre, target for fly animation

  function buildPillNode(sourceLabel, sourceId){
    const text = sourceLabel ? `${sourceLabel} added to report` : 'Added to report';
    const pill = document.createElement('button');
    pill.className = 'rb-add-prompt';
    pill.type = 'button';
    pill.setAttribute('data-rb-action', 'toggle');
    pill.setAttribute('aria-label', text);
    pill.setAttribute('title', text);
    if (sourceId) pill.setAttribute('data-rb-source', sourceId);
    pill.innerHTML = `
      <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.9" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
        <path d="M3 7.5l3 3 5-7"/>
      </svg>
      <span data-rb-add-prompt-label>${escapeHtmlLocal(text)}</span>
    `;
    return pill;
  }

  function startPillFly(pill){
    if (!pill || !pill.parentNode) return;
    clearTimeout(pill._holdTimer);
    pill.classList.remove('is-visible');
    // Compute fly distance from this pill's CURRENT bottom to the chip's centre
    const bottomVal = parseInt(pill.style.bottom, 10) || PILL_BASE_BOTTOM;
    const flyDistance = bottomVal - CHIP_CENTER_BOTTOM;
    pill.style.setProperty('--rb-fly-distance', flyDistance + 'px');
    pill.classList.add('is-flying');
    // Chip "catches" the pill just before arrival
    setTimeout(() => flashChip(), 380);
    pill._doneTimer = setTimeout(() => {
      pill.remove();
      reflowPillStack();
    }, 540);
  }

  /* Make column headers in the drawer click-to-edit. Called after each render
     so newly-rendered <th> cells get the contenteditable affordance. */
  function decorateEditableHeaders(){
    if (!rootEl) return;
    rootEl.querySelectorAll('.rb-item').forEach(article => {
      const itemId = article.dataset.rbId;
      if (!itemId) return;
      const headerRows = article.querySelectorAll('.rb-item-rendered table thead tr');
      if (!headerRows.length) return;
      const lastRow = headerRows[headerRows.length - 1];
      [...lastRow.children].forEach((cell, idx) => {
        cell.setAttribute('contenteditable', 'plaintext-only');
        cell.classList.add('rb-editable-header');
        cell.dataset.rbItemId = itemId;
        cell.dataset.rbColIdx = String(idx);
        cell.title = 'Click to rename column';
      });
    });
  }
  function saveHeaderOverride(itemId, colIdx, newText){
    const item = state.items.find(i => i.id === itemId);
    if (!item) return;
    if (!Array.isArray(item.headerOverrides)) item.headerOverrides = [];

    // Compare against the original (pre-override) header text - if it matches,
    // clear the override so the column reverts to whatever the source produces.
    const tmp = document.createElement('div');
    tmp.innerHTML = item.html;
    const headerRows = tmp.querySelectorAll('table thead tr');
    const lastRow = headerRows[headerRows.length - 1];
    const originalText = ((lastRow?.children[colIdx]?.textContent) || '').trim();

    const trimmed = String(newText || '').trim();
    if (trimmed === '' || trimmed === originalText){
      item.headerOverrides[colIdx] = null;
    } else {
      item.headerOverrides[colIdx] = trimmed;
    }
    item.updatedAt = new Date().toISOString();
    save();
    render(); // by the time blur fires, focus has already moved off the cell
  }

  function reflowPillStack(){
    if (!rootEl) return;
    const pills = rootEl.querySelectorAll('.rb-add-prompt:not(.is-flying)');
    // Newest pill (last added to DOM) sits closest to the chip; stack upward
    [...pills].reverse().forEach((pill, i) => {
      pill.style.bottom = (PILL_BASE_BOTTOM + i * (PILL_HEIGHT + PILL_GAP)) + 'px';
    });
  }

  function showAddPrompt(sourceLabel, sourceId){
    if (!state.minimized) return;
    if (!rootEl) return;

    // Dedupe by sourceId - if a non-flying pill for this source already exists,
    // refresh its label and reset its hold timer instead of stacking a new one.
    if (sourceId){
      const existing = rootEl.querySelector(
        `.rb-add-prompt[data-rb-source="${CSS.escape(sourceId)}"]:not(.is-flying)`
      );
      if (existing){
        const text = sourceLabel ? `${sourceLabel} added to report` : 'Added to report';
        const labelEl = existing.querySelector('[data-rb-add-prompt-label]');
        if (labelEl) labelEl.textContent = text;
        existing.setAttribute('title', text);
        existing.setAttribute('aria-label', text);
        clearTimeout(existing._holdTimer);
        existing._holdTimer = setTimeout(() => startPillFly(existing), 1500);
        return;
      }
    }

    // Fresh pill - append, position, animate
    const pill = buildPillNode(sourceLabel, sourceId);
    rootEl.appendChild(pill);
    reflowPillStack(); // sets bottom for the new (newest) pill at the base position

    void pill.offsetWidth;             // reflow so the entrance animation runs
    pill.classList.add('is-visible');

    pill._holdTimer = setTimeout(() => startPillFly(pill), 1500);
  }

  function hideAddPrompt(){
    if (!rootEl) return;
    rootEl.querySelectorAll('.rb-add-prompt').forEach(pill => {
      clearTimeout(pill._holdTimer);
      clearTimeout(pill._doneTimer);
      pill.remove();
    });
  }
  function dismissOnboarding(){
    state.onboardingSeen = true;
    save();
    render();
  }

  /* ---------- per-item ops ---------- */
  function moveToTop(id){
    const idx = state.items.findIndex(i => i.id === id);
    if (idx <= 0) return;
    const [moved] = state.items.splice(idx, 1);
    state.items.unshift(moved);
    save();
    render();
  }
  function moveToBottom(id){
    const idx = state.items.findIndex(i => i.id === id);
    if (idx < 0 || idx === state.items.length - 1) return;
    const [moved] = state.items.splice(idx, 1);
    state.items.push(moved);
    save();
    render();
  }
  function moveUp(id){
    const idx = state.items.findIndex(i => i.id === id);
    if (idx <= 0) return;
    const tmp = state.items[idx];
    state.items[idx] = state.items[idx - 1];
    state.items[idx - 1] = tmp;
    save();
    render();
  }
  function moveDown(id){
    const idx = state.items.findIndex(i => i.id === id);
    if (idx < 0 || idx >= state.items.length - 1) return;
    const tmp = state.items[idx];
    state.items[idx] = state.items[idx + 1];
    state.items[idx + 1] = tmp;
    save();
    render();
  }
  function refreshItem(id){
    const item = state.items.find(i => i.id === id);
    if (!item) return;
    // Split items have sourceId like "bat-apa::CVLT-3 Indices" - only the
    // prefix is a real DOM id.
    const parentId = item.sourceId.split('::')[0];
    const container = document.getElementById(parentId);
    if (!container || !container.querySelector('.apa-table')){
      if (typeof showToast === 'function') showToast('Source table is empty - nothing to refresh', true);
      return;
    }
    if (typeof buildStandaloneHtml !== 'function') return;
    const fullHtml = applyMeaningfulTitle(buildStandaloneHtml(container), parentId);

    if (item.sourceId.includes('::')){
      // Re-extract just this group from the latest source HTML
      const splits = extractGroupsFromHtml(fullHtml, parentId);
      const groupName = item.sourceId.split('::')[1];
      const split = splits.find(s => s.name === groupName);
      if (!split){
        if (typeof showToast === 'function') showToast(`"${groupName}" no longer in source`, true);
        return;
      }
      item.html = split.html;
      item.title = split.title || split.name;
      item.sourceTool = split.name;
    } else {
      item.html = fullHtml;
      item.title = buildIntelligentTitle(item.sourceId, fullHtml);
    }
    item.updatedAt = new Date().toISOString();
    lastChangedItemId = id;
    save();
    render();
    if (typeof showToast === 'function') showToast('↻ Refreshed from source');
  }
  function toggleColumn(itemId, colIdx){
    const item = state.items.find(i => i.id === itemId);
    if (!item) return;
    const idx = Number(colIdx);
    if (!Number.isFinite(idx)) return;
    if (!Array.isArray(item.hiddenColumns)) item.hiddenColumns = [];
    const i = item.hiddenColumns.indexOf(idx);
    if (i >= 0) item.hiddenColumns.splice(i, 1);
    else        item.hiddenColumns.push(idx);
    save();
    render();
  }
  async function copyItem(id){
    const idx = state.items.findIndex(i => i.id === id);
    if (idx < 0) return;
    const item = state.items[idx];
    const itemHtml = effectiveItemHtml(item, idx);
    const html = `<div style="font-family:'Times New Roman',serif;font-size:11pt;color:#000;">${itemHtml}</div>`;
    const tmp = document.createElement('div');
    tmp.innerHTML = itemHtml;
    const plain = tmp.textContent.replace(/\s+/g, ' ').trim();
    try {
      if (navigator.clipboard && window.ClipboardItem){
        await navigator.clipboard.write([new ClipboardItem({
          'text/html':  new Blob([html],  { type:'text/html' }),
          'text/plain': new Blob([plain], { type:'text/plain' })
        })]);
      } else {
        await navigator.clipboard.writeText(plain);
      }
      if (typeof showToast === 'function') showToast('✓ Table copied to clipboard');
      maybeShowKofiToast();
    } catch(e){
      console.error(e);
      if (typeof showToast === 'function') showToast('Copy failed - try selecting manually', true);
    }
  }

  /* ---------- rendering ---------- */
  function renderItems(){
    const body = rootEl?.querySelector('[data-rb-body]');
    if (!body) return;
    if (!state.items.length){
      body.innerHTML = `
        <div class="rb-empty">
          <svg viewBox="0 0 32 32" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
            <path d="M8 4h12l4 4v20H8z"/>
            <path d="M20 4v4h4"/>
            <line x1="12" y1="14" x2="22" y2="14"/>
            <line x1="12" y1="19" x2="20" y2="19"/>
            <line x1="12" y1="24" x2="18" y2="24"/>
          </svg>
          <div class="rb-empty-title">No tables yet</div>
          <div class="rb-empty-sub">Tables appear here automatically as you enter scores in any tool. Drag to reorder, click <strong>↺ New report</strong> to start fresh.</div>
        </div>`;
      return;
    }
    body.innerHTML = state.items.map((it, i) => {
      const cols = getItemColumns(it);
      const hidden = new Set(it.hiddenColumns || []);
      const colsHtml = cols.length ? `
        <div class="rb-item-menu-section">Columns</div>
        ${cols.map(c => `
          <button class="rb-item-menu-item rb-col-toggle" type="button" role="menuitemcheckbox"
            aria-checked="${!hidden.has(c.idx)}"
            data-rb-col-toggle="${c.idx}"
            data-rb-item-id="${escapeHtmlLocal(it.id)}">
            <span class="rb-col-check ${hidden.has(c.idx) ? '' : 'is-checked'}" aria-hidden="true">
              <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6.5l2 2 4-5"/></svg>
            </span>
            <span class="rb-col-label">${escapeHtmlLocal(c.label)}</span>
          </button>
        `).join('')}
        <div class="rb-item-menu-sep"></div>
      ` : '';
      return `
      <article class="rb-item" data-rb-id="${escapeHtmlLocal(it.id)}">
        <header class="rb-item-header">
          <span class="rb-item-grip" draggable="true" title="Drag to reorder (or focus + Alt+↑/↓)" role="button" tabindex="0" aria-label="Drag to reorder. Keyboard: Alt + Arrow Up or Arrow Down.">
            <svg class="rb-item-grip-icon" viewBox="0 0 14 18" fill="currentColor" aria-hidden="true"><circle cx="5" cy="4" r="1.4"/><circle cx="9" cy="4" r="1.4"/><circle cx="5" cy="9" r="1.4"/><circle cx="9" cy="9" r="1.4"/><circle cx="5" cy="14" r="1.4"/><circle cx="9" cy="14" r="1.4"/></svg>
          </span>
          <span class="rb-item-num">${i + 1}</span>
          <div class="rb-item-meta">
            <span class="rb-item-source">${escapeHtmlLocal(it.sourceTool)}</span>
            <span class="rb-item-time" data-rb-time="${escapeHtmlLocal(it.updatedAt || it.addedAt)}">${formatRelative(it.updatedAt || it.addedAt)}</span>
          </div>
          <div class="rb-item-actions">
            <button class="rb-item-actionbtn rb-item-copy" type="button" data-rb-item-copy="${escapeHtmlLocal(it.id)}" aria-label="Copy this table" title="Copy this table to clipboard">
              <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="4" y="4" width="8" height="8" rx="1"/><path d="M2 9V3a1 1 0 0 1 1-1h6"/></svg>
              <span>Copy</span>
            </button>
            <div class="rb-item-options-wrap">
              <button class="rb-item-actionbtn rb-item-options" type="button" data-rb-item-options="${escapeHtmlLocal(it.id)}" aria-label="More options" title="More options">
                <svg viewBox="0 0 14 14" fill="currentColor" aria-hidden="true"><circle cx="3" cy="7" r="1.3"/><circle cx="7" cy="7" r="1.3"/><circle cx="11" cy="7" r="1.3"/></svg>
                <span>Options</span>
              </button>
              <div class="rb-item-menu" data-rb-item-menu="${escapeHtmlLocal(it.id)}" role="menu">
                ${colsHtml}
                <button class="rb-item-menu-item" data-rb-item-action="up" data-rb-item-id="${escapeHtmlLocal(it.id)}" type="button" role="menuitem"${i === 0 ? ' disabled' : ''}>
                  <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M3 6l4-4 4 4"/><path d="M7 2v10"/></svg>
                  Move up
                </button>
                <button class="rb-item-menu-item" data-rb-item-action="down" data-rb-item-id="${escapeHtmlLocal(it.id)}" type="button" role="menuitem"${i === state.items.length - 1 ? ' disabled' : ''}>
                  <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M7 2v10"/><path d="M3 8l4 4 4-4"/></svg>
                  Move down
                </button>
                <button class="rb-item-menu-item" data-rb-item-action="refresh" data-rb-item-id="${escapeHtmlLocal(it.id)}" type="button" role="menuitem">
                  <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M2 7a5 5 0 0 1 9-3"/><path d="M11 1v3.5h-3"/><path d="M12 7a5 5 0 0 1-9 3"/><path d="M3 13V9.5h3"/></svg>
                  Refresh from source
                </button>
                <div class="rb-item-menu-sep"></div>
                <button class="rb-item-menu-item rb-item-menu-danger" data-rb-item-action="remove" data-rb-item-id="${escapeHtmlLocal(it.id)}" type="button" role="menuitem">
                  <svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M2 4h10"/><path d="M5.5 4V2.5h3V4"/><path d="M3 4l1 8h6l1-8"/></svg>
                  Remove
                </button>
              </div>
            </div>
            <button class="rb-item-remove" type="button" data-rb-remove="${escapeHtmlLocal(it.id)}" aria-label="Remove this table from the report" title="Remove from report">
              <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" aria-hidden="true"><path d="M3 3l6 6M9 3l-6 6"/></svg>
            </button>
          </div>
        </header>
        <div class="rb-item-rendered">${effectiveItemHtml(it, i)}</div>
      </article>
    `;}).join('');
  }

  function render(){
    if (!rootEl) return;
    const drawer = rootEl.querySelector('.rb-drawer');
    const onboarding = rootEl.querySelector('[data-rb-onboarding]');
    rootEl.querySelectorAll('[data-rb-count]').forEach(el => el.textContent = String(state.items.length));
    rootEl.querySelectorAll('[data-rb-count-text]').forEach(el => el.textContent = `${state.items.length} item${state.items.length === 1 ? '' : 's'}`);
    rootEl.querySelectorAll('[data-rb-action="clear"], [data-rb-action="copy"], [data-rb-action="export-word"], [data-rb-action="export-excel"]').forEach(b => b.disabled = !state.items.length);
    rootEl.dataset.state = state.minimized ? 'closed' : 'open';
    rootEl.classList.toggle('is-open', !state.minimized);
    rootEl.classList.toggle('is-maximised', state.maximised && !state.minimized);

    // Chip is ALWAYS visible (only its position shifts via CSS based on data-state)
    // Drawer is only shown when open
    if (drawer) drawer.hidden = state.minimized;

    // Drawer is now a fixed-size floating popover - no inline width override,
    // no body padding push.
    if (drawer) drawer.style.width = '';
    document.body.style.paddingRight = '';
    document.body.classList.remove('rb-docked-active');

    // Onboarding hint: show only on first ever visit, when drawer is closed
    if (onboarding){
      const showHint = !state.onboardingSeen && state.minimized;
      onboarding.hidden = !showHint;
      onboarding.setAttribute('aria-hidden', String(!showHint));
    }

    renderItems();
    decorateEditableHeaders();

    // Auto-scroll to the most recently changed item - only if drawer is open
    if (lastChangedItemId && !state.minimized){
      const target = lastChangedItemId;
      lastChangedItemId = null;
      requestAnimationFrame(() => {
        const el = rootEl?.querySelector(`.rb-item[data-rb-id="${target}"]`);
        el?.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      });
    }
  }

  /* "+ Add to report" buttons removed - adds are automatic via MutationObserver.
     If a saved bundle from before still references these IDs, they're harmless. */

  /* ---------- init ---------- */
  function init(){
    load();
    injectUI();
    bindEvents();
    setupAllObservers();
    render();
    setInterval(() => {
      if (state.minimized) return;
      rootEl?.querySelectorAll('[data-rb-time]').forEach(el => {
        el.textContent = formatRelative(el.dataset.rbTime);
      });
    }, 30000);

    // Time-based Ko-fi prompt - fires after 5 min on site IF the user has
    // never seen the modal (across any tab/session) and no export has
    // already triggered it this page load. Persistence is via localStorage
    // inside maybeShowKofiToast itself.
    if (!hasSeenKofi()){
      setTimeout(() => {
        if (!kofiToastShown && !hasSeenKofi()) maybeShowKofiToast();
      }, 5 * 60 * 1000);
    }

    // Header "Buy me a coffee" button - opens the modal in-page rather than
    // letting the anchor navigate to a new tab. Force-opens regardless of
    // whether the modal has been shown before (the seen flag only suppresses
    // automatic triggers, not deliberate user clicks).
    document.querySelectorAll('.topbar-kofi').forEach(el => {
      el.addEventListener('click', e => {
        e.preventDefault();
        maybeShowKofiToast({ force: true });
      });
    });
  }

  return { init, addOrReplace, remove, clear, clearSilent, setSuppressed, isSuppressed, copyAll, exportWord, exportExcel, open, close, toggle, showKofiPrompt: maybeShowKofiToast };
})();

if (document.readyState === 'loading'){
  document.addEventListener('DOMContentLoaded', () => ReportBundle.init());
} else {
  ReportBundle.init();
}
