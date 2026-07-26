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
  /* Every caller is a clinical display column — score conversions, SD Index
     changes, RCI statistics, predicted scores. Two old behaviours misfired
     there. Values under 0.0001 switched to exponential notation, so a column
     of figures like -1.42 could suddenly contain "-6.67e-5"; and toFixed keeps
     the sign when a small negative rounds away, printing "-0.00", which is not
     a number anyone writes. A value that rounds to zero at the requested
     precision is displayed as zero, unsigned. Anything genuinely needing more
     precision should ask for more decimal places, not a different notation. */
  const out = n.toFixed(dp);
  return /^-0(\.0*)?$/.test(out) ? out.slice(1) : out;
}
function fmtPct(p){
  if (p == null || isNaN(p)) return '-';
  // A percentile rank lives in the OPEN interval (0, 100) — 0 and 100 are both
  // impossible. Rounding the extreme tails to 2dp used to emit "0.00" and
  // "100.00", with two consequences: report narrative printed "0th percentile"
  // / "100th percentile" via reportOrdinal, and toZ (which returns null for
  // v <= 0 || v >= 100) could no longer convert a percentile-typed row back to
  // a standard score at all. Clamp inside the interval instead.
  // Only reachable at the Wechsler floor/ceiling (FSIQ 40 -> 0.0032th,
  // FSIQ 160 -> 99.9968th) and at T = 90; FSIQ 45 / 155 already resolve fine.
  if (p < 0.1)  return Math.max(p, 0.01).toFixed(2);
  if (p > 99.9) return Math.min(p, 99.99).toFixed(2);
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

/* ---------- NAVIGATION ----------
   Single entry point for switching top-level sections. Every path in — sidebar,
   topnav, home cards, footer links, the auth overlay, and the browser's
   back/forward buttons — goes through here, so the active section, the nav
   highlight, the scroll position and the URL can never disagree.

   opts.navItem  the element that was clicked, when there is one (decides which
                 sidebar group to open, since several nav items share a target)
   opts.history  false to change section without touching history (popstate)
   opts.replace  true to replace the current entry instead of pushing one
   opts.scroll   false to leave the scroll position alone                     */
const NAV_DEFAULT_SECTION = 'home';

function isNavigableSection(id){
  const el = id && document.getElementById(id);
  return !!(el && el.classList.contains('section'));
}
function sectionFromHash(){
  const id = (location.hash || '').replace(/^#/, '');
  return isNavigableSection(id) ? id : null;
}
function writeNavHistory(target){
  const url = `#${target}`;
  const state = { section: target };
  // Re-selecting the current section replaces rather than pushes, so repeated
  // clicks on the same nav item don't pile up identical history entries.
  if (location.hash === url) history.replaceState(state, '', url);
  else history.pushState(state, '', url);
}

/* The global reduce rule in design-system.css zeroes animation DURATIONS but
   not delays, and it can't tell JS not to schedule work. Anything timed from
   script has to ask directly. */
const reducedMotionQuery = window.matchMedia
  ? window.matchMedia('(prefers-reduced-motion: reduce)')
  : null;
function prefersReducedMotion(){
  return !!(reducedMotionQuery && reducedMotionQuery.matches);
}

/* Stagger the arriving page's own blocks a beat behind the page itself.
   Delays are set inline rather than in CSS so they can be omitted entirely
   under reduced motion — a delay left in place there would leave these
   elements invisible for its full duration, since `both` holds the from-state
   and the reduce rule only collapses the duration. */
const STAGGER_STEP_MS = 45;
const STAGGER_MAX = 6;
function staggerSectionContent(section){
  if (!section) return;
  // Only for the fallback transition. Where View Transitions are available the
  // crossfade already carries the arrival, and staggering underneath it reads
  // as two competing animations rather than one move.
  if (document.startViewTransition) return;
  const blocks = [...section.children].filter(el => el.nodeType === 1).slice(0, STAGGER_MAX);
  blocks.forEach((el, i) => {
    el.classList.remove('stagger-in');
    el.style.animationDelay = '';
    if (prefersReducedMotion()) return;
    // Force a reflow so re-navigating to the same section replays the animation
    // rather than being ignored as an unchanged class list.
    void el.offsetWidth;
    el.style.animationDelay = `${i * STAGGER_STEP_MS}ms`;
    el.classList.add('stagger-in');
    // Clear on animationend, but also on a timer: animationend never arrives if
    // the element is replaced mid-flight or the tab stops compositing, and the
    // class uses `both`, so a stranded one would pin the element's final frame.
    let cleared = false;
    const clear = () => {
      if (cleared) return;
      cleared = true;
      el.classList.remove('stagger-in');
      el.style.animationDelay = '';
      el.removeEventListener('animationend', clear);
    };
    el.addEventListener('animationend', clear);
    setTimeout(clear, i * STAGGER_STEP_MS + 600);
  });
}

/* Run the page-switch animation around `swap`. The section change itself stays
   synchronous from the caller's point of view when motion is off, and is
   deferred by one short beat when it's on. Guarded by a timeout as well as
   animationend, because an animation on a hidden or replaced element may never
   fire the event — the swap must happen regardless. */
let navAnimToken = 0;
function runPageTransition(swap){
  const main = document.querySelector('.main');
  // Skip the animation when it can't be seen. A hidden tab doesn't composite,
  // so animationend never arrives and setTimeout is throttled to ~1s — waiting
  // on either would stall the swap for a second with nothing to show for it.
  if (!main || prefersReducedMotion() || !main.animate || document.hidden){
    swap();
    return;
  }

  /* Preferred path. The browser snapshots the page, applies the swap, then
     crossfades the two — no dead frame between out and in, and the section's
     height change happens under the fade rather than snapping the layout.
     Styled by the ::view-transition-* rules in styles.css. */
  if (document.startViewTransition){
    ++navAnimToken;                 // cancel any in-flight fallback animation
    main.classList.remove('is-leaving', 'is-entering');
    document.startViewTransition(() => swap());
    return;
  }

  const token = ++navAnimToken;
  main.classList.remove('is-entering');
  main.classList.add('is-leaving');

  let done = false;
  const finish = () => {
    if (done) return;
    done = true;
    // A newer navigation started mid-flight; let that one own the animation.
    if (token !== navAnimToken){ swap(); return; }
    main.classList.remove('is-leaving');
    swap();
    void main.offsetWidth;
    main.classList.add('is-entering');
    const endEnter = e => {
      if (e && e.target !== main) return;
      main.classList.remove('is-entering');
      main.removeEventListener('animationend', endEnter);
    };
    main.addEventListener('animationend', endEnter);
    setTimeout(() => { if (token === navAnimToken) endEnter(); }, 600);
  };
  main.addEventListener('animationend', function out(e){
    if (e.target !== main) return;
    main.removeEventListener('animationend', out);
    finish();
  });
  setTimeout(finish, 200);
}

function navigateTo(target, opts){
  const o = opts || {};
  if (!isNavigableSection(target)) return false;

  // Gate before anything visible changes, so a declined unlock leaves the UI
  // exactly as it was rather than half-navigated.
  if (target === 'custom-tests' && !isCustomTestsUnlocked()){
    if (!requestCustomTestsUnlock()) return false;
  }

  /* Highlight exactly one nav item. Several can share a target — the five
     Change Analysis methods all point at #change-analysis — so lighting every
     match would show all five as active at once. Prefer the one clicked;
     otherwise take the first that points here. (The topnav highlight is driven
     separately by syncTopnav, which observes the active section.) */
  const matching = [...document.querySelectorAll(`.nav-item[data-target="${target}"]`)];
  const chosen = (o.navItem && matching.includes(o.navItem)) ? o.navItem : matching[0];
  document.querySelectorAll('.nav-item').forEach(n => {
    n.classList.toggle('active', n === chosen);
  });

  const groupNav = (chosen && chosen.closest('.nav-group'))
    ? chosen
    : document.querySelector(`.nav-group .nav-item[data-target="${target}"]`);
  if (groupNav && typeof openOnlyNavGroup === 'function'){
    openOnlyNavGroup(groupNav.closest('.nav-group'));
  }

  // Nav highlight and history update immediately; only the section swap itself
  // waits for the outgoing page to clear, so the URL and the sidebar never lag
  // behind the click.
  if (o.history !== false) writeNavHistory(target);

  const alreadyThere = document.getElementById(target).classList.contains('active');
  const swap = () => {
    document.querySelectorAll('.section').forEach(s => {
      s.classList.toggle('active', s.id === target);
    });
    // The document scrolls, not .main — the redesign moved the scroll container
    // to <body>, which left the old `.main.scrollTop = 0` a silent no-op and
    // dropped you into the middle of every page you navigated to.
    if (o.scroll !== false) window.scrollTo({ top: 0, behavior: 'auto' });
    staggerSectionContent(document.getElementById(target));
  };

  // Re-selecting the page you're already on shouldn't flash it away and back.
  if (alreadyThere) swap();
  else runPageTransition(swap);
  return true;
}
window.navigateTo = navigateTo;

document.querySelectorAll('.nav-item').forEach(el => {
  el.addEventListener('click', () => navigateTo(el.dataset.target, { navItem: el }));
});

window.addEventListener('popstate', e => {
  const target = (e.state && e.state.section) || sectionFromHash() || NAV_DEFAULT_SECTION;
  navigateTo(target, { history: false });
});

/* Restore the section named in the URL on load, so a refresh or a shared link
   lands where it should instead of always bouncing to Home. Runs on
   DOMContentLoaded — late enough that the inline script which builds the
   Change Analysis section has already added it. */
document.addEventListener('DOMContentLoaded', () => {
  const target = sectionFromHash();
  const overlay = document.getElementById('auth-overlay');
  const gated = overlay && !overlay.hidden && overlay.classList.contains('is-visible');
  if (target && target !== NAV_DEFAULT_SECTION && !gated){
    navigateTo(target, { history: false });
  }
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
// Briefly flip a copy button to a "Copied ✓" state for clear confirmation.
function flashCopiedButton(btn){
  if (!btn) return;
  if (btn._restore != null){ clearTimeout(btn._timer); btn.innerHTML = btn._restore; }
  btn._restore = btn.innerHTML;
  btn.classList.add('is-copied');
  btn.innerHTML = '<svg viewBox="0 0 16 16" width="14" height="14" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M3 8.5l3.2 3.2L13 4"/></svg> Copied';
  btn._timer = setTimeout(() => {
    btn.innerHTML = btn._restore; btn.classList.remove('is-copied'); btn._restore = null;
  }, 1500);
}
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
    flashCopiedButton(document.querySelector(`[data-copy="${containerId}"]`));
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
  // Bottom border under the last row that actually has CONTENT. Using the last
  // non-empty row (rather than tbody tr:last-child) means the bottom rule never
  // goes missing when a table happens to end with an empty/placeholder row.
  clone.querySelectorAll('.apa-table tbody').forEach(tbody => {
    const rows = [...tbody.querySelectorAll('tr')];
    let lastRow = null;
    for (let k = rows.length - 1; k >= 0; k--){
      if (rows[k].textContent.trim()){ lastRow = rows[k]; break; }
    }
    if (lastRow){
      lastRow.querySelectorAll('td, th').forEach(td => {
        const existing = td.getAttribute('style') || '';
        td.setAttribute('style', existing + 'border-bottom:1.5pt solid #000;padding-bottom:3pt;');
      });
    }
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
  'rci-basic': {name:'Example index score',sd:'15',r:'0.90',t1:'100',t2:'89',isExample:true},
  'rci-practice': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89',isExample:true},
  'rci-srb': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',t1:'100',t2:'89',isExample:true},
  'rci-crawford': {name:'Example index score',m1:'100',sd1:'15',m2:'103',sd2:'15',r:'0.90',n:'100',t1:'100',t2:'89',isExample:true}
};

/* A row is "blank" if the user has entered nothing in it — the seeded
   placeholder rows that sit below an example fall into this category.
   Returns false for example rows so they're never accidentally swept out.
   Used by the autofill loaders to remove a stale placeholder when the
   user picks a test from the family list. */
/* `groupKey` is included so that an as-yet-unfilled row belonging to a
   user-created test group never counts as a stale placeholder — otherwise
   starting a second custom test would sweep away the first one. */
const ROW_BLANK_KEYS = ['name','group','groupKey','raw','score','t1','t2','sd','m1','sd1','m2','sd2','r','rCorrected','n'];
function isRowBlank(r){
  if (!r) return true;
  if (r.isExample) return false;
  return ROW_BLANK_KEYS.every(k => r[k] == null || r[k] === '');
}
function dropFirstBlankRow(rows){
  const idx = rows.findIndex(isRowBlank);
  if (idx >= 0) rows.splice(idx, 1);
}
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
/* Horizontal geometry of the converter curve, shared by the code that DRAWS it
   and the code that maps a click back to a z score. These were duplicated, and
   drifted: the drag handler still held the pre-redesign 640/94/48 while the
   curve was being drawn at 940/70/36, so clicking the gridline at 70 entered
   63.4 — out by up to 6.6 standard-score points. One source of truth so they
   cannot diverge again. Height and vertical padding stay local to drawCurve;
   only the horizontal mapping is shared. */
const CURVE_GEOM = { W: 940, padL: 70, padR: 36, xMin: -3.5, xMax: 3.5 };

function drawCurve(z, scoreType){
  const svg = document.getElementById('conv-curve');
  /* Cockpit layout: viewBox is ~2.6:1 so when the SVG fills the right column
     (~830px wide) it renders ~320px tall — matching the equivalents column's
     natural height for a balanced row with no dead vertical space. */
  const W = CURVE_GEOM.W, H = 360, padL = CURVE_GEOM.padL, padR = CURVE_GEOM.padR;
  /* bottomPad was 124 to leave room for the per-row conversion table beneath
     the baseline (Standard / Scaled / T / z / Percentile / Classification).
     The cockpit layout shows that information in the side column instead,
     so we crop the SVG to ~30px below the baseline. The table-generating
     code (rowsSvg below) is kept but rendered into a hidden <g>. */
  const bottomPad = 30;
  /* topPad was 76 to leave room for the "Average · 50th percentile" callout
     above the curve peak. The callout is dropped in the cockpit layout — the
     side column already shows the percentile and the classification strip
     below shows the band — so we shrink the top margin to just enough for
     the −1 SD / +1 SD labels above the curve. */
  const topPad = 16;
  const curveH = H - topPad - bottomPad;
  const xMin = CURVE_GEOM.xMin, xMax = CURVE_GEOM.xMax;
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
      stroke="${bold ? '#555' : '#999'}" stroke-width="${bold ? 1.6 : 1.0}" opacity="0.55"/>`;
    sdLabelTexts += `<text x="${lx}" y="${(curveY - 12).toFixed(1)}"
      font-family="IBM Plex Sans" font-size="${bold ? 15 : 13}" fill="${bold ? '#444' : '#777'}"
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
  /* Percentage of the distribution in each SD band, DERIVED rather than typed.
     The two outer labels used to read 2.14%, which is P(2 < |Z| < 3) — the
     figure for a band stopping at 3 SD, applied to bands that run to the edge
     of the plot and beyond. The tail past 3 SD was simply dropped, so the six
     labels summed to 99.72% instead of 100%. Computing them from the same
     normal CDF the curve is drawn with keeps them honest and makes them follow
     xMin/xMax automatically if the plot range ever changes. The outer bands are
     deliberately integrated to infinity, not to the plot edge: they label
     "beyond 2 SD", and the sliver past the edge belongs to that band. */
  const bandEdges = [-Infinity, -2, -1, 0, 1, 2, Infinity];
  const bandPcts = bandEdges.slice(0, -1).map((lo, i) => {
    const hi = bandEdges[i + 1];
    const area = (normCDF(hi) - normCDF(lo)) * 100;
    return {
      zMid: [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5][i],
      pct: area.toFixed(2) + '%'
    };
  });
  let pctLabels = '';
  bandPcts.forEach(({ zMid, pct: p }) => {
    const bx = xScale(zMid);
    const curveY = yScale(normPDF(zMid));
    const labelY = curveY + (base - curveY) * 0.45;
    pctLabels += `<text x="${bx.toFixed(1)}" y="${labelY.toFixed(1)}"
      font-family="IBM Plex Mono" font-size="12" fill="#888" text-anchor="middle" opacity="0.85">${p}</text>`;
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
    twoTail += `<text x="${cx-10}" y="${base-8}" font-family="IBM Plex Mono" font-size="11" fill="#9A9285" text-anchor="end">${belowN}% scored below ◂</text>`;
  if (cx < W - padR - 70)
    twoTail += `<text x="${cx+10}" y="${base-8}" font-family="IBM Plex Mono" font-size="11" fill="#9A9285" text-anchor="start">▸ ${aboveN}% scored above</text>`;

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
    <g class="conv-curve-rows" aria-hidden="true">${rowsSvg}</g>
    ${zLine}${zDot}${zInner}
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
  // Geometry comes from CURVE_GEOM so this can never drift from what
  // drawCurve actually renders — see the note there.
  const { W, padL, padR, xMin, xMax } = CURVE_GEOM;
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

/* Rows are grouped under a header by their group KEY, not their display name.
   Autofilled rows use the database family name as both (no groupKey). Rows the
   user creates get a synthetic, stable groupKey so the header stays a single
   group while it is being renamed — and so two tests can share a name (or be
   left unnamed) without silently merging into one block. */
let batteryGroupSeq = 0;
function batteryGroupKeyOf(r){ return (r && (r.groupKey || r.group)) || ''; }
function batteryGroupNameOf(key){
  const row = batteryRows.find(r => batteryGroupKeyOf(r) === key);
  return row ? (row.group || '') : '';
}
function batteryGroupIsCustom(key){
  return batteryRows.some(r => batteryGroupKeyOf(r) === key && r.groupKey);
}
/* Header text for a group: the user's name, the family name with its age band
   stripped, or a placeholder while a new test is still unnamed. */
function batteryGroupLabel(key){
  const name = batteryGroupNameOf(key);
  if (name) return stripAgeRange(name);
  return batteryGroupIsCustom(key) ? 'Untitled test' : '';
}

function batteryAddRow(initial){
  batteryRows.push(initial || { name:'', raw:'', score:'' });
  renderBattery();
}
/* Start a new custom test: its own header, its own score type, ready to name. */
function batteryAddGroup(scoreType){
  const groupKey = `cg-${++batteryGroupSeq}`;
  dropFirstBlankRow(batteryRows);
  batteryRows.push({ name:'', raw:'', score:'', group:'', groupKey, scoreType });
  renderBattery();
  const input = document.querySelector(`#bat-table [data-group-rename="${groupKey}"]`);
  if (input) input.focus();
}
/* Append one more row to an existing group, directly beneath its last row. */
function batteryAddRowToGroup(key){
  let lastIdx = -1;
  batteryRows.forEach((r, i) => { if (batteryGroupKeyOf(r) === key) lastIdx = i; });
  if (lastIdx < 0) return;
  const sibling = batteryRows[lastIdx];
  batteryRows.splice(lastIdx + 1, 0, {
    name:'', raw:'', score:'',
    group: sibling.group,
    groupKey: sibling.groupKey,
    scoreType: sibling.scoreType
  });
  renderBattery();
  const inputs = document.querySelectorAll(`#bat-table tbody input[data-f="name"]`);
  if (inputs[lastIdx + 1]) inputs[lastIdx + 1].focus();
}
function batteryRenameGroup(key, name){
  batteryRows.forEach(r => { if (batteryGroupKeyOf(r) === key) r.group = name; });
}
function batteryRemove(i){ batteryRows.splice(i, 1); renderBattery(); }
function batteryRemoveGroup(key){
  batteryRows = batteryRows.filter(r => batteryGroupKeyOf(r) !== key);
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
// Infer the score type for an INDIVIDUAL subtest. A single test family can mix
// metrics (e.g. D-KEFS Advanced Color-Word Interference has scaled "Net Correct
// Responses" alongside standard-score "...Index" measures). We decide per subtest
// using, in order of reliability:
//   1. the subtest's own normative mean from the database (most reliable) -
//      nearest standard metric: z≈0, scaled≈10, T≈50, standard≈100;
//   2. keywords in the subtest's own name ("...Index", "IQ", "Quotient", etc.);
//   3. the family-name heuristic (legacy fallback).
function inferScoreTypeForSubtest(familyName, subtestName, stats){
  // 1. Data-driven: classify by the normative mean's nearest standard metric.
  let mean = null;
  if (stats){
    if (typeof stats.m1 === 'number') mean = stats.m1;
    else if (typeof stats.m2 === 'number') mean = stats.m2;
  }
  if (mean != null && isFinite(mean)){
    if (mean < 5)  return 'z';        // z-scores  (mean ≈ 0)
    if (mean < 30) return 'scaled';   // scaled    (mean ≈ 10)
    if (mean < 75) return 't';        // T-scores  (mean ≈ 50)
    return 'standard';                // standard/index (mean ≈ 100)
  }
  // 2. Name-driven: the subtest's own name often states its metric.
  const s = (subtestName || '').toLowerCase();
  if (/\bindex\b/.test(s) || /\bindices\b/.test(s) || /\biq\b/.test(s) ||
      /\bquotient\b/.test(s) || /\bcomposite\b/.test(s) || /\bstandard\s+score\b/.test(s)){
    return 'standard';
  }
  if (/\bt[-\s]?score\b/.test(s)) return 't';
  if (/\bscaled\b/.test(s)) return 'scaled';
  // 3. Fall back to the family-level heuristic.
  return inferScoreType(familyName);
}
// Strip age-band suffixes for compact group labels in editable and APA tables.
function stripAgeRange(name){
  if (!name) return name;
  return name
    .replace(/\s*·\s*Ages?\s+[\d–\-–]+\s*$/i, '')
    .replace(/\s*·\s*All\s+Ages\s*$/i, '')
    .replace(/\s*\(all ages\)\s*$/i, '');
}
/* Change-Analysis ONLY: the reliability-comparison qualifier for a family, so
   the clinician sees which comparison the coefficients represent.
   - CVLT-3 ships alternate-form reliability (Standard ↔ Alternate Form; CVLT-3
     Manual, Delis et al., 2017, Table 3.4) — it publishes no same-form retest.
   - RBANS data here is same-form test–retest (Form A → Form A; RBANS Update
     Manual, Randolph 2012, Tables 3.8–3.9), NOT alternate-form.
   Returns '' for families with no special qualifier. Used in the RCI dropdown,
   table group headers, and footnote — NOT on Score Tables or any other page. */
function caReliabilityQualifier(family){
  const f = family || '';
  if (/^CVLT-3\b/i.test(f)) return 'Standard to Alternate Form';
  if (/^RBANS\b/i.test(f)){
    const m = f.match(/\(Form ([BCD])\)/);          // alternate-form variants
    if (m) return `Form A → Form ${m[1]} (alternate form)`;
    return 'Form A → Form A (test-retest)';          // default = same-form retest
  }
  return '';
}
// True for RBANS alternate-form entries — change-analysis-only data (reliable
// change with a different form). Filtered out of Score Tables and other single-
// administration pages so they don't clutter them.
function isAltFormFamily(name){
  return /\(Form [BCD]\)/.test(name || '');
}
// Change-Analysis dropdown: the retest/delay interval the coefficients were
// measured over, shown on each age-band item so the clinician can match it to
// their own reassessment gap. Sources: CVLT-3 Manual Table 3.4; RBANS Update
// Tables 3.8-3.12; WAIS-IV / WISC-V / WMS-IV / D-KEFS / CVLT-C technical manuals.
// D-KEFS Advanced must be tested before plain D-KEFS (more specific prefix).
function caIntervalLabel(family){
  const f = family || '';
  if (/^CVLT-3\b/i.test(f))          return 'median ~20 days';
  if (/^CVLT-C\b/i.test(f))          return 'median ~28 days';
  if (/^D-KEFS Advanced\b/i.test(f)) return 'mean ~28 days';
  if (/^D-KEFS\b/i.test(f))          return 'mean ~25 days';
  if (/^RBANS\b/i.test(f)){
    if (/\(Form [BCD]\)/.test(f)) return '1-7 days';
    if (/12-19/.test(f))          return '14-31 days';
    if (/20-89/.test(f))          return '~9 months';
  }
  if (/^WAIS-IV\b/i.test(f)) return 'mean ~22 days';
  if (/^WISC-V\b/i.test(f))  return 'mean ~26 days';
  if (/^WMS-IV\b/i.test(f))  return 'mean ~23 days';
  return '';
}
// Change-Analysis display rename for a family base name, where the visible
// edition should differ from the shared normDB key. The RBANS data here is the
// RBANS Update (Randolph, 2012), so show "RBANS Update" on this page only; the
// internal "(Form B/C/D)" marker is dropped because the qualifier shows it.
function caFamilyDisplayName(base){
  return (base || '')
    .replace(/\s*\(Form [BCD]\)/, '')
    .replace(/^RBANS\b/, 'RBANS Update');
}
// Group-header display label for the Change-Analysis tables (base name, age
// range stripped, edition rename + reliability qualifier appended when applicable).
function caGroupDisplay(group){
  const base = caFamilyDisplayName(stripAgeRange(group));
  const q = caReliabilityQualifier(group);
  return q ? `${base} · ${q}` : base;
}
function scoreTypeLabel(type){
  return {scaled:'Scaled Score', standard:'Standard Score', t:'T-Score', z:'Z-Score'}[type] || 'Score';
}
function scoreTypeAbbr(type){
  return {scaled:'Scaled', standard:'Standard', t:'T-Score', z:'Z-Score'}[type] || '';
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
  ['bat-prem-score','bat-prem-threshold','bat-prem-link-btn'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.disabled = !enabled;
  });
}
function getBatteryPremorbidComparison(){
  const enabled = document.getElementById('bat-prem-enable')?.checked;
  const scoreEl = document.getElementById('bat-prem-score');
  /* Prefer the unrounded estimate carried by applyPremorbidLink; fall back to
     the box itself when the user has typed a value in by hand. */
  const exact = parseFloat(scoreEl?.dataset.estimateExact);
  const estimate = Number.isFinite(exact) ? exact : parseFloat(scoreEl?.value);
  if (!enabled || isNaN(estimate)) return null;
  const lo  = parseFloat(scoreEl.dataset.lowerBound);
  const hi  = parseFloat(scoreEl.dataset.upperBound);
  const see = parseFloat(scoreEl.dataset.see);
  return {
    estimate,
    lowerBound: Number.isFinite(lo)  ? lo  : null,
    upperBound: Number.isFinite(hi)  ? hi  : null,
    see:        Number.isFinite(see) ? see : null,
    modelLabel: scoreEl.dataset.modelLabel || null,
    ciLabel:    scoreEl.dataset.ciLabel || null
  };
}
/* Asterisk logic — controlled by the "SD / SEE-based CI" toggle.
   SD MODE (bat-prem-threshold = "sd", default):
       *   ≥ 1   SD below premorbid point estimate
       **  ≥ 1.5 SD below premorbid point estimate
       *** ≥ 2   SD below premorbid point estimate
   SEE MODE (bat-prem-threshold = "see", requires autofill SEE):
       *   below estimate − 1.645·SEE   (90% CI lower bound)
       **  below estimate − 1.960·SEE   (95% CI lower bound)
       *** below estimate − 2.576·SEE   (99% CI lower bound)
   SEE mode silently falls back to SD mode when SEE is not available.    */
const PREMORBID_CI_Z = { ninety: 1.645, ninetyFive: 1.960, ninetyNine: 2.576 };
/* The EFFECTIVE premorbid flagging mode, which is not always the selected one:
   SEE mode silently falls back to SD mode when the linked model has no usable
   SEE. Every consumer must ask this function rather than reading the control
   directly, or one of them ends up describing a mode the table is not using —
   which is exactly how the exported APA note came to misstate the asterisks. */
function batteryPremorbidMode(prem){
  const selected = document.getElementById('bat-prem-threshold')?.value === 'see' ? 'see' : 'sd';
  return (selected === 'see' && prem && Number.isFinite(prem.see) && prem.see > 0) ? 'see' : 'sd';
}
function batteryPremorbidStars(ss, prem){
  if (!prem || !Number.isFinite(ss) || !Number.isFinite(prem.estimate)) return '';
  if (batteryPremorbidMode(prem) === 'see'){
    const t90 = prem.estimate - PREMORBID_CI_Z.ninety      * prem.see;
    const t95 = prem.estimate - PREMORBID_CI_Z.ninetyFive  * prem.see;
    const t99 = prem.estimate - PREMORBID_CI_Z.ninetyNine  * prem.see;
    if (ss <= t99) return '***';
    if (ss <= t95) return '**';
    if (ss <= t90) return '*';
    return '';
  }
  const diffSd = (prem.estimate - ss) / 15;
  if (diffSd >= 2)   return '***';
  if (diffSd >= 1.5) return '**';
  if (diffSd >= 1)   return '*';
  return '';
}

/* Read available premorbid model estimates from the Premorbid page's
   results table. Returns an array of {label, fsiq, lo, hi, see}. Rows
   that haven't been computed (showing "-") are skipped. SEE is needed
   to derive multi-tier CI thresholds for the asterisk logic. */
/* Premorbid estimates offered to the Score Tables asterisk logic.

   Reads preState.estimateRows, NOT the rendered table. calcPremorbid caches the
   rows there with `val` unrounded, whereas the table displays fmtIntOrDash(val)
   — already rounded to a whole number. Scraping the cell fed a rounded figure
   into a further calculation and could flip an asterisk tier: a true estimate
   of 104.6 renders as 105, and against a score of 90 that turns 0.973 SD (no
   asterisk) into exactly 1.000 SD (asterisk).

   The CI bounds are still derived the same way the table derives them, so the
   status line beside the link keeps matching what the Estimates tab shows. */
function getPremorbidEstimateOptions(){
  const rows = (typeof preState === 'object' && preState) ? preState.estimateRows : null;
  if (!Array.isArray(rows)) return [];
  const mult = Number.isFinite(preState.ciMult) ? preState.ciMult : 1.645;
  const out = [];
  rows.forEach(r => {
    if (!r || !r.name || !Number.isFinite(r.val)) return;
    const see = Number.isFinite(r.see) ? r.see : null;
    const margin = see != null ? Math.round(mult * see) : null;
    out.push({
      label: r.name,
      fsiq: r.val,                                        // unrounded
      lo: margin != null ? Math.round(r.val) - margin : null,
      hi: margin != null ? Math.round(r.val) + margin : null,
      see
    });
  });
  return out;
}
/* Single source of truth for the #pre-ci control, which emits "0.90" / "0.95".
   FIVE places used to parse this value independently, and one of them compared
   it against '95' — a string the <select> can never produce — so it returned
   the 90% branch unconditionally while the bounds printed beside it were
   computed at 1.96. A genuine 95% interval was labelled "90% CI" on the
   Battery premorbid link. Route everything through here; do not re-parse the
   control anywhere else. Accepts either "0.95" or "95" so a future change to
   the option values cannot resurrect the mismatch. */
function premorbidCi(){
  const el = document.getElementById('pre-ci');
  const n = parseFloat(el ? el.value : '');
  const is95 = Number.isFinite(n) && (n > 1 ? n : n * 100) >= 95;
  return {
    mult:  is95 ? 1.96    : 1.645,   // z multiplier for the bounds
    label: is95 ? '95% CI' : '90% CI', // "90% CI" form, for the battery link
    short: is95 ? '95%'    : '90%'     // "90%" form, for table column headers
  };
}
function getPremorbidCiLevel(){ return premorbidCi().label; }
/* Update the link-status block beside the autofill button. When a link
   is active, the autofill button hides and the status block takes its
   slot in the source row. The status shows all three asterisk-threshold
   CIs (90/95/99%) so the clinician can see what each asterisk means.

   Asterisks are rendered with spaces (* * *) instead of (***) — the
   spaced form avoids the bold-rendering quirk where the middle glyph
   appears misaligned against its neighbours. */
function updatePremorbidLinkStatus(){
  const scoreEl  = document.getElementById('bat-prem-score');
  const statusEl = document.getElementById('bat-prem-link-status');
  const linkWrap = document.getElementById('bat-prem-link-wrap');
  if (!scoreEl || !statusEl){ return; }
  const see   = parseFloat(scoreEl.dataset.see);
  const est   = parseFloat(scoreEl.value);
  const lo    = parseFloat(scoreEl.dataset.lowerBound);
  const hi    = parseFloat(scoreEl.dataset.upperBound);
  const label = scoreEl.dataset.modelLabel;
  const ci    = scoreEl.dataset.ciLabel;
  const clearBtn = '<button type="button" class="bat-prem-link-clear" id="bat-prem-link-clear" aria-label="Unlink premorbid estimate">×</button>';
  const linkActive = (Number.isFinite(see) && see > 0 && Number.isFinite(est) && label) || (Number.isFinite(lo) && label);

  /* When the link is active, hide the autofill button (mutually
     exclusive with the status block in the source row). */
  if (linkWrap) linkWrap.style.display = linkActive ? 'none' : '';

  if (Number.isFinite(see) && see > 0 && Number.isFinite(est) && label){
    const mode = batteryPremorbidMode({ see });
    let t1, t2, t3;
    if (mode === 'see'){
      t1 = Math.round(est - 1.645 * see);
      t2 = Math.round(est - 1.960 * see);
      t3 = Math.round(est - 2.576 * see);
    } else {
      t1 = Math.round(est - 1.0 * 15);
      t2 = Math.round(est - 1.5 * 15);
      t3 = Math.round(est - 2.0 * 15);
    }
    statusEl.innerHTML =
      '<span class="bat-prem-link-status-body">' +
        `<strong>${escapeHtml(label)}</strong> ${Math.round(est)}` +
        ` · <span class="bat-prem-link-tier"><strong>*</strong> ≤${t1}</span>` +
        ` · <span class="bat-prem-link-tier"><strong>* *</strong> ≤${t2}</span>` +
        ` · <span class="bat-prem-link-tier"><strong>* * *</strong> ≤${t3}</span>` +
      '</span>' + clearBtn;
    statusEl.hidden = false;
    document.getElementById('bat-prem-link-clear')?.addEventListener('click', clearPremorbidLink);
    return;
  }
  if (Number.isFinite(lo) && label){
    /* Fallback for models without SEE — show just the single chosen CI. */
    const range = Number.isFinite(hi) ? `${lo}–${hi}` : `≥${lo}`;
    statusEl.innerHTML =
      `<span class="bat-prem-link-text">Using <strong>${escapeHtml(label)}</strong>` +
      (ci ? ` · ${escapeHtml(ci)} ${range}` : ` · lower bound ${lo}`) + '</span>' + clearBtn;
    statusEl.hidden = false;
    document.getElementById('bat-prem-link-clear')?.addEventListener('click', clearPremorbidLink);
    return;
  }
  statusEl.innerHTML = '';
  statusEl.hidden = true;
}
function clearPremorbidLink(){
  const scoreEl = document.getElementById('bat-prem-score');
  if (!scoreEl) return;
  delete scoreEl.dataset.lowerBound;
  delete scoreEl.dataset.upperBound;
  delete scoreEl.dataset.see;
  delete scoreEl.dataset.modelLabel;
  delete scoreEl.dataset.ciLabel;
  updatePremorbidLinkStatus();
  renderBattery();
}
function applyPremorbidLink(option){
  const scoreEl = document.getElementById('bat-prem-score');
  if (!scoreEl || !option) return;
  /* Mark the next legacy-input event as programmatic so the manual-edit
     clear-on-input handler doesn't wipe the data attributes we're about to
     set. */
  scoreEl.dataset.programmaticUpdate = '1';
  /* The visible box shows the rounded estimate, matching the Estimates tab.
     The unrounded value rides along in a data attribute and is what the
     asterisk thresholds actually use — rounding 104.6 to 105 was enough to
     turn 0.973 SD into exactly 1.000 SD and add an asterisk that should not
     be there. Cleared on manual edit alongside the other link metadata. */
  scoreEl.value = String(Math.round(option.fsiq));
  scoreEl.dataset.estimateExact = String(option.fsiq);
  if (option.lo != null) scoreEl.dataset.lowerBound = String(option.lo);
  else delete scoreEl.dataset.lowerBound;
  if (option.hi != null) scoreEl.dataset.upperBound = String(option.hi);
  else delete scoreEl.dataset.upperBound;
  if (option.see != null) scoreEl.dataset.see = String(option.see);
  else delete scoreEl.dataset.see;
  scoreEl.dataset.modelLabel = option.label;
  scoreEl.dataset.ciLabel    = getPremorbidCiLevel();
  /* Mirror to the visible card input so the user sees the value change. */
  const cardScore = document.getElementById('ds-prem-card-score');
  if (cardScore && cardScore.value !== scoreEl.value) cardScore.value = scoreEl.value;
  /* Make sure the checkbox is on so the comparison actually fires. */
  const enable = document.getElementById('bat-prem-enable');
  if (enable && !enable.checked){
    enable.checked = true;
    enable.dispatchEvent(new Event('change', { bubbles: true }));
  }
  delete scoreEl.dataset.programmaticUpdate;
  updatePremorbidLinkStatus();
  renderBattery();
}
/* Build / open / close the popover that lists available premorbid models. */
function openPremorbidLinkPopover(){
  const popover = document.getElementById('bat-prem-link-popover');
  if (!popover) return;
  const options = getPremorbidEstimateOptions();
  if (!options.length){
    popover.innerHTML =
      '<div class="bat-prem-link-empty">No estimates yet. Open the <strong>Premorbid</strong> page and enter inputs to generate a model.</div>';
  } else {
    popover.innerHTML = options.map((o, i) => {
      const range = (o.lo != null && o.hi != null) ? `${o.lo}–${o.hi}` : '—';
      return `<button type="button" class="bat-prem-link-option" data-idx="${i}">` +
        `<span class="bat-prem-link-option-name">${escapeHtml(o.label)}</span>` +
        `<span class="bat-prem-link-option-stats">FSIQ ${Math.round(o.fsiq)} · ${getPremorbidCiLevel()} ${range}</span>` +
      '</button>';
    }).join('');
    popover.querySelectorAll('.bat-prem-link-option').forEach(btn => {
      btn.addEventListener('click', () => {
        const i = parseInt(btn.dataset.idx, 10);
        applyPremorbidLink(options[i]);
        closePremorbidLinkPopover();
      });
    });
  }
  popover.classList.add('show');
  // dismiss on outside click
  setTimeout(() => document.addEventListener('click', premorbidLinkOutsideClick), 0);
}
function closePremorbidLinkPopover(){
  document.getElementById('bat-prem-link-popover')?.classList.remove('show');
  document.removeEventListener('click', premorbidLinkOutsideClick);
}
function premorbidLinkOutsideClick(e){
  const pop = document.getElementById('bat-prem-link-popover');
  const btn = document.getElementById('bat-prem-link-btn');
  if (!pop) return;
  if (pop.contains(e.target) || (btn && btn.contains(e.target))) return;
  closePremorbidLinkPopover();
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
/* ---------- TABLE KEYBOARD NAVIGATION ----------
   One model for every data-entry table in the app.

   Tab follows the shape of the source you are reading from. Patient scores are
   transcribed one document at a time - all the baseline scores from one report,
   then all the retest scores from another - so Tab runs DOWN a score column and
   wraps to the top of the next one. Normative parameters (M, SD, r, N) come from
   a single row of a test manual, so Tab runs ACROSS the row and on to the next.

   Enter / Shift+Enter always step down / up the current column.
   Ctrl (or Cmd) + arrows move freely around the grid.

   Plain arrows are deliberately left alone: Left/Right must stay available for
   fixing a typo mid-word, and Up/Down still drive the number spinners. */
function tknCells(scope, sel){
  return Array.from(scope.querySelectorAll(sel))
    .filter(el => el.offsetParent !== null && !el.disabled && !el.readOnly);
}
function tknFocus(el){
  if (!el) return true;
  el.focus();
  if (typeof el.select === 'function') el.select();
  return true;
}
function setupTableKeyNav(scope, config){
  if (!scope) return;
  /* Delegated on the container so it survives the innerHTML re-render that every
     render* function performs. Guarded so repeated setup calls don't stack. */
  if (scope.dataset.tknWired === '1') return;
  scope.dataset.tknWired = '1';

  const sel   = config.selector;
  const keyOf = config.keyOf;
  const resolve = v => (typeof v === 'function' ? v() : (v || []));

  const columnCells = field => tknCells(scope, sel).filter(el => keyOf(el) === field);
  const groupFor = field => resolve(config.rowGroups).find(g => g.includes(field)) || null;

  function tabTarget(el, back){
    const field = keyOf(el);
    const group = groupFor(field);
    if (group){
      /* Row-wise: flattening the group in DOM order runs across the row and then
         continues onto the next row, which is how a manual's table reads. */
      const cells = tknCells(scope, sel).filter(x => group.includes(keyOf(x)));
      const i = cells.indexOf(el);
      if (i === -1 || cells.length < 2) return null;
      return cells[(i + (back ? -1 : 1) + cells.length) % cells.length];
    }
    const cols = resolve(config.columns);
    const ci = cols.indexOf(field);
    if (ci === -1) return null;
    const cells = columnCells(field);
    const i = cells.indexOf(el);
    if (i === -1) return null;
    if (!back && i < cells.length - 1) return cells[i + 1];
    if (back && i > 0) return cells[i - 1];
    /* Off the end of this column - carry on at the adjacent column, skipping any
       column that is currently empty or hidden (SDI drops SD outside raw mode). */
    const n = cols.length;
    for (let step = 1; step <= n; step++){
      const at = (((ci + (back ? -step : step)) % n) + n) % n;
      const next = columnCells(cols[at]);
      if (next.length) return back ? next[next.length - 1] : next[0];
    }
    return null;
  }

  function stepColumn(el, delta){
    const cells = columnCells(keyOf(el));
    const i = cells.indexOf(el);
    if (i === -1 || cells.length < 2) return null;
    return cells[(i + delta + cells.length) % cells.length];
  }

  function stepGrid(el, dx, dy){
    const grid = Array.from(scope.querySelectorAll('tr'))
      .map(tr => tknCells(tr, sel))
      .filter(cells => cells.length);
    let r = -1, c = -1;
    grid.forEach((row, ri) => { const ci = row.indexOf(el); if (ci !== -1){ r = ri; c = ci; } });
    if (r === -1) return null;
    if (dy){
      const row = grid[r + dy];
      return row ? (row[Math.min(c, row.length - 1)] || null) : null;
    }
    return grid[r][c + dx] || null;
  }

  scope.addEventListener('keydown', e => {
    const el = e.target.closest(sel);
    if (!el || !scope.contains(el)) return;
    const mod = e.ctrlKey || e.metaKey;

    if (e.key === 'Tab' && !e.altKey && !mod){
      const t = tabTarget(el, e.shiftKey);
      if (t){ e.preventDefault(); tknFocus(t); }
      return;
    }
    if (e.key === 'Enter' && !e.altKey && !mod){
      const t = stepColumn(el, e.shiftKey ? -1 : 1);
      if (t){ e.preventDefault(); tknFocus(t); }
      return;
    }
    if (mod && !e.altKey && e.key.indexOf('Arrow') === 0){
      const t = e.key === 'ArrowUp'   ? stepGrid(el,  0, -1)
              : e.key === 'ArrowDown' ? stepGrid(el,  0,  1)
              : e.key === 'ArrowLeft' ? stepGrid(el, -1,  0)
              :                         stepGrid(el,  1,  0);
      if (t){ e.preventDefault(); tknFocus(t); }
    }
  });
}
function getBatteryRowReliability(row){
  if (!row.group) return null;
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  if (!db) return null;
  const family = db[row.group];
  if (!family) return null;
  const entry = family[row.name];
  if (!entry || typeof entry !== 'object') return null;
  const r  = Number.isFinite(entry.rCorrected) ? entry.rCorrected :
             Number.isFinite(entry.r)           ? entry.r           : null;
  const sd = Number.isFinite(entry.sd1)         ? entry.sd1         : null;
  return (r !== null && sd !== null) ? { r, sd } : null;
}
/* Hard limits imposed by the score scale itself, used to keep an interval from
   printing a value the scale cannot produce. Only scaled scores have a
   universal ceiling (Wechsler subtests run 1-19); "standard" and T are generic
   metrics whose usable range varies by instrument, so no ceiling is asserted
   for them, and z is unbounded. */
const BATTERY_SCORE_BOUNDS = { scaled: { min: 1, max: 19 } };

function getBatteryCiHtml(ss, row, level){
  if (!Number.isFinite(ss)) return '';
  const rel = getBatteryRowReliability(row);
  if (!rel) return '';
  const zMult  = level === '90' ? 1.645 : 1.960;
  const sem    = rel.sd * Math.sqrt(1 - rel.r);
  const loRaw  = ss - zMult * sem;
  const hiRaw  = ss + zMult * sem;
  const type   = rowScoreType(row);
  if (type === 'z'){
    return `${Math.round(loRaw * 10) / 10}–${Math.round(hiRaw * 10) / 10}`;
  }
  // The lower end was already floored, but nothing capped the upper end, so a
  // scaled score of 19 printed "17-21" — a value above the top of the scale.
  const bounds = BATTERY_SCORE_BOUNDS[type];
  const lo = Math.max(bounds ? bounds.min : 1, Math.round(loRaw));
  const hi = bounds ? Math.min(bounds.max, Math.round(hiRaw)) : Math.round(hiRaw);
  return `${lo}–${hi}`;
}

function renderBattery(){
  syncBatteryPremorbidControls();
  updateBatteryScoreHeader();
  const cls     = document.getElementById('bat-class').value;
  const ciLevel = document.getElementById('bat-ci-level')?.value || 'off';
  const tbody = document.querySelector('#bat-table tbody');
  tbody.innerHTML = '';
  document.getElementById('bat-table').classList.toggle('ci-enabled', ciLevel !== 'off');
  const ciHead = document.querySelector('#bat-table .bat-ci-head');
  if (ciHead) ciHead.textContent = ciLevel === 'off' ? 'CI' : `${ciLevel}% CI`;
  const repeatBtn = document.getElementById('bat-add-repeat');
  if (repeatBtn && batteryRows.length === 0) repeatBtn.hidden = true;
  let lastGroup = null;
  batteryRows.forEach((r, i) => {
    // Inject a group header when the group changes
    const gKey = batteryGroupKeyOf(r);
    if (gKey && gKey !== lastGroup){
      const ghr = document.createElement('tr');
      ghr.className = 'group-header';
      // A group can hold mixed score types (e.g. scaled subtests + standard indices).
      // Show the shared label when uniform, otherwise "Mixed" (each row shows its own tag).
      const groupTypes = new Set(batteryRows.filter(x => batteryGroupKeyOf(x) === gKey).map(x => rowScoreType(x)));
      const stLabel = groupTypes.size > 1 ? 'Mixed' : scoreTypeLabel([...groupTypes][0] || r.scoreType || inferScoreType(r.group));
      // Custom tests get an editable name; database families keep their fixed one.
      const nameHtml = r.groupKey
        ? `<input class="group-name-input" data-group-rename="${escapeAttr(gKey)}" value="${escapeAttr(r.group)}" placeholder="Name this test" aria-label="Test name" autocomplete="off">`
        : escapeHtml(batteryGroupLabel(gKey));
      ghr.innerHTML = `<td colspan="8">${nameHtml} <span class="type-badge">· ${stLabel}</span>` +
        `<button class="group-add-row" data-add-to-group="${escapeAttr(gKey)}" title="Add a row to this test" aria-label="Add a row to this test">+</button>` +
        `<button class="group-remove" data-rm-group="${escapeAttr(gKey)}" title="Remove group">×</button></td>`;
      tbody.appendChild(ghr);
      lastGroup = gKey;
    } else if (!gKey){
      lastGroup = null;
    }
    const rowType = rowScoreType(r);
    const z = toZ(r.score, rowType);
    const pct = z == null ? '' : fmtPct(normCDF(z) * 100);
    const details = batteryClassificationDetails(r, cls);
    const ss = parseFloat(r.score);
    const ciHtml = ciLevel !== 'off' ? getBatteryCiHtml(ss, r, ciLevel) : '';
    const tr = document.createElement('tr');
    if (gKey) tr.className = 'in-group';
    const abbr = scoreTypeAbbr(rowType);
    const typeTag = abbr ? `<span class="bat-row-type-tag">(${abbr})</span>` : '';
    tr.innerHTML = `
      <td class="row-num">${i+1}${typeTag}</td>
      <td><input type="text" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Subtest name"></td>
      <td class="bat-raw-cell"><input type="number" step="any" data-r="${i}" data-f="raw" value="${escapeAttr(r.raw)}"></td>
      <td><input type="number" step="any" data-r="${i}" data-f="score" value="${escapeAttr(r.score)}"></td>
      <td class="computed bat-ci-cell">${ciHtml}</td>
      <td class="computed">${pct}</td>
      <td class="computed ${details.className}">${details.html}</td>
      <td class="row-actions"><button onclick="batteryRemove(${i})" title="Remove">×</button></td>
    `;
    tbody.appendChild(tr);
  });
  // Wire group-header buttons
  tbody.querySelectorAll('[data-rm-group]').forEach(b => {
    b.addEventListener('click', () => batteryRemoveGroup(b.dataset.rmGroup));
  });
  tbody.querySelectorAll('[data-add-to-group]').forEach(b => {
    b.addEventListener('click', () => batteryAddRowToGroup(b.dataset.addToGroup));
  });
  // Renaming a custom test updates state and the APA preview on every keystroke,
  // but defers the full table rebuild to blur so the caret is not thrown away.
  tbody.querySelectorAll('[data-group-rename]').forEach(inp => {
    inp.addEventListener('input', e => {
      batteryRenameGroup(e.target.dataset.groupRename, e.target.value);
      renderBatteryApa();
    });
    inp.addEventListener('blur', () => renderBattery());
    inp.addEventListener('keydown', e => { if (e.key === 'Enter') e.target.blur(); });
  });
  // In-place updates while typing. Scoped to data-r cells so the group-header
  // name inputs (which carry no row index) don't fall through to here.
  tbody.querySelectorAll('input[data-r]').forEach(inp => {
    inp.addEventListener('input', e => {
      const i = +e.target.dataset.r, f = e.target.dataset.f;
      batteryRows[i][f] = e.target.value;
      const tr = e.target.closest('tr');
      const rowType = rowScoreType(batteryRows[i]);
      const z = toZ(batteryRows[i].score, rowType);
      const cells = tr.querySelectorAll('.computed');
      const ciCell = cells[0]; // bat-ci-cell
      const pctCell = cells[1];
      const clsCell = cells[2];
      pctCell.textContent = z == null ? '' : fmtPct(normCDF(z) * 100);
      const details = batteryClassificationDetails(batteryRows[i], cls);
      clsCell.className = `computed ${details.className}`.trim();
      clsCell.innerHTML = details.html;
      // Refresh unconditionally. The interval depends on the score AND on the
      // row name (the reliability is looked up by subtest name) AND on the
      // row's score type. Gating this on `f === 'score'` meant renaming a row
      // left the previous subtest's interval on screen while renderBatteryApa()
      // below rebuilt the APA preview with the new name — two visible tables
      // disagreeing. Recomputing on every keystroke is cheap; an explicit
      // field list would just reintroduce the same bug the next time a field
      // starts feeding the interval.
      if (ciCell){
        ciCell.innerHTML = ciLevel !== 'off' ? getBatteryCiHtml(parseFloat(batteryRows[i].score), batteryRows[i], ciLevel) : '';
      }
      renderBatteryApa();
    });
  });
    /* Raw and scaled scores are both transcribed a column at a time. */
    setupTableKeyNav(tbody, {
      selector: 'input[data-f]',
      keyOf: el => el.dataset.f,
      columns: ['name', 'raw', 'score']
    });
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
function buildApaTableFromColumns(outId, columns, rows, groupLabelFn, groupDisplayFn){
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
      const groupText = groupDisplayFn ? groupDisplayFn(group) : stripAgeRange(group);
      body += `<tr class="apa-group"><td colspan="${visible.length}">${escapeHtml(groupText)}</td></tr>`;
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

/* ============================================================
   APA TABLE NOTES — single source of truth
   ------------------------------------------------------------
   One canonical phrasing per table, all in one place. Each entry takes
   the render context and returns an array of sentences; apaNoteHtml
   drops the empty ones and joins them.

   Keep them short. These are pasted into a report, where an over-long
   note costs the clinician more to trim than it saves. State what a
   reader cannot infer from the column headers, and nothing else.
   ============================================================ */
const APA_NOTES = {
  'bat': ctx => [
    `Classification follows ${ctx.classification === 'wechsler' ? 'Wechsler conventions' : 'Guilmette et al. (2020)'}.`,
    ctx.mixedTypes ? 'Scores are reported in their native standardised metric.' : '',
    // States the BASIS, not just the level. These intervals use test-retest
    // reliability, so they run wider than manual-published intervals (which
    // use internal consistency); without this line a reader comparing against
    // the manual has no way to know why the numbers differ. Full rationale is
    // in Methods & References.
    ctx.ciLevel && ctx.ciLevel !== 'off'
      ? `Confidence intervals are ${ctx.ciLevel}%, calculated as the obtained score ± z × SEM using test–retest reliability; these are wider than manual-published intervals based on internal consistency.`
      : '',
    // Must follow the EFFECTIVE flagging mode (batteryPremorbidMode). This
    // previously described the SD thresholds unconditionally, so with SEE
    // flagging active the exported note misstated what the asterisks in its
    // own table meant.
    ctx.premorbid != null
      ? (ctx.premorbidMode === 'see'
          ? `Asterisks mark scores falling below the premorbid estimate of ${ctx.premorbid} by more than the model's standard error of estimate: * beyond the 90% bound, ** beyond the 95% bound, *** beyond the 99% bound.`
          : `Asterisks mark scores below the premorbid estimate of ${ctx.premorbid}: * ≥ 1 SD, ** ≥ 1.5 SD, *** ≥ 2 SD.`)
      : ''
  ],
  'sdi': ctx => [
    'SD Δ = (retest − test) ÷ SD.',
    ctx.mixedTypes ? 'Scores are reported in their native standardised metric.' : '',
    `Significance threshold = ${ctx.thresholdLabel}.`,
    '<i>p</i>-values are two-tailed.'
  ],
  'rci': ctx => [
    ctx.methodSentence,
    `Reliable change threshold = ${ctx.thresholdLabel}.`,
    ctx.rSentence,
    ctx.formSentence,
    '<i>p</i>-values are two-tailed.'
  ],
  'pre-estimates': ctx => [
    'FSIQ = Full Scale IQ estimate.',
    `CI = confidence interval based on ${ctx.ciMultiplier} × SEE.`,
    '<i>r</i> = predictor-criterion correlation. SEE = standard error of estimate.'
  ],
  'pre-predict': () => [
    'WAIS-IV indices are predicted from ToPF, education and sex; WMS-IV indices from ToPF and age.',
    'Difference = Achieved − Predicted.',
    'Base rate = estimated % at or below this discrepancy, from a normal model with SD = SEE (negative discrepancies only). These are parametric estimates, not observed standardisation-sample frequencies, and run slightly low against published empirical figures.'
  ],
  'pre-opiepredict': () => [
    'OPIE-4 prorated scores are predicted from age and sex with Vocabulary and/or Matrix Reasoning.',
    'Illustrative only in a UK context; these values should not be quoted as concrete premorbid estimates. The published equations also carry US education, ethnicity and region terms which are not applied, so every patient is scored at the US reference category (12th-grade high-school graduate, not African-American, not resident in the US West). Those categories have no valid UK equivalent.',
    'The three FSIQ rows predict three different prorated criteria, as do the three GAI rows; they are not interchangeable and are not expected to agree.',
    'Difference = Achieved − Predicted.',
    'Base rate = % of the US standardisation sample at or below this discrepancy (ACS Table eA5.12).'
  ]
};
/* Render a registered note as its APA block. Returns '' when every sentence
   is empty, so a table with nothing to qualify carries no note at all. */
function apaNoteHtml(id, ctx){
  const build = APA_NOTES[id];
  if (!build) return '';
  const body = build(ctx || {}).filter(s => s && s.trim()).join(' ');
  return body ? `<div class="apa-note"><strong>Note.</strong> ${body}</div>` : '';
}
/* Fill the on-screen info boxes that mirror an APA note, so the guidance beside
   the interactive table and the note in the generated table stay identical.
   Only context-free notes can be mirrored this way. */
function renderStaticApaNotes(){
  document.querySelectorAll('[data-apa-note]').forEach(el => {
    const build = APA_NOTES[el.dataset.apaNote];
    if (!build) return;
    const body = build({}).filter(s => s && s.trim()).join(' ');
    if (body) el.innerHTML = `<strong>Note.</strong> ${body}`;
  });
}
document.addEventListener('DOMContentLoaded', renderStaticApaNotes);

function renderBatteryApa(){
  const cls      = document.getElementById('bat-class').value;
  const title    = document.getElementById('bat-title').value || 'Test scores';
  const out      = document.getElementById('bat-apa');
  const ciLevel  = document.getElementById('bat-ci-level')?.value || 'off';
  const rawHidden = document.getElementById('bat-table')?.classList.contains('raw-hidden') ?? true;
  const valid    = batteryRows.filter(r => r.name && !r.isExample);
  const completed = valid.filter(r => r.score !== '' && !isNaN(r.score));
  const types    = new Set((completed.length ? completed : valid).map(r => rowScoreType(r)));
  const headerLabel = types.size === 1 ? scoreTypeLabel([...types][0]) : 'Score';

  /* Sync column visibility with the table's current settings before building */
  if (!apaColumnState['bat-apa']) apaColumnState['bat-apa'] = new Set();
  if (rawHidden) apaColumnState['bat-apa'].delete('raw');
  else           apaColumnState['bat-apa'].add('raw');
  if (ciLevel === 'off') apaColumnState['bat-apa'].delete('ci');
  else                   apaColumnState['bat-apa'].add('ci');

  const ciLabel = ciLevel !== 'off' ? `${ciLevel}% CI` : 'CI';
  const columns = [
    { key:'subtest',        label:'Subtest',        num:false, render:r => escapeHtml(r.name) },
    { key:'raw',            label:'Raw Score',       group:'Scores', num:true,  defaultVisible:!rawHidden, render:r => escapeHtml(r.raw || '-') },
    { key:'score',          label:headerLabel,       group:'Scores', num:true,  render:r => escapeHtml(r.score || '') },
    { key:'ci',             label:ciLabel,           group:'Scores', num:true,  defaultVisible:ciLevel !== 'off', render:r => { const ss = parseFloat(r.score); return ciLevel !== 'off' ? getBatteryCiHtml(ss, r, ciLevel) : ''; }},
    { key:'percentile',     label:'Percentile',      group:'Scores', num:true,  render:r => { const z = toZ(r.score, rowScoreType(r)); return z == null ? '' : fmtPct(normCDF(z) * 100); }},
    { key:'classification', label:'Classification',  group:'Interpretation', num:false, render:r => batteryClassificationDetails(r, cls).html }
  ];
  updateApaColumnControls('bat-apa', columns, renderBatteryApa);
  if (valid.length === 0){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Add or select at least one subtest to preview the APA table.</div>';
    return;
  }
  const prem = getBatteryPremorbidComparison();
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(title)}</div>
    ${buildApaTableFromColumns('bat-apa', columns, valid, batteryGroupKeyOf, batteryGroupLabel)}
    ${apaNoteHtml('bat', {
      classification: cls,
      mixedTypes: types.size > 1,
      ciLevel,
      premorbid: prem ? fmt(prem.estimate, 0) : null,
      premorbidMode: prem ? batteryPremorbidMode(prem) : null
    })}
  `;
}

// Auto-fill: append subtest names from a family as a new group
function loadFamilyIntoBattery(family){
  const db = getMergedDB();
  if (!db[family]) return;
  const names = Object.keys(db[family]);
  // A custom family created but not yet populated. Bail BEFORE dropping a blank
  // row, or the user loses a row and sees nothing added.
  if (!names.length){
    showToast(`${family} has no tests yet - add some in the Norms Database`, true);
    return;
  }
  // If the family is already on the table, append only the subtests that aren't
  // there yet, so tests added to a custom family later can still be pulled in.
  const present = new Set(batteryRows.filter(r => batteryGroupKeyOf(r) === family).map(r => r.name));
  const missing = names.filter(n => !present.has(n));
  if (!missing.length){
    showToast(`${family} is already loaded`, true);
    return;
  }
  // Infer the score type PER subtest, not once for the whole family, so mixed
  // families (e.g. scaled subtests + standard-score indices) categorise correctly.
  const newRows = missing.map(name => ({
    name, raw:'', score:'', group:family,
    scoreType: inferScoreTypeForSubtest(family, name, db[family][name])
  }));
  if (present.size){
    // Slot the new subtests in beneath the group's existing rows so they land
    // under the header that is already there rather than starting a second one.
    let lastIdx = -1;
    batteryRows.forEach((r, i) => { if (batteryGroupKeyOf(r) === family) lastIdx = i; });
    batteryRows.splice(lastIdx + 1, 0, ...newRows);
  } else {
    // Sweep out the stale placeholder row that sits below the example before
    // appending the autofilled tests.
    dropFirstBlankRow(batteryRows);
    batteryRows.push(...newRows);
  }
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
  // Opening the list keeps whatever has been typed, rather than resetting the
  // filter. The visible control is the design-system proxy input, which forwards
  // its `input` events here — without the listener below, typing did nothing.
  inp.addEventListener('focus', () => { list.classList.add('show'); filterFamilyListEl(list, inp.value || ''); });
  inp.addEventListener('click', () => { list.classList.add('show'); filterFamilyListEl(list, inp.value || ''); });
  inp.addEventListener('input', () => { list.classList.add('show'); filterFamilyListEl(list, inp.value || ''); });
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
  // Custom families sort alphabetically among ~30 built-ins, so they need a mark
  // to be findable. Styled by .combo-custom-tag in styles.css.
  return isCustom ? '<span class="combo-custom-tag">Custom</span>' : '';
}
function comboFooterHtml(){
  return '<div class="combo-footer"><span class="combo-count">0 selected</span><button class="btn btn-ghost combo-clear" type="button">Clear</button><button class="btn btn-primary combo-add" type="button" disabled>Add selected tests</button></div>';
}
function comboAgeBandNoteHtml(){
  return '<div class="combo-ageband-note"><span class="combo-ageband-note-icon">ℹ</span><span><strong>Age bands</strong>: more age-specific but smaller <em>N</em> (less stable <em>r</em>). <strong>All Ages</strong>: larger <em>N</em>, stronger <em>r</em>. <strong>Greyed time</strong> = the retest interval each <em>r</em> was measured over.</span></div>';
}
function comboCheckboxItemHtml(f, isCustom, indented, groupKey, displayLabel, suffix){
  const cls = 'combo-item combo-check' + (indented ? ' combo-indented' : '');
  let label = displayLabel
    ? escapeHtml(displayLabel)
    : (indented ? ageBandLabel(f) : escapeHtml(f));
  if (suffix) label += ` <span class="combo-interval">· ${escapeHtml(suffix)}</span>`;
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
  // Hide change-analysis-only alternate-form entries from Score Tables.
  const families = Object.keys(db).sort().filter(f => !isAltFormFamily(f));
  // Battery page: collapse age bands to a single entry per family, no
  // age-band note - norms don't affect the resulting table here.
  list.innerHTML = comboFooterHtml() + buildFamilyListHtml(families, { flat: true });
  wireMultiSelectFamilyList(list, families => {
    families.forEach(loadFamilyIntoBattery);
    const inp = document.getElementById('bat-family-input');
    if (inp){ inp.value = ''; inp.focus(); }
  });
}

/* Popover for score type when adding a manual row. The same popover serves
   "+ Add row" (a loose row) and "+ New test" (a row under its own new heading);
   pendingMode records which button opened it. */
(function(){
  const addBtn    = document.getElementById('bat-add');
  const groupBtn  = document.getElementById('bat-add-group');
  const repeatBtn = document.getElementById('bat-add-repeat');
  const popover   = document.getElementById('bat-type-popover');
  if (!addBtn || !popover) return;

  const TYPE_LABELS = { scaled:'Scaled Score', standard:'Standard Score', t:'T-Score', z:'Z-Score' };
  let lastType = null;
  let pendingMode = 'row';

  function setLastType(type){
    lastType = type;
    if (repeatBtn){
      repeatBtn.textContent = `+ ${TYPE_LABELS[type] || type}`;
      repeatBtn.hidden = false;
    }
  }

  function closePopover(){
    popover.classList.remove('is-open');
    addBtn.setAttribute('aria-expanded', 'false');
    if (groupBtn) groupBtn.setAttribute('aria-expanded', 'false');
  }

  function openFor(mode, trigger, e){
    e.stopPropagation();
    const alreadyOpenForThis = popover.classList.contains('is-open') && pendingMode === mode;
    closePopover();
    if (alreadyOpenForThis) return;
    pendingMode = mode;
    popover.classList.add('is-open');
    trigger.setAttribute('aria-expanded', 'true');
  }

  addBtn.addEventListener('click', e => openFor('row', addBtn, e));
  if (groupBtn) groupBtn.addEventListener('click', e => openFor('group', groupBtn, e));

  popover.querySelectorAll('[data-add-type]').forEach(btn => {
    btn.addEventListener('click', () => {
      closePopover();
      const type = btn.dataset.addType;
      if (pendingMode === 'group'){
        batteryAddGroup(type);
      } else {
        setLastType(type);
        batteryAddRow({ name:'', raw:'', score:'', scoreType: type });
      }
    });
  });

  if (repeatBtn){
    repeatBtn.addEventListener('click', () => {
      if (lastType) batteryAddRow({ name:'', raw:'', score:'', scoreType: lastType });
    });
  }

  document.addEventListener('click', e => {
    if (!popover.contains(e.target) && e.target !== addBtn && e.target !== groupBtn){
      closePopover();
    }
  });

})();
document.getElementById('bat-class').addEventListener('change', renderBattery);
document.getElementById('bat-title').addEventListener('input', renderBatteryApa);
document.getElementById('bat-prem-enable').addEventListener('change', renderBattery);
document.getElementById('bat-prem-score').addEventListener('input', e => {
  /* Skip the clear when the value was set by applyPremorbidLink — that
     update is programmatic and should preserve the data attributes. */
  if (e.target.dataset.programmaticUpdate) return;
  /* A manual edit means the user no longer wants the CI-linked anchor.
     Clear the link metadata so the comparison falls back to the point
     estimate they're typing. */
  delete e.target.dataset.lowerBound;
  delete e.target.dataset.upperBound;
  delete e.target.dataset.see;
  delete e.target.dataset.modelLabel;
  delete e.target.dataset.ciLabel;
  delete e.target.dataset.estimateExact;   // else it would shadow the typed value
  updatePremorbidLinkStatus();
  renderBattery();
});
document.getElementById('bat-prem-threshold').addEventListener('change', () => { updatePremorbidLinkStatus(); renderBattery(); });
document.getElementById('bat-ci-level').addEventListener('change', renderBattery);
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
// The divisor is a property of the ROW, not of the table. A battery mixes metrics
// (scaled subtests alongside standard-score indices), and dividing every row by one
// shared SD silently rescales most of them. Rows carry their own .scoreType — set
// from the Add-row popover, or inferred per subtest on auto-fill. The #sdi-type
// select is only the default for rows that never got one. Mirrors rowScoreType().
function sdiRowScoreType(r){
  return (r && r.scoreType) || document.getElementById('sdi-type').value;
}
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
  return (parseFloat(r.t2) - parseFloat(r.t1)) / sdiSdUnit(sdiRowScoreType(r));
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
}
/* The Change-Analysis tables are table-layout:fixed, so a <colgroup> — not the
   cell widths — decides every column. renderSdi owns the SDI one because the
   column COUNT changes with the mode: raw mode adds an SD column, and only
   standardised mode carries the score-type pill in the # column, which needs
   roughly the width Score Tables gives it. Percentages, summing to 100. */
const SDI_COL_WIDTHS = {
  //            #   Subtest  d1  d2  [SD]  SDΔ   p   Sig  ×
  //  11% ≈ 127px: enough for a two-digit row number plus the longest pill,
  //  "(Standard)", on one line. Below that it wraps under the number.
  standardised: [11, 28,     13, 13,       10,  10,  13,  2],
  raw:          [ 3, 28,     12, 12,  12,  10,  10,  11,  2]
};
function syncSdiColgroup(raw){
  const table = document.getElementById('sdi-table');
  if (!table) return;
  const widths = SDI_COL_WIDTHS[raw ? 'raw' : 'standardised'];
  let cg = table.querySelector('colgroup');
  if (!cg){
    cg = document.createElement('colgroup');
    table.insertBefore(cg, table.firstChild);
  }
  while (cg.children.length > widths.length) cg.removeChild(cg.lastChild);
  while (cg.children.length < widths.length) cg.appendChild(document.createElement('col'));
  widths.forEach((w, i) => { cg.children[i].style.width = `${w}%`; });
}
function renderSdi(){
  const raw = sdiMode() === 'raw';
  const cv = parseFloat(document.getElementById('sdi-cv').value);
  document.getElementById('sdi-type-field').classList.toggle('is-hidden', raw);
  document.getElementById('sdi-raw-help').style.display = raw ? 'block' : 'none';
  syncSdiColgroup(raw);
  renderSdiHead();
  const tbody = document.querySelector('#sdi-table tbody');
  tbody.innerHTML = '';
  let lastGroup = null;
  sdiRows.forEach((r, i) => {
    if (r.group && r.group !== lastGroup){
      const ghr = document.createElement('tr');
      ghr.className = 'group-header';
      const colspan = raw ? 9 : 8;
      // A group can hold mixed score types (e.g. scaled subtests + standard indices).
      // Show the shared label when uniform, otherwise "Mixed" (each row shows its own tag).
      // Raw mode divides by the SD typed into the row, so no metric applies.
      let badge = '';
      if (!raw){
        const groupTypes = new Set(sdiRows.filter(x => x.group === r.group).map(x => sdiRowScoreType(x)));
        const stLabel = groupTypes.size > 1 ? 'Mixed' : scoreTypeLabel([...groupTypes][0]);
        badge = ` <span class="type-badge">· ${escapeHtml(stLabel)}</span>`;
      }
      ghr.innerHTML = `<td colspan="${colspan}">${escapeHtml(stripAgeRange(r.group))}${badge}<button class="group-remove" data-rm-sdi-group="${escapeAttr(r.group)}" title="Remove group">×</button></td>`;
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
    const abbr = raw ? '' : scoreTypeAbbr(sdiRowScoreType(r));
    const typeTag = abbr ? `<span class="bat-row-type-tag">(${escapeHtml(abbr)})</span>` : '';
    tr.innerHTML = `
      <td class="row-num">${i+1}${typeTag}</td>
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
  /* SD is a per-row normative value, but it is still entered a column at a time,
     and it only exists in raw mode. */
  setupTableKeyNav(tbody, {
    selector: 'input[data-f]',
    keyOf: el => el.dataset.f,
    columns: () => sdiMode() === 'raw' ? ['name', 't1', 't2', 'sd'] : ['name', 't1', 't2']
  });
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
  const named = sdiRows.filter(r => r.name && !r.isExample);
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
    ${apaNoteHtml('sdi', { thresholdLabel: cvDesc, mixedTypes: !raw && new Set(named.map(sdiRowScoreType)).size > 1 })}
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
  dropFirstBlankRow(sdiRows);
  // Infer the score type PER subtest, not once for the whole family, so mixed
  // families (e.g. scaled subtests + standard-score indices) divide by their own SD.
  subtests.forEach(([name, p]) => {
    sdiRows.push({
      name, t1:'', t2:'', sd: raw ? sdiNormSd(p) : '', group:family,
      scoreType: inferScoreTypeForSubtest(family, name, p)
    });
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
  // Alternate-form entries are reliable-change (RCI) only; keep them out of SDI.
  const families = Object.keys(db).sort().filter(f => !isAltFormFamily(f));
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
(function(){
  const addBtn    = document.getElementById('sdi-add');
  const repeatBtn = document.getElementById('sdi-add-repeat');
  const popover   = document.getElementById('sdi-type-popover');
  const TYPE_LABELS = { standard:'Standard Score', scaled:'Scaled Score', t:'T-Score', z:'Z-Score' };
  let lastType = null;
  // Deliberately does NOT write back to #sdi-type. Doing so re-divided every row
  // already on the table by the newly picked SD. The choice belongs to the row
  // being added and travels on it as .scoreType.
  function setLastType(type){
    lastType = type;
    if (repeatBtn){ repeatBtn.textContent = `+ ${TYPE_LABELS[type] || type}`; repeatBtn.hidden = false; }
  }
  addBtn.addEventListener('click', e => {
    e.stopPropagation();
    const open = popover.classList.toggle('is-open');
    addBtn.setAttribute('aria-expanded', open ? 'true' : 'false');
  });
  popover.querySelectorAll('[data-sdi-type]').forEach(btn => {
    btn.addEventListener('click', () => {
      popover.classList.remove('is-open');
      addBtn.setAttribute('aria-expanded', 'false');
      const type = btn.dataset.sdiType;
      setLastType(type);
      sdiAddRow({ name:'', t1:'', t2:'', sd:'', scoreType: type });
    });
  });
  if (repeatBtn) repeatBtn.addEventListener('click', () => {
    if (lastType) sdiAddRow({ name:'', t1:'', t2:'', sd:'', scoreType: lastType });
  });
  document.addEventListener('click', e => {
    if (!popover.contains(e.target) && e.target !== addBtn){
      popover.classList.remove('is-open');
      addBtn.setAttribute('aria-expanded', 'false');
    }
  });
})();
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
// All four RCI methods share ONE row set, so a person's scores and norms are
// entered once and carry across every method. Each row holds the full superset
// of fields; each method reads the subset it needs. NEVER reassign a method's
// .rows (e.g. `= [...]`) — that breaks the shared link. Mutate in place instead.
const RCI_SHARED_ROWS = [];
const rciState = {
  // useCorrectedR defaults OFF on all four methods. The corrected r is the
  // retest correlation rescaled to the normative sample's variability, so it
  // describes a population with SD 15 / 3 — not the retest sample whose SD
  // (sd1, sd2) is what these formulas multiply it by. Pairing the two mixes
  // populations and understates the error variance by ~20%. WAIS-IV FSIQ:
  //   13.8^2 x (1 - .95) = 9.52   sample SD with the sample's own r   OK
  //   15.0^2 x (1 - .96) = 9.00   normative SD with the corrected r   OK
  //   13.8^2 x (1 - .96) = 7.62   the mix this default used to give   wrong
  // Raw r is also what the source methods specify: Jacobson & Truax take the
  // test-retest correlation and the variance of the initial testing in the
  // study itself; Iverson, McSweeney and Crawford likewise work from the
  // retest sample's own statistics.
  //
  // The user-facing toggle exists on Basic and Practice ONLY. On those two, r
  // is purely an error term, so ticking it widens or narrows the interval and
  // nothing else — a defensible choice if the SD is also set to the normative
  // value. On SRB and Crawford r is a fitted regression SLOPE (r x sd2/sd1),
  // so substituting a population-corrected value moves the predicted score
  // itself and yields a line that was never fitted to anything. Worked example
  // from the database entry with the largest r -> rCorrected gap, RBANS Story
  // Memory Ages 20-89 (m1 11.6, sd1 1.8, m2 12.5, sd2 2.4, r .45, rCorr .80),
  // patient T1 9.8 / T2 10.1:
  //   raw r .45        slope 0.600   predicted 11.42   SEE 2.143   RCI -0.616
  //   corrected r .80  slope 1.067   predicted 10.58   SEE 1.440   RCI -0.333
  // The predicted score moves 0.84 points, the SEE shrinks by 33%, and the
  // slope exceeds 1.0 — i.e. the model now says scores diverge from the mean on
  // retest, which the observed r of .45 flatly contradicts. These two keep the
  // field permanently false; with no checkbox in the markup, nothing sets it.
  'rci-basic':    { rows:RCI_SHARED_ROWS, cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-practice': { rows:RCI_SHARED_ROWS, cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-srb':      { rows:RCI_SHARED_ROWS, cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' },
  'rci-crawford': { rows:RCI_SHARED_ROWS, cv:0.95, useCorrectedR:false, d1:'Date 1', d2:'Date 2', title:'Reliable change analysis' }
};
// Blank rows carry the full superset so any method can read/fill them.
function newRciRow(method){
  return { name:'', sd:'', m1:'', sd1:'', m2:'', sd2:'', r:'', rCorrected:'', n:'', t1:'', t2:'' };
}
function rciAddRow(method){ rciState[method].rows.push(newRciRow(method)); renderRci(method); }
function rciRemove(method, i){ rciState[method].rows.splice(i, 1); renderRci(method); }
function rciRemoveGroup(method, group){
  const rows = rciState[method].rows;
  for (let i = rows.length - 1; i >= 0; i--){ if (rows[i].group === group) rows.splice(i, 1); }
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
/* Session-level "use corrected r" toggle, stored on
   rciState[method].useCorrectedR. OFF by default, and permanently false on SRB/Crawford — see the note on
   rciState. Turning it ON rescales the reliability to the normative population
   without rescaling the SD it is paired with, so it is offered for comparison
   rather than as the default basis. */
/* Which reliability FIELD is in force for a given row — 'rCorrected' only when
   the toggle is on AND that row actually has a corrected value, otherwise 'r'.
   The table binds its reliability cell to this, so the cell shows and edits the
   coefficient the calculation is really using.

   Two things previously went wrong without it. The cell was hard-bound to
   data-f="r", so with the toggle on it displayed the raw coefficient while the
   maths used the corrected one; and because the generic input handler writes
   rciState[m].rows[i][f], overtyping that cell wrote to `r` — a field the
   calculation was ignoring — so the edit appeared to be accepted and changed
   nothing.

   Note this is deliberately per ROW, not per column: within one table some
   rows carry a corrected value and some do not (357 of the 590 database
   entries have none), so a column-level label would misdescribe the mixed
   case. The APA note already words this as "where published". */
function rciReliabilityField(method, row){
  const wantCorrected = !!(rciState[method] && rciState[method].useCorrectedR === true);
  return (wantCorrected && Number.isFinite(parseFloat(row.rCorrected))) ? 'rCorrected' : 'r';
}
/* Render the reliability cell against whichever field is in force. */
function rciReliabilityCell(method, i, row){
  const field = rciReliabilityField(method, row);
  const isCorrected = field === 'rCorrected';
  const title = isCorrected
    ? 'Corrected (attenuation-adjusted) test-retest r. This is the value in use; editing it changes the result. Untick "Use corrected r" to enter a raw coefficient instead.'
    : 'Raw test-retest r. This is the value in use; editing it changes the result.';
  return `<td><input type="number" step="0.01" min="0" max="1" data-m="${method}" data-r="${i}"`
       + ` data-f="${field}"${isCorrected ? ' data-r-corrected="1"' : ''}`
       + ` title="${escapeAttr(title)}" value="${escapeAttr(row[field])}"></td>`;
}
function rciEffectiveR(method, row){
  const wantCorrected = !!(rciState[method] && rciState[method].useCorrectedR === true);
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
  // Crawford & Garthwaite (2007) Eq. 4. The (n-1)/(n-2) factor gives the unbiased
  // residual SD, which is what the t reference distribution on n-2 df requires.
  const see = sd2 * Math.sqrt((1 - rel*rel) * (n - 1) / (n - 2));
  // Crawford & Garthwaite (2007) Eq. 5: standard error of prediction for a new case
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
  /* Reports significance only, NOT a clinical direction.
     This used to map the sign of the statistic straight onto "Reliable
     improvement" / "Reliable decline", which is wrong for any measure where a
     higher score is a worse result — and normDB carries 119 such entries
     (intrusions, perseverations, errors, false positives, repetitions). A
     patient making 16 recall intrusions at T1 and 2 at T2 has plainly
     improved, but the statistic is negative, so the table asserted "Reliable
     decline" in a report.
     Rather than infer valence from a per-measure direction flag the database
     does not have, the label states only that the change exceeded the
     threshold. The signed statistic sits in the adjacent column (fmt(calc.rci,
     2)), so the direction of movement remains visible and the clinician
     applies their own knowledge of which way is better for that measure.
     Uses sig-yes / sig-no, the same neutral pair the SDI table already uses
     for the same question; sig-improve and sig-decline are no longer emitted
     by anything. */
  return { label:'Reliable change', cls:'sig-yes' };
}



// Make keyboard entry match the user's current task.
// If the cursor is in Date 1/Date 2, Tab stays within the patient-score cells.
// If the cursor is in the test/normative fields, Tab stays within those fields.
function setupRciContextualTabbing(tbody){
  if (!tbody) return;

  tbody.querySelectorAll('input[data-f]').forEach(el => {
    const field = el.dataset.f;
    el.tabIndex = 0;
    if (field === 't1' || field === 't2'){
      el.classList.add('patient-score-tab-target');
      el.title = 'Tab moves down this score column, then on to the next column.';
      el.setAttribute('aria-label', field === 't1' ? 'Date 1 score' : 'Date 2 score');
    } else {
      el.classList.add('norm-tab-target');
      // Append rather than replace: some normative cells carry a more specific
      // tooltip set in their own markup (the reliability cell states which
      // coefficient it holds and that editing it changes the result), and
      // overwriting it silently discarded that. Guarded so repeated calls
      // cannot stack the hint.
      const tabHint = 'Tab moves across the normative fields, then on to the next test.';
      const existing = (el.getAttribute('title') || '').trim();
      el.title = existing && !existing.includes(tabHint) ? existing + ' ' + tabHint : tabHint;
    }
  });

  /* Patient scores are read off two separate reports, so Tab runs down one column
     at a time. The normative fields come from one row of a manual, so Tab runs
     across them. See setupTableKeyNav. */
  setupTableKeyNav(tbody, {
    selector: 'input[data-f]',
    keyOf: el => el.dataset.f,
    columns: ['t1', 't2'],
    rowGroups: [['name', 'sd', 'r', 'm1', 'sd1', 'm2', 'sd2', 'n']]
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
      ghr.innerHTML = `<td colspan="${colCount}">${escapeHtml(caGroupDisplay(r.group))}<button class="group-remove" data-rm-rci-group="${escapeAttr(r.group)}" data-method="${method}" title="Remove group">×</button></td>`;
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
        ${rciReliabilityCell(method, i, r)}
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
        ${rciReliabilityCell(method, i, r)}
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
        ${rciReliabilityCell(method, i, r)}
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
        ${rciReliabilityCell(method, i, r)}
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
  const valid = st.rows.filter(r => r.name && !r.isExample);
  const cvLabel = st.cv === 0.95 ? '95%' : '90%';
  const cvLabelZ = st.cv === 0.95 ? '95% (z = 1.96)' : '90% (z = 1.645)';
  const methodSentence = {
    'rci-basic':    'RCI (z) is computed per Jacobson and Truax (1991).',
    'rci-practice': 'RCI (z) is computed per Iverson (2001), adjusted for practice effects.',
    'rci-srb':      'RCI (z) is computed per McSweeney et al. (1993); Ŷ₂ = predicted retest score.',
    'rci-crawford': '<i>t</i>(RB) is the Crawford regression-based reliable-change statistic.'
  }[method];
  /* Crawford's critical value is a t quantile on df = N − 2, so it varies by
     row and cannot be stated as a single number the way the z-based methods
     can. Name the distribution instead of quoting a value that is never used. */
  const thresholdLabel = method === 'rci-crawford'
    ? `${cvLabel} (two-tailed <i>t</i>, df = <i>N</i> − 2)`
    : cvLabelZ;
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
  let rSentence = '';
  {
    const wantCorrected = st.useCorrectedR === true;   // default OFF; always false on SRB/Crawford
    if (!wantCorrected){
      rSentence = 'Raw test-retest <i>r</i> was used.';
    } else {
      const fallbackTests = valid
        .filter(r => !Number.isFinite(parseFloat(r.rCorrected)))
        .map(r => r.name);
      if (fallbackTests.length === 0){
        rSentence = 'Corrected (attenuation-adjusted) test-retest <i>r</i> was used.';
      } else {
        // Naming every fallback test made this note grow without limit — it is
        // the reason the Change Analysis note ran to several hundred characters.
        // Name a few, then count.
        const named = fallbackTests.length <= 3
          ? escapeHtml(fallbackTests.join(', '))
          : `${fallbackTests.length} tests`;
        rSentence = `Corrected (attenuation-adjusted) test-retest <i>r</i> was used where published; raw <i>r</i> for ${named}.`;
      }
    }
  }
  // Change-Analysis: state which reliability comparison + source was applied,
  // per test family present (CVLT-3 = alternate-form, RBANS = same-form retest).
  const formParts = [];
  if (valid.some(r => /^CVLT-3\b/i.test(r.group || ''))){
    formParts.push('CVLT-3 coefficients are alternate-form reliabilities (Delis et al., 2017).');
  }
  const rbansAA  = valid.some(r => /^RBANS\b/i.test(r.group || '') && !isAltFormFamily(r.group || ''));
  const rbansAlt = valid.some(r => /^RBANS\b/i.test(r.group || '') &&  isAltFormFamily(r.group || ''));
  if (rbansAA){
    formParts.push('RBANS Update coefficients are Form A retest reliabilities (Randolph, 2012).');
  }
  if (rbansAlt){
    formParts.push('RBANS Update coefficients are Form A to B/C/D alternate-form reliabilities (Randolph, 2012).');
  }
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">${escapeHtml(st.title)}</div>
    ${buildApaTableFromColumns(outId, columns, valid, r => r.group, caGroupDisplay)}
    ${apaNoteHtml('rci', { methodSentence, thresholdLabel, rSentence, formSentence: formParts.join(' ') })}
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
/* Use corrected r toggle. Present on the Basic and Practice pages only, and
   unticked by default. SRB and Crawford deliberately have no checkbox because
   r is a fitted regression slope there — see the note on rciState. This binds
   by class, so those two simply never get a listener and their useCorrectedR
   stays false. */
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
    } else if (flat && !members.some(isCustom)){
      // Pick the "All Ages" canonical entry if present, else the first.
      const canon = members.find(m => /·\s*All\s+Ages\s*$/i.test(m)) || members[0];
      html += comboCheckboxItemHtml(canon, false, false, groupKey, base);
    } else if (flat){
      // Custom families sharing a base name are separate user-authored tests
      // that can hold entirely different subtests — collapsing them to one
      // canonical entry makes every variant but the first unselectable here.
      members.forEach(f => { html += comboCheckboxItemHtml(f, isCustom(f), false, groupKey); });
    } else {
      const headingText = (opts && typeof opts.headingLabel === 'function') ? opts.headingLabel(base) : base;
      html += `<div class="combo-group-heading" data-group="${escapeAttr(groupKey)}">${escapeHtml(headingText)}</div>`;
      // Wrap age-banded variants in a flex row so they render side-by-side
      // as compact pills instead of stacking - cuts vertical scroll roughly
      // in half on long family lists (e.g. CVLT-3 INDICES + TRIALS, etc.).
      html += `<div class="combo-indented-row" data-group="${escapeAttr(groupKey)}">`;
      members.forEach(f => {
        const suffix = (opts && typeof opts.itemSuffix === 'function') ? opts.itemSuffix(f) : '';
        html += comboCheckboxItemHtml(f, isCustom(f), true, groupKey, null, suffix);
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
  // Change-Analysis dropdown: edition rename + reliability-comparison qualifier
  // on the family heading (e.g. RBANS → "RBANS Update", CVLT-3 → "… Standard to
  // Alternate Form").
  const headingLabel = base => {
    const name = caFamilyDisplayName(base);
    const q = caReliabilityQualifier(base);
    return q ? `${name} · ${q}` : name;
  };
  list.innerHTML = comboFooterHtml() + comboAgeBandNoteHtml() + buildFamilyListHtml(families, { headingLabel, itemSuffix: caIntervalLabel });
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
  // Data is shared across all methods now — don't add the same family twice.
  if (rciState[method].rows.some(r => r.group === family)){
    if (typeof showToast === 'function') showToast(`${family} is already loaded`, true);
    renderRci(method);
    return;
  }
  // Store the FULL superset on every row so each method reads the fields it
  // needs (basic: sd + r; practice/SRB: m1/sd1/m2/sd2/r; Crawford also: n).
  const newRows = subtests.map(([name, p]) => ({
    name, group:family,
    sd:  (p.sd1 != null ? String(p.sd1) : ''),
    m1:  (p.m1  != null ? String(p.m1)  : ''),
    sd1: (p.sd1 != null ? String(p.sd1) : ''),
    m2:  (p.m2  != null ? String(p.m2)  : ''),
    sd2: (p.sd2 != null ? String(p.sd2) : ''),
    r:   (p.r   != null ? String(p.r)   : ''),
    rCorrected: (p.rCorrected != null ? String(p.rCorrected) : ''),
    n:   (p.n   != null ? String(p.n)   : ''),
    t1:'', t2:''
  }));
  // Sweep out the stale placeholder row, then append in place — mutate the
  // SHARED array (never reassign .rows, which would break the shared link).
  dropFirstBlankRow(rciState[method].rows);
  rciState[method].rows.push(...newRows);
  renderRci(method);
  // Toast suppressed - the working-report pill is the single feedback channel
}
function clearMethodRows(method){
  rciState[method].rows.length = 0;
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
/* Rebuilds every consumer of the normative database after a custom test is
   added, deleted or imported. Also runs once at top level during boot.

   Nothing may be called from here that is not defined. This function is a
   TOP-LEVEL init statement, so a throw inside it kills every top-level
   statement that follows it in the file - which is roughly a third of the app,
   including the whole Working Report bundle. It used to end with a call to
   refreshReportWriterOptions(), deleted with the Report Writer while the call
   site stayed, and that is exactly what happened: everything after the
   refreshAll() call at the bottom of this file silently stopped running.
   check.js section 17 now asserts that every project function called anywhere
   in the source actually exists. */
function refreshAll(selectedFamily){
  refreshFamilySelect(selectedFamily);
  renderDbList();
  rebuildAllFamilyLists();
  rebuildBatteryFamilyList();
  rebuildSdiFamilyList();
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
  /* Pure normative entry - every value on a row is read off one row of a test
     manual, so Tab runs across the row and then on to the next. */
  setupTableKeyNav(body, {
    selector: 'input[data-k]',
    keyOf: el => el.dataset.k,
    rowGroups: [['name', 'm1', 'sd1', 'm2', 'sd2', 'r', 'n']]
  });
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
function preCiMult(){ return premorbidCi().mult; }
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
  // Scope to the pre-tabs strip — the topnav dropdown shares data-pre-tab
  // attribute names and would otherwise steal the .active class.
  document.querySelectorAll('#premorbid .pre-tabs .tab').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.pre-tab-content').forEach(c => c.classList.remove('active'));
  const tabBtn = document.querySelector(`#premorbid .pre-tabs .tab[data-pre-tab="${tabId}"]`);
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
  // Initialise on "estimates" — Inputs is now an always-visible aside in the
  // restructured layout, so its tab is hidden and shouldn't be the default.
  // Set active classes directly (avoiding switchPreTab's scroll behaviour
  // which can fire before the section is visible). Scope to .pre-tabs because
  // the topnav dropdown shares the data-pre-tab attribute.
  document.querySelectorAll('#premorbid .pre-tabs .tab').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.pre-tab-content').forEach(c => c.classList.remove('active'));
  const initBtn = document.querySelector('#premorbid .pre-tabs .tab[data-pre-tab="estimates"]');
  if (initBtn) initBtn.classList.add('active');
  const initContent = document.getElementById('pre-estimates');
  if (initContent) initContent.classList.add('active');
  updatePreNav('estimates');
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

  // Wire achieved inputs — shared by index, so a change here autofills the
  // matching OPIE-4 rows too (and vice versa) via setSharedAchieved.
  document.querySelectorAll('[data-pre-ach]').forEach(inp => {
    inp.addEventListener('input', e => {
      setSharedAchieved(e.target.dataset.preAch, e.target.value, e.target);
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
  // CI level comes from premorbidCi() — see its note; do not re-parse #pre-ci.
  // Sex coding differs between models:
  //   TOPF / WAIS-IV manual demographic equations:  Female = 1, Male = 2
  //   OPIE-4 (Holdnack et al., 2013):               Female = 0, Male = 1
  const sexC_topf = sex === 'Female' ? 1 : sex === 'Male' ? 2 : null;
  const sexC_opie = sex === 'Female' ? 0 : sex === 'Male' ? 1 : null;
  const occC = occ ? OCC_CODE[occ] : null;
  const mult = preCiMult();
  const ciPct = premorbidCi().short;

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
  // Computing a value additionally requires age AND sex:
  //   - age must fall in the fitted 16-90 range; the Age3 term is unbounded and
  //     a mistyped age extrapolates hard (age 110 inflates FSIQ by ~18 points).
  //   - sex is NOT optional. OPIE-4 codes Female = 0, so defaulting a blank
  //     field to 0 silently returned the FEMALE equation rather than no answer,
  //     under-stating male patients by up to 5 points with nothing on screen to
  //     signal it. Gate on it the same way age is gated.
  const hasVC = vc != null, hasMR = mr != null;
  const opieAgeOk = age != null && age >= OPIE_AGE_MIN && age <= OPIE_AGE_MAX;
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
    // Compute prediction only once age (in range) and sex are also available
    if (opieAgeOk && sexC_opie != null){
      let pred = c.intercept + (c.age != null ? c.age * age : 0) + c.age3 * Math.pow(age, 3) + c.sex * sexC_opie;
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
    // Round the point estimate and the margin separately so the bounds stay
    // symmetric about the printed estimate. Must match calcOpiePredict's
    // convention exactly or the same model prints different bounds on the
    // Estimates tab and the OPIE-4 tab (it disagreed on ~21% of inputs).
    const margin = (row.see == null) ? null : Math.round(mult*row.see);
    const lo   = (row.val == null || margin == null) ? '-' : String(Math.round(row.val) - margin);
    const hi   = (row.val == null || margin == null) ? '-' : String(Math.round(row.val) + margin);
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
  renderPremorbidForestPlot(rows, mult, ciPct);
}

/* === Forest plot rendering ===
   Modern JAMA-style tabular figure. Five columns: Model · r · SEE ·
   Estimate (CI) · Forest plot. Filled diamonds mark the point estimate
   (sized by inverse SEE); horizontal whiskers span the CI; one thin
   vertical reference at the population mean (FSIQ 100). All-charcoal
   grayscale palette, IBM Plex Sans throughout, tabular figures for
   clean column alignment. */
function renderPremorbidForestPlot(rows, mult, ciPct){
  const wrap = document.getElementById('pre-forest-plot-wrap');
  const svg  = document.getElementById('pre-forest-plot');
  if (!wrap || !svg) return;

  // Figure is always visible — empty rows render with em-dash placeholders
  // so the table scaffold (model names, r, SEE) shows even before any
  // inputs are entered. Gives the user a preview of what's coming.
  wrap.hidden = false;

  const ciLabel = document.getElementById('pre-forest-ci-pct');
  if (ciLabel) ciLabel.textContent = ciPct;

  // ── Geometry ──────────────────────────────────────────────────────
  // Modern JAMA-style tabular figure: numeric columns on the LEFT
  // (Model · r · SEE · Estimate (CI)), forest plot zone on the RIGHT.
  const W = 880;
  const padTop = 42;            // headers + rule + breath
  const padBottom = 52;         // axis ticks, labels, caption
  const rowHeight = 40;         // generous breathing room
  const plotHeight = rows.length * rowHeight;
  const H = padTop + plotHeight + padBottom;

  // Numeric column anchors — right-aligned x positions (except Model).
  // Spacing balanced for even ~45px gaps between visible header edges.
  const COL_MODEL_X       = 0;       // left-aligned, max ~270px
  const COL_R_X           = 312;     // right-aligned
  const COL_SEE_X         = 384;     // right-aligned (72 from R)
  const COL_ESTIMATE_X    = 528;     // right-aligned (144 from SEE — wider content)

  // Plot zone (right side). plotRight is pulled inward from the viewBox
  // edge so the centred "145" tick label has room to render without
  // being clipped by the SVG bounds.
  const plotLeft  = 564;
  const plotRight = W - 22;
  const plotWidth = plotRight - plotLeft;

  const xMin = 55, xMax = 145;
  const xScale = (ss) => plotLeft + ((Math.max(xMin, Math.min(xMax, ss)) - xMin) / (xMax - xMin)) * plotWidth;

  // Precision scaling: SEE ∈ [7.5, 10.5] → half-size ∈ [8, 5]
  // (diamond's half-width = half-height = `s`, so full diameter = 2s)
  function markerSize(see){
    if (!Number.isFinite(see) || see <= 0) return 6;
    const minS = 5, maxS = 8;
    const minSEE = 7.5, maxSEE = 10.5;
    const t = Math.max(0, Math.min(1, (see - minSEE) / (maxSEE - minSEE)));
    return maxS - t * (maxS - minS);
  }

  // ── Palette ───────────────────────────────────────────────────────
  // Soft charcoal grayscale throughout. No accent — the data itself is
  // the focal point.
  const INK      = '#1A1A1A';
  const INK_SOFT = '#2E2E2E';
  const MUTED    = '#6B6660';
  const FAINT    = '#9A968D';
  const RULE     = 'rgba(20,20,20,0.85)';
  const RULE_LN  = 0.6;
  const MEAN_LN  = 'rgba(20,20,20,0.18)';

  const SANS = '"IBM Plex Sans", "Inter", system-ui, -apple-system, "Helvetica Neue", Arial, sans-serif';
  // Helper: TXT(size, weight, fill, anchor, extra)
  const TXT = (size, weight, fill, anchor, extra) =>
    `font-family='${SANS}' font-size='${size}' font-weight='${weight || '400'}' fill='${fill}'`
    + (anchor ? ` text-anchor='${anchor}'` : '')
    + ` font-variant-numeric='tabular-nums' style='font-feature-settings:"tnum" 1'`
    + (extra ? ` ${extra}` : '');

  let out = '';

  // ── Top rule (subtle) ─────────────────────────────────────────────
  out += `<line x1="0" y1="2" x2="${W}" y2="2" stroke="${RULE}" stroke-width="${RULE_LN}"/>`;

  // ── Column headers (uppercase, letterspaced, weight 600, muted) ───
  const headerY = 22;
  const headerAttrs = `letter-spacing='0.10em'`;
  out += `<text x="${COL_MODEL_X}"            y="${headerY}" ${TXT('9.5', '600', MUTED, null, headerAttrs)}>MODEL</text>`;
  out += `<text x="${COL_R_X}"                y="${headerY}" ${TXT('9.5', '600', MUTED, 'end', headerAttrs)}>R</text>`;
  out += `<text x="${COL_SEE_X}"              y="${headerY}" ${TXT('9.5', '600', MUTED, 'end', headerAttrs)}>SEE</text>`;
  out += `<text x="${COL_ESTIMATE_X}"         y="${headerY}" ${TXT('9.5', '600', MUTED, 'end', headerAttrs)}>ESTIMATE (${ciPct} CI)</text>`;
  out += `<text x="${plotLeft + plotWidth/2}" y="${headerY}" ${TXT('9.5', '600', MUTED, 'middle', headerAttrs)}>FOREST PLOT</text>`;

  // Header bottom rule
  out += `<line x1="0" y1="${padTop - 10}" x2="${W}" y2="${padTop - 10}" stroke="${RULE}" stroke-width="${RULE_LN}"/>`;

  // ── Population-mean reference ─────────────────────────────────────
  const meanX = xScale(100);
  out += `<line x1="${meanX}" y1="${padTop - 2}" x2="${meanX}" y2="${padTop + plotHeight + 2}" stroke="${MEAN_LN}" stroke-width="0.7"/>`;

  // ── Data rows ─────────────────────────────────────────────────────
  rows.forEach((row, i) => {
    const y = padTop + rowHeight * (i + 0.5);
    const baseline = y + 4;

    // Model name
    const nameText = row.name.replace(/&/g, '&amp;');
    out += `<text x="${COL_MODEL_X}" y="${baseline}" ${TXT('12', '450', INK)}>${nameText}</text>`;

    // r
    if (row.r != null){
      out += `<text x="${COL_R_X}" y="${baseline}" ${TXT('12', '400', INK_SOFT, 'end')}>${row.r.toFixed(2)}</text>`;
    } else {
      out += `<text x="${COL_R_X}" y="${baseline}" ${TXT('12', '400', FAINT, 'end')}>—</text>`;
    }

    // SEE
    if (row.see != null){
      out += `<text x="${COL_SEE_X}" y="${baseline}" ${TXT('12', '400', INK_SOFT, 'end')}>${row.see.toFixed(2)}</text>`;
    } else {
      out += `<text x="${COL_SEE_X}" y="${baseline}" ${TXT('12', '400', FAINT, 'end')}>—</text>`;
    }

    if (row.val != null && row.see != null){
      const lo = row.val - mult * row.see;
      const hi = row.val + mult * row.see;
      const cx = xScale(row.val);
      const x1 = xScale(lo);
      const x2 = xScale(hi);

      // Estimate (CI) — point estimate in weight-700 charcoal so it
      // reads as the headline number; CI bounds in light muted grey
      // so the eye locks on the answer first.
      out += `<text x="${COL_ESTIMATE_X}" y="${baseline}" ${TXT('12.5', '700', INK, 'end')}>${Math.round(row.val)}`
          +  `<tspan font-weight='400' fill='${MUTED}' font-size='11.5'>  (${Math.round(lo)}–${Math.round(hi)})</tspan></text>`;

      // Whisker (charcoal, slightly thicker for cleaner read)
      out += `<line x1="${x1}" y1="${y}" x2="${x2}" y2="${y}" stroke="${INK}" stroke-width="1" stroke-linecap="round"/>`;
      // Small end caps
      out += `<line x1="${x1}" y1="${y - 3.5}" x2="${x1}" y2="${y + 3.5}" stroke="${INK}" stroke-width="1" stroke-linecap="round"/>`;
      out += `<line x1="${x2}" y1="${y - 3.5}" x2="${x2}" y2="${y + 3.5}" stroke="${INK}" stroke-width="1" stroke-linecap="round"/>`;

      // Filled DIAMOND marker — pointed corners make the point estimate
      // findable at a glance. Sized by precision. White stroke gives
      // crisp separation from the whisker behind it.
      const s = markerSize(row.see);
      const dia = `${cx},${y - s} ${cx + s},${y} ${cx},${y + s} ${cx - s},${y}`;
      out += `<polygon points="${dia}" fill="${INK}" stroke="#FFFFFF" stroke-width="1"/>`;
    } else {
      // Empty state — quiet em-dash. Plot zone stays empty for this row.
      out += `<text x="${COL_ESTIMATE_X}" y="${baseline}" ${TXT('12', '400', FAINT, 'end')}>—</text>`;
    }
  });

  // ── Data-block bottom rule ────────────────────────────────────────
  const dataBottomY = padTop + plotHeight + 6;
  out += `<line x1="0" y1="${dataBottomY}" x2="${W}" y2="${dataBottomY}" stroke="${RULE}" stroke-width="${RULE_LN}"/>`;

  // ── Axis ──────────────────────────────────────────────────────────
  const axisY = dataBottomY;
  const ticks = [55, 70, 85, 100, 115, 130, 145];
  ticks.forEach(t => {
    const x = xScale(t);
    const isMean = (t === 100);
    out += `<line x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY + (isMean ? 6 : 3.5)}" stroke="${INK}" stroke-width="0.7"/>`;
    out += `<text x="${x}" y="${axisY + 18}" ${TXT('10.5', isMean ? '500' : '400', isMean ? INK : MUTED, 'middle')}>${t}</text>`;
  });

  // Axis caption removed per request — the plot's axis labels speak for themselves.

  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = out;
}

// === Predicted vs Actual calculation (ToPF-based) ===
function calcPredict(){
  const topf = preNum('pre-topf');
  const sex  = preStr('pre-sex');
  const edu  = preNum('pre-edu');
  const age  = preNum('pre-age');
  // CI level comes from premorbidCi() — see its note; do not re-parse #pre-ci.
  // ToPF-predicted WAIS-IV indices use the WAIS-IV manual sex coding: F = 1, M = 2.
  const sexC_topf = sex === 'Female' ? 1 : sex === 'Male' ? 2 : null;
  const mult = preCiMult();
  const ciPct = premorbidCi().short;

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
  else if (diff < 0 && Number.isFinite(see) && see > 0){
    /* Past the end of the tabulated range. This used to assert "< 0.01%",
       which is false for seven of the eight columns — one row past the last
       entry the model gives PRI 0.081%, PSI 0.060%, DMI 0.048%, VWMI 0.038%,
       IMI 0.033%, VCI 0.019% and WMI 0.012%. Only FSIQ, with the narrowest
       SEE, actually falls under 0.01% there.
       BASE_RATES is itself round(Phi(d / SEE), 4) — see the note above it in
       data.js — so continuing the same model past the table is exact rather
       than an approximation, and keeps the column internally consistent.
       fmtPctBr rounds to 2dp, so genuinely tiny values would print "0.00%" —
       asserting zero, which is worse than the bound it replaced. Fall back to
       "< 0.01%" only where that is actually true. */
    const p = normCDF(diff / see);
    brEl.textContent = p < 0.0001 ? '< 0.01%' : fmtPctBr(p);
  }
  else brEl.textContent = '-';
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
    // Same rounding convention as the on-screen table (calcPremorbid) and the
    // OPIE-4 tab: round(estimate) ± round(z·SEE). Keeps the report text and the
    // screen in agreement — they previously differed by a point on ~21% of values.
    const margin = (r.see != null) ? Math.round(mult*r.see) : null;
    const lo = (margin != null) ? Math.round(r.val) - margin : '-';
    const hi = (margin != null) ? Math.round(r.val) + margin : '-';
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
    ${apaNoteHtml('pre-estimates', { ciMultiplier: ciPct === '95%' ? '1.96' : '1.645' })}
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
    ${apaNoteHtml('pre-predict')}
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

  const hasVC = vc != null, hasMR = mr != null;

  function predict(c){
    let pred = c.intercept + (c.age != null ? c.age * age : 0) + c.age3 * Math.pow(age, 3) + c.sex * sexC_opie;
    if (c.age6 != null) pred += c.age6 * Math.pow(age, 6);   // VCI only
    if (c.vc != null && hasVC) pred += c.vc * vc;
    if (c.mr != null && hasMR) pred += c.mr * mr;
    return pred;
  }

  const list = [];
  // Age must be inside the fitted 16-90 range (the Age3/Age6 terms are unbounded
  // and extrapolate hard), and sex must be present — Female is coded 0, so a
  // blank field previously produced the female equation silently.
  const canUseAge = age != null && age >= OPIE_AGE_MIN && age <= OPIE_AGE_MAX && sexC_opie != null;
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
  pushModel('GAI_MR', 'GAI - Matrix only', OPIE_PRORATED_GAI.MR, false, true, 'opiePredGAI_MR');

  // Single-index models (VCI from Vocabulary, PRI from Matrix Reasoning)
  pushModel('VCI', 'VCI - Vocab', OPIE_PRORATED_INDEX.VCI, true, false, 'opiePredVCI');
  pushModel('PRI', 'PRI - Matrix', OPIE_PRORATED_INDEX.PRI, false, true, 'opiePredPRI');

  return list;
}

function opieBaseRateFor(rowKey, diff){
  // NB `if (!diff)` treated a difference of exactly 0 as missing data and printed
  // '-', even though 0 is the single most likely outcome and -1 / +1 both resolve
  // to 'common'. Test for absence explicitly instead.
  if (diff == null || !Number.isFinite(diff)) return '-';
  const key = diff > 0 ? `+${diff}` : String(diff);
  const row = OPIE_BASE_RATES[key];
  if (row && row[rowKey] != null) return fmtPctBr(row[rowKey]);
  // Differences < 5 points are not tabulated in Table eA5.12 - always 'common'.
  if (Math.abs(diff) < 5) return 'common';
  // Beyond the range of the published table - for negative discrepancies this means
  // the base rate is smaller than the minimum tabulated value (< 0.1%)
  if (diff < 0) return '< 0.1%';
  return '-';
}

// ----- Shared achieved scores across premorbid methods -----
// The patient's obtained index scores are the same whichever premorbid model
// predicted them, so every achieved-score input (ToPF tab + OPIE-4 tab) reads
// and writes ONE store keyed by index (preState.achieved). Entering a score on
// one method autofills it on the others (FSIQ / VCI / PRI are shared; the three
// FSIQ rows and three GAI rows within the OPIE tab also share a single value).
function opieModelIndex(key){
  if (/^FSIQ_/.test(key)) return 'FSIQ';
  if (/^GAI_/.test(key))  return 'GAI';
  return key; // VCI, PRI
}
// Recompute the OPIE-4 table's difference / base-rate cells in place (no input
// rebuild, so focus is preserved), then refresh its APA output.
function refreshOpieDerived(){
  const tb = preGet('pre-opiepredict-tbody');
  if (tb){
    (preState.opieRows || []).forEach(row => {
      const inp = tb.querySelector(`[data-pre-opie-ach="${row.key}"]`);
      const tr = inp ? inp.closest('tr') : null;
      if (!tr) return;
      const diffCell = tr.querySelector('.diff');
      const brCell = tr.querySelector('.opie-base-rate');
      const raw = preState.achieved[opieModelIndex(row.key)];
      const achVal = (raw != null && raw !== '' && !isNaN(parseFloat(raw))) ? parseFloat(raw) : null;
      if (diffCell) diffCell.className = 'num diff';
      if (achVal == null || row.val == null || !Number.isFinite(row.val)){
        if (diffCell) diffCell.textContent = '-';
        if (brCell) brCell.textContent = '-';
      } else {
        const diff = Math.round(achVal - row.val);
        if (diffCell){ diffCell.textContent = (diff > 0 ? '+' : '') + diff; if (diff < 0) diffCell.classList.add('neg'); else if (diff > 0) diffCell.classList.add('pos'); }
        if (brCell) brCell.textContent = opieBaseRateFor(row.key, diff);
      }
    });
  }
  renderOpiePredictApa();
}
// Single entry point for any achieved-score change: store it, mirror it into
// every other input for that index, and recompute both methods' outputs.
function setSharedAchieved(idx, value, sourceEl){
  preState.achieved[idx] = value;
  document.querySelectorAll(`[data-pre-ach="${idx}"]`).forEach(el => { if (el !== sourceEl) el.value = value; });
  document.querySelectorAll('[data-pre-opie-ach]').forEach(el => {
    if (el !== sourceEl && opieModelIndex(el.dataset.preOpieAch) === idx) el.value = value;
  });
  if (typeof calcPredict === 'function') calcPredict();   // ToPF derived cells + APA
  refreshOpieDerived();                                   // OPIE derived cells + APA
}

function calcOpiePredict(){
  const tbody = preGet('pre-opiepredict-tbody');
  if (!tbody) return;

  // CI level comes from premorbidCi() — see its note; do not re-parse #pre-ci.
  const mult = preCiMult();
  const ciPct = premorbidCi().short;
  document.querySelectorAll('.pre-ci-lo-hdr-3').forEach(e => e.textContent = `Lower ${ciPct}`);
  document.querySelectorAll('.pre-ci-hi-hdr-3').forEach(e => e.textContent = `Upper ${ciPct}`);

  const rows = getOpiePredictions();
  preState.opieRows = rows;
  preState.ciPct = ciPct;

  // Group: FSIQ rows, then GAI rows, then single-index (VCI / PRI) rows
  const fsiq = rows.filter(r => r.key.startsWith('FSIQ_'));
  const gai  = rows.filter(r => r.key.startsWith('GAI_'));
  const idx  = rows.filter(r => r.key === 'VCI' || r.key === 'PRI');

  function modelCell(row){
    const tip = PRE_MODEL_TOOLTIPS[row.tipKey] || '';
    return `<td class="model model-has-tip">${escapeHtml(row.label)}<span class="model-info-dot" data-tooltip="${escapeAttr(tip)}" aria-label="${escapeAttr(row.label + '. ' + tip)}">?</span></td>`;
  }

  function rowHtml(row){
    const ach = preState.achieved[opieModelIndex(row.key)];
    const achVal = (ach != null && ach !== '') ? parseFloat(ach) : null;
    const hasPred = row.val != null && Number.isFinite(row.val);
    // Parametric prediction interval: round the point estimate and the margin
    // separately so the bounds are symmetric — round(predicted) ± round(z·SEE).
    const margin = Math.round(mult * row.see);
    const lo = hasPred ? String(Math.round(row.val) - margin) : '-';
    const hi = hasPred ? String(Math.round(row.val) + margin) : '-';
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
  if (idx.length){
    html += '<tr class="group-row"><td colspan="7">Index predictions (VCI / PRI)</td></tr>';
    html += idx.map(rowHtml).join('');
  }
  tbody.innerHTML = html;

  // Wire the achieved inputs (re-wired each render because rows are rebuilt).
  // Achieved scores are shared by index, so a change autofills the matching
  // rows here AND the ToPF tab via setSharedAchieved.
  tbody.querySelectorAll('[data-pre-opie-ach]').forEach(inp => {
    inp.addEventListener('input', e => {
      setSharedAchieved(opieModelIndex(e.target.dataset.preOpieAch), e.target.value, e.target);
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

  const byGroup = { FSIQ:[], GAI:[], INDEX:[] };
  rows.forEach(r => {
    if (r.key.startsWith('FSIQ_')) byGroup.FSIQ.push(r);
    else if (r.key.startsWith('GAI_')) byGroup.GAI.push(r);
    else if (r.key === 'VCI' || r.key === 'PRI') byGroup.INDEX.push(r);
  });
  const groupLabels = { FSIQ:'Prorated FSIQ', GAI:'Prorated GAI', INDEX:'Index (VCI / PRI)' };

  let body = '';
  ['FSIQ','GAI','INDEX'].forEach(group => {
    if (!byGroup[group].length) return;
    body += `<tr><td colspan="7" style="font-style:italic;padding-top:6px;padding-bottom:2px">${groupLabels[group]}</td></tr>`;
    byGroup[group].forEach(r => {
      const achRaw = preState.achieved[opieModelIndex(r.key)];
      const ach = achRaw != null && achRaw !== '' && !isNaN(parseFloat(achRaw)) ? parseFloat(achRaw) : null;
      const hasPred = r.val != null && Number.isFinite(r.val);
      // Parametric PI: round(predicted) ± round(z × SEE) — symmetric.
      const margin = Math.round(mult * r.see);
      const lo = hasPred ? Math.round(r.val) - margin : '-';
      const hi = hasPred ? Math.round(r.val) + margin : '-';
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
    ${apaNoteHtml('pre-opiepredict')}
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
    const example = sdiMode() === 'raw'
      ? {name:'Example memory score',t1:'42',t2:'36',sd:'8',isExample:true}
      : {name:'Example memory score',t1:'9',t2:'6',scoreType:'scaled',isExample:true};
    sdiRows.push(example); renderSdi(); showToast('Example row added'); return;
  }
  if (method === 'battery'){
    const example = {name:'Example subtest', raw:'25', score:'10', isExample:true};
    batteryRows.push(example); renderBattery(); showToast('Example row added'); return;
  }
  // Single superset example so the row is complete on every method (shared data).
  const rciExample = {name:'Example index score', isExample:true, sd:'15', m1:'100', sd1:'15', m2:'103', sd2:'15', r:'0.90', rCorrected:'', n:'100', t1:'100', t2:'89'};
  if (rciState[method]){ rciState[method].rows.push({...rciExample}); renderRci(method); showToast('Example row added'); }
}

/* ---------- INITIAL POPULATION ---------- */
// Add a starter example row (purely illustrative, excluded from analyses)
// and one blank row to each editable table.
batteryRows = [{name:'Example subtest', raw:'25', score:'10', isExample:true}, {name:'', raw:'', score:''}];
sdiRows = [
  sdiMode() === 'raw'
    ? {name:'Example memory score',t1:'42',t2:'36',sd:'8',isExample:true}
    : {name:'Example memory score',t1:'9',t2:'6',scoreType:'scaled',isExample:true},
  {}
];
// One shared example + blank row across all RCI methods (mutate in place to keep the shared link).
RCI_SHARED_ROWS.length = 0;
RCI_SHARED_ROWS.push({name:'Example index score', isExample:true, sd:'15', m1:'100', sd1:'15', m2:'103', sd2:'15', r:'0.90', rCorrected:'', n:'100', t1:'100', t2:'89'}, {});
['rci-basic', 'rci-practice', 'rci-srb', 'rci-crawford'].forEach(m => renderRci(m));
renderBattery();
renderSdi();
applyCalculatorPolish();
enhanceCalculatorWorkflow();
enhanceApaToolbars();

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

/* Premorbid → Score-Tables auto-link: picker button + popover wiring.
   Uses document-level delegation because the design-system rebuilds the
   ds-prem-card later, replacing the original button. Delegation keeps the
   handler attached regardless of when the markup arrives. */
document.addEventListener('click', e => {
  const btn = e.target.closest('#bat-prem-link-btn');
  if (!btn || btn.disabled) return;
  e.stopPropagation();
  const popover = document.getElementById('bat-prem-link-popover');
  if (popover && popover.classList.contains('show')) closePremorbidLinkPopover();
  else openPremorbidLinkPopover();
});

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

  /* Login/logout drops the user back to Home. Routed through navigateTo so this
     path can't drift from the main one; history is left alone because an auth
     transition isn't somewhere the back button should return to. */
  function activateHomePage(){
    navigateTo('home', { history: false });
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
    battery: 'Score Tables',
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
    'bat-apa':           'Score Tables',
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
    'bat-apa':              'Results',
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
  let dragGroupIds = null; // when dragging a merged-battery group, all member ids
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
    // Automatic "Saved you time today?" Ko-fi prompts are disabled by request.
    // The toast only appears when the user explicitly opens it (e.g. the header
    // "Buy me a coffee" button, which passes { force:true }).
    if (!force) return;
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
    // Deliberately no .apa-table-num fallback — the report strips table numbers,
    // so falling back to one would label the item "Table 1" for no reason.
    const titleEl = container.querySelector('.apa-table-title');
    return (titleEl?.textContent || fallback || 'APA Table').trim();
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

  /* Drop the captured "Table N" line. The working report is pasted into a host
     document that does its own numbering, so a number here is always wrong by
     the time it lands. Stripped at render rather than at capture, so reports
     saved before this change lose the number too (same approach as
     ensureBottomRule below). The tool pages keep their own numbers. */
  function stripTableNumber(html){
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    tmp.querySelectorAll('.apa-table-num').forEach(el => el.remove());
    return tmp.innerHTML;
  }
  /* Ensure every table has a bottom rule under its last row WITH content — even
     for items captured before this fix (their HTML snapshot is stored). Applied
     at render/copy time so existing saved reports are corrected too. */
  function ensureBottomRule(html){
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    tmp.querySelectorAll('.apa-table tbody').forEach(tbody => {
      const rows = [...tbody.querySelectorAll('tr')];
      let last = null;
      for (let k = rows.length - 1; k >= 0; k--){ if (rows[k].textContent.trim()){ last = rows[k]; break; } }
      if (!last) return;
      last.querySelectorAll('td, th').forEach(td => {
        const s = (td.getAttribute('style') || '').replace(/border-bottom\s*:[^;]*;?/gi, '');
        td.setAttribute('style', s + 'border-bottom:1.5pt solid #000;padding-bottom:3pt;');
      });
    });
    return tmp.innerHTML;
  }
  function effectiveItemHtml(item){
    let html = applyHiddenColumns(item.html, item.hiddenColumns || []);
    html = applyHeaderOverrides(html, item.headerOverrides || []);
    html = stripTableNumber(html);
    html = ensureBottomRule(html);
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
  // ---- Merge tables from the same battery into one combined table ----
  let mergeBattery = true; // default ON (toggle in the report toolbar)
  // Map a detected family name (e.g. "CVLT-3 Indices") to its battery in the
  // catalog → full name + the sub-section label ("Indices").
  function catalogBatteryFor(fam){
    if (!fam || typeof REPORT_TEST_CATALOG === 'undefined') return null;
    for (const t of REPORT_TEST_CATALOG){
      const fams = t.families || [];
      const base = fams.find(f => fam === f || fam.startsWith(f + ' ') || fam.startsWith(f));
      if (base || fam === t.name || fam.startsWith(t.name + ' ')){
        const stripFrom = base || t.name;
        let sub = fam.slice(stripFrom.length).replace(/^[\s·•:.\-]+/, '').trim();
        if (!sub) sub = fam.replace(t.name, '').replace(/^[\s·•:.\-]+/, '').trim();
        return { id: t.id, name: t.name, longName: t.longName || t.name, subLabel: sub || t.name };
      }
    }
    return null;
  }
  // Build one combined table from several same-battery tables, matching the
  // mockup: full battery-name title, each source table's header row kept as a
  // labelled sub-section, body rows beneath, one Note at the bottom.
  function buildMergedTableHtml(longName, sections){
    const wrap = document.createElement('div');
    wrap.innerHTML =
      `<div class="apa-table-title" style="font-family:'Times New Roman',serif;font-size:11pt;font-style:italic;color:#000;margin:0 0 8pt 0;line-height:1.4;">${longName}</div>`;
    const table = document.createElement('table');
    table.className = 'apa-table';
    table.setAttribute('style', "border-collapse:collapse;font-family:'Times New Roman',serif;font-size:11pt;color:#000;width:auto;");
    let note = null;

    // Parse each section's source table once.
    const parsed = sections.map(sec => {
      const tmp = document.createElement('div');
      tmp.innerHTML = sec.html;
      const src = tmp.querySelector('.apa-table');
      const headRows = src ? [...src.querySelectorAll('thead tr')] : [];
      const colRow = headRows.length ? headRows[headRows.length - 1] : null;
      const bodyRows = src ? [...src.querySelectorAll('tbody tr')].filter(r => !r.classList.contains('apa-group')) : [];
      const noteEl = tmp.querySelector('.apa-note');
      return { sec, src, headRows, colRow, bodyRows, noteEl, colCount: colRow ? colRow.children.length : 0 };
    }).filter(p => p.src);
    /* Merge the sections' notes instead of keeping only the first. Sections of
       one battery usually share a note, so dedupe on text — but where they
       differ (a section using raw r, say) the caveat has to survive, and the
       old first-note-wins rule silently dropped it. */
    {
      const seen = new Set();
      const bodies = [];
      parsed.forEach(p => {
        if (!p.noteEl) return;
        const key = (p.noteEl.textContent || '').replace(/\s+/g, ' ').trim();
        if (!key || seen.has(key)) return;
        seen.add(key);
        const clone = p.noteEl.cloneNode(true);
        // Drop the leading "Note." label from all but the first — one label
        // introduces the combined note.
        if (bodies.length){
          const label = clone.querySelector('strong');
          if (label) label.remove();
        }
        bodies.push(clone);
      });
      if (bodies.length){
        note = bodies[0];
        bodies.slice(1).forEach(extra => {
          note.appendChild(document.createTextNode(' '));
          while (extra.firstChild) note.appendChild(extra.firstChild);
        });
      }
    }

    // A single shared column-header row only works when every section has the
    // same columns; otherwise fall back to per-section headers (legacy) to avoid
    // misalignment.
    const counts = parsed.map(p => p.colCount);
    const uniform = parsed.length > 0 && counts.every(c => c > 0 && c === counts[0]);
    const tbody = document.createElement('tbody');

    if (uniform){
      // One column-header row at the top (labels from the first section; its
      // first cell already reads "Subtest"). Each sub-section becomes a
      // full-width bold divider row, with data rows beneath.
      const thead = document.createElement('thead');
      thead.appendChild(parsed[0].colRow.cloneNode(true));
      table.appendChild(thead);
      const colspan = counts[0];
      parsed.forEach((p, si) => {
        if (p.sec.subLabel){
          const dr = document.createElement('tr');
          dr.className = 'apa-group';
          const td = document.createElement('td');
          td.setAttribute('colspan', String(colspan));
          td.setAttribute('style', "font-family:'Times New Roman',serif;font-size:11pt;font-weight:bold;color:#000;text-align:left;padding:" + (si > 0 ? '10pt' : '4pt') + " 0 2pt;");
          td.textContent = p.sec.subLabel;
          dr.appendChild(td);
          tbody.appendChild(dr);
        }
        p.bodyRows.forEach(r => tbody.appendChild(r.cloneNode(true)));
      });
    } else {
      // Fallback: per-section header rows, sub-label in the first header cell.
      parsed.forEach((p, si) => {
        p.headRows.forEach((hr, hi) => {
          const tr = hr.cloneNode(true);
          if (hi === 0 && p.sec.subLabel){
            const fc = tr.querySelector('th, td');
            if (fc) fc.textContent = p.sec.subLabel;
          }
          if (hi === 0 && si > 0){
            tr.querySelectorAll('th, td').forEach(c => {
              const s = (c.getAttribute('style') || '').replace(/border-top\s*:[^;]*;?/gi, '');
              c.setAttribute('style', s + 'border-top:1.5pt solid #000;');
            });
          }
          tbody.appendChild(tr);
        });
        p.bodyRows.forEach(r => tbody.appendChild(r.cloneNode(true)));
      });
    }
    table.appendChild(tbody);
    // Strip per-section bottom rules inherited from the source tables (these
    // showed as lines between merged sections); only the final row should carry
    // the closing rule, added below.
    table.querySelectorAll('tbody td, tbody th').forEach(c => {
      const s = (c.getAttribute('style') || '').replace(/border-bottom\s*:[^;]*;?/gi, '');
      c.setAttribute('style', s);
    });
    // Bottom rule under the last row with content
    const rows = [...table.querySelectorAll('tbody tr')];
    for (let k = rows.length - 1; k >= 0; k--){
      if (rows[k].textContent.trim()){
        rows[k].querySelectorAll('td, th').forEach(c => {
          const s = (c.getAttribute('style') || '').replace(/border-bottom\s*:[^;]*;?/gi, '');
          c.setAttribute('style', s + 'border-bottom:1.5pt solid #000;padding-bottom:3pt;');
        });
        break;
      }
    }
    wrap.appendChild(table);
    if (note) wrap.appendChild(note);
    return wrap.innerHTML;
  }
  // Group report items by battery; merge groups of 2+, pass singles through.
  /* Group state.items into ordered render blocks. When merge is on, same-battery
     items (2+) collapse into one block. Each block carries its member ids so the
     edit view can attach group-level controls (reorder / remove / copy). The
     html is processed (hidden cols, overrides, bottom rule, merge) but NOT yet
     renumbered — callers renumber by block position. */
  // The label shown on each merged section header (e.g. "Indices", "Core
  // Subtests"). Prefer the catalog sub-label derived from the table's own group
  // row; but split items keep only data rows (no group row), so detectTestFamily
  // collapses to the bare battery name — in that case fall back to the item's
  // sourceTool/title, which still carries the sub-section, and strip the battery
  // name + age suffix off it.
  function deriveSubLabel(it, info){
    if (info && info.subLabel && info.subLabel !== info.name) return info.subLabel;
    const raw = stripAgeRange(((it && (it.sourceTool || it.title)) || '')).trim();
    if (info && info.name && raw){
      const esc = info.name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      const stripped = raw.replace(new RegExp('^' + esc + '\\b'), '').replace(/^[\s·•:.\-]+/, '').trim();
      if (stripped) return stripped;
    }
    return (info && info.subLabel) || raw || (info ? info.name : '');
  }
  // Default sub-section order within a merged battery: Indices → Core →
  // Supplementary → Process → (anything else). CVLT trials sit just after the
  // indices. Lower number = earlier.
  function subSectionRank(subLabel){
    const s = (subLabel || '').toLowerCase();
    if (/index|indices|composite/.test(s)) return 0;
    if (/core/.test(s))                    return 1;
    if (/trial/.test(s))                   return 1;
    if (/supplement/.test(s))              return 2;
    if (/process/.test(s))                 return 3;
    return 50;
  }
  function computeBlocks(items){
    const list = items || state.items;
    if (!(mergeBattery && list.length > 1)){
      return list.map(it => ({ ids: [it.id], items: [it], isMerged: false, longName: null, html: effectiveItemHtml(it) }));
    }
    const groups = [];
    const byKey = new Map();
    list.forEach(it => {
      const processed = effectiveItemHtml(it);
      const info = catalogBatteryFor(detectTestFamily(processed));
      const key = info ? info.id : null;
      const member = { id: it.id, item: it, html: processed, subLabel: deriveSubLabel(it, info) };
      if (key && byKey.has(key)){
        byKey.get(key).members.push(member);
      } else {
        const g = { key, longName: info ? info.longName : null, members: [member] };
        groups.push(g);
        if (key) byKey.set(key, g);
      }
    });
    return groups.map(g => {
      // Order members: a user-pinned mergeRank wins; otherwise the canonical
      // sub-section order. Stable on insertion order as the tiebreak.
      const ordered = g.members
        .map((m, i) => ({ m, i }))
        .sort((a, b) => {
          const ka = Number.isFinite(a.m.item.mergeRank) ? a.m.item.mergeRank : subSectionRank(a.m.subLabel);
          const kb = Number.isFinite(b.m.item.mergeRank) ? b.m.item.mergeRank : subSectionRank(b.m.subLabel);
          return ka !== kb ? ka - kb : a.i - b.i;
        })
        .map(x => x.m);
      const isMerged = ordered.length > 1;
      return {
        key: g.key,
        ids: ordered.map(m => m.id),
        items: ordered.map(m => m.item),
        isMerged,
        longName: g.longName,
        html: isMerged
          ? buildMergedTableHtml(g.longName, ordered.map(m => ({ html: m.html, subLabel: m.subLabel })))
          : ordered[0].html
      };
    });
  }
  function mergeReportBlocks(items){
    return computeBlocks(items).map(b => b.html);
  }
  function buildReportHtmlBody(){
    if (!state.items.length) return '<p style="color:#888;font-style:italic;">No items in the working report yet.</p>';
    // Word ignores margins on <div>s between pasted tables, so insert blank
    // paragraphs between tables — that is the reliable way to add a visible gap.
    const spacer = '<p style="margin:0;font-family:\'Times New Roman\',serif;font-size:12pt;line-height:12pt;">&nbsp;</p>';
    return mergeReportBlocks(state.items)
      .map(html => `<div style="page-break-inside:avoid;">${html}</div>`)
      .join(spacer + spacer);
  }
  async function copyAll(){
    if (!state.items.length){
      if (typeof showToast === 'function') showToast('Working report is empty', true);
      return;
    }
    const html = `<div style="font-family:'Times New Roman',serif;font-size:11pt;color:#000;">${buildReportHtmlBody()}</div>`;
    // Build the plain-text alternative from the SAME blocks as the HTML, so the
    // two clipboard flavours agree on table count when battery merging is on.
    const blocks = mergeReportBlocks(state.items);
    const plain = blocks.map(blockHtml => {
      const tmp = document.createElement('div');
      tmp.innerHTML = blockHtml;
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
      if (typeof showToast === 'function') showToast(`✓ ${blocks.length} table${blocks.length===1?'':'s'} copied`);
      if (typeof flashCopiedButton === 'function') flashCopiedButton(rootEl && rootEl.querySelector('[data-rb-action="copy"]'));
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
    // Same blocks as the HTML/Word output, so a merged battery exports as one
    // section here too rather than as its individual source tables.
    const sections = mergeReportBlocks(state.items).map((blockHtml, i) => {
      const tmp = document.createElement('div');
      tmp.innerHTML = blockHtml;
      const table = tmp.querySelector('table');
      if (!table) return '';
      const heading = (tmp.querySelector('.apa-table-title')?.textContent || '').trim()
                   || `Table ${i + 1}`;
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
            <button class="btn rb-action-merge" data-rb-action="toggle-merge" type="button" title="Combine tables from the same test battery into one (e.g. all CVLT-3 tables → one table)">
              Merge by battery: On
            </button>
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
        else if (action === 'toggle-merge') {
          mergeBattery = !mergeBattery;
          rootEl.querySelectorAll('[data-rb-action="toggle-merge"]').forEach(b => { b.textContent = 'Merge by battery: ' + (mergeBattery ? 'On' : 'Off'); });
          render();
        }
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
      // Block-level reorder (single item OR merged group) from the Options menu
      const blockAct = e.target.closest('[data-rb-block-action]');
      if (blockAct){
        e.preventDefault();
        e.stopPropagation();
        moveBlockBy(blockAct.dataset.rbAnchor, blockAct.dataset.rbBlockAction);
        rootEl.querySelectorAll('.rb-item-menu.is-open').forEach(m => m.classList.remove('is-open'));
        return;
      }
      // Reorder a table WITHIN its merged group; keep the menu open to chain moves
      const memberAct = e.target.closest('[data-rb-member-action]');
      if (memberAct){
        e.preventDefault();
        e.stopPropagation();
        const mid = memberAct.dataset.rbMemberId;
        moveWithinGroup(mid, memberAct.dataset.rbMemberAction);
        requestAnimationFrame(() => {
          const cards = rootEl?.querySelectorAll('.rb-item-group') || [];
          for (const card of cards){
            if ((card.dataset.rbGroupIds || '').split(',').includes(mid)){
              card.querySelector('.rb-item-menu')?.classList.add('is-open');
              break;
            }
          }
        });
        return;
      }
      // Copy a whole merged group (the combined table)
      const groupCopy = e.target.closest('[data-rb-group-copy]');
      if (groupCopy){
        e.preventDefault();
        e.stopPropagation();
        copyGroup(groupCopy.dataset.rbGroupCopy.split(','));
        return;
      }
      // Remove a whole merged group (× or "Remove all")
      const groupRemove = e.target.closest('[data-rb-group-remove]');
      if (groupRemove){
        e.preventDefault();
        e.stopPropagation();
        removeMany(groupRemove.dataset.rbGroupRemove.split(','));
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
      dragGroupIds = item.dataset.rbGroupIds ? item.dataset.rbGroupIds.split(',') : [dragId];
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
      dragGroupIds = null;
    });
    body.addEventListener('dragover', e => {
      const target = e.target.closest('.rb-item');
      if (!target || !dragId || (dragGroupIds || [dragId]).includes(target.dataset.rbId)) return;
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
      const fromIds = dragGroupIds || (dragId ? [dragId] : null);
      if (!target || !fromIds || fromIds.includes(target.dataset.rbId)){ dragId = null; dragGroupIds = null; return; }
      e.preventDefault();
      const rect = target.getBoundingClientRect();
      const after = (e.clientY - rect.top) > rect.height / 2;
      dragId = null;
      dragGroupIds = null;
      moveBlock(fromIds, target.dataset.rbId, after);
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
      // Block-aware: moves the whole merged group when the card is a group.
      moveBlockBy(id, e.key === 'ArrowUp' ? 'up' : 'down');
      lastChangedItemId = id;
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
    // Re-render the body so it switches between per-item edit cards and the
    // merged, export-faithful preview document.
    renderItems();
    if (!previewMode) decorateEditableHeaders();
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
    const itemHtml = effectiveItemHtml(item);
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
      if (typeof flashCopiedButton === 'function') flashCopiedButton(rootEl && rootEl.querySelector(`[data-rb-item-copy="${id}"]`));
      maybeShowKofiToast();
    } catch(e){
      console.error(e);
      if (typeof showToast === 'function') showToast('Copy failed - try selecting manually', true);
    }
  }

  /* ---------- merged-group (Edit view) operations ----------
     A "block" is one render unit: a single item, or a merged battery group of
     several items. These let group cards reorder / remove / copy as a unit. */
  function blockIdsFor(id){
    const b = computeBlocks().find(bl => bl.ids.includes(id));
    return b ? b.ids : [id];
  }
  // Move a set of item ids so they sit (contiguously) before/after the target
  // block. Used by both group drag and single-item drag, so a drop never lands
  // in the middle of a merged group.
  function moveBlock(fromIds, toId, dropAfter){
    const idSet = new Set(fromIds);
    if (idSet.has(toId)) return;
    const targetIds = blockIdsFor(toId);
    const moved = state.items.filter(i => idSet.has(i.id));
    if (!moved.length) return;
    const rest = state.items.filter(i => !idSet.has(i.id));
    const targetIdxs = targetIds.map(id => rest.findIndex(i => i.id === id)).filter(x => x >= 0);
    if (!targetIdxs.length){ state.items = rest.concat(moved); save(); render(); return; }
    const insertAt = dropAfter ? Math.max(...targetIdxs) + 1 : Math.min(...targetIdxs);
    rest.splice(insertAt, 0, ...moved);
    state.items = rest;
    save();
    render();
  }
  // Reorder one member up/down WITHIN its merged group (swaps it with the
  // adjacent group sibling), so e.g. the indices table can sit above core.
  function moveWithinGroup(memberId, dir){
    const block = computeBlocks().find(b => b.ids.includes(memberId));
    if (!block) return;
    const order = block.ids.slice();                 // current displayed order
    const pos = order.indexOf(memberId);
    const swapPos = dir === 'up' ? pos - 1 : pos + 1;
    if (swapPos < 0 || swapPos >= order.length) return;
    [order[pos], order[swapPos]] = [order[swapPos], order[pos]];
    // Pin the new within-group order via mergeRank so it survives re-render and
    // overrides the canonical default.
    order.forEach((id, idx) => {
      const it = state.items.find(i => i.id === id);
      if (it) it.mergeRank = idx;
    });
    save();
    render();
  }
  // Move the block that contains anchorId up/down one block via the Options menu.
  function moveBlockBy(anchorId, dir){
    const blocks = computeBlocks();
    const idx = blocks.findIndex(b => b.ids.includes(anchorId));
    if (idx < 0) return;
    const swapWith = dir === 'up' ? idx - 1 : idx + 1;
    if (swapWith < 0 || swapWith >= blocks.length) return;
    moveBlock(blocks[idx].ids, blocks[swapWith].ids[0], dir === 'down');
  }
  function removeMany(ids){
    const set = new Set(ids);
    state.items = state.items.filter(i => !set.has(i.id));
    save();
    render();
  }
  async function copyGroup(ids){
    const idKey = ids.join(',');
    const block = computeBlocks().find(b => b.ids.join(',') === idKey)
               || computeBlocks().find(b => ids.every(id => b.ids.includes(id)));
    if (!block) return;
    const inner = block.html;
    const html = `<div style="font-family:'Times New Roman',serif;font-size:11pt;color:#000;">${inner}</div>`;
    const tmp = document.createElement('div');
    tmp.innerHTML = inner;
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
      if (typeof showToast === 'function') showToast('✓ Merged table copied to clipboard');
      if (typeof flashCopiedButton === 'function') flashCopiedButton(rootEl && rootEl.querySelector(`[data-rb-group-copy="${idKey}"]`));
      maybeShowKofiToast();
    } catch(e){
      console.error(e);
      if (typeof showToast === 'function') showToast('Copy failed - try selecting manually', true);
    }
  }

  /* ---------- rendering ---------- */
  const RB_SVG = {
    grip: '<svg class="rb-item-grip-icon" viewBox="0 0 14 18" fill="currentColor" aria-hidden="true"><circle cx="5" cy="4" r="1.4"/><circle cx="9" cy="4" r="1.4"/><circle cx="5" cy="9" r="1.4"/><circle cx="9" cy="9" r="1.4"/><circle cx="5" cy="14" r="1.4"/><circle cx="9" cy="14" r="1.4"/></svg>',
    copy: '<svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="4" y="4" width="8" height="8" rx="1"/><path d="M2 9V3a1 1 0 0 1 1-1h6"/></svg>',
    options: '<svg viewBox="0 0 14 14" fill="currentColor" aria-hidden="true"><circle cx="3" cy="7" r="1.3"/><circle cx="7" cy="7" r="1.3"/><circle cx="11" cy="7" r="1.3"/></svg>',
    up: '<svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M3 6l4-4 4 4"/><path d="M7 2v10"/></svg>',
    down: '<svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M7 2v10"/><path d="M3 8l4 4 4-4"/></svg>',
    refresh: '<svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M2 7a5 5 0 0 1 9-3"/><path d="M11 1v3.5h-3"/><path d="M12 7a5 5 0 0 1-9 3"/><path d="M3 13V9.5h3"/></svg>',
    trash: '<svg viewBox="0 0 14 14" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M2 4h10"/><path d="M5.5 4V2.5h3V4"/><path d="M3 4l1 8h6l1-8"/></svg>',
    close: '<svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" aria-hidden="true"><path d="M3 3l6 6M9 3l-6 6"/></svg>'
  };
  const GRIP_ATTRS = 'draggable="true" title="Drag to reorder (or focus + Alt+↑/↓)" role="button" tabindex="0" aria-label="Drag to reorder. Keyboard: Alt + Arrow Up or Arrow Down."';

  /* A single (non-merged) report item card — full per-item controls. */
  function itemCardHtml(it, num, isFirst, isLast){
    const id = escapeHtmlLocal(it.id);
    const cols = getItemColumns(it);
    const hidden = new Set(it.hiddenColumns || []);
    const colsHtml = cols.length ? `
      <div class="rb-item-menu-section">Columns</div>
      ${cols.map(c => `
        <button class="rb-item-menu-item rb-col-toggle" type="button" role="menuitemcheckbox"
          aria-checked="${!hidden.has(c.idx)}"
          data-rb-col-toggle="${c.idx}"
          data-rb-item-id="${id}">
          <span class="rb-col-check ${hidden.has(c.idx) ? '' : 'is-checked'}" aria-hidden="true">
            <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6.5l2 2 4-5"/></svg>
          </span>
          <span class="rb-col-label">${escapeHtmlLocal(c.label)}</span>
        </button>
      `).join('')}
      <div class="rb-item-menu-sep"></div>
    ` : '';
    return `
    <article class="rb-item" data-rb-id="${id}">
      <header class="rb-item-header">
        <span class="rb-item-grip" ${GRIP_ATTRS}>${RB_SVG.grip}</span>
        <span class="rb-item-num">${num}</span>
        <div class="rb-item-meta">
          <span class="rb-item-source">${escapeHtmlLocal(it.sourceTool)}</span>
          <span class="rb-item-time" data-rb-time="${escapeHtmlLocal(it.updatedAt || it.addedAt)}">${formatRelative(it.updatedAt || it.addedAt)}</span>
        </div>
        <div class="rb-item-actions">
          <button class="rb-item-actionbtn rb-item-copy" type="button" data-rb-item-copy="${id}" aria-label="Copy this table" title="Copy this table to clipboard">${RB_SVG.copy}<span>Copy</span></button>
          <div class="rb-item-options-wrap">
            <button class="rb-item-actionbtn rb-item-options" type="button" data-rb-item-options="${id}" aria-label="More options" title="More options">${RB_SVG.options}<span>Options</span></button>
            <div class="rb-item-menu" data-rb-item-menu="${id}" role="menu">
              ${colsHtml}
              <button class="rb-item-menu-item" data-rb-block-action="up" data-rb-anchor="${id}" type="button" role="menuitem"${isFirst ? ' disabled' : ''}>${RB_SVG.up} Move up</button>
              <button class="rb-item-menu-item" data-rb-block-action="down" data-rb-anchor="${id}" type="button" role="menuitem"${isLast ? ' disabled' : ''}>${RB_SVG.down} Move down</button>
              <button class="rb-item-menu-item" data-rb-item-action="refresh" data-rb-item-id="${id}" type="button" role="menuitem">${RB_SVG.refresh} Refresh from source</button>
              <div class="rb-item-menu-sep"></div>
              <button class="rb-item-menu-item rb-item-menu-danger" data-rb-item-action="remove" data-rb-item-id="${id}" type="button" role="menuitem">${RB_SVG.trash} Remove</button>
            </div>
          </div>
          <button class="rb-item-remove" type="button" data-rb-remove="${id}" aria-label="Remove this table from the report" title="Remove from report">${RB_SVG.close}</button>
        </div>
      </header>
      <div class="rb-item-rendered">${effectiveItemHtml(it)}</div>
    </article>`;
  }

  /* A merged battery group card — group-level controls (one combined table,
     reorder/copy/remove the group, plus per-member remove). */
  function groupCardHtml(block, num, isFirst, isLast){
    const idKey = escapeHtmlLocal(block.ids.join(','));
    const anchor = escapeHtmlLocal(block.ids[0]);
    // Stable menu id keyed to the battery, so it survives within-group reorders
    // (the first member — and thus the anchor — can change as rows move).
    const menuId = escapeHtmlLocal('grp-' + (block.key || block.ids[0]));
    const title = escapeHtmlLocal(block.longName || 'Merged tables');
    const count = block.items.length;
    const members = block.items.map((it, mi) => {
      const mid = escapeHtmlLocal(it.id);
      const label = escapeHtmlLocal(it.sourceTool);
      return `
      <div class="rb-member-row">
        <span class="rb-member-label" title="${label}">${label}</span>
        <span class="rb-member-actions">
          <button class="rb-member-btn" data-rb-member-action="up" data-rb-member-id="${mid}" type="button" title="Move up within group" aria-label="Move ${label} up within group"${mi === 0 ? ' disabled' : ''}>${RB_SVG.up}</button>
          <button class="rb-member-btn" data-rb-member-action="down" data-rb-member-id="${mid}" type="button" title="Move down within group" aria-label="Move ${label} down within group"${mi === count - 1 ? ' disabled' : ''}>${RB_SVG.down}</button>
          <button class="rb-member-btn rb-member-btn-danger" data-rb-item-action="remove" data-rb-item-id="${mid}" type="button" title="Remove from report" aria-label="Remove ${label}">${RB_SVG.close}</button>
        </span>
      </div>`;
    }).join('');
    return `
    <article class="rb-item rb-item-group" data-rb-id="${anchor}" data-rb-group-ids="${idKey}">
      <header class="rb-item-header">
        <span class="rb-item-grip" ${GRIP_ATTRS}>${RB_SVG.grip}</span>
        <span class="rb-item-num">${num}</span>
        <div class="rb-item-meta">
          <span class="rb-item-source">${title}<span class="rb-merge-badge">Merged</span></span>
          <span class="rb-item-time">${count} tables combined</span>
        </div>
        <div class="rb-item-actions">
          <button class="rb-item-actionbtn rb-item-copy" type="button" data-rb-group-copy="${idKey}" aria-label="Copy the merged table" title="Copy the merged table to clipboard">${RB_SVG.copy}<span>Copy</span></button>
          <div class="rb-item-options-wrap">
            <button class="rb-item-actionbtn rb-item-options" type="button" data-rb-item-options="${menuId}" aria-label="More options" title="More options">${RB_SVG.options}<span>Options</span></button>
            <div class="rb-item-menu rb-item-menu-group" data-rb-item-menu="${menuId}" role="menu">
              <button class="rb-item-menu-item" data-rb-block-action="up" data-rb-anchor="${anchor}" type="button" role="menuitem"${isFirst ? ' disabled' : ''}>${RB_SVG.up} Move group up</button>
              <button class="rb-item-menu-item" data-rb-block-action="down" data-rb-anchor="${anchor}" type="button" role="menuitem"${isLast ? ' disabled' : ''}>${RB_SVG.down} Move group down</button>
              <div class="rb-item-menu-sep"></div>
              <div class="rb-item-menu-section">Order tables within this group</div>
              ${members}
              <div class="rb-item-menu-sep"></div>
              <button class="rb-item-menu-item rb-item-menu-danger" data-rb-group-remove="${idKey}" type="button" role="menuitem">${RB_SVG.trash} Remove all (${count})</button>
            </div>
          </div>
          <button class="rb-item-remove" type="button" data-rb-group-remove="${idKey}" aria-label="Remove the whole merged group" title="Remove the whole merged group">${RB_SVG.close}</button>
        </div>
      </header>
      <div class="rb-item-rendered">${block.html}</div>
    </article>`;
  }

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

    const blocks = computeBlocks();

    /* Preview mode renders the EXACT document that gets copied/exported —
       battery-merging applied when the toggle is on — as Word-style pages,
       with no editing chrome. */
    if (previewMode){
      body.innerHTML = blocks.map((b, i) =>
        `<article class="rb-item rb-item-preview"><div class="rb-item-rendered">${b.html}</div></article>`
      ).join('');
      return;
    }

    /* Edit mode: one card per block. Single items keep full per-item controls;
       merged battery groups render the combined table with group-level controls
       (reorder / copy / remove the group, plus per-member remove). */
    body.innerHTML = blocks.map((b, i) =>
      b.isMerged
        ? groupCardHtml(b, i + 1, i === 0, i === blocks.length - 1)
        : itemCardHtml(b.items[0], i + 1, i === 0, i === blocks.length - 1)
    ).join('');
  }

  let rbPrevCount = 0;
  function render(){
    if (!rootEl) return;
    const drawer = rootEl.querySelector('.rb-drawer');
    const onboarding = rootEl.querySelector('[data-rb-onboarding]');
    const countNow = state.items.length;
    rootEl.querySelectorAll('[data-rb-count]').forEach(el => {
      el.textContent = String(countNow);
      // Pop the badge when the count goes UP (a new table was captured).
      if (countNow > rbPrevCount){
        el.classList.remove('rb-count-bump'); void el.offsetWidth; el.classList.add('rb-count-bump');
        setTimeout(() => el.classList.remove('rb-count-bump'), 500);
      }
    });
    rbPrevCount = countNow;
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
