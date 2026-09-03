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
/* 'raw' is a deliberate no-conversion case, not an oversight. A raw score
   carries no metric, so there is no z, no percentile and no classification to
   be had from the number alone — see the metric:'raw' note in data.js. Both
   directions return null so every derived cell downstream stays blank instead
   of showing a figure derived from the wrong scale. Do not give 'raw' a
   fallback: falling through to 'scaled' is the exact bug this closes. */
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
    case 'raw': return null;
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
    case 'raw': return null;
    default: return null;
  }
}
/* Every score type the app understands. 'percentile' is deliberately absent:
   it is an OUTPUT of the converter, not a metric a stored measure sits on, and
   it is not offered by any row-type selector. Used to validate a declared
   metric on a database entry, so an unrecognised value is rejected at the door
   rather than silently ignored downstream. */
const SCORE_METRICS = new Set(['z', 't', 'scaled', 'standard', 'raw']);

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

/* The Norms Database page used to be Custom Tests, where a clinician could type
   new normative entries straight into the scoring database. It was gated behind
   a window.prompt for that reason.

   The typed form is gone and the page is now a read-only view of published
   norms, so there is nothing left to protect: the coefficients in it are
   printed in the manuals they are cited from. The gate also never provided any
   real protection, the password having been a literal in the shipped bundle.

   It had a concrete cost too. window.prompt() is ignored in a sandboxed iframe,
   so the call returned null, navigateTo returned false, and the nav item did
   nothing at all with no error — which is how this was found.

   Import is the only thing that still writes to the database, and it validates
   every entry through ctValidateEntry. */

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
    /* The age field is in the top bar and outlives the section swap, but
       whether anything on the new page reads it does not — so the pip is
       re-evaluated per page. Defined later in the file; this runs on click,
       long after init. */
    if (typeof refreshPatientAgeIndicator === 'function') refreshPatientAgeIndicator();
    if (typeof refreshBatteryAgePrompt === 'function') refreshBatteryAgePrompt();
    /* The popover is anchored to a control on Score Tables, so it cannot
       survive a page change — its anchor goes with it. */
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
function showToast(msg, isError, ms){
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = 'toast show' + (isError ? ' error' : '');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => t.className = 'toast', ms || 2200);
}

/* ---------- FIRST-EXPORT NOTICE ----------------------------------------
   Copy is the boundary where a number crosses from this app into a
   medico-legal document, which makes it the right moment to say that there
   are terms. It is NOT the right moment to state them: at the instant of
   copying, "make sure you are qualified" is un-actionable, the clinician
   having already chosen the method and entered the scores.

   So this fires ONCE, ever, and then never again. A per-copy notice on the
   app's core output path would be dismissed reflexively within a week, and
   a notice everyone dismisses is worth little to the user and probably
   little as evidence that anyone was informed. Once is also why it can be
   non-blocking: it never gates an export.

   Every copy path routes through here rather than calling showToast
   directly, so a new export route cannot silently miss it. */
const EXPORT_NOTICE_KEY = 'paExportNoticeSeen';
function exportToast(msg){
  let first = false;
  try { first = localStorage.getItem(EXPORT_NOTICE_KEY) !== 'true'; }
  catch(e){ first = false; }   // private mode: confirm the copy, skip the notice
  if (!first){ showToast(msg); return; }
  try { localStorage.setItem(EXPORT_NOTICE_KEY, 'true'); } catch(e){}
  showToast(msg + '  Interpretation and verification rest with you - see Privacy & use in the footer.', false, 9000);
  /* The toast is white-space:nowrap, sized for "Table copied". This one is a
     sentence, and on a narrow window a nowrap toast that long runs off the
     screen. Added after showToast, which assigns className wholesale — which
     is also what clears it again on the next ordinary toast. */
  const t = document.getElementById('toast');
  if (t) t.classList.add('is-longform');
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
    exportToast('✓ Table copied - ready to paste into your report');
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
  /* A section that relabels one column carries that label on its group row.
     The rule above has just overwritten it as a bold, left-aligned section
     name — but it is a COLUMN heading, so it has to sit over its own column
     exactly as that column's data cells do, and not read as part of the
     section name. Same 7pt horizontal padding, and centred wherever the column
     is numeric, which is how td.num is styled above. */
  clone.querySelectorAll('.apa-table tbody tr.apa-group td.apa-group-col-label').forEach(td => {
    const centred = td.classList.contains('num') ? 'text-align:center;' : '';
    td.setAttribute('style',
      "padding:5pt 7pt 2.5pt;border:none;font-family:'Times New Roman',serif;color:#333;"
      + "font-style:italic;font-weight:normal;font-size:9.5pt;line-height:1.05;" + centred);
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
  if (val === '' || isNaN(val)){ out.style.display = 'none'; renderConvFullTable(); return; }
  out.style.display = 'block';
  const z = toZ(val, type);
  if (z == null){ out.style.display = 'none'; renderConvFullTable(); return; }
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
  renderConvFullTable();
}

/* ---------- Full conversion table view ----------
   A scannable lookup table, the appendix every manual carries. Every cell
   is computed through the SAME toZ/fromZ/fmt/fmtPct/wechslerDesc/aanDesc
   the one-at-a-time converter uses — no second table of numbers exists to
   drift. Rows are anchored on ONE metric at that metric's real integer
   values (z in 0.25 steps): a merged table anchored on standard scores
   would print scaled "7.3", a score nobody can report, so the fractional
   equivalents stay in the non-anchored columns where they honestly mark
   that two metrics do not map cleanly. Descending, as the manuals print.
   The entered score's row is highlighted and scrolled into view, so the
   input above doubles as "find me in the table". */
const CONV_FULL_ANCHORS = {
  standard: { min: 40,  max: 160, step: 1,    dp: 0 },
  t:        { min: 10,  max: 90,  step: 1,    dp: 0 },
  scaled:   { min: 1,   max: 19,  step: 1,    dp: 0 },
  z:        { min: -4,  max: 4,   step: 0.25, dp: 2 }
};
let convFullAnchor = 'standard';
function renderConvFullTable(){
  const wrap = document.getElementById('conv-full-scroll');
  if (!wrap) return;
  const a = CONV_FULL_ANCHORS[convFullAnchor];
  const type = document.getElementById('conv-type')?.value;
  const val = document.getElementById('conv-value')?.value;
  const zIn = (val === '' || val == null) ? null : toZ(val, type);

  /* The highlighted row is the anchor value nearest the entered score,
     and only if the score actually lands inside the table's range. */
  let currentIdx = -1;
  const steps = Math.round((a.max - a.min) / a.step);
  if (zIn != null){
    const equiv = fromZ(zIn, convFullAnchor);
    const idx = Math.round((equiv - a.min) / a.step);
    if (idx >= 0 && idx <= steps) currentIdx = idx;
  }

  const cols = ['standard', 't', 'scaled', 'z'];
  /* One whole score of a metric spans this much of the z axis. A column
     COARSER than the anchor's row spacing prints each whole score only on
     the row where it actually falls, blank between — the manuals' own
     convention, and the honest one: repeating "21" down five standard-score
     rows would claim each of them converts to scaled 21. A column FINER
     than the anchor never repeats, so it prints its nearest whole score on
     every row. z always prints exact to 2 dp: it is the continuous ruler
     that shows where every score sits. */
  const CONV_Z_UNIT = { standard: 1/15, t: 1/10, scaled: 1/3 };
  const anchorZStep = convFullAnchor === 'z' ? a.step : CONV_Z_UNIT[convFullAnchor];
  /* Compute anchor values from the integer step index rather than by
     repeated addition, so the z column's 0.25 steps stay exact. Rows are
     built ascending and rendered descending. */
  const rowVals = [];
  for (let i = 0; i <= steps; i++){
    const v = a.min + i * a.step;
    rowVals.push({ v, z: toZ(v, convFullAnchor), cells: {} });
  }
  cols.forEach(m => {
    if (m === convFullAnchor){
      rowVals.forEach(r => { r.cells[m] = fmt(r.v, a.dp); });
      return;
    }
    if (m === 'z'){
      rowVals.forEach(r => { r.cells[m] = fmt(r.z, 2); });
      return;
    }
    if (CONV_Z_UNIT[m] > anchorZStep + 1e-9){
      /* Sparse: place each whole score of this metric on its nearest row.
         Spacing exceeds one row by construction, so no two collide. */
      rowVals.forEach(r => { r.cells[m] = ''; });
      const range = CONV_FULL_ANCHORS[m];
      for (let t = range.min; t <= range.max; t++){
        const idx = Math.round((fromZ(toZ(t, m), convFullAnchor) - a.min) / a.step);
        if (idx >= 0 && idx <= steps) rowVals[idx].cells[m] = fmt(t, 0);
      }
    } else {
      rowVals.forEach(r => { r.cells[m] = fmt(fromZ(r.z, m), 0); });
    }
  });
  const rowsHtml = [];
  for (let i = steps; i >= 0; i--){
    const r = rowVals[i];
    const ss = fromZ(r.z, 'standard');
    const cells = cols.map(m =>
      `<td class="conv-full-cell${m === convFullAnchor ? ' is-anchor' : ''}">${r.cells[m]}</td>`
    ).join('');
    rowsHtml.push(`<tr class="conv-full-row${i === currentIdx ? ' is-current' : ''}">
      ${cells}
      <td class="conv-full-cell">${fmtPct(fromZ(r.z, 'percentile'))}</td>
      <td class="conv-full-desc">${wechslerDesc(ss)}</td>
      <td class="conv-full-desc">${aanDesc(ss)}</td>
    </tr>`);
  }
  wrap.innerHTML = `<table class="conv-full-table">
    <thead><tr>
      <th class="conv-full-th${convFullAnchor === 'standard' ? ' is-anchor' : ''}">Standard</th>
      <th class="conv-full-th${convFullAnchor === 't' ? ' is-anchor' : ''}">T</th>
      <th class="conv-full-th${convFullAnchor === 'scaled' ? ' is-anchor' : ''}">Scaled</th>
      <th class="conv-full-th${convFullAnchor === 'z' ? ' is-anchor' : ''}">z</th>
      <th class="conv-full-th">Percentile</th>
      <th class="conv-full-th">Wechsler</th>
      <th class="conv-full-th">AACN</th>
    </tr></thead>
    <tbody>${rowsHtml.join('')}</tbody>
  </table>`;
  const cur = wrap.querySelector('.conv-full-row.is-current');
  if (cur) wrap.scrollTop = Math.max(0, cur.offsetTop - wrap.clientHeight / 2);
}
function convFullSetAnchor(name){
  if (!CONV_FULL_ANCHORS[name]) return;
  convFullAnchor = name;
  document.querySelectorAll('#conv-full-anchor [data-anchor]').forEach(b => {
    const on = b.dataset.anchor === name;
    b.classList.toggle('is-active', on);
    b.setAttribute('aria-selected', String(on));
  });
  renderConvFullTable();
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
  /* The full table follows the input type so the highlighted row lands in
     the anchored column — except percentile, which is an output, not a
     metric the table can anchor on (see SCORE_METRICS). */
  if (CONV_FULL_ANCHORS[this.value]) convFullSetAnchor(this.value);
  renderConverter();
});
document.getElementById('conv-value').addEventListener('input', renderConverter);
document.getElementById('conv-full-anchor')?.addEventListener('click', e => {
  const btn = e.target.closest('[data-anchor]');
  if (btn) convFullSetAnchor(btn.dataset.anchor);
});
/* Opening the disclosure re-renders so the highlighted row is scrolled to
   centre — scrollTop set while the details was closed does not stick. */
document.getElementById('conv-full-details')?.addEventListener('toggle', function(){
  if (this.open) renderConvFullTable();
});

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
  /* -1. What the CLINICIAN types in, which is not always what m1/sd1 are.

     CVLT-C is the case that forced these apart. Its normative statistics here
     are raw (metric:'raw'), which is what Change Analysis needs — but nobody
     transcribes a raw CVLT-C score onto a results table; they transcribe the
     standardised score off the record form. Manual Table A.1 gives List A
     Total Trials 1-5 as a T-score; Table A.2 gives every other index on a
     scale running -5.0 to +5.0, i.e. a z-score.

     BEWARE the manual's wording: Table A.2 is headed "Standard Score
     Equivalents", but that is NOT this app's 'standard' (M 100, SD 15). A
     CVLT-C "standard score" of -1.0 is the 16th percentile; read as M 100 /
     SD 15 it would be z = -6.7. reportedAs holds the real metric, not the
     manual's label for it. */
  if (stats && SCORE_METRICS.has(stats.reportedAs)) return stats.reportedAs;
  /* 0. Declared, and therefore not a guess — it wins over everything below.

     metric:'raw' is the case that forced this to exist: a raw score is not on
     a standardised metric at all, and the mean-based heuristic has no raw
     category, so it picked whichever standard metric the raw mean sat nearest.
     That is how RBANS List Recognition (raw, M 19.6, SD 0.8) came to be read
     as a scaled score. See the metric:'raw' note in data.js.

     The other four values are accepted too, so a custom test can also correct
     a standardised metric the heuristic gets wrong — the Norms Database score
     type column writes this field. normDB itself only ever uses 'raw'. */
  if (stats && SCORE_METRICS.has(stats.metric)) return stats.metric;
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
/* True for families with no retest data at all — every entry lacks m2/sd2/r,
   because none is published. WAIS-IV longest span is the first of these: the
   manual gives normative base rates but no reliability coefficient.

   They must be kept out of the Change Analysis and SD Index dropdowns. Loading
   one there produces a table of rows that can never compute anything, which
   reads as the app being broken rather than as data being absent. */
function isSingleAdministrationFamily(name){
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  const fam = db && db[name];
  if (!fam) return false;
  const entries = Object.values(fam).filter(e => e && typeof e === 'object');
  return entries.length > 0 && entries.every(e => e.singleAdministration);
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
  return {scaled:'Scaled Score', standard:'Standard Score', t:'T-Score', z:'Z-Score', raw:'Raw Score'}[type] || 'Score';
}
function scoreTypeAbbr(type){
  return {scaled:'Scaled', standard:'Standard', t:'T-Score', z:'Z-Score', raw:'Raw'}[type] || '';
}
/* The header must describe the rows that are actually there, not the
   "default for ungrouped rows" dropdown.

   It used to read straight from #bat-type, so autofilling CVLT-C — every row a
   T-score or a z-score — printed a column headed "Scaled Score" over them. The
   APA export already derived its header from the rows (renderBatteryApa's
   headerLabel), so screen and export disagreed too. Same rule in both places
   now: one label when the rows agree, "Score" when they do not, and the
   dropdown only while the table is empty. */
function updateBatteryScoreHeader(){
  const fallback = document.getElementById('bat-type')?.value;
  const scored = batteryRows.filter(r => !r.isExample);
  const types = new Set(scored.map(r => rowScoreType(r)));
  const label = scored.length === 0 ? scoreTypeLabel(fallback)
              : types.size === 1    ? scoreTypeLabel([...types][0])
              : 'Score';
  document.querySelectorAll('#bat-table .bat-score-head').forEach(th => {
    th.innerHTML = `<span class="bat-score-head-label">${escapeHtml(label)}</span>`;
    th.title = types.size > 1
      ? 'This table mixes score types. Each row shows its own metric next to the row number.'
      : `Manual rows use the selected default score type: ${scoreTypeLabel(fallback)}. Auto-filled test families show their inferred score type in the group row.`;
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
/* ---------- PUBLISHED BASE-RATE LOOKUP ----------
   Some measures are raw counts with no metric to convert — WAIS-IV longest
   span is the case this was built for. They cannot be scored by formula, but
   the manual publishes the cumulative percentage of the normative sample at
   each raw score, so the percentile can be read off instead of derived.

   Preferred over the normal approximation not because it is more accurate —
   computing from the printed M/SD lands within about 3 percentile points at
   every span — but because it IS the published figure. In a report that gets
   scrutinised, a number citable to Table C.4 needs no defending. */
function batteryBaseRateEntry(row){
  if (!row || !row.group || !row.name) return null;
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  const e = db && db[row.group] && db[row.group][row.name];
  return (e && e.baseRates && typeof e.baseRates === 'object') ? e : null;
}
/* baseRates[x] is the percentage scoring x OR HIGHER, so P(X >= v) is the
   value at the smallest tabulated span that is at least v. Spans the table
   skips (there is no row for 1) therefore inherit the next one up, which is
   correct: nobody scores between two adjacent whole spans. */
function baseRateAtOrAbove(entry, v){
  const spans = Object.keys(entry.baseRates).map(Number).filter(Number.isFinite);
  const atOrAbove = spans.filter(s => s >= v);
  if (!atOrAbove.length) return 0;               // above the top of the table
  return entry.baseRates[Math.min(...atOrAbove)];
}
/* Percentile RANK, i.e. the percentage scoring below, using the midpoint
   convention for a discrete score: P(X < v) + 0.5 * P(X = v). The column means
   "percentage scoring below" everywhere else in the app, so the published
   "or higher" figure is converted rather than shown raw — a column that
   silently reverses direction for some rows would be a trap. */
function baseRatePercentile(entry, value){
  const v = parseFloat(value);
  if (!Number.isFinite(v) || !Number.isInteger(v)) return null;   // spans are whole numbers
  const pGe = baseRateAtOrAbove(entry, v);
  const pGt = baseRateAtOrAbove(entry, v + 1);
  return (100 - pGe) + 0.5 * (pGe - pGt);
}
/* A BASE-RATE ROW IS ENTERED IN THE RAW COLUMN, and only there. A longest
   span of 6 is a raw count — it has no scaled/standard equivalent to put in
   the Score column, so that box is disabled for these rows and the value is
   typed under Raw Score. (It used to be entered in the Score column, which
   implied a metric score existed; renderBattery migrates any such legacy
   value into the raw field.) */
function batteryBaseRateValue(r){ return r ? r.raw : ''; }

/* THE AGE IS REQUIRED BEFORE A BASE-RATE ROW SCORES. The published lookup is
   per normative age band (Tables C.4–C.5 differ hugely by band — Longest Digit
   Span Backward at a span of 4 is the 22nd percentile at 20-24 and the 59th at
   85-90), and the band lives in the group the row was added from. Requiring
   the patient age (a) stops a table being scored with no age on record, and
   (b) catches the entered age contradicting the band the row was picked from —
   scoring an 85-year-old on the 20-24 table silently would be a misstatement
   on a printed report. Returns 'ok' | 'no-age' | 'out-of-band'. */
function batteryBaseRateAgeState(r){
  if (!batteryBaseRateEntry(r)) return 'ok';
  const band = ageBandRange(r.group);
  if (!band) return 'ok';               // un-banded base-rate group: nothing to hold the age against
  const age = batteryPatientAge();
  if (age == null) return 'no-age';
  return (age >= band.lo && age <= band.hi) ? 'ok' : 'out-of-band';
}

/* One answer for "what percentile is this row", used by the table, the in-place
   update and the APA export so the three cannot drift. */
function batteryRowPercentile(r){
  const entry = batteryBaseRateEntry(r);
  if (entry){
    if (batteryBaseRateAgeState(r) !== 'ok') return null;
    return baseRatePercentile(entry, batteryBaseRateValue(r));
  }
  const z = toZ(r.score, rowScoreType(r));
  return z == null ? null : normCDF(z) * 100;
}

/* ---------- PERCENTILE vs BASE RATE: two quantities, one column ----------
   A base rate and a percentile run in OPPOSITE directions. The manual's 79 for
   a longest span of 6 means 79% scored 6 or more; the percentile rank for the
   same score is 30, because 30% scored below it.

   The first attempt at this made the whole column switch between the two
   depending on what else was in the table. That was wrong twice over. It meant
   a longest-span row silently changed quantity when an unrelated subtest was
   entered — and because the mode ignored seeded example rows while still
   rendering them, a table could show "Base rate" over a row holding a
   percentile. Both were visible on screen.

   So: a base-rate measure ALWAYS reports its published base rate, and the
   section carries its own column heading saying so. The main heading never
   changes. One quantity per section, each labelled where it is read. */
function batteryGroupIsBaseRate(gKey){
  const rows = batteryRows.filter(r => batteryGroupKeyOf(r) === gKey);
  return rows.length > 0 && rows.every(r => batteryBaseRateEntry(r));
}
const BAT_PCT_LABEL = 'Percentile';
const BAT_BASERATE_LABEL = 'Base rate';
/* Base rates are published to 0.5 and legitimately reach 100 — fmtPct must not
   be used, as it clamps the tails into (0.01, 99.99) to keep percentiles inside
   the open interval. 100% of a sample really can score at or above the floor. */
function fmtBaseRate(v){
  if (v == null || isNaN(v)) return '-';
  return Number.isInteger(v) ? String(v) : v.toFixed(1);
}
/* The cell value, and which quantity it is. A base-rate measure always reports
   its published base rate — never converted, whatever else is in the table. */
function batteryRowPctCell(r){
  const entry = batteryBaseRateEntry(r);
  if (entry){
    const v = parseFloat(batteryBaseRateValue(r));
    if (!Number.isFinite(v) || !Number.isInteger(v)) return null;
    /* A value is entered but the age gate refuses. `hint` is a SCREEN-ONLY
       affordance: the APA export reads .text, which stays empty, so a blocked
       row exports as blank rather than as advice to the reader. */
    const ageState = batteryBaseRateAgeState(r);
    if (ageState !== 'ok') return { value: null, text: '', kind: 'baseRate', hint: ageState };
    const br = baseRateAtOrAbove(entry, v);
    return { value: br, text: fmtBaseRate(br), kind: 'baseRate' };
  }
  const p = batteryRowPercentile(r);
  return p == null ? null : { value: p, text: fmtPct(p), kind: 'percentile' };
}
/* The screen text for a gated cell. Short, factual, and pointing at the fix:
   the master age field is in the top bar. */
function batteryAgeHintHtml(state){
  return state === 'out-of-band'
    ? '<span class="bat-age-hint">age outside this band</span>'
    : '<span class="bat-age-hint">add patient age ↑</span>';
}
function batteryClassificationDetails(r, cls){
  /* A base-rate measure has a real percentile but no metric, so its z comes
     BACK from that percentile rather than from the score. The classification
     bands are defined on standard scores, and reading an empirical percentile
     onto them ("performance at the 5th percentile") is ordinary practice. */
  const brEntry = batteryBaseRateEntry(r);
  let z;
  if (brEntry){
    /* Blank while the age gate refuses: the pct cell carries the hint, and a
       classification derived from an unverified band would be the very number
       the gate exists to withhold. */
    if (batteryBaseRateAgeState(r) !== 'ok') return { text:'', html:'', className:'' };
    const pct = baseRatePercentile(brEntry, batteryBaseRateValue(r));
    if (pct == null) return { text:'', html:'', className:'' };
    // clamp inside the open interval: normInv is undefined at 0 and 100, and
    // the table legitimately reaches both ends.
    z = normInv(Math.min(99.99, Math.max(0.01, pct)) / 100);
  } else {
    z = toZ(r.score, rowScoreType(r));
  }
  if (z == null) return { text:'', html:'', className:'' };
  const ss = fromZ(z, 'standard');
  /* On a higher-is-worse measure the score runs the other way: CVLT-C
     Perseverations is scored so that MORE perseverations give a HIGHER
     standardised score (Manual Table A.2 maps z -1.0 to 0-3 perseverations and
     z +5.0 to 45 or more). A child at z +2.0 has made more perseverations than
     98% of their age group, so the merit bands would label a poor result
     "Very Superior".

     The classification therefore describes PERFORMANCE, computed from the
     reflected score — z +2.0 classifies as z -2.0 does, "Extremely Low", which
     is what a clinician would write. The PERCENTILE is deliberately NOT
     reflected: 98th is what z +2.0 gives and what anyone checking the working
     against the manual will calculate, and a number that cannot be reproduced
     is worth less in a report than one that needs a footnote. The APA note
     carries the convention so the pairing does not read as a contradiction.

     Blanking the label was tried first and rejected: the column exists to
     summarise performance, and a dash reads as missing data rather than as a
     deliberate choice. */
  const perfSs = r.higherIsWorse ? fromZ(-z, 'standard') : ss;
  const desc = cls === 'wechsler' ? wechslerDesc(perfSs) : aanDesc(perfSs);
  const prem = getBatteryPremorbidComparison();
  const extreme = perfSs < 70;
  /* No premorbid asterisks on an error measure. The stars mean "this ability
     falls short of the premorbid estimate", and a perseveration count is not
     an ability being predicted — the comparison has no meaning to assert. */
  const stars = r.higherIsWorse ? '' : batteryPremorbidStars(ss, prem);
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
/* Normative SD of each standardised metric — the multiplier in SD x sqrt(1-r).
   Mirrors fromZ's coefficients. 'raw' is deliberately absent: a raw score is
   not on a metric at all, so those rows fall back to the entry's own sd1. */
const BATTERY_METRIC_SD = { z: 1, t: 10, scaled: 3, standard: 15 };

/* The reliability and the SD must describe the SAME population, or SD^2(1-r)
   is not an error variance. Two pairings are coherent:
     A  rCorrected x normative SD  — both describe the normative population
     B  r          x sd1           — both describe the retest sample
   This used to take rCorrected AND sd1, which is neither. rCorrected is the
   coefficient rescaled to the normative sample's variability, so pairing it
   with the retest sample's SD mixes two populations — the same mistake
   rciState rejects by defaulting useCorrectedR to false, and for the same
   reason. On the RBANS Form C/D groups, whose retest samples are impaired
   (m1 88-92, sd1 17-19), it inflated the SEM by up to 3.2 points.

   Pairing A is preferred because it reproduces the SEM the PUBLISHER prints.
   The CVLT-3 Manual gives an SEM column in Tables 3.4 and 3.5, and
   normative SD x sqrt(1 - corrected r) reproduces all 38 of them exactly —
   pinned in check.js section 23.

   Where no corrected r is published (310 entries) the raw r is used with the
   normative SD anyway. That mixes populations, worth roughly 20% of the error
   variance, and is the lesser error: sd1 is the SD of the STORED statistics,
   which are not always in the units of the score on screen.

   THE SD FOLLOWS THE METRIC ON SCREEN, not the metric of the stored figures.
   sd1 is the SD of the STORED statistics, and those are not always in the
   units the clinician typed. CVLT-C is the case that forced this: its norms
   here are raw word counts (Change Analysis needs them that way) but nobody
   transcribes a raw CVLT-C score — they copy the z or T off the record form,
   which is what reportedAs records. Building the SEM from a raw word-count SD
   and printing it against a z-score gave Discriminability an interval of
   +/-7.5 z at age 8, and Perseverations +/-4.3. sd1 is therefore reached for
   only when the row is SHOWN raw, where it is the one SD in the right units.

   RELIABILITY, in preference order:
     rInternal   — internal consistency, where the publisher reports one AND
                   builds its own published intervals from it. See data.js.
     rStability  — a stability coefficient the publisher tabulates in that SAME
                   reliability table, for measures it excludes from internal
                   consistency. RBANS Update Table 3.6 is the first: footnote a
                   marks five of its fourteen rows "based on test-retest".
     rCorrected  — retest r rescaled to the normative sample.
     r           — the raw retest coefficient.
   Internal consistency comes first because a confidence interval asks how
   precisely THIS administration measured the person, which is the question an
   internal-consistency coefficient answers. That preference is confined to
   this function: reliable change asks a different question and must keep the
   retest coefficient — see the rInternal note in data.js.

   rStability is NOT an exception to the test-retest default; it IS that
   default, sourced from the manual's own reliability table rather than from a
   retest study group. It exists as a separate field only so that a stability
   coefficient is never stored under a name the APA note and Methods &
   References then describe as internal consistency. */
/* Does this row's measure actually consult the patient age? Only measures
   whose publisher tabulates reliability by age band do, so the APA note can
   say the age was used without asserting it on tables where nothing read it.
   Either banded field counts — the question is whether an age was read, not
   which basis the coefficient has. */
function batteryRowUsesAgeBand(row){
  if (!row || !row.group) return false;
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  const entry = db && db[row.group] && db[row.group][row.name];
  return !!(entry && typeof entry === 'object' && (entry.rInternalByAge || entry.rStabilityByAge));
}

/* THE PATIENT'S AGE — one value, one master, one mirror.

   There is one patient and therefore one age. It is typed into #patient-age in
   the top bar, which is on screen on every page, and mirrored into #pre-age on
   Premorbid, which keeps a field because Age is a genuine member of that
   demographic block (OPIE-4 reads it alongside Sex). Two boxes that can hold
   different ages for the same patient is a defect waiting to happen, so they
   are kept in step.

   WHY THE MASTER MOVED TO THE TOP BAR. It used to be #bat-age, declared in the
   Score Tables markup and moved by design-system.js into that page's inline
   control bar, revealed only once a CI level was chosen. Three problems, all
   from the same cause — patient-level data living in a page-level control bar:
   it shipped rendering at 0x0 inside a wholesale-hidden panel; once moved, it
   took width from the test-family search and its label had to be truncated to
   "Age"; and being gated behind the CI toggle meant a clinician who never
   turned CI on never saw an age field. In the top bar it is simply always
   there, on every page, under its own name.

   WHY THAT IS SAFE DESPITE THE DIFFERENT RANGES. #pre-age used to declare
   min="16" because ToPF and OPIE-4 are adult instruments, while Score Tables
   needs from age 5 (CVLT-C) and 8 (D-KEFS Advanced). Both premorbid consumers
   bound the value themselves before computing — OPIE-4 on
   OPIE_AGE_MIN/OPIE_AGE_MAX, the ToPF-predicted WMS-IV indices on
   TOPF_AGE_MIN/TOPF_AGE_MAX. So a paediatric age makes those decline, which is
   correct because they are not normed for children, and changes nothing else.
   The min attribute was a claim the code did not rely on, and is relaxed to
   match so the browser stops flagging a legitimate age as invalid.

   ORDER MATTERS HERE. patientAge() returns the first input holding a finite
   value, so the master is listed first: if the two ever drift, the top-bar box
   — the one the clinician can see from every page — is the one that wins.

   The listeners write the sibling's .value and do NOT dispatch an event on it.
   Dispatching would have each input echo the other indefinitely. */
const PATIENT_AGE_INPUTS = ['patient-age', 'pre-age'];

function patientAge(){
  for (let i = 0; i < PATIENT_AGE_INPUTS.length; i++){
    const el = document.getElementById(PATIENT_AGE_INPUTS[i]);
    const v = el ? parseFloat(el.value) : NaN;
    if (Number.isFinite(v)) return v;
  }
  return null;
}

/* Copy the edited box into the other one. Value only — see above. */
function syncPatientAge(sourceId){
  const src = document.getElementById(sourceId);
  if (!src) return;
  PATIENT_AGE_INPUTS.forEach(id => {
    if (id === sourceId) return;
    const el = document.getElementById(id);
    if (el && el.value !== src.value) el.value = src.value;
  });
}

/* Kept as the Score Tables reader so the CI code states what it wants. */
function batteryPatientAge(){ return patientAge(); }

/* IS THE AGE ACTUALLY BEING READ ON THE PAGE IN FRONT OF YOU?

   The field is permanent and on every page, but only measures whose publisher
   tabulates reliability by age band consult it — D-KEFS Advanced, D-KEFS
   original and WAIS-IV. On a table of CVLT-3 and RBANS rows the age changes
   nothing. A permanently visible "Patient age 72" would imply otherwise, so the
   pip beside it lights only where something reads the value.

   Score Tables uses EXACTLY the expression that decides the APA note's ciAge,
   so the screen and the exported table cannot disagree about whether the age
   was read. Premorbid asks whether the age is in range for a model on that page
   — OPIE-4, or the ToPF-predicted WMS-IV indices — since those are the only two
   consumers there.

   This is a UI affordance, not a printed claim: the load-bearing statement is
   still the APA note, which is unchanged. */
/* HOW MANY ROWS ON SCORE TABLES WOULD READ AN AGE.

   Deliberately independent of whether an age is actually set, because two
   opposite features ask this same question and must never disagree:

     - the topbar pip, which lights when an age IS set and something reads it;
     - the missing-age prompt above the table, which appears when one is NOT.

   If these were computed separately, a table could show the prompt and the pip
   at once, or neither. One predicate makes that unrepresentable.

   The CI level is part of the question for the reliability half only: with
   the interval switched off, an age-band COEFFICIENT changes nothing printed.
   Base-rate rows are the other half, and they are CI-independent — their
   scoring itself is gated on the age (batteryBaseRateAgeState), whatever the
   interval setting. */
function batteryCiAgeBandRowCount(){
  const ci = document.getElementById('bat-ci-level');
  if (!ci || ci.value === 'off') return 0;
  if (!Array.isArray(batteryRows)) return 0;
  return batteryRows.filter(r => r && r.name && !r.isExample && batteryRowUsesAgeBand(r)).length;
}
/* Rows whose SCORING requires the age: the published base-rate lookups. */
function batteryBaseRateRowCount(){
  if (!Array.isArray(batteryRows)) return 0;
  return batteryRows.filter(r => r && r.name && !r.isExample && batteryBaseRateEntry(r)).length;
}
function batteryAgeBandRowCount(){
  return batteryCiAgeBandRowCount() + batteryBaseRateRowCount();
}

function patientAgeIsInUse(){
  if (patientAge() === null) return false;
  const active = document.querySelector('.section.active');
  const page = active ? active.id : '';

  if (page === 'battery') return batteryAgeBandRowCount() > 0;

  if (page === 'premorbid'){
    const a = patientAge();
    const inOpie = (typeof OPIE_AGE_MIN === 'number' && typeof OPIE_AGE_MAX === 'number')
                   && a >= OPIE_AGE_MIN && a <= OPIE_AGE_MAX;
    const inTopf = (typeof TOPF_AGE_MIN === 'number' && typeof TOPF_AGE_MAX === 'number')
                   && a >= TOPF_AGE_MIN && a <= TOPF_AGE_MAX;
    return inOpie || inTopf;
  }

  return false;
}

/* Toggled by CLASS, never by the [hidden] attribute. The last element that
   tried [hidden] here was .ds-inline-bar-age, whose own display:flex outranked
   the browser default and left it on screen — see the CSS cascade note in
   CLAUDE.md. A class the stylesheet owns outright cannot lose that fight. */
function refreshPatientAgeIndicator(){
  const field = document.getElementById('patient-age-field');
  if (!field) return;
  field.classList.toggle('is-live', patientAgeIsInUse());
}

/* THE REQUIRED-AGE BAR.

   The predecessor here was a dismissable POPOVER anchored to the Score CI
   toggle: edge-triggered, with a skip button — an offer, because a blank age
   only cost a sharper interval. Two things ended it (owner decision,
   2026-08): base-rate rows made the age load-bearing for SCORING, not just
   interval width, and the anchor stopped making sense once the ask could
   arise with the interval untouched — it drew its arrow at the 90% CI button
   while talking about scoring. It also simply was not effective: a floating
   popover reads as ignorable, and the ask is not ignorable any more.

   What replaced it is this inline bar directly above the table, shown
   whenever any row actually reads the age (batteryAgeBandRowCount() > 0) and
   none is stored. It is REQUIRED in presentation: no skip, no outside-click
   dismissal, it stays until an age is entered. Note this reverses the old
   "the gate must not ask about the table's contents" rule — that rule
   belonged to the offer, whose n === 0 branch existed so a CI-only ask never
   looked arbitrary; a REQUIREMENT shown on a table where nothing reads the
   age would be a false claim, so the bar is gated on actual readers.

   STATE-DRIVEN, NOT EDGE-TRIGGERED. An inline element cannot "re-open"
   annoyingly the way a popover could, so the whole edge/arming machinery
   (batAgePopLastWanted, batAgePopJustOpened, the opening-click guard, the
   fixed-position zoom mathematics) went with the popover. The bar is
   rewritten only when its message changes, so typing into its own age input
   is never clobbered by a table-keystroke re-render.

   It shares its predicate with the topbar residual (.is-wanted): one
   predicate, so bar and residual cannot disagree. */
function batteryAgeWanted(){
  return batteryAgeBandRowCount() > 0 && patientAge() === null;
}

/* body{zoom:0.9} splits measurement into visual and layout px — anything
   positioning from a measured rect divides by this. Survives the popover it
   was written for: the consent card's pill stacking still reads it, and a
   missing definition would fall through that call site's typeof guard to 1,
   silently reintroducing the exact offset bug CLAUDE.md documents. */
function pageZoomFactor(){
  const z = parseFloat(getComputedStyle(document.body).zoom);
  return Number.isFinite(z) && z > 0 ? z : 1;
}

/* Every clause agrees with its count — written out per branch, not assembled
   from fragments; see the note on the old popover text in check.js §26. */
function batteryAgeBarHtml(){
  const br = batteryBaseRateRowCount();
  const bands = batteryCiAgeBandRowCount();
  const parts = [];
  if (br === 1) parts.push('<strong>1 measure</strong> in this table is scored from an age-banded base-rate table and stays unscored until the age is entered.');
  else if (br > 1) parts.push('<strong>' + br + ' measures</strong> in this table are scored from age-banded base-rate tables and stay unscored until the age is entered.');
  if (bands === 1) parts.push('<strong>1 measure</strong> publishes its reliability by age band; the age selects the correct coefficient for its interval.');
  else if (bands > 1) parts.push('<strong>' + bands + ' measures</strong> publish their reliability by age band; the age selects the correct coefficients for their intervals.');
  return '<div class="bat-age-bar-main">'
    + '<div class="bat-age-bar-title">Patient age required</div>'
    + '<div class="bat-age-bar-body">' + parts.join(' ') + '</div>'
    + '</div>'
    + '<div class="bat-age-bar-row">'
    + '<input type="number" id="bat-age-bar-input" class="bat-age-bar-input" min="5" max="110" step="1" placeholder="e.g. 72" aria-label="Patient age">'
    + '<button type="button" class="bat-age-bar-add" id="bat-age-bar-add">Add age</button>'
    + '</div>';
}

/* The age was typed into the bar above the table — so the field that now
   holds it, up in the top bar, is somewhere the clinician was not looking. The pulse says "it landed here", and answers the question the
   popover leaves behind: where did that value actually go?

   Restarted rather than re-added, so a second commit re-runs the animation
   instead of doing nothing because the class is already present. */
let patientAgePulseTimer;
function pulsePatientAgeField(){
  const field = document.getElementById('patient-age-field');
  if (!field) return;
  field.classList.remove('is-pulsing');
  void field.offsetWidth;                       // force reflow to restart it
  field.classList.add('is-pulsing');
  clearTimeout(patientAgePulseTimer);
  patientAgePulseTimer = setTimeout(() => field.classList.remove('is-pulsing'), 2200);
}

/* Writes through to the master and lets its own listener do the rest — the
   sync, the re-render and the pip all hang off that one 'input' event. */
/* Writes through to the master and lets its own listener do the rest — the
   sync, the re-render and the pip all hang off that one 'input' event. */
function commitBatteryAgeBar(){
  const input = document.getElementById('bat-age-bar-input');
  const master = document.getElementById('patient-age');
  if (!input || !master) return;
  const v = parseFloat(input.value);
  if (!Number.isFinite(v)) { input.focus(); return; }
  master.value = String(v);
  master.dispatchEvent(new Event('input', { bubbles: true }));
  /* After the dispatch, so the pulse lands on a field already showing the new
     value and its .is-live pip rather than on a stale one. */
  pulsePatientAgeField();
}

function refreshBatteryAgePrompt(){
  const wanted = batteryAgeWanted();
  /* ONE predicate for the bar and the topbar residual, so the table cannot
     demand an age while the field beside the age box implies nothing reads
     one, or vice versa. */
  const field = document.getElementById('patient-age-field');
  if (field) field.classList.toggle('is-wanted', wanted);
  const bar = document.getElementById('bat-age-required');
  if (!bar) return;
  bar.hidden = !wanted;
  if (wanted){
    const html = batteryAgeBarHtml();
    /* Rewrite only on change: renderBattery() runs on every keystroke in the
       score cells, and clobbering the bar would clear its half-typed age. */
    if (bar.dataset.sig !== html){ bar.dataset.sig = html; bar.innerHTML = html; }
  }
}

document.addEventListener('click', e => {
  if (e.target.closest('#bat-age-bar-add')) commitBatteryAgeBar();
});
document.addEventListener('keydown', e => {
  if (e.key === 'Enter' && e.target && e.target.id === 'bat-age-bar-input'){
    e.preventDefault();
    commitBatteryAgeBar();
  }
});

/* Pick the published coefficient for the patient's age band.
   rInternalByAge is keyed by the LOWER BOUND of each normative band, so the
   lookup takes the greatest key <= age. That keeps the table compact and
   matches how the manual prints it (8, 9, ... 15, 16-18, 19-29, ... 80-90).

   Returns null rather than guessing when the age is absent or outside the
   normed range, and the caller then falls back to the published all-ages
   average. Both figures are the publisher's own, so either path is citable —
   and a blank age must never empty the column, which is why this degrades
   instead of refusing. */
function bandedReliabilityForAge(byAge, ageMax, age){
  if (!byAge || !Number.isFinite(age)) return null;
  if (Number.isFinite(ageMax) && age > ageMax) return null;
  let band = null;
  Object.keys(byAge).forEach(k => {
    const lo = Number(k);
    if (age >= lo && (band === null || lo > band)) band = lo;
  });
  return band === null ? null : byAge[band];
}
function rInternalForAge(entry, age){
  return bandedReliabilityForAge(entry.rInternalByAge, entry.rInternalAgeMax, age);
}
/* Same lookup, for the stability coefficients a manual tabulates by age
   alongside its internal-consistency ones. Kept as its own function rather
   than a flag so that a grep for either field finds every site that reads it. */
function rStabilityForAge(entry, age){
  return bandedReliabilityForAge(entry.rStabilityByAge, entry.rStabilityAgeMax, age);
}

/* THE RANGE-RESTRICTION CORRECTION — OFFERED, NOT IMPOSED.

   Every interval here is SD x sqrt(1 - rxx), which is only a valid error
   variance when both terms describe the same group. 298 entries pair the retest
   study's own r with the NORMATIVE SD of the displayed metric — two
   populations. Allen & Yen / Magnusson repairs that arithmetically:

       rxx(unrestricted) = 1 - (sd1^2 / normSD^2)(1 - r)

   with sd1 the retest sample's SD at first testing. It is a good formula, not a
   guess: over the 267 entries carrying both a raw r and a published rCorrected
   it reproduces the publisher's own figure to a median error of .003.

   IT IS NOT THE DEFAULT, AND MUST NOT BECOME ONE. Every measure it reaches on
   Score Tables is D-KEFS or D-KEFS Advanced, and both manuals make the
   uncorrected pairing themselves, deliberately — D-KEFS Technical Manual p. 19
   fixes SD 3 for all scaled scores and derives its test-retest SEMs from the
   total sample, which Table 2.8 shows on the page. Applying this by default
   would print a coefficient the cited manual does not contain. It is offered
   because the statistical objection is real, and the APA note names whichever
   basis produced the table. See resolveCiReliability and check.js S27.

   Three refusals, each a distinct way of being wrong:

   - NO NORMATIVE SD. Undefined only when the row is displayed raw, and then
     sd1 already IS the SD in force — same sample, nothing to correct.
   - metric:'raw', EVEN WHEN DISPLAYED ON A STANDARDISED METRIC. CVLT-C is the
     case: sd1 is in words recalled while reportedAs puts the row on T or z, so
     the ratio of the two is meaningless.
   - A RESULT OUTSIDE (0, 1) IS NOT A RELIABILITY. Two entries go negative —
     D-KEFS Design Fluency Switching (Ages 20-49) and D-KEFS Advanced Social
     Sorting Total Number of Conceptual Level Responses (Ages 8-18) — both
     retest samples being more variable than the norm group on a coefficient
     already near .22. They keep the published r even with the toggle on.

   NEVER REACHES RELIABLE CHANGE. There r pairs with sd1/sd2, the same sample's
   SDs, so the pairing is already self-consistent; and in McSweeney and Crawford
   r is a fitted regression slope, where a population-corrected value gives a
   line that was never fitted to anything. */
function derivedCorrectedR(entry, normSD){
  if (!entry || !Number.isFinite(normSD) || normSD <= 0) return null;
  if (entry.metric === 'raw') return null;
  if (!Number.isFinite(entry.r) || !Number.isFinite(entry.sd1) || entry.sd1 <= 0) return null;
  const rxx = 1 - (entry.sd1 * entry.sd1) / (normSD * normSD) * (1 - entry.r);
  return (rxx > 0 && rxx < 1) ? rxx : null;
}

/* Score Tables' reliability-basis control, read here so that every consumer —
   the table, the APA note and the Data page — asks one question of one input.
   Absent element means published, which is also the markup's default: a
   missing control must never silently switch the app onto derived numbers. */
function batteryCiCorrectRetest(){
  return document.getElementById('bat-ci-basis')?.value === 'corrected';
}

/* THE ONE PLACE THAT DECIDES WHICH COEFFICIENT AN INTERVAL RESTS ON.

   Score Tables renders the interval; the Data page tells the clinician what it
   rests on. Those were two functions reading the same fields in the same order,
   pinned together by check.js S32 — which caught nothing for a while and then
   had to be strengthened twice, because a mirror is only ever as good as the
   next person's memory. They now share this one resolver, so the two cannot
   drift; S32 still drives both entry points, which is where the remaining
   drift risk lives (the arguments each one passes in).

   Returns the coefficient AND the name of its source, because .90 from a
   split-half table and .90 from a retest study are different claims and this
   project treats an on-screen statement of method as a contract.

   `age` is optional throughout. The Data page has no patient, so it renders
   the blank-age answer, which is also the app's own fallback.

   `correctRetest` is the Score Tables reliability-basis control, and it reaches
   ONLY the last branch — a published coefficient is never overwritten, whatever
   the toggle says. */
function resolveCiReliability(entry, normSD, age, correctRetest){
  if (!entry || typeof entry !== 'object') return null;
  const banded = !!(entry.rInternalByAge || entry.rStabilityByAge);
  const rBand  = rInternalForAge(entry, age);
  const rSBand = rStabilityForAge(entry, age);
  let r = null, basis = null;
  if (Number.isFinite(rBand))                 { r = rBand;            basis = 'internal consistency'; }
  else if (Number.isFinite(rSBand))           { r = rSBand;           basis = 'stability, published'; }
  else if (Number.isFinite(entry.rInternal))  { r = entry.rInternal;  basis = 'internal consistency'; }
  else if (Number.isFinite(entry.rStability)) { r = entry.rStability; basis = 'stability, published'; }
  else if (Number.isFinite(entry.rCorrected)) { r = entry.rCorrected; basis = 'retest, corrected'; }
  /* THE UNCORRECTED RETEST COEFFICIENT, AND WHY IT STAYS UNCORRECTED.

     This is the last resort, and it pairs a coefficient measured on the retest
     STUDY's sample with the NORMATIVE SD of the displayed metric. Those are
     different populations, which normally makes SD sqrt(1 - rxx) an invalid
     error variance — the very mismatch the Score Tables CI fixed in the other
     direction when it stopped pairing rCorrected with sd1.

     APPLYING THE RANGE-RESTRICTION CORRECTION HERE WAS BUILT, TESTED AND
     REJECTED. Allen & Yen / Magnusson gives

         rxx(unrestricted) = 1 - (sd1^2 / normSD^2)(1 - r)

     and it is a good formula: over the 267 entries carrying both a raw r and a
     published rCorrected it reproduces the publisher's own figure to a median
     error of .003. It moves 246 stored entries, 46 of them reachable from
     Score Tables — and every one of those 46 is D-KEFS or D-KEFS Advanced.

     CORRECTING BY DEFAULT WAS REJECTED, because EVERY ONE of those 46 belongs
     to a manual that has already made this pairing itself, deliberately.
     D-KEFS Technical Manual p. 19 says test-retest SEMs "were derived from the
     total sample of cases" and fixes "The standard deviation unit is 3 for all
     D-KEFS scaled scores" — i.e. 3 x sqrt(1 - r) on the UNCORRECTED
     total-sample r. Table 2.8 proves it: Design Fluency's three All Ages SEMs
     are 1.94 / 1.97 / 2.47, and 3 x sqrt(1 - r) gives exactly those, where the
     corrected coefficient gives 1.78 / 1.95 / 2.43. The remaining 21 are
     D-KEFS Advanced Trail Making and Verbal Fluency, whose Table 3.4 rows ARE
     the retest coefficients — that manual's stated choice for its speeded
     tests.

     So correcting here by default would print a number the cited manual does
     not contain, on measures where the publisher's own arithmetic is on the
     page. That is the same ground on which the unrounded Fisher's-z WAIS-IV
     average was declined: the app renders the coefficient it actually used,
     and a clinician cross-checking the manual must find it there.

     WHAT SHIPPED IS THE LABEL, AND THE CHOICE. "retest, uncorrected" says
     plainly that this interval rests on a retest-study coefficient used with a
     normative SD, so the compromise is visible rather than silent; and the
     Score Tables reliability control lets a clinician take the corrected
     reading deliberately, with the APA note saying which they took. Off by
     default, exactly as the RCI pages' corrected-r toggle is.

     Two quite different states share the published r and the label separates
     them: with a normative SD in force the pairing is the publisher's own
     compromise; with none, the row is shown raw and r sits beside its own
     sample's SD, which is simply correct and needs no qualifier — and is why
     the correction is not offered there either. */
  else if (Number.isFinite(entry.r))          { r = entry.r;
    basis = Number.isFinite(normSD) ? 'retest, uncorrected' : 'retest';
    if (correctRetest){
      const derived = derivedCorrectedR(entry, normSD);
      if (derived !== null)                   { r = derived;          basis = 'retest, corrected here'; }
    }
  }
  if (r === null){
    /* No coefficient at all. Two quite different reasons, and the row should
       say which: a base-rate measure is scored by published lookup and never
       had one, whereas a raw RBANS subtest is simply absent from its manual's
       reliability table. Neither prints an interval. */
    if (entry.baseRates) return { r:null, basis:'base rate — no interval', none:true };
    return { r:null, basis:'none published', none:true };
  }
  /* Unconditional, and it has to be: D-KEFS (original) publishes
     rInternalByAge with NO all-ages average, so at a blank age it lands on the
     retest branch above while still being a measure whose printed figure moves
     the moment an age is entered. Suffixing only the internal-consistency
     branches would silence exactly the rows that most need the warning. */
  if (banded) basis += ' · by age';
  return { r, basis, none:false };
}

/* `correctRetest` is optional and defaults to the control's current state, the
   same way getBatteryCiHtml defaults `age`. Passing it explicitly is what lets
   both readings be exercised without a DOM. */
function getBatteryRowReliability(row, type, age, correctRetest){
  if (!row.group) return null;
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  if (!db) return null;
  const family = db[row.group];
  if (!family) return null;
  const entry = family[row.name];
  if (!entry || typeof entry !== 'object') return null;
  /* BATTERY_METRIC_SD has no 'raw' key, so a row shown raw falls through to
     the entry's own SD — the only one in its units — and, for the same reason,
     takes no range-restriction correction. */
  const normSD = BATTERY_METRIC_SD[type];
  const effCorr = correctRetest !== undefined ? correctRetest : batteryCiCorrectRetest();
  const rel = resolveCiReliability(entry, normSD, age, effCorr);
  if (!rel || rel.none) return null;
  const sd = Number.isFinite(normSD) ? normSD
           : (Number.isFinite(entry.sd1) ? entry.sd1 : null);
  return sd !== null ? { r: rel.r, sd, basis: rel.basis } : null;
}
/* DOES ANY ROW IN THIS TABLE REST ON THIS BASIS?

   Asked of the SHIPPED resolver, at the age and control state actually in
   force, rather than re-tested against the entry's fields — a second reading of
   the preference chain is a second thing to keep in step, and the APA note is a
   contract about what produced the numbers above it.

   PREFIX MATCH, because the basis can carry a " · by age" suffix: D-KEFS
   (original) publishes rInternalByAge with no all-ages average, so at a blank
   age such a row is BOTH on a retest coefficient and a measure whose figure
   moves the moment an age is entered. */
function batteryBasisPresent(rows, prefix){
  return rows.some(r =>
    ((getBatteryRowReliability(r, rowScoreType(r), batteryPatientAge()) || {}).basis || '')
      .startsWith(prefix));
}

/* WHY A ROW'S INTERVAL IS BLANK — and the two reasons must stay apart.

   getBatteryRowReliability returns null for both, which is right for the
   arithmetic and useless for the note: a base-rate measure is scored by
   published lookup and never had a coefficient, whereas a raw RBANS subtest is
   simply absent from its manual's reliability table. resolveCiReliability
   already distinguishes them ('base rate — no interval' vs 'none published'),
   so ask it directly rather than re-deriving the difference here.

   Verified over normDB at age 45: 51 base-rate rows and 4 none-published,
   against 209 that do print an interval. */
function batteryBlankCiReasons(rows){
  const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
  const out = { baseRate: false, nonePublished: false };
  if (!db) return out;
  rows.forEach(row => {
    const entry = row.group && db[row.group] ? db[row.group][row.name] : null;
    if (!entry || typeof entry !== 'object') return;
    const type = rowScoreType(row);
    const rel = resolveCiReliability(entry, BATTERY_METRIC_SD[type], batteryPatientAge());
    if (!rel || !rel.none) return;
    if (rel.basis === 'base rate — no interval') out.baseRate = true;
    else out.nonePublished = true;
  });
  return out;
}

/* Hard limits imposed by the score scale itself, used to keep an interval from
   printing a value the scale cannot produce. Only scaled scores have a
   universal ceiling (Wechsler subtests run 1-19); "standard" and T are generic
   metrics whose usable range varies by instrument, so no ceiling is asserted
   for them, and z is unbounded. */
const BATTERY_SCORE_BOUNDS = { scaled: { min: 1, max: 19 } };

/* `age` is optional and defaults to whatever the age box holds; passing it
   explicitly is what lets the interval be exercised without a DOM. */
function getBatteryCiHtml(ss, row, level, age){
  if (!Number.isFinite(ss)) return '';
  // Resolved first: the displayed metric now decides which SD the SEM uses.
  const type   = rowScoreType(row);
  const effAge = age !== undefined ? age : batteryPatientAge();
  const rel    = getBatteryRowReliability(row, type, effAge);
  if (!rel) return '';
  const zMult  = level === '90' ? 1.645 : 1.960;
  const sem    = rel.sd * Math.sqrt(1 - rel.r);
  const loRaw  = ss - zMult * sem;
  const hiRaw  = ss + zMult * sem;
  /* z and raw share one treatment: 1dp, no floor and no ceiling.
     The interval itself is valid for a raw score — sd1 is in raw units, so
     SD x sqrt(1-r) is a raw-unit SEM — but the integer floor of 1 applied
     below would be wrong for the many raw measures whose scale starts at 0
     (intrusions, perseverations, false positives all have normative means
     under 6), and no raw ceiling can be asserted without knowing each
     measure's maximum. */
  if (type === 'z' || type === 'raw'){
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
  /* Base-rate rows are entered in the RAW column (see batteryBaseRateValue).
     Two consequences handled here, before the rows render:
     - a value typed in the Score column while that was the entry route (or
       restored from an older session) is migrated into the raw field, once,
       so nothing already entered stops scoring;
     - the raw column must be VISIBLE while such rows exist, whatever the
       Show Raw toggle says — a hidden entry column is an unusable row. The
       .raw-forced class outranks .raw-hidden in styles.css by specificity. */
  let hasBaseRateRows = false;
  const rebandAge = batteryPatientAge();
  batteryRows.forEach(r => {
    if (!batteryBaseRateEntry(r)) return;
    hasBaseRateRows = true;
    if ((r.raw === '' || r.raw == null) && r.score !== '' && r.score != null){
      r.raw = r.score; r.score = '';
    }
    /* AUTO-REBAND. The dropdown offers one entry per base-rate family and the
       band is a mechanical function of the required age, so when the age and
       the row's band disagree, the row is re-pointed at the sibling band
       containing the age — silently correct beats a complaint the clinician
       has to resolve by re-adding the row. Only where the target band
       actually publishes the measure; otherwise the row keeps its group and
       the out-of-band hint states the refusal. */
    if (rebandAge != null){
      const band = ageBandRange(r.group);
      if (band && (rebandAge < band.lo || rebandAge > band.hi)){
        const target = baseRateGroupForAge(r.group, r.name, rebandAge);
        if (target) r.group = target;
      }
    }
  });
  document.getElementById('bat-table').classList.toggle('raw-forced', hasBaseRateRows);
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
      /* A base-rate group is entered as raw spans and scored by lookup — its
         entries' means sit near the scaled range, so the inferred metric would
         badge it "Scaled Score", over a section whose Score boxes are disabled.
         Name what is actually entered instead. */
      const stLabel = batteryGroupIsBaseRate(gKey) ? 'Raw score'
        : groupTypes.size > 1 ? 'Mixed' : scoreTypeLabel([...groupTypes][0] || r.scoreType || inferScoreType(r.group));
      // Custom tests get an editable name; database families keep their fixed one.
      const nameHtml = r.groupKey
        ? `<input class="group-name-input" data-group-rename="${escapeAttr(gKey)}" value="${escapeAttr(r.group)}" placeholder="Name this test" aria-label="Test name" autocomplete="off">`
        : escapeHtml(batteryGroupLabel(gKey));
      ghr.innerHTML = `<td colspan="8">${nameHtml} <span class="type-badge">· ${stLabel}</span>` +
        `<button class="group-add-row" data-add-to-group="${escapeAttr(gKey)}" title="Add a row to this test" aria-label="Add a row to this test">+</button>` +
        `<button class="group-remove" data-rm-group="${escapeAttr(gKey)}" title="Remove group">×</button></td>`;
      tbody.appendChild(ghr);
      /* A base-rate section reports a different quantity from the rest of the
         table, so it states its own column heading rather than borrowing the
         one in the thead. Only the affected column is relabelled; the others
         are left blank so the row reads as a heading for this section, not as
         a second table header. */
      if (batteryGroupIsBaseRate(gKey)){
        const sub = document.createElement('tr');
        sub.className = 'group-col-header';
        sub.innerHTML = `<td></td><td></td><td class="bat-raw-cell"></td><td></td>` +
          `<td class="bat-ci-cell"></td>` +
          `<td class="group-col-label" title="Percentage of the normative sample obtaining this score or higher, as published (WAIS-IV Manual, Tables C.4–C.5). A HIGH base rate means a common, and therefore lower, score.">${BAT_BASERATE_LABEL}</td>` +
          `<td></td><td></td>`;
        tbody.appendChild(sub);
      }
      lastGroup = gKey;
    } else if (!gKey){
      lastGroup = null;
    }
    const rowType = rowScoreType(r);
    const z = toZ(r.score, rowType);
    const pctCellVal = batteryRowPctCell(r);
    const pct = pctCellVal ? (pctCellVal.hint ? batteryAgeHintHtml(pctCellVal.hint) : pctCellVal.text) : '';
    const details = batteryClassificationDetails(r, cls);
    const ss = parseFloat(r.score);
    const ciHtml = ciLevel !== 'off' ? getBatteryCiHtml(ss, r, ciLevel) : '';
    const tr = document.createElement('tr');
    if (gKey) tr.className = 'in-group';
    /* A base-rate row has no metric score — the raw span IS the entry — so its
       Score box is disabled rather than left inviting a number that has no
       meaning here. Disabled, not removed: the cell keeps the column aligned.
       Its inferred-metric tag is suppressed for the same reason: "(Scaled)"
       beside a disabled Score box would claim a metric the row does not have. */
    const isBr = !!batteryBaseRateEntry(r);
    const abbr = isBr ? '' : scoreTypeAbbr(rowType);
    const typeTag = abbr ? `<span class="bat-row-type-tag">(${abbr})</span>` : '';
    const scoreCell = isBr
      ? `<td class="bat-score-na-cell"><input type="number" disabled class="bat-score-na" title="Raw-score measure: enter the span in the Raw Score column — the base rate is read from the published table." aria-label="Not applicable — raw-score measure"></td>`
      : `<td><input type="number" step="any" data-r="${i}" data-f="score" value="${escapeAttr(r.score)}"></td>`;
    tr.innerHTML = `
      <td class="row-num">${i+1}${typeTag}</td>
      <td><input type="text" data-r="${i}" data-f="name" value="${escapeAttr(r.name)}" placeholder="Subtest name"></td>
      <td class="bat-raw-cell"><input type="number" step="any" data-r="${i}" data-f="raw" value="${escapeAttr(r.raw)}"></td>
      ${scoreCell}
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
      /* Each row's quantity is fixed by the measure, not by what else is in the
         table, so only the edited row needs rewriting — no cross-row refresh,
         and no heading to keep in step. */
      const cellVal = batteryRowPctCell(batteryRows[i]);
      if (cellVal && cellVal.hint) pctCell.innerHTML = batteryAgeHintHtml(cellVal.hint);
      else pctCell.textContent = cellVal ? cellVal.text : '';
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
  /* Whether the age is being read depends on the rows and the CI level, both of
     which can change without the age itself being touched — adding a D-KEFS row
     or switching CI off. Hooked here so every one of those paths updates the
     pip AND the prompt; renderBattery() is what they all end in. */
  refreshPatientAgeIndicator();
  refreshBatteryAgePrompt();
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
      /* A section reporting a different quantity from the rest of the table
         states its own column heading. Driven by the column definitions rather
         than hard-coded here, so this stays generic. */
      const subLabels = visible.map(c => (c.groupLabel ? c.groupLabel(group) : ''));
      if (subLabels.some(Boolean)){
        /* Section name and its column labels share ONE row. Given as separate
           rows they left a visible gap above the first data row, and the label
           sat closer to the data than to the heading it belongs to. The name
           takes the first cell; every other cell carries its own label, which
           is blank for all but the relabelled column. */
        body += `<tr class="apa-group apa-group-labelled">` + visible.map((c, i) => {
          const tdCls = `${c.num ? 'num ' : ''}col-${c.key}`.trim();
          // The name wins the first cell; a groupLabel on column 0 would be
          // shadowed, which is correct — the section has to be named.
          if (i === 0) return `<td class="${tdCls}">${escapeHtml(groupText)}</td>`;
          return `<td class="${tdCls} apa-group-col-label">${escapeHtml(subLabels[i] || '')}</td>`;
        }).join('') + `</tr>`;
      } else {
        body += `<tr class="apa-group"><td colspan="${visible.length}">${escapeHtml(groupText)}</td></tr>`;
      }
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
    /* Raw rows print a score but no percentile or classification, because a raw
       score carries no metric to derive them from. Without this line the blank
       cells read as an oversight rather than as the deliberate refusal they
       are. */
    ctx.hasRaw
      ? 'Raw-score measures are reported as obtained; percentile and classification are not derived, as these measures are not on a standardised metric.'
      : '',
    /* Names the source, because these percentiles come from a published table
       rather than from a metric conversion, and says which direction they run
       in — the manual tabulates "or higher", the column reports "below". */
    /* One sentence, because the quantity no longer depends on what else is in
       the table — a base-rate section always reports a base rate, and says so
       in its own column heading. Defines the direction, which is the opposite
       of the percentile column above it. */
    ctx.hasBaseRates
      ? 'Base rate = percentage of the normative sample obtaining the same score or higher (WAIS-IV Administration and Scoring Manual, Tables C.4–C.5). A higher base rate indicates a more common, and therefore lower, score.'
      : '',
    /* Without this the em-dash in the classification column reads as missing
       data rather than as a deliberate refusal. */
    /* Without this the pairing of a high percentile with a low classification
       reads as a contradiction rather than as the convention it is. */
    ctx.hasHigherIsWorse
      ? 'On error measures (perseverations, intrusions, false positives) a higher score indicates more errors. Percentiles are reported as obtained; classifications describe performance, so a high percentile corresponds to a low classification.'
      : '',
    // States the BASIS, not just the level; without it a reader comparing
    // against the manual has no way to know why the numbers differ. Names the
    // SD as well as the reliability, because both terms are choices and the
    // pairing of them is the methodological claim — see
    // getBatteryRowReliability. Says "retest or alternate-form" rather than
    // "test-retest": CVLT-3 publishes alternate-form coefficients, so the
    // unconditional wording misdescribed every CVLT-3 table. The comparison
    // with the manual is conditional for the same reason — not every
    // publisher reports internal consistency at all. Full rationale is in
    // Methods & References.
    ctx.ciLevel && ctx.ciLevel !== 'off'
      ? `Confidence intervals are ${ctx.ciLevel}%, calculated as the obtained score ± z × SEM, where SEM = SD × √(1 − r) on the normative standard deviation of the reported metric. The reliability r is the coefficient each test's manual uses for its own published intervals: internal consistency where the publisher reports one, otherwise the test–retest coefficient.`
      : '',
    /* THE SENTENCE ABOVE IS FALSE WITHOUT THIS ONE, on any table holding a
       measure whose publisher tabulates reliability BY AGE BAND while no age
       was entered.

       D-KEFS (original) is the case: it publishes rInternalByAge with no
       all-ages average, so at a blank age resolveCiReliability falls through
       to the retest branch — and the sentence above then claims internal
       consistency was used "where the publisher reports one" when the
       publisher does report one and the app did not use it. Verified: D-KEFS
       Tower Total Achievement takes r .44 at a blank age against .72 at 45.

       This is the same defect the comment above records fixing once already
       for CVLT-3, recurring by a different route. The predicate is exact:
       'retest, uncorrected · by age' can only arise when the entry HAS banded
       coefficients (the suffix is added only for those) and the retest r was
       used anyway. batteryBasisPresent asks the shipped resolver at the age
       actually in force, so it cannot disagree with the printed intervals. */
    ctx.ageBandNotSelected
      ? 'One or more measures here publish their reliability by age band. No patient age was entered, so the published all-ages or test–retest coefficient was used for those measures instead; entering an age narrows their intervals.'
      : '',
    /* A BLANK INTERVAL MUST SAY WHY IT IS BLANK. Every other deliberate blank
       in this table already carries a sentence — the raw percentile above and
       the higher-is-worse pairing — on the stated ground that an absent value
       reads as an oversight rather than as the refusal it is. The CI column
       was the one left out: 55 measures reachable from Score Tables print no
       interval, for two quite different reasons that must not be merged.

       Sourced from resolveCiReliability rather than getBatteryRowReliability,
       because the latter collapses both to null and the distinction is the
       whole point of the sentence. */
    ctx.blankCiBaseRate
      ? 'Longest-span measures are scored from a published base-rate table rather than by conversion, and have no reliability coefficient, so no confidence interval is shown for them.'
      : '',
    ctx.blankCiNonePublished
      ? 'Where a measure is absent from its manual\'s reliability table, no coefficient is available and no confidence interval is shown for it.'
      : '',
    /* NO SENTENCE FOR THE UNCORRECTED PAIRING. It used to have one, and it has
       been moved wholesale to Methods & References.

       The test is what an exported table can be MISREAD as without it. Here,
       nothing: the retest coefficient beside a normative SD is what those
       manuals themselves do, so the interval matches the published one and a
       reader checking it finds agreement. It is a caveat about precision, which
       belongs with the other caveats, not a departure needing declaration.

       The derived coefficient below fails that same test, which is why it keeps
       its sentence. */
    /* THE ONLY COEFFICIENT IN THIS APP THAT IS NOT A PUBLISHED NUMBER, so a
       table built on it has to say so — a reader cross-checking the manual will
       not find this value there, and would reasonably conclude the table is
       wrong. Kept in the note for that reason alone, when everything else about
       the CI method moved to Methods & References: this is the one thing the
       exported artefact cannot be left to imply. The derivation itself is on
       the Methods page; what survives here is the fact and its consequence.

       Stated only when a measure in THIS table actually used it — the control
       is a page setting, so it can be on while every measure present takes a
       published coefficient and nothing was in fact corrected. Same rule the
       age sentence below follows. */
    ctx.hasDerivedR
      ? 'Coefficients for measures whose manual publishes no normative-sample reliability have been corrected to the normative sample\'s variability for this table. These are not the values those manuals print, so the intervals concerned differ from the published ones.'
      : '',
    /* Which age band the coefficients were drawn from. Without this the
       interval cannot be reproduced from the manual, because for these
       measures the published reliability changes with age.

       THE TAIL NO LONGER NAMES THE FALLBACK COEFFICIENT, and that is a fix
       rather than a trim. It used to end "otherwise the total-sample retest
       coefficient", which the reliability-basis control can falsify: D-KEFS
       (original) publishes rInternalByAge with no all-ages average, so an
       out-of-range age lands on the retest branch, and with the correction on
       that branch no longer yields the manual's figure — Tower's stored r is
       .44 against a corrected .359. Naming the fallback needed a conditional
       clause to stay true; saying only that a fallback happened needs none, and
       says everything a reader must know to reproduce the row. Which
       coefficient each fallback resolves to is on the Methods page. */
    ctx.ciAge != null
      ? `For measures whose published reliability is tabulated by age, coefficients are those for age ${ctx.ciAge}, or the publisher's all-ages figure where that age falls outside a measure's normed range.`
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
    /* Index mode has no divisor for a raw-score row, so those rows compute
       nothing. Say so, and say where they DO work, rather than leaving a row
       of blanks with no explanation. */
    ctx.hasRawInIndexMode
      ? 'Raw-score measures are not scored in index mode, which divides by the metric’s SD; use raw mode, which divides by the normative SD entered for each measure.'
      : '',
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
    'Base rate = published % at or below this discrepancy (ToPF-UK manual, negative discrepancies only). The manual derives these from a normal model with SD = SEE rather than from observed standardisation-sample frequencies.'
  ],
  /* Performance Validity page — one note under the single combined table.

     EVERY CLAUSE IS CONDITIONAL ON WHAT THE TABLE ACTUALLY HOLDS. The note
     had grown a sentence per feature until it ran to a paragraph on a
     table of one measure, citing five papers the reader could not see the
     scores for. It now cites only the measures present, states the
     independence rule once for whichever pairs are present rather than
     once per pair, and explains the dashes in one clause covering both
     reasons an accuracy pair can be absent (an index published as a base
     rate, and the Effort Scale's AUC).

     The measure names in the table are written out in full, so the note
     carries no abbreviation roster — expanding "EI" for a table that
     never says "EI" was pure length. */
  'pvt': ctx => {
    const sources = [];
    if (ctx.hasEi)   sources.push('Effort Index (Silverberg et al., 2007)');
    if (ctx.hasEs)   sources.push('Effort Scale (Novitski et al., 2012)');
    if (ctx.hasRds)  sources.push('Reliable Digit Span (Greiffenstein et al., 1994)');
    if (ctx.hasDs)   sources.push('Digit Span indices (Iverson & Tulsky, 2003; Axelrod et al., 2006, WAIS-III)');
    if (ctx.hasRey)  sources.push('Rey 15-Item (Boone et al., 2002)');
    /* The source entry follows the basis: with a borrowed cut-off in force
       the table rests on two publications, not one, and the Sources line is
       where a reader looks for that. */
    if (ctx.hasCvlt3) sources.push(ctx.cvlt3Borrowed
      ? `CVLT-3 Forced Choice (base rates Delis et al., 2017, Tables D.13–D.15; cut-off and accuracy ${ctx.cvlt3Cite}, CVLT-II)`
      : 'CVLT-3 Forced Choice (Delis et al., 2017, Tables D.13–D.15)');
    if (ctx.hasTomm) sources.push('TOMM (Tombaugh, 1996; cut-offs Martin et al., 2020)');
    const shared = [];
    if (ctx.bothRbans)     shared.push('the two RBANS indices');
    if (ctx.bothDigitSpan) shared.push('the digit-span indices');
    return [
      /* The mirror drops the citations: the tab strip and the on-page
         references state every measure and source in full. The exported
         note keeps them — the licensed onScreen difference. */
      (ctx.onScreen || !sources.length) ? '' : `Sources: ${sources.join('; ')}.`,
      '"Fail" = score beyond the published cut-off, not a determination of invalidity; probable invalidity is conventionally supported by failure of at least two independent indicators (Larrabee, 2014).',
      ctx.hasDashes
        ? 'Sensitivity and specificity are the published values at the applied cut-off; a dash marks an index published as a base rate or an AUC rather than as a pair.'
        : 'Sensitivity and specificity are the published values at the applied cut-off.',
      shared.length ? `Indices sharing a subtest count as one indicator: ${shared.join(' and ')}.` : '',
      /* The CVLT-3 is the one measure here with no published cut-off, so an
         exported table must say where its threshold came from — otherwise
         the Cut-off column implies the manual printed one. */
      ctx.hasCvlt3 && !ctx.cvlt3Borrowed
        ? 'The CVLT-3 Forced Choice manual publishes base rates by age band, not a cut-off; its threshold is the rarest score whose published base rate is at or below the stated per-test false-positive criterion, so it varies with age.'
        : '',
      ctx.cvlt3Borrowed
        ? `The CVLT-3 manual publishes no cut-off for Forced Choice Recognition; the cut-off and accuracy applied here are CVLT-II figures (${ctx.cvlt3Cite}), the trial being structurally identical across editions. The base rate shown is the CVLT-3 manual's own.`
        : '',
      /* A cut-off derived on one edition of a test does not automatically
         transfer to another, and the exported table is where a reviewer
         asks. Emitted only for the measures actually in the table. */
      ctx.versionCaveats && ctx.versionCaveats.length
        ? `Instrument versions: ${ctx.versionCaveats.join(' ')}`
        : '',
      ctx.esGated
        ? 'The Effort Scale is not computed because its screening gate was not met (Novitski et al., 2012).'
        : '',
      ctx.hasTomm
        ? 'Specificity falls in dementia and severe impairment; traditional TOMM cut-offs should not be interpreted there.'
        : 'Specificity falls in dementia and severe impairment.'
    ];
  },
  'pre-opiepredict': () => [
    'OPIE-4 prorated scores are predicted from age and sex with Vocabulary and/or Matrix Reasoning.',
    '<span class="uk-caution-red">Illustrative only in a UK context as this is derived from US regression equations. The numbers should not be considered to be accurate in a UK context.</span> The published equations also use US education, ethnicity and region terms which are not applied, so every patient is scored at the US reference category (12th-grade high-school graduate, not African-American, not resident in the US West). These estimates would likely run high for patients who left school early and low for graduates.'
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
   Only context-free notes can be mirrored this way.

   `onScreen` is the ONE licensed difference between the two renderings: a note
   may drop a sentence that the surrounding page already states in full. It may
   not add, soften or reword one — the exported note stays the superset. */
function renderStaticApaNotes(){
  document.querySelectorAll('[data-apa-note]').forEach(el => {
    const build = APA_NOTES[el.dataset.apaNote];
    if (!build) return;
    const body = build({ onScreen: true }).filter(s => s && s.trim()).join(' ');
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
    /* A base-rate row's entered value lives in the raw field; the export's
       Score column carries it so the span is printed whether or not the Raw
       column is toggled on. */
    { key:'score',          label:headerLabel,       group:'Scores', num:true,  render:r => escapeHtml(batteryBaseRateEntry(r) ? (r.raw || '') : (r.score || '')) },
    { key:'ci',             label:ciLabel,           group:'Scores', num:true,  defaultVisible:ciLevel !== 'off', render:r => { const ss = parseFloat(r.score); return ciLevel !== 'off' ? getBatteryCiHtml(ss, r, ciLevel) : ''; }},
    { key:'percentile',     label:BAT_PCT_LABEL,     group:'Scores', num:true,
      /* Sections whose rows report a base rate relabel this column for
         themselves; everything else reads under the main heading. */
      groupLabel: g => (batteryGroupIsBaseRate(g) ? BAT_BASERATE_LABEL : ''),
      render:r => { const c = batteryRowPctCell(r); return c ? c.text : ''; }},
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
      hasRaw: valid.some(r => rowScoreType(r) === 'raw' && !batteryBaseRateEntry(r)),
      hasBaseRates: valid.some(r => batteryBaseRateEntry(r) && r.raw !== '' && !isNaN(r.raw)),
      hasHigherIsWorse: valid.some(r => r.higherIsWorse && r.score !== '' && !isNaN(r.score)),
      ciLevel,
      /* Reported ONLY when it actually changed a coefficient. A single field
         silently narrowing or widening every interval in the table is exactly
         the kind of thing a reader must be able to check, so the note names
         the age it rests on — but claiming an age was used on a table where
         no measure reads one would be its own misstatement. */
      ciAge: (ciLevel !== 'off' && batteryPatientAge() !== null
              && valid.some(r => batteryRowUsesAgeBand(r))) ? batteryPatientAge() : null,
      /* A publisher tabulates reliability by age band and the band was not
         selected — no age entered, or one outside the normed range. The suffix
         is added only to banded entries, so this prefix cannot match anything
         else. Without it the CI-method sentence above is false. */
      ageBandNotSelected: ciLevel !== 'off' && batteryBasisPresent(valid, 'retest, uncorrected · by age'),
      /* Why a blank interval is blank. Two reasons, never merged. */
      blankCiBaseRate:      ciLevel !== 'off' && batteryBlankCiReasons(valid).baseRate,
      blankCiNonePublished: ciLevel !== 'off' && batteryBlankCiReasons(valid).nonePublished,
      /* The corrected reading is a departure from every manual cited in this
         table, so it is the one thing here a reader most needs told. Asked the
         same way and for the same reason: the control being on does not mean a
         coefficient was actually corrected. */
      hasDerivedR: ciLevel !== 'off' && batteryBasisPresent(valid, 'retest, corrected here'),
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
    showToast(`${family} has no tests yet - add some in Data`, true);
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
  /* higherIsWorse rides along on the row so the classification cell can see it
     without another database lookup. On these measures a high score means MORE
     errors, so a merit label like "Very Superior" would assert the opposite of
     what the number means — see batteryClassificationDetails. */
  const newRows = missing.map(name => ({
    name, raw:'', score:'', group:family,
    scoreType: inferScoreTypeForSubtest(family, name, db[family][name]),
    higherIsWorse: !!(db[family][name] && db[family][name].higherIsWorse)
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
  /* No popover re-arming here any more. It used to be needed because the
     popover fired once per patient and only a reset could re-enable it; the
     trigger is now simply "interval on, no age stored", so switching CI off
     and on asks again by itself and clearing the table changes nothing about
     whether the question applies. */
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
/* The reliability basis drives the table AND the Data page, which reports what
   the app's intervals rest on — leaving that page on the old reading would make
   it advertise a coefficient the table is no longer using. Guarded because
   renderDbList needs the Data page's markup, which is present in index.html but
   not in every harness that loads app.js. */
document.getElementById('bat-ci-basis').addEventListener('change', () => {
  renderBattery();
  if (typeof renderDbList === 'function' && document.getElementById('db-tbody')) renderDbList();
});
/* The master patient age. Bound here rather than beside the nav code because
   everything it recalculates is on this page or the premorbid one.

   'input', not 'change': waiting for blur leaves a stale age visibly paired
   with fresh numbers. Mirrors into #pre-age and recalculates the premorbid
   page, because the two boxes hold one patient's age — see the note on
   PATIENT_AGE_INPUTS. The element is in the top bar, which is markup, so it
   exists by the time this runs; the optional chaining is belt-and-braces. */
document.getElementById('patient-age')?.addEventListener('input', () => {
  syncPatientAge('patient-age');
  renderBattery();
  if (typeof calcPremorbid === 'function'){ calcPremorbid(); calcPredict(); calcOpiePredict(); }
  refreshPatientAgeIndicator();
  rebuildFamilyListsForAge();
  /* The Data page reports which coefficient an interval rests on, and that
     answer moves with the age for every banded measure. Left un-rendered it
     would keep showing the previous patient's reading — the same staleness the
     reliability-basis control is re-rendered for. Guarded because renderDbList
     needs that page's markup. */
  if (typeof renderDbList === 'function' && document.getElementById('db-tbody')) renderDbList();
});
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
/* Returns undefined for 'raw' (and anything unrecognised) BY DESIGN — a raw
   score has no metric SD to divide by. sdiComputeChange must therefore treat a
   missing unit as "cannot compute" rather than dividing and producing NaN, or
   worse, having a unit invented for it. */
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
  /* Index mode divides by the metric's SD, which only exists for a
     standardised score. Raw-score rows (metric:'raw' in normDB) have no such
     unit, so they are refused here rather than scored against a borrowed one.
     Dividing RBANS List Recognition — a raw score with SD 0.8 — by the scaled
     unit of 3 understated every change by a factor of 3.75, turning a 2-point
     drop (2.5 SD, clearly reliable) into 0.67 SD, "not significant".
     Raw mode still handles these correctly: it uses the row's own SD, which
     autofill populates from sd1. */
  const unit = sdiSdUnit(sdiRowScoreType(r));
  if (!Number.isFinite(unit) || unit <= 0) return null;
  return (parseFloat(r.t2) - parseFloat(r.t1)) / unit;
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
    ${apaNoteHtml('sdi', {
      thresholdLabel: cvDesc,
      mixedTypes: !raw && new Set(named.map(sdiRowScoreType)).size > 1,
      hasRawInIndexMode: !raw && named.some(r => sdiRowScoreType(r) === 'raw')
    })}
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
    /* sdiNormSd returns sd1, which is in whatever units m1/sd1 are — raw, for
       the CVLT-C entries. Prefilling it would be wrong for a measure the
       clinician enters in a DIFFERENT metric: they would be typing T-scores
       into columns divided by a raw-unit SD. Where the entry metric differs
       from the stored one, leave the SD for them to supply. Index mode is
       unaffected: it divides by the entry metric's own unit. */
    const entryMetricDiffers = SCORE_METRICS.has(p && p.reportedAs) && p.reportedAs !== p.metric;
    sdiRows.push({
      name, t1:'', t2:'', group:family,
      sd: (raw && !entryMetricDiffers) ? sdiNormSd(p) : '',
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
  /* Alternate-form entries are reliable-change (RCI) only; keep them out of SDI.
     Single-administration families (base rates, no reliability) are excluded
     for the opposite reason — they have no second administration to compare. */
  const families = Object.keys(db).sort()
    .filter(f => !isAltFormFamily(f) && !isSingleAdministrationFamily(f));
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
/* The four methods share ONE row set, so entering scores on one of them fills
   in the other three as a side effect. The working report therefore collects a
   method only once it has been ACCEPTED, and entering data IS the acceptance —
   this marks it, so the tab the clinician is actually working in collects
   silently, exactly as it did before. The other three offer themselves through
   the "Add to report" control on their own APA toolbar.

   Deliberately NOT called from the view settings: the CI level, the corrected-r
   toggle and the column pickers all re-render the APA table without a score
   being touched, and re-selecting the CI level on a method never opened was
   precisely how a third table used to appear in the report. */
function rciMarkMethodUsed(method){
  if (!rciState[method]) return;
  if (typeof ReportBundle === 'undefined' || !ReportBundle.acceptSource) return;
  ReportBundle.acceptSource(`${method}-apa`, { capture: false });
}
function rciAddRow(method){ rciMarkMethodUsed(method); rciState[method].rows.push(newRciRow(method)); renderRci(method); }
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
      rciMarkMethodUsed(target);
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
      rciMarkMethodUsed(m);
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
//
// That last clause is only true for families scored by converting a metric —
// the age band picks which normative sample the m1/sd1 came from, and Score
// Tables never touches those. It is NOT true for a family scored from a
// published base-rate table, where the band IS the lookup. Collapsing WAIS-IV
// longest span would have scored every patient against All Ages: Longest Digit
// Span Backward at a span of 4 is the 22nd percentile at 20-24 and the 59th at
// 85-90, so up to 37 percentile points would have gone silently wrong.
//
// The second case is a band that is not a band at all but a SEPARATE BATTERY.
// WMS-IV publishes an Adult battery (16-69) and an Older Adult battery (65-90)
// with different subtests — no Designs, no Spatial Addition, no Visual Working
// Memory Index in the older one — and different reliability for the SAME
// measure. Age cannot choose between them, because 65-69 appears in both
// (Logical Memory I is .79 in the adult battery and .83 in the older adult):
// the clinician chose which to administer, so the clinician must choose here.
// Collapsing them showed an 80-year-old the adult battery's measure list, and
// its coefficients — 10 of the 12 shared measures printed a wrong interval.
//
// Declared per entry rather than inferred. A rule like "the two groups hold
// different measures" would have worked today and broken silently the first
// time a family's bands diverged for some other reason.
function familyScoredByAgeBand(name){
  const fam = normDB[name];
  if (!fam) return false;
  return Object.values(fam).some(e => e && typeof e === 'object' && (e.baseRates || e.separateBattery));
}

/* Base-rate families ONLY: every entry scored by published lookup. Distinct
   from familyScoredByAgeBand, which also matches separateBattery (WMS-IV) —
   there the band is a CLINICAL choice of battery (§30) and must stay
   selectable; here it is a mechanical function of the patient age. */
function familyGroupIsBaseRate(name){
  const fam = normDB[name];
  if (!fam) return false;
  const entries = Object.values(fam).filter(e => e && typeof e === 'object');
  return entries.length > 0 && entries.every(e => e.baseRates);
}

/* The banded sibling group that contains this age AND publishes this measure.
   Per MEASURE, not per family: Longest Letter-Number Sequence is published
   only to 65-69 while its siblings run to 85-90, so a 70-year-old has a band
   for the digit spans and none for LNS. */
function baseRateGroupForAge(group, name, age){
  if (age == null) return null;
  const base = familyBaseName(group);
  const db = getMergedDB();
  for (const g of Object.keys(db)){
    if (!hasAgeBandSuffix(g) || familyBaseName(g) !== base) continue;
    const band = ageBandRange(g);
    if (!band || age < band.lo || age > band.hi) continue;
    const e = db[g] && db[g][name];
    if (e && e.baseRates) return g;
  }
  return null;
}

/* ── AGE-BAND FILTERING OF THE FAMILY DROPDOWNS ──────────────────────────────

   With one patient age now on screen everywhere, the dropdowns can stop
   offering every band of every family at once. The rule is deliberately
   conservative, because the band means different things on different pages:
   on Score Tables it is usually a metric conversion that does not change the
   table, but for a base-rate family it IS the lookup; on Change Analysis and
   the SD Index it selects the retest sample the r was measured in.

   THREE CLAUSES, per family:

     1. every band whose span contains the age — PLURAL, because WMS-IV
        publishes overlapping 16-69 and 65-90 bands, so a 67-year-old is
        legitimately in both and picking one would hide a real choice;
     2. plus "· All Ages" wherever the family publishes one, so the clinician
        keeps the larger-N / stronger-r option that comboAgeBandNoteHtml
        describes. That choice is a clinical judgement, not something an age
        should make silently — which is why this is a filter, not a jump;
     3. and if 1 and 2 both come back empty, EVERY member, unchanged.

   CLAUSE 3 IS THE SAFETY RAIL, and it is not hypothetical. Eight of the 28
   numerically banded families publish no All Ages entry at all, and for those
   an age can match nothing: CVLT-C Subtests is banded by the three normative
   sample ages (Age 8 / Age 12 / Age 16), so 88 of the 91 ages from 5 to 95
   match none of them. Without clause 3, a clinician testing a 10-year-old
   would open Quick Add and find CVLT-C simply gone — which reads as the app
   being broken rather than as a filter doing its job, the same failure mode as
   a blank age emptying the CI column. RBANS, RBANS Form B, CVLT-3 and WMS-IV
   are the others.

   THE INVARIANT: an age may quieten a family's list. It may never remove the
   family. check.js §29 asserts that over every family at every age 5-95.

   The escape hatch is real rather than decorative: whenever the list has been
   narrowed it says so and offers "Show all bands", because a silently shorter
   list is indistinguishable from missing data. */
let familyListShowAllBands = false;

/* "· Ages 20-49" -> {lo:20, hi:49};  "· Age 8" -> {lo:8, hi:8}.
   null for "· All Ages" and for names carrying no band, both of which are
   handled by name above rather than by range. */
function ageBandRange(f){
  const m = String(f).match(/·\s*Ages?\s+(\d+)\s*(?:[–-]\s*(\d+))?\s*$/i);
  if (!m) return null;
  const lo = +m[1];
  return { lo, hi: m[2] != null ? +m[2] : lo };
}
function isAllAgesFamily(f){ return /·\s*All\s+Ages\s*$/i.test(f); }

function familyMembersForAge(members, age){
  if (age == null || familyListShowAllBands) return members;
  const keep = members.filter(f => {
    if (isAllAgesFamily(f)) return true;      // clause 2
    const r = ageBandRange(f);
    if (!r) return true;                      // unbanded — not ours to filter
    return age >= r.lo && age <= r.hi;        // clause 1
  });
  return keep.length ? keep : members;        // clause 3
}

/* Shown only when something was actually hidden, so a list that was already
   as short as it can be does not claim to have been filtered. */
function comboAgeFilterNoteHtml(age, narrowed){
  if (age == null) return '';
  if (familyListShowAllBands){
    return '<div class="combo-agefilter-note"><span>Showing <strong>all</strong> age bands.</span>'
         + `<button type="button" class="combo-agefilter-toggle" data-combo-agefilter="filter">Filter to age ${escapeHtml(String(age))}</button></div>`;
  }
  if (!narrowed) return '';
  return `<div class="combo-agefilter-note"><span>Showing bands for <strong>age ${escapeHtml(String(age))}</strong>, plus All Ages where published.</span>`
       + '<button type="button" class="combo-agefilter-toggle" data-combo-agefilter="all">Show all bands</button></div>';
}
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

  /* Read once for the whole list rather than per group, so a list cannot be
     built half against one age and half against another. */
  const age = (typeof patientAge === 'function') ? patientAge() : null;
  let narrowed = false;

  let html = '';
  order.forEach(base => {
    const allMembers = groups[base];
    const members = familyMembersForAge(allMembers, age);
    const groupKey = `grp:${base}`;
    /* A base-rate family collapses to ONE entry on the flat (Score Tables)
       list: its band is a mechanical function of the now-required patient
       age, so fourteen selectable bands were fourteen chances to pick the
       wrong one. The stored value is the band containing the age where one
       is known — with no age yet, the first band stands in and renderBattery
       auto-rebands the rows the moment the age arrives. Not marked as
       "narrowed": nothing selectable was hidden, and the show-all-bands note
       must not claim otherwise. WMS-IV does NOT come this way —
       familyGroupIsBaseRate is false for separateBattery groups, whose band
       choice is the clinician's (§30). */
    if (flat && !allMembers.some(isCustom) && allMembers.length > 0
        && hasAgeBandSuffix(allMembers[0]) && allMembers.every(familyGroupIsBaseRate)){
      const inBand = age != null
        ? allMembers.find(m => { const b = ageBandRange(m); return b && age >= b.lo && age <= b.hi; })
        : null;
      const bandText = inBand ? (inBand.match(/·\s*(.+)$/) || [,''])[1] : '';
      html += comboCheckboxItemHtml(inBand || allMembers[0], false, false, groupKey, base,
        inBand ? bandText : 'band set by age');
      return;
    }
    if (members.length < allMembers.length) narrowed = true;
    if (members.length === 1 && !hasAgeBandSuffix(members[0])){
      html += comboCheckboxItemHtml(members[0], isCustom(members[0]), false, groupKey);
    } else if (flat && !members.some(isCustom) && !members.some(familyScoredByAgeBand)){
      // Pick the "All Ages" canonical entry if present, else the first.
      const canon = members.find(m => /·\s*All\s+Ages\s*$/i.test(m)) || members[0];
      html += comboCheckboxItemHtml(canon, false, false, groupKey, base);
    } else if (flat && members.some(familyScoredByAgeBand)){
      /* Age band changes the score here, so it has to stay selectable even on
         the flat list. Rendered as the indented age-band pills the Change
         Analysis pages use, which is the existing pattern for "the band
         matters". */
      html += `<div class="combo-group-heading" data-group="${escapeAttr(groupKey)}">${escapeHtml(base)}</div>`;
      html += `<div class="combo-indented-row" data-group="${escapeAttr(groupKey)}">`;
      members.forEach(f => { html += comboCheckboxItemHtml(f, isCustom(f), true, groupKey); });
      html += `</div>`;
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
  return comboAgeFilterNoteHtml(age, narrowed) + comboOptionsHtml(html);
}
function rebuildAllFamilyLists(){
  document.querySelectorAll('.rci-family-list').forEach(list => populateFamilyList(list));
}

/* Every dropdown is built from the age at build time, so a changed age has to
   rebuild them. Guarded on the parsed value rather than the keystroke: typing
   "7" then "2" for 72 would otherwise rebuild three lists twice, once against
   the nonsense age 7. */
let lastFamilyListAge;
function rebuildFamilyListsForAge(force){
  const age = (typeof patientAge === 'function') ? patientAge() : null;
  if (!force && age === lastFamilyListAge) return;
  lastFamilyListAge = age;
  try { if (typeof rebuildBatteryFamilyList === 'function') rebuildBatteryFamilyList(); } catch(e){}
  try { if (typeof rebuildSdiFamilyList === 'function') rebuildSdiFamilyList(); } catch(e){}
  try { rebuildAllFamilyLists(); } catch(e){}
}

/* mousedown, not click, and preventDefault: the toggle sits inside an open
   combo list, and letting the press blur the search input would close the list
   out from under the very button being pressed. Same reason the footer buttons
   are given this treatment in wireMultiSelectFamilyList. */
document.addEventListener('mousedown', e => {
  if (e.target.closest('[data-combo-agefilter]')) e.preventDefault();
});
document.addEventListener('click', e => {
  const btn = e.target.closest('[data-combo-agefilter]');
  if (!btn) return;
  e.preventDefault();
  e.stopPropagation();
  familyListShowAllBands = (btn.dataset.comboAgefilter === 'all');
  rebuildFamilyListsForAge(true);
});
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
  /* Every Change Analysis method needs a second administration and a
     reliability coefficient. Families that publish neither cannot be scored
     here at all, so they are not offered. */
  let families = Object.keys(db).sort().filter(f => !isSingleAdministrationFamily(f));
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
  rciMarkMethodUsed(method);
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
/* A corrupt or hand-edited localStorage value must not take the app down.
   getCustom() is reached from getMergedDB(), which every autofill path and the
   Score Tables reliability lookup call, so a throw here is not a local failure.
   JSON.parse also returns non-objects for valid JSON ("null", "[1,2]", a bare
   string), and spreading those into the merged database produced families
   named "0" and "1". Anything that is not a plain object is discarded. */
function getCustom(){
  let parsed;
  try { parsed = JSON.parse(localStorage.getItem('customTests') || '{}'); }
  catch (e){ console.warn('customTests is not valid JSON; ignoring it', e); return {}; }
  if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) return {};
  return parsed;
}
function saveCustom(c){ localStorage.setItem('customTests', JSON.stringify(c)); }

/* The SAME rules the "Add subtest" form enforces (see the ct-add-subtest
   handler), applied to one entry and returned rather than shown as a toast, so
   the import path can reuse them.

   Import used to apply none of these. A hand-edited or shared file could put a
   negative sd1 into the clinical database, which passes Number.isFinite and
   then prints a confidence interval whose lower bound sits ABOVE its upper
   bound; or an r of 1.5, which makes sqrt(1-r) NaN. Export/import is the
   intended way to move a database between machines, so the file is not a
   trusted input.

   Returns { ok:true, entry } with a clean, numeric-typed entry, or
   { ok:false, reason }. */
function ctValidateEntry(raw){
  if (!raw || typeof raw !== 'object' || Array.isArray(raw)) return { ok:false, reason:'not an object' };
  const num = v => (typeof v === 'number' ? v : (typeof v === 'string' && v.trim() !== '' ? Number(v) : NaN));
  const m1 = num(raw.m1), sd1 = num(raw.sd1), m2 = num(raw.m2), sd2 = num(raw.sd2), r = num(raw.r);
  if (![m1, sd1, m2, sd2, r].every(Number.isFinite)) return { ok:false, reason:'m1, sd1, m2, sd2 and r must all be numbers' };
  if (!(r >= 0 && r < 1))   return { ok:false, reason:`r must be in [0, 1), got ${r}` };
  if (!(sd1 > 0 && sd2 > 0)) return { ok:false, reason:'SDs must be positive' };
  const entry = { m1, sd1, m2, sd2, r };
  if (raw.rCorrected != null && raw.rCorrected !== ''){
    const rc = num(raw.rCorrected);
    if (!(Number.isFinite(rc) && rc >= 0 && rc < 1)) return { ok:false, reason:`rCorrected must be in [0, 1), got ${raw.rCorrected}` };
    entry.rCorrected = rc;
  }
  if (raw.n != null && raw.n !== ''){
    const n = num(raw.n);
    if (!(Number.isFinite(n) && Math.floor(n) === n && n >= 3)) return { ok:false, reason:`N must be a whole number >= 3, got ${raw.n}` };
    entry.n = n;
  }
  /* A declared metric overrides the mean-based guess, so an unrecognised value
     must be refused rather than stored and silently ignored. '' / absent means
     "work it out from the mean", which is the default. */
  if (raw.metric != null && raw.metric !== ''){
    if (!SCORE_METRICS.has(raw.metric)){
      return { ok:false, reason:`metric must be one of ${[...SCORE_METRICS].join(', ')}, got "${raw.metric}"` };
    }
    entry.metric = raw.metric;
  }
  return { ok:true, entry };
}

/* Which instrument a group belongs to. D-KEFS Advanced must be tested before
   plain D-KEFS, the latter being a prefix of the former. */
function dbInstrumentOf(family, isCustom){
  if (isCustom) return 'Custom';
  if (/^D-KEFS Advanced/.test(family)) return 'D-KEFS Advanced';
  return family.split(' ')[0];
}

/* THE RELIABILITY THIS ROW'S CONFIDENCE INTERVAL IS ACTUALLY BUILT FROM.

   No longer a mirror of getBatteryRowReliability — both now call
   resolveCiReliability, so the page cannot advertise a coefficient the app
   does not use. What is left to get wrong is the ARGUMENTS, and there is one
   that matters: the range-restriction correction needs the normative SD, which
   depends on the metric the score is displayed in. Score Tables gets that from
   the row (rowScoreType), which for an auto-filled row is
   inferScoreTypeForSubtest's answer — so this page asks the same question of
   the same function rather than inventing a shortcut. check.js S32 drives both
   entry points over every entry and compares.

   No patient age, so this renders the BLANK-AGE answer, which is also the
   app's own fallback. Where a by-age table exists the basis says so, because
   the printed figure will move once an age is entered on Score Tables.

   IT DOES FOLLOW THE RELIABILITY-BASIS CONTROL, though that control lives on
   Score Tables. The whole point of this page is to say what the app's
   intervals rest on; showing the published reading while the table is set to
   corrected would be the exact defect the shared resolver exists to prevent.
   renderDbList is therefore re-run when the control changes. */
function dbReliabilityBasis(entry, family, name, correctRetest, age){
  const type = inferScoreTypeForSubtest(family, name, entry);
  const effCorr = correctRetest !== undefined ? correctRetest : batteryCiCorrectRetest();
  /* THE PATIENT AGE IS PART OF THE ANSWER, and passing `undefined` here was a
     defect rather than a simplification.

     These two columns exist to say what the app's intervals actually rest on.
     Hard-coding a blank age made them state the no-age reading whatever the
     patient age was — so with an age of 45 set, 74 of 209 entries reachable
     from Score Tables showed a different r from the one that produced the
     printed interval, and 13 named a different KIND of coefficient entirely.
     D-KEFS Tower Total Achievement read "retest, uncorrected · by age, r .44"
     while the table built its interval on "internal consistency · by age,
     r .72". That is precisely what CLAUDE.md warns of: the page telling a
     clinician their interval rests on a coefficient it does not.

     The page already follows the reliability-basis control from Score Tables
     for exactly this reason. The age is the same kind of state and is now
     followed the same way. */
  const effAge = age !== undefined ? age : batteryPatientAge();
  const rel = resolveCiReliability(entry, BATTERY_METRIC_SD[type], effAge, effCorr);
  return rel || { r:null, basis:'none published', none:true };
}

/* Sort state for the Norms Database table. Key is a column id, not an index,
   so reordering the columns cannot silently re-point the sort at a different
   one. null means "as the database is ordered", which groups an instrument's
   entries together and is the sensible default for browsing. */
let dbSort = { key: null, dir: 1 };

/* The columns, in order. `get` pulls the value used for BOTH display and
   sorting, so a column can never sort on something other than what it shows. */
const DB_COLUMNS = [
  { key:'instrument', label:'Instrument', text:true,  get:r => r.instrument },
  { key:'category',   label:'Category',   text:true,  get:r => r.category },
  { key:'band',       label:'Band',       text:true,  get:r => r.band },
  { key:'measure',    label:'Measure',    text:true,  get:r => r.name },
  { key:'m1',         label:'M₁',      dp:2, get:r => r.e.m1 },
  { key:'sd1',        label:'SD₁',     dp:2, get:r => r.e.sd1 },
  { key:'m2',         label:'M₂',      dp:2, get:r => r.e.m2 },
  { key:'sd2',        label:'SD₂',     dp:2, get:r => r.e.sd2 },
  { key:'r',          label:'r',       dp:2, get:r => r.e.r },
  { key:'rCorrected', label:'corr. r', dp:2, get:r => r.e.rCorrected,
    title:'Attenuation-corrected test–retest correlation' },
  { key:'n',          label:'N',       dp:0, get:r => r.e.n },
  { key:'ci',         label:'CI r',    dp:2, get:r => r.rel.r, emphasis:true,
    title:'The reliability this measure\'s confidence interval is built from, with no patient age entered' },
  { key:'basis',      label:'Basis',   text:true, get:r => r.rel.basis,
    title:'Where that coefficient comes from' }
];

/* Every entry as a flat row. Instrument, category and age band are derived
   from the group key, which is "<Instrument> <Category> · <Age band>". */
function dbFlatRows(){
  const merged = getMergedDB();
  const custom = getCustom();
  const rows = [];
  Object.keys(merged).sort().forEach(family => {
    const isCustom = !!custom[family];
    const instrument = dbInstrumentOf(family, isCustom);
    const parts = family.split(' · ');
    const base = parts[0];
    const band = parts[1] || '—';
    const category = isCustom ? base : (base.slice(instrument.length).trim() || base);
    Object.entries(merged[family]).forEach(([name, e]) => {
      if (!e || typeof e !== 'object') return;
      rows.push({ family, isCustom, instrument, category, band, name, e, rel: dbReliabilityBasis(e, family, name) });
    });
  });
  return rows;
}

function dbBasisClass(basis){
  if (basis.startsWith('internal')) return 'db-chip-internal';
  if (basis.startsWith('stability')) return 'db-chip-stability';
  if (basis.startsWith('retest')) return 'db-chip-retest';
  return 'db-chip-none';
}

function dbFillSelect(el, placeholder, values){
  const keep = el.value;
  el.innerHTML = '<option value="">' + escapeHtml(placeholder) + '</option>' +
    values.map(v => '<option' + (v === keep ? ' selected' : '') + '>' + escapeHtml(v) + '</option>').join('');
  if (el.value !== keep) el.value = '';
}

function renderDbList(){
  const search = document.getElementById('ct-search').value.toLowerCase().trim();
  const all = dbFlatRows();
  const fInst  = document.getElementById('db-f-inst');
  const fCat   = document.getElementById('db-f-cat');
  const fBand  = document.getElementById('db-f-band');
  const fBasis = document.getElementById('db-f-basis');

  const uniq = (fn, from) => [...new Set(from.map(fn))].sort();
  dbFillSelect(fInst, 'All instruments', uniq(r => r.instrument, all));
  /* Category depends on the instrument, so the dropdown only ever offers
     combinations that exist — otherwise picking one silently empties the
     table and reads as a fault. */
  dbFillSelect(fCat, 'All categories', uniq(r => r.category, all.filter(r => !fInst.value || r.instrument === fInst.value)));
  dbFillSelect(fBand, 'All age bands', uniq(r => r.band, all));
  dbFillSelect(fBasis, 'Any basis', uniq(r => r.rel.basis.replace(' · by age', ''), all));

  let rows = all.filter(r => {
    if (fInst.value && r.instrument !== fInst.value) return false;
    if (fCat.value && r.category !== fCat.value) return false;
    if (fBand.value && r.band !== fBand.value) return false;
    if (fBasis.value && r.rel.basis.replace(' · by age', '') !== fBasis.value) return false;
    if (search){
      const hay = (r.name + ' ' + r.instrument + ' ' + r.category + ' ' + r.band).toLowerCase();
      if (!hay.includes(search)) return false;
    }
    return true;
  });

  if (dbSort.key){
    const col = DB_COLUMNS.find(c => c.key === dbSort.key);
    if (col){
      rows = rows.slice().sort((a, b) => {
        const x = col.get(a), y = col.get(b);
        // Blanks sort last whichever way the column is pointed, so "worst
        // first" never opens with a column of dashes.
        if (x == null && y == null) return 0;
        if (x == null) return 1;
        if (y == null) return -1;
        return (col.text ? String(x).localeCompare(String(y)) : x - y) * dbSort.dir;
      });
    }
  }

  const thead = document.getElementById('db-thead');
  thead.innerHTML = '<tr>' + DB_COLUMNS.map(c => {
    const on = dbSort.key === c.key;
    return '<th class="db-th' + (c.text ? ' db-th-text' : '') + '"' +
      (on ? ' aria-sort="' + (dbSort.dir === 1 ? 'ascending' : 'descending') + '"' : '') +
      (c.title ? ' title="' + escapeAttr(c.title) + '"' : '') +
      ' data-db-sort="' + escapeAttr(c.key) + '" tabindex="0" role="button">' + c.label +
      (on ? '<span class="db-arrow">' + (dbSort.dir === 1 ? '▲' : '▼') + '</span>' : '') + '</th>';
  }).join('') + '<th class="db-th"></th></tr>';

  const sortBy = (key) => {
    if (dbSort.key === key) dbSort.dir = -dbSort.dir;
    else { dbSort.key = key; dbSort.dir = 1; }
    renderDbList();
  };
  thead.querySelectorAll('[data-db-sort]').forEach(th => {
    th.addEventListener('click', () => sortBy(th.dataset.dbSort));
    th.addEventListener('keydown', ev => {
      if (ev.key === 'Enter' || ev.key === ' '){ ev.preventDefault(); sortBy(th.dataset.dbSort); }
    });
  });

  const tbody = document.getElementById('db-tbody');
  tbody.innerHTML = rows.map(r => '<tr>' + DB_COLUMNS.map(c => {
    const v = c.get(r);
    if (c.key === 'basis'){
      const byAge = v.includes(' · by age');
      return '<td class="db-td-text"><span class="db-chip ' + dbBasisClass(v) + '">' +
        escapeHtml(v.replace(' · by age', '')) + '</span>' +
        (byAge ? '<span class="db-byage">by age</span>' : '') + '</td>';
    }
    if (c.text){
      const tag = c.key === 'measure' && r.e.metric
        ? '<span class="db-metric-tag">' + escapeHtml(scoreTypeAbbr(r.e.metric) || r.e.metric) + '</span>' : '';
      const cust = c.key === 'instrument' && r.isCustom ? '<span class="custom-tag">Custom</span>' : '';
      return '<td class="db-td-text' + (c.key === 'measure' ? ' db-td-measure' : ' db-td-ctx') + '">' +
        escapeHtml(String(v)) + tag + cust + '</td>';
    }
    return '<td class="db-td-num' + (c.emphasis ? ' db-td-rel' : '') + '">' +
      (v == null ? '<span class="db-dash">–</span>' : escapeHtml(fmt(v, c.dp))) + '</td>';
  }).join('') +
    '<td class="db-td-text">' + (r.isCustom
      ? '<button class="btn btn-ghost btn-icon" data-del-sub="' + escapeAttr(r.family) + '::' + escapeAttr(r.name) + '" title="Remove this custom measure">×</button>'
      : '') + '</td></tr>').join('') ||
    '<tr><td class="db-td-text" colspan="' + (DB_COLUMNS.length + 1) + '"><span class="db-dash">Nothing matches these filters.</span></td></tr>';

  document.getElementById('db-count').innerHTML =
    'Showing <b>' + rows.length + '</b> of <b>' + all.length + '</b> measures';
  const sortCol = DB_COLUMNS.find(c => c.key === dbSort.key);
  document.getElementById('db-sortnote').textContent = sortCol
    ? 'Sorted by ' + sortCol.label + (dbSort.dir === 1 ? ' ascending' : ' descending')
    : 'Click a column heading to sort';

  /* Deleting a custom measure removes its family too once the last one goes,
     which is why there is no separate "delete family" control any more — the
     flat table has no group header to hang one on. */
  tbody.querySelectorAll('[data-del-sub]').forEach(b => b.addEventListener('click', () => {
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
/* `selectedFamily` is vestigial: it used to re-select the family that had just
   been added in the entry form, and is kept only so the many existing
   refreshAll(x) call sites need no edit. Nothing reads it now. */
function refreshAll(selectedFamily){
  renderDbList();
  rebuildAllFamilyLists();
  rebuildBatteryFamilyList();
  rebuildSdiFamilyList();
}
/* Score types a custom entry may declare. '' means "work it out from the mean",
   which is the old behaviour and stays the default so nothing changes for
   entries already saved.

   This column exists because the mean-based guess has no way to recognise a
   raw score. A user-created test with raw norms — say a recognition total out
   of 20, M 19.6, SD 0.8 — was read as a scaled score, and a genuinely
   below-average 19 printed as the 99.9th percentile, "Very Superior". The
   built-in raw families are tagged in normDB; this is the same escape hatch
   for tests the clinician enters. It also lets any other misclassification be
   corrected by hand rather than argued with. */
document.getElementById('ct-search').addEventListener('input', renderDbList);
['db-f-inst', 'db-f-cat', 'db-f-band', 'db-f-basis'].forEach(id => {
  const el = document.getElementById(id);
  if (el) el.addEventListener('change', renderDbList);
});
{
  const clear = document.getElementById('db-f-clear');
  if (clear) clear.addEventListener('click', () => {
    ['db-f-inst', 'db-f-cat', 'db-f-band', 'db-f-basis'].forEach(id => {
      const el = document.getElementById(id); if (el) el.value = '';
    });
    document.getElementById('ct-search').value = '';
    dbSort = { key: null, dir: 1 };
    renderDbList();
  });
}
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
    /* The try covers PARSING ONLY. It used to wrap refreshAll() too, so when
       refreshAll() threw, a successful import saved the data, announced
       success, then announced "Invalid JSON file" from the catch. Anything
       that is not a parse failure must stay outside this block. */
    let data;
    try {
      data = JSON.parse(reader.result);
    } catch (err) {
      showToast('Invalid JSON file', true);
      return;
    }
    /* typeof null is 'object', and so is an array - both used to get past the
       old guard and produce families keyed "0", "1", ... */
    if (!data || typeof data !== 'object' || Array.isArray(data)){
      showToast('That file is not a custom-test database', true);
      return;
    }

    /* Validate every entry against the same rules the on-screen form enforces,
       and import only what passes. A partially-bad file imports its good rows
       and names what it dropped, rather than being rejected wholesale or - as
       before - accepted wholesale. */
    const c = getCustom();
    let accepted = 0;
    const rejected = [];
    Object.entries(data).forEach(([fam, subs]) => {
      if (!fam || !subs || typeof subs !== 'object' || Array.isArray(subs)){
        rejected.push(`${fam || '(unnamed family)'}: not a family of subtests`);
        return;
      }
      const clean = {};
      Object.entries(subs).forEach(([name, raw]) => {
        const v = ctValidateEntry(raw);
        if (v.ok){ clean[name] = v.entry; accepted++; }
        else rejected.push(`${fam} / ${name}: ${v.reason}`);
      });
      if (Object.keys(clean).length) c[fam] = { ...(c[fam] || {}), ...clean };
    });

    if (rejected.length){
      console.warn('Custom test import skipped ' + rejected.length + ' entr' +
                   (rejected.length === 1 ? 'y' : 'ies') + ':\n  ' + rejected.join('\n  '));
    }
    if (!accepted){
      showToast(rejected.length ? `No valid entries - ${rejected[0]}` : 'File contained no test data', true);
      return;
    }
    saveCustom(c);
    refreshAll();
    showToast(rejected.length
      ? `✓ Imported ${accepted}; skipped ${rejected.length} invalid (see console)`
      : `✓ Imported ${accepted} entr${accepted === 1 ? 'y' : 'ies'}`,
      rejected.length > 0);
  };
  reader.readAsText(file);
  e.target.value = '';
});

/* ============================================================
   08 · PREMORBID ESTIMATE
   Ports the user's reference implementation. All formulas/data
   from WAIS-IV/WMS-IV manuals + Crawford & Allan (1997) + OPIE-4.
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

  /* 3. Crawford & Allan (1997) — The Clinical Neuropsychologist, 11(2), 192-197
     (the paper is 1997; a 2001 date circulated here previously was a citation
     error — vol. 11 of that journal is 1997, and the authors' 2001 paper is the
     Crawford, Millar & Milne clinical-judgement comparison in BJCP 40).
     Derived from 200 healthy adults representative of the adult UK population;
     the same lab sample is described in the authors' companion WAIS-R papers as
     ages 16-83. Gate the FLOOR only: #pre-age accepts from 5, and an adult
     WAIS-R equation must not return an estimate for a child. No ceiling is
     enforced — the age term is linear and shallow (0.18/yr), unlike OPIE's
     cubic, and the sample maximum is attested via the companion papers rather
     than printed in the 1997 brief report itself. */
  let v3 = null;
  const caAgeOk = age != null && age >= CRAWFORD_ALLAN_AGE_MIN;
  if (occC != null && edu != null && caAgeOk){
    v3 = 87.14 - 5.21*occC + 1.78*edu + 0.18*age;
  }
  rows.push({ name:'Crawford & Allan (1997) Demographic', val:v3, see:9.11, r:0.73, tipKey:'crawfordAllan' });

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

  /* Say WHY rows are empty when the age is simply outside a model's norms.
     Only when an age is actually present and out of range — a blank age means
     "not entered yet", which is a different thing and needs no explanation.
     Two branches, because a paediatric age blanks BOTH age-gated models while
     an over-90 age blanks OPIE-4 alone. The old single message called every
     non-OPIE model age-free, which was false: the Crawford & Allan (1997)
     equation carries +0.18 × age. Only the two ToPF models on this tab are
     truly age-free. */
  const rangeNote = document.getElementById('pre-age-range-note');
  if (rangeNote){
    const belowAdult = age != null && age < CRAWFORD_ALLAN_AGE_MIN;
    const aboveOpie  = age != null && age > OPIE_AGE_MAX;
    rangeNote.hidden = !(belowAdult || aboveOpie);
    if (belowAdult){
      rangeNote.innerHTML = `<strong>OPIE-4 and the Crawford &amp; Allan model do not apply at age ${escapeHtml(String(age))}.</strong> `
        + `OPIE-4 is fitted for ages ${OPIE_AGE_MIN}–${OPIE_AGE_MAX}, and the Crawford &amp; Allan (1997) equation `
        + `was derived from an adult sample, so neither produces an estimate. `
        + `The two ToPF models carry no age term and are unaffected.`;
    } else if (aboveOpie){
      rangeNote.innerHTML = `<strong>OPIE-4 does not apply at age ${escapeHtml(String(age))}.</strong> `
        + `Its equations are fitted for ages ${OPIE_AGE_MIN}–${OPIE_AGE_MAX}, so no estimate is produced. `
        + `The two ToPF models carry no age term and are unaffected. The Crawford &amp; Allan estimate is still `
        + `produced, but note an age this far above its adult derivation sample extrapolates its linear age term.`;
    }
  }

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
  /* Age is used by the WMS-IV branch below and NOT by the WAIS-IV one. The
     WAIS-IV demographic equations run on ToPF score, education and sex; the
     WMS-IV ToPF-predicted indices carry an explicit age term. */
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
  /* WMS-IV predictions. Gated on the age being IN RANGE, not merely present.
     This used to test `age != null` alone, which was harmless while #pre-age
     was itself bounded 16-90 by the markup — the input could not hold anything
     else. Now that the box carries the shared patient age and accepts from 5
     (see PATIENT_AGE_INPUTS), presence is no longer a proxy for validity: a
     paediatric age would otherwise be fed straight into an adult equation and
     return a plausible-looking index for a child. */
  const topfAgeOk = age != null && age >= TOPF_AGE_MIN && age <= TOPF_AGE_MAX;
  WMS_COEF.forEach(c => {
    let pred = null;
    if (topf != null && topfAgeOk){
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

  /* NOTHING ENTERED YET → NO TABLE. Its two siblings already do this
     (renderPreEstimatesApa needs a valid estimate, renderPrePredictApa needs an
     achieved score); this one did not, and getOpiePredictions ALWAYS returns its
     eight models — `val:null` where the inputs are absent — so the container
     always held a full .apa-table of dashes.

     That is not cosmetic. The Working Report's observer collects any container
     holding an .apa-table, and the top-bar patient age recalculates the premorbid
     page from EVERY page (the #patient-age listener calls calcOpiePredict). So
     typing an age anywhere in the app mutated this container and put an empty
     OPIE-4 table in the report, for a patient with no OPIE data and possibly on a
     page that has nothing to do with premorbid estimation.

     The test is a finite prediction on at least one model, which needs age in
     range, sex, and Vocabulary or Matrix Reasoning — i.e. the clinician has
     actually run OPIE-4. NOT the presence of an achieved score: unlike the ToPF
     tab, this is the only place OPIE-4 predictions are reported, so requiring one
     would keep the premorbid estimate itself out of the report. */
  if (!rows.some(r => r.val != null && Number.isFinite(r.val))){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Enter age, sex and at least one of Vocabulary or Matrix Reasoning to generate OPIE-4 predictions.</div>';
    return;
  }

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
  /* The mirror half of the shared patient age. Separate from the loop above so
     the premorbid recalcs stay exactly as they were; this only adds the mirror
     into #patient-age and the Score Tables re-render. */
  const ageEl = preGet('pre-age');
  if (ageEl){
    const mirror = () => {
      syncPatientAge('pre-age');
      if (typeof renderBattery === 'function') renderBattery();
      if (typeof refreshPatientAgeIndicator === 'function') refreshPatientAgeIndicator();
      if (typeof rebuildFamilyListsForAge === 'function') rebuildFamilyListsForAge();
    };
    ageEl.addEventListener('input', mirror);
    ageEl.addEventListener('change', mirror);
  }
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

/* ================================================================
   PERFORMANCE VALIDITY (PVT) PAGE

   Scores four validity indicators against their published cut-offs
   (constants in data.js, pinned by check.js §38):
     - RBANS Effort Index      Silverberg et al. (2007)
     - RBANS Effort Scale      Novitski et al. (2012) — GATED, see below
     - Reliable Digit Span     Greiffenstein et al. (1994)
     - TOMM                    Martin et al. (2020) meta-analysis

   Design constraints, in order of importance:
   1. Results are CUT-OFF COMPARISONS, never verdicts. "Fail" means the
      score is beyond the published cut-off — the conventional PVT term —
      and the APA note defines it as exactly that. No output may label an
      examinee or a protocol invalid from a single indicator.
   2. The ES gate is mandatory, not advisory: in intact examinees free
      recall normally exceeds ceiling-limited recognition, so an ungated
      ES over-flags. The gate unmet is rendered as "not computed", never
      as a number.
   3. TOMM PPP/NPP are DERIVED by Bayes from the pinned meta-analytic
      sensitivity/specificity and the selected base rate, exactly as
      Martin et al. built Tables 16-17 from the same values. 57 of the 60
      published cells reproduce at 2 dp (§38); deriving rather than
      storing avoids transcribing the one cell that does not.
   4. One APA export for the whole page (pvt-apa, on the Summary tab):
      the four indicators belong in one report table, and four
      near-identical tables would recreate the Change Analysis bloat
      that consent gating exists to fix.

   The two RBANS inputs (Digit Span, List Recognition) appear on both the
   EI and ES tabs; pvtState is the master and every [data-pvt-field] input
   is a synced view of it, so a value entered on either tab carries to the
   other. Raw entry only — this page never reads normDB.
   ================================================================ */

const pvtState = { ds:'', lrec:'', listRecall:'', storyRecall:'', figRecall:'' };

/* Parse one raw-score field. Returns:  undefined = empty,  null = not a
   whole number in [min, max],  number = usable. The two non-values stay
   distinct so an out-of-range entry reads as "check the value" rather
   than as still-waiting-for-input. */
function pvtInt(v, min, max){
  if (v === '' || v === null || v === undefined) return undefined;
  const n = Number(v);
  if (!Number.isFinite(n) || !Number.isInteger(n) || n < min || n > max) return null;
  return n;
}

/* ---------- Effort Index (Silverberg et al., 2007) ---------- */
function pvtEiWeight(bands, raw){
  for (const b of bands) if (raw >= b.min && raw <= b.max) return b.w;
  return null;
}
function getPvtEi(){
  const ds = pvtInt(pvtState.ds, 0, 16);
  const lr = pvtInt(pvtState.lrec, 0, 20);
  if (ds === undefined && lr === undefined) return { empty: true };
  if (ds === null || lr === null) return { invalid: true };
  if (ds === undefined || lr === undefined) return { partial: true };
  const wDs = pvtEiWeight(PVT_EI_WEIGHTS.digitSpan, ds);
  const wLr = pvtEiWeight(PVT_EI_WEIGHTS.listRecognition, lr);
  const cutKey = document.getElementById('pvt-ei-cutoff')?.value === 'sensitive' ? 'sensitive' : 'standard';
  const cut = PVT_EI_CUTOFFS[cutKey];
  const ei = wDs + wLr;
  return { ds, lr, wDs, wLr, ei, cutKey, cut, fail: ei > cut };
}

/* ---------- Effort Scale (Novitski et al., 2012) ---------- */
function getPvtEs(){
  const ds   = pvtInt(pvtState.ds, 0, 16);
  const lr   = pvtInt(pvtState.lrec, 0, 20);
  const list = pvtInt(pvtState.listRecall, 0, 10);
  const stor = pvtInt(pvtState.storyRecall, 0, 12);
  const fig  = pvtInt(pvtState.figRecall, 0, 20);
  const vals = [ds, lr, list, stor, fig];
  if (vals.every(v => v === undefined)) return { empty: true };
  if (vals.some(v => v === null)) return { invalid: true };
  if (vals.some(v => v === undefined)) return { partial: true };
  const g = PVT_ES.gate;
  const gateMet = ds < g.digitSpanBelow || lr < g.listRecognitionBelow || (ds + lr) < g.combinedBelow;
  if (!gateMet) return { ds, lr, gated: true };
  const es = (lr - (list + stor + fig)) + ds;
  return { ds, lr, list, stor, fig, es, fail: es < PVT_ES.cutoff };
}

/* ---------- Reliable Digit Span (Greiffenstein et al., 1994) ---------- */
function getPvtRds(){
  const f = pvtInt(document.getElementById('pvt-rds-f')?.value, 0, 9);
  const b = pvtInt(document.getElementById('pvt-rds-b')?.value, 0, 8);
  if (f === undefined && b === undefined) return { empty: true };
  if (f === null || b === null) return { invalid: true };
  if (f === undefined || b === undefined) return { partial: true };
  const sum = f + b;
  /* Floor rule: failing at least one trial each of the lowest items is
     recorded as RDS = 3 — with two-number entry that is any sum below 3. */
  const floored = sum < PVT_RDS.floor;
  const rds = floored ? PVT_RDS.floor : sum;
  const conservative = document.getElementById('pvt-rds-cutoff')?.value !== 'traditional';
  const cut = conservative ? PVT_RDS.cutoffConservative : PVT_RDS.cutoffTraditional;
  return { f, b, rds, floored, conservative, cut, fail: rds <= cut };
}

/* ---------- WAIS Digit Span ACSS (Iverson & Tulsky, 2003; Axelrod et al., 2006) ---------- */
function getPvtDs(){
  const acss  = pvtInt(document.getElementById('pvt-ds-acss')?.value, 1, 19);
  const vocab = pvtInt(document.getElementById('pvt-ds-vocab')?.value, 1, 19);
  if (acss === undefined && vocab === undefined) return { empty: true };
  if (acss === null || vocab === null) return { invalid: true };
  if (acss === undefined) return { partial: true };   // vocab alone computes nothing
  const conservative = document.getElementById('pvt-ds-cutoff')?.value !== 'sensitive';
  const cut = conservative ? PVT_DS.cutoffConservative : PVT_DS.cutoffSensitive;
  const s = { acss, conservative, cut, fail: acss <= cut };
  if (vocab !== undefined){
    s.vocab = vocab;
    s.diff = vocab - acss;
    s.diffFail = s.diff >= PVT_DS.vocabDiffCutoff;
  }
  const lsf = pvtInt(document.getElementById('pvt-ds-lsf')?.value, 0, 9);
  const lsb = pvtInt(document.getElementById('pvt-ds-lsb')?.value, 0, 8);
  if (lsf === null || lsb === null) return { invalid: true };
  if (lsf !== undefined){
    s.lsf = lsf;
    /* Iverson & Tulsky limit the forward-span index to under-55s (its base
       rate climbs to 11% in the oldest band). Without an age it is
       WITHHELD, never guessed — the same posture as the base-rate rows on
       Score Tables. */
    const age = (typeof patientAge === 'function') ? patientAge() : null;
    if (age === null) s.lsfState = 'no-age';
    else if (age >= PVT_DS.lsfAgeBelow) s.lsfState = 'over-age';
    else { s.lsfState = 'ok'; s.lsfFail = lsf <= PVT_DS.lsfCutoff; }
  }
  if (lsb !== undefined){
    s.lsb = lsb;
    s.lsbFail = lsb <= PVT_DS.lsbCutoff;
  }
  return s;
}

/* ---------- Rey 15-Item + recognition (Boone et al., 2002) ---------- */
function getPvtRey(){
  const recall = pvtInt(document.getElementById('pvt-rey-recall')?.value, 0, 15);
  const recog  = pvtInt(document.getElementById('pvt-rey-recog')?.value, 0, 15);
  const fp     = pvtInt(document.getElementById('pvt-rey-fp')?.value, 0, 15);
  if (recall === undefined && recog === undefined && fp === undefined) return { empty: true };
  if (recall === null || recog === null || fp === null) return { invalid: true };
  if (recall === undefined) return { partial: true };   // recognition alone computes nothing
  const s = { recall, recallFail: recall < PVT_REY15.recallCutoff };
  /* The combination needs BOTH recognition numbers — a recognition trial
     with an unrecorded false-positive count is half a score. */
  if (recog !== undefined && fp !== undefined){
    s.recog = recog; s.fp = fp;
    s.combo = recall + (recog - fp);
    s.comboFail = s.combo < PVT_REY15.comboCutoff;
  } else if (recog !== undefined || fp !== undefined){
    s.comboPartial = true;
  }
  return s;
}

/* ---------- Rey 15-Item administration ----------
   Displays the stimulus for the protocol's 10 seconds, then the
   recognition page, and writes the two recognition scores back into the
   tab so they cannot be mis-transcribed. Free recall stays a manual
   entry: the examinee draws it on paper.

   MOUNTED ON <body>, NEVER INSIDE THE SECTION. position:fixed resolves
   against the nearest TRANSFORMED ancestor, and staggerSectionContent
   animates panel children on every page entry — an overlay built inside
   the panel would be positioned against a moving box. Same reason the
   consent card lives on body.

   TEST SECURITY: nothing here renders until the clinician clicks through
   an explicit confirmation in-session. There is no deep link to it, it
   is not in the page's static markup, and it never reaches an export. */
const REY15_SHAPE_SVG = {
  circle:        '<circle cx="26" cy="26" r="19"/>',
  square:        '<rect x="7" y="7" width="38" height="38"/>',
  triangle:      '<path d="M26 6 L46 45 L6 45 Z"/>',
  diamond:       '<path d="M26 5 L45 26 L26 47 L7 26 Z"/>',
  pentagon:      '<path d="M26 6 L45 20 L38 43 L14 43 L7 20 Z"/>',
  parallelogram: '<path d="M17 12 L48 12 L35 40 L4 40 Z"/>',
  rule1:         '<line x1="6" y1="26" x2="46" y2="26"/>',
  rule2:         '<line x1="6" y1="19" x2="46" y2="19"/><line x1="6" y1="33" x2="46" y2="33"/>',
  rule3:         '<line x1="6" y1="15" x2="46" y2="15"/><line x1="6" y1="26" x2="46" y2="26"/><line x1="6" y1="37" x2="46" y2="37"/>'
};
function rey15ItemHtml(id){
  if (REY15_SHAPES.indexOf(id) !== -1){
    return `<svg class="rey-glyph" viewBox="0 0 52 52" fill="none" stroke="currentColor" stroke-width="3" stroke-linejoin="round" aria-hidden="true">${REY15_SHAPE_SVG[id]}</svg>`;
  }
  return `<span class="rey-glyph rey-glyph-text">${escapeHtml(id)}</span>`;
}
function rey15IsTarget(id){
  return REY15_RECALL_ROWS.some(row => row.indexOf(id) !== -1);
}

const reyAdmin = { el: null, step: null, selected: null, timer: null };

function reyAdminClose(){
  if (reyAdmin.timer){ clearInterval(reyAdmin.timer); reyAdmin.timer = null; }
  if (reyAdmin.el){ reyAdmin.el.remove(); reyAdmin.el = null; }
  reyAdmin.step = null;
  document.removeEventListener('keydown', reyAdminKey);
  document.getElementById('pvt-rey-administer')?.focus();
}
function reyAdminKey(e){ if (e.key === 'Escape') reyAdminClose(); }

function reyAdminOpen(){
  reyAdminClose();
  reyAdmin.selected = new Set();
  const el = document.createElement('div');
  el.className = 'rey-admin';
  el.setAttribute('role', 'dialog');
  el.setAttribute('aria-modal', 'true');
  el.setAttribute('aria-label', 'Rey 15-Item administration');
  document.body.appendChild(el);
  reyAdmin.el = el;
  el.addEventListener('click', reyAdminClick);
  document.addEventListener('keydown', reyAdminKey);
  reyAdminGo('intro');
}
function reyAdminClick(e){
  const btn = e.target.closest('[data-rey-action]');
  if (btn){ reyAdminAction(btn.dataset.reyAction); return; }
  const item = e.target.closest('[data-rey-item]');
  if (item && reyAdmin.step === 'recognition'){
    const id = item.dataset.reyItem;
    if (reyAdmin.selected.has(id)) reyAdmin.selected.delete(id);
    else reyAdmin.selected.add(id);
    item.classList.toggle('is-picked', reyAdmin.selected.has(id));
    item.setAttribute('aria-pressed', reyAdmin.selected.has(id) ? 'true' : 'false');
    const n = reyAdmin.el.querySelector('[data-rey-picked]');
    if (n) n.textContent = String(reyAdmin.selected.size);
  }
}
function reyAdminAction(action){
  if (action === 'close'){ reyAdminClose(); return; }
  if (action === 'apply'){ reyAdminApply(); return; }
  reyAdminGo(action);
}

/* The 10-second exposure is the protocol, not an animation: reduced
   motion suppresses the countdown's transition, never its duration. */
function reyAdminStartExposure(){
  let left = 10;
  const tick = () => {
    const n = reyAdmin.el?.querySelector('[data-rey-count]');
    if (n) n.textContent = String(left);
    if (left <= 0){
      clearInterval(reyAdmin.timer);
      reyAdmin.timer = null;
      reyAdminGo('draw');
      return;
    }
    left--;
  };
  tick();
  reyAdmin.timer = setInterval(tick, 1000);
}

function reyAdminGo(step){
  if (!reyAdmin.el) return;
  if (reyAdmin.timer){ clearInterval(reyAdmin.timer); reyAdmin.timer = null; }
  reyAdmin.step = step;
  const close = '<button type="button" class="rey-close" data-rey-action="close" aria-label="Close administration">×</button>';
  if (step === 'intro'){
    reyAdmin.el.innerHTML = `${close}<div class="rey-panel">
      <div class="rey-kicker">Rey 15-Item · administration</div>
      <h2 class="rey-title">Ready to begin?</h2>
      <p class="rey-copy">This displays the test stimuli full screen. Start it only with the examinee present, and have paper and a pen ready for their drawing.</p>
      <p class="rey-copy rey-script">Read: “I'm going to show you a page with 15 different things on it for just a short period of time, and I want you to learn as many of the things as you can. When I take away the page, I'll want you to draw as many of the things as you can remember. Keep in mind, there are <strong>15</strong> different things so you will have to learn them very quickly.”</p>
      <div class="rey-actions">
        <button type="button" class="btn" data-rey-action="close">Cancel</button>
        <button type="button" class="btn btn-primary" data-rey-action="exposure">Show the stimulus (10 s)</button>
      </div>
    </div>`;
    return;
  }
  if (step === 'exposure'){
    reyAdmin.el.innerHTML = `<div class="rey-stage">
      <div class="rey-grid rey-grid-3">${REY15_RECALL_ROWS.flat().map(id =>
        `<div class="rey-cell">${rey15ItemHtml(id)}</div>`).join('')}</div>
      <div class="rey-countdown"><span data-rey-count>10</span></div>
    </div>`;
    reyAdminStartExposure();
    return;
  }
  if (step === 'draw'){
    reyAdmin.el.innerHTML = `${close}<div class="rey-panel">
      <div class="rey-kicker">Rey 15-Item · recall</div>
      <h2 class="rey-title">Now draw what you can remember</h2>
      <p class="rey-copy">The examinee draws on paper. Score the free recall yourself and enter it on the tab. This step is not scored on screen.</p>
      <p class="rey-copy rey-script">Then read: “On this page are the 15 things I showed you as well as 15 items that were not on the page. I want you to circle the things you remember from the page I showed you.”</p>
      <div class="rey-actions">
        <button type="button" class="btn" data-rey-action="close">Stop here</button>
        <button type="button" class="btn btn-primary" data-rey-action="recognition">Show the recognition page</button>
      </div>
    </div>`;
    return;
  }
  if (step === 'recognition'){
    reyAdmin.el.innerHTML = `<div class="rey-stage">
      <div class="rey-grid rey-grid-6">${REY15_RECOGNITION_ROWS.flat().map(id =>
        `<button type="button" class="rey-cell rey-pick" data-rey-item="${escapeAttr(id)}" aria-pressed="false">${rey15ItemHtml(id)}</button>`).join('')}</div>
      <div class="rey-stage-foot">
        <span class="rey-picked"><span data-rey-picked>0</span> selected</span>
        <button type="button" class="btn btn-primary" data-rey-action="done">Finish recognition</button>
      </div>
    </div>`;
    return;
  }
  if (step === 'done'){
    const picked = [...reyAdmin.selected];
    const correct = picked.filter(rey15IsTarget).length;
    const fp = picked.length - correct;
    reyAdmin.el.innerHTML = `${close}<div class="rey-panel">
      <div class="rey-kicker">Rey 15-Item · recognition scored</div>
      <h2 class="rey-title">${correct} correct · ${fp} false positive${fp === 1 ? '' : 's'}</h2>
      <p class="rey-copy">Scored from the ${picked.length} item${picked.length === 1 ? '' : 's'} selected. Applying these fills the recognition fields on the tab; enter the free-recall score from the drawing yourself.</p>
      <div class="rey-actions">
        <button type="button" class="btn" data-rey-action="close">Discard</button>
        <button type="button" class="btn btn-primary" data-rey-action="apply">Apply to the tab</button>
      </div>
    </div>`;
  }
}

function reyAdminApply(){
  const picked = [...reyAdmin.selected];
  const correct = picked.filter(rey15IsTarget).length;
  const fp = picked.length - correct;
  const set = (id, v) => {
    const el = document.getElementById(id);
    if (!el) return;
    el.value = String(v);
    el.dispatchEvent(new Event('input', { bubbles: true }));
  };
  set('pvt-rey-recog', correct);
  set('pvt-rey-fp', fp);
  reyAdminClose();
  if (typeof showToast === 'function') showToast('Recognition scores applied');
}

/* ---------- TOMM (Martin et al., 2020) ---------- */
function pvtPPP(sens, spec, br){ return (sens * br) / (sens * br + (1 - spec) * (1 - br)); }
function pvtNPP(sens, spec, br){ return (spec * (1 - br)) / (spec * (1 - br) + (1 - sens) * br); }
function pvtTommCutoffById(id){ return PVT_TOMM_CUTOFFS.find(c => c.id === id) || null; }
function getPvtTomm(){
  const t1  = pvtInt(document.getElementById('pvt-tomm-t1')?.value, 0, 50);
  const t2  = pvtInt(document.getElementById('pvt-tomm-t2')?.value, 0, 50);
  const ret = pvtInt(document.getElementById('pvt-tomm-ret')?.value, 0, 50);
  if (t1 === undefined && t2 === undefined && ret === undefined) return { empty: true };
  if (t1 === null || t2 === null || ret === null) return { invalid: true };
  const t1CutId = document.getElementById('pvt-tomm-t1cut')?.value === 't1-41' ? 't1-41' : 't1-42';
  const lateCut = document.getElementById('pvt-tomm-t2cut')?.value === '49' ? '49' : '45';
  const br = parseFloat(document.getElementById('pvt-tomm-br')?.value) || 0.10;
  const rows = [];
  const addRow = (label, score, cutoff) => {
    if (score === undefined || !cutoff) return;
    rows.push({
      label, score, cutoff,
      fail: score < cutoff.cut,
      ppp: pvtPPP(cutoff.sens, cutoff.spec, br),
      npp: pvtNPP(cutoff.sens, cutoff.spec, br)
    });
  };
  addRow('Trial 1',   t1,  pvtTommCutoffById(t1CutId));
  addRow('Trial 2',   t2,  pvtTommCutoffById('t2-' + lateCut));
  addRow('Retention', ret, pvtTommCutoffById('ret-' + lateCut));
  return { rows, br, anyFail: rows.some(r => r.fail) };
}

/* ---------- result boxes ---------- */
/* ---- CVLT-3 Forced Choice Recognition (Delis et al., 2017, Appendix D) ----
   The one measure here scored from a published base-rate table rather than
   a published cut-off, because the manual prints no cut-off. Everything the
   flag rests on is derived from Tables D.13-D.15 at runtime. */

/* The manual's age band containing this age, or 'All ages' when there is no
   age, or the age falls outside the 16-90 normative range. Both readings are
   published, so both are citable; the caller says which one it used. */
function pvtCvlt3Band(age){
  if (age === null || age === undefined || !isFinite(age)) return 'All ages';
  const b = PVT_CVLT3_BANDS.find(x => age >= x.lo && age <= x.hi);
  return b ? b.key : 'All ages';
}
/* Base rate for a score in a band. Hits read DOWNWARD (percentage scoring
   that many or fewer) and clamp onto the manual's "≤12" row; critical items
   read UPWARD (that many or more) and clamp onto its "≥3" row. */
function pvtCvlt3HitsRate(hits, band){
  const t = PVT_CVLT3_FC_HITS.bands[band];
  if (!t) return null;
  const key = Math.min(Math.max(hits, PVT_CVLT3_FC_HITS.floorRow), PVT_CVLT3_FC_HITS.max);
  const v = t[key];
  return typeof v === 'number' ? v : null;
}
function pvtCvlt3CriticalRate(count, band, which){
  const spec = PVT_CVLT3_CRITICAL[which];
  const t = spec && spec.bands[band];
  if (!t) return null;
  const key = Math.min(Math.max(count, 0), spec.capRow);
  const v = t[key];
  return typeof v === 'number' ? v : null;
}
/* THE DERIVED THRESHOLD. Not a stored cut-off: the rarest score in THIS
   band whose published base rate is still at or below the selected per-test
   false-positive criterion. On hits that is the highest such score, on
   critical items the lowest. It moves with age by construction — 15 hits is
   8.6% overall but 26% at 80-90, so the same score flags in a 30-year-old
   and does not in an 85-year-old, which is the whole reason the band is
   read rather than ignored. Returns null where no score in the table meets
   the criterion, and the measure then reports its base rate without a flag. */
function pvtCvlt3HitsThreshold(band, pct){
  const t = PVT_CVLT3_FC_HITS.bands[band];
  if (!t) return null;
  let best = null;
  for (let h = PVT_CVLT3_FC_HITS.floorRow; h < PVT_CVLT3_FC_HITS.max; h++){
    const r = t[h];
    if (typeof r === 'number' && r <= pct) best = Math.max(best === null ? h : best, h);
  }
  return best;
}
function pvtCvlt3CriticalThreshold(band, pct, which){
  const spec = PVT_CVLT3_CRITICAL[which];
  const t = spec && spec.bands[band];
  if (!t) return null;
  for (let n = 1; n <= spec.capRow; n++){
    const r = t[n];
    if (typeof r === 'number' && r <= pct) return n;
  }
  return null;
}
function pvtCvlt3Basis(){
  const key = document.getElementById('pvt-cvlt3-basis')?.value;
  return PVT_CVLT3_FC_CUTOFFS.find(c => c.key === key) || PVT_CVLT3_FC_CUTOFFS[0];
}
function pvtCvlt3Criterion(){
  const key = document.getElementById('pvt-cvlt3-criterion')?.value;
  return PVT_CVLT3_CRITERIA.find(c => c.key === key) || PVT_CVLT3_CRITERIA[0];
}
function getPvtCvlt3(){
  const hits = pvtInt(document.getElementById('pvt-cvlt3-hits')?.value, 0, PVT_CVLT3_FC_HITS.max);
  const rc   = pvtInt(document.getElementById('pvt-cvlt3-crit-recall')?.value, 0, 16);
  const yc   = pvtInt(document.getElementById('pvt-cvlt3-crit-recog')?.value, 0, 16);
  if (hits === undefined && rc === undefined && yc === undefined) return { empty: true };
  if (hits === null || rc === null || yc === null) return { invalid: true };
  if (hits === undefined) return { partial: true };   // critical items alone score nothing
  const age = (typeof patientAge === 'function') ? patientAge() : null;
  const band = pvtCvlt3Band(age);
  const crit = pvtCvlt3Criterion();
  const s = { hits, age, band, usedAllAges: band === 'All ages', criterion: crit };
  s.hitsRate = pvtCvlt3HitsRate(hits, band);
  /* The base rate is ALWAYS shown — it is the CVLT-3's own published figure
     and the reason to read the age band at all. What the basis selects is
     what the FLAG rests on: the derived threshold, or one of the two
     CVLT-II cut-offs, which carry a published accuracy pair the manual
     does not. Either way the app prints which, and on which instrument. */
  const basis = pvtCvlt3Basis();
  s.basis = basis;
  s.derivedThreshold = pvtCvlt3HitsThreshold(band, crit.pct);
  s.hitsThreshold = basis.cut === null ? s.derivedThreshold : basis.cut;
  s.hitsFail = s.hitsThreshold !== null && hits <= s.hitsThreshold;
  [['recall', rc, 'critRecall'], ['recognition', yc, 'critRecog']].forEach(([which, val, key]) => {
    if (val === undefined) return;
    const th = pvtCvlt3CriticalThreshold(band, crit.pct, which);
    s[key] = {
      which, count: val, label: PVT_CVLT3_CRITICAL[which].label,
      rate: pvtCvlt3CriticalRate(val, band, which),
      threshold: th,
      fail: th !== null && val >= th
    };
  });
  s.anyFail = !!(s.hitsFail || s.critRecall?.fail || s.critRecog?.fail);
  return s;
}

function pvtResultHtml(kind, headline, detail){
  return `<div class="pvt-result-inner is-${kind}"><span class="pvt-result-headline">${headline}</span>${detail ? `<span class="pvt-result-detail">${detail}</span>` : ''}</div>`;
}

/* ---------- the per-index readout ----------
   ONE ROW PER INDEX, each carrying its OWN state.

   This replaces a single card-level state that coloured every line at once:
   a Digit Span card holding four indices printed "ACSS = 8 — Pass" in
   failure red because one span had flagged. The rule that fixes it is that
   VALUES AND LABELS STAY IN INK and a coloured chip beside them carries the
   state, so a card can hold four indices in four different states and each
   reads correctly. Status colour never appears alone — every chip carries
   its word ("Pass", "Flag", "Not evaluated"), so the state survives
   greyscale, colour-blindness and print.

   Each row also plots the score against its cut-off on a slim track. The
   domain is never invented: it is the measure's own possible range (a
   scaled-score metric, a trial out of 50, or the arithmetic minimum and
   maximum of the raw inputs that feed a composite). The shaded zone is the
   side of the boundary that fails, placed on the half-integer between the
   last passing and first failing score, so it is exact for these
   integer-valued scales rather than approximate. */
const PVT_STATE_WORD = { pass: 'Pass', fail: 'Fail', flag: 'Flag', na: 'n/a' };

/* The state glyph. A raised flag for anything beyond its cut-off, a tick
   for a clear index, a dash where nothing was evaluated. The glyph rides
   INSIDE the chip beside its word, so a state is carried three ways —
   shape, word and colour — and reading it never depends on colour alone.

   There is deliberately no score-against-cut-off plot here: what a
   clinician acts on is the cut-off and whether the score crossed it, and
   a track showing HOW FAR it fell asks them to read a distance that
   carries no published meaning. The cut-off itself is printed on every
   row. */
const PVT_STATE_ICON = {
  flag: '<svg class="pvt-flag-icon" viewBox="0 0 12 12" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><line x1="3" y1="1.5" x2="3" y2="10.5"/><path d="M3 2.2 L9.5 2.2 L8 4.6 L9.5 7 L3 7 Z" fill="currentColor" stroke="none"/></svg>',
  pass: '<svg class="pvt-flag-icon" viewBox="0 0 12 12" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><polyline points="2.5,6.4 4.9,8.8 9.5,3.2"/></svg>',
  na:   '<svg class="pvt-flag-icon" viewBox="0 0 12 12" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"><line x1="2.8" y1="6" x2="9.2" y2="6"/></svg>'
};

/* state: 'pass' | 'fail' | 'flag' | 'na'. `word` overrides the default
   chip text (a base-rate index reads "Flag", a withheld one "Not
   evaluated"). */
function pvtIndexRowHtml(o){
  const word = o.word || PVT_STATE_WORD[o.state] || '';
  const chipState = o.state === 'flag' ? 'fail' : o.state;
  const icon = PVT_STATE_ICON[o.state === 'fail' ? 'flag' : o.state] || '';
  const valueHtml = (o.value === undefined || o.value === null || o.value === '')
    ? '<span class="pvt-index-value is-none">—</span>'
    : `<span class="pvt-index-value">${escapeHtml(String(o.value))}</span>`;
  return `<div class="pvt-index">
    <div class="pvt-index-main">
      <div class="pvt-index-top">
        <span class="pvt-index-label">${o.label}</span>
        ${valueHtml}
      </div>
      ${o.meta ? `<div class="pvt-index-meta">${o.meta}</div>` : ''}
    </div>
    <span class="pvt-flag is-${chipState}">${icon}${word}</span>
  </div>`;
}

/* Wrap rows in the card. `note` carries any arithmetic or explanation that
   belongs to the whole measure rather than to one index. */
function pvtReadoutHtml(rows, note){
  const flagged = rows.filter(r => r.state === 'fail' || r.state === 'flag').length;
  const scored  = rows.filter(r => r.state !== 'na').length;
  const head = (rows.length > 1 && scored > 0)
    ? `<div class="pvt-readout-head"><strong>${flagged}</strong> of ${scored} scored ${scored === 1 ? 'index' : 'indices'} beyond ${flagged === 1 ? 'its' : 'their'} cut-off</div>`
    : '';
  return `<div class="pvt-readout">
    ${head}
    ${rows.map(pvtIndexRowHtml).join('')}
    ${note ? `<div class="pvt-readout-note">${note}</div>` : ''}
  </div>`;
}
const PVT_PROMPTS = {
  partial: 'Enter the remaining score(s) to compute this index.',
  invalid: 'Check the entered values — one is outside the possible raw-score range.'
};
function pvtStatusWord(fail){ return fail ? 'Fail' : 'Pass'; }

function renderPvtEi(){
  const out = document.getElementById('pvt-ei-result');
  if (!out) return;
  const s = getPvtEi();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter both raw scores to compute the Effort Index.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.partial); return; }
  const acc = PVT_EI_ACCURACY[s.cutKey];
  out.innerHTML = pvtReadoutHtml([{
    label: 'Effort Index', value: s.ei, state: s.fail ? 'fail' : 'pass',
    meta: `Cut-off &gt; ${s.cut}${s.cutKey === 'sensitive' ? ' (screening)' : ''} · sens. ${acc.sens} · spec. ${acc.spec} · higher = less credible`
  }], `Digit Span ${s.ds} → weight ${s.wDs}; List Recognition ${s.lr} → weight ${s.wLr}; sum ${s.ei}.${
    s.fail ? ' Corroborate with an independent, preferably forced-choice, measure.' : ''}`);
}

function renderPvtEs(){
  const out = document.getElementById('pvt-es-result');
  if (!out) return;
  const s = getPvtEs();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter all five raw scores to evaluate the Effort Scale.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.partial); return; }
  if (s.gated){
    out.innerHTML = pvtReadoutHtml([{
      label: 'Effort Scale', value: null, state: 'na', word: 'Not computed',
      meta: `Gate not met — Digit Span ${s.ds} (&lt; 9), List Recognition ${s.lr} (&lt; 19), sum ${s.ds + s.lr} (&lt; 28)`
    }], 'On a profile this strong the ES over-flags, so no score is reported (Novitski et al., 2012).');
    return;
  }
  /* Domain is the arithmetic range of the inputs that feed the formula:
     List Recognition 0-20 minus total recall 0-42, plus Digit Span 0-16. */
  out.innerHTML = pvtReadoutHtml([{
    label: 'Effort Scale', value: s.es, state: s.fail ? 'fail' : 'pass',
    meta: `Cut-off &lt; ${PVT_ES.cutoff} · ROC AUC ${PVT_ES_ACCURACY.auc} (no published sens/spec) · lower = less credible`
  }], `(${s.lr} − [${s.list} + ${s.stor} + ${s.fig}]) + ${s.ds} = ${s.es}. Gate met.${
    s.fail ? ' Confirm with a stand-alone forced-choice measure.' : ''}`);
}

function renderPvtRds(){
  const out = document.getElementById('pvt-rds-result');
  if (!out) return;
  const s = getPvtRds();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter both span lengths to compute Reliable Digit Span.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.partial); return; }
  const racc = PVT_RDS_ACCURACY[s.conservative ? 'conservative' : 'traditional'];
  out.innerHTML = pvtReadoutHtml([{
    label: 'Reliable Digit Span', value: s.rds, state: s.fail ? 'fail' : 'pass',
    meta: `Cut-off ≤ ${s.cut} (${s.conservative ? 'conservative' : 'traditional'}) · sens. ${racc.sens} · spec. ${racc.spec} · lower = less credible`
  }], `${s.f} forward + ${s.b} backward${s.floored ? ` = ${s.f + s.b}; the floor rule records any RDS below 3 as 3` : ` = ${s.rds}`}.${
    s.fail ? ' Weigh genuine attentional impairment, anxiety and aphasia before concluding.' : ''}`);
}

function renderPvtDs(){
  const out = document.getElementById('pvt-ds-result');
  if (!out) return;
  const s = getPvtDs();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter the Digit Span scaled score to evaluate this index.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', 'Enter the Digit Span scaled score. Vocabulary alone computes nothing.'); return; }
  const dacc = PVT_DS_ACCURACY[s.conservative ? 'conservative' : 'sensitive'];
  const rows = [{
    label: 'Digit Span ACSS', value: s.acss, state: s.fail ? 'fail' : 'pass',
    meta: `Cut-off ≤ ${s.cut}${s.conservative ? ' (conservative)' : ''} · sens. ${dacc.sens} · spec. ${dacc.spec}`
  }];
  if (s.diff !== undefined) rows.push({
    label: 'Vocabulary − Digit Span', value: s.diff,
    state: s.diffFail ? 'flag' : 'pass', word: s.diffFail ? 'Flag' : 'Pass',
    meta: `Cut-off ≥ ${PVT_DS.vocabDiffCutoff} · base rate ${PVT_DS_VOCABDIFF_BASERATES.standardization} standardisation / ${PVT_DS_VOCABDIFF_BASERATES.clinical} clinical`
  });
  if (s.lsf !== undefined){
    if (s.lsfState === 'ok') rows.push({
      label: 'Longest span forward', value: s.lsf,
      state: s.lsfFail ? 'flag' : 'pass', word: s.lsfFail ? 'Flag' : 'Pass',
      meta: `Cut-off ≤ ${PVT_DS.lsfCutoff} (age &lt; ${PVT_DS.lsfAgeBelow}) · base rate ${PVT_DS_SPAN_BASERATES.lsf.standardization} standardisation / ${PVT_DS_SPAN_BASERATES.lsf.clinical} clinical`
    });
    else rows.push({
      label: 'Longest span forward', value: s.lsf, state: 'na',
      word: s.lsfState === 'no-age' ? 'Needs age' : 'n/a',
      meta: s.lsfState === 'no-age'
        ? `Enter the patient age in the top bar. This index applies under age ${PVT_DS.lsfAgeBelow}.`
        : `Not applicable at age ${PVT_DS.lsfAgeBelow}+. Iverson &amp; Tulsky limit this index to under-55s.`
    });
  }
  if (s.lsb !== undefined) rows.push({
    label: 'Longest span backward', value: s.lsb,
    state: s.lsbFail ? 'flag' : 'pass', word: s.lsbFail ? 'Flag' : 'Pass',
    meta: `Cut-off ≤ ${PVT_DS.lsbCutoff} · base rate ${PVT_DS_SPAN_BASERATES.lsb.standardization} standardisation / ${PVT_DS_SPAN_BASERATES.lsb.clinical} clinical`
  });
  const anyFail = s.fail || s.diffFail || s.lsfFail || s.lsbFail;
  out.innerHTML = pvtReadoutHtml(rows, anyFail
    ? 'Base-rate indices state how often a score this extreme occurs in the standardisation and clinical samples; corroborate with an independent, preferably forced-choice, measure.'
    : '');
}

function renderPvtRey(){
  const out = document.getElementById('pvt-rey15-result');
  if (!out) return;
  const s = getPvtRey();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter the free-recall score to evaluate this test.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', 'Enter the free-recall score. The recognition trial alone computes nothing.'); return; }
  const rows = [{
    label: 'Free recall', value: s.recall, state: s.recallFail ? 'fail' : 'pass',
    meta: `Cut-off &lt; ${PVT_REY15.recallCutoff} · sens. ${PVT_REY15_ACCURACY.recall.sens} · spec. ${PVT_REY15_ACCURACY.recall.spec}`
  }];
  let note = '';
  if (s.combo !== undefined){
    rows.push({
      label: 'Combination score', value: s.combo, state: s.comboFail ? 'fail' : 'pass',
      meta: `Cut-off &lt; ${PVT_REY15.comboCutoff} · sens. ${PVT_REY15_ACCURACY.combo.sens} · spec. ${PVT_REY15_ACCURACY.combo.spec} · the more sensitive index`
    });
    note = `${s.recall} + (${s.recog} − ${s.fp}) = ${s.combo}. The combination raises sensitivity from ${PVT_REY15_ACCURACY.recall.sens} to ${PVT_REY15_ACCURACY.combo.sens} at comparable specificity (Boone et al., 2002).`;
  } else if (s.comboPartial){
    rows.push({
      label: 'Combination score', value: null, state: 'na', word: 'Incomplete',
      meta: 'Enter both recognition correct and false positives to compute this index'
    });
  }
  if ((s.recallFail || s.comboFail) && !note) note = 'Corroborate with an independent, preferably forced-choice, measure.';
  out.innerHTML = pvtReadoutHtml(rows, note);
}

function renderPvtTomm(){
  const out = document.getElementById('pvt-tomm-result');
  const power = document.getElementById('pvt-tomm-power');
  if (!out || !power) return;
  const s = getPvtTomm();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter at least one trial score to evaluate the TOMM.', 'The result appears here and the outcome joins the running summary.'); power.innerHTML = ''; return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); power.innerHTML = ''; return; }
  out.innerHTML = pvtReadoutHtml(s.rows.map(r => ({
    label: `TOMM ${r.label}`, value: r.score, state: r.fail ? 'fail' : 'pass',
    meta: `Cut-off &lt; ${r.cutoff.cut} · sens. ${r.cutoff.sensRange || r.cutoff.sens.toFixed(2).replace(/^0/, '')} · spec. ${r.cutoff.specRange || r.cutoff.spec.toFixed(2).replace(/^0/, '')} · PPP ${r.ppp.toFixed(2).replace(/^0/, '')} at this base rate`
  })), s.anyFail
    ? 'At low base rates a single failure has modest positive predictive power. See the table below.'
    : '');
  const brPct = Math.round(s.br * 100);
  power.innerHTML = `
    <div class="pvt-card">
      <div class="pvt-card-kicker">Predictive power at a ${brPct}% base rate of invalidity</div>
      <table class="pvt-table">
        <thead><tr><th>Trial · cut-off</th><th>Sens.</th><th>Spec.</th><th>PPP</th><th>NPP</th></tr></thead>
        <tbody>${s.rows.map(r => `
          <tr><td>${r.cutoff.label.replace('<', '&lt;')}</td>
          <td>${r.cutoff.sensRange || r.cutoff.sens.toFixed(2).replace(/^0/, '')}</td>
          <td>${r.cutoff.specRange || r.cutoff.spec.toFixed(2).replace(/^0/, '')}</td>
          <td>${r.ppp.toFixed(2).replace(/^0/, '')}</td>
          <td>${r.npp.toFixed(2).replace(/^0/, '')}</td></tr>`).join('')}
        </tbody>
      </table>
      <p class="pvt-agg-copy" style="margin-top:10px;margin-bottom:0">PPP = probability a failure is a true positive; NPP = probability a pass is a true negative, computed by Bayes' theorem from the meta-analytic sens/spec and the selected base rate (Martin et al., 2020, Tables 16–17).</p>
    </div>`;
}

/* ---------- summary + APA export ---------- */
/* One row per reportable line. TOMM contributes one row per entered trial
   but counts as ONE indicator; EI and ES share their RBANS subtests and
   also count as one. The independence arithmetic is the point of the
   summary — see Larrabee (2014). */
function getPvtSummaryRows(){
  const rows = [];
  const ei = getPvtEi();
  if (ei.ei !== undefined) rows.push({
    id: 'ei', group: 'rbans', measure: 'RBANS Effort Index',
    score: String(ei.ei), cutoff: `> ${ei.cut}${ei.cutKey === 'sensitive' ? ' (screening)' : ''}`,
    sens: PVT_EI_ACCURACY[ei.cutKey].sens, spec: PVT_EI_ACCURACY[ei.cutKey].spec,
    result: pvtStatusWord(ei.fail), fail: ei.fail
  });
  const es = getPvtEs();
  /* The ES publishes no sensitivity/specificity pair (its discrimination is
     ROC AUC = .91) — the columns print dashes and the note explains. */
  if (es.gated) rows.push({
    id: 'es', group: 'rbans', measure: 'RBANS Effort Scale',
    score: '—', cutoff: `< ${PVT_ES.cutoff}`, sens: '—', spec: '—',
    result: 'Not computed (gate not met)', fail: false, gated: true
  });
  else if (es.es !== undefined) rows.push({
    id: 'es', group: 'rbans', measure: 'RBANS Effort Scale',
    score: String(es.es), cutoff: `< ${PVT_ES.cutoff}`, sens: '—', spec: '—',
    result: pvtStatusWord(es.fail), fail: es.fail
  });
  const rds = getPvtRds();
  if (rds.rds !== undefined){
    const acc = PVT_RDS_ACCURACY[rds.conservative ? 'conservative' : 'traditional'];
    rows.push({
      id: 'rds', group: 'rds', measure: 'Reliable Digit Span',
      score: String(rds.rds), cutoff: `≤ ${rds.cut}${rds.conservative ? ' (conservative)' : ''}`,
      sens: acc.sens, spec: acc.spec,
      result: pvtStatusWord(rds.fail), fail: rds.fail
    });
  }
  /* Same subtest as RDS, so the SAME group — one indicator between them,
     exactly as EI/ES share the RBANS group. */
  const ds = getPvtDs();
  if (ds.acss !== undefined){
    const acc = PVT_DS_ACCURACY[ds.conservative ? 'conservative' : 'sensitive'];
    rows.push({
      id: 'ds', group: 'rds', measure: 'Digit Span scaled score',
      score: String(ds.acss), cutoff: `≤ ${ds.cut}${ds.conservative ? ' (conservative)' : ''}`,
      sens: acc.sens, spec: acc.spec,
      result: pvtStatusWord(ds.fail), fail: ds.fail
    });
    if (ds.diff !== undefined) rows.push({
      id: 'ds-vocab', group: 'rds', measure: 'Vocabulary − Digit Span',
      score: String(ds.diff), cutoff: `≥ ${PVT_DS.vocabDiffCutoff}`,
      sens: '—', spec: '—',
      result: ds.diffFail ? 'Flagged' : 'Not flagged', fail: ds.diffFail
    });
    /* A withheld forward-span index (no age, or 55+) exports NO row — an
       unevaluated index in a report table would be a false claim. */
    if (ds.lsf !== undefined && ds.lsfState === 'ok') rows.push({
      id: 'ds-lsf', group: 'rds', measure: 'Longest span forward',
      score: String(ds.lsf), cutoff: `≤ ${PVT_DS.lsfCutoff} (age < ${PVT_DS.lsfAgeBelow})`,
      sens: '—', spec: '—',
      result: ds.lsfFail ? 'Flagged' : 'Not flagged', fail: ds.lsfFail
    });
    if (ds.lsb !== undefined) rows.push({
      id: 'ds-lsb', group: 'rds', measure: 'Longest span backward',
      score: String(ds.lsb), cutoff: `≤ ${PVT_DS.lsbCutoff}`,
      sens: '—', spec: '—',
      result: ds.lsbFail ? 'Flagged' : 'Not flagged', fail: ds.lsbFail
    });
  }
  /* Stand-alone: shares no subtest with anything above, so its own group. */
  const rey = getPvtRey();
  if (rey.recall !== undefined){
    rows.push({
      id: 'rey-recall', group: 'rey15', measure: 'Rey 15-Item free recall',
      score: String(rey.recall), cutoff: `< ${PVT_REY15.recallCutoff}`,
      sens: PVT_REY15_ACCURACY.recall.sens, spec: PVT_REY15_ACCURACY.recall.spec,
      result: pvtStatusWord(rey.recallFail), fail: rey.recallFail
    });
    if (rey.combo !== undefined) rows.push({
      id: 'rey-combo', group: 'rey15', measure: 'Rey 15-Item combination',
      score: String(rey.combo), cutoff: `< ${PVT_REY15.comboCutoff}`,
      sens: PVT_REY15_ACCURACY.combo.sens, spec: PVT_REY15_ACCURACY.combo.spec,
      result: pvtStatusWord(rey.comboFail), fail: rey.comboFail
    });
  }
  /* CVLT-3: a word list, sharing no subtest with anything above, so its own
     group. Hits and critical items come from ONE administration of the same
     trial, so they sit in that one group together — the same non-independence
     rule that collapses EI/ES and the digit-span indices. The cut-off column
     names the derived threshold and the criterion it came from, so an
     exported table never implies the manual published a cut score. */
  const cv = getPvtCvlt3();
  if (cv.hits !== undefined){
    rows.push({
      id: 'cvlt3', group: 'cvlt3', measure: 'CVLT-3 Forced Choice hits',
      score: `${cv.hits}/${PVT_CVLT3_FC_HITS.max}`,
      cutoff: cv.hitsThreshold === null ? '—'
        : cv.basis.cut === null ? `≤ ${cv.hitsThreshold} (base rate, ages ${cv.band})`
        : `≤ ${cv.hitsThreshold} (CVLT-II)`,
      sens: cv.basis.sens, spec: cv.basis.spec,
      result: pvtStatusWord(cv.hitsFail), fail: cv.hitsFail
    });
    [cv.critRecall, cv.critRecog].forEach((c, i) => {
      if (!c) return;
      rows.push({
        id: 'cvlt3-crit-' + (i === 0 ? 'recall' : 'recog'), group: 'cvlt3',
        measure: 'CVLT-3 ' + c.label.charAt(0).toLowerCase() + c.label.slice(1),
        score: String(c.count),
        cutoff: c.threshold === null ? '—' : `≥ ${c.threshold} (base rate, ages ${cv.band})`,
        sens: '—', spec: '—',
        result: c.fail ? 'Flagged' : 'Not flagged', fail: c.fail
      });
    });
  }
  const tomm = getPvtTomm();
  if (tomm.rows) tomm.rows.forEach(r => rows.push({
    id: 'tomm', group: 'tomm', measure: `TOMM ${r.label}`,
    score: String(r.score), cutoff: `< ${r.cutoff.cut}`,
    sens: r.cutoff.sensRange || r.cutoff.sens.toFixed(2).replace(/^0/, ''),
    spec: r.cutoff.specRange || r.cutoff.spec.toFixed(2).replace(/^0/, ''),
    result: pvtStatusWord(r.fail), fail: r.fail
  }));
  return rows;
}
function pvtIndicatorCounts(rows){
  const groups = [...new Set(rows.filter(r => !r.gated).map(r => r.group))];
  const failed = groups.filter(g => rows.some(r => r.group === g && r.fail));
  return { total: groups.length, failed: failed.length };
}

function pvtFmtRate(r){
  if (r === null || r === undefined) return '—';
  /* A published 0.0 is a real cell here (nobody in the standardisation
     sample), unlike the ToPF table where zero cells are omitted. Print it
     as "< 0.1%" rather than "0%": the sample is finite, so the honest claim
     is that the rate rounds to zero, not that the score cannot occur. */
  if (r === 0) return '&lt; 0.1%';
  /* One decimal always, as the manual prints it — "3%" beside a stored
     3.0 invites a reader to wonder which is the real figure. */
  return r.toFixed(1) + '%';
}
function renderPvtCvlt3(){
  const out = document.getElementById('pvt-cvlt3-result');
  if (!out) return;
  const s = getPvtCvlt3();
  if (s.empty){ out.innerHTML = pvtResultHtml('empty', 'Enter the Forced Choice total hits to look up its published base rate.', 'The result appears here and the outcome joins the running summary.'); return; }
  if (s.invalid){ out.innerHTML = pvtResultHtml('empty', PVT_PROMPTS.invalid); return; }
  if (s.partial){ out.innerHTML = pvtResultHtml('empty', 'Enter the Forced Choice total hits. The critical-item counts are read against it.'); return; }
  const derived = s.basis.cut === null;
  const rows = [{
    label: 'Forced Choice total hits', value: `${s.hits}/${PVT_CVLT3_FC_HITS.max}`,
    state: s.hitsFail ? 'fail' : 'pass',
    meta: `Base rate ${pvtFmtRate(s.hitsRate)} at ages ${s.band} (Table D.13)${
      s.hitsThreshold === null
        ? ' · no score in this band reaches the criterion'
        : derived
          ? ` · flags at ≤ ${s.hitsThreshold}, derived at the ${s.criterion.pct}% criterion`
          : ` · cut-off ≤ ${s.hitsThreshold} · sens. ${s.basis.sens} · spec. ${s.basis.spec} (CVLT-II)`}`
  }];
  [s.critRecall, s.critRecog].forEach(c => {
    if (!c) return;
    rows.push({
      label: c.label, value: String(c.count),
      state: c.fail ? 'fail' : 'pass',
      meta: `Base rate ${pvtFmtRate(c.rate)} at ages ${s.band} (Table ${PVT_CVLT3_CRITICAL[c.which].table})${
        c.threshold === null ? ' · no count in this band reaches the criterion'
                             : ` · flags at ≥ ${c.threshold}`}`
    });
  });
  /* Which instrument the flag came from is the one thing a reader of the
     report must not have to guess, so it is stated whenever the flag rests
     on a borrowed cut-off rather than on the CVLT-3's own table. */
  const basisNote = derived ? '' :
    ` The flag rests on the CVLT-II cut-off ≤ ${s.basis.cut}, not on a CVLT-3 figure; the base rate beside it is the CVLT-3 manual's.`;
  const note = `${s.usedAllAges
      ? (s.age === null ? 'No patient age entered, so the manual&rsquo;s All ages column is used — enter an age in the top bar for the banded reading, which can differ sharply.'
                        : `Age ${s.age} falls outside the CVLT-3 normative range (16&ndash;90), so the All ages column is used.`)
      : `Read against the ${s.band} band for the entered age of ${s.age}.`
    }${basisNote} A poor forced-choice score is a strong indicator of exaggeration, but a perfect or near-perfect one does not rule it out (Delis et al., 2017).`;
  out.innerHTML = pvtReadoutHtml(rows, note);
}

function renderPvtSummary(){
  const host = document.getElementById('pvt-summary-body');
  if (!host) return;
  const rows = getPvtSummaryRows();
  if (rows.length === 0){
    host.innerHTML = '<div class="pvt-result">' + pvtResultHtml('empty', 'Nothing entered yet. Score at least one measure on the other tabs and the summary builds itself. The tab strip above tracks each measure as you go.') + '</div>';
    return;
  }
  const c = pvtIndicatorCounts(rows);
  const countStat = c.total === 0 ? '' :
    `<div class="pvt-count">
      <span class="pvt-count-num${c.failed > 0 ? ' is-fail' : ''}">${c.failed}<span style="color:var(--faint)">/</span>${c.total}</span>
      <span class="pvt-count-label">independent indicator${c.total === 1 ? '' : 's'} beyond the selected cut-off${c.failed === 1 ? '' : 's'}. The two RBANS indices count as one, the digit-span indices as one, the CVLT-3 forced-choice scores as one, and the TOMM trials as one. A single failure is a hypothesis, not a conclusion (Larrabee, 2014).</span>
    </div>`;
  host.innerHTML = `
    <div class="pvt-card">
      <div class="pvt-card-kicker">Indicators scored this session</div>
      ${countStat}
      <table class="pvt-table">
        <thead><tr><th>Measure</th><th>Score</th><th>Cut-off</th><th>Sens.</th><th>Spec.</th><th>Result</th></tr></thead>
        <tbody>${rows.map(r => `<tr><td>${r.measure}</td><td>${r.score}</td><td>${r.cutoff.replace('<', '&lt;')}</td><td>${r.sens}</td><td>${r.spec}</td><td class="${r.fail ? 'pvt-cell-fail' : ''}">${r.result}</td></tr>`).join('')}</tbody>
      </table>
    </div>`;
}

function renderPvtApa(){
  const out = document.getElementById('pvt-apa');
  if (!out) return;
  const rows = getPvtSummaryRows();
  /* EMPTY-STATE GUARD — same contract §35 enforces on the premorbid
     renderers: the working-report observer's only test is the presence of
     .apa-table, so a renderer that always writes one always offers an
     empty table to the report. No data, no table. */
  if (rows.length === 0){
    out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">Enter at least one validity measure to preview.</div>';
    return;
  }
  const ei = getPvtEi(), es = getPvtEs();
  out.innerHTML = `
    <div class="apa-table-num">Table 1</div>
    <div class="apa-table-title">Performance validity indicators</div>
    <table class="apa-table">
      <thead><tr><th>Measure</th><th class="num">Score</th><th>Cut-off</th><th class="num">Sens.</th><th class="num">Spec.</th><th>Result</th></tr></thead>
      <tbody>${rows.map(r => `<tr><td>${escapeHtml(r.measure)}</td><td class="num">${escapeHtml(r.score)}</td><td>${escapeHtml(r.cutoff)}</td><td class="num">${escapeHtml(r.sens)}</td><td class="num">${escapeHtml(r.spec)}</td><td>${escapeHtml(r.result)}</td></tr>`).join('')}</tbody>
    </table>
    ${apaNoteHtml('pvt', {
      hasEi:   rows.some(r => r.id === 'ei'),
      hasEs:   rows.some(r => r.id === 'es'),
      hasRds:  rows.some(r => r.id === 'rds'),
      hasDs:   rows.some(r => r.id.startsWith('ds')),
      hasRey:  rows.some(r => r.id.startsWith('rey')),
      hasCvlt3: rows.some(r => r.group === 'cvlt3'),
      cvlt3Borrowed: rows.some(r => r.group === 'cvlt3') && (typeof getPvtCvlt3 === 'function') && getPvtCvlt3().basis?.cut !== null && getPvtCvlt3().basis !== undefined,
      cvlt3Cite: (getPvtCvlt3().basis?.cut === 15) ? 'Erdodi et al., 2018' : 'Schwartz et al., 2016',
      /* One sentence per measure present that carries a version mismatch,
         drawn from the same PVT_INSTRUMENTS the cards render. The CVLT-3's
         is stated above in its own clause, so it is not repeated here. */
      versionCaveats: (typeof PVT_INSTRUMENTS === 'object' ? Object.keys(PVT_INSTRUMENTS) : [])
        .filter(tab => PVT_INSTRUMENTS[tab].mismatch && tab !== 'cvlt3')
        .filter(tab => rows.some(r => r.id === tab || r.id.startsWith(tab + '-') || (tab === 'ds' && r.id.startsWith('ds'))))
        .map(tab => PVT_INSTRUMENTS[tab].mismatch),
      hasTomm: rows.some(r => r.group === 'tomm'),
      esGated: !!es.gated,
      bothRbans:     rows.some(r => r.id === 'ei')  && rows.some(r => r.id === 'es'),
      bothDigitSpan: rows.some(r => r.id === 'rds') && rows.some(r => r.id.startsWith('ds')),
      /* Any row whose accuracy columns print a dash — a base-rate index or
         the Effort Scale — so the note explains the dashes only when some
         are on the page. */
      hasDashes: rows.some(r => r.sens === '—')
    })}
  `;
}

/* The tab strip's live status chips. One chip per measure, restated from
   the same getPvt* state the result cards render, so strip and card cannot
   disagree. The summary chip carries the independent-indicator count. */
/* An unscored tab shows NO chip at all: a dash in a pill reads as seven
   pieces of dead chrome before anything is entered, so the pill appears
   only once it has a state to report. Scored chips carry the same glyphs
   as the readout rows (shape + word + colour, never colour alone). */
function pvtChip(el, text, kind){
  if (!el) return;
  el.classList.remove('is-pass', 'is-fail', 'is-na');
  el.classList.toggle('is-idle', !kind);
  if (!kind){ el.textContent = ''; return; }
  const icon = PVT_STATE_ICON[kind === 'fail' ? 'flag' : kind === 'na' ? 'na' : 'pass'] || '';
  el.innerHTML = icon + escapeHtml(text);
  el.classList.add('is-' + kind);
}

/* The caution under each result is neutral reference text until the method
   on that tab actually flags — then it takes the amber emphasis, because
   "corroborate before concluding" is advice about a result that now
   exists. Driven from the same states the chips render, in one place. */
function pvtCautionFlag(tab, flagged){
  document.querySelector(`#pvt-${tab} .pvt-caution`)?.classList.toggle('is-flagged', !!flagged);
}
function renderPvtNav(){
  const ei = getPvtEi();
  pvtChip(document.getElementById('pvt-status-ei'),
    ei.ei !== undefined ? (ei.fail ? 'Fail' : 'Pass') : '—',
    ei.ei !== undefined ? (ei.fail ? 'fail' : 'pass') : null);
  const es = getPvtEs();
  pvtChip(document.getElementById('pvt-status-es'),
    es.gated ? 'Gated' : es.es !== undefined ? (es.fail ? 'Fail' : 'Pass') : '—',
    es.gated ? 'na' : es.es !== undefined ? (es.fail ? 'fail' : 'pass') : null);
  const rds = getPvtRds();
  pvtChip(document.getElementById('pvt-status-rds'),
    rds.rds !== undefined ? (rds.fail ? 'Fail' : 'Pass') : '—',
    rds.rds !== undefined ? (rds.fail ? 'fail' : 'pass') : null);
  const dsS = getPvtDs();
  const dsAny = dsS.acss !== undefined;
  const dsFail = dsAny && (dsS.fail || dsS.diffFail || dsS.lsfFail || dsS.lsbFail);
  pvtChip(document.getElementById('pvt-status-ds'),
    dsAny ? (dsFail ? 'Fail' : 'Pass') : '—',
    dsAny ? (dsFail ? 'fail' : 'pass') : null);
  const reyS = getPvtRey();
  const reyAny = reyS.recall !== undefined;
  const reyFail = reyAny && (reyS.recallFail || !!reyS.comboFail);
  pvtChip(document.getElementById('pvt-status-rey15'),
    reyAny ? (reyFail ? 'Fail' : 'Pass') : '—',
    reyAny ? (reyFail ? 'fail' : 'pass') : null);
  const cvS = getPvtCvlt3();
  const cvAny = cvS.hits !== undefined;
  pvtChip(document.getElementById('pvt-status-cvlt3'),
    cvAny ? (cvS.anyFail ? 'Fail' : 'Pass') : '—',
    cvAny ? (cvS.anyFail ? 'fail' : 'pass') : null);
  const tomm = getPvtTomm();
  const tommHas = !!(tomm.rows && tomm.rows.length);
  pvtChip(document.getElementById('pvt-status-tomm'),
    tommHas ? (tomm.anyFail ? 'Fail' : 'Pass') : '—',
    tommHas ? (tomm.anyFail ? 'fail' : 'pass') : null);
  const c = pvtIndicatorCounts(getPvtSummaryRows());
  pvtChip(document.getElementById('pvt-status-summary'),
    c.total > 0 ? `${c.failed}/${c.total}` : '—',
    c.total > 0 ? (c.failed > 0 ? 'fail' : 'pass') : null);
  pvtCautionFlag('ei', ei.ei !== undefined && ei.fail);
  pvtCautionFlag('es', es.es !== undefined && !es.gated && es.fail);
  pvtCautionFlag('rds', rds.rds !== undefined && rds.fail);
  pvtCautionFlag('ds', dsFail);
  pvtCautionFlag('rey15', reyFail);
  pvtCautionFlag('cvlt3', cvAny && cvS.anyFail);
  pvtCautionFlag('tomm', tommHas && tomm.anyFail);
}

/* The published-accuracy line beside each cut-off select — same strings the
   summary table and APA export print, so screen and export cannot differ. */
/* The "Derived on" line under every method's citation, and the mismatch
   warning where the edition a clinician is administering is not the one the
   cut-off came from. Both are written into placeholders in the markup from
   PVT_INSTRUMENTS, so the cards and the About roster cannot drift. */
function pvtInstrumentLineHtml(tab){
  const i = PVT_INSTRUMENTS[tab];
  if (!i) return '';
  return `<span class="pvt-derived"><span class="pvt-derived-label">Derived on</span> ${i.derived}${
    i.unspecified ? ' <span class="pvt-derived-gap">— form/edition not recorded here</span>' : ''}</span>${
    i.mismatch ? `<span class="pvt-derived-warn">${i.mismatch}</span>` : ''}`;
}
function renderPvtInstruments(){
  document.querySelectorAll('#validity [data-pvt-derived]').forEach(el => {
    const html = pvtInstrumentLineHtml(el.dataset.pvtDerived);
    if (el.innerHTML !== html) el.innerHTML = html;
  });
}
function renderPvtAccuracy(){
  const eiEl = document.getElementById('pvt-ei-accuracy');
  if (eiEl){
    const key = document.getElementById('pvt-ei-cutoff')?.value === 'sensitive' ? 'sensitive' : 'standard';
    eiEl.textContent = `Published accuracy at this cut-off: sens. ${PVT_EI_ACCURACY[key].sens} · spec. ${PVT_EI_ACCURACY[key].spec} (Silverberg et al., 2007)`;
  }
  const dsEl = document.getElementById('pvt-ds-accuracy');
  if (dsEl){
    const key = document.getElementById('pvt-ds-cutoff')?.value === 'sensitive' ? 'sensitive' : 'conservative';
    dsEl.textContent = key === 'conservative'
      ? `Published accuracy at this cut-off: sens. ${PVT_DS_ACCURACY.conservative.sens} · spec. ${PVT_DS_ACCURACY.conservative.spec} (Axelrod et al., 2006); base rate 3.8% standardisation / 3.4% clinical (Iverson & Tulsky, 2003)`
      : `Published accuracy at this cut-off: sens. ${PVT_DS_ACCURACY.sensitive.sens} · spec. ${PVT_DS_ACCURACY.sensitive.spec} (Axelrod et al., 2006)`;
  }
  const cvEl = document.getElementById('pvt-cvlt3-accuracy');
  if (cvEl) cvEl.textContent = pvtCvlt3Basis().cite;
  const rdsEl = document.getElementById('pvt-rds-accuracy');
  if (rdsEl){
    const key = document.getElementById('pvt-rds-cutoff')?.value === 'traditional' ? 'traditional' : 'conservative';
    rdsEl.textContent = `Published accuracy at this cut-off: sens. ${PVT_RDS_ACCURACY[key].sens} · spec. ${PVT_RDS_ACCURACY[key].spec} (Schroeder et al., 2012)`;
  }
}

/* ---------- the running rail ----------
   A single element beside every method panel, so the aggregate picture is
   in view while the clinician is still entering scores rather than one
   tab away. It reads the SAME getPvtSummaryRows/pvtIndicatorCounts the
   Summary tab and the APA export read — three surfaces, one derivation,
   which is why a measure can never appear scored here and unscored there. */
/* The rail is grouped by INDEPENDENT INDICATOR, not by measure, because
   the count above it is over indicators: the two RBANS indices are one,
   and the two digit-span measures are one. Listing six flat rows beside a
   "3 of 4" count invited exactly the arithmetic the page keeps warning
   against, so the grouping now shows on screen what the note says in
   words. Rows are buttons: the rail doubles as navigation. */
const PVT_RAIL_GROUPS = [
  { id: 'rbans', label: 'RBANS', members: [
    { tab: 'ei',    label: 'Effort Index', ids: ['ei'] },
    { tab: 'es',    label: 'Effort Scale', ids: ['es'] }
  ]},
  { id: 'rds', label: 'Digit span', members: [
    { tab: 'rds',   label: 'Reliable Digit Span', ids: ['rds'] },
    { tab: 'ds',    label: 'Digit Span indices',  ids: ['ds', 'ds-vocab', 'ds-lsf', 'ds-lsb'] }
  ]},
  { id: 'rey15', label: 'Rey 15-Item', members: [
    { tab: 'rey15', label: 'Recall & recognition', ids: ['rey-recall', 'rey-combo'] }
  ]},
  { id: 'cvlt3', label: 'CVLT-3', members: [
    { tab: 'cvlt3', label: 'Forced choice', ids: ['cvlt3', 'cvlt3-crit-recall', 'cvlt3-crit-recog'] }
  ]},
  { id: 'tomm', label: 'TOMM', members: [
    { tab: 'tomm',  label: 'Forced choice', ids: ['tomm'] }
  ]}
];
const PVT_RAIL_ICON = {
  flag: PVT_STATE_ICON.flag,
  pass: PVT_STATE_ICON.pass,
  gated: PVT_STATE_ICON.na,
  idle: '<svg class="pvt-flag-icon" viewBox="0 0 12 12" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.3"><circle cx="6" cy="6" r="3.4"/></svg>'
};
function renderPvtRail(){
  const rail = document.getElementById('pvt-rail');
  if (!rail) return;
  const rows = getPvtSummaryRows();
  const c = pvtIndicatorCounts(rows);
  const current = document.querySelector('#validity .pvt-method-tab.active')?.dataset.pvtTab;

  /* Status reads as words in a right-hand column, the way a progress panel
     states a value against its label — colour marks a flag rather than
     carrying the meaning on its own. */
  let scored = 0, total = 0;
  const groups = PVT_RAIL_GROUPS.map(g => {
    const members = g.members.map(m => {
      total++;
      const mine = rows.filter(r => m.ids.includes(r.id));
      const gated = mine.some(r => r.gated);
      if (mine.length) scored++;
      const flagged = mine.filter(r => r.fail).length;
      const state = gated ? 'gated' : !mine.length ? 'idle' : flagged ? 'flag' : 'pass';
      const value = gated ? 'gate not met'
        : !mine.length ? 'not scored'
        : flagged ? `${flagged} flagged`
        : 'within cut-offs';
      return `<button type="button" class="pvt-rail-row is-${state}${m.tab === current ? ' is-current' : ''}" data-rail-tab="${m.tab}">
        <span class="pvt-rail-label">${m.label}</span>
        <span class="pvt-rail-value">${value}</span>
      </button>`;
    }).join('');
    return `<div class="pvt-rail-group">
      <div class="pvt-rail-group-head">
        <span class="pvt-rail-group-name">${g.label}</span>
        ${g.members.length > 1 ? '<span class="pvt-rail-group-tag">counts as one</span>' : ''}
      </div>${members}</div>`;
  }).join('');

  /* With nothing scored the section shows the bare 0/6 count alone — the
     workspace card already carries the "enter scores" instruction, and the
     Flagged stat is withheld rather than printed as 0/0. */
  const hint = c.failed >= 2 ? 'Two or more independent failures support probable invalidity (Larrabee, 2014).'
    : c.failed === 1 ? 'A single failure is a hypothesis to corroborate, not a conclusion.'
    : 'Nothing beyond its cut-off so far.';

  /* The collapse state lives on the rail, which survives this re-render;
     the card is rebuilt from it so screen and ARIA cannot drift. */
  const collapsed = rail.classList.contains('is-collapsed');
  rail.innerHTML = `<div class="pvt-rail-card${collapsed ? ' is-collapsed' : ''}">
    <div class="pvt-rail-head">
      <span class="pvt-rail-title">Progress</span>
      <button type="button" class="pvt-rail-toggle" data-rail-toggle aria-label="${collapsed ? 'Expand' : 'Collapse'} summary" aria-expanded="${collapsed ? 'false' : 'true'}">
        <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><polyline points="6,3 11,8 6,13"/></svg>
      </button>
    </div>
    <div class="pvt-rail-body">
      <div class="pvt-rail-section">
        <div class="pvt-rail-kicker">Indicators</div>
        ${c.total > 0 ? `<div class="pvt-rail-stat${c.failed ? ' is-flagged' : ''}">
          <span class="pvt-rail-stat-label">Flagged</span>
          <span class="pvt-rail-stat-value">${c.failed}<span class="pvt-rail-sep">/</span>${c.total}</span>
        </div>` : ''}
        <div class="pvt-rail-stat">
          <span class="pvt-rail-stat-label">Measures scored</span>
          <span class="pvt-rail-stat-value">${scored}<span class="pvt-rail-sep">/</span>${total}</span>
        </div>
        ${c.total > 0 ? `<p class="pvt-rail-hint">${hint}</p>` : ''}
      </div>
      <div class="pvt-rail-section">
        <div class="pvt-rail-kicker">By indicator</div>
        ${groups}
      </div>
      <button type="button" class="pvt-rail-cta" data-rail-tab="summary">Full summary &amp; APA table →</button>
    </div>
  </div>`;
}

function renderPvtAll(){
  renderPvtEi();
  renderPvtEs();
  renderPvtRds();
  renderPvtDs();
  renderPvtRey();
  renderPvtTomm();
  renderPvtCvlt3();
  renderPvtAccuracy();
  renderPvtInstruments();
  renderPvtNav();
  renderPvtRail();
  renderPvtSummary();
  renderPvtApa();
}

function switchPvtTab(name){
  document.querySelectorAll('#validity .pvt-method-tab').forEach(t => {
    const on = t.dataset.pvtTab === name;
    t.classList.toggle('active', on);
    t.setAttribute('aria-selected', on ? 'true' : 'false');
  });
  document.querySelectorAll('#validity .pvt-tab-content').forEach(c =>
    c.classList.toggle('active', c.id === 'pvt-' + name));
  /* The rail would restate the Summary tab's own table beside it. */
  document.querySelector('#validity .pvt-sheet')?.classList.toggle('is-summary', name === 'summary');
  renderPvtRail();
}

function clearPvt(){
  Object.keys(pvtState).forEach(k => { pvtState[k] = ''; });
  document.querySelectorAll('#validity [data-pvt-field]').forEach(inp => { inp.value = ''; });
  ['pvt-rds-f','pvt-rds-b','pvt-ds-acss','pvt-ds-vocab','pvt-ds-lsf','pvt-ds-lsb','pvt-rey-recall','pvt-rey-recog','pvt-rey-fp','pvt-tomm-t1','pvt-tomm-t2','pvt-tomm-ret','pvt-cvlt3-hits','pvt-cvlt3-crit-recall','pvt-cvlt3-crit-recog'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = '';
  });
  renderPvtAll();
}

/* ---------- the About landing tab ----------
   The page lands here rather than on the Effort Index, which was only
   ever first by markup order. Same idiom as the Change Analysis overview:
   one row per measure, click to jump, a Counts-as column carrying the
   independence grouping that the ≥ 2-failure rule depends on. Accuracy
   cells are rendered from the SAME constants the method tabs print
   (PVT_*_ACCURACY, PVT_TOMM_CUTOFFS), so this table cannot drift from
   the pages it describes. */
function renderPvtAboutPanel(){
  const el = document.getElementById('pvt-about-body');
  if (!el) return;
  const t2 = PVT_TOMM_CUTOFFS.find(c => c.id === 't2-45');
  const rows = [
    { tab: 'ei', title: 'Effort Index', cite: 'Silverberg et al. (2007)',
      source: 'RBANS embedded', group: 'RBANS', cut: 'EI &gt; 3',
      sens: PVT_EI_ACCURACY.standard.sens, spec: PVT_EI_ACCURACY.standard.spec,
      desc: 'Weighted sum of the Digit Span and List Recognition raw scores (0–12); higher = less credible.' },
    { tab: 'es', title: 'Effort Scale', cite: 'Novitski et al. (2012)',
      source: 'RBANS embedded', group: 'RBANS', cut: 'ES &lt; 12',
      sens: null, spec: null, acc: `ROC AUC ${PVT_ES_ACCURACY.auc}`,
      desc: 'List Recognition minus the three recall scores, plus Digit Span, all raw; computed only where its gate is met, because an ungated ES over-flags intact examinees.' },
    { tab: 'rds', title: 'Reliable Digit Span', cite: 'Greiffenstein et al. (1994); Schroeder et al. (2012)',
      source: 'WAIS embedded', group: 'Digit span', cut: 'RDS ≤ 6',
      sens: PVT_RDS_ACCURACY.conservative.sens, spec: PVT_RDS_ACCURACY.conservative.spec,
      desc: 'Longest forward span plus longest backward span passed on both trials.' },
    { tab: 'ds', title: 'Digit Span indices', cite: 'Iverson &amp; Tulsky (2003); Axelrod et al. (2006)',
      source: 'WAIS embedded', group: 'Digit span', cut: 'ACSS ≤ 5',
      sens: PVT_DS_ACCURACY.conservative.sens, spec: PVT_DS_ACCURACY.conservative.spec,
      desc: 'Age-corrected scaled-score cut-off, with Vocabulary − Digit Span and longest-span base rates alongside.' },
    { tab: 'rey15', title: 'Rey 15-Item', cite: 'Boone et al. (2002)',
      source: 'Stand-alone', group: 'Rey 15-Item', cut: 'Recall + recognition',
      sens: PVT_REY15_ACCURACY.combo.sens, spec: PVT_REY15_ACCURACY.combo.spec,
      desc: 'Free recall of fifteen over-learned items, with a recognition trial that raises sensitivity.' },
    { tab: 'cvlt3', title: 'CVLT-3 Forced Choice', cite: 'Delis et al. (2017); Erdodi et al. (2018)',
      source: 'Embedded', group: 'CVLT-3', cut: 'Base rate by age',
      sens: null, spec: null, acc: 'no pair published for the base-rate reading',
      desc: 'Total hits on the Forced Choice trial, read against the CVLT-3 manual’s age-banded base rates, with the CVLT-II cut-offs selectable alongside.' },
    { tab: 'tomm', title: 'TOMM', cite: 'Tombaugh (1996); Martin et al. (2020)',
      source: 'Stand-alone', group: 'TOMM', cut: 'Trial 2 &lt; 45',
      sens: t2 ? (t2.sensRange || String(t2.sens).replace('0.', '.')) : '—',
      spec: t2 ? (t2.specRange || String(t2.spec).replace('0.', '.')) : '—',
      desc: 'Fifty-item forced-choice picture recognition; robust to most genuine impairment, though specificity falls in dementia.' }
  ];
  const body = rows.map(r => `<tr class="pvt-overview-row" data-about-tab="${r.tab}" tabindex="0" role="button" aria-label="Open ${r.title}">
    <td class="pvt-overview-measure">
      <div class="pvt-overview-titlerow">
        <span class="pvt-overview-title">${r.title}</span>
        <button type="button" class="pvt-overview-info" data-pvtip="${r.desc}" aria-label="More about ${r.title}" tabindex="0">?</button>
      </div>
      <span class="pvt-overview-cite">${r.cite}</span>
    </td>
    <td class="pvt-overview-cell pvt-overview-derived">${(PVT_INSTRUMENTS[r.tab] || {}).derived || r.source}${
      (PVT_INSTRUMENTS[r.tab] || {}).unspecified ? '<span class="pvt-overview-sub">edition not recorded</span>' : ''}${
      (PVT_INSTRUMENTS[r.tab] || {}).mismatch ? '<span class="pvt-overview-sub is-warn">version caveat, see the measure</span>' : ''}<span class="pvt-overview-sub">${
      (PVT_INSTRUMENTS[r.tab] || {}).kind || r.source}</span></td>
    <td class="pvt-overview-cell">${r.group}</td>
    <td class="pvt-overview-cell pvt-overview-acc">${r.acc
      ? `<span class="pvt-overview-cut">${r.cut}</span>${r.acc}`
      : `<span class="pvt-overview-cut">${r.cut}</span>sens. ${r.sens} · spec. ${r.spec}`}</td>
    <td class="pvt-overview-arrow" aria-hidden="true"><svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><line x1="3" y1="8" x2="13" y2="8"/><polyline points="9,4 13,8 9,12"/></svg></td>
  </tr>`).join('');
  el.innerHTML = `<div class="pvt-overview-wrap">
    <table class="pvt-overview-table">
      <thead><tr>
        <th class="pvt-overview-th is-measure">Measure</th>
        <th class="pvt-overview-th"><span class="pvt-overview-colh" tabindex="0" data-pvtip="The instrument and edition each cut-off was calibrated on. Embedded indices are computed from subtests that also measure genuine ability; stand-alone tests are administered solely to assess performance validity. A cut-off derived on one edition does not automatically transfer to another.">Derived on</span></th>
        <th class="pvt-overview-th"><span class="pvt-overview-colh" tabindex="0" data-pvtip="Measures derived from the same administration of the same instrument share error and are not independent evidence: each named group counts as ONE indicator in the aggregation.">Counts as</span></th>
        <th class="pvt-overview-th"><span class="pvt-overview-colh is-end" tabindex="0" data-pvtip="Published accuracy at the cut-off named in the cell — the same figures printed beside each measure's cut-off selector.">Accuracy at default cut-off</span></th>
        <th aria-hidden="true"></th>
      </tr></thead>
      <tbody>${body}</tbody>
    </table>
    <p class="pvt-overview-foot">Failing <strong>two or more independent</strong> indicators supports probable invalidity (Larrabee, 2014). Measures sharing an instrument count as one between them, which is what the running summary's "counts as one" tags mean. The Summary tab reports Larrabee's published classification accuracy at each failure count. No single index is a verdict.</p>
  </div>`;
  el.querySelectorAll('[data-about-tab]').forEach(row => {
    const go = () => switchPvtTab(row.dataset.aboutTab);
    row.addEventListener('click', e => { if (e.target.closest('.pvt-overview-info')) return; go(); });
    row.addEventListener('keydown', e => {
      if (e.target.closest('.pvt-overview-info')) return;
      if (e.key === 'Enter' || e.key === ' '){ e.preventDefault(); go(); }
    });
  });
}

function setupPvtPage(){
  const root = document.getElementById('validity');
  if (!root) return;
  renderPvtAboutPanel();
  root.querySelectorAll('.pvt-method-tab').forEach(tab => {
    tab.addEventListener('click', () => switchPvtTab(tab.dataset.pvtTab));
  });
  /* The shared RBANS fields appear on both the EI and ES tabs; pvtState is
     the master and every input with the same data-pvt-field mirrors it. */
  root.querySelectorAll('[data-pvt-field]').forEach(inp => {
    inp.addEventListener('input', () => {
      const field = inp.dataset.pvtField;
      pvtState[field] = inp.value;
      root.querySelectorAll(`[data-pvt-field="${field}"]`).forEach(other => {
        if (other !== inp && other.value !== inp.value) other.value = inp.value;
      });
      renderPvtAll();
    });
  });
  ['pvt-rds-f','pvt-rds-b','pvt-ds-acss','pvt-ds-vocab','pvt-ds-lsf','pvt-ds-lsb','pvt-rey-recall','pvt-rey-recog','pvt-rey-fp','pvt-tomm-t1','pvt-tomm-t2','pvt-tomm-ret','pvt-cvlt3-hits','pvt-cvlt3-crit-recall','pvt-cvlt3-crit-recog'].forEach(id => {
    document.getElementById(id)?.addEventListener('input', renderPvtAll);
  });
  /* The forward-span index reads the shared top-bar age, so an age typed
     anywhere must re-evaluate it. Safe against the §35 hazard: the APA
     container only ever holds a table once a score is entered here, so an
     age keystroke alone cannot conjure one. */
  ['patient-age','pre-age'].forEach(id => {
    document.getElementById(id)?.addEventListener('input', renderPvtAll);
  });
  document.getElementById('pvt-rey-administer')?.addEventListener('click', reyAdminOpen);
  /* The rail doubles as navigation — delegated, because it re-renders on
     every keystroke and bound handlers would not survive. */
  document.getElementById('pvt-rail')?.addEventListener('click', e => {
    if (e.target.closest('[data-rail-toggle]')){
      const card = document.querySelector('#pvt-rail .pvt-rail-card');
      const btn = e.target.closest('[data-rail-toggle]');
      const collapsed = card?.classList.toggle('is-collapsed');
      btn.setAttribute('aria-expanded', collapsed ? 'false' : 'true');
      btn.setAttribute('aria-label', collapsed ? 'Expand summary' : 'Collapse summary');
      /* The state must survive the next re-render, which happens on every
         keystroke — so it lives on the rail, not only on the card. */
      document.getElementById('pvt-rail')?.classList.toggle('is-collapsed', !!collapsed);
      return;
    }
    const row = e.target.closest('[data-rail-tab]');
    if (row) switchPvtTab(row.dataset.railTab);
  });
  ['pvt-ei-cutoff','pvt-rds-cutoff','pvt-ds-cutoff','pvt-tomm-t1cut','pvt-tomm-t2cut','pvt-tomm-br','pvt-cvlt3-criterion','pvt-cvlt3-basis'].forEach(id => {
    document.getElementById(id)?.addEventListener('change', renderPvtAll);
  });
  renderPvtAll();
}


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

// Performance Validity page
setupPvtPage();

// Final initialization
refreshAll();

/* ---------- GLOBAL CLEAR — Topbar "New patient" button ----------
   Wipes every tool's session data in one action: battery rows, SDI rows,
   the four RCI methods, premorbid inputs, working report items.

   The auto-add MutationObserver on each APA container is suppressed during
   the clear cascade so the residual table HTML doesn't sneak items back
   into the bundle after we wipe it. */
(function wireGlobalClear(){
  const btn = document.getElementById('topbar-clear-all');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const ok = confirm('Start a new patient?\n\nThis clears every table you\'ve been working on across all tools, including the Working Report, and the patient age. It cannot be undone.');
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
    // Performance Validity page
    try { if (typeof clearPvt === 'function') clearPvt(); } catch(e){}
    /* Premorbid inputs (predictors + demographics), and the master patient age.
       #patient-age is listed explicitly rather than left to the #pre-age mirror:
       the mirror would in fact carry the blank across, but a patient identity
       that survives "New patient" because one listener was reordered is the
       failure this button exists to prevent. Clear both, assert neither. */
    try {
      ['patient-age','pre-topf','pre-vc','pre-mr','pre-sex','pre-occ','pre-edu','pre-age'].forEach(id => {
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
    charts: 'charts',
    validity: 'validity',
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
    charts: 'Score Charts',
    validity: 'Performance Validity',
    premorbid: 'Premorbid Estimation',
    'custom-tests': 'Data',
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
   - Report items are SESSION-ONLY (sessionStorage): a reload keeps them, but
     closing the tab, window or browser clears them — they are patient data.
     Only view preferences persist across sessions, in localStorage.
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
    'pre-opiepredict-apa':  'Premorbid · OPIE-4 Predicted',
    'pvt-apa':              'Performance Validity'
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
    'pre-opiepredict-apa':  'OPIE-4-Predicted vs Achieved',
    'pvt-apa':              'Performance Validity Indicators'
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
    /* Validity tables mention RBANS / Digit Span / TOMM as the indicators'
       host tests, not as a test family whose results these are — titling the
       table "Performance Validity Indicators: RBANS" would misdescribe it. */
    if (parentId && (parentId.startsWith('pre-') || parentId.startsWith('pvt-'))){
      return method || 'APA Table';
    }

    const family = explicitFamily || detectTestFamily(html);
    if (method && family) return `${method}: ${family}`;
    if (method) return method;
    return family || 'APA Table';
  }
  const SOURCE_IDS = Object.keys(SOURCE_LABELS);

  let state = { items: [], minimized: true, drawerWidth: null, onboardingSeen: false, maximised: false, consent: {} };
  /* ---------- consent-gated sources ----------
     The four Change Analysis methods share ONE row set (RCI_SHARED_ROWS), so a
     table entered on one method is already complete on the other three: the
     moment anything on a sibling tab caused a re-render — the CI selector, a
     date, an autofill — its APA table appeared and the observer below collected
     it, putting up to four near-identical tables in the report from one entry.

     So these four do not auto-collect. A method is collected once it has been
     ACCEPTED, which happens either because the clinician entered data in it
     (rciMarkMethodUsed, which is what makes the method they are actually
     working in behave exactly as before) or because they pressed the "Add to
     report" control on that method's own APA toolbar. Consent is per method and
     persists with the bundle, so an accepted method goes on updating silently.

     Nothing else in the app is gated: every other tool owns its own data, so
     its table appearing IS the clinician's intent. */
  const CONSENT_SOURCES = new Set(['rci-basic-apa','rci-practice-apa','rci-srb-apa','rci-crawford-apa']);
  /* Split items are keyed "<sourceId>::<family>" — consent lives on the parent. */
  function consentParent(sourceId){ return String(sourceId || '').split('::')[0]; }
  function isConsentGated(sourceId){ return CONSENT_SOURCES.has(consentParent(sourceId)); }
  function isSourceAccepted(sourceId){
    const parent = consentParent(sourceId);
    if (!CONSENT_SOURCES.has(parent)) return true;
    return state.consent[parent] === true;
  }
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

  /* ---------- persistence ----------

     PATIENT DATA IS SESSION-ONLY. The report items (and the consent map that
     belongs to them) are patient data, and closing the tab, window or browser
     must clear them — so they live in sessionStorage, which survives a reload
     of the same tab and nothing else. Only the non-patient VIEW PREFERENCES
     (drawer minimised, drawer width, onboarding seen) persist in localStorage,
     under their own key, so the drawer does not re-onboard a clinician every
     morning.

     The bundle used to keep everything — items included — in
     localStorage[STORAGE_KEY]. Existing installs therefore have patient data
     sitting on disk, which is exactly what this change exists to stop: load()
     migrates a legacy blob into the session once, then DELETES the
     localStorage copy unconditionally. Merely switching the writes to
     sessionStorage would have left every previously saved report on disk
     indefinitely.

     sessionStorage is per-tab, so two tabs now hold two independent reports
     where localStorage shared one. That is the correct reading of "session":
     each tab is its own patient session. */
  const PREFS_KEY = 'workingReportPrefs_v1';
  function load(){
    try {
      const rawPrefs = localStorage.getItem(PREFS_KEY);
      if (rawPrefs){
        const p = JSON.parse(rawPrefs);
        const w = Number(p.drawerWidth);
        state.minimized = p.minimized !== false;
        state.drawerWidth = (Number.isFinite(w) && w >= 320) ? w : null;
        state.onboardingSeen = p.onboardingSeen === true;
      }
    } catch(e){}
    try {
      /* The session blob first; a legacy localStorage bundle only when the
         session holds nothing (first load after the change, or an old tab). */
      const raw = sessionStorage.getItem(STORAGE_KEY) || localStorage.getItem(STORAGE_KEY);
      if (raw){
        const parsed = JSON.parse(raw);
        const w = Number(parsed.drawerWidth);
        state = {
          items: Array.isArray(parsed.items) ? parsed.items : [],
          minimized: parsed.minimized !== false,
          drawerWidth: (Number.isFinite(w) && w >= 320) ? w : null,
          onboardingSeen: parsed.onboardingSeen === true,
          maximised: parsed.maximised === true,
          consent: (parsed.consent && typeof parsed.consent === 'object') ? parsed.consent : {}
        };
      }
    } catch(e){ state = { items:[], minimized:true, drawerWidth:null, onboardingSeen:false, maximised:false, consent:{} }; }
    /* Purge the legacy on-disk copy whether or not it parsed: patient data
       must not persist past this load, and a corrupt blob is still a blob. */
    try { localStorage.removeItem(STORAGE_KEY); } catch(e){}
    /* A bundle saved before consent existed carries collected RCI tables and no
       consent map. Those tables were accepted by the act of collecting them, so
       grant it — otherwise a report in progress would silently stop updating. */
    state.items.forEach(i => {
      const parent = consentParent(i.sourceId);
      if (CONSENT_SOURCES.has(parent)) state.consent[parent] = true;
    });
    save();   // re-home whatever was loaded: session data to sessionStorage, prefs to localStorage
  }
  function save(){
    try { sessionStorage.setItem(STORAGE_KEY, JSON.stringify(state)); } catch(e){}
    try {
      localStorage.setItem(PREFS_KEY, JSON.stringify({
        minimized: state.minimized,
        drawerWidth: state.drawerWidth,
        onboardingSeen: state.onboardingSeen
      }));
    } catch(e){}
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
     REMOVING the cells at hidden indices.

     This used to set display:none instead, which only works for consumers that
     render CSS. The screen and the Word/HTML output honoured it, but the two
     text-based outputs read the DOM directly and cannot see it:
     exportExcel maps over tr.children, and copyAll's text/plain flavour calls
     textContent — and textContent returns the text of display:none nodes. So a
     column hidden in the working report came back in the Excel export and in
     any plain-text paste. Removing the cell is the only form of "hidden" that
     every consumer agrees on, and it also removes the dependency on Word's
     HTML importer honouring display:none, which is not something this app can
     verify.

     Column indices are counted in GRID columns, not in cells, so a cell that
     spans several columns is narrowed rather than dropped. Two rows in these
     tables need that: the "Scores"/"Interpretation" spanner above the column
     labels, and the full-width group rows that carry a test family name. Under
     the old display:none approach their colspans were left untouched and
     silently over-spanned by one for every hidden column. */
  function applyHiddenColumns(html, hiddenColumns){
    if (!hiddenColumns || !hiddenColumns.length) return html;
    const hidden = new Set(hiddenColumns);
    const tmp = document.createElement('div');
    tmp.innerHTML = html;
    tmp.querySelectorAll('table').forEach(table => {
      table.querySelectorAll('tr').forEach(tr => {
        let col = 0;
        [...tr.children].forEach(cell => {
          const span = Math.max(1, parseInt(cell.getAttribute('colspan') || '1', 10));
          let dropped = 0;
          for (let k = col; k < col + span; k++) if (hidden.has(k)) dropped++;
          if (dropped >= span) cell.remove();
          else if (dropped > 0){
            const left = span - dropped;
            if (left > 1) cell.setAttribute('colspan', String(left));
            else cell.removeAttribute('colspan');
          }
          col += span;
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
  /* hiddenOverride lets a merged block impose ONE column set on all its
     members. Without it each member applied its own, and two members that
     disagreed produced a table whose sections had different column counts —
     a Full Scale IQ row printing "100 | Average" under headings
     "Score | Percentile | Classification", i.e. a value sitting under the
     wrong heading in an exported report. See mergedHiddenColumns(). */
  function effectiveItemHtml(item, hiddenOverride){
    const hidden = hiddenOverride || item.hiddenColumns || [];
    let html = applyHiddenColumns(item.html, hidden);
    html = applyHeaderOverrides(html, item.headerOverrides || []);
    html = stripTableNumber(html);
    html = ensureBottomRule(html);
    return html;
  }
  /* The column set a merged block shows: a column is hidden only where EVERY
     member hides it. Intersection rather than union, deliberately — hiding is
     a presentation preference, so resolving a disagreement by showing the
     column costs the reader a column they did not want, while resolving it the
     other way would silently drop a member's data from the report. */
  function mergedHiddenColumns(items){
    const sets = items.map(it => new Set(it.hiddenColumns || []));
    if (!sets.length) return [];
    return [...sets[0]].filter(idx => sets.every(s => s.has(idx)));
  }

  /* Detect the test family name (CVLT-3, WAIS-IV, etc.) from a captured APA
     table's content. Used to label the "Added to report" pill with the
     actual test rather than the analysis type. */
  const TEST_FAMILY_PATTERNS = [
    /* Order matters: the first match wins, so the more specific edition must
       come before the bare abbreviation. CVLT-C sat missing here while being
       catalogued as its own battery, so a split CVLT-C item resolved to plain
       'CVLT' and never matched anything in REPORT_TEST_CATALOG. */
    'CVLT-3', 'CVLT-C', 'CVLT-II', 'CVLT',
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
    /* FIRST cell of each group row only. A group row may now carry per-column
       labels in its later cells — longest span relabels the percentile column
       to "Base rate" — and scanning every cell picked that up as the test
       family, which broke both the title and the battery merge. Only the first
       cell ever holds the section name. */
    const groupCells = [...tmp.querySelectorAll('table tbody tr.apa-group')]
      .map(tr => tr.children[0]).filter(Boolean);
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
    // those are predicted outcomes, not the test family. Validity sources:
    // same shape - RBANS / Digit Span / TOMM name the indicators' host
    // tests, so the pill would read "RBANS added to report".
    if (parentId && (parentId.startsWith('pre-') || parentId.startsWith('pvt-'))){
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
        /* A section that relabels a column for itself carries those labels on
           its group row (longest span reports a base rate, not a percentile).
           The group row is dropped from the split — its name becomes the
           item's title — so the labels have to be lifted off it here or the
           report loses them, and base rates end up under "Percentile". */
        const cells = [...row.children];
        current = {
          name: ((row.querySelector('td')?.textContent) || '').trim(),
          colLabels: cells.length > 1
            ? cells.map((c, i) => (i === 0 ? '' : c.textContent.trim()))
            : [],
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
      /* Re-emit the section's column labels ahead of its rows. The first cell
         is left blank because the section name is the item's title now. Kept
         as .apa-group so the merge continues to treat it as a divider rather
         than data, and .apa-group-labelled so the merge can find it again. */
      if (section.colLabels && section.colLabels.some(Boolean)){
        const lr = document.createElement('tr');
        /* Deliberately NOT .apa-group. This row has no section name — the name
           became the item's title — so anything hunting for a section divider
           must skip it. Giving it .apa-group made detectTestFamily read
           "Base rate" as the test family and the battery merge stopped
           matching. */
        lr.className = 'apa-col-label-row';
        /* This row is built AFTER buildStandaloneHtml has inlined every style,
           so it gets none of that treatment. Rather than restate the padding
           and alignment and risk drifting from it, each label copies them off
           the first data cell of its own column — which is exactly the column
           it labels, so the two line up by construction. */
        const refRow = section.rows[0];
        section.colLabels.forEach((label, i) => {
          const td = document.createElement('td');
          if (i > 0){
            td.className = 'apa-group-col-label';
            const ref = refRow && refRow.children[i];
            const refStyle = ref ? (ref.getAttribute('style') || '') : '';
            const pad   = /padding\s*:\s*([^;]+)/i.exec(refStyle);
            const align = /text-align\s*:\s*([^;]+)/i.exec(refStyle);
            td.setAttribute('style',
              (pad ? 'padding:' + pad[1].trim() + ';' : '')
              + (align ? 'text-align:' + align[1].trim() + ';' : '')
              + "font-family:'Times New Roman',serif;font-size:9.5pt;font-style:italic;color:#333;");
            td.textContent = label || '';
          }
          lr.appendChild(td);
        });
        cloneTbody.appendChild(lr);
      }
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
    const gone = state.items.find(i => i.id === id);
    state.items = state.items.filter(i => i.id !== id);
    /* Deleting a gated method's last table is a decision about that method, so
       it withdraws consent — otherwise the next keystroke anywhere in the
       shared row set would silently put it straight back. Only when the last
       one goes: a split table losing one family must keep collecting the rest. */
    if (gone){
      const parent = consentParent(gone.sourceId);
      if (CONSENT_SOURCES.has(parent) && !state.items.some(i => consentParent(i.sourceId) === parent)){
        delete state.consent[parent];
      }
    }
    save();
    render();
    refreshConsentControls();
  }
  function clear(){
    if (!state.items.length) return;
    if (!window.confirm('Start a new report?\n\nThis clears all collected tables - the change cannot be undone.')) return;
    state.items = [];
    state.consent = {};
    save();
    render();
    refreshConsentControls();
    if (typeof showToast === 'function') showToast('New report started');
  }
  /* Consent goes with the report it belongs to. A method accepted for the last
     patient must not go on silently collecting for the next one. */
  function clearSilent(){
    state.items = [];
    state.consent = {};
    save();
    render();
    refreshConsentControls();
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
     - EXCEPT for the four consent-gated sources above, which collect nothing
       until accepted. See CONSENT_SOURCES.
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
    const handler = debounce(() => captureSource(sourceId), 350);
    const obs = new MutationObserver(handler);
    obs.observe(container, { childList: true, subtree: true, characterData: true });
    observers.set(sourceId, obs);
  }
  /* The capture itself, factored out of the observer so the "Add to report"
     control can run exactly the same path the moment consent is given —
     otherwise an accepted table would sit out of the report until the next
     stray mutation happened to fire. */
  function captureSource(sourceId, opts){
    const o = opts || {};
    const container = document.getElementById(sourceId);
    if (!container) return;
    {
      // Bail out during a global clear so we don't immediately re-populate
      // the bundle from residual APA HTML during the clear cascade.
      if (_suppressed) return;
      // Always update the row-count tracker first, even if we're about to
      // early-return, so cleared tables don't leave stale counts behind.
      const currentRowCount = countTableRows(container);
      const rowsAdded = o.forceAnnounce === true || currentRowCount > (lastRowCount[sourceId] || 0);
      lastRowCount[sourceId] = currentRowCount;

      /* The offer control appears and disappears with the table, so it is
         refreshed on every pass — including the ones that early-return. */
      refreshConsentControls();
      if (!container.querySelector('.apa-table')) return;
      /* Not accepted yet: render nothing into the report. The control on this
         method's APA toolbar is now offering it. Sits ABOVE the auto-split
         branch so a multi-family table cannot slip past as child items. */
      if (!isSourceAccepted(sourceId)) return;
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
    }
    refreshConsentControls();
  }
  function setupAllObservers(){
    SOURCE_IDS.forEach(id => ensureObserver(id));
  }
  /* Which method is on screen decides which card shows, and switching method
     tabs changes that without touching any APA container — so the capture
     observers above never fire for it. Watch the panels' own class instead,
     the way the top-bar nav sync does: it covers the method tabs, the page nav
     and the URL-restore path without needing a hook in each. */
  function watchMethodPanels(){
    CONSENT_SOURCES.forEach(sourceId => {
      const panel = document.getElementById(sourceId)?.closest('section');
      if (!panel) return;
      new MutationObserver(() => refreshConsentControls())
        .observe(panel, { attributes: true, attributeFilter: ['class'] });
    });
    /* In the consolidated Change Analysis page the method panels live inside
       #change-analysis, so leaving the page entirely leaves their own classes
       untouched. */
    const host = document.getElementById('change-analysis');
    if (host){
      new MutationObserver(() => refreshConsentControls())
        .observe(host, { attributes: true, attributeFilter: ['class'] });
    }
  }

  /* ---------- accepting a consent-gated source ----------
     Two callers, and the difference between them is the whole design:
       - capture:true  — the clinician pressed "Add to report", so collect the
                         table now and announce it exactly as an auto-add does.
       - capture:false — the clinician entered data in this method. Consent is
                         recorded, and the render that entry is about to cause
                         reaches the observer as an ordinary auto-add. Calling
                         capture here instead would re-save and re-render the
                         drawer on every keystroke. */
  function acceptSource(sourceId, opts){
    const parent = consentParent(sourceId);
    if (!CONSENT_SOURCES.has(parent)) return;
    if (state.consent[parent] === true) return;
    state.consent[parent] = true;
    save();
    if (!opts || opts.capture !== false){
      captureSource(parent, { forceAnnounce: true });
    }
    refreshConsentControls();
  }
  function resetConsent(){
    state.consent = {};
    save();
    refreshConsentControls();
  }

  /* ---------- the offer, inline in each gated method panel ----------
     It goes IN THE PAGE, directly under that method's own table.

     Two places it must not go, both established by driving the real UI:
       - the APA toolbar. `styles.css` hides every page's inline APA panel
         outright (`.apa-wrap{display:none!important}`) because the drawer is
         the canonical view of a rendered table, so a control there is on
         screen in the markup and nowhere at all to the clinician.
       - a transient prompt near the chip. It can be missed, and with the
         "+ Add to report" buttons long gone there would then be no way to
         collect that method's table short of typing in it.

     In-flow and state-driven, so it cannot go stale, cannot be missed while
     the method is open, and needs no stacking or drawer-position handling.
     It offers and it never blocks: the calculator above it works untouched
     whether the offer is taken, declined or ignored. */
  /* ONE element, mounted on <body>, never inside the method panel.

     Mounting it in the panel and letting `display:none` on the inactive panels
     hide it was tried and is wrong twice over: `position:fixed` resolves against
     the nearest TRANSFORMED ancestor, and staggerSectionContent animates panel
     children on every page entry, so the card would anchor to a moving box; and
     a card that scrolls with a thirteen-row table is below the fold exactly when
     the clinician needs it. On <body> its position is the viewport, full stop. */
  let offerEl = null;
  function offerHost(){
    if (offerEl && offerEl.isConnected) return offerEl;
    offerEl = document.createElement('div');
    offerEl.id = 'rb-offer-host';
    document.body.appendChild(offerEl);
    return offerEl;
  }
  /* The panel, not the APA container: the container lives inside .apa-wrap,
     which is display:none on every page, so it is never itself on screen. */
  function methodPanelVisible(sourceId){
    const container = document.getElementById(sourceId);
    const panel = container ? container.closest('section') : null;
    return !!(panel && panel.getClientRects().length);
  }
  function offerableSource(){
    /* Only ever one: the method currently open. The other three are behind
       display:none, and stacking cards for pages nobody is looking at would be
       four questions about one decision. */
    for (const sourceId of CONSENT_SOURCES){
      if (!methodPanelVisible(sourceId)) continue;
      if (isSourceAccepted(sourceId)) continue;
      if (state.consent[consentParent(sourceId)] === 'declined') continue;
      const container = document.getElementById(sourceId);
      if (container && container.querySelector('.apa-table')) return sourceId;
    }
    return null;
  }
  function refreshConsentControls(){
    const host = offerHost();
    /* Two things already own this corner and both outrank the offer:

       - the open drawer, where the clinician is looking at the report itself.
       - the first-run onboarding bubble, which explains what the Working Report
         IS. Offering a table to someone who has not yet been told there is a
         report to put it in is the wrong order, and the two cards overlap.

       Both are transient, and the offer is state-driven, so it cannot go stale
       while it waits: closing the drawer or dismissing the hint brings it back. */
    const hintUp = !state.onboardingSeen && state.minimized;
    const sourceId = (state.minimized && !hintUp) ? offerableSource() : null;
    if (!sourceId){
      if (host.firstChild){
        host.innerHTML = '';
        /* Clear the marker too, or a method that goes away and comes back —
           switch tab and switch back — would be blocked by the guard below and
           never re-render. */
        delete host.dataset.rbFor;
        reflowPillStack();
      }
      return;
    }
    /* Rebuilding while it is already offering the same method would restart the
       entrance animation on every keystroke in the table. The pill reflow still
       runs: pills stack off the card's measured height, and a pill can arrive
       (or the card rewrap on resize) while the offer itself is unchanged. */
    if (host.dataset.rbFor === sourceId){ reflowPillStack(); return; }
    host.dataset.rbFor = sourceId;
    const label = escapeHtmlLocal(SOURCE_METHOD_NAMES[sourceId] || 'this method');
    host.innerHTML = `
      <div class="rb-offer-box" role="status">
        <div class="rb-offer-text">
          <strong>Not in the working report</strong>
          <span>These scores were entered on another Change Analysis tab, so the ${label} table has not been collected.</span>
        </div>
        <div class="rb-offer-actions">
          <button class="rb-offer-add" type="button" data-rb-accept="${escapeHtmlLocal(sourceId)}">Add to report</button>
          <button class="rb-offer-skip" type="button" data-rb-decline="${escapeHtmlLocal(sourceId)}">Not this one</button>
        </div>
      </div>`;
    /* Pills stack from the chip upward; with the card there they start above it. */
    reflowPillStack();
  }
  /* Declining is per patient, not permanent: it clears with the bundle, and
     entering data in the method accepts it outright. */
  function declineSource(sourceId){
    const parent = consentParent(sourceId);
    if (!CONSENT_SOURCES.has(parent)) return;
    if (state.consent[parent] === true) return;
    state.consent[parent] = 'declined';
    save();
    refreshConsentControls();
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
      const bodyRows = src
        ? [...src.querySelectorAll('tbody tr')].filter(r =>
            !r.classList.contains('apa-group') && !r.classList.contains('apa-col-label-row'))
        : [];
      const noteEl = tmp.querySelector('.apa-note');
      /* A section that relabels a column for itself — longest span reports a
         base rate, not a percentile — carries that on its group row. The line
         above drops every .apa-group row, so without this the merged table
         would print base rates under the "Percentile" heading with nothing
         saying so, which is the misreading the label exists to prevent. */
      /* Either shape: the unsplit table labels its group row, while a split
         item carries a standalone label row with no section name. */
      const labelRow = src
        ? src.querySelector('tbody tr.apa-col-label-row, tbody tr.apa-group-labelled')
        : null;
      const colLabels = labelRow
        ? [...labelRow.children].map((c, i) => (i === 0 ? '' : c.textContent.trim()))
        : [];
      return { sec, src, headRows, colRow, bodyRows, noteEl, colLabels,
               colCount: colRow ? colRow.children.length : 0 };
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
          const pad = (si > 0 ? '10pt' : '4pt') + ' 0 2pt';
          const nameStyle = "font-family:'Times New Roman',serif;font-size:11pt;font-weight:bold;color:#000;text-align:left;padding:" + pad + ";";
          if (p.colLabels.some(Boolean)){
            /* This section relabels a column, so the divider is built cell by
               cell and carries that label in its own column — same shape the
               unmerged table uses. Inline styles because this HTML is pasted
               into Word, where the stylesheet does not follow it. */
            dr.className = 'apa-group apa-group-labelled';
            for (let ci = 0; ci < colspan; ci++){
              const td = document.createElement('td');
              if (ci === 0){
                td.setAttribute('style', nameStyle);
                td.textContent = p.sec.subLabel;
              } else {
                /* Copy padding and alignment from this column's own data cells
                   so the label sits over the numbers it describes. Falls back
                   to the divider's padding when the section has no rows yet. */
                const ref = p.bodyRows[0] && p.bodyRows[0].children[ci];
                const refStyle = ref ? (ref.getAttribute('style') || '') : '';
                const rPad   = /padding\s*:\s*([^;]+)/i.exec(refStyle);
                const rAlign = /text-align\s*:\s*([^;]+)/i.exec(refStyle);
                /* font-weight:normal is explicit: this cell sits on a row
                   classed .apa-group, whose stylesheet rule bolds it as a
                   section name. It is a column label, not a heading. */
                td.setAttribute('style',
                  'padding:' + (rPad ? rPad[1].trim() : pad) + ';'
                  + (rAlign ? 'text-align:' + rAlign[1].trim() + ';' : '')
                  + "font-family:'Times New Roman',serif;font-size:9.5pt;font-style:italic;font-weight:normal;color:#333;");
                td.textContent = p.colLabels[ci] || '';
              }
              dr.appendChild(td);
            }
          } else {
            const td = document.createElement('td');
            td.setAttribute('colspan', String(colspan));
            td.setAttribute('style', nameStyle);
            td.textContent = p.sec.subLabel;
            dr.appendChild(td);
          }
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
  /* The battery an item may be merged under, or null if it must stand alone.

     Premorbid tables are refused outright. Their content names WAIS-IV and
     WMS-IV because those are the PREDICTED criteria, so detectTestFamily reads
     them as Wechsler tables and the merge would stack a predicted-score table
     into an achieved-score one — the two sets of numbers side by side in one
     column, which is the single worst thing this feature could do. Same reason
     buildIntelligentTitle and pillLabelFor skip family detection for pre-*. */
  function mergeableBattery(item, processedHtml){
    const parentId = (item.sourceId || '').split('::')[0];
    if (parentId.startsWith('pre-')) return null;
    return catalogBatteryFor(detectTestFamily(processedHtml));
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
      const info = mergeableBattery(it, processed);
      /* Key on the battery AND the parent tool. Two tables merge only when they
         are the same instrument AND the same analysis — which is exactly the
         case the split creates, one Score Tables entry arriving as "WAIS-IV
         Indices" + "WAIS-IV Core Subtests".

         Keying on the battery alone would also fuse a Score Tables results
         table with a reliable-change table for the same instrument. Those carry
         different columns and answer different questions, and stacking them
         under one title would misrepresent both. */
      const key = info ? `${(it.sourceId || '').split('::')[0]}::${info.id}` : null;
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
      /* Re-render the members against ONE column set before combining them.
         Each was processed above with its own hiddenColumns, and members that
         disagreed produced sections of different widths in the same table. */
      if (isMerged){
        const common = mergedHiddenColumns(ordered.map(m => m.item));
        ordered.forEach(m => { m.html = effectiveItemHtml(m.item, common); });
      }
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
      if (typeof exportToast === 'function') exportToast(`✓ ${blocks.length} table${blocks.length===1?'':'s'} copied`);
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
        <span class="rb-chip-count is-zero" data-rb-count>0</span>
      </button>
      <div class="rb-onboarding" data-rb-onboarding hidden aria-hidden="true">
        <div class="rb-onboarding-bubble">
          <div class="rb-onboarding-title">Your Working Report</div>
          <div class="rb-onboarding-body">
            <p>Every APA-formatted table you generate across the suite is <strong>auto-saved into one place</strong> as you work.</p>
            <p>Open the panel to <strong>reorder, hide columns, edit titles, and export</strong> the whole bundle to Word, Excel, or your clipboard - paste straight into your report.</p>
            <p>Everything lives in your browser. <strong>Nothing is uploaded, nothing leaves your device</strong> — and the report clears itself when you close the tab.</p>
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
    /* The consent control lives on the calculator page, not in the drawer, so
       it needs its own listener — the handler below returns early for anything
       outside #report-bundle-root. */
    document.addEventListener('click', e => {
      const acceptBtn = e.target.closest('[data-rb-accept]');
      if (acceptBtn){ e.preventDefault(); acceptSource(acceptBtn.dataset.rbAccept); return; }
      const declineBtn = e.target.closest('[data-rb-decline]');
      if (declineBtn){ e.preventDefault(); declineSource(declineBtn.dataset.rbDecline); }
    });
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
          // Forward to the topbar's "New patient" button so it does
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
    /* The consent card occupies the same corner, so pills start above it when
       it is showing. Measured rather than assumed: the card wraps to two or
       three lines depending on the method name and the viewport. Its rect is in
       VISUAL px (body carries zoom:0.9) while `bottom` is applied in LAYOUT px,
       so the height is divided back out — see pageZoomFactor. */
    const card = document.querySelector('#rb-offer-host .rb-offer-box');
    let base = PILL_BASE_BOTTOM;
    if (card){
      const z = (typeof pageZoomFactor === 'function' ? pageZoomFactor() : 1) || 1;
      base += Math.round(card.getBoundingClientRect().height / z) + PILL_GAP;
    }
    // Newest pill (last added to DOM) sits closest to the chip; stack upward
    [...pills].reverse().forEach((pill, i) => {
      pill.style.bottom = (base + i * (PILL_HEIGHT + PILL_GAP)) + 'px';
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
      if (typeof exportToast === 'function') exportToast('✓ Table copied to clipboard');
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
      if (typeof exportToast === 'function') exportToast('✓ Merged table copied to clipboard');
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
      el.classList.toggle('is-zero', countNow === 0);
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
    /* The card hides while the drawer is open, and open/close routes through
       here. Safe from recursion: refreshConsentControls never calls render. */
    refreshConsentControls();

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
    watchMethodPanels();
    render();
    refreshConsentControls();
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

  return { init, addOrReplace, remove, clear, clearSilent, setSuppressed, isSuppressed, copyAll, exportWord, exportExcel, open, close, toggle, showKofiPrompt: maybeShowKofiToast,
           acceptSource, isSourceAccepted, isConsentGated, resetConsent, refreshConsentControls };
})();

if (document.readyState === 'loading'){
  document.addEventListener('DOMContentLoaded', () => ReportBundle.init());
} else {
  ReportBundle.init();
}
