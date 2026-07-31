/* ============================================================================
   tools/check.js — headless numeric regression checks

   Run:  node tools/check.js
   Exit: 0 if everything passes, 1 if anything fails.

   WHY THIS EXISTS
   data.js is a plain script of globals with no DOM dependency, so it can be
   loaded in Node via the `vm` module and its tables checked directly. app.js
   is mostly DOM-bound, so the formulas that matter are RE-IMPLEMENTED here and
   run against data.js's real coefficients. That means this file is a second,
   independent statement of each formula: if app.js and this script ever
   disagree, one of them has drifted and you need to find out which.

   Keep that property. Do not import formulas from app.js — the duplication is
   the point.

   WHAT KIND OF CHECK BELONGS HERE
   Two kinds only:
     1. INVARIANTS  — things that must be true by mathematics or by the shape
                      of the data (monotonicity, positivity, round-trips).
     2. PINNED      — values verified against a published source document, with
                      the source named in the check. Never invent an expected
                      number; if you cannot cite it, derive it in the check and
                      show the arithmetic.

   ADDING A CHECK
   Call check('name', fn) or checkClose('name', actual, expected, tol, 'source').
   Group it under an existing section banner or add a new one.
   ============================================================================ */

'use strict';

const fs = require('fs');
const vm = require('vm');
const path = require('path');

// ---------------------------------------------------------------------------
// Load data.js as a script into a sandbox and pull out the globals we need.
// ---------------------------------------------------------------------------
const ROOT = path.resolve(__dirname, '..');
const sandbox = {};
vm.createContext(sandbox);
vm.runInContext(
  fs.readFileSync(path.join(ROOT, 'data.js'), 'utf8') +
    ';globalThis.__EXPORTS = { TOPF_TO_FSIQ, WAIS_COEF, WMS_COEF,' +
    ' OPIE_PRORATED_FSIQ, OPIE_PRORATED_GAI, OPIE_PRORATED_INDEX,' +
    ' BASE_RATES, OPIE_BASE_RATES, OCC_CODE, normDB,' +
    ' OPIE_AGE_MIN, OPIE_AGE_MAX, PRE_MODEL_TOOLTIPS };',
  sandbox
);
const D = sandbox.__EXPORTS;

// ---------------------------------------------------------------------------
// Source extraction. Several checks pull a function STRAIGHT OUT of the shipped
// file and execute it, rather than re-implementing it, so that they test what
// actually ships. Only safe for pure functions (no DOM, no globals) or ones
// whose dependencies are stubbed at the call site. If a target is refactored
// so it is no longer a top-level function, extraction throws rather than
// silently passing.
// ---------------------------------------------------------------------------
function extractFn(source, name) {
  const start = source.indexOf('function ' + name + '(');
  if (start === -1) throw new Error('could not find function ' + name + ' in app.js');
  let depth = 0, i = source.indexOf('{', start);
  if (i === -1) throw new Error('malformed function ' + name);
  for (let j = i; j < source.length; j++) {
    if (source[j] === '{') depth++;
    else if (source[j] === '}') { depth--; if (depth === 0) return source.slice(start, j + 1); }
  }
  throw new Error('unbalanced braces in ' + name);
}

const APP_SRC = fs.readFileSync(path.join(ROOT, 'app.js'), 'utf8');
const HTML_SRC = fs.readFileSync(path.join(ROOT, 'index.html'), 'utf8');

/* The confidence-interval block of Methods & References, delimited by a pair of
   HTML comments in index.html.

   IT USED TO BE ONE <p>. Sections 24 and 28 both asserted their contract over
   `/<strong>Confidence intervals...<\/strong>[\s\S]*?<\/p>/`, which quietly made
   the boundary of a paragraph load-bearing: naming seven instruments, their
   coefficient tables and their exclusions inside one <p> produced 1,100 words
   of unbroken prose, and splitting it for legibility would have silently
   truncated what the checks could see rather than failing. The sentinels move
   the boundary off the markup's shape, so the block can be paragraphs, a
   heading and a table without weakening either assertion.

   Throws rather than returning '' if a sentinel is missing: a check that
   silently passes over an empty string is worse than no check. */
function methodsCiBlock() {
  const open = HTML_SRC.indexOf('<!-- methods-ci');
  const close = HTML_SRC.indexOf('<!-- /methods-ci -->');
  if (open === -1 || close === -1 || close < open) {
    throw new Error('the methods-ci sentinels are missing from index.html — '
      + 'sections 24 and 28 cannot locate the confidence-interval block');
  }
  /* Start AFTER the opening comment closes, not at it. The comment explains
     what these checks assert, so any instrument name written into it would
     satisfy the very check it describes — a note about a check passing the
     check. Slicing past it means only rendered text can. */
  const afterComment = HTML_SRC.indexOf('-->', open);
  if (afterComment === -1 || afterComment > close) {
    throw new Error('the opening methods-ci sentinel is unterminated');
  }
  return HTML_SRC.slice(afterComment + 3, close);
}

/* Just the exception roster, which is the table listing each instrument whose
   interval is built on internal consistency, what coefficient it uses, and what
   it is not applied to.

   Section 24 asserts over THIS rather than over the whole block, because the
   block also carries a sentence naming every instrument whose published SEM
   table the app reproduces — which is most of them. Slicing the wider region
   would let that sentence satisfy the check on its own, so an instrument could
   lose its roster row, and with it its stated coefficient and exclusions, while
   the check stayed green. Proven by mutation: deleting the WISC-V row passes
   against the block and fails against the roster. */
function methodsCiRoster() {
  const open = HTML_SRC.indexOf('<table class="methods-table" id="methods-ci-roster"');
  if (open === -1) throw new Error('the methods-ci-roster table is missing from index.html');
  const close = HTML_SRC.indexOf('</table>', open);
  if (close === -1) throw new Error('the methods-ci-roster table is unclosed');
  return HTML_SRC.slice(open, close);
}

/* Pull a top-level `const NAME = ...;` declaration straight out of the shipped
   file, so a sandbox gets the REAL value rather than a copy that can drift.
   Same contract as extractFn: throws rather than silently passing. */
function extractConst(source, name) {
  const re = new RegExp('^const\\s+' + name + '\\s*=[^\\n]*(?:\\n(?!const |function |let |var )[^\\n]*)*?;\\s*$', 'm');
  const m = source.match(re);
  if (!m) throw new Error('could not find const ' + name + ' in app.js');
  return m[0];
}

// ---------------------------------------------------------------------------
// Tiny harness
// ---------------------------------------------------------------------------
let passed = 0;
const failures = [];
let section = '';

function heading(name) {
  section = name;
  console.log('\n' + name);
  console.log('-'.repeat(name.length));
}

function check(name, fn) {
  let ok, detail = '';
  try {
    const r = fn();
    if (r === true || r === undefined) ok = true;
    else { ok = false; detail = String(r); }
  } catch (e) {
    ok = false;
    detail = e && e.message ? e.message : String(e);
  }
  if (ok) { passed++; console.log('  PASS  ' + name); }
  else { failures.push({ section, name, detail }); console.log('  FAIL  ' + name + (detail ? '  -> ' + detail : '')); }
}

// Floating-point comparison. Default tolerance is 5e-3, loose enough for
// values quoted to 2dp in a manual, tight enough to catch a real change.
function checkClose(name, actual, expected, tol, source) {
  tol = (tol == null) ? 5e-3 : tol;
  check(name + (source ? '  [' + source + ']' : ''), () => {
    if (!Number.isFinite(actual)) return 'got ' + actual;
    const d = Math.abs(actual - expected);
    return d <= tol ? true : 'got ' + actual.toFixed(6) + ', expected ' + expected + ' (diff ' + d.toExponential(2) + ', tol ' + tol + ')';
  });
}

// ---------------------------------------------------------------------------
// Statistical primitives, re-implemented independently of app.js.
// High-accuracy normal CDF via the incomplete gamma function — deliberately a
// different algorithm from app.js's erf approximation, so agreement between
// them is meaningful rather than circular.
// ---------------------------------------------------------------------------
function lgamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}
function gammp(a, x) {
  if (x < a + 1) {                              // series expansion
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 1; n < 500; n++) {
      ap++; del *= x / ap; sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-17) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
  }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d; // continued fraction
  for (let i = 1; i < 500; i++) {
    const an = -i * (i - a);
    b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
    c = b + an / c;         if (Math.abs(c) < 1e-300) c = 1e-300;
    d = 1 / d; const del = d * c; h *= del;
    if (Math.abs(del - 1) < 1e-17) break;
  }
  return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
}
function normCDF(z) {
  // The incomplete-gamma route returns NaN at the infinities; app.js's erf
  // route handles them, so match it rather than diverging on the edge cases.
  if (z === Infinity) return 1;
  if (z === -Infinity) return 0;
  const t = Math.abs(z) / Math.SQRT2;
  const e = gammp(0.5, t * t);
  return 0.5 * (1 + (z >= 0 ? e : -e));
}

// SEM and change-score helpers. These mirror app.js but are written out here
// so a drift in either shows up as a failure.
const semOf      = (sd, r)            => sd * Math.sqrt(1 - r);
const sdiffJT    = (sd, r)            => Math.SQRT2 * semOf(sd, r);                 // Jacobson & Truax (1991)
const sdiffIv    = (sd1, sd2, r)      => Math.sqrt(semOf(sd1, r) ** 2 + semOf(sd2, r) ** 2); // Iverson (2001)
const sdDiffObs  = (sd1, sd2, r)      => Math.sqrt(sd1 * sd1 + sd2 * sd2 - 2 * r * sd1 * sd2); // observed difference SD

// OPIE-4 prediction, mirroring app.js getOpiePredictions.predict
function opiePredict(c, vc, mr, age, sexCode) {
  let p = c.intercept + (c.age != null ? c.age * age : 0) + c.age3 * Math.pow(age, 3) + c.sex * sexCode;
  if (c.age6 != null) p += c.age6 * Math.pow(age, 6);
  if (c.vc != null && vc != null) p += c.vc * vc;
  if (c.mr != null && mr != null) p += c.mr * mr;
  return p;
}

function eachNormEntry(fn) {
  for (const group in D.normDB) {
    for (const name in D.normDB[group]) {
      const e = D.normDB[group][name];
      if (e && typeof e === 'object') fn(e, group, name);
    }
  }
}

/* ==========================================================================
   1. Statistical primitives
   ========================================================================== */
heading('1. Statistical primitives');

check('normCDF(0) is exactly 0.5', () => normCDF(0) === 0.5 || 'got ' + normCDF(0));
checkClose('normCDF(-1.959964) = 0.025', normCDF(-1.959964), 0.025, 1e-6, 'standard normal');
checkClose('normCDF(-1.644854) = 0.05',  normCDF(-1.644854), 0.05,  1e-6, 'standard normal');
checkClose('normCDF(1) - normCDF(-1) = 0.682689', normCDF(1) - normCDF(-1), 0.6826895, 1e-6, 'standard normal');
check('normCDF is monotone increasing over [-5, 5]', () => {
  let prev = -Infinity;
  for (let z = -5; z <= 5; z += 0.01) {
    const v = normCDF(z);
    if (v < prev) return 'decreased at z=' + z.toFixed(2);
    prev = v;
  }
  return true;
});
check('normCDF is symmetric: F(-z) = 1 - F(z)', () => {
  for (const z of [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 4]) {
    if (Math.abs(normCDF(-z) - (1 - normCDF(z))) > 1e-12) return 'asymmetric at z=' + z;
  }
  return true;
});

/* ==========================================================================
   2. Score-type conversions (round trips)
   ========================================================================== */
heading('2. Score-type conversions');

const SCALES = { standard: { m: 100, sd: 15 }, t: { m: 50, sd: 10 }, scaled: { m: 10, sd: 3 }, z: { m: 0, sd: 1 } };
check('z -> scale -> z round-trips for every scale', () => {
  for (const key in SCALES) {
    const { m, sd } = SCALES[key];
    for (const z of [-3, -1.5, -0.4, 0, 0.7, 2, 3]) {
      const back = ((z * sd + m) - m) / sd;
      if (Math.abs(back - z) > 1e-12) return key + ' failed at z=' + z;
    }
  }
  return true;
});
check('the curve drag maps a click back to the z it was drawn at', () => {
  // drawCurve and the drag handler used to hold separate copies of the plot
  // geometry and had drifted (640/94/48 against 940/70/36), so clicking the
  // gridline at 70 entered 63.4 — out by up to 6.6 standard-score points. Both
  // now read CURVE_GEOM; this asserts the mapping is its own inverse.
  const m = APP_SRC.match(/const CURVE_GEOM = \{[^}]*\}/);
  if (!m) return 'CURVE_GEOM not found in app.js — did the shared constant go?';
  const c = {}; vm.createContext(c);
  vm.runInContext(m[0] + ';globalThis.G = CURVE_GEOM;', c);
  const G = c.G;
  const span = G.W - G.padL - G.padR;
  const xOf = (z) => G.padL + ((z - G.xMin) / (G.xMax - G.xMin)) * span;
  const cssWidth = 830;                                   // rendered width; cancels out
  const clientToZ = (z) => {
    const localX = (xOf(z) * cssWidth / G.W) * (G.W / cssWidth);   // app.js line 1075
    return G.xMin + (localX - G.padL) / span * (G.xMax - G.xMin);
  };
  for (const z of [-3, -2, -1, 0, 1, 2, 3]) {
    const err = Math.abs(clientToZ(z) - z) * 15;           // in standard-score points
    if (err > 1e-9) return 'click at z=' + z + ' returns ' + clientToZ(z).toFixed(4) +
                           ' (' + err.toFixed(2) + ' standard-score points out)';
  }
  return true;
});
check('percentile -> z -> percentile round-trips', () => {
  // invert normCDF by bisection; a drift here means the curve is wrong
  const inv = (p) => {
    let lo = -6, hi = 6;
    for (let i = 0; i < 200; i++) {
      const mid = (lo + hi) / 2;
      if (normCDF(mid) < p) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2;
  };
  for (const p of [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]) {
    if (Math.abs(normCDF(inv(p)) - p) > 1e-9) return 'failed at p=' + p;
  }
  return true;
});

/* ==========================================================================
   3. normDB structural integrity
   ========================================================================== */
heading('3. normDB structural integrity');

/* Two shapes of entry now live in normDB. A RETEST entry carries the full
   m1/sd1/m2/sd2/r set and feeds Change Analysis. A SINGLE-ADMINISTRATION entry
   carries only m1/sd1 plus a published base-rate table, because no reliability
   coefficient exists for it - WAIS-IV longest span is the first of these.

   The retest invariants below stay exactly as strict as they were; they simply
   no longer apply to entries that never claimed to have a second
   administration. eachRetestEntry keeps that distinction in one place. */
function eachRetestEntry(fn) {
  eachNormEntry((e, g, n) => { if (!e.singleAdministration) fn(e, g, n); });
}

check('116 groups and 671 entries present', () => {
  const groups = Object.keys(D.normDB).length;
  let n = 0; eachNormEntry(() => n++);
  return (groups === 116 && n === 671) || 'got ' + groups + ' groups / ' + n + ' entries';
});
check('every retest entry has m1, sd1, m2, sd2 and r', () => {
  const bad = [];
  eachRetestEntry((e, g, n) => {
    for (const f of ['m1', 'sd1', 'm2', 'sd2', 'r']) {
      if (!Number.isFinite(e[f])) bad.push(g + ' / ' + n + ' missing ' + f);
    }
  });
  return bad.length === 0 || bad.length + ' problems, first: ' + bad[0];
});
check('every single-administration entry has m1, sd1 and a published basis, and NO r', () => {
  /* singleAdministration means one thing only: no second testing, so nothing
     here may claim a retest coefficient.

     It USED to also require baseRates, which conflated "no retest data" with
     "scored by base-rate lookup". RBANS · All Ages broke that: those entries
     are scored by ordinary metric conversion and carry the reliability the
     RBANS Update tabulates by age, with no base rates anywhere. The invariant
     that actually matters is that the entry has SOME published basis for what
     it puts on screen — a base-rate table, or a reliability coefficient — and
     is not simply an empty row.

     The exception is a raw entry with neither. RBANS Line Orientation, Picture
     Naming, List Recall and List Recognition are absent from Table 3.6, so
     they derive nothing at all: no percentile, no classification, no interval.
     They are allowed through because they still let a clinician record the raw
     score, but they must be tagged raw, so this cannot become a way to smuggle
     in a standardised measure with nothing behind it. */
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (!e.singleAdministration) return;
    if (!Number.isFinite(e.m1) || !Number.isFinite(e.sd1)) bad.push(g + ' / ' + n + ' missing m1/sd1');
    const hasRates = !!(e.baseRates && Object.keys(e.baseRates).length);
    const hasRel = ['rInternal', 'rInternalByAge', 'rStability', 'rStabilityByAge']
      .some((f) => e[f] != null);
    if (!hasRates && !hasRel && e.metric !== 'raw') {
      bad.push(g + ' / ' + n + ' has neither base rates nor a published reliability, and is not raw');
    }
    /* Carrying an r here would be a claim the manual does not make, and would
       put the family back into the Change Analysis dropdowns. */
    for (const f of ['m2', 'sd2', 'r', 'rCorrected']) {
      if (e[f] != null) bad.push(g + ' / ' + n + ' unexpectedly carries ' + f);
    }
  });
  return bad.length === 0 || bad.join('; ');
});
check('all correlations lie in [0, 1)', () => {
  const bad = [];
  eachRetestEntry((e, g, n) => {
    if (!(e.r >= 0 && e.r < 1)) bad.push(g + ' / ' + n + ' r=' + e.r);
    if (Number.isFinite(e.rCorrected) && !(e.rCorrected >= 0 && e.rCorrected < 1)) bad.push(g + ' / ' + n + ' rCorrected=' + e.rCorrected);
  });
  return bad.length === 0 || bad.join('; ');
});
check('all standard deviations are strictly positive', () => {
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (!(e.sd1 > 0)) bad.push(g + ' / ' + n + ' sd1');
    if (!e.singleAdministration && !(e.sd2 > 0)) bad.push(g + ' / ' + n + ' sd2');
  });
  return bad.length === 0 || bad.join('; ');
});
check('difference variance sd1^2 + sd2^2 - 2r*sd1*sd2 is positive everywhere', () => {
  // Bounded below by (sd1-sd2)^2, so a failure means corrupt data, not maths.
  const bad = [];
  eachRetestEntry((e, g, n) => {
    const v = e.sd1 ** 2 + e.sd2 ** 2 - 2 * e.r * e.sd1 * e.sd2;
    if (!(v > 0)) bad.push(g + ' / ' + n + ' v=' + v);
  });
  return bad.length === 0 || bad.join('; ');
});
check('rCorrected present on exactly 267 of the retest entries', () => {
  let n = 0; eachNormEntry(e => { if (Number.isFinite(e.rCorrected)) n++; });
  return n === 267 || 'got ' + n;
});
check('D-KEFS and CVLT-C carry no rCorrected, but WISC-V now does', () => {
  /* WISC-V USED TO BE ON THIS LIST, AND THAT WAS AN ERROR OF FACT, not a
     property of the source. Its Technical Manual Table 4.7 prints a corrected
     r for all 34 rows; the database simply never captured it. Nothing scored
     differently — rInternal outranks rCorrected in the CI chain and reliable
     change defaults to the raw r — but the corrected-r option on the Basic RCI
     page was silently unusable for WISC-V alone.

     D-KEFS and CVLT-C are genuinely without it, so the check keeps biting for
     them, and now asserts the WISC-V position instead of assuming it. */
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (/^(D-KEFS|CVLT-C)\b/.test(g) && Number.isFinite(e.rCorrected)) bad.push(g + ' / ' + n + ' unexpectedly has one');
  });
  let wisc = 0, wiscTotal = 0;
  eachNormEntry((e, g) => {
    if (!/^WISC-V\b/.test(g)) return;
    wiscTotal++;
    if (Number.isFinite(e.rCorrected)) wisc++;
  });
  if (wisc !== wiscTotal) bad.push(wiscTotal - wisc + ' WISC-V entries lack the corrected r Table 4.7 publishes');
  return bad.length === 0 || bad.slice(0, 3).join('; ');
});
check('group keys use the U+00B7 middle dot as separator', () => {
  const bad = Object.keys(D.normDB).filter(k => k.includes(' - ') || !k.includes('·'));
  return bad.length === 0 || 'first offender: ' + bad[0];
});

/* ==========================================================================
   4. WAIS-IV normDB pinned to the Technical Manual, Table 4.5
   Source: WAIS-IV Technical and Interpretive Manual, Table 4.5 (Stability
   Coefficients). All 702 WAIS-IV fields were verified against it; these are
   spot anchors that would catch a bad edit.
   ========================================================================== */
heading('4. WAIS-IV normDB vs Technical Manual Table 4.5');

const FSIQ_ALL = D.normDB['WAIS-IV Indices · All Ages']['Full Scale IQ'];
const SRC45 = 'WAIS-IV Tech Manual Table 4.5';
checkClose('FSIQ All Ages m1 = 99.7',        FSIQ_ALL.m1,         99.7, 1e-9, SRC45);
checkClose('FSIQ All Ages sd1 = 13.8',       FSIQ_ALL.sd1,        13.8, 1e-9, SRC45);
checkClose('FSIQ All Ages m2 = 104.0',       FSIQ_ALL.m2,        104.0, 1e-9, SRC45);
checkClose('FSIQ All Ages sd2 = 15.0',       FSIQ_ALL.sd2,        15.0, 1e-9, SRC45);
checkClose('FSIQ All Ages r = 0.95',         FSIQ_ALL.r,          0.95, 1e-9, SRC45);
checkClose('FSIQ All Ages rCorrected = 0.96', FSIQ_ALL.rCorrected, 0.96, 1e-9, SRC45);

const VCI_ALL = D.normDB['WAIS-IV Indices · All Ages']['Verbal Comprehension Index'];
checkClose('VCI All Ages sd1 = 14.4', VCI_ALL.sd1, 14.4, 1e-9, SRC45);
checkClose('VCI All Ages r = 0.95',   VCI_ALL.r,   0.95, 1e-9, SRC45);

check('all five WAIS-IV age bands plus All Ages are present for Indices', () => {
  const want = ['All Ages', 'Ages 16-29', 'Ages 30-54', 'Ages 55-69', 'Ages 70-90'];
  const missing = want.filter(b => !D.normDB['WAIS-IV Indices · ' + b]);
  return missing.length === 0 || 'missing: ' + missing.join(', ');
});

check('rCorrected is the Allen & Yen correction toward the normative SD', () => {
  // r_c = (r*k)/sqrt(1 - r^2 + r^2 k^2), k = SDnorm/sd1.
  // Confirms the stored corrected values sit on the normative metric, which is
  // WHY they must not be paired with sd1 in a SEM. Tolerance is loose because
  // the published values are rounded to 2dp and the relation is sensitive.
  const bad = [];
  for (const grp of ['Indices', 'Core Subtests']) {
    const key = 'WAIS-IV ' + grp + ' · All Ages';
    const S = grp === 'Indices' ? 15 : 3;
    for (const n in D.normDB[key]) {
      const e = D.normDB[key][n];
      const k = S / e.sd1;
      const pred = (e.r * k) / Math.sqrt(1 - e.r * e.r + e.r * e.r * k * k);
      if (Math.abs(pred - e.rCorrected) > 0.04) bad.push(n + ' pred ' + pred.toFixed(3) + ' vs ' + e.rCorrected);
    }
  }
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   5. OPIE-4 coefficients pinned to Table eA5.8
   Source: Holdnack, Schoenberg, Lange & Iverson (2013), Table eA5.8.
   All 40 coefficient slots were verified against it.
   ========================================================================== */
heading('5. OPIE-4 coefficients vs Table eA5.8');

const EA58 = 'Holdnack et al. 2013 Table eA5.8';
const OPIE_PUBLISHED = {
  'FSIQ.VC_MR': { t: D.OPIE_PRORATED_FSIQ.VC_MR,  intercept: 65.77827122, vc: 0.646258435, mr: 1.182068623, age: -0.197692558, age3: 3.73292e-05, sex: 1.955504838 },
  'FSIQ.VC':    { t: D.OPIE_PRORATED_FSIQ.VC,     intercept: 86.63733022, vc: 0.825479066,                  age: -0.355783733, age3: 3.73292e-05, sex: 2.795447219 },
  'FSIQ.MR':    { t: D.OPIE_PRORATED_FSIQ.MR,     intercept: 62.02281403,                 mr: 1.719384768,                     age3: 2.75723e-05, sex: 1.509479923 },
  'GAI.VC_MR':  { t: D.OPIE_PRORATED_GAI.VC_MR,   intercept: 60.14203956, vc: 0.763136717, mr: 1.127062322, age: -0.246247784, age3: 4.16209e-05, sex: 4.708926488 },
  'GAI.VC':     { t: D.OPIE_PRORATED_GAI.VC,      intercept: 79.65445374, vc: 0.921039566,                  age: -0.378906405, age3: 3.99793e-05, sex: 5.001045997 },
  'GAI.MR':     { t: D.OPIE_PRORATED_GAI.MR,      intercept: 56.74797323,                 mr: 1.722933286,                     age3: 2.81088e-05, sex: 3.784565031 },
  'VCI':        { t: D.OPIE_PRORATED_INDEX.VCI,   intercept: 76.77718804, vc: 1.046716258,                  age: -0.600717374, age3: 8.88487e-05, sex: 4.196941127, age6: -4.42948e-11 },
  'PRI':        { t: D.OPIE_PRORATED_INDEX.PRI,   intercept: 60.63732634,                 mr: 1.858062305,                     age3: 2.87764e-05, sex: 4.162137929 },
};
for (const label in OPIE_PUBLISHED) {
  const spec = OPIE_PUBLISHED[label];
  check('OPIE ' + label + ' matches published coefficients  [' + EA58 + ']', () => {
    const bad = [];
    for (const f in spec) {
      if (f === 't') continue;
      const got = spec.t[f];
      if (!Number.isFinite(got) || Math.abs(got - spec[f]) > 1e-12) bad.push(f + ': got ' + got + ', want ' + spec[f]);
    }
    return bad.length === 0 || bad.join('; ');
  });
}
check('FSIQ VC_MR and VC share the same age3 (present in the source table)', () =>
  D.OPIE_PRORATED_FSIQ.VC_MR.age3 === D.OPIE_PRORATED_FSIQ.VC.age3 || 'they differ');
check('VCI is the only branch carrying an age6 term', () => {
  const all = { ...D.OPIE_PRORATED_FSIQ, ...D.OPIE_PRORATED_GAI, ...D.OPIE_PRORATED_INDEX };
  const withAge6 = Object.keys(all).filter(k => all[k].age6 != null);
  return (withAge6.length === 1 && all.VCI.age6 != null) || 'got ' + withAge6.join(',');
});
check('OPIE fitted age range is 16-90', () => (D.OPIE_AGE_MIN === 16 && D.OPIE_AGE_MAX === 90) || 'got ' + D.OPIE_AGE_MIN + '-' + D.OPIE_AGE_MAX);

/* ==========================================================================
   6. OPIE-4 predictions
   Pinned to values computed from the published coefficients for a worked case.
   ========================================================================== */
heading('6. OPIE-4 predictions (45yo male, Vocab 38, Matrix 16)');

const CASE = 'worked case, coefficients per ' + EA58;
const P = (c) => opiePredict(c, 38, 16, 45, 1);
checkClose('FSIQ Vocab+Matrix = 105.71', P(D.OPIE_PRORATED_FSIQ.VC_MR),  105.71, 5e-3, CASE);
checkClose('FSIQ Vocab only   = 108.19', P(D.OPIE_PRORATED_FSIQ.VC),     108.19, 5e-3, CASE);
checkClose('FSIQ Matrix only  =  93.55', P(D.OPIE_PRORATED_FSIQ.MR),      93.55, 5e-3, CASE);
checkClose('GAI  Vocab+Matrix = 104.59', P(D.OPIE_PRORATED_GAI.VC_MR),   104.59, 5e-3, CASE);
checkClose('GAI  Vocab only   = 106.25', P(D.OPIE_PRORATED_GAI.VC),      106.25, 5e-3, CASE);
checkClose('GAI  Matrix only  =  90.66', P(D.OPIE_PRORATED_GAI.MR),       90.66, 5e-3, CASE);
checkClose('VCI  Vocab        = 101.45', P(D.OPIE_PRORATED_INDEX.VCI),   101.45, 5e-3, CASE);
checkClose('PRI  Matrix       =  97.15', P(D.OPIE_PRORATED_INDEX.PRI),    97.15, 5e-3, CASE);

check('sex term is non-zero on every branch (blank sex must not be silently female)', () => {
  // Female = 0 in OPIE-4, so a zero default would return the female equation.
  // app.js gates on sex being present; this guards the coefficient side of it.
  const all = { ...D.OPIE_PRORATED_FSIQ, ...D.OPIE_PRORATED_GAI, ...D.OPIE_PRORATED_INDEX };
  const zero = Object.keys(all).filter(k => !(Math.abs(all[k].sex) > 0));
  return zero.length === 0 || 'zero sex coefficient on: ' + zero.join(',');
});

check('VCI and PRI branches self-calibrate to 100 at their own normative raw', () => {
  // Solve each single-subtest branch for the raw score giving 100, then confirm
  // it round-trips. Catches an intercept or slope edit on either branch.
  for (const age of [20, 45, 75]) {
    const cv = D.OPIE_PRORATED_INDEX.VCI, cm = D.OPIE_PRORATED_INDEX.PRI;
    const v = (100 - opiePredict(cv, 0, null, age, 0.5)) / cv.vc;
    const m = (100 - opiePredict(cm, null, 0, age, 0.5)) / cm.mr;
    if (Math.abs(opiePredict(cv, v, null, age, 0.5) - 100) > 1e-9) return 'VCI failed at age ' + age;
    if (Math.abs(opiePredict(cm, null, m, age, 0.5) - 100) > 1e-9) return 'PRI failed at age ' + age;
  }
  return true;
});

/* ==========================================================================
   7. Reliable change
   The methods here must agree with what each page prints on screen and cites
   in its APA note. If you change a formula, change the on-screen text too.
   ========================================================================== */
heading('7. Reliable change indices');

// RCI defaults to the RAW r on all four methods: the corrected r is scaled to
// the normative population and must not be paired with the retest sample's SD.
const eF = FSIQ_ALL;
checkClose('J&T 90% FSIQ threshold = 7.18 pts',
  1.645 * sdiffJT(eF.sd1, eF.r), 7.18, 5e-3, 'Jacobson & Truax 1991, raw r');
checkClose('Iverson 90% FSIQ threshold = 7.50 pts',
  1.645 * sdiffIv(eF.sd1, eF.sd2, eF.r), 7.50, 5e-3, 'Iverson 2001, raw r');
checkClose('J&T 95% FSIQ threshold = 8.55 pts',
  1.96 * sdiffJT(eF.sd1, eF.r), 8.55, 5e-3, 'Jacobson & Truax 1991, raw r');

check('pairing sd1 with rCorrected understates the error variance', () => {
  // Guards the reasoning behind the raw-r default. Both coherent pairings land
  // near 9; the mix lands near 7.6. If this ever stops holding, the data changed.
  const coherentSample = eF.sd1 ** 2 * (1 - eF.r);
  const coherentNorm   = 15 ** 2 * (1 - eF.rCorrected);
  const mixed          = eF.sd1 ** 2 * (1 - eF.rCorrected);
  if (!(mixed < coherentSample && mixed < coherentNorm)) return 'mix no longer lowest';
  if (Math.abs(coherentSample - coherentNorm) > 1.0) return 'the two coherent pairings diverged: ' + coherentSample.toFixed(2) + ' vs ' + coherentNorm.toFixed(2);
  return true;
});

check('Iverson Sdiff equals the observed difference SD when sd1 = sd2', () => {
  // Algebraic identity; the two diverge by exactly r*(sd1-sd2)^2 in variance.
  for (const [sd, r] of [[15, 0.9], [3, 0.75], [10, 0.6]]) {
    if (Math.abs(sdiffIv(sd, sd, r) - sdDiffObs(sd, sd, r)) > 1e-12) return 'failed at sd=' + sd;
  }
  return true;
});

check('Iverson Sdiff never exceeds the observed difference SD', () => {
  // (sd1^2+sd2^2)(1-r) <= sd1^2+sd2^2-2r*sd1*sd2 by AM-GM, with equality at sd1=sd2.
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (sdiffIv(e.sd1, e.sd2, e.r) > sdDiffObs(e.sd1, e.sd2, e.r) + 1e-9) bad.push(g + ' / ' + n);
  });
  return bad.length === 0 || bad.join('; ');
});

check('the outcome label reports significance only, never a direction', () => {
  // Extracted from app.js so this tests the shipped function. normDB carries
  // 119 higher-is-worse measures (intrusions, perseverations, errors, false
  // positives, repetitions); mapping the sign of the statistic onto
  // "improvement" / "decline" asserts the wrong clinical conclusion for every
  // one of them. The signed statistic is displayed separately, so direction
  // remains visible without the app interpreting it.
  const c = {};
  vm.createContext(c);
  vm.runInContext(
    ['erf', 'normCDF', 'tInv', 'rciOutcome'].map((n) => extractFn(APP_SRC, n)).join('\n') +
      ';globalThis.__O = rciOutcome;',
    c
  );
  const outcome = c.__O;
  const banned = /improve|decline|deteriorat|worse|better|gain|loss/i;
  for (const cv of [0.95, 0.90]) {
    for (const r of [-5, -3.2, -2.5, -1.96, -1.645, -1, 0, 1, 1.645, 1.96, 2.5, 3.2, 5]) {
      const o = outcome(r, cv);
      if (banned.test(o.label)) return 'directional label at rci=' + r + ': ' + o.label;
      // and the verdict must depend only on magnitude, not sign
      const mirror = outcome(-r, cv);
      if (o.label !== mirror.label) return 'asymmetric at rci=' + r + ': ' + o.label + ' vs ' + mirror.label;
    }
  }
  return true;
});
check('the outcome still separates significant from non-significant', () => {
  const c = {};
  vm.createContext(c);
  vm.runInContext(
    ['erf', 'normCDF', 'tInv', 'rciOutcome'].map((n) => extractFn(APP_SRC, n)).join('\n') +
      ';globalThis.__O = rciOutcome;',
    c
  );
  // Neutral wording must not collapse the two verdicts into one.
  const sig = c.__O(2.5, 0.95), ns = c.__O(1.0, 0.95);
  if (sig.label === ns.label) return 'both verdicts read "' + sig.label + '"';
  if (sig.cls === ns.cls) return 'both verdicts share the class ' + sig.cls;
  // and the z threshold must still bite in the right place
  if (c.__O(1.95, 0.95).label === sig.label) return '1.95 was called significant at 95%';
  if (c.__O(1.97, 0.95).label !== sig.label) return '1.97 was not called significant at 95%';
  return true;
});
check('the Crawford threshold is labelled as t, not as a fixed z', () => {
  // Crawford is the only method here that uses a t critical value, on
  // df = N - 2. Quoting "1.96" on its dropdown named a number that is never
  // applied: with the smallest norms shipped (N = 25) the real 95% threshold
  // is 2.069, so a statistic of 2.02 read as significant against the label
  // while the app correctly called it non-significant.
  const html = fs.readFileSync(path.join(ROOT, 'index.html'), 'utf8');
  const m = html.match(/<select class="rci-cv" data-target="rci-crawford"[^>]*>([\s\S]*?)<\/select>/);
  if (!m) return 'Crawford significance dropdown not found';
  if (/1\.96|1\.645/.test(m[1])) return 'dropdown still quotes a fixed z value: ' + m[1].replace(/\s+/g, ' ').trim();
  if (!/t\(N−2\)|t\(N-2\)/.test(m[1])) return 'dropdown no longer names the t distribution';
  // The SRB dropdown genuinely IS z-based, so its 1.96 must stay.
  const srb = html.match(/<select class="rci-cv" data-target="rci-srb"[^>]*>([\s\S]*?)<\/select>/);
  if (!srb || !/1\.96/.test(srb[1])) return 'the SRB dropdown should still state 1.96 — it really does use z';
  return true;
});
check('the gap the Crawford label hid is real across the shipped norms', () => {
  // Guards the reasoning: if tInv ever collapsed to the normal value this
  // would stop being true and the relabel would look unmotivated.
  const c = {};
  vm.createContext(c);
  vm.runInContext(APP_SRC.split('\n').slice(0, 105).join('\n') + ';globalThis.t = tInv;', c);
  const ns = [];
  eachNormEntry((e) => { if (Number.isFinite(e.n)) ns.push(e.n); });
  if (!ns.length) return 'no entries carry an N';
  const smallest = Math.min(...ns);
  const gap = c.t(0.975, smallest - 2) - 1.96;
  if (!(gap > 0.05)) return 'at the smallest shipped N (' + smallest + ') the gap is only ' + gap.toFixed(4);
  // and it must shrink toward zero as N grows
  if (!(c.t(0.975, 296) - 1.96 < gap)) return 'the gap does not shrink with N';
  return true;
});
/* Only retest entries reach the regression methods; the single-administration
   families are filtered out of the Change Analysis dropdowns precisely because
   they have no r to fit a slope from. */
check('McSweeney SRB slope and SEE are computable for every retest entry', () => {
  const bad = [];
  eachRetestEntry((e, g, n) => {
    const slope = e.r * (e.sd2 / e.sd1);
    const see = e.sd2 * Math.sqrt(1 - e.r * e.r);
    if (!Number.isFinite(slope) || !Number.isFinite(see) || see <= 0) bad.push(g + ' / ' + n);
  });
  return bad.length === 0 || bad.join('; ');
});

/* Crawford's residual SD and standard error of prediction were pinned by
   nothing at all, which is how the "View formula" block on that page came to
   print McSweeney's SEE — SD2 sqrt(1 - r^2), with no (N-1)/(N-2) — above a
   statistic computed with the factor. Nobody saw it, the block having been
   display:none since the five methods were consolidated onto one page, and the
   whole block is now deleted. These two checks make the arithmetic itself
   assertable, so the next description of it can be checked against something.

   Both formulas are re-derived here rather than imported, per this file's rule.
   Eq. 4: the residual SD carries (N-1)/(N-2) because r and sd2 are sample
   statistics on N-1 df while the t reference distribution wants N-2. Eq. 5:
   predicting a NEW case adds 1 (the case's own error) and 1/N + the leverage
   term (uncertainty in the fitted line). */
check('Crawford SEE and SEpred reproduce Crawford & Garthwaite (2007) Eqs. 4-5', () => {
  const c = {};
  vm.createContext(c);
  vm.runInContext(
    APP_SRC.split('\n').slice(0, 105).join('\n')
    + ';' + extractFn(APP_SRC, 'rciEffectiveR')
    + ';' + extractFn(APP_SRC, 'calcCrawfordRow')
    + ';var rciState = {};'
    + ';globalThis.__C = calcCrawfordRow;', c);
  const bad = [];
  // Drive it over the shipped norms, at a baseline deliberately off the mean so
  // the leverage term cannot be zero and pass by cancellation.
  eachRetestEntry((e, g, n) => {
    if (!Number.isFinite(e.n) || e.n < 3) return;
    const x1 = e.m1 + e.sd1, x2 = e.m2;
    const got = c.__C({ m1:e.m1, sd1:e.sd1, m2:e.m2, sd2:e.sd2, r:e.r, n:e.n, t1:x1, t2:x2 }, 'rci-crawford');
    if (!got) { bad.push(g + ' / ' + n + ': no result'); return; }
    const slope = e.r * (e.sd2 / e.sd1);
    const predicted = (e.m2 - slope * e.m1) + slope * x1;
    const see = e.sd2 * Math.sqrt((1 - e.r * e.r) * (e.n - 1) / (e.n - 2));
    const sePred = see * Math.sqrt(1 + 1 / e.n + ((x1 - e.m1) ** 2) / ((e.n - 1) * e.sd1 * e.sd1));
    const t = (x2 - predicted) / sePred;
    const off = [];
    if (Math.abs(got.see - see) > 1e-12) off.push('SEE ' + got.see + ' vs ' + see);
    if (Math.abs(got.sePred - sePred) > 1e-12) off.push('SEpred ' + got.sePred + ' vs ' + sePred);
    if (got.df !== e.n - 2) off.push('df ' + got.df + ' vs ' + (e.n - 2));
    if (Math.abs(got.rci - t) > 1e-12) off.push('t ' + got.rci + ' vs ' + t);
    if (off.length) bad.push(g + ' / ' + n + ': ' + off.join(', '));
  });
  return bad.length === 0 || bad.length + ' entries differ, first: ' + bad[0];
});

check('the (N−1)/(N−2) factor is what separates Crawford SEE from McSweeney SEE', () => {
  // Guards the reasoning above: if the factor were negligible, describing the
  // two SEEs with one formula would be harmless and this pinning unmotivated.
  // It is not — it always inflates, and most at the smallest shipped N.
  const ns = [];
  eachRetestEntry((e) => { if (Number.isFinite(e.n) && e.n >= 3) ns.push(e.n); });
  if (!ns.length) return 'no entries carry an N';
  const smallest = Math.min(...ns);
  const infl = Math.sqrt((smallest - 1) / (smallest - 2)) - 1;
  if (!(infl > 0.01)) return 'at the smallest shipped N (' + smallest + ') the factor inflates SEE by only ' + (infl * 100).toFixed(2) + '%';
  const biggest = Math.max(...ns);
  if (!(Math.sqrt((biggest - 1) / (biggest - 2)) - 1 < infl)) return 'the factor does not shrink as N grows';
  return true;
});

check('the dead per-method formula blocks are gone from the Change Analysis page', () => {
  // They were display:none from the moment the five methods were consolidated,
  // so their text could drift from the code unseen — and did, on Crawford.
  // The Score Charts disclosure is live and must stay.
  const html = fs.readFileSync(path.join(ROOT, 'index.html'), 'utf8');
  const all = [...html.matchAll(/<details class="formula-disclosure"/g)].length;
  if (all !== 1) return html.match(/View formula/) ? 'a per-method "View formula" block is back' : 'expected exactly the Score Charts disclosure, found ' + all;
  if (!/#charts[\s\S]{0,80}|id="charts"/.test(html)) return 'Score Charts section not found';
  if (/#change-analysis \.formula-disclosure/.test(html)) return 'the page-scoped CSS for the deleted blocks is back';
  return true;
});

/* ==========================================================================
   8. Base rates
   ========================================================================== */
heading('8. Base rates');

check('BASE_RATES has 40 discrepancy rows and 298 populated cells', () => {
  const rows = Object.keys(D.BASE_RATES).length;
  let cells = 0;
  for (const d in D.BASE_RATES) cells += Object.keys(D.BASE_RATES[d]).length;
  return (rows === 40 && cells === 298) || 'got ' + rows + ' rows / ' + cells + ' cells';
});

check('BASE_RATES cells all reproduce round(Phi(d / SEE), 4)', () => {
  // These are parametric normal values, NOT observed frequencies. They are
  // labelled as such throughout the UI. If this check fails, either the SEEs
  // changed or someone substituted real empirical figures — in which case the
  // labels in app.js / index.html / data.js must change too.
  const see = {};
  D.WAIS_COEF.forEach(c => { see[c.idx] = c.see; });
  D.WMS_COEF.forEach(c => { see[c.idx] = c.see; });
  const bad = [];
  for (const d in D.BASE_RATES) {
    for (const k in D.BASE_RATES[d]) {
      const model = Math.round(normCDF(Number(d) / see[k]) * 10000) / 10000;
      if (Math.abs(D.BASE_RATES[d][k] - model) > 1.01e-4) bad.push(d + '/' + k + ' table ' + D.BASE_RATES[d][k] + ' vs model ' + model);
    }
  }
  return bad.length === 0 || bad.length + ' cells differ, first: ' + bad[0];
});

check('BASE_RATES decrease monotonically as the discrepancy grows more negative', () => {
  const idxs = new Set();
  for (const d in D.BASE_RATES) Object.keys(D.BASE_RATES[d]).forEach(k => idxs.add(k));
  const bad = [];
  for (const k of idxs) {
    let prev = Infinity;
    for (let d = -1; d >= -40; d--) {
      const row = D.BASE_RATES[String(d)];
      if (!row || row[k] == null) continue;
      if (row[k] > prev) bad.push(k + ' rises at d=' + d);
      prev = row[k];
    }
  }
  return bad.length === 0 || bad.join('; ');
});

check('past the tabulated range the base rate is computed, not floored at 0.01%', () => {
  // "< 0.01%" was asserted for any discrepancy beyond the table. That is false
  // for seven of the eight columns: one row past the last entry the model gives
  // PRI 0.081%, PSI 0.060%, DMI 0.048%, VWMI 0.038%, IMI 0.033%, VCI 0.019%,
  // WMI 0.012%. Only FSIQ, with the narrowest SEE, is genuinely under 0.01%.
  const see = {};
  D.WAIS_COEF.forEach((c) => { see[c.idx] = c.see; });
  D.WMS_COEF.forEach((c) => { see[c.idx] = c.see; });
  const lastRow = {};
  for (const d in D.BASE_RATES) {
    for (const k in D.BASE_RATES[d]) {
      const n = Number(d);
      if (lastRow[k] == null || n < lastRow[k]) lastRow[k] = n;
    }
  }
  const overstated = [];
  for (const k in see) {
    const d = lastRow[k] - 1;
    const p = normCDF(d / see[k]) * 100;
    if (p >= 0.01) overstated.push(k + ' at d=' + d + ' is ' + p.toFixed(3) + '%');
  }
  // The point of the check: most columns DO exceed the old bound, so a floor
  // would be a false claim. If this ever stops holding the fix is unmotivated.
  if (overstated.length < 5) return 'expected most columns to exceed 0.01%, got ' + overstated.length;
  // And FSIQ must remain the exception where the bound is honest.
  const fsiq = normCDF((lastRow.FSIQ - 1) / see.FSIQ) * 100;
  if (fsiq >= 0.01) return 'FSIQ at d=' + (lastRow.FSIQ - 1) + ' is ' + fsiq.toFixed(4) + '%, no longer under the bound';
  return true;
});
check('the base-rate cell can never print a rounded-to-zero percentage', () => {
  // fmtPctBr rounds to 2dp, so a computed tail would show "0.00%" — asserting
  // zero, which is worse than the bound it replaced. Below 0.01% the display
  // must fall back to "< 0.01%", which is then true.
  const render = (p) => (p < 0.0001 ? '< 0.01%' : (p * 100).toFixed(2) + '%');
  for (const d of [-33, -39, -41, -50, -60, -80, -100]) {
    for (const s of [8.441, 9.277, 10.617, 12.038, 12.165, 12.368, 12.42]) {
      const out = render(normCDF(d / s));
      if (out === '0.00%') return 'd=' + d + ' SEE=' + s + ' printed 0.00%';
    }
  }
  return true;
});
check('fmt never emits a negative zero or scientific notation', () => {
  // fmt feeds score conversions, SD Index changes, RCI statistics and predicted
  // scores. It used to print "-0.00" and, below 0.0001, values like "-6.67e-5"
  // in the middle of a two-decimal column.
  const f = (() => {
    const c = {}; vm.createContext(c);
    vm.runInContext(extractFn(APP_SRC, 'fmt') + ';globalThis.__F = fmt;', c);
    return c.__F;
  })();
  const cases = [[-0.00066667, 2], [-0.000066667, 2], [0.00004, 2], [-0.004, 2], [-0.4, 0], [0, 2]];
  for (const [v, dp] of cases) {
    const out = f(v, dp);
    if (/e[-+]/i.test(out)) return 'fmt(' + v + ',' + dp + ') = ' + out + ' — scientific notation';
    if (/^-0(\.0*)?$/.test(out)) return 'fmt(' + v + ',' + dp + ') = ' + out + ' — negative zero';
  }
  // and it must still keep the sign where the value does not round away
  for (const [v, dp, want] of [[-0.006, 2, '-0.01'], [-0.05, 1, '-0.1'], [-1.42, 2, '-1.42']]) {
    if (f(v, dp) !== want) return 'fmt(' + v + ',' + dp + ') = ' + f(v, dp) + ', want ' + want;
  }
  return true;
});
check('the curve band labels are derived and sum to 100%', () => {
  // The two outer labels were hard-coded 2.14% — P(2<|Z|<3) — applied to bands
  // running to the edge of the plot, so the six summed to 99.72%.
  if (/pct:\s*'2\.14%'/.test(APP_SRC)) return 'the hard-coded 2.14% label is back';
  const edges = [-Infinity, -2, -1, 0, 1, 2, Infinity];
  let total = 0;
  for (let i = 0; i < edges.length - 1; i++) {
    total += Number(((normCDF(edges[i + 1]) - normCDF(edges[i])) * 100).toFixed(2));
  }
  return Math.abs(total - 100) < 0.005 || 'labels would sum to ' + total.toFixed(2) + '%';
});
check('OPIE_BASE_RATES values sit on an empirical count grid (unlike BASE_RATES)', () => {
  // Genuinely empirical: values fall on ~1/1020 and ~1/2040 steps. This is the
  // contrast that proves BASE_RATES is parametric.
  const vals = new Set();
  for (const d in D.OPIE_BASE_RATES) for (const k in D.OPIE_BASE_RATES[d]) vals.add(D.OPIE_BASE_RATES[d][k]);
  const smallest = Math.min(...vals);
  return (smallest > 0 && smallest < 0.002) || 'smallest value ' + smallest + ' looks wrong for a count grid';
});

/* ==========================================================================
   9. Lookup and coefficient tables
   ========================================================================== */
heading('9. Lookup and coefficient tables');

check('TOPF_TO_FSIQ has 71 entries for raw scores 0-70', () => D.TOPF_TO_FSIQ.length === 71 || 'got ' + D.TOPF_TO_FSIQ.length);
check('TOPF_TO_FSIQ is monotone non-decreasing', () => {
  for (let i = 1; i < D.TOPF_TO_FSIQ.length; i++) {
    if (D.TOPF_TO_FSIQ[i] < D.TOPF_TO_FSIQ[i - 1]) return 'drops at raw ' + i;
  }
  return true;
});
check('TOPF raw 48 maps to FSIQ 100', () => D.TOPF_TO_FSIQ[48] === 100 || 'got ' + D.TOPF_TO_FSIQ[48]);

check('ToPF + Demographics returns ~100 for an average examinee', () => {
  // Calibration sanity on a published equation. ToPF raw 48 is the FSIQ-100
  // anchor; US mean education is ~12-13 years; sex coded mid.
  const c = D.WAIS_COEF[0];
  const v = c.intercept + c.b1 * 48 + c.b2 * 48 * 48 + c.b3 * Math.pow(48, 3) + c.edu * 12.5 + c.sex * 1.5;
  return Math.abs(v - 100) < 1.5 || 'got ' + v.toFixed(2);
});

check('WAIS_COEF and WMS_COEF carry the expected indices and finite SEEs', () => {
  const wais = D.WAIS_COEF.map(c => c.idx).join(',');
  const wms = D.WMS_COEF.map(c => c.idx).join(',');
  if (wais !== 'FSIQ,VCI,PRI,WMI,PSI') return 'WAIS_COEF idx = ' + wais;
  if (wms !== 'IMI,DMI,VWMI') return 'WMS_COEF idx = ' + wms;
  const bad = [...D.WAIS_COEF, ...D.WMS_COEF].filter(c => !(c.see > 0));
  return bad.length === 0 || 'bad see on ' + bad.map(c => c.idx).join(',');
});

check('Crawford & Allan (2001) returns a plausible value for a mid-range case', () => {
  // 87.14 - 5.21*occ + 1.78*edu + 0.18*age  (app.js calcPremorbid, row 3).
  // Occupation class 3 (Skilled), 12 years education, age 45:
  //   87.14 - 15.63 + 21.36 + 8.10 = 100.97
  // Arithmetic shown so the expected value is auditable rather than asserted.
  const v = 87.14 - 5.21 * 3 + 1.78 * 12 + 0.18 * 45;
  return Math.abs(v - 100.97) < 0.005 || 'got ' + v.toFixed(2);
});

check('OCC_CODE maps the five occupational classes to 1-5', () => {
  const vals = Object.values(D.OCC_CODE).sort();
  return vals.join(',') === '1,2,3,4,5' || 'got ' + vals.join(',');
});

/* ==========================================================================
   10. Percentile display at the distribution tails
   A percentile rank lives in the OPEN interval (0, 100). Rounding the extreme
   tails to 2dp once produced "0.00" and "100.00", which broke conversion back
   to a standard score (toZ rejects v <= 0 || v >= 100) and, while the Report
   Writer existed, reached its narrative as "0th percentile" / "100th
   percentile". That page has since been removed along with reportOrdinal, but
   fmtPct still feeds the Score Converter, Score Tables and the battery
   percentile column, so the clamp still matters. Extracted from app.js.
   ========================================================================== */
heading('10. Percentile display at the tails');

// (extractFn and APP_SRC are defined near the top of this file.)
const appFns = {};
vm.createContext(appFns);
vm.runInContext(
  extractFn(APP_SRC, 'fmtPct') + '\n;globalThis.__F = { fmtPct };',
  appFns
);
const fmtPctCheck = appFns.__F.fmtPct;
const pctOf = (z) => normCDF(z) * 100;

check('fmtPct never emits an impossible rank (0.00 or 100.00)', () => {
  // Sweep the whole reachable z range, well past the Wechsler limits.
  for (let z = -6; z <= 6; z += 0.001) {
    const s = fmtPctCheck(pctOf(z));
    if (s === '0.00' || s === '100.00') return 'emitted ' + s + ' at z=' + z.toFixed(3);
    const v = parseFloat(s);
    if (!(v > 0 && v < 100)) return 'outside the open interval: ' + s + ' at z=' + z.toFixed(3);
  }
  return true;
});


check('tail percentiles still convert back to a standard score', () => {
  // toZ returns null for v <= 0 || v >= 100, so a clamped value must stay inside.
  for (const z of [-6, -4, -3.6667, 0, 3.6667, 4, 6]) {
    const v = parseFloat(fmtPctCheck(pctOf(z)));
    if (!(v > 0 && v < 100)) return 'z=' + z + ' produced ' + v + ', which toZ would reject';
  }
  return true;
});



/* ==========================================================================
   11. Effect-size calculator (app-effectsize-page.js)
   As in section 10, functions are extracted from the shipped file and
   executed, not re-implemented. All expectations are either algebraic
   identities, internal-consistency requirements (two branches of the same
   file that must agree), or round trips against the file's own output.
   ========================================================================== */
heading('11. Effect-size calculator');

const ES_SRC = fs.readFileSync(path.join(ROOT, 'app-effectsize-page.js'), 'utf8');
const esFns = {};
vm.createContext(esFns);
// ciMultiplier calls tInv, which lives in app.js and is a plain global there —
// app.js loads before this module in index.html. Supply app.js's primitives
// block (lines 1-105) so the extracted code runs as it does in the page.
vm.runInContext(
  APP_SRC.split('\n').slice(0, 105).join('\n') + '\n' +
    ['normPdf', 'normCdf', 'statToD', 'dToStat', 'descD', 'ciMultiplier', 'ciString', 'sdFrom', 'targetShares']
      .map((n) => extractFn(ES_SRC, n)).join('\n') +
    '\n;globalThis.__ES = { statToD, dToStat, descD, ciMultiplier, ciString, sdFrom, targetShares };',
  esFns
);
const ES = esFns.__ES;

check('zstat and t branches agree (both d = 2·stat/√N)', () => {
  // t ≈ z for any reasonable N, so identical inputs must give identical d.
  // The zstat branch previously returned z/√N — Rosenthal's r, half of d.
  for (const [v, N] of [[2.5, 100], [1.96, 64], [3.0, 200]]) {
    const dt = ES.statToD('t', v, N), dz = ES.statToD('zstat', v, N);
    if (Math.abs(dt - dz) > 1e-12) return 'v=' + v + ' N=' + N + ': t->' + dt + ' vs zstat->' + dz;
  }
  return true;
});
check('every aux statistic round-trips through dToStat', () => {
  for (const [type, v, a] of [['zstat', 2.5, 100], ['t', 2.5, 100], ['g', 0.8, 52], ['md', 7.5, 15]]) {
    const back = ES.dToStat(type, ES.statToD(type, v, a), a);
    if (Math.abs(back - v) > 1e-9) return type + ' round trip gave ' + back + ', want ' + v;
  }
  return true;
});
check('g conversion uses plain Hedges J, matching compute()\'s output grid', () => {
  // g = J·d with J = 1 − 3/(4N − 9); the file previously used √J here while
  // the output grid used J, an internal contradiction.
  const N = 52, J = 1 - 3 / (4 * N - 9);
  const d = ES.statToD('g', 1, N);
  return Math.abs(d - 1 / J) < 1e-9 || 'statToD(g,1,52) = ' + d + ', want 1/J = ' + (1 / J);
});
check('CI-upper input round-trips the file\'s own displayed CI', () => {
  // Feed ciString's own upper bound straight back into sdFrom and require the
  // original SD out. This was pinned to the literal 105.88, which silently
  // assumed the old fixed 1.96 multiplier and would have masked the switch to
  // a t value; a round trip stays valid whatever multiplier is in use.
  for (const [m, sd, n] of [[100, 15, 25], [100, 15, 5], [50, 10, 12], [105, 15, 60]]) {
    const upper = Number(ES.ciString(m, sd, n).split(', ')[1]);
    const back = ES.sdFrom('ciu', upper, n, m);
    // ciString prints to 2dp, so the bound carries up to half a unit in the
    // last place. Inverting multiplies that by sqrt(n)/k, which is the exact
    // amount the round trip can lose — derive the tolerance from it rather
    // than picking a number that happens to pass.
    const tol = 1.5 * 0.005 * Math.sqrt(n) / ES.ciMultiplier(n);
    if (Math.abs(back - sd) > tol) {
      return 'n=' + n + ': upper ' + upper + ' mapped back to SD ' + back.toFixed(4) +
             ', want ' + sd + ' (tol ' + tol.toFixed(4) + ')';
    }
  }
  return true;
});
check('the CI multiplier is the t value for n-1 df, not a flat 1.96', () => {
  const want = { 5: 2.776, 10: 2.262, 20: 2.093, 30: 2.045, 60: 2.001 };  // published
  for (const n in want) {
    const got = ES.ciMultiplier(Number(n));
    if (Math.abs(got - want[n]) > 5e-3) return 'n=' + n + ' gave ' + got.toFixed(4) + ', want ' + want[n];
  }
  if (Math.abs(ES.ciMultiplier(100000) - 1.96) > 1e-3) return 'does not converge to 1.96 at large n';
  return true;
});
check('CI-upper without n returns null rather than a silent SE', () =>
  ES.sdFrom('ciu', 105.88, null, 100) === null || 'returned a value with no n');
check('descD uses Cohen/Sawilowsky anchors as bin floors', () => {
  const want = [[0.005, 'Negligible'], [0.05, 'Very Small'], [0.30, 'Small'], [0.65, 'Medium'],
                [1.00, 'Large'], [1.50, 'Very Large'], [2.50, 'Huge'],
                [0.20, 'Small'], [0.50, 'Medium'], [0.80, 'Large'], [-0.65, 'Medium']];
  for (const [d, label] of want) {
    const got = ES.descD(d).label;
    if (got !== label) return 'descD(' + d + ') = ' + got + ', want ' + label;
  }
  return true;
});
check('descD agrees with the slider/classifyD bands at the shared anchors', () => {
  // The other two classifiers are coarser (they collapse <0.2 and >=1.2), but
  // within 0.2–1.2 all three must give the same word.
  const slider = (d) => { const a = Math.abs(d);
    return a >= 1.2 ? 'Very Large' : a >= 0.8 ? 'Large' : a >= 0.5 ? 'Medium' : a >= 0.2 ? 'Small' : 'Negligible'; };
  for (const d of [0.2, 0.35, 0.5, 0.65, 0.8, 1.0, 1.19]) {
    if (ES.descD(d).label !== slider(d)) return 'disagree at d=' + d + ': ' + ES.descD(d).label + ' vs ' + slider(d);
  }
  return true;
});
check('target shares are weighted by group size (screening scenario)', () => {
  // n1=40 M=70 SD=10 vs n2=400 M=100 SD=10, cut 85 (midpoint, equal densities):
  // membership must follow the 40:400 prior -> 1/11 and 10/11.
  const ts = ES.targetShares(70, 10, 40, 100, 10, 400, 85);
  if (Math.abs(ts.share1 - 1 / 11) > 1e-9) return 'share1 = ' + ts.share1 + ', want ' + (1 / 11);
  if (Math.abs(ts.share1 + ts.share2 - 1) > 1e-12) return 'shares do not sum to 1';
  return true;
});
check('likelihood ratio stays density-based (priors-free) after the share fix', () => {
  // At the midpoint the densities are equal, so the LR must be exactly 1
  // regardless of the 40:400 group sizes. Deriving it from the weighted
  // shares would turn it into posterior odds and break this.
  const ts = ES.targetShares(70, 10, 40, 100, 10, 400, 85);
  return Math.abs(ts.ratio - 1) < 1e-9 || 'ratio = ' + ts.ratio;
});
check('missing n falls back to equal weighting, flagged as such', () => {
  const ts = ES.targetShares(70, 10, null, 100, 10, 400, 85);
  return (Math.abs(ts.share1 - 0.5) < 1e-9 && ts.weightedByN === false) ||
    'share1 = ' + ts.share1 + ', weightedByN = ' + ts.weightedByN;
});

/* ==========================================================================
   12. Premorbid CI level and battery flagging mode
   Both of these are DOM-reading functions, so they are extracted from app.js
   and run against a minimal document stub. Both previously described a
   setting other than the one actually in force.
   ========================================================================== */
heading('12. Premorbid CI level and battery flagging mode');

// Minimal stub: just enough for premorbidCi / batteryPremorbidMode to read a
// control value. If either grows a real DOM dependency, this fails loudly.
function withControls(values, fnNames, body) {
  const ctx = {
    document: {
      getElementById: (id) => (id in values ? { value: values[id] } : null)
    }
  };
  vm.createContext(ctx);
  vm.runInContext(
    fnNames.map((n) => extractFn(APP_SRC, n)).join('\n') +
      '\n;globalThis.__P = { ' + fnNames.join(', ') + ' };',
    ctx
  );
  return body(ctx.__P);
}

check('premorbidCi label, multiplier and header all agree', () => {
  const cases = [
    ['0.90', '90% CI', '90%', 1.645],
    ['0.95', '95% CI', '95%', 1.96],
    // Accepts the bare-percent encoding too, so a change to the <select>
    // option values cannot resurrect the original mismatch.
    ['90',   '90% CI', '90%', 1.645],
    ['95',   '95% CI', '95%', 1.96],
  ];
  for (const [val, label, short, mult] of cases) {
    const got = withControls({ 'pre-ci': val }, ['premorbidCi'], (P) => P.premorbidCi());
    if (got.label !== label || got.short !== short || Math.abs(got.mult - mult) > 1e-12) {
      return 'value ' + JSON.stringify(val) + ' gave ' + JSON.stringify(got);
    }
  }
  return true;
});

check('premorbidCi defaults to 90% when the control is absent or unparseable', () => {
  for (const stub of [{}, { 'pre-ci': '' }, { 'pre-ci': 'abc' }]) {
    const got = withControls(stub, ['premorbidCi'], (P) => P.premorbidCi());
    if (got.label !== '90% CI' || Math.abs(got.mult - 1.645) > 1e-12) {
      return JSON.stringify(stub) + ' gave ' + JSON.stringify(got);
    }
  }
  return true;
});

check('batteryPremorbidMode falls back to sd when SEE mode has no usable SEE', () => {
  const run = (ctrl, prem) =>
    withControls({ 'bat-prem-threshold': ctrl }, ['batteryPremorbidMode'], (P) => P.batteryPremorbidMode(prem));
  const cases = [
    ['sd',  { estimate: 110, see: 9.87 }, 'sd'],
    ['see', { estimate: 110, see: 9.87 }, 'see'],
    ['see', { estimate: 110 },            'sd'],   // no SEE at all
    ['see', { estimate: 110, see: 0 },    'sd'],   // SEE present but unusable
    ['see', null,                         'sd'],
  ];
  for (const [ctrl, prem, want] of cases) {
    const got = run(ctrl, prem);
    if (got !== want) return 'control=' + ctrl + ' prem=' + JSON.stringify(prem) + ' -> ' + got + ', want ' + want;
  }
  return true;
});

check('the premorbid anchor is read unrounded, not scraped from the table', () => {
  // getPremorbidEstimateOptions used to parse the rendered cell, which
  // fmtIntOrDash had already rounded. A true estimate of 104.6 became 105, and
  // against a score of 90 that turns 0.973 SD into exactly 1.000 SD — an
  // asterisk that should not be there. It must read preState.estimateRows.
  const src = extractFn(APP_SRC, 'getPremorbidEstimateOptions');
  if (/querySelectorAll|textContent/.test(src)) {
    return 'still reads the DOM rather than preState.estimateRows';
  }
  const c = { preState: { estimateRows: [
        { name: 'ToPF Raw Score', val: 104.6, see: 9.867 },
        { name: 'No estimate yet', val: null,  see: 9.11 }
      ], ciMult: 1.645 } };
  vm.createContext(c);
  vm.runInContext(src + ';globalThis.__G = getPremorbidEstimateOptions;', c);
  const opts = c.__G();
  if (opts.length !== 1) return 'expected 1 usable option, got ' + opts.length;
  if (Math.abs(opts[0].fsiq - 104.6) > 1e-9) return 'fsiq came back as ' + opts[0].fsiq + ', want 104.6 unrounded';
  // and the tier that the rounding used to flip must now resolve correctly
  const diffSd = (opts[0].fsiq - 90) / 15;
  if (!(diffSd < 1)) return 'diffSd = ' + diffSd.toFixed(4) + ' — the rounding artefact is back';
  return true;
});

check('the battery APA note describes the flagging mode actually in force', () => {
  // APA_NOTES is an object literal; brace-match it the same way as a function.
  const start = APP_SRC.indexOf('const APA_NOTES = {');
  if (start === -1) return 'could not locate APA_NOTES in app.js';
  let depth = 0, src = null;
  for (let j = APP_SRC.indexOf('{', start); j < APP_SRC.length; j++) {
    if (APP_SRC[j] === '{') depth++;
    else if (APP_SRC[j] === '}') { depth--; if (depth === 0) { src = APP_SRC.slice(start, j + 1); break; } }
  }
  if (!src) return 'unbalanced braces in APA_NOTES';
  const c = {}; vm.createContext(c);
  vm.runInContext(src + ';globalThis.__N = APA_NOTES;', c);
  const noteFor = (mode) => c.__N.bat({ classification: 'wechsler', mixedTypes: false,
    ciLevel: 'off', premorbid: '110', premorbidMode: mode }).filter(Boolean).join(' ');
  const sd = noteFor('sd'), see = noteFor('see');
  if (!/1 SD/.test(sd) || /standard error/.test(sd)) return 'sd-mode note is wrong: ' + sd;
  if (!/standard error/.test(see) || /1 SD/.test(see)) return 'see-mode note is wrong: ' + see;
  if (sd === see) return 'the note does not vary with the mode';
  return true;
});

/* ==========================================================================
   13. Score Tables confidence intervals
   getBatteryCiHtml is extracted from app.js and run against a stubbed
   getMergedDB (returning the real normDB) and a minimal document.
   ========================================================================== */
heading('13. Score Tables confidence intervals');

/* Exposes BOTH the renderer and the reliability lookup. The renderer rounds
   its bounds for display, so a SEM recovered from the printed string carries
   that rounding — fine for the bounds checks here, useless for pinning a SEM
   against a published figure. Section 23 uses the lookup directly. */
const batteryCtx = (() => {
  const boundsSrc = APP_SRC.match(/const BATTERY_SCORE_BOUNDS = \{[^;]*;/);
  if (!boundsSrc) throw new Error('BATTERY_SCORE_BOUNDS not found in app.js');
  /* The normative SD per metric — the multiplier in the SEM. Extracted rather
     than restated here so the checks cannot drift from the shipped table. */
  const metricSdSrc = APP_SRC.match(/const BATTERY_METRIC_SD = \{[^;]*;/);
  if (!metricSdSrc) throw new Error('BATTERY_METRIC_SD not found in app.js');
  const ctx = {
    normDB: D.normDB,
    getMergedDB: () => D.normDB,
    /* #bat-type answers 'scaled' for rowScoreType; #bat-ci-basis is ABSENT on
       purpose, so batteryCiCorrectRetest has to fall back to published. The
       checks pass the flag explicitly when they want the other reading. */
    document: { getElementById: (id) => (id === 'bat-ci-basis' ? null : { value: 'scaled' }) }
  };
  vm.createContext(ctx);
  vm.runInContext(
    boundsSrc[0] + '\n' + metricSdSrc[0] + '\n' +
      (APP_SRC.match(/const PATIENT_AGE_INPUTS = \[[^\]]*\];/) || [''])[0] + '\n' +
      ['patientAge', 'bandedReliabilityForAge', 'rInternalForAge', 'rStabilityForAge',
       'derivedCorrectedR', 'batteryCiCorrectRetest', 'resolveCiReliability',
       'batteryPatientAge', 'getBatteryRowReliability', 'rowScoreType', 'getBatteryCiHtml']
        .map((n) => extractFn(APP_SRC, n)).join('\n') +
      '\n;globalThis.__B = getBatteryCiHtml;'
      + '\n;globalThis.__REL = getBatteryRowReliability;'
      + '\n;globalThis.__AGE = rInternalForAge;'
      + '\n;globalThis.__STAB = rStabilityForAge;'
      + '\n;globalThis.__CORR = batteryCiCorrectRetest;',
    ctx
  );
  return ctx;
})();
const batteryCi  = batteryCtx.__B;
const batteryRel = batteryCtx.__REL;
const rInternalForAge = batteryCtx.__AGE;
const rStabilityForAge = batteryCtx.__STAB;
const batteryCiCorrectRetest = batteryCtx.__CORR;

const SCALED_ROW   = { name: 'Block Design', group: 'WAIS-IV Core Subtests · All Ages', scoreType: 'scaled' };
const STANDARD_ROW = { name: 'Full Scale IQ', group: 'WAIS-IV Indices · All Ages',      scoreType: 'standard' };

check('scaled-score intervals never exceed the 1-19 scale limits', () => {
  // A score of 19 previously produced "17-21"; 21 is not a possible scaled score.
  for (let ss = 1; ss <= 19; ss++) {
    const out = batteryCi(ss, SCALED_ROW, '90');
    const [lo, hi] = out.split('–').map(Number);
    if (!(lo >= 1 && hi <= 19)) return 'score ' + ss + ' gave ' + out;
  }
  return true;
});
check('the scaled-score ceiling actually binds at the top of the scale', () => {
  // Guards against the cap being silently removed: at 19 the upper end must
  // be pinned to 19 rather than floating above it.
  const out = batteryCi(19, SCALED_ROW, '90');
  return out.endsWith('–19') || 'score 19 gave ' + out;
});
check('standard-metric intervals are NOT capped at 19', () => {
  // The cap is a property of the scaled-score scale only. Applying it to
  // index scores would be catastrophic, so pin that it does not happen.
  const out = batteryCi(160, STANDARD_ROW, '90');
  const hi = Number(out.split('–')[1]);
  return hi > 19 || 'FSIQ 160 gave ' + out;
});
check('the interval depends on the row NAME, not just the score', () => {
  // The reliability is looked up by subtest name, which is why the on-screen
  // cell must refresh on a rename. Two subtests in the same family with
  // different reliabilities must give different intervals at the same score.
  const mr = batteryCi(10, { ...SCALED_ROW, name: 'Matrix Reasoning' }, '90');
  const vc = batteryCi(10, { ...SCALED_ROW, name: 'Vocabulary' }, '90');
  return mr !== vc || 'Matrix Reasoning and Vocabulary both gave ' + mr;
});
check('an unknown subtest name yields no interval rather than a wrong one', () => {
  const out = batteryCi(10, { ...SCALED_ROW, name: 'Not A Real Subtest' }, '90');
  return out === '' || 'got ' + JSON.stringify(out);
});

/* ==========================================================================
   14. Reliable-change reliability cell
   The displayed coefficient must be the one the calculation uses, and the
   cell must edit that same field. Extracted from app.js and run against a
   stubbed rciState.
   ========================================================================== */
heading('14. Reliable-change reliability cell');

function rciCellFns(toggleOn) {
  const ctx = { rciState: { 'rci-basic': { useCorrectedR: toggleOn } } };
  vm.createContext(ctx);
  vm.runInContext(
    ['escapeHtml', 'escapeAttr', 'rciReliabilityField', 'rciReliabilityCell', 'rciEffectiveR']
      .map((n) => extractFn(APP_SRC, n)).join('\n') +
      '\n;globalThis.__R = { rciReliabilityField, rciReliabilityCell, rciEffectiveR };',
    ctx
  );
  return ctx.__R;
}
const cellValue = (html) => (html.match(/value="([^"]*)"/) || [])[1];
const cellField = (html) => (html.match(/data-f="([^"]*)"/) || [])[1];

check('the displayed coefficient is always the one used in the calculation', () => {
  const rows = [
    ['corrected present and differing', { r: '0.95', rCorrected: '0.96' }],
    ['no corrected value',              { r: '0.83', rCorrected: '' }],
    ['corrected field absent',          { r: '0.70' }],
  ];
  for (const toggle of [false, true]) {
    const F = rciCellFns(toggle);
    for (const [label, row] of rows) {
      const shown = parseFloat(cellValue(F.rciReliabilityCell('rci-basic', 0, row)));
      const used = F.rciEffectiveR('rci-basic', row).value;
      if (!(Math.abs(shown - used) < 1e-12)) {
        return 'toggle=' + toggle + ' ' + label + ': shows ' + shown + ' but uses ' + used;
      }
    }
  }
  return true;
});

check('the cell edits the field the calculation reads', () => {
  // The generic RCI input handler writes rciState[m].rows[i][dataset.f], so
  // data-f must name the field in force or an edit is silently discarded.
  const on = rciCellFns(true), off = rciCellFns(false);
  const withCorr = { r: '0.95', rCorrected: '0.96' };
  const without  = { r: '0.83', rCorrected: '' };
  if (cellField(on.rciReliabilityCell('rci-basic', 0, withCorr)) !== 'rCorrected') return 'toggle on, corrected present: should bind rCorrected';
  if (cellField(on.rciReliabilityCell('rci-basic', 0, without))  !== 'r')          return 'toggle on, none present: should fall back to r';
  if (cellField(off.rciReliabilityCell('rci-basic', 0, withCorr)) !== 'r')         return 'toggle off: should always bind r';
  return true;
});

check('the field decision is per row, not per column', () => {
  // Within one table some rows carry a corrected value and some do not, so a
  // column-level treatment would misdescribe the mixed case.
  const F = rciCellFns(true);
  const a = F.rciReliabilityField('rci-basic', { r: '0.9', rCorrected: '0.92' });
  const b = F.rciReliabilityField('rci-basic', { r: '0.9', rCorrected: '' });
  return (a === 'rCorrected' && b === 'r') || 'got ' + a + ' and ' + b;
});

/* ==========================================================================
   15. Documentation contracts
   Text that must stay in step with the code.
   ========================================================================== */
heading('15. Documentation contracts');

check('every OPIE tooltip warns it is illustrative only in a UK context', () => {
  const keys = Object.keys(D.PRE_MODEL_TOOLTIPS).filter(k => /^opie/i.test(k));
  const missing = keys.filter(k => !/ILLUSTRATIVE ONLY/i.test(D.PRE_MODEL_TOOLTIPS[k]));
  return missing.length === 0 || 'missing the warning: ' + missing.join(', ');
});

check('the OPIE tooltips state that sex and an age range are required', () => {
  const t = D.PRE_MODEL_TOOLTIPS.opieDefault || '';
  return (/sex/i.test(t) && /16-90|16–90/.test(t)) || 'opieDefault no longer states the input requirements';
});

check('every adjustment column on the Change Analysis overview explains itself', () => {
  /* The overview table asserts a tick pattern against four column headings.
     A reader who does not know what "regression to mean" is cannot read the
     pattern at all, so each heading carries a hover explainer. Two ways that
     silently rots: a heading gains a key with no entry (the attribute renders
     the literal "undefined"), or a column is added to ADJUST_BY_METHOD's flag
     array without a heading, which slides every tick one column left. */
  const heads = [...HTML_SRC.matchAll(/makeColHead\('([a-z]+)'/g)].map(m => m[1]);
  const block = HTML_SRC.match(/const DIMENSION_EXPLAINERS = \{[\s\S]*?\n  \};/);
  if (!block) return 'DIMENSION_EXPLAINERS is no longer declared in index.html';
  const defined = [...block[0].matchAll(/'([a-z]+)':\s*'/g)].map(m => m[1]);
  const bad = [];

  const undefinedKeys = heads.filter(k => !defined.includes(k));
  if (undefinedKeys.length) bad.push('headings with no explainer: ' + undefinedKeys.join(', '));
  const unused = defined.filter(k => !heads.includes(k));
  if (unused.length) bad.push('explainers no heading uses: ' + unused.join(', '));

  /* One heading per tick, checked against the flags the rows actually emit. */
  const flags = HTML_SRC.match(/'rci-crawford':\s*\[([01,\s]+)\]/);
  if (!flags) return 'could not locate ADJUST_BY_METHOD in index.html';
  const width = flags[1].split(',').length;
  if (heads.length !== width) {
    bad.push(heads.length + ' headings against ' + width + ' tick columns');
  }
  return bad.length === 0 || bad.join('; ');
});

check('the OPIE UK caveat is exported once and shown on screen once', () => {
  /* The OPIE tab used to carry the same warning twice: the .caution-box at the
     top of the tab, and again in the info-box mirroring the APA note directly
     below the table. Two statements of one caveat read as two caveats, and the
     second one is skimmed.

     Which copy goes is not arbitrary. The exported note is the only one that
     leaves the app, so it keeps the caveat whatever else happens; the mirror
     drops it because the page around it already says it, at greater length.
     Three things therefore have to hold together, and each fails silently on
     its own: the export keeps it, the mirror does not repeat it, and the
     caution-box that justifies the mirror's silence is still in the markup. */
  const c = {};
  vm.createContext(c);
  vm.runInContext(APP_SRC.match(/const APA_NOTES = \{[\s\S]*?\n\};/)[0] + ';globalThis.__N = APA_NOTES;', c);
  const join = (ctx) => c.__N['pre-opiepredict'](ctx).filter(s => s && s.trim()).join(' ');
  const exported = join({}), onScreen = join({ onScreen: true });
  const CAVEAT = /US reference category/i;
  const bad = [];

  if (!CAVEAT.test(exported)) bad.push('the exported note has lost the UK caveat');
  if (CAVEAT.test(onScreen)) bad.push('the on-screen mirror still repeats the caution box');
  if (!/<div class="caution-box"[^>]*>\s*<strong>Illustrative only/.test(HTML_SRC)) {
    bad.push('the caution-box is gone from the OPIE tab, so nothing on screen carries the caveat');
  }
  /* The flag has to reach the note, or the split above is inert and the
     mirror silently reverts to printing the caveat twice. */
  if (!/onScreen:\s*true/.test(extractFn(APP_SRC, 'renderStaticApaNotes'))) {
    bad.push('renderStaticApaNotes no longer marks its render as on-screen');
  }
  /* The mirror may drop a sentence the page already states; it may not say
     anything the export does not. */
  const extra = onScreen.split(/(?<=\.)\s+/).filter(s => s && !exported.includes(s));
  if (extra.length) bad.push('the mirror says something the export does not: ' + extra.join(' | '));

  return bad.length === 0 || bad.join('; ');
});

check('the CI method lives on the Methods page, and the note stays short', () => {
  /* THE DIVISION OF LABOUR BETWEEN THE TWO, which is a judgement that will not
     survive on its own. The APA note is pasted into a report and read by people
     who do not have this app; the Methods page is read by whoever wants the
     reasoning. Detail migrates back into the note one well-meant sentence at a
     time unless something objects.

     The test for what MUST stay in the note is what an exported table could be
     misread as without it. The uncorrected pairing fails that test — it is what
     those manuals themselves do, so the interval matches the published one and
     a reader checking it finds agreement. The DERIVED coefficient passes it: a
     reader cross-checking the manual will not find that value and would
     reasonably conclude the table is wrong, so the note must declare it however
     terse the rest becomes. That one sentence is the load-bearing assertion
     here; the length bound only stops the rest growing back around it. */
  const c = {};
  vm.createContext(c);
  vm.runInContext(APP_SRC.match(/const APA_NOTES = \{[\s\S]*?\n\};/)[0] + ';globalThis.__N = APA_NOTES;', c);
  const note = (ctx) => c.__N.bat({ classification: 'wechsler', ciLevel: '95', ...ctx }).filter(Boolean).join(' ');
  const bad = [];

  /* A corrected table must say the coefficients are not the published ones.
     Both halves: that they were corrected, and what that costs the reader. */
  const derived = note({ hasDerivedR: true });
  if (!/corrected to the normative sample/i.test(derived)) {
    bad.push('a corrected table no longer says its coefficients were corrected');
  }
  if (!/not the values those manuals print|differ from the published/i.test(derived)) {
    bad.push('a corrected table no longer warns that its intervals differ from the published ones');
  }
  if (/corrected/i.test(note({}))) {
    bad.push('an uncorrected table claims a correction it did not make');
  }

  /* And the note must stay a note. 900 characters is roughly twice the current
     default-state length — loose enough not to fire on ordinary rewording,
     tight enough that a migrated paragraph trips it. */
  const plain = note({ ciAge: 45 }).replace(/<[^>]+>/g, '');
  if (plain.length > 900) bad.push('the CI note is ' + plain.length + ' characters; the method belongs on the Methods page');

  /* The Methods page has to have actually received what left the note, or the
     reasoning is simply gone. These are the two claims that moved. */
  const about = HTML_SRC.slice(HTML_SRC.indexOf('Methods &amp; conventions'));
  [[/retest, uncorrected/, 'the retest-uncorrected label is not explained on the Methods page'],
   [/Table 2\.8/, 'the Methods page no longer cites the D-KEFS evidence for the uncorrected pairing'],
   [/reliability<\/strong> control|<strong>reliability<\/strong>/i, 'the Methods page does not describe the reliability control'],
   [/median error of \.003/, 'the Methods page no longer states how the correction was validated']]
    .forEach(([re, msg]) => { if (!re.test(about)) bad.push(msg); });

  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   16. Wiring — does the app actually start, and is the corrected-r toggle
   confined to the two methods where it is defensible?

   These are not maths checks. They exist because removing the Report Writer
   took the whole init tail with it (it sat immediately after
   setupReportWriter()), leaving a dozen functions defined but never called.
   Every calculator loaded with no example rows, the premorbid page never
   initialised, autofill was never wired, and the "Clear all tables" button
   had no handler. All 102 numeric checks still passed, because none of them
   ask whether the app boots. A user found it by opening the page.
   ========================================================================== */
heading('16. Wiring');

/* styles.css, for the few checks that assert a layout invariant the markup
   and the renderer both depend on. design-system.css is deliberately NOT read:
   nothing asserted here has a rule in it. */
const CSS_SRC = fs.readFileSync(path.join(ROOT, 'styles.css'), 'utf8');
const DATA_SRC = fs.readFileSync(path.join(ROOT, 'data.js'), 'utf8');

/* Bare top-level calls, i.e. at column 0 and not inside any function. Anything
   here that app.js defines but never calls is dead init. */
const INIT_CALLS = [
  'renderBattery', 'renderSdi', 'applyCalculatorPolish',
  'enhanceCalculatorWorkflow', 'enhanceApaToolbars', 'buildDescCarousels',
  'renderConverter', 'setupPreTabs', 'buildPredictTable',
  'setupPremorbidListeners', 'calcPremorbid', 'calcPredict',
  'calcOpiePredict', 'wireBatteryAutofill', 'wireSdiAutofill', 'refreshAll'
];

check('every init function app.js defines is also invoked at top level', () => {
  const dead = INIT_CALLS.filter(fn => {
    const defined = new RegExp('^function\\s+' + fn + '\\s*\\(', 'm').test(APP_SRC);
    const called  = new RegExp('^' + fn + '\\(\\);', 'm').test(APP_SRC);
    return defined && !called;
  });
  return dead.length === 0
    || 'defined but never called at top level: ' + dead.join(', ');
});

check('the calculator tables are seeded with their example rows on load', () => {
  const missing = [];
  if (!/^batteryRows = \[/m.test(APP_SRC)) missing.push('batteryRows');
  if (!/^sdiRows = \[/m.test(APP_SRC)) missing.push('sdiRows');
  if (!/^RCI_SHARED_ROWS\.push\(/m.test(APP_SRC)) missing.push('RCI_SHARED_ROWS');
  if (!/forEach\(m => renderRci\(m\)\)/.test(APP_SRC)) missing.push('initial renderRci');
  return missing.length === 0 || 'no top-level seeding for: ' + missing.join(', ');
});

check('the "Clear all tables" button has a handler bound', () => {
  const inMarkup = /id="topbar-clear-all"/.test(HTML_SRC);
  const bound = /getElementById\('topbar-clear-all'\)/.test(APP_SRC);
  if (inMarkup && !bound) return 'button is in index.html but nothing binds it';
  return true;
});

/* The Score Charts page lives in its own script (app-viz-page.js), so section
   16's INIT_CALLS cannot see it. Its wiring is four separate registrations
   that must all exist or the page is silently unreachable or stale. */
check('the Score Charts page is fully wired: section, nav, script tag, precache', () => {
  const missing = [];
  if (!/<section class="section" id="charts">/.test(HTML_SRC)) missing.push('section#charts');
  if (!/class="nav-item" data-target="charts"/.test(HTML_SRC)) missing.push('sidebar nav entry');
  if (!/data-target="charts" data-bucket="charts"/.test(HTML_SRC)) missing.push('topnav entry');
  if (!/<script src="app-viz-page\.js\?v=/.test(HTML_SRC)) missing.push('script tag');
  const SW_SRC = fs.readFileSync(path.join(ROOT, 'service-worker.js'), 'utf8');
  if (!SW_SRC.includes("'./app-viz-page.js'")) missing.push('service-worker precache');
  return missing.length === 0 || 'missing: ' + missing.join(', ');
});

/* r is an error term in Jacobson & Truax and Iverson, but a fitted regression
   slope (r x sd2/sd1) in McSweeney/SRB and Crawford & Garthwaite. Offering the
   population-corrected r on the latter two changes the predicted score and
   produces a regression line that was never fitted — so the checkbox must not
   reappear on those pages. */
check('the corrected-r toggle appears only on the Basic and Practice methods', () => {
  const targets = [...HTML_SRC.matchAll(/class="rci-use-corrected-r"\s+data-target="([^"]+)"/g)]
    .map(m => m[1]).sort();
  const expected = ['rci-basic', 'rci-practice'];
  return JSON.stringify(targets) === JSON.stringify(expected)
    || 'toggle targets are [' + targets.join(', ') + '], expected [' + expected.join(', ') + ']';
});

check('SRB and Crawford treat r as a regression slope, so cannot use corrected r', () => {
  const bad = [];
  for (const fn of ['calcSrbRow', 'calcCrawfordRow']) {
    const src = extractFn(APP_SRC, fn);
    if (!/const slope = rel \* \(sd2 \/ sd1\)/.test(src)) bad.push(fn + ': slope no longer r*sd2/sd1');
  }
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   17. Every project function that is called actually exists

   Section 16 asks whether the init functions are CALLED. It cannot tell you
   whether the call SUCCEEDS, and that is a different failure with the same
   blast radius.

   Removing the Report Writer deleted refreshReportWriterOptions() but left
   the call to it at the end of refreshAll(). refreshAll() is itself a
   top-level init statement, so it threw during boot and every top-level
   statement after it stopped running - about a third of app.js, including
   the entire Working Report bundle, the "Clear all tables" handler and the
   Score Converter's Distribution tab. All 107 checks passed throughout,
   because a source scan sees the call and is satisfied.

   This check closes that gap without a DOM: collect every name BOUND
   anywhere in the shipped scripts (function declarations, const/let/var,
   and parameter lists), then assert that every bare call to a name in the
   project's own naming families resolves to one of them.

   Bare calls only - `foo(` and not `obj.foo(` - so DOM and builtin methods
   are not candidates. The prefix list is deliberately conservative: it
   omits get/is/has/show/add/open/close, which collide with globals such as
   getComputedStyle and getSelection.
   ========================================================================== */
heading('17. Called functions exist');

const PROJECT_SCRIPTS = ['app.js', 'design-system.js', 'app-effectsize-page.js', 'app-viz-page.js', 'data.js'];

/* Comments are not call sites. Without this, a comment recording WHY a
   function was removed - which is exactly the documentation this outage
   deserves - would fail the check that the removal is complete. Line comments
   are only stripped when the `//` is not preceded by a colon, so URLs survive. */
function stripComments(src) {
  return src
    .replace(/\/\*[\s\S]*?\*\//g, ' ')
    .replace(/(^|[^:])\/\/[^\n]*/g, '$1');
}

const ALL_SRC = stripComments(
  PROJECT_SCRIPTS.map(f => fs.readFileSync(path.join(ROOT, f), 'utf8')).join('\n;\n')
);

// Names the project binds, gathered generously so a miss means "bound nowhere".
function collectBoundNames(src) {
  const bound = new Set();
  for (const m of src.matchAll(/\bfunction\s+([A-Za-z_$][\w$]*)/g)) bound.add(m[1]);
  for (const m of src.matchAll(/\b(?:const|let|var)\s+([A-Za-z_$][\w$]*)/g)) bound.add(m[1]);
  // Parameter lists of both `function (a, b)` and `(a, b) =>`, plus single
  // `x =>`. Params can shadow or supply a callable, so they count as bound.
  for (const m of src.matchAll(/\(([^()]*)\)\s*(?:=>|\{)/g)) {
    for (const p of m[1].split(',')) {
      const id = p.trim().replace(/[=.].*$/, '').replace(/^\.\.\./, '').trim();
      if (/^[A-Za-z_$][\w$]*$/.test(id)) bound.add(id);
    }
  }
  for (const m of src.matchAll(/(?:^|[^.\w$])([A-Za-z_$][\w$]*)\s*=>/g)) bound.add(m[1]);
  return bound;
}

// Function-name families used in this codebase. See CLAUDE.md, "Navigating app.js".
const PROJECT_PREFIXES =
  /^(refresh|render|calc|setup|wire|rebuild|populate|infer|apa|rci|sdi|bat|battery|opie|premorbid|update|sync|enhance|draw|fmt|load|clear|ensure|extract|split|merge|build)[A-Z_]/;

// Globals that match a project prefix but are supplied by the platform.
const PLATFORM_GLOBALS = new Set(['clearTimeout', 'clearInterval', 'clearImmediate']);

check('every project function that is called is defined somewhere', () => {
  const bound = collectBoundNames(ALL_SRC);
  const missing = new Map();
  // A bare call: not preceded by a dot, and not itself a declaration keyword.
  for (const m of ALL_SRC.matchAll(/(^|[^.\w$])([A-Za-z_$][\w$]*)\s*\(/g)) {
    const name = m[2];
    if (!PROJECT_PREFIXES.test(name)) continue;
    if (PLATFORM_GLOBALS.has(name)) continue;
    if (bound.has(name)) continue;
    if (!missing.has(name)) {
      const line = ALL_SRC.slice(0, m.index).split('\n').length;
      missing.set(name, line);
    }
  }
  return missing.size === 0
    || 'called but never defined: ' +
       [...missing].map(([n, l]) => n + '() (first call at combined line ' + l + ')').join(', ');
});

/* The specific shape of the outage: refreshAll() runs at top level, so
   anything it calls must resolve or the rest of the file dies with it. */
check('refreshAll() calls nothing that is undefined', () => {
  const bound = collectBoundNames(ALL_SRC);
  const src = stripComments(extractFn(APP_SRC, 'refreshAll'));
  const bad = [];
  for (const m of src.matchAll(/(^|[^.\w$])([A-Za-z_$][\w$]*)\s*\(/g)) {
    const name = m[2];
    if (name === 'refreshAll' || PLATFORM_GLOBALS.has(name)) continue;
    if (!bound.has(name)) bad.push(name);
  }
  return bad.length === 0 || 'refreshAll() calls undefined: ' + bad.join(', ');
});

check('no Report Writer leftovers remain in the shipped scripts', () => {
  const hits = [];
  if (/refreshReportWriterOptions|setupReportWriter|\brwAuto[A-Za-z]*/.test(ALL_SRC)) {
    for (const m of ALL_SRC.matchAll(/refreshReportWriterOptions|setupReportWriter|\brwAuto[A-Za-z]*/g)) {
      hits.push(m[0]);
    }
  }
  return hits.length === 0 || 'leftover references: ' + [...new Set(hits)].join(', ');
});

/* ==========================================================================
   18. Raw-score measures are never given a standardised metric

   normDB holds 63 measures that are not on a standardised scale. Before
   metric:'raw' existed, inferScoreTypeForSubtest guessed the score type from
   the normative mean, which has no raw category, so it picked whichever
   standard metric the raw mean sat nearest. Score Tables then converted with
   that metric's constants and printed a percentile and a Wechsler
   classification.

   Worst case, verified in the running page: RBANS List Recognition (raw,
   M 19.6, SD 0.8, out of 20). A raw 19 is z = (19-19.6)/0.8 = -0.75, the 23rd
   percentile. Read as a scaled score it is z = (19-10)/3 = +3.00, the 99.9th
   percentile, "Very Superior" — 56.3 standard-score points out, sign inverted.
   ========================================================================== */
heading('18. Raw-score measures carry no derived scores');

/* CVLT-C declares it in the family name, so those three groups are raw
   throughout. RBANS is NOT: it is a genuinely mixed family, and treating the
   whole of it as raw was an over-generalisation that made eight perfectly good
   scaled measures per age band stop producing a percentile.

   What settles it without the manual: the eight scaled measures cluster
   between 9.6 and 11.6 despite having raw maxima of 12, 16, 20, 24, 40 and 89.
   Raw scores on scales that different cannot all land near 10 — that only
   happens on a shared metric. Coding is the clearest case, raw max 89 with a
   mean of 10.8. The four listed below break the cluster and are raw. */
const RAW_FAMILIES_WHOLE = [
  'CVLT-C Subtests (Raw Scores) · Age 8',
  'CVLT-C Subtests (Raw Scores) · Age 12',
  'CVLT-C Subtests (Raw Scores) · Age 16',
];
const RBANS_RAW_SUBTESTS = ['Line Orientation', 'Picture Naming', 'List Recognition', 'List Recall'];
const RBANS_GROUPS = ['RBANS Subtests · Ages 12-19', 'RBANS Subtests · Ages 20-89',
                      'RBANS Subtests · All Ages'];

check('the CVLT-C raw-score families are tagged on every entry', () => {
  const bad = [];
  for (const fam of RAW_FAMILIES_WHOLE) {
    const g = D.normDB[fam];
    if (!g) { bad.push(fam + ' missing from normDB'); continue; }
    for (const name in g) {
      if (g[name] && typeof g[name] === 'object' && g[name].metric !== 'raw') bad.push(fam + ' / ' + name);
    }
  }
  return bad.length === 0 || bad.length + ' untagged, first: ' + bad[0];
});

check('RBANS is tagged as the mixed family it is: 4 raw per band, 8 scaled', () => {
  const bad = [];
  for (const fam of RBANS_GROUPS) {
    const g = D.normDB[fam];
    if (!g) { bad.push(fam + ' missing'); continue; }
    for (const name in g) {
      const e = g[name];
      if (!e || typeof e !== 'object') continue;
      const shouldBeRaw = RBANS_RAW_SUBTESTS.includes(name);
      const isRaw = e.metric === 'raw';
      if (shouldBeRaw && !isRaw) bad.push(fam + ' / ' + name + ' should be raw');
      if (!shouldBeRaw && isRaw) bad.push(fam + ' / ' + name + ' is tagged raw but is on a scaled metric');
    }
  }
  return bad.length === 0 || bad.join('; ');
});

/* The cluster argument, stated as an invariant so it cannot silently rot: on a
   shared standardised metric the untagged RBANS means must sit near 10, and
   they must do so regardless of how different their raw ceilings are. */
check('every untagged RBANS subtest has a mean consistent with a scaled score', () => {
  const bad = [];
  for (const fam of RBANS_GROUPS) {
    for (const name in D.normDB[fam]) {
      const e = D.normDB[fam][name];
      if (!e || typeof e !== 'object' || e.metric === 'raw') continue;
      if (Math.abs(e.m1 - 10) > 3) bad.push(fam + ' / ' + name + ' m1=' + e.m1);
      if (!(e.sd1 >= 1.5 && e.sd1 <= 4.5)) bad.push(fam + ' / ' + name + ' sd1=' + e.sd1);
    }
  }
  return bad.length === 0 || 'not scaled-like: ' + bad.join('; ');
});

/* The WAIS-IV longest-span families are raw too, but they are raw WITH a
   published base-rate table, so they score properly rather than going blank.
   Counted separately to keep the "raw and unscoreable" total honest. */
const SPAN_GROUPS = Object.keys(D.normDB).filter(g => /^WAIS-IV Longest Span/.test(g));

check('exactly 102 entries are tagged raw, and none outside the known groups', () => {
  const fams = new Set([...RAW_FAMILIES_WHOLE, ...RBANS_GROUPS, ...SPAN_GROUPS]);
  let n = 0; const stray = [];
  eachNormEntry((e, g, name) => {
    if (e.metric !== 'raw') return;
    n++;
    if (!fams.has(g)) stray.push(g + ' / ' + name);
  });
  if (stray.length) return 'tagged raw outside the known families: ' + stray.join(', ');
  /* 12 RBANS now, not 8: the · All Ages group repeats the same four raw
     subtests. They are the four absent from RBANS Update Table 3.6, which is
     the manual's own confirmation that its reliability table covers the eight
     SCALED subtests only. */
  return n === 102 || 'got ' + n + ' tagged raw, expected 102 (39 CVLT-C + 12 RBANS + 51 longest span)';
});

/* No standardised metric in this app has an SD below 1, so an entry whose sd1
   is under 1 cannot be a standardised score. This is the evidence the RBANS
   families are raw, stated as an invariant rather than taken on trust. */
check('every entry with sd1 < 1 is tagged raw', () => {
  const bad = [];
  eachNormEntry((e, g, name) => {
    if (e.sd1 < 1 && e.metric !== 'raw') bad.push(g + ' / ' + name + ' sd1=' + e.sd1);
  });
  return bad.length === 0 || bad.join('; ');
});

{
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(
    extractConst(APP_SRC, 'SCORE_METRICS') + '\n' +
    ['toZ', 'fromZ', 'normInv', 'normCDF', 'erf', 'inferScoreTypeForSubtest', 'inferScoreType',
     'scoreTypeLabel', 'scoreTypeAbbr', 'sdiSdUnit']
      .map(n => extractFn(APP_SRC, n)).join('\n') +
    ';globalThis.__E={toZ,fromZ,inferScoreTypeForSubtest,scoreTypeLabel,scoreTypeAbbr,sdiSdUnit};',
    ctx
  );
  const A = ctx.__E;

  /* Two fields, two questions. metric says what m1/sd1 are; reportedAs says
     what the clinician types in. Where a measure has both — CVLT-C, whose
     norms here are raw but whose record form gives a T or z score —
     reportedAs must win, because the entry pages are asking the second
     question. Where there is no reportedAs, metric answers both. */
  check('a raw-normed entry resolves to its REPORTED metric, else to raw', () => {
    const bad = [];
    eachNormEntry((e, g, name) => {
      if (e.metric !== 'raw') return;
      const t = A.inferScoreTypeForSubtest(g, name, e);
      const want = e.reportedAs || 'raw';
      if (t !== want) bad.push(g + ' / ' + name + ' -> ' + t + ', expected ' + want);
    });
    return bad.length === 0 || bad.length + ' wrong, first: ' + bad[0];
  });

  check('reportedAs takes precedence over metric, not the other way round', () => {
    const probe = { m1: 42, sd1: 9, m2: 48, sd2: 9, r: 0.7, metric: 'raw', reportedAs: 't' };
    const got = A.inferScoreTypeForSubtest('Some Family', 'Some Measure', probe);
    if (got !== 't') return 'got ' + got + '; metric is overriding reportedAs';
    // and with no reportedAs it must fall back to raw
    const bare = { m1: 42, sd1: 9, m2: 48, sd2: 9, r: 0.7, metric: 'raw' };
    const got2 = A.inferScoreTypeForSubtest('Some Family', 'Some Measure', bare);
    return got2 === 'raw' || 'without reportedAs it gave ' + got2 + ', expected raw';
  });

  check('toZ and fromZ refuse the raw metric in both directions', () => {
    const problems = [];
    for (const v of [0, 1, 5, 19, 42, 92.76, -3]) {
      if (A.toZ(v, 'raw') !== null) problems.push('toZ(' + v + ') = ' + A.toZ(v, 'raw'));
    }
    for (const z of [-3, -0.75, 0, 1, 3]) {
      if (A.fromZ(z, 'raw') !== null) problems.push('fromZ(' + z + ') = ' + A.fromZ(z, 'raw'));
    }
    return problems.length === 0 || problems.join('; ');
  });

  /* The specific number this closes. Pinned to normDB rather than typed in,
     so it moves if the data does. */
  check('RBANS List Recognition raw 19 yields no percentile, not the 99.9th', () => {
    const e = D.normDB['RBANS Subtests · Ages 20-89']['List Recognition'];
    const t = A.inferScoreTypeForSubtest('RBANS Subtests · Ages 20-89', 'List Recognition', e);
    if (t !== 'raw') return 'typed as ' + t + ', not raw';
    if (A.toZ(19, t) !== null) return 'toZ returned ' + A.toZ(19, t);
    // and show the arithmetic the tag prevents
    const wrongZ = (19 - 10) / 3;                  // read as a scaled score
    const trueZ  = (19 - e.m1) / e.sd1;            // against its own norms
    const gap = Math.abs(15 * (wrongZ - trueZ));
    return gap > 50 || 'expected the averted error to exceed 50 SS points, got ' + gap.toFixed(1);
  });

  check('sdiSdUnit has no divisor for raw, so index mode cannot score it', () => {
    return A.sdiSdUnit('raw') === undefined || 'got ' + A.sdiSdUnit('raw');
  });

  check('the raw metric has its own visible label, not the "Score" fallback', () => {
    if (A.scoreTypeLabel('raw') !== 'Raw Score') return 'label is ' + A.scoreTypeLabel('raw');
    if (A.scoreTypeAbbr('raw') !== 'Raw') return 'abbr is ' + A.scoreTypeAbbr('raw');
    return true;
  });
}

/* ---- CVLT-C: what the clinician types is not what m1/sd1 are ----
   Source: CVLT-C Manual, Table A.1 (T-score equivalents for List A Total
   Trials 1-5, by age group) and Table A.2 (equivalents for every other index,
   on a scale running -5.0 to +5.0).

   Table A.2 is headed "Standard Score Equivalents", but a scale of -5.0 to
   +5.0 in 0.5 steps is a z-score, NOT this app's 'standard' (M 100, SD 15).
   A CVLT-C "standard score" of -1.0 is the 16th percentile; read as M 100 /
   SD 15 it would be z = -6.7, off the bottom of the scale. reportedAs must
   therefore hold 'z', never 'standard', for these. */
check('CVLT-C List A Total is a T-score, every other index a z-score', () => {
  const bad = [];
  eachNormEntry((e, g, name) => {
    if (!/^CVLT-C/.test(g)) return;
    const want = name === 'List A Trials 1-5 Total' ? 't' : 'z';
    if (e.reportedAs !== want) bad.push(g + ' / ' + name + ' reportedAs=' + e.reportedAs + ', expected ' + want);
  });
  return bad.length === 0 || bad.length + ' wrong, first: ' + bad[0];
});

check('no CVLT-C index is typed \'standard\' — the manual\'s wording is a trap', () => {
  const bad = [];
  eachNormEntry((e, g, name) => {
    if (/^CVLT-C/.test(g) && e.reportedAs === 'standard') bad.push(g + ' / ' + name);
  });
  return bad.length === 0
    || 'read the manual\'s "Standard Score" as M 100 / SD 15: ' + bad.join(', ');
});

check('reportedAs never conflicts with what m1/sd1 actually are', () => {
  const bad = [];
  eachNormEntry((e, g, name) => {
    if (!e.reportedAs) return;
    if (!e.metric) bad.push(g + ' / ' + name + ' declares reportedAs but not what m1/sd1 are');
  });
  return bad.length === 0 || bad.join('; ');
});

/* CVLT-C Manual Table A.2: Perseverations z -1.0 maps to 0-3 and z +5.0 to
   45 or more, so a HIGHER score means MORE errors. Recognition Hits and
   Discriminability run the normal way and must NOT be flagged. */
check('the four higher-is-worse CVLT-C measures are flagged, and only those', () => {
  const WORSE = ['Perseverations', 'Free-Recall Intrusions', 'Cued-Recall Intrusions', 'False Positives'];
  const bad = [];
  eachNormEntry((e, g, name) => {
    if (!/^CVLT-C/.test(g)) return;
    const want = WORSE.includes(name);
    if (want && !e.higherIsWorse) bad.push(g + ' / ' + name + ' not flagged');
    if (!want && e.higherIsWorse) bad.push(g + ' / ' + name + ' wrongly flagged');
  });
  return bad.length === 0 || bad.join('; ');
});

/* The classification describes PERFORMANCE, so on an error measure it is
   computed from the reflected score. The percentile is deliberately left
   unreflected — 98th is what z +2.0 gives and what a reader checking against
   the manual will calculate. */
check('a higher-is-worse row classifies on the REFLECTED score', () => {
  const src = extractFn(APP_SRC, 'batteryClassificationDetails');
  if (!/r\.higherIsWorse/.test(src)) return 'the classification cell ignores higherIsWorse';
  if (!/fromZ\(-z, 'standard'\)/.test(src)) return 'the score is no longer reflected for the label';
  // the descriptor must be built from the reflected value, not the raw one
  if (!/wechslerDesc\(perfSs\)/.test(src) || !/aanDesc\(perfSs\)/.test(src)) {
    return 'the descriptor is built from the unreflected score, so an error measure reads as a merit';
  }
  return true;
});

check('the percentile column is NOT reflected for error measures', () => {
  // the percentile renderer must know nothing about direction
  const bat = extractFn(APP_SRC, 'renderBattery');
  const apa = extractFn(APP_SRC, 'renderBatteryApa');
  const bad = [];
  if (/higherIsWorse/.test(bat.match(/const pct = [^\n]*/)?.[0] || '')) bad.push('renderBattery');
  if (/higherIsWorse[^\n]*normCDF|normCDF[^\n]*higherIsWorse/.test(apa)) bad.push('renderBatteryApa');
  return bad.length === 0
    || 'the percentile is being reflected in ' + bad.join(', ') + '; it must stay reproducible from z';
});

check('no premorbid asterisks are asserted on an error measure', () => {
  const src = extractFn(APP_SRC, 'batteryClassificationDetails');
  return /r\.higherIsWorse \? '' : batteryPremorbidStars/.test(src)
    || 'error measures are being compared against a premorbid estimate, which has no meaning for an error count';
});

check('the APA note explains the high-percentile / low-classification pairing', () => {
  const notes = APP_SRC.slice(APP_SRC.indexOf('const APA_NOTES'), APP_SRC.indexOf('const APA_NOTES') + 4000);
  if (!/hasHigherIsWorse/.test(notes)) return 'the note never mentions error measures';
  if (!/classification/i.test(notes.match(/hasHigherIsWorse[\s\S]{0,400}/)?.[0] || '')) {
    return 'the note does not say how to read the classification on these rows';
  }
  return true;
});

check('autofill carries higherIsWorse onto the row', () => {
  const src = extractFn(APP_SRC, 'loadFamilyIntoBattery');
  return /higherIsWorse:/.test(src)
    || 'loadFamilyIntoBattery drops the flag, so the classification cell can never see it';
});

/* The editable table and the APA export must agree about what the score
   column holds. The header used to come from the #bat-type dropdown, which
   describes the default for MANUAL rows, so an autofilled CVLT-C table was
   headed "Scaled Score" over T-scores and z-scores while the export said
   "Score". */
check('the score column header describes the rows, not the default dropdown', () => {
  const src = extractFn(APP_SRC, 'updateBatteryScoreHeader');
  if (!/batteryRows/.test(src)) return 'the header still ignores the rows';
  if (!/types\.size === 1/.test(src)) return 'the header no longer collapses a uniform table to one label';
  if (!/'Score'/.test(src)) return 'a mixed table no longer falls back to the generic "Score"';
  return true;
});

check('screen and APA export derive that header the same way', () => {
  const apa = extractFn(APP_SRC, 'renderBatteryApa');
  const head = extractFn(APP_SRC, 'updateBatteryScoreHeader');
  const rule = /types\.size === 1[\s\S]{0,90}scoreTypeLabel/;
  return (rule.test(apa) && rule.test(head))
    || 'the two headers use different rules, so the table and its export can disagree';
});

check('SDI raw mode does not prefill a raw SD for a standardised entry metric', () => {
  const src = extractFn(APP_SRC, 'loadFamilyIntoSdi');
  return /entryMetricDiffers/.test(src)
    || 'raw mode prefills sd1 regardless, so a T-score row divides by a raw-unit SD';
});

check('sdiComputeChange rejects a row whose metric has no SD unit', () => {
  const src = extractFn(APP_SRC, 'sdiComputeChange');
  return /Number\.isFinite\(unit\)/.test(src)
    || 'sdiComputeChange no longer guards a missing SD unit — a raw row would divide by undefined';
});

/* ==========================================================================
   19. Imported custom tests are validated like typed ones

   The "Add subtest" form checks that m1/sd1/m2/sd2/r are all present and
   numeric, that r is in [0,1), that the SDs are positive and that N (if given)
   is a whole number >= 3. Import applied NONE of that and merged the file
   straight into the clinical database, so a negative sd1 - which passes
   Number.isFinite - produced a confidence interval with its lower bound above
   its upper bound, and an r of 1.5 made sqrt(1-r) NaN.

   ctValidateEntry holds those rules once.

   THE TYPED FORM IS GONE. The "Add a test family" and "Add subtests" UI was
   removed once the page became a norms browser, so the four checks that drove
   ctEntryRowHtml, ctReadRows, CT_METRIC_OPTIONS and the add-subtest handler
   went with it — they tested markup that no longer exists.

   ctValidateEntry itself stays and matters MORE than before, not less: import
   is now the only way custom data enters the clinical database, so it is the
   sole gatekeeper. Everything below still pins it.
   ========================================================================== */
heading('19. Custom-test entry validation');

{
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(
    extractConst(APP_SRC, 'SCORE_METRICS') + '\n' +
    extractFn(APP_SRC, 'ctValidateEntry') + ';globalThis.__V = ctValidateEntry;', ctx);
  const V = ctx.__V;
  const good = { m1: 100, sd1: 15, m2: 103, sd2: 15, r: 0.9, n: 100 };

  check('a well-formed entry is accepted and returned numeric-typed', () => {
    const v = V(good);
    if (!v.ok) return 'rejected: ' + v.reason;
    const e = v.entry;
    const bad = ['m1','sd1','m2','sd2','r','n'].filter(k => typeof e[k] !== 'number');
    return bad.length === 0 || 'not coerced to number: ' + bad.join(', ');
  });

  check('numeric strings are accepted and coerced, as the form does', () => {
    const v = V({ m1:'100', sd1:'15', m2:'103', sd2:'15', r:'0.9' });
    return (v.ok && v.entry.sd1 === 15 && v.entry.r === 0.9) || 'got ' + JSON.stringify(v);
  });

  check('the reliability bounds match the form: r must be in [0, 1)', () => {
    const cases = [[-0.1, false], [0, true], [0.5, true], [0.999, true], [1, false], [1.5, false]];
    const bad = cases.filter(([r, want]) => V({ ...good, r }).ok !== want).map(([r]) => 'r=' + r);
    return bad.length === 0 || 'wrong verdict for ' + bad.join(', ');
  });

  check('negative and zero SDs are rejected, in both sd1 and sd2', () => {
    const bad = [];
    for (const sd of [-15, 0]) {
      if (V({ ...good, sd1: sd }).ok) bad.push('sd1=' + sd);
      if (V({ ...good, sd2: sd }).ok) bad.push('sd2=' + sd);
    }
    return bad.length === 0 || 'accepted ' + bad.join(', ');
  });

  check('N must be a whole number >= 3 when present, and may be omitted', () => {
    const bad = [];
    if (!V({ m1:100, sd1:15, m2:103, sd2:15, r:0.9 }).ok) bad.push('omitted N rejected');
    for (const [n, want] of [[2, false], [3, true], [100, true], [10.5, false], ['abc', false]]) {
      if (V({ ...good, n }).ok !== want) bad.push('n=' + n);
    }
    return bad.length === 0 || bad.join(', ');
  });

  check('a missing required field is rejected, not defaulted', () => {
    const bad = [];
    for (const k of ['m1', 'sd1', 'm2', 'sd2', 'r']) {
      const e = { ...good }; delete e[k];
      if (V(e).ok) bad.push('accepted without ' + k);
    }
    return bad.length === 0 || bad.join(', ');
  });

  check('rCorrected, when supplied, is held to the same [0, 1) bound as r', () => {
    const bad = [];
    if (!V({ ...good, rCorrected: 0.93 }).ok) bad.push('valid rCorrected rejected');
    if (V({ ...good, rCorrected: 1.2 }).ok) bad.push('rCorrected 1.2 accepted');
    if (V({ ...good, rCorrected: -0.1 }).ok) bad.push('rCorrected -0.1 accepted');
    if (V({ ...good, rCorrected: '' }).entry?.rCorrected !== undefined) bad.push('empty rCorrected became a value');
    return bad.length === 0 || bad.join(', ');
  });

  /* An imported file must not smuggle in a metric the app cannot resolve, and
     every metric the UI can produce must survive an export/import round trip. */
  check('every metric the score-type column offers survives a round trip', () => {
    const bad = [];
    for (const m of ['z', 't', 'scaled', 'standard', 'raw']) {
      if (V({ ...good, metric: m }).entry?.metric !== m) bad.push(m + ' was dropped');
    }
    return bad.length === 0 || bad.join(', ');
  });

  check('an unresolvable metric is refused, not stored and ignored', () => {
    const bad = [];
    for (const m of ['percentile', 'Scaled', 'sten', 'raw ', 0, true]) {
      if (V({ ...good, metric: m }).ok) bad.push(JSON.stringify(m));
    }
    if (V(good).entry.metric !== undefined) bad.push('invented a metric when none was given');
    if (V({ ...good, metric: '' }).entry.metric !== undefined) bad.push('empty string became a declaration');
    return bad.length === 0 || 'accepted ' + bad.join(', ');
  });

  check('non-objects are rejected rather than spread into the database', () => {
    const bad = [null, undefined, 42, 'x', [1, 2], []].filter(v => V(v).ok);
    return bad.length === 0 || 'accepted: ' + JSON.stringify(bad);
  });
}

check('the import handler validates every entry before saving', () => {
  const start = APP_SRC.indexOf("getElementById('ct-import')");
  if (start === -1) return 'import handler not found';
  const src = APP_SRC.slice(start, start + 2600);
  if (!/ctValidateEntry\(/.test(src)) return 'import no longer calls ctValidateEntry';
  if (!/Array\.isArray\(data\)/.test(src)) return 'import no longer rejects a top-level array';
  return true;
});

/* The bug where a successful import reported failure: refreshAll() sat inside
   the try, so anything it threw was caught by the "Invalid JSON file" handler. */
check('the import try block covers parsing only, not refreshAll()', () => {
  const start = APP_SRC.indexOf("getElementById('ct-import')");
  const src = APP_SRC.slice(start, start + 2600);
  const tryStart = src.indexOf('try {');
  const tryEnd = src.indexOf('catch', tryStart);
  if (tryStart === -1 || tryEnd === -1) return 'could not locate the try/catch';
  const guarded = src.slice(tryStart, tryEnd);
  return !/refreshAll\s*\(/.test(guarded)
    || 'refreshAll() is inside the try again - a successful import will report "Invalid JSON file"';
});

/* ---- the Norms Database score-type column ----
   The last route into the raw-score bug: a user-created test with raw norms
   was guessed from the mean, because the entry form had nowhere to say what
   the score type was. A custom "Recognition total" (M 19.6, SD 0.8) at 19 -
   genuinely the 23rd percentile - printed as the 99.9th, "Very Superior". */

check('a declared metric overrides the mean-based guess', () => {
  const src = extractFn(APP_SRC, 'inferScoreTypeForSubtest');
  if (!/SCORE_METRICS\.has\(stats\.metric\)/.test(src)) return 'declared metric is no longer honoured';
  // and it must come BEFORE the mean heuristic, or the guess wins
  const declaredAt = src.indexOf('SCORE_METRICS.has(stats.metric)');
  const meanAt = src.indexOf("return 'z'");
  return (declaredAt !== -1 && declaredAt < meanAt)
    || 'the declared metric is checked after the mean heuristic, so the guess wins';
});

check('getCustom survives malformed localStorage instead of throwing', () => {
  const src = extractFn(APP_SRC, 'getCustom');
  if (!/try\s*\{/.test(src)) return 'getCustom no longer guards JSON.parse';
  if (!/Array\.isArray/.test(src)) return 'getCustom no longer rejects a parsed array';
  return true;
});

/* ==========================================================================
   20. Working Report: a hidden column is hidden in EVERY output

   applyHiddenColumns used to set display:none. That only works for consumers
   that render CSS. The screen and the Word/HTML output honoured it; the two
   text-based outputs could not see it, because they read the DOM directly:
   exportExcel maps over tr.children, and copyAll's text/plain flavour uses
   textContent, which returns the text of display:none nodes. A column hidden
   in the working report therefore reappeared in the Excel export and in any
   plain-text paste — a silent disagreement between the table on screen and
   the table that reaches the report.

   These are structural: applyHiddenColumns is DOM-bound, so it is verified in
   the browser. What is asserted here is that the mechanism cannot regress to a
   CSS-only one, and that every export still routes through it.
   ========================================================================== */
heading('20. Working Report hidden columns');

check('hidden columns are REMOVED, not just styled out of sight', () => {
  const src = extractFn(APP_SRC, 'applyHiddenColumns');
  if (/display\s*:\s*none/i.test(src)){
    return 'applyHiddenColumns sets display:none again — the Excel export and '
         + 'plain-text paste cannot see that, so the column comes back';
  }
  if (!/\.remove\(\)/.test(src)) return 'applyHiddenColumns no longer removes any cell';
  return true;
});

check('cells spanning several columns are narrowed, not dropped wholesale', () => {
  const src = extractFn(APP_SRC, 'applyHiddenColumns');
  if (!/colspan/i.test(src)) return 'applyHiddenColumns ignores colspan';
  if (!/setAttribute\('colspan'/.test(src)) return 'a partly-hidden spanner is never narrowed';
  if (!/dropped >= span/.test(src)) return 'a fully-hidden spanner is not removed';
  return true;
});

/* Every output must go through effectiveItemHtml, which is what applies the
   filter. A new export that read item.html directly would bypass it.
   Named rather than counted: a count tells you the number moved, not whether
   an OUTPUT lost its filtering, which is the only thing this check is for. */
check('every output path renders through effectiveItemHtml', () => {
  const missing = ['computeBlocks', 'copyItem', 'itemCardHtml'].filter(fn => {
    let src;
    try { src = extractFn(APP_SRC, fn); } catch (e) { return true; }
    return !/effectiveItemHtml\(/.test(src);
  });
  return missing.length === 0
    || 'these render item HTML without the hidden-column filter: ' + missing.join(', ');
});

check('the CSV export builds its rows from the filtered block, not the raw item', () => {
  const start = APP_SRC.indexOf('function exportExcel(');
  if (start === -1) return 'exportExcel not found';
  const src = APP_SRC.slice(start, start + 1600);
  if (!/mergeReportBlocks\(/.test(src)) return 'exportExcel no longer uses mergeReportBlocks';
  if (/item\.html/.test(src)) return 'exportExcel reads item.html directly, bypassing the column filter';
  return true;
});

check('the plain-text clipboard flavour uses the same blocks as the HTML one', () => {
  const start = APP_SRC.indexOf('async function copyAll(');
  if (start === -1) return 'copyAll not found';
  const src = APP_SRC.slice(start, start + 1400);
  if (!/mergeReportBlocks\(/.test(src)) return 'copyAll no longer uses mergeReportBlocks';
  if (/item\.html/.test(src)) return 'copyAll reads item.html directly, bypassing the column filter';
  return true;
});

/* ==========================================================================
   21. Working Report: Merge by battery

   The toolbar button existed and did nothing for want of REPORT_TEST_CATALOG,
   which was referenced but had never been defined. The catalog now exists, so
   these guard the two ways it can silently stop working again: a battery whose
   `name` no longer matches what detectTestFamily() produces, and a normDB
   family that no catalog entry claims.
   ========================================================================== */
heading('21. Merge by battery');

const CATALOG = (() => {
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(extractConst(
    fs.readFileSync(path.join(ROOT, 'data.js'), 'utf8'), 'REPORT_TEST_CATALOG')
    + ';globalThis.__C = REPORT_TEST_CATALOG;', ctx);
  return ctx.__C;
})();

check('REPORT_TEST_CATALOG exists and every entry is complete', () => {
  if (!Array.isArray(CATALOG) || !CATALOG.length) return 'catalog missing or empty';
  const bad = [];
  CATALOG.forEach(t => {
    for (const f of ['id', 'name', 'longName']) {
      if (!t[f] || typeof t[f] !== 'string') bad.push((t.id || '?') + ' missing ' + f);
    }
    if (!Array.isArray(t.families) || !t.families.length) bad.push((t.id || '?') + ' has no families');
  });
  return bad.length === 0 || bad.join('; ');
});

check('catalog ids are unique, so two batteries cannot collide', () => {
  const ids = CATALOG.map(t => t.id);
  const dupes = ids.filter((id, i) => ids.indexOf(id) !== i);
  return dupes.length === 0 || 'duplicate ids: ' + [...new Set(dupes)].join(', ');
});

/* detectTestFamily() can only ever return a TEST_FAMILY_PATTERNS entry, so a
   catalog name outside that list is unreachable — the battery would simply
   never merge, silently, which is the exact failure this feature just had. */
check('every catalog battery name is one detectTestFamily can actually produce', () => {
  const m = APP_SRC.match(/const TEST_FAMILY_PATTERNS = \[([\s\S]*?)\];/);
  if (!m) return 'TEST_FAMILY_PATTERNS not found in app.js';
  const patterns = [...m[1].matchAll(/'([^']+)'/g)].map(x => x[1]);
  const unreachable = CATALOG.filter(t => !patterns.includes(t.name)).map(t => t.name);
  return unreachable.length === 0
    || 'not in TEST_FAMILY_PATTERNS, so never detected: ' + unreachable.join(', ');
});

check('every normDB family is claimed by exactly one battery', () => {
  const bad = [];
  Object.keys(D.normDB).forEach(group => {
    const hits = CATALOG.filter(t =>
      t.families.some(f => group === f || group.startsWith(f + ' ')) ||
      group.startsWith(t.name + ' '));
    if (hits.length === 0) bad.push(group + ' -> no battery');
    else if (hits.length > 1) bad.push(group + ' -> ' + hits.map(h => h.id).join(' AND '));
  });
  return bad.length === 0 || bad.length + ' problems, first: ' + bad[0];
});

/* Premorbid tables name WAIS-IV / WMS-IV as the PREDICTED criterion, so
   family detection reads them as Wechsler tables. Merging one into an
   achieved-score table would put predicted and obtained scores in one column. */
check('premorbid tables are refused a battery, so they can never merge in', () => {
  const src = extractFn(APP_SRC, 'mergeableBattery');
  return /startsWith\('pre-'\)/.test(src) && /return null/.test(src)
    || 'mergeableBattery no longer refuses pre-* sources';
});

check('the merge key includes the parent tool, not the battery alone', () => {
  const src = extractFn(APP_SRC, 'computeBlocks');
  return /sourceId \|\| ''\)\.split\('::'\)\[0\]/.test(src)
    || 'blocks are keyed on the battery alone; a results table could fuse with '
     + 'a reliable-change table for the same instrument';
});

{
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(extractFn(APP_SRC, 'mergedHiddenColumns') + ';globalThis.__M = mergedHiddenColumns;', ctx);
  const M = ctx.__M;

  /* Members that disagreed about a hidden column produced a table whose
     sections had different column counts — a Full Scale IQ row printing
     "100 | Average" under "Score | Percentile | Classification". */
  check('a merged block hides a column only where every member hides it', () => {
    const cases = [
      [[{ hiddenColumns: [2] }, { hiddenColumns: [] }],      [],     'one hides, one does not'],
      [[{ hiddenColumns: [2] }, { hiddenColumns: [2] }],     [2],    'both hide the same'],
      [[{ hiddenColumns: [1, 2] }, { hiddenColumns: [2, 3] }], [2],  'partial overlap'],
      [[{ hiddenColumns: [] }, { hiddenColumns: [] }],       [],     'neither hides'],
      [[{ hiddenColumns: [1] }, {}],                          [],    'member with no list at all'],
    ];
    const bad = cases.filter(([items, want]) =>
      JSON.stringify(M(items).sort()) !== JSON.stringify(want.sort()));
    return bad.length === 0
      || 'wrong for: ' + bad.map(c => c[2] + ' (got ' + JSON.stringify(M(c[0])) + ')').join('; ');
  });

  check('merged members are re-rendered against that common column set', () => {
    const src = extractFn(APP_SRC, 'computeBlocks');
    if (!/mergedHiddenColumns\(/.test(src)) return 'computeBlocks never computes a common set';
    return /effectiveItemHtml\(m\.item, common\)/.test(src)
      || 'the common set is computed but the members are not re-rendered with it';
  });
}

/* ==========================================================================
   22. WAIS-IV longest span — published base rates

   Source: WAIS-IV Administration and Scoring Manual (GB), Tables C.4 and C.5,
   cumulative percentages of the normative sample obtaining each raw longest
   span, by age group.

   These are the first entries scored by LOOKUP rather than by converting a
   metric, and the first with no retest data at all. The direction of the
   tabulated value is the thing most likely to be got wrong by a later editor,
   so it is re-derived here from the data rather than asserted.
   ========================================================================== */
heading('22. WAIS-IV longest span base rates');

const SPAN_MEASURES = ['Longest Digit Span Forward', 'Longest Digit Span Backward',
                       'Longest Digit Span Sequence', 'Longest Letter-Number Sequence'];

check('14 age bands, 51 measure-bands, Letter-Number stopping at 65-69', () => {
  if (SPAN_GROUPS.length !== 14) return 'got ' + SPAN_GROUPS.length + ' bands, expected 14';
  let n = 0; const lns = [];
  SPAN_GROUPS.forEach(g => {
    Object.keys(D.normDB[g]).forEach(m => {
      n++;
      if (!SPAN_MEASURES.includes(m)) lns.push(g + ' / ' + m + ' is not a longest-span measure');
      if (m === 'Longest Letter-Number Sequence') lns.push(g);
    });
  });
  if (n !== 51) return 'got ' + n + ' measure-bands, expected 51';
  /* LLNS is published for 16-17 through 65-69 only — 9 bands. Any more would
     mean invented data; any fewer, a transcription drop. */
  const lnsBands = SPAN_GROUPS.filter(g => D.normDB[g]['Longest Letter-Number Sequence']);
  return lnsBands.length === 9
    || 'Letter-Number Sequence present in ' + lnsBands.length + ' bands, expected 9';
});

/* THE DIRECTION CHECK. baseRates[x] is the percentage scoring x OR HIGHER.
   For a non-negative whole-number score, E[X] = the sum of P(X >= x) over
   x >= 1, so reconstructing the mean that way must reproduce the printed m1.
   If a later editor flips the table to "or lower", this fails loudly. */
check('base rates run "or higher": reconstructed means match the published m1', () => {
  const bad = [];
  let worst = 0;
  SPAN_GROUPS.forEach(g => {
    Object.entries(D.normDB[g]).forEach(([name, e]) => {
      const spans = Object.keys(e.baseRates).map(Number).sort((a, b) => a - b);
      const top = Math.max(...spans);
      let est = 0;
      for (let x = 1; x <= top; x++) {
        let p = e.baseRates[x];
        if (p == null) {                       // the table has no row for 1
          const above = spans.filter(s => s >= x);
          p = above.length ? e.baseRates[Math.min(...above)] : 0;
        }
        est += p / 100;
      }
      const diff = Math.abs(est - e.m1);
      if (diff > worst) worst = diff;
      if (diff > 0.08) bad.push(g + ' / ' + name + ' m1=' + e.m1 + ' reconstructed=' + est.toFixed(3));
    });
  });
  return bad.length === 0
    || bad.length + ' mismatched (worst ' + worst.toFixed(3) + '), first: ' + bad[0];
});

check('base rates are monotone, bounded, and reach 100% at the bottom', () => {
  const bad = [];
  SPAN_GROUPS.forEach(g => {
    Object.entries(D.normDB[g]).forEach(([name, e]) => {
      const spans = Object.keys(e.baseRates).map(Number).sort((a, b) => a - b);
      spans.forEach(s => {
        const v = e.baseRates[s];
        if (!(v >= 0 && v <= 100)) bad.push(g + ' / ' + name + ' ' + s + ' -> ' + v);
      });
      for (let i = 1; i < spans.length; i++) {
        if (e.baseRates[spans[i]] > e.baseRates[spans[i - 1]] + 1e-9) {
          bad.push(g + ' / ' + name + ' rises from ' + spans[i - 1] + ' to ' + spans[i]);
        }
      }
      // everyone scores at least the lowest tabulated span
      if (e.baseRates[spans[0]] !== 100) bad.push(g + ' / ' + name + ' bottom = ' + e.baseRates[spans[0]] + '%, expected 100');
    });
  });
  return bad.length === 0 || bad.length + ' problems, first: ' + bad[0];
});

check('the published median sits where the cumulative percentage crosses 50%', () => {
  const bad = [];
  SPAN_GROUPS.forEach(g => {
    Object.entries(D.normDB[g]).forEach(([name, e]) => {
      if (!Number.isFinite(e.median)) return;
      const spans = Object.keys(e.baseRates).map(Number);
      const atOrAbove50 = spans.filter(s => e.baseRates[s] >= 50);
      if (!atOrAbove50.length) return;
      const crossing = Math.max(...atOrAbove50);
      if (Math.abs(crossing - e.median) > 0.5) {
        bad.push(g + ' / ' + name + ' median ' + e.median + ' vs crossing ' + crossing);
      }
    });
  });
  return bad.length === 0 || bad.join('; ');
});

{
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(
    ['baseRateAtOrAbove', 'baseRatePercentile'].map(n => extractFn(APP_SRC, n)).join('\n') +
    ';globalThis.__B = { baseRateAtOrAbove, baseRatePercentile };', ctx);
  const B = ctx.__B;
  const LDSF = D.normDB['WAIS-IV Longest Span (Process) · All Ages']['Longest Digit Span Forward'];

  /* Worked from the published figures: P(>=6) = 77.3 and P(>=7) = 53.9, so
     P(<6) = 22.7 and P(=6) = 23.4, giving 22.7 + 11.7 = 34.4. */
  checkClose('LDSF All Ages, span 6 -> 34.4th percentile',
    B.baseRatePercentile(LDSF, 6), 34.4, 0.05, 'WAIS-IV Manual Table C.4');

  check('the percentile converts "or higher" into "percentage below"', () => {
    // must be the complement direction: a higher span => a higher percentile
    const p5 = B.baseRatePercentile(LDSF, 5), p8 = B.baseRatePercentile(LDSF, 8);
    return p8 > p5 || 'span 8 gives ' + p8 + ' but span 5 gives ' + p5 + ' — direction inverted';
  });

  check('spans the table skips inherit the next tabulated value', () => {
    // there is no row for 1; P(X >= 1) must still be 100
    return B.baseRateAtOrAbove(LDSF, 1) === 100 || 'got ' + B.baseRateAtOrAbove(LDSF, 1);
  });

  check('a span above the top of the table yields 0% at or above', () => {
    return B.baseRateAtOrAbove(LDSF, 10) === 0 || 'got ' + B.baseRateAtOrAbove(LDSF, 10);
  });

  check('non-whole and non-numeric spans are refused rather than guessed', () => {
    const bad = [6.5, 'abc', '', null].filter(v => B.baseRatePercentile(LDSF, v) !== null);
    return bad.length === 0 || 'accepted ' + JSON.stringify(bad);
  });
}

check('a base-rate row is scored by lookup, not by its metric', () => {
  const src = extractFn(APP_SRC, 'batteryRowPercentile');
  if (!/batteryBaseRateEntry\(/.test(src)) return 'batteryRowPercentile ignores base rates';
  const brAt = src.indexOf('batteryBaseRateEntry');
  const tzAt = src.indexOf('toZ(');
  return (brAt !== -1 && brAt < tzAt)
    || 'the metric conversion is tried first, so a raw span would return null before the lookup runs';
});

/* All three display sites go through batteryRowPctCell, which is the single
   place that decides percentile-vs-base-rate and formats accordingly. Named
   rather than counted, so a refactor that moves calls around cannot pass while
   a site quietly computes its own. */
check('screen, in-place edit and APA export all render via batteryRowPctCell', () => {
  const missing = ['renderBattery', 'renderBatteryApa'].filter(fn => {
    let src; try { src = extractFn(APP_SRC, fn); } catch (e) { return true; }
    return !/batteryRowPctCell\(/.test(src);
  });
  /* The in-place handler lives inside renderBattery, so it is covered above,
     but assert its own path explicitly since it is the one that can go stale
     while the full render stays correct. */
  const rb = extractFn(APP_SRC, 'renderBattery');
  if (!/pctCell\.textContent = cellVal/.test(rb)) missing.push('in-place update');
  return missing.length === 0
    || 'these do not use batteryRowPctCell: ' + missing.join(', ');
});

check('batteryRowPctCell is the only thing that picks between the two quantities', () => {
  const src = extractFn(APP_SRC, 'batteryRowPctCell');
  if (!/baseRateAtOrAbove\(/.test(src)) return 'it no longer reports the published base rate';
  if (!/batteryRowPercentile\(/.test(src)) return 'it no longer falls back to the percentile rank';
  return true;
});

/* A base-rate measure must report a base rate REGARDLESS of what else is in the
   table. An earlier version switched the whole column between the two
   quantities depending on the other rows, which meant a longest-span value
   silently changed meaning when an unrelated subtest was entered — and, because
   the switch ignored seeded example rows while still rendering them, a table
   could print "Base rate" above a cell holding a percentile. */
check('a base-rate row reports a base rate whatever else is in the table', () => {
  const src = extractFn(APP_SRC, 'batteryRowPctCell');
  if (/mode/.test(src)) {
    return 'batteryRowPctCell takes a mode again, so a base-rate row can be converted '
         + 'depending on its neighbours';
  }
  // the base-rate branch must come first and return unconditionally
  const brAt = src.indexOf('batteryBaseRateEntry');
  const pctAt = src.indexOf('batteryRowPercentile');
  return (brAt !== -1 && brAt < pctAt)
    || 'the percentile path is reached before the base-rate one';
});

/* Because the quantity is per-measure, the label has to be per-section too. */
check('a base-rate section carries its own column heading, on screen and in export', () => {
  const bad = [];
  const rb = extractFn(APP_SRC, 'renderBattery');
  if (!/batteryGroupIsBaseRate\(gKey\)/.test(rb)) bad.push('editable table emits no section heading');
  if (!/group-col-header/.test(rb)) bad.push('editable table heading row has no class');
  const apa = extractFn(APP_SRC, 'buildApaTableFromColumns');
  if (!/groupLabel/.test(apa)) bad.push('APA builder ignores per-section column labels');
  if (!/apa-group-col-label/.test(apa)) bad.push('APA section column label has no class');
  /* The label shares the section-name row. Split across two rows it left a gap
     above the first data row and sat nearer the data than the heading. */
  if (!/apa-group apa-group-labelled/.test(apa)) {
    bad.push('the section label is no longer on the section-name row');
  }
  // and a section with no labels must keep the plain full-width name row
  if (!/colspan="\$\{visible\.length\}"/.test(apa)) {
    bad.push('an unlabelled section no longer spans the table width');
  }
  if (!/groupLabel: g => \(batteryGroupIsBaseRate\(g\)/.test(APP_SRC)) {
    bad.push('the percentile column defines no groupLabel, so the export loses the distinction');
  }
  return bad.length === 0 || bad.join('; ');
});

check('the main column heading is fixed and never switches quantity', () => {
  if (/batteryPctColumnMode|updateBatteryPctHeader/.test(APP_SRC)) {
    return 'the toggling heading is back; a longest-span row would change quantity '
         + 'when an unrelated subtest is entered';
  }
  return /const BAT_PCT_LABEL = 'Percentile'/.test(APP_SRC)
    || 'the fixed heading constant is gone';
});

/* The example row is seeded with isExample and is still rendered. Any rule that
   decides how to LABEL a column must therefore account for it, or the heading
   and the cells disagree — which is exactly how this broke the first time. */
check('seeded example rows cannot desynchronise a heading from its cells', () => {
  // with a per-row quantity there is no table-wide decision to get wrong
  const src = extractFn(APP_SRC, 'batteryRowPctCell');
  return !/isExample/.test(src)
    || 'batteryRowPctCell inspects isExample, which means the quantity is being '
     + 'decided table-wide again';
});

check('base rates are not run through fmtPct, which clamps 100 away', () => {
  const src = extractFn(APP_SRC, 'batteryRowPctCell');
  const brBranch = src.slice(0, src.indexOf('batteryRowPercentile'));
  return !/fmtPct\(/.test(brBranch)
    || 'the base-rate branch uses fmtPct, which would print a published 100% as 99.99';
});

{
  const ctx = {};
  vm.createContext(ctx);
  vm.runInContext(extractFn(APP_SRC, 'fmtBaseRate') + ';globalThis.__F = fmtBaseRate;', ctx);
  const F = ctx.__F;
  check('fmtBaseRate keeps the published precision, including a flat 100', () => {
    const cases = [[100, '100'], [99.5, '99.5'], [79, '79'], [61.5, '61.5'], [3, '3'], [0, '0']];
    const bad = cases.filter(([v, want]) => F(v) !== want).map(([v]) => v + ' -> ' + F(v));
    return bad.length === 0 || bad.join(', ');
  });
}

/* Score Tables collapses age-banded families to one entry, on the grounds that
   the band does not change the resulting table. True when the score comes from
   a metric conversion; false when it comes from an age-banded lookup. */
check('age bands stay selectable for families scored from a base-rate table', () => {
  if (!/function familyScoredByAgeBand/.test(APP_SRC)) return 'the exemption no longer exists';
  const src = extractFn(APP_SRC, 'buildFamilyListHtml');
  if (!/!members\.some\(familyScoredByAgeBand\)/.test(src)) {
    return 'the flat branch no longer exempts base-rate families, so every patient '
         + 'would be scored against All Ages';
  }
  return /flat && members\.some\(familyScoredByAgeBand\)/.test(src)
    || 'base-rate families are exempted from flattening but never rendered';
});

/* The size of what that exemption prevents, pinned to the shipped data. */
check('the age band moves a longest-span percentile by enough to matter', () => {
  const bands = Object.keys(D.normDB).filter(g => /^WAIS-IV Longest Span/.test(g));
  const pct = (e, v) => {
    const spans = Object.keys(e.baseRates).map(Number);
    const ge = x => { const a = spans.filter(s => s >= x); return a.length ? e.baseRates[Math.min(...a)] : 0; };
    return (100 - ge(v)) + 0.5 * (ge(v) - ge(v + 1));
  };
  const vals = bands
    .map(b => D.normDB[b]['Longest Digit Span Backward'])
    .filter(Boolean).map(e => pct(e, 4));
  const spread = Math.max(...vals) - Math.min(...vals);
  return spread > 20
    || 'spread is only ' + spread.toFixed(1) + ' points; if the data changed this much, '
     + 'revisit whether age bands still need to be selectable';
});

/* The merge rebuilds the table from each section's source and drops every
   .apa-group row on the way through, so a section that relabels a column has to
   have that label carried across explicitly. Otherwise a merged report prints
   base rates under the "Percentile" heading with nothing saying so. */
/* The report splits a table into one item per family BEFORE any merging, and
   the split drops the group row (its name becomes the item title). The column
   labels have to be lifted off it first, or they never reach the report at all
   and the merge has nothing left to carry. */
check('splitting a table into report items preserves a section\'s column label', () => {
  const src = extractFn(APP_SRC, 'extractGroupsFromHtml');
  const bad = [];
  if (!/colLabels/.test(src)) bad.push('the split does not lift the labels off the group row');
  if (!/apa-group-labelled/.test(src)) bad.push('the split re-emits no labelled row for the merge to find');
  if (!/section\.colLabels/.test(src)) bad.push('the captured labels are never written back');
  /* The re-emitted row must NOT be an .apa-group row. It has no section name,
     so anything hunting for a section divider has to skip it — when it did
     carry .apa-group, detectTestFamily read "Base rate" as the test family and
     the battery merge stopped matching. */
  // the ROW, not the cell — .apa-group-col-label on the <td> is correct
  if (/\blr\.className\s*=\s*'apa-group/.test(src)) {
    bad.push('the re-emitted label row is classed as a section divider again');
  }
  if (!/apa-col-label-row/.test(src)) bad.push('the re-emitted row has no distinguishing class');
  return bad.length === 0 || bad.join('; ');
});

/* A column label that does not line up with the numbers beneath it is worse
   than no label — it reads as belonging to the neighbouring column. Alignment
   is therefore never restated; it is taken from the column itself. */
check('a section column label inherits its column\'s alignment, never its own', () => {
  const css = fs.readFileSync(path.join(ROOT, 'styles.css'), 'utf8');
  const bad = [];
  for (const sel of ['apa-group-col-label', 'apa-col-label-row']) {
    const m = new RegExp('\\.apa-table tr[^{]*' + sel + '[^{]*\\{([^}]*)\\}').exec(css);
    if (!m) { bad.push(sel + ' rule missing'); continue; }
    // strip comments — the rule explains WHY it sets no alignment
    const decls = m[1].replace(/\/\*[\s\S]*?\*\//g, '');
    if (/text-align/.test(decls)) bad.push(sel + ' sets its own text-align');
  }
  return bad.length === 0 || bad.join('; ');
});

check('labels built after style inlining copy padding and alignment from their column', () => {
  const bad = [];
  const split = extractFn(APP_SRC, 'extractGroupsFromHtml');
  if (!/text-align\\s\*:\\s\*\(\[\^;\]\+\)/.test(split) && !/text-align/.test(split)) {
    bad.push('the split label does not copy alignment from a reference cell');
  }
  if (!/section\.rows\[0\]/.test(split)) bad.push('the split label has no reference cell');
  const merged = extractFn(APP_SRC, 'buildMergedTableHtml');
  if (!/p\.bodyRows\[0\]/.test(merged)) bad.push('the merged label has no reference cell');
  /* buildStandaloneHtml rewrites every group-row cell as a bold left-aligned
     section name; the label needs restoring after that or it inherits it. */
  const inliner = extractFn(APP_SRC, 'buildStandaloneHtml');
  if (!/apa-group-col-label/.test(inliner)) {
    bad.push('the Word inliner leaves the label styled as a section name');
  }
  return bad.length === 0 || bad.join('; ');
});

check('detectTestFamily reads only a group row\'s FIRST cell', () => {
  const src = extractFn(APP_SRC, 'detectTestFamily');
  if (/querySelectorAll\('table tbody tr\.apa-group td'\)/.test(src)) {
    return 'it scans every cell of a group row again, so a per-column label such '
         + 'as "Base rate" will be read as the test family';
  }
  return /tr\.children\[0\]/.test(src)
    || 'it no longer takes the first cell of each group row';
});

check('merging a battery preserves a section\'s own column label', () => {
  const src = extractFn(APP_SRC, 'buildMergedTableHtml');
  const bad = [];
  if (!/apa-group-labelled/.test(src)) bad.push('the merge never looks for a relabelled section');
  if (!/colLabels/.test(src)) bad.push('the merge does not carry the labels across');
  if (!/p\.colLabels\.some\(Boolean\)/.test(src)) bad.push('the merged divider never switches to per-cell form');
  /* Word does not receive the stylesheet, so the merged output has to inline
     its styling — the rest of this function already does. */
  if (!/font-style:italic/.test(src)) bad.push('the carried label has no inline styling for Word');
  return bad.length === 0 || bad.join('; ');
});

check('single-administration families are kept out of Change Analysis and SDI', () => {
  if (!/function isSingleAdministrationFamily/.test(APP_SRC)) return 'the filter no longer exists';
  const rci = extractFn(APP_SRC, 'populateFamilyList');
  const sdi = extractFn(APP_SRC, 'rebuildSdiFamilyList');
  const bad = [];
  if (!/isSingleAdministrationFamily/.test(rci)) bad.push('Change Analysis');
  if (!/isSingleAdministrationFamily/.test(sdi)) bad.push('SD Index');
  return bad.length === 0
    || bad.join(' and ') + ' would offer a family that can never compute anything';
});

/* ==========================================================================
   23. Score Tables SEM — the SD paired with the reliability

   SEM = SD x sqrt(1 - r) is only an error variance when both terms describe
   the SAME population. rCorrected is rescaled to the normative sample's
   variability, so it belongs with the normative SD (15 / 10 / 3 / 1); the raw
   r belongs with sd1. getBatteryRowReliability used to take rCorrected AND
   sd1, which is neither.

   PINNED SOURCE: CVLT-3 Manual (US, 2017), Table 3.4, "Alternate Form
   Reliability, Standard and Alternate Forms", pp. 38-39. That table publishes
   an SEM column alongside r and corrected r, which lets the pairing be tested
   against the publisher's own arithmetic rather than asserted.
   ========================================================================== */
heading('23. Score Tables SEM — the SD paired with the reliability');

/* CVLT-3 Manual Table 3.4, SEM column. Ages 16-44 then Ages 45-90. */
const CVLT3_SEM = {
  'CVLT-3 Trials · Ages 16-44': {
    'Trial 1': 2.10, 'Trial 2': 2.06, 'Trial 3': 1.64, 'Trial 4': 1.94, 'Trial 5': 2.12,
    'List B Correct': 1.97, 'Short Delay Free Recall': 1.62, 'Short Delay Cued Recall': 1.99,
    'Long Delay Free Recall': 1.87, 'Long Delay Cued Recall': 1.90, 'Recognition': 2.14,
    'Recognition False Positive': 1.97, 'Recognition Discrimination': 1.72,
    'Discrimination Nonparametric': 1.77, 'Total Intrusions': 2.60, 'Total Repetitions': 1.92
  },
  'CVLT-3 Trials · Ages 45-90': {
    'Trial 1': 1.99, 'Trial 2': 1.92, 'Trial 3': 1.82, 'Trial 4': 1.72, 'Trial 5': 1.62,
    'List B Correct': 2.06, 'Short Delay Free Recall': 1.62, 'Short Delay Cued Recall': 1.67,
    'Long Delay Free Recall': 1.62, 'Long Delay Cued Recall': 1.77, 'Recognition': 1.85,
    'Recognition False Positive': 2.01, 'Recognition Discrimination': 2.03,
    'Discrimination Nonparametric': 1.94, 'Total Intrusions': 1.97, 'Total Repetitions': 1.64
  },
  'CVLT-3 Indices · Ages 16-44': {
    'T1-5 Correct': 7.50, 'Delayed Recall Correct': 7.19, 'Total Recall Correct': 6.87
  },
  'CVLT-3 Indices · Ages 45-90': {
    'T1-5 Correct': 6.71, 'Delayed Recall Correct': 6.71, 'Total Recall Correct': 6.18
  }
};

check('the published SEM equals normative SD x sqrt(1 - corrected r)', () => {
  /* Derives the pairing rule from the source rather than assuming it. All 38
     published SEMs must fall out of the normative SD (15 for the indices,
     3 for the scaled trials) and the entry's own corrected r. */
  const bad = [];
  let n = 0;
  Object.entries(CVLT3_SEM).forEach(([group, tab]) => {
    const normSd = /Indices/.test(group) ? 15 : 3;
    Object.entries(tab).forEach(([name, published]) => {
      const e = D.normDB[group] && D.normDB[group][name];
      if (!e) { bad.push(group + ' / ' + name + ' is missing from normDB'); return; }
      n++;
      const derived = normSd * Math.sqrt(1 - e.rCorrected);
      // Published to 2dp, so half a unit in the last place is the tolerance.
      if (Math.abs(derived - published) > 0.005) {
        bad.push(name + ' [' + group + '] derived ' + derived.toFixed(3) + ', published ' + published);
      }
    });
  });
  if (n !== 38) return 'checked ' + n + ' measures, expected 38';
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the SHIPPED function reproduces the published CVLT-3 SEM', () => {
  /* The rule above is arithmetic; this asserts app.js actually implements it.
     Reads getBatteryRowReliability directly — the rendered interval is rounded
     for display, which would swamp a 2dp comparison. */
  const bad = [];
  Object.entries(CVLT3_SEM).forEach(([group, tab]) => {
    const type = /Indices/.test(group) ? 'standard' : 'scaled';
    Object.entries(tab).forEach(([name, published]) => {
      const rel = batteryRel({ name, group, scoreType: type }, type);
      if (!rel) { bad.push(name + ' [' + group + '] yielded no reliability'); return; }
      const sem = rel.sd * Math.sqrt(1 - rel.r);
      if (Math.abs(sem - published) > 0.005) {
        bad.push(name + ' [' + group + '] gave SEM ' + sem.toFixed(3) + ', published ' + published);
      }
    });
  });
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('rCorrected is never paired with sd1', () => {
  /* The defect this section exists for. RBANS Form C Delayed Memory is the
     widest gap in the database: its retest sample is impaired, so sd1 is 18.7
     against a normative 15, and rCorrected is .24. Pairing them gave a SEM of
     16.30 where 13.08 is correct — 6.3 standard-score points at each end of a
     95% interval. */
  const row = { name: 'Delayed Memory', group: 'RBANS Indices (Form C) · All Ages', scoreType: 'standard' };
  const e = D.normDB[row.group][row.name];
  const rel = batteryRel(row, 'standard');
  if (!rel) return 'no reliability returned for a row that has one';
  if (rel.sd === e.sd1) return 'the SD is still sd1 (' + e.sd1 + ')';
  return rel.sd === 15 || 'the SD is ' + rel.sd + ', expected the normative 15';
});

check('a raw-metric row still takes its SD from sd1', () => {
  /* Raw entries have no normative SD, so sd1 is the only one in their units.
     RBANS List Recognition is raw (see data.js metric:'raw'). */
  const row = { name: 'List Recognition', group: 'RBANS Subtests · Ages 20-89', scoreType: 'raw' };
  const e = D.normDB[row.group][row.name];
  if (e.metric !== 'raw') return 'List Recognition is no longer tagged raw — pick another raw row';
  const rel = batteryRel(row, 'raw');
  if (!rel) return 'a raw row yielded no reliability at all';
  return rel.sd === e.sd1 || 'raw row used SD ' + rel.sd + ', expected sd1 ' + e.sd1;
});

check('an entry entered on a standardised metric never uses its raw sd1', () => {
  /* Was a pinned known defect; now the thing it warned about. All 39 CVLT-C
     entries store raw statistics but are entered as z or T (reportedAs), so
     their SD must come from the entry metric, never from the raw sd1. At age 8
     the old behaviour gave Discriminability an interval of +/-7.5 on the z
     scale and Perseverations +/-4.3. */
  const bad = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    Object.entries(tab).forEach(([name, e]) => {
      if (!e || typeof e !== 'object') return;
      if (!(e.metric === 'raw' && e.reportedAs && e.reportedAs !== 'raw')) return;
      const rel = batteryRel({ name, group, scoreType: e.reportedAs }, e.reportedAs);
      if (!rel) { bad.push(name + ' [' + group + '] yielded no reliability'); return; }
      const want = { z: 1, t: 10, scaled: 3, standard: 15 }[e.reportedAs];
      if (rel.sd !== want) {
        bad.push(name + ' [' + group + '] entered as ' + e.reportedAs
               + ' but used SD ' + rel.sd + (rel.sd === e.sd1 ? ' (its raw sd1)' : '') + ', expected ' + want);
      }
    });
  });
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

/* ==========================================================================
   24. rInternal — internal-consistency reliability, CI only

   PINNED SOURCE: CVLT-C Manual (US), Chapter 6, pp. 85-88.
     Table 6.5  odd/even split-half (Spearman-Brown) and coefficient alpha for
                List A Trials 1-5 Total, ages 5-16, WITH the SEM for each.
     p. 87      states SEM = SD sqrt(1 - rxx) and directs the reader to add and
                subtract from the child's standard T score.
     pp. 87-88  worked example: an 8-year-old with T = 45 has a 95% interval of
                38-52 and a 90% interval of 39-51.
   ========================================================================== */
heading('24. rInternal — internal-consistency reliability, CI only');

/* CVLT-C Manual Table 6.5, odd/even row: rxx and the SEM printed beside it. */
const CVLTC_T65 = {
  5:  [0.88, 3.46], 6:  [0.85, 3.87], 7:  [0.87, 3.60], 8:  [0.87, 3.60],
  9:  [0.91, 3.00], 10: [0.88, 3.46], 11: [0.87, 3.60], 12: [0.89, 3.31],
  13: [0.88, 3.46], 14: [0.89, 3.31], 15: [0.86, 3.74], 16: [0.84, 4.00]
};

check('Table 6.5 SEMs are computed on the T metric, not the raw SD', () => {
  /* The finding the whole change rests on. If these SEMs used the raw SD they
     would not reproduce from a constant; they reproduce from SD = 10 in every
     one of the 12 age bands.

     The manual TRUNCATES rather than rounds — 3.6055 prints as 3.60, 3.3166 as
     3.31 — so the assertion is that the printed value sits in [derived - 0.01,
     derived]. Stated as an interval rather than by truncating in JS, because
     sqrt(0.09) lands a hair under 0.3 in binary floating point and Math.floor
     would turn an exact 3.00 into 2.99. */
  const bad = [];
  const EPS = 1e-9;
  Object.entries(CVLTC_T65).forEach(([age, [rxx, printed]]) => {
    const derived = 10 * Math.sqrt(1 - rxx);
    if (!(printed <= derived + EPS && printed > derived - 0.01 - EPS)) {
      bad.push('age ' + age + ': 10*sqrt(1-' + rxx + ') = ' + derived.toFixed(4)
             + ', manual prints ' + printed);
    }
  });
  return bad.length === 0 || bad.join('; ');
});

check("the manual's worked example reproduces exactly (age 8, T 45)", () => {
  /* CVLT-C Manual pp. 87-88. Uses the shipped renderer, so this is an
     end-to-end assertion: data, reliability choice, SD choice and rounding. */
  const row = { name: 'List A Trials 1-5 Total',
                group: 'CVLT-C Subtests (Raw Scores) · Age 8', scoreType: 't' };
  const got95 = batteryCi(45, row, '95');
  const got90 = batteryCi(45, row, '90');
  const bad = [];
  if (got95 !== '38–52') bad.push('95% gave ' + got95 + ', manual prints 38-52');
  if (got90 !== '39–51') bad.push('90% gave ' + got90 + ', manual prints 39-51');
  return bad.length === 0 || bad.join('; ');
});

check('the three List A Trials 1-5 entries carry the published split-half r', () => {
  const want = { 8: 0.87, 12: 0.89, 16: 0.84 };   // Table 6.5, odd/even
  const bad = [];
  Object.entries(want).forEach(([age, rxx]) => {
    const e = D.normDB['CVLT-C Subtests (Raw Scores) · Age ' + age]['List A Trials 1-5 Total'];
    if (e.rInternal !== rxx) bad.push('age ' + age + ' has rInternal ' + e.rInternal + ', Table 6.5 gives ' + rxx);
    /* The retest coefficient must survive alongside it — Change Analysis
       still needs it, and it is a different quantity. */
    if (!Number.isFinite(e.r)) bad.push('age ' + age + ' lost its retest r');
  });
  return bad.length === 0 || bad.join('; ');
});

check('rInternal is used ONLY for the confidence interval, never by RCI', () => {
  /* The load-bearing separation. With the retest r, SD sqrt(2(1-r)) is the
     spread of the difference distribution — the quantity reliable change is
     measured against. An internal-consistency coefficient strips out real
     between-session fluctuation, narrows the interval and overcalls change;
     on SRB and Crawford it also moves a fitted regression slope. */
  const RCI_FNS = ['calcBasicRow', 'calcPracticeRow', 'calcSrbRow', 'calcCrawfordRow', 'rciEffectiveR'];
  const bad = RCI_FNS.filter((n) => /rInternal/.test(extractFn(APP_SRC, n)));
  if (bad.length) return bad.join(', ') + ' read rInternal; reliable change must use the retest coefficient';
  /* And it must actually be reachable from the CI path, or the field is inert.
     The preference chain moved into resolveCiReliability when Score Tables and
     the Data page were made to share one resolver, so that is where the read
     now has to be — and getBatteryRowReliability has to still call it. */
  if (!/rInternal/.test(extractFn(APP_SRC, 'resolveCiReliability'))) {
    return 'resolveCiReliability never reads rInternal, so the field does nothing';
  }
  return /resolveCiReliability\(/.test(extractFn(APP_SRC, 'getBatteryRowReliability'))
    || 'getBatteryRowReliability no longer calls resolveCiReliability, so the CI path has left the chain';
});

check('rInternal appears only where a publisher derives its own intervals from it', () => {
  /* Not a general-purpose field. Adding it to a measure whose manual does not
     publish an internal-consistency interval would silently change that row's
     basis, so the roster is pinned.

     Seven sources so far:
       CVLT-C List A Trials 1-5 Total   3 entries, no by-age table
       D-KEFS Advanced CWIT/Tower/SST/RISK  26 entries, All Ages groups only
       D-KEFS (original)                13 entries, All Ages groups only
       WAIS-IV Technical Manual Table 4.1   21 entries, All Ages groups only
       RBANS Update Table 3.6            9 entries, All Ages groups only
       WMS-IV Technical Manual Table 3.1 30 entries, on the two BATTERY groups
       WISC-V Technical Manual Table 4.1 29 entries, All Ages groups only

     RBANS is the first source whose reliability table is MIXED-BASIS within
     itself: five of its fourteen rows carry footnote a, "estimates based on
     test-retest", and are held in rStability* instead. They must not appear
     here, and section 29 asserts that separately.

     Keyed on EITHER field. D-KEFS original carries rInternalByAge with no
     rInternal — its manual publishes no all-ages average — so a roster that
     tested rInternal alone would let 11 entries change the basis of a printed
     interval without ever being counted. See section 27. */
  const carriers = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    Object.entries(tab).forEach(([name, e]) => {
      if (e && typeof e === 'object' && (Number.isFinite(e.rInternal) || e.rInternalByAge)) {
        carriers.push(group + ' / ' + name);
      }
    });
  });
  const cvltc = carriers.filter((c) => /^CVLT-C .* \/ List A Trials 1-5 Total$/.test(c));
  const dkefs = carriers.filter((c) => /^D-KEFS Advanced (Colour-Word Interference|Tower|Social Sorting|Risk-Reward Decision) · All Ages \//.test(c));
  const orig  = carriers.filter((c) => /^D-KEFS (?!Advanced)(Trail Making Test|Verbal Fluency|Sorting Test|Word Context Test|Tower Test|Word Proverb Test|Twenty Questions Test) · All Ages \//.test(c));
  const wais  = carriers.filter((c) => /^WAIS-IV (Core Subtests|Indices|Process Scores|Supplementary Subtests) · All Ages \//.test(c));
  const rbans = carriers.filter((c) => /^RBANS (Subtests|Indices) · All Ages \//.test(c));
  const wms   = carriers.filter((c) => /^WMS-IV (Subtests|Indices) · Ages (16-69|65-90) \//.test(c));
  const wisc  = carriers.filter((c) => /^WISC-V (Subtests|Indices|Process Scores) · All Ages \//.test(c));
  const stray = carriers.filter((c) => !cvltc.includes(c) && !dkefs.includes(c) && !orig.includes(c) && !wais.includes(c) && !rbans.includes(c) && !wms.includes(c) && !wisc.includes(c));
  const bad = [];
  if (cvltc.length !== 3) bad.push('CVLT-C carriers: ' + cvltc.length + ', expected 3');
  if (dkefs.length !== 26) bad.push('D-KEFS Advanced carriers: ' + dkefs.length + ', expected 26');
  if (orig.length !== 13) bad.push('D-KEFS original carriers: ' + orig.length + ', expected 13');
  if (wais.length !== 21) bad.push('WAIS-IV carriers: ' + wais.length + ', expected 21');
  if (rbans.length !== 9) bad.push('RBANS carriers: ' + rbans.length + ', expected 9');
  /* 30, not 32: Table 3.1 gives 32 measure-battery cells (20 adult rows + 12
     older adult), and VPA II Word Recall carries footnote b in BOTH batteries,
     so its 2 are held in rStability and counted by section 30 instead. */
  if (wms.length !== 30) bad.push('WMS-IV carriers: ' + wms.length + ', expected 30');
  /* 19, not 22: Coding, Symbol Search and Cancellation are the three the
     manual names as improper for split-half, and are held in rStability. */
  /* 29 = 13 core subtests (16 less the three speeded) + 11 indices (6 primary
     + 5 ancillary) + 5 process scores (7 less the two Cancellation ones).
     The five excluded are held in rStability and counted by section 29. */
  if (wisc.length !== 29) bad.push('WISC-V carriers: ' + wisc.length + ', expected 29');
  if (stray.length) bad.push('undocumented: ' + stray.slice(0, 3).join(', '));
  return bad.length === 0 || bad.join('; ');
});

check('Methods & References names every instrument whose CI uses internal consistency', () => {
  /* IT WENT STALE ONCE, SILENTLY. The paragraph said internal consistency
     applied to "the CVLT-C List A Trials 1-5 Total score alone" while 42
     entries across 12 families were already using it — the sentence was not
     updated when D-KEFS Advanced was added, and nothing failed.

     On-screen text is a contract (see CLAUDE.md), and this one describes the
     basis of every printed interval, so it has to track the data. Checked at
     INSTRUMENT level rather than per measure, because the roster reasonably
     groups measures by test. */
  const para = methodsCiRoster();
  const instruments = new Set();
  Object.entries(D.normDB).forEach(([group, tab]) => {
    Object.values(tab).forEach((e) => {
      if (!e || typeof e !== 'object') return;
      if (!Number.isFinite(e.rInternal) && !e.rInternalByAge) return;
      instruments.add(group.startsWith('D-KEFS Advanced') ? 'D-KEFS Advanced' : group.split(' ')[0]);
    });
  });
  const missing = [...instruments].filter((i) => !para.includes(i));
  return missing.length === 0
    || missing.join(', ') + ' use an internal-consistency coefficient but are not named in Methods & References';
});

check('D-KEFS Advanced TMT and VFT stay on stability coefficients', () => {
  /* Their manual says split-half and alpha do not give accurate estimates for
     SPEEDED measures and uses stability coefficients for these two tests
     instead. Adding rInternal here would substitute a coefficient the
     publisher deliberately did not use. */
  const bad = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    if (!/^D-KEFS Advanced (Trail Making|Verbal Fluency)/.test(group)) return;
    Object.entries(tab).forEach(([name, e]) => {
      if (e && typeof e === 'object' && (Number.isFinite(e.rInternal) || e.rInternalByAge)) {
        bad.push(group + ' / ' + name);
      }
    });
  });
  return bad.length === 0
    || bad.slice(0, 3).join('; ') + ' carry an internal-consistency coefficient; these tests are speeded';
});

check('a by-age table lives on an All Ages group, or on a separate battery', () => {
  /* THE RULE: a normative-band lookup may not sit on a group whose own band
     came from a RETEST study, because the two banding schemes are unrelated.
     D-KEFS Advanced is the case it was written for — its groups are the retest
     bands 8-18/19-59/60-90, nothing like Table 3.4's normative bands. Score
     Tables collapses those families to All Ages, and only Change Analysis and
     the SD Index see the banded ones, where the retest r is wanted anyway.

     THE EXCEPTION, and it is principled rather than a carve-out: WMS-IV's two
     groups are not retest bands of one norm set, they are two BATTERIES. The
     manual norms an Adult battery over 16-69 and an Older Adult battery over
     65-90, Table 3.1 prints a column block for each, and the retest study used
     those same two ranges. So the group's range IS the battery's range, and
     the by-age keys are that battery's own normative sub-bands, lying wholly
     inside it. There is no mismatch to protect against.

     Marked with separateBattery rather than inferred from the group name, so
     the exemption is declared by the data and greppable. Anything else banded
     that acquires a by-age table still fails, which is the point. */
  const bad = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    if (/· All Ages$/.test(group)) return;
    Object.entries(tab).forEach(([name, e]) => {
      if (!e || typeof e !== 'object') return;
      if (!e.rInternalByAge && !e.rStabilityByAge) return;
      if (e.separateBattery) return;
      bad.push(group + ' / ' + name);
    });
  });
  if (bad.length) return bad.slice(0, 3).join('; ') + ' — a retest-banded group cannot carry a normative-band table';
  /* ...and separateBattery must not spread. It is a claim about how an
     instrument is published, not a way round the rule above. */
  const marked = new Set();
  Object.entries(D.normDB).forEach(([group, tab]) => {
    Object.values(tab).forEach((e) => { if (e && typeof e === 'object' && e.separateBattery) marked.add(group); });
  });
  const stray = [...marked].filter((g) => !/^WMS-IV (Subtests|Indices) · Ages (16-69|65-90)$/.test(g));
  if (stray.length) return 'separateBattery has escaped WMS-IV: ' + stray.join(', ');
  if (marked.size !== 4) return 'expected 4 separateBattery groups, found ' + marked.size;
  return true;
});

/* ==========================================================================
   25. Age-banded reliability (rInternalByAge)

   PINNED SOURCE: D-KEFS Advanced Manual, Table 3.4 (reliability coefficients
   by age group and overall sample). Its bands are the NORMATIVE ones — single
   years 8-15, then 16-18, 19-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-90 —
   which is why the table cannot be folded into the retest study's bands.
   ========================================================================== */
heading('25. Age-banded reliability (rInternalByAge)');

const CWIT_ALL = 'D-KEFS Advanced Colour-Word Interference · All Ages';
const SST_ALL  = 'D-KEFS Advanced Social Sorting · All Ages';

check('spot values match Table 3.4 as printed', () => {
  /* Transcription guard. Taken from three different tests and both ends of
     the age range, so a systematic off-by-one column shift would surface. */
  const want = [
    [CWIT_ALL, 'Inhibition Net Correct Responses',            8,  0.69],
    [CWIT_ALL, 'Inhibition Net Correct Responses',            80, 0.95],
    [CWIT_ALL, 'Inhibition/Switching Net Correct Responses',  10, 0.90],
    [SST_ALL,  'Total Number of Conceptual Level Responses',  15, 0.79],
    [SST_ALL,  'Percent Conceptual Level Responses',          80, 0.96],
    ['D-KEFS Advanced Tower · All Ages', 'Adjusted Mean Unproductive Responses', 14, 0.22]
  ];
  const bad = [];
  want.forEach(([g, n, age, rxx]) => {
    const got = D.normDB[g][n].rInternalByAge[age];
    if (got !== rxx) bad.push(n + ' age ' + age + ': stored ' + got + ', Table 3.4 prints ' + rxx);
  });
  return bad.length === 0 || bad.join('; ');
});

check('the lookup takes the band containing the age, not the nearest key', () => {
  /* Keys are band LOWER bounds, so 12 -> the 12 band but 17 -> the 16-18 band
     and 45 -> the 40-49 band. Getting this wrong would silently read a
     neighbouring band for every patient not born on a band boundary. */
  const e = D.normDB[CWIT_ALL]['Inhibition Net Correct Responses'];
  const cases = [
    [8, e.rInternalByAge[8]], [12, e.rInternalByAge[12]],
    [16, e.rInternalByAge[16]], [17, e.rInternalByAge[16]], [18, e.rInternalByAge[16]],
    [19, e.rInternalByAge[19]], [29, e.rInternalByAge[19]],
    [45, e.rInternalByAge[40]], [59, e.rInternalByAge[50]], [90, e.rInternalByAge[80]]
  ];
  const bad = [];
  cases.forEach(([age, expected]) => {
    const got = rInternalForAge(e, age);
    if (got !== expected) bad.push('age ' + age + ' gave ' + got + ', expected ' + expected);
  });
  return bad.length === 0 || bad.join('; ');
});

check('an absent or out-of-range age falls back to the published average', () => {
  /* The interval must never empty because the age box is blank, and must not
     extrapolate past the normed range. Both fall back to rInternal, which is
     the manual's own all-ages figure — a citable answer either way. */
  const row = { name: 'Inhibition Net Correct Responses', group: CWIT_ALL, scoreType: 'scaled' };
  const e = D.normDB[CWIT_ALL]['Inhibition Net Correct Responses'];
  const bad = [];
  [null, undefined, NaN, 4, 91, 200].forEach((age) => {
    if (rInternalForAge(e, age) !== null) bad.push('rInternalForAge did not decline age ' + age);
    const rel = batteryRel(row, 'scaled', age);
    if (!rel) { bad.push('no interval at age ' + age); return; }
    if (rel.r !== e.rInternal) bad.push('age ' + age + ' used r ' + rel.r + ', expected the average ' + e.rInternal);
  });
  // ...and a valid age must actually change it, or none of this does anything.
  const inBand = batteryRel(row, 'scaled', 8);
  if (inBand.r !== e.rInternalByAge[8]) bad.push('age 8 did not pick up the band coefficient');
  if (inBand.r === e.rInternal) bad.push('age 8 and the average are indistinguishable — pick a sharper probe');
  return bad.length === 0 || bad.join('; ');
});

check('age changes the printed interval on an age-banded measure', () => {
  /* End to end through the shipped renderer. Social Sorting Conceptual Level
     runs .79 at age 15 and .95 at 9, which must not print the same interval. */
  const row = { name: 'Total Number of Conceptual Level Responses', group: SST_ALL, scoreType: 'scaled' };
  const young = batteryCi(10, row, '95', 9);
  const older = batteryCi(10, row, '95', 15);
  if (!young || !older) return 'one of the intervals did not render';
  return young !== older
    || 'age 9 and age 15 both gave ' + young + ', so the age is not reaching the renderer';
});

check('the APA note states the age whenever one changed a coefficient', () => {
  /* A single field silently driving every interval in the table is only
     defensible if the table says which age it rests on. */
  const c = {};
  vm.createContext(c);
  vm.runInContext(APP_SRC.match(/const APA_NOTES = \{[\s\S]*?\n\};/)[0] + ';globalThis.__N = APA_NOTES;', c);
  const note = (ciAge) => c.__N.bat({ classification: 'wechsler', ciLevel: '95', ciAge }).filter(Boolean).join(' ');
  const bad = [];
  if (!/age 34/.test(note(34))) bad.push('the note does not name the age it used');
  if (/age \d/.test(note(null))) bad.push('the note claims an age on a table where none was used');
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   26. Shared patient age

   One patient, one age, two inputs (#bat-age, #pre-age).

   TWO premorbid models read age: OPIE-4 (calcPremorbid) and the ToPF-predicted
   WMS-IV indices (calcPredict). Both are adult-only. Before the age was shared,
   #pre-age was hard-bounded 16-90 by the markup, so "age is present" was a
   proxy for "age is valid" and the WMS-IV branch tested only `age != null`.
   That proxy is gone: the box now accepts from 5. Every age consumer must
   therefore bound the value itself, or a paediatric age typed on Score Tables
   returns a plausible-looking adult index.

   These checks exist because that is exactly what went wrong while writing
   this section — the WMS-IV age term was missed on a first read.
   ========================================================================== */
heading('26. Shared patient age');

check('the ToPF-predicted WMS-IV equations bound the age they use', () => {
  /* WMS_COEF entries carry an explicit age term. Presence alone is not enough. */
  const src = extractFn(APP_SRC, 'calcPredict');
  const body = src.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/[^\n]*/g, '');
  if (!/\bage\b/.test(body)) return true;          // fine if it stops using age
  const bad = [];
  if (!/TOPF_AGE_MIN/.test(body) || !/TOPF_AGE_MAX/.test(body)) {
    bad.push('age is used but not bounded by TOPF_AGE_MIN/TOPF_AGE_MAX');
  }
  if (/if \(topf != null && age != null\)/.test(body)) {
    bad.push('the WMS-IV branch still gates on presence alone');
  }
  return bad.length === 0 || bad.join('; ');
});

check('the ToPF age bounds are their own constants, not borrowed from OPIE', () => {
  /* Same numbers today, different publishers. Sharing the constant would let a
     revision to one model silently move the other model's range. */
  const bad = [];
  if (!/const TOPF_AGE_MIN\s*=\s*\d+/.test(DATA_SRC)) bad.push('TOPF_AGE_MIN is not defined in data.js');
  if (!/const TOPF_AGE_MAX\s*=\s*\d+/.test(DATA_SRC)) bad.push('TOPF_AGE_MAX is not defined in data.js');
  const predict = extractFn(APP_SRC, 'calcPredict');
  if (/OPIE_AGE_(MIN|MAX)/.test(predict)) bad.push('calcPredict bounds ToPF ages with the OPIE constants');
  return bad.length === 0 || bad.join('; ');
});

check('every use of age in calcPremorbid is range-gated', () => {
  /* Not merely "is age present" — an adult-only equation fed age 9 would
     return a plausible-looking FSIQ for a child. */
  const src = extractFn(APP_SRC, 'calcPremorbid');
  if (!/OPIE_AGE_MIN/.test(src) || !/OPIE_AGE_MAX/.test(src)) {
    return 'calcPremorbid no longer bounds age by OPIE_AGE_MIN/OPIE_AGE_MAX';
  }
  // The prediction itself must sit behind that gate.
  return /opieAgeOk\s*&&/.test(src)
    || 'the OPIE-4 prediction is no longer guarded by opieAgeOk';
});

check('the two age inputs are declared as one shared value', () => {
  const bad = [];
  if (!/const PATIENT_AGE_INPUTS = \[[^\]]*'bat-age'[^\]]*'pre-age'[^\]]*\]/.test(APP_SRC)) {
    bad.push('PATIENT_AGE_INPUTS no longer names both inputs');
  }
  const sync = extractFn(APP_SRC, 'syncPatientAge');
  /* Dispatching on the sibling would make the two inputs echo each other
     without end, so the mirror must only assign .value. */
  if (/dispatchEvent/.test(sync)) bad.push('syncPatientAge dispatches an event on the sibling input — infinite echo');
  if (!/if \(id === sourceId\) return;/.test(sync)) bad.push('syncPatientAge does not skip the edited input');
  return bad.length === 0 || bad.join('; ');
});

check('#pre-age accepts the full shared range, not just adults', () => {
  /* It holds the shared age now. Leaving min="16" would make the browser
     flag a legitimate paediatric age as invalid on a page that no longer
     owns the value. */
  const m = HTML_SRC.match(/<input[^>]*id="pre-age"[^>]*>/);
  if (!m) return 'the #pre-age input is gone';
  const min = m[0].match(/min="(\d+)"/);
  if (!min) return 'no min attribute at all — was it removed by accident?';
  return Number(min[1]) <= 5
    || 'min="' + min[1] + '" still excludes the ages Score Tables is normed for';
});

check('the age input is reachable, not stranded in the hidden legacy panel', () => {
  /* IT SHIPPED UNREACHABLE ONCE. #bat-age is declared inside
     .bat-premorbid-block, which buildBatteryInlineBar hides wholesale with
     ds-legacy-hidden — so the field rendered at 0x0 and no user could ever
     type in it. Every check passed, because they all called the reliability
     code directly and never asked whether the input could be reached.

     It stays in the markup (app.js wires its listener at init and
     design-system.js is deferred, so the element must pre-exist) and is MOVED
     into the control bar. These assertions pin that arrangement. */
  const DS_SRC = fs.readFileSync(path.join(ROOT, 'design-system.js'), 'utf8');
  const bad = [];
  const fn = DS_SRC.slice(DS_SRC.indexOf('function buildBatteryInlineBar'));
  if (!/ds-inline-bar-age/.test(fn)) bad.push('the control bar has no slot for the age field');
  if (!/appendChild\(ageInput\)/.test(fn)) bad.push('the age input is never moved into the bar');
  /* Rebuilding rather than moving would drop the listener app.js already
     attached, leaving a box that looks right and changes nothing. */
  if (/createElement\(['"]input['"]\)[^;]*bat-age/.test(fn)) bad.push('the age input is rebuilt rather than moved');
  if (!/ageSlot\.hidden\s*=/.test(fn)) bad.push('the age slot is never revealed with the CI toggle');
  /* And setting .hidden must actually hide it. The slot is also a
     .ds-inline-bar-section, whose display:flex outranks the browser default
     for [hidden] — so without an explicit rule the box stayed on screen with
     the CI toggle off. It did exactly that before this rule was added. */
  const DS_CSS = fs.readFileSync(path.join(ROOT, 'design-system.css'), 'utf8');
  if (!/\.ds-inline-bar-age\[hidden\]\s*\{[^}]*display:\s*none\s*!important/.test(DS_CSS)) {
    bad.push('nothing overrides .ds-inline-bar-section display:flex, so [hidden] will not hide the age box');
  }
  /* The family dropdown is pinned to the search box's width (left:0/right:0),
     and the search sits in the flexible section, so ANY section added to this
     bar narrows the dropdown. Adding the age box clipped its header row and
     wrapped the family names. The floor is what stops that recurring. */
  if (!/\.ds-inline-bar-search \.combo-list\{[^}]*min-width:/.test(DS_CSS.replace(/\/\*[\s\S]*?\*\//g, ''))) {
    bad.push('the family dropdown has no min-width, so a crowded bar will clip it again');
  }
  /* And it must still be declared in the markup, or the move has nothing to
     move and app.js has nothing to bind. */
  if (!/id="bat-age"/.test(HTML_SRC)) bad.push('#bat-age is no longer declared in index.html');
  return bad.length === 0 || bad.join('; ');
});

check('an out-of-range age is explained rather than left blank', () => {
  const src = extractFn(APP_SRC, 'calcPremorbid');
  const bad = [];
  if (!/pre-age-range-note/.test(src)) bad.push('calcPremorbid never touches the explanation box');
  if (!/OPIE_AGE_MIN\s*\|\|\s*age\s*>\s*OPIE_AGE_MAX/.test(src)) bad.push('the box is not driven by the OPIE range');
  if (!/<div class="caution-box" id="pre-age-range-note"/.test(HTML_SRC)) bad.push('the box is missing from the markup');
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   27. D-KEFS (original) internal consistency

   PINNED SOURCE: D-KEFS Technical Manual, Chapter 2 "Evidence of Reliability",
   pp. 18-19 and the per-test tables that follow.

     p. 18   "For some D-KEFS measures, the internal consistency coefficients
              were used; for other measures the test-retest correlations were
              employed."
     p. 19   SEM = SD sqrt(1 - rxx), and "The standard deviation unit is
              3 for all D-KEFS scaled scores."  CI = observed +/- z(SEM).
     p. 19   internal-consistency SEMs are tabulated BY AGE; test-retest SEMs
              "were derived from the total sample of cases", the by-band retest
              samples being too small.

   So the manual runs two regimes and picks per test. Tables 2.1/2.4/2.12/2.18/
   2.21/2.24 give the coefficients; 2.3/2.6/2.14/2.20/2.23/2.26 give the SEM
   they produce. Design Fluency has no coefficient table at all -- "Item
   interdependence precluded the use of internal consistency procedures" -- and
   its Table 2.8 is All Ages only, from the retest r.

   ONE MEASURE IS AN EXCEPTION. Twenty Questions Initial Abstraction is the only
   place in the manual where the coefficient table and its SEM table cannot both
   be right. The stored values are Table 2.15 verbatim, and the discrepancy is
   pinned as a fact rather than smoothed away. Full reasoning in data.js.
   ========================================================================== */
heading('27. D-KEFS (original) internal consistency');

/* The SEM column of each internal-consistency test's SEM table, keyed by the
   age band's LOWER bound to match rInternalByAge. */
const DKEFS_SEM = [
  ['D-KEFS Trail Making Test · All Ages', 'Combined Number + Letter',
   { 8:1.41, 9:1.59, 10:1.96, 11:1.92, 12:1.69, 13:1.66, 14:1.38, 15:1.59, 16:1.66, 20:1.41, 30:1.41, 40:1.52, 50:1.31, 60:1.33, 70:1.89, 80:1.43 }],  // Tables 2.1/2.3
  ['D-KEFS Verbal Fluency · All Ages', 'Letter Fluency',
   { 8:1.69, 9:1.62, 10:1.33, 11:1.48, 12:1.43, 13:1.31, 14:1.36, 15:1.41, 16:1.36, 20:1.16, 30:0.97, 40:1.45, 50:0.97, 60:1.16, 70:1.08, 80:1.13 }],  // Tables 2.4/2.6
  ['D-KEFS Verbal Fluency · All Ages', 'Category Fluency',
   { 8:1.89, 9:1.50, 10:1.62, 11:1.94, 12:1.59, 13:1.76, 14:1.71, 15:2.06, 16:1.89, 20:1.87, 30:1.48, 40:1.82, 50:1.85, 60:1.80, 70:1.78, 80:1.48 }],  // Tables 2.4/2.6
  ['D-KEFS Verbal Fluency · All Ages', 'Category Switching',
   { 8:2.37, 9:2.06, 10:1.99, 11:1.85, 12:1.85, 13:1.85, 14:2.03, 15:2.25, 16:2.15, 20:2.27, 30:1.71, 40:1.69, 50:2.23, 60:2.13, 70:2.03, 80:2.01 }],  // Tables 2.4/2.6
  ['D-KEFS Verbal Fluency · All Ages', 'Switching Accuracy',
   { 8:2.06, 9:1.57, 10:1.80, 11:1.48, 12:1.78, 13:2.03, 14:1.55, 15:1.82, 16:2.06, 20:1.92, 30:1.59, 40:1.85, 50:2.06, 60:2.11, 70:1.80, 80:1.62 }],  // Tables 2.4/2.6
  ['D-KEFS Sorting Test · All Ages', 'Free Sorting Confirmed Sorts',
   { 8:1.92, 9:1.94, 10:1.33, 11:1.64, 12:1.85, 13:1.57, 14:1.26, 15:2.01, 16:1.59, 20:1.41, 30:1.26, 40:1.31, 50:1.13, 60:1.31, 70:1.31, 80:1.43 }],  // Tables 2.12/2.14
  ['D-KEFS Sorting Test · All Ages', 'Free Sorting Description Total Score',
   { 8:1.85, 9:1.80, 10:1.43, 11:1.57, 12:1.80, 13:1.64, 14:1.33, 15:2.01, 16:1.55, 20:1.43, 30:1.24, 40:1.33, 50:1.19, 60:1.36, 70:1.28, 80:1.45 }],  // Tables 2.12/2.14
  ['D-KEFS Sorting Test · All Ages', 'Sort Recognition Total Description Score',
   { 8:1.52, 9:1.62, 10:1.85, 11:1.59, 12:1.73, 13:1.85, 14:1.59, 15:1.59, 16:1.52, 20:1.50, 30:1.45, 40:1.36, 50:1.52, 60:1.31, 70:1.38, 80:1.64 }],  // Tables 2.12/2.14
  ['D-KEFS Word Context Test · All Ages', 'Total First Trial Consistently Correct',
   { 8:2.01, 9:2.08, 10:2.08, 11:1.92, 12:2.08, 13:2.18, 14:1.71, 15:1.62, 16:2.11, 20:1.69, 30:1.73, 40:2.06, 50:1.59, 60:1.71, 70:1.52, 80:1.69 }],  // Tables 2.18/2.20
  ['D-KEFS Tower Test · All Ages', 'Total Achievement Score',
   { 8:1.99, 9:1.62, 10:1.21, 11:1.87, 12:1.87, 13:2.01, 14:2.27, 15:1.89, 16:1.89, 20:1.85, 30:1.59, 40:1.59, 50:1.99, 60:1.59, 70:1.41, 80:1.87 }],  // Tables 2.21/2.23
  ['D-KEFS Word Proverb Test · All Ages', 'Total Achievement Score',
   { 16:1.69, 20:1.62, 30:1.33, 40:1.48, 50:1.43, 60:1.31, 70:1.36, 80:1.41 }],  // Tables 2.24/2.26
  /* Twenty Questions' OTHER column. It reproduces cleanly, which is what makes
     its neighbour's failure a property of that one column rather than of the
     table or the test — see TWENTY_Q_IA_SEM below. */
  ['D-KEFS Twenty Questions Test · All Ages', 'Total Weighted Achievement',
   { 8:2.25, 9:2.15, 10:2.20, 11:2.30, 12:2.11, 13:2.35, 14:2.06, 15:2.37, 16:2.85, 20:2.13, 30:2.37, 40:2.45, 50:2.58, 60:2.18, 70:2.01, 80:2.15 }],  // Tables 2.15/2.17
];

/* Table 2.17, Initial Abstraction Score. Deliberately NOT in DKEFS_SEM: this
   column does not follow from Table 2.15, and the check below pins that fact
   rather than papering over it. */
const TWENTY_Q_IA_SEM = { 8:1.11, 9:1.59, 10:1.45, 11:1.33, 12:1.31, 13:1.33, 14:1.36, 15:1.24,
                          16:1.50, 20:1.24, 30:1.64, 40:1.76, 50:1.21, 60:1.19, 70:1.21, 80:1.26 };

check('the stored coefficient is the one the published SEM was built from', () => {
  /* The finding the whole change rests on, and a transcription guard on all
     168 coefficients at once.

     Stated as "the rxx implied by the printed SEM agrees with the stored rxx",
     because the manual computes from UNROUNDED coefficients and prints both
     quantities to 2dp. Two roundings of the same number can differ by half a
     unit in the last place at each end; the worst case here is 0.0064, so the
     tolerance is 0.007. Asserting equality of the SEM instead would fail on
     seven cells for pure typography. */
  const bad = [];
  let n = 0;
  DKEFS_SEM.forEach(([group, name, tab]) => {
    const e = D.normDB[group] && D.normDB[group][name];
    if (!e) { bad.push(group + ' / ' + name + ' is missing from normDB'); return; }
    if (!e.rInternalByAge) { bad.push(group + ' / ' + name + ' carries no rInternalByAge'); return; }
    Object.entries(tab).forEach(([age, sem]) => {
      n++;
      const stored = e.rInternalByAge[age];
      if (!Number.isFinite(stored)) { bad.push(name + ' has no band ' + age); return; }
      const implied = 1 - Math.pow(sem / 3, 2);
      if (Math.abs(implied - stored) > 0.007) {
        bad.push(name + ' [' + group.replace('D-KEFS ', '').replace(' · All Ages', '') + '] age ' + age
               + ': stored ' + stored + ', SEM ' + sem + ' implies ' + implied.toFixed(3));
      }
    });
  });
  if (n !== 184) return 'checked ' + n + ' cells, expected 184';
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('Twenty Questions Initial Abstraction is Table 2.15, verbatim', () => {
  /* The one measure where the manual contradicts itself, so the stored values
     are pinned digit for digit against the coefficient table rather than being
     allowed to drift toward whatever makes Table 2.17 reconcile. */
  const T215 = { 8:0.85, 9:0.72, 10:0.76, 11:0.83, 12:0.82, 13:0.81, 14:0.80, 15:0.87,
                 16:0.74, 20:0.85, 30:0.77, 40:0.75, 50:0.86, 60:0.85, 70:0.87, 80:0.77 };
  const e = D.normDB['D-KEFS Twenty Questions Test · All Ages']['Initial Abstraction Score'];
  const bad = [];
  Object.entries(T215).forEach(([age, rxx]) => {
    if (e.rInternalByAge[age] !== rxx) {
      bad.push('age ' + age + ': stored ' + e.rInternalByAge[age] + ', Table 2.15 prints ' + rxx);
    }
  });
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('...and Table 2.17 still fails to reproduce from it', () => {
  /* Pins the DEFECT, not the fix. Two ways this check earns its place:

     If someone "corrects" a stored coefficient to make the SEM table work —
     .75 to .66 at 40-49 is the tempting one — the count drops and this fails.

     If a corrected printing ever makes the two tables agree, this also fails,
     which is the signal to re-read the note in data.js and reconsider. A check
     that quietly passed either way would record nothing.

     Contrast is the point: the sibling column in the same table reproduces 16
     of 16 through the check above, so this is a property of one column. */
  const e = D.normDB['D-KEFS Twenty Questions Test · All Ages']['Initial Abstraction Score'];
  let off = 0;
  Object.entries(TWENTY_Q_IA_SEM).forEach(([age, sem]) => {
    if (Math.abs((1 - Math.pow(sem / 3, 2)) - e.rInternalByAge[age]) > 0.007) off++;
  });
  if (off !== 12) {
    return 'Table 2.17 now disagrees in ' + off + ' of 16 bands, not the documented 12 — '
         + 'either a coefficient was edited to suit it, or the source changed; see data.js';
  }
  /* The worked example the documentation quotes, so the note cannot go stale. */
  const at40 = 3 * Math.sqrt(1 - e.rInternalByAge[40]);
  if (Math.abs(at40 - 1.50) > 0.005) return 'the 40-49 example no longer gives 1.50, it gives ' + at40.toFixed(3);
  return Math.abs(TWENTY_Q_IA_SEM[40] - 1.76) < 1e-9
    || 'the 40-49 printed SEM is no longer 1.76';
});

check('the retest r could NOT have produced those SEMs', () => {
  /* Falsification, and the reason this is not merely a plausible reading of
     the manual's prose. If the by-age SEM tables were built on stability
     instead, the stored retest r would reproduce them. It reproduces 4 of 168
     -- and only by coincidence, at bands where the two coefficients happen to
     coincide. Anything above a handful means the wrong coefficient is stored. */
  let hits = 0, n = 0;
  DKEFS_SEM.forEach(([group, name, tab]) => {
    const e = D.normDB[group][name];
    Object.values(tab).forEach((sem) => {
      n++;
      if (Math.abs((1 - Math.pow(sem / 3, 2)) - e.r) <= 0.007) hits++;
    });
  });
  return hits <= 8 || 'the retest r reproduces ' + hits + ' of ' + n
    + ' published SEMs — the by-age tables may be stability after all, re-read Chapter 2';
});

check('D-KEFS original carries no bare rInternal, because none is published', () => {
  /* Unlike D-KEFS Advanced, not one of the eight coefficient tables prints an
     all-ages average. There is therefore nothing citable to put in rInternal,
     and averaging the bands would be inventing a number. The consequence is
     deliberate: a blank age falls through to the retest r, which IS the
     manual's own second regime — see the next check. */
  const bad = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    if (!/^D-KEFS (?!Advanced)/.test(group)) return;
    Object.entries(tab).forEach(([name, e]) => {
      if (e && typeof e === 'object' && Number.isFinite(e.rInternal)) bad.push(group + ' / ' + name);
    });
  });
  return bad.length === 0
    || bad.slice(0, 3).join('; ') + ' — the D-KEFS manual publishes no all-ages coefficient to cite';
});

check('a blank age falls back to the retest r, and an age changes the interval', () => {
  /* Both halves matter. The fallback must be the retest r (the manual's
     total-sample regime), and a valid age must actually reach the renderer, or
     the 168 coefficients are inert. Tower is the sharpest probe: retest r .44
     against .84 at age 10. */
  const group = 'D-KEFS Tower Test · All Ages';
  const name = 'Total Achievement Score';
  const e = D.normDB[group][name];
  const row = { name, group, scoreType: 'scaled' };
  const bad = [];
  [null, undefined, NaN, 4, 90, 200].forEach((age) => {
    const rel = batteryRel(row, 'scaled', age);
    if (!rel) { bad.push('no interval at age ' + age); return; }
    if (rel.r !== e.r) bad.push('age ' + age + ' used r ' + rel.r + ', expected the retest ' + e.r);
  });
  const inBand = batteryRel(row, 'scaled', 10);
  if (inBand.r !== e.rInternalByAge[10]) bad.push('age 10 did not pick up the band coefficient');
  const young = batteryCi(10, row, '95', 10);
  const blank = batteryCi(10, row, '95', null);
  if (young === blank) bad.push('age 10 and a blank age both printed ' + young);
  return bad.length === 0 || bad.join('; ');
});

check('the measures the manual excludes carry nothing', () => {
  /* Four exclusions, four different reasons. Each would be an easy and silent
     mistake to make, so each is pinned.

     Colour-Word Interference  its only coefficient table (2.9) is for the
       Combined Colour Naming + Word Reading COMPOSITE, which is not a measure
       in this database. Attaching it to Colour Naming or Word Reading would
       give a condition the composite's reliability.
     Design Fluency  "Item interdependence precluded the use of internal
       consistency procedures, and reliability was investigated with test-retest
       procedures."
     Trail Making, the five conditions  Table 2.1 is the composite alone.

     Twenty Questions is NOT on this list — both its measures carry Table 2.15.
     Its Initial Abstraction column has its own problem, pinned separately
     above. */
  const bad = [];
  const carries = (group, name) => {
    const e = D.normDB[group] && D.normDB[group][name];
    return !!(e && (Number.isFinite(e.rInternal) || e.rInternalByAge));
  };
  ['Colour Naming', 'Word Reading', 'Inhibition', 'Inhibition/Switching'].forEach((n) => {
    if (carries('D-KEFS Colour-Word Interference · All Ages', n)) bad.push('CWIT ' + n);
  });
  ['Filled Dots', 'Empty Dots', 'Switching'].forEach((n) => {
    if (carries('D-KEFS Design Fluency · All Ages', n)) bad.push('Design Fluency ' + n);
  });
  ['Visual Scanning', 'Number Sequencing', 'Letter Sequencing', 'Switching', 'Motor Speed'].forEach((n) => {
    if (carries('D-KEFS Trail Making Test · All Ages', n)) bad.push('TMT ' + n);
  });
  return bad.length === 0
    || bad.slice(0, 4).join('; ') + ' — the manual publishes no usable internal-consistency coefficient for these';
});

check('the two D-KEFS manuals disagree, and both are honoured', () => {
  /* The trap this section exists to stop. D-KEFS Advanced rejects internal
     consistency for Trail Making and Verbal Fluency because they are SPEEDED;
     the original manual publishes it for all four Verbal Fluency measures and
     for the Trail Making composite, and rejects Design Fluency for ITEM
     INTERDEPENDENCE instead. Same two test names, opposite answers. Anyone
     "harmonising" the two families would break one of them, so assert the
     inversion directly. */
  const bad = [];
  const has = (g, n) => {
    const e = D.normDB[g] && D.normDB[g][n];
    return !!(e && (Number.isFinite(e.rInternal) || e.rInternalByAge));
  };
  if (has('D-KEFS Advanced Verbal Fluency · All Ages', 'Letter Fluency Total Correct')) {
    bad.push('Advanced Verbal Fluency gained a coefficient; that manual calls it speeded');
  }
  if (!has('D-KEFS Verbal Fluency · All Ages', 'Letter Fluency')) {
    bad.push('original Verbal Fluency lost its coefficient; Table 2.4 publishes one');
  }
  if (!has('D-KEFS Trail Making Test · All Ages', 'Combined Number + Letter')) {
    bad.push('the original TMT composite lost its coefficient; Table 2.1 publishes one');
  }
  if (has('D-KEFS Design Fluency · All Ages', 'Filled Dots')) {
    bad.push('Design Fluency gained one; its manual says item interdependence precluded it');
  }
  return bad.length === 0 || bad.join('; ');
});

check('the manual pairs its normative SD with the UNCORRECTED retest r', () => {
  /* PINNED SOURCE: D-KEFS Technical Manual Table 2.8, Design Fluency, All Ages
     — SEM 1.94 / 1.97 / 2.47 for Filled Dots, Empty Dots and Switching.

     WHY THIS CHECK EXISTS, AND IT IS NOT ABOUT DESIGN FLUENCY. Every measure
     in this app whose interval rests on a bare retest r pairs that r with the
     NORMATIVE SD of its metric — 3 for a D-KEFS scaled score. r describes the
     retest study's sample, the SD describes the norm group, and pairing terms
     from two populations normally makes SD sqrt(1 - rxx) an invalid error
     variance. That argument is sound, and the Allen & Yen / Magnusson
     correction that repairs it,

         rxx = 1 - (sd1^2 / normSD^2)(1 - r)

     reproduces the published rCorrected to a median error of .003 across the
     267 entries carrying both. Applying it in resolveCiReliability was built
     and tested; it moves 246 entries, 46 of them reachable from Score Tables
     — 25 D-KEFS and 21 D-KEFS Advanced, and nothing else. Both counts are
     asserted now: 246 by section 32, and the 46 with its cost breakdown by the
     check pinning the Methods & References wording.

     IT WAS REJECTED BECAUSE THE PUBLISHER GOT THERE FIRST. The manual states
     the pairing (p. 19: test-retest SEMs "were derived from the total sample
     of cases"; "The standard deviation unit is 3 for all D-KEFS scaled
     scores") and Table 2.8 is that arithmetic on the page. Correcting would
     print a coefficient the cited manual does not contain, on the one family
     where its own working is checkable — the same ground on which the
     unrounded Fisher's-z WAIS-IV average was declined.

     It is OFFERED, though — the reliability-basis control on Score Tables
     applies it deliberately, with the APA note saying so. What must never
     happen is it becoming the default, or leaking into a caller that did not
     ask, so this drives the SHIPPED renderer in both states.

     Four assertions, and the middle two are the load-bearing ones:
       1. the uncorrected r reproduces Table 2.8 exactly;
       2. the DEFAULT reading — no flag passed at all — still gives it, so the
          published interval is what an untouched app prints;
       3. the corrected reading does NOT reproduce Table 2.8, so the two
          readings genuinely separate and the choice is a real one;
       4. a missing control reads as published, so a harness or a page whose
          control failed to build cannot land on derived coefficients.

     If a corrected printing of the manual ever makes the corrected value fit,
     3 fails, which is the signal to re-read this note rather than assume a
     regression. */
  const T28 = { 'Filled Dots': 1.94, 'Empty Dots': 1.97, 'Switching': 2.47 };
  const group = 'D-KEFS Design Fluency · All Ages';
  const bad = [];
  let corrHits = 0;
  Object.entries(T28).forEach(([name, sem]) => {
    const e = D.normDB[group][name];
    const row = { name, group, scoreType: 'scaled' };
    const at = (rel) => Math.round(3 * Math.sqrt(1 - rel.r) * 100) / 100;

    const plain = Math.round(3 * Math.sqrt(1 - e.r) * 100) / 100;
    if (plain !== sem) bad.push(name + ': 3 sqrt(1 - r) = ' + plain + ', Table 2.8 prints ' + sem);

    // 2. the default, with no flag passed — what an untouched app prints
    const dflt = batteryRel(row, 'scaled', undefined);
    if (at(dflt) !== sem) {
      bad.push(name + ': the DEFAULT reading gives SEM ' + at(dflt) + ', Table 2.8 prints ' + sem
        + ' — the correction has become the default');
    }
    if (dflt.basis !== 'retest, uncorrected') bad.push(name + ': default basis is "' + dflt.basis + '"');
    if (batteryRel(row, 'scaled', undefined, false).r !== dflt.r) {
      bad.push(name + ': explicitly asking for published differs from the default');
    }

    // 3. and the corrected reading must be a DIFFERENT number
    const corr = batteryRel(row, 'scaled', undefined, true);
    if (corr.basis !== 'retest, corrected here') bad.push(name + ': corrected basis is "' + corr.basis + '"');
    if (at(corr) === sem) corrHits++;
  });
  if (corrHits) {
    bad.push(corrHits + ' of 3 also reproduce from the range-corrected coefficient — '
      + 'the two readings no longer separate, so the argument in this check needs re-checking');
  }
  /* 4. The control's own default, read straight from the markup rather than
        from the function, because that is the value a browser actually loads. */
  if (!/id="bat-ci-basis"[^>]*value="published"/.test(HTML_SRC)) {
    bad.push('#bat-ci-basis does not default to published in index.html');
  }
  if (batteryCiCorrectRetest() !== false) {
    bad.push('batteryCiCorrectRetest returns true with no control present; a missing control must read as published');
  }
  /* A published coefficient must never be overwritten, whatever the toggle
     says. Cheapest proof: a measure with internal consistency is unmoved. */
  const wais = { name: 'Vocabulary', group: 'WAIS-IV Core Subtests · All Ages', scoreType: 'scaled' };
  if (batteryRel(wais, 'scaled', undefined, true).r !== batteryRel(wais, 'scaled', undefined, false).r) {
    bad.push('the toggle moved a published internal-consistency coefficient; it must reach only the retest branch');
  }
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   28. WAIS-IV internal consistency — Technical Manual Table 4.1

   PINNED SOURCE: WAIS-IV Technical and Interpretive Manual (GB), chapter 4
   "Reliability and Errors of Measurement", pp. 42-43.
     p. 42     "Reliability coefficients were obtained utilizing the split-half
               and the Cronbach's coefficient alpha methods... Because Symbol
               Search, Coding, and Cancellation are Processing Speed subtests,
               the split-half coefficient is not a proper reliability estimate.
               Therefore, test-retest stability coefficients were used as the
               reliability estimates for these subtests... corrected for the
               normative sample's variability (Allen & Yen, 1979; Magnusson,
               1967)."
     Table 4.1 24 measures x 13 NORMATIVE age bands (16-17, 18-19, 20-24,
               25-29, 30-34, 35-44, 45-54, 55-64, 65-69, 70-74, 75-79, 80-84,
               85-90) plus an overall average via Fisher's z.
     Table 4.3 the SEMs those coefficients produce, same 24 x 13 shape, plus an
               overall average SEM. 300 published cells; LN, FW and CA print a
               dash above 65-69, being normed only that far.
     Table 4.3 note. "The reliability coefficients shown in Table 4.1 and the
               population standard deviations (i.e., 3 for the subtests and 15
               for the composite scores) were used to compute the SEMs."
               a "The average SEMs were calculated by averaging the squared
               SEMs for each age group and obtaining the square root of the
               result."

   Table 4.3 was obtained after the coefficients were stored, and it closes what
   this section previously recorded as unproven: that the manual builds its own
   SEMs from Table 4.1. It is now arithmetic, not the publisher's word for it —
   all 300 banded cells below, exactly, to the printed 2dp.

   It settles the SD question too, and independently of section 23. That section
   derived "pair the coefficient with the NORMATIVE SD, not sd1" from the CVLT-3
   manual's SEM column; this note states the same rule outright for a different
   publisher, and the 300 cells only reproduce that way. The pairing 6de0d81
   replaced — sd1 against a normatively-scaled coefficient — reproduces 29.
   ========================================================================== */
heading('28. WAIS-IV internal consistency — Table 4.1');

const WAIS_ALL = {
  core: 'WAIS-IV Core Subtests · All Ages',
  supp: 'WAIS-IV Supplementary Subtests · All Ages',
  proc: 'WAIS-IV Process Scores · All Ages',
  idx:  'WAIS-IV Indices · All Ages'
};

/* Table 4.1, the three SPEEDED rows only, as printed: one value per normative
   band. These are stability coefficients, not internal consistency. */
const WAIS_SPEEDED_T41 = [
  [WAIS_ALL.core, 'Symbol Search', 0.81,
    { 16:0.81, 18:0.81, 20:0.81, 25:0.81, 30:0.73, 35:0.73, 45:0.73, 55:0.81, 65:0.81, 70:0.86, 75:0.86, 80:0.86, 85:0.86 }],
  [WAIS_ALL.core, 'Coding', 0.86,
    { 16:0.85, 18:0.85, 20:0.85, 25:0.85, 30:0.84, 35:0.84, 45:0.84, 55:0.89, 65:0.89, 70:0.86, 75:0.86, 80:0.86, 85:0.86 }],
  [WAIS_ALL.supp, 'Cancellation', 0.78,
    { 16:0.81, 18:0.81, 20:0.81, 25:0.81, 30:0.71, 35:0.71, 45:0.71, 55:0.80, 65:0.80 }]
];

/* Normative band lower bound -> the retest study band this database uses. */
const WAIS_RETEST_BAND = { 16:'16-29', 18:'16-29', 20:'16-29', 25:'16-29',
                           30:'30-54', 35:'30-54', 45:'30-54',
                           55:'55-69', 65:'55-69',
                           70:'70-90', 75:'70-90', 80:'70-90', 85:'70-90' };

check('the speeded three already hold Table 4.1, in rCorrected', () => {
  /* THE LOAD-BEARING CHECK, and the reason those three need no new field.
     The manual gives them a corrected stability coefficient, which is exactly
     what rCorrected stores. Table 4.1 repeats one retest-study value across
     every normative band it spans, so all 35 published cells must reproduce
     from the banded groups — and the overall average from All Ages.

     This is also the transcription guard for the whole section: the band
     mapping below is what turns Table 4.1's 13 normative bands into this
     database's four, and if it were wrong these cells would not line up. */
  const bad = [];
  let n = 0;
  WAIS_SPEEDED_T41.forEach(([allGroup, name, overall, byBand]) => {
    const base = allGroup.replace(' · All Ages', '');
    Object.entries(byBand).forEach(([lo, rxx]) => {
      const g = base + ' · Ages ' + WAIS_RETEST_BAND[lo];
      const e = D.normDB[g] && D.normDB[g][name];
      if (!e) { bad.push('no entry at ' + g + ' / ' + name); return; }
      n++;
      if (e.rCorrected !== rxx) {
        bad.push(name + ' band ' + lo + ' [' + g + ']: rCorrected ' + e.rCorrected + ', Table 4.1 prints ' + rxx);
      }
    });
    n++;
    const eA = D.normDB[allGroup][name];
    if (eA.rCorrected !== overall) bad.push(name + ' All Ages: rCorrected ' + eA.rCorrected + ', Table 4.1 prints ' + overall);
  });
  if (n !== 38) bad.push('expected 38 published cells (35 banded + 3 averages), tested ' + n);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the raw retest r could NOT have produced those cells', () => {
  /* Falsification, and the reason the match above is not a coincidence of two
     similar numbers. rCorrected hits 35 of 35; the uncorrected r hits none.
     If this ever starts hitting, the wrong field is being compared. */
  let hits = 0, n = 0;
  WAIS_SPEEDED_T41.forEach(([allGroup, name, overall, byBand]) => {
    const base = allGroup.replace(' · All Ages', '');
    Object.entries(byBand).forEach(([lo, rxx]) => {
      const e = D.normDB[base + ' · Ages ' + WAIS_RETEST_BAND[lo]][name];
      n++; if (e.r === rxx) hits++;
    });
    n++; if (D.normDB[allGroup][name].r === overall) hits++;
  });
  return hits <= 2 || 'the uncorrected r reproduces ' + hits + ' of ' + n
    + ' cells — check which field Table 4.1 actually holds';
});

check('the speeded three carry no internal-consistency field', () => {
  /* Their Table 4.1 value is a stability coefficient. Storing it in rInternal
     would change no number but would label it internal consistency, which the
     Methods & References paragraph then asserts on screen. D-KEFS Advanced
     excludes TMT and VFT for the same reason — see section 24. */
  const bad = [];
  WAIS_SPEEDED_T41.forEach(([group, name]) => {
    const e = D.normDB[group][name];
    if (Number.isFinite(e.rInternal) || e.rInternalByAge) bad.push(group + ' / ' + name);
  });
  return bad.length === 0
    || bad.join('; ') + ' — the WAIS-IV manual calls the split-half coefficient improper for these';
});

check('spot values match Table 4.1 as printed', () => {
  /* Transcription guard over the 21 stored measures. Drawn from all four
     families and both ends of the age range, so a column shift would surface. */
  const want = [
    [WAIS_ALL.core, 'Block Design',               16, 0.88],
    [WAIS_ALL.core, 'Vocabulary',                 85, 0.96],
    [WAIS_ALL.core, 'Digit Span',                 25, 0.94],
    [WAIS_ALL.supp, 'Letter-Number Sequencing',   65, 0.88],
    [WAIS_ALL.proc, 'Digit Span Sequencing',      85, 0.92],
    [WAIS_ALL.idx,  'Perceptual Reasoning Index', 80, 0.92],
    [WAIS_ALL.idx,  'Full Scale IQ',              16, 0.97]
  ];
  const bad = [];
  want.forEach(([g, n, age, rxx]) => {
    const got = D.normDB[g][n].rInternalByAge[age];
    if (got !== rxx) bad.push(n + ' band ' + age + ': stored ' + got + ', Table 4.1 prints ' + rxx);
  });
  return bad.length === 0 || bad.join('; ');
});

/* Table 4.3 as printed. Rows are [family, population SD, normDB name, the 13
   normative bands in order, overall average SEM]; null is the dash the manual
   prints above 65-69 for the three measures normed only that far. Read off the
   supplied table rather than retyped. */
const WAIS_T43_BANDS = [16, 18, 20, 25, 30, 35, 45, 55, 65, 70, 75, 80, 85];
const WAIS_T43 = [
  ['core', 3, 'Block Design', [1.04, 1.08, 1.2, 0.95, 0.9, 0.99, 0.95, 1.04, 1.08, 0.99, 1.27, 1.34, 1.12], 1.08],
  ['core', 3, 'Similarities', [1.31, 1.16, 1.16, 1.12, 1.08, 1.04, 1.04, 1.08, 1.04, 0.95, 1.12, 0.9, 0.9], 1.07],
  ['core', 3, 'Digit Span', [0.99, 0.85, 0.9, 0.73, 0.73, 0.73, 0.73, 0.85, 0.79, 0.73, 0.79, 0.85, 0.85], 0.81],
  ['core', 3, 'Matrix Reasoning', [1.04, 1.08, 1.04, 0.9, 0.9, 0.95, 0.95, 0.95, 0.9, 0.9, 0.73, 1.12, 0.85], 0.95],
  ['core', 3, 'Vocabulary', [0.79, 0.79, 0.73, 0.79, 0.79, 0.73, 0.73, 0.73, 0.67, 0.67, 0.73, 0.73, 0.6], 0.73],
  ['core', 3, 'Arithmetic', [0.99, 1.04, 1.2, 0.99, 0.95, 0.99, 0.9, 1.04, 0.99, 1.2, 0.95, 0.99, 1.12], 1.03],
  ['core', 3, 'Symbol Search', [1.31, 1.31, 1.31, 1.31, 1.56, 1.56, 1.56, 1.31, 1.31, 1.12, 1.12, 1.12, 1.12], 1.32],
  ['core', 3, 'Visual Puzzles', [0.95, 0.99, 0.95, 0.9, 0.95, 1.04, 0.85, 0.99, 0.85, 0.99, 0.99, 1.27, 1.41], 1.02],
  ['core', 3, 'Information', [0.99, 0.9, 0.9, 0.9, 0.9, 0.85, 0.73, 0.67, 0.73, 0.73, 0.73, 0.73, 0.6], 0.8],
  ['core', 3, 'Coding', [1.16, 1.16, 1.16, 1.16, 1.2, 1.2, 1.2, 0.99, 0.99, 1.12, 1.12, 1.12, 1.12], 1.13],
  ['supp', 3, 'Letter-Number Sequencing', [0.95, 0.95, 1.16, 1.04, 0.9, 1.12, 1.04, 1.08, 1.04, null, null, null, null], 1.03],
  ['supp', 3, 'Figure Weights', [0.95, 0.85, 1.04, 0.9, 0.85, 0.99, 0.9, 0.99, 0.95, null, null, null, null], 0.94],
  ['supp', 3, 'Comprehension', [1.27, 1.08, 1.08, 1.04, 1.16, 0.99, 1.16, 1.08, 1.08, 1.27, 1.08, 1.08, 0.95], 1.11],
  ['supp', 3, 'Cancellation', [1.31, 1.31, 1.31, 1.31, 1.62, 1.62, 1.62, 1.34, 1.34, null, null, null, null], 1.43],
  ['supp', 3, 'Picture Completion', [1.34, 1.2, 1.27, 1.27, 1.12, 1.2, 1.24, 1.12, 1.16, 0.99, 1.12, 1.24, 1.27], 1.2],
  ['proc', 3, 'Block Design No Time Bonus', [1.08, 1.12, 1.31, 1.04, 0.99, 1.08, 1.08, 1.12, 1.16, 1.04, 1.27, 1.34, 1.12], 1.14],
  ['proc', 3, 'Digit Span Forward', [1.44, 1.47, 1.34, 1.16, 1.44, 1.2, 1.24, 1.44, 1.37, 1.04, 1.31, 1.31, 1.41], 1.33],
  ['proc', 3, 'Digit Span Backward', [1.37, 1.34, 1.34, 1.2, 1.2, 1.12, 1.12, 1.27, 1.41, 1.37, 1.34, 1.44, 1.27], 1.3],
  ['proc', 3, 'Digit Span Sequencing', [1.56, 1.37, 1.27, 1.31, 1.34, 1.24, 1.27, 1.37, 1.31, 1.12, 1.12, 1.2, 0.85], 1.27],
  ['idx', 15, 'Verbal Comprehension Index', [3.67, 3, 3, 3, 3, 3, 2.6, 2.6, 2.6, 2.6, 3, 2.6, 2.12], 2.85],
  ['idx', 15, 'Perceptual Reasoning Index', [3.35, 3.67, 3.67, 3, 3, 3.35, 3, 3.35, 3.35, 3.35, 3.67, 4.24, 3.97], 3.48],
  ['idx', 15, 'Working Memory Index', [3.97, 3.67, 4.24, 3.35, 3.35, 3.35, 3.35, 3.67, 3.67, 3.97, 3.35, 3.67, 3.97], 3.67],
  ['idx', 15, 'Processing Speed Index', [5.2, 4.74, 4.74, 4.74, 5.41, 5.41, 5.41, 4.5, 4.5, 4.5, 4.24, 4.24, 4.24], 4.78],
  ['idx', 15, 'Full Scale IQ', [2.6, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12, 2.12], 2.16]
];

/* The coefficient the app would use for a given measure at a given band: the
   stored internal-consistency lookup, except on the speeded three, where the
   manual's Table 4.1 value is the corrected stability coefficient rCorrected
   already holds. Deliberately reads the same fields the renderer reads. */
function waisCoefficientAt(fam, name, lo) {
  const all = D.normDB[WAIS_ALL[fam]][name];
  if (all.rInternalByAge && Number.isFinite(all.rInternalByAge[lo])) return all.rInternalByAge[lo];
  const banded = D.normDB[WAIS_ALL[fam].replace('· All Ages', '· Ages ' + WAIS_RETEST_BAND[lo])];
  return banded && banded[name] ? banded[name].rCorrected : null;
}
const round2 = (v) => Math.round(v * 100) / 100;

check('Table 4.3 reproduces from the stored coefficients, all 300 cells', () => {
  /* THE CLAIM THIS SECTION USED TO RECORD AS UNPROVEN. The manual's own SEM
     table is now arithmetic rather than the publisher's word: every printed
     cell equals populationSD * sqrt(1 - the coefficient this app stores),
     to the printed 2dp, with nothing to spare and no tolerance allowed.

     It covers the speeded three as well, which is a second and independent
     confirmation of the 38/38 rCorrected finding above — this time through a
     different table. */
  const bad = [];
  let n = 0;
  WAIS_T43.forEach(([fam, sd, name, cells]) => {
    cells.forEach((printed, i) => {
      if (printed === null) return;
      const lo = WAIS_T43_BANDS[i];
      const rxx = waisCoefficientAt(fam, name, lo);
      if (rxx == null) { bad.push(name + ' band ' + lo + ': no coefficient stored'); return; }
      n++;
      const der = round2(sd * Math.sqrt(1 - rxx));
      if (der !== printed) bad.push(name + ' band ' + lo + ': rxx ' + rxx + ' gives ' + der + ', Table 4.3 prints ' + printed);
    });
  });
  if (n !== 300) bad.push('expected 300 published cells, tested ' + n);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('no other SD or coefficient could have produced Table 4.3', () => {
  /* Falsification, four ways, because a table of SEMs around 1.0 is exactly
     the kind of thing that half-fits several models.

     The sd1 row is the one that matters beyond this section: it is the pairing
     6de0d81 replaced, and Table 4.3's note independently states the rule
     section 23 derived from a different publisher — "the population standard
     deviations (i.e., 3 for the subtests and 15 for the composite scores) were
     used". If that row ever climbs, the SD the renderer picks has drifted. */
  const tally = { rawR: 0, corrected: 0, sd1: 0 };
  WAIS_T43.forEach(([fam, sd, name, cells]) => {
    cells.forEach((printed, i) => {
      if (printed === null) return;
      const lo = WAIS_T43_BANDS[i];
      const rxx = waisCoefficientAt(fam, name, lo);
      const g = D.normDB[WAIS_ALL[fam].replace('· All Ages', '· Ages ' + WAIS_RETEST_BAND[lo])];
      const e = g && g[name];
      if (!e) return;
      if (Number.isFinite(e.r) && round2(sd * Math.sqrt(1 - e.r)) === printed) tally.rawR++;
      if (Number.isFinite(e.rCorrected) && round2(sd * Math.sqrt(1 - e.rCorrected)) === printed) tally.corrected++;
      if (Number.isFinite(e.sd1) && rxx != null && round2(e.sd1 * Math.sqrt(1 - rxx)) === printed) tally.sd1++;
    });
  });
  const bad = [];
  if (tally.rawR > 10) bad.push('the retest r reproduces ' + tally.rawR + ' of 300 — it should reproduce almost none');
  // rCorrected legitimately accounts for the 38 speeded cells; well above that
  // would mean the two fields have converged and the pin has stopped biting.
  if (tally.corrected > 60) bad.push('rCorrected reproduces ' + tally.corrected + ' of 300, beyond the 38 speeded cells it owns');
  if (tally.sd1 > 45) bad.push('sd1 reproduces ' + tally.sd1 + ' of 300 — the retest SD should not fit the manual\'s SEMs');
  return bad.length === 0 || bad.join('; ');
});

check('the manual averages its two tables by two different methods', () => {
  /* THE REASON NO STORAGE CHOICE RECONCILES TABLES 4.1 AND 4.3, and the thing
     to read before anyone tries to make them agree.

     Table 4.3 footnote a: "The average SEMs were calculated by averaging the
     squared SEMs for each age group and obtaining the square root of the
     result." Expand that. RMS of the band SEMs
       = sqrt( mean( SD^2 (1 - r_i) ) )
       = SD * sqrt( 1 - mean(r_i) )
     — algebraically identical to running the ARITHMETIC MEAN of the band
     coefficients through the SEM formula.

     But Table 4.1's average column is a FISHER-Z average; the check below
     reproduces all 24 stored averages that way and only 18 by the arithmetic
     mean. So the two tables average the same 13 numbers by different methods,
     and Table 4.3's average CANNOT be derived from Table 4.1's average. This
     is a property of the source, not a defect in this app, and no choice of
     what to store makes both columns come out right.

     Pinned as a fact, the way section 27 pins the Twenty Questions
     contradiction rather than smoothing it. The discrimination is decisive —
     20 of 24 against 5 of 24 — so a check that merely passed would record
     nothing. If a corrected printing ever makes the two agree, this fails,
     which is the signal to re-read the decision in data.js. */
  const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
  const tanh = (z) => (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  const bad = [];
  let arith = 0, fisher = 0, stored = 0, n = 0;
  WAIS_T43.forEach(([fam, sd, name, cells, avg]) => {
    n++;
    const bands = [];
    cells.forEach((printed, i) => {
      if (printed === null) return;
      const rxx = waisCoefficientAt(fam, name, WAIS_T43_BANDS[i]);
      if (rxx != null) bands.push(rxx);
    });
    const mean = bands.reduce((a, r) => a + r, 0) / bands.length;
    const fz = tanh(bands.reduce((a, r) => a + atanh(r), 0) / bands.length);
    const all = D.normDB[WAIS_ALL[fam]][name];
    const overall = Number.isFinite(all.rInternal) ? all.rInternal : all.rCorrected;
    // Tolerance only on the arithmetic route, because that is the one claimed
    // to be exact: the residue is the manual computing from unrounded
    // coefficients, damped by sqrt(13) across bands. Max observed 0.0067.
    if (Math.abs(sd * Math.sqrt(1 - mean) - avg) <= 0.007) arith++;
    if (round2(sd * Math.sqrt(1 - fz)) === avg) fisher++;
    if (round2(sd * Math.sqrt(1 - overall)) === avg) stored++;
  });
  if (arith !== n) bad.push('the arithmetic-mean (RMS) reading fits only ' + arith + ' of ' + n + ' — re-read footnote a');
  if (fisher > 10) bad.push('the Fisher-z reading now fits ' + fisher + ' of ' + n + ' — the two tables have converged, and the decision recorded in data.js needs re-reading');
  if (stored > 10) bad.push('the stored average now reproduces ' + stored + ' of ' + n + ' — the divergence this section documents has gone');
  return bad.length === 0 || bad.join('; ');
});

check('the published overall average is the Fisher-z average of its own bands', () => {
  /* THE BLANK-AGE FALLBACK, WHICH NOTHING PREVIOUSLY PINNED. rInternal is what
     a patient with no age entered is scored on, and until Table 4.3 arrived
     there was no way to hold it: the spot check above only reads
     rInternalByAge, and "an average is present" only tests that it is a number.
     Verified by mutation — moving Vocabulary's average .94 -> .95 passed all
     238 checks.

     Table 4.3 fixes the 13 band coefficients exactly, and Table 4.1's average
     column is computed with Fisher's z, so the average is now derivable from
     values this file already pins rather than taken on trust.

     Note WHICH average, because it discriminates: Fisher's z reproduces all 24
     stored values to the printed 2dp, the arithmetic mean only 18. That 18 is
     also why the sibling check below refuses to assert the average DIFFERS
     from the mean of its bands — over coefficients this close together the two
     agree more often than not. */
  const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
  const tanh = (z) => (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  const bad = [];
  let n = 0, arith = 0;
  WAIS_T43.forEach(([fam, , name]) => {
    const all = D.normDB[WAIS_ALL[fam]][name];
    const stored = Number.isFinite(all.rInternal) ? all.rInternal : all.rCorrected;
    const bands = [];
    WAIS_T43_BANDS.forEach((lo) => {
      const rxx = waisCoefficientAt(fam, name, lo);
      // The speeded three repeat one retest value per span, so de-duplicating
      // is wrong: Table 4.1 prints a value in every band and averages all 13.
      if (rxx != null && !(all.rInternalAgeMax && lo > all.rInternalAgeMax)) bands.push(rxx);
    });
    if (!bands.length) { bad.push(name + ': no band coefficients to average'); return; }
    n++;
    const fisher = round2(tanh(bands.reduce((a, r) => a + atanh(r), 0) / bands.length));
    if (fisher !== stored) bad.push(name + ': Fisher-z average of its bands is ' + fisher + ', stored ' + stored);
    if (round2(bands.reduce((a, r) => a + r, 0) / bands.length) === stored) arith++;
  });
  if (n !== 24) bad.push('expected 24 measures, averaged ' + n);
  if (arith > 20) bad.push('the arithmetic mean now fits ' + arith + ' of ' + n + ' — the two readings have converged, and this check no longer shows the manual uses Fisher\'s z');
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the blank-age fallback departs from Table 4.3 in exactly two printed margins', () => {
  /* What the decision costs, measured the way section 27 insists this kind of
     discrepancy is judged: on the margin actually printed, not on the gap
     between the two statistics.

     Over 24 measures x both offered CI levels, the app's blank-age fallback and
     the manual's average SEM round to the same half-width in 46 of 48 cases.
     The two that differ are knife-edge:
       FSIQ at 90%   1.645 * 2.1213 = 3.49 -> +/-3;  1.645 * 2.16 = 3.55 -> +/-4
       DSB  at 95%   1.960 * 1.2728 = 2.49 -> +/-2;  1.960 * 1.30 = 2.55 -> +/-3

     WHAT ACTUALLY DRIVES THAT GAP, because the obvious answer is wrong. It is
     not mainly the averaging method. Decomposing FSIQ:
       arithmetic mean of the bands  .97923  -> 2.1617   (= Table 4.3's 2.16)
       Fisher-z average, unrounded   .97936  -> 2.1547
       the stored published average  .98     -> 2.1213   (what the app prints)
     Method costs 0.007; ROUNDING THE PUBLISHED AVERAGE TO 2DP costs 0.033,
     about five times more. On VCI it is 0.044 against 0.19.

     THE THIRD OPTION, WEIGHED AND DECLINED — do not "discover" it and switch
     without reading this. Using the Fisher-z average of the bands UNROUNDED
     would match the manual's printed margin on 24 of 24, both CI levels, and
     invents nothing: the bands are pinned above and Fisher's z is Table 4.1's
     own stated method. It was declined because the app renders the coefficient
     it actually used, so the reliability cell would print .979 where the manual
     prints .98 — a clinician cross-checking Table 4.1 would find a number that
     is not there, and that display contract is load-bearing (see the r vs
     rCorrected note in CLAUDE.md). Reviewed and kept as-is, 2026-07-31.

     Entering an age — the primary path, and the one Table 4.3 verifies exactly
     — sidesteps the whole question.

     Pinned as a fact so it cannot drift unnoticed in either direction: if a
     third measure joins them the fallback needs re-arguing, and if these two
     ever agree the note above is stale. */
  const bad = [];
  const seen = [];
  WAIS_T43.forEach(([fam, sd, name, , avg]) => {
    const all = D.normDB[WAIS_ALL[fam]][name];
    const overall = Number.isFinite(all.rInternal) ? all.rInternal : all.rCorrected;
    const appSem = sd * Math.sqrt(1 - overall);
    [['90', 1.645], ['95', 1.96]].forEach(([label, z]) => {
      if (Math.round(z * appSem) !== Math.round(z * avg)) seen.push(name + ' @' + label);
    });
  });
  const want = ['Digit Span Backward @95', 'Full Scale IQ @90'];
  const got = seen.slice().sort();
  if (got.join(' | ') !== want.join(' | ')) {
    bad.push('expected exactly [' + want.join(', ') + '], got [' + got.join(', ') + ']');
  }
  return bad.length === 0 || bad.join('; ');
});

check('every stored band key is a Table 4.1 band, and an average is present', () => {
  /* Shape, over all 21 measures at once. The keys must be the 13 NORMATIVE
     lower bounds — not this database's retest bands, which is the mistake the
     "All Ages only" rule in section 24 exists to prevent — and rInternal must
     be a published overall average, never a mean of the bands. */
  const BANDS = [16, 18, 20, 25, 30, 35, 45, 55, 65, 70, 75, 80, 85];
  const bad = [];
  let seen = 0;
  Object.values(WAIS_ALL).forEach((group) => {
    Object.entries(D.normDB[group]).forEach(([name, e]) => {
      if (!e.rInternalByAge) return;
      seen++;
      const keys = Object.keys(e.rInternalByAge).map(Number).sort((a, b) => a - b);
      const stray = keys.filter((k) => !BANDS.includes(k));
      if (stray.length) bad.push(name + ' has non-Table-4.1 band(s) ' + stray.join(','));
      // Bands must be a PREFIX of the full list: the only gap the manual has is
      // a truncation at the top (LN, FW stop at 65-69), never a hole.
      if (keys.join() !== BANDS.slice(0, keys.length).join()) bad.push(name + ' bands are not a prefix of Table 4.1: ' + keys.join());
      if (!Number.isFinite(e.rInternal)) bad.push(name + ' has no published overall average');
      /* No "rInternal must differ from the mean of its bands" assertion here.
         It was tried and it fires on good data: Table 4.1's averages are
         computed with Fisher's z, but over 13 coefficients this close together
         the two agree to 2dp often enough that Block Design, Vocabulary,
         Picture Completion and the WMI all land on the arithmetic mean exactly.
         The guard against an invented average is the spot check above and the
         roster in section 24, not arithmetic. */
    });
  });
  if (seen !== 21) bad.push('expected 21 measures carrying the field, found ' + seen);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('LN and FW stop at 69, where the manual prints a dash', () => {
  /* Both are normed to 69, so Table 4.1 has no band above 65-69 for them.
     Without rInternalAgeMax 69 the lookup would take the greatest key <= age
     and silently score a 90-year-old on the 65-69 coefficient. The fallback is
     the published overall average, which is computed over ages 16-69 — the
     only sample there is, and therefore still citable. */
  const bad = [];
  [['Letter-Number Sequencing', 0.88], ['Figure Weights', 0.9]].forEach(([name, avg]) => {
    const e = D.normDB[WAIS_ALL.supp][name];
    if (e.rInternalAgeMax !== 69) bad.push(name + ' has rInternalAgeMax ' + e.rInternalAgeMax + ', expected 69');
    if (rInternalForAge(e, 70) !== null) bad.push(name + ' still returned a band at age 70');
    if (rInternalForAge(e, 65) !== 0.88 && name === 'Letter-Number Sequencing') bad.push('LN lost its 65-69 band');
    const rel = batteryRel({ name, group: WAIS_ALL.supp, scoreType: 'scaled' }, 'scaled', 80);
    if (!rel) bad.push(name + ' produced no interval at age 80');
    else if (rel.r !== avg) bad.push(name + ' at age 80 used r ' + rel.r + ', expected the average ' + avg);
  });
  // ...and a measure normed to 90 must NOT be truncated.
  const vc = D.normDB[WAIS_ALL.core].Vocabulary;
  if (vc.rInternalAgeMax !== 90) bad.push('Vocabulary was truncated to ' + vc.rInternalAgeMax);
  if (rInternalForAge(vc, 88) !== 0.96) bad.push('Vocabulary at 88 did not read the 85-90 band');
  return bad.length === 0 || bad.join('; ');
});

check('a blank age falls back to the published average, and an age moves the interval', () => {
  /* The optional-age contract, end to end through the shipped renderer. A
     blank age must never empty the column, and the fallback must be the
     manual's own overall average — not the retest r, which WAIS-IV publishes
     an alternative to. PRI is the sharpest probe: .95 average against .92 at
     80-84, on an SD-15 metric where a coefficient change is most visible. */
  const row = { name: 'Perceptual Reasoning Index', group: WAIS_ALL.idx, scoreType: 'standard' };
  const e = D.normDB[WAIS_ALL.idx]['Perceptual Reasoning Index'];
  const bad = [];
  [null, undefined, NaN, 4, 91, 200].forEach((age) => {
    const rel = batteryRel(row, 'standard', age);
    if (!rel) { bad.push('no interval at age ' + age); return; }
    if (rel.r !== e.rInternal) bad.push('age ' + age + ' used r ' + rel.r + ', expected the average ' + e.rInternal);
  });
  const banded = batteryRel(row, 'standard', 82);
  if (banded.r !== 0.92) bad.push('age 82 did not pick up the 80-84 band');
  const blankCi = batteryCi(100, row, '95', null);
  const oldCi   = batteryCi(100, row, '95', 82);
  if (!blankCi || !oldCi) bad.push('one of the intervals did not render');
  else if (blankCi === oldCi) bad.push('a blank age and age 82 both printed ' + blankCi);
  return bad.length === 0 || bad.join('; ');
});

check('Methods & References states the WAIS-IV position, both halves of it', () => {
  /* The instrument-level naming check in section 24 cannot catch this one:
     "WAIS-IV" was already in that paragraph before this change, in the
     OPPOSITE sense — it named the manual only as a user of test-retest
     coefficients. The roster grew and the naming check stayed green. So assert
     the substance: the paragraph must say WAIS-IV uses internal consistency,
     and must still name the three subtests for which it does not. */
  const para = methodsCiBlock();
  const bad = [];
  if (!/WAIS[-‑–]IV/.test(para)) bad.push('WAIS-IV is not named at all');
  if (!/Table 4\.1/.test(para)) bad.push('the paragraph does not cite Table 4.1');
  ['Symbol Search', 'Coding', 'Cancellation'].forEach((n) => {
    if (!para.includes(n)) bad.push(n + ' is no longer named as an exclusion');
  });
  return bad.length === 0 || bad.join('; ');
});


/* ==========================================================================
   29. RBANS Update — Tables 3.6 and 3.7

   PINNED SOURCE: RBANS Update Manual (Randolph 2012), p. 42.
     Table 3.6  "Reliability Coefficients of RBANS Subtest and Index Scores by
                Age Group". 14 measures x 9 normative bands (12-13, 14-15,
                16-19, 20-39, 40-49, 50-59, 60-69, 70-79, 80-89) plus an
                overall average.
                footnote a  "Reliability estimates based on test-retest."
                            — on Figure Copy, Semantic Fluency, Coding,
                            Story Recall and Figure Recall ONLY.
                footnote b  "Average reliability coefficients were calculated
                            using Fisher's z transformation."
     Table 3.7  the SEMs, same shape. Note: "SEMs expressed in scaled score
                (subtest) and standard score (index) units."
                footnote a  the average SEM is the RMS across bands.

   This is the fifth internal-consistency source and the first whose table is
   MIXED-BASIS within itself, which is why rStability* exists. It is also the
   first family that needed new · All Ages groups: every other RBANS group
   holds a retest study, and Score Tables was falling through to whichever was
   listed first — scoring an 80-year-old on 55 adolescents.
   ========================================================================== */
heading('29. RBANS Update — Tables 3.6 and 3.7');

/* [family, population SD, normDB name, footnote a?, Table 3.6 bands, its
   average, Table 3.7 bands, its average]. Read off the supplied tables. */
const RBANS_BANDS = [12, 14, 16, 20, 40, 50, 60, 70, 80];
const RBANS_ALL = { sub: 'RBANS Subtests · All Ages', idx: 'RBANS Indices · All Ages' };
const RBANS_T = [
  ['sub', 3, "List Learning", false,
    [0.91, 0.88, 0.8, 0.82, 0.88, 0.85, 0.8, 0.86, 0.84], 0.85,
    [0.9, 1.04, 1.34, 1.27, 1.04, 1.16, 1.34, 1.12, 1.2], 1.17],
  ['sub', 3, "Story Memory", false,
    [0.87, 0.55, 0.79, 0.71, 0.73, 0.82, 0.79, 0.8, 0.84], 0.78,
    [1.08, 2.01, 1.37, 1.62, 1.56, 1.27, 1.37, 1.34, 1.2], 1.45],
  ['sub', 3, "Figure Copy", true,
    [0.42, 0.42, 0.42, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54], 0.5,
    [2.28, 2.28, 2.28, 2.03, 2.03, 2.03, 2.03, 2.03, 2.03], 2.12],
  ['sub', 3, "Semantic Fluency", true,
    [0.65, 0.65, 0.65, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52], 0.57,
    [1.77, 1.77, 1.77, 2.08, 2.08, 2.08, 2.08, 2.08, 2.08], 1.98],
  ['sub', 3, "Digit Span", false,
    [0.71, 0.86, 0.85, 0.84, 0.83, 0.85, 0.76, 0.86, 0.83], 0.83,
    [1.62, 1.12, 1.16, 1.2, 1.24, 1.16, 1.47, 1.12, 1.24], 1.27],
  ['sub', 3, "Coding", true,
    [0.76, 0.76, 0.76, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83], 0.81,
    [1.47, 1.47, 1.47, 1.24, 1.24, 1.24, 1.24, 1.24, 1.24], 1.32],
  ['sub', 3, "Story Recall", true,
    [0.45, 0.45, 0.45, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72], 0.54,
    [2.22, 2.22, 2.22, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59], 1.82],
  ['sub', 3, "Figure Recall", true,
    [0.71, 0.71, 0.71, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55], 0.59,
    [1.62, 1.62, 1.62, 2.01, 2.01, 2.01, 2.01, 2.01, 2.01], 1.89],
  ['idx', 15, "Immediate Memory", false,
    [0.93, 0.81, 0.86, 0.84, 0.88, 0.89, 0.85, 0.89, 0.9], 0.88,
    [3.97, 6.54, 5.61, 6, 5.2, 4.97, 5.81, 4.97, 4.74], 5.36],
  ['idx', 15, "Visuospatial/Constructional", false,
    [0.64, 0.67, 0.53, 0.77, 0.81, 0.82, 0.84, 0.81, 0.78], 0.75,
    [9, 8.62, 10.28, 7.19, 6.54, 6.36, 6, 6.54, 7.04], 7.63],
  ['idx', 15, "Attention", false,
    [0.81, 0.82, 0.81, 0.84, 0.84, 0.85, 0.83, 0.88, 0.85], 0.84,
    [6.54, 6.36, 6.54, 6, 6, 5.81, 6.18, 5.2, 5.81], 6.06],
  ['idx', 15, "Language", false,
    [0.79, 0.8, 0.74, 0.75, 0.76, 0.87, 0.85, 0.81, 0.83], 0.8,
    [6.87, 6.71, 7.65, 7.5, 7.35, 5.41, 5.81, 6.54, 6.18], 6.71],
  ['idx', 15, "Delayed Memory", false,
    [0.85, 0.85, 0.84, 0.84, 0.83, 0.84, 0.85, 0.83, 0.81], 0.84,
    [5.81, 5.81, 6, 6, 6.18, 6, 5.81, 6.18, 6.54], 6.04],
  ['idx', 15, "Total Scale", false,
    [0.92, 0.91, 0.9, 0.92, 0.94, 0.95, 0.93, 0.93, 0.94], 0.93,
    [4.24, 4.5, 4.74, 4.24, 3.67, 3.35, 3.97, 3.97, 3.67], 4.06]
];
const rb2 = (v) => Math.round(v * 100) / 100;

check('Table 3.7 reproduces from the stored coefficients, all 126 cells', () => {
  /* The bar this project sets for setting the retest default aside: the
     publisher must report the coefficient AND derive its own published
     intervals from it. Table 3.7 is that derivation, and it reproduces
     exactly at the printed 2dp on the SDs the manual's own note states. */
  const bad = [];
  let n = 0;
  RBANS_T.forEach(([fam, sd, name, retest, coef, , sems]) => {
    const e = D.normDB[RBANS_ALL[fam]][name];
    const byAge = retest ? e.rStabilityByAge : e.rInternalByAge;
    if (!byAge) { bad.push(name + ' carries no by-age lookup'); return; }
    RBANS_BANDS.forEach((lo, i) => {
      n++;
      if (byAge[lo] !== coef[i]) { bad.push(name + ' @' + lo + ': stored ' + byAge[lo] + ', Table 3.6 prints ' + coef[i]); return; }
      const der = rb2(sd * Math.sqrt(1 - byAge[lo]));
      if (der !== sems[i]) bad.push(name + ' @' + lo + ': rxx ' + byAge[lo] + ' gives SEM ' + der + ', Table 3.7 prints ' + sems[i]);
    });
  });
  if (n !== 126) bad.push('expected 126 cells, tested ' + n);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the five footnote-a measures are stability, and are never called internal consistency', () => {
  /* THE REASON THERE ARE TWO FIELDS. Table 3.6 marks these five "estimates
     based on test-retest". Storing them in rInternal would change no number
     while labelling a stability coefficient as internal consistency, which
     Methods & References and the APA note then assert on screen — exactly the
     mislabelling WAIS-IV Table 4.1 avoids for its three speeded subtests, and
     D-KEFS Advanced for TMT and VFT.

     Asserted in both directions, so neither field can drift into the other. */
  const bad = [];
  RBANS_T.forEach(([fam, , name, retest]) => {
    const e = D.normDB[RBANS_ALL[fam]][name];
    if (retest) {
      if (!e.rStabilityByAge) bad.push(name + ' is footnote a but has no rStabilityByAge');
      if (Number.isFinite(e.rInternal) || e.rInternalByAge) bad.push(name + ' is footnote a yet carries an internal-consistency field');
    } else {
      if (!e.rInternalByAge) bad.push(name + ' is unmarked but has no rInternalByAge');
      if (Number.isFinite(e.rStability) || e.rStabilityByAge) bad.push(name + ' is unmarked yet carries a stability field');
    }
  });
  const stab = [];
  Object.entries(D.normDB).forEach(([g, tab]) => Object.entries(tab).forEach(([n, e]) => {
    if (e && typeof e === 'object' && (Number.isFinite(e.rStability) || e.rStabilityByAge)) stab.push(g + ' / ' + n);
  }));
  /* Seven carriers now, across two manuals: RBANS Update Table 3.6's five
     footnote-a subtests, WMS-IV Table 3.1's one footnote-b measure in each of
     its two batteries, and WISC-V's three Processing Speed subtests. All are
     stability coefficients printed inside an otherwise internal-consistency
     reliability table, which is exactly what the field is for. Sections 30 and
     31 pin the WMS-IV and WISC-V carriers. */
  const rbans = stab.filter((c) => /^RBANS Subtests · All Ages \//.test(c));
  const wms = stab.filter((c) => /^WMS-IV Subtests · Ages (16-69|65-90) \/ Verbal Paired Associates II - Word Recall$/.test(c));
  const wisc = stab.filter((c) => /^WISC-V (Subtests|Process Scores) · All Ages \/ (Coding|Symbol Search|Cancellation|Cancellation Random|Cancellation Structured)$/.test(c));
  if (rbans.length !== 5) bad.push('RBANS rStability carriers: ' + rbans.length + ', expected 5');
  if (wms.length !== 2) bad.push('WMS-IV rStability carriers: ' + wms.length + ', expected 2');
  /* 5 for WISC-V: the three speeded subtests the manual names, plus the two
     Cancellation process scores derived from the same speeded task. */
  if (wisc.length !== 5) bad.push('WISC-V rStability carriers: ' + wisc.length + ', expected 5');
  if (stab.length !== 12) bad.push('rStability appears on ' + stab.length + ' entries, expected 12 (RBANS 5 + WMS-IV 2 + WISC-V 5)');
  return bad.length === 0 || bad.join('; ');
});

check('the adult stability values are this database\'s own rCorrected, 5 of 5', () => {
  /* THE TRANSCRIPTION PROOF, and the same one the WAIS-IV speeded three give.
     Table 3.6's footnote-a rows are corrected stability coefficients, so for
     the adult bands they must equal what the retest group already stores. The
     raw r must NOT fit, or the comparison is measuring nothing.

     The adolescent bands match only 2 of 5. They are a different retest sample
     and Table 3.6 is stored as printed rather than reconciled to Table 3.8 —
     recorded here so the mismatch reads as known rather than as an error. */
  let corr = 0, raw = 0, adoCorr = 0;
  RBANS_T.forEach(([fam, , name, retest, coef]) => {
    if (!retest) return;
    const adult = D.normDB['RBANS Subtests · Ages 20-89'][name];
    const ado = D.normDB['RBANS Subtests · Ages 12-19'][name];
    if (adult.rCorrected === coef[3]) corr++;
    if (adult.r === coef[3]) raw++;
    if (ado.rCorrected === coef[0]) adoCorr++;
  });
  const bad = [];
  if (corr !== 5) bad.push('adult bands match stored rCorrected ' + corr + '/5, expected 5');
  if (raw > 1) bad.push('the raw retest r reproduces ' + raw + '/5 — the comparison has stopped discriminating');
  if (adoCorr !== 2) bad.push('adolescent bands match rCorrected ' + adoCorr + '/5, expected 2 — re-read the note in data.js');
  return bad.length === 0 || bad.join('; ');
});

check('the four raw subtests carry no coefficient and print no interval', () => {
  /* They appear nowhere in Table 3.6 — the manual publishes reliability for
     its eight SCALED subtests only, which is its own confirmation of the
     raw/scaled split section 18 asserts from the data shape.

     Their m1/sd1 are the adult retest descriptives, carried over so the row
     stays selectable and declares its metric. Nothing may derive from them:
     a raw row has no percentile, no classification, and now no interval. If
     an interval ever appears here it means a coefficient was invented. */
  const bad = [];
  RBANS_RAW_SUBTESTS.forEach((name) => {
    const e = D.normDB[RBANS_ALL.sub][name];
    if (!e) { bad.push(name + ' missing from the All Ages group'); return; }
    if (e.metric !== 'raw') bad.push(name + ' is not tagged raw');
    ['rInternal', 'rInternalByAge', 'rStability', 'rStabilityByAge', 'r', 'rCorrected'].forEach((f) => {
      if (e[f] != null) bad.push(name + ' carries ' + f + ', which Table 3.6 does not publish');
    });
    const row = { name, group: RBANS_ALL.sub, scoreType: 'raw' };
    [null, 8, 45, 88].forEach((age) => {
      if (batteryCi(e.m1, row, '95', age)) bad.push(name + ' printed an interval at age ' + age);
    });
  });
  if (RBANS_T.some(([, , n]) => RBANS_RAW_SUBTESTS.includes(n))) bad.push('a raw subtest has appeared in Table 3.6');
  return bad.length === 0 || bad.join('; ');
});

check('the age band drives the interval, and a blank age falls back to the average', () => {
  /* End to end through the shipped renderer, on both bases.

     Immediate Memory is the sharpest probe and the reason these groups exist:
     before them Score Tables read the Ages 12-19 retest group for everyone and
     printed 85-115 at an index of 100. Table 3.6 gives .90 at 80-89. */
  const bad = [];
  const im = { name: 'Immediate Memory', group: RBANS_ALL.idx, scoreType: 'standard' };
  const e = D.normDB[RBANS_ALL.idx]['Immediate Memory'];
  if (batteryRel(im, 'standard', 82).r !== 0.9) bad.push('Immediate Memory at 82 did not read the 80-89 band');
  if (batteryRel(im, 'standard', 12).r !== 0.93) bad.push('Immediate Memory at 12 did not read the 12-13 band');
  [null, undefined, NaN, 5, 95].forEach((age) => {
    const rel = batteryRel(im, 'standard', age);
    if (!rel) { bad.push('no interval at age ' + age); return; }
    if (rel.r !== e.rInternal) bad.push('age ' + age + ' used r ' + rel.r + ', expected the published average ' + e.rInternal);
  });
  // the same, on a stability measure, so the second field is proven live
  const cod = { name: 'Coding', group: RBANS_ALL.sub, scoreType: 'scaled' };
  if (batteryRel(cod, 'scaled', 14).r !== 0.76) bad.push('Coding at 14 did not read its adolescent band');
  if (batteryRel(cod, 'scaled', 45).r !== 0.83) bad.push('Coding at 45 did not read its adult band');
  if (batteryRel(cod, 'scaled', null).r !== 0.81) bad.push('Coding with a blank age did not fall back to its published average');
  // and the defect these groups fix: the old and new answers must differ
  const wasCi = batteryCi(100, { name: 'Immediate Memory', group: 'RBANS Indices · Ages 12-19', scoreType: 'standard' }, '95', 82);
  const nowCi = batteryCi(100, im, '95', 82);
  if (wasCi === nowCi) bad.push('the All Ages group prints the same interval as the adolescent retest group — the fix is inert');
  return bad.length === 0 || bad.join('; ');
});

check('the All Ages groups stay out of Change Analysis and the SD Index', () => {
  /* They hold no second testing, so a row loaded there could never compute.
     isSingleAdministrationFamily is what keeps them out; the retest groups
     must keep answering false, or reliable change loses RBANS entirely. */
  const ctx = { getMergedDB: () => D.normDB, normDB: D.normDB };
  vm.createContext(ctx);
  vm.runInContext(extractFn(APP_SRC, 'isSingleAdministrationFamily') + ';globalThis.__F = isSingleAdministrationFamily;', ctx);
  const f = ctx.__F;
  const bad = [];
  Object.values(RBANS_ALL).forEach((g) => { if (!f(g)) bad.push(g + ' is not recognised as single-administration'); });
  ['RBANS Subtests · Ages 12-19', 'RBANS Subtests · Ages 20-89', 'RBANS Indices · Ages 12-19', 'RBANS Indices · Ages 20-89']
    .forEach((g) => { if (f(g)) bad.push(g + ' was wrongly excluded from Change Analysis'); });
  return bad.length === 0 || bad.join('; ');
});

check('RBANS averages its two tables the way WAIS-IV does, by two different methods', () => {
  /* Footnote b on Table 3.6 says Fisher's z; footnote a on Table 3.7 says the
     RMS of the band SEMs, which is algebraically SD sqrt(1 - ARITHMETIC mean).
     The same inconsistency section 28 records for WAIS-IV Tables 4.1 and 4.3,
     in a second manual — so it is a habit of the publisher, not a one-off, and
     the stored average stays the published COEFFICIENT either way. */
  const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
  const tanh = (z) => (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  let fisher = 0, arith = 0, rms = 0;
  RBANS_T.forEach(([, sd, , , coef, cAvg, sems, sAvg]) => {
    if (rb2(tanh(coef.reduce((a, r) => a + atanh(r), 0) / coef.length)) === cAvg) fisher++;
    if (rb2(coef.reduce((a, r) => a + r, 0) / coef.length) === cAvg) arith++;
    if (rb2(Math.sqrt(sems.reduce((a, s) => a + s * s, 0) / sems.length)) === sAvg) rms++;
  });
  const bad = [];
  if (fisher < 12) bad.push('Fisher z reproduces only ' + fisher + '/14 of the Table 3.6 averages');
  if (fisher <= arith) bad.push('Fisher z (' + fisher + ') no longer beats the arithmetic mean (' + arith + ') — re-read footnote b');
  if (rms !== 14) bad.push('the RMS reading fits ' + rms + '/14 of the Table 3.7 averages');
  return bad.length === 0 || bad.join('; ');
});

check('the stored averages are Table 3.6\'s own average column', () => {
  /* The blank-age fallback, pinned to the printed value rather than derived —
     the manual prints this column, unlike D-KEFS original. */
  const bad = [];
  RBANS_T.forEach(([fam, , name, retest, , cAvg]) => {
    const e = D.normDB[RBANS_ALL[fam]][name];
    const got = retest ? e.rStability : e.rInternal;
    if (got !== cAvg) bad.push(name + ': stored ' + got + ', Table 3.6 prints ' + cAvg);
    const max = retest ? e.rStabilityAgeMax : e.rInternalAgeMax;
    if (max !== 89) bad.push(name + ' has age max ' + max + ', expected 89 (the table stops at 80-89)');
  });
  return bad.length === 0 || bad.join('; ');
});


/* ==========================================================================
   30. WMS-IV — Technical Manual Tables 3.1 and 3.3

   PINNED SOURCE: WMS-IV Technical and Interpretive Manual (GB), chapter 3
   "Evidence of Reliability", pp. 44-46.
     p. 44      "Internal consistency reliability coefficients were calculated
                utilizing the split-half and alpha methods. In addition,
                stability coefficients were used on subtests for which the
                internal consistencies were not appropriate."
                "Internal consistency measures are used for all subtests
                excluding the Verbal Paired Associates Word Recall Scaled
                Score and recognition memory measures."
     Table 3.1  reliability by age group, printed as two battery blocks —
                Adult 16-69 (9 bands) and Older Adult 65-90 (5 bands).
                footnote a  averages via Fisher's z
                footnote b  "Reliability is based on test-retest correlation."
                            — on VPA II Word Recall only.
     Table 3.3  the SEMs, same two-block shape.
                footnote a  the average SEM is the RMS across bands.

   THE RECOGNITION MEMORY MEASURES ARE EXCLUDED, AND MUST STAY EXCLUDED. The
   manual reports their reliability as a DECISION-CONSISTENCY percentage — the
   percent agreement of impaired/not-impaired classification at a 10th
   percentile cut, adopted because those tasks are cumulative percentages with
   skewed distributions and attenuated retest correlations. A percent agreement
   is not a correlation and cannot enter SEM = SD sqrt(1 - rxx). They appear
   neither in Table 3.1 nor in normDB, so there is nothing to remove; the check
   below exists so that adding one later cannot pass unnoticed.
   ========================================================================== */
heading('30. WMS-IV — Tables 3.1 and 3.3');

/* [family, battery, population SD, normDB name, footnote b?, Table 3.1 bands,
   its average, Table 3.3 bands, its average]. */
const WMS_BANDS = { adult: [16, 18, 20, 25, 30, 35, 45, 55, 65], older: [65, 70, 75, 80, 85] };
const WMS_GROUP = (fam, bat) => (fam === 'idx' ? 'WMS-IV Indices' : 'WMS-IV Subtests') +
  (bat === 'adult' ? ' · Ages 16-69' : ' · Ages 65-90');
const WMS_T = [
  ['sub', 'adult', 3, "Logical Memory I", false,
    [0.8, 0.84, 0.79, 0.87, 0.82, 0.86, 0.77, 0.81, 0.79], 0.82,
    [1.34, 1.2, 1.37, 1.08, 1.27, 1.12, 1.44, 1.31, 1.37], 1.28],
  ['sub', 'older', 3, "Logical Memory I", false,
    [0.83, 0.88, 0.87, 0.84, 0.88], 0.86,
    [1.24, 1.04, 1.08, 1.2, 1.04], 1.12],
  ['sub', 'adult', 3, "Logical Memory II", false,
    [0.8, 0.87, 0.85, 0.9, 0.81, 0.9, 0.87, 0.84, 0.8], 0.85,
    [1.34, 1.08, 1.16, 0.95, 1.31, 0.95, 1.08, 1.2, 1.34], 1.17],
  ['sub', 'older', 3, "Logical Memory II", false,
    [0.85, 0.87, 0.88, 0.87, 0.89], 0.87,
    [1.16, 1.08, 1.04, 1.08, 0.99], 1.07],
  ['sub', 'adult', 3, "Verbal Paired Associates I", false,
    [0.93, 0.93, 0.94, 0.95, 0.94, 0.94, 0.94, 0.94, 0.93], 0.94,
    [0.79, 0.79, 0.73, 0.67, 0.73, 0.73, 0.73, 0.73, 0.79], 0.74],
  ['sub', 'older', 3, "Verbal Paired Associates I", false,
    [0.94, 0.92, 0.92, 0.92, 0.93], 0.93,
    [0.73, 0.85, 0.85, 0.85, 0.79], 0.82],
  ['sub', 'adult', 3, "Verbal Paired Associates II", false,
    [0.87, 0.84, 0.84, 0.89, 0.82, 0.84, 0.85, 0.85, 0.82], 0.85,
    [1.08, 1.2, 1.2, 0.99, 1.27, 1.2, 1.16, 1.16, 1.27], 1.17],
  ['sub', 'older', 3, "Verbal Paired Associates II", false,
    [0.82, 0.71, 0.71, 0.78, 0.68], 0.74,
    [1.27, 1.62, 1.62, 1.41, 1.7], 1.53],
  ['sub', 'adult', 3, "Verbal Paired Associates II - Word Recall", true,
    [0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76], 0.76,
    [1.47, 1.47, 1.47, 1.47, 1.47, 1.47, 1.47, 1.47, 1.47], 1.47],
  ['sub', 'older', 3, "Verbal Paired Associates II - Word Recall", true,
    [0.76, 0.76, 0.76, 0.76, 0.76], 0.76,
    [1.47, 1.47, 1.47, 1.47, 1.47], 1.47],
  ['sub', 'adult', 3, "Designs I", false,
    [0.84, 0.9, 0.87, 0.83, 0.88, 0.85, 0.83, 0.82, 0.82], 0.85,
    [1.2, 0.95, 1.08, 1.24, 1.04, 1.16, 1.24, 1.27, 1.27], 1.17],
  ['sub', 'adult', 3, "Designs I - Content", false,
    [0.71, 0.77, 0.75, 0.88, 0.77, 0.76, 0.75, 0.81, 0.66], 0.77,
    [1.62, 1.44, 1.5, 1.04, 1.44, 1.47, 1.5, 1.31, 1.75], 1.46],
  ['sub', 'adult', 3, "Designs I - Spatial", false,
    [0.76, 0.78, 0.7, 0.78, 0.83, 0.73, 0.76, 0.71, 0.77], 0.76,
    [1.47, 1.41, 1.64, 1.41, 1.24, 1.56, 1.47, 1.62, 1.44], 1.48],
  ['sub', 'adult', 3, "Designs II", false,
    [0.88, 0.9, 0.87, 0.84, 0.87, 0.81, 0.8, 0.83, 0.82], 0.85,
    [1.04, 0.95, 1.08, 1.2, 1.08, 1.31, 1.34, 1.24, 1.27], 1.17],
  ['sub', 'adult', 3, "Designs II - Content", false,
    [0.81, 0.76, 0.79, 0.84, 0.79, 0.77, 0.7, 0.75, 0.71], 0.77,
    [1.31, 1.47, 1.37, 1.2, 1.37, 1.44, 1.64, 1.5, 1.62], 1.44],
  ['sub', 'adult', 3, "Designs II - Spatial", false,
    [0.74, 0.82, 0.69, 0.73, 0.75, 0.7, 0.68, 0.81, 0.67], 0.74,
    [1.53, 1.27, 1.67, 1.56, 1.5, 1.64, 1.7, 1.31, 1.72], 1.55],
  ['sub', 'adult', 3, "Visual Reproduction I", false,
    [0.92, 0.94, 0.88, 0.95, 0.93, 0.92, 0.93, 0.96, 0.92], 0.93,
    [0.85, 0.73, 1.04, 0.67, 0.79, 0.85, 0.79, 0.6, 0.85], 0.81],
  ['sub', 'older', 3, "Visual Reproduction I", false,
    [0.92, 0.94, 0.92, 0.92, 0.93], 0.93,
    [0.85, 0.73, 0.85, 0.85, 0.79], 0.82],
  ['sub', 'adult', 3, "Visual Reproduction II", false,
    [0.97, 0.98, 0.97, 0.97, 0.96, 0.97, 0.97, 0.98, 0.98], 0.97,
    [0.52, 0.42, 0.52, 0.52, 0.6, 0.52, 0.52, 0.42, 0.42], 0.5],
  ['sub', 'older', 3, "Visual Reproduction II", false,
    [0.95, 0.96, 0.96, 0.97, 0.96], 0.96,
    [0.67, 0.6, 0.6, 0.52, 0.6], 0.6],
  ['sub', 'adult', 3, "Spatial Addition", false,
    [0.89, 0.89, 0.92, 0.92, 0.89, 0.91, 0.92, 0.93, 0.9], 0.91,
    [0.99, 0.99, 0.85, 0.85, 0.99, 0.9, 0.85, 0.79, 0.95], 0.91],
  ['sub', 'adult', 3, "Symbol Span", false,
    [0.81, 0.88, 0.89, 0.9, 0.92, 0.86, 0.89, 0.88, 0.85], 0.88,
    [1.31, 1.04, 0.99, 0.95, 0.85, 1.12, 0.99, 1.04, 1.16], 1.06],
  ['sub', 'older', 3, "Symbol Span", false,
    [0.88, 0.87, 0.81, 0.84, 0.76], 0.84,
    [1.04, 1.08, 1.31, 1.2, 1.47], 1.23],
  ['idx', 'adult', 15, "Auditory Memory Index", false,
    [0.94, 0.95, 0.95, 0.97, 0.94, 0.96, 0.95, 0.95, 0.94], 0.95,
    [3.67, 3.35, 3.35, 2.6, 3.67, 3, 3.35, 3.35, 3.67], 3.35],
  ['idx', 'older', 15, "Auditory Memory Index", false,
    [0.95, 0.94, 0.95, 0.95, 0.95], 0.95,
    [3.35, 3.67, 3.35, 3.35, 3.35], 3.42],
  ['idx', 'adult', 15, "Visual Memory Index", false,
    [0.96, 0.97, 0.96, 0.96, 0.96, 0.96, 0.95, 0.96, 0.95], 0.96,
    [3, 2.6, 3, 3, 3, 3, 3.35, 3, 3.35], 3.04],
  ['idx', 'older', 15, "Visual Memory Index", false,
    [0.96, 0.97, 0.96, 0.97, 0.97], 0.97,
    [3, 2.6, 3, 2.6, 2.6], 2.77],
  ['idx', 'adult', 15, "Visual Working Memory Index", false,
    [0.89, 0.92, 0.94, 0.94, 0.94, 0.91, 0.93, 0.93, 0.92], 0.93,
    [4.97, 4.24, 3.67, 3.67, 3.67, 4.5, 3.97, 3.97, 4.24], 4.12],
  ['idx', 'adult', 15, "Immediate Memory Index", false,
    [0.93, 0.95, 0.95, 0.96, 0.95, 0.95, 0.94, 0.94, 0.93], 0.95,
    [3.97, 3.35, 3.35, 3, 3.35, 3.35, 3.67, 3.67, 3.97], 3.53],
  ['idx', 'older', 15, "Immediate Memory Index", false,
    [0.94, 0.95, 0.94, 0.94, 0.96], 0.95,
    [3.67, 3.35, 3.67, 3.67, 3], 3.48],
  ['idx', 'adult', 15, "Delayed Memory Index", false,
    [0.93, 0.95, 0.95, 0.96, 0.92, 0.94, 0.94, 0.94, 0.92], 0.94,
    [3.97, 3.35, 3.35, 3, 4.24, 3.67, 3.67, 3.67, 4.24], 3.71],
  ['idx', 'older', 15, "Delayed Memory Index", false,
    [0.92, 0.91, 0.91, 0.93, 0.92], 0.92,
    [4.24, 4.5, 4.5, 3.97, 4.24], 4.29]
];
const wms2 = (v) => Math.round(v * 100) / 100;

check('Table 3.3 reproduces from the stored coefficients, all 240 cells', () => {
  /* The bar, met arithmetically rather than on the manual's word — which is
     the lesson section 28 had to learn the hard way when Table 4.3 turned up
     and settled what had been recorded as unproven. Every printed SEM equals
     populationSD * sqrt(1 - the stored coefficient), at the printed 2dp. */
  const bad = [];
  let n = 0;
  WMS_T.forEach(([fam, bat, sd, name, retest, coef, , sems]) => {
    const e = D.normDB[WMS_GROUP(fam, bat)][name];
    const byAge = retest ? e.rStabilityByAge : e.rInternalByAge;
    if (!byAge) { bad.push(name + ' (' + bat + ') carries no by-age lookup'); return; }
    WMS_BANDS[bat].forEach((lo, i) => {
      n++;
      if (byAge[lo] !== coef[i]) { bad.push(name + ' ' + bat + ' @' + lo + ': stored ' + byAge[lo] + ', Table 3.1 prints ' + coef[i]); return; }
      const der = wms2(sd * Math.sqrt(1 - byAge[lo]));
      if (der !== sems[i]) bad.push(name + ' ' + bat + ' @' + lo + ': rxx ' + byAge[lo] + ' gives ' + der + ', Table 3.3 prints ' + sems[i]);
    });
  });
  if (n !== 240) bad.push('expected 240 published cells, tested ' + n);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the two batteries are kept apart, and 65-69 proves why', () => {
  /* WMS-IV is the reason separateBattery exists. Its two groups are not two
     bands of one norm set — they are two batteries with different subtests,
     and 65-69 is normed in BOTH with different coefficients. Age therefore
     cannot choose between them; the clinician did, when they decided which
     battery to administer, so the band has to stay selectable on Score Tables.

     Asserted on the data, not on the prose: the overlap must exist and must
     disagree, or the whole argument for separateBattery evaporates. */
  const bad = [];
  let overlap = 0, differ = 0;
  WMS_T.filter(([, bat]) => bat === 'older').forEach(([fam, , , name, retest, coef]) => {
    const adult = WMS_T.find(([f, b, , n]) => f === fam && b === 'adult' && n === name);
    if (!adult) { bad.push(name + ' is in the older battery but not the adult one'); return; }
    overlap++;
    // index 8 is the adult 65-69 band; index 0 is the older adult 65-69 band
    if (adult[5][8] !== coef[0]) differ++;
    if (retest !== adult[4]) bad.push(name + ' has a different reliability basis in each battery');
  });
  if (overlap !== 12) bad.push('expected 12 measures in both batteries, found ' + overlap);
  if (differ < 6) bad.push('only ' + differ + ' of ' + overlap + ' differ at 65-69 — the batteries have converged and separateBattery needs re-arguing');
  // the adult battery must hold the measures the older one drops
  ['Designs I', 'Designs II', 'Spatial Addition'].forEach((n) => {
    if (D.normDB['WMS-IV Subtests · Ages 65-90'][n]) bad.push(n + ' has appeared in the Older Adult battery');
  });
  if (D.normDB['WMS-IV Indices · Ages 65-90']['Visual Working Memory Index']) bad.push('VWMI has appeared in the Older Adult battery');
  return bad.length === 0 || bad.join('; ');
});

check('Score Tables keeps the WMS-IV battery selectable', () => {
  /* Without this the flat dropdown collapses the family to whichever group is
     listed first — Ages 16-69 — and an 80-year-old is shown the Adult
     battery's measure list and its coefficients. 10 of the 12 shared measures
     printed a different interval. Drives the shipped predicate. */
  const ctx = { normDB: D.normDB, getMergedDB: () => D.normDB };
  vm.createContext(ctx);
  vm.runInContext(extractFn(APP_SRC, 'familyScoredByAgeBand') + ';globalThis.__F = familyScoredByAgeBand;', ctx);
  const f = ctx.__F;
  const bad = [];
  ['WMS-IV Subtests · Ages 16-69', 'WMS-IV Subtests · Ages 65-90',
   'WMS-IV Indices · Ages 16-69', 'WMS-IV Indices · Ages 65-90']
    .forEach((g) => { if (!f(g)) bad.push(g + ' is not band-selectable'); });
  /* ...and the families that SHOULD collapse still do. RBANS is the pointed
     one: it went the other way, to All Ages groups, in the same body of work. */
  ['RBANS Subtests · Ages 12-19', 'RBANS Subtests · All Ages', 'WAIS-IV Core Subtests · All Ages',
   'CVLT-3 Indices · Ages 16-44'].forEach((g) => { if (f(g)) bad.push(g + ' was wrongly made band-selectable'); });
  // the base-rate case the predicate was written for must still hold
  if (!f('WAIS-IV Longest Span (Process) · Ages 16-17')) bad.push('the base-rate exemption has been lost');
  return bad.length === 0 || bad.join('; ');
});

check('the footnote-b measure is stability in both batteries, and matches rCorrected', () => {
  /* VPA II Word Recall is the one row Table 3.1 marks "based on test-retest",
     because free recall has no consistent item count to correlate across
     examinees. It is a CORRECTED stability coefficient — the manual says so
     (Allen & Yen 1979; Magnusson 1967) — so it must equal what this database
     already stores as rCorrected, and must not equal the raw r. Same
     transcription proof the WAIS-IV speeded three and the RBANS footnote-a
     five give. */
  const bad = [];
  let corr = 0, raw = 0;
  [['WMS-IV Subtests · Ages 16-69', 'adult'], ['WMS-IV Subtests · Ages 65-90', 'older']].forEach(([g, bat]) => {
    const e = D.normDB[g]['Verbal Paired Associates II - Word Recall'];
    const row = WMS_T.find(([, b, , n]) => b === bat && n === 'Verbal Paired Associates II - Word Recall');
    if (Number.isFinite(e.rInternal) || e.rInternalByAge) bad.push(g + ' calls the footnote-b measure internal consistency');
    if (!e.rStabilityByAge) bad.push(g + ' has no rStabilityByAge on the footnote-b measure');
    if (e.rCorrected === row[6]) corr++;
    if (e.r === row[6]) raw++;
    // the manual prints one value across every band of the battery
    if (new Set(row[5]).size !== 1) bad.push(bat + ': Table 3.1 no longer prints a constant across the battery');
  });
  if (corr !== 2) bad.push('matches stored rCorrected ' + corr + '/2, expected 2');
  if (raw > 0) bad.push('the raw r reproduces ' + raw + '/2 — the comparison has stopped discriminating');
  return bad.length === 0 || bad.join('; ');
});

check('no recognition memory measure has acquired a reliability coefficient', () => {
  /* Their published reliability is a decision-consistency PERCENTAGE, not a
     correlation, so it can never be fed to SEM = SD sqrt(1 - rxx). None is in
     normDB today. This check is here so that adding one — which would look
     like ordinary data entry — cannot pass silently. */
  const bad = [];
  Object.entries(D.normDB).forEach(([group, tab]) => {
    if (!/^WMS-IV /.test(group)) return;
    Object.entries(tab).forEach(([name, e]) => {
      if (!/Recognition/i.test(name)) return;
      bad.push(group + ' / ' + name + ' — recognition memory reliability is a percent agreement, not a coefficient');
    });
  });
  return bad.length === 0 || bad.join('; ');
});

check('WMS-IV averages its two tables by two different methods, as WAIS-IV and RBANS do', () => {
  /* Third manual, same pattern: Table 3.1 footnote a says Fisher's z, Table
     3.3 footnote a says the RMS of the band SEMs — which is algebraically
     SD sqrt(1 - the ARITHMETIC mean). So the average SEM cannot be derived
     from the average coefficient, and the stored average stays the published
     COEFFICIENT, per the decision recorded in data.js for WAIS-IV. */
  const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
  const tanh = (z) => (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  let fisher = 0, arith = 0, rms = 0, fromAvg = 0;
  WMS_T.forEach(([, , sd, , , coef, cAvg, sems, sAvg]) => {
    if (wms2(tanh(coef.reduce((a, r) => a + atanh(r), 0) / coef.length)) === cAvg) fisher++;
    if (wms2(coef.reduce((a, r) => a + r, 0) / coef.length) === cAvg) arith++;
    if (wms2(Math.sqrt(sems.reduce((a, s) => a + s * s, 0) / sems.length)) === sAvg) rms++;
    if (wms2(sd * Math.sqrt(1 - cAvg)) === sAvg) fromAvg++;
  });
  const bad = [];
  if (fisher !== 32) bad.push('Fisher z reproduces ' + fisher + '/32 of the Table 3.1 averages');
  if (fisher <= arith) bad.push('Fisher z (' + fisher + ') no longer beats the arithmetic mean (' + arith + ')');
  if (rms !== 32) bad.push('the RMS reading fits ' + rms + '/32 of the Table 3.3 averages');
  if (fromAvg > 12) bad.push('the average coefficient now reproduces ' + fromAvg + '/32 of the average SEMs — the two tables have converged');
  return bad.length === 0 || bad.join('; ');
});

check('an age inside the battery moves the interval; outside it falls back', () => {
  /* End to end through the shipped renderer, on both bases and both batteries.
     rInternalAgeMax is the battery's own upper bound — 69 for Adult, 90 for
     Older Adult — so an age outside it takes the published average rather
     than silently re-reading the topmost band. */
  const bad = [];
  const adultRow = { name: 'Logical Memory I', group: 'WMS-IV Subtests · Ages 16-69', scoreType: 'scaled' };
  const olderRow = { name: 'Logical Memory I', group: 'WMS-IV Subtests · Ages 65-90', scoreType: 'scaled' };
  const ea = D.normDB[adultRow.group][adultRow.name];
  const eo = D.normDB[olderRow.group][olderRow.name];
  if (batteryRel(adultRow, 'scaled', 25).r !== 0.87) bad.push('adult LM I at 25 did not read the 25-29 band');
  if (batteryRel(olderRow, 'scaled', 82).r !== 0.84) bad.push('older LM I at 82 did not read the 80-84 band');
  // the same age, different battery, different coefficient — the 65-69 overlap
  if (batteryRel(adultRow, 'scaled', 67).r === batteryRel(olderRow, 'scaled', 67).r) {
    bad.push('at 67 both batteries returned the same coefficient — the overlap is not being honoured');
  }
  // past the Adult battery's ceiling, fall back rather than re-read 65-69
  [null, 5, 75].forEach((age) => {
    const rel = batteryRel(adultRow, 'scaled', age);
    if (!rel) { bad.push('adult LM I gave no interval at age ' + age); return; }
    if (rel.r !== ea.rInternal) bad.push('adult LM I at age ' + age + ' used ' + rel.r + ', expected the average ' + ea.rInternal);
  });
  [null, 40, 95].forEach((age) => {
    const rel = batteryRel(olderRow, 'scaled', age);
    if (!rel) { bad.push('older LM I gave no interval at age ' + age); return; }
    if (rel.r !== eo.rInternal) bad.push('older LM I at age ' + age + ' used ' + rel.r + ', expected the average ' + eo.rInternal);
  });
  /* The claim worth asserting is that the BATTERY CHOICE changes the answer,
     not that it changes it for any one measure. Asserting the latter on LM I
     was wrong and this check caught it: at 82 the adult battery falls back to
     its average .82 and the older battery reads .84, and 3*sqrt(1-r) rounds
     both to the same margin. Coincidence, not a fault.

     So count across all 12 measures the two batteries share. Four differ, and
     the number is fixed by the pinned tables plus the fallback rule, so a
     change means one of those moved: VPA I, VPA II, VMI and DMI. */
  const shared = [];
  ['WMS-IV Subtests', 'WMS-IV Indices'].forEach((base) => {
    const sd = base.endsWith('Indices') ? 15 : 3;
    Object.keys(D.normDB[base + ' · Ages 65-90']).forEach((name) => {
      const a = batteryCi(sd === 15 ? 100 : 10, { name, group: base + ' · Ages 16-69', scoreType: sd === 15 ? 'standard' : 'scaled' }, '95', 82);
      const o = batteryCi(sd === 15 ? 100 : 10, { name, group: base + ' · Ages 65-90', scoreType: sd === 15 ? 'standard' : 'scaled' }, '95', 82);
      if (!a || !o) { bad.push(name + ' produced no interval in one of the batteries'); return; }
      if (a !== o) shared.push(name);
    });
  });
  if (shared.length !== 4) {
    bad.push('at 82 the batteries differ on ' + shared.length + ' of 12 measures, expected 4 (' + shared.join(', ') + ')');
  }
  return bad.length === 0 || bad.join('; ');
});


/* ==========================================================================
   31. WISC-V — Technical and Interpretive Manual Tables 4.1 and 4.4

   PINNED SOURCE: WISC-V Technical and Interpretive Manual, "Reliability and
   Errors of Measurement", pp. 56-62.
     p. 56      "Reliability coefficients were obtained utilizing the split-half
                method... The split-half coefficient is not a proper reliability
                estimate for Coding, Symbol Search, Cancellation, Naming Speed
                Literacy, Naming Speed Quantity, Immediate Symbol Translation,
                or Delayed Symbol Translation. Therefore, test-retest
                coefficients were used as the reliability estimates for these
                subtests... corrected for the normative sample's variability."
     Table 4.1  reliability by SINGLE YEAR of age, 6 to 16, plus an overall
                average via Fisher's z (footnote a).
     Table 4.4  the SEMs. Its note: "The reliability coefficients shown in
                Table 4.1 and the population standard deviations (i.e., 3 for
                the scaled scores and 15 for standard scores) were used to
                compute the SEMs." Average SEM is the RMS (footnote a).

   THE FOURTH PUBLISHER TO STATE THE SD RULE outright, after CVLT-3 (section
   23), WAIS-IV Table 4.3 and WMS-IV. Pair the coefficient with the NORMATIVE
   SD, never with sd1 — which is what 6de0d81 fixed.

   WISC-V IS THE ONLY FAMILY HERE WITH NO rCorrected AT ALL, so the usual
   transcription proof — matching the manual's stability rows against this
   database's stored rCorrected — is unavailable. Table 4.4 does that work
   instead: 242 cells cannot all reproduce from mistyped coefficients. The
   manual's own worked example is driven through the shipped renderer below,
   which is the same standard CVLT-C is held to in section 24.
   ========================================================================== */
heading('31. WISC-V — Tables 4.1 and 4.4');

/* [family, population SD, normDB name, stability?, Table 4.1 ages 6-16, its
   average, Table 4.4 ages 6-16, its average]. */
const WISC_AGES = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
const WISC_ALL = { sub: 'WISC-V Subtests · All Ages', idx: 'WISC-V Indices · All Ages' };
const WISC_T = [
  ['sub', 3, "Similarities", false,
    [0.89, 0.87, 0.85, 0.88, 0.88, 0.81, 0.87, 0.89, 0.85, 0.85, 0.87], 0.87,
    [0.99, 1.08, 1.16, 1.04, 1.04, 1.31, 1.08, 0.99, 1.16, 1.16, 1.08], 1.1],
  ['sub', 3, "Vocabulary", false,
    [0.83, 0.86, 0.83, 0.87, 0.87, 0.87, 0.91, 0.86, 0.88, 0.89, 0.9], 0.87,
    [1.24, 1.12, 1.24, 1.08, 1.08, 1.08, 0.9, 1.12, 1.04, 0.99, 0.95], 1.08],
  ['sub', 3, "Information", false,
    [0.82, 0.86, 0.81, 0.86, 0.82, 0.81, 0.89, 0.88, 0.85, 0.88, 0.9], 0.86,
    [1.27, 1.12, 1.31, 1.12, 1.27, 1.31, 0.99, 1.04, 1.16, 1.04, 0.95], 1.15],
  ['sub', 3, "Comprehension", false,
    [0.76, 0.86, 0.8, 0.84, 0.82, 0.79, 0.87, 0.82, 0.88, 0.83, 0.8], 0.83,
    [1.47, 1.12, 1.34, 1.2, 1.27, 1.37, 1.08, 1.27, 1.04, 1.24, 1.34], 1.26],
  ['sub', 3, "Block Design", false,
    [0.84, 0.86, 0.88, 0.83, 0.81, 0.83, 0.85, 0.82, 0.8, 0.86, 0.85], 0.84,
    [1.2, 1.12, 1.04, 1.24, 1.31, 1.24, 1.16, 1.27, 1.34, 1.12, 1.16], 1.2],
  ['sub', 3, "Visual Puzzles", false,
    [0.89, 0.9, 0.87, 0.9, 0.89, 0.88, 0.9, 0.89, 0.89, 0.92, 0.9], 0.89,
    [0.99, 0.95, 1.08, 0.95, 0.99, 1.04, 0.95, 0.99, 0.99, 0.85, 0.95], 0.98],
  ['sub', 3, "Matrix Reasoning", false,
    [0.89, 0.88, 0.89, 0.87, 0.85, 0.82, 0.9, 0.84, 0.84, 0.86, 0.86], 0.87,
    [0.99, 1.04, 0.99, 1.08, 1.16, 1.27, 0.95, 1.2, 1.2, 1.12, 1.12], 1.11],
  ['sub', 3, "Figure Weights", false,
    [0.91, 0.94, 0.94, 0.94, 0.96, 0.94, 0.95, 0.95, 0.94, 0.94, 0.93], 0.94,
    [0.9, 0.73, 0.73, 0.73, 0.6, 0.73, 0.67, 0.67, 0.73, 0.73, 0.79], 0.73],
  ['sub', 3, "Picture Concepts", false,
    [0.88, 0.82, 0.83, 0.82, 0.85, 0.85, 0.83, 0.8, 0.81, 0.8, 0.78], 0.83,
    [1.04, 1.27, 1.24, 1.27, 1.16, 1.16, 1.24, 1.34, 1.31, 1.34, 1.41], 1.26],
  ['sub', 3, "Arithmetic", false,
    [0.88, 0.9, 0.88, 0.9, 0.92, 0.9, 0.91, 0.91, 0.92, 0.89, 0.92], 0.9,
    [1.04, 0.95, 1.04, 0.95, 0.85, 0.95, 0.9, 0.9, 0.85, 0.99, 0.85], 0.94],
  ['sub', 3, "Digit Span", false,
    [0.92, 0.92, 0.9, 0.89, 0.9, 0.91, 0.93, 0.92, 0.92, 0.92, 0.92], 0.91,
    [0.85, 0.85, 0.95, 0.99, 0.95, 0.9, 0.79, 0.85, 0.85, 0.85, 0.85], 0.88],
  ['sub', 3, "Picture Span", false,
    [0.86, 0.82, 0.87, 0.87, 0.83, 0.84, 0.86, 0.85, 0.83, 0.84, 0.84], 0.85,
    [1.12, 1.27, 1.08, 1.08, 1.24, 1.2, 1.12, 1.16, 1.24, 1.2, 1.2], 1.18],
  ['sub', 3, "Letter-Number Sequencing", false,
    [0.93, 0.9, 0.83, 0.87, 0.8, 0.82, 0.85, 0.86, 0.89, 0.86, 0.82], 0.86,
    [0.79, 0.95, 1.24, 1.08, 1.34, 1.27, 1.16, 1.12, 0.99, 1.12, 1.27], 1.13],
  ['sub', 3, "Coding", true,
    [0.78, 0.78, 0.79, 0.79, 0.81, 0.81, 0.82, 0.82, 0.86, 0.86, 0.86], 0.82,
    [1.41, 1.41, 1.37, 1.37, 1.31, 1.31, 1.27, 1.27, 1.12, 1.12, 1.12], 1.28],
  ['sub', 3, "Symbol Search", true,
    [0.83, 0.83, 0.8, 0.8, 0.79, 0.79, 0.67, 0.67, 0.87, 0.87, 0.87], 0.81,
    [1.24, 1.24, 1.34, 1.34, 1.37, 1.37, 1.72, 1.72, 1.08, 1.08, 1.08], 1.34],
  ['sub', 3, "Cancellation", true,
    [0.8, 0.8, 0.83, 0.83, 0.84, 0.84, 0.81, 0.81, 0.81, 0.81, 0.81], 0.82,
    [1.34, 1.34, 1.24, 1.24, 1.2, 1.2, 1.31, 1.31, 1.31, 1.31, 1.31], 1.28],
  ['idx', 15, "Verbal Comprehension Index", false,
    [0.91, 0.92, 0.9, 0.93, 0.93, 0.9, 0.94, 0.93, 0.92, 0.92, 0.93], 0.92,
    [4.5, 4.24, 4.74, 3.97, 3.97, 4.74, 3.67, 3.97, 4.24, 4.24, 3.97], 4.22],
  ['idx', 15, "Visuospatial Index", false,
    [0.91, 0.92, 0.92, 0.91, 0.91, 0.91, 0.93, 0.91, 0.9, 0.93, 0.92], 0.92,
    [4.5, 4.24, 4.24, 4.5, 4.5, 4.5, 3.97, 4.5, 4.74, 3.97, 4.24], 4.36],
  ['idx', 15, "Fluid Reasoning Index", false,
    [0.93, 0.94, 0.94, 0.93, 0.93, 0.92, 0.95, 0.93, 0.93, 0.93, 0.93], 0.93,
    [3.97, 3.67, 3.67, 3.97, 3.97, 4.24, 3.35, 3.97, 3.97, 3.97, 3.97], 3.89],
  ['idx', 15, "Working Memory Index", false,
    [0.92, 0.91, 0.92, 0.92, 0.91, 0.92, 0.93, 0.92, 0.92, 0.92, 0.92], 0.92,
    [4.24, 4.5, 4.24, 4.24, 4.5, 4.24, 3.97, 4.24, 4.24, 4.24, 4.24], 4.26],
  ['idx', 15, "Processing Speed Index", false,
    [0.88, 0.88, 0.86, 0.87, 0.87, 0.88, 0.84, 0.84, 0.91, 0.91, 0.92], 0.88,
    [5.2, 5.2, 5.61, 5.41, 5.41, 5.2, 6, 6, 4.5, 4.5, 4.24], 5.24],
  ['idx', 15, "Full Scale IQ", false,
    [0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.96, 0.96, 0.97, 0.97], 0.96,
    [3, 3, 3, 3, 3, 3, 2.6, 3, 3, 2.6, 2.6], 2.9]
];
const wc2 = (v) => Math.round(v * 100) / 100;

check('Table 4.4 reproduces from the stored coefficients, all 242 cells', () => {
  const bad = [];
  let n = 0;
  WISC_T.forEach(([fam, sd, name, stab, coef, , sems]) => {
    const e = D.normDB[WISC_ALL[fam]][name];
    const byAge = stab ? e.rStabilityByAge : e.rInternalByAge;
    if (!byAge) { bad.push(name + ' carries no by-age lookup'); return; }
    WISC_AGES.forEach((a, i) => {
      n++;
      if (byAge[a] !== coef[i]) { bad.push(name + ' @' + a + ': stored ' + byAge[a] + ', Table 4.1 prints ' + coef[i]); return; }
      const der = wc2(sd * Math.sqrt(1 - byAge[a]));
      if (der !== sems[i]) bad.push(name + ' @' + a + ': rxx ' + byAge[a] + ' gives ' + der + ', Table 4.4 prints ' + sems[i]);
    });
  });
  if (n !== 242) bad.push('expected 242 cells, tested ' + n);
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check("the manual's worked example reproduces through the shipped renderer", () => {
  /* p. 62, and the strongest evidence available for this family because
     WISC-V stores no rCorrected to cross-check against:

       "If a 6-year-old child obtained an FSIQ score of 108, the practitioner
        can be 95% confident that the child's true FSIQ score falls in the
        range of 102-114 (because the 95% confidence interval is 108 +/- 1.96
        SEM, where the SEM is 3.00), and 90% confident that the child's true
        FSIQ score is in the range of 103-113 (108 +/- 1.65 SEM)."

     It also settles the CONVENTION. This manual documents two methods: its own
     Tables A.2-A.7 build intervals around the ESTIMATED TRUE score using the
     SEE (Dudek 1979), while the passage above gives the observed-score form
     for those who "prefer to calculate confidence intervals... in the most
     parsimonious manner". This app uses the latter and says so on screen, so
     the worked example is the right thing to reproduce — and it does, exactly,
     at both offered levels. */
  const row = { name: 'Full Scale IQ', group: WISC_ALL.idx, scoreType: 'standard' };
  const bad = [];
  const sem = 15 * Math.sqrt(1 - D.normDB[WISC_ALL.idx]['Full Scale IQ'].rInternalByAge[6]);
  if (wc2(sem) !== 3) bad.push('the age-6 FSIQ SEM is ' + wc2(sem) + ', the manual says 3.00');
  const strip = (h) => String(h).replace(/<[^>]*>/g, '').replace(/[–—]/g, '-').trim();
  const got95 = strip(batteryCi(108, row, '95', 6));
  const got90 = strip(batteryCi(108, row, '90', 6));
  if (got95 !== '102-114') bad.push('95% gave ' + got95 + ', the manual prints 102-114');
  if (got90 !== '103-113') bad.push('90% gave ' + got90 + ', the manual prints 103-113');
  return bad.length === 0 || bad.join('; ');
});

check('the three speeded subtests are stability, and are never called internal consistency', () => {
  /* The manual names them itself, and names four complementary subtests
     alongside them that this app does not hold. Asserted in both directions:
     these three must carry rStability and no internal-consistency field, and
     no OTHER WISC-V measure may carry rStability. */
  const bad = [];
  const SPEEDED = ['Coding', 'Symbol Search', 'Cancellation'];
  Object.entries(D.normDB[WISC_ALL.sub]).forEach(([name, e]) => {
    const speeded = SPEEDED.includes(name);
    if (speeded) {
      if (!e.rStabilityByAge) bad.push(name + ' is speeded but has no rStabilityByAge');
      if (Number.isFinite(e.rInternal) || e.rInternalByAge) bad.push(name + ' is speeded yet carries an internal-consistency field');
    } else if (Number.isFinite(e.rStability) || e.rStabilityByAge) {
      bad.push(name + ' is not speeded yet carries a stability field');
    }
  });
  /* PSI is a hybrid and is included anyway, exactly as WAIS-IV's is: p. 61
     concedes its average "is based on test-retest reliabilities", but the
     manual computes and labels every composite coefficient as internal
     consistency and publishes no other reliability for them. */
  const psi = D.normDB[WISC_ALL.idx]['Processing Speed Index'];
  if (!psi.rInternalByAge) bad.push('PSI lost its internal-consistency coefficient');
  return bad.length === 0 || bad.join('; ');
});

check('the age lookup is by single year, and stops at the normed range', () => {
  /* WISC-V is normed 6-16, so rInternalAgeMax is 16 and the keys are single
     years rather than bands — the finest lookup in the database. An age
     outside the range must fall back to the published average rather than
     re-reading age 16, which for an adult mistyped into a child's record
     would otherwise pass silently. */
  const bad = [];
  const row = { name: 'Vocabulary', group: WISC_ALL.sub, scoreType: 'scaled' };
  const e = D.normDB[WISC_ALL.sub].Vocabulary;
  if (Object.keys(e.rInternalByAge).length !== 11) bad.push('Vocabulary has ' + Object.keys(e.rInternalByAge).length + ' keys, expected 11 single years');
  if (batteryRel(row, 'scaled', 6).r !== 0.83) bad.push('age 6 did not read its own year');
  if (batteryRel(row, 'scaled', 16).r !== 0.9) bad.push('age 16 did not read its own year');
  [null, 5, 17, 40].forEach((age) => {
    const rel = batteryRel(row, 'scaled', age);
    if (!rel) { bad.push('no interval at age ' + age); return; }
    if (rel.r !== e.rInternal) bad.push('age ' + age + ' used ' + rel.r + ', expected the published average ' + e.rInternal);
  });
  // a single year must actually move the printed interval somewhere
  const moved = Object.keys(D.normDB[WISC_ALL.sub]).filter((name) => {
    const r = { name, group: WISC_ALL.sub, scoreType: 'scaled' };
    return batteryCi(10, r, '95', 6) !== batteryCi(10, r, '95', 16);
  });
  if (!moved.length) bad.push('no subtest interval differs between age 6 and age 16');
  return bad.length === 0 || bad.join('; ');
});

check('WISC-V splits its two average columns the way the other manuals do', () => {
  /* Fourth manual with the pattern recorded in section 28: Table 4.1's average
     is Fisher's z, Table 4.4's is the RMS of the band SEMs — which is
     SD sqrt(1 - the ARITHMETIC mean). The two cannot be reconciled, and the
     stored average stays the published COEFFICIENT. */
  const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
  const tanh = (z) => (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  let fisher = 0, arith = 0, rms = 0, fromAvg = 0;
  WISC_T.forEach(([, sd, , , coef, cAvg, sems, sAvg]) => {
    if (wc2(tanh(coef.reduce((a, r) => a + atanh(r), 0) / coef.length)) === cAvg) fisher++;
    if (wc2(coef.reduce((a, r) => a + r, 0) / coef.length) === cAvg) arith++;
    if (wc2(Math.sqrt(sems.reduce((a, s) => a + s * s, 0) / sems.length)) === sAvg) rms++;
    if (wc2(sd * Math.sqrt(1 - cAvg)) === sAvg) fromAvg++;
  });
  const bad = [];
  if (fisher !== 22) bad.push('Fisher z reproduces ' + fisher + '/22 of the Table 4.1 averages');
  if (fisher <= arith) bad.push('Fisher z (' + fisher + ') no longer beats the arithmetic mean (' + arith + ')');
  if (rms !== 22) bad.push('the RMS reading fits ' + rms + '/22 of the Table 4.4 averages');
  if (fromAvg > 8) bad.push('the average coefficient now reproduces ' + fromAvg + '/22 of the average SEMs');
  return bad.length === 0 || bad.join('; ');
});

check('Table 4.7 — the retest study, pinned for every measure it covers', () => {
  /* PINNED SOURCE: WISC-V Technical and Interpretive Manual, Table 4.7,
     "Stability Coefficients of Subtest, Process, and Composite Scores", all
     ages. 34 rows: first- and second-testing means and SDs, n, r12, and a
     corrected r.

     Transcribed from a photograph of the table rather than a spreadsheet, so
     it was checked two ways before being stored. Cohen's (1996) Formula 10.4
     reproduces the printed Standard Difference from the means and SDs on 34 of
     34 rows, which no misread digit survives; and the 22 measures already in
     normDB matched their m1/sd1/m2/sd2/n/r exactly, which also proves the r12
     column was not confused with the corrected one — comparing stored r
     against the corrected column matches 1 of 22, against r12 22 of 22. */
  const NEWROWS = [
  ['proc', 'Block Design No Time Bonus', 9.6, 2.8, 10.8, 3.2, 207, 0.78, 0.81],
  ['proc', 'Block Design Partial', 9.8, 3.1, 10.8, 3.1, 208, 0.84, 0.84],
  ['proc', 'Digit Span Forward', 9.9, 2.7, 10.2, 3.0, 213, 0.78, 0.82],
  ['proc', 'Digit Span Backward', 9.9, 2.8, 10.2, 2.9, 201, 0.72, 0.76],
  ['proc', 'Digit Span Sequencing', 9.5, 2.7, 10.1, 2.9, 200, 0.74, 0.79],
  ['proc', 'Cancellation Random', 9.9, 2.9, 11.1, 2.9, 200, 0.78, 0.81],
  ['proc', 'Cancellation Structured', 9.9, 2.8, 10.9, 3.1, 209, 0.78, 0.82],
  ['idx', 'Quantitative Reasoning Index', 99.2, 13.4, 102.4, 13.7, 216, 0.76, 0.81],
  ['idx', 'Auditory Working Memory Index', 98.7, 13.3, 100.9, 14.6, 217, 0.85, 0.88],
  ['idx', 'Nonverbal Index', 98.5, 14.0, 105.5, 13.8, 216, 0.86, 0.88],
  ['idx', 'General Ability Index', 98.0, 13.7, 103.6, 13.3, 213, 0.89, 0.91],
  ['idx', 'Cognitive Proficiency Index', 99.3, 14.0, 105.5, 14.8, 212, 0.84, 0.86],
  ];
  const BACKFILL = [
  ['sub', 'Similarities', 0.88],
  ['sub', 'Vocabulary', 0.9],
  ['sub', 'Information', 0.88],
  ['sub', 'Comprehension', 0.83],
  ['sub', 'Block Design', 0.81],
  ['sub', 'Visual Puzzles', 0.8],
  ['sub', 'Matrix Reasoning', 0.78],
  ['sub', 'Figure Weights', 0.82],
  ['sub', 'Picture Concepts', 0.71],
  ['sub', 'Arithmetic', 0.84],
  ['sub', 'Digit Span', 0.82],
  ['sub', 'Picture Span', 0.8],
  ['sub', 'Letter-Number Sequencing', 0.82],
  ['sub', 'Coding', 0.81],
  ['sub', 'Symbol Search', 0.8],
  ['sub', 'Cancellation', 0.82],
  ['idx', 'Verbal Comprehension Index', 0.94],
  ['idx', 'Visuospatial Index', 0.84],
  ['idx', 'Fluid Reasoning Index', 0.75],
  ['idx', 'Working Memory Index', 0.82],
  ['idx', 'Processing Speed Index', 0.83],
  ['idx', 'Full Scale IQ', 0.92],
  ];
  const G = { sub:'WISC-V Subtests · All Ages', idx:'WISC-V Indices · All Ages',
              proc:'WISC-V Process Scores · All Ages' };
  const bad = [];
  NEWROWS.forEach(([g, name, m1, sd1, m2, sd2, n, r, rc]) => {
    const e = D.normDB[G[g]][name];
    if (!e) { bad.push(name + ' is missing'); return; }
    const want = { m1, sd1, m2, sd2, n, r, rCorrected: rc };
    Object.entries(want).forEach(([f, v]) => {
      if (e[f] !== v) bad.push(name + ' ' + f + ': stored ' + e[f] + ', Table 4.7 prints ' + v);
    });
  });
  BACKFILL.forEach(([g, name, rc]) => {
    const e = D.normDB[G[g]][name];
    if (e.rCorrected !== rc) bad.push(name + ' rCorrected: stored ' + e.rCorrected + ', Table 4.7 prints ' + rc);
  });
  if (NEWROWS.length + BACKFILL.length !== 34) bad.push('expected all 34 rows of Table 4.7, have ' + (NEWROWS.length + BACKFILL.length));
  /* The checksum that validated the photograph, kept live: if a mean or SD is
     ever edited, its row's Standard Difference stops reconciling. */
  const SD_DIFF = { 'Block Design No Time Bonus':0.40, 'Block Design Partial':0.32,
    'Digit Span Forward':0.11, 'Digit Span Backward':0.11, 'Digit Span Sequencing':0.21,
    'Cancellation Random':0.41, 'Cancellation Structured':0.34,
    'Quantitative Reasoning Index':0.24, 'Auditory Working Memory Index':0.16,
    'Nonverbal Index':0.50, 'General Ability Index':0.41, 'Cognitive Proficiency Index':0.43 };
  Object.entries(SD_DIFF).forEach(([name, want]) => {
    const g = name.includes('Index') ? 'idx' : 'proc';
    const e = D.normDB[G[g]][name];
    const d = (e.m2 - e.m1) / Math.sqrt((e.sd1 ** 2 + e.sd2 ** 2) / 2);
    if (Math.round(d * 100) / 100 !== want) bad.push(name + ': standard difference ' + d.toFixed(3) + ', Table 4.7 prints ' + want);
  });
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the two Cancellation process scores are stability, the other five are not', () => {
  /* The manual names Cancellation among the subtests for which "the split-half
     coefficient is not a proper reliability estimate", and its two process
     scores inherit that. Confirmed from the data as well as the prose: a
     stability coefficient is broadcast from the retest study's coarse bands,
     so it repeats across the 11 single-year bands, where an internal-
     consistency row varies freely. */
  const P = D.normDB['WISC-V Process Scores · All Ages'];
  const distinct = (name) => {
    const e = P[name];
    const t = e.rInternalByAge || e.rStabilityByAge;
    return new Set(Object.values(t)).size;
  };
  const bad = [];
  ['Cancellation Random', 'Cancellation Structured'].forEach((n) => {
    if (!P[n].rStabilityByAge) bad.push(n + ' is not held as a stability coefficient');
    if (Number.isFinite(P[n].rInternal) || P[n].rInternalByAge) bad.push(n + ' is labelled internal consistency');
    if (distinct(n) > 4) bad.push(n + ' varies across ' + distinct(n) + ' values — that is not a broadcast stability coefficient');
  });
  ['Block Design No Time Bonus', 'Block Design Partial', 'Digit Span Forward',
   'Digit Span Backward', 'Digit Span Sequencing'].forEach((n) => {
    if (!P[n].rInternalByAge) bad.push(n + ' lost its internal-consistency lookup');
    if (Number.isFinite(P[n].rStability) || P[n].rStabilityByAge) bad.push(n + ' is labelled stability');
    if (distinct(n) < 6) bad.push(n + ' varies across only ' + distinct(n) + ' values — check whether it is really internal consistency');
  });
  return bad.length === 0 || bad.join('; ');
});

check('the stored averages are Table 4.1\'s own average column', () => {
  const bad = [];
  WISC_T.forEach(([fam, , name, stab, , cAvg]) => {
    const e = D.normDB[WISC_ALL[fam]][name];
    const got = stab ? e.rStability : e.rInternal;
    if (got !== cAvg) bad.push(name + ': stored ' + got + ', Table 4.1 prints ' + cAvg);
    const max = stab ? e.rStabilityAgeMax : e.rInternalAgeMax;
    if (max !== 16) bad.push(name + ' has age max ' + max + ', expected 16');
  });
  return bad.length === 0 || bad.join('; ');
});


/* ==========================================================================
   32. Norms Database view

   The page that lets a clinician SEE what the app scores from. Until now it
   listed all 115 groups flat and alphabetically — D-KEFS alone ran to 36
   consecutive groups — and its columns stopped at the retest r, so the 131
   entries whose interval is built from a published internal-consistency or
   stability coefficient showed no sign of it.

   The risk a reliability column introduces is that it drifts from the renderer
   and starts advertising a coefficient the app does not use. That is what the
   first check exists to prevent.
   ========================================================================== */
heading('32. Norms Database view');

const dbCtx = (() => {
  /* No #bat-ci-basis here, which is the point: batteryCiCorrectRetest must
     read "published" from a missing control, so a harness (or a page that
     failed to build the bar) can never land on derived coefficients by
     accident. The checks below pass the flag explicitly to drive both states. */
  const ctx = { normDB: D.normDB, getMergedDB: () => D.normDB,
                document: { getElementById: () => null } };
  vm.createContext(ctx);
  vm.runInContext(
    (APP_SRC.match(/const BATTERY_METRIC_SD = \{[^;]*;/) || [''])[0] + '\n' +
    (APP_SRC.match(/const SCORE_METRICS = [^;]*;/) || [''])[0] + '\n' +
    ['inferScoreType', 'inferScoreTypeForSubtest', 'bandedReliabilityForAge',
     'rInternalForAge', 'rStabilityForAge', 'derivedCorrectedR',
     'batteryCiCorrectRetest', 'resolveCiReliability',
     'dbReliabilityBasis', 'dbInstrumentOf'].map((n) => extractFn(APP_SRC, n)).join('\n') +
    ';globalThis.__BASIS = dbReliabilityBasis; globalThis.__INST = dbInstrumentOf;'
    + 'globalThis.inferScoreTypeForSubtest = inferScoreTypeForSubtest;',
    ctx
  );
  return ctx;
})();
const dbBasis = dbCtx.__BASIS;
const dbInstrument = dbCtx.__INST;
const dbType = dbCtx.inferScoreTypeForSubtest;

check('the reliability control costs what Methods & References says it costs', () => {
  /* Methods & References quantifies the control so a clinician can decide
     whether to touch it with the size of the change in view: "46 of the
     measures reachable from Score Tables, every one of them D-KEFS or D-KEFS
     Advanced: at the 95% level 9 intervals widen and 4 narrow, 33 are unchanged
     after rounding, and the largest single change is 2 scaled-score points."

     Those five figures were unpinned. Section 32 already asserts the control
     moves 246 stored ENTRIES, but that is the whole database; the sentence on
     screen is about the far smaller set a clinician can actually select, and
     nothing tied the two. A new family entering the database on a bare retest r
     would move the printed counts without failing anything.

     Derived here rather than transcribed: the reachable set is reconstructed
     the way buildFamilyListHtml picks it in flat mode (one canonical All Ages
     group per family base, except where the age band is itself the lookup), and
     the half-widths come from the SHIPPED resolver in both states. So this
     re-runs the arithmetic behind the sentence rather than restating it. */
  const Z95 = 1.96;
  const SD = { standard: 15, t: 10, scaled: 3, z: 1 };
  const hasBandSuffix = (f) => /·\s*[^·]*$/.test(f) && /·/.test(f);
  const baseOf = (f) => (f.includes('·') ? f.slice(0, f.lastIndexOf('·')).trim() : f);

  const groups = {};
  const order = [];
  Object.keys(D.normDB).forEach((f) => {
    const base = hasBandSuffix(f) ? baseOf(f) : f;
    if (!groups[base]) { groups[base] = []; order.push(base); }
    groups[base].push(f);
  });
  const reachable = [];
  order.forEach((base) => {
    const members = groups[base];
    /* familyScoredByAgeBand: the band IS the lookup, so every band stays
       selectable instead of collapsing to one canonical entry. */
    const banded = members.some((f) => Object.values(D.normDB[f])
      .some((e) => e && typeof e === 'object' && (e.baseRates || e.separateBattery)));
    if (members.length === 1 && !hasBandSuffix(members[0])) reachable.push(members[0]);
    else if (banded) members.forEach((f) => reachable.push(f));
    else reachable.push(members.find((m) => /·\s*All\s+Ages\s*$/i.test(m)) || members[0]);
  });

  let moved = 0, wider = 0, narrower = 0, same = 0, maxShift = 0;
  const families = new Set();
  reachable.forEach((group) => {
    Object.entries(D.normDB[group]).forEach(([name, e]) => {
      if (!e || typeof e !== 'object') return;
      const normSD = SD[dbType(group, name, e)];
      const pub = dbCtx.resolveCiReliability(e, normSD, undefined, false);
      const cor = dbCtx.resolveCiReliability(e, normSD, undefined, true);
      if (!pub || pub.none || !cor || cor.none || pub.r === cor.r) return;
      moved++;
      families.add(dbInstrument(group));
      const sd = Number.isFinite(normSD) ? normSD : e.sd1;
      const hp = Math.round(Z95 * sd * Math.sqrt(1 - pub.r));
      const hc = Math.round(Z95 * sd * Math.sqrt(1 - cor.r));
      if (hc > hp) wider++; else if (hc < hp) narrower++; else same++;
      maxShift = Math.max(maxShift, Math.abs(hc - hp));
    });
  });

  const bad = [];
  const got = { moved, wider, narrower, same, maxShift };
  const want = { moved: 46, wider: 9, narrower: 4, same: 33, maxShift: 2 };
  Object.keys(want).forEach((k) => {
    if (got[k] !== want[k]) bad.push(k + ' is ' + got[k] + ', Methods & References says ' + want[k]);
  });
  const stray = [...families].filter((f) => !/^D-KEFS/.test(f));
  if (stray.length) bad.push('the control also moves ' + stray.join(', ')
    + '; Methods & References says every moved measure is D-KEFS or D-KEFS Advanced');

  /* And the sentence must still be on the page saying it. Asserted against the
     rendered figures above, so editing one without the other fails. */
  const block = methodsCiBlock().replace(/<[^>]+>/g, '');
  [[String(want.moved) + ' of the measures reachable', 'the count of moved measures'],
   [String(want.wider) + ' intervals widen', 'the count that widen'],
   [String(want.narrower) + ' narrow', 'the count that narrow'],
   [String(want.same) + ' are unchanged', 'the count left unchanged'],
   [String(want.maxShift) + ' scaled-score points', 'the largest single change']]
    .forEach(([needle, what]) => {
      if (!block.includes(needle)) bad.push(what + ' is no longer stated on the Methods page');
    });

  return bad.length === 0 || bad.join('; ');
});

check('the reliability column shows what the app would actually use', () => {
  /* THE LOAD-BEARING CHECK. Score Tables renders the interval; this page tells
     the clinician what it rests on. If they diverge, the page advertises a
     coefficient the report does not use.

     The two used to be separate functions reading the same fields in the same
     order, pinned together here — a mirror, and this check had to be
     strengthened twice because a mirror is only as good as the next person's
     memory. They now share resolveCiReliability, so the ORDER cannot drift.
     What can still drift is the ARGUMENTS each entry point passes in, which is
     what this now covers: driven over every entry at a blank age, through both
     public functions, comparing the coefficient AND the basis string.

     The score type is the app's own inferScoreTypeForSubtest rather than a
     shortcut restated here, because that is what an auto-filled Score Tables
     row carries and what dbReliabilityBasis asks for. Restating it would let
     the two disagree about the metric while agreeing about everything else.

     DRIVEN IN BOTH STATES OF THE RELIABILITY-BASIS CONTROL. That control lives
     on Score Tables but this page has to follow it: showing the published
     reading while the table is set to corrected is exactly the defect the
     shared resolver exists to prevent, and it would be invisible — the page
     would look right, on its own terms, while contradicting the report. */
  const bad = [];
  let n = 0, withR = 0, moved = 0;
  Object.entries(D.normDB).forEach(([group, tab]) => {
    Object.entries(tab).forEach(([name, e]) => {
      if (!e || typeof e !== 'object') return;
      n++;
      const type = dbType(group, name, e);
      [false, true].forEach((corr) => {
        const shown = dbBasis(e, group, name, corr);
        const used = batteryRel({ name, group, scoreType: type }, type, undefined, corr);
        const at = corr ? ' [corrected]' : '';
        if (!used) {
          if (!shown.none) bad.push(group + ' / ' + name + at + ': table shows ' + shown.r + ' but the app builds no interval');
          return;
        }
        if (!corr) withR++;
        if (shown.r !== used.r) bad.push(group + ' / ' + name + at + ': table shows ' + shown.r + ', the renderer uses ' + used.r);
        if (shown.basis !== used.basis) bad.push(group + ' / ' + name + at + ': table says "' + shown.basis + '", the renderer says "' + used.basis + '"');
      });
      if (dbBasis(e, group, name, true).r !== dbBasis(e, group, name, false).r) moved++;
    });
  });
  if (n !== 671) bad.push('expected 671 entries, walked ' + n);
  if (withR < 550) bad.push('only ' + withR + ' entries produced an interval — the comparison has stopped covering the database');
  /* And the control must actually do something, or the agreement above is
     vacuous — two functions can hardly disagree about a flag neither reads.
     From the 298 entries on the retest pairing:

       -36  CVLT-C, sd1 in raw words while the row displays T or z
       - 2  corrected value outside (0, 1), so not a reliability
       -14  sd1 exactly equals the normative SD, so the correction is identity
       ---
        246 actually move

     The 14 are worth knowing about before treating a shortfall as a bug: a
     D-KEFS retest sample whose SD is a round 3.0 is unrestricted by
     definition, and there is nothing for the formula to do. */
  if (moved !== 246) bad.push('the basis control moves ' + moved + ' entries, expected 246');
  /* And the shared resolver must actually be shared. Both entry points calling
     it is what makes the order impossible to drift; if one grows its own copy
     of the chain the check above would still pass on the day it was written. */
  if (!/resolveCiReliability\(/.test(extractFn(APP_SRC, 'dbReliabilityBasis'))) {
    bad.push('dbReliabilityBasis no longer calls resolveCiReliability — the page has its own copy of the chain again');
  }
  if (!/resolveCiReliability\(/.test(extractFn(APP_SRC, 'getBatteryRowReliability'))) {
    bad.push('getBatteryRowReliability no longer calls resolveCiReliability');
  }
  return bad.length === 0 || bad.slice(0, 4).join('; ');
});

check('the basis names the right source, and says when age will move it', () => {
  /* The figure alone is not enough: .90 from a manual's split-half table and
     .90 from a retest study are different claims, and CLAUDE.md treats an
     on-screen statement of method as a contract. Also asserts the "by age"
     suffix, without which a clinician would not know the printed value changes
     the moment an age is entered on Score Tables. */
  const bad = [];
  const cases = [
    ['WAIS-IV Indices · All Ages', 'Full Scale IQ', 'internal consistency · by age'],
    ['RBANS Subtests · All Ages', 'Coding', 'stability, published · by age'],
    ['WMS-IV Subtests · Ages 16-69', 'Verbal Paired Associates II - Word Recall', 'stability, published · by age'],
    ['WISC-V Subtests · All Ages', 'Symbol Search', 'stability, published · by age'],
    ['CVLT-3 Indices · Ages 16-44', 'T1-5 Correct', 'retest, corrected'],
    ['WISC-V Indices · All Ages', 'Fluid Reasoning Index', 'internal consistency · by age'],
    ['RBANS Subtests · All Ages', 'List Recognition', 'none published'],
    ['WAIS-IV Longest Span (Process) · Ages 16-17', 'Longest Digit Span Forward', 'base rate — no interval'],
    /* THE TWO RETEST LABELS, WHICH MUST NOT COLLAPSE INTO EACH OTHER. Both
       rest on the retest study's own correlation, but they are different
       claims about the SD it is paired with, and the whole point of the pair
       is that a clinician can tell them apart.

       "retest, uncorrected" — a normative SD is in force, so the coefficient
       and the SD describe different populations. That is the publisher's own
       pairing for these measures (D-KEFS Technical Manual p. 19), not a defect
       to be silently repaired; see resolveCiReliability.
       "retest" — the row is displayed raw, so sd1 IS the SD in force and the
       coefficient sits beside its own sample's SD. No qualifier is warranted,
       and adding one would read as a fault where there is none.

       CVLT-C lands on "uncorrected" and belongs there, which is not obvious:
       metric:'raw' makes it look like the second case, but reportedAs puts the
       row on a T or z metric, so the SD in force is the normative 10 or 1
       while r was measured on raw words recalled. Same mismatch, and the one
       place where no correction could repair it even in principle. */
    ['D-KEFS Colour-Word Interference · All Ages', 'Colour Naming', 'retest, uncorrected'],
    ['D-KEFS Tower Test · All Ages', 'Total Achievement Score', 'retest, uncorrected · by age'],
    ['CVLT-C Subtests (Raw Scores) · Age 8', 'Perseverations', 'retest, uncorrected'],
    ['RBANS Subtests · Ages 20-89', 'List Recognition', 'retest']
  ];
  cases.forEach(([g, n, want]) => {
    const e = D.normDB[g] && D.normDB[g][n];
    if (!e) { bad.push('missing fixture ' + g + ' / ' + n); return; }
    const got = dbBasis(e, g, n).basis;
    if (got !== want) bad.push(g + ' / ' + n + ': basis "' + got + '", expected "' + want + '"');
  });
  /* The two no-interval reasons must stay distinguishable. A base-rate measure
     never had a coefficient; a raw RBANS subtest is absent from its manual's
     reliability table. Both print nothing, for quite different reasons. */
  const kinds = new Set();
  Object.entries(D.normDB).forEach(([g, tab]) => Object.entries(tab).forEach(([n, e]) => {
    if (e && typeof e === 'object') { const b = dbBasis(e, g, n); if (b.none) kinds.add(b.basis); }
  }));
  if (kinds.size !== 2) bad.push('expected 2 distinct no-interval reasons, found ' + [...kinds].join(' / '));
  /* THE WHOLE DISTRIBUTION, so that a change of basis anywhere in the database
     has to be acknowledged here rather than sliding through unnoticed. These
     are counts of what the app builds its intervals from today:

       internal consistency        118  (3 bare + 115 by age)
       stability, published         12
       retest, corrected           180  a published rCorrected
       retest, uncorrected         298  (285 + 13 by age) — the retest study's
                                        own r used with a normative SD, which
                                        is these manuals' own pairing
       retest                        8  raw display, r beside its own sample SD
       none published / base rate   55  no interval at all

     PLAIN "retest" IS EXACTLY THE FOUR RAW RBANS SUBTESTS IN EACH OF THE TWO
     RETEST BANDS. Nothing else in the database is displayed raw, so if that 8
     moves, either a family has been tagged metric:'raw' or one has lost the
     tag — both of which change what a percentile means, not just a label. */
  const tally = {};
  Object.entries(D.normDB).forEach(([g, tab]) => Object.entries(tab).forEach(([n, e]) => {
    if (!e || typeof e !== 'object') return;
    const b = dbBasis(e, g, n).basis;
    tally[b] = (tally[b] || 0) + 1;
  }));
  const want = {
    'internal consistency': 3, 'internal consistency · by age': 115,
    'stability, published · by age': 12, 'retest, corrected': 180,
    'retest, uncorrected': 285, 'retest, uncorrected · by age': 13,
    'retest': 8, 'none published': 4, 'base rate — no interval': 51
  };
  Object.entries(want).forEach(([k, v]) => {
    if (tally[k] !== v) bad.push('basis "' + k + '": ' + (tally[k] || 0) + ' entries, expected ' + v);
  });
  Object.keys(tally).forEach((k) => { if (!(k in want)) bad.push('unaccounted basis "' + k + '"'); });
  return bad.length === 0 || bad.join('; ');
});

check('every group resolves to exactly one instrument', () => {
  /* The Instrument column and its filter are built from the data at render
     time, so a new family cannot be left out of either. D-KEFS Advanced must
     be tested before plain D-KEFS, the latter being a prefix of the former —
     the same ordering trap caIntervalLabel documents. */
  const bad = [];
  const tabs = {};
  Object.keys(D.normDB).forEach((g) => {
    const inst = dbInstrument(g, false);
    tabs[inst] = (tabs[inst] || 0) + 1;
  });
  const want = { 'CVLT-3':4, 'CVLT-C':3, 'D-KEFS':36, 'D-KEFS Advanced':23, 'RBANS':9, 'WAIS-IV':34, 'WISC-V':3, 'WMS-IV':4 };
  Object.entries(want).forEach(([k, v]) => {
    if (tabs[k] !== v) bad.push(k + ': ' + (tabs[k] || 0) + ' groups, expected ' + v);
  });
  const stray = Object.keys(tabs).filter((k) => !(k in want));
  if (stray.length) bad.push('unexpected instrument(s): ' + stray.join(', '));
  if (dbInstrument('anything at all', true) !== 'Custom') bad.push('a custom family is not labelled Custom');
  if (dbInstrument('D-KEFS Advanced Tower · All Ages', false) === 'D-KEFS') bad.push('D-KEFS Advanced is being swallowed by the D-KEFS prefix');
  return bad.length === 0 || bad.join('; ');
});

check('every column shows and sorts on the same value', () => {
  /* THE GRID IS GONE, AND WITH IT THE HAZARD THE OLD CHECK GUARDED. The
     grouped list laid each row out as a CSS grid whose column count was
     declared in two places; a mismatch slid every cell one place left, under
     the wrong heading. This page is a real <table>, which cannot do that.

     What CAN go wrong now is a column sorting on something other than what it
     displays — a "CI r" heading ordered by the retest r would be invisible and
     completely wrong. DB_COLUMNS gives each column ONE `get`, used for both,
     so the two cannot diverge; this pins that shape rather than the values. */
  const src = (APP_SRC.match(/const DB_COLUMNS = \[[\s\S]*?\n\];/) || [''])[0];
  if (!src) return 'DB_COLUMNS is gone';
  const bad = [];
  const keys = [...src.matchAll(/\{ key:'([a-zA-Z0-9]+)'/g)].map((m) => m[1]);
  const gets = (src.match(/get:r =>/g) || []).length;
  if (keys.length !== 13) bad.push('expected 13 columns, found ' + keys.length);
  if (gets !== keys.length) bad.push(keys.length + ' columns but ' + gets + ' getters — a column has no single source');
  ['instrument', 'category', 'band', 'measure', 'ci', 'basis'].forEach((k) => {
    if (!keys.includes(k)) bad.push('the ' + k + ' column is gone');
  });
  if (new Set(keys).size !== keys.length) bad.push('two columns share a key, so sorting one would point at the other');
  if (!/col\.get\(a\)[\s\S]{0,40}col\.get\(b\)/.test(APP_SRC)) {
    bad.push('the sort no longer reads the column getter, so it can order by something other than what is shown');
  }
  /* AND THE TWO RELIABILITY COLUMNS MUST READ dbReliabilityBasis, not a raw
     field. Found by mutation: pointing the CI r getter at r.e.r passed
     everything, because the column then showed and sorted on the same thing —
     the retest coefficient — under a heading promising the one the interval is
     built from. The neighbouring check proves dbReliabilityBasis is right; this
     proves the column actually asks it. */
  if (!/key:'ci'[^}]*get:r => r\.rel\.r/.test(src)) {
    bad.push('the CI r column no longer reads dbReliabilityBasis, so it can print a coefficient the interval does not use');
  }
  if (!/key:'basis'[^}]*get:r => r\.rel\.basis/.test(src)) {
    bad.push('the Basis column no longer reads dbReliabilityBasis');
  }
  return bad.length === 0 || bad.join('; ');
});

check('the typed entry form is gone, and took none of its neighbours with it', () => {
  /* THE DELETION GUARD. refreshAll() is a top-level init statement that used to
     call refreshFamilySelect(); leaving that call behind while removing the
     callee is the exact failure CLAUDE.md records for
     refreshReportWriterOptions — it throws at boot and every top-level
     statement after it silently stops running.

     Asserted both ways: the add-form machinery must be absent, and the things
     that survive it must still be present and wired. */
  const gone = ['refreshFamilySelect', 'makeFamilyName', 'CT_METRIC_OPTIONS', 'ctEntryRowHtml',
                'ctRenumberRows', 'ctRenderRows', 'ctReadRows', 'ctInitEntryRows'];
  const bad = [];
  gone.forEach((n) => { if (new RegExp('\\b' + n + '\\b').test(APP_SRC)) bad.push(n + ' is still referenced'); });
  ['getCustom', 'saveCustom', 'getMergedDB', 'renderDbList', 'refreshAll', 'ctValidateEntry']
    .forEach((n) => { if (!new RegExp('^function ' + n + '\\s*\\(', 'm').test(APP_SRC)) bad.push(n + ' was removed and must not have been'); });
  if (!/^refreshAll\(\);/m.test(APP_SRC)) bad.push('refreshAll is no longer called at top level');
  /* Import is now the ONLY route into the clinical database, so its validation
     and its own handler have to survive. */
  if (!/ct-import/.test(APP_SRC)) bad.push('the import handler is gone, leaving no way to load custom tests');
  if (!/ctValidateEntry\(/.test(APP_SRC)) bad.push('import no longer validates what it saves');
  ['ct-add-family', 'ct-add-subtest', 'ct-entry-body', 'ct-family-select']
    .forEach((id) => { if (HTML_SRC.includes('id="' + id + '"')) bad.push('markup for ' + id + ' is still in index.html'); });
  ['ct-search', 'ct-export', 'ct-import', 'db-list', 'db-thead', 'db-tbody',
   'db-filters', 'db-f-inst', 'db-f-cat', 'db-f-band', 'db-f-basis', 'db-count']
    .forEach((id) => { if (!HTML_SRC.includes('id="' + id + '"')) bad.push('markup for ' + id + ' is missing'); });
  return bad.length === 0 || bad.join('; ');
});

/* ==========================================================================
   33. The single-file bundle

   Psychometric_Assistant.html is committed, and it is the copy that gets
   emailed around, so it is the copy a clinician is most likely to be
   reading. Nothing had ever checked it.

   bundle.ps1 is local-only (.gitignore excludes *.ps1), so it cannot be
   tested here. What CAN be tested is its output, and the output had two
   faults that had been committed and shipped:

     1. It appends " - Bundled" to the document <title> so the standalone
        file is distinguishable from the served page. It did that by
        matching <title> as text, which also hit the ten <title> tags that
        live INSIDE JavaScript strings - nine SVG tooltips in
        app-viz-page.js and the Word export's document title in app.js. So
        hovering a Score Charts row read "...95% CI 92-108 - Bundled", and
        every exported Word document was titled "Assessment Report -
        Bundled". On-screen text is a contract; this was not text anyone
        wrote.

     2. The em dash in that marker was written as CP1252 and re-read as
        UTF-8, so the tab title rendered the three-character mojibake
        instead of a dash.

   The check that catches both, and anything else of the kind, is the
   strongest one available: every inlined source must appear in the bundle
   VERBATIM. Only the BOM is allowed to differ (app.js and styles.css carry
   one; the bundler strips it). That also enforces step 2 of CLAUDE.md's
   shipping list mechanically - a source edited without a rebuild now fails
   here instead of shipping a stale bundle.

   The pre-existing mojibake in styles.css COMMENTS is deliberately not
   caught: it is present in the source too, so the bundler is being
   faithful, and faithfulness is the only thing this section judges.
   ========================================================================== */
heading('33. The single-file bundle');

const BUNDLE_PATH = path.join(ROOT, 'Psychometric_Assistant.html');
const BUNDLE_SRC = fs.readFileSync(BUNDLE_PATH, 'utf8');
const INLINED_SOURCES = ['data.js', 'app.js', 'design-system.js', 'app-effectsize-page.js',
                         'app-viz-page.js', 'styles.css', 'design-system.css'];
// U+00E2 U+20AC U+201D - an em dash encoded CP1252 and decoded UTF-8.
const MANGLED_EM_DASH = 'â€”';

check('every inlined source appears in the bundle verbatim', () => {
  const bad = [];
  INLINED_SOURCES.forEach((f) => {
    const src = fs.readFileSync(path.join(ROOT, f), 'utf8').replace(/^﻿/, '');
    if (!BUNDLE_SRC.includes(src)) bad.push(f + ' does not appear verbatim - the bundle is stale, ' +
                                                'or bundle.ps1 rewrote it on the way in (re-run ./bundle.ps1)');
  });
  return bad.length === 0 || bad.join('; ');
});

check('the "Bundled" marker is on the document title and nowhere else', () => {
  /* One marker, in the <head>. Ten is what the text-matching bug produced;
     zero is fine too, should the marker ever be dropped - what must never
     come back is a marker inside a JS string, where it reaches an SVG
     tooltip or an exported Word document. */
  const total = (BUNDLE_SRC.match(/Bundled/g) || []).length;
  if (total === 0) return true;
  const bad = [];
  if (total > 1) bad.push(total + ' occurrences of "Bundled" - only the document title may carry it');
  const head = BUNDLE_SRC.slice(0, BUNDLE_SRC.indexOf('</head>'));
  const titled = /<title>[^<]*Bundled[^<]*<\/title>/.test(head);
  if (!titled) bad.push('the marker is not on the document title');
  return bad.length === 0 || bad.join('; ');
});

check('the bundle introduces no mojibake of its own', () => {
  /* Scoped to what the bundler ADDS. styles.css ships 138 mangled em dashes
     inside its own comments; those are the source's, not the bundler's, and
     the verbatim check above is what holds them still. */
  const bad = [];
  const marker = BUNDLE_SRC.match(/<title>[^<]*Bundled[^<]*<\/title>/);
  if (marker && marker[0].includes(MANGLED_EM_DASH)) {
    bad.push('the title marker\'s dash is CP1252 read as UTF-8 - bundle.ps1 must write UTF-8');
  }
  const count = (s) => s.split(MANGLED_EM_DASH).length - 1;
  const fromSources = INLINED_SOURCES.concat(['index.html'])
    .reduce((n, f) => n + count(fs.readFileSync(path.join(ROOT, f), 'utf8')), 0);
  const inBundle = count(BUNDLE_SRC);
  if (inBundle > fromSources) {
    bad.push(inBundle + ' mangled dashes in the bundle against ' + fromSources +
             ' in its sources - bundle.ps1 added ' + (inBundle - fromSources));
  }
  return bad.length === 0 || bad.join('; ');
});

// ---------------------------------------------------------------------------
// Summary
// ---------------------------------------------------------------------------
console.log('\n' + '='.repeat(60));
if (failures.length === 0) {
  console.log('ALL ' + passed + ' CHECKS PASSED');
  process.exit(0);
}
console.log(failures.length + ' FAILED, ' + passed + ' passed\n');
failures.forEach(f => console.log('  [' + f.section + '] ' + f.name + '\n      ' + f.detail));
process.exit(1);
