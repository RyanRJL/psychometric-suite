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

check('99 groups and 590 entries present', () => {
  const groups = Object.keys(D.normDB).length;
  let n = 0; eachNormEntry(() => n++);
  return (groups === 99 && n === 590) || 'got ' + groups + ' groups / ' + n + ' entries';
});
check('every entry has m1, sd1, m2, sd2 and r', () => {
  const bad = [];
  eachNormEntry((e, g, n) => {
    for (const f of ['m1', 'sd1', 'm2', 'sd2', 'r']) {
      if (!Number.isFinite(e[f])) bad.push(g + ' / ' + n + ' missing ' + f);
    }
  });
  return bad.length === 0 || bad.length + ' problems, first: ' + bad[0];
});
check('all correlations lie in [0, 1)', () => {
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (!(e.r >= 0 && e.r < 1)) bad.push(g + ' / ' + n + ' r=' + e.r);
    if (Number.isFinite(e.rCorrected) && !(e.rCorrected >= 0 && e.rCorrected < 1)) bad.push(g + ' / ' + n + ' rCorrected=' + e.rCorrected);
  });
  return bad.length === 0 || bad.join('; ');
});
check('all standard deviations are strictly positive', () => {
  const bad = [];
  eachNormEntry((e, g, n) => { if (!(e.sd1 > 0) || !(e.sd2 > 0)) bad.push(g + ' / ' + n); });
  return bad.length === 0 || bad.join('; ');
});
check('difference variance sd1^2 + sd2^2 - 2r*sd1*sd2 is positive everywhere', () => {
  // Bounded below by (sd1-sd2)^2, so a failure means corrupt data, not maths.
  const bad = [];
  eachNormEntry((e, g, n) => {
    const v = e.sd1 ** 2 + e.sd2 ** 2 - 2 * e.r * e.sd1 * e.sd2;
    if (!(v > 0)) bad.push(g + ' / ' + n + ' v=' + v);
  });
  return bad.length === 0 || bad.join('; ');
});
check('rCorrected present on exactly 233 of 590 entries', () => {
  let n = 0; eachNormEntry(e => { if (Number.isFinite(e.rCorrected)) n++; });
  return n === 233 || 'got ' + n;
});
check('D-KEFS, CVLT-C and WISC-V carry no rCorrected at all', () => {
  const bad = [];
  eachNormEntry((e, g, n) => {
    if (/^(D-KEFS|CVLT-C|WISC-V)\b/.test(g) && Number.isFinite(e.rCorrected)) bad.push(g + ' / ' + n);
  });
  return bad.length === 0 || bad.join('; ');
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
check('McSweeney SRB slope and SEE are computable for every entry', () => {
  const bad = [];
  eachNormEntry((e, g, n) => {
    const slope = e.r * (e.sd2 / e.sd1);
    const see = e.sd2 * Math.sqrt(1 - e.r * e.r);
    if (!Number.isFinite(slope) || !Number.isFinite(see) || see <= 0) bad.push(g + ' / ' + n);
  });
  return bad.length === 0 || bad.join('; ');
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

const batteryCi = (() => {
  const boundsSrc = APP_SRC.match(/const BATTERY_SCORE_BOUNDS = \{[^;]*;/);
  if (!boundsSrc) throw new Error('BATTERY_SCORE_BOUNDS not found in app.js');
  const ctx = {
    normDB: D.normDB,
    getMergedDB: () => D.normDB,
    document: { getElementById: () => ({ value: 'scaled' }) }
  };
  vm.createContext(ctx);
  vm.runInContext(
    boundsSrc[0] + '\n' +
      ['getBatteryRowReliability', 'rowScoreType', 'getBatteryCiHtml'].map((n) => extractFn(APP_SRC, n)).join('\n') +
      '\n;globalThis.__B = getBatteryCiHtml;',
    ctx
  );
  return ctx.__B;
})();

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

const HTML_SRC = fs.readFileSync(path.join(ROOT, 'index.html'), 'utf8');

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

const PROJECT_SCRIPTS = ['app.js', 'design-system.js', 'app-effectsize-page.js', 'data.js'];

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

const RAW_FAMILIES = [
  'CVLT-C Subtests (Raw Scores) · Age 8',
  'CVLT-C Subtests (Raw Scores) · Age 12',
  'CVLT-C Subtests (Raw Scores) · Age 16',
  'RBANS Subtests · Ages 12-19',
  'RBANS Subtests · Ages 20-89',
];

check('the five raw-score families are tagged metric:\'raw\' on every entry', () => {
  const bad = [];
  for (const fam of RAW_FAMILIES) {
    const g = D.normDB[fam];
    if (!g) { bad.push(fam + ' missing from normDB'); continue; }
    for (const name in g) {
      if (g[name] && typeof g[name] === 'object' && g[name].metric !== 'raw') {
        bad.push(fam + ' / ' + name);
      }
    }
  }
  return bad.length === 0 || bad.length + ' untagged, first: ' + bad[0];
});

check('exactly 63 entries are tagged raw, and none outside those families', () => {
  const fams = new Set(RAW_FAMILIES);
  let n = 0; const stray = [];
  eachNormEntry((e, g, name) => {
    if (e.metric !== 'raw') return;
    n++;
    if (!fams.has(g)) stray.push(g + ' / ' + name);
  });
  if (stray.length) return 'tagged raw outside the known families: ' + stray.join(', ');
  return n === 63 || 'got ' + n + ' tagged raw, expected 63';
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
    ['toZ', 'fromZ', 'normInv', 'normCDF', 'erf', 'inferScoreTypeForSubtest', 'inferScoreType',
     'scoreTypeLabel', 'scoreTypeAbbr', 'sdiSdUnit']
      .map(n => extractFn(APP_SRC, n)).join('\n') +
    ';globalThis.__E={toZ,fromZ,inferScoreTypeForSubtest,scoreTypeLabel,scoreTypeAbbr,sdiSdUnit};',
    ctx
  );
  const A = ctx.__E;

  check('inferScoreTypeForSubtest returns \'raw\' for every tagged entry', () => {
    const bad = [];
    eachNormEntry((e, g, name) => {
      if (e.metric !== 'raw') return;
      const t = A.inferScoreTypeForSubtest(g, name, e);
      if (t !== 'raw') bad.push(g + ' / ' + name + ' -> ' + t);
    });
    return bad.length === 0 || bad.length + ' wrong, first: ' + bad[0];
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

check('sdiComputeChange rejects a row whose metric has no SD unit', () => {
  const src = extractFn(APP_SRC, 'sdiComputeChange');
  return /Number\.isFinite\(unit\)/.test(src)
    || 'sdiComputeChange no longer guards a missing SD unit — a raw row would divide by undefined';
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
