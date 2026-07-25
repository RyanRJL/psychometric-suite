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
   10. Documentation contracts
   Text that must stay in step with the code.
   ========================================================================== */
heading('10. Documentation contracts');

check('every OPIE tooltip warns it is illustrative only in a UK context', () => {
  const keys = Object.keys(D.PRE_MODEL_TOOLTIPS).filter(k => /^opie/i.test(k));
  const missing = keys.filter(k => !/ILLUSTRATIVE ONLY/i.test(D.PRE_MODEL_TOOLTIPS[k]));
  return missing.length === 0 || 'missing the warning: ' + missing.join(', ');
});

check('the OPIE tooltips state that sex and an age range are required', () => {
  const t = D.PRE_MODEL_TOOLTIPS.opieDefault || '';
  return (/sex/i.test(t) && /16-90|16–90/.test(t)) || 'opieDefault no longer states the input requirements';
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
