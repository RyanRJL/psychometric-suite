/* =====================================================================
   Profile Analysis - page module (#profile)

   Crawford, Garthwaite & Gault (2007), "Estimating the percentage of the
   population with abnormally low scores (or abnormally large score
   differences) on standardized neuropsychological test batteries: A
   generic method with applications." Neuropsychology, 21, 419-430.

   THE QUESTION. A battery raises one a single row cannot answer: this
   patient has three index scores below the 5th percentile - how unusual
   is that? By definition 5% of the population falls below the 5th
   percentile on any ONE measure. Across the four correlated WAIS-IV
   indices, 13.8% show at least one. Reading rows independently overcalls
   impairment by enough to change a conclusion, which is the paper's own
   argument for the method.

   WHY THIS IS A PAGE AND NOT A VIEW. Every other calculator here is a
   view of data entered elsewhere, and Score Charts is deliberately so.
   This one is not, and the difference is real: the method needs the four
   Index scores specifically, on their own metric, as a complete set.
   Score Tables holds whatever subtests were administered, in any mixture,
   and a profile computed from a partial or subtest-level table would be
   answering a different question from the one the matrix describes. The
   inputs here are therefore typed, and the page carries no other data.

   THE ARITHMETIC IS NOT HERE. profileAbnormality() and the matrix live in
   app.js and data.js, verified in check.js sections 42-43 against the
   paper's Appendix, the exact binomial and all 12 rows of its Table 1.
   This file decides what to ask for and how to say the answer.
   ===================================================================== */
(function(){
  'use strict';

  const section = document.getElementById('profile');
  if (!section) return;

  /* The four WAIS-IV Indices, in the manual's own order. Only Indices:
     the stored matrix covers subtests and process scores too, but a
     profile over a battery the clinician did not administer as a set is
     not a question this page is entitled to answer. */
  const PROF_INDICES = [
    { key:'VCI', label:'Verbal Comprehension', short:'VCI' },
    { key:'PRI', label:'Perceptual Reasoning', short:'PRI' },
    { key:'WMI', label:'Working Memory',       short:'WMI' },
    { key:'PSI', label:'Processing Speed',     short:'PSI' }
  ];

  /* The paper's own criteria, its default in bold on screen. `z` is the
     standard normal deviate for an abnormally LOW score; the paper gives
     -1.645, -1.282 and -1.0 explicitly and the other two follow from the
     same definition. `pct` is the percentile each names, which is what a
     clinician recognises. */
  const PROF_CRITERIA = [
    { id:'sd1',   label:'below 1 SD (15.9th percentile)', z:-1.0,    pct:'15.9th' },
    { id:'p10',   label:'below the 10th percentile',      z:-1.282,  pct:'10th' },
    { id:'p5',    label:'below the 5th percentile',       z:-1.645,  pct:'5th', isDefault:true },
    { id:'p2',    label:'below the 2nd percentile',       z:-2.054,  pct:'2nd' },
    { id:'p1',    label:'below the 1st percentile',       z:-2.326,  pct:'1st' }
  ];
  /* Differences are two-tailed - a large difference in either direction is
     the finding - so this is the |z| beyond which a difference is called
     abnormal, not a low-score criterion. The paper uses 1.960 throughout. */
  const PROF_DIFF_Z = 1.960;
  const PROF_TRIALS = 200000;

  const profState = { criterion:'p5', scores:{} };

  /* The simulation depends on the criterion and nothing else - scores do
     not enter it - so a keystroke in a score box must not re-run it.
     Cached per criterion, which makes typing free. */
  const profCache = {};

  function profEl(id){ return document.getElementById(id); }
  function profNum(v){
    if (v === null || v === undefined || v === '') return null;
    const n = parseFloat(v);
    return Number.isFinite(n) ? n : null;
  }
  function profCriterion(){
    return PROF_CRITERIA.find(c => c.id === profState.criterion) || PROF_CRITERIA[2];
  }

  /* The four-Index submatrix, read out of the shipped WAIS4_INTERCORR at
     the moment it is needed rather than copied here. One source: a second
     literal would be a second thing to keep in step with Table 5.1. */
  function profMatrix(){
    if (typeof WAIS4_INTERCORR === 'undefined') return null;
    const M = WAIS4_INTERCORR;
    const pos = {};
    M.order.forEach((n, i) => { pos[n] = i; });
    const keys = PROF_INDICES.map(x => x.key);
    const out = [];
    for (let i = 0; i < keys.length; i++){
      const row = [];
      for (let j = 0; j < keys.length; j++){
        if (i === j){ row.push(1); continue; }
        const a = keys[i], b = keys[j];
        const v = pos[a] > pos[b] ? M.r[a + '|' + b] : M.r[b + '|' + a];
        if (!Number.isFinite(v)) return null;
        row.push(v);
      }
      out.push(row);
    }
    return out;
  }

  function profSimulate(){
    const c = profCriterion();
    if (profCache[c.id]) return profCache[c.id];
    const R = profMatrix();
    if (!R || typeof profileAbnormality !== 'function') return null;
    const res = profileAbnormality(R, { lowZ:c.z, diffZ:PROF_DIFF_Z, trials:PROF_TRIALS });
    if (res) profCache[c.id] = res;
    return res;
  }

  /* The patient's own counts, on exactly the definitions the simulation
     used. Anything computed differently here would report a count against
     a percentage that answers a different question. */
  function profCounts(){
    const R = profMatrix();
    if (!R) return null;
    const c = profCriterion();
    const z = [], present = [];
    PROF_INDICES.forEach((idx, i) => {
      const v = profNum(profState.scores[idx.key]);
      present[i] = v !== null;
      z[i] = v === null ? null : (v - 100) / 15;
    });
    if (!present.every(Boolean)) return { complete:false };

    let low = 0;
    const lowWhich = [];
    z.forEach((v, i) => { if (v < c.z){ low++; lowWhich.push(PROF_INDICES[i].short); } });

    let pair = 0;
    const pairWhich = [];
    for (let i = 0; i < 4; i++){
      for (let j = i + 1; j < 4; j++){
        const t = PROF_DIFF_Z * Math.sqrt(2 - 2 * R[i][j]);
        if (Math.abs(z[i] - z[j]) > t){
          pair++;
          pairWhich.push(PROF_INDICES[i].short + '–' + PROF_INDICES[j].short);
        }
      }
    }

    const means = (typeof profileMatrixMeans === 'function') ? profileMatrixMeans(R) : null;
    let dev = 0;
    const devWhich = [];
    if (means){
      const mean = z.reduce((a, b) => a + b, 0) / 4;
      z.forEach((v, i) => {
        const t = PROF_DIFF_Z * Math.sqrt(1 + means.grandMean - 2 * means.rowMean[i]);
        if (Math.abs(mean - v) > t){ dev++; devWhich.push(PROF_INDICES[i].short); }
      });
    }
    return { complete:true, low, lowWhich, pair, pairWhich, dev, devWhich };
  }

  /* "j or more" for a count of zero is the whole population, so it is not
     a finding and must not be printed as one - "100% of people show 0 or
     more" is true and useless. Zero gets its own sentence. */
  function profPctFor(series, j){
    if (!series || j < 1) return null;
    const v = series[j - 1];
    return Number.isFinite(v) ? v : null;
  }
  function profFmtPct(v){
    if (v === null) return '—';
    if (v > 0 && v < 0.01) return '< 0.01%';
    return (v >= 10 ? v.toFixed(1) : v.toFixed(2)) + '%';
  }

  function profInputsHtml(){
    return PROF_INDICES.map(idx =>
      '<div class="field">'
      + '<label>' + idx.label + ' <span class="hint">' + idx.short + '</span></label>'
      + '<input type="number" class="prof-score-input" data-prof-index="' + idx.key + '" '
      + 'min="40" max="160" step="1" inputmode="numeric" placeholder="e.g. 100" '
      + 'aria-label="' + idx.label + ' Index score">'
      + '</div>'
    ).join('');
  }

  function profCriterionHtml(){
    return PROF_CRITERIA.map(c =>
      '<option value="' + c.id + '"' + (c.id === profState.criterion ? ' selected' : '') + '>'
      + c.label + '</option>'
    ).join('');
  }

  /* One card per question. The patient's count is the headline; the
     percentage under it is what makes the count mean something. With no
     scores entered the card still states the population figure for j = 1,
     because that is the number a clinician most often wants and it needs
     no patient data at all. */
  function profResultCard(title, tip, count, series, unit){
    const pct = count === null ? profPctFor(series, 1) : profPctFor(series, count);
    let value, desc;
    if (count === null){
      value = profFmtPct(pct);
      desc = 'of the population show one or more. Enter all four Index scores for this patient’s own count.';
    } else if (count === 0){
      value = '0';
      desc = 'None of this patient’s ' + unit + ' meet the criterion.';
    } else {
      value = String(count);
      desc = profFmtPct(pct) + ' of the population show ' + count + ' or more.';
    }
    return '<div class="prof-card" data-tooltip="' + escapeHtml(tip) + '">'
      + '<div class="prof-card-label">' + escapeHtml(title) + '</div>'
      + '<div class="prof-card-value">' + value + '</div>'
      + '<div class="prof-card-desc">' + desc + '</div>'
      + '</div>';
  }

  function profWhichHtml(which){
    if (!which || !which.length) return '';
    return '<div class="prof-which">' + escapeHtml(which.join(', ')) + '</div>';
  }

  /* The full distribution, so the clinician can see how steeply the
     percentages fall rather than only the one cell their patient landed
     in. The patient's own column is marked, not extracted. */
  function profTableHtml(res, counts){
    if (!res) return '';
    const rows = [
      { label:'Abnormally low Index scores', series:res.lowScores, n:4, got:counts && counts.complete ? counts.low : null },
      { label:'Abnormally large pairwise differences', series:res.pairwise, n:6, got:counts && counts.complete ? counts.pair : null },
      { label:'Abnormally large deviations from own mean', series:res.deviations, n:4, got:counts && counts.complete ? counts.dev : null }
    ];
    const maxJ = 6;
    let head = '<tr><th class="prof-th-label">Percentage of the population showing…</th>';
    for (let j = 1; j <= maxJ; j++) head += '<th>' + j + ' or more</th>';
    head += '</tr>';
    let body = '';
    for (const r of rows){
      body += '<tr><td class="prof-td-label">' + escapeHtml(r.label) + '</td>';
      for (let j = 1; j <= maxJ; j++){
        if (j > r.n){ body += '<td class="prof-td-na">—</td>'; continue; }
        const cls = (r.got !== null && r.got === j) ? ' class="prof-td-hit"' : '';
        body += '<td' + cls + '>' + profFmtPct(profPctFor(r.series, j)) + '</td>';
      }
      body += '</tr>';
    }
    return '<table class="prof-table"><thead>' + head + '</thead><tbody>' + body + '</tbody></table>';
  }

  /* THE MONTE CARLO ERROR IS SHOWN, not hidden behind decimal places. The
     engine reports sqrt(p(1-p)/n), the exact sampling error of a
     proportion; printing more digits than that supports would imply a
     precision the method does not have. */
  function profPrecisionHtml(res){
    if (!res || typeof profileAbnormalityStdErr !== 'function') return '';
    const se = profileAbnormalityStdErr(res.lowScores[0], res.trials);
    if (!Number.isFinite(se)) return '';
    return '<p class="prof-precision">Percentages are Monte Carlo estimates from '
      + res.trials.toLocaleString() + ' simulated cases, so they carry a sampling error of about '
      + '±' + se.toFixed(2) + ' percentage points at the largest value in this table. '
      + 'The simulation is seeded, so the same profile always returns the same figures.</p>';
  }

  function renderProfile(){
    const res = profSimulate();
    const counts = profCounts();
    const c = profCriterion();
    const out = profEl('prof-results');
    if (!out) return;

    if (!res){
      out.innerHTML = '<div class="prof-empty">The WAIS-IV intercorrelation matrix is unavailable, so no profile can be computed.</div>';
      renderProfileApa();
      return;
    }
    const complete = counts && counts.complete;
    let html = '<div class="prof-cards">'
      + profResultCard('Abnormally low Index scores',
          'Index scores ' + c.label + '. Across four correlated Indices, showing one is far more common than the criterion alone suggests.',
          complete ? counts.low : null, res.lowScores, 'Index scores')
      + profResultCard('Abnormal pairwise differences',
          'Differences between any two of the four Indices larger than 95% of the population shows, regardless of direction.',
          complete ? counts.pair : null, res.pairwise, 'Index pairs')
      + profResultCard('Abnormal deviations from own mean',
          'Index scores differing from this patient’s own four-Index mean by more than 95% of the population does.',
          complete ? counts.dev : null, res.deviations, 'Index scores')
      + '</div>';
    if (complete && (counts.lowWhich.length || counts.pairWhich.length || counts.devWhich.length)){
      html += '<div class="prof-which-row">'
        + '<div><span class="prof-which-label">Low</span>' + (counts.lowWhich.length ? profWhichHtml(counts.lowWhich) : '<div class="prof-which prof-which-none">none</div>') + '</div>'
        + '<div><span class="prof-which-label">Differences</span>' + (counts.pairWhich.length ? profWhichHtml(counts.pairWhich) : '<div class="prof-which prof-which-none">none</div>') + '</div>'
        + '<div><span class="prof-which-label">Deviations</span>' + (counts.devWhich.length ? profWhichHtml(counts.devWhich) : '<div class="prof-which prof-which-none">none</div>') + '</div>'
        + '</div>';
    }
    html += profTableHtml(res, counts) + profPrecisionHtml(res);
    out.innerHTML = html;
    renderProfileApa();
  }

  /* THE APA TABLE IS THE PATIENT'S, NOT THE POPULATION'S, and it is
     emitted only once all four scores are present.

     The Working Report collects any container that holds an .apa-table, so
     a renderer that writes its shell unconditionally offers an empty table
     from every page - the defect CLAUDE.md records for renderOpiePredictApa.
     A population table needs no patient data at all, so emitting one would
     put a table in the report for a patient who has not been assessed. */
  function renderProfileApa(){
    const out = profEl('prof-apa');
    if (!out) return;
    const counts = profCounts();
    const res = profSimulate();
    if (!counts || !counts.complete || !res){
      out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">'
        + 'Enter all four Index scores to generate the APA table.</div>';
      return;
    }
    const c = profCriterion();
    const rows = [
      ['Index scores ' + c.label, counts.low, profPctFor(res.lowScores, counts.low)],
      ['Abnormally large pairwise differences', counts.pair, profPctFor(res.pairwise, counts.pair)],
      ['Abnormally large deviations from own mean', counts.dev, profPctFor(res.deviations, counts.dev)]
    ];
    let body = '';
    for (const [label, n, pct] of rows){
      body += '<tr><td>' + escapeHtml(label) + '</td><td class="num">' + n + '</td>'
        + '<td class="num">' + (n === 0 ? '—' : profFmtPct(pct)) + '</td></tr>';
    }
    const scoreLine = PROF_INDICES
      .map(i => i.short + ' ' + profNum(profState.scores[i.key]))
      .join(', ');
    out.innerHTML =
      '<div class="apa-table-num">Table 1</div>'
      + '<div class="apa-table-title">Number of abnormal WAIS-IV Index findings and their base rates</div>'
      + '<table class="apa-table"><thead><tr>'
      + '<th>Finding</th><th class="num">Number observed</th><th class="num">Base rate</th>'
      + '</tr></thead><tbody>' + body + '</tbody></table>'
      + (typeof apaNoteHtml === 'function'
          ? apaNoteHtml('prof', { criterion:c.label, pct:c.pct, trials:res.trials, scores:scoreLine })
          : '');
  }

  /* Wiring. Inputs are read into profState and re-rendered in place: the
     results block is a sibling of the inputs, never their parent, so a
     re-render cannot take the focused box out from under the clinician. */
  function setupProfile(){
    const inputs = profEl('prof-inputs');
    if (inputs) inputs.innerHTML = profInputsHtml();
    const sel = profEl('prof-criterion');
    if (sel) sel.innerHTML = profCriterionHtml();

    section.addEventListener('input', e => {
      const t = e.target;
      if (t && t.dataset && t.dataset.profIndex){
        profState.scores[t.dataset.profIndex] = t.value;
        renderProfile();
      }
    });
    section.addEventListener('change', e => {
      if (e.target && e.target.id === 'prof-criterion'){
        profState.criterion = e.target.value;
        renderProfile();
      }
    });
    const clear = profEl('prof-clear');
    if (clear){
      clear.addEventListener('click', () => {
        profState.scores = {};
        section.querySelectorAll('.prof-score-input').forEach(el => { el.value = ''; });
        renderProfile();
      });
    }
    renderProfile();
  }

  setupProfile();
  window.renderProfile = renderProfile;
})();
