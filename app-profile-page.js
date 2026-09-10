/* =====================================================================
   Profile Analysis - page module (#profile)

   Crawford, Garthwaite & Gault (2007), Neuropsychology, 21, 419-430.

   THE QUESTION. A battery raises one a single row cannot answer: this
   patient has three scores below the 5th percentile - how unusual is
   that? By definition 5% of the population falls below the 5th
   percentile on any ONE measure. Across four correlated WAIS-IV Indices
   13.8% show at least one; across ten core subtests, 24.9%. Reading rows
   independently overcalls impairment by enough to change a conclusion.

   ANY SET OF MEASURES, NOT A FIXED BATTERY. An earlier version offered
   three canned sets on the belief that the method needed a complete
   battery. It does not, and the paper is explicit three times over: the
   program "requires the user to specify the number of tests in the
   battery (up to a maximum of 20 tests)... and enter the correlation
   between tests in the form of a lower triangular matrix"; the
   conclusion speaks of abnormality "from among the overall set of tests
   administered"; and Crawford's own supplementary programs state that
   "the methods can be applied when only a subset of the Index scores
   have been administered". The requirement is only that R covers the
   measures actually in hand.

   So the clinician picks what was administered and R is the sub-matrix
   over exactly those. The one thing that cannot be done is quote a
   ten-measure figure for seven measures - which is why every percentage
   on screen is recomputed for the current selection and nothing is
   carried over.

   THE ARITHMETIC IS NOT HERE. profileAbnormality() and the matrix live
   in app.js and data.js, verified in check.js sections 42-43 against the
   paper's Appendix, the exact binomial and all 12 rows of its Table 1.
   ===================================================================== */
(function(){
  'use strict';

  const section = document.getElementById('profile');
  if (!section) return;

  /* ---------- what may be profiled together ----------

     A profile must never hold a measure and a piece of itself. Two
     measures conflict when their COMPONENT SETS intersect, which catches
     three distinct cases with one rule:

       - a composite and its own subtest   (VCI with Vocabulary)
       - two composites sharing subtests   (FSIQ with VCI)
       - a subtest and its process scores  (Digit Span with Digit Span
                                            Forward, which is part of it)

     Crawford's own programs analyse Index scores OR subtests for exactly
     this reason. Expansion is recursive: WMI holds Digit Span, which
     itself holds the three Digit Span process scores. */
  const PROF_COMPOSED_OF = {
    VCI:  ['SI', 'VC', 'IN', 'CO'],
    PRI:  ['BD', 'MR', 'VP', 'FW', 'PCm'],
    WMI:  ['DS', 'AR', 'LN'],
    PSI:  ['SS', 'CD', 'CA'],
    FSIQ: ['BD', 'SI', 'DS', 'MR', 'VC', 'AR', 'SS', 'VP', 'IN', 'CD'],
    DS:   ['DSF', 'DSB', 'DSS']
  };
  /* Block Design No Time Bonus is the SAME administration rescored - the
     matrix puts them at r = .97 - so they are one measure for this
     purpose even though neither contains the other. */
  const PROF_ALIAS = { BDN: 'BD' };

  function profComponents(key, seen){
    const k = PROF_ALIAS[key] || key;
    const out = seen || new Set();
    const parts = PROF_COMPOSED_OF[k];
    if (!parts){ out.add(k); return out; }
    parts.forEach(p => profComponents(p, out));
    return out;
  }
  function profConflicts(a, b){
    if (a === b) return true;
    const A = profComponents(a), B = profComponents(b);
    for (const x of A) if (B.has(x)) return true;
    return false;
  }

  /* ---------- the roster on screen ----------
     Grouped the way a clinician reads a record form, and each group's
     members are exactly the matrix's own keys. */
  const PROF_GROUPS = [
    { id:'composites', label:'Composites',    keys:['FSIQ', 'VCI', 'PRI', 'WMI', 'PSI'] },
    { id:'core',       label:'Core subtests', keys:['BD', 'SI', 'DS', 'MR', 'VC', 'AR', 'SS', 'VP', 'IN', 'CD'] },
    { id:'supp',       label:'Supplementary', keys:['LN', 'FW', 'CO', 'CA', 'PCm'] },
    { id:'process',    label:'Process scores', keys:['BDN', 'DSF', 'DSB', 'DSS'] }
  ];
  const PROF_COMPOSITES = ['FSIQ', 'VCI', 'PRI', 'WMI', 'PSI'];

  /* Composites are M 100 / SD 15; everything else is a scaled score,
     M 10 / SD 3. Per measure, not per page: a single conversion would
     score a scaled 8 as z = -6.1 and call every subtest catastrophic. */
  function profMetric(key){
    return PROF_COMPOSITES.indexOf(key) !== -1
      ? { mean:100, sd:15, label:'Index', min:40, max:160 }
      : { mean:10,  sd:3,  label:'scaled', min:1,  max:19 };
  }

  const PROF_CRITERIA = [
    { id:'sd1', label:'below 1 SD (15.9th percentile)', z:-1.0,   pct:'15.9th' },
    { id:'p10', label:'below the 10th percentile',      z:-1.282, pct:'10th' },
    { id:'p5',  label:'below the 5th percentile',       z:-1.645, pct:'5th' },
    { id:'p2',  label:'below the 2nd percentile',       z:-2.054, pct:'2nd' },
    { id:'p1',  label:'below the 1st percentile',       z:-2.326, pct:'1st' }
  ];
  const PROF_DIFF_Z = 1.960;   // two-tailed; a large difference either way is the finding
  const PROF_TRIALS = 200000;
  const PROF_MIN = 2;          // one measure is not a profile

  const profState = { criterion:'p5', selected:[], scores:{}, pulled:{} };

  /* Keyed on the SELECTION and the criterion. Scores do not enter the
     simulation, so typing a score must never re-run it - and the run
     costs ~200ms at 15 measures, comparing 105 pairs per simulated case. */
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
  function profMatrixObj(){
    return (typeof WAIS4_INTERCORR !== 'undefined') ? WAIS4_INTERCORR : null;
  }
  function profLabel(key){
    const M = profMatrixObj();
    return (M && M.labels && M.labels[key]) || key;
  }
  function profRestricted(){
    const M = profMatrixObj();
    const list = (M && M.restrictedTo16_69) || [];
    return profState.selected.filter(k => list.indexOf(k) !== -1);
  }
  /* Any selected measure that is not a composite is a scaled score, and
     the paper's own limitation applies to those - see profCaveatHtml. */
  function profHasScaled(){
    return profState.selected.some(k => PROF_COMPOSITES.indexOf(k) === -1);
  }
  function profDisabledBy(key){
    for (const s of profState.selected){
      if (s !== key && profConflicts(s, key)) return s;
    }
    return null;
  }

  /* The sub-matrix over the CURRENT selection, read out of the shipped
     matrix at the moment it is needed. */
  function profMatrix(){
    const M = profMatrixObj();
    if (!M) return null;
    const pos = {};
    M.order.forEach((n, i) => { pos[n] = i; });
    const keys = profState.selected;
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
    if (profState.selected.length < PROF_MIN) return null;
    const key = profState.selected.join(',') + '|' + profCriterion().id;
    if (profCache[key]) return profCache[key];
    const R = profMatrix();
    if (!R || typeof profileAbnormality !== 'function') return null;
    const c = profCriterion();
    const res = profileAbnormality(R, { lowZ:c.z, diffZ:PROF_DIFF_Z, trials:PROF_TRIALS });
    if (res) profCache[key] = res;
    return res;
  }

  /* The patient's own counts, on exactly the definitions the simulation
     used - anything computed differently here would report a count
     against a percentage answering a different question. */
  function profCounts(){
    const keys = profState.selected;
    if (keys.length < PROF_MIN) return { complete:false, entered:0, total:keys.length };
    const R = profMatrix();
    if (!R) return null;
    const c = profCriterion();
    const z = [], present = [];
    keys.forEach((k, i) => {
      const v = profNum(profState.scores[k]);
      const m = profMetric(k);
      present[i] = v !== null;
      z[i] = v === null ? null : (v - m.mean) / m.sd;
    });
    const entered = present.filter(Boolean).length;
    if (entered < keys.length) return { complete:false, entered, total:keys.length };

    let low = 0; const lowWhich = [];
    z.forEach((v, i) => { if (v < c.z){ low++; lowWhich.push(profLabel(keys[i])); } });

    let pair = 0; const pairWhich = [];
    for (let i = 0; i < keys.length; i++){
      for (let j = i + 1; j < keys.length; j++){
        const t = PROF_DIFF_Z * Math.sqrt(2 - 2 * R[i][j]);
        if (Math.abs(z[i] - z[j]) > t){
          pair++;
          pairWhich.push(profLabel(keys[i]) + ' – ' + profLabel(keys[j]));
        }
      }
    }

    const means = (typeof profileMatrixMeans === 'function') ? profileMatrixMeans(R) : null;
    let dev = 0; const devWhich = [];
    if (means){
      const mean = z.reduce((a, b) => a + b, 0) / keys.length;
      z.forEach((v, i) => {
        const t = PROF_DIFF_Z * Math.sqrt(1 + means.grandMean - 2 * means.rowMean[i]);
        if (Math.abs(mean - v) > t){ dev++; devWhich.push(profLabel(keys[i])); }
      });
    }
    return { complete:true, entered, total:keys.length, low, lowWhich, pair, pairWhich, dev, devWhich };
  }

  /* ---------- pulling scores from Score Tables ----------

     PREFILL, NEVER BINDING. Score Tables holds whatever was administered
     in whatever mixture; this page needs a coherent selection. So a
     pulled value fills the box, selects the measure and is labelled as
     pulled - and the clinician can change or deselect it. Nothing here
     writes back.

     Matched on the measure NAME, which is the same string in both files:
     24 of the 28 WAIS-IV names in normDB are Table 5.1 labels verbatim.
     The four that are not are the Longest Span base-rate measures, which
     have no correlations and correctly cannot be profiled. */
  function profScoreTableRows(){
    if (typeof batteryRows === 'undefined' || !Array.isArray(batteryRows)) return {};
    const M = profMatrixObj();
    if (!M) return {};
    const byName = {};
    M.order.forEach(k => { byName[M.labels[k]] = k; });
    const out = {};
    batteryRows.forEach(r => {
      if (!r || r.isExample) return;
      const key = byName[r.name];
      if (!key) return;
      const v = profNum(r.score);
      if (v === null) return;
      /* A group key naming another instrument cannot be a WAIS-IV
         measure however the row is titled. */
      const group = (typeof batteryGroupKeyOf === 'function') ? batteryGroupKeyOf(r) : '';
      if (group && !/WAIS-IV/.test(group)) return;
      out[key] = v;
    });
    return out;
  }

  function profPull(){
    const found = profScoreTableRows();
    const keys = Object.keys(found);
    if (!keys.length) return 0;
    let added = 0;
    keys.forEach(k => {
      if (profState.selected.indexOf(k) !== -1){
        profState.scores[k] = String(found[k]);
        profState.pulled[k] = true;
        return;
      }
      /* Skip anything that would conflict with what is already chosen -
         Score Tables may legitimately hold both an Index and its
         subtests, and this page may not. */
      if (profDisabledBy(k)) return;
      profState.selected.push(k);
      profState.scores[k] = String(found[k]);
      profState.pulled[k] = true;
      added++;
    });
    profSortSelection();
    return added;
  }
  function profSortSelection(){
    const M = profMatrixObj();
    if (!M) return;
    const pos = {};
    M.order.forEach((n, i) => { pos[n] = i; });
    profState.selected.sort((a, b) => pos[a] - pos[b]);
  }

  /* ---------- rendering ---------- */
  function profFmtPct(v){
    if (v === null || v === undefined) return '—';
    if (v > 0 && v < 0.01) return '< 0.01%';
    return (v >= 10 ? v.toFixed(1) : v.toFixed(2)) + '%';
  }
  function profPctFor(series, j){
    if (!series || j < 1) return null;
    const v = series[j - 1];
    return Number.isFinite(v) ? v : null;
  }

  function profPickerHtml(){
    const found = profScoreTableRows();
    return PROF_GROUPS.map(g => {
      const rows = g.keys.map(k => {
        const on = profState.selected.indexOf(k) !== -1;
        const blockedBy = on ? null : profDisabledBy(k);
        const m = profMetric(k);
        const cls = 'prof-pick' + (on ? ' is-on' : '') + (blockedBy ? ' is-blocked' : '');
        const title = blockedBy
          ? 'Cannot be profiled alongside ' + profLabel(blockedBy) + ' — one contains the other.'
          : '';
        return '<label class="' + cls + '"' + (title ? ' title="' + escapeHtml(title) + '"' : '') + '>'
          + '<input type="checkbox" class="prof-pick-box" data-prof-key="' + k + '"'
          + (on ? ' checked' : '') + (blockedBy ? ' disabled' : '') + '>'
          + '<span class="prof-pick-name">' + escapeHtml(profLabel(k)) + '</span>'
          + (found[k] !== undefined && !on ? '<span class="prof-pick-avail" title="On Score Tables">•</span>' : '')
          + '<input type="number" class="prof-pick-score" data-prof-score="' + k + '" '
          + 'min="' + m.min + '" max="' + m.max + '" step="1" inputmode="numeric" '
          + 'value="' + escapeHtml(profState.scores[k] === undefined ? '' : profState.scores[k]) + '" '
          + (on ? '' : 'disabled ') + 'aria-label="' + escapeHtml(profLabel(k)) + ' score">'
          + (profState.pulled[k] && on ? '<span class="prof-pick-pulled" title="Pulled from Score Tables">↩</span>' : '')
          + '</label>';
      }).join('');
      return '<div class="prof-group"><div class="prof-group-head">' + escapeHtml(g.label) + '</div>'
        + '<div class="prof-group-rows">' + rows + '</div></div>';
    }).join('');
  }

  function profCardHtml(title, tip, count, series, unit){
    const pct = count === null ? profPctFor(series, 1) : profPctFor(series, count);
    let value, desc;
    if (count === null){
      value = profFmtPct(pct);
      desc = 'of the population show one or more.';
    } else if (count === 0){
      value = '0';
      desc = 'None of these ' + unit + ' meet the criterion.';
    } else {
      value = String(count);
      desc = profFmtPct(pct) + ' of the population show ' + count + ' or more.';
    }
    return '<div class="prof-card" data-tooltip="' + escapeHtml(tip) + '">'
      + '<div class="prof-card-label">' + escapeHtml(title) + '</div>'
      + '<div class="prof-card-value">' + value + '</div>'
      + '<div class="prof-card-desc">' + desc + '</div></div>';
  }

  function profWhichHtml(which){
    if (!which || !which.length) return '<div class="prof-which prof-which-none">none</div>';
    return '<div class="prof-which">' + escapeHtml(which.join(', ')) + '</div>';
  }

  function profTableHtml(res, counts){
    if (!res) return '';
    const k = profState.selected.length;
    const pairs = k * (k - 1) / 2;
    const rows = [
      { label:'Abnormally low scores', series:res.lowScores, n:k, got:counts && counts.complete ? counts.low : null },
      { label:'Abnormally large pairwise differences', series:res.pairwise, n:pairs, got:counts && counts.complete ? counts.pair : null },
      { label:'Abnormally large deviations from own mean', series:res.deviations, n:k, got:counts && counts.complete ? counts.dev : null }
    ];
    /* Capped: 15 measures give 105 possible pairs and a 105-column table
       is unreadable. The cards carry this patient's own count whatever
       j they landed on, so the table is context, not the answer. */
    const CAP = 8;
    const maxJ = Math.min(CAP, Math.max.apply(null, rows.map(r => r.n)));
    let head = '<tr><th class="prof-th-label">Percentage of the population showing…</th>';
    for (let j = 1; j <= maxJ; j++) head += '<th>' + j + '+</th>';
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

  /* THE PAPER'S OWN LIMITATION, and it lands on every subtest profile.
     Multivariate normality assumes CONTINUOUS scores: "if the tests in a
     battery have a limited number of possible raw scores... or scaled
     scores, the accuracy of the estimates will suffer. For example, in
     contrast to Wechsler Index scores, Wechsler subtest scaled scores
     have a limited range of scores (the distances between scaled score
     points represent one third of a standard deviation)."

     So the author says a subtest-level profile is less accurate than an
     Index-level one. On-screen text is a contract here; the page cannot
     offer subtests and stay silent about it. */
  function profCaveatHtml(){
    const bits = [];
    if (profHasScaled()){
      bits.push('<strong>Scaled scores are coarse.</strong> The method assumes continuous scores. '
        + 'Crawford, Garthwaite &amp; Gault note that with a limited range of scaled scores — a scaled point is a third of a standard deviation — '
        + 'the accuracy of these estimates suffers, in contrast to Index scores. Treat a subtest-level profile as the rougher reading.');
    }
    const r = profRestricted();
    if (r.length){
      const names = r.map(profLabel).join(', ');
      bits.push('<strong>Normed to age 69.</strong> ' + escapeHtml(names)
        + (r.length === 1 ? ' is' : ' are') + ' normed only for ages 16:0–69:11, so the correlations used for '
        + (r.length === 1 ? 'it' : 'them') + ' come from that sample.');
    }
    return bits.length ? '<div class="prof-caveat">' + bits.map(b => '<p>' + b + '</p>').join('') + '</div>' : '';
  }

  function profPrecisionHtml(res){
    if (!res || typeof profileAbnormalityStdErr !== 'function') return '';
    const se = profileAbnormalityStdErr(res.lowScores[0], res.trials);
    if (!Number.isFinite(se)) return '';
    return '<p class="prof-precision">Monte Carlo estimates from ' + res.trials.toLocaleString()
      + ' simulated cases — sampling error about ±' + se.toFixed(2)
      + ' percentage points at the largest value here. Seeded, so the same profile always returns the same figures.</p>';
  }

  function renderProfile(){
    const picker = profEl('prof-picker');
    if (picker) picker.innerHTML = profPickerHtml();

    const out = profEl('prof-results');
    if (!out) return;
    const k = profState.selected.length;
    const c = profCriterion();

    if (k < PROF_MIN){
      out.innerHTML = '<div class="prof-empty">Choose at least two measures. '
        + 'A profile is about how findings accumulate across a battery, so one measure has nothing to accumulate over.</div>';
      renderProfileApa();
      return;
    }
    const res = profSimulate();
    const counts = profCounts();
    if (!res){
      out.innerHTML = '<div class="prof-empty">These measures do not form a usable correlation matrix, so no profile can be computed.</div>';
      renderProfileApa();
      return;
    }
    const complete = counts && counts.complete;
    let html = '<div class="prof-selected-line">' + k + ' measures · '
      + escapeHtml(c.label) + (complete ? '' : ' · ' + counts.entered + ' of ' + k + ' scored') + '</div>';
    html += '<div class="prof-cards">'
      + profCardHtml('Abnormally low scores',
          'Scores ' + c.label + '. Across ' + k + ' correlated measures, showing one is far more common than the criterion alone suggests.',
          complete ? counts.low : null, res.lowScores, 'measures')
      + profCardHtml('Abnormal pairwise differences',
          'Differences between any two of the ' + k + ' measures larger than 95% of the population shows, regardless of direction.',
          complete ? counts.pair : null, res.pairwise, 'pairs')
      + profCardHtml('Abnormal deviations from own mean',
          'Scores differing from this patient’s own mean across the set by more than 95% of the population does.',
          complete ? counts.dev : null, res.deviations, 'measures')
      + '</div>';
    if (complete){
      html += '<div class="prof-which-row">'
        + '<div><span class="prof-which-label">Low</span>' + profWhichHtml(counts.lowWhich) + '</div>'
        + '<div><span class="prof-which-label">Differences</span>' + profWhichHtml(counts.pairWhich) + '</div>'
        + '<div><span class="prof-which-label">Deviations</span>' + profWhichHtml(counts.devWhich) + '</div>'
        + '</div>';
    }
    html += profCaveatHtml();
    html += '<details class="prof-details"><summary>Full distribution and precision</summary>'
      + '<div class="prof-details-body">' + profTableHtml(res, counts) + profPrecisionHtml(res) + '</div></details>';
    out.innerHTML = html;
    renderProfileApa();
  }

  /* THE APA TABLE IS THE PATIENT'S, and is emitted only once every chosen
     measure has a score. The population percentages need no patient data
     at all, so an unguarded renderer would file a full-looking table for
     someone who has not been assessed. */
  function renderProfileApa(){
    const out = profEl('prof-apa');
    if (!out) return;
    const counts = profCounts();
    const res = profSimulate();
    if (!counts || !counts.complete || !res){
      out.innerHTML = '<div style="color:var(--faint);font-style:italic;font-family:var(--sans);font-size:13px">'
        + 'Choose at least two measures and score them all to generate the APA table.</div>';
      return;
    }
    const c = profCriterion();
    const rows = [
      ['Scores ' + c.label, counts.low, profPctFor(res.lowScores, counts.low)],
      ['Abnormally large pairwise differences', counts.pair, profPctFor(res.pairwise, counts.pair)],
      ['Abnormally large deviations from own mean', counts.dev, profPctFor(res.deviations, counts.dev)]
    ];
    let body = '';
    for (const [label, n, pct] of rows){
      body += '<tr><td>' + escapeHtml(label) + '</td><td class="num">' + n + '</td>'
        + '<td class="num">' + (n === 0 ? '—' : profFmtPct(pct)) + '</td></tr>';
    }
    const scoreLine = profState.selected
      .map(k => profLabel(k) + ' ' + profNum(profState.scores[k]))
      .join(', ');
    out.innerHTML =
      '<div class="apa-table-num">Table 1</div>'
      + '<div class="apa-table-title">Number of abnormal WAIS-IV findings and their base rates</div>'
      + '<table class="apa-table"><thead><tr>'
      + '<th>Finding</th><th class="num">Number observed</th><th class="num">Base rate</th>'
      + '</tr></thead><tbody>' + body + '</tbody></table>'
      + (typeof apaNoteHtml === 'function'
          ? apaNoteHtml('prof', { criterion:c.label, pct:c.pct, trials:res.trials, scores:scoreLine,
                                  k:profState.selected.length,
                                  metric:profHasScaled() ? 'scaled scores among them' : 'Index scores',
                                  coarse:profHasScaled(),
                                  restricted:profRestricted().map(profLabel).join(', ') })
          : '');
  }

  /* ---------- wiring ---------- */
  function setupProfile(){
    const sel = profEl('prof-criterion');
    if (sel){
      sel.innerHTML = PROF_CRITERIA.map(c =>
        '<option value="' + c.id + '"' + (c.id === profState.criterion ? ' selected' : '') + '>'
        + escapeHtml(c.label) + '</option>').join('');
    }
    const compo = profEl('prof-composition');
    if (compo){
      compo.innerHTML = 'VCI <em>is</em> ' + (PROF_COMPOSED_OF.VCI || []).map(profLabel).map(escapeHtml).join(' + ');
    }

    section.addEventListener('change', e => {
      const t = e.target;
      if (!t) return;
      if (t.dataset && t.dataset.profKey){
        const key = t.dataset.profKey;
        const at = profState.selected.indexOf(key);
        if (t.checked && at === -1) profState.selected.push(key);
        else if (!t.checked && at !== -1){
          profState.selected.splice(at, 1);
          delete profState.pulled[key];
        }
        profSortSelection();
        renderProfile();
      } else if (t.id === 'prof-criterion'){
        profState.criterion = t.value;
        renderProfile();
      }
    });
    /* Scores update in place: the picker is re-rendered only when the
       SELECTION changes, so typing never takes the focused box away. */
    section.addEventListener('input', e => {
      const t = e.target;
      if (t && t.dataset && t.dataset.profScore){
        profState.scores[t.dataset.profScore] = t.value;
        delete profState.pulled[t.dataset.profScore];
        renderResultsOnly();
      }
    });

    const pull = profEl('prof-pull');
    if (pull){
      pull.addEventListener('click', () => {
        const n = profPull();
        renderProfile();
        if (typeof showToast === 'function'){
          showToast(n ? n + ' measure' + (n === 1 ? '' : 's') + ' pulled from Score Tables'
                      : 'No new WAIS-IV scores on Score Tables to pull', !n);
        }
      });
    }
    const clear = profEl('prof-clear');
    if (clear){
      clear.addEventListener('click', () => {
        profState.selected = [];
        profState.scores = {};
        profState.pulled = {};
        renderProfile();
      });
    }
    renderProfile();
  }

  /* Results without touching the picker - see the input handler above. */
  function renderResultsOnly(){
    const out = profEl('prof-results');
    if (!out) return;
    const picker = profEl('prof-picker');
    const keep = picker ? picker.innerHTML : null;
    renderProfile();
    if (picker && keep !== null) picker.innerHTML = keep;
  }

  setupProfile();
  window.renderProfile = renderProfile;
})();
