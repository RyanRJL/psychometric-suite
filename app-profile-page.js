/* =====================================================================
   Profile Analysis - page module (#profile)

   Crawford, Garthwaite & Gault (2007), Neuropsychology, 21, 419-430.

   THE QUESTION, AND IT IS THE WHOLE PAGE. A battery raises one a single
   row cannot answer: this patient has three scores below the 5th
   percentile - how unusual is that? By definition 5% of the population
   falls below the 5th percentile on any ONE measure. Across four
   correlated WAIS-IV Indices 13.8% show at least one; across ten core
   subtests, 24.9%. Reading rows independently overcalls impairment by
   enough to change a conclusion.

   So the answer is a COUNT paired with a BASE RATE - "2 of 4 are
   abnormally low; 4.4% of the population show 2 or more" - and the page
   is built around that pairing. Each of the three questions gets a card
   carrying the count, the sentence, and the population distribution of
   that count as a bar per possible value, with this patient's own bar
   marked and everything at or above it shaded. The shaded part IS the
   base rate, so the two halves are one picture rather than two numbers.

   An earlier version had it the other way round: a tall list of measures
   with a marker each took most of the page and the answer was a strip of
   small figures at the top. The measures are the INPUT. The count and
   its base rate are the OUTPUT, and the output is what this page exists
   to print.

   THREE INSTRUMENTS, ONE PER TAB. WAIS-IV (Table 5.1), WMS-IV Adult
   (Table 4.1) and WMS-IV Older Adult (Table 4.2). The two WMS-IV
   batteries are separate entries rather than one instrument with an age
   switch: they are different normative samples with different measure
   lists and different coefficients, ages 65-69 are normed in both, and
   the clinician already chose between them when they chose what to
   administer. The group key on the Score Tables row carries that choice,
   so the tab is read from the data rather than asked for again. Only
   instruments Score Tables actually holds a score for get a tab, and a
   lone instrument gets none at all.

   WMS-IV also bends the level rule that WAIS-IV set. Its five indices
   are not a partition of its subtests - AMI/VMI/VWMI split them by
   modality and IMI/DMI split the SAME subtests by delay - so AMI and IMI
   share Logical Memory I and can never be profiled together. It
   therefore gets two index levels where WAIS-IV needs one. Nothing about
   the rule changed; an instrument that indexes its subtests twice simply
   has two conflict-free ways to read them.

   NOTHING IS TYPED HERE. Every score is read from Score Tables, live -
   see profScoreTableRowsFor. The page holds no score of its own, so there
   is no second place a WAIS-IV or WMS-IV score can live and no way for
   the two to disagree. The measures that went in are a single strip of chips: a
   chip turns red when that score met the criterion, so the red chips are
   exactly the ones the headline counted. Clicking a chip drops it out.

   ANY SET OF MEASURES, NOT A FIXED BATTERY. The paper is explicit three
   times over: the program "requires the user to specify the number of
   tests in the battery (up to a maximum of 20 tests)... and enter the
   correlation between tests in the form of a lower triangular matrix";
   the conclusion speaks of abnormality "from among the overall set of
   tests administered"; and Crawford's own supplementary programs state
   that "the methods can be applied when only a subset of the Index
   scores have been administered". The requirement is only that R covers
   the measures in hand. What cannot be done is quote a ten-measure
   figure for seven measures, so profSimulate keys its cache on the
   selection itself and nothing is carried over.

   THE ARITHMETIC IS NOT HERE. profileAbnormality() and the matrices live
   in app.js and data.js, verified in check.js sections 42, 43 and 48 -
   the engine against the paper's Appendix, the exact binomial and all 12
   rows of its Table 1; each matrix against its own printed page.
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
    wais4: {
      VCI:  ['SI', 'VC', 'IN', 'CO'],
      PRI:  ['BD', 'MR', 'VP', 'FW', 'PCm'],
      WMI:  ['DS', 'AR', 'LN'],
      PSI:  ['SS', 'CD', 'CA'],
      FSIQ: ['BD', 'SI', 'DS', 'MR', 'VC', 'AR', 'SS', 'VP', 'IN', 'CD'],
      DS:   ['DSF', 'DSB', 'DSS']
    },
    /* WMS-IV, ADULT BATTERY. Its five indices are NOT a partition of its
       subtests - they cross. AMI and VMI split them by modality, IMI and DMI
       split the SAME subtests by delay, so AMI and IMI both hold Logical
       Memory I and can never share a profile. That is not a quirk of this
       page: it is the part-whole rule meeting an instrument that indexes its
       subtests twice, and it is why WMS-IV gets two index levels where
       WAIS-IV needs one.

       Composition is not taken on trust. The manual's own Mean row is a sum
       of scaled scores, so it counts the subtests in each index: adult AMI
       40.0 and IMI 40.1 are four subtests each, VWMI 20.0 is two. And the
       shaded upper triangle prints a corrected coefficient for exactly those
       pairs where a subtest sits inside an index - LM I has one against AMI
       and IMI and against nothing else, SA has one against VWMI alone - which
       is the manual stating the membership cell by cell. check.js drives
       both rather than restating this list. */
    wms4: {
      AMI:  ['LM1', 'LM2', 'VPA1', 'VPA2'],
      VMI:  ['DE1', 'DE2', 'VR1', 'VR2'],
      VWMI: ['SA', 'SSP'],
      IMI:  ['LM1', 'VPA1', 'DE1', 'VR1'],
      DMI:  ['LM2', 'VPA2', 'DE2', 'VR2'],
      /* Process scores sit inside the subtest they are scored from, exactly
         as the three Digit Span process scores do on WAIS-IV. */
      VPA2: ['VPAWR'],
      DE1:  ['DE1C', 'DE1S'],
      DE2:  ['DE2C', 'DE2S']
    },
    /* OLDER ADULT BATTERY. Designs, Spatial Addition and the VWMI are not
       administered, so VMI is the two Visual Reproduction subtests and Symbol
       Span belongs to no index at all - which the manual confirms by printing
       no corrected cell for it. */
    wms4o: {
      AMI:  ['LM1', 'LM2', 'VPA1', 'VPA2'],
      VMI:  ['VR1', 'VR2'],
      IMI:  ['LM1', 'VPA1', 'VR1'],
      DMI:  ['LM2', 'VPA2', 'VR2'],
      VPA2: ['VPAWR']
    }
  };
  /* Block Design No Time Bonus is the SAME administration rescored - the
     matrix puts them at r = .97 - so they are one measure for this
     purpose even though neither contains the other. WMS-IV has no such
     pair: every process score there is a part, not a rescoring. */
  const PROF_ALIAS = { wais4: { BDN: 'BD' }, wms4: {}, wms4o: {} };

  function profComponents(key, seen){
    const inst = profState.instrument;
    const k = (PROF_ALIAS[inst] || {})[key] || key;
    const out = seen || new Set();
    const parts = (PROF_COMPOSED_OF[inst] || {})[k];
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

  /* ---------- the roster on screen: THE PART-WHOLE RULE MEANS YOU ARE
     ALWAYS CHOOSING A LEVEL ----------

     Indices and their own subtests can never share a profile. An earlier
     version listed all four groups at once and greyed out whatever the
     current selection blocked - which, with the four Indices ticked, was
     ten of fourteen rows, each carrying a line of explanation. The page
     ran past the fold before the answer was reached.

     So the page asks which level is being profiled and lists only that.
     Each group below is INTERNALLY CONFLICT-FREE by construction, so a
     conflict cannot arise on screen at all. The rule is still enforced
     underneath (profDisabledBy, applied on every sync); it simply has
     nothing left to block. check.js asserts the conflict-freedom rather
     than trusting the list.

     FULL SCALE IQ IS DELIBERATELY NOT OFFERED. It contains ten of the
     fifteen subtests and all four Indices, so it has no level: every
     profile it could legally join (FSIQ with the five supplementary
     subtests, say) is one no clinician would run. It is named in the
     footer when Score Tables holds it, so its absence reads as a rule
     rather than as missing data. */
  const PROF_INSTRUMENTS = [
    { id:'wais4', label:'WAIS-IV', name:'WAIS-IV',
      /* The group key is what decides the instrument, never the row title:
         a custom test called "Vocabulary" is not a WAIS-IV subtest. */
      groupRe: /WAIS-IV/,
      composites: ['FSIQ', 'VCI', 'PRI', 'WMI', 'PSI'],
      matrix: () => (typeof WAIS4_INTERCORR !== 'undefined') ? WAIS4_INTERCORR : null,
      levels: [
        { id:'w4-indices',  label:'Indices',        keys:['VCI', 'PRI', 'WMI', 'PSI'] },
        { id:'w4-subtests', label:'Subtests',
          keys:['BD', 'SI', 'DS', 'MR', 'VC', 'AR', 'SS', 'VP', 'IN', 'CD', 'LN', 'FW', 'CO', 'CA', 'PCm'] },
        { id:'w4-process',  label:'Process scores', keys:['BDN', 'DSF', 'DSB', 'DSS'] }
      ] },
    /* TWO BATTERIES, TWO ENTRIES, because they are two normative samples with
       two matrices and two measure lists - the separation normDB already
       makes. Ages 65-69 are normed in both, so the age cannot choose: the
       battery is whichever group the clinician picked on Score Tables, which
       is what the group key holds. */
    { id:'wms4', label:'WMS-IV · Adult battery', name:'WMS-IV',
      groupRe: /^WMS-IV .*Ages 16-69/,
      composites: ['AMI', 'VMI', 'VWMI', 'IMI', 'DMI'],
      matrix: () => (typeof WMS4_INTERCORR !== 'undefined') ? WMS4_INTERCORR.adult : null,
      /* THE INDICES NEED TWO LEVELS. AMI/VMI/VWMI split the battery by
         modality and IMI/DMI split the same subtests by delay, so the five
         cannot be profiled together - AMI and IMI share Logical Memory I and
         Verbal Paired Associates I. Each list below is internally
         conflict-free, which is the invariant check.js asserts rather than
         trusting the lists.

         VWMI appears in both because it is orthogonal to each split: Spatial
         Addition and Symbol Span are in no other index, so it conflicts with
         nothing on either level and leaving it out would cost a measure for
         no reason. */
      levels: [
        { id:'wm-mod',  label:'Indices by modality',  keys:['AMI', 'VMI', 'VWMI'] },
        { id:'wm-time', label:'Indices by delay',     keys:['IMI', 'DMI', 'VWMI'] },
        { id:'wm-sub',  label:'Subtests',
          keys:['LM1', 'LM2', 'VPA1', 'VPA2', 'DE1', 'DE2', 'VR1', 'VR2', 'SA', 'SSP'] },
        { id:'wm-proc', label:'Process scores',
          keys:['VPAWR', 'DE1C', 'DE1S', 'DE2C', 'DE2S'] }
      ] },
    { id:'wms4o', label:'WMS-IV · Older Adult battery', name:'WMS-IV',
      groupRe: /^WMS-IV .*Ages 65-90/,
      composites: ['AMI', 'VMI', 'IMI', 'DMI'],
      matrix: () => (typeof WMS4_INTERCORR !== 'undefined') ? WMS4_INTERCORR.older : null,
      levels: [
        { id:'wo-mod',  label:'Indices by modality',  keys:['AMI', 'VMI'] },
        { id:'wo-time', label:'Indices by delay',     keys:['IMI', 'DMI'] },
        { id:'wo-sub',  label:'Subtests',
          keys:['LM1', 'LM2', 'VPA1', 'VPA2', 'VR1', 'VR2', 'SSP'] },
        /* One process score is published for this battery, so this level can
           never reach PROF_MIN. It is listed anyway: profLevelsHtml hides a
           level with nothing scored and disables one holding a single
           measure, which reads as "one is not a profile" rather than as an
           omission. */
        { id:'wo-proc', label:'Process scores', keys:['VPAWR'] }
      ] }
  ];
  function profInstrument(){
    return PROF_INSTRUMENTS.find(x => x.id === profState.instrument) || PROF_INSTRUMENTS[0];
  }
  function profLevels(){ return profInstrument().levels; }

  /* Composites are M 100 / SD 15; everything else is a scaled score,
     M 10 / SD 3. Per measure, not per page: a single conversion would
     score a scaled 8 as z = -6.1 and call every subtest catastrophic. */
  function profIsComposite(key){
    return profInstrument().composites.indexOf(key) !== -1;
  }
  function profMetric(key){
    return profIsComposite(key)
      ? { mean:100, sd:15, label:'Index', min:40, max:160 }
      : { mean:10,  sd:3,  label:'scaled', min:1,  max:19 };
  }

  /* THE FIVE ROWS OF THE PAPER'S OWN TABLES 2 AND 3, in its own order:
     "a range of (increasingly stringent) definitions of abnormality were
     applied, ranging from a score that was more than one standard deviation
     below the mean (i.e., below the 15.9th percentile) to a score that fell
     below the 1st percentile." The default is 5%, which is the two statements
     the paper makes in both table notes: "We define abnormality as a score
     falling below the 5th percentile, and for this reason we have presented
     these results in bold."

     `diffZ` IS THE SAME CRITERION, TWO-TAILED, and it follows the selector
     rather than sitting at 1.960 for every row - the paper's program does
     exactly this: "For simplicity, the criterion selected by the user to
     define an abnormally low score is also used to define abnormally large
     pairwise differences and abnormally large deviations from individuals'
     mean scores. For example, if an abnormally low score is defined as a
     score falling below the 5th percentile, then an abnormally large pairwise
     difference (or deviation) is defined as a difference (or deviation) that
     is exceeded by less than 5% of the normal population, regardless of sign."

     So `diffZ` is the deviate cutting the SAME tail area in two halves,
     Phi-inverse(1 - p/2) where p = Phi(z) is the low-score tail. It is a
     stored constant rather than a derivation because each value is a standard
     normal deviate anyone can check against a table, and because the pairing
     is the thing under test - deriving it here and in check.js would be one
     formula agreeing with itself. The arithmetic, to 4 dp:

       z -1.0    p 15.8655%  ->  1.4101   (tail 84.1% of the population)
       z -1.282  p  9.9921%  ->  1.6452                    90%
       z -1.645  p  4.9985%  ->  1.9601                    95%   <- the paper's
       z -2.054  p  1.9988%  ->  2.3266                    98%
       z -2.326  p  1.0009%  ->  2.5755                    99%

     `diffPct` is that tail as the percentage the difference must EXCEED, which
     is what the cards and the exported note say in words. */
  const PROF_CRITERIA = [
    { id:'sd1', label:'below 1 SD (15.9th percentile)', z:-1.0,   pct:'15.9th', diffZ:1.410, diffPct:'84.1%' },
    { id:'p10', label:'below the 10th percentile',      z:-1.282, pct:'10th',   diffZ:1.645, diffPct:'90%' },
    { id:'p5',  label:'below the 5th percentile',       z:-1.645, pct:'5th',    diffZ:1.960, diffPct:'95%' },
    { id:'p2',  label:'below the 2nd percentile',       z:-2.054, pct:'2nd',    diffZ:2.326, diffPct:'98%' },
    { id:'p1',  label:'below the 1st percentile',       z:-2.326, pct:'1st',    diffZ:2.576, diffPct:'99%' }
  ];
  const PROF_TRIALS = 200000;
  const PROF_MIN = 2;          // one measure is not a profile

  /* `scores` is a CACHE of what Score Tables held at the last sync, not a
     second store: profPull rewrites it wholesale and nothing else ever
     assigns to it. `excluded` is per level and clears when the level
     changes, the levels having no measures in common. */
  const profState = { criterion:'p5', instrument:'wais4', level:'w4-indices', selected:[], scores:{}, excluded:{} };

  /* Keyed on the SELECTION and the criterion. Scores do not enter the
     simulation, so a score changing on Score Tables must never re-run it
     - and the run costs ~200ms at 15 measures, comparing 105 pairs per
     simulated case. */
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
    try { return profInstrument().matrix() || null; } catch (e){ return null; }
  }
  function profLabel(key){
    const M = profMatrixObj();
    return (M && M.labels && M.labels[key]) || key;
  }
  /* An Index is universally known by its acronym and its full name is
     three words long, so the chip carries whichever reads shorter. */
  function profShortLabel(key){
    return profIsComposite(key) ? key : profLabel(key);
  }
  function profGroup(id){
    const L = profLevels();
    return L.find(g => g.id === id) || L[0];
  }
  function profRestricted(){
    const M = profMatrixObj();
    const list = (M && M.restrictedTo16_69) || [];
    return profState.selected.filter(k => list.indexOf(k) !== -1);
  }
  /* Any selected measure that is not a composite is a scaled score, and
     the paper's own limitation applies to those - see profCaveatHtml. */
  function profHasScaled(){
    return profState.selected.some(k => !profIsComposite(k));
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
    const key = profState.instrument + '|' + profState.selected.join(',') + '|' + profCriterion().id;
    if (profCache[key]) return profCache[key];
    const R = profMatrix();
    if (!R || typeof profileAbnormality !== 'function') return null;
    const c = profCriterion();
    const res = profileAbnormality(R, { lowZ:c.z, diffZ:c.diffZ, trials:PROF_TRIALS });
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
    /* Selection is derived from scored rows, so this cannot normally
       fire - it is the guard for a score vanishing from Score Tables
       between a sync and a render. */
    if (entered < keys.length) return { complete:false, entered, total:keys.length };

    let low = 0; const lowKeys = [];
    z.forEach((v, i) => { if (v < c.z){ low++; lowKeys.push(keys[i]); } });

    let pair = 0;
    for (let i = 0; i < keys.length; i++){
      for (let j = i + 1; j < keys.length; j++){
        const t = c.diffZ * Math.sqrt(2 - 2 * R[i][j]);
        if (Math.abs(z[i] - z[j]) > t) pair++;
      }
    }

    const means = (typeof profileMatrixMeans === 'function') ? profileMatrixMeans(R) : null;
    let dev = 0;
    if (means){
      const mean = z.reduce((a, b) => a + b, 0) / keys.length;
      z.forEach((v, i) => {
        const t = c.diffZ * Math.sqrt(1 + means.grandMean - 2 * means.rowMean[i]);
        if (Math.abs(mean - v) > t) dev++;
      });
    }
    return { complete:true, entered, total:keys.length, low, lowKeys, pair, dev,
             k:keys.length, pairs:keys.length * (keys.length - 1) / 2 };
  }

  /* ---------- the scores, which come from Score Tables and nowhere else ----------

     Matched on the measure NAME, which is the same string in both files:
     24 of the 28 WAIS-IV names in normDB are Table 5.1 labels verbatim.
     The four that are not are the Longest Span base-rate measures, which
     have no correlations and correctly cannot be profiled. WMS-IV is
     matched the same way, and its matrix carries normDB's spelling rather
     than the manual's column heads ("Designs I - Content" for the printed
     "DE I Content") for exactly this reason: the name is the join. */
  function profScoreTableRowsFor(inst){
    if (typeof batteryRows === 'undefined' || !Array.isArray(batteryRows)) return {};
    const M = inst.matrix ? inst.matrix() : null;
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
      /* THE GROUP KEY DECIDES, NOT THE ROW TITLE. A row named "Logical
         Memory I" is a WMS-IV Adult measure only if it was taken from the
         Ages 16-69 group; the Older Adult battery carries the same measure
         names against a different normative sample and a different matrix,
         and a custom test may carry any name at all. A row with no group
         key cannot be placed and is refused rather than guessed at. */
      const group = (typeof batteryGroupKeyOf === 'function') ? batteryGroupKeyOf(r) : '';
      if (!group || !inst.groupRe.test(group)) return;
      out[key] = v;
    });
    return out;
  }
  function profScoreTableRows(){
    return profScoreTableRowsFor(profInstrument());
  }
  /* Which instruments Score Tables actually holds something for, in the
     registry's own order. The tabs are built from this, so an instrument
     with nothing entered never appears. */
  function profInstrumentsFound(){
    return PROF_INSTRUMENTS.map(inst => ({ inst, found:profScoreTableRowsFor(inst) }))
      .filter(x => Object.keys(x.found).length > 0);
  }

  /* Rebuild the selection from Score Tables. Returns the number of
     measures now in the profile. This is the ONLY writer of
     profState.selected and profState.scores.

     The conflict guard is kept even though a level is conflict-free by
     construction: it is the thing that makes that construction safe to
     change, and it costs nothing. Nothing here writes back - this page
     is a reader of Score Tables, never an editor of it. */
  function profPull(){
    /* THE INSTRUMENT FIRST, and only when the current one has nothing:
       a clinician who has chosen a tab keeps it while it holds a score,
       exactly as they keep a level. Landing on the page with only WMS-IV
       entered must not show an empty WAIS-IV. */
    if (!Object.keys(profScoreTableRows()).length){
      const alt = profInstrumentsFound()[0];
      if (alt && alt.inst.id !== profState.instrument) profSetInstrument(alt.inst.id);
    }
    const found = profScoreTableRows();
    /* Fall to the first level that can actually form a profile, so
       landing on the page with only subtests entered does not show an
       empty Indices level. Only when the current level has nothing. */
    if (profAvailable(profState.level, found).length === 0){
      const L = profLevels();
      const alt = L.find(g => profAvailable(g.id, found).length >= PROF_MIN)
               || L.find(g => profAvailable(g.id, found).length > 0);
      if (alt && alt.id !== profState.level){ profState.level = alt.id; profState.excluded = {}; }
    }
    profState.selected = [];
    profState.scores = {};
    profAvailable(profState.level, found).forEach(k => {
      if (profState.excluded[k]) return;
      if (profDisabledBy(k)) return;
      profState.selected.push(k);
      profState.scores[k] = String(found[k]);
    });
    return profState.selected.length;
  }
  /* Switching instrument resets the level too: a level id belongs to one
     instrument, and every exclusion was a judgement about a measure in the
     one being left. */
  function profSetInstrument(id){
    if (id === profState.instrument) return;
    profState.instrument = id;
    profState.level = profLevels()[0].id;
    profState.excluded = {};
  }

  /* The measures at a level that Score Tables actually holds a score
     for, in the matrix's own order. */
  function profAvailable(levelId, found){
    const src = found || profScoreTableRows();
    return profGroup(levelId).keys.filter(k => src[k] !== undefined);
  }

  /* ---------- rendering ---------- */
  function profFmtPct(v){
    if (v === null || v === undefined) return '—';
    if (v > 0 && v < 0.01) return '< 0.01%';
    return (v >= 10 ? v.toFixed(1) : v.toFixed(2)) + '%';
  }
  /* series[j-1] is the percentage showing j OR MORE. j = 0 is refused:
     "100% show 0 or more" is true, useless, and reads as a finding. */
  function profPctFor(series, j){
    if (!series || j < 1) return null;
    const v = series[j - 1];
    return Number.isFinite(v) ? v : null;
  }

  /* The point distribution, differenced out of the cumulative one the
     engine returns, so the chart needs no second pass over the trials:
     P(0) = 100 - P(>=1); P(j) = P(>=j) - P(>=j+1). */
  function profPointDist(series, n){
    const out = [];
    out[0] = 100 - (series[0] === undefined ? 0 : series[0]);
    for (let j = 1; j <= n; j++){
      const hi = series[j] === undefined ? 0 : series[j];
      out[j] = (series[j - 1] === undefined ? 0 : series[j - 1]) - hi;
    }
    return out;
  }
  /* Fifteen measures give 105 possible pairs and almost all of that
     range carries no mass, so the axis stops where the mass does - but
     never before this patient's own count, which must always be on
     screen for the shading to mean anything. */
  function profDistTrim(dist, here){
    let last = 0;
    for (let j = 0; j < dist.length; j++) if (dist[j] >= 0.05) last = j;
    return Math.max(last, here, 2);
  }

  function profDistHtml(series, n, here, unit){
    if (!series) return '';
    const dist = profPointDist(series, n);
    const top = profDistTrim(dist, here);
    let max = 0;
    for (let j = 0; j <= top; j++) if (dist[j] > max) max = dist[j];
    const every = top <= 12 ? 1 : Math.ceil((top + 1) / 8);
    let bars = '', axis = '';
    for (let j = 0; j <= top; j++){
      const h = max > 0 ? Math.max(2, Math.round(dist[j] / max * 100)) : 2;
      /* The patient's own bar solid, everything above it shaded: the
         shaded region is the base rate the sentence quotes. Nothing is
         shaded at a count of zero, which is quoted no base rate. */
      const cls = j === here ? 'prof-bar prof-here'
                : (here > 0 && j > here ? 'prof-bar prof-tail' : 'prof-bar');
      bars += '<div class="' + cls + '" style="height:' + h + '%" title="'
        + escapeHtml('exactly ' + j + ' — ' + profFmtPct(dist[j]) + ' of the population') + '"></div>';
      /* A CARET, BECAUSE HEIGHT ALONE CANNOT FIND THIS BAR. Most people show
         none, so P(0) sets the scale and every bar that matters is a hairline
         beside it - at four Indices the patient's own bar was 2px. Rescaling
         to make it visible would misstate the distribution, so the marker goes
         under the axis instead, where it costs the shape nothing. */
      axis += '<span' + (j === here ? ' class="prof-here"' : '') + '>'
        + (j === here || j % every === 0 ? j : '') + '</span>';
    }
    return '<div class="prof-dist"><div class="prof-dist-title">Population — how many '
      + escapeHtml(unit) + ' are abnormal</div>'
      + '<div class="prof-bars">' + bars + '</div>'
      + '<div class="prof-xaxis">' + axis + '</div></div>';
  }

  /* THE CLAIM ITSELF, AT FULL WIDTH. The per-count chart above shows the
     shape; this shows the sentence. One bar for the whole population, split at
     this patient's count, with everything at or beyond it shaded - so a base
     rate of 4.37% is 4.37% of the bar, read directly rather than inferred from
     a stack of columns whose scale is set by P(0). */
  function profTailBarHtml(pct){
    if (!Number.isFinite(pct)) return '';
    const w = Math.max(0.6, Math.min(100, pct));   // a hairline still has to be visible
    return '<div class="prof-tailbar" title="' + escapeHtml(profFmtPct(pct)
      + ' of the population') + '"><span class="prof-tailbar-fill" style="width:' + w + '%"></span></div>';
  }

  /* THE CARD IS THE PAGE: a count, the sentence that reports it, and the
     distribution the count sits in. */
  function profCardHtml(title, tip, count, denom, series, n, unit, unitOne){
    if (count === null){
      return '<div class="prof-card prof-quiet" data-tooltip="' + escapeHtml(tip) + '">'
        + '<div class="prof-card-label">' + escapeHtml(title) + '</div>'
        + '<div class="prof-card-count"><span class="prof-card-n">—</span></div>'
        + '<div class="prof-card-lede">Not enough scored measures for a profile.</div></div>';
    }
    const pct = profPctFor(series, count);
    const se = (pct !== null && typeof profileAbnormalityStdErr === 'function')
      ? profileAbnormalityStdErr(pct, PROF_TRIALS) : null;
    const lede = count === 0
      ? 'None of these ' + unit + ' meet the criterion.'
      : count + ' of ' + denom + ' ' + (count === 1 ? unitOne : unit)
        + (count === 1 ? ' is' : ' are') + ' abnormal.';
    const read = count === 0
      ? '<div class="prof-read prof-read-none">No base rate is quoted for a count of zero.</div>'
      : profTailBarHtml(pct)
        + '<div class="prof-read"><span><b>' + profFmtPct(pct)
        + '</b> of the population show ' + count + ' or more'
        + (Number.isFinite(se) ? ' <span class="prof-se">± ' + se.toFixed(2) + '</span>' : '')
        + '</span></div>';
    return '<div class="prof-card' + (count === 0 ? ' prof-quiet' : '') + '" data-tooltip="'
      + escapeHtml(tip) + '">'
      + '<div class="prof-card-label">' + escapeHtml(title) + '</div>'
      + '<div class="prof-card-count"><span class="prof-card-n">' + count + '</span>'
      + '<span class="prof-card-of">of ' + denom + '</span></div>'
      + '<div class="prof-card-lede">' + escapeHtml(lede) + '</div>'
      + profDistHtml(series, n, count, unit) + read + '</div>';
  }

  /* A chip carries the measure and its score, and turns red when that
     score met the criterion - so the red chips ARE the ones the headline
     counted, and the count can be checked by eye. */
  function profChipsHtml(counts){
    const found = profScoreTableRows();
    const avail = profAvailable(profState.level, found);
    const lowSet = new Set((counts && counts.complete && counts.lowKeys) || []);
    if (!avail.length) return '';
    return avail.map(k => {
      const on = !profState.excluded[k];
      const low = on && lowSet.has(k);
      const cls = 'prof-chip' + (on ? (low ? ' prof-low' : '') : ' prof-off');
      return '<button type="button" class="' + cls + '" data-prof-chip="' + k + '" title="'
        + escapeHtml(on ? 'Click to leave ' + profLabel(k) + ' out of the profile'
                        : 'Click to put ' + profLabel(k) + ' back in') + '">'
        + escapeHtml(profShortLabel(k)) + ' <b>' + escapeHtml(String(found[k])) + '</b></button>';
    }).join('');
  }

  /* ONE TAB PER INSTRUMENT SCORED, and none at all when only one is - a
     single tab is a label pretending to be a control. The battery is named
     on the tab (Adult / Older Adult) because the two are different normative
     samples that share every measure name, so the tab is the only thing on
     screen that says which was used. */
  function profInstrumentsHtml(){
    const list = profInstrumentsFound();
    if (list.length < 2) return '';
    return '<div class="prof-insts">' + list.map(x => {
      const on = x.inst.id === profState.instrument;
      return '<button type="button" class="prof-inst' + (on ? ' prof-on' : '') + '"'
        + ' data-prof-inst="' + x.inst.id + '">' + escapeHtml(x.inst.label)
        + ' (' + Object.keys(x.found).length + ')</button>';
    }).join('') + '</div>';
  }

  function profLevelsHtml(){
    const found = profScoreTableRows();
    return profLevels().map(g => {
      const n = profAvailable(g.id, found).length;
      if (!n) return '';
      const on = g.id === profState.level;
      return '<button type="button" class="prof-lvl' + (on ? ' prof-on' : '') + '"'
        + (n < PROF_MIN ? ' disabled title="Only one is scored on Score Tables, and one measure is not a profile"' : '')
        + ' data-prof-level="' + g.id + '">' + escapeHtml(g.label) + ' (' + n + ')</button>';
    }).join('');
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
      bits.push('<strong>Scaled scores are coarse.</strong> The method assumes continuous scores; '
        + 'Crawford, Garthwaite &amp; Gault note that a limited range of scaled scores costs accuracy, '
        + 'in contrast to Index scores, a scaled point being a third of a standard deviation.');
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

  /* The footer says what went in and, when Score Tables holds it, why
     Full Scale IQ is not among the choices - an absence with no stated
     reason reads as missing data. */
  function profFootHtml(found){
    const avail = profAvailable(profState.level, found);
    const bits = [profState.selected.length + ' of ' + avail.length + ' included · read from Score Tables, '
      + 'nothing is typed on this page'];
    if (profLevels().length > 1){
      bits.push('an index and its own subtests cannot share a profile, so each level is profiled separately');
    }
    /* An absence with no stated reason reads as missing data, so each
       measure the instrument deliberately withholds says why, and only
       while Score Tables actually holds it. */
    if (profState.instrument === 'wais4' && found.FSIQ !== undefined){
      bits.push('Full Scale IQ contains every other measure, so it cannot join a profile');
    }
    return '<div class="prof-foot"><span>' + escapeHtml(bits.join(' · ')) + '</span>'
      + '<a class="prof-foot-link" data-prof-goto="battery">Edit on Score Tables →</a></div>';
  }

  function renderProfile(){
    const out = profEl('prof-results');
    if (!out) return;
    profPull();
    const found = profScoreTableRows();
    const c = profCriterion();

    /* Nothing to profile: the honest state is a pointer at where scores
       live, not an empty set of cards. */
    if (!Object.keys(found).length){
      const names = [...new Set(PROF_INSTRUMENTS.map(x => x.name))].join(' or ');
      out.innerHTML = '<div class="prof-empty">'
        + '<div class="prof-empty-h">No ' + escapeHtml(names) + ' scores yet</div>'
        + '<p class="prof-empty-p">A profile is built from scores already entered on Score Tables. '
        + 'Add a ' + escapeHtml(names) + ' measure there and it appears here.</p>'
        + '<button class="btn" type="button" data-prof-goto="battery">Go to Score Tables</button></div>';
      renderProfileApa();
      return;
    }

    const counts = profCounts();
    const res = profSimulate();
    const k = profState.selected.length;
    const complete = !!(counts && counts.complete && res);

    let html = '<div class="prof-strip">'
      + profInstrumentsHtml()
      + '<div class="prof-strip-head">'
        + '<span class="prof-lvls">' + profLevelsHtml() + '</span>'
        + '<span class="prof-crit">Abnormal means <select id="prof-criterion" '
        + 'aria-label="Criterion for an abnormally low score"></select></span>'
      + '</div>'
      + '<div class="prof-chips">' + profChipsHtml(counts) + '</div>'
      + profFootHtml(found) + '</div>';

    if (!complete){
      html += '<div class="prof-empty prof-empty-inline">'
        + (k < PROF_MIN
            ? 'Choose at least two measures. A profile is about how findings accumulate across a battery, '
              + 'so one measure has nothing to accumulate over.'
            : 'These measures do not form a usable correlation matrix, so no profile can be computed.')
        + '</div>';
    } else {
      html += '<div class="prof-cards">'
        + profCardHtml('Abnormally low scores',
            'Scores ' + c.label + '. Across ' + k + ' correlated measures, showing one is far more common than the criterion alone suggests.',
            counts.low, counts.k, res.lowScores, counts.k, 'scores', 'score')
        + profCardHtml('Abnormal pairwise differences',
            'Differences between any two of the ' + k + ' measures larger than ' + c.diffPct + ' of the population shows, regardless of direction.',
            counts.pair, counts.pairs, res.pairwise, counts.pairs, 'pairs', 'pair')
        + profCardHtml('Abnormal deviations from own mean',
            'Scores differing from this patient’s own mean across the set by more than ' + c.diffPct + ' of the population does.',
            counts.dev, counts.k, res.deviations, counts.k, 'scores', 'score')
        + '</div>';
      html += profCaveatHtml();
      html += '<details class="prof-details"><summary>Full distribution and precision</summary>'
        + '<div class="prof-details-body">' + profTableHtml(res, counts) + profPrecisionHtml(res) + '</div></details>';
    }
    out.innerHTML = html;

    const sel = profEl('prof-criterion');
    if (sel){
      sel.innerHTML = PROF_CRITERIA.map(x =>
        '<option value="' + x.id + '"' + (x.id === profState.criterion ? ' selected' : '') + '>'
        + escapeHtml(x.label) + '</option>').join('');
    }
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
        + 'Score at least two compatible WAIS-IV measures on Score Tables to generate the APA table.</div>';
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
    const inst = profInstrument();
    const M = profMatrixObj();
    const scoreLine = profState.selected
      .map(k => profLabel(k) + ' ' + profNum(profState.scores[k]))
      .join(', ');
    out.innerHTML =
      '<div class="apa-table-num">Table 1</div>'
      + '<div class="apa-table-title">Number of abnormal ' + escapeHtml(inst.name)
        + ' findings and their base rates</div>'
      + '<table class="apa-table"><thead><tr>'
      + '<th>Finding</th><th class="num">Number observed</th><th class="num">Base rate</th>'
      + '</tr></thead><tbody>' + body + '</tbody></table>'
      + (typeof apaNoteHtml === 'function'
          ? apaNoteHtml('prof', { criterion:c.label, pct:c.pct, diffPct:c.diffPct, trials:res.trials, scores:scoreLine,
                                  k:profState.selected.length,
                                  metric:profHasScaled() ? 'scaled scores among them' : 'Index scores',
                                  scoreKind:profHasScaled() ? 'Scores' : 'Index scores',
                                  coarse:profHasScaled(),
                                  restricted:profRestricted().map(profLabel).join(', '),
                                  /* The SHORT form - see `citation` in data.js.
                                     Falls back to the full title rather than to
                                     nothing, so a matrix without one still names
                                     its table. */
                                  matrixSource:(M && (M.citation || M.source)) || '' })
          : '');
  }

  /* Score Tables is the only source of scores here, so a score typed
     there is the only thing that can change what this page shows. Called
     from renderBattery. Coalesced across a burst of keystrokes: the
     simulation is cached on the SELECTION, so a score-only change is a
     cache hit, but the keystroke that first scores a new measure changes
     the selection and costs a run. */
  let profSyncTimer = null;
  function profileScoresChanged(){
    if (profSyncTimer) clearTimeout(profSyncTimer);
    profSyncTimer = setTimeout(() => { profSyncTimer = null; renderProfile(); }, 180);
  }

  /* ---------- wiring ---------- */
  function setupProfile(){
    section.addEventListener('click', e => {
      const t = e.target.closest
        ? e.target.closest('[data-prof-inst],[data-prof-level],[data-prof-chip],[data-prof-goto]') : null;
      if (!t) return;
      if (t.dataset.profInst){
        if (t.dataset.profInst === profState.instrument) return;
        profSetInstrument(t.dataset.profInst);
        renderProfile();
      } else if (t.dataset.profLevel){
        if (t.dataset.profLevel === profState.level) return;
        profState.level = t.dataset.profLevel;
        profState.excluded = {};   // the levels share no measures
        renderProfile();
      } else if (t.dataset.profChip){
        const k = t.dataset.profChip;
        if (profState.excluded[k]) delete profState.excluded[k];
        else profState.excluded[k] = true;
        renderProfile();
      } else if (t.dataset.profGoto && typeof navigateTo === 'function'){
        navigateTo(t.dataset.profGoto);
      }
    });
    section.addEventListener('change', e => {
      if (e.target && e.target.id === 'prof-criterion'){
        profState.criterion = e.target.value;
        renderProfile();
      }
    });
    renderProfile();
  }

  setupProfile();
  window.renderProfile = renderProfile;
  window.profileScoresChanged = profileScoresChanged;
})();
