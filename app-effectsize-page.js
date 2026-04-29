/* =====================================================================
   Effect Size Tools — page module (#effectsize)
   ===================================================================== */
(function(){
  'use strict';

  // ----- Standard normal helpers -----
  function normPdf(x){ return Math.exp(-0.5*x*x) / Math.sqrt(2*Math.PI); }
  // Abramowitz & Stegun 26.2.17
  function normCdf(x){
    const sign = x < 0 ? -1 : 1;
    const ax = Math.abs(x) / Math.sqrt(2);
    const t = 1 / (1 + 0.3275911 * ax);
    const y = 1 - (((((1.061405429*t - 1.453152027)*t) + 1.421413741)*t - 0.284496736)*t + 0.254829592)*t * Math.exp(-ax*ax);
    return 0.5 * (1 + sign * y);
  }

  const $ = id => document.getElementById(id);
  const els = {};
  ['es-stat-type','es-stat-value','es-stat-aux','es-stat-aux-field','es-stat-aux-label',
   'es-g1-n','es-g1-mean','es-g1-disp','es-g1-disp-val','es-g1-disp-label','es-g1-ci',
   'es-g2-n','es-g2-mean','es-g2-disp','es-g2-disp-val','es-g2-disp-label','es-g2-ci',
   'es-pooled-n','es-pooled-mean','es-pooled-sd','es-pooled-se','es-pooled-md','es-clear','es-reset-all','es-apply-groups',
   'es-out-d','es-out-d-desc','es-out-g','es-out-g-desc','es-out-r','es-out-r2',
   'es-out-f','es-out-z','es-out-or','es-out-ovl','es-out-u3','es-out-cles','es-out-nnt',
   'es-out-r-desc','es-out-r2-desc','es-out-f-desc','es-out-similar','es-method-warning',
   'es-curve','es-curve-caption',
   'es-target','es-tgt-g1-above','es-tgt-g2-above','es-tgt-share-g1','es-tgt-share-g2','es-tgt-ratio','es-target-explain-toggle','es-target-explain',
   'es-cl-summary','es-cl-u3line','es-d-slider','es-d-slider-val'
  ].forEach(id => els[id] = $(id));

  if (!els['es-stat-type']) return;
  const esState = { statTouched: false, statAutoFilled: false, source: 'stat' };
  const SIMILAR_EFFECTS = [
    { label: 'Heavy smokers vs never (lung cancer)', d: 2.60 },
    { label: 'Smokers vs never (lung cancer)', d: 1.75 },
    { label: 'UK male vs female height', d: 1.37 },
    { label: 'Former vs never smoker (lung cancer)', d: 1.10 },
    { label: 'Clozapine vs placebo (schizophrenia)', d: 0.89 },
    { label: 'CBT vs usual care (social anxiety)', d: 0.44 },
    { label: 'CBT vs usual care (depression)', d: 0.43 },
    { label: 'Paracetamol vs placebo (headache)', d: 0.26 }
  ];

  // Aux field per statistic (some conversions need N or pooled SD)
  const AUX_FOR = {
    g:     { label: 'Total N (n1 + n2)', hint: 'N' },
    t:     { label: 'Total N (df + 2)',  hint: 'N' },
    md:    { label: 'Pooled SD',         hint: 'SD' },
    zstat: { label: 'Total N',           hint: 'N' }
  };
  function refreshAuxField(){
    const t = els['es-stat-type'].value;
    const cfg = AUX_FOR[t];
    if (cfg){
      els['es-stat-aux-field'].style.display = '';
      els['es-stat-aux-label'].textContent = cfg.label;
      els['es-stat-aux'].placeholder = cfg.hint;
    } else {
      els['es-stat-aux-field'].style.display = 'none';
    }
  }

  // Convert any statistic -> Cohen's d
  function statToD(type, val, aux){
    if (val === '' || val === null || isNaN(val)) return null;
    const v = Number(val);
    const a = (aux === '' || aux === null || isNaN(aux)) ? null : Number(aux);
    switch(type){
      case 'd':  return v;
      case 'g':
        if (a == null || a < 4) return v;
        return v / Math.sqrt(1 - 3/(4*a - 9));
      case 'r':  if (Math.abs(v) >= 1) return null; return 2*v / Math.sqrt(1 - v*v);
      case 'r2': if (v < 0 || v >= 1) return null; return 2*Math.sqrt(v) / Math.sqrt(1 - v);
      case 'f':  return 2*v;
      case 'z': {
        // Fisher's z -> r -> d
        const r = Math.tanh(v);
        if (Math.abs(r) >= 1) return null;
        return 2*r / Math.sqrt(1 - r*r);
      }
      case 'or':
        if (v <= 0) return null;
        return Math.log(v) * Math.sqrt(3) / Math.PI;
      case 't':
        // Requires equal-group assumption when only total N is supplied.
        // With total N=a (df=a-2), d = 2t/sqrt(a)
        if (a == null || a <= 0) return null;
        return (2*v) / Math.sqrt(a);
      case 'zstat': if (a == null || a <= 0) return null; return v / Math.sqrt(a);
      case 'md': if (a == null || a <= 0) return null; return v / a;
      default: return null;
    }
  }

  function readGroupData(){
    const num = id => { const v = els[id].value; return v === '' ? null : Number(v); };
    const n1 = num('es-g1-n'),  n2 = num('es-g2-n');
    const m1 = num('es-g1-mean'), m2 = num('es-g2-mean');
    const dv1 = num('es-g1-disp-val'), dv2 = num('es-g2-disp-val');
    const dt1 = els['es-g1-disp'].value, dt2 = els['es-g2-disp'].value;

    function sdFrom(disp, val, n, mean){
      if (val == null) return null;
      if (disp === 'sd')  return val;
      if (disp === 'se')  return n != null && n > 0 ? val * Math.sqrt(n) : null;
      if (disp === 'ciu') return mean != null ? (val - mean) / 1.96 : null;
      return null;
    }
    const sd1 = sdFrom(dt1, dv1, n1, m1);
    const sd2 = sdFrom(dt2, dv2, n2, m2);

    const out = { n1, n2, m1, m2, sd1, sd2, dt1, dt2 };
    if (n1 && n2 && sd1 != null && sd2 != null && n1 > 1 && n2 > 1){
      out.pooledN  = n1 + n2;
      if (m1 != null && m2 != null){
        out.pooledMean = ((n1 * m1) + (n2 * m2)) / (n1 + n2);
      }
      out.pooledSD = Math.sqrt(((n1-1)*sd1*sd1 + (n2-1)*sd2*sd2) / (n1 + n2 - 2));
      out.pooledSE = out.pooledSD * Math.sqrt(1/n1 + 1/n2);
      if (m1 != null && m2 != null){
        out.md = m1 - m2;
        out.dFromGroups = (m1 - m2) / out.pooledSD;
      }
    }
    return out;
  }

  function fmtPct(v, digits){
    if (v == null || isNaN(v) || !isFinite(v)) return '—';
    return (v*100).toFixed(digits != null ? digits : 2) + '%';
  }
  function ciString(mean, sd, n){
    if (mean == null || sd == null || n == null || n < 2) return '—';
    const se = sd / Math.sqrt(n);
    return (mean - 1.96*se).toFixed(2) + ', ' + (mean + 1.96*se).toFixed(2);
  }
  function descD(d){
    const a = Math.abs(d);
    if (a <= 0.01) return { label: 'Very Small', mag: 0 };
    if (a <= 0.20) return { label: 'Small',      mag: 1 };
    if (a <= 0.50) return { label: 'Medium',     mag: 2 };
    if (a <= 0.80) return { label: 'Large',      mag: 3 };
    if (a <= 1.20) return { label: 'Very Large', mag: 4 };
    if (a <= 2.00) return { label: 'Huge',       mag: 5 };
    return                 { label: 'Huge',      mag: 6 };
  }
  function descR(r){
    const a = Math.abs(r);
    if (a < 0.10) return { label: 'Very Small', mag: 0 };
    if (a < 0.30) return { label: 'Small',      mag: 1 };
    if (a < 0.50) return { label: 'Medium',     mag: 2 };
    return               { label: 'Large',      mag: 3 };
  }
  function descR2(r2){
    if (r2 == null || r2 < 0) return null;
    if (r2 < 0.01) return { label: 'Very Small', mag: 0 };
    if (r2 < 0.09) return { label: 'Small',      mag: 1 };
    if (r2 < 0.25) return { label: 'Medium',     mag: 2 };
    return                { label: 'Large',      mag: 3 };
  }
  function descF(f){
    const a = Math.abs(f);
    if (a < 0.10) return { label: 'Very Small', mag: 0 };
    if (a < 0.25) return { label: 'Small',      mag: 1 };
    if (a < 0.40) return { label: 'Medium',     mag: 2 };
    return               { label: 'Large',      mag: 3 };
  }
  function applyDesc(el, descObj){
    if (!descObj){ el.textContent = '—'; el.removeAttribute('data-mag'); return; }
    el.textContent = descObj.label;
    el.setAttribute('data-mag', String(descObj.mag));
  }

  function setMethodMessage(type, html){
    const el = els['es-method-warning'];
    if (!el) return;
    el.className = 'es-method-warning';
    if (!html){
      el.innerHTML = '';
      return;
    }
    el.classList.add('is-visible');
    if (type === 'error') el.classList.add('is-error');
    el.innerHTML = html;
  }

  function similarEffect(d){
    if (d == null || isNaN(d) || !isFinite(d)) return '—';
    const ad = Math.abs(d);
    if (ad > 3) return '—';
    const closest = SIMILAR_EFFECTS.reduce((best, item) => {
      const gap = Math.abs(item.d - ad);
      return !best || gap < best.gap ? { ...item, gap } : best;
    }, null);
    return closest ? closest.label : '—';
  }

  function compute(){
    refreshAuxField();
    const grp = readGroupData();

    els['es-g1-ci'].textContent = ciString(grp.m1, grp.sd1, grp.n1);
    els['es-g2-ci'].textContent = ciString(grp.m2, grp.sd2, grp.n2);
    els['es-pooled-n'].textContent  = grp.pooledN  != null ? grp.pooledN : '—';
    els['es-pooled-mean'].textContent = grp.pooledMean != null ? grp.pooledMean.toFixed(3) : '—';
    els['es-pooled-sd'].textContent = grp.pooledSD != null ? grp.pooledSD.toFixed(3) : '—';
    els['es-pooled-se'].textContent = grp.pooledSE != null ? grp.pooledSE.toFixed(3) : '—';
    els['es-pooled-md'].textContent = grp.md != null ? grp.md.toFixed(3) : '—';

    const dispLabel = v => v === 'sd' ? 'Standard Deviation' : v === 'se' ? 'Standard Error' : '95% CI (Upper)';
    els['es-g1-disp-label'].textContent = dispLabel(els['es-g1-disp'].value);
    els['es-g2-disp-label'].textContent = dispLabel(els['es-g2-disp'].value);

    const dGroup = grp.dFromGroups != null ? grp.dFromGroups : null;
    const hasGroup = dGroup != null && !isNaN(dGroup) && isFinite(dGroup);
    if (hasGroup && !esState.statTouched){
      els['es-stat-type'].value = 'd';
      els['es-stat-value'].value = dGroup.toFixed(3);
      els['es-stat-aux'].value = '';
      esState.statAutoFilled = true;
      refreshAuxField();
    } else if (!hasGroup && esState.statAutoFilled && !esState.statTouched){
      els['es-stat-value'].value = '';
      esState.statAutoFilled = false;
    }

    const dStat = statToD(els['es-stat-type'].value, els['es-stat-value'].value, els['es-stat-aux'].value);
    const hasStat = dStat != null && !isNaN(dStat) && isFinite(dStat);
    let d = null;
    if (esState.source === 'groups'){
      d = hasGroup ? dGroup : null;
    } else {
      d = hasStat ? dStat : null;
    }
    setMethodMessage('', '');

    if (d == null || isNaN(d) || !isFinite(d)){
      ['es-out-d','es-out-g','es-out-r','es-out-r2','es-out-f','es-out-z','es-out-or',
       'es-out-ovl','es-out-u3','es-out-cles','es-out-nnt','es-out-similar'].forEach(k => els[k].textContent = '—');
      ['es-out-d-desc','es-out-g-desc','es-out-r-desc','es-out-r2-desc','es-out-f-desc']
        .forEach(k => applyDesc(els[k], null));
      drawCurve(0, true);
      computeTarget(grp, null);
      renderCommonLanguage(null);
      if (els['es-d-slider-val']) els['es-d-slider-val'].textContent = '—';
      return;
    }

    const ad = Math.abs(d);
    const r  = d / Math.sqrt(d*d + 4);
    const r2 = (d*d) / (d*d + 4);
    const f  = d / 2;
    const z  = 0.5 * Math.log((1 + r) / (1 - r));
    const or = Math.exp(d * Math.PI / Math.sqrt(3));
    const ovl = 2 * (1 - normCdf(ad / 2));
    const u3 = 1 - normCdf(-d);
    const cles = normCdf(d / Math.sqrt(2));
    const nntDenom = 2*cles - 1;
    const nnt = nntDenom !== 0 ? 1 / nntDenom : Infinity;

    let gVal;
    if (grp.pooledN != null){
      gVal = d * (1 - 3/(4*grp.pooledN - 9));
    } else if (els['es-stat-type'].value === 'g'){
      gVal = Number(els['es-stat-value'].value);
    } else {
      gVal = d;
    }

    els['es-out-d'].textContent  = d.toFixed(3);
    applyDesc(els['es-out-d-desc'], descD(d));
    els['es-out-g'].textContent  = isFinite(gVal) ? gVal.toFixed(3) : '—';
    applyDesc(els['es-out-g-desc'], isFinite(gVal) ? descD(gVal) : null);
    els['es-out-r'].textContent  = r.toFixed(3);
    applyDesc(els['es-out-r-desc'], descR(r));
    els['es-out-r2'].textContent = (r2*100).toFixed(2) + '%';
    applyDesc(els['es-out-r2-desc'], descR2(r2));
    els['es-out-f'].textContent  = f.toFixed(3);
    applyDesc(els['es-out-f-desc'], descF(f));
    els['es-out-z'].textContent  = z.toFixed(3);
    els['es-out-or'].textContent = or > 999 ? or.toExponential(2) : or.toFixed(3);
    els['es-out-ovl'].textContent = (ovl*100).toFixed(2) + '%';
    els['es-out-u3'].textContent  = (u3*100).toFixed(2) + '%';
    els['es-out-cles'].textContent = (cles*100).toFixed(2) + '%';
    els['es-out-nnt'].textContent = isFinite(nnt) && Math.abs(nnt) < 1e4 ? Math.abs(nnt).toFixed(2) : '—';
    els['es-out-similar'].textContent = similarEffect(d);

    drawCurve(d, false);
    computeTarget(grp, d);
    renderCommonLanguage({ d, cles, u3, ovl, r });
    if (els['es-d-slider']) els['es-d-slider'].value = String(Math.max(-3, Math.min(3, d)));
    if (els['es-d-slider-val']) els['es-d-slider-val'].textContent = d.toFixed(1);
  }

  function renderCommonLanguage(payload){
    const summaryEl = els['es-cl-summary'];
    const u3LineEl = els['es-cl-u3line'];
    if (!summaryEl || !u3LineEl) return;

    if (!payload){
      summaryEl.textContent = 'Enter a valid effect size to generate a plain-English interpretation.';
      u3LineEl.textContent = 'The average person in Group 1 is above about — of Group 2 (Cohen\'s U₃).';
      return;
    }

    const d = payload.d;
    const direction = d >= 0 ? 'Group 1' : 'Group 2';
    const other = d >= 0 ? 'Group 2' : 'Group 1';
    const clesPct = (payload.cles * 100);
    const u3Pct = (payload.u3 * 100);
    summaryEl.innerHTML = `At this effect size (<strong>d = ${d.toFixed(3)}</strong>), a randomly selected person from <strong>${direction}</strong> is likely to score higher than a randomly selected person from <strong>${other}</strong> about <strong>${clesPct.toFixed(1)}%</strong> of the time.`;
    u3LineEl.innerHTML = `The average person in <strong>${direction}</strong> is above about <strong>${u3Pct.toFixed(1)}%</strong> of <strong>${other}</strong> (Cohen's U\u2083).`;
  }

  function drawCurve(d, blank){
    const W = 480, H = 220, padL = 18, padR = 18, padT = 10, padB = 24;
    const xMin = -5, xMax = 5;
    const xRange = xMax - xMin;
    const dShow = Math.max(-4, Math.min(4, d || 0));
    const yMax = normPdf(0);
    const innerW = W - padL - padR, innerH = H - padT - padB;
    const xS = x => padL + (x - xMin)/xRange * innerW;
    const yS = y => padT + innerH - (y / (yMax * 1.05)) * innerH;

    function pathFor(mu){
      let s = '';
      for (let i = 0; i <= 200; i++){
        const x = xMin + (xRange * i / 200);
        const y = normPdf(x - mu);
        s += (i===0?'M':'L') + xS(x).toFixed(2) + ',' + yS(y).toFixed(2) + ' ';
      }
      return s;
    }
    function areaFor(mu){
      let s = 'M' + xS(xMin).toFixed(2) + ',' + yS(0).toFixed(2) + ' ';
      for (let i = 0; i <= 200; i++){
        const x = xMin + (xRange * i / 200);
        const y = normPdf(x - mu);
        s += 'L' + xS(x).toFixed(2) + ',' + yS(y).toFixed(2) + ' ';
      }
      s += 'L' + xS(xMax).toFixed(2) + ',' + yS(0).toFixed(2) + ' Z';
      return s;
    }

    let ticks = '';
    for (let t = -4; t <= 4; t++){
      ticks += '<line x1="' + xS(t) + '" y1="' + (padT+innerH) + '" x2="' + xS(t) + '" y2="' + (padT+innerH+4) + '" stroke="#A8A29E" stroke-width="1"/>';
      ticks += '<text x="' + xS(t) + '" y="' + (padT+innerH+18) + '" text-anchor="middle" font-family="IBM Plex Mono,monospace" font-size="10" fill="#6B6B6B">' + t + '</text>';
    }
    const meanG2 = '<line x1="' + xS(0) + '" y1="' + yS(yMax) + '" x2="' + xS(0) + '" y2="' + yS(0) + '" stroke="#7a8a9a" stroke-width="1" stroke-dasharray="3,3"/>';
    const meanG1 = '<line x1="' + xS(dShow) + '" y1="' + yS(yMax) + '" x2="' + xS(dShow) + '" y2="' + yS(0) + '" stroke="#3d7550" stroke-width="1" stroke-dasharray="3,3"/>';

    const svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg">'
      + '<line x1="' + padL + '" y1="' + (padT+innerH) + '" x2="' + (W-padR) + '" y2="' + (padT+innerH) + '" stroke="#1A1A1A" stroke-width="1"/>'
      + ticks
      + '<path d="' + areaFor(0) + '" fill="#7a8a9a" fill-opacity="0.18"/>'
      + '<path d="' + areaFor(dShow) + '" fill="#3d7550" fill-opacity="0.18"/>'
      + meanG2 + meanG1
      + '<path d="' + pathFor(0) + '" fill="none" stroke="#7a8a9a" stroke-width="1.6"/>'
      + '<path d="' + pathFor(dShow) + '" fill="none" stroke="#3d7550" stroke-width="1.6"/>'
      + '</svg>';
    els['es-curve'].innerHTML = svg;
    if (blank){
      els['es-curve-caption'].innerHTML = 'Enter a value to visualise the effect.';
    } else {
      els['es-curve-caption'].innerHTML = 'Curves shifted by Cohen’s <em>d</em> = <strong>' + d.toFixed(3) + '</strong>. Greater separation indicates a larger effect.';
    }
  }

  function computeTarget(grp, d){
    const t = els['es-target'].value === '' ? null : Number(els['es-target'].value);
    const can = grp.m1 != null && grp.m2 != null && grp.sd1 != null && grp.sd2 != null && grp.sd1 > 0 && grp.sd2 > 0 && t != null;
    if (!can){
      ['es-tgt-g1-above','es-tgt-g2-above','es-tgt-share-g1','es-tgt-share-g2','es-tgt-ratio']
        .forEach(k => els[k].textContent = '—');
      return;
    }
    const z1 = (t - grp.m1)/grp.sd1, z2 = (t - grp.m2)/grp.sd2;
    const above1 = 1 - normCdf(z1), above2 = 1 - normCdf(z2);
    const dens1 = normPdf(z1) / grp.sd1, dens2 = normPdf(z2) / grp.sd2;
    const total = dens1 + dens2;
    const share1 = total > 0 ? dens1/total : null, share2 = total > 0 ? dens2/total : null;
    const ratio = share1 != null && share2 != null && share1 > 0 && share2 > 0
      ? Math.min(share1/share2, share2/share1) : null;
    els['es-tgt-g1-above'].textContent = fmtPct(above1);
    els['es-tgt-g2-above'].textContent = fmtPct(above2);
    els['es-tgt-share-g1'].textContent = fmtPct(share1);
    els['es-tgt-share-g2'].textContent = fmtPct(share2);
    els['es-tgt-ratio'].textContent = ratio != null ? ratio.toFixed(3) : '—';
  }

  const watched = ['es-stat-type','es-stat-value','es-stat-aux',
    'es-g1-n','es-g1-mean','es-g1-disp','es-g1-disp-val',
    'es-g2-n','es-g2-mean','es-g2-disp','es-g2-disp-val',
    'es-target'];
  watched.forEach(id => {
    const el = els[id];
    if (!el) return;
    el.addEventListener('input', () => {
      if (id === 'es-stat-type' || id === 'es-stat-value' || id === 'es-stat-aux'){
        esState.statTouched = true;
        esState.statAutoFilled = false;
        esState.source = 'stat';
      }
      compute();
    });
    el.addEventListener('change', () => {
      if (id === 'es-stat-type' || id === 'es-stat-value' || id === 'es-stat-aux'){
        esState.statTouched = true;
        esState.statAutoFilled = false;
        esState.source = 'stat';
      }
      compute();
    });
  });
  function switchEffectMode(mode){
    const tab = document.querySelector('#effectsize .es-tab[data-mode="' + mode + '"]');
    if (!tab) return;
    document.querySelectorAll('#effectsize .es-tab').forEach(b => b.classList.toggle('is-active', b === tab));
    document.querySelectorAll('#effectsize .es-mode-pane').forEach(p => p.classList.toggle('is-active', p.dataset.pane === mode));
  }
  function switchWorkspaceView(view){
    const tab = document.querySelector('#effectsize .es-view-tab[data-view="' + view + '"]');
    if (!tab) return;
    document.querySelectorAll('#effectsize .es-view-tab').forEach(b => b.classList.toggle('is-active', b === tab));
    document.querySelectorAll('#effectsize .es-view-pane').forEach(p => p.classList.toggle('is-active', p.dataset.viewPane === view));
  }

  function clearGroupData(){
    ['es-g1-n','es-g1-mean','es-g1-disp-val','es-g2-n','es-g2-mean','es-g2-disp-val','es-target']
      .forEach(id => { if (els[id]) els[id].value = ''; });
    compute();
  }

  function clearAllInputs(){
    ['es-stat-value','es-stat-aux','es-g1-n','es-g1-mean','es-g1-disp-val',
     'es-g2-n','es-g2-mean','es-g2-disp-val','es-target']
      .forEach(id => { if (els[id]) els[id].value = ''; });
    els['es-stat-type'].value = 'd';
    els['es-g1-disp'].value = 'sd';
    els['es-g2-disp'].value = 'sd';
    esState.statTouched = false;
    esState.statAutoFilled = false;
    compute();
  }

  els['es-clear'].addEventListener('click', clearGroupData);
  els['es-reset-all'].addEventListener('click', clearAllInputs);
  els['es-method-warning'].addEventListener('click', e => {
    const action = e.target.closest('[data-es-action]')?.dataset.esAction;
    if (!action) return;
    if (action === 'open-groups') switchEffectMode('groups');
    if (action === 'clear-groups') clearGroupData();
    if (action === 'clear-all') clearAllInputs();
  });
  if (els['es-apply-groups']){
    els['es-apply-groups'].addEventListener('click', () => {
      esState.source = 'groups';
      switchEffectMode('groups');
      compute();
    });
  }
  if (els['es-target-explain-toggle'] && els['es-target-explain']){
    els['es-target-explain-toggle'].addEventListener('click', () => {
      const isHidden = els['es-target-explain'].hasAttribute('hidden');
      if (isHidden) els['es-target-explain'].removeAttribute('hidden');
      else els['es-target-explain'].setAttribute('hidden', '');
      els['es-target-explain-toggle'].setAttribute('aria-expanded', String(isHidden));
    });
  }
  if (els['es-d-slider']){
    const onSlide = () => {
      const d = Number(els['es-d-slider'].value);
      els['es-stat-type'].value = 'd';
      els['es-stat-value'].value = d.toFixed(1);
      els['es-stat-aux'].value = '';
      esState.statTouched = true;
      esState.statAutoFilled = false;
      esState.source = 'stat';
      if (els['es-d-slider-val']) els['es-d-slider-val'].textContent = d.toFixed(1);
      switchEffectMode('stat');
      compute();
    };
    els['es-d-slider'].addEventListener('input', onSlide);
    els['es-d-slider'].addEventListener('change', onSlide);
  }

  // Tab wiring
  document.querySelectorAll('#effectsize .es-tab').forEach(btn => {
    btn.addEventListener('click', () => {
      switchEffectMode(btn.dataset.mode);
      compute();
    });
  });
  document.querySelectorAll('#effectsize .es-view-tab').forEach(btn => {
    btn.addEventListener('click', () => switchWorkspaceView(btn.dataset.view));
  });

  refreshAuxField();
  compute();
})();
