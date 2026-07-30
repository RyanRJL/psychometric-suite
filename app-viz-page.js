/* =====================================================================
   Score Charts - page module (#charts)

   One chart per measure, drawn in the measure's NATIVE metric - no score
   is converted for display. The page is a VIEW of the Score Tables state:
   every number a card prints (percentile, base rate, classification, CI)
   is produced by the same functions the table itself calls, so screen and
   chart cannot disagree. Settings (classification scheme, CI level,
   premorbid comparison) are read live from the Score Tables controls for
   the same reason - this page deliberately has none of its own.

   Chart membership follows the table's own machinery:
     - standardised rows (standard / t / scaled / z) -> banded axis card
     - base-rate rows (batteryBaseRateEntry)          -> published step chart
     - raw rows (rowScoreType 'raw')                  -> no chart; score + CI
       only, because toZ('raw') is null by design and a raw score has no
       position on any normative axis (see the metric:'raw' note in data.js)
     - higherIsWorse rows keep the axis position of the score AS OBTAINED,
       but the classification bands are reflected, matching the table's
       convention: classification describes performance, percentile is
       reported as obtained.
   ===================================================================== */
(function(){
  'use strict';

  const section = document.getElementById('charts');
  if (!section) return;
  const grid = document.getElementById('viz-grid');
  const emptyEl = document.getElementById('viz-empty');
  const settingsEl = document.getElementById('viz-settings-line');
  if (!grid) return;

  /* Native axis per metric. Ranges cover the usable span of each scale
     (Wechsler FSIQ floor/ceiling 40-160, T 10-90, scaled 1-19, z +-4). */
  const VIZ_AXIS = {
    standard: { min: 40, max: 160, ticks: [40, 70, 100, 130, 160] },
    t:        { min: 10, max: 90,  ticks: [10, 30, 50, 70, 90] },
    scaled:   { min: 1,  max: 19,  ticks: [1, 4, 7, 10, 13, 16, 19] },
    z:        { min: -4, max: 4,   ticks: [-4, -2, 0, 2, 4] }
  };
  /* Classification band boundaries on the standard-score metric; converted
     to each card's native metric via the same fromZ used everywhere else. */
  const VIZ_BAND_EDGES_SS = [70, 80, 90, 110, 120, 130];

  const SVG_W = 560, PAD_L = 10, PAD_R = 10;

  function vizX(v, axis){
    const t = (v - axis.min) / (axis.max - axis.min);
    return PAD_L + Math.min(1, Math.max(0, t)) * (SVG_W - PAD_L - PAD_R);
  }
  function vizFmtNum(v){
    return (Math.round(v * 10) / 10).toString();
  }
  function vizLevels(cls){
    return cls === 'wechsler' ? WECHSLER_LEVELS : AACN_LEVELS;
  }

  /* The CI whisker draws exactly the interval the table PRINTS - parsed
     back from getBatteryCiHtml's own output rather than recomputed, so a
     change to the table's rounding or bounds convention propagates here
     without a second site to keep in step. Bounds are joined with U+2013,
     which cannot collide with a negative sign. */
  function vizParseCi(ciText){
    if (!ciText) return null;
    const parts = ciText.split('–').map(parseFloat);
    if (parts.length !== 2 || parts.some(isNaN)) return null;
    return { lo: parts[0], hi: parts[1] };
  }

  /* ---------- standardised-metric card: banded axis + dot + whisker ---- */
  function vizAxisSvg(r, type, cls, ciText){
    const axis = VIZ_AXIS[type];
    const score = parseFloat(r.score);
    const levels = vizLevels(cls);
    const bandTop = 22, bandH = 40, midY = bandTop + bandH / 2;

    // Band edges in the native metric, clipped to the axis range.
    const edges = [axis.min]
      .concat(VIZ_BAND_EDGES_SS.map(e => fromZ((e - 100) / 15, type))
        .filter(e => e > axis.min && e < axis.max))
      .concat([axis.max]);

    let svg = '';
    for (let i = 0; i < edges.length - 1; i++){
      /* On an error measure the classification runs against the score, so
         the band colouring reflects - the high end of the axis wears the
         low-performance colour. Boundaries are symmetric about the metric
         midline, so reversing the colour order IS the reflection. */
      const ci = r.higherIsWorse ? (edges.length - 2 - i) : i;
      const x0 = vizX(edges[i], axis), x1 = vizX(edges[i + 1], axis);
      const level = levels[ci] || {};
      svg += `<rect class="viz-band" x="${x0.toFixed(1)}" y="${bandTop}" width="${(x1 - x0).toFixed(1)}" height="${bandH}" fill="${DESC_COLOURS[ci]}">` +
             `<title>${escapeHtml(level.label || '')} (standard-score ${escapeHtml(level.range || '')})</title></rect>`;
    }

    // Axis ticks and labels.
    for (const tval of axis.ticks){
      const x = vizX(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${bandTop + bandH}" x2="${x}" y2="${bandTop + bandH + 4}"/>`;
      svg += `<text class="viz-tick-label" x="${x}" y="${bandTop + bandH + 16}" text-anchor="middle">${tval}</text>`;
    }

    // CI whisker - the printed interval, clamped into the drawable range.
    const ci = vizParseCi(ciText);
    if (ci){
      const x0 = vizX(ci.lo, axis).toFixed(1), x1 = vizX(ci.hi, axis).toFixed(1);
      svg += `<line class="viz-whisker" x1="${x0}" y1="${midY}" x2="${x1}" y2="${midY}"/>` +
             `<line class="viz-whisker" x1="${x0}" y1="${midY - 5}" x2="${x0}" y2="${midY + 5}"/>` +
             `<line class="viz-whisker" x1="${x1}" y1="${midY - 5}" x2="${x1}" y2="${midY + 5}"/>`;
    }

    // Score dot, ringed in the surface colour; value labelled above.
    const sx = vizX(score, axis);
    const offScale = score < axis.min || score > axis.max;
    svg += `<circle class="viz-dot" cx="${sx.toFixed(1)}" cy="${midY}" r="6"/>`;
    const lx = Math.min(SVG_W - PAD_R - 6, Math.max(PAD_L + 6, sx));
    svg += `<text class="viz-dot-label" x="${lx.toFixed(1)}" y="${bandTop - 7}" text-anchor="middle">${vizFmtNum(score)}${offScale ? ' (off scale)' : ''}</text>`;

    return `<svg class="viz-svg" viewBox="0 0 ${SVG_W} 96" role="img" aria-label="${escapeAttr(r.name)}: score ${vizFmtNum(score)} on the ${escapeAttr(scoreTypeLabel(type).toLowerCase())} scale">${svg}</svg>`;
  }

  /* ---------- base-rate card: the published cumulative table as steps ---- */
  function vizBaseRateSvg(r, entry){
    const spans = Object.keys(entry.baseRates).map(Number)
      .filter(Number.isFinite).sort((a, b) => a - b);
    if (!spans.length) return '';
    const v = parseFloat(r.score);
    const x0dom = spans[0], x1dom = spans[spans.length - 1] + 1;
    const top = 16, bottom = 116, H = 150;
    const xOf = s => PAD_L + 24 + (s - x0dom) / (x1dom - x0dom) * (SVG_W - PAD_L - PAD_R - 24);
    const yOf = p => bottom - (p / 100) * (bottom - top);

    let svg = '';
    for (const g of [0, 25, 50, 75, 100]){
      svg += `<line class="viz-grid-line" x1="${xOf(x0dom).toFixed(1)}" y1="${yOf(g).toFixed(1)}" x2="${xOf(x1dom).toFixed(1)}" y2="${yOf(g).toFixed(1)}"/>`;
      if (g % 50 === 0) svg += `<text class="viz-tick-label" x="${(xOf(x0dom) - 6).toFixed(1)}" y="${(yOf(g) + 3).toFixed(1)}" text-anchor="end">${g}</text>`;
    }

    // Descending step path: the base rate holds from each span to the next.
    let d = '';
    for (let i = 0; i < spans.length; i++){
      const x = xOf(spans[i]), y = yOf(entry.baseRates[spans[i]]);
      d += (i === 0 ? `M ${x.toFixed(1)} ${y.toFixed(1)}` : ` H ${x.toFixed(1)} V ${y.toFixed(1)}`);
    }
    d += ` H ${xOf(x1dom).toFixed(1)}`;
    svg += `<path class="viz-step" d="${d}"/>`;

    const labelEvery = spans.length > 12 ? 2 : 1;
    spans.forEach((s, i) => {
      if (i % labelEvery) return;
      svg += `<text class="viz-tick-label" x="${xOf(s).toFixed(1)}" y="${bottom + 16}" text-anchor="middle">${s}</text>`;
    });

    // Patient marker on the published step for their span.
    if (Number.isFinite(v) && Number.isInteger(v)){
      const br = baseRateAtOrAbove(entry, v);
      const px = xOf(Math.min(Math.max(v, x0dom), x1dom));
      svg += `<line class="viz-marker-line" x1="${px.toFixed(1)}" y1="${top}" x2="${px.toFixed(1)}" y2="${bottom}"/>` +
             `<circle class="viz-dot" cx="${px.toFixed(1)}" cy="${yOf(br).toFixed(1)}" r="6"/>` +
             `<text class="viz-dot-label" x="${Math.min(SVG_W - PAD_R - 10, px + 9).toFixed(1)}" y="${(yOf(br) - 9).toFixed(1)}">${fmtBaseRate(br)}</text>`;
    }

    return `<svg class="viz-svg viz-svg-tall" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="${escapeAttr(r.name)}: percentage of the normative sample obtaining each span or higher">${svg}</svg>`;
  }

  /* ---------- card assembly ------------------------------------------- */
  function vizFacts(parts){
    return parts.filter(Boolean).map(p =>
      `<span class="viz-fact"><span class="viz-fact-label">${p[0]}</span> ${p[1]}</span>`
    ).join('');
  }

  function vizCard(r, cls, ciLevel){
    const type = rowScoreType(r);
    const brEntry = batteryBaseRateEntry(r);
    const score = parseFloat(r.score);
    const ciText = ciLevel !== 'off' ? getBatteryCiHtml(score, r, ciLevel) : '';
    const pctCell = batteryRowPctCell(r);
    const details = batteryClassificationDetails(r, cls);

    let body, note = '';
    if (brEntry){
      body = vizBaseRateSvg(r, brEntry);
      note = 'Chart shows the published percentage of the normative sample obtaining each span or higher (WAIS-IV Manual, Tables C.4–C.5). A HIGH base rate means a common, and therefore lower, score.';
    } else if (VIZ_AXIS[type]){
      body = vizAxisSvg(r, type, cls, ciText);
      if (r.higherIsWorse){
        note = 'Error measure: a higher score means more errors. The percentile is reported as obtained; the classification bands describe performance, so they run reversed on this card.';
      }
    } else {
      /* raw (or unrecognised) - no normative axis exists, so nothing is
         drawn. The score and its raw-unit CI are still real numbers. */
      body = `<div class="viz-noplot">Not plotted — raw-score measure with no standardised metric. Percentile and classification are not derived.</div>`;
    }

    const facts = vizFacts([
      [scoreTypeLabel(type), Number.isFinite(score) ? escapeHtml(String(r.score)) : '–'],
      ciText ? [`${ciLevel}% CI`, ciText] : null,
      pctCell ? [pctCell.kind === 'baseRate' ? BAT_BASERATE_LABEL : BAT_PCT_LABEL, pctCell.text] : null,
      details.text ? ['Classification', details.html] : null
    ]);

    return `<div class="viz-card">
      <div class="viz-card-name">${escapeHtml(r.name)}</div>
      ${body}
      <div class="viz-card-facts">${facts}</div>
      ${note ? `<div class="viz-card-note">${note}</div>` : ''}
    </div>`;
  }

  /* ---------- page render --------------------------------------------- */
  function renderVizSettingsLine(cls, ciLevel){
    if (!settingsEl) return;
    const prem = getBatteryPremorbidComparison();
    const bits = [
      cls === 'wechsler' ? 'Wechsler classification' : 'Guilmette et al. (2020) classification',
      ciLevel === 'off' ? 'confidence intervals off' : `${ciLevel}% confidence intervals`,
      prem ? 'premorbid comparison on (asterisks as on Score Tables)' : 'premorbid comparison off'
    ];
    settingsEl.textContent = `Following Score Tables settings: ${bits.join(' · ')}. Change them on the Score Tables page.`;
  }

  function renderVizPage(){
    const cls = document.getElementById('bat-class')?.value || 'wechsler';
    const ciLevel = document.getElementById('bat-ci-level')?.value || 'off';
    renderVizSettingsLine(cls, ciLevel);

    /* Same membership rule as the APA export: named rows with a numeric
       score, example rows excluded - a chart of a seeded example could be
       mistaken for patient data. */
    const rows = batteryRows.filter(r =>
      r.name && !r.isExample && r.score !== '' && !isNaN(parseFloat(r.score)));

    if (!rows.length){
      grid.innerHTML = '';
      if (emptyEl) emptyEl.hidden = false;
      return;
    }
    if (emptyEl) emptyEl.hidden = true;

    let html = '', lastGroup = null;
    for (const r of rows){
      const gKey = batteryGroupKeyOf(r);
      if (gKey && gKey !== lastGroup){
        const label = batteryGroupLabel(gKey);
        if (label) html += `<div class="viz-group-head">${escapeHtml(label)}</div>`;
        lastGroup = gKey;
      } else if (!gKey && lastGroup !== null){
        lastGroup = null;
        html += `<div class="viz-group-head">Ungrouped measures</div>`;
      }
      html += vizCard(r, cls, ciLevel);
    }
    grid.innerHTML = html;
  }

  /* Re-render whenever the page becomes visible, so it always reflects the
     current Score Tables state. Same pattern as syncTopnav: observe the
     section's class attribute rather than hooking every navigation path. */
  new MutationObserver(() => {
    if (section.classList.contains('active')) renderVizPage();
  }).observe(section, { attributes: true, attributeFilter: ['class'] });

  if (section.classList.contains('active')) renderVizPage();
})();
