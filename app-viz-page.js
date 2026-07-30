/* =====================================================================
   Score Charts - page module (#charts)

   ONE CHART PER TEST: each battery group (test family) gets a single
   card, and every trial/subtest in it is a ROW of that card's chart -
   dot + CI whisker on a shared native-metric axis, with the row's
   classification bands drawn as a strip behind it. Rows within one
   instrument share a metric, which is what makes the shared axis safe;
   a family that mixes metrics gets one sub-panel per metric inside the
   same card, never two metrics on one axis.

   The page is a VIEW of the Score Tables state: every number a row
   prints (percentile, base rate, classification, CI) is produced by the
   same functions the table itself calls, so screen and chart cannot
   disagree. Settings (classification scheme, CI level, premorbid
   comparison) are read live from the Score Tables controls - this page
   deliberately has none of its own.

   Row types, following the table's own membership machinery:
     - standard / t / scaled / z rows -> banded strip, dot, printed-CI
       whisker (the whisker is parsed back from getBatteryCiHtml's own
       output, so the table's rounding convention propagates here with
       no second site to keep in step)
     - base-rate rows (batteryBaseRateEntry) -> a step sparkline of the
       published cumulative table with the patient's span marked
     - higherIsWorse rows keep the axis position of the score AS
       OBTAINED, but that row's band strip is reversed - the table's
       convention: percentile as obtained, classification (and so the
       bands) describes performance
     - raw rows (rowScoreType 'raw') are listed, not plotted: toZ('raw')
       is null by design and a raw score has no position on any
       normative axis (see the metric:'raw' note in data.js)
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
     to each panel's native metric via the same fromZ used everywhere else. */
  const VIZ_BAND_EDGES_SS = [70, 80, 90, 110, 120, 130];

  /* One shared column layout for every row chart, so panels align. */
  const SVG_W = 940;
  const COL = { nameEnd: 205, scoreMid: 228, plotX: 262, plotW: 478, rightX: 754 };
  const HEADER_H = 20, ROW_H = 28, AXIS_H = 24;

  function vizX(v, axis){
    const t = (v - axis.min) / (axis.max - axis.min);
    return COL.plotX + Math.min(1, Math.max(0, t)) * COL.plotW;
  }
  function vizLevels(cls){
    return cls === 'wechsler' ? WECHSLER_LEVELS : AACN_LEVELS;
  }
  function vizTruncate(name){
    return name.length > 28 ? name.slice(0, 27) + '…' : name;
  }
  /* Bounds are joined with U+2013, which cannot collide with a negative
     sign, so the printed interval parses back unambiguously. */
  function vizParseCi(ciText){
    if (!ciText) return null;
    const parts = ciText.split('–').map(parseFloat);
    if (parts.length !== 2 || parts.some(isNaN)) return null;
    return { lo: parts[0], hi: parts[1] };
  }
  function vizRowFacts(r, type, ciLevel, ciText, pctCell, details){
    const bits = [`${scoreTypeLabel(type)} ${r.score}`];
    if (ciText) bits.push(`${ciLevel}% CI ${ciText}`);
    if (pctCell) bits.push(`${pctCell.kind === 'baseRate' ? BAT_BASERATE_LABEL : BAT_PCT_LABEL} ${pctCell.text}`);
    if (details && details.text) bits.push(details.text);
    return bits.join(' · ');
  }

  /* Column headings inside the SVG so they scale with the columns. */
  function vizHeaderSvg(pctLabel){
    return `<text class="viz-col-head" x="${COL.scoreMid}" y="13" text-anchor="middle">Score</text>` +
           `<text class="viz-col-head" x="${COL.rightX}" y="13">${pctLabel} · Classification</text>`;
  }

  /* ---------- metric panel: one row per subtest on a shared axis ------- */
  function vizMetricPanelSvg(rows, type, cls, ciLevel){
    const axis = VIZ_AXIS[type];
    const levels = vizLevels(cls);
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;
    const anyReversed = rows.some(r => r.higherIsWorse);

    // Band edges in the native metric, clipped to the axis range. The
    // boundaries are symmetric about the metric midline, so a reversed
    // (error-measure) row reverses the colour order only.
    const edges = [axis.min]
      .concat(VIZ_BAND_EDGES_SS.map(e => fromZ((e - 100) / 15, type))
        .filter(e => e > axis.min && e < axis.max))
      .concat([axis.max]);

    let svg = vizHeaderSvg(BAT_PCT_LABEL);

    rows.forEach((r, i) => {
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const score = parseFloat(r.score);
      const ciText = ciLevel !== 'off' ? getBatteryCiHtml(score, r, ciLevel) : '';
      const pctCell = batteryRowPctCell(r);
      const details = batteryClassificationDetails(r, cls);

      let row = `<title>${escapeHtml(r.name)}: ${escapeHtml(vizRowFacts(r, type, ciLevel, ciText, pctCell, details))}${r.higherIsWorse ? ' (error measure - bands reversed)' : ''}</title>`;
      row += `<text class="viz-row-name" x="${COL.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name))}${r.higherIsWorse ? ' ↓' : ''}</text>`;
      row += `<text class="viz-row-score" x="${COL.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${escapeHtml(String(r.score))}</text>`;

      for (let b = 0; b < edges.length - 1; b++){
        const ci = r.higherIsWorse ? (edges.length - 2 - b) : b;
        const x0 = vizX(edges[b], axis), x1 = vizX(edges[b + 1], axis);
        const level = levels[ci] || {};
        row += `<rect class="viz-band" x="${x0.toFixed(1)}" y="${rowY + 3}" width="${(x1 - x0).toFixed(1)}" height="${ROW_H - 6}" fill="${DESC_COLOURS[ci]}">` +
               `<title>${escapeHtml(level.label || '')} (standard-score ${escapeHtml(level.range || '')})</title></rect>`;
      }

      const ci = vizParseCi(ciText);
      if (ci){
        const x0 = vizX(ci.lo, axis).toFixed(1), x1 = vizX(ci.hi, axis).toFixed(1);
        row += `<line class="viz-whisker" x1="${x0}" y1="${midY}" x2="${x1}" y2="${midY}"/>` +
               `<line class="viz-whisker" x1="${x0}" y1="${midY - 5}" x2="${x0}" y2="${midY + 5}"/>` +
               `<line class="viz-whisker" x1="${x1}" y1="${midY - 5}" x2="${x1}" y2="${midY + 5}"/>`;
      }

      row += `<circle class="viz-dot" cx="${vizX(score, axis).toFixed(1)}" cy="${midY}" r="5.5"/>`;

      const right = [pctCell ? pctCell.text : '–', details && details.text ? details.text : '–'].join(' · ');
      row += `<text class="viz-row-fact" x="${COL.rightX}" y="${midY + 4.5}">${escapeHtml(right)}</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    // Shared axis under the last row.
    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL.plotX}" y1="${axisY + 2}" x2="${COL.plotX + COL.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizX(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
    }

    return { svg: `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="${escapeAttr(scoreTypeLabel(type))} measures">${svg}</svg>`, anyReversed };
  }

  /* ---------- base-rate panel: one published step sparkline per span --- */
  function vizBaseRatePanelSvg(rows){
    // Shared span domain across the family so the sparklines align.
    let lo = Infinity, hi = -Infinity;
    const entries = rows.map(r => batteryBaseRateEntry(r));
    entries.forEach(e => {
      Object.keys(e.baseRates).map(Number).filter(Number.isFinite).forEach(s => {
        if (s < lo) lo = s;
        if (s > hi) hi = s;
      });
    });
    if (!isFinite(lo)) return { svg: '' };
    const x0dom = lo, x1dom = hi + 1;
    const xOf = s => COL.plotX + (s - x0dom) / (x1dom - x0dom) * COL.plotW;
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;

    let svg = vizHeaderSvg(BAT_BASERATE_LABEL);

    rows.forEach((r, i) => {
      const entry = entries[i];
      const rowY = HEADER_H + i * ROW_H;
      const yOf = p => rowY + (ROW_H - 3) - (p / 100) * (ROW_H - 6);
      const v = parseFloat(r.score);
      const cls = document.getElementById('bat-class')?.value || 'wechsler';
      const details = batteryClassificationDetails(r, cls);
      const spans = Object.keys(entry.baseRates).map(Number).filter(Number.isFinite).sort((a, b) => a - b);

      const hasScore = Number.isFinite(v) && Number.isInteger(v);
      const br = hasScore ? baseRateAtOrAbove(entry, v) : null;
      let row = `<title>${escapeHtml(r.name)}: span ${escapeHtml(String(r.score))}, ${BAT_BASERATE_LABEL.toLowerCase()} ${hasScore ? fmtBaseRate(br) : '–'}% obtained this span or higher</title>`;
      row += `<text class="viz-row-name" x="${COL.nameEnd}" y="${rowY + ROW_H / 2 + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name))}</text>`;
      row += `<text class="viz-row-score" x="${COL.scoreMid}" y="${rowY + ROW_H / 2 + 4.5}" text-anchor="middle">${escapeHtml(String(r.score))}</text>`;

      let d = '';
      spans.forEach((s, k) => {
        const x = xOf(s), y = yOf(entry.baseRates[s]);
        d += (k === 0 ? `M ${x.toFixed(1)} ${y.toFixed(1)}` : ` H ${x.toFixed(1)} V ${y.toFixed(1)}`);
      });
      d += ` H ${xOf(x1dom).toFixed(1)}`;
      row += `<path class="viz-step" d="${d}"/>`;

      if (hasScore){
        const px = xOf(Math.min(Math.max(v, x0dom), x1dom));
        row += `<circle class="viz-dot" cx="${px.toFixed(1)}" cy="${yOf(br).toFixed(1)}" r="5.5"/>`;
      }

      const right = [hasScore ? fmtBaseRate(br) : '–', details && details.text ? details.text : '–'].join(' · ');
      row += `<text class="viz-row-fact" x="${COL.rightX}" y="${rowY + ROW_H / 2 + 4.5}">${escapeHtml(right)}</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL.plotX}" y1="${axisY + 2}" x2="${COL.plotX + COL.plotW}" y2="${axisY + 2}"/>`;
    const step = (x1dom - x0dom) > 14 ? 2 : 1;
    for (let s = x0dom; s <= x1dom - 1; s += step){
      const x = xOf(s).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${s}</text>`;
    }

    return { svg: `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="Published base rates by span">${svg}</svg>` };
  }

  /* ---------- family card: all of a test's trials in one place --------- */
  function vizFamilyCard(title, rows, cls, ciLevel){
    // Partition the family's rows: base-rate lookups, plottable metrics
    // (kept apart per metric - never two metrics on one axis), raw rest.
    const brRows = [], rawRows = [], byType = new Map();
    for (const r of rows){
      if (batteryBaseRateEntry(r)) { brRows.push(r); continue; }
      const type = rowScoreType(r);
      if (VIZ_AXIS[type]){
        if (!byType.has(type)) byType.set(type, []);
        byType.get(type).push(r);
      } else {
        rawRows.push(r);
      }
    }

    let body = '', anyReversed = false;
    const panelCount = byType.size + (brRows.length ? 1 : 0);
    for (const [type, tRows] of byType){
      if (panelCount > 1) body += `<div class="viz-sub-caption">${escapeHtml(scoreTypeLabel(type))}s</div>`;
      const panel = vizMetricPanelSvg(tRows, type, cls, ciLevel);
      body += panel.svg;
      anyReversed = anyReversed || panel.anyReversed;
    }
    if (brRows.length){
      if (panelCount > 1) body += `<div class="viz-sub-caption">Longest span (published base rates)</div>`;
      body += vizBaseRatePanelSvg(brRows).svg;
    }
    if (rawRows.length){
      body += rawRows.map(r => {
        const score = parseFloat(r.score);
        const ciText = ciLevel !== 'off' ? getBatteryCiHtml(score, r, ciLevel) : '';
        return `<div class="viz-raw-row"><b>${escapeHtml(r.name)}</b> — raw ${escapeHtml(String(r.score))}` +
               `${ciText ? `, ${ciLevel}% CI ${ciText}` : ''} · not plotted (no standardised metric; percentile and classification are not derived)</div>`;
      }).join('');
    }

    const notes = [];
    if (anyReversed) notes.push('↓ error measure: a higher score means more errors. Percentiles are as obtained; the classification bands describe performance, so that row’s bands run reversed.');
    if (brRows.length) notes.push('Step lines show the published percentage of the normative sample obtaining each span or higher (WAIS-IV Manual, Tables C.4–C.5). A HIGH base rate means a common, and therefore lower, score.');

    return `<div class="viz-card">
      <div class="viz-card-name">${escapeHtml(title)}</div>
      ${body}
      ${notes.length ? `<div class="viz-card-note">${notes.join('<br>')}</div>` : ''}
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

    // One card per test family, in table order; ungrouped rows form one card.
    const families = new Map();
    for (const r of rows){
      const gKey = batteryGroupKeyOf(r);
      if (!families.has(gKey)) families.set(gKey, []);
      families.get(gKey).push(r);
    }
    let html = '';
    for (const [gKey, fRows] of families){
      const title = gKey ? (batteryGroupLabel(gKey) || 'Untitled test') : 'Ungrouped measures';
      html += vizFamilyCard(title, fRows, cls, ciLevel);
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
