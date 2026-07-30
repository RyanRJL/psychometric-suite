/* =====================================================================
   Score Charts - page module (#charts)

   ONE CHART PER TEST: each battery group (test family) gets a single
   card, and every trial/subtest in it is a ROW of that card's chart -
   dot + CI whisker on a shared native-metric axis, with the row's
   classification bands drawn as a strip behind it. Rows within one
   instrument share a metric, which is what makes the shared axis safe;
   a family that mixes metrics gets one sub-panel per metric inside the
   same card, never two metrics on one axis.

   Four blocks, one per source table: Score Tables (score against
   classification bands, with a view selector - native metric, percentile,
   standard score, or raw against the measure's own norms), Premorbid
   (predicted vs achieved against the prediction interval), Change
   Analysis (both testings against the selected method's reliable-change
   interval) and the SD Index (change in SD units against its
   significance band).

   The page is a VIEW of those tables' state: every number a row prints
   (percentile, base rate, classification, CI, RCI, outcome) is produced
   by the same functions the tables themselves call, so screen and chart
   cannot disagree. Settings - classification scheme, CI level, premorbid
   comparison, each method's confidence level, the SDI mode - are read
   live from those pages' own controls. This page deliberately has no
   settings of its own; its only two controls - how to show the Score
   Tables axis, and which RCI method to draw - change no state anywhere.

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
  /* Change rows print a PAIR of scores ("92→85") where a score card prints
     one, so they get their own geometry - the name ends earlier and the
     score column sits wider. Sharing COL made long names collide with the
     pair. */
  const COL_CHANGE = { nameEnd: 186, scoreMid: 232, plotX: 286, plotW: 454, rightX: 754 };
  const HEADER_H = 20, ROW_H = 28, AXIS_H = 24;

  function vizX(v, axis){
    const t = (v - axis.min) / (axis.max - axis.min);
    return COL.plotX + Math.min(1, Math.max(0, t)) * COL.plotW;
  }
  /* A value outside the axis range is CLAMPED to the edge by vizX, so
     without a marker a score 20 SD out is drawn indistinguishably from one
     just past the end. The caret says "this continues beyond the axis";
     the row's own figures give the actual value. */
  function vizDotSvg(value, axis, xFn, midY, cls){
    const off = value < axis.min ? -1 : value > axis.max ? 1 : 0;
    const cx = xFn(value, axis);
    let out = `<circle class="${cls || 'viz-dot'}${off ? ' viz-dot-off' : ''}" cx="${cx.toFixed(1)}" cy="${midY}" r="5.5"/>`;
    if (off) out += `<text class="viz-offscale" x="${(cx + off * 9).toFixed(1)}" y="${midY + 3.5}" text-anchor="${off > 0 ? 'start' : 'end'}">${off > 0 ? '\u25B8' : '\u25C2'}</text>`;
    return out;
  }

  function vizLevels(cls){
    return cls === 'wechsler' ? WECHSLER_LEVELS : AACN_LEVELS;
  }
  function vizTruncate(name, max){
    const n = max || 28;
    return name.length > n ? name.slice(0, n - 1) + '…' : name;
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

      row += vizDotSvg(score, axis, vizX, midY);

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

  /* ================== VIEW MODES for the Score Tables block =============
     'native' keeps every row on its entry metric (the default, and the only
     view in which nothing is converted). The other three answer "show me
     this a different way", each with its own honesty constraint:

       'percentile'  every plottable row on one percentile axis
       'standard'    every plottable row on the standard-score metric
       'raw'         raw scores against their own normative mean, in that
                     measure's own SD units - no percentile, no bands

     A base-rate row joins a converted axis through the percentile its
     PUBLISHED table gives (midpoint convention), which is exactly the route
     batteryClassificationDetails already takes to classify one. A raw row
     never joins a standardised axis, in either converted view. */
  const VIZ_SCORE_VIEWS = [
    ['native',     'Native metric'],
    ['percentile', 'Percentile'],
    ['standard',   'Standard score'],
    ['raw',        'Raw vs norms']
  ];
  let vizScoreView = 'native';
  const VIZ_SCORE_AXIS_NOTE = {
    native: '',
    percentile: 'Axis is the percentile rank, converted from each row’s own metric so the whole test shares one scale. A percentile scale is deliberately non-linear: it compresses the tails, so two clearly different low scores can sit close together, and confidence intervals become asymmetric. Scores are shown as entered.',
    standard: 'Axis is the standard-score metric (M 100, SD 15), converted from each row’s own metric so the whole test shares one scale. Scores are shown as entered.',
    raw: 'Axis is the raw score’s distance from that measure’s own normative mean, in that measure’s own SD units, with ±1 SD shaded — descriptive only. No percentile or classification is derived, because the shape of a raw-score distribution is not known from its mean and SD alone. Only measures whose stored norms are in raw units can appear.'
  };

  const VIZ_CONV_AXIS = {
    percentile: { min: 0, max: 100, ticks: [1, 10, 25, 50, 75, 90, 99] },
    standard:   { min: 40, max: 160, ticks: [40, 70, 100, 130, 160] }
  };
  function vizConvValue(z, view){
    return view === 'percentile' ? normCDF(z) * 100 : fromZ(z, 'standard');
  }
  /* Axis position for a row in a converted view, or null if it has none. */
  function vizConvertedPos(r, view){
    const brEntry = batteryBaseRateEntry(r);
    let z;
    if (brEntry){
      const pct = baseRatePercentile(brEntry, r.score);
      if (pct == null) return null;
      z = normInv(Math.min(99.99, Math.max(0.01, pct)) / 100);
    } else {
      z = toZ(r.score, rowScoreType(r));
    }
    return z == null ? null : vizConvValue(z, view);
  }

  function vizConvertedPanelSvg(rows, view, cls, ciLevel){
    const axis = VIZ_CONV_AXIS[view];
    const levels = vizLevels(cls);
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;
    const anyReversed = rows.some(r => r.higherIsWorse);
    const mixed = new Set(rows.map(r => rowScoreType(r))).size > 1;

    // Band edges: the standard-score boundaries, expressed in this view.
    const edges = [axis.min]
      .concat(VIZ_BAND_EDGES_SS.map(e => vizConvValue((e - 100) / 15, view)))
      .concat([axis.max]);

    let svg = `<text class="viz-col-head" x="${COL.scoreMid}" y="13" text-anchor="middle">Score</text>` +
              `<text class="viz-col-head" x="${COL.rightX}" y="13">${view === 'percentile' ? 'Percentile' : 'Standard score'} · Classification</text>`;

    rows.forEach((r, i) => {
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const type = rowScoreType(r);
      const pos = vizConvertedPos(r, view);
      const pctCell = batteryRowPctCell(r);
      const details = batteryClassificationDetails(r, cls);
      const ciText = ciLevel !== 'off' ? getBatteryCiHtml(parseFloat(r.score), r, ciLevel) : '';

      /* Truncate the NAME, then append the tag - truncating the joined
         string cut the tag itself in half ("(T-Scor…"), which reads as a
         corrupted metric rather than a shortened name. */
      const tag = mixed ? ` (${scoreTypeAbbr(type)})` : '';
      const shownName = vizTruncate(r.name, mixed ? 20 : 28) + tag;
      let row = `<title>${escapeHtml(r.name)}: ${escapeHtml(vizRowFacts(r, type, ciLevel, ciText, pctCell, details))}${r.higherIsWorse ? ' (error measure - bands reversed)' : ''}</title>`;
      row += `<text class="viz-row-name" x="${COL.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(shownName)}${r.higherIsWorse ? ' ↓' : ''}</text>`;
      row += `<text class="viz-row-score" x="${COL.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${escapeHtml(String(r.score))}</text>`;

      for (let b = 0; b < edges.length - 1; b++){
        const bi = r.higherIsWorse ? (edges.length - 2 - b) : b;
        const x0 = vizX(edges[b], axis), x1 = vizX(edges[b + 1], axis);
        const level = levels[bi] || {};
        row += `<rect class="viz-band" x="${x0.toFixed(1)}" y="${rowY + 3}" width="${Math.max(0, x1 - x0).toFixed(1)}" height="${ROW_H - 6}" fill="${DESC_COLOURS[bi]}">` +
               `<title>${escapeHtml(level.label || '')} (standard-score ${escapeHtml(level.range || '')})</title></rect>`;
      }

      /* The interval is converted bound-by-bound from the one the table
         prints, so it stays the same interval - asymmetric on a percentile
         axis because the axis is non-linear, not because it moved. */
      const ci = vizParseCi(ciText);
      if (ci && !batteryBaseRateEntry(r)){
        const zLo = toZ(ci.lo, type), zHi = toZ(ci.hi, type);
        if (zLo != null && zHi != null){
          const x0 = vizX(vizConvValue(zLo, view), axis).toFixed(1);
          const x1 = vizX(vizConvValue(zHi, view), axis).toFixed(1);
          row += `<line class="viz-whisker" x1="${x0}" y1="${midY}" x2="${x1}" y2="${midY}"/>` +
                 `<line class="viz-whisker" x1="${x0}" y1="${midY - 5}" x2="${x0}" y2="${midY + 5}"/>` +
                 `<line class="viz-whisker" x1="${x1}" y1="${midY - 5}" x2="${x1}" y2="${midY + 5}"/>`;
        }
      }

      row += vizDotSvg(pos, axis, vizX, midY);
      const shown = view === 'percentile' ? (pctCell ? pctCell.text : fmtPct(pos)) : String(Math.round(pos));
      row += `<text class="viz-row-fact" x="${COL.rightX}" y="${midY + 4.5}">${escapeHtml(shown)} · ${escapeHtml(details.text || '–')}</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL.plotX}" y1="${axisY + 2}" x2="${COL.plotX + COL.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizX(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
    }
    return { svg: `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="${view === 'percentile' ? 'Percentile' : 'Standard score'} view">${svg}</svg>`, anyReversed, mixed };
  }

  /* ---------- raw view: raw score against its own normative mean -------
     Which cell holds the raw score depends on the entry metric, and the
     two cases are not interchangeable. Where rowScoreType is 'raw' the
     Score column IS the raw score (RBANS List Recognition). Where a raw
     measure is reported on a standardised metric (every CVLT-C index,
     reportedAs z or t) the Score column holds that standardised score and
     the raw count sits in the Raw column. */
  function vizRawValueOf(r){
    const v = rowScoreType(r) === 'raw' ? r.score : r.raw;
    return (v === undefined || v === null || v === '' || isNaN(parseFloat(v))) ? null : parseFloat(v);
  }
  /* Normative mean/SD in RAW units - only where the stored statistics are
     declared raw. For a standardised measure m1/sd1 are on that metric, so
     a typed raw score cannot be placed against them at all. */
  function vizRawNormOf(r){
    if (!r.group || !r.name) return null;
    const db = typeof getMergedDB === 'function' ? getMergedDB() : null;
    const e = db && db[r.group] && db[r.group][r.name];
    if (!e || e.metric !== 'raw') return null;
    if (!Number.isFinite(e.m1) || !Number.isFinite(e.sd1) || e.sd1 <= 0) return null;
    const v = vizRawValueOf(r);
    return v === null ? null : { m1: e.m1, sd1: e.sd1, value: v };
  }

  function vizRawNormPanelSvg(rows){
    const axis = { min: -3.5, max: 3.5, ticks: [-3, -2, -1, 0, 1, 2, 3] };
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;
    let svg = `<text class="viz-col-head" x="${COL.scoreMid}" y="13" text-anchor="middle">Raw</text>` +
              `<text class="viz-col-head" x="${COL.rightX}" y="13">Normative M (SD) · Distance</text>`;

    rows.forEach((r, i) => {
      const n = vizRawNormOf(r);
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const sd = (n.value - n.m1) / n.sd1;

      let row = `<title>${escapeHtml(r.name)}: raw ${n.value} against normative M ${fmt(n.m1, 2)} SD ${fmt(n.sd1, 2)} — ${fmt(sd, 2)} SD from the mean (descriptive; no percentile derived)</title>`;
      row += `<text class="viz-row-name" x="${COL.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name))}</text>`;
      row += `<text class="viz-row-score" x="${COL.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${n.value}</text>`;
      // A band of +-1 SD around the mean, as an orientation aid only.
      const b0 = vizX(-1, axis), b1 = vizX(1, axis);
      row += `<rect class="viz-rci-band" x="${b0.toFixed(1)}" y="${rowY + 3}" width="${(b1 - b0).toFixed(1)}" height="${ROW_H - 6}"/>`;
      const zx = vizX(0, axis).toFixed(1);
      row += `<line class="viz-rci-expected" x1="${zx}" y1="${rowY + 3}" x2="${zx}" y2="${rowY + ROW_H - 3}"/>`;
      row += vizDotSvg(sd, axis, vizX, midY);
      row += `<text class="viz-row-fact" x="${COL.rightX}" y="${midY + 4.5}">${fmt(n.m1, 1)} (${fmt(n.sd1, 1)}) · ${fmt(sd, 2)} SD</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL.plotX}" y1="${axisY + 2}" x2="${COL.plotX + COL.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizX(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval > 0 ? '+' + tval : tval}</text>`;
    }
    return `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="Raw score against normative mean in SD units">${svg}</svg>`;
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
    const notes = [];

    if (vizScoreView === 'raw'){
      /* Raw view is its own thing: a raw score has no standardised axis, so
         the only honest common axis is distance from that measure's OWN
         normative mean, in its own SD units. Bands and percentiles are
         deliberately absent - see vizRawNormPanelSvg. */
      const plot = [], listed = [];
      for (const r of rows){
        (vizRawNormOf(r) ? plot : listed).push(r);
      }
      if (plot.length) body += vizRawNormPanelSvg(plot);
      /* One shared explanation per reason, rather than repeating a long
         sentence on every row - eight identical paragraphs read as an
         error list rather than as a note. */
      const noRaw = listed.filter(r => vizRawValueOf(r) === null);
      const noNorm = listed.filter(r => vizRawValueOf(r) !== null);
      if (noNorm.length){
        body += `<div class="viz-raw-row">Not plotted — held on a standardised metric, so a raw score cannot be placed against ` +
          `this measure’s normative mean and SD: ` +
          noNorm.map(r => `<b>${escapeHtml(r.name)}</b> (raw ${escapeHtml(String(vizRawValueOf(r)))})`).join(', ') + `.</div>`;
      }
      if (noRaw.length){
        body += `<div class="viz-raw-row">No raw score entered on Score Tables: ` +
          noRaw.map(r => `<b>${escapeHtml(r.name)}</b>`).join(', ') + `.</div>`;
      }
    } else if (vizScoreView !== 'native'){
      /* Converted views put every plottable row of the test on ONE axis,
         which is the point of asking for them. Base-rate rows join via the
         percentile their published table gives, exactly as the table's own
         classification does. */
      const conv = rows.filter(r => vizConvertedPos(r, vizScoreView) !== null);
      const excluded = rows.filter(r => vizConvertedPos(r, vizScoreView) === null);
      if (conv.length){
        const panel = vizConvertedPanelSvg(conv, vizScoreView, cls, ciLevel);
        body += panel.svg;
        anyReversed = panel.anyReversed;
        if (panel.mixed) notes.push('This test mixes score metrics, so each row carries its own metric tag. The Score column shows the score as entered; only the axis position is converted.');
      }
      if (excluded.length){
        body += `<div class="viz-raw-row">Not plotted on a converted axis — raw-score measures with no standardised metric: ` +
          excluded.map(r => `<b>${escapeHtml(r.name)}</b>`).join(', ') + `.</div>`;
      }
    } else {
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
      if (brRows.length) notes.push('Step lines show the published percentage of the normative sample obtaining each span or higher (WAIS-IV Manual, Tables C.4–C.5). A HIGH base rate means a common, and therefore lower, score.');
    }

    if (anyReversed) notes.unshift('↓ error measure: a higher score means more errors. Percentiles are as obtained; the classification bands describe performance, so that row’s bands run reversed.');

    return `<div class="viz-card">
      <div class="viz-card-name">${escapeHtml(title)}</div>
      ${body}
      ${notes.length ? `<div class="viz-card-note">${notes.join('<br>')}</div>` : ''}
    </div>`;
  }

  /* ================== CHANGE ANALYSIS (RCI methods) ====================
     One card per test group from the shared Change Analysis rows, one row
     per trial: Date 1 (open circle) and Date 2 (filled dot) in score
     units, against the selected method's reliable-change interval - the
     band |t2 - expected| < crit x SE, which is exactly the region where
     the method's table prints "No reliable change". Interval and outcome
     come from the same calc*Row functions and rciOutcome the tables call,
     at the confidence level set on that method's page. The method tabs
     only choose which method's results to DRAW - they change no state. */
  const VIZ_RCI_METHODS = [
    ['rci-basic',    'Basic (Jacobson & Truax)'],
    ['rci-practice', 'Practice-Adjusted (Iverson)'],
    ['rci-srb',      'SRB (McSweeney)'],
    ['rci-crawford', 'Crawford & Garthwaite']
  ];
  let vizRciMethod = 'rci-basic';

  function vizRciCalc(r, method){
    const calc = method === 'rci-basic'    ? calcBasicRow(r, method)
               : method === 'rci-practice' ? calcPracticeRow(r, method)
               : method === 'rci-srb'      ? calcSrbRow(r, method)
               :                             calcCrawfordRow(r, method);
    if (!calc) return null;
    const cv = rciState[method].cv;
    /* Same critical value rciOutcome uses: t-quantile for Crawford, z else. */
    const crit = (calc.df != null && isFinite(calc.df) && calc.df > 0)
      ? tInv(cv === 0.95 ? 0.975 : 0.95, calc.df)
      : (cv === 0.95 ? 1.96 : 1.645);
    const se = method === 'rci-basic' ? calc.sed
             : method === 'rci-practice' ? calc.sdiff
             : method === 'rci-srb' ? calc.see
             : calc.sePred;
    const expected = method === 'rci-basic' ? parseFloat(r.t1)
                   : method === 'rci-practice' ? parseFloat(r.t1) + (parseFloat(r.m2) - parseFloat(r.m1))
                   : calc.predicted;
    return { calc, expected, half: crit * se, outcome: rciOutcome(calc.rci, cv, calc.df) };
  }

  /* Data-driven axis (change rows carry no declared metric): pad the data
     range and pick a clean tick step. */
  function vizNiceAxis(lo, hi){
    if (!(hi > lo)){ lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.08;
    lo -= pad; hi += pad;
    const target = (hi - lo) / 5;
    const mag = Math.pow(10, Math.floor(Math.log10(target)));
    const step = [1, 2, 2.5, 5, 10].map(s => s * mag).find(s => s >= target) || 10 * mag;
    const ticks = [];
    for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9; t += step) ticks.push(Math.round(t * 100) / 100);
    return { min: lo, max: hi, ticks };
  }
  function vizXd(v, axis){
    const t = (v - axis.min) / (axis.max - axis.min);
    return COL_CHANGE.plotX + Math.min(1, Math.max(0, t)) * COL_CHANGE.plotW;
  }

  function vizRciPanelSvg(rows, method){
    const results = rows.map(r => vizRciCalc(r, method));
    let lo = Infinity, hi = -Infinity;
    rows.forEach((r, i) => {
      const res = results[i];
      if (!res) return;
      for (const v of [parseFloat(r.t1), parseFloat(r.t2), res.expected - res.half, res.expected + res.half]){
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    });
    if (!isFinite(lo)) return '';
    const axis = vizNiceAxis(lo, hi);
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;

    let svg = `<text class="viz-col-head" x="${COL_CHANGE.scoreMid}" y="13" text-anchor="middle">Scores</text>` +
              `<text class="viz-col-head" x="${COL_CHANGE.rightX}" y="13">RCI · Outcome</text>`;

    rows.forEach((r, i) => {
      const res = results[i];
      if (!res) return;
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);

      let row = `<title>${escapeHtml(r.name)}: ${escapeHtml(String(r.t1))} → ${escapeHtml(String(r.t2))}, RCI ${fmt(res.calc.rci, 2)}, p ${fmtP(res.calc.p)} — ${res.outcome.label}</title>`;
      row += `<text class="viz-row-name" x="${COL_CHANGE.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name, 24))}</text>`;
      row += `<text class="viz-row-score" x="${COL_CHANGE.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${escapeHtml(String(r.t1))} → ${escapeHtml(String(r.t2))}</text>`;

      const bx0 = vizXd(res.expected - res.half, axis), bx1 = vizXd(res.expected + res.half, axis);
      row += `<rect class="viz-rci-band" x="${bx0.toFixed(1)}" y="${rowY + 3}" width="${(bx1 - bx0).toFixed(1)}" height="${ROW_H - 6}"/>`;
      const ex = vizXd(res.expected, axis).toFixed(1);
      row += `<line class="viz-rci-expected" x1="${ex}" y1="${rowY + 3}" x2="${ex}" y2="${rowY + ROW_H - 3}"/>`;

      const x1 = vizXd(t1, axis), x2 = vizXd(t2, axis);
      row += `<line class="viz-link-line" x1="${x1.toFixed(1)}" y1="${midY}" x2="${x2.toFixed(1)}" y2="${midY}"/>`;
      row += `<circle class="viz-dot-open" cx="${x1.toFixed(1)}" cy="${midY}" r="5"/>`;
      row += `<circle class="viz-dot" cx="${x2.toFixed(1)}" cy="${midY}" r="5.5"/>`;

      row += `<text class="viz-row-fact" x="${COL_CHANGE.rightX}" y="${midY + 4.5}">${fmt(res.calc.rci, 2)} · ${escapeHtml(res.outcome.label)}</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL_CHANGE.plotX}" y1="${axisY + 2}" x2="${COL_CHANGE.plotX + COL_CHANGE.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizXd(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
    }
    return `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="Reliable change">${svg}</svg>`;
  }

  function vizChangeBlockHtml(){
    const all = RCI_SHARED_ROWS.filter(r => r.name && !r.isExample);
    if (!all.length) return '';
    const method = vizRciMethod;
    const st = rciState[method];
    const chartable = all.filter(r => vizRciCalc(r, method));

    const tabs = `<div class="viz-method-tabs">` + VIZ_RCI_METHODS.map(([m, label]) =>
      `<button type="button" class="viz-method-tab${m === method ? ' is-active' : ''}" data-viz-method="${m}">${label}</button>`
    ).join('') + `</div>`;

    let cards = '';
    if (!chartable.length){
      cards = `<div class="viz-empty"><p>Rows are entered on Change Analysis, but none has the values this method needs yet (scores at both dates plus the normative parameters). Complete them on the ${escapeHtml(VIZ_RCI_METHODS.find(x => x[0] === method)[1])} page.</p></div>`;
    } else {
      const families = new Map();
      for (const r of chartable){
        const g = r.group || '';
        if (!families.has(g)) families.set(g, []);
        families.get(g).push(r);
      }
      for (const [g, fRows] of families){
        const title = g ? caGroupDisplay(g) : 'Ungrouped measures';
        cards += `<div class="viz-card"><div class="viz-card-name">${escapeHtml(title)}</div>` +
          vizRciPanelSvg(fRows, method) +
          `</div>`;
      }
    }
    const cvLabel = st.cv === 0.95 ? '95%' : '90%';
    const legend = `<p class="viz-legend"><span class="viz-legend-open">○</span> ${escapeHtml(st.d1 || 'Date 1')} · <span class="viz-legend-filled">●</span> ${escapeHtml(st.d2 || 'Date 2')} · shaded band: no reliable change at the ${cvLabel} level set on that method's page. Outcomes state significance only, not direction.</p>`;
    return `<div class="viz-section-head">Change Analysis</div>${tabs}${legend}${cards}`;
  }

  /* ================== PREMORBID: predicted vs achieved ==================
     The Estimates tab already draws its own model-comparison forest plot,
     so this block charts the question the other two tabs ask: does the
     achieved score fall below what the model predicted, and how unusual
     is that shortfall.

     Predicted values, intervals, differences and base rates are READ FROM
     THE CELLS THE TABLE PRINTS rather than recomputed. Both tabs round the
     estimate and the margin separately (round(estimate) +- round(z x SEE))
     precisely so the bounds stay symmetric about the printed value, and
     the two tabs disagreed on ~21% of inputs before that convention was
     shared - so re-deriving them here would create a third site to keep in
     step. The DOM is the single source of truth for what the clinician is
     looking at. */
  function vizPreCellNum(id){
    const el = document.getElementById(id);
    if (!el) return null;
    const v = parseFloat((el.textContent || '').replace(/[+]/g, ''));
    return Number.isFinite(v) ? v : null;
  }
  function vizPreAchieved(idx){
    const a = preState.achieved[idx];
    const v = (a === undefined || a === null || a === '') ? NaN : parseFloat(a);
    return Number.isFinite(v) ? v : null;
  }

  /* ToPF-predicted WAIS-IV / WMS-IV rows, straight off the predict table. */
  function vizTopfPredictRows(){
    const out = [];
    for (const c of [].concat(WAIS_COEF, WMS_COEF)){
      const pred = vizPreCellNum('pred-' + c.idx);
      if (pred == null) continue;
      out.push({
        label: c.label || c.idx,
        predicted: pred,
        lo: vizPreCellNum('pred-' + c.idx + '-lo'),
        hi: vizPreCellNum('pred-' + c.idx + '-hi'),
        achieved: vizPreAchieved(c.idx),
        diffText: (document.getElementById('diff-' + c.idx) || {}).textContent || '',
        brText: (document.getElementById('br-' + c.idx) || {}).textContent || ''
      });
    }
    return out;
  }
  /* OPIE-4 rows. This table is rebuilt wholesale on every render and its
     cells carry no ids, so the values come from the cached row objects the
     render itself stored (preState.opieRows) with the same rounding
     convention applied. */
  function vizOpiePredictRows(){
    const rows = preState.opieRows || [];
    const mult = preCiMult();
    const out = [];
    for (const r of rows){
      if (r.val == null || !Number.isFinite(r.val)) continue;
      const margin = Math.round(mult * r.see);
      const idx = opieModelIndex(r.key);
      const ach = vizPreAchieved(idx);
      const diff = ach == null ? null : Math.round(ach - r.val);
      out.push({
        label: r.label,
        predicted: r.val,
        lo: Math.round(r.val) - margin,
        hi: Math.round(r.val) + margin,
        achieved: ach,
        diffText: diff == null ? '' : (diff > 0 ? '+' : '') + diff,
        brText: diff == null ? '' : opieBaseRateFor(r.key, diff)
      });
    }
    return out;
  }

  function vizPrePanelSvg(rows){
    // These are all standard-score indices, so the standard axis applies
    // directly; a reference line marks the population mean of 100.
    const axis = VIZ_AXIS.standard;
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;
    let svg = `<text class="viz-col-head" x="${COL_CHANGE.scoreMid}" y="13" text-anchor="middle">Pred → Ach</text>` +
              `<text class="viz-col-head" x="${COL_CHANGE.rightX}" y="13">Difference · Base rate</text>`;

    const refX = vizXd(100, axis).toFixed(1);
    svg += `<line class="viz-rci-expected" x1="${refX}" y1="${HEADER_H}" x2="${refX}" y2="${HEADER_H + rows.length * ROW_H}"/>`;

    rows.forEach((row, i) => {
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const pRound = Math.round(row.predicted);

      let g = `<title>${escapeHtml(row.label)}: predicted ${pRound}${row.lo != null ? ` (${row.lo}–${row.hi})` : ''}` +
              `${row.achieved != null ? `, achieved ${row.achieved}, difference ${row.diffText}${row.brText ? `, base rate ${row.brText}` : ''}` : ', no achieved score entered'}</title>`;
      g += `<text class="viz-row-name" x="${COL_CHANGE.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(row.label, 26))}</text>`;
      g += `<text class="viz-row-score" x="${COL_CHANGE.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${pRound}${row.achieved != null ? ` → ${row.achieved}` : ''}</text>`;

      if (row.lo != null && row.hi != null){
        const x0 = vizXd(row.lo, axis).toFixed(1), x1 = vizXd(row.hi, axis).toFixed(1);
        g += `<rect class="viz-rci-band" x="${x0}" y="${rowY + 3}" width="${(x1 - x0).toFixed(1)}" height="${ROW_H - 6}"/>` +
             `<line class="viz-whisker" x1="${x0}" y1="${midY}" x2="${x1}" y2="${midY}"/>` +
             `<line class="viz-whisker" x1="${x0}" y1="${midY - 5}" x2="${x0}" y2="${midY + 5}"/>` +
             `<line class="viz-whisker" x1="${x1}" y1="${midY - 5}" x2="${x1}" y2="${midY + 5}"/>`;
      }
      g += `<circle class="viz-dot-open" cx="${vizXd(row.predicted, axis).toFixed(1)}" cy="${midY}" r="5"/>`;
      if (row.achieved != null){
        g += `<circle class="viz-dot" cx="${vizXd(row.achieved, axis).toFixed(1)}" cy="${midY}" r="5.5"/>`;
      }
      const right = row.achieved != null
        ? `${row.diffText}${row.brText && row.brText !== '-' ? ` · ${row.brText}` : ''}`
        : 'awaiting achieved score';
      g += `<text class="viz-row-fact" x="${COL_CHANGE.rightX}" y="${midY + 4.5}">${escapeHtml(right)}</text>`;
      svg += `<g class="viz-row">${g}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL_CHANGE.plotX}" y1="${axisY + 2}" x2="${COL_CHANGE.plotX + COL_CHANGE.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizXd(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
    }
    return `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="Predicted versus achieved">${svg}</svg>`;
  }

  function vizPremorbidBlockHtml(){
    const topf = vizTopfPredictRows();
    const opie = vizOpiePredictRows();
    if (!topf.length && !opie.length) return '';
    const ciPct = preState.ciPct || premorbidCi().short;
    let cards = '';
    if (topf.length){
      cards += `<div class="viz-card"><div class="viz-card-name">ToPF-predicted WAIS-IV / WMS-IV</div>` +
        vizPrePanelSvg(topf) +
        `<div class="viz-card-note">Base rates are the percentage of the normative sample with a discrepancy at least this large. BASE_RATES is a parametric normal model, round(Φ(d / SEE)), not observed frequencies.</div></div>`;
    }
    if (opie.length){
      cards += `<div class="viz-card"><div class="viz-card-name">OPIE-4-predicted WAIS-IV</div>` +
        vizPrePanelSvg(opie) +
        `<div class="viz-card-note">OPIE-4 is illustrative only for UK use: the coefficients reproduce Holdnack et al. (2013) Table eA5.8 exactly, but the published equations also carry US education, ethnicity and region terms that are not applied, so every patient is scored at the US reference category. Its base rates are empirical.</div></div>`;
    }
    const legend = `<p class="viz-legend"><span class="viz-legend-open">○</span> predicted, with its ${escapeHtml(ciPct)} prediction interval shaded · <span class="viz-legend-filled">●</span> achieved · thin line marks the population mean of 100. Enter achieved scores on the premorbid page; rows without one show the prediction alone.</p>`;
    return `<div class="viz-section-head">Premorbid — predicted vs achieved</div>${legend}${cards}`;
  }

  /* ================== SD INDEX ========================================= */
  function vizSdiBlockHtml(){
    const rows = sdiRows.filter(r => r.name && !r.isExample);
    if (!rows.length) return '';
    const cv = parseFloat(document.getElementById('sdi-cv')?.value) || 0.95;
    const crit = cv === 0.90 ? 1.645 : cv === 0.95 ? 1.96 : cv;
    const axis = { min: -4, max: 4, ticks: [-4, -3, -2, -1, 0, 1, 2, 3, 4] };
    const chartable = rows.filter(r => sdiComputeChange(r) !== null);
    if (!chartable.length){
      return `<div class="viz-section-head">Standard Deviation Index</div>` +
        `<div class="viz-empty"><p>Rows are entered on the SD Index page, but none can compute yet (both scores, and an SD in raw mode). Complete them there.</p></div>`;
    }

    const families = new Map();
    for (const r of chartable){
      const g = r.group || '';
      if (!families.has(g)) families.set(g, []);
      families.get(g).push(r);
    }
    let cards = '';
    for (const [g, fRows] of families){
      const H = HEADER_H + fRows.length * ROW_H + AXIS_H;
      let svg = `<text class="viz-col-head" x="${COL_CHANGE.scoreMid}" y="13" text-anchor="middle">Scores</text>` +
                `<text class="viz-col-head" x="${COL_CHANGE.rightX}" y="13">SD Δ · Significance</text>`;
      fRows.forEach((r, i) => {
        const change = sdiComputeChange(r);
        const rowY = HEADER_H + i * ROW_H;
        const midY = rowY + ROW_H / 2;
        const sig = sdiCvHit(change, cv);
        const label = sig ? 'Significant change' : 'No significant change';
        let row = `<title>${escapeHtml(r.name)}: ${escapeHtml(String(r.t1))} → ${escapeHtml(String(r.t2))}, SD Δ ${fmt(change, 2)} — ${label}</title>`;
        row += `<text class="viz-row-name" x="${COL_CHANGE.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name, 24))}</text>`;
        row += `<text class="viz-row-score" x="${COL_CHANGE.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${escapeHtml(String(r.t1))} → ${escapeHtml(String(r.t2))}</text>`;
        const bx0 = vizXd(-crit, axis), bx1 = vizXd(crit, axis);
        row += `<rect class="viz-rci-band" x="${bx0.toFixed(1)}" y="${rowY + 3}" width="${(bx1 - bx0).toFixed(1)}" height="${ROW_H - 6}"/>`;
        const zx = vizXd(0, axis).toFixed(1);
        row += `<line class="viz-rci-expected" x1="${zx}" y1="${rowY + 3}" x2="${zx}" y2="${rowY + ROW_H - 3}"/>`;
        row += vizDotSvg(change, axis, vizXd, midY);
        row += `<text class="viz-row-fact" x="${COL_CHANGE.rightX}" y="${midY + 4.5}">${fmt(change, 2)} · ${label}</text>`;
        svg += `<g class="viz-row">${row}</g>`;
      });
      const axisY = HEADER_H + fRows.length * ROW_H;
      svg += `<line class="viz-tick" x1="${COL_CHANGE.plotX}" y1="${axisY + 2}" x2="${COL_CHANGE.plotX + COL_CHANGE.plotW}" y2="${axisY + 2}"/>`;
      for (const tval of axis.ticks){
        const x = vizXd(tval, axis).toFixed(1);
        svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
               `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
      }
      const title = g ? stripAgeRange(g) : 'Ungrouped measures';
      cards += `<div class="viz-card"><div class="viz-card-name">${escapeHtml(title)}</div>` +
        `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="SD Index change">${svg}</svg></div>`;
    }
    const cvLabel = cv === 0.95 ? '95% (±1.96 SD)' : '90% (±1.645 SD)';
    const legend = `<p class="viz-legend">Dot: change between testings in SD units, on the SD Index page's own divisor per row. Shaded band: no significant change at ${cvLabel}, as set on the SD Index page.</p>`;
    return `<div class="viz-section-head">Standard Deviation Index</div>${legend}${cards}`;
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

    let html = '';

    if (rows.length){
      // One card per test family, in table order; ungrouped rows form one card.
      const families = new Map();
      for (const r of rows){
        const gKey = batteryGroupKeyOf(r);
        if (!families.has(gKey)) families.set(gKey, []);
        families.get(gKey).push(r);
      }
      html += `<div class="viz-section-head">Score Tables</div>`;
      html += `<div class="viz-method-tabs" role="group" aria-label="How to show the scores">` +
        VIZ_SCORE_VIEWS.map(([v, label]) =>
          `<button type="button" class="viz-method-tab${v === vizScoreView ? ' is-active' : ''}" data-viz-view="${v}">${label}</button>`
        ).join('') + `</div>`;
      /* The axis caption describes the whole block, so it is stated once
         here rather than repeated on every card. */
      if (VIZ_SCORE_AXIS_NOTE[vizScoreView]){
        html += `<p class="viz-legend">${VIZ_SCORE_AXIS_NOTE[vizScoreView]}</p>`;
      }
      for (const [gKey, fRows] of families){
        const title = gKey ? (batteryGroupLabel(gKey) || 'Untitled test') : 'Ungrouped measures';
        html += vizFamilyCard(title, fRows, cls, ciLevel);
      }
    }

    html += vizPremorbidBlockHtml();
    html += vizChangeBlockHtml();
    html += vizSdiBlockHtml();

    if (!html){
      grid.innerHTML = '';
      if (emptyEl) emptyEl.hidden = false;
      return;
    }
    if (emptyEl) emptyEl.hidden = true;
    grid.innerHTML = html;
  }

  /* The method tabs pick which method's results to draw. Delegated, because
     the grid is replaced wholesale on every render. */
  grid.addEventListener('click', e => {
    const method = e.target.closest('[data-viz-method]');
    if (method){ vizRciMethod = method.dataset.vizMethod; renderVizPage(); return; }
    const view = e.target.closest('[data-viz-view]');
    if (view){ vizScoreView = view.dataset.vizView; renderVizPage(); }
  });

  /* Re-render whenever the page becomes visible, so it always reflects the
     current Score Tables state. Same pattern as syncTopnav: observe the
     section's class attribute rather than hooking every navigation path. */
  new MutationObserver(() => {
    if (section.classList.contains('active')) renderVizPage();
  }).observe(section, { attributes: true, attributeFilter: ['class'] });

  if (section.classList.contains('active')) renderVizPage();
})();
