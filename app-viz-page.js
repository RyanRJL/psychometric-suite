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
   standard score, or the raw scores themselves), Premorbid
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

  /* ---------- geometry: author the chart at its REAL pixel width --------
     The SVG used to be authored at a fixed 940 units and scaled by the
     browser to whatever the stage was - 635px here, a factor of 0.68. Two
     things followed, both bad: every label rendered at 68% of its intended
     size, and the chart's HEIGHT became a by-product of its width rather
     than something that could fill the space.

     Authoring at the measured width makes the scale 1:1, so a 13px label is
     13px and a 28px row really is 28px. Columns are fractions of that width
     rather than fixed offsets, and the plot band gets the larger share it
     was always short of (51% before, 56% now). */
  const COL_FRAC        = { nameEnd: 0.20,  scoreMid: 0.235, plotX: 0.27, plotW: 0.56, rightX: 0.845 };
  const COL_FRAC_CHANGE = { nameEnd: 0.175, scoreMid: 0.225, plotX: 0.28, plotW: 0.55, rightX: 0.845 };
  const SVG_W_MIN = 480;

  let SVG_W = 940;
  let COL        = { nameEnd: 205, scoreMid: 228, plotX: 262, plotW: 478, rightX: 754 };
  /* Change rows print a PAIR of scores ("92→85") where a score card prints
     one, so they get their own geometry - the name ends earlier and the
     score column sits wider. Sharing COL made long names collide with the
     pair. */
  let COL_CHANGE = { nameEnd: 186, scoreMid: 232, plotX: 286, plotW: 454, rightX: 754 };

  /* Stage width is derived from the container, not read off the stage, so
     the geometry is known BEFORE the first render rather than after it. */
  function vizSetGeometry(){
    const railW = 248, gap = 16, cardPad = 34;   // must track the CSS
    const avail = (grid.clientWidth || 900) - railW - gap - cardPad;
    SVG_W = Math.max(SVG_W_MIN, Math.round(avail));
    const px = f => Math.round(f * SVG_W);
    COL        = { nameEnd: px(COL_FRAC.nameEnd), scoreMid: px(COL_FRAC.scoreMid),
                   plotX: px(COL_FRAC.plotX), plotW: px(COL_FRAC.plotW), rightX: px(COL_FRAC.rightX) };
    COL_CHANGE = { nameEnd: px(COL_FRAC_CHANGE.nameEnd), scoreMid: px(COL_FRAC_CHANGE.scoreMid),
                   plotX: px(COL_FRAC_CHANGE.plotX), plotW: px(COL_FRAC_CHANGE.plotW), rightX: px(COL_FRAC_CHANGE.rightX) };
  }
  const HEADER_H = 20, AXIS_H = 24;
  /* Row PITCH is variable, and this is what makes a tall chart fit.
     Scaling the whole SVG down (max-height + letterbox) shrinks the TEXT
     with it, which is why the first attempt refused to shrink far enough to
     be useful. The rendered scale is set by WIDTH - card width / SVG_W - so
     reducing the row pitch in user units makes the drawing shorter while
     every label stays exactly the same size on screen. */
  const ROW_H_DEFAULT = 28, ROW_H_MIN = 17, ROW_H_MAX = 46;
  let ROW_H = ROW_H_DEFAULT;

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
  /* The name column is now a real pixel width, so the budget scales with it
     (~0.52 of the 13px font's em width per character). */
  function vizTruncate(name, max){
    const n = max || Math.max(12, Math.floor(COL.nameEnd / 6.8));
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

  /* ================== EXPORT: copy to clipboard / save as PNG ===========
     An SVG rendered through <img> gets NO stylesheet - the document's CSS
     does not cross that boundary - so every mark would rasterise as black
     fill on transparent. The clone therefore carries its computed styles
     inline, and two CSS features that SVG-in-<img> cannot honour are
     resolved here instead:

       - text-transform, which is a CSS property with no SVG equivalent, so
         the text content itself is upper-cased in the clone;
       - webfonts, which an <img>-embedded SVG may not load at all, so an
         explicit family stack is set and the raster falls back to a local
         face. Numbers and labels stay legible; the glyphs may differ from
         the screen.

     A background rectangle is painted first: a transparent PNG dropped into
     a report reads as a hole wherever the document is not white. */
  const VIZ_EXPORT_SCALE = 2;   // 2x for a crisp figure in a printed report
  const VIZ_INLINE_PROPS = [
    'fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-linecap',
    'stroke-linejoin', 'stroke-dasharray', 'opacity', 'font-size',
    'font-family', 'font-weight', 'letter-spacing', 'text-anchor'
  ];

  function vizSurfaceColour(){
    const probe = getComputedStyle(document.body).getPropertyValue('--ds-surface');
    return (probe || '').trim() || '#FDFAF6';
  }

  function vizSvgToString(svg, title){
    const clone = svg.cloneNode(true);
    const src = svg.querySelectorAll('*');
    const dst = clone.querySelectorAll('*');
    for (let i = 0; i < src.length; i++){
      const cs = getComputedStyle(src[i]);
      let decl = '';
      for (const prop of VIZ_INLINE_PROPS){
        const v = cs.getPropertyValue(prop);
        if (v) decl += `${prop}:${v};`;
      }
      dst[i].setAttribute('style', decl);
      if (cs.getPropertyValue('text-transform') === 'uppercase' && dst[i].tagName === 'text'){
        dst[i].textContent = (dst[i].textContent || '').toUpperCase();
      }
    }

    const vb = (svg.getAttribute('viewBox') || '').split(/\s+/).map(Number);
    const w = vb[2] || svg.clientWidth || 940;
    const h = vb[3] || 200;
    /* The exported figure names itself. A chart pasted into a report with no
       test name on it cannot be identified once it leaves the app. */
    const titleH = title ? 26 : 0;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('width', w);
    clone.setAttribute('height', h + titleH);
    clone.setAttribute('viewBox', `0 ${-titleH} ${w} ${h + titleH}`);
    if (title){
      const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      t.setAttribute('x', '0');
      t.setAttribute('y', String(-titleH + 14));
      t.setAttribute('style', `font-family:${VIZ_FONT_STACK};font-size:14px;font-weight:600;fill:#1A1410;`);
      t.textContent = title;
      clone.insertBefore(t, clone.firstChild);
    }
    const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    bg.setAttribute('x', '0');
    bg.setAttribute('y', String(-titleH));
    bg.setAttribute('width', String(w));
    bg.setAttribute('height', String(h + titleH));
    bg.setAttribute('fill', vizSurfaceColour());
    clone.insertBefore(bg, clone.firstChild);

    return { xml: new XMLSerializer().serializeToString(clone), w, h: h + titleH };
  }
  const VIZ_FONT_STACK = "Inter, system-ui, -apple-system, 'Segoe UI', sans-serif";

  function vizSvgToPngBlob(svg, title){
    return new Promise((resolve, reject) => {
      const { xml, w, h } = vizSvgToString(svg, title);
      const url = URL.createObjectURL(new Blob([xml], { type: 'image/svg+xml;charset=utf-8' }));
      const img = new Image();
      img.onload = () => {
        try {
          const canvas = document.createElement('canvas');
          canvas.width = Math.round(w * VIZ_EXPORT_SCALE);
          canvas.height = Math.round(h * VIZ_EXPORT_SCALE);
          const ctx = canvas.getContext('2d');
          ctx.fillStyle = vizSurfaceColour();
          ctx.fillRect(0, 0, canvas.width, canvas.height);
          ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
          URL.revokeObjectURL(url);
          canvas.toBlob(b => b ? resolve(b) : reject(new Error('PNG encoding failed')), 'image/png');
        } catch (e){ URL.revokeObjectURL(url); reject(e); }
      };
      img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('SVG could not be rasterised')); };
      img.src = url;
    });
  }

  function vizFileName(title){
    return (title || 'chart').replace(/[^\w\s-]/g, '').trim().replace(/\s+/g, '-').toLowerCase() + '.png';
  }

  async function vizCopyChart(svg, title){
    try {
      const blob = await vizSvgToPngBlob(svg, title);
      /* Clipboard image write needs a secure context and the ClipboardItem
         API. Opening the bundle straight off disk qualifies in current
         Chrome and Firefox, but an older browser will not have it - so the
         failure is reported with the alternative rather than silently
         doing nothing. */
      if (!navigator.clipboard || typeof ClipboardItem === 'undefined'){
        throw new Error('no-clipboard-api');
      }
      await navigator.clipboard.write([new ClipboardItem({ 'image/png': blob })]);
      if (typeof showToast === 'function') showToast('✓ Chart copied to clipboard');
    } catch (e){
      if (typeof showToast === 'function'){
        showToast(e && e.message === 'no-clipboard-api'
          ? 'This browser cannot copy images — use Save PNG instead'
          : 'Could not copy the chart — use Save PNG instead', true);
      }
    }
  }

  async function vizSaveChart(svg, title){
    try {
      const blob = await vizSvgToPngBlob(svg, title);
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = vizFileName(title);
      document.body.appendChild(a);
      a.click();
      a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
      if (typeof showToast === 'function') showToast('✓ Chart saved as PNG');
    } catch (e){
      if (typeof showToast === 'function') showToast('Could not save the chart as PNG', true);
    }
  }

  const VIZ_ICON_COPY = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"></rect><path d="M5 15H4a2 2 0 01-2-2V4a2 2 0 012-2h9a2 2 0 012 2v1"></path></svg>';
  const VIZ_ICON_SAVE = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" aria-hidden="true"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>';

  /* ---------- page furniture, borrowed from the rest of the app --------
     .panel is the app's card and is global, so it is used as-is rather
     than re-invented. (The in-page h2.block-title heading is no longer
     needed: with one source pane on screen at a time, the source tabs
     name what is being shown.) The tab and control classes are standalone because the
     app's own tab components are page-scoped and would arrive unstyled
     outside their section - see the note at the top of the CSS block. */
  function vizTab(attr, value, label, active){
    return `<button type="button" class="viz-tab" role="tab" aria-selected="${active ? 'true' : 'false'}" ${attr}="${escapeAttr(value)}">${escapeHtml(label)}</button>`;
  }
  /* Controls and the display toggle share one panel. They were two stacked
     blocks, which cost ~100px of height on a page whose whole point is
     fitting in one window. */
  function vizControls(kicker, tabs, note, extra){
    return `<div class="panel">
      <div class="viz-controls-body">
        <div class="viz-controls-main">
          <div class="panel-kicker">${escapeHtml(kicker)}</div>
          <div class="viz-tabs" role="tablist" aria-label="${escapeAttr(kicker)}">${tabs}</div>
        </div>
        ${extra || ''}
      </div>
      ${note ? `<p class="viz-axis-note">${note}</p>` : ''}
    </div>`;
  }
  function vizDisplayToggle(){
    return `<div class="viz-controls-main viz-controls-right">
      <div class="panel-kicker">Show</div>
      <div class="viz-tabs" role="tablist" aria-label="How many charts to show">
        ${vizTab('data-viz-display', 'single', 'One at a time', vizSingle)}
        ${vizTab('data-viz-display', 'grid', 'All charts', !vizSingle)}
      </div>
    </div>`;
  }

  /* Card shell: the app's own .panel, a head shaped like .apa-toolbar
     (micro-label left, actions right), and the app's .btn for the actions. */
  function vizPanel(title, body, note){
    return `<div class="panel" data-chart-title="${escapeAttr(title)}">
      <div class="viz-card-head">
        <span class="viz-card-title">${escapeHtml(title)}</span>
        <div class="viz-card-actions">
          <button type="button" class="btn btn-icon" data-viz-copy title="Copy this chart to the clipboard as a PNG">${VIZ_ICON_COPY}Copy</button>
          <button type="button" class="btn btn-icon" data-viz-save title="Save this chart as a PNG file">${VIZ_ICON_SAVE}PNG</button>
        </div>
      </div>
      ${body}
      ${note ? `<div class="viz-card-note">${note}</div>` : ''}
    </div>`;
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
    ['raw',        'Raw scores']
  ];
  let vizScoreView = 'native';
  const VIZ_SCORE_AXIS_NOTE = {
    native: '',
    percentile: 'Axis is the percentile rank, converted from each row’s own metric so the whole test shares one scale. A percentile scale is deliberately non-linear: it compresses the tails, so two clearly different low scores can sit close together, and confidence intervals become asymmetric. Scores are shown as entered.',
    standard: 'Axis is the standard-score metric (M 100, SD 15), converted from each row’s own metric so the whole test shares one scale. Scores are shown as entered.',
    raw: 'The raw scores as entered, so their spread is visible directly. Nothing is derived from them — no norms, no classification bands, no percentile. The axis spans the values entered for that test, not each measure’s possible range, because no raw maximum is held for any measure; where a test’s measures are counted on different scales their positions are not comparable with each other.'
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

  /* ---------- raw view: the raw scores themselves ----------------------
     Deliberately NOT a normative comparison. This view exists so the
     clinician can see the raw numbers and how they sit relative to each
     other - the spread of what was actually recorded - with nothing
     derived from them: no norms, no bands, no percentile, no
     classification.

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

  function vizRawSpreadPanelSvg(rows){
    const vals = rows.map(vizRawValueOf);
    /* Axis spans the values actually entered. It cannot span each measure's
       full possible range, because the app holds no raw maximum for any
       measure - so the axis is the observed spread, not "out of N". */
    const axis = vizNiceAxis(Math.min(...vals), Math.max(...vals));
    const H = HEADER_H + rows.length * ROW_H + AXIS_H;
    let svg = `<text class="viz-col-head" x="${COL.scoreMid}" y="13" text-anchor="middle">Raw</text>` +
              `<text class="viz-col-head" x="${COL.rightX}" y="13">Score as entered</text>`;

    rows.forEach((r, i) => {
      const v = vals[i];
      const rowY = HEADER_H + i * ROW_H;
      const midY = rowY + ROW_H / 2;
      const type = rowScoreType(r);
      /* The score column of a raw-metric row already holds this number, so
         repeating it on the right would say the same thing twice. */
      const alsoScored = type !== 'raw' && r.score !== '' && !isNaN(parseFloat(r.score));

      let row = `<title>${escapeHtml(r.name)}: raw ${v}${alsoScored ? `, ${scoreTypeLabel(type).toLowerCase()} ${r.score}` : ''}</title>`;
      row += `<text class="viz-row-name" x="${COL.nameEnd}" y="${midY + 4.5}" text-anchor="end">${escapeHtml(vizTruncate(r.name))}</text>`;
      row += `<text class="viz-row-score" x="${COL.scoreMid}" y="${midY + 4.5}" text-anchor="middle">${v}</text>`;
      // A hairline across the plot gives each dot a baseline to sit on, so
      // the eye reads the row's position without a band implying a range.
      row += `<line class="viz-raw-rule" x1="${COL.plotX}" y1="${midY}" x2="${COL.plotX + COL.plotW}" y2="${midY}"/>`;
      row += vizDotSvg(v, axis, vizX, midY);
      row += `<text class="viz-row-fact" x="${COL.rightX}" y="${midY + 4.5}">${alsoScored ? escapeHtml(`${r.score} (${scoreTypeAbbr(type)})`) : '–'}</text>`;
      svg += `<g class="viz-row">${row}</g>`;
    });

    const axisY = HEADER_H + rows.length * ROW_H;
    svg += `<line class="viz-tick" x1="${COL.plotX}" y1="${axisY + 2}" x2="${COL.plotX + COL.plotW}" y2="${axisY + 2}"/>`;
    for (const tval of axis.ticks){
      const x = vizX(tval, axis).toFixed(1);
      svg += `<line class="viz-tick" x1="${x}" y1="${axisY + 2}" x2="${x}" y2="${axisY + 6}"/>` +
             `<text class="viz-tick-label" x="${x}" y="${axisY + 18}" text-anchor="middle">${tval}</text>`;
    }
    return `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="Raw scores as entered">${svg}</svg>`;
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
      /* Every row that HAS a raw score is plotted, whatever metric it is
         reported on - the point of this view is to see the recorded
         numbers, so nothing is filtered on normative grounds. */
      const plot = rows.filter(r => vizRawValueOf(r) !== null);
      const listed = rows.filter(r => vizRawValueOf(r) === null);
      if (plot.length) body += vizRawSpreadPanelSvg(plot);
      /* One shared line rather than a paragraph per row - repeating the
         same sentence eight times reads as an error list, not a note. */
      if (listed.length){
        body += `<div class="viz-raw-row">No raw score entered on Score Tables: ` +
          listed.map(r => `<b>${escapeHtml(r.name)}</b>`).join(', ') + `.</div>`;
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
        /* One line, not one box per measure: four dashed boxes repeating the
           same sentence cost ~160px on RBANS and read as an error list. Each
           measure keeps its score and interval; the reason is stated once. */
        body += `<div class="viz-raw-row">Not plotted — no standardised metric, so no percentile or classification is derived: ` +
          rawRows.map(r => {
            const ciText = ciLevel !== 'off' ? getBatteryCiHtml(parseFloat(r.score), r, ciLevel) : '';
            return `<b>${escapeHtml(r.name)}</b> ${escapeHtml(String(r.score))}${ciText ? ` (${ciText})` : ''}`;
          }).join(', ') + `.</div>`;
      }
      if (brRows.length) notes.push('Step lines show the published percentage of the normative sample obtaining each span or higher (WAIS-IV Manual, Tables C.4–C.5). A HIGH base rate means a common, and therefore lower, score.');
    }

    if (anyReversed) notes.unshift('↓ error measure: a higher score means more errors. Percentiles are as obtained; the classification bands describe performance, so that row’s bands run reversed.');

    return vizPanel(title, body, notes.length ? notes.join('<br>') : '');
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

    const ctrl = { kicker:'Method', tabs: VIZ_RCI_METHODS.map(([m, label]) =>
      vizTab('data-viz-method', m, label, m === method)).join(''), note:'' };

    const cards = [];
    let empty = '';
    if (!chartable.length){
      empty = `<p class="viz-empty-line">Rows are entered on Change Analysis, but none has the values the ${escapeHtml(VIZ_RCI_METHODS.find(x => x[0] === method)[1])} method needs yet — scores at both dates, plus the normative parameters.</p>`;
    } else {
      const families = new Map();
      for (const r of chartable){
        const g = r.group || '';
        if (!families.has(g)) families.set(g, []);
        families.get(g).push(r);
      }
      for (const [g, fRows] of families){
        const title = g ? caGroupDisplay(g) : 'Ungrouped measures';
        cards.push({ title, html: vizPanel(title, vizRciPanelSvg(fRows, method), '') });
      }
    }
    const cvLabel = st.cv === 0.95 ? '95%' : '90%';
    const legend = `<p class="viz-legend"><span class="viz-legend-open">○</span> ${escapeHtml(st.d1 || 'Date 1')} · <span class="viz-legend-filled">●</span> ${escapeHtml(st.d2 || 'Date 2')} · shaded band: no reliable change at the ${cvLabel} level set on that method's page. Outcomes state significance only, not direction.</p>`;
    return { key:'change', label:'Change Analysis', ctrl, legend, cards, empty };
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
    const cards = [];
    if (topf.length){
      cards.push({ title:'ToPF-predicted WAIS-IV / WMS-IV', html: vizPanel('ToPF-predicted WAIS-IV / WMS-IV', vizPrePanelSvg(topf),
        'Base rates are the percentage of the normative sample with a discrepancy at least this large. BASE_RATES is a parametric normal model, round(Φ(d / SEE)), not observed frequencies.') });
    }
    if (opie.length){
      cards.push({ title:'OPIE-4-predicted WAIS-IV', html: vizPanel('OPIE-4-predicted WAIS-IV', vizPrePanelSvg(opie),
        'OPIE-4 is illustrative only for UK use: the coefficients reproduce Holdnack et al. (2013) Table eA5.8 exactly, but the published equations also carry US education, ethnicity and region terms that are not applied, so every patient is scored at the US reference category. Its base rates are empirical.') });
    }
    const legend = `<p class="viz-legend"><span class="viz-legend-open">○</span> predicted, with its ${escapeHtml(ciPct)} prediction interval shaded · <span class="viz-legend-filled">●</span> achieved · thin line marks the population mean of 100. Enter achieved scores on the premorbid page; rows without one show the prediction alone.</p>`;
    return { key:'premorbid', label:'Premorbid', ctrl:null, legend, cards, empty:'' };
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
      return { key:'sdi', label:'SD Index', ctrl:null, legend:'', cards:[],
        empty:'<p class="viz-empty-line">Rows are entered on the SD Index page, but none can compute yet — each needs both scores, and an SD in raw mode.</p>' };
    }

    const families = new Map();
    for (const r of chartable){
      const g = r.group || '';
      if (!families.has(g)) families.set(g, []);
      families.get(g).push(r);
    }
    const cards = [];
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
      cards.push({ title, html: vizPanel(title, `<svg class="viz-svg" viewBox="0 0 ${SVG_W} ${H}" role="img" aria-label="SD Index change">${svg}</svg>`, '') });
    }
    const cvLabel = cv === 0.95 ? '95% (±1.96 SD)' : '90% (±1.645 SD)';
    const legend = `<p class="viz-legend">Dot: change between testings in SD units, on the SD Index page's own divisor per row. Shaded band: no significant change at ${cvLabel}, as set on the SD Index page.</p>`;
    return { key:'sdi', label:'SD Index', ctrl:null, legend, cards, empty:'' };
  }

  /* ---------- page render --------------------------------------------- */
  /* What the table settings are, stated where they are useful rather than in
     a heading: these are Score Tables' settings, not this page's, and the
     page has no way to change them. */
  function vizSettingsText(cls, ciLevel){
    const prem = getBatteryPremorbidComparison();
    return [
      cls === 'wechsler' ? 'Wechsler classification' : 'Guilmette et al. (2020)',
      ciLevel === 'off' ? 'no confidence intervals' : `${ciLevel}% confidence intervals`,
      prem ? 'premorbid comparison on' : 'premorbid comparison off'
    ].join(' · ');
  }

  function vizScoreTablesBlock(cls, ciLevel){
    /* Same membership rule as the APA export: named rows with a numeric
       score, example rows excluded - a chart of a seeded example could be
       mistaken for patient data. */
    const rows = batteryRows.filter(r =>
      r.name && !r.isExample && r.score !== '' && !isNaN(parseFloat(r.score)));
    if (!rows.length) return null;

    // One card per test family, in table order; ungrouped rows form one card.
    const families = new Map();
    for (const r of rows){
      const gKey = batteryGroupKeyOf(r);
      if (!families.has(gKey)) families.set(gKey, []);
      families.get(gKey).push(r);
    }
    /* The axis caption describes the whole block, so it is stated once in
       the controls panel rather than repeated on every card. */
    const ctrl = { kicker:'Show scores as',
      tabs: VIZ_SCORE_VIEWS.map(([v, label]) => vizTab('data-viz-view', v, label, v === vizScoreView)).join(''),
      note: VIZ_SCORE_AXIS_NOTE[vizScoreView] };
    const cards = [];
    for (const [gKey, fRows] of families){
      const title = gKey ? (batteryGroupLabel(gKey) || 'Untitled test') : 'Ungrouped measures';
      cards.push({ title, html: vizFamilyCard(title, fRows, cls, ciLevel) });
    }
    return { key:'battery', label:'Score Tables', ctrl, legend:'', cards, empty:'' };
  }

  /* ---------- one source at a time -------------------------------------
     The rest of the app does not present a long scrolling page: premorbid
     shows one pane of four with a Back/Next bar and a "Step 2 of 4"
     counter, and Change Analysis is split across five nav entries. This
     page now follows that, with the four sources as panes.

     Within a pane the cards can be shown all at once (the profile across a
     battery is the comparison a clinician actually reads) or one at a time
     for a close look. Both are kept because they answer different
     questions; neither is a substitute for the other. */
  let vizSource = null;
  /* Single is the DEFAULT. The rest of this site does not present a long
     scrolling page, and "all charts at once" is the deliberate overview
     rather than the resting state. */
  let vizSingle = true;
  const vizCardIndex = {};

  /* ---------- shell: control rail + chart stage ------------------------
     Two columns, as on Score Tables (.bat-panels-grid) and premorbid (its
     always-visible Inputs aside): everything you SET lives in a rail on the
     left, and the chart occupies the stage on the right.

     The page itself does not scroll. The layout is given an explicit height
     and each column scrolls INSIDE it, so the header, the controls and the
     chart's own navigation stay put while you move between charts - which
     is the behaviour that makes it read as an app rather than a document. */
  function vizRailHtml(blocks, active, settingsText){
    let out = '<aside class="viz-rail">';

    out += '<div class="viz-rail-group"><div class="panel-kicker">Source</div><div class="viz-rail-list" role="tablist" aria-label="Chart source">';
    for (const b of blocks){
      out += `<button type="button" class="viz-rail-item" role="tab" aria-selected="${b.key === active.key}" data-viz-source="${b.key}">` +
             `<span class="viz-rail-label">${escapeHtml(b.label)}</span>` +
             `<span class="viz-rail-count">${b.cards.length || ''}</span></button>`;
    }
    out += '</div></div>';

    if (active.ctrl && active.ctrl.tabs){
      out += `<div class="viz-rail-group"><div class="panel-kicker">${escapeHtml(active.ctrl.kicker)}</div>` +
             `<div class="viz-rail-list" role="tablist" aria-label="${escapeAttr(active.ctrl.kicker)}">${active.ctrl.tabs}</div></div>`;
    }

    if (active.cards.length > 1){
      out += '<div class="viz-rail-group"><div class="panel-kicker">Charts</div><div class="viz-rail-list" role="tablist" aria-label="Charts">';
      const idx = vizCardIndex[active.key] || 0;
      active.cards.forEach((c, i) => {
        out += `<button type="button" class="viz-rail-item viz-rail-chart" role="tab" aria-selected="${!vizSingle ? 'false' : i === idx}" data-viz-goto="${i}">` +
               `<span class="viz-rail-num">${i + 1}</span><span class="viz-rail-label">${escapeHtml(c.title)}</span></button>`;
      });
      out += '</div></div>';
      out += `<div class="viz-rail-group"><div class="panel-kicker">Show</div>
        <div class="viz-rail-list viz-rail-row" role="tablist" aria-label="How many charts to show">
          ${vizTab('data-viz-display', 'single', 'One', vizSingle)}
          ${vizTab('data-viz-display', 'grid', 'All', !vizSingle)}
        </div></div>`;
    }

    if (active.ctrl && active.ctrl.note){
      out += `<p class="viz-axis-note">${active.ctrl.note}</p>`;
    }
    /* The legend used to float above the stage, orphaned from the chart it
       described; the rail had 59% of its height empty. */
    if (active.legend){
      out += `<div class="viz-rail-group viz-rail-foot">${active.legend}</div>`;
    }
    if (settingsText){
      out += `<div class="viz-rail-group viz-rail-foot">
        <div class="panel-kicker">Table settings</div>
        <p class="viz-axis-note">${escapeHtml(settingsText)}. Change these on the Score Tables page.</p>
      </div>`;
    }
    return out + '</aside>';
  }

  function vizStageHtml(block){
    const n = block.cards.length;
    if (block.empty) return `<div class="viz-stage">${block.empty}</div>`;

    let out = '<div class="viz-stage">';

    if (vizSingle && n >= 1){
      const idx = Math.min(Math.max(vizCardIndex[block.key] || 0, 0), n - 1);
      vizCardIndex[block.key] = idx;
      out += `<div class="viz-cards viz-cards-single">${block.cards[idx].html}</div>`;
      if (n > 1){
        out += `<div class="viz-gallery-nav">
          <button type="button" class="btn" data-viz-card="prev"${idx === 0 ? ' disabled' : ''}>← Previous</button>
          <div class="viz-step-count">${idx + 1} of ${n}</div>
          <button type="button" class="btn" data-viz-card="next"${idx === n - 1 ? ' disabled' : ''}>Next →</button>
        </div>`;
      }
    } else {
      out += `<div class="viz-cards">${block.cards.map(c => c.html).join('')}</div>`;
    }
    return out + '</div>';
  }

  /* ---------- size the layout to the window ---------------------------
     The page is given an explicit height and its two columns scroll inside
     it, so the document itself never scrolls. Then the chart's ROW PITCH is
     tightened until the stage needs no scrollbar either.

     Pitch, not scale, is the lever that matters: the rendered scale is set
     by WIDTH (stage width / SVG_W), so shortening the row pitch in user
     units makes the drawing shorter while every label stays exactly the
     same size on screen. Letterboxing via max-height would shrink the text
     with it, which is why it is only the last resort.

     Trap worth keeping: `zoom:0.9` on <body> means getBoundingClientRect
     returns VISUAL pixels while CSS lengths and scrollHeight are LAYOUT
     pixels, and mixing them sizes everything ~11% wrong. */
  let vizFitPass = 0;

  function vizFitToWindow(){
    const layout = grid.querySelector('.viz-layout');
    const stage = grid.querySelector('.viz-stage');
    if (!layout || !stage) return;

    const zoom = parseFloat(getComputedStyle(document.body).zoom) || 1;
    const topVisual = layout.getBoundingClientRect().top;
    /* Room between the top of the layout and the bottom of the window, less
       the footer strip that sits below it. */
    const footer = document.querySelector('.site-footer, footer');
    const footerH = footer ? footer.getBoundingClientRect().height : 0;
    let availLayout = (window.innerHeight - topVisual - footerH - 16) / zoom;
    if (availLayout < 200) return;
    layout.style.height = Math.floor(availLayout) + 'px';

    /* Whatever still sits below the layout - the disclosure block, margins,
       the footer's own spacing - is measured rather than itemised: set the
       height, see what the document has left over, take that off. An
       enumerated version of this was wrong twice, once by 128px flat. */
    const spill = document.body.scrollHeight - document.body.clientHeight;
    if (spill > 0){
      availLayout = Math.max(200, availLayout - spill);
      layout.style.height = Math.floor(availLayout) + 'px';
    }

    if (!vizSingle){ grid.style.removeProperty('--viz-fit-h'); return; }
    const svg = stage.querySelector('.viz-cards-single .viz-svg');
    if (!svg) return;

    const rect = svg.getBoundingClientRect();
    const overflow = stage.scrollHeight - stage.clientHeight;
    const svgLayoutH = rect.height / zoom;

    /* SPARE ROOM, measured from the stage's actual content rather than from
       scrollHeight. When the stage does not overflow, scrollHeight equals
       clientHeight, so a scrollHeight-based figure can only ever say "it
       fits" - which is why the pitch shrank when it had to but never grew,
       leaving a five-row chart in a fifth of the stage. Summing the
       children's heights measures the free space in both directions. */
    const gap = parseFloat(getComputedStyle(stage).rowGap) || 0;
    const kids = [...stage.children];
    const contentLayoutH = kids.reduce((a, el) => a + el.getBoundingClientRect().height / zoom, 0)
                         + gap * Math.max(0, kids.length - 1);
    const free = stage.clientHeight - contentLayoutH;
    const availForSvg = svgLayoutH + free - 6;

    const vb = (svg.getAttribute('viewBox') || '').split(/\s+/).map(Number);
    const rows = Math.max(1, Math.round(((vb[3] || 0) - HEADER_H - AXIS_H) / ROW_H));
    const unitsPerLayoutPx = SVG_W / (rect.width / zoom);
    const availUnits = availForSvg * unitsPerLayoutPx;

    let want = Math.floor((availUnits - HEADER_H - AXIS_H) / rows);
    /* Grows as well as shrinks: a five-row chart used to sit in a quarter
       of the stage because the pitch could only ever come down. */
    want = Math.max(ROW_H_MIN, Math.min(ROW_H_MAX, want));

    if (want !== ROW_H && vizFitPass < 2){
      ROW_H = want; vizFitPass++;
      renderVizPage();
      return;
    }
    /* Pitch is as tight as it goes; letterbox whatever is left rather than
       leave the stage scrolling. */
    if (overflow > 0 && availForSvg > 90){
      grid.style.setProperty('--viz-fit-h', Math.floor(availForSvg) + 'px');
    } else {
      grid.style.removeProperty('--viz-fit-h');
    }
  }

  function renderVizPage(){
    vizSetGeometry();
    const cls = document.getElementById('bat-class')?.value || 'wechsler';
    const ciLevel = document.getElementById('bat-ci-level')?.value || 'off';
    const settingsText = vizSettingsText(cls, ciLevel);

    const blocks = [
      vizScoreTablesBlock(cls, ciLevel),
      vizPremorbidBlockHtml(),
      vizChangeBlockHtml(),
      vizSdiBlockHtml()
    ].filter(b => b && (b.cards.length || b.empty));

    if (!blocks.length){
      grid.innerHTML = '';
      if (emptyEl) emptyEl.hidden = false;
      return;
    }
    if (emptyEl) emptyEl.hidden = true;

    /* A source can disappear between renders - clearing the SDI table, say -
       so a remembered selection has to be re-checked, not trusted. */
    let active = blocks.find(b => b.key === vizSource) || blocks[0];
    vizSource = active.key;

    grid.innerHTML = `<div class="viz-layout">${vizRailHtml(blocks, active, settingsText)}${vizStageHtml(active)}</div>`;
    requestAnimationFrame(vizFitToWindow);
  }

  /* The available height changes with the window, not just with a render. */
  window.addEventListener('resize', () => {
    if (section.classList.contains('active')){ ROW_H = ROW_H_DEFAULT; vizFitPass = 0; renderVizPage(); }
  });

  /* One delegated handler for the whole page, because the content is
     replaced wholesale on every render. */
  grid.addEventListener('click', e => {
    /* Any deliberate change restarts the fit from the natural pitch. */
    if (e.target.closest('[data-viz-source],[data-viz-display],[data-viz-card],[data-viz-method],[data-viz-view],[data-viz-goto]')){
      ROW_H = ROW_H_DEFAULT; vizFitPass = 0;
    }
    const src = e.target.closest('[data-viz-source]');
    if (src){ vizSource = src.dataset.vizSource; renderVizPage(); grid.scrollIntoView({block:'start', behavior:'auto'}); return; }
    const disp = e.target.closest('[data-viz-display]');
    if (disp){ vizSingle = disp.dataset.vizDisplay === 'single'; renderVizPage(); return; }
    const goto = e.target.closest('[data-viz-goto]');
    if (goto){ vizCardIndex[vizSource] = +goto.dataset.vizGoto; vizSingle = true; renderVizPage(); return; }
    const card = e.target.closest('[data-viz-card]');
    if (card){ vizStepCard(card.dataset.vizCard === 'next' ? 1 : -1); return; }
    const method = e.target.closest('[data-viz-method]');
    if (method){ vizRciMethod = method.dataset.vizMethod; renderVizPage(); return; }
    const view = e.target.closest('[data-viz-view]');
    if (view){ vizScoreView = view.dataset.vizView; renderVizPage(); return; }

    const copy = e.target.closest('[data-viz-copy]');
    const save = e.target.closest('[data-viz-save]');
    const btn = copy || save;
    if (!btn) return;
    const panel = btn.closest('.panel');
    const svg = panel && panel.querySelector('svg.viz-svg');
    if (!svg){
      if (typeof showToast === 'function') showToast('Nothing to export on this card', true);
      return;
    }
    const title = (panel.querySelector('.viz-card-title') || {}).textContent || 'Chart';
    (copy ? vizCopyChart : vizSaveChart)(svg, title);
  });

  function vizStepCard(delta){
    const cur = vizCardIndex[vizSource] || 0;
    vizCardIndex[vizSource] = Math.max(0, cur + delta);
    ROW_H = ROW_H_DEFAULT; vizFitPass = 0;   // the next chart re-fits from scratch
    renderVizPage();
  }

  /* Arrow keys page through the charts in single mode. Guarded so it never
     steals a keystroke from a field - every other page in this app is a
     data-entry page, and the tables have their own arrow handling. */
  document.addEventListener('keydown', e => {
    if (!section.classList.contains('active') || !vizSingle) return;
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    const t = document.activeElement;
    if (t && /^(INPUT|SELECT|TEXTAREA)$/.test(t.tagName)) return;
    if (t && t.isContentEditable) return;
    if (e.key === 'ArrowRight'){ e.preventDefault(); vizStepCard(1); }
    else if (e.key === 'ArrowLeft'){ e.preventDefault(); vizStepCard(-1); }
  });

  /* Re-render whenever the page becomes visible, so it always reflects the
     current Score Tables state. Same pattern as syncTopnav: observe the
     section's class attribute rather than hooking every navigation path. */
  new MutationObserver(() => {
    if (section.classList.contains('active')) renderVizPage();
  }).observe(section, { attributes: true, attributeFilter: ['class'] });

  if (section.classList.contains('active')) renderVizPage();
})();
