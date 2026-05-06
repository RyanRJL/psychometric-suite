/* =============================================================================
   PSYCHOMETRIC ASSISTANT — DESIGN SYSTEM JS (practice build)
   Small enhancements that turn the SPA into something that feels
   organisation-built rather than solo-coded:
     1. Dynamic <title> per page         (Stripe / Linear convention)
     2. Standardise scattered button microcopy
     3. Mark page as ready (kills FOUC during stylesheet load)
   ========================================================================== */

(function(){
  'use strict';

  /* --- 1. Dynamic page titles ----------------------------------------------
     Listens for nav-item clicks and updates the document title to:
       "<Page name> — Psychometric Assistant"
     Falls back gracefully if the page label can't be resolved. */
  const SITE = 'Psychometric Assistant';
  const TITLE_MAP = {
    'home':           'Home',
    'converter':      'Score Converter',
    'battery':        'Neuropsych Report Tables',
    'report-writer':  'Report Writer',
    'effectsize':     'Effect Sizes',
    'change-analysis':'Change Analysis',
    'rci-basic':      'Standard Deviation Index',
    'rci-practice':   'Simple Reliable Change',
    'rci-srb':        'Practice-Adjusted',
    'rci-mcsweeney':  'McSweeney Regression-Based',
    'rci-crawford':   'Crawford Regression-Based',
    'premorbid':      'Premorbid Estimation',
    'about':          'Methods & References',
    'custom-tests':   'Norms Database'
  };

  function setTitleForTarget(target){
    const label = TITLE_MAP[target] || target;
    document.title = `${label} — ${SITE}`;
  }

  function syncTitleFromActiveSection(){
    const active = document.querySelector('section.section.active');
    if (active && active.id) setTitleForTarget(active.id);
  }

  document.addEventListener('click', e => {
    const navBtn = e.target.closest('[data-target]');
    if (!navBtn) return;
    const target = navBtn.dataset.target;
    if (target) setTitleForTarget(target);
  }, true);

  /* Run once on load and on hash change so reloads land on the right title */
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', syncTitleFromActiveSection);
  } else {
    syncTitleFromActiveSection();
  }
  window.addEventListener('hashchange', syncTitleFromActiveSection);

  /* --- 2. Microcopy standardisation ----------------------------------------
     Swap the small handful of inconsistent button labels. Sentence case,
     no period. Cheap to maintain — small explicit map. */
  const COPY_REPLACEMENTS = [
    [/^\s*\+\s*Add another subtest\s*$/i, '+ Add row'],
    [/^\s*\+\s*Add test row\s*$/i,        '+ Add row'],
    [/^\s*Clear all rows\s*$/i,           'Clear all'],
    [/^\s*\+\s*Load example row\s*$/i,    'Load example']
  ];

  function applyMicrocopy(root){
    if (!root) return;
    const candidates = root.querySelectorAll('button, .btn, .add-row-btn');
    candidates.forEach(el => {
      // Only operate on buttons whose visible text matches the patterns.
      // Don't touch elements with rich child structure (icons + text).
      const txt = el.textContent.trim();
      for (const [re, rep] of COPY_REPLACEMENTS){
        if (re.test(txt)){
          // Preserve any leading icon by replacing only the trailing text node
          // when present; otherwise just set textContent.
          const lastText = Array.from(el.childNodes).reverse().find(n => n.nodeType === 3);
          if (lastText) lastText.textContent = ' ' + rep.replace(/^\+\s*/, '');
          else el.textContent = rep;
          break;
        }
      }
    });
  }

  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', () => applyMicrocopy(document));
  } else {
    applyMicrocopy(document);
  }

  /* --- 3. Ready flag --------------------------------------------------------
     Adds a class to <html> once everything has loaded. Useful hook for
     fade-in transitions or for selectively suppressing animations until
     the page has settled. */
  window.addEventListener('load', () => {
    document.documentElement.classList.add('ds-ready');
  });

  /* --- 4b. Inline control bar — Change Analysis pages ----------------------
     For each Change Analysis sub-method (sdi, rci-basic, rci-practice,
     rci-srb, rci-crawford), build an inline control bar that sits directly
     above the data table and contains:
       LEFT  — quick-add autofill search (over Subtest + norms columns)
       RIGHT — 90/95% segmented toggle (over t(RB)/p result columns)
     The bar is just a UI shell — it programmatically updates the original
     hidden controls so all existing app logic keeps working. */
  const CHANGE_METHODS = [
    { id: 'sdi',           hasThreshold: true },
    { id: 'rci-basic',     hasThreshold: true },
    { id: 'rci-practice',  hasThreshold: true },
    { id: 'rci-srb',       hasThreshold: true },
    { id: 'rci-crawford',  hasThreshold: true }
  ];

  /* Compact button labels for threshold options. Maps the underlying
     <option value> to a short visible label. SDI has 5 options
     (two confidence intervals + three SD thresholds); the others have 2. */
  const THRESHOLD_LABEL = {
    '0.90': '90% CI · 1.645',
    '0.95': '95% CI · 1.96',
    '1':    '1σ',
    '1.5':  '1.5σ',
    '2':    '2σ'
  };

  /* Score-type → group-row label. Used by updateScoresHeader() to
     re-label the "Scores" column-group cell when the dropdown changes. */
  const SCORE_LABEL = {
    standard: 'Standard scores',
    scaled:   'Scaled scores',
    t:        'T-scores',
    z:        'Z-scores',
    raw:      'Raw scores'
  };
  const SCORE_LABEL_VALUES = Object.values(SCORE_LABEL);
  function updateScoresHeader(methodId){
    const table = document.getElementById(methodId + '-table');
    if (!table) return;
    const scoreSelect = methodId === 'sdi'
      ? document.getElementById('sdi-type')
      : document.querySelector('.rci-score-type[data-target="' + methodId + '"]');
    if (!scoreSelect) return;
    const label = SCORE_LABEL[scoreSelect.value] || 'Scores';
    const groupRow = table.querySelector('thead .table-group-row');
    if (!groupRow) return;
    Array.from(groupRow.children).forEach(th => {
      const text = (th.textContent || '').trim();
      if (text === 'Scores' || SCORE_LABEL_VALUES.includes(text)){
        th.textContent = label;
      }
    });
  }
  function thresholdLabel(opt){
    return THRESHOLD_LABEL[opt.value] || opt.textContent.trim();
  }

  function buildInlineControlBar(method){
    const section = document.getElementById(method.id);
    if (!section) return;

    /* Hide legacy panels — runs every call so it works even after the bar
       is already built. Idempotent. The hide can fail on the first call
       (consolidation hasn't wrapped panels yet) and succeed on a retry. */
    const setupGrid = section.querySelector('.change-setup-grid');
    if (setupGrid) setupGrid.classList.add('ds-legacy-hidden');
    section.querySelectorAll(':scope > .panel').forEach(p => p.classList.add('ds-legacy-hidden'));

    const table = section.querySelector('table.input-table');
    if (!table) return;
    if (section.querySelector('.ds-inline-bar')) return;     /* already built */

    /* Locate the legacy controls we'll remote-control. Crawford & co use
       class-based selectors; SDI uses ID-based. Try both. */
    const sigSelect    = section.querySelector('.rci-cv, #sdi-cv');
    const scoreSelect  = section.querySelector('.rci-score-type, #sdi-type');
    const familyInput  = section.querySelector('.combo-input.rci-family-input, .combo-input');
    const familyList   = section.querySelector('.combo-list.rci-family-list, .combo-list');

    /* Build the bar — three sections: LEFT autofill, MIDDLE score type,
       RIGHT threshold toggle. Each section omitted if its source isn't found.
       The .combo class on the search wrapper is REQUIRED — there is a
       global document mousedown listener (app.js ~2134) that closes any
       open .combo-list when the click target's ancestor chain doesn't
       include .combo. Without this class, every checkbox click would
       collapse the dropdown before the change event registered. */
    const bar = document.createElement('div');
    bar.className = 'ds-inline-bar';
    bar.innerHTML = `
      <div class="ds-inline-bar-section ds-inline-bar-left">
        ${familyInput ? `
          <span class="ds-inline-bar-label">Quick add</span>
          <div class="ds-inline-bar-search combo">
            <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
              <circle cx="7" cy="7" r="4.5"/>
              <path d="M11 11l3 3"/>
            </svg>
            <input type="text" class="ds-inline-bar-input" placeholder="Search a test family (WAIS-IV, WMS-IV, RBANS…)" autocomplete="off">
          </div>
        ` : ''}
      </div>
      ${scoreSelect ? `
        <div class="ds-inline-bar-section ds-inline-bar-middle">
          <span class="ds-inline-bar-label">Score type</span>
          <select class="ds-inline-bar-select" aria-label="Input score type">
            ${[...scoreSelect.options].map(opt => `<option value="${opt.value}"${opt.value===scoreSelect.value?' selected':''}>${opt.textContent}</option>`).join('')}
          </select>
        </div>
      ` : ''}
      <div class="ds-inline-bar-section ds-inline-bar-right">
        ${method.hasThreshold && sigSelect ? `
          <span class="ds-inline-bar-label">Threshold</span>
          <div class="ds-inline-bar-toggle" role="radiogroup" aria-label="Significance threshold">
            ${[...sigSelect.options].map(opt => `
              <button type="button" class="ds-inline-bar-toggle-btn" data-cv="${opt.value}" role="radio" title="${opt.textContent.trim().replace(/"/g,'&quot;')}">${thresholdLabel(opt)}</button>
            `).join('')}
          </div>
        ` : ''}
      </div>
    `;
    table.parentNode.insertBefore(bar, table);

    /* Wire the autofill input → original (hidden) input */
    if (familyInput && familyList){
      const newInput = bar.querySelector('.ds-inline-bar-input');
      const dropdownHost = bar.querySelector('.ds-inline-bar-search');
      if (dropdownHost && familyList){
        dropdownHost.appendChild(familyList);   /* re-parent the dropdown */
      }
      if (newInput){
        const sync = (v) => {
          familyInput.value = v;
          familyInput.dispatchEvent(new Event('input',  { bubbles:true }));
          familyInput.dispatchEvent(new Event('change', { bubbles:true }));
        };
        newInput.addEventListener('focus', () => {
          familyInput.dispatchEvent(new Event('focus', { bubbles:true }));
          familyList.classList.add('show');
        });
        newInput.addEventListener('blur',  () => setTimeout(() => {
          if (!familyList.matches(':hover')) familyList.classList.remove('show');
        }, 180));
        newInput.addEventListener('input', () => sync(newInput.value));
        newInput.addEventListener('keydown', e => {
          familyInput.dispatchEvent(new KeyboardEvent('keydown', {
            key: e.key, bubbles: true
          }));
          if (e.key === 'Escape') familyList.classList.remove('show');
        });
      }
    }

    /* Wire score-type dropdown → original (hidden) <select>.
       Also sync the table's group-row "Scores" header to reflect the
       chosen score type each time it changes. */
    if (scoreSelect){
      const newSelect = bar.querySelector('.ds-inline-bar-select');
      if (newSelect){
        newSelect.addEventListener('change', () => {
          scoreSelect.value = newSelect.value;
          scoreSelect.dispatchEvent(new Event('change', { bubbles:true }));
          updateScoresHeader(method.id);
        });
        scoreSelect.addEventListener('change', () => {
          if (newSelect.value !== scoreSelect.value) newSelect.value = scoreSelect.value;
          updateScoresHeader(method.id);
        });
      }
      /* Set the initial header text once on bar build */
      updateScoresHeader(method.id);
    }

    /* Wire the segmented threshold toggle → original (hidden) <select> */
    if (sigSelect && method.hasThreshold){
      const buttons = bar.querySelectorAll('.ds-inline-bar-toggle-btn');
      const refresh = () => {
        const v = sigSelect.value;
        buttons.forEach(b => {
          const active = b.dataset.cv === v;
          b.classList.toggle('is-active', active);
          b.setAttribute('aria-checked', active ? 'true' : 'false');
        });
      };
      buttons.forEach(b => {
        b.addEventListener('click', () => {
          sigSelect.value = b.dataset.cv;
          sigSelect.dispatchEvent(new Event('change', { bubbles:true }));
          refresh();
        });
      });
      sigSelect.addEventListener('change', refresh);
      refresh();
    }

    /* (legacy panel hiding already done at top of function) */
  }

  function buildAllChangeBars(){
    CHANGE_METHODS.forEach(buildInlineControlBar);
    buildBatteryInlineBar();
  }

  /* --- Battery / Neuropsych Tables inline bar -----------------------------
     Same pattern as Change Analysis but four sections:
       LEFT          Quick add search
       MIDDLE-1      Classification (Wechsler / AACN — segmented toggle)
       MIDDLE-2      Score Type (Scaled / Standard / T / Z — dropdown)
       RIGHT         Premorbid (checkbox; score input appears when ticked)
     Hides the original Configuration + Premorbid Comparison side panels
     and the legacy bat-table-head autofill row. */
  function buildBatteryInlineBar(){
    const section = document.getElementById('battery');
    if (!section) return;

    /* Hide legacy chrome (always — idempotent) */
    section.querySelectorAll('.bat-config-panel, .bat-premorbid-block, .bat-side-panels, .bat-table-head').forEach(el => {
      el.classList.add('ds-legacy-hidden');
    });

    if (section.querySelector('.ds-inline-bar-battery')) return;          /* already built */

    const familyInput   = document.getElementById('bat-family-input');
    const familyList    = document.getElementById('bat-family-list');
    const classSelect   = document.getElementById('bat-class');
    const typeSelect    = document.getElementById('bat-type');
    const premCheckbox  = document.getElementById('bat-prem-enable');
    const premScoreInput= document.getElementById('bat-prem-score');
    const table         = document.getElementById('bat-table');
    if (!table) return;

    const bar = document.createElement('div');
    bar.className = 'ds-inline-bar ds-inline-bar-battery';
    bar.innerHTML = `
      <div class="ds-inline-bar-section ds-inline-bar-left">
        ${familyInput ? `
          <span class="ds-inline-bar-label">Quick add</span>
          <div class="ds-inline-bar-search combo">
            <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
              <circle cx="7" cy="7" r="4.5"/>
              <path d="M11 11l3 3"/>
            </svg>
            <input type="text" class="ds-inline-bar-input" placeholder="Search a test family…" autocomplete="off">
          </div>
        ` : ''}
      </div>
      ${typeSelect ? `
        <div class="ds-inline-bar-section ds-inline-bar-middle">
          <span class="ds-inline-bar-label">Score type</span>
          <select class="ds-inline-bar-select" aria-label="Score type">
            ${[...typeSelect.options].map(opt => `<option value="${opt.value}"${opt.value===typeSelect.value?' selected':''}>${opt.textContent}</option>`).join('')}
          </select>
        </div>
      ` : ''}
      ${classSelect ? `
        <div class="ds-inline-bar-section ds-inline-bar-right">
          <span class="ds-inline-bar-label">Class</span>
          <div class="ds-inline-bar-toggle" role="radiogroup" aria-label="Classification system">
            ${[...classSelect.options].map(opt => `
              <button type="button" class="ds-inline-bar-toggle-btn" data-class="${opt.value}" role="radio">${opt.textContent.trim()}</button>
            `).join('')}
          </div>
        </div>
      ` : ''}
    `;

    /* Premorbid card — sits ABOVE the inline bar. Contains a header eyebrow,
       description, checkbox to enable, and a score input that appears when
       enabled. More descriptive than a single bar slot allowed. */
    const premCard = document.createElement('div');
    premCard.className = 'ds-prem-card' + (premCheckbox && premCheckbox.checked ? ' is-enabled' : '');
    if (premCheckbox && premScoreInput){
      premCard.innerHTML = `
        <span class="ds-inline-bar-label ds-prem-card-eyebrow">Premorbid comparison</span>
        <div class="ds-prem-card-row">
          <label class="ds-prem-card-toggle">
            <input type="checkbox" class="ds-prem-card-checkbox"${premCheckbox.checked ? ' checked' : ''}>
            <span class="ds-prem-card-mark" aria-hidden="true">
              <svg viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="2.5,6.5 5,9 9.5,3.5"/></svg>
            </span>
            <span class="ds-prem-card-desc">
              <span class="ds-prem-card-desc-line">Flag subtest scores meaningfully below an estimated premorbid ability.</span>
              <span class="ds-prem-card-desc-line">Adds <strong>*</strong> when ≥1 SD below, <strong>**</strong> when ≥1.5 SD, <strong>***</strong> when ≥2 SD.</span>
            </span>
          </label>
          <div class="ds-prem-card-input-wrap">
            <label class="ds-prem-card-input-label" for="ds-prem-card-score">Premorbid score</label>
            <input type="number" id="ds-prem-card-score" class="ds-prem-card-score" placeholder="e.g. 110" step="any"${premCheckbox.checked ? '' : ' disabled'} aria-label="Premorbid standard score">
          </div>
        </div>
      `;
    }
    /* Insert: premorbid card first (above), then inline bar, then table */
    const tableBlock = section.querySelector('.bat-table-block') || table;
    if (premCheckbox && premScoreInput) tableBlock.parentNode.insertBefore(premCard, tableBlock);
    tableBlock.parentNode.insertBefore(bar, tableBlock);

    /* Wire Quick add — same pattern as Change Analysis */
    if (familyInput && familyList){
      const newInput = bar.querySelector('.ds-inline-bar-input');
      const dropdownHost = bar.querySelector('.ds-inline-bar-search');
      if (dropdownHost && familyList) dropdownHost.appendChild(familyList);
      if (newInput){
        newInput.addEventListener('focus', () => {
          familyInput.dispatchEvent(new Event('focus', { bubbles:true }));
          familyList.classList.add('show');
        });
        newInput.addEventListener('blur', () => setTimeout(() => {
          if (!familyList.matches(':hover')) familyList.classList.remove('show');
        }, 180));
        newInput.addEventListener('input', () => {
          familyInput.value = newInput.value;
          familyInput.dispatchEvent(new Event('input',  { bubbles:true }));
          familyInput.dispatchEvent(new Event('change', { bubbles:true }));
        });
        newInput.addEventListener('keydown', e => {
          familyInput.dispatchEvent(new KeyboardEvent('keydown', { key: e.key, bubbles: true }));
          if (e.key === 'Escape') familyList.classList.remove('show');
        });
      }
    }

    /* Wire Classification toggle — `data-class` attribute is unique to the
       class buttons (threshold buttons use `data-cv`), so a bare attribute
       selector is the safest match regardless of where the section lives. */
    if (classSelect){
      const buttons = bar.querySelectorAll('.ds-inline-bar-toggle-btn[data-class]');
      const refresh = () => {
        const v = classSelect.value;
        buttons.forEach(b => {
          const active = b.dataset.class === v;
          b.classList.toggle('is-active', active);
          b.setAttribute('aria-checked', active ? 'true' : 'false');
        });
      };
      buttons.forEach(b => {
        b.addEventListener('click', () => {
          classSelect.value = b.dataset.class;
          classSelect.dispatchEvent(new Event('change', { bubbles:true }));
          refresh();
        });
      });
      classSelect.addEventListener('change', refresh);
      refresh();
    }

    /* Wire Score Type dropdown */
    if (typeSelect){
      const newSelect = bar.querySelector('.ds-inline-bar-middle .ds-inline-bar-select');
      if (newSelect){
        newSelect.addEventListener('change', () => {
          typeSelect.value = newSelect.value;
          typeSelect.dispatchEvent(new Event('change', { bubbles:true }));
        });
        typeSelect.addEventListener('change', () => {
          if (newSelect.value !== typeSelect.value) newSelect.value = typeSelect.value;
        });
      }
    }

    /* Wire Premorbid CARD checkbox + score input — card lives above the bar.
       Checkbox controls whether the score input is enabled and the card
       picks up its "is-enabled" treatment (lighter on rest, terracotta-tinted
       when active). */
    if (premCheckbox && premScoreInput){
      const newCheckbox = premCard.querySelector('.ds-prem-card-checkbox');
      const newScore    = premCard.querySelector('.ds-prem-card-score');
      if (newCheckbox){
        const syncFromOriginal = () => {
          newCheckbox.checked = premCheckbox.checked;
          if (newScore){
            newScore.disabled = !premCheckbox.checked;
            if (newScore.value !== premScoreInput.value) newScore.value = premScoreInput.value;
          }
          premCard.classList.toggle('is-enabled', premCheckbox.checked);
        };
        newCheckbox.addEventListener('change', () => {
          premCheckbox.checked = newCheckbox.checked;
          premCheckbox.dispatchEvent(new Event('change', { bubbles:true }));
          syncFromOriginal();
          if (newCheckbox.checked && newScore) setTimeout(() => newScore.focus(), 50);
        });
        premCheckbox.addEventListener('change', syncFromOriginal);
        syncFromOriginal();
      }
      if (newScore){
        newScore.addEventListener('input', () => {
          premScoreInput.value = newScore.value;
          premScoreInput.dispatchEvent(new Event('input', { bubbles:true }));
          premScoreInput.dispatchEvent(new Event('change', { bubbles:true }));
        });
        premScoreInput.addEventListener('input', () => {
          if (newScore.value !== premScoreInput.value) newScore.value = premScoreInput.value;
        });
      }
    }
  }

  /* Use a MutationObserver on the change-analysis container so the bar
     gets built whenever the SPA finishes consolidating its panels, AND
     gets re-built if the user switches between sub-methods that may
     re-render. Belt-and-braces: also fire on a few timeouts for the case
     where the observer attaches AFTER the relevant mutations fire. */
  function startChangeObserver(){
    /* Run immediately in case panels are already in their final state */
    buildAllChangeBars();

    const watch = document.body;
    if (!watch || typeof MutationObserver === 'undefined'){
      /* Fallback: timer-based retries */
      [100, 300, 600, 1200, 2400].forEach(ms => setTimeout(buildAllChangeBars, ms));
      return;
    }

    let scheduled = false;
    const observer = new MutationObserver(() => {
      if (scheduled) return;
      scheduled = true;
      requestAnimationFrame(() => {
        scheduled = false;
        buildAllChangeBars();
      });
    });
    observer.observe(watch, { childList: true, subtree: true });

    /* Belt-and-braces timer fires for cases the observer might miss */
    [100, 300, 600, 1200, 2400, 5000].forEach(ms => setTimeout(buildAllChangeBars, ms));
  }

  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', startChangeObserver);
  } else {
    startChangeObserver();
  }

  /* --- 4e. Tag Simple RCI norm cells -------------------------------------
     The Simple Reliable Change table (#rci-basic-table) is skipped by the
     legacy setupNormsLockToggle, so its SD and r columns never get the
     data-norm-cell="true" attribute. Tag them ourselves so the borderless
     CSS rules in section 34 apply consistently. */
  function tagSimpleRciNorms(){
    const table = document.getElementById('rci-basic-table');
    if (!table || !table.tBodies[0]) return;
    /* Columns: 0=#, 1=Subtest, 2=SD, 3=r, 4=Date1, 5=Date2, 6=RCI, 7=p, 8=Outcome, 9=actions
       The two norm cells are columns 2 and 3 (zero-indexed). */
    const NORM_COLS = [2, 3];
    /* Tag header row */
    table.querySelectorAll('thead tr').forEach(row => {
      Array.from(row.children).forEach((cell, i) => {
        if (NORM_COLS.includes(i)) cell.dataset.normCell = 'true';
      });
    });
    /* Tag body rows */
    table.tBodies[0].querySelectorAll('tr').forEach(row => {
      Array.from(row.children).forEach((cell, i) => {
        if (NORM_COLS.includes(i)) cell.dataset.normCell = 'true';
      });
    });
  }
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', tagSimpleRciNorms);
  } else {
    tagSimpleRciNorms();
  }
  /* Re-tag whenever rows are added/removed */
  if (typeof MutationObserver !== 'undefined'){
    const t = document.getElementById('rci-basic-table');
    if (t && t.tBodies[0]){
      new MutationObserver(tagSimpleRciNorms).observe(t.tBodies[0], { childList: true, subtree: true });
    } else {
      /* Try again after the SPA renders */
      [200, 600, 1200, 2400].forEach(ms => setTimeout(tagSimpleRciNorms, ms));
    }
  }

  /* --- 4d. Working Report chip — green orbit on nav-tab click --------------
     Visual feedback: a green light sweeps around the chip's circumference
     each time the user clicks a top-nav tab. One-shot 1.4s animation,
     re-triggerable by removing/re-adding the class with a forced reflow. */
  document.addEventListener('click', e => {
    const navTab = e.target.closest('.topnav-item, .topnav-drop-item');
    if (!navTab) return;
    const chip = document.querySelector('.rb-chip');
    if (!chip) return;
    chip.classList.remove('rb-chip-orbit');
    void chip.offsetWidth;                         /* force reflow to restart animation */
    chip.classList.add('rb-chip-orbit');
    setTimeout(() => chip.classList.remove('rb-chip-orbit'), 2500);  /* matches CSS 2.4s + tiny buffer */
  }, true);

  /* --- 4c. Remove the deprecated "Lock norms" feature ----------------------
     The original app injects a Lock/Unlock pill into each Change Analysis
     table's "Norms" group cell, defaults the table to "locked" (norms
     faded + inputs disabled), and toggles via the pill. We've removed this
     feature, so:
       1. Strip the .norms-locked class from every input-table
       2. Re-enable any norm-cell inputs that were `disabled`
       3. Remove the Lock/Unlock pill if it's already in the DOM
     CSS in section 33 hides any visual remnants. */
  function clearNormsLock(){
    document.querySelectorAll('.input-table.norms-locked').forEach(t => {
      t.classList.remove('norms-locked');
    });
    document.querySelectorAll('.input-table [data-norm-cell="true"] input[disabled]').forEach(inp => {
      inp.disabled = false;
      inp.removeAttribute('disabled');
    });
    document.querySelectorAll('.norms-toggle-pill').forEach(p => p.remove());
  }
  /* Run on DOM ready and again on every mutation, so newly added rows /
     re-rendered tables don't reapply the lock. */
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', clearNormsLock);
  } else {
    clearNormsLock();
  }
  if (typeof MutationObserver !== 'undefined'){
    let nlScheduled = false;
    new MutationObserver(() => {
      if (nlScheduled) return;
      nlScheduled = true;
      requestAnimationFrame(() => {
        nlScheduled = false;
        clearNormsLock();
      });
    }).observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ['class', 'disabled'] });
  }
  /* Belt-and-braces timer fires too */
  [200, 600, 1200].forEach(ms => setTimeout(clearNormsLock, ms));

  /* --- 4. Chip auto-fade on scroll ------------------------------------------
     While the user is actively scrolling (i.e. reading / working), the
     Working Report chip recedes (opacity + slight scale-down) so it isn't
     a hot spot on top of content. When scrolling stops for ~600ms or the
     user scrolls back up, it returns to full presence. The drawer being
     open suspends this entirely — chip stays fully visible. */
  let scrollTimer = null;
  let lastScrollY = window.scrollY;
  const RECEDE_CLASS = 'ds-chip-recede';
  const RECEDE_DELAY = 600;

  function onScroll(){
    const drawerOpen = document.querySelector('.rb-root.is-open');
    if (drawerOpen){
      document.body.classList.remove(RECEDE_CLASS);
      return;
    }
    const dy = window.scrollY - lastScrollY;
    lastScrollY = window.scrollY;
    // Only recede on downward scroll — upward scroll restores immediately
    if (dy > 2) document.body.classList.add(RECEDE_CLASS);
    else if (dy < -2) document.body.classList.remove(RECEDE_CLASS);
    // Reset timer — when scrolling stops, restore the chip
    clearTimeout(scrollTimer);
    scrollTimer = setTimeout(() => {
      document.body.classList.remove(RECEDE_CLASS);
    }, RECEDE_DELAY);
  }
  // Listen on the scrollable element used by the SPA — `<main>` in this app
  function attachChipFade(){
    const mainEl = document.querySelector('.main');
    const target = mainEl || window;
    target.addEventListener('scroll', onScroll, { passive: true });
  }
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', attachChipFade);
  } else {
    attachChipFade();
  }

  // ============================================================================
  // PREMORBID — 2-column layout (sticky inputs LEFT, tabbed tables RIGHT)
  //
  // Restructures the existing markup at runtime so this is a pure additive
  // overlay (no HTML edits required, easy to revert by deleting this block).
  // Also builds a minimalistic model-availability tick list inside the aside
  // that updates live as the user fills in fields.
  // ============================================================================

  function restructurePremorbid(){
    const sec = document.getElementById('premorbid');
    if (!sec || sec.classList.contains('ds-prem-restructured')) return;

    const tabs        = sec.querySelector('.pre-tabs');
    const inputsTab   = sec.querySelector('#pre-inputs');
    const inputsPanel = inputsTab && inputsTab.querySelector('.premorbid-input-panel');
    const estimates   = sec.querySelector('#pre-estimates');
    const predict     = sec.querySelector('#pre-predict');
    const opiePredict = sec.querySelector('#pre-opiepredict');
    const tabNav      = sec.querySelector('.pre-tab-nav');

    if (!tabs || !inputsPanel || !estimates) return;

    const layout = document.createElement('div');
    layout.className = 'ds-prem-layout';
    const aside = document.createElement('aside');
    aside.className = 'ds-prem-aside';
    const main  = document.createElement('div');
    main.className = 'ds-prem-main';
    layout.appendChild(aside);
    layout.appendChild(main);

    tabs.parentNode.insertBefore(layout, tabs);

    aside.appendChild(inputsPanel);

    // Shorten field labels for the compact aside layout. Map by input id so
    // we don't depend on the source label text. Idempotent.
    const labelMap = {
      'pre-topf':  'ToPF Raw',
      'pre-vc':    'Vocabulary',
      'pre-mr':    'Matrix Reasoning',
      'pre-sex':   'Sex',
      'pre-occ':   'Occupation',
      'pre-edu':   'Education (yrs)',
      'pre-age':   'Age',
      'pre-ci':    'Confidence',
      'pre-title': 'Title'
    };
    Object.keys(labelMap).forEach(id => {
      const ctrl = aside.querySelector('#' + id);
      if (!ctrl) return;
      const wrap = ctrl.closest('.field');
      if (!wrap) return;
      const lbl = wrap.querySelector('label');
      if (!lbl) return;
      const sub = lbl.querySelector('.pre-label-sub');
      lbl.textContent = labelMap[id];
      if (sub) lbl.appendChild(sub);
    });
    if (inputsTab) inputsTab.style.display = 'none';

    main.appendChild(tabs);
    if (estimates)   main.appendChild(estimates);
    if (predict)     main.appendChild(predict);
    if (opiePredict) main.appendChild(opiePredict);
    if (tabNav)      main.appendChild(tabNav);

    const inputsTabBtn = tabs.querySelector('[data-pre-tab="inputs"]');
    const estimatesBtn = tabs.querySelector('[data-pre-tab="estimates"]');
    if (inputsTabBtn && inputsTabBtn.classList.contains('active') && estimatesBtn){
      estimatesBtn.click();
    }
    if (inputsTabBtn && inputsTabBtn.parentNode){
      inputsTabBtn.parentNode.removeChild(inputsTabBtn);
    }

    const ticks = document.createElement('div');
    ticks.className = 'ds-prem-ticks';
    ticks.id = 'ds-prem-ticks';
    ticks.setAttribute('aria-live', 'polite');
    ticks.innerHTML =
      '<div class="ds-prem-ticks-head">'
      + '<div class="ds-prem-ticks-eyebrow">Predictions</div>'
      + '<div class="ds-prem-ticks-count" id="ds-prem-ticks-count">0 / 4</div>'
      + '</div>'
      + '<ul class="ds-prem-ticks-list">'
      +   '<li class="ds-prem-ticks-row" data-model="topf-raw">'
      +     '<span class="ds-prem-ticks-mark" aria-hidden="true"></span>'
      +     '<span class="ds-prem-ticks-name">ToPF raw score</span>'
      +   '</li>'
      +   '<li class="ds-prem-ticks-row" data-model="topf-demo">'
      +     '<span class="ds-prem-ticks-mark" aria-hidden="true"></span>'
      +     '<span class="ds-prem-ticks-name">ToPF + Demographics</span>'
      +   '</li>'
      +   '<li class="ds-prem-ticks-row" data-model="crawford">'
      +     '<span class="ds-prem-ticks-mark" aria-hidden="true"></span>'
      +     '<span class="ds-prem-ticks-name">Crawford &amp; Allan</span>'
      +   '</li>'
      +   '<li class="ds-prem-ticks-row" data-model="opie4">'
      +     '<span class="ds-prem-ticks-mark" aria-hidden="true"></span>'
      +     '<span class="ds-prem-ticks-name">OPIE-4 prorated FSIQ</span>'
      +   '</li>'
      + '</ul>';
    aside.appendChild(ticks);

    sec.classList.add('ds-prem-restructured');
  }

  function premHas(id){
    const el = document.getElementById(id);
    if (!el) return false;
    return (el.value || '').trim() !== '';
  }
  function updatePremorbidTicks(){
    const list = document.getElementById('ds-prem-ticks');
    if (!list) return;
    const has = {
      topf: premHas('pre-topf'), vc:  premHas('pre-vc'),  mr:  premHas('pre-mr'),
      sex:  premHas('pre-sex'),  occ: premHas('pre-occ'), edu: premHas('pre-edu'),
      age:  premHas('pre-age')
    };
    const ready = {
      'topf-raw' : has.topf,
      'topf-demo': has.topf && has.sex && has.occ && has.edu,
      'crawford' : has.sex && has.occ && has.edu && has.age,
      'opie4'    : has.age && (has.vc || has.mr)
    };
    let n = 0;
    list.querySelectorAll('.ds-prem-ticks-row').forEach(row => {
      const k = row.getAttribute('data-model');
      const r = !!ready[k];
      row.classList.toggle('is-ready', r);
      if (r) n++;
    });
    const counter = document.getElementById('ds-prem-ticks-count');
    if (counter){
      counter.textContent = n + ' / 4';
      counter.classList.toggle('is-all-ready', n === 4);
    }
  }
  function bindPremorbidTicks(){
    const ids = ['pre-topf','pre-vc','pre-mr','pre-sex','pre-occ','pre-edu','pre-age'];
    ids.forEach(id => {
      const el = document.getElementById(id);
      if (!el || el.dataset.dsPremBound === '1') return;
      el.dataset.dsPremBound = '1';
      el.addEventListener('input',  updatePremorbidTicks);
      el.addEventListener('change', updatePremorbidTicks);
    });
    updatePremorbidTicks();
  }

  function initPremorbidLayout(){
    restructurePremorbid();
    bindPremorbidTicks();
  }

  // ============================================================================
  // COMBO DROPDOWNS — clamp height to available viewport space
  // The legacy CSS used a fixed 330–520px max-height. On smaller viewports
  // (or when the input sits low on the page) this clipped the bottom of the
  // dropdown below the viewport, forcing the user to scroll the WINDOW.
  // Solution: when a combo-list opens, measure its top position and set its
  // max-height to (viewport − top − margin). The inner .combo-options gets a
  // matching cap so its overflow-y still triggers.
  // ============================================================================
  function clampComboList(list){
    if (!list || !list.classList.contains('show')) return;
    const rect = list.getBoundingClientRect();
    const margin = 16;
    const available = Math.max(220, window.innerHeight - rect.top - margin);
    // setProperty with 'important' so we beat the CSS `!important` cap
    list.style.setProperty('max-height', available + 'px', 'important');
    const opts = list.querySelector('.combo-options');
    if (opts){
      const footer = list.querySelector('.combo-footer');
      const note   = list.querySelector('.combo-ageband-note');
      const reserve = (footer ? footer.offsetHeight : 0)
                    + (note   ? note.offsetHeight   : 0)
                    + 12;
      const optMax = Math.max(160, available - reserve);
      opts.style.setProperty('max-height', optMax + 'px', 'important');
    }
  }
  // Watch every combo-list for class changes (the legacy code toggles `show`)
  function observeComboLists(){
    document.querySelectorAll('.combo-list').forEach(list => {
      if (list.dataset.dsComboObs === '1') return;
      list.dataset.dsComboObs = '1';
      const obs = new MutationObserver(() => {
        if (list.classList.contains('show')) clampComboList(list);
      });
      obs.observe(list, { attributes: true, attributeFilter: ['class'] });
    });
  }
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', observeComboLists);
  } else {
    observeComboLists();
  }
  // Re-bind when new dropdowns are added later (SPA / dynamic re-render)
  setTimeout(observeComboLists, 500);
  setTimeout(observeComboLists, 1500);
  // Re-clamp on window resize (open dropdown follows the new viewport)
  window.addEventListener('resize', () => {
    document.querySelectorAll('.combo-list.show').forEach(clampComboList);
  });

  // Topnav dropdown shortcuts: items with data-pre-tab="..." should activate
  // the matching tab inside the Premorbid section after navigation.
  document.addEventListener('click', e => {
    const item = e.target.closest('[data-pre-tab][data-target="premorbid"]');
    if (!item) return;
    const want = item.dataset.preTab;
    // Defer so the SPA's section-switch logic has run.
    setTimeout(() => {
      const tabBtn = document.querySelector(
        '#premorbid .pre-tabs [data-pre-tab="' + want + '"]'
      );
      if (tabBtn) tabBtn.click();
    }, 0);
  }, true);
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', initPremorbidLayout);
  } else {
    initPremorbidLayout();
  }

})();
