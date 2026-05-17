/* =============================================================================
   PSYCHOMETRIC ASSISTANT - DESIGN SYSTEM JS (practice build)
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
       "<Page name> - Psychometric Assistant"
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
    document.title = `${label} - ${SITE}`;
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
     no period. Cheap to maintain - small explicit map. */
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

  /* --- 4b. Inline control bar - Change Analysis pages ----------------------
     For each Change Analysis sub-method (sdi, rci-basic, rci-practice,
     rci-srb, rci-crawford), build an inline control bar that sits directly
     above the data table and contains:
       LEFT  - quick-add autofill search (over Subtest + norms columns)
       RIGHT - 90/95% segmented toggle (over t(RB)/p result columns)
     The bar is just a UI shell - it programmatically updates the original
     hidden controls so all existing app logic keeps working. */
  const CHANGE_METHODS = [
    { id: 'sdi',           hasThreshold: true },
    { id: 'rci-basic',     hasThreshold: true, hasCorrCorrToggle: true },
    { id: 'rci-practice',  hasThreshold: true, hasCorrCorrToggle: true },
    { id: 'rci-srb',       hasThreshold: true, hasCorrCorrToggle: true },
    { id: 'rci-crawford',  hasThreshold: true, hasCorrCorrToggle: true }
  ];

  /* Compact button labels for threshold options. Maps the underlying
     <option value> to a short visible label. SDI has 5 options
     (two confidence intervals + three SD thresholds); the others have 2. */
  const THRESHOLD_LABEL = {
    '0.90': '90% CI',
    '0.95': '95% CI',
    '1':    '1 SD',
    '1.5':  '1.5 SD',
    '2':    '2 SD'
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

    /* Hide legacy panels - runs every call so it works even after the bar
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
    const corrCheckbox = section.querySelector('.rci-use-corrected-r');

    /* Build the bar - three sections: LEFT autofill, MIDDLE score type,
       RIGHT threshold toggle. Each section omitted if its source isn't found.
       The .combo class on the search wrapper is REQUIRED - there is a
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
        ${method.hasCorrCorrToggle && corrCheckbox ? `
          <label class="ds-inline-bar-corr-toggle" title="Use the corrected (attenuation-adjusted) test–retest correlation when available. Falls back to the raw r when corrected r is not provided for a given test.">
            <input type="checkbox" class="ds-inline-bar-corr-input"${corrCheckbox.checked ? ' checked' : ''}>
            <span>Corrected <em>r</em></span>
          </label>
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

    /* Wire the corrected-r toggle → original (hidden) checkbox.
       Two-way sync so calling code can flip the legacy box and the visible
       toggle stays in lockstep. */
    if (corrCheckbox && method.hasCorrCorrToggle){
      const newCorr = bar.querySelector('.ds-inline-bar-corr-input');
      if (newCorr){
        newCorr.checked = !!corrCheckbox.checked;
        newCorr.addEventListener('change', () => {
          corrCheckbox.checked = newCorr.checked;
          corrCheckbox.dispatchEvent(new Event('change', { bubbles:true }));
        });
        corrCheckbox.addEventListener('change', () => {
          if (newCorr.checked !== corrCheckbox.checked) newCorr.checked = corrCheckbox.checked;
        });
      }
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
       MIDDLE-1      Classification (Wechsler / AACN - segmented toggle)
       MIDDLE-2      Score Type (Scaled / Standard / T / Z - dropdown)
       RIGHT         Premorbid (checkbox; score input appears when ticked)
     Hides the original Configuration + Premorbid Comparison side panels
     and the legacy bat-table-head autofill row. */
  function buildBatteryInlineBar(){
    const section = document.getElementById('battery');
    if (!section) return;

    /* Hide legacy chrome (always - idempotent) */
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

    /* Premorbid card - sits ABOVE the inline bar. Contains a header eyebrow,
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

    /* Wire Quick add - same pattern as Change Analysis */
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

    /* Wire Classification toggle - `data-class` attribute is unique to the
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

    /* Wire Premorbid CARD checkbox + score input - card lives above the bar.
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

  /* --- 4d. Working Report chip - green orbit on nav-tab click --------------
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
     open suspends this entirely - chip stays fully visible. */
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
    // Only recede on downward scroll - upward scroll restores immediately
    if (dy > 2) document.body.classList.add(RECEDE_CLASS);
    else if (dy < -2) document.body.classList.remove(RECEDE_CLASS);
    // Reset timer - when scrolling stops, restore the chip
    clearTimeout(scrollTimer);
    scrollTimer = setTimeout(() => {
      document.body.classList.remove(RECEDE_CLASS);
    }, RECEDE_DELAY);
  }
  // Listen on the scrollable element used by the SPA - `<main>` in this app
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
  // PREMORBID - 2-column layout (sticky inputs LEFT, tabbed tables RIGHT)
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

    /* Replace the Confidence Interval <select> with a segmented toggle, matching
       the threshold pills used on the Change Analysis bars. The original select
       stays in place (visually hidden) so all downstream JS that reads its
       value continues to work. */
    const ciSel = aside.querySelector('#pre-ci');
    if (ciSel && !aside.querySelector('.ds-prem-ci-toggle')){
      const wrap = ciSel.closest('.field');
      if (wrap){
        const toggle = document.createElement('div');
        toggle.className = 'ds-prem-ci-toggle ds-inline-bar-toggle';
        toggle.setAttribute('role', 'radiogroup');
        toggle.setAttribute('aria-label', 'Confidence interval');
        Array.from(ciSel.options).forEach(opt => {
          const btn = document.createElement('button');
          btn.type = 'button';
          btn.className = 'ds-inline-bar-toggle-btn';
          btn.setAttribute('role', 'radio');
          btn.dataset.val = opt.value;
          btn.textContent = (opt.value === '0.95') ? '95%' : '90%';
          btn.title = opt.textContent.trim();
          toggle.appendChild(btn);
        });
        // Insert toggle after the label, before the select
        ciSel.insertAdjacentElement('beforebegin', toggle);
        // Hide the original select but keep it functional for state
        ciSel.style.display = 'none';
        const refresh = () => {
          const v = ciSel.value;
          toggle.querySelectorAll('button').forEach(b => {
            const active = b.dataset.val === v;
            b.classList.toggle('is-active', active);
            b.setAttribute('aria-checked', active ? 'true' : 'false');
          });
        };
        toggle.querySelectorAll('button').forEach(b => {
          b.addEventListener('click', () => {
            ciSel.value = b.dataset.val;
            ciSel.dispatchEvent(new Event('change', { bubbles:true }));
            refresh();
          });
        });
        ciSel.addEventListener('change', refresh);
        refresh();
      }
    }

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
  // COMBO DROPDOWNS - clamp height to available viewport space
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

  // ============================================================================
  // REPORT WRITER · auto-generated narrative
  //
  // The page becomes a single A4-style canvas of automatically generated
  // descriptive prose, sourced from the user's scores entered elsewhere on
  // the site. Two settings live in a slim toolbar at the top: Reference and
  // Descriptor system. Everything else is derived.
  //
  // Sentence rules (from MEMORY.md + recent feedback):
  //   - Pronoun-led ("Her performance on X was Y") not impersonal/passive
  //   - Long-form score labels ("standard score of 100", not "SS=100")
  //   - Single canonical phrasing per type; cycle templates for variety
  //   - Minimal connectors (no "Furthermore," "In addition,")
  //   - Sentences shouldn't start with bare names; pronouns are fine
  // ============================================================================

  /* Domain order for the generated report (clinical convention). */
  const RWAUTO_DOMAIN_ORDER = [
    'validity',
    'premorbid',
    'intellectual',
    'attention_working_memory',
    'processing_speed',
    'verbal_learning_memory',
    'visual_learning_memory',
    'language',
    'visuospatial',
    'executive',
    'mood'
  ];
  /* Per-domain prose label for opener sentences. */
  const RWAUTO_DOMAIN_LABEL = {
    validity:                'performance validity',
    premorbid:               'premorbid intellectual functioning',
    intellectual:            'general intellectual functioning',
    attention_working_memory:'attention and working memory',
    processing_speed:        'processing speed',
    verbal_learning_memory:  'verbal learning and memory',
    visual_learning_memory:  'visual learning and memory',
    language:                'language',
    visuospatial:            'visuospatial and constructional ability',
    executive:               'executive functioning',
    mood:                    'mood and symptom reporting'
  };
  const RWAUTO_DOMAIN_HEADING = {
    validity:                'Performance Validity',
    premorbid:               'Premorbid Functioning',
    intellectual:            'General Intellectual Functioning',
    attention_working_memory:'Attention and Working Memory',
    processing_speed:        'Processing Speed',
    verbal_learning_memory:  'Verbal Learning and Memory',
    visual_learning_memory:  'Visual Learning and Memory',
    language:                'Language',
    visuospatial:            'Visuospatial and Constructional Ability',
    executive:               'Executive Functioning',
    mood:                    'Mood and Symptom Reporting'
  };

  /* Family / subtest → domain mapping. Order matters for the regex-based
     classifier — first match wins. */
  const RWAUTO_DOMAIN_RULES = [
    { re: /tomm|word memory|reliable digit|medical symptom validity|msvt/i, d: 'validity' },
    { re: /topf|opie|wtar|crawford.*allan|nart/i,                            d: 'premorbid' },

    // Memory — verbal vs visual
    { re: /(cvlt|verbal paired|word list|hopkins verbal|rey auditory|hvlt)/i,                     d: 'verbal_learning_memory' },
    { re: /logical memory/i,                                                                       d: 'verbal_learning_memory' },
    { re: /(rey complex figure.*recall|visual reproduction|brief visuospatial|bvmt)/i,             d: 'visual_learning_memory' },
    { re: /designs|family pictures|spatial recall/i,                                               d: 'visual_learning_memory' },
    { re: /(immediate|delayed|auditory|visual).*memory.*index/i,                                   d: 'verbal_learning_memory' },
    { re: /\b(imi|dmi|ammi|gmi|vmi)\b/i,                                                           d: 'verbal_learning_memory' },

    // Attention / working memory / processing speed
    { re: /digit span|letter[- ]number|spatial span|symbol span|arithmetic/i,                       d: 'attention_working_memory' },
    { re: /\b(wmi|working memory index)\b/i,                                                        d: 'attention_working_memory' },
    { re: /(coding|symbol search|cancellation)/i,                                                   d: 'processing_speed' },
    { re: /trail making|trails|dkefs.*trail|tmt/i,                                                  d: 'processing_speed' },
    { re: /\b(psi|processing speed)\b/i,                                                            d: 'processing_speed' },

    // D-KEFS Color-Word baseline conditions → Processing Speed.
    // (Must precede the language `naming` rule so "Color Naming" doesn't
    //  get swept into Language, and precede the executive `color[- ]word`
    //  rule so the baseline conditions don't get swept into Executive.)
    { re: /\b(word reading|colou?r naming)\b/i,                                                     d: 'processing_speed' },

    // D-KEFS Color-Word interference conditions (Inhibition, Inhibition/Switching) → Executive
    { re: /(color[- ]word|stroop|inhibition|set[- ]shifting)/i,                                     d: 'executive' },

    // Language
    { re: /(verbal fluency|category fluency|letter fluency|fas|cowat|naming|boston)/i,              d: 'language' },
    { re: /vocabulary|similarities|comprehension|information/i,                                     d: 'language' },
    { re: /\bvci\b/i,                                                                               d: 'language' },

    // Visuospatial
    { re: /(block design|matrix reasoning|visual puzzles|figure weights|picture concepts|picture completion)/i, d: 'visuospatial' },
    { re: /(rey complex figure copy|figure copy|object assembly|judgement of line)/i,              d: 'visuospatial' },
    { re: /\bpri\b/i,                                                                               d: 'visuospatial' },

    // Other executive measures (after the Color-Word split above).
    // D-KEFS Twenty Questions: Total Weighted Achievement, Initial Abstraction Score
    { re: /(tower|sorting|design fluency|wcst|twenty questions|weighted achievement|initial abstraction)/i, d: 'executive' },

    // Mood
    { re: /\b(hads|phq|gad|bdi|bai|dass|pcl|geriatric depression)\b|mood|anxiety|depression/i,      d: 'mood' },

    // Intellectual / IQ-level
    { re: /(fsiq|gai|full[- ]?scale|ratio iq|gen.*ability|wais|wisc|wppsi|rbans.*total|kbit)/i,     d: 'intellectual' }
  ];

  function rwAutoDetectDomain(text){
    const t = String(text || '');
    if (!t.trim()) return null;
    for (const r of RWAUTO_DOMAIN_RULES){
      if (r.re.test(t)) return r.d;
    }
    return null;
  }

  /* Get the active reference (mirrors the legacy reportReference) for
     pronoun selection. The legacy function is in app.js scope. */
  function rwAutoReference(){
    try { if (typeof reportReference === 'function') return reportReference(); } catch(e){}
    // Sensible fallback
    return { firstPerson: false, name: 'the patient', subject: 'they', possessive: 'their' };
  }
  function rwAutoPossessive(cap){
    const ref = rwAutoReference();
    // ref.possessive arrives capitalized ("Her"); lowercase first so the
    // cap=false branch returns "her", not "Her".
    const p = (ref.firstPerson ? 'your' : String(ref.possessive || 'their')).toLowerCase();
    return cap ? p.charAt(0).toUpperCase() + p.slice(1) : p;
  }
  function rwAutoSubject(cap){
    const ref = rwAutoReference();
    const s = (ref.firstPerson ? 'you' : String(ref.subject || 'they')).toLowerCase();
    return cap ? s.charAt(0).toUpperCase() + s.slice(1) : s;
  }
  function rwAutoSubjectVerb(){
    // first person uses base verb; third person uses third-person singular if needed
    // (callers will phrase their own verb forms; this is mainly a helper hook)
    return rwAutoReference().firstPerson ? 'you' : 'they';
  }

  /* Standard-score → classification (using the legacy engine for system + label). */
  function rwAutoDescriptor(ss){
    try { return (typeof reportDescriptorFor === 'function') ? reportDescriptorFor(ss) : ''; } catch(e){ return ''; }
  }
  /* Convert a row's value in any score type to a standard-score equivalent. */
  function rwAutoToStandard(score, scoreType){
    const n = parseFloat(score);
    if (!Number.isFinite(n)) return null;
    switch ((scoreType || 'standard').toLowerCase()){
      case 'standard': return n;
      case 'scaled':   return 100 + (n - 10) * 5;        // SS = 100 + (Sc - 10) × 5
      case 't':        return 100 + (n - 50) * 1.5;      // SS = 100 + (T - 50) × 1.5
      case 'z':        return 100 + n * 15;
      default:         return n;
    }
  }
  function rwAutoScoreLabel(scoreType){
    const t = String(scoreType || 'standard').toLowerCase();
    return ({ standard: 'standard score', scaled: 'scaled score', t: 'T-score', z: 'z-score', percentile: 'percentile' })[t] || 'standard score';
  }

  /* Format a list of strings with Oxford comma. */
  function rwAutoListJoin(items){
    const arr = items.filter(Boolean);
    if (arr.length === 0) return '';
    if (arr.length === 1) return arr[0];
    if (arr.length === 2) return arr[0] + ' and ' + arr[1];
    return arr.slice(0, -1).join(', ') + ', and ' + arr[arr.length - 1];
  }

  function rwAutoEsc(s){
    return String(s == null ? '' : s).replace(/[&<>"']/g, c =>
      ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  }

  /* Extract a clean brand name from a battery row's `group` string.
     "WAIS-IV Core Subtests · All Ages"   → "WAIS-IV"
     "WMS-IV Subtests · Ages 16-69"        → "WMS-IV"
     "D-KEFS Trail Making Test · All Ages" → "D-KEFS"
     Falls back to the part before " · " stripped of common qualifiers. */
  const RWAUTO_BRAND_PATTERNS = [
    /\bCVLT-3\b/i, /\bCVLT-II\b/i, /\bCVLT\b/i,
    /\bD-KEFS\b/i,
    /\bRBANS\b/i,
    /\bWAIS-V\b/i, /\bWAIS-IV\b/i, /\bWAIS-III\b/i, /\bWAIS\b/i,
    /\bWMS-V\b/i, /\bWMS-IV\b/i, /\bWMS-III\b/i, /\bWMS\b/i,
    /\bWISC-V\b/i, /\bWISC-IV\b/i, /\bWISC\b/i,
    /\bTOMM\b/i,
    /\bToPF\b/i,
    /\bOPIE-4\b/i, /\bOPIE-3\b/i, /\bOPIE\b/i
  ];
  function rwAutoExtractBrand(group){
    const g = String(group || '');
    for (const p of RWAUTO_BRAND_PATTERNS){
      const m = g.match(p);
      if (m) return m[0].toUpperCase().includes('WAIS') || m[0].toUpperCase().includes('WMS') || m[0].toUpperCase().includes('WISC') || m[0].toUpperCase().includes('CVLT') || m[0].toUpperCase().includes('OPIE') || m[0].toUpperCase().includes('RBANS')
        ? m[0].toUpperCase()
        : (m[0] === 'D-KEFS' ? 'D-KEFS' : (m[0] === 'TOMM' ? 'TOMM' : (m[0] === 'ToPF' ? 'ToPF' : m[0])));
    }
    // Fallback: strip " · ..." suffix and common qualifiers
    return g.replace(/\s*·.*$/, '')
            .replace(/\s+(Core )?Subtests?\b.*$/i, '')
            .replace(/\s+Tests?\b.*$/i, '')
            .replace(/\s+Indices?\b.*$/i, '')
            .replace(/\s+Trials?\b.*$/i, '')
            .trim();
  }

  /* Ordinal suffix for percentile (1st, 2nd, 3rd, 11th, 21st, etc.) */
  function rwAutoOrdinal(n){
    const r = Math.round(n);
    const last = r % 100;
    if (last >= 11 && last <= 13) return r + 'th';
    switch (r % 10){
      case 1: return r + 'st';
      case 2: return r + 'nd';
      case 3: return r + 'rd';
      default: return r + 'th';
    }
  }
  /* Title-case version of the score-type label ("scaled score" → "Scaled Score"). */
  function rwAutoScoreLabelCap(scoreType){
    const t = String(scoreType || 'standard').toLowerCase();
    return ({ standard: 'Standard Score', scaled: 'Scaled Score', t: 'T-Score', z: 'Z-Score', percentile: 'Percentile' })[t] || 'Standard Score';
  }
  /* Convert a standard-score equivalent to a percentile (1-99 clamp). */
  function rwAutoSsToPercentile(ss){
    if (!Number.isFinite(ss)) return null;
    let pct;
    try { pct = (typeof normCDF === 'function') ? normCDF((ss - 100) / 15) * 100 : null; }
    catch(e){ pct = null; }
    if (pct == null) return null;
    return Math.max(1, Math.min(99, Math.round(pct)));
  }

  /* User-controlled toggles for what appears inside per-test parentheticals.
     Persisted to localStorage so reloads remember the preference. */
  const RWAUTO_TOGGLES_KEY = 'rwAutoToggles_v1';
  let rwAutoShowStd = true;
  let rwAutoShowPct = false;
  function rwAutoLoadToggles(){
    try {
      const raw = localStorage.getItem(RWAUTO_TOGGLES_KEY);
      if (raw){
        const p = JSON.parse(raw);
        if (typeof p.std === 'boolean') rwAutoShowStd = p.std;
        if (typeof p.pct === 'boolean') rwAutoShowPct = p.pct;
      }
    } catch(e){}
  }
  function rwAutoSaveToggles(){
    try { localStorage.setItem(RWAUTO_TOGGLES_KEY, JSON.stringify({ std: rwAutoShowStd, pct: rwAutoShowPct })); } catch(e){}
  }

  /* Build the parenthetical for a single row, honouring the toggles.
     Returns the FULL string including surrounding " (…)" — or empty string
     when both toggles are off (no parens at all).
     `mode` controls capitalisation of the score label:
        'cap'    → "Scaled Score 10" (used in collapse + CVLT narrative)
        'lower'  → "scaled score of 10" (used in standard battery sentences) */
  function rwAutoFmtBracket(row, mode){
    const score = row.score;
    const ss = row.ss;
    const pct = rwAutoSsToPercentile(ss);
    const parts = [];
    if (rwAutoShowStd){
      if (mode === 'cap'){
        parts.push(rwAutoScoreLabelCap(row.scoreType) + ' ' + score);
      } else {
        parts.push(rwAutoScoreLabel(row.scoreType) + ' of ' + score);
      }
    }
    if (rwAutoShowPct && pct != null){
      parts.push(rwAutoOrdinal(pct) + ' percentile');
    }
    if (!parts.length) return '';
    return ' (' + parts.join(', ') + ')';
  }
  /* Compact "score summary" string for a single row, no surrounding parens.
     Honours both bracket toggles. Used inside multi-part bracket sentences
     (e.g. CVLT stable-delay combined sentence, recognition "up from LDFR"). */
  function rwAutoScoreSummary(row){
    if (!row) return '';
    const parts = [];
    if (rwAutoShowStd){
      parts.push(rwAutoScoreLabelCap(row.scoreType) + ' ' + row.score);
    }
    if (rwAutoShowPct){
      const p = rwAutoSsToPercentile(row.ss);
      if (p != null) parts.push(rwAutoOrdinal(p) + ' percentile');
    }
    return parts.join(', ');
  }

  /* Plural-aware bracket for a list of same-descriptor rows (different scores).
     Used by rwAutoBatterySentences when scores differ across the cluster.
     Rule 6 amend 2: when both bracket toggles are on AND multiple scores are
     listed, percentiles are introduced with "at the …" and pluralised at the
     end ("at the 2nd, 5th, and 5th percentiles") rather than the previous
     comma-joined "percentiles 2nd, 5th, and 5th" form. */
  function rwAutoFmtBracketList(items){
    if (!items.length) return '';
    const scores = items.map(r => r.score);
    const parts = [];
    if (rwAutoShowStd){
      const label = rwAutoScoreLabel(items[0].scoreType);
      parts.push(label + (scores.length > 1 ? 's of ' : ' of ') + rwAutoListJoin(scores));
    }
    if (rwAutoShowPct){
      const pcts = items.map(r => rwAutoSsToPercentile(r.ss)).filter(p => p != null);
      if (pcts.length){
        if (pcts.length > 1){
          parts.push('at the ' + rwAutoListJoin(pcts.map(rwAutoOrdinal)) + ' percentiles');
        } else {
          parts.push('at the ' + rwAutoOrdinal(pcts[0]) + ' percentile');
        }
      }
    }
    if (!parts.length) return '';
    return ' (' + parts.join(', ') + ')';
  }

  /* Subject token for sentence index `i` within a domain section.
     Even index (0, 2, 4…) → full possessive name ("Mrs Doe's").
     Odd index (1, 3, 5…)  → possessive pronoun ("Her").
     For first-person reference, always "Your" (no rotation needed). */
  function rwAutoFullPossessive(){
    const ref = rwAutoReference();
    if (ref.firstPerson) return 'Your';
    if (ref.name){
      try { return (typeof reportPossessiveName === 'function') ? reportPossessiveName(ref.name) : (ref.name + "'s"); }
      catch(e){ return ref.name + "'s"; }
    }
    return rwAutoPossessive(true);
  }
  function rwAutoSubjectFor(i){
    const ref = rwAutoReference();
    if (ref.firstPerson) return { possessive: 'Your', possLower: 'your' };
    return (i % 2 === 0)
      ? { possessive: rwAutoFullPossessive(),     possLower: rwAutoFullPossessive() }
      : { possessive: rwAutoPossessive(true),     possLower: rwAutoPossessive(false) };
  }

  /* Read the user's selected CI from the Premorbid page. */
  function rwAutoReadCi(){
    const sel = document.getElementById('pre-ci');
    const val = sel ? parseFloat(sel.value) : 0.90;
    return val === 0.95 ? { z: 1.96, label: '95%' } : { z: 1.645, label: '90%' };
  }

  /* Strip the help-icon "?" glyph that the legacy table appends to model names. */
  function rwAutoCleanModelName(s){
    return String(s || '')
      .replace(/[?¿]+\s*$/g, '')         // trailing "?" tooltips
      .replace(/\s+/g, ' ')
      .trim();
  }

  /* Map a model row to its display name + canonical sort order. */
  const RWAUTO_PREMODEL = {
    'topf-demo':  { order: 1, label: 'the Test of Premorbid Functioning (ToPF) with demographic adjustment' },
    'topf-raw':   { order: 2, label: 'the ToPF (raw score)' },
    'crawford':   { order: 3, label: 'the Crawford & Allan (2001) demographic predictor' },
    'opie4':      { order: 4, label: 'the OPIE-4 (Vocabulary and/or Matrix)' }
  };
  function rwAutoIdentifyPreModel(rawName){
    const n = rwAutoCleanModelName(rawName).toLowerCase();
    if (n.startsWith('topf') && /demographic|demo/.test(n))    return 'topf-demo';
    if (n.startsWith('topf') && /raw/.test(n))                  return 'topf-raw';
    if (n.startsWith('crawford'))                               return 'crawford';
    if (n.startsWith('opie'))                                   return 'opie4';
    return null;
  }

  /* ----- Read each data source ----- */

  /* Battery rows → cognitive measures, grouped by detected domain. */
  function rwAutoCollectBattery(){
    let rows = [];
    try { rows = (typeof batteryRows !== 'undefined' && Array.isArray(batteryRows)) ? batteryRows : []; } catch(e){}
    return rows.filter(r => r && r.name && r.score !== '' && r.score != null).map(r => {
      const ss = rwAutoToStandard(r.score, r.scoreType);
      const desc = ss != null ? rwAutoDescriptor(ss) : '';
      const familyText = r.group || '';
      const subtestText = r.name || '';
      const brand = rwAutoExtractBrand(familyText);
      const domain = rwAutoDetectDomain(subtestText) || rwAutoDetectDomain(familyText) || 'intellectual';
      return {
        kind: 'battery',
        family: familyText,
        brand: brand || familyText,
        subtest: subtestText,
        score: r.score,
        scoreType: r.scoreType || 'standard',
        ss,
        descriptor: desc,
        domain
      };
    });
  }

  /* Change rows from SDI + RCI methods, mapped to the test family they
     describe. Each gets a "change finding" string we'll embed under the
     matching cognitive domain. */
  function rwAutoCollectChange(){
    const out = [];
    try {
      if (typeof sdiRows !== 'undefined' && Array.isArray(sdiRows)){
        sdiRows.forEach(r => {
          if (!r || !r.name) return;
          const t1 = parseFloat(r.t1), t2 = parseFloat(r.t2);
          if (!Number.isFinite(t1) || !Number.isFinite(t2)) return;
          const sd = parseFloat(r.sd);
          if (!Number.isFinite(sd) || sd <= 0) return;
          const delta = (t2 - t1) / sd;
          const family = r.group || '';
          const domain = rwAutoDetectDomain(r.name) || rwAutoDetectDomain(family) || null;
          if (!domain) return;
          out.push({ kind: 'sdi', name: r.name, family, t1, t2, delta, domain });
        });
      }
    } catch(e){}
    try {
      if (typeof rciState !== 'undefined' && rciState){
        Object.entries(rciState).forEach(([method, state]) => {
          (state.rows || []).forEach(r => {
            if (!r || !r.name) return;
            const family = r.group || '';
            const domain = rwAutoDetectDomain(r.name) || rwAutoDetectDomain(family) || null;
            if (!domain) return;
            // Compute display values via the legacy calc helpers. Calc fns
            // return the test statistic on the `rci` field (a z for basic /
            // practice / SRB; a t-stat for Crawford). They don't compute a
            // pass/fail flag — derive it from p vs the method's chosen cv.
            let stat = null, p = null, sig = null;
            try {
              const calc = (method === 'rci-basic'    && typeof calcBasicRow === 'function')        ? calcBasicRow(r, method)
                         : (method === 'rci-practice' && typeof calcPracticeRow === 'function')     ? calcPracticeRow(r, method)
                         : (method === 'rci-srb'      && typeof calcSrbRow === 'function')          ? calcSrbRow(r, method)
                         : (method === 'rci-crawford' && typeof calcCrawfordRow === 'function')     ? calcCrawfordRow(r, method)
                         : null;
              if (calc){
                stat = calc.rci;
                p = calc.p;
                const cv = (state && Number.isFinite(parseFloat(state.cv))) ? parseFloat(state.cv) : 0.95;
                sig = (p != null && Number.isFinite(p)) ? (p < (1 - cv)) : null;
              }
            } catch(e){}
            out.push({ kind: method, name: r.name, family, stat, p, sig, domain });
          });
        });
      }
    } catch(e){}
    return out;
  }

  /* Read premorbid inputs + computed estimates from the Premorbid page.
     Returns FSIQ + CI bounds per model. Estimates are sorted in canonical
     report order: ToPF demographic → ToPF raw → Crawford → OPIE-4. */
  function rwAutoCollectPremorbid(){
    const ids = ['pre-topf','pre-vc','pre-mr','pre-sex','pre-occ','pre-edu','pre-age'];
    const filled = ids.filter(id => {
      const el = document.getElementById(id);
      return el && (el.value != null && String(el.value).trim() !== '');
    });
    if (!filled.length) return null;

    const estimates = [];
    document.querySelectorAll('#pre-results-table tbody tr').forEach(tr => {
      const tds = tr.querySelectorAll('td');
      if (tds.length < 4) return;
      const rawModel = rwAutoCleanModelName(tds[0]?.textContent || '');
      const fsiq = parseFloat((tds[1]?.textContent || '').replace(/[^\d.\-]/g, ''));
      const lo   = parseFloat((tds[2]?.textContent || '').replace(/[^\d.\-]/g, ''));
      const hi   = parseFloat((tds[3]?.textContent || '').replace(/[^\d.\-]/g, ''));
      if (!rawModel || !Number.isFinite(fsiq)) return;
      const id = rwAutoIdentifyPreModel(rawModel);
      const meta = id ? RWAUTO_PREMODEL[id] : null;
      estimates.push({
        id, rawModel,
        label: meta ? meta.label : ('the ' + rawModel),
        order: meta ? meta.order : 99,
        fsiq: Math.round(fsiq),
        lo: Number.isFinite(lo) ? Math.round(lo) : null,
        hi: Number.isFinite(hi) ? Math.round(hi) : null
      });
    });
    estimates.sort((a, b) => a.order - b.order);

    // OPIE-4 label is dynamic — depends on which subtests the user supplied
    // for the calculation (Vocabulary, Matrix Reasoning, or both).
    const vc = document.getElementById('pre-vc');
    const mr = document.getElementById('pre-mr');
    const hasVc = vc && String(vc.value || '').trim() !== '';
    const hasMr = mr && String(mr.value || '').trim() !== '';
    estimates.forEach(e => {
      if (e.id !== 'opie4') return;
      let inputs;
      if (hasVc && hasMr)      inputs = 'Vocabulary and Matrix Reasoning';
      else if (hasVc)           inputs = 'Vocabulary';
      else if (hasMr)           inputs = 'Matrix Reasoning';
      else                      inputs = 'Vocabulary and/or Matrix Reasoning';
      e.label = 'the OPIE-4 (' + inputs + ')';
    });

    return { hasInputs: true, estimates };
  }

  /* ----- Sentence builders ----- */

  /* First-paragraph opener for a domain section.
     If subtest count > 4, collapse to "tests from the [brand]" rather than
     listing every subtest by name. */
  function rwAutoBrandOpenerFirst(domainId, brand, rows, i){
    const subtests = Array.from(new Set(rows.map(r => r.subtest).filter(Boolean)));
    const ability = RWAUTO_DOMAIN_LABEL[domainId] || 'this domain';
    const subj = rwAutoSubjectFor(i);
    const using = (subtests.length > 4)
      ? 'tests from the ' + brand
      : rwAutoListJoin(subtests) + ' from the ' + brand;
    return subj.possessive + ' ' + ability + ' was assessed using ' + using + '.';
  }

  /* Subsequent-paragraph opener within a domain.
     Rule 8 amend: when the subject cycle lands on the FULL-NAME position, the
     bare subject-verb construction "Mrs Doe also undertook…" violates rule 2.
     Use the same "[Subject's] [domain] was also assessed using …" form as the
     first opener instead. Pronoun position keeps "She also undertook…". */
  function rwAutoBrandOpenerNext(domainId, brand, rows, i){
    const subtests = Array.from(new Set(rows.map(r => r.subtest).filter(Boolean)));
    const tail = (subtests.length > 4)
      ? 'tests from the ' + brand
      : rwAutoListJoin(subtests) + ' from the ' + brand;
    const ref = rwAutoReference();
    if (ref.firstPerson){
      return 'You also undertook ' + tail + '.';
    }
    // Even cycle index → full-name position. Use the possessive-name +
    // domain-ability form so the sentence isn't a bare subject-verb.
    if (i % 2 === 0){
      const possName = rwAutoFullPossessive();
      const ability = RWAUTO_DOMAIN_LABEL[domainId] || 'this domain';
      return possName + ' ' + ability + ' was also assessed using ' + tail + '.';
    }
    // Odd cycle index → pronoun position (unchanged: "She also undertook…").
    return rwAutoSubject(true) + ' also undertook ' + tail + '.';
  }

  /* Is the brand a CVLT variant? */
  function rwAutoIsCvltBrand(brand){
    return /\bCVLT(-3|-II)?\b/i.test(String(brand || ''));
  }

  /* CVLT-style verbal-learning narrative — replaces the generic per-test
     listing for verbal_learning_memory + CVLT brand combinations. Reports:
       1. Learning curve (first → last immediate trial)
       2. Short-delay recall (SDFR / SDCR)
       3. Long-delay recall (LDFR / LDCR)
       4. Recognition assistance (slightly | notably | plain) */
  function rwAutoCvltNarrative(rows, startIdx){
    const find = (re) => rows.find(r => re.test(r.subtest || ''));
    const trialNum = s => {
      // Match "Trial 3", "Trial3", "T3" (followed by space/end), but NOT
      // composites like "T1-5 Correct" where the digit is followed by "-".
      const m = String(s || '').match(/^Trial\s*([1-5])(?!\d)/i);
      return m ? parseInt(m[1], 10) : null;
    };
    // Strict learning-curve trigger: only "Trial 1" through "Trial 5". Excludes
    // composite measures like "T1-5 Correct" (which is a Standard Score, M=100,
    // not a per-trial scaled score). Engine bug fix per rule 18 brief.
    const isImmediateTrial = r => /^Trial\s+[1-5]\b/i.test(r.subtest || '');
    const trials = rows.filter(isImmediateTrial)
                       .sort((a, b) => (trialNum(a.subtest) || 0) - (trialNum(b.subtest) || 0));
    const sdfr = find(/short delay free recall|sdfr/i);
    const sdcr = find(/short delay cued recall|sdcr/i);
    const ldfr = find(/long delay free recall|ldfr/i);
    const ldcr = find(/long delay cued recall|ldcr/i);
    // Main recognition = the simple Hits / Total Recognition score. The
    // discrimination indices (parametric / nonparametric) and false-positive
    // count are auxiliary measures.
    const recog = find(/^recognition$|^recognition\s+hits?\b|^total recognition\b/i);

    // Identify which rows are "main" so the rest can flow into auxiliary or
    // composite paragraphs. Anything not in mainRows is split by score type:
    //   scoreType = 'standard'  → composite (Total Recall Correct, etc.)
    //   anything else           → auxiliary (intrusions, repetitions, etc.)
    // Rule 18: composite measures get their own paragraph after auxiliary.
    const mainRows = new Set();
    if (sdfr) mainRows.add(sdfr);
    if (sdcr) mainRows.add(sdcr);
    if (ldfr) mainRows.add(ldfr);
    if (ldcr) mainRows.add(ldcr);
    if (recog) mainRows.add(recog);
    trials.forEach(t => mainRows.add(t));
    const remainder = rows.filter(r => !mainRows.has(r));
    const composites = remainder.filter(r => (r.scoreType || '').toLowerCase() === 'standard');
    const compositeSet = new Set(composites);
    const auxiliary = remainder.filter(r => !compositeSet.has(r));

    const sentences = [];
    let i = startIdx;
    // CVLT narrative uses capitalised "Scaled Score" form for parens
    const paren = (r) => rwAutoFmtBracket(r, 'cap');

    // 1. Learning curve — Rule 10 (ranged) OR Rule 13 (stable when first/last
    //    fall in the same descriptor band)
    if (trials.length >= 2){
      const first = trials[0], last = trials[trials.length - 1];
      const subj = rwAutoSubjectFor(i);
      const stable = first.descriptor && last.descriptor && first.descriptor === last.descriptor;
      if (stable){
        sentences.push(
          subj.possessive + ' performance was stable across the five learning trials, ' +
          first.descriptor.toLowerCase() + ' on both ' + first.subtest + paren(first) +
          ' and ' + last.subtest + paren(last) + '.'
        );
      } else {
        sentences.push(
          subj.possessive + ' performance ranged from ' + first.descriptor.toLowerCase() +
          ' on ' + first.subtest + paren(first) + ' to ' +
          last.descriptor.toLowerCase() + ' by ' + last.subtest + paren(last) + '.'
        );
      }
      i++;
    }

    // 2/3. Delays — Rule 10 separate sentences OR Rule 14 combined
    //      (when SDFR and LDFR fall in the same descriptor band).
    const shortPrimary = sdfr || sdcr;
    const longPrimary  = ldfr || ldcr;
    const stableDelays = sdfr && ldfr && sdfr.descriptor && sdfr.descriptor === ldfr.descriptor;

    if (stableDelays){
      const subj = rwAutoSubjectFor(i);
      // Combined bracket: "Scaled Score X for SDFR, Scaled Score Y for LDFR"
      const sdfrSum = rwAutoScoreSummary(sdfr);
      const ldfrSum = rwAutoScoreSummary(ldfr);
      const combined = (sdfrSum && ldfrSum)
        ? ' (' + sdfrSum + ' for ' + sdfr.subtest + ', ' + ldfrSum + ' for ' + ldfr.subtest + ')'
        : '';
      sentences.push(
        subj.possessive + ' performance after both short and long delays was ' +
        sdfr.descriptor.toLowerCase() + combined + '.'
      );
      i++;
    } else {
      // Rule 10 separate sentences
      if (shortPrimary){
        const subj = rwAutoSubjectFor(i);
        sentences.push(
          'After a short delay, ' + subj.possLower + ' performance on ' +
          shortPrimary.subtest + ' was ' + shortPrimary.descriptor.toLowerCase() +
          paren(shortPrimary) + '.'
        );
        i++;
      }
      if (longPrimary){
        const subj = rwAutoSubjectFor(i);
        sentences.push(
          'Finally, after a long delay, ' + subj.possLower + ' performance on ' +
          longPrimary.subtest + ' was ' + longPrimary.descriptor.toLowerCase() +
          paren(longPrimary) + '.'
        );
        i++;
      }
    }

    // 4. Recognition — Rule 16 (assisted, with "up from LDFR") or Rule 15
    //    (delta < 3 → "not significantly assisted"). LDFR is needed for the
    //    delta calc; without it, falls back to a plain recognition report.
    if (recog){
      const subj = rwAutoSubjectFor(i);
      const diff = longPrimary
        ? parseFloat(recog.score) - parseFloat(longPrimary.score)
        : NaN;
      let phrase = null;
      let includeUpFrom = false;
      if (Number.isFinite(diff)){
        if (diff >= 6){      phrase = 'notably assisted by the recognition trial'; includeUpFrom = true; }
        else if (diff >= 3){ phrase = 'slightly assisted by the recognition trial'; includeUpFrom = true; }
        else                 phrase = 'not significantly assisted by the recognition trial';
      } else {
        // No LDFR available — fall back to a plain recognition report
        sentences.push(subj.possessive + ' recognition performance was ' + recog.descriptor.toLowerCase() + paren(recog) + '.');
        i++;
      }

      if (phrase){
        // Rule 16: include "up from {LDFR summary}" when assistance fired AND
        // at least one bracket toggle is on. When both toggles are off, the
        // entire parenthetical (including up-from) is omitted.
        const recogSum = rwAutoScoreSummary(recog);
        const ldfrSum  = (includeUpFrom && longPrimary) ? rwAutoScoreSummary(longPrimary) : '';
        let bracket = '';
        if (recogSum){
          const tail = ldfrSum
            ? recogSum + ', up from ' + ldfrSum + ' on long delay free recall'
            : recogSum;
          bracket = ' (' + tail + ')';
        }
        sentences.push(subj.possessive + ' performance was ' + phrase + bracket + '.');
        i++;
      }
    }

    // Return main narrative + auxiliary + composite rows. Composer renders
    // them as three separate paragraphs (main → auxiliary → composites).
    return { sentences, nextIdx: i, auxiliary, composites };
  }

  /* Per-descriptor-group result sentences with subject cycle and identical-
     score collapse. When all rows in a descriptor cluster share the SAME
     score, prose collapses to:
        "[Subject's] performance on X, Y and Z was all in the average range
         (Scaled Score 10, 50th percentile)."
     Otherwise the existing list form is used. */
  function rwAutoBatterySentences(rows, startIdx){
    if (!rows.length) return { sentences: [], nextIdx: startIdx };
    const groups = {};
    rows.forEach(r => {
      if (!r.descriptor || r.ss == null) return;
      const key = (r.descriptor || '').toLowerCase() + '|' + (r.scoreType || 'standard');
      if (!groups[key]) groups[key] = { descriptor: r.descriptor, scoreType: r.scoreType, items: [] };
      groups[key].items.push(r);
    });
    const sentences = [];
    let i = startIdx;
    Object.keys(groups).forEach(k => {
      const g = groups[k];
      const labels = g.items.map(r => r.subtest);
      const scores = g.items.map(r => r.score);
      const labelList = rwAutoListJoin(labels);
      const subj = rwAutoSubjectFor(i);
      const descLower = g.descriptor.toLowerCase();

      const allSameScore = scores.length > 1 && scores.every(s => String(s) === String(scores[0]));

      if (allSameScore){
        // Collapse to singular: capitalised form, percentile in brackets per
        // toggle. Rule 6 amend 1: two-item form uses "were both in the X range";
        // three or more items keep the existing "was all in the X range" form.
        const bracket = rwAutoFmtBracket(g.items[0], 'cap');
        const verbPhrase = (g.items.length === 2)
          ? ' were both in the '
          : ' was all in the ';
        sentences.push(subj.possessive + ' performance on ' + labelList +
                       verbPhrase + descLower + ' range' + bracket + '.');
      } else if (g.items.length === 1){
        const bracket = rwAutoFmtBracket(g.items[0], 'lower');
        sentences.push(subj.possessive + ' performance on ' + labelList +
                       ' was ' + descLower + bracket + '.');
      } else {
        // Multiple different scores in same descriptor cluster
        const bracket = rwAutoFmtBracketList(g.items);
        sentences.push(subj.possessive + ' performance on ' + labelList +
                       ' was ' + descLower + bracket + '.');
      }
      i++;
    });
    return { sentences, nextIdx: i };
  }

  /* Rule 19 — change-analysis paragraph(s) for one domain.
     One paragraph per (battery × method) combination. Subject cycle resets at
     the top of each paragraph (opener at idx=0 → full possessive name). The
     stat / p bracket is always shown regardless of rule 12 toggles. */

  /* Human-readable label for each method, used in the opener after "using a/an". */
  const RWAUTO_METHOD_LABEL = {
    'sdi':          'SD discrepancy analysis',
    'rci-basic':    'Simple Reliable Change Index',
    'rci-practice': 'Practice-Adjusted Reliable Change Index',
    'rci-srb':      'McSweeney Regression-Based change analysis',
    'rci-crawford': 'Crawford & Garthwaite regression-based change analysis'
  };
  /* Stat label used inside the per-test bracket. */
  const RWAUTO_STAT_LABEL = {
    'sdi':          'SD-Δ',
    'rci-basic':    'z',
    'rci-practice': 'z',
    'rci-srb':      'z',
    'rci-crawford': 't(RB)'
  };
  /* Pick "a"/"an" by first letter of the type label. */
  function rwAutoIndefArticle(s){
    const c = String(s || '').trim().charAt(0).toLowerCase();
    return /[aeiou]/.test(c) ? 'an' : 'a';
  }
  /* Date labels for a given method. RCI methods carry these on rciState; SDI's
     labels live on the date-header inputs. Falls back to "Date 1"/"Date 2". */
  function rwAutoMethodDates(kind){
    if (kind === 'sdi'){
      const d1 = document.getElementById('sdi-d1-head')?.value;
      const d2 = document.getElementById('sdi-d2-head')?.value;
      return { d1: d1 || 'Date 1', d2: d2 || 'Date 2' };
    }
    try {
      if (typeof rciState !== 'undefined' && rciState[kind]){
        return { d1: rciState[kind].d1 || 'Date 1', d2: rciState[kind].d2 || 'Date 2' };
      }
    } catch(e){}
    return { d1: 'Date 1', d2: 'Date 2' };
  }
  /* Format the (stat, p) bracket. Always shown — not governed by rule 12. */
  function rwAutoChangeBracket(c){
    let stat, p;
    if (c.kind === 'sdi'){
      stat = c.delta;
      // Approximate p from |delta| via standard normal (consistent with
      // basic-RCI two-tailed convention)
      try { p = (typeof normCDF === 'function') ? 2 * (1 - normCDF(Math.abs(c.delta))) : null; }
      catch(e){ p = null; }
    } else {
      stat = c.stat;
      p = c.p;
    }
    const statStr = (stat != null && Number.isFinite(stat))
      ? (RWAUTO_STAT_LABEL[c.kind] || 'stat') + ' = ' + stat.toFixed(2)
      : null;
    const pStr = (p != null && Number.isFinite(p))
      ? (p < 0.001 ? '<i>p</i> < .001' : '<i>p</i> = ' + p.toFixed(3).replace(/^0\./, '.'))
      : null;
    const parts = [];
    if (statStr) parts.push(statStr);
    if (pStr) parts.push(pStr);
    return parts.length ? ' (' + parts.join(', ') + ')' : '';
  }
  /* Decide outcome wording. For SDI uses |delta| ≥ 1.96 as the reliable
     threshold (standard 95% z-cutoff); for RCI methods uses the calc's `sig`. */
  function rwAutoChangeOutcome(c){
    let sig, stat;
    if (c.kind === 'sdi'){
      stat = c.delta;
      sig = Number.isFinite(c.delta) && Math.abs(c.delta) >= 1.96;
    } else {
      stat = c.stat;
      sig = c.sig === true || c.sig === 'sig';
    }
    if (!sig) return 'did not reliably change';
    return (Number.isFinite(stat) && stat < 0) ? 'reliably declined' : 'reliably improved';
  }

  function rwAutoChangeOpener(brand, kind, i){
    const subj = rwAutoSubjectFor(i);
    const dates = rwAutoMethodDates(kind);
    const methodLabel = RWAUTO_METHOD_LABEL[kind] || 'change analysis';
    const article = rwAutoIndefArticle(methodLabel);
    return subj.possessive + ' performance on the ' + brand + ' in ' + dates.d1 +
           ' was compared against ' + rwAutoPossessive(false) +
           ' current performance (' + dates.d2 + ') using ' +
           article + ' ' + methodLabel + '.';
  }

  function rwAutoChangeSentence(c, i){
    const subj = rwAutoSubjectFor(i);
    const dates = rwAutoMethodDates(c.kind);
    const dirPhrase = rwAutoChangeOutcome(c);
    const bracket = rwAutoChangeBracket(c);
    return subj.possessive + ' performance on ' + c.name + ' ' + dirPhrase +
           ' since ' + dates.d1 + bracket + '.';
  }

  /* ── Rule 19 CVLT cluster helpers ─────────────────────────── */

  function rwAutoIsReliable(c){
    if (c.kind === 'sdi') return Number.isFinite(c.delta) && Math.abs(c.delta) >= 1.96;
    return c.sig === true || c.sig === 'sig';
  }

  function rwAutoChangeDir(c){
    if (!rwAutoIsReliable(c)) return null;
    const v = c.kind === 'sdi' ? c.delta : c.stat;
    return (Number.isFinite(v) && v < 0) ? 'decline' : 'improvement';
  }

  function cvltNameList(items){
    const ns = items.map(c => c.name);
    if (ns.length === 1) return ns[0];
    if (ns.length === 2) return ns[0] + ' and ' + ns[1];
    return ns.slice(0,-1).join(', ') + ', and ' + ns[ns.length-1];
  }

  /* "t(RB) range: A.AA to B.BB" across items (sorted ascending). */
  function cvltStatRangeStr(items, statLabel){
    const vals = items.map(c => c.kind === 'sdi' ? c.delta : c.stat).filter(v => Number.isFinite(v));
    if (!vals.length) return '';
    vals.sort((a,b) => a - b);
    return statLabel + ' range: ' + vals[0].toFixed(2) + ' to ' + vals[vals.length-1].toFixed(2);
  }

  /* p-value string, value only, no italic prefix: '.003' or '< .001'. */
  function cvltPStr(p){
    if (p == null || !Number.isFinite(p)) return '';
    if (p < 0.001) return '< .001';
    return p.toFixed(3).replace(/^0\./, '.');
  }

  /* Paired stat bracket for same-direction delay sentences:
     " (t(RB) = X and Y, p = P1 and P2)" */
  function cvltPairedBracket(m1, m2, sl){
    const v1 = m1.kind === 'sdi' ? m1.delta : m1.stat;
    const v2 = m2.kind === 'sdi' ? m2.delta : m2.stat;
    const p1 = cvltPStr(m1.p);
    const p2 = cvltPStr(m2.p);
    let out = ' (';
    if (Number.isFinite(v1) && Number.isFinite(v2)){
      out += sl + ' = ' + v1.toFixed(2) + ' and ' + v2.toFixed(2);
      if (p1 && p2) out += ', <i>p</i> = ' + p1 + ' and ' + p2;
    }
    return out + ')';
  }

  /* CVLT-specific cluster-based change paragraph.
     Implements Rule 19 with Learning / Short Delay / Long Delay /
     Recognition / Auxiliary cluster logic. */
  function rwAutoChangeParagraphCvlt(g){
    const kind  = g.kind;
    const brand = g.brand;
    const dates = rwAutoMethodDates(kind);
    const sl    = RWAUTO_STAT_LABEL[kind] || 't(RB)';

    const byName = {};
    g.items.forEach(c => { byName[c.name] = c; });

    const CLUSTERS = [
      { id:'learning',    label:'Learning Trials', size:'five',  tests:['Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'] },
      { id:'short_delay', label:'Short Delay',      size:'two',   tests:['Short Delay Free Recall','Short Delay Cued Recall'] },
      { id:'long_delay',  label:'Long Delay',       size:'two',   tests:['Long Delay Free Recall','Long Delay Cued Recall'] },
      { id:'recognition', label:'Recognition',      size:'four',  tests:['Recognition','Recognition False Positive','Recognition Discrimination','Discrimination Nonparametric'] },
      { id:'auxiliary',   label:'Auxiliary',        size:'three', tests:['List B Correct','Total Intrusions','Total Repetitions'] },
    ];

    let idx = 0;
    const sents = [];
    const skipped = [];

    sents.push(rwAutoChangeOpener(brand, kind, idx));
    idx++;

    CLUSTERS.forEach(cl => {
      const present = cl.tests.map(t => byName[t]).filter(Boolean);
      if (!present.length) return;

      const changed   = present.filter(c => rwAutoIsReliable(c));
      if (!changed.length){ skipped.push(cl.label); return; }

      const declines  = changed.filter(c => rwAutoChangeDir(c) === 'decline');
      const improves  = changed.filter(c => rwAutoChangeDir(c) === 'improvement');
      const dom       = declines.length >= improves.length ? 'decline' : 'improvement';
      const majority  = dom === 'decline' ? declines : improves;
      const minority  = dom === 'decline' ? improves : declines;
      const nullItems = present.filter(c => !rwAutoIsReliable(c));

      if (cl.id === 'learning'){
        if (changed.length >= 2){
          // Most pronounced = highest |stat| in majority direction
          const top = [...majority].sort((a,b) =>
            Math.abs(b.kind === 'sdi' ? b.delta : b.stat) -
            Math.abs(a.kind === 'sdi' ? a.delta : a.stat))[0];
          // Count sentence — does NOT advance subject cycle
          sents.push(majority.length + ' of the ' + cl.size +
            ' learning trials showed reliable ' + dom + ' since ' + dates.d1 +
            ', most pronounced in ' + top.name + rwAutoChangeBracket(top) + '.');
          // However: null items first, no statistics
          if (nullItems.length){
            const subj = rwAutoSubjectFor(idx); idx++;
            sents.push('However, ' + subj.possessive + ' performance on ' +
              cvltNameList(nullItems) + ' did not reliably change.');
          }
          // However: opposite direction, with statistics
          if (minority.length){
            const oppVerb = dom === 'decline' ? 'reliably improved' : 'reliably declined';
            const bkt = minority.length === 1
              ? rwAutoChangeBracket(minority[0])
              : ' (' + cvltStatRangeStr(minority, sl) + ')';
            const subj = rwAutoSubjectFor(idx); idx++;
            sents.push('However, ' + subj.possessive + ' performance on ' +
              cvltNameList(minority) + ' ' + oppVerb + ' since ' + dates.d1 + bkt + '.');
          }
        } else {
          sents.push(rwAutoChangeSentence(changed[0], idx)); idx++;
        }

      } else if (cl.id === 'short_delay' || cl.id === 'long_delay'){
        if (present.length < 2){
          sents.push(rwAutoChangeSentence(changed[0], idx)); idx++;
        } else {
          const [m1, m2] = present;
          const d1 = rwAutoChangeDir(m1);
          const d2 = rwAutoChangeDir(m2);
          if (d1 !== null && d2 !== null && d1 === d2){
            // Both same direction — combined sentence
            const verb = d1 === 'decline' ? 'reliably declined' : 'reliably improved';
            const subj = rwAutoSubjectFor(idx); idx++;
            sents.push(subj.possessive + ' performance on ' + m1.name + ' and ' + m2.name +
              ' ' + verb + ' since ' + dates.d1 + cvltPairedBracket(m1, m2, sl) + '.');
          } else {
            // Mixed or one null — individual sentences
            [m1, m2].forEach(m => {
              if (rwAutoIsReliable(m)){
                sents.push(rwAutoChangeSentence(m, idx)); idx++;
              } else {
                const subj = rwAutoSubjectFor(idx); idx++;
                sents.push(subj.possessive + ' performance on ' + m.name +
                  ' did not reliably change.');
              }
            });
          }
        }

      } else {
        // Recognition or Auxiliary
        const clLabel = cl.id === 'recognition' ? 'recognition' : 'auxiliary';
        if (changed.length >= 2){
          const rangeStr = cvltStatRangeStr(majority, sl);
          // Count sentence — does NOT advance subject cycle
          sents.push(majority.length + ' of the ' + cl.size +
            ' ' + clLabel + ' measures showed reliable ' + dom + ' since ' +
            dates.d1 + ' (' + rangeStr + ').');
          if (nullItems.length){
            const subj = rwAutoSubjectFor(idx); idx++;
            sents.push('However, ' + subj.possessive + ' performance on ' +
              cvltNameList(nullItems) + ' did not reliably change.');
          }
          if (minority.length){
            const oppVerb = dom === 'decline' ? 'reliably improved' : 'reliably declined';
            const bkt = minority.length === 1
              ? rwAutoChangeBracket(minority[0])
              : ' (' + cvltStatRangeStr(minority, sl) + ')';
            const subj = rwAutoSubjectFor(idx); idx++;
            sents.push('However, ' + subj.possessive + ' performance on ' +
              cvltNameList(minority) + ' ' + oppVerb + ' since ' + dates.d1 + bkt + '.');
          }
        } else {
          sents.push(rwAutoChangeSentence(changed[0], idx)); idx++;
        }
      }
    });

    // Trailing sentence: clusters present but with zero reliable change
    if (skipped.length === 1){
      sents.push('No reliable change was observed in the ' + skipped[0] + ' measures.');
    } else if (skipped.length >= 2){
      const last = skipped[skipped.length-1];
      const rest = skipped.slice(0,-1);
      sents.push('No reliable change was observed in the ' +
        rest.join(', ') + ', and ' + last + ' measures.');
    }

    return '<p>' + sents.join(' ') + '</p>';
  }

  /* Build the change paragraphs for one domain's worth of change rows. Returns
     an array of <p>…</p> HTML strings (zero or more), one per (brand, kind).
     CVLT brands get the cluster-based Rule 19 treatment; others get the generic
     per-test sentences. */
  function rwAutoChangeParagraphs(changes){
    if (!changes.length) return [];
    const order = [];
    const groups = new Map();
    changes.forEach(c => {
      const brand = rwAutoExtractBrand(c.family) || 'this assessment';
      const key = brand + '|' + c.kind;
      if (!groups.has(key)){
        groups.set(key, { brand, kind: c.kind, items: [] });
        order.push(key);
      }
      groups.get(key).items.push(c);
    });
    const paragraphs = [];
    order.forEach(key => {
      const g = groups.get(key);
      if (rwAutoIsCvltBrand(g.brand)){
        paragraphs.push(rwAutoChangeParagraphCvlt(g));
      } else {
        let idx = 0;
        const opener = rwAutoChangeOpener(g.brand, g.kind, idx);
        idx++;
        const sentences = [opener];
        g.items.forEach(c => {
          sentences.push(rwAutoChangeSentence(c, idx));
          idx++;
        });
        paragraphs.push('<p>' + sentences.join(' ') + '</p>');
      }
    });
    return paragraphs;
  }

  /* Validity section — Rule 17.
     Single general statement; no per-test scores listed.
     Pass: "[Subject's] performance on all validity measures was within
            acceptable limits, and the findings below are likely a reliable
            reflection of [their] abilities."
     Fail: "[Subject's] performance on validity testing fell below acceptable
            limits, and the findings below should be interpreted with caution." */
  function rwAutoValiditySection(rows){
    if (!rows.length) return '';
    const failed = rows.find(r => r.ss != null && r.ss < 45);
    const subj0 = rwAutoSubjectFor(0); // even index → full possessive name
    const possLow = rwAutoPossessive(false); // pronoun "her"/"his"/"your"
    const heading = '<h3>' + RWAUTO_DOMAIN_HEADING.validity + '</h3>';
    if (failed){
      return heading +
        '<p>' + subj0.possessive + ' performance on validity testing fell below acceptable limits, and the findings below should be interpreted with caution.</p>';
    }
    return heading +
      '<p>' + subj0.possessive + ' performance on all validity measures was within acceptable limits, and the findings below are likely a reliable reflection of ' + possLow + ' abilities.</p>';
  }

  /* Premorbid section in prescribed format:
       "Mrs Doe's baseline cognitive ability was estimated using the Test
        of Premorbid Functioning (ToPF) with demographic adjustment
        (FSIQ ≈ 102, 90% CI 88-116), as well as the ToPF raw score only
        (FSIQ ≈ 97, 90% CI 81-113), the Crawford & Allan (2001) demographic
        predictor (FSIQ ≈ 99, 90% CI 84-114), and the OPIE-4 prorated FSIQ
        (FSIQ ≈ 113, 90% CI 98-128). All of these models converged on the
        prediction that Mrs Doe's ability was likely to fall within the
        average range."
     The convergence sentence is only emitted when every estimate maps to
     the SAME descriptor band (strict). */
  function rwAutoPremorbidSection(prem){
    if (!prem || !prem.hasInputs) return '';
    if (!prem.estimates || !prem.estimates.length){
      // Inputs entered but no estimates computed yet — minimal one-liner.
      const subj = rwAutoSubjectFor(0);
      return '<h3>' + RWAUTO_DOMAIN_HEADING.premorbid + '</h3>' +
        '<p>' + subj.possessive + ' baseline cognitive ability was estimated using the available demographic and behavioural predictors.</p>';
    }
    const ci = rwAutoReadCi();
    const subj0 = rwAutoSubjectFor(0);
    const possLowerName = rwAutoFullPossessive(); // for second-sentence convergence ("Mrs Doe's ability...")
    // Build "(FSIQ ≈ X, 90% CI X-X)" parenthetical per model
    const modelClauses = prem.estimates.map(e => {
      const ciText = (e.lo != null && e.hi != null)
        ? ', ' + ci.label + ' CI ' + e.lo + '-' + e.hi
        : '';
      return e.label + ' (FSIQ ≈ ' + e.fsiq + ciText + ')';
    });
    // Join with "as well as" for the second clause, then commas + "and" for the rest
    let modelsText;
    if (modelClauses.length === 1){
      modelsText = modelClauses[0];
    } else if (modelClauses.length === 2){
      modelsText = modelClauses[0] + ', as well as ' + modelClauses[1];
    } else {
      modelsText = modelClauses[0] + ', as well as ' + modelClauses.slice(1, -1).join(', ') + ', and ' + modelClauses[modelClauses.length - 1];
    }

    const sent1 = subj0.possessive + ' baseline cognitive ability was estimated using ' + modelsText + '.';

    // Strict convergence check: do all estimates produce the SAME descriptor?
    // Only emit a convergence sentence when there are 2+ models. Pluralisation
    // adapts:  ≥3 models → "All of these models...";  exactly 2 → "Both of these models...".
    const descs = prem.estimates.map(e => rwAutoDescriptor(e.fsiq).toLowerCase()).filter(Boolean);
    const allSame = descs.length === prem.estimates.length && descs.every(d => d === descs[0]);
    let sent2 = '';
    if (allSame && descs[0] && prem.estimates.length >= 2){
      const lead = prem.estimates.length === 2 ? 'Both of these models' : 'All of these models';
      sent2 = ' ' + lead + ' converged on the prediction that ' + possLowerName + ' ability was likely to fall within the ' + descs[0] + ' range.';
    }

    return '<h3>' + RWAUTO_DOMAIN_HEADING.premorbid + '</h3><p>' + sent1 + sent2 + '</p>';
  }

  /* ----- Composer ----- */

  function rwAutoCompose(){
    const battery = rwAutoCollectBattery();
    const changes = rwAutoCollectChange();
    const prem    = rwAutoCollectPremorbid();

    // Group battery by domain
    const byDomain = {};
    battery.forEach(r => {
      const d = r.domain || 'intellectual';
      (byDomain[d] = byDomain[d] || []).push(r);
    });
    const changesByDomain = {};
    changes.forEach(c => {
      const d = c.domain;
      (changesByDomain[d] = changesByDomain[d] || []).push(c);
    });

    const blocks = [];
    let any = false;

    // 1. Validity
    if (byDomain.validity && byDomain.validity.length){
      blocks.push(rwAutoValiditySection(byDomain.validity));
      any = true;
    }
    // 2. Premorbid
    if (prem && prem.hasInputs){
      blocks.push(rwAutoPremorbidSection(prem));
      any = true;
    }

    // 3. Cognitive domains (in convention order, skipping validity/premorbid).
    //    Within a domain, rows are grouped by BRAND (e.g. WAIS-IV, WMS-IV).
    //    Each brand becomes its own paragraph:
    //      Para 1: "[Subject's] [domain] was assessed using X, Y from the [brand1]."
    //               + per-test result sentences for brand1 rows.
    //      Para 2: "She also undertook A, B from the [brand2]."
    //               + per-test result sentences for brand2 rows.
    //    Sentence index continues across brands so name/pronoun cycle stays consistent.
    RWAUTO_DOMAIN_ORDER.forEach(d => {
      if (d === 'validity' || d === 'premorbid') return;
      const rows = byDomain[d] || [];
      const ch   = changesByDomain[d] || [];
      if (!rows.length && !ch.length) return;
      any = true;
      const heading = '<h3>' + RWAUTO_DOMAIN_HEADING[d] + '</h3>';

      // Group rows by brand, preserving first-seen order
      const brandOrder = [];
      const byBrand = {};
      rows.forEach(r => {
        const b = r.brand || 'this assessment';
        if (!byBrand[b]){ byBrand[b] = []; brandOrder.push(b); }
        byBrand[b].push(r);
      });

      const paragraphs = [];
      let idx = 0;
      brandOrder.forEach((brand, bi) => {
        const brandRows = byBrand[brand];
        const opener = (bi === 0)
          ? rwAutoBrandOpenerFirst(d, brand, brandRows, idx)
          : rwAutoBrandOpenerNext(d, brand, brandRows, idx);
        idx++;
        // CVLT verbal-learning gets its own structured narrative (learning
        // curve → short delay → long delay → recognition). Other brands use
        // the generic per-test sentence builder.
        const useCvlt = (d === 'verbal_learning_memory') && rwAutoIsCvltBrand(brand);
        const body = useCvlt
          ? rwAutoCvltNarrative(brandRows, idx)
          : rwAutoBatterySentences(brandRows, idx);
        idx = body.nextIdx;
        paragraphs.push('<p>' + [opener].concat(body.sentences).join(' ') + '</p>');

        // Auxiliary measures (CVLT only) — render in their own paragraph
        // after the main narrative. Uses the standard per-test sentences.
        if (body.auxiliary && body.auxiliary.length){
          const aux = rwAutoBatterySentences(body.auxiliary, idx);
          idx = aux.nextIdx;
          if (aux.sentences.length){
            paragraphs.push('<p>' + aux.sentences.join(' ') + '</p>');
          }
        }
        // Composite measures (CVLT only) — Standard Score scale (T1-5 Correct,
        // Total Recall Correct, Delayed Recall Correct). Rule 18: rendered as
        // a third paragraph after the auxiliary paragraph.
        if (body.composites && body.composites.length){
          const comp = rwAutoBatterySentences(body.composites, idx);
          idx = comp.nextIdx;
          if (comp.sentences.length){
            paragraphs.push('<p>' + comp.sentences.join(' ') + '</p>');
          }
        }
      });

      // Change rows: rule 19 — one paragraph per (battery × method), each with
      // its own opener and a fresh subject cycle. Appended at end of domain.
      if (ch.length){
        rwAutoChangeParagraphs(ch).forEach(p => paragraphs.push(p));
      }

      blocks.push(heading + paragraphs.join(''));
    });

    if (!any){
      return '<div class="rwauto-empty">' +
        '<div class="rwauto-empty-title">No scores entered yet</div>' +
        '<div class="rwauto-empty-sub">Enter scores in <strong>Neuropsych Tables</strong>, <strong>Premorbid Estimate</strong>, or <strong>Change Analysis</strong>. The descriptive narrative will build itself here.</div>' +
      '</div>';
    }

    return blocks.join('');
  }

  /* ----- Render + wire ----- */

  let rwAutoBound = false;
  let rwAutoRenderQueued = false;
  function rwAutoRender(){
    const paper = document.getElementById('rwauto-paper');
    if (!paper) return;
    paper.innerHTML = rwAutoCompose();
  }
  function rwAutoQueueRender(){
    if (rwAutoRenderQueued) return;
    rwAutoRenderQueued = true;
    setTimeout(() => { rwAutoRenderQueued = false; rwAutoRender(); }, 200);
  }

  function rwAutoBindToolbar(){
    const refSel  = document.getElementById('rwauto-reference');
    const descSel = document.getElementById('rwauto-descriptor');
    const legacyRef  = document.getElementById('rw-reference');
    const legacyDesc = document.getElementById('rw-descriptor-system');
    if (!refSel || !descSel) return;

    if (legacyRef){
      refSel.innerHTML = legacyRef.innerHTML;
      refSel.value = legacyRef.value;
      refSel.addEventListener('change', () => {
        legacyRef.value = refSel.value;
        legacyRef.dispatchEvent(new Event('change', { bubbles: true }));
        rwAutoQueueRender();
      });
      legacyRef.addEventListener('change', () => {
        refSel.value = legacyRef.value;
        rwAutoQueueRender();
      });
    }
    if (legacyDesc){
      descSel.innerHTML = legacyDesc.innerHTML;
      descSel.value = legacyDesc.value;
      descSel.addEventListener('change', () => {
        legacyDesc.value = descSel.value;
        legacyDesc.dispatchEvent(new Event('change', { bubbles: true }));
        rwAutoQueueRender();
      });
      legacyDesc.addEventListener('change', () => {
        descSel.value = legacyDesc.value;
        rwAutoQueueRender();
      });
    }

    // Bracket toggles (Standard score / Percentile)
    const stdCb = document.getElementById('rwauto-show-std');
    const pctCb = document.getElementById('rwauto-show-pct');
    if (stdCb){
      stdCb.checked = rwAutoShowStd;
      stdCb.addEventListener('change', () => {
        rwAutoShowStd = stdCb.checked;
        rwAutoSaveToggles();
        rwAutoQueueRender();
      });
    }
    if (pctCb){
      pctCb.checked = rwAutoShowPct;
      pctCb.addEventListener('change', () => {
        rwAutoShowPct = pctCb.checked;
        rwAutoSaveToggles();
        rwAutoQueueRender();
      });
    }
  }

  function rwAutoBindCopy(){
    const btn = document.getElementById('rwauto-copy');
    const paper = document.getElementById('rwauto-paper');
    if (!btn || !paper) return;
    btn.addEventListener('click', () => {
      const html = paper.innerHTML;
      const text = paper.innerText;
      try {
        if (navigator.clipboard && window.ClipboardItem){
          const item = new ClipboardItem({
            'text/html': new Blob([html], { type: 'text/html' }),
            'text/plain': new Blob([text], { type: 'text/plain' })
          });
          navigator.clipboard.write([item]).then(() => {
            if (typeof showToast === 'function') showToast('✓ Report copied');
          });
        } else if (navigator.clipboard){
          navigator.clipboard.writeText(text).then(() => {
            if (typeof showToast === 'function') showToast('✓ Report copied');
          });
        }
      } catch(e){
        if (typeof showToast === 'function') showToast('Copy failed - try selecting manually', true);
      }
    });
  }

  function rwAutoInit(){
    if (!document.getElementById('rwauto-shell')) return;
    rwAutoLoadToggles();
    rwAutoBindToolbar();
    rwAutoBindCopy();
    rwAutoRender();

    if (rwAutoBound) return;
    rwAutoBound = true;

    // Re-render on any score change anywhere on the site
    document.addEventListener('input', e => {
      if (!e.target || !e.target.matches) return;
      if (e.target.matches('#battery input, #battery select, #sdi input, #sdi select, #rci-basic input, #rci-basic select, #rci-practice input, #rci-practice select, #rci-srb input, #rci-srb select, #rci-crawford input, #rci-crawford select, #premorbid input, #premorbid select')){
        rwAutoQueueRender();
      }
    });
    document.addEventListener('change', e => {
      if (!e.target || !e.target.matches) return;
      if (e.target.matches('#battery input, #battery select, #sdi input, #sdi select, #rci-basic input, #rci-basic select, #rci-practice input, #rci-practice select, #rci-srb input, #rci-srb select, #rci-crawford input, #rci-crawford select, #premorbid input, #premorbid select')){
        rwAutoQueueRender();
      }
    });
    // Re-render when navigating back to Report Writer
    document.querySelectorAll('[data-target="report-writer"]').forEach(btn => {
      btn.addEventListener('click', () => setTimeout(rwAutoRender, 80));
    });
  }
  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', rwAutoInit);
  } else {
    rwAutoInit();
  }

})();
