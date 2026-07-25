# Psychometric Suite

A clinical neuropsychology calculator suite used to produce real patient reports,
including medico-legal ones. **Numbers that appear on screen end up in documents that
get scrutinised.** Accuracy outranks speed, features and tidiness.

Static HTML/CSS/JS. **No build step, no npm, no framework, no bundler.** Files are
opened directly in a browser or served locally. That is deliberate — the tool must work
by double-clicking a single file — so do not introduce a toolchain.

---

## Shipping a change

Run these in order. Skipping step 2 or 3 means you are testing stale code.

1. **Edit** the source files (`app.js`, `data.js`, `index.html`, `styles.css`, …).
2. **Regenerate the bundle** — `Psychometric_Assistant.html` is committed and must be
   rebuilt or it goes out of sync with source:
   ```bash
   ./bundle.ps1
   ```
3. **Bust the caches** — two separate mechanisms, both needed:
   - the `?v=…` query strings on the `<script>`/`<link>` tags in `index.html`
   - `CACHE_VERSION` in `service-worker.js`
4. **Verify the maths** — must pass before committing:
   ```bash
   node tools/check.js
   ```
5. **Verify in the browser** if the change is visible. See the trap below.
6. **Commit.** Include `Psychometric_Assistant.html` in the same commit as its sources.

`bundle.ps1` inlines every local CSS and JS file into one self-contained
`Psychometric_Assistant.html`, stripping the `?v=` strings as it goes. Google Fonts stay
as external links.

### Dev server

`.claude/launch.json` defines a PowerShell static server (`.claude/serve.ps1`) on port
3131. If another session already holds that port, add a second entry on a different port
rather than killing theirs.

### What is not committed

`.gitignore` excludes `.claude/`, `*.ps1`, `fix_encoding*`, `temp_opie_*`, `*.ipynb`.
So `bundle.ps1` and everything under `.claude/` — including `launch.json` and
`settings.json` — are **local only**. Do not put anything load-bearing there.

---

## Traps that will cost you an hour

### The service worker serves stale files

Bumping `CACHE_VERSION` alone does **not** refresh an already-open tab. The old worker
stays in control and keeps serving its cached `index.html`. Symptom: your edit is
provably in the file and provably served by the dev server, but the page does not show
it, and you start doubting a correct change.

Before testing in the browser, always:

```js
for (const r of await navigator.serviceWorker.getRegistrations()) await r.unregister();
for (const k of await caches.keys()) await caches.delete(k);
```

then navigate again. This bites roughly every second session.

### The CSS cascade will silently eat your rule

Two stylesheets load in order: `styles.css` (9.5k lines, ~1700 `!important`) then
`design-system.css` (4.1k lines, ~590 `!important`), which is a later restyling layer
that deliberately overrides the first. About **70 class selectors are defined in both
files**, and `.input-table` alone is declared ~86 times within `styles.css`.

The consequence: a correct, more-specific rule can still lose, because page-scoped
selectors like `#premorbid #pre-opiepredict .info-box` (specificity 2-1-0, with
`!important`) outrank `.info-box.is-caution` (0-2-0).

**Rule: when adding styled UI, create a new standalone class. Never add a modifier to a
class that `design-system.css` or a page-scoped block already touches.** A standalone
name sidesteps the whole cascade and needs no `!important` of its own. `.caution-box` in
`styles.css` exists precisely because `.info-box.is-caution` was tried first and lost —
its comment records the reasoning.

If a style mysteriously does not apply, enumerate the competing rules rather than
guessing:

```js
[...document.styleSheets].flatMap(s => [...s.cssRules])
  .filter(r => r.selectorText && r.selectorText.includes('your-class'))
  .map(r => r.cssText)
```

### On-screen text is a contract

Every calculator **prints its own formula** in a `<details class="formula-disclosure">`
block and **cites its method** in the APA note, section header, references list, method
picker and tooltip. Those are load-bearing claims, not decoration.

If you change a formula you must change the text — and if the text names a published
method, changing the formula means the page now implements a *different* method and must
be re-branded everywhere. Grep the citation before touching any calculation. A recent
change to the practice-adjusted RCI had to be scoped narrowly for exactly this reason:
the page is branded Iverson (2001) in six places.

---

## Layout

| File | Lines | Role |
|---|---|---|
| `data.js` | 1.1k | All constants, coefficients, lookup tables and `normDB`. Plain script, defines globals. Loads first. |
| `app.js` | 9.3k | Every calculator. 291 top-level functions, 12 numbered section banners. |
| `design-system.js` | 3.0k | Page titles, microcopy, FOUC handling, and the **report writer** (62 `rwAuto*` functions that turn scores into narrative sentences). |
| `app-effectsize-page.js` | 0.7k | The effect-size calculator, self-contained. |
| `index.html` | 6.5k | All pages in one document; sections shown/hidden by nav. |
| `styles.css` | 9.5k | Original stylesheet. |
| `design-system.css` | 4.1k | Later restyling layer that overrides the above. |
| `service-worker.js` | 79 | Cache-first PWA worker. |
| `tools/check.js` | — | Headless numeric regression checks. |

### Navigating `app.js`

Section banners (`/* ====`) mark the boundaries:

| Line | Section |
|---|---|
| ~1117 | Battery table |
| ~1872 | APA table notes — single source of truth |
| ~2252 | SDI |
| ~2608 | RCI calculators (Basic / Practice / SRB / Crawford) |
| ~3096 | Per-method autofill from the normative database |
| ~3284 | Custom tests database management |
| ~3566 | Premorbid estimation |
| ~4697 | Report writer — descriptive narrative builder |
| ~7231 | Working report bundle |

Function-name families make grep effective: `calc*` computes, `render*` draws,
`get*` derives, `fmt*` formats, and page prefixes are `pre*` (premorbid), `rci*`,
`sdi*`, `bat*` (battery), `opie*`. To find any calculator, grep the prefix.

Statistical primitives sit at the top: `erf` (~7), `normCDF` (~16), `tInv` (~73).
General formatters follow at ~132–144 (`fmt`, `fmtPct`, `fmtP`). Note that
`fmtIntOrDash` and `fmtPctBr` are **not** there — they are premorbid-local, at ~3595.

Pages hold state in module-level objects: `preState`, `rciState`, `batteryRows`,
`sdiRows`. Note that tables are sometimes **updated in place** rather than re-rendered,
specifically so a focused input does not lose focus mid-typing — see `refreshOpieDerived`
for the pattern. Follow it when touching a table that contains inputs.

---

## `normDB`

The reliability database: 99 groups, 590 entries, seven instrument families (D-KEFS,
CVLT-3, CVLT-C, RBANS, WAIS-IV, WISC-V, WMS-IV).

Group keys are `"<Instrument> <Category> · <Age band>"` — the separator is **U+00B7
MIDDLE DOT**, not a hyphen. Entries look like:

```js
"Full Scale IQ": { m1:99.7, sd1:13.8, m2:104, sd2:15, r:0.95, rCorrected:0.96, n:298 }
```

- `m1`/`sd1`, `m2`/`sd2` — mean and SD at first and second testing
- `r` — the observed test-retest correlation **in that retest sample**
- `rCorrected` — `r` rescaled to the **normative population's** variability
  (Allen & Yen / Magnusson)
- `n` — retest sample size, present only where published

**`r` and `rCorrected` describe different populations, and this matters.** Every SEM in
this app is `SD × √(1 − reliability)`, which is only a valid error variance when both
terms describe the same group. `rCorrected` belongs with a normative SD (15 / 3 / 10 / 1);
`r` belongs with `sd1`/`sd2`. Pairing `sd1` with `rCorrected` understates the error
variance by ~20% and is what the RCI pages used to do. All four RCI methods now default
to raw `r`. `tools/check.js` guards this.

Field coverage: `m1`, `sd1`, `m2`, `sd2` and `r` are present on **all 590** entries;
`rCorrected` on only **233** — D-KEFS, CVLT-C and WISC-V have none at all, so any feature
depending on it must degrade gracefully.

Users can add custom tests; `getMergedDB()` merges those over `normDB`.

---

## Domain conventions

**Sex coding differs between instruments and is a live footgun.**
- ToPF / WAIS-IV demographic equations: Female = 1, Male = 2
- OPIE-4: **Female = 0**, Male = 1

Because OPIE-4 codes Female as 0, defaulting a blank sex field to 0 silently returns the
*female* equation rather than no answer. OPIE-4 therefore gates on sex being present, the
same way it gates on age. Never default a coded categorical to 0 without checking what 0
means.

**Score types** are `standard` (M 100, SD 15), `t` (50/10), `scaled` (10/3), `z` (0/1).
`rowScoreType(row)` resolves per-row overrides against the page default.

**Confidence intervals** round the estimate and the margin separately —
`round(estimate) ± round(z × SEE)` — so bounds stay symmetric about the printed value.
All display sites must use the same convention or the screen and the exported report
disagree; the premorbid Estimates tab, its APA export and the OPIE-4 tab are three
separate sites that must agree.

The user-facing CI selector offers **90% (z 1.645) and 95% (z 1.960) only**. Separately,
premorbid-vs-achieved significance flagging uses a three-tier `*`/`**`/`***` scheme at
1.645 / 1.960 / 2.576 — see `PREMORBID_CI_Z` (~1348). Do not conflate the two.

**APA notes** live in the `APA_NOTES` object (~line 1883, under its banner at ~1872) and
render via `apaNoteHtml`.
That is the single source of truth for the note text under every exported table — change
methodology there, not in the markup.

**OPIE-4 is labelled illustrative-only for UK use.** The coefficients reproduce
Holdnack et al. (2013) Table eA5.8 exactly, but the published equations also carry US
education, ethnicity and region terms that are not applied, so every patient is scored at
the US reference category. Those US education categories have no valid UK equivalent —
the dummies encode how unusual an attainment level is within the US population, not years
of schooling, so they cannot be mapped across. See the long comment above
`OPIE_PRORATED_FSIQ` in `data.js` before touching anything OPIE-related.

**`BASE_RATES` is a parametric normal model, not observed frequencies.** Every one of its
298 cells equals `round(Φ(d / SEE), 4)`. It is labelled as such everywhere it appears.
`OPIE_BASE_RATES`, by contrast, *is* empirical — its values sit on a count grid. If you
ever replace `BASE_RATES` with real published frequencies, update the labels in `app.js`,
`index.html` and `data.js` too.

---

## Verifying calculations

`node tools/check.js` runs 67 headless checks: statistical primitives, score-conversion
round trips, `normDB` structural integrity, WAIS-IV values pinned to Technical Manual
Table 4.5, OPIE-4 coefficients pinned to Table eA5.8, worked OPIE predictions, reliable-
change thresholds, base-rate reconstruction and monotonicity, and documentation contracts.

It loads `data.js` through Node's `vm` module and **re-implements the formulas
independently** rather than importing them from `app.js`. That duplication is
intentional: if the two ever disagree, one has drifted and you need to find out which.
Do not "simplify" it by importing from `app.js`.

Only add checks that are either **invariants** (true by mathematics or data shape) or
**pinned to a cited source**. Never invent an expected number — if you cannot cite it,
derive it in the check and show the arithmetic.

For anything DOM-bound, drive the real page instead:

```js
// after clearing the service worker
calcPremorbid(); calcOpiePredict();
[...document.querySelectorAll('#pre-opiepredict-tbody tr')].map(tr => tr.textContent)
```

---

## Working style

- **Verify before reporting.** A false alarm costs more than a missed issue. Validate
  your test method against a known-good case first — if a calibration test flags an
  equation, run it against a known-correct equation in the same file to prove the method
  works before concluding the code is wrong.
- **Check whether an apparent defect is project-wide** before calling it a bug in one
  place. A finding that "r and SEE are inconsistent" dissolved once every model in the
  app — including pairs taken straight from published manuals — showed the same pattern.
- **State plainly what is proven versus inferred**, and say when a source could not be
  obtained.
- **Commit each completed change**, with the reasoning and the numbers in the message.
  These commits are the audit trail for clinical decisions.
- **Never adjust a published coefficient to make a result look right.** If a published
  model behaves oddly, document it and surface the caveat to the clinician.
