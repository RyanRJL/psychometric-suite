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

Two stylesheets load in order: `styles.css` (9.1k lines, ~1600 `!important`) then
`design-system.css` (4.1k lines, ~580 `!important`), which is a later restyling layer
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

### A large deletion will take its neighbours with it

`app.js` ends feature code and top-level init code in the same flat file, with no
marker between them. Removing the Report Writer as one contiguous 2306-line block
also removed the app's entire init tail, because that tail happened to start on
the line after `setupReportWriter();`. Exactly one line of the deletion belonged
to the feature.

The result: a dozen functions still defined, none of them called. No error, no
console warning, and all 102 numeric checks still green — because none of them
asked whether the app boots. Every calculator loaded with no example rows, the
premorbid page never initialised, autofill was never wired, and "Clear all
tables" was a button with no handler.

**When deleting a block, read the last twenty lines of it as carefully as the
first, and diff the list of top-level statements before and after.** `check.js`
section 16 now guards this, but it only knows about the calls listed in
`INIT_CALLS` — add to that list when you add init code.

It bit a second time, in a subtler form. The same deletion removed
`refreshReportWriterOptions()` but left the **call** to it at the end of
`refreshAll()` — and `refreshAll()` is itself a top-level init statement, so it
threw during boot and every top-level statement after it stopped running. That
is 2,350 lines: the whole Working Report bundle, "Clear all tables", the Score
Converter's Distribution tab, and the top-bar nav sync. No console error was
surfaced, and all 107 checks passed, because section 16 verifies that init
functions are *called* — it cannot tell whether the call *succeeds*.

`check.js` section 17 closes that: every bare call to a name in the project's
own naming families must resolve to something bound in the four shipped
scripts. **If a `typeof SomeConst` throws `ReferenceError` instead of returning
`"undefined"`, that const's declaration never ran** — a top-level statement
above it threw. That probe is the fastest way to find the boundary.

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
| `app.js` | 7.1k | Every calculator. 206 top-level functions, 11 section banners. |
| `design-system.js` | 1.1k | Page titles, microcopy, FOUC handling, the inline control bars. |
| `app-effectsize-page.js` | 0.8k | The effect-size calculator, self-contained. |
| `app-viz-page.js` | 0.3k | Score Charts page — one SVG chart per test family, one row per trial/subtest, self-contained. Draws only what the table's own functions return; has no settings of its own. |
| `index.html` | 6.3k | All pages in one document; sections shown/hidden by nav. |
| `styles.css` | 9.1k | Original stylesheet. |
| `design-system.css` | 4.1k | Later restyling layer that overrides the above. |
| `service-worker.js` | 79 | Cache-first PWA worker. |
| `tools/check.js` | — | Headless numeric regression checks. |

The **Report Writer** — a page that generated descriptive narrative prose from
entered scores, ~2.3k lines in `app.js` plus ~2.0k in `design-system.js` — was
removed in July 2026. If you find a stray `rwAuto*` reference or `.rw-` CSS rule,
it is a leftover and can go.

**`REPORT_TEST_CATALOG`** (`data.js`) is what lets the Working Report's "Merge by
battery" button work. The report auto-splits a Score Tables table into one item per
test family, and the merge puts same-instrument items back together as one table with
a labelled sub-section each. The catalog was referenced by the merge code but had never
been defined, so `catalogBatteryFor()` always returned `null` and the button did nothing
— silently, with no error. Two invariants keep it honest, both in `check.js` §21: every
catalog `name` must appear in `TEST_FAMILY_PATTERNS` (`app.js`), because that is the only
thing `detectTestFamily()` can return; and every `normDB` group must be claimed by
exactly one battery. ToPF and OPIE are deliberately **not** catalogued — their tables
name WAIS-IV/WMS-IV as the *predicted* criterion, and merging one into an achieved-score
table would put predicted and obtained scores in the same column. `mergeableBattery()`
refuses `pre-*` sources outright for the same reason.

Do not confuse the Report Writer with the **Working Report bundle**
(`app.js`, "WORKING REPORT BUNDLE v2"), which is live and collects APA tables from
every calculator into a drawer.

### Navigating `app.js`

Section banners (`/* ====`) mark the boundaries:

| Line | Section |
|---|---|
| ~1125 | Battery table (Score Tables page) |
| ~1923 | APA table notes — single source of truth |
| ~2317 | SDI |
| ~2673 | RCI calculators (Basic / Practice / SRB / Crawford) |
| ~3213 | Per-method autofill from the normative database |
| ~3401 | Custom tests database management |
| ~3683 | Premorbid estimation |
| ~4811 | Auth overlay |
| ~4957 | Top-bar navigation bucket sync |
| ~5029 | Score Converter view-mode tabs |
| ~5048 | Working report bundle |

Line numbers drift with every edit — treat them as a starting point and grep to
confirm. Function-name families make grep effective: `calc*` computes, `render*`
draws, `get*` derives, `fmt*` formats, and page prefixes are `pre*` (premorbid),
`rci*`, `sdi*`, `bat*` (battery), `opie*`. To find any calculator, grep the prefix.

Statistical primitives sit at the top: `erf` (~7), `normCDF` (~16), `tInv` (~73).
General formatters follow at ~132–152 (`fmt`, `fmtPct`, `fmtP`). Note that
`fmtIntOrDash` and `fmtPctBr` are **not** there — they are premorbid-local, at ~3709.

Pages hold state in module-level objects: `preState`, `rciState`, `batteryRows`,
`sdiRows`. Note that tables are sometimes **updated in place** rather than re-rendered,
specifically so a focused input does not lose focus mid-typing — see `refreshOpieDerived`
for the pattern. Follow it when touching a table that contains inputs.

---

## `normDB`

The reliability database: 113 groups, 641 entries, seven instrument families (D-KEFS,
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
- `metric` — present only as `'raw'`, on 63 entries. See below.

### `metric:'raw'` — the entries that carry no scale

Most entries are on a standardised metric, but `normDB` does not record which one.
Everything except these 63 entries has its score type **inferred from the normative
mean** by `inferScoreTypeForSubtest` (nearest of z≈0, scaled≈10, T≈50, standard≈100).

That heuristic has no raw category, so for a raw score it silently returned whichever
standard metric the raw mean sat nearest. `metric:'raw'` exists to stop it guessing:

| Group | Entries |
|---|---|
| `CVLT-C Subtests (Raw Scores) · Age 8 / 12 / 16` | 39 (all of them) |
| `RBANS Subtests · Ages 12-19 / Ages 20-89` | 4 of 12 per band |

CVLT-C declares it in the family name, so those three groups are raw throughout.

**RBANS is a genuinely MIXED family — do not tag it wholesale.** Only these four are
raw: `Line Orientation`, `Picture Naming`, `List Recognition`, `List Recall`. The other
eight are scaled scores and must keep producing a percentile. Tagging the whole family
was an early mistake here and it silently killed the derived scores on eight good
measures per band.

What settles it without the manual: **the eight scaled measures cluster between 9.6 and
11.6 despite raw maxima of 12, 16, 20, 24, 40 and 89.** Raw scores on scales that
different cannot all land near 10 — that only happens on a shared metric. `Coding` is
the clearest case: raw max 89, mean 10.8. The four raw ones break the cluster —
`Line Orientation` sits at 80% of a 20-point scale, `List Recognition` at 98% with
SD 0.8, `Picture Naming` at 98% with SD 0.4, and `List Recall` at 6.2, which no
normative sample would average on a scaled metric.

Two invariants in `check.js` §18 hold that shape: any entry with `sd1 < 1` must be
tagged (**no standardised metric in this app has an SD below 1**), and every *untagged*
RBANS subtest must have `m1` within 3 of 10 and `sd1` between 1.5 and 4.5.

### `reportedAs` — what the clinician types, which is not always what `m1`/`sd1` are

`metric` describes the **stored statistics**. `reportedAs` describes the **entry
metric**. For most of the database they are the same thing and only `metric` exists.
CVLT-C is where they part company: its norms here are raw, which is what Change
Analysis needs, but nobody transcribes a raw CVLT-C score onto a results table — they
copy the standardised score off the record form.

| Measure | `metric` | `reportedAs` | Source |
|---|---|---|---|
| `List A Trials 1-5 Total` | `raw` | `t` | Manual Table A.1 |
| every other CVLT-C index | `raw` | `z` | Manual Table A.2 |

**The manual's wording is a trap.** Table A.2 is headed "*Standard Score* Equivalents",
but its scale runs −5.0 to +5.0 in 0.5 steps — that is a **z-score**, not this app's
`standard` (M 100, SD 15). A CVLT-C "standard score" of −1.0 is the 16th percentile;
read as M 100 / SD 15 it becomes z = −6.7. `check.js` §18 asserts that no CVLT-C entry
is ever typed `standard`.

`inferScoreTypeForSubtest` returns `reportedAs` first, then `metric`, then the mean
heuristic. Score Tables and the SD Index ask the entry question, so they get
`reportedAs`; Change Analysis reads `m1`/`sd1` directly and is unaffected. Because
these standardised scores are already age-corrected, the age band chosen on autofill
does not change a Score Tables percentile — it only matters for Change Analysis.

### `rInternal` — the one place internal consistency beats test–retest

`r` is a **stability** coefficient (same people, retested). `rInternal` is a
**consistency** coefficient (one sitting, split in half). Different questions, not
interchangeable.

It exists because the app's default is test–retest for every measure, and that default
is right almost everywhere (see the CI paragraphs in Methods & References). The bar for
setting it aside is deliberately high: **the publisher must both report an
internal-consistency coefficient AND derive its own published intervals from it.**

Three manuals clear it, and `check.js` §24 pins the roster so a fourth cannot appear
undocumented:

| Source | Entries | Fields |
|---|---|---|
| CVLT-C `List A Trials 1-5 Total`, ages 8 / 12 / 16 | 3 | `rInternal` only |
| D-KEFS Advanced — CWIT, Tower, Social Sorting, Risk-Reward | 26 | `rInternal` + `rInternalByAge` |
| D-KEFS (original) — VFT, Sorting, Word Context, Tower, Proverb, TMT composite | 11 | `rInternalByAge` **only** |

That roster check keys on **either** field. D-KEFS original carries no `rInternal` at
all, so a check testing `rInternal` alone would let 11 entries change the basis of a
printed interval without ever being counted.

- CVLT-C Manual Table 6.5 gives odd/even split-half by age — .87 / .89 / .84 for the
  three bands here — alongside an SEM for each.
- Those SEMs are on the **T metric**: every one of the 12 age bands equals
  `10 × √(1 − rxx)`. The manual truncates rather than rounds (3.6055 → 3.60), so assert
  the printed value lies in `[derived − 0.01, derived]`, not equality after
  `Math.floor` — `√0.09` lands a hair under 0.3 and would truncate an exact 3.00 to 2.99.
- p. 87 states `SEM = SD √(1 − rxx)` and says to add and subtract from *the child's
  standard T score*. Its worked example — an 8-year-old at T 45, 95% → 38–52, 90% →
  39–51 — reproduces exactly, and §24 drives the real renderer to prove it.

**`rInternal` must never reach Change Analysis or the RCI pages.** With the retest `r`,
`SD√(2(1−r))` *is* the spread of the difference distribution, which is the quantity
reliable change is measured against; an internal-consistency coefficient strips out real
between-session fluctuation, narrows the interval and overcalls change. On SRB and
Crawford it would also move a fitted regression slope. §24 asserts that none of
`calcBasicRow`, `calcPracticeRow`, `calcSrbRow`, `calcCrawfordRow` or `rciEffectiveR`
mentions the field, and that `getBatteryRowReliability` does — so it cannot become inert
either.

Note this is why CVLT-C keeps **both** coefficients on the same entry. Deleting `r`
because "`rInternal` is better" would break Change Analysis for that measure.

#### `rInternalByAge` — when the published coefficient changes with age

D-KEFS Advanced publishes reliability **by normative age band**, so one number per
measure would be a compromise. 26 entries carry a lookup instead:

```js
rInternal: 0.86,                                    // published Average — the fallback
rInternalByAge: { 8:0.69, 9:0.75, …, 80:0.95 },     // Table 3.4, keyed by band LOWER bound
rInternalAgeMax: 90
```

Keys are the band's **lower bound**; `rInternalForAge` takes the greatest key ≤ age. So
17 reads the 16-18 band and 45 reads 40-49. Same shape as `baseRates` — a lookup living
on an entry rather than a separate group.

The age is read through `batteryPatientAge()`, which since 616a70a delegates to
`patientAge()` — **one shared value** behind `#bat-age` on Score Tables and `#pre-age` on
Premorbid, pinned by `check.js` §26. It is **optional**. Blank, or outside
`rInternalAgeMax`, falls back to `rInternal` where one is published — the manual's own
all-ages figure, so both paths are citable — and to the retest `r` where none is, which is
the D-KEFS original case below. **Never make it required**: a blank age silently emptying
the CI column would read as the app being broken.

Three constraints, all in `check.js` §25:

- **`rInternalByAge` only on `· All Ages` groups.** The banded D-KEFS Advanced groups
  (8-18 / 19-59 / 60-90) hold the *retest study's* bands, which are not the normative
  bands the reliability table uses — the same mismatch as `baseRates` vs `WAIS-IV Process
  Scores`. It works because Score Tables collapses these families to `All Ages`
  (`buildFamilyListHtml`, no `baseRates` so no `familyScoredByAgeBand`), while Change
  Analysis and the SD Index see only the banded ones, where the retest `r` is wanted.
- **TMT and VFT get none of this.** Their manual states split-half and alpha are not
  accurate for **speeded** measures and uses stability coefficients for those two tests.
  Their Table 3.4 rows *are* the retest coefficients, and match `normDB`'s stored `r`
  exactly at ages 8-18 (.47/.47/.62/.77/.76/.40/.65/.75/.74).
- **The APA note names the age**, but only when a measure in the table actually read one
  (`batteryRowUsesAgeBand`). One field silently driving every interval has to be visible
  on the exported table, and claiming an age on a table where nothing used it would be
  its own misstatement.

#### D-KEFS (original) — `rInternalByAge` with no `rInternal`

The D-KEFS Technical Manual, Chapter 2, states the rule outright (p. 18): "For some
D-KEFS measures, the internal consistency coefficients were used; for other measures the
test–retest correlations were employed." It prints `SEM = SD √(1 − rxx)` with "The
standard deviation unit is **3** for all D-KEFS scaled scores", and
`CI = observed ± z(SEM)` at 90% and 95%.

So it runs **two regimes**, chosen per test: internal consistency tabulated *by age*, or
the retest `r` on the *total sample* where no coefficient table exists — the by-band
retest samples being "too small".

Which regime feeds the SEM is settled arithmetically, not by reading the prose:
`3 × √(1 − rxx)` on the internal-consistency coefficient reproduces **190 of the 216**
published SEM cells exactly, and all but one measure to within 0.006 (the manual computes
from unrounded coefficients and prints both to 2 dp). The retest `r` reproduces **2 of
200**. `check.js` §27 pins all 168 cells for the measures actually used, and asserts the
falsification too.

**No `rInternal` fallback exists, and none may be invented.** Unlike Advanced, not one of
the eight tables prints an all-ages average. A blank or out-of-range age therefore falls
through `rInternalForAge` to the retest `r` — which is precisely the manual's own second
regime, confirmed by Design Fluency Table 2.8, whose three All Ages SEMs all reproduce as
`3 × √(1 − r)`. Both paths are the publisher's own figures, so the column stays citable
either way. Averaging the bands to manufacture an `rInternal` would be inventing a number.

Four exclusions, four different reasons — each an easy silent mistake, each pinned:

- **Colour-Word Interference, all four.** Its only coefficient table (2.9) is for the
  *Combined Colour Naming + Word Reading composite*, which is not a measure in `normDB`.
  Attaching it to Colour Naming or Word Reading would give a condition the composite's
  reliability.
- **Design Fluency, all three.** "Item interdependence precluded the use of internal
  consistency procedures."
- **Trail Making, five of six.** Table 2.1 is the composite alone, so only
  `Combined Number + Letter` carries it.
- **Twenty Questions, both — held deliberately.** Table 2.17's SEM column does not follow
  from Table 2.15's coefficients: 15 of 16 bands disagree, by up to .094 and in both
  directions. At 40–49, rxx .75 implies SEM 1.50 where the manual prints 1.76. No other
  coefficient column, row shift, constant SD or uncorrected split-half fits either, so one
  of the two tables is misprinted and the source cannot say which. The bar is therefore not
  demonstrably met. `Total Weighted Achievement` reproduces 16/16 from its own coefficients
  and is held only because it shares the suspect table — revisit both together if another
  printing settles it.

Word Proverb is banded 16–19 upward, the Proverb Test not being administered under 16, so
a child's age finds no band and correctly falls back.

**The two D-KEFS manuals reach opposite conclusions on the same two test names.** Advanced
excludes Trail Making and Verbal Fluency for *speededness*; the original publishes internal
consistency for all four Verbal Fluency measures and the TMT composite, and excludes Design
Fluency for *item interdependence* instead. Do not "harmonise" them — §27 asserts the
inversion in both directions.

Reliability method is a **per-manual question, not a per-app policy.** Four sources, four
different answers: D-KEFS Advanced rejects internal consistency for *speededness*; D-KEFS
original accepts it for those same tests and rejects Design Fluency for *item
interdependence*; CVLT-3 and CVLT-C reject it for *item interdependence* on a word-list
task; CVLT-C nonetheless publishes split-half for its one non-interdependent score. Read
the manual before assuming any of them generalises.

### `baseRates` — measures scored by lookup rather than by conversion

**WAIS-IV Longest Span (Process)** is the first family scored from a published
table instead of a formula, and the first with **no retest data at all**. Source:
WAIS-IV Administration and Scoring Manual (GB), Tables C.4 and C.5 — 14 age bands
× Longest Digit Span Forward / Backward / Sequence, plus Longest Letter-Number
Sequence for the nine bands 16-17 to 65-69, which is all the manual publishes.

`baseRates[x]` is **the percentage scoring x or higher**. That is derived, not
assumed: for a whole-number score E[X] = Σ P(X ≥ x), and reconstructing the mean
that way reproduces every printed `m1` across all 51 measure-bands to within
0.055. `check.js` §22 re-runs it, so flipping the table to "or lower" fails loudly.

**A base-rate measure ALWAYS reports its base rate**, whatever else is in the
table, and **its section carries its own column heading**. The main heading never
changes.

An earlier version switched the whole column between the two quantities depending
on the other rows. That was wrong twice over: a longest-span value silently
changed meaning when an unrelated subtest was entered, and — because the switch
skipped seeded `isExample` rows while still rendering them — a table could print
"Base rate" above a cell holding a percentile. Both were visible on screen. One
quantity per measure, labelled where it is read.

`batteryGroupIsBaseRate()` drives both the editable table (a `.group-col-header`
row after the group header) and the APA export (a `.apa-group-cols` row, via the
column's `groupLabel`), so screen and export are labelled identically.

The percentile rank — used for the classification, and for any non-base-rate row
— uses the **midpoint convention** for a discrete score: `P(X < v) + ½·P(X = v)`,
so 21 + 8.75 = 29.75 for that span of 6. Not `100 − P(X ≥ v)`, which sounds
simpler but runs 6–16 points below the percentile the same measure's published M
and SD would give, and returns exactly 0 at the bottom of the scale — a
percentile the app treats as impossible (see `fmtPct`).

`fmtBaseRate`, not `fmtPct`, formats a base rate: `fmtPct` clamps into
(0.01, 99.99) to keep percentiles inside the open interval, and a published base
rate legitimately reaches 100.

Three things follow, each with a check:

- **`singleAdministration:true`** marks entries with no `m2`/`sd2`/`r`, because
  none is published. They are filtered out of Change Analysis and the SD Index —
  loading one there gives rows that can never compute, which reads as the app
  being broken rather than as data being absent.
- **The age bands are the normative ones** (16-17, 18-19, …, 85-90) and
  deliberately do **not** match `WAIS-IV Process Scores`, which uses the retest
  study's bands (16-29, 30-54, 55-69, 70-90). They cannot be merged.
- **Score Tables must not collapse these age bands.** Its flat dropdown normally
  shows one entry per family, on the grounds that the band does not change the
  table — true for a metric conversion, false when the band *is* the lookup.
  Longest Digit Span Backward at a span of 4 is the **22nd** percentile at 20-24
  and the **59th** at 85-90. `familyScoredByAgeBand()` exempts them.

Note the table was not adopted for accuracy: computing from the printed M/SD
lands within ~3 percentile points at every span. It was adopted because it is the
**published, citable** figure, which is worth more in a report than a close
approximation.

### `higherIsWorse` — measures where a high score is a bad result

On CVLT-C error measures the standardised score runs with the error count: Table A.2
maps z −1.0 to 0–3 perseverations and z +5.0 to 45 or more. A child at z +2.0 has made
more perseverations than 98% of their age group, and would otherwise have been labelled
"Very Superior".

Flagged on four measures: `Perseverations`, `Free-Recall Intrusions`,
`Cued-Recall Intrusions`, `False Positives`. `Recognition Hits` and `Discriminability`
run the normal way and must **not** be flagged.

The two columns are treated differently, on purpose:

- **Classification describes performance**, so it is computed from the **reflected**
  score. z +2.0 classifies as z −2.0 does. A child with more perseverations than 98% of
  their age group reads "Borderline", which is what a clinician would write. It works
  both ways: `Cued-Recall Intrusions` at z −1.0 — fewer intrusions than average — reads
  "High Average".
- **The percentile is NOT reflected.** 98th is what z +2.0 gives and what anyone
  checking the working against the manual will calculate. A number that cannot be
  reproduced is worth less in a report than one that needs a footnote.

So an error row pairs a **high percentile with a low classification**. The APA note
carries that convention so the pairing does not read as a contradiction. Blanking the
classification was tried first and rejected: the column exists to summarise performance,
and a dash reads as missing data rather than as a deliberate choice.

**No premorbid asterisks** are placed on these rows. The stars mean "this ability falls
short of the premorbid estimate", and an error count is not an ability being predicted,
so there is no comparison to assert.

`toZ`/`fromZ` return `null` for `'raw'` in both directions, so percentile and
classification go blank rather than being invented. **Do not give `'raw'` a fallback** —
falling through to `scaled` is the bug this closes: RBANS `List Recognition` at a raw 19
is the 23rd percentile, but read as a scaled score it printed the 99.9th, "Very
Superior" — 56 standard-score points out, with the sign inverted.

The confidence interval is still shown for raw rows and is valid (`sd1` is in raw units,
so `SD × √(1−r)` is a raw-unit SEM), but takes the z-score treatment — 1 dp, no floor,
no ceiling — because many raw measures start at 0. SDI **index** mode refuses raw rows
(no metric SD to divide by); SDI **raw** mode scores them correctly from the row's own SD.

Beware the near-miss: `D-KEFS Sorting Test · Ages 8-19` has SD ≈ 1.6 against a scaled
unit of 3 and looks raw. It is not — the Ages 20-49 band of the same test shows SD ≈ 3.2,
so it is a scaled score in a range-restricted retest sample. Check a sibling age band
before tagging anything new.

**`r` and `rCorrected` describe different populations, and this matters.** Every SEM in
this app is `SD × √(1 − reliability)`, which is only a valid error variance when both
terms describe the same group. `rCorrected` belongs with a normative SD (15 / 3 / 10 / 1);
`r` belongs with `sd1`/`sd2`. Pairing `sd1` with `rCorrected` understates the error
variance by ~20% and is what the RCI pages used to do. All four RCI methods now default
to raw `r`, and the reliability cell renders whichever field is actually in force, so
what is displayed is what is used. `tools/check.js` guards both.

The **Score Tables CI** made the same mistake for longer, and independently: it took
`rCorrected` but multiplied by `sd1`. It now takes the SD from **the metric the score is
displayed in**, reaching for `sd1` only when the row is shown raw. The rule is not a
convention picked here — the CVLT-3 manual publishes an SEM column in Tables 3.4/3.5,
and `normativeSD × √(1 − rCorrected)` reproduces **all 38** of them exactly, which is
what `check.js` §23 pins. The old pairing was wrong in both directions and by up to 6.3
standard-score points at each end (RBANS Form C Delayed Memory, whose retest sample is
impaired: `sd1` 18.7 against a normative 15).

Note `r` does not play the same role in all four methods. In Jacobson & Truax and
Iverson it is purely an error term. In McSweeney and Crawford it is a **fitted
regression slope** (`slope = r × sd2/sd1`), so substituting a population-corrected
value there changes the predicted score, not just the interval — it yields a
regression line that was never fitted to anything. The corrected-`r` option is
therefore defensible only on the Basic page, and then only if the SD is also set to
the normative value.

**Change outcomes report significance only** — "Reliable change" / "No reliable
change", never "improvement" or "decline". `normDB` holds 119 higher-is-worse
measures (intrusions, perseverations, errors, false positives, repetitions) and no
score-direction flag, so inferring valence from the sign of the statistic asserts the
wrong clinical conclusion for all of them. The signed statistic is displayed
alongside, so direction stays visible without the app interpreting it.

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
1.645 / 1.960 / 2.576 — see `PREMORBID_CI_Z` (~1356). Do not conflate the two.

**APA notes** live in the `APA_NOTES` object (~line 1934, under its banner at ~1923) and
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

`node tools/check.js` runs 193 headless checks: statistical primitives, score-conversion
round trips, `normDB` structural integrity, WAIS-IV values pinned to Technical Manual
Table 4.5, OPIE-4 coefficients pinned to Table eA5.8, worked OPIE predictions, reliable-
change thresholds and direction-neutral outcome labels, base-rate reconstruction and
monotonicity, percentile-tail clamping, the effect-size calculator, Score Tables
confidence intervals, documentation contracts, wiring (§16–17) and the raw-score
metric (§18).

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
