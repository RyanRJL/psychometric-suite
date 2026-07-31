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
| `data.js` | 1.8k | All constants, coefficients, lookup tables and `normDB`. Plain script, defines globals. Loads first. |
| `app.js` | 8.2k | Every calculator. 218 top-level functions, 11 section banners. |
| `design-system.js` | 1.1k | Page titles, microcopy, FOUC handling, the inline control bars. |
| `app-effectsize-page.js` | 0.8k | The effect-size calculator, self-contained. |
| `app-viz-page.js` | 1.4k | Score Charts page — one SVG chart per test family, one row per trial/subtest, self-contained. Draws only what the table's own functions return; has no settings of its own. |
| `index.html` | 6.4k | All pages in one document; sections shown/hidden by nav. |
| `styles.css` | 9.2k | Original stylesheet. |
| `design-system.css` | 4.5k | Later restyling layer that overrides the above. |
| `service-worker.js` | 80 | Cache-first PWA worker. |
| `tools/check.js` | 4.7k | Headless numeric regression checks. Larger than every source file but `app.js`, `index.html` and the two stylesheets — the pinned source tables account for most of it. |

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
| ~1180 | Battery table (Score Tables page) |
| ~2377 | APA table notes — single source of truth |
| ~2847 | SDI |
| ~3232 | RCI calculators (Basic / Practice / SRB / Crawford) |
| ~3797 | Per-method autofill from the normative database |
| ~4024 | Norms Database (custom-test storage, import, the browser) |
| ~4453 | Premorbid estimation |
| ~5771 | Auth overlay |
| ~5917 | Top-bar navigation bucket sync |
| ~5991 | Score Converter view-mode tabs |
| ~6010 | Working report bundle |

Line numbers drift with every edit — treat them as a starting point and grep to
confirm. Function-name families make grep effective: `calc*` computes, `render*`
draws, `get*` derives, `fmt*` formats, and page prefixes are `pre*` (premorbid),
`rci*`, `sdi*`, `bat*` (battery), `opie*`. To find any calculator, grep the prefix.

Statistical primitives sit at the top: `erf` (~7), `normCDF` (~16), `tInv` (~73).
General formatters follow at ~147–175 (`fmt`, `fmtPct`, `fmtP`). Note that
`fmtIntOrDash` and `fmtPctBr` are **not** there — they are premorbid-local, at ~4479.

Pages hold state in module-level objects: `preState`, `rciState`, `batteryRows`,
`sdiRows`. Note that tables are sometimes **updated in place** rather than re-rendered,
specifically so a focused input does not lose focus mid-typing — see `refreshOpieDerived`
for the pattern. Follow it when touching a table that contains inputs.

---

## `normDB`

The reliability database: 116 groups, 671 entries, seven instrument families (D-KEFS,
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

Seven manuals clear it, and `check.js` §24 pins the roster so an eighth cannot appear
undocumented:

| Source | Entries | Fields |
|---|---|---|
| CVLT-C `List A Trials 1-5 Total`, ages 8 / 12 / 16 | 3 | `rInternal` only |
| D-KEFS Advanced — CWIT, Tower, Social Sorting, Risk-Reward | 26 | `rInternal` + `rInternalByAge` |
| D-KEFS (original) — VFT, Sorting, 20Q, Word Context, Tower, Proverb, TMT composite | 13 | `rInternalByAge` **only** |
| WAIS-IV — Table 4.1, all but the three speeded subtests | 21 | `rInternal` + `rInternalByAge` |
| RBANS Update — Table 3.6, the nine rows without footnote a | 9 | `rInternal` + `rInternalByAge` |
| WMS-IV — Table 3.1, both batteries, all but VPA II Word Recall | 30 | `rInternal` + `rInternalByAge` |
| WISC-V — Table 4.1, all but the three speeded subtests and the two Cancellation process scores | 29 | `rInternal` + `rInternalByAge` |

That roster check keys on **either** field. D-KEFS original carries no `rInternal` at
all, so a check testing `rInternal` alone would let 13 entries change the basis of a
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
`patientAge()` — **one shared value**. The master is **`#patient-age` in the top bar**,
mirrored into `#pre-age` on Premorbid; `PATIENT_AGE_INPUTS` lists the master **first**,
because `patientAge()` returns the first input holding a finite value. Pinned by
`check.js` §26. It is **optional**. Blank, or outside `rInternalAgeMax`, falls back to
`rInternal` where one is published — the manual's own all-ages figure, so both paths are
citable — and to the retest `r` where none is, which is the D-KEFS original case below.
**Never make it required**: a blank age silently emptying the CI column would read as the
app being broken.

It used to be `#bat-age`, declared in the Score Tables markup and moved by
`design-system.js` into that page's inline control bar. That cost three fixes in
succession, all from one cause — **patient-level data living in a page-level control
bar**: it shipped rendering at 0×0 inside a wholesale-hidden panel (`d6e80e8`); once
moved, it took width from the test-family search and its label had to be truncated to
"Age" (`5cea561`); and being revealed only with the CI toggle, a clinician who never
turned CI on never saw an age field at all. The inline bar now holds **view settings
only** — Show Raw, Score CI, Classification. If you are tempted to put another
patient-level field on a page, put it in the top bar instead.

A `.ds-patient-field.is-live` pip beside the field lights only where a measure actually
reads the age — `patientAgeIsInUse()`, which on Score Tables uses **the same expression
that decides the APA note's `ciAge`**, so screen and export cannot disagree. It is
refreshed from `renderBattery()` (rows and CI level change without the age being touched)
and from `navigateTo()` (the page changes with neither changing). The pip is toggled by a
**class**, never `[hidden]` — the element it replaced used the attribute and stayed on
screen, because its own `display:flex` outranked the browser default.

#### The missing-age popover

`#bat-age-pop` is a glass popover anchored to the **Score CI toggle** — the control that
creates the moment a blank age starts costing a sharper interval. It carries its own age
input, so the ask is also the fastest way to answer it; the value is written through to
`#patient-age`, which stays the master.

**Edge-triggered, not state-triggered.** `renderBattery()` runs on every keystroke in the
table, so opening whenever the condition holds would re-open it continuously while
someone types scores. It fires on the *transition* into "CI on, age-band measures present,
no age" — which also fixes an ordering problem a click handler on the CI button would
have: switching CI on with D-KEFS already loaded fires, **and** autofilling D-KEFS while CI
is already on fires. A click handler catches only the first, and would fire on an RBANS
table where the age changes nothing.

**Once per patient**, re-armed by an explicit reset — never by clearing the age alone, or
someone who skipped and then blanked the field is asked twice about the same person.
**There are two reset controls and conflating them is how this shipped broken**: the
table's own "Clear all" (`#bat-clear` → `clearBattery`) empties the table, and the topbar
"New patient" (`#topbar-clear-all`) empties every table *and* the age. Only the topbar one
re-armed at first, so clearing the table and autofilling again got silence. Both re-arm
now; removing a single row must not, or a table being edited down re-asks mid-edit.
`check.js` §26 pins exactly **two** runtime `batAgePopArmed = true` sites, by name.

**The opening click must not also dismiss it.** The popover opens from inside
`renderBattery()`, which usually runs *during* a click — "Add selected tests", the CI
toggle, a row edit. That same click keeps bubbling to the document-level outside-click
handler, which finds a popover that is now open with a target outside it, and closes it:
OPEN then CLOSE, nothing on screen, no error. `batAgePopJustOpened` is set on open and
cleared on a `setTimeout(…, 0)`, so it survives exactly the synchronous dispatch of the
opening click. An earlier version special-cased the CI toggle *by selector*, which fixed
one route and left autofill — the commonest — broken; the guard has to be about **when**,
not about which element.

**That bug was invisible to every check and every browser probe**, because they all called
`renderBattery()` directly and no click was ever in flight. It only exists on the real UI
path. When verifying anything that opens from a render, drive the actual controls —
`.combo-add`, `#bat-clear`, `[data-bat-ci]` — not the functions behind them.

The pip and the popover are opposite halves of one question and both route through
**`batteryAgeBandRowCount()`**, so a table lighting the pip while asking for an age is
unrepresentable. That count excludes `isExample` rows: seeded rows are not the clinician's
data and must not be counted into a claim about "measures in this table".

**It offers, it never scolds, and it never blocks.** No backdrop, nothing dimmed, Escape
and outside-click dismiss, and the skip reads *"Skip — use all-ages"* rather than
"Dismiss" because that names the alternative actually being chosen. Primary tint, never
`--ds-warning`. §26 pins all of it, including that every clause agrees with the count —
the first version read "1 measure in this table **publish their** reliability", which is
what assembling the sentence from shared fragments produces.

Committing an age from the popover pulses the topbar field for ~2s
(`.ds-patient-field.is-pulsing`, three 0.7s cycles). The popover sits over the table, so
the field the value lands in is somewhere the clinician was not looking — the pulse
answers "where did that go?". **Only on commit from the popover**: typing directly into
the field must not pulse it, because a field that animates as you type in it reads as an
error. §26 pins both halves, the restart-on-reflow, and the reduced-motion fallback —
motion *is* the mechanism here, so suppressing it has to leave a held ring rather than
nothing.

A transient popover can be missed, so **`.ds-patient-field.is-wanted`** is the residual:
a dashed tint on the topbar field whenever an age would be read and none is set. Purely
state-driven, so it cannot go stale, and mutually exclusive with `.is-live` by
construction. That trace is why a purely transient prompt was not enough.

##### `body{zoom:0.9}` splits measurement in two

`styles.css` line 33 applies a deliberate global 10% downscale. Anything positioning an
element from a measured rect has to respect it, and mixing the halves is silent:

| | space |
|---|---|
| `getBoundingClientRect()` | **visual** px, already scaled |
| `offsetWidth` / `offsetHeight` | **layout** px — reports 320 for a box occupying 282 |
| `style.top` / `style.left` | **layout** px — the browser applies the zoom |

The first version of `positionBatteryAgePop` read rects and `offsetWidth` together and
wrote the result straight back, putting the popover **19px above its anchor and 58px off
centre**. Nothing threw, no check noticed, and re-positioning did not shift it — the error
is systematic, not a layout-timing race, which is the tell. It was found by measuring the
rendered page. `pageZoomFactor()` reads the factor off the element rather than hardcoding
`0.9`, and §26 fails if `offsetWidth` returns, if the divide disappears, or if that `0.9`
leaves `styles.css` without this being re-derived.

The topbar button is **"New patient"**, not "Clear all tables": it clears every table
*and* the age, because once age is a header-level property of the patient, two controls
that each clear half of one is how a previous age survives onto the next person's report.

#### The age also filters the family dropdowns

`familyMembersForAge` (`app.js`, beside the age-band helpers) narrows each family in
Quick Add, Change Analysis and the SD Index to **the band(s) containing the age, plus
`· All Ages` where one is published**. `buildFamilyListHtml` applies it per family and
`rebuildFamilyListsForAge` rebuilds all three lists when the parsed age changes.

Three clauses, and the third is the one that matters:

1. **Plural, not singular.** WMS-IV publishes overlapping `Ages 16-69` and `Ages 65-90`,
   so a 67-year-old is in both and picking one would hide a real choice.
2. **`· All Ages` always survives.** `comboAgeBandNoteHtml` tells the clinician All Ages
   has the larger *N* and stronger *r*; filtering it away would delete the option the app
   itself recommends weighing. The band choice is a clinical judgement, not an age lookup
   — which is why this filters rather than auto-selects.
3. **If 1 and 2 are both empty, show every member unchanged.** Eight of the 28
   numerically banded families publish no All Ages entry, and `CVLT-C Subtests (Raw
   Scores)` is banded by the three normative *sample ages* (`Age 8` / `Age 12` /
   `Age 16`) — so **88 of the 91 ages from 5 to 95 match no band**. Without this clause a
   clinician testing a 10-year-old would open Quick Add and find CVLT-C gone, which reads
   as the app being broken rather than as a filter working. RBANS, RBANS Form B, CVLT-3
   and WMS-IV are the others.

**The invariant: an age may quieten a family's list; it may never remove the family.**
`check.js` §33 asserts it over every family at every age 5–95, and it drives the shipped
functions rather than re-implementing them, because the thing under test is the
interaction between the rule and the actual shape of `normDB` — a second copy of the rule
would agree with itself while both drifted from the database.

Whenever the list has been narrowed it says so and offers **"Show all bands"**
(`familyListShowAllBands`), because a silently shorter dropdown is indistinguishable from
missing data. That note must **not** carry `.combo-item`: `filterFamilyListEl` hides those
by search text, so it would vanish the moment anyone typed — exactly when a shortened list
is most confusing. §33 pins that too.

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

##### Table 2.8 is also why the range-restriction correction is offered but never the default

**Do not "discover" this and switch the default.** The reasoning is sound right up to the
last step, which is exactly what makes it worth recording.

Every interval in this app is `SD × √(1 − rxx)`, valid only when both terms describe the
same group. 298 entries fall through to the retest study's own `r` paired with a
**normative** SD — two populations. `rCorrected` is the fix where a manual prints one, but
D-KEFS, D-KEFS Advanced and CVLT-C print none at all. The Allen & Yen / Magnusson
correction repairs it arithmetically:

```
rxx(unrestricted) = 1 − (sd1² ÷ normSD²)(1 − r)
```

and it is a **good formula, not a guess**: over the 267 entries carrying both a raw `r` and
a published `rCorrected` it reproduces the publisher's own figure to a median error of
**.003**, 193 of them inside .005. It moves **246 stored entries, 46 of them reachable from
Score Tables** — 25 D-KEFS, 21 D-KEFS Advanced, and nothing else.

Making it the default was rejected because **every one of those 46 belongs to a manual that
has already made this pairing itself, deliberately.** D-KEFS Technical Manual p. 19 says
test–retest SEMs "were derived from the total sample of cases" and fixes "The standard
deviation unit is 3 for all D-KEFS scaled scores" — i.e. `3 × √(1 − r)` on the
**uncorrected** total-sample `r`. Table 2.8 is that arithmetic on the page: Design Fluency's
three All Ages SEMs are 1.94 / 1.97 / 2.47, and `3 × √(1 − r)` gives exactly those, where
the corrected coefficient gives 1.78 / 1.95 / 2.43. The remaining 21 are D-KEFS Advanced
Trail Making and Verbal Fluency, whose Table 3.4 rows *are* the retest coefficients — that
manual's stated choice for its speeded tests.

So correcting by default would print a coefficient the cited manual does not contain, on the
one family whose own working is checkable. **Same ground as the declined unrounded Fisher's-z
WAIS-IV average**: the app renders the coefficient it actually used, and a clinician
cross-checking the manual must find it there.

**What shipped is the label and the choice.** `retest, uncorrected` on the Data page names
the compromise rather than hiding it, and the **Score Tables reliability control**
(`#bat-ci-basis`, Published / Corrected) lets a clinician take the corrected reading
deliberately. Off by default, exactly as the RCI pages' corrected-`r` toggle is, and for the
same reason. The APA note carries a different sentence in each state, so an exported table
always says which basis produced it. The Data page **follows the control** — showing the
published reading while the table is set to corrected would be precisely the drift the shared
resolver exists to prevent.

What the choice costs, so it can be made with the size in view: at 95% it moves **9 of the 46
printed margins wider and 4 narrower; 33 are unchanged** once rounded. The largest single
shift is **2 scaled-score points**, on D-KEFS Advanced Multitasking Index.

Three counts to know before treating a shortfall as a bug. Of the 298 entries on the retest
pairing, **36 CVLT-C** are refused (`sd1` in raw words while the row displays T or z), **2**
produce a value outside (0, 1) — D-KEFS Design Fluency Switching (Ages 20-49) and D-KEFS
Advanced Social Sorting Total Number of Conceptual Level Responses (Ages 8-18), both retest
samples being more variable than the norm group on a coefficient near .22 — and **14** have
`sd1` exactly equal to the normative SD, where the correction is the identity because the
sample was not restricted at all.

`check.js` §27 drives the shipped renderer in **both** states and asserts four things: the
uncorrected `r` reproduces Table 2.8; the **default** reading still does, so an untouched app
prints the published interval; the corrected reading does **not**, so the two genuinely
separate; and a **missing control reads as published**, so no harness can land on derived
coefficients by accident. It also proves the toggle cannot move a published coefficient.
§32 drives the Data page and the table in both states and pins the 246. If a corrected
printing of the manual ever makes the two readings agree, §27 fails — that is the signal to
re-read this section rather than assume a regression.

Four exclusions, four different reasons — each an easy silent mistake, each pinned:

- **Colour-Word Interference, all four.** Its only coefficient table (2.9) is for the
  *Combined Colour Naming + Word Reading composite*, which is not a measure in `normDB`.
  Attaching it to Colour Naming or Word Reading would give a condition the composite's
  reliability.
- **Design Fluency, all three.** "Item interdependence precluded the use of internal
  consistency procedures."
- **Trail Making, five of six.** Table 2.1 is the composite alone, so only
  `Combined Number + Letter` carries it.
Twenty Questions is **not** an exclusion — both its measures carry Table 2.15 — but its
`Initial Abstraction Score` is the one place in the manual where two published tables
cannot both be right. See below before touching it.

Word Proverb is banded 16–19 upward, the Proverb Test not being administered under 16, so
a child's age finds no band and correctly falls back.

##### Twenty Questions Initial Abstraction — where the manual contradicts itself

**The stored values are Table 2.15, verbatim. Do not "fix" them.**

Everywhere else this manual's coefficient tables and SEM tables satisfy its own formula.
Here they cannot both be right:

```
ages 40–49   Table 2.15 prints rxx .75
             3 × √(1 − .75) = 1.50
             Table 2.17 prints SEM 1.76, which needs rxx .656
```

12 of 16 bands disagree, by up to .094, **in both directions** — so it is not a transform
anyone has missed. Ruled out by test: every other published coefficient column, every row
shift, any constant SD other than 3, and the Spearman-Brown-uncorrected split-half. The
column beside it in the same table reproduces 16/16. Nothing inside the document says
which table is the misprint — there is no worked example or third statistic for this test
to triangulate on, which is what settled the equivalent question for CVLT-C.

**Why 2.15 and not 2.17.** The app stores coefficients and derives the SEM, so there is no
field for a published SEM. Honouring 2.17 would mean storing `.656`, a number the manual
never prints — the one thing this project forbids. 2.15 is the publisher's own coefficient
through the publisher's own formula, so every stored digit is citable.

**What it costs.** Three bands print one scaled-score point away from Table 2.17: ages 11
and 12 (±2 here, ±3 there) and 80–89 (±3 here, ±2 there). Note the bands with the *largest*
coefficient gaps — 40–49 at .094, 30–39 at .069 — are **not** among them; both readings
round to ±3 there. Severity in the table and visibility on screen are nearly unrelated, so
judge this kind of discrepancy on the printed margin, not the coefficient.

`check.js` §27 pins the discrepancy itself. If someone edits a coefficient to make 2.17
reconcile, it fails; if a corrected printing ever makes the tables agree, it *also* fails,
which is the signal to re-read this section.

**The two D-KEFS manuals reach opposite conclusions on the same two test names.** Advanced
excludes Trail Making and Verbal Fluency for *speededness*; the original publishes internal
consistency for all four Verbal Fluency measures and the TMT composite, and excludes Design
Fluency for *item interdependence* instead. Do not "harmonise" them — §27 asserts the
inversion in both directions.

Reliability method is a **per-manual question, not a per-app policy.** Five sources, five
different answers: D-KEFS Advanced rejects internal consistency for *speededness*; D-KEFS
original accepts it for those same tests and rejects Design Fluency for *item
interdependence*; CVLT-3 and CVLT-C reject it for *item interdependence* on a word-list
task; CVLT-C nonetheless publishes split-half for its one non-interdependent score; WAIS-IV
accepts it everywhere except its three Processing Speed subtests, again for *speededness*.
Read the manual before assuming any of them generalises.

#### WAIS-IV — Table 4.1, and why the speeded three needed no change

Technical Manual (GB) ch. 4, p. 42: "Because Symbol Search, Coding, and Cancellation are
Processing Speed subtests, the split-half coefficient is not a proper reliability
estimate. Therefore, test-retest stability coefficients were used… corrected for the
normative sample's variability." So Table 4.1 is a **mixed-basis** table: 21 rows of
split-half/alpha, three of corrected stability.

Those three keep `rCorrected`, because that field **already holds Table 4.1's own value** —
all 38 published cells match it to the digit and not one matches the raw `r`, no other
measure exceeding 3 of 13 (`check.js` §28). Storing it in `rInternal` would change no
number while labelling a stability coefficient as internal consistency, which the Methods
& References paragraph then asserts on screen.

That 38/38 match is also the **transcription proof** for the whole family. Table 4.1
repeats one retest-study value across each span of normative bands it covers — Symbol
Search .81 at 16-17/18-19/20-24/25-29, .73 at 30-34/35-44/45-54 — and those spans
reconstruct `normDB`'s 16-29 / 30-54 / 55-69 / 70-90 groups exactly.

- **`rInternalAgeMax` is 69 for Letter-Number Sequencing and Figure Weights.** The manual
  prints a dash above 65-69 for both, they being normed to 69. Without it the "greatest
  key ≤ age" lookup would score a 90-year-old on the 65-69 coefficient.
- **PSI and FSIQ are hybrids and are included anyway.** The manual concedes the PSI
  average "is based on test-retest subtest reliabilities", but it computes and labels
  every composite coefficient as internal consistency and publishes no other reliability
  for them. Excluding PSI while keeping FSIQ — which also draws on Processing Speed
  subtests — would be incoherent.
- **`WAIS-IV Longest Span` gets none of this** — base-rate scored, no retest data, absent
  from Table 4.1.

The §24 naming check could not catch a stale paragraph here: **"WAIS-IV" was already in
that paragraph in the opposite sense**, naming the manual only as a user of test–retest
coefficients, so the roster grew and the check stayed green. §28 pins the substance of the
wording instead.

##### Table 4.3 — the SEM table, which settles what used to be taken on trust

This section long carried a "NOT verified" note: that the manual builds its own SEMs from
Table 4.1 rested on its stated method, Table 4.3 being unavailable. It is available now,
and **all 300 published cells equal `populationSD × √(1 − rxx)` from the stored
coefficients, exactly, at the printed 2 dp** (`check.js` §28). Nothing needed changing.

Three things follow, and the second and third matter more than the confirmation.

- **A second publisher states the SD rule.** Table 4.3's note: the SEMs use "the
  reliability coefficients shown in Table 4.1 and the population standard deviations
  (i.e., 3 for the subtests and 15 for the composite scores)". That is exactly what §23
  derived independently from CVLT-3, and what `6de0d81` fixed. Over the 300 cells: this
  app's pairing **300**, retest `r` 5, `rCorrected` 46 (the 38 speeded cells it owns, plus
  8 coincidences), and **`sd1` against a normatively-scaled coefficient — the old pairing
  — 29**.
- **The blank-age fallback was unpinned, and now is not.** `rInternal` is what a patient
  with no age is scored on, and nothing held it — the spot checks read `rInternalByAge`
  only. Proven by mutation: moving Vocabulary's average .94 → .95 passed all 238 checks.
  Because Table 4.3 fixes the 13 bands exactly and Table 4.1 averages with **Fisher's z**,
  the average is now derivable from pinned values — Fisher's z reproduces all 24 stored
  averages at 2 dp where the plain arithmetic mean manages 18. **That 18 is also why the
  neighbouring check refuses to assert the average *differs* from the mean of its bands.**
- **The two tables average the same 13 numbers by different methods, so nothing
  reconciles them.** Table 4.3's footnote defines the average SEM as the RMS of the band
  SEMs, and `RMS = √(mean(SD²(1−rᵢ))) = SD√(1 − mean(rᵢ))` — the **arithmetic** mean of the
  coefficients, which reproduces 20 of 24 against **5** for the **Fisher's z** that Table
  4.1's own average column uses. That is a property of the source. The blank-age interval
  therefore keeps coming from Table 4.1's average *coefficient*: it is what an unknown band
  needs and it is published, where inverting the average SEM would mean storing
  `1 − (2.16/15)² = .9793` for FSIQ — the Twenty Questions rule. §28 pins the inconsistency
  itself, so it fails if a corrected printing ever makes the two agree.

It costs **two printed margins in 48** (24 measures × the two CI levels), both on a
rounding knife-edge: FSIQ at 90% gives ±3 against ±4 (3.49 vs 3.55) and Digit Span
Backward at 95% ±2 against ±3 (2.49 vs 2.55). §28 pins the pair, so a third joining them
fails the check.

**The gap is mostly rounding, not method** — the obvious reading is wrong, and it is worth
knowing before re-opening this. Decomposing FSIQ: the arithmetic mean of the bands
(.97923) gives 2.1617, which *is* Table 4.3's printed 2.16; Fisher's z unrounded (.97936)
gives 2.1547; the stored published average (.98) gives 2.1213. **Method costs 0.007;
rounding the published average to 2 dp costs 0.033** — about five times more. On VCI it is
0.044 against 0.19.

**A third option was weighed and declined — do not "discover" it and switch.** Using the
Fisher-z average of the bands **unrounded** matches the manual's printed margin on 24 of
24 at both CI levels, and invents nothing: the bands are pinned exactly by Table 4.3 and
Fisher's z is Table 4.1's own stated method. It was declined because the app renders the
coefficient it actually used, so the reliability cell would print `.979` where the manual
prints `.98`, and a clinician cross-checking Table 4.1 would find a number that is not
there. Reviewed and kept as-is, 2026-07-31. Entering an age sidesteps the whole question,
and that path is exact.

#### `rStability` — the same table, the other basis

**RBANS Update Table 3.6 is the first published reliability table in this database that
is mixed-basis *within itself*.** Its footnote a marks five of fourteen rows "Reliability
estimates based on test–retest": Figure Copy, Semantic Fluency, Coding, Story Recall and
Figure Recall. The other nine are internal consistency.

Those five cannot take `rInternal*` — that field's name is asserted on screen by the APA
note and by Methods & References, so storing a stability coefficient there would state
something false, exactly as it would for WAIS-IV's three speeded subtests. But they also
should not fall back to the retest *group's* coefficient, because Table 3.6 publishes a
value for them **by normative age band** and derives its own printed SEM from it.

Hence `rStability` / `rStabilityByAge` / `rStabilityAgeMax`, read by
`getBatteryRowReliability` immediately after the internal-consistency fields. **It is not
an exception to the test–retest default — it *is* that default**, sourced from the
manual's own reliability table rather than from a retest study group. Five entries carry
it and `check.js` §29 asserts it has not escaped RBANS.

Corroborated: the adult bands equal this database's stored `rCorrected` **5 of 5**,
against 1 of 5 for the raw `r` — the same transcription proof the WAIS-IV speeded three
give. The adolescent bands match only 2 of 5, those being a different retest sample;
Table 3.6 is stored as printed rather than reconciled to Table 3.8, and §29 pins the 2 so
the mismatch reads as known rather than as an error.

#### RBANS needed new `· All Ages` groups, and that fixed a live defect

Every other RBANS group holds a retest study banded the way that study sampled (12-19,
20-89). Score Tables shows one entry per instrument and picks the `· All Ages` group —
and with none present, `buildFamilyListHtml` fell through to **whichever group was listed
first**. Nothing chose it; it was object order. So every RBANS patient was scored on 55
adolescents, differing from the 20-89 study on **12 of 18 measures**. Immediate Memory at
an index of 100 printed 85–115 where Table 3.6 gives 90–110 at age 80.

`RBANS Subtests · All Ages` and `RBANS Indices · All Ages` now carry the normative
metrics (10/3 and 100/15, which are definitions, not data) and Table 3.6's coefficients.
Every entry is `singleAdministration:true`, so `isSingleAdministrationFamily` keeps them
out of Change Analysis and the SD Index, which go on using the retest groups — the same
separation §25 enforces for D-KEFS Advanced.

**The four raw subtests appear nowhere in Table 3.6** — the manual publishes reliability
for its eight *scaled* subtests only, which is its own confirmation of the raw/scaled
split §18 asserts from the data shape. They therefore print **no interval at all** in the
All Ages group. That loss is deliberate: the interval they used to show came from 55
adolescents whatever the patient's age, and `List Recognition` alone runs `r` .70 there
against .27 in adults, so the honest adult interval is nearly twice the width that was
printed. Their `m1`/`sd1` are the adult retest descriptives, carried over only so the row
stays selectable and declares its metric; §29 asserts nothing derives from them.

Note this made `singleAdministration` mean less than it used to. §3 required base rates
on every such entry, conflating "no retest data" with "scored by base-rate lookup". These
entries are scored by ordinary conversion, so the check now asks for m1/sd1, no retest
fields, and **some** published basis — base rates, a reliability, or the raw tag.

#### WMS-IV — two batteries, and the first band that must stay selectable

**WMS-IV's age-band groups are separate norm sets, not separate retest samples.** The
Adult battery is normed 16-69 (15 subtests, 5 indices); the Older Adult battery is normed
65-90 and drops Designs, Spatial Addition and VWMI, leaving 8 and 4. Ages **65-69 are
normed in both**, with different coefficients — Logical Memory I is `.79` adult against
`.83` older adult — so age cannot pick the battery. The clinician picked it when they
decided what to administer.

So this family gets the opposite treatment to RBANS: no `· All Ages` group, and the band
stays **selectable**. `familyScoredByAgeBand()` now returns true for `baseRates` **or**
`separateBattery`, the latter declared per entry rather than inferred — a rule like "the
two groups hold different measures" would have worked today and broken silently later.
Collapsing them had been showing an 80-year-old the adult measure list and its
coefficients, wrong on **10 of the 12 shared measures**. §30 asserts the 65-69 overlap
both *exists* and *disagrees*, since that is the entire argument for the exemption.

Tables 3.1 and 3.3, pp. 44-46, verified arithmetically: **all 240 published SEM cells**
equal `SD √(1−rxx)` at the printed 2 dp, on 3 and 15. §25's "All Ages only" rule now
reads "All Ages, or a separate battery", and `separateBattery` is pinned to WMS-IV's four
groups so it cannot spread.

- **VPA II Word Recall carries footnote b** — free recall has no consistent item count,
  so internal consistency is inapplicable. It takes `rStability`, and Table 3.1's `.76`
  equals the stored `rCorrected` in **both** batteries and neither raw `r`.
- **The recognition memory measures are absent, and must stay absent.** Their published
  reliability is a **decision-consistency percentage** — percent agreement of impaired vs
  not-impaired at a 10th-percentile cut — because those tasks are cumulative percentages
  with skewed distributions. A percent agreement is not a correlation and cannot enter
  `SD √(1−rxx)`. They are in neither Table 3.1 nor `normDB`; §30 fails if one appears.
- **`rInternalAgeMax` is the battery's ceiling**, 69 and 90, so an age past it takes the
  published average rather than re-reading the top band.

**Three manuals now split their two average columns the same way** — WAIS-IV 4.1/4.3,
RBANS 3.6/3.7 and WMS-IV 3.1/3.3 all average coefficients by Fisher's z and SEMs by RMS,
which is `SD √(1 − arithmetic mean)`. It is a publisher habit, not a one-off, and the
stored average stays the published *coefficient* in every case.

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
`rCorrected` on **267** — D-KEFS and CVLT-C have none at all, so any feature depending on it
must degrade gracefully.

**WISC-V used to be on that list, and that was an error of fact rather than a property of the
source.** Its Table 4.7 prints a corrected `r` for all 34 rows; the database had simply never
captured it. Nothing scored differently — `rInternal` outranks `rCorrected` in the CI chain
and reliable change defaults to the raw `r` — but the corrected-`r` option on the Basic RCI
page was silently unusable for WISC-V alone. All 22 existing entries were backfilled from
Table 4.7 and §3 now asserts the position instead of assuming it.

WISC-V is still **verified differently from the rest** in one respect: Tables 4.1/4.4 carry
the transcription proof on their own (242 cells), and §31 drives the shipped renderer through
the manual's worked example — a 6-year-old at FSIQ 108 gives 102–114 at 95% and 103–113 at
90%, exactly as p. 62 prints. It also has the finest lookup in the database: **single year of
age**, 6 to 16, not bands.

Table 4.7 also filled the family out. It was holding 22 of the 34 measures the manual
publishes; the 7 process scores now live in `WISC-V Process Scores · All Ages`, mirroring
WAIS-IV, and the 5 ancillary indices sit with the primary ones. **`Cancellation Random` and
`Cancellation Structured` take `rStability`**, the manual naming Cancellation among the
subtests for which split-half is improper — confirmed from the data too, since a stability
coefficient broadcast from the retest study's coarse bands repeats across single years (3 and
2 distinct values across 11 bands, against 7–9 for the internal-consistency rows).

That table came from a photograph rather than a spreadsheet, so it was checked two ways before
being stored, and both checks are now live in §31: Cohen's Formula 10.4 reproduces the printed
Standard Difference from the means and SDs on 34 of 34 rows, and the 22 measures already held
matched their `m1`/`sd1`/`m2`/`sd2`/`n`/`r` exactly — which also proves the `r₁₂` column was
not confused with the corrected one, stored `r` matching `r₁₂` 22 of 22 and the corrected
column 1 of 22.

Users can add custom tests; `getMergedDB()` merges those over `normDB`.

### The Norms Database page

The only view of `normDB` a clinician has. It was a flat alphabetical list of all
115 groups, fully expanded — D-KEFS alone ran to **36 consecutive groups** — with columns
stopping at the retest `r`.

It is now **one table over every entry**. Instrument, category and age band come out of the
group key (`"<Instrument> <Category> · <Age band>"`) and become real columns, so they can
be filtered and sorted rather than scrolled past. That is the point of the shape: it is the
only arrangement that answers a question cutting *across* instruments — which measures are
weakest, which are still on a raw retest coefficient, which move with the patient's age.

Two columns say **which reliability the interval is actually built from and what it rests
on**. 131 entries are scored on a coefficient the old columns never showed.
`dbReliabilityBasis()` and `getBatteryRowReliability()` were two functions reading the same
fields in the same order — **if they drift, the page tells a clinician their interval rests
on a coefficient it does not** — and the mirror had to be strengthened twice, because a
mirror is only as good as the next person's memory. They now **share one resolver,
`resolveCiReliability(entry, normSD, age)`**, so the preference order cannot drift at all.
What can still drift is the *arguments*: the correction-relevant one is the normative SD,
which depends on the displayed metric, so the page asks `inferScoreTypeForSubtest` — the
same function autofill uses to set `row.scoreType`. §32 drives both entry points over all
671 entries and compares the coefficient **and** the basis string.

Internal consistency and published stability **never co-occur on an entry** — 0 of 659, a
measure being one or the other within its manual's reliability table. What does co-occur is
the published coefficient and the retest study's: 112 entries carry internal consistency
alongside a retest `r`, 5 carry stability alongside one. Both are visible, `r`/`corr. r`
being the retest study and `CI r` the published figure, with the basis chip saying which
kind it is. A pair of Internal/Stability columns was mocked up and declined: one of the two
is blank on every row.

Two rows print no interval, for different reasons, and the basis distinguishes them:
`base rate — no interval` (scored by published lookup, never had a coefficient) and
`none published` (the four raw RBANS subtests, absent from Table 3.6).

**The two retest labels are a distinction, not a synonym pair**, and §32 pins the whole
distribution so a change of basis anywhere has to be acknowledged rather than slide through:

| Basis | Entries | What it means |
|---|---|---|
| `internal consistency` (+ `· by age`) | 3 + 115 | published split-half / alpha |
| `stability, published · by age` | 12 | the manual's own reliability table, stability rows |
| `retest, corrected` | 180 | a published `rCorrected` |
| `retest, uncorrected` (+ `· by age`) | 285 + 13 | the retest study's own `r`, used with a **normative** SD |
| `retest` | 8 | raw display — `r` beside its own sample's `sd1` |
| `none published` / `base rate` | 4 + 51 | no interval at all |

`retest, uncorrected` is the honest name for a real compromise: `r` describes the retest
study's sample and the SD describes the norm group. **It is not a defect to be repaired by
default — see the D-KEFS note above.** With the Score Tables reliability control set to
Corrected, 246 of those entries move to a sixth basis, `retest, corrected here`, and this
page follows the control rather than showing a stale reading. Plain `retest` is exactly the four raw RBANS subtests in each
of the two retest bands; nothing else in the database is displayed raw, so if that 8 moves,
a family has gained or lost `metric:'raw'`, which changes what a percentile means.
CVLT-C lands on `uncorrected` and belongs there, which is not obvious: `metric:'raw'` makes
it look like the second case, but `reportedAs` puts the row on T or z, so the SD in force is
the normative 10 or 1 while `r` was measured on raw words recalled.

**`DB_COLUMNS` gives each column one `get`, used for both display and sorting**, so a
column cannot order by something other than what it shows. §32 pins that, and pins the two
reliability columns to `dbReliabilityBasis` specifically — found by mutation that pointing
`CI r` at the raw retest `r` otherwise passed everything, because the column then showed
and sorted on the same wrong number.

It is a real `<table>`, deliberately. The grouped list laid each row out as a CSS grid whose
column count was declared in two places, and a mismatch slid every cell one place left under
the wrong heading. A table cannot do that, so that whole class of fault is gone rather than
guarded.

**The typed entry form was removed.** Adding a family and typing its norms row by row is
gone; import remains, so `ctValidateEntry` is now the *sole* gatekeeper into the clinical
database rather than one of two — §19 keeps pinning it. Custom measures are deleted per row;
a family disappears with its last one, which is why there is no separate "delete family"
control. The deletion was the exact hazard recorded above: `refreshAll()` is a top-level
init statement that called `refreshFamilySelect()`, so removing the callee without the call
would have thrown at boot and killed every statement after it. §17 catches that.

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
1.645 / 1.960 / 2.576 — see `PREMORBID_CI_Z` (~1470). Do not conflate the two.

**APA notes** live in the `APA_NOTES` object (~line 2388, under its banner at ~2377) and
render via `apaNoteHtml`.
That is the single source of truth for the note text under every exported table — change
methodology there, not in the markup.

Each note is rendered **twice**: into the exported table, and into the on-screen `info-box`
that mirrors it (`renderStaticApaNotes`, via `data-apa-note`). `ctx.onScreen` — passed only
by the mirror — is the **one licensed difference** between the two. A note may use it to
**drop** a sentence the surrounding page already states in full; it may not add, soften or
reword one, so the exported note stays the superset. One note uses it: `pre-opiepredict`
suppresses the UK caveat on screen, because the `.caution-box` at the top of that tab says
the same thing at greater length and two statements of one caveat read as two caveats. The
exported copy is the one that leaves the app, so it keeps the caveat unconditionally.
`check.js` §15 pins all four halves — export keeps it, mirror drops it, caution-box still
exists, and the flag still reaches the note (without which the split is inert).

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

`node tools/check.js` runs 290 headless checks: statistical primitives, score-conversion
round trips, `normDB` structural integrity, WAIS-IV values pinned to Technical Manual
Tables 4.5 (§4) and 4.1/4.3 (§28), RBANS Update Tables 3.6/3.7 (§29), WMS-IV Tables 3.1/3.3 (§30), WISC-V Tables 4.1/4.4 (§31),
OPIE-4 coefficients
pinned to Table eA5.8, worked OPIE
predictions, reliable-change thresholds and direction-neutral outcome labels, base-rate
reconstruction and monotonicity, percentile-tail clamping, the effect-size calculator,
Score Tables confidence intervals, documentation contracts, wiring (§16–17), the
raw-score metric (§18), the Norms Database view (§32) and age-band filtering of the
family dropdowns (§33).

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
