# Report Writer · Phrasing Iteration Brief

> Drop this file (and `phrasing-engine.html`) into a new Claude conversation.
> Use the conversation as a phrasing playground — no code commits, just prose
> iteration. Once you're happy with the wording, paste the finalised
> templates back into the original development conversation.

---

## What this engine does

The Report Writer page on a clinical neuropsych website auto-generates a
descriptive textual report from scores entered elsewhere on the site
(Neuropsych Tables, Change Analysis, Premorbid Estimate). The clinician
can then copy the prose into their Word report.

The page is **read-only narrative** — no test picker, no score entry. Just
a Reference selector (Mrs Doe / Mr Doe / You), a Descriptor system selector
(Wechsler / AACN), bracket-content toggles (Standard score, Percentile),
and the rendered report below.

## Accumulated writing rules

These are the user's explicit preferences. **Do not violate.**

1. **Pronoun-led, never passive.** Always lead with "Her", "Mrs Doe's",
   "His", "Your", etc. Never "Performance on X was Y" or "Estimates
   yielded…". Always something like "Her performance on X was Y".

2. **Sentences should not start with bare names.** Pronouns at the start
   are fine ("Her performance…"). Possessive names are fine ("Mrs Doe's
   performance…"). But *not* "Mrs Doe scored…" or similar bare-name
   subject-verb starts.

3. **Subject cycle within each domain section.** First sentence (the
   opener) uses the full possessive name ("Mrs Doe's"). Subsequent
   sentences alternate: pronoun → name → pronoun → name. So second
   sentence is "Her", third is "Mrs Doe's", fourth is "Her", etc.

4. **Long-form score labels.** "standard score of 100", not "SS=100".
   "scaled score of 10", not "Sc=10". "T-score of 50", not "T=50".
   When all-same-score collapse fires, use Title Case: "Scaled Score 10".

5. **Minimal connectors.** No "Furthermore,", "In addition,", "Moreover,".
   Sentences should flow as plain statements.

6. **Identical-score collapse.** When 2+ subtests in the same descriptor
   cluster have *identical scores*, collapse to:
   *"Mrs Doe's performance on X, Y, and Z was all in the average range
   (Scaled Score 10, 50th percentile)."*
   Use Title Case for the score label and singular form. Include
   percentile when the toggle is on.

7. **Brand stripping.** Battery `group` strings like
   *"WAIS-IV Core Subtests · All Ages"* or *"D-KEFS Trail Making Test ·
   All Ages"* should display as just *"WAIS-IV"* / *"D-KEFS"*.

8. **Per-brand paragraph split.** Within a domain, each test brand gets
   its own paragraph. First brand's paragraph opens with
   *"Mrs Doe's [domain] was assessed using X, Y from the [brand]."*.
   Subsequent brands open with
   *"She also undertook A, B from the [brand]."*. The subject cycle
   continues across paragraphs (so the opener of brand 2 might be
   *"Mrs Doe also undertook…"* if the cycle index lands on name).

9. **Subtest collapse when > 4.** If a brand has more than 4 subtests in
   the same domain, collapse the opener to *"using tests from the
   [brand]"* rather than listing every subtest by name.

10. **CVLT-specific narrative.** For CVLT verbal-learning brands in the
    *Verbal Learning and Memory* domain, replace the generic per-test
    listing with a structured narrative:
    - **Learning curve**: *"Her performance ranged from \[descriptor]
      on Trial 1 (Scaled Score X, Yth percentile) to \[descriptor] by
      Trial 5 (Scaled Score X, Yth percentile)."*
    - **Short delay**: *"After a short delay, Mrs Doe's performance on
      Short Delay Free Recall was X (Scaled Score …, Yth percentile)."*
    - **Long delay**: *"Finally, after a long delay, her performance on
      Long Delay Free Recall was X (Scaled Score …)."*
    - **Recognition assistance**:
      - 6+ scaled-score increase from LDFR → "notably assisted by the
        recognition trial"
      - 3–5 increase → "slightly assisted"
      - <3 → fall back to plain recognition score report
    - **Auxiliary measures** (List B Correct, Recognition False Positive,
      Recognition Discrimination, Discrimination Nonparametric, Total
      Intrusions, Total Repetitions) go in a *separate paragraph* after
      the main narrative, using the standard per-test sentence builder.

11. **Premorbid format** — fixed structure:
    *"Mrs Doe's baseline cognitive ability was estimated using the
    Test of Premorbid Functioning (ToPF) with demographic adjustment
    (FSIQ ≈ 102, 90% CI 88–116), as well as the ToPF (raw score)
    (FSIQ ≈ 100, 90% CI 85–115), the Crawford & Allan (2001)
    demographic predictor (FSIQ ≈ 99, 90% CI 84–114), and the OPIE-4
    (Vocabulary and Matrix Reasoning) (FSIQ ≈ 105, 90% CI 90–120)."*
    - Order: ToPF demographic → ToPF raw → Crawford → OPIE-4
    - OPIE-4 label is dynamic: *(Vocabulary)*, *(Matrix Reasoning)*,
      or *(Vocabulary and Matrix Reasoning)* depending on inputs
    - CI label uses the user's selected interval (90% or 95%)
    - Convergence sentence emitted only when 2+ models AND every
      estimate produces the same descriptor band:
      - 2 models → "Both of these models converged on the prediction
        that Mrs Doe's ability was likely to fall within the
        \[descriptor] range."
      - 3+ models → "All of these models converged…"
      - 1 model → omit the second sentence

12. **Bracket content toggles** — what's in the parenthetical after each
    test result is user-controlled:
    - Standard score on, Percentile off → *"(scaled score of 10)"*
    - Standard score off, Percentile on → *"(50th percentile)"*
    - Both on → *"(scaled score of 10, 50th percentile)"*
    - Both off → no parens at all → *"…was average."*
    - Premorbid CIs and change-analysis statistics ignore these toggles.

## Section ordering (clinical convention)

Domains render in this order; empty domains are skipped:

1. Performance Validity (TOMM, embedded indicators)
2. Premorbid Functioning (ToPF, OPIE-4)
3. General Intellectual Functioning (WAIS, WISC, RBANS)
4. Attention and Working Memory (Digit Span, Symbol Span, etc.)
5. Processing Speed (Coding, Symbol Search, Trail Making, Color Naming, Word Reading)
6. Verbal Learning and Memory (CVLT, WMS Logical Memory)
7. Visual Learning and Memory (BVMT, WMS Visual Reproduction)
8. Language (Verbal Fluency, Vocabulary, Naming)
9. Visuospatial (Block Design, Matrix Reasoning, Figure Copy)
10. Executive Functioning (Color-Word Inhibition, Tower, Sorting, WCST)
11. Mood and Symptom Reporting (HADS, PHQ, BDI)

Change-analysis findings (RCI / SDI) are embedded into the matching
domain as a separate paragraph at the end of that domain.

## Test-family → domain mapping (current rules, regex order matters)

```
TOMM, Word Memory Test, Reliable Digit, MSVT  → Performance Validity
ToPF, OPIE, WTAR, Crawford & Allan, NART       → Premorbid Functioning

CVLT-3, CVLT-II, Verbal Paired, Word List,
  Hopkins Verbal, Rey Auditory, HVLT,
  Logical Memory                                → Verbal Learning and Memory
Rey Complex Figure Recall, Visual Reproduction,
  BVMT, Designs, Family Pictures,
  Spatial Recall                                → Visual Learning and Memory
*memory.*index, IMI, DMI, AMMI, GMI, VMI        → Verbal Learning and Memory

Digit Span, Letter-Number, Spatial Span,
  Symbol Span, Arithmetic, WMI                  → Attention and Working Memory
Coding, Symbol Search, Cancellation,
  Trail Making, TMT, PSI                        → Processing Speed

D-KEFS Color-Word Word Reading, Color Naming    → Processing Speed
D-KEFS Color-Word Inhibition,
  Inhibition/Switching                          → Executive Functioning

Verbal Fluency, Category/Letter Fluency,
  FAS, COWAT, Naming, Boston, Vocabulary,
  Similarities, Comprehension, Information,
  VCI                                           → Language
Block Design, Matrix Reasoning, Visual Puzzles,
  Figure Weights, Picture Concepts,
  Rey Complex Figure Copy, PRI                  → Visuospatial
Tower, Sorting, Design Fluency, WCST            → Executive Functioning
HADS, PHQ, GAD, BDI, BAI, mood, anxiety,
  depression                                    → Mood and Symptom Reporting

FSIQ, GAI, Full Scale, Ratio IQ, WAIS, WISC,
  RBANS Total, KBIT                             → General Intellectual Functioning
```

## Current sample output (working baseline)

Given:
- Reference: Mrs Doe (she/her)
- Descriptor system: Wechsler
- Brackets: Standard score on, Percentile off
- Battery rows: WAIS-IV (Digit Span 10, Symbol Search 12, Coding 7, Vocabulary 13)
                + WMS-IV (Logical Memory I 10, Visual Reproduction I 11)
                + CVLT-3 full set (Trial 1–5, List B, SDFR, SDCR, LDFR, LDCR,
                                   Recognition, Recognition False Positive,
                                   Recognition Discrimination, Discrimination
                                   Nonparametric, Total Intrusions,
                                   Total Repetitions)
                + D-KEFS Color-Word (Word Reading 12, Color Naming 10,
                                     Inhibition 8, Inhibition/Switching 9)
- Premorbid: 4 estimates all in average range

The engine produces (selected sections):

```
Premorbid Functioning

Mrs Doe's baseline cognitive ability was estimated using the Test of
Premorbid Functioning (ToPF) with demographic adjustment (FSIQ ≈ 102,
90% CI 88–116), as well as the ToPF (raw score) (FSIQ ≈ 100, 90% CI
85–115), the Crawford & Allan (2001) demographic predictor (FSIQ ≈ 99,
90% CI 84–114), and the OPIE-4 (Vocabulary and Matrix Reasoning)
(FSIQ ≈ 105, 90% CI 90–120). All of these models converged on the
prediction that Mrs Doe's ability was likely to fall within the average
range.

Attention and Working Memory

Mrs Doe's attention and working memory was assessed using Digit Span
from the WAIS-IV. Her performance on Digit Span was average (scaled
score of 10).

Processing Speed

Mrs Doe's processing speed was assessed using Symbol Search and Coding
from the WAIS-IV. Her performance on Symbol Search was high average
(scaled score of 12). Mrs Doe's performance on Coding was low average
(scaled score of 7).

She also undertook Color Naming and Word Reading from the D-KEFS. Her
performance on Color Naming was average (scaled score of 10). Mrs Doe's
performance on Word Reading was high average (scaled score of 12).

Verbal Learning and Memory

Mrs Doe's verbal learning and memory was assessed using tests from the
CVLT-3. Her performance ranged from superior on Trial 1 (Scaled Score
15, 95th percentile) to average by Trial 5 (Scaled Score 10, 50th
percentile). After a short delay, Mrs Doe's performance on Short Delay
Free Recall was average (Scaled Score 8, 25th percentile). Finally,
after a long delay, her performance on Long Delay Free Recall was low
average (Scaled Score 6, 9th percentile). Mrs Doe's performance was
slightly assisted by the recognition trial (Scaled Score 11, 63rd
percentile).

Her performance on List B Correct and Total Intrusions was low average
(scaled scores of 7 and 6). Mrs Doe's performance on Recognition False
Positive was extremely low (scaled score of 2). Her performance on
Recognition Discrimination and Discrimination Nonparametric was
borderline (scaled scores of 4 and 5). Mrs Doe's performance on Total
Repetitions was average (scaled score of 8).

Executive Functioning

Mrs Doe's executive functioning was assessed using Inhibition and
Inhibition/Switching from the D-KEFS. Her performance on Inhibition and
Inhibition/Switching was average (scaled scores of 8 and 9).
```

## Areas the user wants to iterate on

The user has flagged that the page is the most complex and wants tighter
language polish. Specific areas worth simulating:

- Wording for **mixed-descriptor clusters** (e.g. some average, some low
  average) — should they collapse differently, or is the current grouped
  form fine?
- Wording for **change-over-time** sentences when statistically reliable
  vs. not — does the current "*On retest, her performance on X showed
  reliable decline (t(RB) = …)*" feel right?
- Wording for **performance validity** when failed — should it be more
  direct ("This invalidates the assessment") or more measured?
- Wording for the **learning curve** when scores are stable across trials
  (currently still says "ranged from X to X" which reads oddly when the
  start and end descriptors are identical).
- Wording for **mood/behavioural inventories** — currently routed to a
  Mood section but no specific narrative template exists yet.
- **Domain transitions** — should there be any connecting wording across
  domains, or do clean breaks via headings work?

## How to iterate in the new conversation

1. Open a new Claude conversation.
2. Attach this file plus `phrasing-engine.html` (a runnable simulator).
3. Tell the new Claude:
   > *"This is a phrasing-iteration playground for a clinical neuropsych
   > report writer. The brief explains the rules and current state. The
   > HTML file simulates the engine. Help me iterate on phrasings — I'll
   > give you scenarios and feedback, you draft variations and update
   > the templates. Don't change the structural rules or the engine code
   > unless I explicitly ask. Focus on word-level polish."*
4. Paste a scenario like *"the patient has CVLT scores all 10 across all
   trials — the current learning-curve sentence reads oddly. Suggest 3
   alternatives."*
5. Iterate on each variation; pick the one you like.
6. When done, copy the finalised template strings into a single message
   and bring back to the original conversation. I'll patch the engine.

## Caveats

- The new Claude won't know about the dev environment or have access to
  the actual codebase. Keep iteration to **template strings only**.
- Structural changes (e.g. "add a new sentence type that fires when X")
  still need to come back here.
- The HTML simulator runs the current engine in a browser — open it,
  edit the seed data, see the rendered output. Useful for what-if checks.
