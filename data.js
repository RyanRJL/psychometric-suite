/* ============================================================
   DATA · constants, lookup tables, coefficients and tooltip strings
   Extracted from Psychometric Calculators · Clinical Suite
   ============================================================ */

let normDB = {}; // populated below + custom

/* ============================================================
   REPORT_TEST_CATALOG — batteries the Working Report may merge

   The Working Report auto-splits a Score Tables table into one item per test
   family, so a battery entered in one go arrives as several items ("WAIS-IV
   Indices", "WAIS-IV Core Subtests", …). "Merge by battery" puts them back
   together as one table with a labelled sub-section per family, which is how a
   report presents a battery.

   This list is what tells the merge which families belong to the same
   instrument. It was referenced by the code but had never existed — the merge
   silently did nothing for want of it.

   Each entry:
     id        stable key the merge groups on; never shown
     name      the abbreviation detectTestFamily() finds in the table text.
               MUST match a TEST_FAMILY_PATTERNS entry in app.js or nothing
               will ever resolve to this battery.
     longName  the merged table's title, so a medico-legal reader gets the
               instrument in full rather than an abbreviation
     families  prefixes that map onto this battery. Kept separate from `name`
               so an instrument whose normDB groups do not start with the
               abbreviation can still be attached.

   Only the seven instrument families in normDB are listed. ToPF and OPIE are
   deliberately absent: they are premorbid PREDICTORS, and their tables name
   WAIS-IV / WMS-IV as the predicted criterion. Cataloguing them would let a
   predicted-score table merge into an achieved-score one. The merge also
   refuses pre-* sources outright — belt and braces, since that particular
   confusion would put predicted and obtained scores in one column.
   ============================================================ */
const REPORT_TEST_CATALOG = [
  { id:'wais-iv', name:'WAIS-IV', families:['WAIS-IV'],
    longName:'Wechsler Adult Intelligence Scale – Fourth Edition (WAIS-IV)' },
  { id:'wms-iv',  name:'WMS-IV',  families:['WMS-IV'],
    longName:'Wechsler Memory Scale – Fourth Edition (WMS-IV)' },
  { id:'wisc-v',  name:'WISC-V',  families:['WISC-V'],
    longName:'Wechsler Intelligence Scale for Children – Fifth Edition (WISC-V)' },
  { id:'cvlt-3',  name:'CVLT-3',  families:['CVLT-3'],
    longName:'California Verbal Learning Test – Third Edition (CVLT-3)' },
  { id:'cvlt-c',  name:'CVLT-C',  families:['CVLT-C'],
    longName:'California Verbal Learning Test – Children’s Version (CVLT-C)' },
  { id:'d-kefs',  name:'D-KEFS',  families:['D-KEFS'],
    longName:'Delis–Kaplan Executive Function System (D-KEFS)' },
  { id:'rbans',   name:'RBANS',   families:['RBANS'],
    longName:'Repeatable Battery for the Assessment of Neuropsychological Status (RBANS Update)' },
];

// ToPF Raw (0-70) → estimated FSIQ
const TOPF_TO_FSIQ = [
  42,45,48,51,54,57,59,62,64,66,
  68,70,72,74,75,77,78,80,81,82,
  83,85,86,86,87,88,89,90,90,91,
  92,92,93,93,94,94,94,95,95,96,
  96,96,97,97,97,98,98,99,100,101,
  101,102,102,103,104,105,106,107,108,109,
  110,111,113,114,116,117,119,121,123,125,127
];

// WAIS-IV regression coefficients (ToPF + Demographics)
// Formula: Intercept + B1·ToPF + B2·ToPF² + B3·ToPF³ + Edu·Education + Sex·SexCode
// SexCode: Female = 1, Male = 2
const WAIS_COEF = [
  {idx:'FSIQ', label:'Full Scale IQ',              intercept:29.991, b1:2.09426,   b2:-0.0404559, b3:0.000340705, edu:1.4617126, sex:4.925, see:8.441},
  {idx:'VCI',  label:'Verbal Comprehension Index', intercept:52.873, b1:0.7367445, b2:-0.0120053, b3:0.000152838, edu:1.3280933, sex:3.491, see:9.277},
  {idx:'PRI',  label:'Perceptual Reasoning Index', intercept:20.695, b1:3.0900198, b2:-0.0642406, b3:0.000487938, edu:1.4928918, sex:5.734, see:12.368},
  {idx:'WMI',  label:'Working Memory Index',       intercept:42.779, b1:1.3486859, b2:-0.0177809, b3:0.000132589, edu:0.8494128, sex:7.112, see:10.617},
  {idx:'PSI',  label:'Processing Speed Index',     intercept:52.309, b1:2.0864576, b2:-0.0450319, b3:0.0003621,   edu:0.8079689, sex:0,     see:12.038}
];

// WMS-IV regression coefficients (ToPF + Age)
const WMS_COEF = [
  {idx:'IMI',  label:'Immediate Memory Index',      intercept:84.77,  b1:0.5426, age:-0.1779, see:12.032},
  {idx:'DMI',  label:'Delayed Memory Index',        intercept:84.438, b1:0.5598, age:-0.1904, see:12.42},
  {idx:'VWMI', label:'Visual Working Memory Index', intercept:87.994, b1:0.5265, age:-0.1783, see:12.165}
];

/* ============================================================================
   OPIE-4 (Holdnack, Schoenberg, Lange & Iverson, 2013; Table eA5.8)

   ⚠ ILLUSTRATIVE ONLY IN A UK CONTEXT. Do not quote these numbers as concrete
   premorbid estimates in a UK report. Reason below.

   The coefficients in the three tables that follow have each been checked
   against Table eA5.8 and reproduce it exactly — intercepts, subtest weights,
   Age, Age³, Age⁶ and Sex. Nothing here is mistranscribed.

   What is NOT applied: the published equations also carry Education, Ethnicity
   (African-American) and Region (West) terms, all dummy-coded (1 = belongs to
   group, 0 = does not). Those terms are deliberately omitted, which fixes every
   prediction at the US reference category:

       12th-grade high-school graduate · not African-American · not US West

   Omitting them is not a rounding detail. In the source table the education
   block spans up to 18.9 points (Matrix-only FSIQ, K–7th grade −11.57 through
   Bachelor's +7.13) and the African-American term reaches −8.31.

   THAT 18.9 DOES NOT RECONCILE WITH THE ENDPOINTS BESIDE IT: 7.13 − (−11.57)
   is 18.70. It cannot be settled from this project, because the education,
   ethnicity and region coefficients are stored nowhere here — only the applied
   terms are, and check.js §5 pins those alone. Table eA5.8 itself is needed.
   Two candidates: the span is right and the parenthetical names the wrong top
   category (a 17+/graduate dummy at ≈ +7.33 would give exactly 18.90, and US
   education dummy sets in this literature usually do run past Bachelor's), or
   the endpoints are right and 18.9 is a slip for 18.70. The first is likelier
   on shape, which is not evidence.

   So the caution box on the OPIE-4 tab no longer states a figure at all: it
   gives the direction of the bias and says the tool cannot quantify it, which
   is true whichever way the above resolves. DO NOT PUT A NUMBER BACK ON SCREEN
   until Table eA5.8 has been re-read — an unciteable figure in a warning is
   worse than no figure, this being a page whose whole job is to say that these
   numbers should not be quoted.

   Why they are not mapped to UK equivalents: the dummies encode how unusual a
   given attainment level is *within the US population they were fitted on*, not
   years of schooling, and that does not transfer.
     · The US reference (12th grade, leaving at 18) has no clean UK counterpart —
       matched by years it is A-levels, matched by population position it is
       GCSE/O-level. One whole category apart, worth ~3 points.
     · UK school leaving age was 14, then 15 (1947), then 16 only from 1972, so
       "no qualifications" runs from ~11% of UK 16–24s to ~45% of the over-75s.
       The US "K–7th grade" coefficient of −11.57 is large precisely because that
       group is a small, heavily selected ~5% of Americans. Applying it to a UK
       80-year-old — for whom leaving school at 14 was normative — would
       under-estimate premorbid ability by ~11 points and hide real decline.
     · UK Level 4+ includes HNC/HND, which is sub-bachelor's, so it is a broader
       and less selected group than the US Bachelor's category.

   Doing this properly would mean percentile-equating the UK cohort-specific
   qualification distribution against the US 2008 distribution — a novel
   adaptation that would need validating before it informed a report.

   For the UK demographic angle use Crawford & Allan (1997), already in this
   file's model set and UK-normed. Treat OPIE-4 as a subtest-based estimate.

   Sex coding: Female = 0, Male = 1 (per the source table). Because Female is
   coded 0, a blank Sex field would silently return the female equation — so
   app.js requires Sex before computing any OPIE-4 value.

   Fitted age range 16–90 (OPIE_AGE_MIN/MAX); the Age³ and Age⁶ terms are
   unbounded and extrapolate wildly outside it.
   ============================================================================ */
const OPIE_AGE_MIN = 16;
const OPIE_AGE_MAX = 90;

/* Age range for the ToPF-predicted WMS-IV equations, which carry an explicit
   age term (WMS_COEF[].age). Numerically the same as the OPIE range, but a
   SEPARATE constant on purpose: same number, different provenance, and one
   publisher revising its norms must not silently move the other model's bound.

   The values are the range the app has always asserted for these models — the
   #pre-age input was hard-bounded 16-90 from the outset, which was the only
   thing keeping an out-of-range age out of these equations before the age
   field became shared with Score Tables. They are the adult ToPF/WMS-IV
   standardisation range and should be confirmed against the ToPF manual if a
   paediatric or very elderly case ever needs to be defended in a report. */
const TOPF_AGE_MIN = 16;
const TOPF_AGE_MAX = 90;

/* Floor for the Crawford & Allan (1997) demographic equation — The Clinical
   Neuropsychologist, 11(2), 192-197 (long miscited here as 2001; vol. 11 is
   1997). Derived from 200 healthy adults representative of the adult UK
   population; the same lab sample is described in the authors' companion
   WAIS-R papers as ages 16-83. A floor only, and again a SEPARATE constant:
   the number matches OPIE_AGE_MIN today, but the provenance is a different
   publication and the two must be able to move independently. No ceiling
   constant exists on purpose — the equation's age term is linear and shallow
   (0.18/yr), and the sample maximum comes from companion papers rather than
   the 1997 brief report itself, so the app warns above the sample rather than
   refusing (see the range note in calcPremorbid). */
const CRAWFORD_ALLAN_AGE_MIN = 16;

// OPIE-4 prorated FSIQ regression coefficients (Table eA5.8). Verified against
// source. See the block above before using these numbers in a UK context.
const OPIE_PRORATED_FSIQ = {
  // SEE from the chapter regression tables (Holdnack et al., 2013); used for the
  // parametric prediction interval predicted ± round(z × SEE).
  VC_MR: { intercept:65.77827122, vc:0.646258435, mr:1.182068623, age:-0.197692558, age3:0.0000373292, sex:1.955504838, r:0.80, see:8.44 },
  VC:    { intercept:86.63733022, vc:0.825479066,                 age:-0.355783733, age3:0.0000373292, sex:2.795447219, r:0.71, see:9.61 },
  MR:    { intercept:62.02281403,                mr:1.719384768,                    age3:0.0000275723, sex:1.509479923, r:0.66, see:10.10 }
};

// OPIE-4 prorated GAI regression coefficients (Table eA5.8). Verified against
// source, including the MR-only branch — an earlier comment here claimed that
// column was truncated and omitted; it is present and complete in the table, and
// the coefficients below match it. R/SEE from Table 5.19 (Prorated GAI),
// post-Age step. See the block above OPIE_PRORATED_FSIQ before using in the UK.
const OPIE_PRORATED_GAI = {
  VC_MR: { intercept:60.14203956, vc:0.763136717, mr:1.127062322, age:-0.246247784, age3:0.0000416209, sex:4.708926488, r:0.79, see:8.75 },
  VC:    { intercept:79.65445374, vc:0.921039566,                 age:-0.378906405, age3:0.0000399793, sex:5.001045997, r:0.75, see:9.10 },
  MR:    { intercept:56.74797323,                mr:1.722933286,                    age3:0.0000281088, sex:3.784565031, r:0.66, see:10.66 }
};

// OPIE-4 prorated VCI / PRI (single-index) regression coefficients (Table eA5.8).
// Verified against source. VCI predicted from Vocabulary; PRI from Matrix
// Reasoning. NB: VCI is the only equation with an Age⁶ term; PRI has no linear
// Age term (Age³ only). See the block above OPIE_PRORATED_FSIQ before using in
// the UK — the omitted education block reaches −6.56 on PRI and −4.44 on VCI.
const OPIE_PRORATED_INDEX = {
  VCI: { intercept:76.77718804, vc:1.046716258, age:-0.600717374, age3:0.0000888487, age6:-4.42948e-11, sex:4.196941127, r:0.74, see:8.46 },
  PRI: { intercept:60.63732634,                mr:1.858062305,                        age3:0.0000287764,                  sex:4.162137929, r:0.62, see:12.34 }
};

const OCC_CODE = { 'Professional':1, 'Intermediate':2, 'Skilled':3, 'Semi Skilled':4, 'Unskilled':5 };

const PRE_MODEL_TOOLTIPS = {
  topfRaw: 'Assumes the ToPF word-reading raw score is a resistant estimate of premorbid ability. Input required: ToPF raw score only. The FSIQ estimate is returned from the ToPF raw-score look-up table; CI uses the model SEE.',
  topfDemo: 'Assumes premorbid FSIQ is best estimated by combining ToPF performance with demographic predictors. Inputs required: ToPF raw score, years of education, and sex. Uses the cubic ToPF + education + sex regression equation; CI uses the model SEE.',
  crawfordAllan: 'Demographic-only estimate from Crawford & Allan (1997), intended for UK-normed demographic prediction. Inputs required: occupation class, years of education, and age (16+, the equation being derived from an adult UK sample). Does not use the ToPF score.',
  opieDefault: 'ILLUSTRATIVE ONLY in a UK context - do not quote as a concrete premorbid estimate. OPIE-4 estimate of prorated FSIQ. Inputs required: age (16-90), sex, plus Vocabulary raw and/or Matrix Reasoning raw. The equation switches automatically to the matching single- or two-subtest model. The published equation also carries US education, ethnicity and region terms which are NOT applied, so every patient is scored as a US 12th-grade high-school graduate, not African-American, not resident in the US West. Those US categories have no valid UK equivalent, so no substitute adjustment is made. Expect the estimate to run high for patients who left school early and low for graduates. CI uses the branch-specific SEE.',
  opieVCMR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4 two-subtest estimate of prorated FSIQ (Holdnack et al., 2013; Table eA5.8). Inputs required: Vocabulary raw, Matrix Reasoning raw, age (16-90) and sex. Predicts a prorated FSIQ that excludes Vocabulary and Matrix Reasoning, removing part-whole correlation inflation present in the standard-FSIQ equations. US education, ethnicity and region terms are not applied - see the Vocab and/or Matrix model note.',
  opieVC: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4 single-subtest estimate of prorated FSIQ (Vocabulary branch). Inputs required: Vocabulary raw, age (16-90) and sex. US education, ethnicity and region terms are not applied; the omitted education block spans -6.68 to +5.50 points on this branch. CI uses the Vocabulary-branch SEE.',
  opieMR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4 single-subtest estimate of prorated FSIQ (Matrix Reasoning branch). Inputs required: Matrix Reasoning raw, age (16-90) and sex. The published equation uses only Age3 (no linear age term) for this branch. US education, ethnicity and region terms are not applied; the omitted education block spans -11.57 to +7.37 points on this branch, the widest of any model here. CI uses the Matrix-branch SEE.',
  predictWais: 'ToPF-predicted WAIS-IV index model. Inputs required: ToPF raw score, years of education, and sex. Difference is Achieved − Predicted; base rates are shown for negative discrepancies only.',
  predictWms: 'ToPF-predicted WMS-IV index model. Inputs required: ToPF raw score and age. Difference is Achieved − Predicted; base rates are shown for negative discrepancies only.',
  opiePredFSIQ_VCMR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated FSIQ (two-subtest). Compare against the patient\'s prorated FSIQ calculated excluding Vocabulary and Matrix Reasoning per WAIS-IV manual. Note the three FSIQ rows predict three DIFFERENT prorated criteria, so they are not interchangeable and are not expected to agree.',
  opiePredFSIQ_VC: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated FSIQ (Vocabulary only). Compare against the patient\'s prorated FSIQ calculated excluding Vocabulary per WAIS-IV manual - not against the ordinary FSIQ.',
  opiePredFSIQ_MR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated FSIQ (Matrix Reasoning only). Compare against the patient\'s prorated FSIQ calculated excluding Matrix Reasoning per WAIS-IV manual - not against the ordinary FSIQ. This branch carries the widest omitted education block (-11.57 to +7.37).',
  opiePredGAI_VCMR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated GAI (two-subtest). Compare against the patient\'s prorated GAI calculated excluding Vocabulary and Matrix Reasoning per WAIS-IV manual.',
  opiePredGAI_VC: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated GAI (Vocabulary only). Compare against the patient\'s prorated GAI calculated excluding Vocabulary per WAIS-IV manual.',
  opiePredGAI_MR: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted prorated GAI (Matrix Reasoning only). Compare against the patient\'s prorated GAI calculated excluding Matrix Reasoning per WAIS-IV manual. Omitted education block spans -9.41 to +8.10 on this branch.',
  opiePredVCI: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted VCI (Vocabulary only). Compare against the patient\'s achieved VCI. The VCI equation is the only one carrying an Age⁶ term.',
  opiePredPRI: 'ILLUSTRATIVE ONLY in a UK context. OPIE-4-predicted PRI (Matrix Reasoning only). Compare against the patient\'s achieved PRI. This equation uses Age³ only, with no linear age term.'
};

/* Discrepancy (negative integer) → proportion at or below it. WAIS-IV indices
   (FSIQ, VCI, PRI, WMI, PSI) and WMS-IV indices (IMI, DMI, VWMI).

   SOURCE: ToPF-UK manual (Wechsler, 2011), the ToPF-predicted vs obtained
   discrepancy base rates. Used as published.

   THE PUBLISHER DERIVED THEM PARAMETRICALLY, and the app says so wherever they
   appear — but it is the publisher's arithmetic, not ours, so every printed
   figure is citable. Each cell is round(Φ(d / SEE), 4): not approximately,
   exactly. Give each column its unrounded SEE and all 298 stored cells
   reproduce with no residual. A predicted-difference table has to be built this
   way — a standardisation sample of ~1,000 yields no observed frequency at all
   40 discrepancy points, let alone a monotone one.

   THAT ARITHMETIC IS ALSO THE TRANSCRIPTION PROOF, and it is what distinguishes
   this block from a downstream re-computation wearing a citation. Six of the
   eight columns fit exactly on the SEE printed in WAIS_COEF / WMS_COEF above.
   The other two fit only at MORE precision than the printed coefficient carries:

     IMI    fits only for SEE ∈ [12.03159, 12.03173] — printed value 12.032
     VWMI   fits only for SEE ∈ [12.16512, 12.16564] — printed value 12.165

   Both bands EXCLUDE the printed SEE — no value within rounding distance of it
   fits every cell. Anyone rebuilding this table from the published coefficients
   misses 4 IMI cells and 2 VWMI cells. Whoever produced it held the unrounded
   regression output, which is the publisher. check.js §8 pins the exclusion, so
   a "tidy-up" to the rounded SEE cannot pass.

   Note it is parametric where OPIE_BASE_RATES below is empirical — that one sits
   on a ~1/1020 count grid. Across the decisive −5 to −20 band the two published
   tables differ by ~10% relatively on models of almost identical SEE: d = −15
   gives 3.78% here against ~4.3% there. That is two publishers answering the
   question differently, not a defect in either, and both are used as published.

   Cells the manual prints as 0.00 are deliberately NOT stored. The lookup falls
   through to the same model continued past the table (see renderPreRow), which
   prints "< 0.01%" — true, where a stored 0.00 would assert a base rate of zero. */
const BASE_RATES = {
  '-1':{FSIQ:.4528,VCI:.4571,PRI:.4678,WMI:.4625,PSI:.4669,IMI:.4669,DMI:.4679,VWMI:.4672},
  '-2':{FSIQ:.4064,VCI:.4147,PRI:.4358,WMI:.4253,PSI:.434, IMI:.434, DMI:.436, VWMI:.4347},
  '-3':{FSIQ:.3611,VCI:.3732,PRI:.4042,WMI:.3888,PSI:.4016,IMI:.4015,DMI:.4046,VWMI:.4026},
  '-4':{FSIQ:.3178,VCI:.3332,PRI:.3732,WMI:.3532,PSI:.3698,IMI:.3698,DMI:.3737,VWMI:.3712},
  '-5':{FSIQ:.2768,VCI:.295, PRI:.343, WMI:.3188,PSI:.3389,IMI:.3389,DMI:.3436,VWMI:.3405},
  '-6':{FSIQ:.2386,VCI:.2589,PRI:.3138,WMI:.286, PSI:.3091,IMI:.309, DMI:.3145,VWMI:.3109},
  '-7':{FSIQ:.2035,VCI:.2253,PRI:.2857,WMI:.2548,PSI:.2805,IMI:.2804,DMI:.2865,VWMI:.2825},
  '-8':{FSIQ:.1716,VCI:.1942,PRI:.2589,WMI:.2256,PSI:.2532,IMI:.2531,DMI:.2597,VWMI:.2554},
  '-9':{FSIQ:.1432,VCI:.166, PRI:.2334,WMI:.1983,PSI:.2273,IMI:.2272,DMI:.2343,VWMI:.2297},
  '-10':{FSIQ:.1181,VCI:.1405,PRI:.2094,WMI:.1731,PSI:.2031,IMI:.2029,DMI:.2104,VWMI:.2055},
  '-11':{FSIQ:.0963,VCI:.1179,PRI:.1869,WMI:.1501,PSI:.1804,IMI:.1803,DMI:.1879,VWMI:.1829},
  '-12':{FSIQ:.0776,VCI:.0979,PRI:.166, WMI:.1292,PSI:.1594,IMI:.1593,DMI:.167, VWMI:.162},
  '-13':{FSIQ:.0618,VCI:.0806,PRI:.1466,WMI:.1104,PSI:.1401,IMI:.14,  DMI:.1476,VWMI:.1426},
  '-14':{FSIQ:.0486,VCI:.0656,PRI:.1288,WMI:.0936,PSI:.1224,IMI:.1223,DMI:.1298,VWMI:.1249},
  '-15':{FSIQ:.0378,VCI:.0529,PRI:.1126,WMI:.0789,PSI:.1064,IMI:.1063,DMI:.1136,VWMI:.1088},
  '-16':{FSIQ:.029, VCI:.0423,PRI:.0979,WMI:.0659,PSI:.0919,IMI:.0918,DMI:.0988,VWMI:.0942},
  '-17':{FSIQ:.022, VCI:.0334,PRI:.0846,WMI:.0547,PSI:.0789,IMI:.0788,DMI:.0855,VWMI:.0811},
  '-18':{FSIQ:.0165,VCI:.0262,PRI:.0728,WMI:.045, PSI:.0674,IMI:.0673,DMI:.0736,VWMI:.0695},
  '-19':{FSIQ:.0122,VCI:.0203,PRI:.0622,WMI:.0368,PSI:.0572,IMI:.0571,DMI:.063, VWMI:.0592},
  '-20':{FSIQ:.0089,VCI:.0155,PRI:.0529,WMI:.0298,PSI:.0483,IMI:.0482,DMI:.0537,VWMI:.0501},
  '-21':{FSIQ:.0064,VCI:.0118,PRI:.0448,WMI:.024, PSI:.0405,IMI:.0405,DMI:.0454,VWMI:.0422},
  '-22':{FSIQ:.0046,VCI:.0089,PRI:.0376,WMI:.0191,PSI:.0338,IMI:.0337,DMI:.0383,VWMI:.0353},
  '-23':{FSIQ:.0032,VCI:.0066,PRI:.0315,WMI:.0151,PSI:.028, IMI:.028, DMI:.032, VWMI:.0293},
  '-24':{FSIQ:.0022,VCI:.0048,PRI:.0262,WMI:.0119,PSI:.0231,IMI:.023, DMI:.0267,VWMI:.0243},
  '-25':{FSIQ:.0015,VCI:.0035,PRI:.0216,WMI:.0093,PSI:.0189,IMI:.0189,DMI:.0221,VWMI:.0199},
  '-26':{FSIQ:.001, VCI:.0025,PRI:.0178,WMI:.0072,PSI:.0154,IMI:.0153,DMI:.0182,VWMI:.0163},
  '-27':{FSIQ:.0007,VCI:.0018,PRI:.0145,WMI:.0055,PSI:.0125,IMI:.0124,DMI:.0149,VWMI:.0132},
  '-28':{FSIQ:.0005,VCI:.0013,PRI:.0118,WMI:.0042,PSI:.01,  IMI:.01,  DMI:.0121,VWMI:.0107},
  '-29':{FSIQ:.0003,VCI:.0009,PRI:.0095,WMI:.0032,PSI:.008, IMI:.008, DMI:.0098,VWMI:.0086},
  '-30':{FSIQ:.0002,VCI:.0006,PRI:.0076,WMI:.0024,PSI:.0063,IMI:.0063,DMI:.0079,VWMI:.0068},
  '-31':{FSIQ:.0001,VCI:.0004,PRI:.0061,WMI:.0018,PSI:.005, IMI:.005, DMI:.0063,VWMI:.0054},
  '-32':{FSIQ:.0001,VCI:.0003,PRI:.0048,WMI:.0013,PSI:.0039,IMI:.0039,DMI:.005, VWMI:.0043},
  /* VCI is published down to -36 and FSIQ only to -32; below those the manual
     prints 0.00, which is not stored. The four VCI cells here were missing while
     the table's provenance was unknown and the fallback covered them; -35 and
     -36 printed "< 0.01%" against the manual's 0.01%. */
  '-33':{          VCI:.0002, PRI:.0038,WMI:.0009,PSI:.0031,IMI:.003, DMI:.0039,VWMI:.0033},
  '-34':{          VCI:.0001, PRI:.003, WMI:.0007,PSI:.0024,IMI:.0024,DMI:.0031,VWMI:.0026},
  '-35':{          VCI:.0001, PRI:.0023,WMI:.0005,PSI:.0018,IMI:.0018,DMI:.0024,VWMI:.002},
  '-36':{          VCI:.0001, PRI:.0018,WMI:.0003,PSI:.0014,IMI:.0014,DMI:.0019,VWMI:.0015},
  '-37':{                     PRI:.0014,WMI:.0002,PSI:.0011,IMI:.0011,DMI:.0014,VWMI:.0012},
  '-38':{                     PRI:.0011,WMI:.0002,PSI:.0008,IMI:.0008,DMI:.0011,VWMI:.0009},
  '-39':{                                                   IMI:.0006,DMI:.0008,VWMI:.0007},
  '-40':{                                                   IMI:.0004,DMI:.0006,VWMI:.0005}
};

// OPIE-4 base rates from WAIS-IV / WMS-IV / ACS Table eA5.12 (prorated actual vs predicted).
// Values are stored as proportions for display via fmtPctBr.
const OPIE_BASE_RATES = {
  '-37':{FSIQ_VC:0.001},
  '-36':{FSIQ_VC:0.002},
  '-35':{FSIQ_VC:0.0049,GAI_VC:0.001},
  '-34':{FSIQ_VC:0.0049,GAI_VC:0.002,PRI:0.0005},
  '-33':{FSIQ_VC:0.0049,GAI_VC:0.002,PRI:0.0015},
  '-32':{FSIQ_VC:0.0049,GAI_VC:0.002,PRI:0.0015},
  '-31':{FSIQ_VC:0.0049,GAI_VC:0.0029,PRI:0.0025},
  '-30':{FSIQ_VC:0.0049,GAI_VC:0.0029,PRI:0.003},
  '-29':{FSIQ_VC:0.0049,GAI_VC:0.0029,PRI:0.003},
  '-28':{FSIQ_VC:0.0049,GAI_VC:0.0029,GAI_MR:0.001,PRI:0.0034},
  '-27':{FSIQ_VC:0.0049,FSIQ_MR:0.0029,GAI_VC_MR:0.001,GAI_VC:0.0039,GAI_MR:0.002,VCI:0.001,PRI:0.0044},
  '-26':{FSIQ_VC:0.0069,FSIQ_MR:0.0039,GAI_VC_MR:0.002,GAI_VC:0.0059,GAI_MR:0.0039,VCI:0.0049,PRI:0.0064},
  '-25':{FSIQ_VC_MR:0.001,FSIQ_VC:0.0079,FSIQ_MR:0.0039,GAI_VC_MR:0.002,GAI_VC:0.0059,GAI_MR:0.0088,VCI:0.0049,PRI:0.0098},
  '-24':{FSIQ_VC_MR:0.002,FSIQ_VC:0.0098,FSIQ_MR:0.0049,GAI_VC_MR:0.002,GAI_VC:0.0088,GAI_MR:0.0088,VCI:0.0049,PRI:0.0113},
  '-23':{FSIQ_VC_MR:0.002,FSIQ_VC:0.0108,FSIQ_MR:0.0098,GAI_VC_MR:0.0039,GAI_VC:0.0088,GAI_MR:0.0138,VCI:0.0049,PRI:0.0153},
  '-22':{FSIQ_VC_MR:0.0049,FSIQ_VC:0.0108,FSIQ_MR:0.0118,GAI_VC_MR:0.0059,GAI_VC:0.0108,GAI_MR:0.0138,VCI:0.0049,PRI:0.0202},
  '-21':{FSIQ_VC_MR:0.0059,FSIQ_VC:0.0177,FSIQ_MR:0.0157,GAI_VC_MR:0.0069,GAI_VC:0.0147,GAI_MR:0.0177,VCI:0.0059,PRI:0.0276},
  '-20':{FSIQ_VC_MR:0.0098,FSIQ_VC:0.0216,FSIQ_MR:0.0216,GAI_VC_MR:0.0098,GAI_VC:0.0196,GAI_MR:0.0216,VCI:0.0088,PRI:0.033},
  '-19':{FSIQ_VC_MR:0.0138,FSIQ_VC:0.0295,FSIQ_MR:0.0265,GAI_VC_MR:0.0167,GAI_VC:0.0236,GAI_MR:0.0255,VCI:0.0118,PRI:0.0389},
  '-18':{FSIQ_VC_MR:0.0187,FSIQ_VC:0.0373,FSIQ_MR:0.0432,GAI_VC_MR:0.0216,GAI_VC:0.0314,GAI_MR:0.0324,VCI:0.0187,PRI:0.0507},
  '-17':{FSIQ_VC_MR:0.0295,FSIQ_VC:0.0462,FSIQ_MR:0.0501,GAI_VC_MR:0.0314,GAI_VC:0.0354,GAI_MR:0.0452,VCI:0.0246,PRI:0.064},
  '-16':{FSIQ_VC_MR:0.0354,FSIQ_VC:0.057,FSIQ_MR:0.0619,GAI_VC_MR:0.0393,GAI_VC:0.0452,GAI_MR:0.0629,VCI:0.0314,PRI:0.0763},
  '-15':{FSIQ_VC_MR:0.0432,FSIQ_VC:0.0668,FSIQ_MR:0.0678,GAI_VC_MR:0.0472,GAI_VC:0.053,GAI_MR:0.0747,VCI:0.0413,PRI:0.0891},
  '-14':{FSIQ_VC_MR:0.0511,FSIQ_VC:0.0796,FSIQ_MR:0.0806,GAI_VC_MR:0.058,GAI_VC:0.0629,GAI_MR:0.0943,VCI:0.0511,PRI:0.1068},
  '-13':{FSIQ_VC_MR:0.0639,FSIQ_VC:0.0972,FSIQ_MR:0.0992,GAI_VC_MR:0.0697,GAI_VC:0.0756,GAI_MR:0.1169,VCI:0.0668,PRI:0.1216},
  '-12':{FSIQ_VC_MR:0.0815,FSIQ_VC:0.113,FSIQ_MR:0.1198,GAI_VC_MR:0.0855,GAI_VC:0.0914,GAI_MR:0.1297,VCI:0.0796,PRI:0.1437},
  '-11':{FSIQ_VC_MR:0.0982,FSIQ_VC:0.1277,FSIQ_MR:0.1395,GAI_VC_MR:0.1051,GAI_VC:0.1081,GAI_MR:0.1611,VCI:0.1012,PRI:0.1698},
  '-10':{FSIQ_VC_MR:0.1306,FSIQ_VC:0.1483,FSIQ_MR:0.1739,GAI_VC_MR:0.1306,GAI_VC:0.1297,GAI_MR:0.1896,VCI:0.1238,PRI:0.2013},
  '-9':{FSIQ_VC_MR:0.1591,FSIQ_VC:0.1709,FSIQ_MR:0.2004,GAI_VC_MR:0.1591,GAI_VC:0.1582,GAI_MR:0.2279,VCI:0.1473,PRI:0.2343},
  '-8':{FSIQ_VC_MR:0.1906,FSIQ_VC:0.2073,FSIQ_MR:0.2358,GAI_VC_MR:0.1945,GAI_VC:0.1916,GAI_MR:0.2525,VCI:0.1758,PRI:0.2648},
  '-7':{FSIQ_VC_MR:0.223,FSIQ_VC:0.2387,FSIQ_MR:0.2692,GAI_VC_MR:0.2279,GAI_VC:0.2269,GAI_MR:0.2898,VCI:0.2092,PRI:0.2933},
  '-6':{FSIQ_VC_MR:0.2495,FSIQ_VC:0.2741,FSIQ_MR:0.3016,GAI_VC_MR:0.2603,GAI_VC:0.2593,GAI_MR:0.3193,VCI:0.2544,PRI:0.3209},
  '-5':{FSIQ_VC_MR:0.3055,FSIQ_VC:0.3104,FSIQ_MR:0.336,GAI_VC_MR:0.3075,GAI_VC:0.3006,GAI_MR:0.3536,VCI:0.2976,PRI:0.3484},
  '-4':{FSIQ_VC_MR:0.3438,FSIQ_VC:0.3448,FSIQ_MR:0.3644,GAI_VC_MR:0.3566,GAI_VC:0.3438},
  '-3':{FSIQ_VC_MR:0.3811,FSIQ_VC:0.3772,FSIQ_MR:0.3998,GAI_VC_MR:0.3949,GAI_VC:0.389},
  '-2':{FSIQ_VC_MR:0.4322,FSIQ_VC:0.4293,FSIQ_MR:0.446,GAI_VC_MR:0.4381,GAI_VC:0.4381},
  '-1':{FSIQ_VC_MR:0.4882,FSIQ_VC:0.4725,FSIQ_MR:0.4971,GAI_VC_MR:0.4882,GAI_VC:0.4813},
  '+1':{FSIQ_VC_MR:0.4686,FSIQ_VC:0.4872,FSIQ_MR:0.4725,GAI_VC_MR:0.4617,GAI_VC:0.4656},
  '+2':{FSIQ_VC_MR:0.4273,FSIQ_VC:0.4489,FSIQ_MR:0.4303,GAI_VC_MR:0.4194,GAI_VC:0.4106},
  '+3':{FSIQ_VC_MR:0.3821,FSIQ_VC:0.3969,FSIQ_MR:0.388,GAI_VC_MR:0.3802,GAI_VC:0.3752},
  '+4':{FSIQ_VC_MR:0.3399,FSIQ_VC:0.3418,FSIQ_MR:0.3507,GAI_VC_MR:0.335,GAI_VC:0.336},
  '+5':{FSIQ_VC_MR:0.3065,FSIQ_VC:0.3075,FSIQ_MR:0.3114,GAI_VC_MR:0.2917,GAI_VC:0.2849,GAI_MR:0.333,VCI:0.2898,PRI:0.3376},
  '+6':{FSIQ_VC_MR:0.2603,FSIQ_VC:0.2721,FSIQ_MR:0.2809,GAI_VC_MR:0.2593,GAI_VC:0.2574,GAI_MR:0.2917,VCI:0.2475,PRI:0.3056},
  '+7':{FSIQ_VC_MR:0.2151,FSIQ_VC:0.2446,FSIQ_MR:0.2485,GAI_VC_MR:0.2191,GAI_VC:0.2308,GAI_MR:0.2593,VCI:0.2102,PRI:0.2781},
  '+8':{FSIQ_VC_MR:0.1857,FSIQ_VC:0.2033,FSIQ_MR:0.222,GAI_VC_MR:0.1906,GAI_VC:0.2033,GAI_MR:0.2279,VCI:0.1768,PRI:0.2451},
  '+9':{FSIQ_VC_MR:0.1532,FSIQ_VC:0.1827,FSIQ_MR:0.1916,GAI_VC_MR:0.1532,GAI_VC:0.1768,GAI_MR:0.2083,VCI:0.1375,PRI:0.2165},
  '+10':{FSIQ_VC_MR:0.1257,FSIQ_VC:0.1552,FSIQ_MR:0.168,GAI_VC_MR:0.1267,GAI_VC:0.1483,GAI_MR:0.1876,VCI:0.1159,PRI:0.1929},
  '+11':{FSIQ_VC_MR:0.0992,FSIQ_VC:0.1306,FSIQ_MR:0.1424,GAI_VC_MR:0.1022,GAI_VC:0.1277,GAI_MR:0.1572,VCI:0.0972,PRI:0.1688},
  '+12':{FSIQ_VC_MR:0.0776,FSIQ_VC:0.1081,FSIQ_MR:0.1257,GAI_VC_MR:0.0796,GAI_VC:0.1022,GAI_MR:0.1316,VCI:0.0756,PRI:0.1471},
  '+13':{FSIQ_VC_MR:0.0629,FSIQ_VC:0.0923,FSIQ_MR:0.111,GAI_VC_MR:0.0648,GAI_VC:0.0835,GAI_MR:0.1218,VCI:0.0648,PRI:0.1284},
  '+14':{FSIQ_VC_MR:0.0501,FSIQ_VC:0.0707,FSIQ_MR:0.0953,GAI_VC_MR:0.0521,GAI_VC:0.0697,GAI_MR:0.1031,VCI:0.0521,PRI:0.1112},
  '+15':{FSIQ_VC_MR:0.0373,FSIQ_VC:0.0599,FSIQ_MR:0.0806,GAI_VC_MR:0.0452,GAI_VC:0.0589,GAI_MR:0.0855,VCI:0.0452,PRI:0.0969},
  '+16':{FSIQ_VC_MR:0.0334,FSIQ_VC:0.0511,FSIQ_MR:0.0678,GAI_VC_MR:0.0363,GAI_VC:0.0472,GAI_MR:0.0727,VCI:0.0413,PRI:0.0856},
  '+17':{FSIQ_VC_MR:0.0265,FSIQ_VC:0.0413,FSIQ_MR:0.054,GAI_VC_MR:0.0314,GAI_VC:0.0363,GAI_MR:0.0619,VCI:0.0334,PRI:0.0773},
  '+18':{FSIQ_VC_MR:0.0236,FSIQ_VC:0.0373,FSIQ_MR:0.0452,GAI_VC_MR:0.0275,GAI_VC:0.0255,GAI_MR:0.055,VCI:0.0285,PRI:0.0689},
  '+19':{FSIQ_VC_MR:0.0196,FSIQ_VC:0.0265,FSIQ_MR:0.0403,GAI_VC_MR:0.0196,GAI_VC:0.0187,GAI_MR:0.0481,VCI:0.0226,PRI:0.0625},
  '+20':{FSIQ_VC_MR:0.0147,FSIQ_VC:0.0196,FSIQ_MR:0.0344,GAI_VC_MR:0.0147,GAI_VC:0.0177,GAI_MR:0.0393,VCI:0.0167,PRI:0.0507},
  '+21':{FSIQ_VC_MR:0.0088,FSIQ_VC:0.0147,FSIQ_MR:0.0265,GAI_VC_MR:0.0128,GAI_VC:0.0138,GAI_MR:0.0324,VCI:0.0118,PRI:0.0389},
  '+22':{FSIQ_VC_MR:0.0059,FSIQ_VC:0.0108,FSIQ_MR:0.0216,GAI_VC_MR:0.0088,GAI_VC:0.0079,GAI_MR:0.0285,VCI:0.0118,PRI:0.033},
  '+23':{FSIQ_VC_MR:0.0039,FSIQ_VC:0.0108,FSIQ_MR:0.0177,GAI_VC_MR:0.0088,GAI_VC:0.0079,GAI_MR:0.0255,VCI:0.0088,PRI:0.0285},
  '+24':{FSIQ_VC_MR:0.0029,FSIQ_VC:0.0108,FSIQ_MR:0.0138,GAI_VC_MR:0.0088,GAI_VC:0.0059,GAI_MR:0.0226,VCI:0.0069,PRI:0.0241},
  '+25':{FSIQ_VC_MR:0.0029,FSIQ_VC:0.0059,FSIQ_MR:0.0118,GAI_VC_MR:0.0039,GAI_VC:0.0039,GAI_MR:0.0177,VCI:0.0069,PRI:0.0182},
  '+26':{FSIQ_VC_MR:0.0029,FSIQ_VC:0.0059,FSIQ_MR:0.0079,GAI_VC_MR:0.0039,GAI_VC:0.0039,GAI_MR:0.0147,VCI:0.0059,PRI:0.0153},
  '+27':{FSIQ_VC_MR:0.0029,FSIQ_VC:0.0049,FSIQ_MR:0.0069,GAI_VC_MR:0.0039,GAI_VC:0.0029,GAI_MR:0.0108,VCI:0.0049,PRI:0.0123},
  '+28':{FSIQ_VC_MR:0.002,FSIQ_VC:0.0039,FSIQ_MR:0.0049,GAI_VC_MR:0.0039,GAI_VC:0.0029,GAI_MR:0.0098,VCI:0.0049,PRI:0.0103},
  '+29':{FSIQ_VC_MR:0.002,FSIQ_VC:0.002,FSIQ_MR:0.0039,GAI_VC_MR:0.0039,GAI_VC:0.0029,GAI_MR:0.0079,VCI:0.0029,PRI:0.0079},
  '+30':{FSIQ_VC_MR:0.002,FSIQ_VC:0.001,FSIQ_MR:0.0029,GAI_VC_MR:0.0029,GAI_VC:0.0029,GAI_MR:0.0079,VCI:0.002,PRI:0.0054},
  '+31':{FSIQ_VC_MR:0.001,FSIQ_MR:0.002,GAI_VC_MR:0.0029,GAI_VC:0.002,GAI_MR:0.0069,VCI:0.002,PRI:0.003},
  '+32':{FSIQ_MR:0.002,GAI_VC_MR:0.0029,GAI_MR:0.0059,VCI:0.002,PRI:0.0025},
  '+33':{FSIQ_MR:0.001,GAI_VC_MR:0.002,GAI_MR:0.0049,VCI:0.002,PRI:0.0025},
  '+34':{GAI_VC_MR:0.002,GAI_MR:0.0039,VCI:0.002,PRI:0.0025},
  '+35':{GAI_VC_MR:0.002,GAI_MR:0.0039,VCI:0.001,PRI:0.0025},
  '+36':{GAI_VC_MR:0.002,GAI_MR:0.0039,PRI:0.002},
  '+37':{GAI_VC_MR:0.002,GAI_MR:0.0039,PRI:0.002},
  '+38':{GAI_VC_MR:0.001,GAI_MR:0.0039,PRI:0.002},
  '+39':{GAI_MR:0.0039,PRI:0.0015},
  '+40':{GAI_MR:0.0039,PRI:0.0015}
};

/* ============================================================
   NORMATIVE DATABASE
   Reliability parameters from published manuals/literature.
   m1/sd1 = first administration, m2/sd2 = retest (or alternate form)
   r = raw reliability, rCorrected = corrected reliability, n = sample size
   Reliability TYPE varies by test (same-form test–retest vs alternate-form).
   e.g. CVLT-3 = alternate-form (Standard ↔ Alternate Form), Manual Table 3.4.
        RBANS  = same-form test–retest (Form A → Form A), Manual Tables 3.8–3.9.

   metric:'raw'
   -------------
   Present on 63 entries across five groups. Marks a measure that is NOT on a
   standardised metric, so no percentile or classification can be derived from
   the score alone.

   Without it the score type had to be guessed from the normative mean
   (inferScoreTypeForSubtest), which has no concept of a raw score and typed
   these as scaled/T/standard/z. On Score Tables that produced badly wrong
   output in both directions: RBANS List Recognition (M 19.6, SD 0.8, raw out
   of 20) at a raw 19 is the 23rd percentile, but was read as a scaled score
   and printed as the 99.9th, "Very Superior" — 56 standard-score points out,
   with the sign inverted.

   The five groups, and why they are certainly raw:
     CVLT-C Subtests (Raw Scores) · Age 8 / 12 / 16   — declared in the name
     RBANS Subtests · Ages 12-19 / 20-89              — Picture Naming has
       SD 0.4 and List Recognition SD 0.8; no standardised metric in this app
       has an SD below 1, so these cannot be scaled scores.

   Consumers must degrade gracefully: toZ returns null for 'raw', so the
   percentile and classification cells stay blank rather than being invented.
   The reliability parameters themselves are unaffected and remain correct —
   the Change Analysis pages use m1/sd1/m2/sd2/r directly and never needed a
   metric, which is why they were never wrong.

   ------------------------------------------------------------
   rInternal — an INTERNAL-CONSISTENCY coefficient, for the CI only

   `r` is a stability coefficient: the same people retested. `rInternal` is a
   consistency coefficient: one administration, split in half. They answer
   different questions and are NOT interchangeable.

   The field exists only where a publisher does two things: reports an
   internal-consistency coefficient, AND derives its own published confidence
   intervals from it. Four manuals clear that bar, and each for a different
   reason — reliability method is a per-manual question, never a policy this
   app applies across instruments:

     CVLT-C          List A Trials 1-5 Total only, 3 entries. rInternal, no
                     by-age table. Details below.
     D-KEFS Advanced 26 entries on the All Ages groups, rInternal (published
                     Average) + rInternalByAge. TMT and VFT excluded, because
                     that manual rejects internal consistency for SPEEDED
                     measures.
     D-KEFS          13 entries on the All Ages groups, rInternalByAge ONLY —
                     it publishes no all-ages average. See below.
     WAIS-IV         21 entries on the All Ages groups, rInternal +
                     rInternalByAge, from Technical Manual Table 4.1. Symbol
                     Search, Coding and Cancellation excluded — that manual
                     rejects internal consistency for SPEEDED measures too, and
                     those three already hold its chosen coefficient. See below.

     CVLT-C Manual, Table 6.5 — odd/even split-half with Spearman-Brown,
     by age. Age 8 .87, Age 12 .89, Age 16 .84. The same table prints the
     SEM, and every one of its 12 age bands equals 10 x sqrt(1 - rxx), i.e.
     it is computed on the T metric, not on the raw SD. The manual states
     the formula (SEM = SD sqrt(1 - rxx), p. 87) and instructs the reader to
     add and subtract from "the child's standard T score".

     Its worked example settles it. An 8-year-old with T = 45 is given a 95%
     interval of 38-52 and a 90% interval of 39-51. From rxx .87 and SD 10:
     SEM 3.606, so 45 +/- 7 and 45 +/- 6. Both reproduce exactly.

   ------------------------------------------------------------
   D-KEFS (original) — rInternalByAge with NO rInternal

   D-KEFS Technical Manual, Chapter 2 (Evidence of Reliability), pp. 18-19:
   "For some D-KEFS measures, the internal consistency coefficients were used;
   for other measures the test-retest correlations were employed." It prints
   SEM = SD sqrt(1 - rxx) with "The standard deviation unit is 3 for all
   D-KEFS scaled scores", and CI = observed +/- z(SEM) at 90% and 95%.

   So the manual runs TWO regimes, and which applies is a per-test fact:
     - internal consistency, tabulated BY AGE (Tables 2.1/2.4/2.9/2.12/2.15/
       2.18/2.21/2.24, with the matching SEM tables);
     - the retest r on the TOTAL sample, where no internal-consistency table
       exists, because "the samples for the test-retest studies were too small
       when analyzed by age band".

   Proof of which coefficient feeds the SEM, over all 216 published cells:
   3 x sqrt(1 - internal consistency) reproduces 190 exactly and all but one
   measure to within 0.006 of the printed value (last-digit rounding — the
   manual computes from unrounded coefficients). The retest r reproduces
   2 of 200. Not ambiguous.

   NO rInternal FALLBACK EXISTS, and none may be invented. Unlike D-KEFS
   Advanced, not one of the eight tables prints an all-ages average. A blank
   or out-of-range age therefore falls through rInternalForAge to the retest
   `r`, which is exactly the manual's own second regime — confirmed by Design
   Fluency Table 2.8, whose three All Ages SEMs all reproduce as
   3 x sqrt(1 - r). Both paths are the publisher's own figures, so the column
   stays citable whether or not an age is entered.

   WHICH MEASURES, and why the others get nothing:
     Verbal Fluency (4), Sorting (3), Word Context, Tower, Word Proverb  — carry it.
     Trail Making — ONLY "Combined Number + Letter". Table 2.1 is the composite
       score alone; the five conditions have no internal-consistency table.
     Colour-Word Interference — NONE. Its only such table (2.9) is for the
       Combined Colour Naming + Word Reading composite, which is not a measure
       in this database. Attaching it to Colour Naming or Word Reading would
       give a condition the composite's reliability.
     Design Fluency — NONE. "Item interdependence precluded the use of internal
       consistency procedures, and reliability was investigated with test-retest
       procedures."
     Twenty Questions — BOTH measures carry it, from Table 2.15, but Initial
       Abstraction comes with a caveat the other 12 measures do not. Read on
       before touching its numbers.

   Word Proverb is banded 16-19 upward — the Proverb Test is not administered
   under 16 — so a child's age finds no band and correctly falls back.

   ------------------------------------------------------------
   Twenty Questions Initial Abstraction — the manual contradicts itself

   THE STORED VALUES ARE TABLE 2.15, VERBATIM. Do not "fix" them.

   Everywhere else in this manual, Table 2.15-style coefficients and their SEM
   table satisfy the manual's own formula, SEM = 3 sqrt(1 - rxx). For this one
   measure they do not:

     ages 40-49   Table 2.15 prints rxx .75
                  3 sqrt(1 - .75) = 1.50
                  Table 2.17 prints SEM 1.76, which needs rxx .656

   12 of the 16 bands disagree, by up to .094, in BOTH directions — so it is
   not a transform anyone has failed to spot. Ruled out by test: every other
   published coefficient column, every row shift, any constant SD other than 3,
   and the Spearman-Brown-uncorrected split-half. The column beside it in the
   same table, Total Weighted Achievement, reproduces 16 of 16. One of the two
   tables is misprinted and nothing inside the document says which — there is
   no worked example or third statistic for this test to triangulate on.

   WHY 2.15 AND NOT 2.17. This app stores coefficients and derives the SEM, so
   there is no field to put a published SEM in. Honouring 2.17 would mean
   storing .656 — a number the manual never prints — which is the one thing
   the project forbids. 2.15 is the publisher's own printed coefficient put
   through the publisher's own printed formula, so every stored digit is
   citable.

   WHAT IT COSTS. Three bands print a margin one scaled-score point away from
   Table 2.17's: ages 11 and 12 (+/-2 here, +/-3 there) and 80-89 (+/-3 here,
   +/-2 there). Note the bands with the LARGEST coefficient gaps — 40-49 at
   .094 and 30-39 at .069 — are not among them: both readings round to +/-3.
   Severity in the table and visibility on screen are nearly unrelated.

   check.js section 27 pins the discrepancy itself, asserting that Table 2.17
   still fails to reproduce. If a corrected printing ever makes it reconcile,
   that check fires and tells you to re-read this note.
   ------------------------------------------------------------

   Do NOT confuse this with D-KEFS Advanced, whose manual reaches the opposite
   conclusion on the same two test names: it excludes Trail Making and Verbal
   Fluency for speededness, while this manual publishes internal consistency
   for all four Verbal Fluency measures and excludes Design Fluency instead.
   ------------------------------------------------------------

   ------------------------------------------------------------
   WAIS-IV — Table 4.1, and the three subtests that are NOT in it

   WAIS-IV Technical and Interpretive Manual (GB), chapter 4, "Reliability and
   Errors of Measurement", p. 42:

     "Reliability coefficients were obtained utilizing the split-half and the
     Cronbach's coefficient alpha methods... Because Symbol Search, Coding, and
     Cancellation are Processing Speed subtests, the split-half coefficient is
     not a proper reliability estimate. Therefore, test-retest stability
     coefficients were used as the reliability estimates for these subtests...
     The stability coefficient is the correlation between the scores on the
     first and second testings corrected for the normative sample's variability
     (Allen & Yen, 1979; Magnusson, 1967)."

   Table 4.1 is that chapter's reliability table: 24 rows (15 subtests, 4
   process scores, 5 composites) x 13 NORMATIVE age bands plus an overall
   average computed with Fisher's z. 21 of the 24 rows are internal
   consistency and are stored here as rInternal + rInternalByAge.

   THE THREE SPEEDED SUBTESTS NEED NO NEW FIELD — they already hold Table 4.1's
   own coefficient. Table 4.1 carries the corrected stability coefficient for
   Symbol Search, Coding and Cancellation, and that is exactly what rCorrected
   holds: all 38 published cells — 35 banded plus the three overall averages —
   match rCorrected to the digit, and not one matches the raw r. No other
   measure exceeds 3 of 13. Giving them rInternal would therefore change
   nothing numerically while falsely labelling a stability coefficient as
   internal consistency — so they are excluded, as D-KEFS Advanced excludes
   TMT and VFT for the same reason. check.js section 28 pins the 38 cells.

   That match is also what proves the transcription: Table 4.1 propagates one
   retest-study value across the normative bands it spans (Symbol Search .81 at
   16-17/18-19/20-24/25-29, .73 at 30-34/35-44/45-54), and those spans
   reconstruct this database's 16-29 / 30-54 / 55-69 / 70-90 groups exactly.

   PSI AND FSIQ ARE HYBRIDS, and the manual says so (p. 43): the PSI average
   "is based on test-retest subtest reliabilities, which tend to be lower than
   split-half or alpha reliabilities". The manual nonetheless computes and
   labels every composite coefficient as internal consistency (p. 42, "The
   composite score internal consistency reliability coefficients were
   calculated with the formula recommended by Guilford (1954)"), and Table 4.1
   is the only reliability table it publishes for them. They are stored as the
   other composites are; excluding PSI while keeping FSIQ — which also draws on
   Processing Speed subtests — would be incoherent.

   LETTER-NUMBER SEQUENCING AND FIGURE WEIGHTS STOP AT 69, because the manual
   prints a dash above the 65-69 band for both (they are normed to 69 only, as
   is Cancellation). Their rInternalAgeMax is 69, not 90, so an age of 70+
   falls back to the published overall average rather than silently re-reading
   the 65-69 band. That average is computed over ages 16-69, which is the only
   sample there is, so it stays citable.

   TABLE 4.3 CONFIRMS THE LOT, ARITHMETICALLY. This paragraph used to record
   that the manual builds its own SEMs from Table 4.1 as unproven — the claim
   rested on the manual's stated method, Table 4.3 not being available. It is
   available now, and every one of its 300 published cells equals
   populationSD * sqrt(1 - the coefficient stored here), exactly, at the
   printed 2dp. check.js section 28 pins all 300.

   Three things fall out of it, beyond the confirmation itself.

   The SD question is settled by a second publisher. Table 4.3's note says the
   SEMs use "the reliability coefficients shown in Table 4.1 and the population
   standard deviations (i.e., 3 for the subtests and 15 for the composite
   scores)" — which is the rule section 23 derived independently from the
   CVLT-3 manual's SEM column, and the pairing 6de0d81 replaced. Reproductions
   over the 300 cells: this app's pairing 300, the retest r 5, rCorrected 46
   (the 38 speeded cells it legitimately owns, plus 8 coincidences), and sd1
   against the normatively-scaled coefficient — the old pairing — 29.

   The speeded three are confirmed twice over, through a different table: their
   Table 4.3 cells reproduce from rCorrected, as their Table 4.1 cells already
   did 38/38.

   The blank-age fallback is now pinned, and was not before. rInternal is what
   a patient with no age entered is scored on, and nothing held it: the spot
   checks read rInternalByAge only. Confirmed by mutation — moving Vocabulary's
   average .94 -> .95 passed every check in the file. Table 4.3 fixes the 13
   band coefficients exactly, and Table 4.1's average column is computed with
   Fisher's z, so the average is now derivable from pinned values: Fisher's z
   over the bands reproduces all 24 stored averages at 2dp, where the plain
   arithmetic mean manages 18.

   TABLE 4.1 AND TABLE 4.3 AVERAGE THE SAME 13 NUMBERS BY DIFFERENT METHODS,
   so nothing reconciles their two average columns. Table 4.3's footnote a
   defines the average SEM as the RMS of the band SEMs, and
     RMS of SEMs = sqrt(mean(SD^2 (1 - r_i))) = SD * sqrt(1 - mean(r_i))
   which is the ARITHMETIC mean of the band coefficients through the SEM
   formula — it reproduces 20 of 24, against 5 for the Fisher-z average that
   Table 4.1's own column uses. A property of the source, not of this app.

   The blank-age interval therefore keeps coming from Table 4.1's average
   COEFFICIENT: that is what a patient of unknown band needs, and it is
   published. Backing a coefficient out of Table 4.3's average SEM instead
   would mean storing 1 - (2.16/15)^2 = .9793 for FSIQ, a number the manual
   never prints — the Twenty Questions rule.

   The cost is two printed margins in 48 (24 measures x the two offered CI
   levels), both on a rounding knife-edge: FSIQ at 90% prints +/-3 here against
   +/-4 from the average SEM (3.49 vs 3.55), and Digit Span Backward at 95%
   +/-2 against +/-3 (2.49 vs 2.55).

   AND THE GAP IS MOSTLY ROUNDING, NOT METHOD — the obvious reading is wrong.
   Decomposing FSIQ:
     arithmetic mean of the bands  .97923  -> 2.1617   (= Table 4.3's 2.16)
     Fisher-z average, unrounded   .97936  -> 2.1547
     the stored published average  .98     -> 2.1213   (what the app prints)
   The averaging method costs 0.007; rounding the published average to 2dp
   costs 0.033, roughly five times more. On VCI it is 0.044 against 0.19.

   THE THIRD OPTION, WEIGHED AND DECLINED. Using the Fisher-z average of the
   bands UNROUNDED would match the manual's printed margin on 24 of 24 at both
   CI levels, and would invent nothing — the bands are pinned exactly by Table
   4.3, and Fisher's z is Table 4.1's own stated method. It was declined because
   this app renders the coefficient it actually used, so the reliability cell
   would print .979 where the manual prints .98, and a clinician checking
   against Table 4.1 would find a number that is not there. Reviewed and kept
   as-is, 2026-07-31. Entering an age avoids the question entirely, and that
   path is exact.

   LEAVE WAIS-IV Longest Span (Process) ALONE. It is scored by published base
   rate, has no retest data at all, and appears nowhere in Table 4.1.
   ------------------------------------------------------------

   WHERE IT MUST NOT BE USED: Change Analysis and the RCI pages. There the
   reliability is not merely an error term — with the retest r, SD sqrt(2(1-r))
   IS the spread of the difference distribution, which is the quantity reliable
   change is measured against. An internal-consistency coefficient strips out
   real between-session fluctuation, narrows the interval and overcalls change.
   On SRB and Crawford it is worse still, r being a fitted regression slope.
   check.js section 24 asserts that no RCI path reads this field.
   ------------------------------------------------------------
   ============================================================ */
﻿normDB = {
  "CVLT-3 Indices · Ages 16-44": {
    "T1-5 Correct": { m1:97.9, sd1:15.8, m2:98.5, sd2:14, r:0.77, rCorrected:0.75, n:91 },
    "Delayed Recall Correct": { m1:97.9, sd1:15.3, m2:98.4, sd2:15.7, r:0.78, rCorrected:0.77, n:91 },
    "Total Recall Correct": { m1:97.7, sd1:16, m2:98.7, sd2:14.8, r:0.82, rCorrected:0.79, n:91 }
  },
  "CVLT-3 Indices · Ages 45-90": {
    "T1-5 Correct": { m1:102.6, sd1:14.3, m2:101.3, sd2:14.6, r:0.77, rCorrected:0.8, n:122 },
    "Delayed Recall Correct": { m1:100.9, sd1:15.4, m2:101.1, sd2:15.7, r:0.81, rCorrected:0.8, n:122 },
    "Total Recall Correct": { m1:101.7, sd1:15.1, m2:101.2, sd2:15.5, r:0.84, rCorrected:0.83, n:122 }
  },
  "CVLT-3 Trials · Ages 16-44": {
    "Trial 1": { m1:9.6, sd1:3, m2:9.7, sd2:3, r:0.49, rCorrected:0.51, n:91 },
    "Trial 2": { m1:9.8, sd1:3.1, m2:10, sd2:2.6, r:0.56, rCorrected:0.53, n:91 },
    "Trial 3": { m1:9.7, sd1:3, m2:9.7, sd2:3, r:0.69, rCorrected:0.7, n:91 },
    "Trial 4": { m1:9.8, sd1:3, m2:10, sd2:3, r:0.57, rCorrected:0.58, n:91 },
    "Trial 5": { m1:9.3, sd1:3.1, m2:9.1, sd2:2.9, r:0.53, rCorrected:0.5, n:91 },
    "List B Correct": { m1:9.5, sd1:2.7, m2:9.8, sd2:2.8, r:0.45, rCorrected:0.57, n:91 },
    "Short Delay Free Recall": { m1:9.8, sd1:3, m2:9.6, sd2:2.9, r:0.7, rCorrected:0.71, n:91 },
    "Short Delay Cued Recall": { m1:9.8, sd1:3.2, m2:9.6, sd2:3.2, r:0.61, rCorrected:0.56, n:91 },
    "Long Delay Free Recall": { m1:9.6, sd1:3.2, m2:9.5, sd2:3.3, r:0.65, rCorrected:0.61, n:91 },
    "Long Delay Cued Recall": { m1:9.6, sd1:3.2, m2:9.6, sd2:2.8, r:0.65, rCorrected:0.6, n:91 },
    "Recognition": { m1:10, sd1:3.1, m2:9.8, sd2:2.9, r:0.51, rCorrected:0.49, n:91 },
    "Recognition False Positive": { m1:9.5, sd1:2.8, m2:9.1, sd2:2.9, r:0.5, rCorrected:0.57, n:91 },
    "Recognition Discrimination": { m1:9.6, sd1:2.8, m2:9.3, sd2:2.7, r:0.61, rCorrected:0.67, n:91 },
    "Discrimination Nonparametric": { m1:9.8, sd1:2.7, m2:9.4, sd2:2.6, r:0.57, rCorrected:0.65, n:91 },
    "Total Intrusions": { m1:9.8, sd1:3.2, m2:9.7, sd2:2.7, r:0.33, rCorrected:0.25, n:91 },
    "Total Repetitions": { m1:10.2, sd1:3.3, m2:10.1, sd2:2.9, r:0.65, rCorrected:0.59, n:91 }
  },
  "CVLT-3 Trials · Ages 45-90": {
    "Trial 1": { m1:10.3, sd1:2.6, m2:10, sd2:3, r:0.42, rCorrected:0.56, n:122 },
    "Trial 2": { m1:10.5, sd1:2.9, m2:10.3, sd2:3.2, r:0.55, rCorrected:0.59, n:122 },
    "Trial 3": { m1:10.7, sd1:3, m2:10.5, sd2:2.8, r:0.62, rCorrected:0.63, n:122 },
    "Trial 4": { m1:10.5, sd1:3, m2:10.3, sd2:2.6, r:0.67, rCorrected:0.67, n:122 },
    "Trial 5": { m1:10.5, sd1:2.8, m2:10.2, sd2:2.9, r:0.67, rCorrected:0.71, n:122 },
    "List B Correct": { m1:10.1, sd1:2.9, m2:10.1, sd2:3, r:0.48, rCorrected:0.53, n:122 },
    "Short Delay Free Recall": { m1:10.3, sd1:3.2, m2:10.5, sd2:3.2, r:0.74, rCorrected:0.71, n:122 },
    "Short Delay Cued Recall": { m1:10.1, sd1:3.1, m2:10.2, sd2:3.1, r:0.7, rCorrected:0.69, n:122 },
    "Long Delay Free Recall": { m1:10.2, sd1:3.1, m2:10.2, sd2:3.1, r:0.73, rCorrected:0.71, n:122 },
    "Long Delay Cued Recall": { m1:10.1, sd1:3.2, m2:9.9, sd2:3, r:0.69, rCorrected:0.65, n:122 },
    "Recognition": { m1:10.4, sd1:3, m2:10.5, sd2:3, r:0.61, rCorrected:0.62, n:122 },
    "Recognition False Positive": { m1:10.1, sd1:3.2, m2:10.3, sd2:3.2, r:0.59, rCorrected:0.55, n:122 },
    "Recognition Discrimination": { m1:10.2, sd1:3.5, m2:10.3, sd2:3.2, r:0.65, rCorrected:0.54, n:122 },
    "Discrimination Nonparametric": { m1:10.2, sd1:3.4, m2:10.3, sd2:3.3, r:0.68, rCorrected:0.58, n:122 },
    "Total Intrusions": { m1:9.4, sd1:3, m2:9.7, sd2:3.3, r:0.56, rCorrected:0.57, n:122 },
    "Total Repetitions": { m1:9.7, sd1:3.1, m2:10, sd2:3.4, r:0.72, rCorrected:0.7, n:122 }
  },
  "CVLT-C Subtests (Raw Scores) · Age 8": {
    "List A Trials 1-5 Total": { m1:42.47, sd1:9.55, m2:48.32, sd2:9.09, r:0.73, n:35, metric:'raw', reportedAs:'t', rInternal:0.87 },
    "List B Free-Recall Trial": { m1:5.26, sd1:1.71, m2:5.39, sd2:2.16, r:0.59, n:35, metric:'raw', reportedAs:'z' },
    "Short-Delay Free Recall": { m1:8.5, sd1:2.04, m2:9.91, sd2:2.37, r:0.4, n:35, metric:'raw', reportedAs:'z' },
    "Short-Delay Cued Recall": { m1:8.35, sd1:2.07, m2:9.51, sd2:2.79, r:0.75, n:35, metric:'raw', reportedAs:'z' },
    "Long-Delay Free Recall": { m1:8.94, sd1:2.3, m2:10.14, sd2:2.59, r:0.59, n:35, metric:'raw', reportedAs:'z' },
    "Long-Delay Cued Recall": { m1:8.42, sd1:2.21, m2:10.15, sd2:3.09, r:0.69, n:35, metric:'raw', reportedAs:'z' },
    "Semantic Cluster Ratio": { m1:1.44, sd1:0.44, m2:1.83, sd2:0.63, r:0.56, n:35, metric:'raw', reportedAs:'z' },
    "Perseverations": { m1:7.35, sd1:6.94, m2:7.54, sd2:6.39, r:0.9, n:35, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Free-Recall Intrusions": { m1:5.11, sd1:5.55, m2:4.97, sd2:5.99, r:0.74, n:35, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Cued-Recall Intrusions": { m1:2.83, sd1:4.33, m2:2.31, sd2:3.31, r:0.59, n:35, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Recognition Hits": { m1:13.67, sd1:1.45, m2:13.84, sd2:1.33, r:0.38, n:35, metric:'raw', reportedAs:'z' },
    "Discriminability": { m1:92.76, sd1:5.69, m2:94.52, sd2:5.9, r:0.55, n:35, metric:'raw', reportedAs:'z' },
    "False Positives": { m1:1.82, sd1:2.38, m2:1.16, sd2:2, r:0.62, n:35, metric:'raw', reportedAs:'z', higherIsWorse:true }
  },
  "CVLT-C Subtests (Raw Scores) · Age 12": {
    "List A Trials 1-5 Total": { m1:50.64, sd1:7.19, m2:56.47, sd2:8.92, r:0.73, n:40, metric:'raw', reportedAs:'t', rInternal:0.89 },
    "List B Free-Recall Trial": { m1:6.24, sd1:1.58, m2:6.82, sd2:1.57, r:0.26, n:40, metric:'raw', reportedAs:'z' },
    "Short-Delay Free Recall": { m1:9.89, sd1:2.44, m2:12.24, sd2:2.62, r:0.77, n:40, metric:'raw', reportedAs:'z' },
    "Short-Delay Cued Recall": { m1:11.13, sd1:2.21, m2:12.15, sd2:3.14, r:0.49, n:40, metric:'raw', reportedAs:'z' },
    "Long-Delay Free Recall": { m1:10.77, sd1:2.32, m2:12.3, sd2:2.54, r:0.62, n:40, metric:'raw', reportedAs:'z' },
    "Long-Delay Cued Recall": { m1:11.26, sd1:2.18, m2:12.74, sd2:2.53, r:0.69, n:40, metric:'raw', reportedAs:'z' },
    "Semantic Cluster Ratio": { m1:1.49, sd1:0.43, m2:1.83, sd2:0.6, r:0.58, n:40, metric:'raw', reportedAs:'z' },
    "Perseverations": { m1:4.55, sd1:4.32, m2:5.6, sd2:5.7, r:0.32, n:40, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Free-Recall Intrusions": { m1:0.97, sd1:1.44, m2:1.64, sd2:2.67, r:0.56, n:40, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Cued-Recall Intrusions": { m1:0.69, sd1:1.14, m2:0.56, sd2:1.02, r:0.17, n:40, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Recognition Hits": { m1:14.21, sd1:1.02, m2:14.67, sd2:0.82, r:0.24, n:40, metric:'raw', reportedAs:'z' },
    "Discriminability": { m1:96.71, sd1:3.53, m2:98.06, sd2:3.25, r:0.37, n:40, metric:'raw', reportedAs:'z' },
    "False Positives": { m1:0.55, sd1:0.86, m2:0.55, sd2:1.03, r:0.35, n:40, metric:'raw', reportedAs:'z', higherIsWorse:true }
  },
  "CVLT-C Subtests (Raw Scores) · Age 16": {
    "List A Trials 1-5 Total": { m1:53.53, sd1:6.15, m2:62.94, sd2:10.94, r:0.61, n:31, metric:'raw', reportedAs:'t', rInternal:0.84 },
    "List B Free-Recall Trial": { m1:6.67, sd1:1.7, m2:7.6, sd2:2.19, r:0.66, n:31, metric:'raw', reportedAs:'z' },
    "Short-Delay Free Recall": { m1:11.61, sd1:2.09, m2:13.52, sd2:1.91, r:0.48, n:31, metric:'raw', reportedAs:'z' },
    "Short-Delay Cued Recall": { m1:12, sd1:1.52, m2:13.71, sd2:1.59, r:0.59, n:31, metric:'raw', reportedAs:'z' },
    "Long-Delay Free Recall": { m1:11.9, sd1:1.9, m2:13.5, sd2:1.87, r:0.6, n:31, metric:'raw', reportedAs:'z' },
    "Long-Delay Cued Recall": { m1:12.57, sd1:1.72, m2:14, sd2:1.58, r:0.59, n:31, metric:'raw', reportedAs:'z' },
    "Semantic Cluster Ratio": { m1:1.55, sd1:0.53, m2:2.3, sd2:0.66, r:0.53, n:31, metric:'raw', reportedAs:'z' },
    "Perseverations": { m1:3.84, sd1:3.03, m2:5.42, sd2:5.76, r:0.31, n:31, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Free-Recall Intrusions": { m1:2.28, sd1:4.22, m2:2.19, sd2:4.13, r:0.85, n:31, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Cued-Recall Intrusions": { m1:0.66, sd1:1.68, m2:0.84, sd2:1.59, r:0.74, n:31, metric:'raw', reportedAs:'z', higherIsWorse:true },
    "Recognition Hits": { m1:14.59, sd1:0.82, m2:14.79, sd2:0.62, r:0.8, n:31, metric:'raw', reportedAs:'z' },
    "Discriminability": { m1:97.2, sd1:4.17, m2:98.9, sd2:2.23, r:0.78, n:31, metric:'raw', reportedAs:'z' },
    "False Positives": { m1:0.68, sd1:1.51, m2:0.35, sd2:0.91, r:0.78, n:31, metric:'raw', reportedAs:'z', higherIsWorse:true }
  },
  "D-KEFS Colour-Word Interference · Ages 20-49": {
    "Colour Naming": { m1:9.63, sd1:3.15, m2:10.6, sd2:2.58, r:0.86, n:35 },
    "Word Reading": { m1:9.57, sd1:2.93, m2:10.17, sd2:2.15, r:0.49, n:35 },
    "Inhibition": { m1:10.11, sd1:2.63, m2:11.29, sd2:2.19, r:0.71, n:35 },
    "Inhibition/Switching": { m1:10, sd1:2.36, m2:11.09, sd2:2.15, r:0.52, n:35 }
  },
  "D-KEFS Colour-Word Interference · Ages 50-89": {
    "Colour Naming": { m1:9.63, sd1:3.12, m2:10.16, sd2:3.37, r:0.56, n:38 },
    "Word Reading": { m1:9.95, sd1:3, m2:10.4, sd2:2.67, r:0.56, n:38 },
    "Inhibition": { m1:10.43, sd1:3.07, m2:10.97, sd2:3.48, r:0.5, n:37 },
    "Inhibition/Switching": { m1:10.43, sd1:2.68, m2:10.92, sd2:3.48, r:0.57, n:37 }
  },
  "D-KEFS Colour-Word Interference · Ages 8-19": {
    "Colour Naming": { m1:9.96, sd1:2.43, m2:11.04, sd2:2.76, r:0.79, n:28 },
    "Word Reading": { m1:10.04, sd1:2.82, m2:10.04, sd2:3.6, r:0.77, n:28 },
    "Inhibition": { m1:10.07, sd1:3.01, m2:11.54, sd2:2.78, r:0.9, n:28 },
    "Inhibition/Switching": { m1:9.75, sd1:2.94, m2:11.57, sd2:3.25, r:0.8, n:28 }
  },
  "D-KEFS Colour-Word Interference · All Ages": {
    "Colour Naming": { m1:9.72, sd1:2.93, m2:10.55, sd2:2.94, r:0.76, n:101 },
    "Word Reading": { m1:9.84, sd1:2.91, m2:10.22, sd2:2.78, r:0.62, n:101 },
    "Inhibition": { m1:10.22, sd1:2.88, m2:11.24, sd2:2.87, r:0.75, n:100 },
    "Inhibition/Switching": { m1:10.09, sd1:2.64, m2:11.16, sd2:2.99, r:0.65, n:100 }
  },
  "D-KEFS Design Fluency · Ages 20-49": {
    "Filled Dots": { m1:9.37, sd1:2.76, m2:11.89, sd2:3.39, r:0.62, n:35 },
    "Empty Dots": { m1:9.37, sd1:2.68, m2:11.14, sd2:2.79, r:0.73, n:35 },
    "Switching": { m1:9.83, sd1:3.44, m2:11.49, sd2:2.83, r:0.22, n:35 }
  },
  "D-KEFS Design Fluency · Ages 50-89": {
    "Filled Dots": { m1:10.24, sd1:2.72, m2:11.84, sd2:3.06, r:0.43, n:38 },
    "Empty Dots": { m1:10, sd1:2.97, m2:11.08, sd2:2.91, r:0.49, n:38 },
    "Switching": { m1:11.03, sd1:2.59, m2:11.16, sd2:3.33, r:0.58, n:38 }
  },
  "D-KEFS Design Fluency · Ages 8-19": {
    "Filled Dots": { m1:10.21, sd1:2.74, m2:11.75, sd2:3.19, r:0.66, n:28 },
    "Empty Dots": { m1:9.64, sd1:3.38, m2:11.39, sd2:3.1, r:0.43, n:28 },
    "Switching": { m1:9.64, sd1:2.56, m2:11.86, sd2:2.81, r:0.13, n:28 }
  },
  "D-KEFS Design Fluency · All Ages": {
    "Filled Dots": { m1:9.93, sd1:2.74, m2:11.83, sd2:3.19, r:0.58, n:101 },
    "Empty Dots": { m1:9.68, sd1:2.98, m2:11.19, sd2:2.89, r:0.57, n:101 },
    "Switching": { m1:10.23, sd1:2.95, m2:11.47, sd2:3.01, r:0.32, n:101 }
  },
  "D-KEFS Sorting Test · Ages 20-49": {
    "Free Sorting Confirmed Sorts": { m1:9.67, sd1:3.24, m2:11.33, sd2:2.56, r:0.51, n:35 },
    "Free Sorting Description Total Score": { m1:9.45, sd1:3.23, m2:11.3, sd2:2.57, r:0.46, n:35 },
    "Sort Recognition Total Description Score": { m1:9.91, sd1:2.83, m2:10.91, sd2:3.48, r:0.55, n:35 }
  },
  "D-KEFS Sorting Test · Ages 50-89": {
    "Free Sorting Confirmed Sorts": { m1:10.9, sd1:2.98, m2:10.97, sd2:2.44, r:0.62, n:38 },
    "Free Sorting Description Total Score": { m1:10.77, sd1:2.82, m2:10.7, sd2:2.44, r:0.63, n:38 },
    "Sort Recognition Total Description Score": { m1:10.67, sd1:3.02, m2:10.83, sd2:3.03, r:0.73, n:38 }
  },
  "D-KEFS Sorting Test · Ages 8-19": {
    "Free Sorting Confirmed Sorts": { m1:10.22, sd1:1.63, m2:11.67, sd2:2.34, r:0.49, n:28 },
    "Free Sorting Description Total Score": { m1:10.19, sd1:1.57, m2:11.85, sd2:2.48, r:0.67, n:28 },
    "Sort Recognition Total Description Score": { m1:10.22, sd1:2.95, m2:11.81, sd2:2.77, r:0.56, n:28 }
  },
  "D-KEFS Sorting Test · All Ages": {
    "Free Sorting Confirmed Sorts": { m1:10.24, sd1:2.77, m2:11.31, sd2:2.44, r:0.51, n:101, rInternalByAge:{8:0.59, 9:0.58, 10:0.8, 11:0.7, 12:0.62, 13:0.73, 14:0.82, 15:0.55, 16:0.72, 20:0.78, 30:0.82, 40:0.81, 50:0.86, 60:0.81, 70:0.81, 80:0.77}, rInternalAgeMax:89 },
    "Free Sorting Description Total Score": { m1:10.11, sd1:2.72, m2:11.27, sd2:2.51, r:0.5, n:101, rInternalByAge:{8:0.62, 9:0.64, 10:0.77, 11:0.73, 12:0.64, 13:0.7, 14:0.8, 15:0.55, 16:0.73, 20:0.77, 30:0.83, 40:0.8, 50:0.84, 60:0.8, 70:0.82, 80:0.77}, rInternalAgeMax:89 },
    "Sort Recognition Total Description Score": { m1:10.26, sd1:2.92, m2:11.16, sd2:3.13, r:0.6, n:101, rInternalByAge:{8:0.74, 9:0.71, 10:0.62, 11:0.72, 12:0.67, 13:0.62, 14:0.72, 15:0.72, 16:0.74, 20:0.75, 30:0.77, 40:0.8, 50:0.74, 60:0.81, 70:0.79, 80:0.7}, rInternalAgeMax:89 }
  },
  "D-KEFS Tower Test · Ages 20-49": {
    "Total Achievement Score": { m1:10.33, sd1:3.03, m2:11.1, sd2:3.04, r:0.41, n:30 }
  },
  "D-KEFS Tower Test · Ages 50-89": {
    "Total Achievement Score": { m1:9.79, sd1:3.45, m2:11.89, sd2:2.96, r:0.38, n:28 }
  },
  "D-KEFS Tower Test · Ages 8-19": {
    "Total Achievement Score": { m1:11, sd1:3.14, m2:12.08, sd2:2.8, r:0.51, n:25 }
  },
  "D-KEFS Tower Test · All Ages": {
    "Total Achievement Score": { m1:10.35, sd1:3.21, m2:11.66, sd2:2.94, r:0.44, n:83, rInternalByAge:{8:0.56, 9:0.71, 10:0.84, 11:0.61, 12:0.61, 13:0.55, 14:0.43, 15:0.6, 16:0.6, 20:0.62, 30:0.72, 40:0.72, 50:0.56, 60:0.72, 70:0.78, 80:0.61}, rInternalAgeMax:89 }
  },
  "D-KEFS Trail Making Test · Ages 20-49": {
    "Visual Scanning": { m1:10.14, sd1:2.87, m2:10.86, sd2:2.68, r:0.55, n:35 },
    "Number Sequencing": { m1:9.57, sd1:3.18, m2:10.91, sd2:2.22, r:0.54, n:35 },
    "Letter Sequencing": { m1:9.69, sd1:3.32, m2:10.63, sd2:2.29, r:0.48, n:35 },
    "Switching": { m1:9.63, sd1:2.8, m2:10.97, sd2:1.87, r:0.36, n:35 },
    "Motor Speed": { m1:9.77, sd1:3.61, m2:10.51, sd2:3.06, r:0.73, n:35 },
    "Combined Number + Letter": { m1:9.54, sd1:3.37, m2:10.83, sd2:2.26, r:0.64, n:35 }
  },
  "D-KEFS Trail Making Test · Ages 50-89": {
    "Visual Scanning": { m1:10.53, sd1:2.66, m2:10.9, sd2:2.7, r:0.63, n:38 },
    "Number Sequencing": { m1:9.7, sd1:2.92, m2:11.54, sd2:2.21, r:0.37, n:37 },
    "Letter Sequencing": { m1:9.9, sd1:3.25, m2:10.76, sd2:3.04, r:0.7, n:38 },
    "Switching": { m1:10.33, sd1:2.99, m2:10.61, sd2:3.25, r:0.55, n:36 },
    "Motor Speed": { m1:10.42, sd1:2.47, m2:10.45, sd2:3.17, r:0.74, n:38 },
    "Combined Number + Letter": { m1:9.66, sd1:3.23, m2:11.13, sd2:2.92, r:0.6, n:37 }
  },
  "D-KEFS Trail Making Test · Ages 8-19": {
    "Visual Scanning": { m1:9.29, sd1:2.84, m2:11.29, sd2:1.65, r:0.5, n:28 },
    "Number Sequencing": { m1:10.44, sd1:2.85, m2:11.44, sd2:1.99, r:0.77, n:27 },
    "Letter Sequencing": { m1:9.44, sd1:3.19, m2:11.15, sd2:2.69, r:0.57, n:27 },
    "Switching": { m1:9.36, sd1:2.93, m2:10.5, sd2:3.05, r:0.2, n:28 },
    "Motor Speed": { m1:10.32, sd1:2.71, m2:10.68, sd2:2.55, r:0.82, n:28 },
    "Combined Number + Letter": { m1:9.7, sd1:3.54, m2:11.22, sd2:3.02, r:0.78, n:26 }
  },
  "D-KEFS Trail Making Test · All Ages": {
    "Visual Scanning": { m1:10.05, sd1:2.8, m2:10.99, sd2:2.43, r:0.56, n:101 },
    "Number Sequencing": { m1:9.86, sd1:2.99, m2:11.29, sd2:2.15, r:0.59, n:99 },
    "Letter Sequencing": { m1:9.7, sd1:3.23, m2:10.82, sd2:2.68, r:0.59, n:100 },
    "Switching": { m1:9.81, sd1:2.91, m2:10.71, sd2:2.75, r:0.38, n:99 },
    "Motor Speed": { m1:10.17, sd1:2.96, m2:10.54, sd2:2.95, r:0.77, n:101 },
    "Combined Number + Letter": { m1:9.63, sd1:3.33, m2:11.05, sd2:2.71, r:0.66, n:98, rInternalByAge:{8:0.78, 9:0.72, 10:0.57, 11:0.59, 12:0.68, 13:0.69, 14:0.79, 15:0.72, 16:0.69, 20:0.78, 30:0.78, 40:0.74, 50:0.81, 60:0.8, 70:0.6, 80:0.77}, rInternalAgeMax:89 }
  },
  "D-KEFS Twenty Questions Test · Ages 20-49": {
    "Total Weighted Achievement": { m1:10.61, sd1:2.73, m2:10.33, sd2:2.91, r:0.19, n:35 },
    "Initial Abstraction Score": { m1:10.21, sd1:2, m2:9.74, sd2:2.25, r:0.24, n:35 }
  },
  "D-KEFS Twenty Questions Test · Ages 50-89": {
    "Total Weighted Achievement": { m1:9.68, sd1:3.43, m2:10.18, sd2:3.05, r:0.39, n:38 },
    "Initial Abstraction Score": { m1:9.36, sd1:2.63, m2:9.62, sd2:2.37, r:0.42, n:38 }
  },
  "D-KEFS Twenty Questions Test · Ages 8-19": {
    "Total Weighted Achievement": { m1:9.36, sd1:2.84, m2:10.16, sd2:2.87, r:0.06, n:28 },
    "Initial Abstraction Score": { m1:9.64, sd1:2.34, m2:9.61, sd2:2.63, r:0.62, n:28 }
  },
  "D-KEFS Twenty Questions Test · All Ages": {
    "Total Weighted Achievement": { m1:9.92, sd1:3.05, m2:10.23, sd2:2.92, r:0.24, n:101, rInternalByAge:{8:0.44, 9:0.48, 10:0.46, 11:0.41, 12:0.51, 13:0.39, 14:0.53, 15:0.37, 16:0.1, 20:0.5, 30:0.37, 40:0.33, 50:0.26, 60:0.47, 70:0.55, 80:0.48}, rInternalAgeMax:89 },
    /* Initial Abstraction: Table 2.15 verbatim. Its Table 2.17 SEM column does
       NOT follow from these — see the Twenty Questions note above normDB and
       check.js section 27. Do not "correct" .75 to .66 to make 2.17 reconcile. */
    "Initial Abstraction Score": { m1:9.73, sd1:2.35, m2:9.66, sd2:2.38, r:0.43, n:101, rInternalByAge:{8:0.85, 9:0.72, 10:0.76, 11:0.83, 12:0.82, 13:0.81, 14:0.8, 15:0.87, 16:0.74, 20:0.85, 30:0.77, 40:0.75, 50:0.86, 60:0.85, 70:0.87, 80:0.77}, rInternalAgeMax:89 }
  },
  "D-KEFS Verbal Fluency · Ages 20-49": {
    "Letter Fluency": { m1:9.39, sd1:3.11, m2:9.91, sd2:3.64, r:0.76, n:35 },
    "Category Fluency": { m1:10.06, sd1:3.3, m2:10.52, sd2:3.42, r:0.81, n:35 },
    "Category Switching": { m1:9.39, sd1:3.63, m2:8.24, sd2:3.66, r:0.49, n:35 },
    "Switching Accuracy": { m1:10.36, sd1:3, m2:9.3, sd2:3.88, r:0.24, n:35 }
  },
  "D-KEFS Verbal Fluency · Ages 50-89": {
    "Letter Fluency": { m1:9.71, sd1:3.56, m2:9.97, sd2:3.88, r:0.88, n:38 },
    "Category Fluency": { m1:10.11, sd1:3.68, m2:10.24, sd2:3.66, r:0.82, n:38 },
    "Category Switching": { m1:10.71, sd1:3.44, m2:10.79, sd2:3.86, r:0.51, n:38 },
    "Switching Accuracy": { m1:10.55, sd1:3.13, m2:11.26, sd2:3.55, r:0.39, n:38 }
  },
  "D-KEFS Verbal Fluency · Ages 8-19": {
    "Letter Fluency": { m1:9.75, sd1:2.61, m2:10.5, sd2:2.86, r:0.67, n:28 },
    "Category Fluency": { m1:9.18, sd1:2.52, m2:10.14, sd2:3, r:0.7, n:28 },
    "Category Switching": { m1:9.25, sd1:2.86, m2:10.54, sd2:2.32, r:0.65, n:28 },
    "Switching Accuracy": { m1:10.29, sd1:2.77, m2:11, sd2:3.15, r:0.53, n:28 }
  },
  "D-KEFS Verbal Fluency · All Ages": {
    "Letter Fluency": { m1:9.62, sd1:3.14, m2:10.1, sd2:3.51, r:0.8, n:101, rInternalByAge:{8:0.68, 9:0.71, 10:0.8, 11:0.76, 12:0.77, 13:0.81, 14:0.8, 15:0.78, 16:0.8, 20:0.85, 30:0.9, 40:0.77, 50:0.9, 60:0.85, 70:0.87, 80:0.86}, rInternalAgeMax:89 },
    "Category Fluency": { m1:9.83, sd1:3.25, m2:10.3, sd2:3.38, r:0.79, n:101, rInternalByAge:{8:0.6, 9:0.75, 10:0.71, 11:0.58, 12:0.72, 13:0.66, 14:0.68, 15:0.53, 16:0.6, 20:0.61, 30:0.76, 40:0.63, 50:0.62, 60:0.64, 70:0.65, 80:0.76}, rInternalAgeMax:89 },
    "Category Switching": { m1:9.86, sd1:3.39, m2:9.87, sd2:3.58, r:0.52, n:101, rInternalByAge:{8:0.37, 9:0.53, 10:0.56, 11:0.62, 12:0.62, 13:0.62, 14:0.54, 15:0.44, 16:0.48, 20:0.43, 30:0.68, 40:0.68, 50:0.45, 60:0.5, 70:0.54, 80:0.55}, rInternalAgeMax:89 },
    "Switching Accuracy": { m1:10.41, sd1:2.96, m2:10.54, sd2:3.63, r:0.36, n:101, rInternalByAge:{8:0.53, 9:0.73, 10:0.64, 11:0.76, 12:0.65, 13:0.54, 14:0.73, 15:0.63, 16:0.53, 20:0.59, 30:0.72, 40:0.62, 50:0.53, 60:0.51, 70:0.64, 80:0.71}, rInternalAgeMax:89 }
  },
  "D-KEFS Word Context Test · Ages 20-49": {
    "Total First Trial Consistently Correct": { m1:10.26, sd1:3.08, m2:11.97, sd2:2.94, r:0.73, n:35 }
  },
  "D-KEFS Word Context Test · Ages 50-89": {
    "Total First Trial Consistently Correct": { m1:9.42, sd1:2.66, m2:10.61, sd2:3.53, r:0.78, n:38 }
  },
  "D-KEFS Word Context Test · Ages 8-19": {
    "Total First Trial Consistently Correct": { m1:10.57, sd1:2.74, m2:12.25, sd2:3.1, r:0.58, n:28 }
  },
  "D-KEFS Word Context Test · All Ages": {
    "Total First Trial Consistently Correct": { m1:10.03, sd1:2.85, m2:11.54, sd2:3.27, r:0.7, n:101, rInternalByAge:{8:0.55, 9:0.52, 10:0.52, 11:0.59, 12:0.52, 13:0.47, 14:0.68, 15:0.71, 16:0.51, 20:0.68, 30:0.67, 40:0.53, 50:0.72, 60:0.68, 70:0.74, 80:0.68}, rInternalAgeMax:89 }
  },
  "D-KEFS Word Proverb Test · Ages 20-49": {
    "Total Achievement Score: Free Inquiry": { m1:10.13, sd1:2.62, m2:10.56, sd2:2.26, r:0.66, n:35 }
  },
  "D-KEFS Word Proverb Test · Ages 50-89": {
    "Total Achievement Score: Free Inquiry": { m1:9.46, sd1:3.4, m2:10.38, sd2:3.73, r:0.81, n:38 }
  },
  "D-KEFS Word Proverb Test · Ages 16-19": {
    "Total Achievement Score: Free Inquiry": { m1:9.8, sd1:3.29, m2:11.3, sd2:2.5, r:0.9, n:28 }
  },
  "D-KEFS Word Proverb Test · All Ages": {
    "Total Achievement Score": { m1:9.77, sd1:3.07, m2:10.57, sd2:3.04, r:0.76, n:101, rInternalByAge:{16:0.68, 20:0.71, 30:0.8, 40:0.76, 50:0.77, 60:0.81, 70:0.8, 80:0.78}, rInternalAgeMax:89 }
  },
  "D-KEFS Advanced Trail Making · All Ages": {
    "Number Sequencing Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.8, m2:10.4, sd2:2.8, r:0.53, n:224 },
    "Letter Sequencing Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.9, m2:10.4, sd2:3, r:0.52, n:224 },
    "Number–Letter Switching Mean Correct (With Errors) Connection Time": { m1:10, sd1:3.1, m2:11.1, sd2:3.1, r:0.64, n:224 },
    "Switching–Distraction Mean Correct (With Errors) Connection Time": { m1:9.9, sd1:3.1, m2:11.1, sd2:3.1, r:0.69, n:224 },
    "Switching–Working Memory Mean Correct (With Errors) Connection Time": { m1:10.1, sd1:2.9, m2:10.7, sd2:3.3, r:0.66, n:224 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.2, m2:10.7, sd2:2.8, r:0.48, n:224 },
    "Number–Letter Switching Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:3, m2:11.1, sd2:3, r:0.66, n:224 },
    "Switching–Distraction Mean Pure Response Speed (Correct or Error Connections)": { m1:9.9, sd1:3.1, m2:11, sd2:3.2, r:0.7, n:224 },
    "Switching–Working Memory Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:2.9, m2:10.6, sd2:3.3, r:0.7, n:224 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:102.3, sd1:14.3, m2:102.5, sd2:15.1, r:0.67, n:224 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:99.6, sd1:15.8, m2:105.2, sd2:16.4, r:0.79, n:224 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:99.8, sd1:15.4, m2:105, sd2:16.2, r:0.81, n:224 },
    "Multitasking Index": { m1:99.7, sd1:16.4, m2:105.5, sd2:16.4, r:0.73, n:224 }
  },
  "D-KEFS Advanced Trail Making · Ages 8-18": {
    "Number Sequencing Mean Correct (With Errors) Connection Time": { m1:10.1, sd1:3, m2:10.6, sd2:2.9, r:0.47, n:91 },
    "Letter Sequencing Mean Correct (With Errors) Connection Time": { m1:10.2, sd1:2.8, m2:10.6, sd2:2.7, r:0.47, n:91 },
    "Number–Letter Switching Mean Correct (With Errors) Connection Time": { m1:9.8, sd1:3.2, m2:11.5, sd2:3.3, r:0.62, n:91 },
    "Switching–Distraction Mean Correct (With Errors) Connection Time": { m1:9.8, sd1:3.4, m2:11.5, sd2:3.3, r:0.77, n:91 },
    "Switching–Working Memory Mean Correct (With Errors) Connection Time": { m1:9.6, sd1:3.1, m2:10.6, sd2:3.4, r:0.76, n:91 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.3, m2:10.8, sd2:2.9, r:0.4, n:91 },
    "Number–Letter Switching Mean Pure Response Speed (Correct or Error Connections)": { m1:9.7, sd1:3.1, m2:11.5, sd2:3.3, r:0.65, n:91 },
    "Switching–Distraction Mean Pure Response Speed (Correct or Error Connections)": { m1:9.9, sd1:3.2, m2:11.5, sd2:3.2, r:0.75, n:91 },
    "Switching–Working Memory Mean Pure Response Speed (Correct or Error Connections)": { m1:9.5, sd1:2.9, m2:10.6, sd2:3.3, r:0.74, n:91 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:100.5, sd1:14.1, m2:103.4, sd2:14.6, r:0.64, n:91 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:98.2, sd1:16.7, m2:106.7, sd2:17.5, r:0.83, n:91 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:98.3, sd1:15.7, m2:106.8, sd2:16.7, r:0.84, n:91 },
    "Multitasking Index": { m1:98.4, sd1:16.9, m2:106.8, sd2:17, r:0.77, n:91 }
  },
  "D-KEFS Advanced Trail Making · Ages 19-59": {
    "Number Sequencing Mean Correct (With Errors) Connection Time": { m1:11, sd1:2.5, m2:10.8, sd2:2.6, r:0.39, n:67 },
    "Letter Sequencing Mean Correct (With Errors) Connection Time": { m1:11, sd1:3, m2:10.7, sd2:3.1, r:0.46, n:67 },
    "Number–Letter Switching Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:2.9, m2:11.2, sd2:2.5, r:0.59, n:67 },
    "Switching–Distraction Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.8, m2:11.3, sd2:2.8, r:0.48, n:67 },
    "Switching–Working Memory Mean Correct (With Errors) Connection Time": { m1:10.5, sd1:2.5, m2:10.7, sd2:2.9, r:0.44, n:67 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10.1, sd1:2.8, m2:10.8, sd2:2.7, r:0.52, n:67 },
    "Number–Letter Switching Mean Pure Response Speed (Correct or Error Connections)": { m1:10.4, sd1:2.9, m2:11.2, sd2:2.6, r:0.59, n:67 },
    "Switching–Distraction Mean Pure Response Speed (Correct or Error Connections)": { m1:10.5, sd1:2.7, m2:11.2, sd2:2.8, r:0.55, n:67 },
    "Switching–Working Memory Mean Pure Response Speed (Correct or Error Connections)": { m1:10.5, sd1:2.4, m2:10.7, sd2:2.9, r:0.54, n:67 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:105.4, sd1:12.9, m2:104.2, sd2:15.1, r:0.55, n:67 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:102.3, sd1:13.4, m2:106, sd2:14, r:0.67, n:67 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:102.9, sd1:13.1, m2:106, sd2:13.7, r:0.69, n:67 },
    "Multitasking Index": { m1:101.9, sd1:14.1, m2:106.2, sd2:14.2, r:0.63, n:67 }
  },
  "D-KEFS Advanced Trail Making · Ages 60-90": {
    "Number Sequencing Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:2.7, m2:9.8, sd2:2.8, r:0.7, n:66 },
    "Letter Sequencing Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:3.2, m2:10, sd2:3.3, r:0.62, n:66 },
    "Number–Letter Switching Mean Correct (With Errors) Connection Time": { m1:10, sd1:3.3, m2:10.4, sd2:3.5, r:0.7, n:66 },
    "Switching–Distraction Mean Correct (With Errors) Connection Time": { m1:9.5, sd1:3.1, m2:10.4, sd2:2.9, r:0.76, n:66 },
    "Switching–Working Memory Mean Correct (With Errors) Connection Time": { m1:10.2, sd1:3.1, m2:10.7, sd2:3.5, r:0.73, n:66 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.4, m2:10.5, sd2:2.9, r:0.52, n:66 },
    "Number–Letter Switching Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:3, m2:10.5, sd2:3.1, r:0.72, n:66 },
    "Switching–Distraction Mean Pure Response Speed (Correct or Error Connections)": { m1:9.3, sd1:3.2, m2:10, sd2:3.4, r:0.77, n:66 },
    "Switching–Working Memory Mean Pure Response Speed (Correct or Error Connections)": { m1:10.1, sd1:3.4, m2:10.5, sd2:3.6, r:0.79, n:66 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:101.6, sd1:15.4, m2:99.5, sd2:15.4, r:0.78, n:66 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:99, sd1:16.8, m2:102.4, sd2:17.1, r:0.83, n:66 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:98.7, sd1:16.8, m2:101.5, sd2:17.6, r:0.87, n:66 },
    "Multitasking Index": { m1:99.1, sd1:17.9, m2:103, sd2:17.5, r:0.76, n:66 }
  },
  "D-KEFS Advanced Verbal Fluency · All Ages": {
    "Letter Fluency Total Correct": { m1:10.1, sd1:3.3, m2:10.9, sd2:3.3, r:0.81, n:224 },
    "Category Fluency Total Correct": { m1:10.3, sd1:3.1, m2:10.3, sd2:3.1, r:0.79, n:224 },
    "Switching Fluency Total Correct": { m1:10.1, sd1:3.1, m2:10.6, sd2:3.2, r:0.79, n:224 },
    "Switching Fluency Total Accurate Switches": { m1:10.1, sd1:3.2, m2:10.6, sd2:3.1, r:0.78, n:224 },
    "Total Set-Loss Errors": { m1:9.8, sd1:3, m2:10, sd2:2.9, r:0.19, n:224 },
    "Total Repetitions": { m1:9.8, sd1:2.9, m2:9.9, sd2:2.7, r:0.44, n:224 },
    "Switching Fluency Total Correct/Switching Accuracy Index": { m1:100.4, sd1:15.9, m2:102.9, sd2:16.2, r:0.8, n:224 },
    "Combined Letter and Category Fluency Total Correct Index": { m1:100.8, sd1:16.4, m2:103.4, sd2:16.4, r:0.84, n:224 },
    "Combined Letter, Category, and Switching Fluency Total Correct Index": { m1:100.5, sd1:16.3, m2:103.2, sd2:16.2, r:0.87, n:224 }
  },
  "D-KEFS Advanced Verbal Fluency · Ages 8-18": {
    "Letter Fluency Total Correct": { m1:9.6, sd1:3, m2:10.7, sd2:3.2, r:0.76, n:91 },
    "Category Fluency Total Correct": { m1:9.9, sd1:3.2, m2:9.9, sd2:2.8, r:0.8, n:91 },
    "Switching Fluency Total Correct": { m1:9.6, sd1:3.1, m2:10.3, sd2:3.1, r:0.66, n:91 },
    "Switching Fluency Total Accurate Switches": { m1:9.7, sd1:3.2, m2:10.5, sd2:3.1, r:0.68, n:91 },
    "Total Set-Loss Errors": { m1:9.9, sd1:2.9, m2:10.1, sd2:3, r:0.33, n:91 },
    "Total Repetitions": { m1:10.1, sd1:2.8, m2:10, sd2:2.6, r:0.41, n:91 },
    "Switching Fluency Total Correct/Switching Accuracy Index": { m1:98.2, sd1:15.8, m2:101.7, sd2:15.7, r:0.67, n:91 },
    "Combined Letter and Category Fluency Total Correct Index": { m1:98.6, sd1:15.9, m2:101.5, sd2:15.1, r:0.81, n:91 },
    "Combined Letter, Category, and Switching Fluency Total Correct Index": { m1:97.4, sd1:15.1, m2:100.8, sd2:14.7, r:0.78, n:91 }
  },
  "D-KEFS Advanced Verbal Fluency · Ages 19-59": {
    "Letter Fluency Total Correct": { m1:10.9, sd1:3.4, m2:11.7, sd2:3.3, r:0.78, n:67 },
    "Category Fluency Total Correct": { m1:10.9, sd1:2.9, m2:10.8, sd2:3.4, r:0.71, n:67 },
    "Switching Fluency Total Correct": { m1:11, sd1:3.1, m2:11.2, sd2:3.2, r:0.81, n:67 },
    "Switching Fluency Total Accurate Switches": { m1:10.9, sd1:2.9, m2:11.2, sd2:3.2, r:0.76, n:67 },
    "Total Set-Loss Errors": { m1:10, sd1:2.8, m2:10.4, sd2:2.8, r:0, n:67 },
    "Total Repetitions": { m1:9.5, sd1:2.7, m2:9.8, sd2:3, r:0.5, n:67 },
    "Switching Fluency Total Correct/Switching Accuracy Index": { m1:104.8, sd1:15.3, m2:106.1, sd2:16.5, r:0.81, n:67 },
    "Combined Letter and Category Fluency Total Correct Index": { m1:104.8, sd1:16.1, m2:106.9, sd2:17, r:0.78, n:67 },
    "Combined Letter, Category, and Switching Fluency Total Correct Index": { m1:105.1, sd1:16.4, m2:106.9, sd2:16.4, r:0.85, n:67 }
  },
  "D-KEFS Advanced Verbal Fluency · Ages 60-90": {
    "Letter Fluency Total Correct": { m1:10, sd1:3.5, m2:10.5, sd2:3.5, r:0.88, n:66 },
    "Category Fluency Total Correct": { m1:10.1, sd1:3, m2:10.5, sd2:3.1, r:0.85, n:66 },
    "Switching Fluency Total Correct": { m1:9.9, sd1:3.1, m2:10.4, sd2:3.4, r:0.86, n:66 },
    "Switching Fluency Total Accurate Switches": { m1:9.8, sd1:3.2, m2:10.2, sd2:3.1, r:0.87, n:66 },
    "Total Set-Loss Errors": { m1:9.6, sd1:3.4, m2:9.5, sd2:3, r:0.23, n:66 },
    "Total Repetitions": { m1:9.7, sd1:3.3, m2:9.8, sd2:2.6, r:0.41, n:66 },
    "Switching Fluency Total Correct/Switching Accuracy Index": { m1:99, sd1:16, m2:101.3, sd2:16.5, r:0.87, n:66 },
    "Combined Letter and Category Fluency Total Correct Index": { m1:99.9, sd1:17, m2:102.6, sd2:17.1, r:0.91, n:66 },
    "Combined Letter, Category, and Switching Fluency Total Correct Index": { m1:100, sd1:16.9, m2:102.7, sd2:17.5, r:0.94, n:66 }
  },
  "D-KEFS Advanced Colour-Word Interference · All Ages": {
    "Colour Identification Net Correct Responses": { m1:10.2, sd1:3, m2:10.5, sd2:3.2, r:0.75, n:224, rInternal:0.88, rInternalByAge:{8:0.8, 9:0.85, 10:0.77, 11:0.85, 12:0.79, 13:0.89, 14:0.84, 15:0.84, 16:0.85, 19:0.9, 30:0.9, 40:0.91, 50:0.94, 60:0.93, 70:0.92, 80:0.92}, rInternalAgeMax:90 },
    "Word Identification Net Correct Responses": { m1:10.2, sd1:2.8, m2:10.4, sd2:3, r:0.78, n:224, rInternal:0.87, rInternalByAge:{8:0.84, 9:0.89, 10:0.72, 11:0.88, 12:0.84, 13:0.84, 14:0.75, 15:0.84, 16:0.78, 19:0.88, 30:0.9, 40:0.87, 50:0.94, 60:0.9, 70:0.92, 80:0.96}, rInternalAgeMax:90 },
    "Inhibition Net Correct Responses": { m1:10, sd1:3.2, m2:10.7, sd2:3.2, r:0.81, n:224, rInternal:0.86, rInternalByAge:{8:0.69, 9:0.75, 10:0.78, 11:0.77, 12:0.82, 13:0.85, 14:0.84, 15:0.87, 16:0.83, 19:0.87, 30:0.82, 40:0.9, 50:0.92, 60:0.94, 70:0.94, 80:0.95}, rInternalAgeMax:90 },
    "Inhibition/Switching Net Correct Responses": { m1:9.8, sd1:3, m2:11, sd2:3, r:0.72, n:224, rInternal:0.92, rInternalByAge:{8:0.92, 9:0.92, 10:0.9, 11:0.91, 12:0.9, 13:0.93, 14:0.89, 15:0.91, 16:0.93, 19:0.94, 30:0.91, 40:0.93, 50:0.94, 60:0.93, 70:0.93, 80:0.91}, rInternalAgeMax:90 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:9.7, sd1:3.2, m2:10.4, sd2:3, r:0.6, n:224, rInternal:0.81, rInternalByAge:{8:0.86, 9:0.93, 10:0.76, 11:0.81, 12:0.85, 13:0.76, 14:0.81, 15:0.73, 16:0.73, 19:0.52, 30:0.79, 40:0.65, 50:0.85, 60:0.83, 70:0.9, 80:0.86}, rInternalAgeMax:90 },
    "Combined Colour and Word Identification Net Correct Responses Index": { m1:101.1, sd1:14.6, m2:102.4, sd2:15.7, r:0.81, n:224, rInternal:0.93, rInternalByAge:{8:0.9, 9:0.92, 10:0.85, 11:0.93, 12:0.9, 13:0.93, 14:0.88, 15:0.91, 16:0.9, 19:0.94, 30:0.95, 40:0.94, 50:0.97, 60:0.95, 70:0.96, 80:0.97}, rInternalAgeMax:90 },
    "Multitasking Index": { m1:99.6, sd1:15.8, m2:105, sd2:16.2, r:0.82, n:224, rInternal:0.92, rInternalByAge:{8:0.83, 9:0.88, 10:0.89, 11:0.87, 12:0.9, 13:0.93, 14:0.91, 15:0.92, 16:0.92, 19:0.94, 30:0.9, 40:0.94, 50:0.95, 60:0.95, 70:0.95, 80:0.95}, rInternalAgeMax:90 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:100.9, sd1:14.1, m2:104.8, sd2:15.6, r:0.84, n:224, rInternal:0.96, rInternalByAge:{8:0.95, 9:0.96, 10:0.89, 11:0.92, 12:0.94, 13:0.95, 14:0.94, 15:0.97, 16:0.97, 19:0.97, 30:0.97, 40:0.95, 50:0.98, 60:0.98, 70:0.99, 80:0.98}, rInternalAgeMax:90 }
  },
  "D-KEFS Advanced Colour-Word Interference · Ages 8-18": {
    "Colour Identification Net Correct Responses": { m1:10.3, sd1:3, m2:10.4, sd2:3.4, r:0.73, n:91 },
    "Word Identification Net Correct Responses": { m1:10.3, sd1:2.8, m2:10.4, sd2:3.1, r:0.75, n:91 },
    "Inhibition Net Correct Responses": { m1:10.1, sd1:3.4, m2:11, sd2:3.4, r:0.8, n:91 },
    "Inhibition/Switching Net Correct Responses": { m1:9.4, sd1:3.4, m2:10.8, sd2:3.2, r:0.79, n:91 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:9.3, sd1:3.4, m2:10.1, sd2:3.2, r:0.61, n:91 },
    "Combined Colour and Word Identification Net Correct Responses Index": { m1:101.6, sd1:14.4, m2:102.3, sd2:17, r:0.81, n:91 },
    "Multitasking Index": { m1:98.4, sd1:17.3, m2:105.3, sd2:17.4, r:0.84, n:91 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:100.2, sd1:15.3, m2:105.5, sd2:16.7, r:0.8, n:91 }
  },
  "D-KEFS Advanced Colour-Word Interference · Ages 19-59": {
    "Colour Identification Net Correct Responses": { m1:10.7, sd1:3.2, m2:11, sd2:3.3, r:0.8, n:67 },
    "Word Identification Net Correct Responses": { m1:10.4, sd1:3.2, m2:10.3, sd2:3.3, r:0.82, n:67 },
    "Inhibition Net Correct Responses": { m1:10.2, sd1:2.7, m2:10.7, sd2:3, r:0.77, n:67 },
    "Inhibition/Switching Net Correct Responses": { m1:10.5, sd1:2.7, m2:11.6, sd2:2.8, r:0.62, n:67 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:10, sd1:3.1, m2:10.3, sd2:3, r:0.7, n:67 },
    "Combined Colour and Word Identification Net Correct Responses Index": { m1:102.8, sd1:16.6, m2:103.3, sd2:17, r:0.84, n:67 },
    "Multitasking Index": { m1:102.2, sd1:13.4, m2:106.8, sd2:15.1, r:0.74, n:67 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:102.9, sd1:12.9, m2:106.6, sd2:14, r:0.79, n:67 }
  },
  "D-KEFS Advanced Colour-Word Interference · Ages 60-90": {
    "Colour Identification Net Correct Responses": { m1:9.7, sd1:2.6, m2:10.2, sd2:2.6, r:0.72, n:66 },
    "Word Identification Net Correct Responses": { m1:9.9, sd1:2.3, m2:10.5, sd2:2.3, r:0.77, n:66 },
    "Inhibition Net Correct Responses": { m1:9.7, sd1:3.2, m2:10.2, sd2:3, r:0.84, n:66 },
    "Inhibition/Switching Net Correct Responses": { m1:9.8, sd1:2.9, m2:10.7, sd2:3, r:0.72, n:66 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:10, sd1:3, m2:10.8, sd2:2.8, r:0.47, n:66 },
    "Combined Colour and Word Identification Net Correct Responses Index": { m1:98.6, sd1:12.4, m2:101.6, sd2:12.2, r:0.79, n:66 },
    "Multitasking Index": { m1:98.6, sd1:15.8, m2:102.8, sd2:15.5, r:0.85, n:66 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:99.8, sd1:13.8, m2:101.7, sd2:15.5, r:0.9, n:66 }
  },
  "D-KEFS Advanced Tower · All Ages": {
    "Global Performance Score": { m1:9.8, sd1:3.2, m2:11.2, sd2:3.3, r:0.56, n:224, rInternal:0.8, rInternalByAge:{8:0.83, 9:0.8, 10:0.83, 11:0.81, 12:0.78, 13:0.8, 14:0.74, 15:0.77, 16:0.76, 19:0.85, 30:0.77, 40:0.81, 50:0.81, 60:0.8, 70:0.78, 80:0.83}, rInternalAgeMax:90 },
    "Adjusted Mean Pure Response Time": { m1:10.1, sd1:3.2, m2:11.8, sd2:3.2, r:0.75, n:224, rInternal:0.9, rInternalByAge:{8:0.86, 9:0.88, 10:0.9, 11:0.89, 12:0.92, 13:0.87, 14:0.91, 15:0.9, 16:0.9, 19:0.88, 30:0.89, 40:0.9, 50:0.91, 60:0.88, 70:0.92, 80:0.92}, rInternalAgeMax:90 },
    "Adjusted Mean Unproductive Responses": { m1:9.9, sd1:3, m2:11.3, sd2:2.9, r:0.51, n:224, rInternal:0.43, rInternalByAge:{8:0.32, 9:0.45, 10:0.57, 11:0.36, 12:0.3, 13:0.36, 14:0.22, 15:0.41, 16:0.41, 19:0.42, 30:0.48, 40:0.47, 50:0.47, 60:0.54, 70:0.42, 80:0.55}, rInternalAgeMax:90 }
  },
  "D-KEFS Advanced Tower · Ages 8-18": {
    "Global Performance Score": { m1:9.3, sd1:3.1, m2:11.1, sd2:3.2, r:0.42, n:91 },
    "Adjusted Mean Pure Response Time": { m1:9.9, sd1:3.5, m2:12.8, sd2:2.9, r:0.65, n:91 },
    "Adjusted Mean Unproductive Responses": { m1:9.5, sd1:3, m2:11.2, sd2:2.9, r:0.39, n:91 }
  },
  "D-KEFS Advanced Tower · Ages 19-59": {
    "Global Performance Score": { m1:10.7, sd1:2.9, m2:12, sd2:3, r:0.54, n:67 },
    "Adjusted Mean Pure Response Time": { m1:10.8, sd1:3.2, m2:11.9, sd2:3.2, r:0.72, n:67 },
    "Adjusted Mean Unproductive Responses": { m1:10.4, sd1:3, m2:11.9, sd2:2.6, r:0.5, n:67 }
  },
  "D-KEFS Advanced Tower · Ages 60-90": {
    "Global Performance Score": { m1:9.6, sd1:3.5, m2:10.3, sd2:3.5, r:0.7, n:66 },
    "Adjusted Mean Pure Response Time": { m1:9.6, sd1:2.8, m2:10.4, sd2:3, r:0.84, n:66 },
    "Adjusted Mean Unproductive Responses": { m1:10, sd1:2.7, m2:10.8, sd2:3.2, r:0.61, n:66 }
  },
  "D-KEFS Advanced Social Sorting · All Ages": {
    "Global Performance Index": { m1:98.2, sd1:16.1, m2:105.8, sd2:19.5, r:0.59, n:224, rInternal:0.59, rInternalByAge:{8:0.56, 9:0.56, 10:0.56, 11:0.56, 12:0.56, 13:0.56, 14:0.56, 15:0.56, 16:0.56, 19:0.48, 30:0.48, 40:0.48, 50:0.48, 60:0.71, 70:0.71, 80:0.71}, rInternalAgeMax:90 },
    "Total Number of Perseverative Responses": { m1:10, sd1:3.2, m2:12.1, sd2:3.5, r:0.53, n:224, rInternal:0.95, rInternalByAge:{8:0.91, 9:0.95, 10:0.92, 11:0.94, 12:0.93, 13:0.93, 14:0.93, 15:0.86, 16:0.92, 19:0.95, 30:0.97, 40:0.96, 50:0.95, 60:0.97, 70:0.97, 80:0.98}, rInternalAgeMax:90 },
    "Percent Perseverative Responses": { m1:9.9, sd1:3.3, m2:12, sd2:3.5, r:0.46, n:224, rInternal:0.92, rInternalByAge:{8:0.88, 9:0.94, 10:0.89, 11:0.91, 12:0.9, 13:0.89, 14:0.88, 15:0.73, 16:0.85, 19:0.91, 30:0.95, 40:0.94, 50:0.93, 60:0.96, 70:0.96, 80:0.98}, rInternalAgeMax:90 },
    "Total Number of Perseverative Errors": { m1:9.8, sd1:3.2, m2:12.1, sd2:3.6, r:0.52, n:224, rInternal:0.92, rInternalByAge:{8:0.88, 9:0.94, 10:0.9, 11:0.92, 12:0.91, 13:0.91, 14:0.91, 15:0.84, 16:0.89, 19:0.93, 30:0.95, 40:0.94, 50:0.93, 60:0.95, 70:0.94, 80:0.97}, rInternalAgeMax:90 },
    "Percent Perseverative Errors": { m1:9.9, sd1:3.3, m2:12.2, sd2:3.7, r:0.47, n:224, rInternal:0.89, rInternalByAge:{8:0.84, 9:0.92, 10:0.85, 11:0.88, 12:0.87, 13:0.86, 14:0.85, 15:0.69, 16:0.81, 19:0.87, 30:0.93, 40:0.92, 50:0.91, 60:0.94, 70:0.94, 80:0.96}, rInternalAgeMax:90 },
    "Total Number of Errors": { m1:9.7, sd1:3.2, m2:11.6, sd2:3.8, r:0.6, n:224, rInternal:0.95, rInternalByAge:{8:0.94, 9:0.94, 10:0.94, 11:0.94, 12:0.94, 13:0.95, 14:0.95, 15:0.92, 16:0.94, 19:0.94, 30:0.97, 40:0.96, 50:0.96, 60:0.96, 70:0.91, 80:0.94}, rInternalAgeMax:90 },
    "Percent Correct Responses": { m1:9.8, sd1:3.2, m2:11.8, sd2:3.7, r:0.59, n:224, rInternal:0.91, rInternalByAge:{8:0.92, 9:0.92, 10:0.9, 11:0.9, 12:0.9, 13:0.91, 14:0.89, 15:0.82, 16:0.89, 19:0.88, 30:0.94, 40:0.94, 50:0.93, 60:0.94, 70:0.9, 80:0.93}, rInternalAgeMax:90 },
    "Total Number of Nonperseverative Errors": { m1:9.7, sd1:3.5, m2:10.7, sd2:3.6, r:0.54, n:224, rInternal:0.88, rInternalByAge:{8:0.91, 9:0.86, 10:0.9, 11:0.86, 12:0.86, 13:0.91, 14:0.85, 15:0.84, 16:0.92, 19:0.84, 30:0.89, 40:0.9, 50:0.89, 60:0.88, 70:0.87, 80:0.88}, rInternalAgeMax:90 },
    "Percent Nonperseverative Errors": { m1:9.7, sd1:3.4, m2:10.6, sd2:3.5, r:0.52, n:224, rInternal:0.83, rInternalByAge:{8:0.89, 9:0.82, 10:0.87, 11:0.79, 12:0.78, 13:0.87, 14:0.75, 15:0.71, 16:0.87, 19:0.72, 30:0.82, 40:0.86, 50:0.86, 60:0.85, 70:0.87, 80:0.87}, rInternalAgeMax:90 },
    "Total Number of Conceptual Level Responses": { m1:9.7, sd1:3.1, m2:10.4, sd2:3.3, r:0.34, n:224, rInternal:0.9, rInternalByAge:{8:0.91, 9:0.92, 10:0.87, 11:0.88, 12:0.87, 13:0.87, 14:0.89, 15:0.79, 16:0.87, 19:0.87, 30:0.92, 40:0.92, 50:0.92, 60:0.94, 70:0.94, 80:0.95}, rInternalAgeMax:90 },
    "Percent Conceptual Level Responses": { m1:9.7, sd1:3.2, m2:11.6, sd2:3.8, r:0.6, n:224, rInternal:0.95, rInternalByAge:{8:0.95, 9:0.95, 10:0.94, 11:0.93, 12:0.94, 13:0.94, 14:0.94, 15:0.9, 16:0.93, 19:0.92, 30:0.96, 40:0.96, 50:0.96, 60:0.97, 70:0.95, 80:0.96}, rInternalAgeMax:90 }
  },
  "D-KEFS Advanced Social Sorting · Ages 8-18": {
    "Global Performance Index": { m1:96.5, sd1:14.8, m2:106.1, sd2:19.1, r:0.56, n:91 },
    "Total Number of Perseverative Responses": { m1:9.6, sd1:3.1, m2:12.2, sd2:3.3, r:0.52, n:91 },
    "Percent Perseverative Responses": { m1:9.8, sd1:3, m2:12.3, sd2:3.1, r:0.48, n:91 },
    "Total Number of Perseverative Errors": { m1:9.5, sd1:3.1, m2:12.2, sd2:3.4, r:0.54, n:91 },
    "Percent Perseverative Errors": { m1:9.9, sd1:3, m2:12.5, sd2:3.4, r:0.48, n:91 },
    "Total Number of Errors": { m1:9.9, sd1:3, m2:12.2, sd2:3.7, r:0.6, n:91 },
    "Percent Correct Responses": { m1:9.5, sd1:3.2, m2:12, sd2:3.8, r:0.57, n:91 },
    "Total Number of Nonperseverative Errors": { m1:9.7, sd1:3.1, m2:11.1, sd2:3.5, r:0.52, n:91 },
    "Percent Nonperseverative Errors": { m1:9.7, sd1:3.2, m2:10.9, sd2:3.6, r:0.51, n:91 },
    "Total Number of Conceptual Level Responses": { m1:9.5, sd1:3.5, m2:10.3, sd2:3.5, r:0.23, n:91 },
    "Percent Conceptual Level Responses": { m1:9.3, sd1:3.3, m2:11.7, sd2:3.7, r:0.6, n:91 }
  },
  "D-KEFS Advanced Social Sorting · Ages 19-59": {
    "Global Performance Index": { m1:100.2, sd1:16.7, m2:107.8, sd2:19.8, r:0.48, n:67 },
    "Total Number of Perseverative Responses": { m1:10, sd1:3.2, m2:12.5, sd2:3.4, r:0.41, n:67 },
    "Percent Perseverative Responses": { m1:9.9, sd1:3.1, m2:12.3, sd2:3.3, r:0.29, n:67 },
    "Total Number of Perseverative Errors": { m1:9.9, sd1:3.4, m2:12.6, sd2:3.6, r:0.37, n:67 },
    "Percent Perseverative Errors": { m1:9.9, sd1:3.3, m2:12.7, sd2:3.7, r:0.33, n:67 },
    "Total Number of Errors": { m1:10, sd1:3.6, m2:12.1, sd2:3.9, r:0.44, n:67 },
    "Percent Correct Responses": { m1:10, sd1:3.3, m2:12.2, sd2:3.7, r:0.41, n:67 },
    "Total Number of Nonperseverative Errors": { m1:10.2, sd1:3.6, m2:11.2, sd2:3.4, r:0.46, n:67 },
    "Percent Nonperseverative Errors": { m1:10.1, sd1:3.4, m2:10.9, sd2:3.1, r:0.45, n:67 },
    "Total Number of Conceptual Level Responses": { m1:9.7, sd1:2.5, m2:10.2, sd2:2.9, r:0.13, n:67 },
    "Percent Conceptual Level Responses": { m1:10, sd1:3.2, m2:11.9, sd2:3.7, r:0.43, n:67 }
  },
  "D-KEFS Advanced Social Sorting · Ages 60-90": {
    "Global Performance Index": { m1:98.3, sd1:17.1, m2:103.3, sd2:19.9, r:0.71, n:66 },
    "Total Number of Perseverative Responses": { m1:10.5, sd1:3.3, m2:11.6, sd2:3.6, r:0.63, n:66 },
    "Percent Perseverative Responses": { m1:10, sd1:3.8, m2:11.3, sd2:4.1, r:0.58, n:66 },
    "Total Number of Perseverative Errors": { m1:10.1, sd1:3.3, m2:11.3, sd2:3.8, r:0.63, n:66 },
    "Percent Perseverative Errors": { m1:10, sd1:3.7, m2:11.2, sd2:4.1, r:0.57, n:66 },
    "Total Number of Errors": { m1:9.2, sd1:2.9, m2:10.3, sd2:3.7, r:0.72, n:66 },
    "Percent Correct Responses": { m1:10, sd1:3.1, m2:11, sd2:3.6, r:0.74, n:66 },
    "Total Number of Nonperseverative Errors": { m1:9.2, sd1:3.8, m2:9.6, sd2:3.6, r:0.64, n:66 },
    "Percent Nonperseverative Errors": { m1:9.4, sd1:3.7, m2:9.8, sd2:3.5, r:0.58, n:66 },
    "Total Number of Conceptual Level Responses": { m1:9.8, sd1:3.2, m2:10.5, sd2:3.5, r:0.61, n:66 },
    "Percent Conceptual Level Responses": { m1:9.9, sd1:3.3, m2:10.9, sd2:3.9, r:0.73, n:66 }
  },
  "D-KEFS Advanced Risk-Reward Decision · All Ages": {
    "Total Net Earnings": { m1:9.9, sd1:3.3, m2:12.6, sd2:3.7, r:0.56, n:224, rInternal:0.82, rInternalByAge:{19:0.74, 30:0.87, 40:0.82, 50:0.79, 60:0.84, 70:0.83, 80:0.85}, rInternalAgeMax:90 },
    "Net Earnings Races 1-20": { m1:9.8, sd1:3.3, m2:12.9, sd2:3.4, r:0.39, n:224, rInternal:0.45, rInternalByAge:{19:0.46, 30:0.67, 40:0.55, 50:0.42, 60:0.37, 70:0.33, 80:0.28}, rInternalAgeMax:90 },
    "Net Earnings Races 21-40": { m1:9.9, sd1:3.3, m2:10.8, sd2:2.4, r:0.55, n:224, rInternal:0.72, rInternalByAge:{19:0.64, 30:0.78, 40:0.64, 50:0.7, 60:0.75, 70:0.74, 80:0.76}, rInternalAgeMax:90 },
    "Net Earnings Races 41-60": { m1:10, sd1:3.2, m2:11.2, sd2:2.9, r:0.53, n:224, rInternal:0.62, rInternalByAge:{19:0.46, 30:0.65, 40:0.6, 50:0.57, 60:0.64, 70:0.63, 80:0.72}, rInternalAgeMax:90 }
  },
  "D-KEFS Advanced Risk-Reward Decision · Ages 19-59": {
    "Total Net Earnings": { m1:9.8, sd1:3.1, m2:13.3, sd2:3.4, r:0.49, n:67 },
    "Net Earnings Races 1-20": { m1:9.4, sd1:3.1, m2:13.5, sd2:2.9, r:0.33, n:67 },
    "Net Earnings Races 21-40": { m1:9.7, sd1:3.4, m2:10.8, sd2:2, r:0.48, n:67 },
    "Net Earnings Races 41-60": { m1:10.1, sd1:2.8, m2:11.4, sd2:2.7, r:0.43, n:67 }
  },
  "D-KEFS Advanced Risk-Reward Decision · Ages 60-90": {
    "Total Net Earnings": { m1:9.9, sd1:3.5, m2:12, sd2:3.9, r:0.63, n:66 },
    "Net Earnings Races 1-20": { m1:10.2, sd1:3.4, m2:12.4, sd2:3.8, r:0.45, n:66 },
    "Net Earnings Races 21-40": { m1:10.1, sd1:3.3, m2:10.9, sd2:2.8, r:0.61, n:66 },
    "Net Earnings Races 41-60": { m1:9.8, sd1:3.5, m2:10.9, sd2:3, r:0.62, n:66 }
  },
  "RBANS Indices · Ages 12-19": {
    "Immediate Memory": { m1:99.5, sd1:14.5, m2:119.8, sd2:16.5, r:0.73, rCorrected:0.75, n:55 },
    "Visuospatial/Constructional": { m1:99.7, sd1:13.4, m2:98.7, sd2:13.4, r:0.53, rCorrected:0.63, n:55 },
    "Attention": { m1:100.5, sd1:14.6, m2:103.6, sd2:17.3, r:0.69, rCorrected:0.71, n:55 },
    "Language": { m1:100.1, sd1:15.9, m2:104, sd2:15.1, r:0.79, rCorrected:0.76, n:55 },
    "Delayed Memory": { m1:100.6, sd1:12.2, m2:110.1, sd2:16.3, r:0.7, rCorrected:0.8, n:55 },
    "Total Scale": { m1:100.8, sd1:13.3, m2:110.4, sd2:14.8, r:0.81, rCorrected:0.85, n:55 }
  },
  "RBANS Indices · Ages 20-89": {
    "Immediate Memory": { m1:109.3, sd1:12.3, m2:110.6, sd2:14.4, r:0.62, rCorrected:0.75, n:40 },
    "Visuospatial/Constructional": { m1:97.8, sd1:14.2, m2:109.3, sd2:14.5, r:0.65, rCorrected:0.68, n:40 },
    "Attention": { m1:103.8, sd1:15.7, m2:105.7, sd2:16.6, r:0.77, rCorrected:0.75, n:40 },
    "Language": { m1:107.8, sd1:13.4, m2:105.4, sd2:14.8, r:0.64, rCorrected:0.71, n:40 },
    "Delayed Memory": { m1:108, sd1:13.8, m2:110.1, sd2:12.2, r:0.77, rCorrected:0.8, n:40 },
    "Total Scale": { m1:106.7, sd1:13.9, m2:110.6, sd2:13.2, r:0.81, rCorrected:0.84, n:40 }
  },
  "RBANS Subtests · Ages 12-19": {
    "List Learning": { m1:10.2, sd1:3.1, m2:13.4, sd2:3.2, r:0.68, rCorrected:0.66, n:55 },
    "Story Memory": { m1:9.8, sd1:2.7, m2:13, sd2:2.6, r:0.65, rCorrected:0.72, n:55 },
    "Figure Copy": { m1:10, sd1:2.7, m2:10, sd2:2.5, r:0.46, rCorrected:0.57, n:55 },
    "Line Orientation": { m1:16.7, sd1:3, m2:16.9, sd2:2.9, r:0.72, n:55, metric:'raw' },
    "Picture Naming": { m1:9.1, sd1:1, m2:9.2, sd2:0.9, r:0.73, n:55, metric:'raw' },
    "Semantic Fluency": { m1:10, sd1:3.1, m2:10.8, sd2:3, r:0.67, rCorrected:0.65, n:55 },
    "Digit Span": { m1:9.8, sd1:2.7, m2:10, sd2:3.2, r:0.59, rCorrected:0.67, n:55 },
    "Coding": { m1:10.2, sd1:2.7, m2:11.2, sd2:3.1, r:0.75, rCorrected:0.79, n:55 },
    "List Recall": { m1:6.8, sd1:1.8, m2:8.2, sd2:1.8, r:0.66, n:55, metric:'raw' },
    "List Recognition": { m1:19.9, sd1:0.5, m2:19.9, sd2:0.4, r:0.7, n:55, metric:'raw' },
    "Story Recall": { m1:10, sd1:3, m2:11.7, sd2:3.1, r:0.48, rCorrected:0.49, n:55 },
    "Figure Recall": { m1:10.2, sd1:2.5, m2:11.2, sd2:3.1, r:0.58, rCorrected:0.71, n:55 }
  },
  "RBANS Subtests · Ages 20-89": {
    "List Learning": { m1:11.5, sd1:2.9, m2:11.2, sd2:3.3, r:0.49, rCorrected:0.52, n:40 },
    "Story Memory": { m1:11.6, sd1:1.8, m2:12.5, sd2:2.4, r:0.45, rCorrected:0.8, n:40 },
    "Figure Copy": { m1:9.6, sd1:2.8, m2:11.9, sd2:2.6, r:0.47, rCorrected:0.54, n:40 },
    "Line Orientation": { m1:16, sd1:3.4, m2:16.4, sd2:3.7, r:0.49, n:40, metric:'raw' },
    "Picture Naming": { m1:9.8, sd1:0.4, m2:9.7, sd2:0.5, r:0.5, n:40, metric:'raw' },
    "Semantic Fluency": { m1:11.1, sd1:2.9, m2:11.2, sd2:3.3, r:0.49, rCorrected:0.52, n:40 },
    "Digit Span": { m1:10.4, sd1:3.5, m2:10.1, sd2:3.7, r:0.73, rCorrected:0.63, n:40 },
    "Coding": { m1:10.8, sd1:2.5, m2:11.7, sd2:2.8, r:0.76, rCorrected:0.83, n:40 },
    "List Recall": { m1:6.2, sd1:2.4, m2:5.8, sd2:2.7, r:0.6, n:40, metric:'raw' },
    "List Recognition": { m1:19.6, sd1:0.8, m2:19.8, sd2:0.5, r:0.27, n:40, metric:'raw' },
    "Story Recall": { m1:11.6, sd1:2.3, m2:11.6, sd2:2.3, r:0.52, rCorrected:0.72, n:40 },
    "Figure Recall": { m1:10.4, sd1:3, m2:11.5, sd2:3, r:0.55, rCorrected:0.55, n:40 }
  },
  /* ---------------------------------------------------------------
     RBANS · All Ages — the Score Tables groups, and the only RBANS
     groups carrying a reliability the publisher tabulates by age.

     WHY THESE GROUPS EXIST. Every other RBANS group holds a RETEST
     study, banded the way that study sampled (12-19, 20-89). Score
     Tables shows one entry per instrument and picks the "· All Ages"
     group; with none present it fell through to whichever group was
     listed FIRST, so every patient — an 80-year-old included — was
     scored on 55 adolescents. Immediate Memory printed 85-115 where
     Table 3.6 gives 90-110 at that age.

     SOURCE: RBANS Update Manual (Randolph 2012), Tables 3.6 and 3.7,
     p. 42. 14 measures x 9 NORMATIVE age bands. It clears the bar this
     project sets: the manual publishes the coefficients AND derives its
     own printed SEMs from them — all 126 cells of Table 3.7 equal
     SD sqrt(1 - rxx) on the SDs its note states ("SEMs expressed in
     scaled score (subtest) and standard score (index) units", so 3 and
     15). check.js section 29 pins every one.

     TWO BASES IN ONE TABLE, WHICH IS WHY THERE ARE TWO FIELDS.
     Table 3.6 footnote a marks five subtests "Reliability estimates
     based on test-retest": Figure Copy, Semantic Fluency, Coding,
     Story Recall, Figure Recall. Those take rStabilityByAge. The other
     nine are internal consistency and take rInternalByAge. Storing the
     five as rInternal would label a stability coefficient as internal
     consistency, which the APA note and Methods & References then
     assert on screen — the same mislabelling WAIS-IV Table 4.1 avoids
     for its three speeded subtests.

     The five are corroborated: their ADULT values equal this database's
     stored rCorrected 5/5, against 1/5 for the raw r — the same
     transcription proof the WAIS-IV speeded three give. The adolescent
     columns match only 2/5; those bands are a different retest sample,
     and Table 3.6 is taken as printed rather than reconciled to
     Table 3.8.

     THE FOUR RAW SUBTESTS APPEAR NOWHERE IN TABLE 3.6, which is an
     independent confirmation of the raw/scaled split recorded above:
     the manual publishes reliability for its eight SCALED subtests
     only. They therefore carry no coefficient and print no interval.
     That loss is deliberate — the interval they used to show came from
     55 adolescents whatever the patient's age, and List Recognition
     alone runs r .70 there against .27 in adults, so the honest adult
     interval is nearly twice the width that was printed. m1/sd1 are the
     adult retest descriptives, carried over only so the row stays
     selectable and declares its metric; nothing on screen derives from
     them, a raw row having no percentile, no classification and now no
     interval. Pinned in section 29.

     singleAdministration:true throughout: these groups hold no second
     testing, so isSingleAdministrationFamily keeps them out of Change
     Analysis and the SD Index, which go on using the retest groups.
     --------------------------------------------------------------- */
  "RBANS Indices · All Ages": {
    "Immediate Memory": { m1:100, sd1:15, singleAdministration:true, rInternal:0.88, rInternalAgeMax:89, rInternalByAge:{12:0.93, 14:0.81, 16:0.86, 20:0.84, 40:0.88, 50:0.89, 60:0.85, 70:0.89, 80:0.9} },
    "Visuospatial/Constructional": { m1:100, sd1:15, singleAdministration:true, rInternal:0.75, rInternalAgeMax:89, rInternalByAge:{12:0.64, 14:0.67, 16:0.53, 20:0.77, 40:0.81, 50:0.82, 60:0.84, 70:0.81, 80:0.78} },
    "Attention": { m1:100, sd1:15, singleAdministration:true, rInternal:0.84, rInternalAgeMax:89, rInternalByAge:{12:0.81, 14:0.82, 16:0.81, 20:0.84, 40:0.84, 50:0.85, 60:0.83, 70:0.88, 80:0.85} },
    "Language": { m1:100, sd1:15, singleAdministration:true, rInternal:0.8, rInternalAgeMax:89, rInternalByAge:{12:0.79, 14:0.8, 16:0.74, 20:0.75, 40:0.76, 50:0.87, 60:0.85, 70:0.81, 80:0.83} },
    "Delayed Memory": { m1:100, sd1:15, singleAdministration:true, rInternal:0.84, rInternalAgeMax:89, rInternalByAge:{12:0.85, 14:0.85, 16:0.84, 20:0.84, 40:0.83, 50:0.84, 60:0.85, 70:0.83, 80:0.81} },
    "Total Scale": { m1:100, sd1:15, singleAdministration:true, rInternal:0.93, rInternalAgeMax:89, rInternalByAge:{12:0.92, 14:0.91, 16:0.9, 20:0.92, 40:0.94, 50:0.95, 60:0.93, 70:0.93, 80:0.94} }
  },
  "RBANS Subtests · All Ages": {
    "List Learning": { m1:10, sd1:3, singleAdministration:true, rInternal:0.85, rInternalAgeMax:89, rInternalByAge:{12:0.91, 14:0.88, 16:0.8, 20:0.82, 40:0.88, 50:0.85, 60:0.8, 70:0.86, 80:0.84} },
    "Story Memory": { m1:10, sd1:3, singleAdministration:true, rInternal:0.78, rInternalAgeMax:89, rInternalByAge:{12:0.87, 14:0.55, 16:0.79, 20:0.71, 40:0.73, 50:0.82, 60:0.79, 70:0.8, 80:0.84} },
    "Figure Copy": { m1:10, sd1:3, singleAdministration:true, rStability:0.5, rStabilityAgeMax:89, rStabilityByAge:{12:0.42, 14:0.42, 16:0.42, 20:0.54, 40:0.54, 50:0.54, 60:0.54, 70:0.54, 80:0.54} },
    "Line Orientation": { m1:16, sd1:3.4, metric:'raw', singleAdministration:true },
    "Picture Naming": { m1:9.8, sd1:0.4, metric:'raw', singleAdministration:true },
    "Semantic Fluency": { m1:10, sd1:3, singleAdministration:true, rStability:0.57, rStabilityAgeMax:89, rStabilityByAge:{12:0.65, 14:0.65, 16:0.65, 20:0.52, 40:0.52, 50:0.52, 60:0.52, 70:0.52, 80:0.52} },
    "Digit Span": { m1:10, sd1:3, singleAdministration:true, rInternal:0.83, rInternalAgeMax:89, rInternalByAge:{12:0.71, 14:0.86, 16:0.85, 20:0.84, 40:0.83, 50:0.85, 60:0.76, 70:0.86, 80:0.83} },
    "Coding": { m1:10, sd1:3, singleAdministration:true, rStability:0.81, rStabilityAgeMax:89, rStabilityByAge:{12:0.76, 14:0.76, 16:0.76, 20:0.83, 40:0.83, 50:0.83, 60:0.83, 70:0.83, 80:0.83} },
    "List Recall": { m1:6.2, sd1:2.4, metric:'raw', singleAdministration:true },
    "List Recognition": { m1:19.6, sd1:0.8, metric:'raw', singleAdministration:true },
    "Story Recall": { m1:10, sd1:3, singleAdministration:true, rStability:0.54, rStabilityAgeMax:89, rStabilityByAge:{12:0.45, 14:0.45, 16:0.45, 20:0.72, 40:0.72, 50:0.72, 60:0.72, 70:0.72, 80:0.72} },
    "Figure Recall": { m1:10, sd1:3, singleAdministration:true, rStability:0.59, rStabilityAgeMax:89, rStabilityByAge:{12:0.71, 14:0.71, 16:0.71, 20:0.55, 40:0.55, 50:0.55, 60:0.55, 70:0.55, 80:0.55} }
  },
  // RBANS alternate-form (index-level only). m1/sd1 = Form A, m2/sd2 = alt form.
  // Change-Analysis only (filtered out of Score Tables). Randolph 2012:
  // A→B Table 3.10 (ages 20-89), A→C Table 3.11 (all ages), A→D Table 3.12 (all ages).
  "RBANS Indices (Form B) · Ages 20-89": {
    "Immediate Memory": { m1:105.5, sd1:12.9, m2:105.7, sd2:13.2, r:0.56, rCorrected:0.68, n:100 },
    "Visuospatial/Constructional": { m1:99.5, sd1:13.2, m2:100.5, sd2:12.9, r:0.54, rCorrected:0.65, n:100 },
    "Attention": { m1:106.5, sd1:14.4, m2:108.7, sd2:13.4, r:0.78, rCorrected:0.8, n:100 },
    "Language": { m1:105.6, sd1:14.2, m2:101.9, sd2:14.9, r:0.39, rCorrected:0.46, n:100 },
    "Delayed Memory": { m1:105.2, sd1:14, m2:100.7, sd2:15.8, r:0.59, rCorrected:0.64, n:100 },
    "Total Scale": { m1:106.2, sd1:13.8, m2:104.8, sd2:13, r:0.79, rCorrected:0.82, n:100 }
  },
  "RBANS Indices (Form C) · All Ages": {
    "Immediate Memory": { m1:90.7, sd1:17.1, m2:87.3, sd2:14.8, r:0.75, rCorrected:0.61, n:135 },
    "Visuospatial/Constructional": { m1:88.6, sd1:18.6, m2:92.9, sd2:19, r:0.77, rCorrected:0.57, n:135 },
    "Attention": { m1:94.4, sd1:15.7, m2:87.2, sd2:18.2, r:0.59, rCorrected:0.53, n:135 },
    "Language": { m1:96.5, sd1:16.9, m2:95.6, sd2:18.9, r:0.8, rCorrected:0.7, n:135 },
    "Delayed Memory": { m1:88.3, sd1:18.7, m2:91.6, sd2:16.6, r:0.58, rCorrected:0.24, n:135 },
    "Total Scale": { m1:88.6, sd1:16.5, m2:87.3, sd2:17.5, r:0.84, rCorrected:0.75, n:135 }
  },
  "RBANS Indices (Form D) · All Ages": {
    "Immediate Memory": { m1:91.9, sd1:18.4, m2:96.4, sd2:16.1, r:0.75, rCorrected:0.64, n:146 },
    "Visuospatial/Constructional": { m1:89.6, sd1:19.3, m2:90.9, sd2:20.8, r:0.75, rCorrected:0.61, n:146 },
    "Attention": { m1:94.3, sd1:16.7, m2:97.9, sd2:17.1, r:0.68, rCorrected:0.63, n:146 },
    "Language": { m1:94.9, sd1:16.6, m2:96.1, sd2:19.3, r:0.82, rCorrected:0.78, n:146 },
    "Delayed Memory": { m1:91.5, sd1:18.6, m2:93.9, sd2:23.2, r:0.71, rCorrected:0.58, n:146 },
    "Total Scale": { m1:89.2, sd1:19.3, m2:92.9, sd2:21, r:0.87, rCorrected:0.8, n:146 }
  },
  "WAIS-IV Core Subtests · Ages 16-29": {
    "Block Design": { m1:10.1, sd1:3, m2:11.3, sd2:2.9, r:0.81, rCorrected:0.81 },
    "Similarities": { m1:10.5, sd1:3.2, m2:10.9, sd2:3, r:0.82, rCorrected:0.8 },
    "Digit Span": { m1:10.1, sd1:2.8, m2:10.7, sd2:2.9, r:0.71, rCorrected:0.75 },
    "Matrix Reasoning": { m1:10.5, sd1:3.2, m2:10.7, sd2:3, r:0.7, rCorrected:0.66 },
    "Vocabulary": { m1:10.3, sd1:3.1, m2:10.5, sd2:3.2, r:0.91, rCorrected:0.9 },
    "Arithmetic": { m1:10, sd1:2.9, m2:10.5, sd2:2.9, r:0.84, rCorrected:0.85 },
    "Symbol Search": { m1:10.4, sd1:3.3, m2:11.6, sd2:3.4, r:0.84, rCorrected:0.81 },
    "Visual Puzzles": { m1:9.9, sd1:2.6, m2:11, sd2:3, r:0.74, rCorrected:0.8 },
    "Information": { m1:10, sd1:2.9, m2:10.7, sd2:3.1, r:0.9, rCorrected:0.91 },
    "Coding": { m1:9.8, sd1:2.8, m2:10.3, sd2:2.8, r:0.83, rCorrected:0.85 }
  },
  "WAIS-IV Core Subtests · Ages 30-54": {
    "Block Design": { m1:10.4, sd1:2.9, m2:11.4, sd2:3.1, r:0.8, rCorrected:0.81 },
    "Similarities": { m1:9.1, sd1:2.5, m2:9.7, sd2:2.6, r:0.85, rCorrected:0.9 },
    "Digit Span": { m1:10, sd1:3.2, m2:10.5, sd2:3.1, r:0.81, rCorrected:0.78 },
    "Matrix Reasoning": { m1:10.2, sd1:3.5, m2:10.3, sd2:3.4, r:0.85, rCorrected:0.8 },
    "Vocabulary": { m1:9.4, sd1:3.1, m2:9.5, sd2:3, r:0.9, rCorrected:0.89 },
    "Arithmetic": { m1:9.8, sd1:2.5, m2:10.1, sd2:3.1, r:0.76, rCorrected:0.83 },
    "Symbol Search": { m1:10.3, sd1:3.1, m2:11.4, sd2:3.4, r:0.75, rCorrected:0.73 },
    "Visual Puzzles": { m1:10.3, sd1:3.1, m2:11.4, sd2:3, r:0.7, rCorrected:0.68 },
    "Information": { m1:9.6, sd1:2.9, m2:10.4, sd2:2.8, r:0.86, rCorrected:0.87 },
    "Coding": { m1:10.3, sd1:2.9, m2:11.2, sd2:2.7, r:0.83, rCorrected:0.84 }
  },
  "WAIS-IV Core Subtests · Ages 55-69": {
    "Block Design": { m1:10.5, sd1:3.1, m2:11.1, sd2:2.9, r:0.77, rCorrected:0.75 },
    "Similarities": { m1:10.2, sd1:2.5, m2:10.9, sd2:2.5, r:0.81, rCorrected:0.87 },
    "Digit Span": { m1:10, sd1:2.9, m2:10.8, sd2:3.1, r:0.89, rCorrected:0.9 },
    "Matrix Reasoning": { m1:9.9, sd1:3, m2:10.7, sd2:3.3, r:0.72, rCorrected:0.72 },
    "Vocabulary": { m1:10.3, sd1:3.2, m2:10.5, sd2:3, r:0.88, rCorrected:0.86 },
    "Arithmetic": { m1:9.6, sd1:2.9, m2:10.4, sd2:2.7, r:0.8, rCorrected:0.81 },
    "Symbol Search": { m1:10, sd1:2.9, m2:11.1, sd2:3.2, r:0.8, rCorrected:0.81 },
    "Visual Puzzles": { m1:10.4, sd1:2.9, m2:10.8, sd2:3, r:0.73, rCorrected:0.75 },
    "Information": { m1:10.5, sd1:3.1, m2:11.4, sd2:3.4, r:0.92, rCorrected:0.91 },
    "Coding": { m1:10.1, sd1:2.7, m2:10.7, sd2:3.1, r:0.86, rCorrected:0.89 }
  },
  "WAIS-IV Core Subtests · Ages 70-90": {
    "Block Design": { m1:10, sd1:2.6, m2:10.4, sd2:2.5, r:0.79, rCorrected:0.84 },
    "Similarities": { m1:9.8, sd1:2.6, m2:10.3, sd2:3, r:0.84, rCorrected:0.88 },
    "Digit Span": { m1:9.8, sd1:2.9, m2:10.5, sd2:2.9, r:0.84, rCorrected:0.85 },
    "Matrix Reasoning": { m1:10, sd1:2.8, m2:10.3, sd2:2.8, r:0.73, rCorrected:0.76 },
    "Vocabulary": { m1:9.7, sd1:2.8, m2:9.7, sd2:2.8, r:0.91, rCorrected:0.92 },
    "Arithmetic": { m1:10.1, sd1:2.7, m2:10.5, sd2:2.9, r:0.8, rCorrected:0.84 },
    "Symbol Search": { m1:9.8, sd1:2.5, m2:10, sd2:2.9, r:0.8, rCorrected:0.86 },
    "Visual Puzzles": { m1:9.4, sd1:2.5, m2:10.4, sd2:3, r:0.57, rCorrected:0.7 },
    "Information": { m1:9.3, sd1:3.2, m2:9.9, sd2:3.2, r:0.93, rCorrected:0.92 },
    "Coding": { m1:9.7, sd1:2.6, m2:10.3, sd2:2.6, r:0.81, rCorrected:0.86 }
  },
  /* rInternal / rInternalByAge below are WAIS-IV Technical Manual Table 4.1,
     keyed by each normative band's LOWER bound. Symbol Search and Coding are
     deliberately absent — Table 4.1 gives them a corrected stability
     coefficient, which rCorrected already holds. See the rInternal note above. */
  "WAIS-IV Core Subtests · All Ages": {
    "Block Design": { m1:10.2, sd1:2.9, m2:11, sd2:2.8, r:0.79, rCorrected:0.8, n:298, rInternal:0.87, rInternalAgeMax:90, rInternalByAge:{16:0.88, 18:0.87, 20:0.84, 25:0.9, 30:0.91, 35:0.89, 45:0.9, 55:0.88, 65:0.87, 70:0.89, 75:0.82, 80:0.8, 85:0.86} },
    "Similarities": { m1:9.9, sd1:2.8, m2:10.4, sd2:2.8, r:0.83, rCorrected:0.87, n:298, rInternal:0.87, rInternalAgeMax:90, rInternalByAge:{16:0.81, 18:0.85, 20:0.85, 25:0.86, 30:0.87, 35:0.88, 45:0.88, 55:0.87, 65:0.88, 70:0.9, 75:0.86, 80:0.91, 85:0.91} },
    "Digit Span": { m1:10, sd1:2.9, m2:10.6, sd2:3, r:0.82, rCorrected:0.83, n:298, rInternal:0.93, rInternalAgeMax:90, rInternalByAge:{16:0.89, 18:0.92, 20:0.91, 25:0.94, 30:0.94, 35:0.94, 45:0.94, 55:0.92, 65:0.93, 70:0.94, 75:0.93, 80:0.92, 85:0.92} },
    "Matrix Reasoning": { m1:10.1, sd1:3.1, m2:10.5, sd2:3.1, r:0.76, rCorrected:0.74, n:298, rInternal:0.9, rInternalAgeMax:90, rInternalByAge:{16:0.88, 18:0.87, 20:0.88, 25:0.91, 30:0.91, 35:0.9, 45:0.9, 55:0.9, 65:0.91, 70:0.91, 75:0.94, 80:0.86, 85:0.92} },
    "Vocabulary": { m1:9.9, sd1:3, m2:10, sd2:3, r:0.9, rCorrected:0.89, n:298, rInternal:0.94, rInternalAgeMax:90, rInternalByAge:{16:0.93, 18:0.93, 20:0.94, 25:0.93, 30:0.93, 35:0.94, 45:0.94, 55:0.94, 65:0.95, 70:0.95, 75:0.94, 80:0.94, 85:0.96} },
    "Arithmetic": { m1:9.9, sd1:2.8, m2:10.4, sd2:2.9, r:0.8, rCorrected:0.83, n:298, rInternal:0.88, rInternalAgeMax:90, rInternalByAge:{16:0.89, 18:0.88, 20:0.84, 25:0.89, 30:0.9, 35:0.89, 45:0.91, 55:0.88, 65:0.89, 70:0.84, 75:0.9, 80:0.89, 85:0.86} },
    "Symbol Search": { m1:10.1, sd1:2.9, m2:11, sd2:3.3, r:0.8, rCorrected:0.81, n:298 },
    "Visual Puzzles": { m1:10, sd1:2.8, m2:10.9, sd2:3, r:0.69, rCorrected:0.74, n:298, rInternal:0.89, rInternalAgeMax:90, rInternalByAge:{16:0.9, 18:0.89, 20:0.9, 25:0.91, 30:0.9, 35:0.88, 45:0.92, 55:0.89, 65:0.92, 70:0.89, 75:0.89, 80:0.82, 85:0.78} },
    "Information": { m1:9.8, sd1:3, m2:10.5, sd2:3.2, r:0.91, rCorrected:0.9, n:298, rInternal:0.93, rInternalAgeMax:90, rInternalByAge:{16:0.89, 18:0.91, 20:0.91, 25:0.91, 30:0.91, 35:0.92, 45:0.94, 55:0.95, 65:0.94, 70:0.94, 75:0.94, 80:0.94, 85:0.96} },
    "Coding": { m1:10, sd1:2.7, m2:10.6, sd2:2.8, r:0.83, rCorrected:0.86, n:298 }
  },
  "WAIS-IV Indices · Ages 16-29": {
    "Verbal Comprehension Index": { m1:101.4, sd1:14.9, m2:103.6, sd2:15.4, r:0.95, rCorrected:0.95 },
    "Perceptual Reasoning Index": { m1:100.8, sd1:13.9, m2:105.4, sd2:14.4, r:0.84, rCorrected:0.86 },
    "Working Memory Index": { m1:100.5, sd1:14.7, m2:103.2, sd2:14.5, r:0.82, rCorrected:0.83 },
    "Processing Speed Index": { m1:100.6, sd1:15.1, m2:105.4, sd2:15.8, r:0.87, rCorrected:0.87 },
    "Full Scale IQ": { m1:101, sd1:13.8, m2:105.4, sd2:14.9, r:0.94, rCorrected:0.95 }
  },
  "WAIS-IV Indices · Ages 30-54": {
    "Verbal Comprehension Index": { m1:96.4, sd1:13.9, m2:99.1, sd2:13.6, r:0.95, rCorrected:0.96 },
    "Perceptual Reasoning Index": { m1:101.2, sd1:15.4, m2:105.7, sd2:15.9, r:0.88, rCorrected:0.87 },
    "Working Memory Index": { m1:99.5, sd1:14.6, m2:101.4, sd2:15.7, r:0.84, rCorrected:0.85 },
    "Processing Speed Index": { m1:102, sd1:13.7, m2:107.3, sd2:14.7, r:0.76, rCorrected:0.8 },
    "Full Scale IQ": { m1:99.5, sd1:14.5, m2:103.9, sd2:15.9, r:0.96, rCorrected:0.96 }
  },
  "WAIS-IV Indices · Ages 55-69": {
    "Verbal Comprehension Index": { m1:101.5, sd1:14.3, m2:104.8, sd2:14.9, r:0.94, rCorrected:0.95 },
    "Perceptual Reasoning Index": { m1:101.1, sd1:13.8, m2:104.7, sd2:14.4, r:0.87, rCorrected:0.89 },
    "Working Memory Index": { m1:98.7, sd1:13.7, m2:103, sd2:14.1, r:0.9, rCorrected:0.92 },
    "Processing Speed Index": { m1:100.3, sd1:13.2, m2:105.3, sd2:14.9, r:0.89, rCorrected:0.91 },
    "Full Scale IQ": { m1:100.6, sd1:14.5, m2:105.5, sd2:15.4, r:0.96, rCorrected:0.96 }
  },
  "WAIS-IV Indices · Ages 70-90": {
    "Verbal Comprehension Index": { m1:97.8, sd1:14, m2:99.9, sd2:15.3, r:0.95, rCorrected:0.96 },
    "Perceptual Reasoning Index": { m1:98.7, sd1:12.4, m2:101.9, sd2:12.7, r:0.8, rCorrected:0.86 },
    "Working Memory Index": { m1:99.4, sd1:13.4, m2:102.6, sd2:14.8, r:0.89, rCorrected:0.91 },
    "Processing Speed Index": { m1:98.5, sd1:12, m2:101.2, sd2:13.8, r:0.82, rCorrected:0.88 },
    "Full Scale IQ": { m1:98.1, sd1:12.7, m2:101.6, sd2:14, r:0.94, rCorrected:0.96 }
  },
  /* Table 4.1 again. PSI and FSIQ draw on the three speeded subtests, so their
     composite coefficient is a hybrid — the manual labels it internal
     consistency and publishes no other reliability for them. See above. */
  "WAIS-IV Indices · All Ages": {
    "Verbal Comprehension Index": { m1:99.3, sd1:14.4, m2:101.8, sd2:15, r:0.95, rCorrected:0.96, n:298, rInternal:0.96, rInternalAgeMax:90, rInternalByAge:{16:0.94, 18:0.96, 20:0.96, 25:0.96, 30:0.96, 35:0.96, 45:0.97, 55:0.97, 65:0.97, 70:0.97, 75:0.96, 80:0.97, 85:0.98} },
    "Perceptual Reasoning Index": { m1:100.4, sd1:13.8, m2:104.3, sd2:14.3, r:0.85, rCorrected:0.87, n:298, rInternal:0.95, rInternalAgeMax:90, rInternalByAge:{16:0.95, 18:0.94, 20:0.94, 25:0.96, 30:0.96, 35:0.95, 45:0.96, 55:0.95, 65:0.95, 70:0.95, 75:0.94, 80:0.92, 85:0.93} },
    "Working Memory Index": { m1:99.5, sd1:14, m2:102.6, sd2:14.7, r:0.87, rCorrected:0.88, n:298, rInternal:0.94, rInternalAgeMax:90, rInternalByAge:{16:0.93, 18:0.94, 20:0.92, 25:0.95, 30:0.95, 35:0.95, 45:0.95, 55:0.94, 65:0.94, 70:0.93, 75:0.95, 80:0.94, 85:0.93} },
    "Processing Speed Index": { m1:100.2, sd1:13.5, m2:104.6, sd2:14.9, r:0.84, rCorrected:0.87, n:298, rInternal:0.9, rInternalAgeMax:90, rInternalByAge:{16:0.88, 18:0.9, 20:0.9, 25:0.9, 30:0.87, 35:0.87, 45:0.87, 55:0.91, 65:0.91, 70:0.91, 75:0.92, 80:0.92, 85:0.92} },
    "Full Scale IQ": { m1:99.7, sd1:13.8, m2:104, sd2:15, r:0.95, rCorrected:0.96, n:298, rInternal:0.98, rInternalAgeMax:90, rInternalByAge:{16:0.97, 18:0.98, 20:0.98, 25:0.98, 30:0.98, 35:0.98, 45:0.98, 55:0.98, 65:0.98, 70:0.98, 75:0.98, 80:0.98, 85:0.98} }
  },
  "WAIS-IV Process Scores · Ages 16-29": {
    "Block Design No Time Bonus": { m1:10.3, sd1:2.7, m2:11.3, sd2:2.4, r:0.77, rCorrected:0.81 },
    "Digit Span Forward": { m1:10, sd1:2.7, m2:10.4, sd2:3.3, r:0.67, rCorrected:0.73 },
    "Digit Span Backward": { m1:10.4, sd1:2.5, m2:11.2, sd2:3, r:0.51, rCorrected:0.66 },
    "Digit Span Sequencing": { m1:9.9, sd1:3.2, m2:10.5, sd2:3.2, r:0.65, rCorrected:0.6 }
  },
  "WAIS-IV Process Scores · Ages 30-54": {
    "Block Design No Time Bonus": { m1:10.6, sd1:2.9, m2:11.4, sd2:2.9, r:0.7, rCorrected:0.72 },
    "Digit Span Forward": { m1:10.1, sd1:3, m2:10.2, sd2:3, r:0.71, rCorrected:0.71 },
    "Digit Span Backward": { m1:10, sd1:3.2, m2:10.4, sd2:3.3, r:0.73, rCorrected:0.69 },
    "Digit Span Sequencing": { m1:10.2, sd1:3, m2:10.5, sd2:2.7, r:0.7, rCorrected:0.7 }
  },
  "WAIS-IV Process Scores · Ages 55-69": {
    "Block Design No Time Bonus": { m1:10.3, sd1:3.3, m2:11.1, sd2:3.3, r:0.75, rCorrected:0.7 },
    "Digit Span Forward": { m1:9.8, sd1:2.7, m2:10.2, sd2:2.9, r:0.76, rCorrected:0.81 },
    "Digit Span Backward": { m1:10.3, sd1:3, m2:11.1, sd2:3.2, r:0.77, rCorrected:0.77 },
    "Digit Span Sequencing": { m1:9.8, sd1:2.2, m2:10.6, sd2:2.5, r:0.71, rCorrected:0.84 }
  },
  "WAIS-IV Process Scores · Ages 70-90": {
    "Block Design No Time Bonus": { m1:10, sd1:2.6, m2:10.4, sd2:2.5, r:0.8, rCorrected:0.85 },
    "Digit Span Forward": { m1:9.9, sd1:2.9, m2:10, sd2:2.9, r:0.81, rCorrected:0.82 },
    "Digit Span Backward": { m1:10, sd1:3, m2:10.3, sd2:3, r:0.71, rCorrected:0.71 },
    "Digit Span Sequencing": { m1:9.6, sd1:3.1, m2:10.7, sd2:2.8, r:0.72, rCorrected:0.7 }
  },
  "WAIS-IV Process Scores · All Ages": {
    "Block Design No Time Bonus": { m1:10.3, sd1:2.9, m2:11, sd2:2.8, r:0.76, rCorrected:0.78, n:298, rInternal:0.86, rInternalAgeMax:90, rInternalByAge:{16:0.87, 18:0.86, 20:0.81, 25:0.88, 30:0.89, 35:0.87, 45:0.87, 55:0.86, 65:0.85, 70:0.88, 75:0.82, 80:0.8, 85:0.86} },
    "Digit Span Forward": { m1:9.9, sd1:2.8, m2:10.2, sd2:3, r:0.74, rCorrected:0.77, n:298, rInternal:0.81, rInternalAgeMax:90, rInternalByAge:{16:0.77, 18:0.76, 20:0.8, 25:0.85, 30:0.77, 35:0.84, 45:0.83, 55:0.77, 65:0.79, 70:0.88, 75:0.81, 80:0.81, 85:0.78} },
    "Digit Span Backward": { m1:10.2, sd1:2.9, m2:10.7, sd2:3.1, r:0.69, rCorrected:0.71, n:298, rInternal:0.82, rInternalAgeMax:90, rInternalByAge:{16:0.79, 18:0.8, 20:0.8, 25:0.84, 30:0.84, 35:0.86, 45:0.86, 55:0.82, 65:0.78, 70:0.79, 75:0.8, 80:0.77, 85:0.82} },
    "Digit Span Sequencing": { m1:9.9, sd1:2.9, m2:10.6, sd2:2.8, r:0.7, rCorrected:0.72, n:298, rInternal:0.83, rInternalAgeMax:90, rInternalByAge:{16:0.73, 18:0.79, 20:0.82, 25:0.81, 30:0.8, 35:0.83, 45:0.82, 55:0.79, 65:0.81, 70:0.86, 75:0.86, 80:0.84, 85:0.92} }
  },
  "WAIS-IV Supplementary Subtests · Ages 16-29": {
    "Letter-Number Sequencing": { m1:10.1, sd1:2.5, m2:10.8, sd2:3.4, r:0.76, rCorrected:0.83 },
    "Figure Weights": { m1:10.4, sd1:3, m2:11.4, sd2:3.4, r:0.76, rCorrected:0.76 },
    "Comprehension": { m1:10.2, sd1:3.2, m2:10.2, sd2:2.9, r:0.86, rCorrected:0.84 },
    "Cancellation": { m1:9.7, sd1:2.9, m2:10.8, sd2:3.3, r:0.8, rCorrected:0.81 },
    "Picture Completion": { m1:10, sd1:2.7, m2:12.4, sd2:3.2, r:0.74, rCorrected:0.79 }
  },
  "WAIS-IV Supplementary Subtests · Ages 30-54": {
    "Letter-Number Sequencing": { m1:10, sd1:2.8, m2:10.5, sd2:2.9, r:0.81, rCorrected:0.83 },
    "Figure Weights": { m1:9.7, sd1:2.7, m2:10.7, sd2:3.2, r:0.77, rCorrected:0.81 },
    "Comprehension": { m1:9.3, sd1:2.8, m2:9.5, sd2:2.8, r:0.89, rCorrected:0.9 },
    "Cancellation": { m1:10.8, sd1:2.8, m2:11.3, sd2:2.8, r:0.67, rCorrected:0.71 },
    "Picture Completion": { m1:10.1, sd1:3, m2:12.4, sd2:3.3, r:0.68, rCorrected:0.68 }
  },
  "WAIS-IV Supplementary Subtests · Ages 55-69": {
    "Letter-Number Sequencing": { m1:10, sd1:2.8, m2:10.1, sd2:2.8, r:0.7, rCorrected:0.74 },
    "Figure Weights": { m1:9.7, sd1:3.2, m2:10.3, sd2:2.9, r:0.76, rCorrected:0.73 },
    "Comprehension": { m1:10.4, sd1:2.9, m2:10.4, sd2:2.9, r:0.84, rCorrected:0.85 },
    "Cancellation": { m1:10.1, sd1:2.6, m2:10.3, sd2:2.7, r:0.74, rCorrected:0.8 },
    "Picture Completion": { m1:10, sd1:2.8, m2:11.8, sd2:3.6, r:0.78, rCorrected:0.81 }
  },
  "WAIS-IV Supplementary Subtests · Ages 70-90": {
    "Comprehension": { m1:10.2, sd1:2.9, m2:10.6, sd2:3.1, r:0.85, rCorrected:0.86 },
    "Picture Completion": { m1:9.4, sd1:3, m2:10.6, sd2:2.9, r:0.77, rCorrected:0.77 }
  },
  /* Letter-Number Sequencing and Figure Weights carry rInternalAgeMax 69, not
     90: Table 4.1 prints a dash above the 65-69 band for both, those subtests
     being normed to 69. Cancellation, normed to 69 as well, is one of the three
     speeded exclusions and takes no field at all. */
  "WAIS-IV Supplementary Subtests · All Ages": {
    "Letter-Number Sequencing": { m1:10.1, sd1:2.7, m2:10.5, sd2:3.1, r:0.76, rCorrected:0.8, n:298, rInternal:0.88, rInternalAgeMax:69, rInternalByAge:{16:0.9, 18:0.9, 20:0.85, 25:0.88, 30:0.91, 35:0.86, 45:0.88, 55:0.87, 65:0.88} },
    "Figure Weights": { m1:10, sd1:3, m2:10.8, sd2:3.2, r:0.76, rCorrected:0.77, n:298, rInternal:0.9, rInternalAgeMax:69, rInternalByAge:{16:0.9, 18:0.92, 20:0.88, 25:0.91, 30:0.92, 35:0.89, 45:0.91, 55:0.89, 65:0.9} },
    "Comprehension": { m1:10, sd1:3, m2:10.2, sd2:2.9, r:0.86, rCorrected:0.86, n:298, rInternal:0.87, rInternalAgeMax:90, rInternalByAge:{16:0.82, 18:0.87, 20:0.87, 25:0.88, 30:0.85, 35:0.89, 45:0.85, 55:0.87, 65:0.87, 70:0.82, 75:0.87, 80:0.87, 85:0.9} },
    "Cancellation": { m1:10.2, sd1:2.8, m2:10.8, sd2:3, r:0.74, rCorrected:0.78, n:298 },
    "Picture Completion": { m1:9.9, sd1:2.9, m2:11.8, sd2:3.3, r:0.74, rCorrected:0.77, n:298, rInternal:0.84, rInternalAgeMax:90, rInternalByAge:{16:0.8, 18:0.84, 20:0.82, 25:0.82, 30:0.86, 35:0.84, 45:0.83, 55:0.86, 65:0.85, 70:0.89, 75:0.86, 80:0.83, 85:0.82} }
  },
  /* ---------------------------------------------------------------
     WISC-V — Technical and Interpretive Manual, Tables 4.1 and 4.4,
     pp. 56-62. The finest age lookup in this database: reliability by
     SINGLE YEAR of age, 6 to 16, rather than by band.

     Verified arithmetically, which matters more here than anywhere
     else. WISC-V is the ONLY family in normDB with no rCorrected on
     any entry, so the transcription proof used for WAIS-IV, RBANS and
     WMS-IV — matching the manual's stability rows against this
     database's stored rCorrected — is unavailable. Table 4.4 does that
     work instead: all 242 published SEM cells equal SD sqrt(1 - rxx)
     from Table 4.1 at the printed 2dp, on 3 and 15.

     AND THE MANUAL'S WORKED EXAMPLE REPRODUCES THROUGH THE SHIPPED
     RENDERER (p. 62): a 6-year-old with an FSIQ of 108 gets 102-114 at
     95% and 103-113 at 90%. check.js section 31 drives the real
     getBatteryCiHtml to prove it, the same standard CVLT-C is held to.

     THREE SUBTESTS ARE STABILITY, and the manual names them (p. 56):
     "The split-half coefficient is not a proper reliability estimate
     for Coding, Symbol Search, Cancellation, Naming Speed Literacy,
     Naming Speed Quantity, Immediate Symbol Translation, or Delayed
     Symbol Translation. Therefore, test-retest coefficients were
     used... corrected for the normative sample's variability." Only
     the first three are held here; the other four are complementary
     subtests this app does not carry. They take rStability.

     PSI is a hybrid and is included anyway, exactly as WAIS-IV's is:
     p. 61 concedes its average "is based on test-retest
     reliabilities", but the manual computes and labels every composite
     coefficient as internal consistency and publishes no other
     reliability for them.

     THE FOURTH PUBLISHER TO STATE THE SD RULE OUTRIGHT. Table 4.4's
     note: "The reliability coefficients shown in Table 4.1 and the
     population standard deviations (i.e., 3 for the scaled scores and
     15 for standard scores) were used to compute the SEMs." After
     CVLT-3, WAIS-IV Table 4.3 and WMS-IV. Pair the coefficient with
     the NORMATIVE SD, never with sd1.

     ON THE INTERVAL CONVENTION. This manual documents both methods and
     the app matches one of them deliberately: its Tables A.2-A.7 build
     intervals around the ESTIMATED TRUE score using the SEE (Dudek
     1979), while p. 62 gives the observed-score form for those who
     "prefer to calculate confidence intervals... in the most
     parsimonious manner". This app uses the observed-score form and
     says so in Methods & References, which is why the worked example
     above is the right thing to reproduce.

     rInternalAgeMax is 16, the top of the normed range, so an age
     outside 6-16 takes the published average rather than silently
     re-reading age 16.
     --------------------------------------------------------------- */
  "WISC-V Indices · All Ages": {
    "Verbal Comprehension Index": { m1:98.5, sd1:12.8, m2:101.6, sd2:13, r:0.91, rCorrected:0.94, n:215, rInternal:0.92, rInternalAgeMax:16, rInternalByAge:{6:0.91, 7:0.92, 8:0.9, 9:0.93, 10:0.93, 11:0.9, 12:0.94, 13:0.93, 14:0.92, 15:0.92, 16:0.93} },
    "Visuospatial Index": { m1:98.6, sd1:14.7, m2:105.3, sd2:15.1, r:0.84, rCorrected:0.84, n:217, rInternal:0.92, rInternalAgeMax:16, rInternalByAge:{6:0.91, 7:0.92, 8:0.92, 9:0.91, 10:0.91, 11:0.91, 12:0.93, 13:0.91, 14:0.9, 15:0.93, 16:0.92} },
    "Fluid Reasoning Index": { m1:98.7, sd1:13.6, m2:103.6, sd2:12.9, r:0.68, rCorrected:0.75, n:217, rInternal:0.93, rInternalAgeMax:16, rInternalByAge:{6:0.93, 7:0.94, 8:0.94, 9:0.93, 10:0.93, 11:0.92, 12:0.95, 13:0.93, 14:0.93, 15:0.93, 16:0.93} },
    "Working Memory Index": { m1:98.5, sd1:13.8, m2:100.9, sd2:13.8, r:0.79, rCorrected:0.82, n:217, rInternal:0.92, rInternalAgeMax:16, rInternalByAge:{6:0.92, 7:0.91, 8:0.92, 9:0.92, 10:0.91, 11:0.92, 12:0.93, 13:0.92, 14:0.92, 15:0.92, 16:0.92} },
    "Processing Speed Index": { m1:100.3, sd1:14.3, m2:108.2, sd2:16, r:0.81, rCorrected:0.83, n:213, rInternal:0.88, rInternalAgeMax:16, rInternalByAge:{6:0.88, 7:0.88, 8:0.86, 9:0.87, 10:0.87, 11:0.88, 12:0.84, 13:0.84, 14:0.91, 15:0.91, 16:0.92} },
    "Full Scale IQ": { m1:98.3, sd1:13.7, m2:104.3, sd2:13.8, r:0.91, rCorrected:0.92, n:212, rInternal:0.96, rInternalAgeMax:16, rInternalByAge:{6:0.96, 7:0.96, 8:0.96, 9:0.96, 10:0.96, 11:0.96, 12:0.97, 13:0.96, 14:0.96, 15:0.97, 16:0.97} },
    "Quantitative Reasoning Index": { m1:99.2, sd1:13.4, m2:102.4, sd2:13.7, r:0.76, rCorrected:0.81, n:216, rInternal:0.95, rInternalAgeMax:16, rInternalByAge:{6:0.92, 7:0.94, 8:0.94, 9:0.94, 10:0.96, 11:0.95, 12:0.96, 13:0.95, 14:0.96, 15:0.94, 16:0.95} },
    "Auditory Working Memory Index": { m1:98.7, sd1:13.3, m2:100.9, sd2:14.6, r:0.85, rCorrected:0.88, n:217, rInternal:0.93, rInternalAgeMax:16, rInternalByAge:{6:0.96, 7:0.94, 8:0.92, 9:0.93, 10:0.91, 11:0.92, 12:0.94, 13:0.93, 14:0.95, 15:0.93, 16:0.92} },
    "Nonverbal Index": { m1:98.5, sd1:14, m2:105.5, sd2:13.8, r:0.86, rCorrected:0.88, n:216, rInternal:0.95, rInternalAgeMax:16, rInternalByAge:{6:0.95, 7:0.95, 8:0.95, 9:0.95, 10:0.95, 11:0.95, 12:0.96, 13:0.95, 14:0.96, 15:0.96, 16:0.96} },
    "General Ability Index": { m1:98, sd1:13.7, m2:103.6, sd2:13.3, r:0.89, rCorrected:0.91, n:213, rInternal:0.96, rInternalAgeMax:16, rInternalByAge:{6:0.95, 7:0.96, 8:0.96, 9:0.96, 10:0.96, 11:0.95, 12:0.97, 13:0.96, 14:0.95, 15:0.96, 16:0.96} },
    "Cognitive Proficiency Index": { m1:99.3, sd1:14, m2:105.5, sd2:14.8, r:0.84, rCorrected:0.86, n:212, rInternal:0.93, rInternalAgeMax:16, rInternalByAge:{6:0.92, 7:0.92, 8:0.92, 9:0.92, 10:0.92, 11:0.92, 12:0.92, 13:0.91, 14:0.94, 15:0.94, 16:0.94} }
  },
  /* ---------------------------------------------------------------
     WISC-V Process Scores. Mirrors WAIS-IV Process Scores, which holds
     the adult equivalents of the same measures, so the two sit
     consistently beside each other in every dropdown.

     Retest descriptives and coefficients: Technical and Interpretive
     Manual Table 4.7 (all ages). Reliability by single year of age:
     Table 4.1, verified against the SEMs of Table 4.4 — all 132 cells
     for these twelve measures equal SD sqrt(1 - rxx) at 2dp, on 3 and
     15, so they are admitted on the same footing as the rest of the
     family rather than on retest data alone.

     CANCELLATION RANDOM AND STRUCTURED TAKE rStability. The manual
     names Cancellation among the subtests for which "the split-half
     coefficient is not a proper reliability estimate", and its two
     process scores carry the signature: a stability coefficient is
     broadcast from the retest study's coarse bands, so it repeats
     across single years. CAr shows 3 distinct values across the 11
     bands and CAs shows 2, against 7-9 for the five internal-
     consistency rows here and 4-5 for the known-stability Coding,
     Symbol Search and Cancellation.
     --------------------------------------------------------------- */
  "WISC-V Process Scores · All Ages": {
    "Block Design No Time Bonus": { m1:9.6, sd1:2.8, m2:10.8, sd2:3.2, r:0.78, rCorrected:0.81, n:207, rInternal:0.83, rInternalAgeMax:16, rInternalByAge:{6:0.84, 7:0.86, 8:0.87, 9:0.82, 10:0.79, 11:0.82, 12:0.84, 13:0.8, 14:0.78, 15:0.85, 16:0.82} },
    "Block Design Partial": { m1:9.8, sd1:3.1, m2:10.8, sd2:3.1, r:0.84, rCorrected:0.84, n:208, rInternal:0.88, rInternalAgeMax:16, rInternalByAge:{6:0.89, 7:0.91, 8:0.93, 9:0.89, 10:0.87, 11:0.87, 12:0.88, 13:0.85, 14:0.84, 15:0.87, 16:0.88} },
    "Digit Span Forward": { m1:9.9, sd1:2.7, m2:10.2, sd2:3, r:0.78, rCorrected:0.82, n:213, rInternal:0.81, rInternalAgeMax:16, rInternalByAge:{6:0.8, 7:0.76, 8:0.85, 9:0.75, 10:0.76, 11:0.82, 12:0.82, 13:0.85, 14:0.79, 15:0.84, 16:0.83} },
    "Digit Span Backward": { m1:9.9, sd1:2.8, m2:10.2, sd2:2.9, r:0.72, rCorrected:0.76, n:201, rInternal:0.8, rInternalAgeMax:16, rInternalByAge:{6:0.84, 7:0.81, 8:0.79, 9:0.78, 10:0.75, 11:0.83, 12:0.82, 13:0.78, 14:0.8, 15:0.8, 16:0.82} },
    "Digit Span Sequencing": { m1:9.5, sd1:2.7, m2:10.1, sd2:2.9, r:0.74, rCorrected:0.79, n:200, rInternal:0.82, rInternalAgeMax:16, rInternalByAge:{6:0.9, 7:0.81, 8:0.78, 9:0.79, 10:0.83, 11:0.76, 12:0.84, 13:0.77, 14:0.78, 15:0.83, 16:0.85} },
    "Cancellation Random": { m1:9.9, sd1:2.9, m2:11.1, sd2:2.9, r:0.78, rCorrected:0.81, n:200, rStability:0.81, rStabilityAgeMax:16, rStabilityByAge:{6:0.8, 7:0.8, 8:0.82, 9:0.82, 10:0.81, 11:0.81, 12:0.8, 13:0.8, 14:0.81, 15:0.81, 16:0.81} },
    "Cancellation Structured": { m1:9.9, sd1:2.8, m2:10.9, sd2:3.1, r:0.78, rCorrected:0.82, n:209, rStability:0.82, rStabilityAgeMax:16, rStabilityByAge:{6:0.8, 7:0.8, 8:0.82, 9:0.82, 10:0.82, 11:0.82, 12:0.82, 13:0.82, 14:0.82, 15:0.82, 16:0.82} }
  },
  "WISC-V Subtests · All Ages": {
    "Similarities": { m1:9.8, sd1:2.5, m2:10.6, sd2:2.5, r:0.82, rCorrected:0.88, n:213, rInternal:0.87, rInternalAgeMax:16, rInternalByAge:{6:0.89, 7:0.87, 8:0.85, 9:0.88, 10:0.88, 11:0.81, 12:0.87, 13:0.89, 14:0.85, 15:0.85, 16:0.87} },
    "Vocabulary": { m1:9.6, sd1:2.8, m2:10, sd2:2.8, r:0.89, rCorrected:0.9, n:217, rInternal:0.87, rInternalAgeMax:16, rInternalByAge:{6:0.83, 7:0.86, 8:0.83, 9:0.87, 10:0.87, 11:0.87, 12:0.91, 13:0.86, 14:0.88, 15:0.89, 16:0.9} },
    "Information": { m1:9.7, sd1:2.7, m2:10.3, sd2:2.7, r:0.85, rCorrected:0.88, n:218, rInternal:0.86, rInternalAgeMax:16, rInternalByAge:{6:0.82, 7:0.86, 8:0.81, 9:0.86, 10:0.82, 11:0.81, 12:0.89, 13:0.88, 14:0.85, 15:0.88, 16:0.9} },
    "Comprehension": { m1:10, sd1:2.9, m2:10.2, sd2:2.8, r:0.81, rCorrected:0.83, n:214, rInternal:0.83, rInternalAgeMax:16, rInternalByAge:{6:0.76, 7:0.86, 8:0.8, 9:0.84, 10:0.82, 11:0.79, 12:0.87, 13:0.82, 14:0.88, 15:0.83, 16:0.8} },
    "Block Design": { m1:9.6, sd1:2.9, m2:10.8, sd2:3.1, r:0.79, rCorrected:0.81, n:208, rInternal:0.84, rInternalAgeMax:16, rInternalByAge:{6:0.84, 7:0.86, 8:0.88, 9:0.83, 10:0.81, 11:0.83, 12:0.85, 13:0.82, 14:0.8, 15:0.86, 16:0.85} },
    "Visual Puzzles": { m1:9.9, sd1:2.8, m2:11, sd2:2.9, r:0.78, rCorrected:0.8, n:210, rInternal:0.89, rInternalAgeMax:16, rInternalByAge:{6:0.89, 7:0.9, 8:0.87, 9:0.9, 10:0.89, 11:0.88, 12:0.9, 13:0.89, 14:0.89, 15:0.92, 16:0.9} },
    "Matrix Reasoning": { m1:9.6, sd1:2.4, m2:10.6, sd2:2.6, r:0.65, rCorrected:0.78, n:202, rInternal:0.87, rInternalAgeMax:16, rInternalByAge:{6:0.89, 7:0.88, 8:0.89, 9:0.87, 10:0.85, 11:0.82, 12:0.9, 13:0.84, 14:0.84, 15:0.86, 16:0.86} },
    "Figure Weights": { m1:10, sd1:2.6, m2:10.5, sd2:2.6, r:0.76, rCorrected:0.82, n:204, rInternal:0.94, rInternalAgeMax:16, rInternalByAge:{6:0.91, 7:0.94, 8:0.94, 9:0.94, 10:0.96, 11:0.94, 12:0.95, 13:0.95, 14:0.94, 15:0.94, 16:0.93} },
    "Picture Concepts": { m1:9.8, sd1:2.7, m2:10.7, sd2:2.9, r:0.63, rCorrected:0.71, n:203, rInternal:0.83, rInternalAgeMax:16, rInternalByAge:{6:0.88, 7:0.82, 8:0.83, 9:0.82, 10:0.85, 11:0.85, 12:0.83, 13:0.8, 14:0.81, 15:0.8, 16:0.78} },
    "Arithmetic": { m1:9.8, sd1:2.5, m2:10.2, sd2:2.6, r:0.75, rCorrected:0.84, n:205, rInternal:0.9, rInternalAgeMax:16, rInternalByAge:{6:0.88, 7:0.9, 8:0.88, 9:0.9, 10:0.92, 11:0.9, 12:0.91, 13:0.91, 14:0.92, 15:0.89, 16:0.92} },
    "Digit Span": { m1:9.8, sd1:2.8, m2:10.1, sd2:3, r:0.79, rCorrected:0.82, n:214, rInternal:0.91, rInternalAgeMax:16, rInternalByAge:{6:0.92, 7:0.92, 8:0.9, 9:0.89, 10:0.9, 11:0.91, 12:0.93, 13:0.92, 14:0.92, 15:0.92, 16:0.92} },
    "Picture Span": { m1:9.7, sd1:2.5, m2:10.1, sd2:2.6, r:0.72, rCorrected:0.8, n:208, rInternal:0.85, rInternalAgeMax:16, rInternalByAge:{6:0.86, 7:0.82, 8:0.87, 9:0.87, 10:0.83, 11:0.84, 12:0.86, 13:0.85, 14:0.83, 15:0.84, 16:0.84} },
    "Letter-Number Sequencing": { m1:9.8, sd1:2.7, m2:10.2, sd2:2.8, r:0.77, rCorrected:0.82, n:212, rInternal:0.86, rInternalAgeMax:16, rInternalByAge:{6:0.93, 7:0.9, 8:0.83, 9:0.87, 10:0.8, 11:0.82, 12:0.85, 13:0.86, 14:0.89, 15:0.86, 16:0.82} },
    "Coding": { m1:10, sd1:2.9, m2:11.3, sd2:3.1, r:0.79, rCorrected:0.81, n:216, rStability:0.82, rStabilityAgeMax:16, rStabilityByAge:{6:0.78, 7:0.78, 8:0.79, 9:0.79, 10:0.81, 11:0.81, 12:0.82, 13:0.82, 14:0.86, 15:0.86, 16:0.86} },
    "Symbol Search": { m1:10, sd1:2.7, m2:11.5, sd2:3.2, r:0.76, rCorrected:0.8, n:209, rStability:0.81, rStabilityAgeMax:16, rStabilityByAge:{6:0.83, 7:0.83, 8:0.8, 9:0.8, 10:0.79, 11:0.79, 12:0.67, 13:0.67, 14:0.87, 15:0.87, 16:0.87} },
    "Cancellation": { m1:9.8, sd1:2.9, m2:11.1, sd2:3.2, r:0.79, rCorrected:0.82, n:209, rStability:0.82, rStabilityAgeMax:16, rStabilityByAge:{6:0.8, 7:0.8, 8:0.83, 9:0.83, 10:0.84, 11:0.84, 12:0.81, 13:0.81, 14:0.81, 15:0.81, 16:0.81} }
  },

  /* ---- WAIS-IV Longest Span, process level -------------------------------
     Source: WAIS-IV Administration and Scoring Manual (GB), Tables C.4 and
     C.5 — cumulative percentages of the normative sample obtaining each raw
     longest span, by age group.

     baseRates[x] is the percentage scoring x OR HIGHER. That direction is not
     an assumption: for a whole-number score E[X] = the sum of P(X >= x), and
     reconstructing the mean that way reproduces every printed mean across all
     51 measure-bands to within 0.055, which is rounding. check.js re-runs it.

     singleAdministration marks these as having NO retest data — there are no
     reliability coefficients published for longest span, so they carry no
     m2/sd2/r and must stay out of the Change Analysis pages.

     The age bands are the NORMATIVE ones and deliberately do not match the
     "WAIS-IV Process Scores" families, which use the retest study bands
     (16-29, 30-54, 55-69, 70-90). They cannot be merged for that reason.
     Letter-Number Sequence stops at 65-69 in the published table.
     -------------------------------------------------------------------- */
  "WAIS-IV Longest Span (Process) · Ages 16-17": {
    "Longest Digit Span Forward": { m1:6.5, sd1:1.2, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:5.5,8:19,7:50.5,6:77,5:97.5,4:99.5,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.8, sd1:1.2, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3.5,7:9,6:22.5,5:53,4:89,3:99,2:99.5,0:100} },
    "Longest Digit Span Sequence": { m1:5.9, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:1.5,8:10,7:26.5,6:71.5,5:85.5,4:98.5,3:99.5,2:99.5,0:100} },
    "Longest Letter-Number Sequence": { m1:5.3, sd1:1.2, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:4.5,7:14.5,6:44,5:76.5,4:96,3:99,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 18-19": {
    "Longest Digit Span Forward": { m1:6.8, sd1:1.2, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:9,8:26.5,7:64.5,6:84.5,5:98,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:5, sd1:1.4, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3.5,7:15.5,6:34.5,5:57.5,4:90.5,3:97.5,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:6, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:5,8:13,7:28.5,6:73,5:88.5,4:97.5,3:99,2:100,0:100} },
    "Longest Letter-Number Sequence": { m1:5.6, sd1:1.2, median:6,
      metric:'raw', singleAdministration:true, baseRates:{8:5,7:21,6:54,5:84.5,4:95.5,3:98.5,2:99.5,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 20-24": {
    "Longest Digit Span Forward": { m1:6.9, sd1:1.3, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:11,8:33.5,7:63.5,6:82,5:97.5,4:98.5,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:5.1, sd1:1.3, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3.5,7:15,6:32,5:64.5,4:92,3:99,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:6, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:7,8:12,7:27.5,6:74,5:87.5,4:97,3:98.5,2:99.5,0:100} },
    "Longest Letter-Number Sequence": { m1:5.6, sd1:1.2, median:6,
      metric:'raw', singleAdministration:true, baseRates:{8:4.5,7:21.5,6:52,5:83.5,4:97,3:99,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 25-29": {
    "Longest Digit Span Forward": { m1:6.9, sd1:1.3, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:17.5,8:27,7:62,6:86,5:97,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.9, sd1:1.5, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:7,7:18.5,6:32,5:56,4:83.5,3:97.5,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:6, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:2.5,8:10.5,7:30,6:72.5,5:86.5,4:97.5,3:99.5,2:100,0:100} },
    "Longest Letter-Number Sequence": { m1:5.5, sd1:1.2, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:6,7:19.5,6:49,5:84,4:95,3:99.5,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 30-34": {
    "Longest Digit Span Forward": { m1:6.8, sd1:1.4, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:11,8:30,7:60.5,6:82.5,5:97.5,4:99,3:99.5,2:99.5,0:100} },
    "Longest Digit Span Backward": { m1:5.1, sd1:1.5, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:6.5,7:18.5,6:32.5,5:63,4:90,3:97,2:99.5,0:100} },
    "Longest Digit Span Sequence": { m1:6, sd1:1.5, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:5.5,8:11,7:32,6:68.5,5:86.5,4:96,3:98,2:99,0:100} },
    "Longest Letter-Number Sequence": { m1:5.5, sd1:1.3, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:8,7:23,6:49,5:81,4:97,3:98.5,2:99,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 35-44": {
    "Longest Digit Span Forward": { m1:6.8, sd1:1.4, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:14,8:32,7:61.5,6:81,5:95.5,4:99.5,3:99.5,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.9, sd1:1.5, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:6.5,7:15.5,6:36,5:56.5,4:85.5,3:96.5,2:99,0:100} },
    "Longest Digit Span Sequence": { m1:5.9, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:3.5,8:11,7:25,6:70.5,5:83,4:96.5,3:98.5,2:99.5,0:100} },
    "Longest Letter-Number Sequence": { m1:5.4, sd1:1.2, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:5.5,7:15,6:44.5,5:77,4:95.5,3:99,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 45-54": {
    "Longest Digit Span Forward": { m1:6.8, sd1:1.4, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:14,8:30,7:61.5,6:79,5:96,4:99.5,3:99.5,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.8, sd1:1.5, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3,7:16.5,6:31,5:55.5,4:82,3:96,2:99.5,0:100} },
    "Longest Digit Span Sequence": { m1:5.7, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:2,8:7.5,7:23,6:66,5:82,4:95.5,3:97.5,2:99.5,0:100} },
    "Longest Letter-Number Sequence": { m1:5.4, sd1:1.3, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:6,7:19,6:46,5:76,4:93,3:97.5,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 55-64": {
    "Longest Digit Span Forward": { m1:6.5, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:7,8:20,7:48,6:75.5,5:94.5,4:100,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.5, sd1:1.3, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:2,7:8.5,6:23.5,5:45,4:79.5,3:97,2:99.5,0:100} },
    "Longest Digit Span Sequence": { m1:5.6, sd1:1.2, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:1.5,8:4.5,7:15,6:63,5:81.5,4:96.5,3:98,2:100,0:100} },
    "Longest Letter-Number Sequence": { m1:5.2, sd1:1.1, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:2.5,7:10.5,6:40.5,5:74.5,4:92.5,3:99.5,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 65-69": {
    "Longest Digit Span Forward": { m1:6.5, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:6,8:26,7:49.5,6:70,5:95.5,4:100,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.5, sd1:1.3, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:0.5,7:10.5,6:21.5,5:46,4:76.5,3:96,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:5.5, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:0.5,8:4.5,7:15.5,6:59,5:78,4:92.5,3:96,2:99.5,0:100} },
    "Longest Letter-Number Sequence": { m1:5, sd1:1.2, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3,7:10,6:34,5:68,4:88.5,3:99.5,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 70-74": {
    "Longest Digit Span Forward": { m1:6.2, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:8,8:17,7:41,6:68,5:90,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.3, sd1:1.3, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:2,7:7,6:12,5:43,4:74,3:95,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:5.2, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:0,8:4,7:8,6:52,5:74,4:86,3:93,2:100,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 75-79": {
    "Longest Digit Span Forward": { m1:6.3, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:6,8:18,7:41,6:69,5:95,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.4, sd1:1.3, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:1,7:7,6:20,5:45,4:76,3:95,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:5, sd1:1.5, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:0,8:3,7:8,6:51,5:68,4:84,3:92,2:99,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 80-84": {
    "Longest Digit Span Forward": { m1:6.1, sd1:1.3, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:5,8:15,7:35,6:66,5:93,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.1, sd1:1.1, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:0,7:1,6:9,5:33,4:77,3:90,2:99,0:100} },
    "Longest Digit Span Sequence": { m1:4.7, sd1:1.4, median:5,
      metric:'raw', singleAdministration:true, baseRates:{9:0,8:0,7:3,6:36,5:59,4:84,3:91,2:99,0:100} }
  },
  "WAIS-IV Longest Span (Process) · Ages 85-90": {
    "Longest Digit Span Forward": { m1:6, sd1:1.2, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:3,8:11,7:26,6:63,5:94,4:99,3:100,2:100,0:100} },
    "Longest Digit Span Backward": { m1:3.8, sd1:1.2, median:4,
      metric:'raw', singleAdministration:true, baseRates:{8:1,7:2,6:8,5:27,4:55,3:89,2:100,0:100} },
    "Longest Digit Span Sequence": { m1:4.3, sd1:1.6, median:4.5,
      metric:'raw', singleAdministration:true, baseRates:{9:0,8:0,7:2,6:31,5:50,4:67,3:80,2:98,0:100} }
  },
  "WAIS-IV Longest Span (Process) · All Ages": {
    "Longest Digit Span Forward": { m1:6.6, sd1:1.3, median:7,
      metric:'raw', singleAdministration:true, baseRates:{9:9.6,8:25,7:53.9,6:77.3,5:95.9,4:99.3,3:99.9,2:100,0:100} },
    "Longest Digit Span Backward": { m1:4.7, sd1:1.4, median:5,
      metric:'raw', singleAdministration:true, baseRates:{8:3.5,7:12.4,6:26.4,5:51.9,4:82.7,3:96.4,2:99.7,0:100} },
    "Longest Digit Span Sequence": { m1:5.7, sd1:1.4, median:6,
      metric:'raw', singleAdministration:true, baseRates:{9:2.6,8:8,7:21.2,6:63.9,5:80.4,4:93.5,3:96.6,2:99.5,0:100} }
  },
  /* ---------------------------------------------------------------
     WMS-IV — TWO BATTERIES, NOT TWO AGE BANDS. This is the only family
     here whose age-band groups are separate norm sets rather than
     separate retest samples, which is why its entries carry
     separateBattery:true and why Score Tables keeps the band
     SELECTABLE instead of collapsing to one canonical group.

     The Adult battery is normed 16-69 and has 15 subtests and 5
     indices; the Older Adult battery is normed 65-90 and drops
     Designs, Spatial Addition and the Visual Working Memory Index,
     leaving 8 and 4. Ages 65-69 are normed in BOTH, with different
     coefficients — Logical Memory I is .79 in the adult battery and
     .83 in the older adult — so age cannot choose between them. The
     clinician chose when they decided which battery to administer.
     Collapsing them showed an 80-year-old the adult measure list and
     its coefficients; 10 of the 12 shared measures printed a wrong
     interval. check.js section 30 asserts the overlap exists AND
     disagrees, because that is the whole argument for the exemption.

     SOURCE: WMS-IV Technical and Interpretive Manual (GB), chapter 3,
     Tables 3.1 and 3.3, pp. 44-46. It clears the bar arithmetically,
     not on stated method: all 240 published SEM cells in Table 3.3
     equal SD sqrt(1 - rxx) from Table 3.1 at the printed 2dp, on 3
     for subtests and 15 for indices.

     ONE MEASURE IS STABILITY, AND CARRIES footnote b. Verbal Paired
     Associates II Word Recall: "an examinee states as many items as
     he or she can freely recall... there is not a consistent number
     of items to evaluate across individuals", so internal consistency
     is not appropriate. It takes rStability, the field RBANS Update
     Table 3.6 introduced. Corroborated the usual way — Table 3.1's
     .76 equals this database's stored rCorrected in BOTH batteries
     and neither raw r.

     THE RECOGNITION MEMORY MEASURES ARE DELIBERATELY ABSENT, and must
     stay absent. The manual reports their reliability as a DECISION-
     CONSISTENCY percentage: percent agreement of impaired vs not
     impaired at a 10th-percentile cut, used because those tasks are
     cumulative percentages with skewed distributions whose retest
     correlations are attenuated by range restriction. A percent
     agreement is not a correlation and cannot enter
     SEM = SD sqrt(1 - rxx). They are in neither Table 3.1 nor normDB;
     section 30 fails if one ever appears.

     rInternalAgeMax is the BATTERY's ceiling — 69 for Adult, 90 for
     Older Adult — so an age past it takes the published average
     rather than silently re-reading the topmost band.
     --------------------------------------------------------------- */
  "WMS-IV Indices · Ages 16-69": {
    "Auditory Memory Index": { m1:100.1, sd1:14.1, m2:111.6, sd2:14.4, r:0.81, rCorrected:0.83, n:168, separateBattery:true, rInternal:0.95, rInternalAgeMax:69, rInternalByAge:{16:0.94, 18:0.95, 20:0.95, 25:0.97, 30:0.94, 35:0.96, 45:0.95, 55:0.95, 65:0.94} },
    "Visual Memory Index": { m1:100, sd1:14.8, m2:112.1, sd2:16.6, r:0.8, rCorrected:0.81, n:144, separateBattery:true, rInternal:0.96, rInternalAgeMax:69, rInternalByAge:{16:0.96, 18:0.97, 20:0.96, 25:0.96, 30:0.96, 35:0.96, 45:0.95, 55:0.96, 65:0.95} },
    "Visual Working Memory Index": { m1:99.5, sd1:14.4, m2:103.8, sd2:15.6, r:0.82, rCorrected:0.83, n:171, separateBattery:true, rInternal:0.93, rInternalAgeMax:69, rInternalByAge:{16:0.89, 18:0.92, 20:0.94, 25:0.94, 30:0.94, 35:0.91, 45:0.93, 55:0.93, 65:0.92} },
    "Immediate Memory Index": { m1:99.9, sd1:14.9, m2:112.3, sd2:15.6, r:0.81, rCorrected:0.81, n:154, separateBattery:true, rInternal:0.95, rInternalAgeMax:69, rInternalByAge:{16:0.93, 18:0.95, 20:0.95, 25:0.96, 30:0.95, 35:0.95, 45:0.94, 55:0.94, 65:0.93} },
    "Delayed Memory Index": { m1:100.4, sd1:13.9, m2:114.1, sd2:15, r:0.79, rCorrected:0.82, n:150, separateBattery:true, rInternal:0.94, rInternalAgeMax:69, rInternalByAge:{16:0.93, 18:0.95, 20:0.95, 25:0.96, 30:0.92, 35:0.94, 45:0.94, 55:0.94, 65:0.92} }
  },
  "WMS-IV Indices · Ages 65-90": {
    "Auditory Memory Index": { m1:101.5, sd1:12.9, m2:112.1, sd2:14.5, r:0.82, rCorrected:0.87, n:69, separateBattery:true, rInternal:0.95, rInternalAgeMax:90, rInternalByAge:{65:0.95, 70:0.94, 75:0.95, 80:0.95, 85:0.95} },
    "Visual Memory Index": { m1:101.6, sd1:14.6, m2:112.6, sd2:17, r:0.79, rCorrected:0.8, n:70, separateBattery:true, rInternal:0.97, rInternalAgeMax:90, rInternalByAge:{65:0.96, 70:0.97, 75:0.96, 80:0.97, 85:0.97} },
    "Immediate Memory Index": { m1:101.5, sd1:13.9, m2:113.9, sd2:14.2, r:0.84, rCorrected:0.86, n:69, separateBattery:true, rInternal:0.95, rInternalAgeMax:90, rInternalByAge:{65:0.94, 70:0.95, 75:0.94, 80:0.94, 85:0.96} },
    "Delayed Memory Index": { m1:101.7, sd1:13.1, m2:112.7, sd2:14.4, r:0.8, rCorrected:0.85, n:71, separateBattery:true, rInternal:0.92, rInternalAgeMax:90, rInternalByAge:{65:0.92, 70:0.91, 75:0.91, 80:0.93, 85:0.92} }
  },
  "WMS-IV Subtests · Ages 16-69": {
    "Logical Memory I": { m1:10.3, sd1:2.9, m2:12.2, sd2:2.6, r:0.72, rCorrected:0.74, n:173, separateBattery:true, rInternal:0.82, rInternalAgeMax:69, rInternalByAge:{16:0.8, 18:0.84, 20:0.79, 25:0.87, 30:0.82, 35:0.86, 45:0.77, 55:0.81, 65:0.79} },
    "Logical Memory II": { m1:10.3, sd1:2.8, m2:12.6, sd2:2.9, r:0.67, rCorrected:0.71, n:172, separateBattery:true, rInternal:0.85, rInternalAgeMax:69, rInternalByAge:{16:0.8, 18:0.87, 20:0.85, 25:0.9, 30:0.81, 35:0.9, 45:0.87, 55:0.84, 65:0.8} },
    "Verbal Paired Associates I": { m1:9.8, sd1:3.1, m2:12.1, sd2:3.4, r:0.76, rCorrected:0.74, n:171, separateBattery:true, rInternal:0.94, rInternalAgeMax:69, rInternalByAge:{16:0.93, 18:0.93, 20:0.94, 25:0.95, 30:0.94, 35:0.94, 45:0.94, 55:0.94, 65:0.93} },
    "Verbal Paired Associates II": { m1:9.8, sd1:3, m2:10.8, sd2:2.7, r:0.76, rCorrected:0.76, n:170, separateBattery:true, rInternal:0.85, rInternalAgeMax:69, rInternalByAge:{16:0.87, 18:0.84, 20:0.84, 25:0.89, 30:0.82, 35:0.84, 45:0.85, 55:0.85, 65:0.82} },
    "Verbal Paired Associates II - Word Recall": { m1:10, sd1:2.9, m2:10.9, sd2:3.2, r:0.74, rCorrected:0.76, n:170, separateBattery:true, rStability:0.76, rStabilityAgeMax:69, rStabilityByAge:{16:0.76, 18:0.76, 20:0.76, 25:0.76, 30:0.76, 35:0.76, 45:0.76, 55:0.76, 65:0.76} },
    "Designs I": { m1:10, sd1:2.9, m2:11.1, sd2:3.4, r:0.73, rCorrected:0.75, n:157, separateBattery:true, rInternal:0.85, rInternalAgeMax:69, rInternalByAge:{16:0.84, 18:0.9, 20:0.87, 25:0.83, 30:0.88, 35:0.85, 45:0.83, 55:0.82, 65:0.82} },
    "Designs I - Content": { m1:10.1, sd1:3, m2:11.2, sd2:3.3, r:0.64, rCorrected:0.64, n:157, separateBattery:true, rInternal:0.77, rInternalAgeMax:69, rInternalByAge:{16:0.71, 18:0.77, 20:0.75, 25:0.88, 30:0.77, 35:0.76, 45:0.75, 55:0.81, 65:0.66} },
    "Designs I - Spatial": { m1:10.2, sd1:2.9, m2:10.9, sd2:3.1, r:0.56, rCorrected:0.59, n:157, separateBattery:true, rInternal:0.76, rInternalAgeMax:69, rInternalByAge:{16:0.76, 18:0.78, 20:0.7, 25:0.78, 30:0.83, 35:0.73, 45:0.76, 55:0.71, 65:0.77} },
    "Designs II": { m1:10.2, sd1:2.7, m2:11.9, sd2:3.2, r:0.72, rCorrected:0.77, n:151, separateBattery:true, rInternal:0.85, rInternalAgeMax:69, rInternalByAge:{16:0.88, 18:0.9, 20:0.87, 25:0.84, 30:0.87, 35:0.81, 45:0.8, 55:0.83, 65:0.82} },
    "Designs II - Content": { m1:10.3, sd1:3, m2:11.5, sd2:3.4, r:0.64, rCorrected:0.64, n:151, separateBattery:true, rInternal:0.77, rInternalAgeMax:69, rInternalByAge:{16:0.81, 18:0.76, 20:0.79, 25:0.84, 30:0.79, 35:0.77, 45:0.7, 55:0.75, 65:0.71} },
    "Designs II - Spatial": { m1:10.3, sd1:2.7, m2:11.6, sd2:2.6, r:0.5, rCorrected:0.6, n:151, separateBattery:true, rInternal:0.74, rInternalAgeMax:69, rInternalByAge:{16:0.74, 18:0.82, 20:0.69, 25:0.73, 30:0.75, 35:0.7, 45:0.68, 55:0.81, 65:0.67} },
    "Visual Reproduction I": { m1:10, sd1:2.8, m2:11.9, sd2:2.8, r:0.62, rCorrected:0.67, n:170, separateBattery:true, rInternal:0.93, rInternalAgeMax:69, rInternalByAge:{16:0.92, 18:0.94, 20:0.88, 25:0.95, 30:0.93, 35:0.92, 45:0.93, 55:0.96, 65:0.92} },
    "Visual Reproduction II": { m1:10.1, sd1:2.8, m2:12.9, sd2:3, r:0.59, rCorrected:0.64, n:169, separateBattery:true, rInternal:0.97, rInternalAgeMax:69, rInternalByAge:{16:0.97, 18:0.98, 20:0.97, 25:0.97, 30:0.96, 35:0.97, 45:0.97, 55:0.98, 65:0.98} },
    "Spatial Addition": { m1:9.9, sd1:2.8, m2:10.7, sd2:3, r:0.74, rCorrected:0.77, n:172, separateBattery:true, rInternal:0.91, rInternalAgeMax:69, rInternalByAge:{16:0.89, 18:0.89, 20:0.92, 25:0.92, 30:0.89, 35:0.91, 45:0.92, 55:0.93, 65:0.9} },
    "Symbol Span": { m1:10, sd1:3, m2:10.6, sd2:3.1, r:0.72, rCorrected:0.72, n:172, separateBattery:true, rInternal:0.88, rInternalAgeMax:69, rInternalByAge:{16:0.81, 18:0.88, 20:0.89, 25:0.9, 30:0.92, 35:0.86, 45:0.89, 55:0.88, 65:0.85} }
  },
  "WMS-IV Subtests · Ages 65-90": {
    "Logical Memory I": { m1:10, sd1:2.9, m2:12, sd2:3.3, r:0.77, rCorrected:0.79, n:69, separateBattery:true, rInternal:0.86, rInternalAgeMax:90, rInternalByAge:{65:0.83, 70:0.88, 75:0.87, 80:0.84, 85:0.88} },
    "Logical Memory II": { m1:10, sd1:2.7, m2:12.1, sd2:2.8, r:0.71, rCorrected:0.77, n:71, separateBattery:true, rInternal:0.87, rInternalAgeMax:90, rInternalByAge:{65:0.85, 70:0.87, 75:0.88, 80:0.87, 85:0.89} },
    "Verbal Paired Associates I": { m1:10.4, sd1:2.8, m2:12.1, sd2:2.9, r:0.76, rCorrected:0.79, n:71, separateBattery:true, rInternal:0.93, rInternalAgeMax:90, rInternalByAge:{65:0.94, 70:0.92, 75:0.92, 80:0.92, 85:0.93} },
    "Verbal Paired Associates II": { m1:10.4, sd1:2.7, m2:11.5, sd2:2.7, r:0.77, rCorrected:0.81, n:71, separateBattery:true, rInternal:0.74, rInternalAgeMax:90, rInternalByAge:{65:0.82, 70:0.71, 75:0.71, 80:0.78, 85:0.68} },
    "Verbal Paired Associates II - Word Recall": { m1:10.5, sd1:2.8, m2:11.7, sd2:2.7, r:0.72, rCorrected:0.76, n:71, separateBattery:true, rStability:0.76, rStabilityAgeMax:90, rStabilityByAge:{65:0.76, 70:0.76, 75:0.76, 80:0.76, 85:0.76} },
    "Visual Reproduction I": { m1:10.2, sd1:3.1, m2:12, sd2:3.1, r:0.79, rCorrected:0.78, n:71, separateBattery:true, rInternal:0.93, rInternalAgeMax:90, rInternalByAge:{65:0.92, 70:0.94, 75:0.92, 80:0.92, 85:0.93} },
    "Visual Reproduction II": { m1:10.5, sd1:2.8, m2:12.3, sd2:3.1, r:0.64, rCorrected:0.69, n:71, separateBattery:true, rInternal:0.96, rInternalAgeMax:90, rInternalByAge:{65:0.95, 70:0.96, 75:0.96, 80:0.97, 85:0.96} },
    "Symbol Span": { m1:10.1, sd1:2.8, m2:10.7, sd2:3, r:0.69, rCorrected:0.73, n:69, separateBattery:true, rInternal:0.84, rInternalAgeMax:90, rInternalByAge:{65:0.88, 70:0.87, 75:0.81, 80:0.84, 85:0.76} }
  }
};

/* ============================================================================
   PERFORMANCE VALIDITY (PVT) REFERENCE DATA

   Sources (primary; the on-screen references list carries full citations):
   - RBANS Effort Index:  Silverberg, Wertheimer & Fichtenberg (2007), TCN 21(5), 841-854.
   - RBANS Effort Scale:  Novitski, Steele, Karantzoulis & Randolph (2012), ACN 27(2), 190-195.
   - Reliable Digit Span: Greiffenstein, Baker & Gola (1994), Psych. Assessment 6(3), 218-224;
                          Meyers & Volbrecht (1998); Schroeder et al. (2012).
   - TOMM:                Martin, Schroeder, Olsen, Maloy, Boettcher, Ernst & Okut (2020),
                          TCN 34(1), 88-119 (meta-analysis); Tombaugh (1996).
   - Aggregation:         Larrabee (2014), ACN 29(4), 364-373.

   These pages take RAW scores directly and never read normDB. Cut-offs and
   weights are stored verbatim from the cited papers; every table below is
   pinned in tools/check.js §38. PVT results are reported as cut-off
   comparisons, never as verdicts — no output may label an examinee
   "malingering" or the profile "invalid" from a single indicator.
   ========================================================================= */

/* Silverberg et al. (2007) Table 2 — raw-to-weighted conversion for the two
   Effort Index components. Bands are [min, max] inclusive on the raw score;
   note the deliberate gaps (Digit Span never yields a weight of 1 or 4).
   EI = weighted Digit Span + weighted List Recognition, range 0-12. */
const PVT_EI_WEIGHTS = {
  digitSpan: [            // RBANS Digit Span raw, 0-16
    { min: 8, max: 16, w: 0 },
    { min: 7, max: 7,  w: 2 },
    { min: 6, max: 6,  w: 3 },
    { min: 5, max: 5,  w: 5 },
    { min: 0, max: 4,  w: 6 }
  ],
  listRecognition: [      // RBANS List Recognition raw, 0-20
    { min: 18, max: 20, w: 0 },
    { min: 17, max: 17, w: 1 },
    { min: 15, max: 16, w: 2 },
    { min: 13, max: 14, w: 3 },
    { min: 11, max: 12, w: 4 },
    { min: 10, max: 10, w: 5 },
    { min: 0,  max: 9,  w: 6 }
  ]
};
/* EI > 3: standard threshold (~94% specificity in the derivation sample);
   EI > 1: the authors' optimal screening cut-off for post-acute mild TBI. */
const PVT_EI_CUTOFFS = { standard: 3, sensitive: 1 };

/* Novitski et al. (2012). ES = (List Recognition - [List Recall + Story
   Recall + Figure Recall]) + Digit Span, all raw. LOWER is more suspicious.
   The gate is mandatory: computed only where Digit Span < 9 OR List
   Recognition < 19 OR their sum < 28 — in intact examinees free recall
   normally far exceeds ceiling-limited recognition, so an ungated ES
   over-flags. Cut-off: ES < 12 once the gate is met. */
const PVT_ES = {
  gate: { digitSpanBelow: 9, listRecognitionBelow: 19, combinedBelow: 28 },
  cutoff: 12
};

/* Greiffenstein et al. (1994); Meyers & Volbrecht (1998). RDS = longest
   forward string passed on BOTH trials + longest backward string passed on
   BOTH trials (classic variant: Forward + Backward only, no WAIS-IV
   Sequencing). Floor rule: failing at least one trial each of 3 forward and
   2 backward is recorded as RDS = 3.
   <= 7: traditional cut-off. <= 6: conservative — Schroeder et al. (2012)
   found <= 7 had inadequate specificity across clinical groups while <= 6
   held specificity > .90; prefer <= 6 in genuinely impaired populations. */
const PVT_RDS = { floor: 3, cutoffTraditional: 7, cutoffConservative: 6 };

/* Martin et al. (2020) meta-analytic TOMM cut-offs — weighted-mean
   specificity/sensitivity for neurocognitive/psychiatric samples (the most
   clinically generalisable subgroup). Scores are correct responses out of 50;
   "cut" means scores BELOW it fail. sens/spec are the point values the
   meta-analysis uses for its predictive-power tables (Tables 16-17), from
   which the app derives PPP/NPP by Bayes at the chosen base rate.
   specRange/sensRange, where present, are the ranges quoted for the cut-off
   summary table. Trial 1 < 41 row: Denning (2012). */
const PVT_TOMM_CUTOFFS = [
  { id: 't1-41',  trial: 'trial1',    label: 'Trial 1 < 41',   cut: 41, sens: 0.66, spec: 0.93, note: 'More conservative Trial 1 cut-off; meta-analytic weighted-mean values at < 41 (cf. Denning, 2012).' },
  { id: 't1-42',  trial: 'trial1',    label: 'Trial 1 < 42',   cut: 42, sens: 0.69, spec: 0.91, note: 'Meta-analytic optimum for abbreviated/screening use.' },
  { id: 't2-45',  trial: 'trial2',    label: 'Trial 2 < 45',   cut: 45, sens: 0.45, spec: 0.97, specRange: '.96–.98', sensRange: '.45–.55', note: 'Traditional cut-off; very high specificity across settings.' },
  { id: 't2-49',  trial: 'trial2',    label: 'Trial 2 < 49',   cut: 49, sens: 0.63, spec: 0.95, specRange: '.91–.97', sensRange: '.59–.70', note: 'Liberal; only where significant impairment is not expected.' },
  { id: 'ret-45', trial: 'retention', label: 'Retention < 45', cut: 45, sens: 0.55, spec: 0.98, note: 'Traditional cut-off; very high specificity across settings.' },
  { id: 'ret-49', trial: 'retention', label: 'Retention < 49', cut: 49, sens: 0.70, spec: 0.93, note: 'Liberal; only where significant impairment is not expected.' }
];
/* Base-rate options for the PPP/NPP display, as in Martin et al. Tables 16-17. */
const PVT_BASE_RATES = [0.10, 0.20, 0.30, 0.40, 0.50];

/* Published classification accuracy for the selectable cut-offs, shown on
   screen and in the APA export. Strings, not numbers, because the sources
   publish RANGES (across samples or statistical methods) and the app renders
   what the source prints:
   - EI (Silverberg et al., 2007): at > 3, specificity .94 in the mixed
     clinical derivation sample rising to 1.00 in mild TBI and controls
     (Tables 1 and 3); sensitivity .46-.71 across the clinical, simulated-
     naive and simulated-coached malingering groups (Table 3). At > 1,
     specificity .75 (derivation) to .96 (controls); sensitivity .67-.92.
   - RDS (Schroeder et al., 2012): global rates by weighted average and
     Bayesian method — at <= 6, sensitivity .30/.35 and specificity .96/.97;
     at <= 7, sensitivity .48/.58 and specificity .82/.85.
   - ES (Novitski et al., 2012): NO published sensitivity/specificity pair
     exists; discrimination is published as ROC AUC = .91 against amnestic
     patients (vs .61 for the EI in the same comparison). Do not invent a
     pair — the APA table prints dashes and the note explains.
   TOMM accuracy lives on PVT_TOMM_CUTOFFS above. Pinned by check.js §38. */
const PVT_EI_ACCURACY = {
  standard:  { sens: '.46–.71', spec: '.94–1.00' },
  sensitive: { sens: '.67–.92', spec: '.75–.96' }
};
const PVT_RDS_ACCURACY = {
  conservative: { sens: '.30–.35', spec: '.96–.97' },
  traditional:  { sens: '.48–.58', spec: '.82–.85' }
};
const PVT_ES_ACCURACY = { auc: '.91' };

/* WAIS Digit Span Age-Corrected Scaled Score (Iverson & Tulsky, 2003;
   Axelrod, Fichtenberg, Millis & Wertheimer, 2006). BOTH SOURCES ARE
   WAIS-III: Iverson & Tulsky tabulate base rates in the WAIS-III
   standardization sample (N = 2,450) and six clinical groups; Axelrod et al.
   give diagnostic efficiency per cut-off (Table 3) against probable
   malingerers vs moderate/severe TBI, cross-validated on non-litigating
   mild TBI. The page says so and tells the clinician to document the
   edition administered.
   - ACSS <= 5 (conservative, default): standardization base rate 3.8%,
     combined clinical 3.4% (Iverson & Tulsky Tables 1/4 — their suspicion
     guideline is "scaled score of 5, 4, or less"); Axelrod Table 3 at
     <= 5: sensitivity .36, specificity .97.
   - ACSS <= 7 (Axelrod's optimum): sensitivity .75, specificity .69 (TBI),
     .77 on the mild-TBI cross-validation.
   - Vocabulary - Digit Span difference >= 5: a second Iverson & Tulsky
     suspicion index, interpreted against BASE RATES, not sens/spec —
     7.1% of the standardization sample and 2.8% of the combined clinical
     sample score >= 5 (Tables 3/6).
   DS ACSS and RDS come from the SAME subtest: the summary counts them as
   one indicator (the digit-span group), exactly as EI/ES share theirs. */
const PVT_DS = { cutoffConservative: 5, cutoffSensitive: 7, vocabDiffCutoff: 5 };
const PVT_DS_ACCURACY = {
  conservative: { sens: '.36', spec: '.97' },
  sensitive:    { sens: '.75', spec: '.69–.77' }
};
const PVT_DS_VOCABDIFF_BASERATES = { standardization: '7.1%', clinical: '2.8%' };

/* Larrabee (2014), combined clinical sample, 6 PVTs + 1 SVT — classification
   accuracy by number of failures. Percentages as published. */
const PVT_AGGREGATION = [
  { threshold: '≥ 2 of 7 failures', spec: 88.9, sens: 97.6, correct: 92.6 },
  { threshold: '≥ 3 of 7 failures', spec: 96.3, sens: 87.8, correct: 92.6 },
  { threshold: '≥ 4 of 7 failures', spec: 100,  sens: 63.4, correct: 84.2 }
];
