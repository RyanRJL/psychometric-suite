/* ============================================================
   DATA · constants, lookup tables, coefficients and tooltip strings
   Extracted from Psychometric Calculators · Clinical Suite
   ============================================================ */

let normDB = {}; // populated below + custom

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

// OPIE-4 R and SEE (prorated FSIQ - Schoenberg et al., 2011 Table 5.16, post-Age step)
// Sex term included in equations; full-sample SEE used as a conservative CI estimate.
const OPIE_STATS = {
  noAge:   { MR:{r:0.57,see:11.98}, VC:{r:0.70,see:10.52}, MR_VC:{r:0.77,see:9.34} },
  withAge: { MR:{r:0.66,see:10.96}, VC:{r:0.71,see:10.27}, MR_VC:{r:0.80,see:8.80} }
};

// OPIE-4 prorated FSIQ regression coefficients (Table eA5.8)
// UK adaptation: US education, region, and ethnicity terms omitted.
// Sex coding: Female = 0, Male = 1 (per OPIE-4 source)
const OPIE_PRORATED_FSIQ = {
  VC_MR: { intercept:65.77827122, vc:0.646258435, mr:1.182068623, age:-0.197692558, age3:0.0000373292, sex:1.955504838, r:0.80, see:8.80 },
  VC:    { intercept:86.63733022, vc:0.825479066,                 age:-0.355783733, age3:0.0000373292, sex:2.795447219, r:0.71, see:10.27 },
  MR:    { intercept:62.02281403,                mr:1.719384768,                    age3:0.0000275723, sex:1.509479923, r:0.66, see:10.96 }
};

// OPIE-4 prorated GAI regression coefficients (Table eA5.8)
// MR-only branch omitted because the source column is partially truncated in the available source.
// R/SEE from Table 5.19 (Prorated GAI), post-Age step.
const OPIE_PRORATED_GAI = {
  VC_MR: { intercept:60.14203956, vc:0.763136717, mr:1.127062322, age:-0.246247784, age3:0.0000416209, sex:4.708926488, r:0.79, see:9.25 },
  VC:    { intercept:79.65445374, vc:0.921039566,                 age:-0.378906405, age3:0.0000399793, sex:5.001045997, r:0.75, see:9.86 }
};

const OCC_CODE = { 'Professional':1, 'Intermediate':2, 'Skilled':3, 'Semi Skilled':4, 'Unskilled':5 };

const PRE_MODEL_TOOLTIPS = {
  topfRaw: 'Assumes the ToPF word-reading raw score is a resistant estimate of premorbid ability. Input required: ToPF raw score only. The FSIQ estimate is returned from the ToPF raw-score look-up table; CI uses the model SEE.',
  topfDemo: 'Assumes premorbid FSIQ is best estimated by combining ToPF performance with demographic predictors. Inputs required: ToPF raw score, years of education, and sex. Uses the cubic ToPF + education + sex regression equation; CI uses the model SEE.',
  crawfordAllan: 'Demographic-only estimate from Crawford & Allan (2001), intended for UK-normed demographic prediction. Inputs required: occupation class, years of education, and age. Does not use the ToPF score.',
  opieDefault: 'OPIE-4 estimate of prorated FSIQ. Inputs required: age plus Vocabulary raw and/or Matrix Reasoning raw (sex optional). Adapted for UK use: US education, region and ethnicity terms are omitted, so US-specific demographic adjustments are not applied. The equation automatically switches to the matching single- or two-subtest model. CI uses the branch-specific SEE.',
  opieVCMR: 'OPIE-4 two-subtest estimate of prorated FSIQ (Schoenberg et al., 2011; Table eA5.8). Inputs required: Vocabulary raw, Matrix Reasoning raw, and age (sex optional). Predicts a prorated FSIQ that excludes Vocabulary and Matrix Reasoning, removing part-whole correlation inflation present in the standard-FSIQ equations. UK adaptation: US education, region and ethnicity terms omitted.',
  opieVC: 'OPIE-4 single-subtest estimate of prorated FSIQ (Vocabulary branch). Inputs required: Vocabulary raw and age (sex optional). UK adaptation: US education, region and ethnicity terms omitted; CI uses the Vocabulary-branch SEE.',
  opieMR: 'OPIE-4 single-subtest estimate of prorated FSIQ (Matrix Reasoning branch). Inputs required: Matrix Reasoning raw and age (sex optional). The published equation uses only Age³ (no linear age term) for this branch. UK adaptation: US education, region and ethnicity terms omitted; CI uses the Matrix-branch SEE.',
  predictWais: 'ToPF-predicted WAIS-IV index model. Inputs required: ToPF raw score, years of education, and sex. Difference is Achieved − Predicted; base rates are shown for negative discrepancies only.',
  predictWms: 'ToPF-predicted WMS-IV index model. Inputs required: ToPF raw score and age. Difference is Achieved − Predicted; base rates are shown for negative discrepancies only.',
  opiePredFSIQ_VCMR: 'OPIE-4-predicted prorated FSIQ (two-subtest). Compare against the patient\'s prorated FSIQ calculated excluding Vocabulary and Matrix Reasoning per WAIS-IV manual.',
  opiePredFSIQ_VC: 'OPIE-4-predicted prorated FSIQ (Vocabulary only). Compare against the patient\'s prorated FSIQ calculated excluding Vocabulary per WAIS-IV manual.',
  opiePredFSIQ_MR: 'OPIE-4-predicted prorated FSIQ (Matrix Reasoning only). Compare against the patient\'s prorated FSIQ calculated excluding Matrix Reasoning per WAIS-IV manual.',
  opiePredGAI_VCMR: 'OPIE-4-predicted prorated GAI (two-subtest). Compare against the patient\'s prorated GAI calculated excluding Vocabulary and Matrix Reasoning per WAIS-IV manual.',
  opiePredGAI_VC: 'OPIE-4-predicted prorated GAI (Vocabulary only). Compare against the patient\'s prorated GAI calculated excluding Vocabulary per WAIS-IV manual.'
};

// Base rates: discrepancy (negative integer) → proportion of standardisation sample
// WAIS-IV indices (FSIQ, VCI, PRI, WMI, PSI) and WMS-IV indices (IMI, DMI, VWMI)
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
  '-33':{                     PRI:.0038,WMI:.0009,PSI:.0031,IMI:.003, DMI:.0039,VWMI:.0033},
  '-34':{                     PRI:.003, WMI:.0007,PSI:.0024,IMI:.0024,DMI:.0031,VWMI:.0026},
  '-35':{                     PRI:.0023,WMI:.0005,PSI:.0018,IMI:.0018,DMI:.0024,VWMI:.002},
  '-36':{                     PRI:.0018,WMI:.0003,PSI:.0014,IMI:.0014,DMI:.0019,VWMI:.0015},
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
  '-34':{FSIQ_VC:0.0049,GAI_VC:0.002},
  '-33':{FSIQ_VC:0.0049,GAI_VC:0.002},
  '-32':{FSIQ_VC:0.0049,GAI_VC:0.002},
  '-31':{FSIQ_VC:0.0049,GAI_VC:0.0029},
  '-30':{FSIQ_VC:0.0049,GAI_VC:0.0029},
  '-29':{FSIQ_VC:0.0049,GAI_VC:0.0029},
  '-28':{FSIQ_VC:0.0049,GAI_VC:0.0029},
  '-27':{FSIQ_MR:0.0029,FSIQ_VC:0.0049,GAI_VC:0.0039,GAI_VC_MR:0.001},
  '-26':{FSIQ_MR:0.0039,FSIQ_VC:0.0069,GAI_VC:0.0059,GAI_VC_MR:0.002},
  '-25':{FSIQ_MR:0.0039,FSIQ_VC:0.0079,FSIQ_VC_MR:0.001,GAI_VC:0.0059,GAI_VC_MR:0.002},
  '-24':{FSIQ_MR:0.0049,FSIQ_VC:0.0098,FSIQ_VC_MR:0.002,GAI_VC:0.0088,GAI_VC_MR:0.002},
  '-23':{FSIQ_MR:0.0098,FSIQ_VC:0.0108,FSIQ_VC_MR:0.002,GAI_VC:0.0088,GAI_VC_MR:0.0039},
  '-22':{FSIQ_MR:0.0118,FSIQ_VC:0.0108,FSIQ_VC_MR:0.0049,GAI_VC:0.0108,GAI_VC_MR:0.0059},
  '-21':{FSIQ_MR:0.0157,FSIQ_VC:0.0177,FSIQ_VC_MR:0.0059,GAI_VC:0.0147,GAI_VC_MR:0.0069},
  '-20':{FSIQ_MR:0.0216,FSIQ_VC:0.0216,FSIQ_VC_MR:0.0098,GAI_VC:0.0196,GAI_VC_MR:0.0098},
  '-19':{FSIQ_MR:0.0265,FSIQ_VC:0.0295,FSIQ_VC_MR:0.0138,GAI_VC:0.0236,GAI_VC_MR:0.0167},
  '-18':{FSIQ_MR:0.0432,FSIQ_VC:0.0373,FSIQ_VC_MR:0.0187,GAI_VC:0.0314,GAI_VC_MR:0.0216},
  '-17':{FSIQ_MR:0.0501,FSIQ_VC:0.0462,FSIQ_VC_MR:0.0295,GAI_VC:0.0354,GAI_VC_MR:0.0314},
  '-16':{FSIQ_MR:0.0619,FSIQ_VC:0.057,FSIQ_VC_MR:0.0354,GAI_VC:0.0452,GAI_VC_MR:0.0393},
  '-15':{FSIQ_MR:0.0678,FSIQ_VC:0.0668,FSIQ_VC_MR:0.0432,GAI_VC:0.053,GAI_VC_MR:0.0472},
  '-14':{FSIQ_MR:0.0806,FSIQ_VC:0.0796,FSIQ_VC_MR:0.0511,GAI_VC:0.0629,GAI_VC_MR:0.058},
  '-13':{FSIQ_MR:0.0992,FSIQ_VC:0.0972,FSIQ_VC_MR:0.0639,GAI_VC:0.0756,GAI_VC_MR:0.0697},
  '-12':{FSIQ_MR:0.1198,FSIQ_VC:0.113,FSIQ_VC_MR:0.0815,GAI_VC:0.0914,GAI_VC_MR:0.0855},
  '-11':{FSIQ_MR:0.1395,FSIQ_VC:0.1277,FSIQ_VC_MR:0.0982,GAI_VC:0.1081,GAI_VC_MR:0.1051},
  '-10':{FSIQ_MR:0.1739,FSIQ_VC:0.1483,FSIQ_VC_MR:0.1306,GAI_VC:0.1297,GAI_VC_MR:0.1306},
  '-9':{FSIQ_MR:0.2004,FSIQ_VC:0.1709,FSIQ_VC_MR:0.1591,GAI_VC:0.1582,GAI_VC_MR:0.1591},
  '-8':{FSIQ_MR:0.2358,FSIQ_VC:0.2073,FSIQ_VC_MR:0.1906,GAI_VC:0.1916,GAI_VC_MR:0.1945},
  '-7':{FSIQ_MR:0.2692,FSIQ_VC:0.2387,FSIQ_VC_MR:0.223,GAI_VC:0.2269,GAI_VC_MR:0.2279},
  '-6':{FSIQ_MR:0.3016,FSIQ_VC:0.2741,FSIQ_VC_MR:0.2495,GAI_VC:0.2593,GAI_VC_MR:0.2603},
  '-5':{FSIQ_MR:0.336,FSIQ_VC:0.3104,FSIQ_VC_MR:0.3055,GAI_VC:0.3006,GAI_VC_MR:0.3075},
  '-4':{FSIQ_MR:0.3644,FSIQ_VC:0.3448,FSIQ_VC_MR:0.3438,GAI_VC:0.3438,GAI_VC_MR:0.3566},
  '-3':{FSIQ_MR:0.3998,FSIQ_VC:0.3772,FSIQ_VC_MR:0.3811,GAI_VC:0.389,GAI_VC_MR:0.3949},
  '-2':{FSIQ_MR:0.446,FSIQ_VC:0.4293,FSIQ_VC_MR:0.4322,GAI_VC:0.4381,GAI_VC_MR:0.4381},
  '-1':{FSIQ_MR:0.4971,FSIQ_VC:0.4725,FSIQ_VC_MR:0.4882,GAI_VC:0.4813,GAI_VC_MR:0.4882},
  '+1':{FSIQ_MR:0.4725,FSIQ_VC:0.4872,FSIQ_VC_MR:0.4686,GAI_VC:0.4656,GAI_VC_MR:0.4617},
  '+2':{FSIQ_MR:0.4303,FSIQ_VC:0.4489,FSIQ_VC_MR:0.4273,GAI_VC:0.4106,GAI_VC_MR:0.4194},
  '+3':{FSIQ_MR:0.388,FSIQ_VC:0.3969,FSIQ_VC_MR:0.3821,GAI_VC:0.3752,GAI_VC_MR:0.3802},
  '+4':{FSIQ_MR:0.3507,FSIQ_VC:0.3418,FSIQ_VC_MR:0.3399,GAI_VC:0.336,GAI_VC_MR:0.335},
  '+5':{FSIQ_MR:0.3114,FSIQ_VC:0.3075,FSIQ_VC_MR:0.3065,GAI_VC:0.2849,GAI_VC_MR:0.2917},
  '+6':{FSIQ_MR:0.2809,FSIQ_VC:0.2721,FSIQ_VC_MR:0.2603,GAI_VC:0.2574,GAI_VC_MR:0.2593},
  '+7':{FSIQ_MR:0.2485,FSIQ_VC:0.2446,FSIQ_VC_MR:0.2151,GAI_VC:0.2308,GAI_VC_MR:0.2191},
  '+8':{FSIQ_MR:0.222,FSIQ_VC:0.2033,FSIQ_VC_MR:0.1857,GAI_VC:0.2033,GAI_VC_MR:0.1906},
  '+9':{FSIQ_MR:0.1916,FSIQ_VC:0.1827,FSIQ_VC_MR:0.1532,GAI_VC:0.1768,GAI_VC_MR:0.1532},
  '+10':{FSIQ_MR:0.168,FSIQ_VC:0.1552,FSIQ_VC_MR:0.1257,GAI_VC:0.1483,GAI_VC_MR:0.1267},
  '+11':{FSIQ_MR:0.1424,FSIQ_VC:0.1306,FSIQ_VC_MR:0.0992,GAI_VC:0.1277,GAI_VC_MR:0.1022},
  '+12':{FSIQ_MR:0.1257,FSIQ_VC:0.1081,FSIQ_VC_MR:0.0776,GAI_VC:0.1022,GAI_VC_MR:0.0796},
  '+13':{FSIQ_MR:0.111,FSIQ_VC:0.0923,FSIQ_VC_MR:0.0629,GAI_VC:0.0835,GAI_VC_MR:0.0648},
  '+14':{FSIQ_MR:0.0953,FSIQ_VC:0.0707,FSIQ_VC_MR:0.0501,GAI_VC:0.0697,GAI_VC_MR:0.0521},
  '+15':{FSIQ_MR:0.0806,FSIQ_VC:0.0599,FSIQ_VC_MR:0.0373,GAI_VC:0.0589,GAI_VC_MR:0.0452},
  '+16':{FSIQ_MR:0.0678,FSIQ_VC:0.0511,FSIQ_VC_MR:0.0334,GAI_VC:0.0472,GAI_VC_MR:0.0363},
  '+17':{FSIQ_MR:0.054,FSIQ_VC:0.0413,FSIQ_VC_MR:0.0265,GAI_VC:0.0363,GAI_VC_MR:0.0314},
  '+18':{FSIQ_MR:0.0452,FSIQ_VC:0.0373,FSIQ_VC_MR:0.0236,GAI_VC:0.0255,GAI_VC_MR:0.0275},
  '+19':{FSIQ_MR:0.0403,FSIQ_VC:0.0265,FSIQ_VC_MR:0.0196,GAI_VC:0.0187,GAI_VC_MR:0.0196},
  '+20':{FSIQ_MR:0.0344,FSIQ_VC:0.0196,FSIQ_VC_MR:0.0147,GAI_VC:0.0177,GAI_VC_MR:0.0147},
  '+21':{FSIQ_MR:0.0265,FSIQ_VC:0.0147,FSIQ_VC_MR:0.0088,GAI_VC:0.0138,GAI_VC_MR:0.0128},
  '+22':{FSIQ_MR:0.0216,FSIQ_VC:0.0108,FSIQ_VC_MR:0.0059,GAI_VC:0.0079,GAI_VC_MR:0.0088},
  '+23':{FSIQ_MR:0.0177,FSIQ_VC:0.0108,FSIQ_VC_MR:0.0039,GAI_VC:0.0079,GAI_VC_MR:0.0088},
  '+24':{FSIQ_MR:0.0138,FSIQ_VC:0.0108,FSIQ_VC_MR:0.0029,GAI_VC:0.0059,GAI_VC_MR:0.0088},
  '+25':{FSIQ_MR:0.0118,FSIQ_VC:0.0059,FSIQ_VC_MR:0.0029,GAI_VC:0.0039,GAI_VC_MR:0.0039},
  '+26':{FSIQ_MR:0.0079,FSIQ_VC:0.0059,FSIQ_VC_MR:0.0029,GAI_VC:0.0039,GAI_VC_MR:0.0039},
  '+27':{FSIQ_MR:0.0069,FSIQ_VC:0.0049,FSIQ_VC_MR:0.0029,GAI_VC:0.0029,GAI_VC_MR:0.0039},
  '+28':{FSIQ_MR:0.0049,FSIQ_VC:0.0039,FSIQ_VC_MR:0.002,GAI_VC:0.0029,GAI_VC_MR:0.0039},
  '+29':{FSIQ_MR:0.0039,FSIQ_VC:0.002,FSIQ_VC_MR:0.002,GAI_VC:0.0029,GAI_VC_MR:0.0039},
  '+30':{FSIQ_MR:0.0029,FSIQ_VC:0.001,FSIQ_VC_MR:0.002,GAI_VC:0.0029,GAI_VC_MR:0.0029},
  '+31':{FSIQ_MR:0.002,FSIQ_VC_MR:0.001,GAI_VC:0.002,GAI_VC_MR:0.0029},
  '+32':{FSIQ_MR:0.002,GAI_VC_MR:0.0029},
  '+33':{FSIQ_MR:0.001,GAI_VC_MR:0.002},
  '+34':{GAI_VC_MR:0.002},
  '+35':{GAI_VC_MR:0.002},
  '+36':{GAI_VC_MR:0.002},
  '+37':{GAI_VC_MR:0.002},
  '+38':{GAI_VC_MR:0.001}
};

/* ============================================================
   NORMATIVE DATABASE
   Test–retest parameters from published manuals/literature
   m1/sd1 = Test, m2/sd2 = Retest, r = corrected reliability
   ============================================================ */
﻿normDB = {
  "CVLT-3 Indices · Ages 16-44": {
    "T1-5 Correct": { m1:97.9, sd1:15.8, m2:98.5, sd2:14, r:0.77, rCorrected:0.75, n:100 },
    "Delayed Recall Correct": { m1:97.9, sd1:15.3, m2:98.4, sd2:15.7, r:0.78, rCorrected:0.77, n:100 },
    "Total Recall Correct": { m1:97.7, sd1:16, m2:98.7, sd2:14.8, r:0.82, rCorrected:0.79, n:100 }
  },
  "CVLT-3 Indices · Ages 45-90": {
    "T1-5 Correct": { m1:102.6, sd1:14.3, m2:101.3, sd2:14.6, r:0.77, rCorrected:0.8, n:100 },
    "Delayed Recall Correct": { m1:100.9, sd1:15.4, m2:101.1, sd2:15.7, r:0.81, rCorrected:0.8, n:100 },
    "Total Recall Correct": { m1:101.7, sd1:15.1, m2:101.2, sd2:15.5, r:0.84, rCorrected:0.83, n:100 }
  },
  "CVLT-3 Trials · Ages 16-44": {
    "Trial 1": { m1:9.6, sd1:3, m2:9.7, sd2:3, r:0.49, rCorrected:0.51, n:100 },
    "Trial 2": { m1:9.8, sd1:3.1, m2:10, sd2:2.6, r:0.56, rCorrected:0.53, n:100 },
    "Trial 3": { m1:9.7, sd1:3, m2:9.7, sd2:3, r:0.69, rCorrected:0.7, n:100 },
    "Trial 4": { m1:9.8, sd1:3, m2:10, sd2:3, r:0.57, rCorrected:0.58, n:100 },
    "Trial 5": { m1:9.3, sd1:3.1, m2:9.1, sd2:2.9, r:0.53, rCorrected:0.5, n:100 },
    "List B Correct": { m1:9.5, sd1:2.7, m2:9.8, sd2:2.8, r:0.45, rCorrected:0.57, n:100 },
    "Short Delay Free Recall": { m1:9.8, sd1:3, m2:9.6, sd2:2.9, r:0.7, rCorrected:0.71, n:100 },
    "Short Delay Cued Recall": { m1:9.8, sd1:3.2, m2:9.6, sd2:3.2, r:0.61, rCorrected:0.56, n:100 },
    "Long Delay Free Recall": { m1:9.6, sd1:3.2, m2:9.5, sd2:3.3, r:0.65, rCorrected:0.61, n:100 },
    "Long Delay Cued Recall": { m1:9.6, sd1:3.2, m2:9.6, sd2:2.8, r:0.65, rCorrected:0.6, n:100 },
    "Recognition": { m1:10, sd1:3.1, m2:9.8, sd2:2.9, r:0.51, rCorrected:0.49, n:100 },
    "Recognition False Positive": { m1:9.5, sd1:2.8, m2:9.1, sd2:2.9, r:0.5, rCorrected:0.57, n:100 },
    "Recognition Discrimination": { m1:9.6, sd1:2.8, m2:9.3, sd2:2.7, r:0.61, rCorrected:0.67, n:100 },
    "Discrimination Nonparametric": { m1:9.8, sd1:2.7, m2:9.4, sd2:2.6, r:0.57, rCorrected:0.65, n:100 },
    "Total Intrusions": { m1:9.8, sd1:3.2, m2:9.7, sd2:2.7, r:0.33, rCorrected:0.25, n:100 },
    "Total Repetitions": { m1:10.2, sd1:3.3, m2:10.1, sd2:2.9, r:0.65, rCorrected:0.59, n:100 }
  },
  "CVLT-3 Trials · Ages 45-90": {
    "Trial 1": { m1:10.3, sd1:2.6, m2:10, sd2:3, r:0.42, rCorrected:0.56, n:100 },
    "Trial 2": { m1:10.5, sd1:2.9, m2:10.3, sd2:3.2, r:0.55, rCorrected:0.59, n:100 },
    "Trial 3": { m1:10.7, sd1:3, m2:10.5, sd2:2.8, r:0.62, rCorrected:0.63, n:100 },
    "Trial 4": { m1:10.5, sd1:3, m2:10.3, sd2:2.6, r:0.67, rCorrected:0.67, n:100 },
    "Trial 5": { m1:10.5, sd1:2.8, m2:10.2, sd2:2.9, r:0.67, rCorrected:0.71, n:100 },
    "List B Correct": { m1:10.1, sd1:2.9, m2:10.1, sd2:3, r:0.48, rCorrected:0.53, n:100 },
    "Short Delay Free Recall": { m1:10.3, sd1:3.2, m2:10.5, sd2:3.2, r:0.74, rCorrected:0.71, n:100 },
    "Short Delay Cued Recall": { m1:10.1, sd1:3.1, m2:10.2, sd2:3.1, r:0.7, rCorrected:0.69, n:100 },
    "Long Delay Free Recall": { m1:10.2, sd1:3.1, m2:10.2, sd2:3.1, r:0.73, rCorrected:0.71, n:100 },
    "Long Delay Cued Recall": { m1:10.1, sd1:3.2, m2:9.9, sd2:3, r:0.69, rCorrected:0.65, n:100 },
    "Recognition": { m1:10.4, sd1:3, m2:10.5, sd2:3, r:0.61, rCorrected:0.62, n:100 },
    "Recognition False Positive": { m1:10.1, sd1:3.2, m2:10.3, sd2:3.2, r:0.59, rCorrected:0.55, n:100 },
    "Recognition Discrimination": { m1:10.2, sd1:3.5, m2:10.3, sd2:3.2, r:0.65, rCorrected:0.54, n:100 },
    "Discrimination Nonparametric": { m1:10.2, sd1:3.4, m2:10.3, sd2:3.3, r:0.68, rCorrected:0.58, n:100 },
    "Total Intrusions": { m1:9.4, sd1:3, m2:9.7, sd2:3.3, r:0.56, rCorrected:0.57, n:100 },
    "Total Repetitions": { m1:9.7, sd1:3.1, m2:10, sd2:3.4, r:0.72, rCorrected:0.7, n:100 }
  },
  "CVLT-C Subtests (Raw Scores) · Age 8": {
    "List A Trials 1-5 Total": { m1:42.47, sd1:9.55, m2:48.32, sd2:9.09, r:0.73, n:35 },
    "List B Free-Recall Trial": { m1:5.26, sd1:1.71, m2:5.39, sd2:2.16, r:0.59, n:35 },
    "Short-Delay Free Recall": { m1:8.5, sd1:2.04, m2:9.91, sd2:2.37, r:0.4, n:35 },
    "Short-Delay Cued Recall": { m1:8.35, sd1:2.07, m2:9.51, sd2:2.79, r:0.75, n:35 },
    "Long-Delay Free Recall": { m1:8.94, sd1:2.3, m2:10.14, sd2:2.59, r:0.59, n:35 },
    "Long-Delay Cued Recall": { m1:8.42, sd1:2.21, m2:10.15, sd2:3.09, r:0.69, n:35 },
    "Semantic Cluster Ratio": { m1:1.44, sd1:0.44, m2:1.83, sd2:0.63, r:0.56, n:35 },
    "Perseverations": { m1:7.35, sd1:6.94, m2:7.54, sd2:6.39, r:0.9, n:35 },
    "Free-Recall Intrusions": { m1:5.11, sd1:5.55, m2:4.97, sd2:5.99, r:0.74, n:35 },
    "Cued-Recall Intrusions": { m1:2.83, sd1:4.33, m2:2.31, sd2:3.31, r:0.59, n:35 },
    "Recognition Hits": { m1:13.67, sd1:1.45, m2:13.84, sd2:1.33, r:0.38, n:35 },
    "Discriminability": { m1:92.76, sd1:5.69, m2:94.52, sd2:5.9, r:0.55, n:35 },
    "False Positives": { m1:1.82, sd1:2.38, m2:1.16, sd2:2, r:0.62, n:35 }
  },
  "CVLT-C Subtests (Raw Scores) · Age 12": {
    "List A Trials 1-5 Total": { m1:50.64, sd1:7.19, m2:56.47, sd2:8.92, r:0.73, n:40 },
    "List B Free-Recall Trial": { m1:6.24, sd1:1.58, m2:6.82, sd2:1.57, r:0.26, n:40 },
    "Short-Delay Free Recall": { m1:9.89, sd1:2.44, m2:12.24, sd2:2.62, r:0.77, n:40 },
    "Short-Delay Cued Recall": { m1:11.13, sd1:2.21, m2:12.15, sd2:3.14, r:0.49, n:40 },
    "Long-Delay Free Recall": { m1:10.77, sd1:2.32, m2:12.3, sd2:2.54, r:0.62, n:40 },
    "Long-Delay Cued Recall": { m1:11.26, sd1:2.18, m2:12.74, sd2:2.53, r:0.69, n:40 },
    "Semantic Cluster Ratio": { m1:1.49, sd1:0.43, m2:1.83, sd2:0.6, r:0.58, n:40 },
    "Perseverations": { m1:4.55, sd1:4.32, m2:5.6, sd2:5.7, r:0.32, n:40 },
    "Free-Recall Intrusions": { m1:0.97, sd1:1.44, m2:1.64, sd2:2.67, r:0.56, n:40 },
    "Cued-Recall Intrusions": { m1:0.69, sd1:1.14, m2:0.56, sd2:1.02, r:0.17, n:40 },
    "Recognition Hits": { m1:14.21, sd1:1.02, m2:14.67, sd2:0.82, r:0.24, n:40 },
    "Discriminability": { m1:96.71, sd1:3.53, m2:98.06, sd2:3.25, r:0.37, n:40 },
    "False Positives": { m1:0.55, sd1:0.86, m2:0.55, sd2:1.03, r:0.35, n:40 }
  },
  "CVLT-C Subtests (Raw Scores) · Age 16": {
    "List A Trials 1-5 Total": { m1:53.53, sd1:6.15, m2:62.94, sd2:10.94, r:0.61, n:31 },
    "List B Free-Recall Trial": { m1:6.67, sd1:1.7, m2:7.6, sd2:2.19, r:0.66, n:31 },
    "Short-Delay Free Recall": { m1:11.61, sd1:2.09, m2:13.52, sd2:1.91, r:0.48, n:31 },
    "Short-Delay Cued Recall": { m1:12, sd1:1.52, m2:13.71, sd2:1.59, r:0.59, n:31 },
    "Long-Delay Free Recall": { m1:11.9, sd1:1.9, m2:13.5, sd2:1.87, r:0.6, n:31 },
    "Long-Delay Cued Recall": { m1:12.57, sd1:1.72, m2:14, sd2:1.58, r:0.59, n:31 },
    "Semantic Cluster Ratio": { m1:1.55, sd1:0.53, m2:2.3, sd2:0.66, r:0.53, n:31 },
    "Perseverations": { m1:3.84, sd1:3.03, m2:5.42, sd2:5.76, r:0.31, n:31 },
    "Free-Recall Intrusions": { m1:2.28, sd1:4.22, m2:2.19, sd2:4.13, r:0.85, n:31 },
    "Cued-Recall Intrusions": { m1:0.66, sd1:1.68, m2:0.84, sd2:1.59, r:0.74, n:31 },
    "Recognition Hits": { m1:14.59, sd1:0.82, m2:14.79, sd2:0.62, r:0.8, n:31 },
    "Discriminability": { m1:97.2, sd1:4.17, m2:98.9, sd2:2.23, r:0.78, n:31 },
    "False Positives": { m1:0.68, sd1:1.51, m2:0.35, sd2:0.91, r:0.78, n:31 }
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
    "Free Sorting Confirmed Sorts": { m1:10.24, sd1:2.77, m2:11.31, sd2:2.44, r:0.51, n:101 },
    "Free Sorting Description Total Score": { m1:10.11, sd1:2.72, m2:11.27, sd2:2.51, r:0.5, n:101 },
    "Sort Recognition Total Description Score": { m1:10.26, sd1:2.92, m2:11.16, sd2:3.13, r:0.6, n:101 }
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
    "Total Achievement Score": { m1:10.35, sd1:3.21, m2:11.66, sd2:2.94, r:0.44, n:83 }
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
    "Combined Number + Letter": { m1:9.63, sd1:3.33, m2:11.05, sd2:2.71, r:0.66, n:98 }
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
    "Total Weighted Achievement": { m1:9.92, sd1:3.05, m2:10.23, sd2:2.92, r:0.24, n:101 },
    "Initial Abstraction Score": { m1:9.73, sd1:2.35, m2:9.66, sd2:2.38, r:0.43, n:101 }
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
    "Letter Fluency": { m1:9.62, sd1:3.14, m2:10.1, sd2:3.51, r:0.8, n:101 },
    "Category Fluency": { m1:9.83, sd1:3.25, m2:10.3, sd2:3.38, r:0.79, n:101 },
    "Category Switching": { m1:9.86, sd1:3.39, m2:9.87, sd2:3.58, r:0.52, n:101 },
    "Switching Accuracy": { m1:10.41, sd1:2.96, m2:10.54, sd2:3.63, r:0.36, n:101 }
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
    "Total First Trial Consistently Correct": { m1:10.03, sd1:2.85, m2:11.54, sd2:3.27, r:0.7, n:101 }
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
    "Total Achievement Score": { m1:9.77, sd1:3.07, m2:10.57, sd2:3.04, r:0.76, n:101 }
  },
  "D-KEFS Advanced Trail Making · All Ages": {
    "NS Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.8, m2:10.4, sd2:2.8, r:0.53, n:224 },
    "LS Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.9, m2:10.4, sd2:3, r:0.52, n:224 },
    "SW Mean Correct (With Errors) Connection Time": { m1:10, sd1:3.1, m2:11.1, sd2:3.1, r:0.64, n:224 },
    "SWD Mean Correct (With Errors) Connection Time": { m1:9.9, sd1:3.1, m2:11.1, sd2:3.1, r:0.69, n:224 },
    "SWM Mean Correct (With Errors) Connection Time": { m1:10.1, sd1:2.9, m2:10.7, sd2:3.3, r:0.66, n:224 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.2, m2:10.7, sd2:2.8, r:0.48, n:224 },
    "SW Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:3, m2:11.1, sd2:3, r:0.66, n:224 },
    "SWD Mean Pure Response Speed (Correct or Error Connections)": { m1:9.9, sd1:3.1, m2:11, sd2:3.2, r:0.7, n:224 },
    "SWM Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:2.9, m2:10.6, sd2:3.3, r:0.7, n:224 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:102.3, sd1:14.3, m2:102.5, sd2:15.1, r:0.67, n:224 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:99.6, sd1:15.8, m2:105.2, sd2:16.4, r:0.79, n:224 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:99.8, sd1:15.4, m2:105, sd2:16.2, r:0.81, n:224 },
    "Multitasking Index": { m1:99.7, sd1:16.4, m2:105.5, sd2:16.4, r:0.73, n:224 }
  },
  "D-KEFS Advanced Trail Making · Ages 8-18": {
    "NS Mean Correct (With Errors) Connection Time": { m1:10.1, sd1:3, m2:10.6, sd2:2.9, r:0.47, n:91 },
    "LS Mean Correct (With Errors) Connection Time": { m1:10.2, sd1:2.8, m2:10.6, sd2:2.7, r:0.47, n:91 },
    "SW Mean Correct (With Errors) Connection Time": { m1:9.8, sd1:3.2, m2:11.5, sd2:3.3, r:0.62, n:91 },
    "SWD Mean Correct (With Errors) Connection Time": { m1:9.8, sd1:3.4, m2:11.5, sd2:3.3, r:0.77, n:91 },
    "SWM Mean Correct (With Errors) Connection Time": { m1:9.6, sd1:3.1, m2:10.6, sd2:3.4, r:0.76, n:91 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.3, m2:10.8, sd2:2.9, r:0.4, n:91 },
    "SW Mean Pure Response Speed (Correct or Error Connections)": { m1:9.7, sd1:3.1, m2:11.5, sd2:3.3, r:0.65, n:91 },
    "SWD Mean Pure Response Speed (Correct or Error Connections)": { m1:9.9, sd1:3.2, m2:11.5, sd2:3.2, r:0.75, n:91 },
    "SWM Mean Pure Response Speed (Correct or Error Connections)": { m1:9.5, sd1:2.9, m2:10.6, sd2:3.3, r:0.74, n:91 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:100.5, sd1:14.1, m2:103.4, sd2:14.6, r:0.64, n:91 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:98.2, sd1:16.7, m2:106.7, sd2:17.5, r:0.83, n:91 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:98.3, sd1:15.7, m2:106.8, sd2:16.7, r:0.84, n:91 },
    "Multitasking Index": { m1:98.4, sd1:16.9, m2:106.8, sd2:17, r:0.77, n:91 }
  },
  "D-KEFS Advanced Trail Making · Ages 19-59": {
    "NS Mean Correct (With Errors) Connection Time": { m1:11, sd1:2.5, m2:10.8, sd2:2.6, r:0.39, n:67 },
    "LS Mean Correct (With Errors) Connection Time": { m1:11, sd1:3, m2:10.7, sd2:3.1, r:0.46, n:67 },
    "SW Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:2.9, m2:11.2, sd2:2.5, r:0.59, n:67 },
    "SWD Mean Correct (With Errors) Connection Time": { m1:10.4, sd1:2.8, m2:11.3, sd2:2.8, r:0.48, n:67 },
    "SWM Mean Correct (With Errors) Connection Time": { m1:10.5, sd1:2.5, m2:10.7, sd2:2.9, r:0.44, n:67 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10.1, sd1:2.8, m2:10.8, sd2:2.7, r:0.52, n:67 },
    "SW Mean Pure Response Speed (Correct or Error Connections)": { m1:10.4, sd1:2.9, m2:11.2, sd2:2.6, r:0.59, n:67 },
    "SWD Mean Pure Response Speed (Correct or Error Connections)": { m1:10.5, sd1:2.7, m2:11.2, sd2:2.8, r:0.55, n:67 },
    "SWM Mean Pure Response Speed (Correct or Error Connections)": { m1:10.5, sd1:2.4, m2:10.7, sd2:2.9, r:0.54, n:67 },
    "Combined Sequencing Mean Correct (With Errors) Connection Time Index": { m1:105.4, sd1:12.9, m2:104.2, sd2:15.1, r:0.55, n:67 },
    "Combined Switching Mean Correct (With Errors) Connection Time Index": { m1:102.3, sd1:13.4, m2:106, sd2:14, r:0.67, n:67 },
    "Combined Switching Mean Pure Response Speed (Correct or Error Connections) Idx": { m1:102.9, sd1:13.1, m2:106, sd2:13.7, r:0.69, n:67 },
    "Multitasking Index": { m1:101.9, sd1:14.1, m2:106.2, sd2:14.2, r:0.63, n:67 }
  },
  "D-KEFS Advanced Trail Making · Ages 60-90": {
    "NS Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:2.7, m2:9.8, sd2:2.8, r:0.7, n:66 },
    "LS Mean Correct (With Errors) Connection Time": { m1:10.3, sd1:3.2, m2:10, sd2:3.3, r:0.62, n:66 },
    "SW Mean Correct (With Errors) Connection Time": { m1:10, sd1:3.3, m2:10.4, sd2:3.5, r:0.7, n:66 },
    "SWD Mean Correct (With Errors) Connection Time": { m1:9.5, sd1:3.1, m2:10.4, sd2:2.9, r:0.76, n:66 },
    "SWM Mean Correct (With Errors) Connection Time": { m1:10.2, sd1:3.1, m2:10.7, sd2:3.5, r:0.73, n:66 },
    "Combined Switching Total Active Errors (Set-Loss + Sequencing)": { m1:10, sd1:3.4, m2:10.5, sd2:2.9, r:0.52, n:66 },
    "SW Mean Pure Response Speed (Correct or Error Connections)": { m1:10, sd1:3, m2:10.5, sd2:3.1, r:0.72, n:66 },
    "SWD Mean Pure Response Speed (Correct or Error Connections)": { m1:9.3, sd1:3.2, m2:10, sd2:3.4, r:0.77, n:66 },
    "SWM Mean Pure Response Speed (Correct or Error Connections)": { m1:10.1, sd1:3.4, m2:10.5, sd2:3.6, r:0.79, n:66 },
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
  "D-KEFS Advanced Color-Word Interference · All Ages": {
    "Color Identification Net Correct Responses": { m1:10.2, sd1:3, m2:10.5, sd2:3.2, r:0.75, n:224 },
    "Word Identification Net Correct Responses": { m1:10.2, sd1:2.8, m2:10.4, sd2:3, r:0.78, n:224 },
    "Inhibition Net Correct Responses": { m1:10, sd1:3.2, m2:10.7, sd2:3.2, r:0.81, n:224 },
    "Inhibition/Switching Net Correct Responses": { m1:9.8, sd1:3, m2:11, sd2:3, r:0.72, n:224 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:9.7, sd1:3.2, m2:10.4, sd2:3, r:0.6, n:224 },
    "Combined Color and Word Identification Net Correct Responses Index": { m1:101.1, sd1:14.6, m2:102.4, sd2:15.7, r:0.81, n:224 },
    "Multitasking Index": { m1:99.6, sd1:15.8, m2:105, sd2:16.2, r:0.82, n:224 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:100.9, sd1:14.1, m2:104.8, sd2:15.6, r:0.84, n:224 }
  },
  "D-KEFS Advanced Color-Word Interference · Ages 8-18": {
    "Color Identification Net Correct Responses": { m1:10.3, sd1:3, m2:10.4, sd2:3.4, r:0.73, n:91 },
    "Word Identification Net Correct Responses": { m1:10.3, sd1:2.8, m2:10.4, sd2:3.1, r:0.75, n:91 },
    "Inhibition Net Correct Responses": { m1:10.1, sd1:3.4, m2:11, sd2:3.4, r:0.8, n:91 },
    "Inhibition/Switching Net Correct Responses": { m1:9.4, sd1:3.4, m2:10.8, sd2:3.2, r:0.79, n:91 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:9.3, sd1:3.4, m2:10.1, sd2:3.2, r:0.61, n:91 },
    "Combined Color and Word Identification Net Correct Responses Index": { m1:101.6, sd1:14.4, m2:102.3, sd2:17, r:0.81, n:91 },
    "Multitasking Index": { m1:98.4, sd1:17.3, m2:105.3, sd2:17.4, r:0.84, n:91 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:100.2, sd1:15.3, m2:105.5, sd2:16.7, r:0.8, n:91 }
  },
  "D-KEFS Advanced Color-Word Interference · Ages 19-59": {
    "Color Identification Net Correct Responses": { m1:10.7, sd1:3.2, m2:11, sd2:3.3, r:0.8, n:67 },
    "Word Identification Net Correct Responses": { m1:10.4, sd1:3.2, m2:10.3, sd2:3.3, r:0.82, n:67 },
    "Inhibition Net Correct Responses": { m1:10.2, sd1:2.7, m2:10.7, sd2:3, r:0.77, n:67 },
    "Inhibition/Switching Net Correct Responses": { m1:10.5, sd1:2.7, m2:11.6, sd2:2.8, r:0.62, n:67 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:10, sd1:3.1, m2:10.3, sd2:3, r:0.7, n:67 },
    "Combined Color and Word Identification Net Correct Responses Index": { m1:102.8, sd1:16.6, m2:103.3, sd2:17, r:0.84, n:67 },
    "Multitasking Index": { m1:102.2, sd1:13.4, m2:106.8, sd2:15.1, r:0.74, n:67 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:102.9, sd1:12.9, m2:106.6, sd2:14, r:0.79, n:67 }
  },
  "D-KEFS Advanced Color-Word Interference · Ages 60-90": {
    "Color Identification Net Correct Responses": { m1:9.7, sd1:2.6, m2:10.2, sd2:2.6, r:0.72, n:66 },
    "Word Identification Net Correct Responses": { m1:9.9, sd1:2.3, m2:10.5, sd2:2.3, r:0.77, n:66 },
    "Inhibition Net Correct Responses": { m1:9.7, sd1:3.2, m2:10.2, sd2:3, r:0.84, n:66 },
    "Inhibition/Switching Net Correct Responses": { m1:9.8, sd1:2.9, m2:10.7, sd2:3, r:0.72, n:66 },
    "Combined Inhibition and Inhibition/Switching Total Errors": { m1:10, sd1:3, m2:10.8, sd2:2.8, r:0.47, n:66 },
    "Combined Color and Word Identification Net Correct Responses Index": { m1:98.6, sd1:12.4, m2:101.6, sd2:12.2, r:0.79, n:66 },
    "Multitasking Index": { m1:98.6, sd1:15.8, m2:102.8, sd2:15.5, r:0.85, n:66 },
    "Combined Inhib. and Inhib./Switching Mean Pure Response Time Index": { m1:99.8, sd1:13.8, m2:101.7, sd2:15.5, r:0.9, n:66 }
  },
  "D-KEFS Advanced Tower · All Ages": {
    "Global Performance Score": { m1:9.8, sd1:3.2, m2:11.2, sd2:3.3, r:0.56, n:224 },
    "Adjusted Mean Pure Response Time": { m1:10.1, sd1:3.2, m2:11.8, sd2:3.2, r:0.75, n:224 },
    "Adjusted Mean Unproductive Responses": { m1:9.9, sd1:3, m2:11.3, sd2:2.9, r:0.51, n:224 }
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
    "Global Performance Index": { m1:98.2, sd1:16.1, m2:105.8, sd2:19.5, r:0.59, n:224 },
    "Total Number of Perseverative Responses": { m1:10, sd1:3.2, m2:12.1, sd2:3.5, r:0.53, n:224 },
    "Percent Perseverative Responses": { m1:9.9, sd1:3.3, m2:12, sd2:3.5, r:0.46, n:224 },
    "Total Number of Perseverative Errors": { m1:9.8, sd1:3.2, m2:12.1, sd2:3.6, r:0.52, n:224 },
    "Percent Perseverative Errors": { m1:9.9, sd1:3.3, m2:12.2, sd2:3.7, r:0.47, n:224 },
    "Total Number of Errors": { m1:9.7, sd1:3.2, m2:11.6, sd2:3.8, r:0.6, n:224 },
    "Percent Correct Responses": { m1:9.8, sd1:3.2, m2:11.8, sd2:3.7, r:0.59, n:224 },
    "Total Number of Nonperseverative Errors": { m1:9.7, sd1:3.5, m2:10.7, sd2:3.6, r:0.54, n:224 },
    "Percent Nonperseverative Errors": { m1:9.7, sd1:3.4, m2:10.6, sd2:3.5, r:0.52, n:224 },
    "Total Number of Conceptual Level Responses": { m1:9.7, sd1:3.1, m2:10.4, sd2:3.3, r:0.34, n:224 },
    "Percent Conceptual Level Responses": { m1:9.7, sd1:3.2, m2:11.6, sd2:3.8, r:0.6, n:224 }
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
    "Total Net Earnings": { m1:9.9, sd1:3.3, m2:12.6, sd2:3.7, r:0.56, n:224 },
    "Net Earnings Races 1-20": { m1:9.8, sd1:3.3, m2:12.9, sd2:3.4, r:0.39, n:224 },
    "Net Earnings Races 21-40": { m1:9.9, sd1:3.3, m2:10.8, sd2:2.4, r:0.55, n:224 },
    "Net Earnings Races 41-60": { m1:10, sd1:3.2, m2:11.2, sd2:2.9, r:0.53, n:224 }
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
    "Line Orientation": { m1:16.7, sd1:3, m2:16.9, sd2:2.9, r:0.72, n:55 },
    "Picture Naming": { m1:9.1, sd1:1, m2:9.2, sd2:0.9, r:0.73, n:55 },
    "Semantic Fluency": { m1:10, sd1:3.1, m2:10.8, sd2:3, r:0.67, rCorrected:0.65, n:55 },
    "Digit Span": { m1:9.8, sd1:2.7, m2:10, sd2:3.2, r:0.59, rCorrected:0.67, n:55 },
    "Coding": { m1:10.2, sd1:2.7, m2:11.2, sd2:3.1, r:0.75, rCorrected:0.79, n:55 },
    "List Recall": { m1:6.8, sd1:1.8, m2:8.2, sd2:1.8, r:0.66, n:55 },
    "List Recognition": { m1:19.9, sd1:0.5, m2:19.9, sd2:0.4, r:0.7, n:55 },
    "Story Recall": { m1:10, sd1:3, m2:11.7, sd2:3.1, r:0.48, rCorrected:0.49, n:55 },
    "Figure Recall": { m1:10.2, sd1:2.5, m2:11.2, sd2:3.1, r:0.58, rCorrected:0.71, n:55 }
  },
  "RBANS Subtests · Ages 20-89": {
    "List Learning": { m1:11.5, sd1:2.9, m2:11.2, sd2:3.3, r:0.49, rCorrected:0.52, n:40 },
    "Story Memory": { m1:11.6, sd1:1.8, m2:12.5, sd2:2.4, r:0.45, rCorrected:0.8, n:40 },
    "Figure Copy": { m1:9.6, sd1:2.8, m2:11.9, sd2:2.6, r:0.47, rCorrected:0.54, n:40 },
    "Line Orientation": { m1:16, sd1:3.4, m2:16.4, sd2:3.7, r:0.49, n:40 },
    "Picture Naming": { m1:9.8, sd1:0.4, m2:9.7, sd2:0.5, r:0.5, n:40 },
    "Semantic Fluency": { m1:11.1, sd1:2.9, m2:11.2, sd2:3.3, r:0.49, rCorrected:0.52, n:40 },
    "Digit Span": { m1:10.4, sd1:3.5, m2:10.1, sd2:3.7, r:0.73, rCorrected:0.63, n:40 },
    "Coding": { m1:10.8, sd1:2.5, m2:11.7, sd2:2.8, r:0.76, rCorrected:0.83, n:40 },
    "List Recall": { m1:6.2, sd1:2.4, m2:5.8, sd2:2.7, r:0.6, n:40 },
    "List Recognition": { m1:19.6, sd1:0.8, m2:19.8, sd2:0.5, r:0.27, n:40 },
    "Story Recall": { m1:11.6, sd1:2.3, m2:11.6, sd2:2.3, r:0.52, rCorrected:0.72, n:40 },
    "Figure Recall": { m1:10.4, sd1:3, m2:11.5, sd2:3, r:0.55, rCorrected:0.55, n:40 }
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
  "WAIS-IV Core Subtests · All Ages": {
    "Block Design": { m1:10.2, sd1:2.9, m2:11, sd2:2.8, r:0.79, rCorrected:0.8, n:298 },
    "Similarities": { m1:9.9, sd1:2.8, m2:10.4, sd2:2.8, r:0.83, rCorrected:0.87, n:298 },
    "Digit Span": { m1:10, sd1:2.9, m2:10.6, sd2:3, r:0.82, rCorrected:0.83, n:298 },
    "Matrix Reasoning": { m1:10.1, sd1:3.1, m2:10.5, sd2:3.1, r:0.76, rCorrected:0.74, n:298 },
    "Vocabulary": { m1:9.9, sd1:3, m2:10, sd2:3, r:0.9, rCorrected:0.89, n:298 },
    "Arithmetic": { m1:9.9, sd1:2.8, m2:10.4, sd2:2.9, r:0.8, rCorrected:0.83, n:298 },
    "Symbol Search": { m1:10.1, sd1:2.9, m2:11, sd2:3.3, r:0.8, rCorrected:0.81, n:298 },
    "Visual Puzzles": { m1:10, sd1:2.8, m2:10.9, sd2:3, r:0.69, rCorrected:0.74, n:298 },
    "Information": { m1:9.8, sd1:3, m2:10.5, sd2:3.2, r:0.91, rCorrected:0.9, n:298 },
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
  "WAIS-IV Indices · All Ages": {
    "Verbal Comprehension Index": { m1:99.3, sd1:14.4, m2:101.8, sd2:15, r:0.95, rCorrected:0.96, n:298 },
    "Perceptual Reasoning Index": { m1:100.4, sd1:13.8, m2:104.3, sd2:14.3, r:0.85, rCorrected:0.87, n:298 },
    "Working Memory Index": { m1:99.5, sd1:14, m2:102.6, sd2:14.7, r:0.87, rCorrected:0.88, n:298 },
    "Processing Speed Index": { m1:100.2, sd1:13.5, m2:104.6, sd2:14.9, r:0.84, rCorrected:0.87, n:298 },
    "Full Scale IQ": { m1:99.7, sd1:13.8, m2:104, sd2:15, r:0.95, rCorrected:0.96, n:298 }
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
    "Block Design No Time Bonus": { m1:10.3, sd1:2.9, m2:11, sd2:2.8, r:0.76, rCorrected:0.78, n:298 },
    "Digit Span Forward": { m1:9.9, sd1:2.8, m2:10.2, sd2:3, r:0.74, rCorrected:0.77, n:298 },
    "Digit Span Backward": { m1:10.2, sd1:2.9, m2:10.7, sd2:3.1, r:0.69, rCorrected:0.71, n:298 },
    "Digit Span Sequencing": { m1:9.9, sd1:2.9, m2:10.6, sd2:2.8, r:0.7, rCorrected:0.72, n:298 }
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
  "WAIS-IV Supplementary Subtests · All Ages": {
    "Letter-Number Sequencing": { m1:10.1, sd1:2.7, m2:10.5, sd2:3.1, r:0.76, rCorrected:0.8, n:298 },
    "Figure Weights": { m1:10, sd1:3, m2:10.8, sd2:3.2, r:0.76, rCorrected:0.77, n:298 },
    "Comprehension": { m1:10, sd1:3, m2:10.2, sd2:2.9, r:0.86, rCorrected:0.86, n:298 },
    "Cancellation": { m1:10.2, sd1:2.8, m2:10.8, sd2:3, r:0.74, rCorrected:0.78, n:298 },
    "Picture Completion": { m1:9.9, sd1:2.9, m2:11.8, sd2:3.3, r:0.74, rCorrected:0.77, n:298 }
  },
  "WISC-V Indices · All Ages": {
    "Verbal Comprehension Index": { m1:98.5, sd1:12.8, m2:101.6, sd2:13, r:0.91, n:215 },
    "Visuospatial Index": { m1:98.6, sd1:14.7, m2:105.3, sd2:15.1, r:0.84, n:217 },
    "Fluid Reasoning Index": { m1:98.7, sd1:13.6, m2:103.6, sd2:12.9, r:0.68, n:217 },
    "Working Memory Index": { m1:98.5, sd1:13.8, m2:100.9, sd2:13.8, r:0.79, n:217 },
    "Processing Speed Index": { m1:100.3, sd1:14.3, m2:108.2, sd2:16, r:0.81, n:213 },
    "Full Scale IQ": { m1:98.3, sd1:13.7, m2:104.3, sd2:13.8, r:0.91, n:212 }
  },
  "WISC-V Subtests · All Ages": {
    "Similarities": { m1:9.8, sd1:2.5, m2:10.6, sd2:2.5, r:0.82, n:213 },
    "Vocabulary": { m1:9.6, sd1:2.8, m2:10, sd2:2.8, r:0.89, n:217 },
    "Information": { m1:9.7, sd1:2.7, m2:10.3, sd2:2.7, r:0.85, n:218 },
    "Comprehension": { m1:10, sd1:2.9, m2:10.2, sd2:2.8, r:0.81, n:214 },
    "Block Design": { m1:9.6, sd1:2.9, m2:10.8, sd2:3.1, r:0.79, n:208 },
    "Visual Puzzles": { m1:9.9, sd1:2.8, m2:11, sd2:2.9, r:0.78, n:210 },
    "Matrix Reasoning": { m1:9.6, sd1:2.4, m2:10.6, sd2:2.6, r:0.65, n:202 },
    "Figure Weights": { m1:10, sd1:2.6, m2:10.5, sd2:2.6, r:0.76, n:204 },
    "Picture Concepts": { m1:9.8, sd1:2.7, m2:10.7, sd2:2.9, r:0.63, n:203 },
    "Arithmetic": { m1:9.8, sd1:2.5, m2:10.2, sd2:2.6, r:0.75, n:205 },
    "Digit Span": { m1:9.8, sd1:2.8, m2:10.1, sd2:3, r:0.79, n:214 },
    "Picture Span": { m1:9.7, sd1:2.5, m2:10.1, sd2:2.6, r:0.72, n:208 },
    "Letter-Number Sequencing": { m1:9.8, sd1:2.7, m2:10.2, sd2:2.8, r:0.77, n:212 },
    "Coding": { m1:10, sd1:2.9, m2:11.3, sd2:3.1, r:0.79, n:216 },
    "Symbol Search": { m1:10, sd1:2.7, m2:11.5, sd2:3.2, r:0.76, n:209 },
    "Cancellation": { m1:9.8, sd1:2.9, m2:11.1, sd2:3.2, r:0.79, n:209 }
  },
  "WMS-IV Indices · Ages 16-69": {
    "Auditory Memory Index": { m1:100.1, sd1:14.1, m2:111.6, sd2:14.4, r:0.81, rCorrected:0.83, n:168 },
    "Visual Memory Index": { m1:100, sd1:14.8, m2:112.1, sd2:16.6, r:0.8, rCorrected:0.81, n:144 },
    "Visual Working Memory Index": { m1:99.5, sd1:14.4, m2:103.8, sd2:15.6, r:0.82, rCorrected:0.83, n:171 },
    "Immediate Memory Index": { m1:99.9, sd1:14.9, m2:112.3, sd2:15.6, r:0.81, rCorrected:0.81, n:154 },
    "Delayed Memory Index": { m1:100.4, sd1:13.9, m2:114.1, sd2:15, r:0.79, rCorrected:0.82, n:150 }
  },
  "WMS-IV Indices · Ages 65-90": {
    "Auditory Memory Index": { m1:101.5, sd1:12.9, m2:112.1, sd2:14.5, r:0.82, rCorrected:0.87, n:69 },
    "Visual Memory Index": { m1:101.6, sd1:14.6, m2:112.6, sd2:17, r:0.79, rCorrected:0.8, n:70 },
    "Immediate Memory Index": { m1:101.5, sd1:13.9, m2:113.9, sd2:14.2, r:0.84, rCorrected:0.86, n:69 },
    "Delayed Memory Index": { m1:101.7, sd1:13.1, m2:112.7, sd2:14.4, r:0.8, rCorrected:0.85, n:71 }
  },
  "WMS-IV Subtests · Ages 16-69": {
    "Logical Memory I": { m1:10.3, sd1:2.9, m2:12.2, sd2:2.6, r:0.72, rCorrected:0.74, n:173 },
    "Logical Memory II": { m1:10.3, sd1:2.8, m2:12.6, sd2:2.9, r:0.67, rCorrected:0.71, n:172 },
    "Verbal Paired Associates I": { m1:9.8, sd1:3.1, m2:12.1, sd2:3.4, r:0.76, rCorrected:0.74, n:171 },
    "Verbal Paired Associates II": { m1:9.8, sd1:3, m2:10.8, sd2:2.7, r:0.76, rCorrected:0.76, n:170 },
    "Verbal Paired Associates II - Word Recall": { m1:10, sd1:2.9, m2:10.9, sd2:3.2, r:0.74, rCorrected:0.76, n:170 },
    "Designs I": { m1:10, sd1:2.9, m2:11.1, sd2:3.4, r:0.73, rCorrected:0.75, n:157 },
    "Designs I - Content": { m1:10.1, sd1:3, m2:11.2, sd2:3.3, r:0.64, rCorrected:0.64, n:157 },
    "Designs I - Spatial": { m1:10.2, sd1:2.9, m2:10.9, sd2:3.1, r:0.56, rCorrected:0.59, n:157 },
    "Designs II": { m1:10.2, sd1:2.7, m2:11.9, sd2:3.2, r:0.72, rCorrected:0.77, n:151 },
    "Designs II - Content": { m1:10.3, sd1:3, m2:11.5, sd2:3.4, r:0.64, rCorrected:0.64, n:151 },
    "Designs II - Spatial": { m1:10.3, sd1:2.7, m2:11.6, sd2:2.6, r:0.5, rCorrected:0.6, n:151 },
    "Visual Reproduction I": { m1:10, sd1:2.8, m2:11.9, sd2:2.8, r:0.62, rCorrected:0.67, n:170 },
    "Visual Reproduction II": { m1:10.1, sd1:2.8, m2:12.9, sd2:3, r:0.59, rCorrected:0.64, n:169 },
    "Spatial Addition": { m1:9.9, sd1:2.8, m2:10.7, sd2:3, r:0.74, rCorrected:0.77, n:172 },
    "Symbol Span": { m1:10, sd1:3, m2:10.6, sd2:3.1, r:0.72, rCorrected:0.72, n:172 }
  },
  "WMS-IV Subtests · Ages 65-90": {
    "Logical Memory I": { m1:10, sd1:2.9, m2:12, sd2:3.3, r:0.77, rCorrected:0.79, n:69 },
    "Logical Memory II": { m1:10, sd1:2.7, m2:12.1, sd2:2.8, r:0.71, rCorrected:0.77, n:71 },
    "Verbal Paired Associates I": { m1:10.4, sd1:2.8, m2:12.1, sd2:2.9, r:0.76, rCorrected:0.79, n:71 },
    "Verbal Paired Associates II": { m1:10.4, sd1:2.7, m2:11.5, sd2:2.7, r:0.77, rCorrected:0.81, n:71 },
    "Verbal Paired Associates II - Word Recall": { m1:10.5, sd1:2.8, m2:11.7, sd2:2.7, r:0.72, rCorrected:0.76, n:71 },
    "Visual Reproduction I": { m1:10.2, sd1:3.1, m2:12, sd2:3.1, r:0.79, rCorrected:0.78, n:71 },
    "Visual Reproduction II": { m1:10.5, sd1:2.8, m2:12.3, sd2:3.1, r:0.64, rCorrected:0.69, n:71 },
    "Symbol Span": { m1:10.1, sd1:2.8, m2:10.7, sd2:3, r:0.69, rCorrected:0.73, n:69 }
  }
};
