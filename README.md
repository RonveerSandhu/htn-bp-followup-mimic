# Predicting the Absence of Recorded Blood Pressure Follow-up in Adults With Hypertension

Analysis code for:

> Sandhu R, Asghar S, Abu Hashish M, Ye J. *Predicting the Absence of Recorded Blood Pressure Follow-up in Adults With Hypertension: Development and Internal Validation.* (under review)

This repository contains the code used to build the study cohorts, fit the models, and produce every table and figure in the paper and supplement. **It contains no patient data.**

## What the study does

Using MIMIC-IV v2.2 (one academic medical center in Boston), we:

1. **Predict** which adults with hypertension will have no recorded blood pressure measurement within 90 days of an elevated reading (systolic > 140 mmHg), using gradient-boosted trees and two logistic regression comparators.
2. **Evaluate** the model at fixed outreach capacity (top 5–30% of predicted risk), with decision curve analysis.
3. **Audit** performance by race and ethnicity, insurance, sex, and age.
4. **Test** whether recorded contact patterns (trajectory classes and gaps in contact) are associated with cardiovascular hospitalization, using a day-730 landmark design.

## Data access

MIMIC-IV is not included and cannot be shared here. To run the code you need:

1. A credentialed PhysioNet account and a signed data use agreement: https://physionet.org/content/mimiciv/2.2/
2. **MIMIC-IV v2.2** – `hosp` module
3. **MIMIC-IV-ED v2.2** – `ed` module (used only to identify emergency visits)

The ICU module is not used.

## How to run

The notebooks were written for Google Colab with Google Drive.

1. Copy the MIMIC-IV `hosp` and `ed` folders to your Drive, and this repository to a Drive folder.
2. Open each notebook and set the paths in its first cell (`REPO_DIR`, `MIMIC_HOSP`, `MIMIC_ED`, `DATA_DIR`).
3. Run the notebooks **in order, 01 → 08**. Each one reads the files saved by the earlier ones.
4. Python packages: `pip install -r requirements.txt`. R packages (called through `rpy2`): `Rscript install_r_packages.R`.

Patient-level intermediate files (parquet, fitted models) are written to `DATA_DIR`, which must sit **outside** the repository; `.gitignore` blocks them. Only aggregate tables are written to `results/`.

Approximate run times on a standard Colab runtime: notebooks 01–03, 5–10 minutes each; 04, 30–45 minutes; 05, one to several hours (the 5- and 6-class trajectory models); 06 and 06b, 1.5–3 hours each; 07, a few minutes; 08, 10–20 minutes. Notebook 05 saves each trajectory model as it finishes, so a rerun after a disconnect picks up where it stopped (and a rerun after all models are fitted takes a few minutes).

**Results summary.** Each notebook ends with a block that prints the numbers it contributes to the paper.

**Repeatability.** Patients are sorted by ID before every step that depends on row order (the train/validation/test split, cross-validation folds, and multiple imputation), and every random step is seeded, so rerunning the notebooks gives the same results.

## Repository structure

```
htn-bp-followup-mimic/
├── README.md
├── LICENSE
├── requirements.txt
├── install_r_packages.R
├── .gitignore
├── src/
│   ├── htn_utils.py            # MIMIC connection, file IO, code lists (Table S1), study constants
│   └── landmark_models.py      # class order, covariate matrix, splines, Rubin's rules, power (notebooks 06–08)
├── notebooks/
│   ├── 01_cohort.ipynb
│   ├── 02_contact_and_landmark.ipynb
│   ├── 03_prediction_sample.ipynb
│   ├── 04_prediction_models.ipynb
│   ├── 05_trajectories.ipynb
│   ├── 06_association.ipynb
│   ├── 06b_sensitivity.ipynb
│   ├── 07_multimorbidity.ipynb
│   └── 08_final_tables_figures.ipynb
└── results/                    # aggregate tables and figures only (no patient-level rows)
```

## Where each result comes from

| Notebook | What it does | Paper items |
|---|---|---|
| `01_cohort` | Hypertension cohort, index date, exclusions, database coverage rule, death capture window, covariates (Elixhauser, eGFR, HbA1c, BMI), outcomes | Methods: data source and study population; Table S1 |
| `02_contact_and_landmark` | Unique contact days per 91-day quarter, longest gap in contact, day-730 landmark population and day-730-information sample | Figure 1A (association branch); Table S9 (association rows) |
| `03_prediction_sample` | First systolic reading > 140 mmHg, 90-day coverage and survival rules, 62 features cut at the prediction date | Figure 1A (prediction branch); Table S9 (prediction rows) |
| `04_prediction_models` | 60/20/20 split, tuning, isotonic recalibration, test-set evaluation, nested cross-validation, outreach at fixed capacity, decision curve, subgroup audit, group-specific thresholds, alternative targets, Shapley values | Table 2A–B; Figure 2A; Figure S1; Tables S5–S6 |
| `05_trajectories` | Group-based trajectory models (1–6 classes, R `lcmm`), class selection (4 classes selected by the rule; 5 classes kept as a sensitivity analysis), rule-based class naming, posterior probabilities | Table S2; data for Figure 3A |
| `06_association` | Pre-index covariates, multiple imputation (m = 30), cause-specific Cox, Fine–Gray, doubly robust, balance, gap spline and cutoffs, Uno's C, Gray's tests, proportional hazards, minimum detectable hazard ratios | Table 4; Table S7; data for Figure 3B–D and Figure S2 |
| `06b_sensitivity` | All Table S3 sensitivity analyses (other landmark samples, censoring weights, missing-data approaches), outcome components, heart-failure checks | Table S3; Table S8; data for Figure S3 |
| `07_multimorbidity` | Elixhauser co-occurrence network, Louvain communities, phenotype assignment, χ² / Cramér's V, interaction test | Supplementary Methods and Results: multimorbidity |
| `08_final_tables_figures` | Prediction-sample baseline table, subgroup confidence intervals, denominator reconciliation, every figure, all tables in one workbook (`paper_tables.xlsx`) | Table 1; Table 3; Table S4; Table S9; Figures 1–3 and S1–S3 |

## Software

- Python 3.13: pandas 2.2, DuckDB 1.3, lifelines 0.30, scikit-learn 1.6, XGBoost 3.4, statsmodels 0.15, SHAP 0.52, comorbidipy 0.5
- R 4.6: lcmm 2.2.2, cmprsk 2.2-12, survival 3.8

Random seeds are set in each notebook (trajectory models use seed 42).

## Notes

- **Recorded contact only.** MIMIC-IV covers one health system. A patient seen elsewhere looks the same as one who was not seen, so all outcomes describe *recorded* measurements and contact.
- **Dates are shifted** in MIMIC-IV; calendar time is available only as three-year groups, which is why validation uses a random split rather than a temporal one.

## Citation

If you use this code, please cite the paper above and MIMIC-IV:

> Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data.* 2023;10:1.

## License

Code is released under the MIT License (see `LICENSE`). MIMIC-IV data are governed by the PhysioNet Credentialed Health Data License.

## Contact

Ronveer Sandhu, Georgia Institute of Technology
