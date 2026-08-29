# CeCe Pipeline

The pipeline combines co-occurrence network analysis,causal inference (NOTEARS + average causal effect estimation), genome-scale metabolic modeling (MICOM) and machine learning (Random Forest + SHAP) into a composite per-species **CECE score**, which a genetic algorithm then uses to select optimal transplant consortia.

---

## 1. Data Sources

| Dataset | What it provides | Source | Reference |
|---|---|---|---|
| **HMP2 / IBDMDB** (Inflammatory Bowel Disease Multi'omics Database, part of the Integrative Human Microbiome Project) | WGS metagenomic species-level taxonomic profiles + longitudinal metadata for CD, UC and non-IBD subjects; used to build co-occurrence networks and the causal DAG | Portal: [ibdmdb.org](https://ibdmdb.org) (`HMP2_taxonomic_profiles.tsv.gz`, `HMP2_metadata.csv`). Raw shotgun sequencing reads are deposited in NCBI SRA under BioProject accessions **PRJNA398089, PRJNA400072, PRJNA389280** | Lloyd-Price et al. 2019, *Nature* 569:655–662, doi:10.1038/s41586-019-1237-9 |
| **Podlesny et al. 2022 FMT cohort** | 14 pooled FMT trials across 5 disease groups (IBD, rCDI, MetS, MDR, ICI) with donor/recipient/post-FMT abundance and clinical metadata; used to derive engraftment labels and train the ML model | *Cell Reports Medicine* 3(8):100711. Sequencing data deposited in the European Nucleotide Archive (ENA); per-sample accessions listed in the paper's key resources table | Podlesny et al. 2022, doi:10.1016/j.xcrm.2022.100711 |
| **Smillie et al. 2018 engraftment AUC table** | Genome-level engraftment-predictability (AUC) scores per species, merged in as an auxiliary feature (`smillie_engraftment_prob`) | Supplementary data of the published study | Smillie et al. 2018, *Cell Host & Microbe* 23(2):229–240, doi:10.1016/j.chom.2018.01.003 |
| **AGORA v1.03** genome-scale metabolic reconstructions | 818 strain-level GEM `.xml` files used to build per-donor microbial communities for MICOM flux-balance simulation | [Virtual Metabolic Human (VMH)](https://vmh.life/#downloadview) - AGORA v1.03, "without mucins" release | Magnúsdóttir et al. 2017, *Nat Biotechnol* 35:81–89, doi:10.1038/nbt.3703 |
| **Western diet exchange medium** | qiime-formatted flux medium (`western_diet_gut.qza`) used as the growth medium for all MICOM simulations | Zenodo record **3755182** | zenodo.org/records/3755182 |

---

## 2. Pipeline Structure

| # | Notebook | Runs on | Stage |
|---|---|---|---|
| 01 | `01_data_and_EDA.ipynb` | local | Data loading, filtering, transforms, EDA |
| 02 | `02_engraftment_labels.ipynb` | local | Engraftment label derivation, train/test split |
| 03a | `03a_fastspar_prep.ipynb` | local | FastSpar input prep- HMP2 (CD/UC/Healthy) |
| 03b | `03b_fastspar_prep.ipynb` | local | FastSpar input prep- Podlesny (5 disease groups) |
| 04a | `04a_fastspar_COLAB.ipynb` | Colab | SparCC/FastSpar- HMP2 networks |
| 04b | `04b_fastspar_COLAB.ipynb` | Colab | SparCC/FastSpar- Podlesny networks |
| 05 | `05_networks.ipynb` | local | Network construction, keystone scoring |
| 06 | `06_NOTEARS_causal_DAG_COLAB.ipynb` | Colab | Causal DAG (NOTEARS) |
| 07 | `07_causal_effects_CIs.ipynb` | local | Average causal effect (ACE) estimation |
| 08 | `08_micom_COLAB.ipynb` | Colab | Metabolic modeling (MICOM/AGORA) |
| 09 | `09_ML_SHAP_Ablation.ipynb` | local | ML model, SHAP, ablation study |
| 10 | `10_cece_score_and_GA.ipynb` | local | CECE score + genetic algorithm |
| 11 | `11_CECE_additions.ipynb` | local | Score calibration, DRCI feature test |
| 12 | `12_validation.ipynb` | local | Final validation & summary |

---

## 3. Notebook-by-Notebook

### NB01 (Data and EDA)
Loads HMP2 taxonomic profiles (1,638 samples × 1,479 taxa, 130 patients: CD 750 / UC 459 / non-IBD 429). Filters to species level (1,479 → 572 species), applies per-condition prevalence/abundance filtering and computes relative-abundance and CLR transforms. Cross-references against the Podlesny species list, then builds patient-level time-lagged (Xₜ → Xₜ₊₁) sample pairs (≥3 timepoints/patient, max 4-week gap) for the causal-DAG stage.

**Result:** 87 "network species" common to CD ∩ UC ∩ Healthy after filtering.

### NB02 (Engraftment Labels)
Loads Podlesny metadata (1,416 samples, 263 cases) and abundance (860 species). Derives binary engraftment labels (donor abundance > 0.1% and detected post-FMT, preferring ≥28-day follow-up samples), splits 80/20 by case stratified on disease group (210 / 53 cases) and merges the Smillie AUC table as an auxiliary feature.

**Result:** 10,951 species×case records across 203 cases (52.33% engrafted); 8,697 train / 2,254 test; 272 species evaluated overall.

### NB03a / NB03b (FastSpar Prep)
Cleans and reformats count tables for FastSpar input: 03a for the three HMP2 network conditions (CD/UC/Healthy), 03b for the five Podlesny recipient disease groups (IBD, rCDI, MetS, MDR, ICI). Removes negative values, drops all-zero columns, and validates integer/non-negative format.

### NB04a / NB04b (FastSpar-Colab)
Compiles FastSpar (C++, [scwatts/fastspar](https://github.com/scwatts/fastspar)) from source and runs SparCC correlation with 100-bootstrap significance testing: 04a for the HMP2 conditions, 04b for the Podlesny disease groups.

### NB05 (Networks)
Symmetrizes and cleans the correlation/p-value matrices, builds NetworkX graphs (edge if \|r\| > 0.15 and p < 0.05) and scores keystone species by weighted centrality (0.5·betweenness + 0.3·degree + 0.2·eigenvector) with greedy-modularity community detection.

**Result:**

| Condition | Nodes | Edges | Mutualism | Competition | Modularity |
|---|---|---|---|---|---|
| CD | 86 | 636 | 353 | 283 | 0.244 |
| UC | 86 | 834 | 431 | 403 | 0.233 |
| Healthy | 87 | 1,040 | 519 | 521 | 0.193 |

193 edges are shared across all three conditions.

### NB06 (NOTEARS Causal DAG- Colab)
Fits linear NOTEARS (temporal form, Xₜ → Xₜ₊₁) on the union of the top-15 keystones per condition (30 species total), separately for CD and UC. Threshold-sweeps edge weight (0.10–0.30) to find the smallest value giving an acyclic graph, breaks any remaining cycles by removing the lowest-weight edge in each cycle and checks edge stability via 20-resample bootstrap. A sensitivity check (subsampling CD to UC's sample size, 5 seeds) shows the gap is at least partly a sample-size artifact: CD's edge count rises from 22 (at its full n=568) to a mean of 33.4 (range 29-39 across 5 seeds) once subsampled down to UC's n=347, close to UC's own edge count at full sample size. 

**Result:** acyclic, bootstrap-annotated causal DAGs saved for CD and UC.

### NB07 (Causal Effects & CIs)
Estimates each species' average causal effect (ACE) of donor abundance on engraftment via DAG-parent-adjusted logistic regression with 300-resample bootstrap confidence intervals.

**Result:**

| DAG | Species estimated | Significant (uncorrected p<0.05) | FDR q<0.05 | Bonferroni |
|---|---|---|---|---|
| CD | 131 / 264 | 33 | 25 | 0 |
| UC | 37 / 148 | 19 | 18 | 0 |

CD-ACE vs. UC-ACE is only weakly correlated across the 37 species estimated in both (Spearman ρ = 0.085, p = 0.62);causal drivers of engraftment look largely condition-specific.

### NB08 (MICOM-Colab)
Builds per-donor microbial communities from 818 AGORA v1.03 GEM files (genus/species-matched) and runs MICOM cooperative-tradeoff flux balance analysis (fraction = 0.5, Western-diet medium) per donor sample. 73 donors were selected (capped at 15/disease group, 13 for ICI); 67 simulated successfully.

**Result:** 1,100 growth-rate records (118 species) and 56,703 exchange-flux records (308 metabolites: 10,692 secretion / 46,011 consumption fluxes). SCFA production scores (butyrate, propionate, acetate, lactate, succinate, formate) extracted via keyword-matched exchange reactions.

### NB09 (ML, SHAP & Ablation)
Assembles a 23-feature table (causal ACE, ecological keystone-conflict / niche-overlap, MICOM resource-accessibility / crossfeeding / SCFA, disease-group, donor abundance) and trains Random Forest, XGBoost and Logistic Regression.

**Result- test set performance (2,254 records, 52.44% engrafted):**

| Model | AUC-ROC | Precision | Recall | F1 |
|---|---|---|---|---|
| **Random Forest (primary)** | **0.7773** | 0.7235 | 0.7174 | 0.7205 |
| XGBoost | 0.7728 | 0.7207 | 0.7030 | 0.7118 |
| Logistic Regression | 0.7227 | 0.6662 | 0.7462 | 0.7039 |

Bootstrap CI: RF vs. XGBoost overlaps zero (not distinguishable); RF vs. LogReg excludes zero (RF is meaningfully better). **Ablation:** removing `donor_abundance` costs the most (ΔAUC = −0.021); the other feature groups (ecological, causal, metabolic, SCFA, disease) each move full-model AUC by ≤0.002 alone, but all show standalone signal above chance in isolation (isolated-group AUCs 0.70–0.74). SHAP (TreeExplainer) computed for feature-level interpretability.

### NB10 (CECE Score & Genetic Algorithm)
Builds a composite **CECE score** per species from six min-max-scaled components:

| Component | Weight |
|---|---|
| Historical engraftment rate (Podlesny-observed × Smillie AUC blend) | 30% |
| Causal ACE | 20% |
| Keystone conflict (inverted) | 15% |
| Resource accessibility | 15% |
| Niche overlap | 10% |
| Crossfeeding synergy | 10% |

**Result:** 201 species scored, range 0.174–0.654. Top species: *Blautia wexlerae* (0.654), *Anaerostipes hadrus* (0.646), *Faecalibacterium prausnitzii* (0.612), *Ruminococcus bromii* (0.589). A genetic algorithm (population 100, 60 generations, 5-20 species per consortium) selects the donor-specific consortium maximizing RF-predicted engraftment probability, run across 20 donors. **Universal engrafters** (species chosen for ≥50% of donors): *F. prausnitzii* (70%), *A. hadrus* (55%), *Dorea longicatena* (50%), *B. wexlerae* (50%) — 4 species meet the bar. An eco-aware GA variant (adds a pairwise network-compatibility bonus) is compared to the baseline via paired Wilcoxon signed-rank tests.

### NB11 (CECE Additions)
Calibrates CECE score weights against observed test-set engraftment rate (Spearman-optimized) as a comparison to the manual NB10 weights. Tests a **Donor-Recipient Compatibility Index (DRCI)**, an ACE-weighted donor-abundance × recipient-absence mismatch score, as a candidate 24th ML feature.

**Result:** DRCI ranks #1 of 24 features by SHAP importance when included, but adding it slightly hurts held-out performance:

| Model | AUROC | AUPRC |
|---|---|---|
| Without DRCI | 0.7773 (95% CI 0.759–0.797) | 0.7916 |
| With DRCI | 0.7715 (95% CI 0.753–0.791) | 0.7853 |

DRCI is reported but **not adopted** in the final scoring/model. A SHAP-weighted GA fitness variant is also explored as an alternative to the manually-set CECE weights.

### NB12 (Validation)
Final held-out evaluation of the NB09 Random Forest model, plus project-wide summary statistics and limitation checks.

**Result:** test AUC-ROC 0.7773, AUPRC 0.7916, Precision 0.7235, Recall 0.7174, F1 0.7205 (2,254 records, 272 species).

**CECE score validity**: Spearman r=0.625, Pearson r=0.619, both p<0.0001, comparing the CECE score to observed test-set engraftment rate across the 138 species with ≥3 test observations. A calibration curve is plotted alongside.

**Per-disease-group AUC** (existing pooled model applied to each subgroup, no retraining; bootstrap 95% CI, 1,000 resamples): confirms pooling isn't masking uneven subgroup performance:

| Disease group | n | AUC | 95% CI |
|---|---|---|---|
| MDR | 284 | 0.672 | [0.607, 0.731] |
| IBD | 310 | 0.695 | [0.637, 0.756] |
| rCDI | 306 | 0.736 | [0.681, 0.793] |
| MetS | 1,175 | 0.813 | [0.791, 0.836] |
| ICI | 179 | 0.840 | [0.776, 0.895] |

Range: 0.672 (MDR) - 0.840 (ICI), against an overall pooled AUC of 0.7773.

**CD vs UC within IBD** (subtype recovered from Podlesny's `Disease` metadata column; same pooled model, no retraining):

| Subtype | n | AUC | 95% CI |
|---|---|---|---|
| UC | 227 | 0.675 | [0.604, 0.747] |
| CD | 83 | 0.783 | [0.670, 0.879] |

The CD/UC gap within IBD (0.783 vs 0.675) is real and not fully explained by the diagnostics tested (see README.md, Limitations).


---

## 4. Results Summary

| Metric | Value |
|---|---|
| Podlesny FMT cases used | 203 (of 263 total; 210/53 train/test split) |
| Engraftment records | 10,951 (8,697 train / 2,254 test), 52.33% engrafted |
| Species evaluated | 272 |
| HMP2 network species (CD ∩ UC ∩ Healthy) | 87 (from 572 species-level taxa, 1,638 samples, 130 patients) |
| Co-occurrence network edges | CD 636, UC 834, Healthy 1,040 (193 shared) |
| NOTEARS keystone species | 30 |
| Species with ACE estimated | CD 131 (25 FDR-sig.), UC 37 (18 FDR-sig.) |
| MICOM donors simulated | 67 of 73 attempted (818 AGORA models) |
| Final ML model | Random Forest: AUC-ROC 0.7773, AUPRC 0.7916, F1 0.7205 |
| CECE score | 201 species scored (range 0.174–0.654) |
| Universal engrafters | 4 species (≥50% of donors); led by *F. prausnitzii* (70%) and *A. hadrus* (55%) |

---

## 5. Software Environment

See `requirements.txt` for pinned Python dependencies. Two components fall
outside pip:
- **FastSpar** (Notebooks 04a/04b)- compiled from source on Colab, not distributed via PyPI.
- **GLPK solver** — system package required by `cobra`/`micom`'s LP/QP steps in Notebook 08.
