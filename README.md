<h1 align="center"><strong>CeCe</strong></h1>

<p align="center">
  <em>
    A Multi-Layer Computational Framework for Predictive Synthetic
    Faecal Microbiota Transplantation Design
  </em>
</p>

<p align="center">
  <sub>
    Integrating Causal Inference, Ecological Network Analysis,
    Constraint-based Metabolic Modelling and Machine Learning
  </sub>
</p>

<p align="center">
  <sub>
    Undergraduate Dissertation · B.Sc. (Hons with Research) Microbiology · 2026
  </sub>
  <br>
  <sub>
    Built on a personal laptop · Local Python + Google Colab
  </sub>
</p>



## 1. About

Fecal microbiota transplantation (FMT) can transfer hundreds of microbial species from a donor into a recipient, yet successful microbial establishment is not determined by donor composition alone. A species that is abundant in a donor may still fail to establish if its ecological relationships, metabolic requirements or interactions with the recipient microbiome are unfavourable.

*CeCe* explores a recipient-aware approach to this problem.

Rather than asking simply:

> Which microbes does the donor contain?

the framework asks:

> *Which donor-derived microbes are most compatible with the recipient environment and can those microbes be assembled into a coherent synthetic consortium?*

## 2. Conceptual Framework

The CeCe workflow integrates multiple biological layers:

<p align="center">
  <img 
    src="figures/FMT_WORKFLOW.png" 
    alt="CECE-FMT conceptual framework showing the workflow from FMT microbiome data through ecological, causal and metabolic analysis to machine learning, compatibility scoring and synthetic consortium optimization"
    width="100%"
  >
</p>

The framework combines ecological, temporal, causal and metabolic information with machine learning to investigate donor-derived microbial establishment within recipient context.
















---

## What makes this more than "stitched-together tools"

Each of the four layers above is a well-established method on its own. What's actually being tested here is whether *chaining* them — DAG → causal effect → ML feature → GA fitness — produces something a single layer can't:

- **The layers are genuinely chained, not parallel.** The causal DAG's output *is* an input to the ACE estimation, which *is* a feature the Random Forest sees, which *is* what the genetic algorithm optimises against. Break any one link and the downstream ones change.
- **A disease-specific causal DAG was built to *test*, not assume, transferability.** Rather than reusing Crohn's-derived causal effects for ulcerative colitis, a UC-native DAG was fit on its own longitudinal data. The result was honestly inconclusive (CD–UC ACE correlation r = 0.085, 95% CI crossing zero) — reported as a limitation, not smoothed over.
- **Ecological networks are disease-native, not borrowed.** Co-occurrence networks were rebuilt separately for each of the five FMT indications in the Podlesny cohort, rather than applying HMP2's IBD-only networks across unrelated diseases.
- **The composite CECE score is calibrated, not hand-weighted.** Its five feature-group weights are fit by Nelder-Mead optimisation against real observed engraftment correlation, not chosen by intuition.
- **The GA's fitness function was extended and *statistically validated*.** Adding a pairwise ecological-compatibility term produced a measurable, significant shift toward mutualistic consortia (Wilcoxon signed-rank on the compatibility score, *p* = 0.0045, n = 20 donors, paired/seeded comparison) — not just added on faith that it would help.
- **A causal invariance test** checks whether species' causal effects hold up across all five disease groups, not just whether the final classifier is accurate.

---

## Results at a glance

Held-out test set, Random Forest (selected over XGBoost and Logistic Regression via paired bootstrap comparison):

| Metric | Value |
|---|---|
| Pooled AUC-ROC | **0.777** |
| CD (n = 83) | 0.778 |
| UC (n = 227) | 0.675 |
| rCDI | 0.735 |
| MetS | 0.813 |
| MDR | 0.672 |
| ICI | 0.842 |
| IBD (pooled CD+UC) | 0.694 |
| CECE score vs. observed engraftment (Spearman) | *r* = 0.619 |

> The CD/UC gap above is real, investigated, and not yet fully closed — see [Limitations](#limitations-read-before-you-cite-a-number).

---

## Pipeline architecture



## Notebook guide

| # | Notebook | Environment | Purpose | Key outputs |
|---|---|---|---|---|
| 01 | `data_and_EDA` | 💻 Laptop | Load HMP2 + Podlesny + Smillie data; CLR transform, PCA, diversity, FDR-corrected differential abundance | `Fig1_EDA.png`, `*_network_CLR.csv` |
| 02 | `engraftment_labels` | 💻 Laptop | Derive donor-strain engraftment labels; case-level stratified train/test split (locked before labelling) | `engraftment_labels.csv`, `podlesny_wide.pkl` |
| 03a | `fastspar_prep` | 💻 Laptop | Validate/repair integer count matrices — CD / UC / Healthy | `*_network_counts_FIXED.tsv` |
| 03b | `fastspar_prep` | 💻 Laptop | Same, native to all 5 Podlesny disease groups | `{group}_pod_network_counts.tsv` |
| 04a | `fastspar_COLAB` | ☁️ Colab | Compile FastSpar; bootstrap correlations + *p*-values — CD / UC / Healthy | `cor_*.csv`, `pvals_*.csv` |
| 04b | `fastspar_COLAB` | ☁️ Colab | Same, across the 5 Podlesny groups | per-group `cor`/`pval` CSVs |
| 05 | `networks` | 💻 Laptop | Build co-occurrence networks + keystone species, HMP2 *and* native Podlesny groups | `*_network.gexf`, `keystones_*.csv` |
| 06 | `NOTEARS_causal_DAG` | ☁️ Colab | Learn a causal DAG over keystone species from HMP2 longitudinal data (CD-native + UC-native) | `dag_CD.csv`, `dag_UC.csv` (bootstrapped) |
| 07 | `causal_effects_CIs` | 💻 Laptop | Backdoor-adjusted average causal effect (ACE) per species; Bonferroni + FDR correction | `causal_effects.csv`, `causal_effects_UC.csv` |
| 08 | `micom_COLAB` | ☁️ Colab | Constraint-based (FBA) metabolic modelling of donor communities, all 5 disease groups | `micom_growth/exchange/scfa_results.csv` |
| 09 | `ML_SHAP_Ablation` | 💻 Laptop | Train Random Forest (vs. XGBoost / LogReg); real SHAP explainability; feature-group ablation | `cece_rf_model.pkl`, `Fig_SHAP.png` |
| 10 | `cece_score_and_GA` | 💻 Laptop | Composite CECE score; genetic algorithm for de novo consortium design + ecological-compatibility extension | `cece_final_scores.csv`, `ga_optimal_consortia.csv` |
| 11 | `CECE_additions` | 💻 Laptop | Calibrated CECE weighting, Donor-Recipient Compatibility Index (DRCI), SHAP-justified GA weights | `Table_GA_fitness_comparison.csv` |
| 12 | `validation` | 💻 Laptop | Held-out evaluation, per-disease-group breakdown, causal invariance test | `Table1_model_performance.csv`, `Fig_ROC.png` |

---

## Repository structure

```
cece-fmt/
├── notebooks/              # 01–12, run roughly top to bottom (see table above)
├── cece/                   # shared helpers imported by multiple notebooks
│   └── preprocessing.py     # preprocess() — used by NB01 and NB03b
├── data/
│   ├── raw/                 # NOT committed — see Data sources below
│   └── processed/           # NOT committed — regenerated by NB01–08
├── results/
│   ├── networks/            # co-occurrence networks, keystones (.gexf, .csv)
│   └── ...                  # causal effects, MICOM outputs, model artifacts
├── figures/
│   └── final/                # publication-ready figures
├── assets/
│   └── cece_banner.svg
├── environment.yml
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Getting started

```bash
git clone https://github.com/<your-username>/cece-fmt.git
cd cece-fmt

# conda environment (matches the notebooks' kernel name, cece-main)
conda create -n cece-main python=3.10
conda activate cece-main
pip install -r requirements.txt
```

Notebooks `04a`, `04b`, `06`, and `08` are written for **Google Colab** (they compile FastSpar from source / need Colab's compute for MICOM and NOTEARS) — upload them to Colab and mount Google Drive as directed in each notebook's first cell. Everything else runs locally.

Run order follows the numbering: `01 → 02 → {03a, 03b} → {04a, 04b} → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12`. The [pipeline architecture](#pipeline-architecture) diagram above shows which stages can run in parallel.

---

## Data sources

No raw data is committed to this repository — only code. All three source datasets are publicly available:

| Dataset | Used for | Source |
|---|---|---|
| HMP2 / IBDMDB | CD/UC/Healthy longitudinal profiles (causal DAG + ACE layer) | [ibdmdb.org](https://ibdmdb.org) — Lloyd-Price et al., 2019 |
| Podlesny et al. FMT cohort | Donor/recipient profiles across 5 disease indications | Podlesny et al., 2022, *Cell Reports Medicine* |
| Smillie et al. FMT cohort | Additional engraftment validation data | Smillie et al., 2018, *Cell Host & Microbe* |

To reproduce the pipeline, download these from their original sources into `data/raw/` following the instructions at the top of `01_data_and_EDA.ipynb`.

---

## Limitations (read before you cite a number)

- **The CD vs. UC AUC gap (0.778 vs. 0.675) is real and not fully explained.** Feature coverage, dominant-feature signal, and disease-group identity were all tested as explanations and ruled out or only partially adopted — documented in NB12 rather than hidden.
- **The causal layer (DAG + ACE) is HMP2-native only.** Podlesny's FMT samples are cross-sectional, not a genuine longitudinal series, so they aren't suitable for NOTEARS; causal effects for rCDI/MetS/MDR/ICI fall back to the nearest available HMP2 estimate rather than being disease-native.
- **FastSpar network edges use an uncorrected bootstrap *p*-value threshold** (standard in this literature for exploratory network construction, but not FDR-corrected — unlike the ACE estimation in NB07, which is).
- **Some UC-specific causal effect estimates rest on small per-species samples** (several below n=20); their bootstrap CIs are wide and are reported as such.
- **Everything here ran on a personal laptop**, with Google Colab used only for the four genuinely compute-heavy stages (FastSpar bootstrapping, NOTEARS structure learning, MICOM's flux balance simulations).

---

## Citation

```bibtex
@unpublished{cece2026,
  title  = {CECE: A Multi-Layer Computational Framework for Predictive
            Synthetic Faecal Microbiota Transplantation Design Integrating
            Causal Inference, Ecological Network Analysis, Constraint-Based
            Metabolic Modelling and Machine Learning},
  author = {Rachna},  % add your full name / affiliation
  year   = {2026},
  note   = {Undergraduate dissertation. Abstract submitted to BIOINFO/GIW,
            ISCB-Asia 2026, Seoul}
}
```

---

---

## License

Not yet chosen. If you intend this repository to be reused (e.g. alongside a grad-school application or a paper), a permissive license such as MIT is the common choice for academic code — see [choosealicense.com](https://choosealicense.com) — but that's your call to make, not a default to assume.

---

<p align="center"><em>Rachna · BSc Microbiology (Hons) · built on a laptop, not a cluster</em></p>
