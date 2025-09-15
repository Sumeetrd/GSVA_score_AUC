***AIM*** <br>
The goal is to learn a single, biologically sensible and reproducible cutoff on a disease-severity GSVA score that separates inflamed vs. uninflamed biology and, when possible, translates that into Responder vs. Non-responder calls for treated samples. The pipeline prefers supervised evidence (ROCs from curated contrasts) and gracefully falls back to an unsupervised valley/median if labels are thin while leaving a full paper trail of decisions and plots.

<br>**Example of input metadata file:** <br>
| Sample.ID | Condition                    | PatientID | Timepoint |
| --------- | ---------------------------- | --------- | --------- |
| S1        | Treated\_CD\_Uninflamed      | P001      | Baseline  |
| S2        | Treated\_CD\_Inflamed        | P002      | Week6     |
| S3        | Treated\_Responder           | P003      | Week8     |
| S4        | Treated\_non\_responders     | P004      | Week8     |
| S5        | Not\_treated\_CD\_Uninflamed | P005      | Baseline  |
| S6        | Not\_treated\_CD\_Inflamed   | P006      | Baseline  |
| S7        | Treated\_CD\_Uninflamed      | P007      | Week6     |
| S8        | Treated\_CD\_Inflamed        | P008      | Week6     |

**GSVA SCORE file:**
| pathway                       | S1    | S2   | S3    | S4   | S5    | S6   | S7    | S8   |
| ----------------------------- | ----- | ---- | ----- | ---- | ----- | ---- | ----- | ---- |
| Pathway_1 | -0.15 | 0.22 | -0.18 | 0.30 | -0.25 | 0.27 | -0.12 | 0.24 |
| Another\_pathway              | 0.05  | 0.08 | 0.01  | 0.11 | 0.02  | 0.10 | 0.04  | 0.09 |


# Disease Severity Cutoff (GSVA) – Colon (Usti) Example

**What it does**
- Loads metadata (`Sample.ID`, `Condition`) and GSVA scores for a chosen pathway.
- Learns a **cutoff** using best supervised contrast (Youden on ROC); falls back to **unsupervised** valley/median.
- Writes ranked contrasts + publishes **density/ROC** plots per contrast.
- Produces an **all-groups density** with arrows pointing to each group’s peak and the learned cutoff.

**Run**
```bash
Rscript scripts/disease_severity_cutoff.R \
  # (or use env vars to override defaults)
  # META_PATH=data/metadata.csv \
  # SCORES_PATH=data/GSVA_Score.csv \
  # PATHWAY_NAME=Pathway_1 \
  # OUTDIR=outputs

