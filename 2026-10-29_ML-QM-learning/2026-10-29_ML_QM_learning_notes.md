# Learning notes, 29 Oct 2026: ML, QM and computational modeling

## 1. Supervised ML pipeline (`ml_classification_demo.py`)
Binary classifier (tumor = 1, normal = 0) on 60 samples x 2,000 genes (**simulated** expression data; set `USE_REAL_DATA = True` for a real VST matrix).

| Step | Purpose |
|---|---|
| StandardScaler | Same scale for every gene (needed for L1 penalty) |
| SelectKBest(f_classif, k=500) | Keep the 500 most informative genes |
| LogisticRegression(penalty="l1", C=0.5) | LASSO: sparse, interpretable gene signature |
| RandomForestClassifier(500 trees) | Non-linear, robust to noise |
| class_weight="balanced" | Corrects the 50:10 class imbalance |
| StratifiedKFold(5) | Honest performance estimate, same class ratio per fold |

- Feature selection sits **inside** the pipeline so it is redone per fold (no data leakage).
- Simulated results: AUC LASSO 0.81, random forest 0.93. Confusion matrix (LASSO): sensitivity 92%, specificity 60%.
- Accuracy alone misleads with imbalance (always predicting "tumor" gives 83%), so also report AUC, sensitivity, specificity.
- With 60 samples, results are noisy. Repeat CV and validate on an independent cohort.

Concepts: overfitting, L1 vs L2, logistic regression (classification despite its name), ROC/AUC, precision/recall, cross-validation, hyperparameter tuning with nested CV.

## 2. Molecular docking and MD
- Docking score is an estimated binding free energy in **kcal/mol**. More negative = tighter predicted binding. About 1.4 kcal/mol per 10x affinity.
- Limits: rigid receptor, approximate scoring, ignores solvent dynamics, not comparable across proteins.
- Validation: redock the co-crystal ligand (RMSD under 2 A), dock controls, run MD, estimate MM-PBSA, then SPR/ITC.
- Lipinski's rule checks the ligand's drug-likeness, not binding.
- MD workflow: force field (AMBER/CHARMM, GAFF/CGenFF for ligand), solvate (TIP3P), add ions, minimize, NVT, NPT, production.
- MD metrics: RMSD (stability), RMSF (per-residue flexibility), Rg (compactness), H-bond persistence, MM-PBSA/GBSA (binding free energy).

## 3. Quantum mechanics basics
- MM treats atoms as springs; QM solves for electrons. DFT uses electron density.
- Choices: functional (B3LYP, M06-2X), basis set (def2-SVP, def2-TZVP), D3 dispersion, implicit solvent (PCM/SMD).
- Minimum = 0 imaginary frequencies. Transition state = 1 imaginary frequency, confirmed by IRC.
- ML interatomic potentials (MLIPs) are trained on QM energies/forces to give near-DFT accuracy at low cost.
- Tools: ORCA/Psi4/PySCF compute; Avogadro/Chemcraft build and visualize.
- ORCA input example:
```
! B3LYP D3BJ def2-SVP Opt Freq
* xyz 0 1
O  0.000  0.000  0.000
H  0.000  0.760  0.590
H  0.000 -0.760  0.590
*
```
Run with `orca water.inp > water.out`. Check for "ORCA TERMINATED NORMALLY" and zero imaginary frequencies.

## 4. Bioinformatics terms
FDR = false discovery rate (Benjamini-Hochberg). RMA = robust multi-array average (background correction, quantile normalization, summarization). limma = linear models for microarray data. DESeq2 = negative binomial model for RNA-seq counts.

## 5. GenAI and LLM concepts
LLM = transformer predicting next tokens; can hallucinate, so verify output. RAG = retrieve documents to ground answers. Agentic workflow = LLM plans, calls tools, reads results, decides next step.

## Files in this folder
- `ml_classification_demo.py`: pipeline plus four-panel plot (PCA, ROC, confusion matrix, top LASSO genes)
- `ml_summary_simulated.png`: output of the script on simulated data
- `lilly_interview_cheat_sheet.png`: one-page revision sheet
