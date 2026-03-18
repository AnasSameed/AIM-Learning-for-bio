# 🌿 Plant Species — Exploratory Data Analysis

A beginner-to-intermediate Python data analysis project exploring morphological
patterns in plant species data — built as part of a biology-to-data-science portfolio.

---

## 📌 Project Overview

This project performs a full Exploratory Data Analysis (EDA) on plant species data,
simulating a real-world botanical field study. It covers species distribution,
trait comparisons, correlation analysis, and seasonal collection patterns.

**This is exactly the kind of analysis used in:**
- AgriTech companies studying crop traits
- Botanical research institutes
- Life sciences data analyst roles

---

## 📊 What the Analysis Covers

| Analysis | Description |
|----------|-------------|
| Species distribution | Which species were sampled most? |
| Sepal & petal morphology | How do measurements differ across species? |
| Scatter plots | Relationship between petal length and width |
| Boxplots | Variability in sepal width per species |
| Correlation heatmap | Which traits are most related? |
| Regional sampling | Where were most plants collected? |
| Monthly trend | When was data collected across the year? |
| Mean trait comparison | Side-by-side species comparison |

---

## 🗂️ Files

```
plant-species-eda/
├── plant_species_eda.py       ← Python script (run directly)
├── plant_species_eda.ipynb    ← Jupyter Notebook (recommended)
├── plant_species_eda.png      ← Output chart (auto-generated)
├── requirements.txt           ← Dependencies
└── README.md                  ← This file
```

---

## ⚙️ How to Run

### Option A — Jupyter Notebook (recommended)
```bash
pip install -r requirements.txt
jupyter notebook plant_species_eda.ipynb
```

### Option B — Python script
```bash
pip install -r requirements.txt
python plant_species_eda.py
```

---

## 📦 Requirements

```
pandas
numpy
matplotlib
scikit-learn
jupyter
pytorch
KNN

## 🔍 Key Findings

1. **Iris setosa** is clearly separable from other species due to much smaller petal dimensions
2. **Petal length ↔ petal width** have the strongest correlation (r ≈ 0.96)
3. **Sepal width** shows the most overlap across species — least useful for classification
4. **Iris virginica** has the largest average petal measurements

---

## 🧠 What I Learned

- Data cleaning and inspection with `pandas`
- Creating multi-panel visualisations with `matplotlib`
- Understanding correlation and what it means biologically
- Structuring a data analysis project clearly for others to read

---

## 🚀 Next Steps

This EDA is **Week 1** of a 30-day portfolio building project. Next:
- **Week 2:** Crop yield analysis with Indian agricultural data
- **Week 3:** Machine learning — plant disease classifier
- **Week 4:** Bioinformatics — DNA sequence analysis

---

## 👤 Author

**Anas Sameed**  
MSc Botany | Transitioning into Data Science & Bioinformatics  
📍 Bengaluru, India

---

*Dataset: Iris (public domain, UCI ML Repository)*  
*All simulated field study columns (region, month) are randomly generated for demonstration.*
