# pip install xgboost matplotlib scikit-learn

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score
from xgboost import XGBClassifier

# ============================================================
# 1. LOAD & SPLIT DATA
# ============================================================

data = load_breast_cancer()
X, y = data.data, data.target
feature_names = data.feature_names

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

print(f"Train size: {X_train.shape}, Test size: {X_test.shape}")
print(f"Classes: {data.target_names}")  # ['malignant' 'benign']

# =============================================================
# 2. RANDOM FOREST
# =============================================================

print("\n" + "="*50)
print("RANDOM FOREST")
print("="*50)

rf = RandomForestClassifier(
    n_estimators=100,     # number of trees
    max_depth=None,       # grow trees fully
    max_features='sqrt',  # random subset of features per split
    random_state=42,
    n_jobs=-1             # use all CPU cores
)

rf.fit(X_train, y_train)

y_pred_rf = rf.predict(X_test)
y_proba_rf = rf.predict_proba(X_test)[:, 1]

print(classification_report(y_test, y_pred_rf, target_names=data.target_names))
print(f"RF ROC-AUC: {roc_auc_score(y_test, y_proba_rf):.3f}")

# =============================================================
# 3. FEATURE IMPORTANCE
# =============================================================

print("\n" + "="*50)
print("FEATURE IMPORTANCE (Top 10)")
print("="*50)

importances = pd.Series(rf.feature_importances_, index=feature_names)
importances = importances.sort_values(ascending=False)

print(importances.head(10).round(4))

# Plot
importances.head(10).plot(kind='bar', figsize=(10, 5), color='steelblue')
plt.title("Top 10 Feature Importances — Random Forest")
plt.ylabel("Importance Score")
plt.tight_layout()
plt.savefig("feature_importance.png", dpi=150)
plt.show()
print("Saved: feature_importance.png")

# =============================================================
# 4. XGBOOST
# =============================================================

print("\n" + "="*50)
print("XGBOOST")
print("="*50)

xgb = XGBClassifier(
    n_estimators=100,
    learning_rate=0.1,      # how much each tree contributes
    max_depth=4,
    subsample=0.8,          # fraction of samples per tree
    colsample_bytree=0.8,   # fraction of features per tree
    eval_metric='logloss',
    random_state=42,
    verbosity=0
)

xgb.fit(X_train, y_train)

y_pred_xgb = xgb.predict(X_test)
y_proba_xgb = xgb.predict_proba(X_test)[:, 1]

print(classification_report(y_test, y_pred_xgb, target_names=data.target_names))
print(f"XGBoost ROC-AUC: {roc_auc_score(y_test, y_proba_xgb):.3f}")

# =============================================================
# 5. SIDE-BY-SIDE COMPARISON
# =============================================================

print("\n" + "="*50)
print("MODEL COMPARISON")
print("="*50)

results = {
    "Random Forest": (y_pred_rf, y_proba_rf),
    "XGBoost":       (y_pred_xgb, y_proba_xgb),
}

print(f"{'Model':<20} {'Accuracy':<12} {'ROC-AUC'}")
print("-" * 42)
for name, (pred, proba) in results.items():
    acc = accuracy_score(y_test, pred)
    auc = roc_auc_score(y_test, proba)
    print(f"{name:<20} {acc:<12.3f} {auc:.3f}")

# =============================================================
# 6. KEY TAKEAWAYS (printed as reminder)
# =============================================================

print("""
KEY TAKEAWAYS:
- Random Forest   → many trees in parallel, vote → low variance
- XGBoost         → trees in sequence, each fixes errors → high accuracy
- Feature Import. → which biomarkers drive prediction (great for bio/health ML)
- Use ROC-AUC     → better than accuracy for imbalanced medical datasets

GitHub commit → Day07_Ensemble_Methods/
""")
