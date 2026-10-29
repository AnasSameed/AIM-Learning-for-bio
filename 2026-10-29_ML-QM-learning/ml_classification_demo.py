import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, ConfusionMatrixDisplay

USE_REAL_DATA = False
if USE_REAL_DATA:
    X = pd.read_csv("vst_matrix.csv", index_col=0)
    lab = pd.read_csv("labels.csv", index_col=0)
    y = (lab.loc[X.index, "group"] == "Primary Tumor").astype(int).values
    tag = "TCGA-PRAD"
else:
    rng = np.random.default_rng(42)
    y = np.array([1]*50 + [0]*10)                        # 50 tumor, 10 normal
    Xa = rng.normal(8, 2, size=(60, 2000))               # log-expression-like noise
    Xa[y == 1, :60] += 1.5                               # 60 genes up in tumor
    Xa[y == 1, 60:120] -= 1.5                            # 60 genes down in tumor
    X = pd.DataFrame(Xa, columns=[f"gene_{i}" for i in range(2000)])
    tag = "SIMULATED data"

cv = StratifiedKFold(5, shuffle=True, random_state=42)
k = min(500, X.shape[1])
models = {
 "LASSO (L1-logistic)": Pipeline([("s",StandardScaler()),("f",SelectKBest(f_classif,k=k)),
    ("c",LogisticRegression(penalty="l1",solver="liblinear",C=0.5,class_weight="balanced"))]),
 "Random forest": Pipeline([("f",SelectKBest(f_classif,k=k)),
    ("c",RandomForestClassifier(n_estimators=500,class_weight="balanced",random_state=42))]),
}
fig, ax = plt.subplots(2, 2, figsize=(12, 10))
# 1 PCA
Z = PCA(2).fit_transform(StandardScaler().fit_transform(X))
for v,l,c in [(1,"Tumor","tab:red"),(0,"Normal","tab:blue")]:
    ax[0,0].scatter(Z[y==v,0],Z[y==v,1],c=c,label=l,edgecolor="k")
ax[0,0].set_title(f"PCA ({tag})"); ax[0,0].set_xlabel("PC1"); ax[0,0].set_ylabel("PC2"); ax[0,0].legend()
# 2 ROC
preds = {}
for n,m in models.items():
    p = cross_val_predict(m,X,y,cv=cv,method="predict_proba")[:,1]; preds[n]=p
    fpr,tpr,_ = roc_curve(y,p); ax[0,1].plot(fpr,tpr,label=f"{n} AUC={roc_auc_score(y,p):.2f}")
ax[0,1].plot([0,1],[0,1],"k--"); ax[0,1].set_title("5-fold CV ROC"); ax[0,1].set_xlabel("FPR"); ax[0,1].set_ylabel("TPR"); ax[0,1].legend()
# 3 confusion matrix (LASSO)
pl = (preds["LASSO (L1-logistic)"]>0.5).astype(int)
ConfusionMatrixDisplay(confusion_matrix(y,pl),display_labels=["Normal","Tumor"]).plot(ax=ax[1,0],colorbar=False)
ax[1,0].set_title("LASSO confusion matrix (CV predictions)")
# 4 top genes
m = models["LASSO (L1-logistic)"].fit(X,y)
genes = X.columns[m.named_steps["f"].get_support()]
coef = pd.Series(m.named_steps["c"].coef_[0],index=genes)
top = coef.reindex(coef.abs().sort_values(ascending=False).head(10).index)[::-1]
ax[1,1].barh(top.index,top.values,color=["tab:red" if v>0 else "tab:blue" for v in top.values])
ax[1,1].set_title("Top 10 LASSO genes (red = tumor-associated)")
plt.tight_layout(); plt.savefig("ml_summary.png",dpi=150)
print({n:round(roc_auc_score(y,p),3) for n,p in preds.items()})
