import numpy as np
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix
import matplotlib.pyplot as plt
from sklearn.tree import plot_tree

# ─────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────
iris = load_iris()
X = iris.data          # 4 features: sepal/petal length & width
y = iris.target        # 3 classes: setosa, versicolor, virginica

print("Dataset shape:", X.shape)
print("Classes:", iris.target_names)

# ─────────────────────────────────────────
# 2. SPLIT DATA
# ─────────────────────────────────────────
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)
print(f"\nTrain size: {len(X_train)} | Test size: {len(X_test)}")

# ─────────────────────────────────────────
# 3. TRAIN DECISION TREE
# ─────────────────────────────────────────
# max_depth=3 keeps it simple and readable
model = DecisionTreeClassifier(max_depth=3, random_state=42)
model.fit(X_train, y_train)

print("\nModel trained!")

# ─────────────────────────────────────────
# 4. EVALUATE
# ─────────────────────────────────────────
y_pred = model.predict(X_test)
acc = accuracy_score(y_test, y_pred)
print(f"Accuracy: {acc:.2%}")

print("\nConfusion Matrix:")
cm = confusion_matrix(y_test, y_pred)
print(cm)
# Rows = actual class, Columns = predicted class
# Perfect model = numbers only on diagonal

# ─────────────────────────────────────────
# 5. SEE THE TREE RULES (human readable!)
# ─────────────────────────────────────────
print("\nDecision Tree Rules:")
print(export_text(model, feature_names=iris.feature_names))
# This prints exactly how the model makes decisions
# e.g. "if petal length <= 2.45 → setosa"

# ─────────────────────────────────────────
# 6. FEATURE IMPORTANCE
# ─────────────────────────────────────────
print("Feature Importances:")
for name, importance in zip(iris.feature_names, model.feature_importances_):
    bar = "█" * int(importance * 30)
    print(f"  {name:<26} {bar} {importance:.3f}")

# ─────────────────────────────────────────
# 7. VISUALIZE THE TREE
# ─────────────────────────────────────────
plt.figure(figsize=(14, 6))
plot_tree(
    model,
    feature_names=iris.feature_names,
    class_names=iris.target_names,
    filled=True,          # color nodes by class
    rounded=True,
    fontsize=10
)
plt.title("Decision Tree — Iris Dataset", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('day7_decision_tree.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nSaved: day7_decision_tree.png")

# ─────────────────────────────────────────
# 8. PREDICT A NEW FLOWER
# ─────────────────────────────────────────
new_flower = [[5.1, 3.5, 1.4, 0.2]]   # sepal_l, sepal_w, petal_l, petal_w
prediction = model.predict(new_flower)
proba = model.predict_proba(new_flower)

print(f"\nNew flower prediction: {iris.target_names[prediction[0]]}")
print(f"Confidence: {proba.max():.1%}")

print("\nDone! Concepts covered: train/test split, Decision Tree, confusion matrix, feature importance")
