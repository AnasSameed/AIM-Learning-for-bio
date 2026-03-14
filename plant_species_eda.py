# ============================================================
#  Plant Species EDA — Exploratory Data Analysis
#  Author  : [Your Name]
#  Dataset : GBIF / Kaggle Iris + custom plant species data
#  Goal    : Understand species distribution, traits, and
#            patterns in plant data using Python
# ============================================================

# ---------- 1. IMPORTS ----------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings("ignore")

# Set a clean style for all plots
plt.rcParams.update({
    "figure.facecolor": "#FAFAF8",
    "axes.facecolor":   "#FAFAF8",
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.3,
    "grid.linestyle":   "--",
    "font.family":      "DejaVu Sans",
    "font.size":        11,
})

COLORS = ["#2D6A4F", "#52B788", "#95D5B2", "#D8F3DC",
          "#1B4332", "#40916C", "#74C69D", "#B7E4C7"]

print("=" * 55)
print("  Plant Species EDA — Starting Analysis")
print("=" * 55)


# ---------- 2. LOAD DATA ----------
# We use the classic Iris dataset (public domain, no download needed)
# and extend it to simulate a real plant species study.

from sklearn.datasets import load_iris  # pre-installed, always available

def load_plant_data():
    """Load Iris dataset and rename columns for botanical clarity."""
    iris = load_iris(as_frame=True)
    df = iris.frame.copy()

    # Rename to be more botanically meaningful
    df.columns = [
        "sepal_length_cm", "sepal_width_cm",
        "petal_length_cm", "petal_width_cm",
        "species_code"
    ]

    # Map numeric codes to actual species names
    species_map = {
        0: "Iris setosa",
        1: "Iris versicolor",
        2: "Iris virginica"
    }
    df["species"] = df["species_code"].map(species_map)

    # Add a simulated 'region' column (realistic for a field study)
    np.random.seed(42)
    regions = ["Northern Plains", "Eastern Ghats", "Western Coast"]
    df["region"] = np.random.choice(regions, size=len(df),
                                     p=[0.4, 0.35, 0.25])

    # Add a simulated 'month_collected' column
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    df["month_collected"] = np.random.choice(months, size=len(df))

    df.drop(columns=["species_code"], inplace=True)
    return df


df = load_plant_data()
print(f"\n✔  Dataset loaded: {df.shape[0]} samples, {df.shape[1]} columns\n")


# ---------- 3. BASIC INSPECTION ----------
print("--- Shape ---")
print(f"  Rows: {df.shape[0]}  |  Columns: {df.shape[1]}\n")

print("--- Column types ---")
print(df.dtypes.to_string(), "\n")

print("--- First 5 rows ---")
print(df.head().to_string(), "\n")

print("--- Missing values ---")
print(df.isnull().sum().to_string(), "\n")

print("--- Summary statistics ---")
print(df.describe().round(2).to_string(), "\n")

print("--- Species counts ---")
print(df["species"].value_counts().to_string(), "\n")


# ---------- 4. VISUALISATION ----------

fig = plt.figure(figsize=(18, 14))
fig.suptitle("Plant Species — Exploratory Data Analysis",
             fontsize=20, fontweight="bold", y=0.98, color="#1B4332")
gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)


# ── Plot 1: Species distribution (bar chart) ──────────────────
ax1 = fig.add_subplot(gs[0, 0])
counts = df["species"].value_counts()
bars = ax1.bar(counts.index, counts.values,
               color=COLORS[:3], edgecolor="white", linewidth=0.8)
ax1.set_title("Species Distribution", fontweight="bold", color="#1B4332")
ax1.set_xlabel("Species")
ax1.set_ylabel("Count")
ax1.set_xticklabels(counts.index, rotation=15, ha="right", fontsize=9)
for bar, val in zip(bars, counts.values):
    ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
             str(val), ha="center", va="bottom", fontsize=10, fontweight="bold")


# ── Plot 2: Sepal length distribution per species (histogram) ──
ax2 = fig.add_subplot(gs[0, 1])
for i, (sp, grp) in enumerate(df.groupby("species")):
    ax2.hist(grp["sepal_length_cm"], bins=12, alpha=0.7,
             label=sp, color=COLORS[i], edgecolor="white")
ax2.set_title("Sepal Length Distribution", fontweight="bold", color="#1B4332")
ax2.set_xlabel("Sepal Length (cm)")
ax2.set_ylabel("Frequency")
ax2.legend(fontsize=8)


# ── Plot 3: Petal length vs petal width (scatter) ─────────────
ax3 = fig.add_subplot(gs[0, 2])
for i, (sp, grp) in enumerate(df.groupby("species")):
    ax3.scatter(grp["petal_length_cm"], grp["petal_width_cm"],
                label=sp, color=COLORS[i], alpha=0.75, s=40, edgecolors="white")
ax3.set_title("Petal Length vs Width", fontweight="bold", color="#1B4332")
ax3.set_xlabel("Petal Length (cm)")
ax3.set_ylabel("Petal Width (cm)")
ax3.legend(fontsize=8)


# ── Plot 4: Boxplot — sepal width by species ──────────────────
ax4 = fig.add_subplot(gs[1, 0])
species_list = df["species"].unique()
data_by_species = [df[df["species"] == sp]["sepal_width_cm"].values
                   for sp in species_list]
bp = ax4.boxplot(data_by_species, patch_artist=True, notch=False,
                 medianprops=dict(color="#1B4332", linewidth=2))
for patch, color in zip(bp["boxes"], COLORS):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax4.set_title("Sepal Width by Species", fontweight="bold", color="#1B4332")
ax4.set_ylabel("Sepal Width (cm)")
ax4.set_xticklabels([sp.split()[1] for sp in species_list],
                    rotation=10, fontsize=9)


# ── Plot 5: Collection by region (horizontal bar) ─────────────
ax5 = fig.add_subplot(gs[1, 1])
region_counts = df["region"].value_counts()
ax5.barh(region_counts.index, region_counts.values,
         color=COLORS[1:4], edgecolor="white")
ax5.set_title("Samples by Region", fontweight="bold", color="#1B4332")
ax5.set_xlabel("Count")
for i, v in enumerate(region_counts.values):
    ax5.text(v + 0.3, i, str(v), va="center", fontsize=10, fontweight="bold")


# ── Plot 6: Correlation heatmap ───────────────────────────────
ax6 = fig.add_subplot(gs[1, 2])
numeric_cols = ["sepal_length_cm", "sepal_width_cm",
                "petal_length_cm", "petal_width_cm"]
corr = df[numeric_cols].corr()
short_names = ["Sep.L", "Sep.W", "Pet.L", "Pet.W"]
im = ax6.imshow(corr.values, cmap="Greens", vmin=-1, vmax=1)
ax6.set_xticks(range(4))
ax6.set_yticks(range(4))
ax6.set_xticklabels(short_names, fontsize=9)
ax6.set_yticklabels(short_names, fontsize=9)
ax6.set_title("Correlation Heatmap", fontweight="bold", color="#1B4332")
for i in range(4):
    for j in range(4):
        ax6.text(j, i, f"{corr.values[i, j]:.2f}",
                 ha="center", va="center", fontsize=9,
                 color="white" if abs(corr.values[i, j]) > 0.6 else "#1B4332")
plt.colorbar(im, ax=ax6, shrink=0.8)


# ── Plot 7: Monthly collection trend ──────────────────────────
ax7 = fig.add_subplot(gs[2, :2])
month_order = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]
monthly = df.groupby(["month_collected", "species"]).size().unstack(fill_value=0)
monthly = monthly.reindex(month_order)
bottom = np.zeros(12)
for i, col in enumerate(monthly.columns):
    ax7.bar(monthly.index, monthly[col], bottom=bottom,
            label=col, color=COLORS[i], alpha=0.85, edgecolor="white")
    bottom += monthly[col].values
ax7.set_title("Collection Activity by Month & Species",
              fontweight="bold", color="#1B4332")
ax7.set_xlabel("Month")
ax7.set_ylabel("Samples Collected")
ax7.legend(loc="upper right", fontsize=9)


# ── Plot 8: Mean measurements per species (grouped bar) ───────
ax8 = fig.add_subplot(gs[2, 2])
means = df.groupby("species")[numeric_cols].mean()
x = np.arange(len(numeric_cols))
width = 0.25
short_labels = ["Sep.L", "Sep.W", "Pet.L", "Pet.W"]
for i, (sp, row) in enumerate(means.iterrows()):
    ax8.bar(x + i * width, row.values, width,
            label=sp.split()[1], color=COLORS[i], alpha=0.85, edgecolor="white")
ax8.set_title("Mean Measurements by Species",
              fontweight="bold", color="#1B4332")
ax8.set_ylabel("cm")
ax8.set_xticks(x + width)
ax8.set_xticklabels(short_labels, fontsize=9)
ax8.legend(fontsize=8)


plt.savefig("plant_species_eda.png", dpi=150, bbox_inches="tight",
            facecolor="#FAFAF8")
print("✔  Plot saved → plant_species_eda.png")
plt.show()


# ---------- 5. KEY FINDINGS ----------
print("\n" + "=" * 55)
print("  KEY FINDINGS")
print("=" * 55)

# Largest species
top_species = df["species"].value_counts().idxmax()
print(f"\n1. Most sampled species : {top_species}")

# Highest mean petal length
mean_petal = df.groupby("species")["petal_length_cm"].mean()
longest_petal = mean_petal.idxmax()
print(f"2. Longest petals       : {longest_petal} "
      f"({mean_petal[longest_petal]:.2f} cm avg)")

# Best correlated pair
corr_vals = corr.abs().unstack()
corr_vals = corr_vals[corr_vals < 1].sort_values(ascending=False)
best_pair  = corr_vals.index[0]
best_score = corr_vals.iloc[0]
print(f"3. Strongest correlation: {best_pair[0]} ↔ {best_pair[1]} "
      f"(r = {best_score:.2f})")

# Top region
top_region = df["region"].value_counts().idxmax()
print(f"4. Most samples from    : {top_region}")

print("\n✔  Analysis complete!")
print("=" * 55)
