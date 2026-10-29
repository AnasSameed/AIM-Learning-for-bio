"""
Toy asynchronous Boolean network simulator (learning example).
Illustrates the method used in my FOXP3+ Treg thesis: nodes are ON/OFF,
logical rules, random-order asynchronous updates, knockouts, heatmap.
NOTE: simplified 8-node toy model, NOT the full 135-node thesis network.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# node: rule(state) -> bool
RULES = {
    "IL2":    lambda s: s["IL2"],                                  # external input, held
    "IL6":    lambda s: s["IL6"],                                  # external input, held
    "STAT5":  lambda s: s["IL2"] and not s["STAT3"],
    "STAT3":  lambda s: s["IL6"],
    "PI3K":   lambda s: s["IL2"] and not s["PTEN"],
    "PTEN":   lambda s: not s["AKT"],
    "AKT":    lambda s: s["PI3K"],
    "FOXP3":  lambda s: s["STAT5"] and not s["AKT"] and not s["STAT3"],
}
NODES = list(RULES)

def simulate(init, steps=500, ko=(), seed=0):
    rng = np.random.default_rng(seed)
    s = dict(init)
    for k in ko: s[k] = False
    trace = []
    for _ in range(steps):
        n = NODES[rng.integers(len(NODES))]          # asynchronous: one random node per step
        if n not in ko:
            s[n] = bool(RULES[n](s))
        trace.append(s.copy())
    return pd.DataFrame(trace).astype(int)

def mean_foxp3(init, ko=(), reps=50, steps=500):
    # average FOXP3 over the last 100 steps and over replicates
    return np.mean([simulate(init, steps, ko, seed=r)["FOXP3"].iloc[-100:].mean() for r in range(reps)])

healthy  = {n: False for n in NODES} | {"IL2": True,  "IL6": False}
diseased = {n: False for n in NODES} | {"IL2": True,  "IL6": True}

rows = []
for cond, init in [("Healthy (IL2)", healthy), ("Inflamed (IL2+IL6)", diseased)]:
    for ko in [()] + [(n,) for n in NODES if n not in ("IL2", "IL6")]:
        rows.append({"condition": cond, "knockout": ko[0] if ko else "none",
                     "FOXP3_activity": mean_foxp3(init, ko)})
df = pd.DataFrame(rows)
print(df.pivot(index="knockout", columns="condition", values="FOXP3_activity").round(2))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
tr = simulate(healthy, 200)
ax[0].plot(tr["FOXP3"].rolling(10).mean(), label="Healthy")
ax[0].plot(simulate(diseased, 200)["FOXP3"].rolling(10).mean(), label="Inflamed")
ax[0].set(title="FOXP3 activity over time (asynchronous)", xlabel="step", ylabel="rolling mean ON")
ax[0].legend()
sns.heatmap(df.pivot(index="knockout", columns="condition", values="FOXP3_activity"),
            annot=True, cmap="viridis", ax=ax[1])
ax[1].set_title("Knockout perturbation matrix (FOXP3)")
plt.tight_layout(); plt.savefig("boolean_network_toy.png", dpi=150)
