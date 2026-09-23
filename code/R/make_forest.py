import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

models = [
    "1. Crude\n(unadjusted)",
    "2. DAG-adjusted\nregression",
    "3. IPTW\nonly",
    "4. IPTW + adjusted\n(doubly robust)",
    "5. PS-quintile stratified\n(Mantel-Haenszel)",
]
OR   = [1.783, 1.458, 1.293, 1.297, 1.414]
lo   = [1.750, 1.378, 1.224, 1.225, 1.387]
hi   = [1.816, 1.542, 1.366, 1.372, 1.442]

fig, ax = plt.subplots(figsize=(11, 5.6))
y = np.arange(len(models))[::-1]

colors = ["#8a8f98", "#8a8f98", "#c1611a", "#c1611a", "#1f3864"]
for i, yi in enumerate(y):
    ax.plot([lo[i], hi[i]], [yi, yi], color=colors[i], lw=3, solid_capstyle="round", zorder=2)
    ax.scatter([OR[i]], [yi], color=colors[i], s=170, zorder=3, edgecolor="white", linewidth=1.2)
    ax.text(hi[i] + 0.025, yi, f"{OR[i]:.2f}  ({lo[i]:.2f}–{hi[i]:.2f})", va="center", fontsize=13, color="#1a2233", fontweight="bold")

ax.axvline(1.0, color="#5b6472", lw=1.3, linestyle="--", zorder=1)
ax.text(1.0, len(models)-0.35, "no association", ha="center", fontsize=10.5, color="#5b6472", style="italic")

ax.set_yticks(y)
ax.set_yticklabels(models, fontsize=13, color="#1a2233")
ax.set_xlim(0.95, 2.05)
ax.set_xlabel("Odds ratio for preterm birth (smokers vs. non-smokers)", fontsize=12.5, color="#1a2233")
ax.spines[["top","right","left"]].set_visible(False)
ax.spines["bottom"].set_color("#c9ccd1")
ax.tick_params(axis="x", colors="#5b6472", labelsize=11)
ax.set_ylim(-0.6, len(models)-0.4)
ax.grid(axis="x", color="#e7e9ec", lw=0.8, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig("fig_forest.png", dpi=230, bbox_inches="tight", facecolor="#FAF8F4")
print("saved")
