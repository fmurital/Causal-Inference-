import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

fig, ax = plt.subplots(figsize=(13, 8.2))
ax.set_xlim(0, 12)
ax.set_ylim(-1.3, 6.6)
ax.axis("off")

# ---- exposure and outcome nodes ----
T = (1.4, 0.2)
Y = (10.6, 0.2)

confounders = [
    "Maternal age",
    "Race / ethnicity",
    "Education",
    "Prenatal care\ntiming",
    "WIC\nparticipation",
    "Parity\n(nulliparous)",
    "Prior preterm\nbirth",
    "Pre-pregnancy\ndiabetes",
    "Pre-pregnancy\nhypertension",
    "Pre-pregnancy\nBMI",
    "Insurance /\npayer type",
]
n = len(confounders)
xs = np.linspace(2.4, 9.6, n)
ys = [5.6 if i % 2 == 0 else 4.85 for i in range(n)]

def node(ax, xy, text, w, h, fc, ec, fontsize=10.5, fontweight="normal", zorder=3):
    x, y = xy
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle="round,pad=0.02,rounding_size=0.08",
                          linewidth=1.6, edgecolor=ec, facecolor=fc, zorder=zorder)
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
             fontweight=fontweight, zorder=zorder+1, color="#1a1a1a")

def arrow(ax, p1, p2, color="#8a8a8a", lw=1.15, alpha=0.85, rad=0.0, zorder=1):
    a = FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=13,
                         linewidth=lw, color=color, alpha=alpha,
                         connectionstyle=f"arc3,rad={rad}", zorder=zorder,
                         shrinkA=10, shrinkB=10)
    ax.add_patch(a)

# confounder boxes + arrows to T and Y
cw, ch = 1.55, 0.62
for (x, y), label in zip(zip(xs, ys), confounders):
    node(ax, (x, y), label, cw, ch, fc="#eef2f7", ec="#6b7f99", fontsize=9.3)
    rad_t = 0.06 if y > 5.2 else 0.03
    rad_y = -0.06 if y > 5.2 else -0.03
    arrow(ax, (x - cw/4, y - ch/2), (T[0] + 0.35, T[1] + 0.55), color="#9aa5b1", lw=0.9, rad=rad_t*3)
    arrow(ax, (x + cw/4, y - ch/2), (Y[0] - 0.35, Y[1] + 0.55), color="#9aa5b1", lw=0.9, rad=-rad_t*3)

# exposure and outcome nodes (drawn after so arrows sit behind them)
node(ax, T, "Maternal smoking\nduring pregnancy\n(exposure)", 2.15, 1.05, fc="#fde3d0", ec="#c1611a", fontsize=11, fontweight="bold", zorder=5)
node(ax, Y, "Preterm birth\n< 37 weeks\n(outcome)", 2.15, 1.05, fc="#d9ecdd", ec="#237a3e", fontsize=11, fontweight="bold", zorder=5)

# direct causal-effect arrow
arrow(ax, (T[0] + 1.08, T[1] - 0.15), (Y[0] - 1.08, Y[1] - 0.15),
      color="#c1611a", lw=2.6, alpha=1.0, rad=-0.18, zorder=4)
ax.text((T[0]+Y[0])/2, -0.95, "Causal effect of interest\n(odds ratio, adjusted / weighted)",
        ha="center", va="center", fontsize=10.5, color="#c1611a", fontweight="bold")

# legend / caption
ax.text(6.0, 6.25,
        "Directed acyclic graph: maternal smoking and preterm birth (2023 NVSS singleton births)",
        ha="center", va="center", fontsize=13.5, fontweight="bold", color="#20242b")
ax.text(6.0, -1.2,
        "Gray nodes: measured pre-treatment confounders (common causes of exposure and outcome) — adjusted for via the propensity score.\n"
        "Excluded from adjustment: gestational weight gain and other variables that occur after smoking status is determined (mediators / potential colliders).",
        ha="center", va="center", fontsize=8.8, color="#555555", style="italic")

plt.tight_layout()
plt.savefig("fig_dag.png", dpi=230, bbox_inches="tight", facecolor="white")
print("saved fig_dag.png")
