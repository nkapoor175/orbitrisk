"""
Purpose : Extra figures for the IEEE report: pipeline diagram, simulated rejection rates of the
          Wilcoxon test, and the distribution of exact Gaussian Pc.
Reads   : data/processed/type_error_sim.csv, data/processed/pc_exact.csv
Writes  : report/figures/report_pipeline.png, report_type_error.png, report_pc_hist.png
Status  : DONE
Run from the repo root:  python report/make_report_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 9})

# 1. pipeline
steps = ["Space-Track\nTLE history\n(300 objects)", "Clean\n(manoeuvre\nflags)", "SGP4 residuals\nR / T / N\n(145,467)",
         "Error model\nσ(t)=a·t^b\nGaussian / t", "Screening\nk-d tree\n(8,172 conj.)", "Pc: exact,\nMC, IS +\nBeta intervals"]
fig, ax = plt.subplots(figsize=(7.0, 1.5))
w, h, gap = 1.0, 0.9, 0.22
for i, s in enumerate(steps):
    x = i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, 0), w, h, boxstyle="round,pad=0.02", fc="#e8eef6", ec="#14213d", lw=0.8))
    ax.text(x + w / 2, h / 2, s, ha="center", va="center", fontsize=7.5)
    if i < len(steps) - 1:
        ax.annotate("", xy=(x + w + gap - 0.02, h / 2), xytext=(x + w + 0.02, h / 2),
                    arrowprops=dict(arrowstyle="->", color="#14213d", lw=0.9))
ax.set_xlim(-0.05, len(steps) * (w + gap) - gap + 0.05)
ax.set_ylim(-0.05, h + 0.05)
ax.axis("off")
fig.savefig(OUT / "report_pipeline.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# 2. rejection rates
sim = pd.read_csv(PROC / "type_error_sim.csv")
fig, ax = plt.subplots(figsize=(3.4, 2.4))
colors = {"radial": "#2f6fb0", "along": "#c8641e", "cross": "#4a7d4a"}
for comp, col in colors.items():
    for truth, ls, lab in [("gaussian", "--", "Type I"), ("student_t", "-", "power")]:
        s = sim[(sim.direction == comp) & (sim.truth == truth)].sort_values("n")
        ax.plot(s.n, s.reject_rate, ls, color=col, marker="o", ms=3, lw=1, label=f"{lab}, {comp}")
ax.axhline(0.05, color="black", lw=0.7, ls=":")
ax.text(1000, 0.06, "nominal α = 0.05", fontsize=7, ha="right", va="bottom")
ax.set_xscale("log")
ax.set_xlabel("sample size n")
ax.set_ylabel("rejection rate")
ax.set_ylim(0, 1.02)
ax.legend(fontsize=6, ncol=2, loc="center right")
fig.savefig(OUT / "report_type_error.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# 3. Pc distribution
pc = pd.read_csv(PROC / "pc_exact.csv").pc
pos = pc[pc > 0]
fig, ax = plt.subplots(figsize=(3.4, 2.3))
ax.hist(np.log10(pos), bins=np.arange(-60, -1, 1), color="#2f6fb0", alpha=0.85)
for v, lab in [(-6, "10⁻⁶"), (-5, "10⁻⁵"), (-4, "10⁻⁴")]:
    ax.axvline(v, color="#c8641e", lw=0.8, ls="--")
ax.set_xlabel("log₁₀ Pc (exact Gaussian)")
ax.set_ylabel("conjunctions")
ax.set_xlim(-60, -2)
fig.savefig(OUT / "report_pc_hist.png", dpi=300, bbox_inches="tight")
plt.close(fig)
# 4. convergence, sized for one IEEE column (the 3-panel figures/mc_convergence.png is too wide)
conv = pd.read_csv(PROC / "mc_convergence.csv")
exact = pd.read_csv(PROC / "pc_exact.csv").set_index("conj_id").pc
style = {("mc", "gaussian"): ("#2f6fb0", "o", "MC, Gaussian"), ("mc", "student_t"): ("#c8641e", "^", "MC, Student-t"),
         ("is", "gaussian"): ("#2f6fb0", "s", "IS, Gaussian"), ("is", "student_t"): ("#c8641e", "D", "IS, Student-t")}
fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
for ax, cid in zip(axes, (183, 1198)):
    for (method, model), (col, mk, lab) in style.items():
        s = conv[(conv.conj_id == cid) & (conv.method == method) & (conv.error_model == model)].sort_values("n_samples")
        ax.errorbar(s.n_samples, s.pc, yerr=[s.pc - s.ci_low, s.ci_high - s.pc], color=col, marker=mk, ms=3.5,
                    lw=1, ls="-" if method == "is" else "--", capsize=2, label=lab)
    ax.axhline(exact[cid], color="black", lw=0.8, ls=":", label="exact (Gaussian)")
    ax.set_xscale("log")
    # zoom on the converged range; the wide intervals at 10^3 samples run off the top
    settled = conv[(conv.conj_id == cid) & (conv.n_samples >= 100_000)]
    ax.set_ylim(0, 1.6 * settled.ci_high.max())
    ax.set_title(f"conjunction {cid}", fontsize=8)
    ax.set_ylabel("Pc")
axes[0].legend(fontsize=6, ncol=2, loc="upper right")
axes[1].set_xlabel("number of samples")
fig.tight_layout()
fig.savefig(OUT / "report_convergence.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# 5. Gaussian vs Student-t Pc per conjunction, and 6. ROC of Gaussian Pc against Student-t labels
res = pd.read_csv(PROC / "pc_results.csv")
g = res[res.method == "exact"].set_index("conj_id").pc
t = res[(res.method == "is") & (res.error_model == "student_t")].set_index("conj_id").pc.reindex(g.index)
show = (g > 1e-9) | (t > 1e-9)
fig, ax = plt.subplots(figsize=(4.2, 3.6))
ax.loglog(g[show].clip(lower=1e-9), t[show].clip(lower=1e-9), "o", ms=2.5, color="#2f6fb0", alpha=0.6)
lims = [1e-9, 1e-3]
ax.plot(lims, lims, color="black", lw=0.8, ls=":")
ax.axvline(1e-5, color="#c8641e", lw=0.8, ls="--")
ax.axhline(1e-5, color="#c8641e", lw=0.8, ls="--")
for cid in (183, 1198, 1605):
    ax.annotate(str(cid), (g[cid], t[cid]), textcoords="offset points", xytext=(4, -10), fontsize=7)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel("Pc, Gaussian error model")
ax.set_ylabel("Pc, Student-t error model")
fig.savefig(OUT / "report_gauss_vs_t.png", dpi=300, bbox_inches="tight")
plt.close(fig)

label = (t >= 1e-5).to_numpy()
score = g.to_numpy()
order = np.argsort(-score, kind="stable")
tpr = np.concatenate([[0], np.cumsum(label[order]) / label.sum()])
fpr = np.concatenate([[0], np.cumsum(~label[order]) / (~label).sum()])
auc = np.trapezoid(tpr, fpr)
fig, ax = plt.subplots(figsize=(3.4, 3.2))
ax.plot(fpr, tpr, color="#2f6fb0", lw=1.5)
ax.plot([0, 1], [0, 1], color="grey", lw=0.7, ls=":")
pred = score >= 1e-5
ax.plot((pred & ~label).sum() / (~label).sum(), (pred & label).sum() / label.sum(), "o", color="#c8641e", ms=5)
ax.annotate("Gaussian alert at 10⁻⁵", ((pred & ~label).sum() / (~label).sum(), (pred & label).sum() / label.sum()),
            textcoords="offset points", xytext=(8, -4), fontsize=7)
ax.set_xlabel("false-alarm rate")
ax.set_ylabel("detection rate")
ax.set_title(f"AUC = {auc:.3f}", fontsize=9)
fig.savefig(OUT / "report_roc.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"ROC AUC recomputed from pc_results.csv: {auc:.4f}")

print(f"figures written to {OUT}; Pc = 0 for {(pc == 0).sum()} of {len(pc)} conjunctions; "
      f"min positive log10 Pc {np.log10(pos).min():.1f}")
