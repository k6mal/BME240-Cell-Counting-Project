"""
Figures for the README / progress report.

Usage:  python reference/python/make_figures.py --data data/BBBC039 --results results --out figures
Writes: figures/fig2_pipeline_steps.png, figures/fig3_test_results.png
"""
import argparse, json, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage import io
from skimage.segmentation import find_boundaries

sys.path.insert(0, os.path.dirname(__file__))
import cellcount as cc  # noqa: E402
from tune_bbbc039 import read_split  # noqa: E402

BLUE, GREY, ORANGE = "#1f5f99", "#8a8a8a", "#c46a1b"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False})


def letter(ax, s):
    ax.text(-0.02, 1.08, s, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")


def stretch(x):
    lo, hi = np.percentile(x, [1, 99.5])
    return np.clip((x - lo) / (hi - lo + 1e-12), 0, 1)


def fig_steps(data, best, stem, out):
    img = tifffile.imread(os.path.join(data, "images", stem + ".tif"))
    gt = cc.load_bbbc039_mask(io.imread(os.path.join(data, "masks", stem + ".png")))
    n, lab, cents, st = cc.count_cells(img, best, return_steps=True)
    m = cc.match_detections(cents, gt)
    crop = (slice(120, 380), slice(200, 546))
    panels = [("1  Normalized image", stretch(st["normalized"]), "gray"),
              ("2–3  " + ("Top-hat background removal" if best.use_tophat else "No background step")
               + ("" if best.filter_type == "none" else f" + {best.filter_type} filter"), stretch(st["preprocessed"]), "gray"),
              ("4  Threshold (binary mask)", st["threshold"], "gray"),
              ("5  Clean-up", st["cleaned"], "gray"),
              ("6  Watershed split", st["split"], "gray")]
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.4))
    for ax, (t, im, cm) in zip(axs.flat, panels):
        ax.imshow(im[crop], cmap=cm, interpolation="nearest"); ax.set_title(t, loc="left")
    ov = np.dstack([stretch(st["normalized"])] * 3)
    ov[find_boundaries(lab, mode="inner")] = [1.0, 0.8, 0.0]
    ax = axs.flat[5]; ax.imshow(ov[crop])
    c = cents[(cents[:, 0] >= crop[0].start) & (cents[:, 0] < crop[0].stop) &
              (cents[:, 1] >= crop[1].start) & (cents[:, 1] < crop[1].stop)]
    ax.plot(c[:, 1] - crop[1].start, c[:, 0] - crop[0].start, "+", color="#00c8ff", ms=4, mew=0.8)
    ax.set_title(f"7  Count: {n} detected / {m['n_true']} annotated", loc="left")
    for a in axs.flat:
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    fig.text(0.01, 0.005, f"Validation image {stem[:14]} (zoomed region). Whole image: precision {m['precision']:.2f}, "
             f"recall {m['recall']:.2f}, F1 {m['f1']:.2f}.", fontsize=6.5, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(out, "fig2_pipeline_steps.png"), dpi=200)
    plt.close(fig)


def fig_results(results, out):
    per = pd.read_csv(os.path.join(results, "test_per_image.csv"))
    summ = pd.read_csv(os.path.join(results, "test_summary.csv")).set_index("method")
    tuned = per[per.method == "Tuned pipeline"].set_index("image")
    base = per[per.method.str.startswith("Baseline")].set_index("image").loc[tuned.index]

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.8], hspace=0.55, wspace=0.35)
    ax = fig.add_subplot(gs[0, 0])
    lim = max(per.n_true.max(), per.n_auto.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="k", lw=0.7, ls="--", zorder=0)
    ax.scatter(base.n_true, base.n_auto, s=12, color=GREY, lw=0, label="Baseline: Otsu + count")
    ax.scatter(tuned.n_true, tuned.n_auto, s=14, color=BLUE, lw=0, label="Tuned pipeline")
    ax.text(lim * 0.97, lim * 0.93, "y = x", fontsize=6.5, ha="right")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Annotated nuclei per image"); ax.set_ylabel("Detected nuclei per image")
    ax.set_title("Counts on the 50 held-out test images", loc="left")
    ax.legend(frameon=False, loc="upper left", handletextpad=0.2, borderaxespad=0.2)
    letter(ax, "a")

    ax = fig.add_subplot(gs[0, 1])
    mean_c = (tuned.n_auto + tuned.n_true) / 2; d = tuned.n_auto - tuned.n_true
    b, sd = d.mean(), d.std()
    ax.axhline(0, color=GREY, lw=0.6, ls=":")
    for y, ls in [(b, "-"), (b - 1.96 * sd, "--"), (b + 1.96 * sd, "--")]:
        ax.axhline(y, color="k", lw=0.7, ls=ls)
    ax.scatter(mean_c, d, s=14, color=BLUE, lw=0, zorder=3)
    ax.text(0.01, b + 1.5, f"bias {b:+.1f}", transform=ax.get_yaxis_transform(), fontsize=6.5, va="bottom")
    ax.text(0.01, b + 1.96 * sd + 1.5, f"+1.96 SD = {b + 1.96 * sd:+.1f}", transform=ax.get_yaxis_transform(), fontsize=6.5, va="bottom")
    ax.text(0.01, b - 1.96 * sd - 1.5, f"\u22121.96 SD = {b - 1.96 * sd:+.1f}", transform=ax.get_yaxis_transform(), fontsize=6.5, va="top")
    ax.set_xlabel("Mean of detected and annotated count"); ax.set_ylabel("Detected \u2212 annotated (nuclei)")
    ax.set_title("Tuned pipeline undercounts by ~10 nuclei/image", loc="left"); ax.margins(y=0.12)
    letter(ax, "b")

    ax = fig.add_subplot(gs[1, :])
    order = [m for m in summ.index if m != "Tuned pipeline"]
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="k", lw=0.7)
    for yi, m in zip(y, order):
        v = summ.loc[m, "delta_f1_vs_tuned"]; p = summ.loc[m, "p_f1"]
        col = GREY if m.startswith("Baseline") else ORANGE
        ax.plot([0, v], [yi, yi], color=col, lw=1.4); ax.plot(v, yi, "o", color=col, ms=4.5)
        ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.2f}"
        ax.text(0.004, yi, f"\u0394F1 = {v:+.3f}  ({ptxt})", fontsize=6.5, va="center", ha="left")
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("Baseline (Otsu + count only)", "Baseline: Otsu + count only") for m in order])
    ax.set_xlim(-0.14, 0.075); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("Change in mean per-image F1 relative to the tuned pipeline (Wilcoxon signed-rank, n = 50)")
    ax.set_title("Ablation: watershed contributes most; background removal has no measurable effect", loc="left")
    letter(ax, "c")
    fig.savefig(os.path.join(out, "fig3_test_results.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/BBBC039")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    best = cc.Params(**json.load(open(os.path.join(a.results, "best_params.json")))["params"])
    split = read_split(os.path.join(a.data, "metadata"))
    val = sorted(s for s, v in split.items() if v == "validation")
    fig_steps(a.data, best, val[3], a.out)
    if os.path.exists(os.path.join(a.results, "test_summary.csv")):
        fig_results(a.results, a.out)


if __name__ == "__main__":
    main()
