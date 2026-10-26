"""
One-time evaluation on the 50 BBBC039 test images with the locked parameters.

Compares the tuned pipeline against (a) the baseline (Otsu threshold + connected
components, nothing else) and (b) ablations that remove one step at a time.
Paired Wilcoxon signed-rank tests are run on per-image F1 and absolute count error.

Usage:  python reference/python/evaluate_test.py --data data/BBBC039 --out results
Writes: results/test_per_image.csv, results/test_summary.csv
"""
import argparse, json, os, sys
from dataclasses import replace

import numpy as np
import pandas as pd
import tifffile
from scipy.stats import wilcoxon
from skimage import io

sys.path.insert(0, os.path.dirname(__file__))
import cellcount as cc  # noqa: E402
from tune_bbbc039 import read_split  # noqa: E402


def variants(best: cc.Params) -> dict:
    v = {"Tuned pipeline": best,
         "Baseline (Otsu + count only)": replace(best, **cc.BASELINE)}
    if best.use_tophat:
        v["Without background removal"] = replace(best, use_tophat=False)
    if best.filter_type != "none":
        v["Without noise filter"] = replace(best, filter_type="none")
    if best.use_cleanup:
        v["Without clean-up"] = replace(best, use_cleanup=False)
    if best.use_watershed:
        v["Without watershed"] = replace(best, use_watershed=False)
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/BBBC039")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()

    best = cc.Params(**json.load(open(os.path.join(a.out, "best_params.json")))["params"])
    split = read_split(os.path.join(a.data, "metadata"))
    test = sorted(s for s, v in split.items() if v == "test")
    V = variants(best)

    rows = []
    for stem in test:
        img = tifffile.imread(os.path.join(a.data, "images", stem + ".tif"))
        gt = cc.load_bbbc039_mask(io.imread(os.path.join(a.data, "masks", stem + ".png")))
        for name, p in V.items():
            _, _, cents = cc.count_cells(img, p)
            m = cc.match_detections(cents, gt)
            m.update(image=stem, method=name)
            rows.append(m)
    per = pd.DataFrame(rows)
    per["error"] = per.n_auto - per.n_true
    per["abs_error"] = per.error.abs()
    per["pct_error"] = 100 * per.abs_error / per.n_true.where(per.n_true >= 5)
    per.to_csv(os.path.join(a.out, "test_per_image.csv"), index=False)

    full = per[per.method == "Tuned pipeline"].set_index("image")
    summ = []
    for name in V:
        d = per[per.method == name].set_index("image").loc[full.index]
        row = dict(method=name, n_images=len(d),
                   total_true=int(d.n_true.sum()), total_auto=int(d.n_auto.sum()),
                   precision=d.tp.sum() / d.n_auto.sum(), recall=d.tp.sum() / d.n_true.sum(),
                   mean_f1=d.f1.mean(), count_mae=d.abs_error.mean(),
                   mape_pct=d.pct_error.mean(), bias=d.error.mean(),
                   loa_low=d.error.mean() - 1.96 * d.error.std(),
                   loa_high=d.error.mean() + 1.96 * d.error.std(),
                   r2=np.corrcoef(d.n_auto, d.n_true)[0, 1] ** 2)
        if name != "Tuned pipeline":
            row["delta_f1_vs_tuned"] = (d.f1 - full.f1).mean()
            row["p_f1"] = wilcoxon(d.f1, full.f1).pvalue
            row["p_abs_error"] = wilcoxon(d.abs_error, full.abs_error).pvalue
        summ.append(row)
    summ = pd.DataFrame(summ)
    summ["pooled_f1"] = 2 * summ.precision * summ.recall / (summ.precision + summ.recall)
    summ.to_csv(os.path.join(a.out, "test_summary.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
