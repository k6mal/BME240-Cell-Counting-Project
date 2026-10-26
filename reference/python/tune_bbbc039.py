"""
Grid search on the BBBC039 training + validation images (test images are never touched).

Usage:  python reference/python/tune_bbbc039.py --data data/BBBC039 --out results
Writes: results/tuning_per_image.csv, results/tuning_grid.csv, results/best_params.json
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")   # one thread per worker process
import argparse, itertools, json, sys
from dataclasses import replace

import numpy as np
import pandas as pd
import tifffile
import multiprocessing as mp
from functools import partial
from skimage import io

sys.path.insert(0, os.path.dirname(__file__))
import cellcount as cc  # noqa: E402

TOPHAT = [False, True]
FILTERS = [("none", None), ("gaussian", 1.0), ("gaussian", 2.0), ("median", 3), ("median", 5)]
THRESH = [("otsu", None), ("adaptive", 0.02), ("adaptive", 0.05)]
MIN_AREA = [0.10, 0.25]
# watershed settings as (h-minima depth, distance-map smoothing sigma); None = watershed off
H_VALUES = [None, (0.5, 2.0), (1.0, 1.0), (1.0, 2.0), (1.0, 3.0), (2.0, 2.0)]


def read_split(meta_dir):
    split = {}
    for s in ("training", "validation", "test"):
        for line in open(os.path.join(meta_dir, f"{s}.txt")):
            if line.strip():
                split[line.strip().replace(".png", "")] = s
    return split


def calibrate_radius(data, stems):
    areas = []
    for st in stems:
        gt = cc.load_bbbc039_mask(io.imread(os.path.join(data, "masks", st + ".png")))
        areas.extend(np.bincount(gt.ravel())[1:].tolist())
    return float(np.sqrt(np.median(areas) / np.pi))


def make_params(r, tophat, filt, thr, min_area, h):
    p = cc.Params(cell_radius=r, use_tophat=tophat, filter_type=filt[0],
                  thresh_method=thr[0], min_area_frac=min_area,
                  use_watershed=h is not None,
                  h_min=h[0] if h is not None else 1.0, dist_sigma=h[1] if h is not None else 2.0)
    if filt[0] == "gaussian":
        p = replace(p, sigma=filt[1])
    if filt[0] == "median":
        p = replace(p, median_size=filt[1])
    if thr[0] == "adaptive":
        p = replace(p, adapt_offset=thr[1])
    return p


def config_name(tophat, filt, thr, min_area, h):
    f = "none" if filt[0] == "none" else f"{filt[0]}{filt[1]:g}"
    t = "otsu" if thr[0] == "otsu" else f"adaptive{thr[1]:g}"
    return f"tophat{int(tophat)}|{f}|{t}|minA{min_area:g}|h{'off' if h is None else f'{h[0]:g}s{h[1]:g}'}"


def run_image(data, stem, r):
    """Evaluate every configuration on one image, re-using shared intermediate steps."""
    img = tifffile.imread(os.path.join(data, "images", stem + ".tif"))
    gt = cc.load_bbbc039_mask(io.imread(os.path.join(data, "masks", stem + ".png")))
    rows = []
    for tophat, filt in itertools.product(TOPHAT, FILTERS):
        I = cc.preprocess(img, make_params(r, tophat, filt, THRESH[0], 0.25, None))
        for thr in THRESH:
            bw0 = cc.threshold(I, make_params(r, tophat, filt, thr, 0.25, None))
            for min_area in MIN_AREA:
                bw1 = cc.cleanup(bw0, make_params(r, tophat, filt, thr, min_area, None))
                for h in H_VALUES:
                    p = make_params(r, tophat, filt, thr, min_area, h)
                    _, cents, _ = cc.count_objects(cc.split_touching(bw1, p))
                    m = cc.match_detections(cents, gt)
                    m.update(image=stem, config=config_name(tophat, filt, thr, min_area, h))
                    rows.append(m)
    return rows


def _run(data, r, stem):
    return run_image(data, stem, r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/BBBC039")
    ap.add_argument("--out", default="results")
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    split = read_split(os.path.join(a.data, "metadata"))
    tune = [s for s, v in split.items() if v in ("training", "validation")]
    r = calibrate_radius(a.data, [s for s, v in split.items() if v == "training"])
    print(f"calibrated radius r = {r:.2f} px; tuning on {len(tune)} images", flush=True)

    n_jobs = os.cpu_count() - 1 if a.jobs < 1 else a.jobs
    with mp.Pool(n_jobs) as pool:
        out = []
        for i, rows in enumerate(pool.imap_unordered(partial(_run, a.data, r), tune), 1):
            out.append(rows)
            if i % 10 == 0:
                print(f"  {i}/{len(tune)} images", flush=True)
    per = pd.DataFrame([row for rows in out for row in rows])
    per["abs_err"] = (per.n_auto - per.n_true).abs()
    per.to_csv(os.path.join(a.out, "tuning_per_image.csv"), index=False)

    grid = (per.groupby("config")
               .agg(mean_f1=("f1", "mean"), mean_precision=("precision", "mean"),
                    mean_recall=("recall", "mean"), count_mae=("abs_err", "mean"))
               .reset_index()
               .sort_values(["mean_f1", "count_mae"], ascending=[False, True]))
    grid.to_csv(os.path.join(a.out, "tuning_grid.csv"), index=False)

    best = grid.iloc[0].config
    tophat, f, t, mA, h = best.split("|")
    filt = ("none", None) if f == "none" else (f.rstrip("0123456789."), float(f.lstrip("gaussianmedian")))
    thr = ("otsu", None) if t == "otsu" else ("adaptive", float(t.replace("adaptive", "")))
    hv = None if h == "hoff" else tuple(float(x) for x in h[1:].split("s"))
    p = make_params(r, tophat == "tophat1", filt, thr, float(mA[4:]), hv)
    if p.filter_type == "median":
        p = replace(p, median_size=int(p.median_size))
    json.dump(dict(config=best, params=p.to_dict()), open(os.path.join(a.out, "best_params.json"), "w"), indent=2)
    print(grid.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
