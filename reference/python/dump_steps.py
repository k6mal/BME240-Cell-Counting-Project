"""
Save the intermediate masks of the tuned pipeline for a few validation images,
so matlab/parityCheck.m can compare MATLAB and Python step by step.

Usage:  python reference/python/dump_steps.py --data data/BBBC039 --results results --n 3
Writes: results/parity_steps/<image>_{threshold,cleaned,split}.png
"""
import argparse, json, os, sys

import numpy as np
import tifffile
from skimage import io

sys.path.insert(0, os.path.dirname(__file__))
import cellcount as cc  # noqa: E402
from tune_bbbc039 import read_split  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/BBBC039")
    ap.add_argument("--results", default="results")
    ap.add_argument("--n", type=int, default=3)
    a = ap.parse_args()
    out = os.path.join(a.results, "parity_steps")
    os.makedirs(out, exist_ok=True)
    best = cc.Params(**json.load(open(os.path.join(a.results, "best_params.json")))["params"])
    split = read_split(os.path.join(a.data, "metadata"))
    val = sorted(s for s, v in split.items() if v == "validation")
    for stem in val[:a.n]:
        img = tifffile.imread(os.path.join(a.data, "images", stem + ".tif"))
        n, _, _, st = cc.count_cells(img, best, return_steps=True)
        for k in ("threshold", "cleaned", "split"):
            io.imsave(os.path.join(out, f"{stem}_{k}.png"), (st[k] * 255).astype(np.uint8),
                      check_contrast=False)
        print(f"{stem}: {n} nuclei")


if __name__ == "__main__":
    main()
