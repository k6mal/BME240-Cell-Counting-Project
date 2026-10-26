"""
cellcount.py - Python reference implementation of the nucleus-counting pipeline.

It mirrors matlab/countCells.m step for step so results from the two languages
can be compared:

    MATLAB                         Python (this file)
    -----------------------------  ---------------------------------------------
    mat2gray                        normalize()          (min-max to [0, 1])
    imtophat(I, strel('disk',R))    morphology.white_tophat(I, disk(R))
    imgaussfilt(I, sigma)           filters.gaussian(I, sigma)
    medfilt2(I, [k k], 'symmetric') ndimage.median_filter(I, k, mode='reflect')
    graythresh + imbinarize         filters.threshold_otsu
    imboxfilt (local mean) + c      ndimage.uniform_filter + c
    imfill(BW, 'holes')             ndimage.binary_fill_holes
    imopen(BW, strel('disk',1))     ndimage.binary_opening(BW, disk(1))
    bwareaopen(BW, P)               remove_small()  (keeps objects with >= P px)
    bwdist -> imgaussfilt -> imextendedmin(h)  distance_transform_edt -> gaussian -> h_maxima
    watershed -> BW(L==0) = 0       segmentation.watershed(..., watershed_line=True)
    bwconncomp(BW, 4)               measure.label(BW, connectivity=1)
    regionprops Centroid            measure.regionprops centroid
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, replace

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    cell_radius: float = 14.0          # typical nucleus radius r (px), calibrated from data
    use_tophat: bool = True
    tophat_factor: float = 3.0         # top-hat disk radius = factor * r
    filter_type: str = "gaussian"      # 'none' | 'gaussian' | 'median'
    sigma: float = 1.0                 # Gaussian std (px)
    median_size: int = 3               # median window (px, odd)
    thresh_method: str = "otsu"        # 'otsu' | 'adaptive'
    adapt_window_factor: float = 4.0   # adaptive window width = 2*floor(f*r)+1
    adapt_offset: float = 0.05         # pixel is cell if I > local_mean + offset
    open_radius: int = 1
    min_area_frac: float = 0.25        # remove objects < frac * pi r^2
    use_cleanup: bool = True
    use_watershed: bool = True
    h_min: float = 1.0                 # h-minima depth (px)
    dist_sigma: float = 2.0            # Gaussian smoothing of the distance map (px); 0 = off

    def to_dict(self):
        return asdict(self)


BASELINE = dict(use_tophat=False, filter_type="none", thresh_method="otsu",
                use_cleanup=False, use_watershed=False)


# --------------------------------------------------------------------------
# Pipeline steps
# --------------------------------------------------------------------------
def normalize(img: np.ndarray) -> np.ndarray:
    """Grayscale, float, min-max scaled to [0, 1] (MATLAB mat2gray)."""
    img = np.asarray(img, dtype=np.float64)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=2)
    lo, hi = img.min(), img.max()
    return (img - lo) / (hi - lo) if hi > lo else np.zeros_like(img)


def remove_small(bw: np.ndarray, min_px: int) -> np.ndarray:
    """Keep 8-connected objects with at least min_px pixels (MATLAB bwareaopen)."""
    lab, n = ndi.label(bw, structure=np.ones((3, 3)))
    if n == 0:
        return bw
    sizes = np.bincount(lab.ravel())
    keep = sizes >= min_px
    keep[0] = False
    return keep[lab]


def preprocess(img: np.ndarray, p: Params) -> np.ndarray:
    """Steps 1-3: normalize, top-hat background removal, noise filter."""
    I = normalize(img)
    if p.use_tophat:
        rad = int(round(p.tophat_factor * p.cell_radius))
        I = morphology.white_tophat(I, morphology.disk(rad, decomposition="sequence"))
        I = normalize(I)
    if p.filter_type == "gaussian":
        I = filters.gaussian(I, sigma=p.sigma, mode="nearest", truncate=2.0)
    elif p.filter_type == "median":
        I = ndi.median_filter(I, size=p.median_size, mode="reflect")
    return I


def threshold(I: np.ndarray, p: Params) -> np.ndarray:
    """Step 4: global Otsu or local-mean adaptive threshold."""
    if p.thresh_method == "otsu":
        return I > filters.threshold_otsu(I)
    if p.thresh_method == "adaptive":
        w = 2 * int(np.floor(p.adapt_window_factor * p.cell_radius)) + 1
        local_mean = ndi.uniform_filter(I, size=w, mode="nearest")
        return I > local_mean + p.adapt_offset
    raise ValueError(p.thresh_method)


def cleanup(bw: np.ndarray, p: Params) -> np.ndarray:
    """Step 5: fill holes, opening, remove objects smaller than frac * pi r^2."""
    if not p.use_cleanup:
        return bw
    bw = ndi.binary_fill_holes(bw)
    bw = ndi.binary_opening(bw, structure=morphology.disk(p.open_radius))
    return remove_small(bw, int(round(p.min_area_frac * np.pi * p.cell_radius ** 2)))


def split_touching(bw: np.ndarray, p: Params) -> np.ndarray:
    """Step 6: distance-transform watershed with h-minima suppression."""
    if not p.use_watershed or not bw.any():
        return bw
    dist = ndi.distance_transform_edt(bw)
    if p.dist_sigma > 0:
        # The raw Euclidean distance map has pixel-scale ridges that create
        # spurious peaks and over-split single nuclei; light smoothing removes them.
        dist = filters.gaussian(dist, sigma=p.dist_sigma, mode="nearest", truncate=2.0)
    markers = measure.label(morphology.h_maxima(dist, p.h_min), connectivity=2)
    lab = segmentation.watershed(-dist, markers, mask=bw, watershed_line=True)
    out = lab > 0
    if p.use_cleanup:   # drop slivers left behind by the watershed lines
        out = remove_small(out, int(round(p.min_area_frac * np.pi * p.cell_radius ** 2)))
    return out


def count_objects(bw: np.ndarray):
    """Step 7: 4-connected component labeling; count = number of objects."""
    lab = measure.label(bw, connectivity=1)
    props = measure.regionprops(lab)
    cents = np.array([pr.centroid for pr in props]).reshape(-1, 2)  # (row, col)
    areas = np.array([pr.area for pr in props])
    return lab, cents, areas


def count_cells(img: np.ndarray, p: Params = Params(), return_steps: bool = False):
    """Full pipeline. Returns (count, label_image, centroids[, steps])."""
    I = preprocess(img, p)
    bw0 = threshold(I, p)
    bw1 = cleanup(bw0, p)
    bw2 = split_touching(bw1, p)
    lab, cents, _ = count_objects(bw2)
    if return_steps:
        steps = dict(normalized=normalize(img), preprocessed=I,
                     threshold=bw0, cleaned=bw1, split=bw2)
        return lab.max(), lab, cents, steps
    return lab.max(), lab, cents


# --------------------------------------------------------------------------
# Ground truth and scoring
# --------------------------------------------------------------------------
def load_bbbc039_mask(mask: np.ndarray) -> np.ndarray:
    """Decode a BBBC039 mask PNG into an instance label image.

    The masks store a small 'map colouring' in the red channel (values 1, 2, 3):
    touching nuclei always get different values, so every nucleus is a
    4-connected component of a single value.
    """
    r = mask[..., 0] if mask.ndim == 3 else mask
    out = np.zeros(r.shape, dtype=np.int32)
    offset = 0
    for v in np.unique(r):
        if v == 0:
            continue
        lab = measure.label(r == v, connectivity=1)
        lab[lab > 0] += offset
        out = np.where(lab > 0, lab, out)
        offset = out.max()
    return out


def match_detections(centroids: np.ndarray, gt_labels: np.ndarray) -> dict:
    """One-to-one matching: a detection is a true positive if its centroid lies
    inside an annotated nucleus that has not already been matched."""
    n_gt = int(gt_labels.max())
    matched = set()
    tp = 0
    H, W = gt_labels.shape
    for r, c in centroids:
        g = gt_labels[min(int(round(r)), H - 1), min(int(round(c)), W - 1)]
        if g > 0 and g not in matched:
            matched.add(g)
            tp += 1
    n_det = len(centroids)
    fp, fn = n_det - tp, n_gt - tp
    prec = tp / n_det if n_det else 0.0
    rec = tp / n_gt if n_gt else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return dict(n_true=n_gt, n_auto=n_det, tp=tp, fp=fp, fn=fn,
                precision=prec, recall=rec, f1=f1)
