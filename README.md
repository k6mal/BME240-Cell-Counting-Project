# BME 240: Automated Nucleus Counting

Counts fluorescently stained nuclei using classical image processing, and checks the counts against the hand annotations in the BBBC039 dataset. The MATLAB code in `matlab/` is the course deliverable. All results below were produced in MATLAB R2026a. The Python code in `reference/` is an independent cross-check of the same steps: test F1 0.935 in MATLAB, 0.934 in Python (`results/matlab/parity.md`).

**Status:** The MATLAB pipeline, tuning, test, figures and unit tests have been run on BBBC039 in MATLAB R2026a, and all match the Python reference within the agreed tolerances. Not started: the robustness tests on BBBC005 and the manual counts.

## Pipeline

1. Normalize the image to [0, 1].
2. Remove the background with a top-hat filter (disk radius 3r).
3. Noise filter (Gaussian or median). None was selected.
4. Adaptive threshold: a pixel is a nucleus if it is brighter than the local mean (113 px window) plus 0.02.
5. Fill holes, apply an opening (1 px disk), and remove objects smaller than 63 px.
6. Watershed on the distance map, smoothed with σ = 2 px, with h-minima depth 0.5 px.
7. Label the separate objects (`bwconncomp`). The count is the number of objects.

r = 14.1 px is the median nucleus radius in the training annotations. The settings above were chosen from 360 combinations tested on the 150 training and validation images (`runTuning.m`, `results/matlab/tuning_grid.csv`). MATLAB and Python chose the same configuration.

## Results (50 test images, scored once)

| Method | Precision | Recall | F1 | Count error |
|---|---|---|---|---|
| Full pipeline | 0.987 | 0.901 | 0.935 | 10.1 (9.5%) |
| Baseline: Otsu + count | 0.791 | 0.836 | 0.819 | 11.8 (11.8%) |
| Without watershed | 0.985 | 0.787 | 0.874 | 23.0 (19.6%) |
| Without top-hat | 0.986 | 0.901 | 0.935 | 10.0 (9.5%) |
| Without clean-up | 0.960 | 0.917 | 0.931 | 6.0 (6.2%) |

F1 is the mean per image. A detection counts as correct if its centroid falls inside an annotated nucleus that hasn't already been matched. Count error is the mean absolute error per image. Full per-image results, Bland–Altman limits, R² and p values from paired Wilcoxon signed-rank tests are in `results/matlab/` (Python's are in `results/`).

- Watershed matters most. Without it, touching nuclei are counted as one.
- The baseline's count error looks reasonable only because its false detections and its merged nuclei cancel out. Its F1 is much lower (p < 0.001).
- The top-hat and noise filters did not help on BBBC039, which has even lighting and little noise.
- The pipeline undercounts by about 10 nuclei per image (95% limits of agreement −25.3 to +5.4), mostly by missing small or dim nuclei.

![Results](figures/fig3_test_results.png)

![Pipeline steps](figures/fig2_pipeline_steps.png)

## Running it

Download the data (BBBC039, about 80 MB, CC0):

```bash
bash scripts/download_data.sh
```

MATLAB R2020b or newer (run in R2026a) with the Image Processing Toolbox. The Statistics Toolbox (`signrank`) and Parallel Computing Toolbox (`parfor`) are optional; without them `signrankTest.m` and a plain loop are used. Outputs go to `results/matlab/` and `figures/`.

```matlab
cd matlab
runtests('tests')    % unit tests (the BBBC039 tests need the data)
runDemo              % one image, step-by-step figure
runTuning            % 360-configuration grid on training + validation (13 min on 8 workers)
runBBBC039Test       % test set, once: tuned pipeline vs. baseline vs. ablations
makeFigures          % fig1-fig3
compareStrelDisk     % top-hat disk: strel default (N = 4) vs. exact (N = 0)
parityCheck          % MATLAB vs. Python, writes results/matlab/parity.md
```

Python reference:

```bash
pip install -r reference/python/requirements.txt
python reference/python/tune_bbbc039.py  --data data/BBBC039 --out results
python reference/python/evaluate_test.py --data data/BBBC039 --out results
python reference/python/make_figures.py  --data data/BBBC039 --results results --out figures
python reference/python/dump_steps.py    --data data/BBBC039 --results results   # masks for parityCheck
```

## Notes

- Watershed on the raw distance map split many single nuclei in two. Smoothing the map first fixed this.
- The masks store touching nuclei as different values (1–3) in the red channel. Each nucleus is one 4-connected region of a single value (`loadBBBC039Mask.m`).
- One test image (`IXMtest_A12_s7`) fails: it has 69 annotated nuclei but only 20 are detected, because its nuclei are small and broken into fragments.
- Changes to make the MATLAB watershed (`splitTouching.m`) match the Python reference. The algorithm is unchanged.
  - Markers are the tops of peaks at least h high, as in skimage `h_maxima`, built with `imreconstruct`. `imextendedmin` returns the whole plateau within h of each peak. Without clean-up, one plateau covers several noise specks, and every speck it touches survives.
  - Flooding stays inside the foreground, and components with no marker are dropped, as in skimage `watershed(..., mask=bw)`. The old code flooded the whole image and kept them. On synthetic images without clean-up it counted 984–1,630 objects where Python counted 84–93. With these changes it counts 84–94, and 0–1 apart on each image.
  - `watershed` keeps its default 8-connectivity. `watershed(D, 4)` treats a marker whose pixels touch only diagonally as two minima, which over-splits.
- The opening uses `strel('disk', 1, 0)` (the exact 3×3 cross, as in Python). The top-hat disk keeps strel's default approximation (N = 4, 5,877 px vs. 5,525 px for the exact disk). On 5 validation images it gives the same counts and F1 as N = 0, with 99.99% of mask pixels equal, and is about 8× faster (`compareStrelDisk.m`).
- MATLAB vs. Python on the test set: detected counts are identical on 34 of 50 images and never differ by more than 2. On 3 validation images the counts are identical and the masks agree on ≥ 99.9% of pixels. The differences come from the top-hat disk (identical threshold masks with the top-hat off), from scipy's `binary_opening` trimming nuclei that touch the image edge (MATLAB's `imopen` does not), and from where the watershed lines fall (`results/matlab/parity.md`).
- p values come from MATLAB's `signrank`, which uses a normal approximation for n > 15, whereas scipy is exact up to n = 50. This makes the very small p values differ (10⁻⁹ vs. 10⁻¹⁴), but every conclusion is the same.
- Python rounds centroids half-to-even and MATLAB rounds half away from zero, so a few centroids that lie exactly on a pixel edge are scored against a neighbouring pixel.

## Data

Ljosa et al., *Nature Methods* 9, 637 (2012). Caicedo et al., *Cytometry A* 95, 952 (2019). https://bbbc.broadinstitute.org/BBBC039
