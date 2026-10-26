# BME 240: Automated Nucleus Counting

Counts fluorescently stained nuclei using classical image processing, and checks the counts against the hand annotations in the BBBC039 dataset. The MATLAB code is the course deliverable. The Python version implements the same steps and produced the results below.

**Status:** The Python pipeline is tested on BBBC039. The MATLAB code has not been run yet. Not started: the robustness tests on BBBC005 and the manual counts.

## Pipeline

1. Normalize the image to [0, 1].
2. Remove the background with a top-hat filter (disk radius 3r).
3. Noise filter (Gaussian or median). None was selected.
4. Adaptive threshold: a pixel is a nucleus if it is brighter than the local mean (113 px window) plus 0.02.
5. Fill holes, apply an opening (1 px disk), and remove objects smaller than 63 px.
6. Watershed on the distance map, smoothed with σ = 2 px, with h-minima depth 0.5 px.
7. Label the separate objects (`bwconncomp`). The count is the number of objects.

r = 14.1 px is the median nucleus radius in the training annotations. The settings above were chosen from 360 combinations tested on the 150 training and validation images (`results/tuning_grid.csv`).

## Results (50 test images, scored once)

| Method | Precision | Recall | F1 | Count error |
|---|---|---|---|---|
| Full pipeline | 0.987 | 0.899 | 0.934 | 10.2 (9.6%) |
| Baseline: Otsu + count | 0.788 | 0.836 | 0.817 | 12.5 (12.2%) |
| Without watershed | 0.985 | 0.785 | 0.873 | 23.2 (19.9%) |
| Without top-hat | 0.986 | 0.899 | 0.933 | 10.2 (9.6%) |
| Without clean-up | 0.959 | 0.917 | 0.930 | 5.9 (6.2%) |

F1 is the mean per image. A detection counts as correct if its centroid falls inside an annotated nucleus that hasn't already been matched. Count error is the mean absolute error per image. Full per-image results, and p values from paired Wilcoxon signed-rank tests, are in `results/`.

- Watershed matters most. Without it, touching nuclei are counted as one.
- The baseline's count error looks reasonable only because its false detections and its merged nuclei cancel out. Its F1 is much lower (p < 0.001).
- The top-hat and noise filters did not help on BBBC039, which has even lighting and little noise.
- The pipeline undercounts by about 10 nuclei per image, mostly by missing small or dim nuclei.

![Results](figures/fig3_test_results.png)

![Pipeline steps](figures/fig2_pipeline_steps.png)

## Running it

Download the data (BBBC039, about 80 MB, CC0):

```bash
bash scripts/download_data.sh
```

MATLAB (needs the Image Processing Toolbox; `signrank` needs the Statistics Toolbox):

```matlab
cd matlab
runDemo          % one image, step-by-step figure
runBBBC039Test   % test set: full pipeline vs. baseline vs. ablations
```

Python:

```bash
pip install -r requirements.txt
python python/tune_bbbc039.py  --data data/BBBC039 --out results
python python/evaluate_test.py --data data/BBBC039 --out results
python python/make_figures.py  --data data/BBBC039 --results results --out figures
```

## Notes

- Watershed on the raw distance map split many single nuclei in two. Smoothing the map first fixed this.
- The masks store touching nuclei as different values (1–3) in the red channel. Each nucleus is one 4-connected region of a single value (`loadBBBC039Mask.m`).
- One test image (`IXMtest_A12_s7`) fails: it has 69 annotated nuclei but only 19 are detected, because its nuclei are small and broken into fragments.
- MATLAB and scikit-image implement Otsu thresholding and the watershed slightly differently, so the MATLAB results may not match exactly.

## Data

Ljosa et al., *Nature Methods* 9, 637 (2012). Caicedo et al., *Cytometry A* 95, 952 (2019). https://bbbc.broadinstitute.org/BBBC039
