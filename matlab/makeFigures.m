function makeFigures()
%MAKEFIGURES  README / report figures (port of reference/python/make_figures.py).
%   fig1_bbbc039_dataset.png : raw image, annotation, nuclei per image, area histogram
%   fig2_pipeline_steps.png  : the pipeline steps on one validation image
%   fig3_test_results.png    : test counts, Bland-Altman, ablation (needs runBBBC039Test)
%   All figures: 300 dpi, 8 pt text. Also writes results/matlab/bbbc039_dataset_summary.csv.

P = projectPaths();
FS = 8;
C.blue = [31 95 153] / 255;  C.grey = [138 138 138] / 255;  C.orange = [196 106 27] / 255;

figDataset(P, FS, C);
figSteps(P, FS);
if exist(fullfile(P.results, 'test_summary.csv'), 'file')
    figResults(P, FS, C);
else
    fprintf('results/matlab/test_summary.csv not found: run runBBBC039Test first for fig3.\n');
end
end

% -------------------------------------------------------------------------
function figDataset(P, FS, C)
splits = {'training', 'validation', 'test'};
image = {};  split = {};  n_nuclei = [];  median_area_px = [];  trainAreas = [];
for s = 1:3
    stems = readSplit(P.data, splits{s});
    for i = 1:numel(stems)
        G = loadBBBC039Mask(fullfile(P.data, 'masks', [stems{i} '.png']));
        a = accumarray(double(G(G > 0)), 1);  a = a(a > 0);
        image{end+1, 1} = stems{i};  split{end+1, 1} = splits{s};          %#ok<AGROW>
        n_nuclei(end+1, 1) = numel(a);  median_area_px(end+1, 1) = median(a); %#ok<AGROW>
        if s == 1, trainAreas = [trainAreas; a]; end                         %#ok<AGROW>
    end
end
T = table(image, split, n_nuclei, median_area_px);
writetable(T, fullfile(P.results, 'bbbc039_dataset_summary.csv'));
r = sqrt(median(trainAreas) / pi);

stem = T.image{1};                                   % first training image
I = imread(fullfile(P.data, 'images', [stem '.tif']));
G = loadBBBC039Mask(fullfile(P.data, 'masks', [stem '.png']));
In = stretch(normalizeImage(I));

f = newFigure(5.4);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);  imshow(In, 'Parent', ax);
title(ax, sprintf('a  Raw image (%s)', stem(1:14)), 'FontSize', FS, 'FontWeight', 'normal', 'Interpreter', 'none');
ax = nexttile(tl);
B = boundarymask(G);
ov = repmat(In, 1, 1, 3);
ov(cat(3, B, B, B)) = reshape(repmat([1 0.8 0], nnz(B), 1), [], 1);
imshow(ov, 'Parent', ax);
title(ax, sprintf('b  Annotation: %d nuclei', max(G(:))), 'FontSize', FS, 'FontWeight', 'normal');

rng(0);                                              % fixed jitter
ax = nexttile(tl);  hold(ax, 'on');
for s = 1:3
    n = T.n_nuclei(strcmp(T.split, splits{s}));
    x = s + 0.25 * (rand(size(n)) - 0.5);
    plot(ax, x, n, 'o', 'MarkerSize', 3, 'MarkerFaceColor', C.blue, 'MarkerEdgeColor', 'none');
    plot(ax, s + [-0.2 0.2], median(n) * [1 1], 'k-', 'LineWidth', 1);
    lbl{s} = sprintf('%s\\newline%d images\\newline%d nuclei', splits{s}, numel(n), sum(n)); %#ok<AGROW>
end
set(ax, 'XTick', 1:3, 'XTickLabel', lbl, 'XTickLabelRotation', 0, 'XLim', [0.5 3.5], 'FontSize', FS, 'Box', 'off');
ylabel(ax, 'Annotated nuclei per image');
title(ax, 'c  Nuclei per image by split', 'FontSize', FS, 'FontWeight', 'normal');

ax = nexttile(tl);
histogram(ax, trainAreas, 0:50:3000, 'FaceColor', C.grey, 'EdgeColor', 'none');
hold(ax, 'on');
xline(ax, median(trainAreas), 'k--', sprintf('median %d px  (r = %.2f px)', round(median(trainAreas)), r), ...
      'FontSize', FS, 'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
set(ax, 'FontSize', FS, 'Box', 'off');
xlabel(ax, 'Nucleus area (px), training masks');  ylabel(ax, 'Nuclei');
title(ax, 'd  Nucleus size', 'FontSize', FS, 'FontWeight', 'normal');
exportgraphics(f, fullfile(P.figures, 'fig1_bbbc039_dataset.png'), 'Resolution', 300);
close(f);
end

% -------------------------------------------------------------------------
function figSteps(P, FS)
p = defaultParams();
val = sort(readSplit(P.data, 'validation'));
stem = val{4};                                       % same image as the Python figure
I = imread(fullfile(P.data, 'images', [stem '.tif']));
G = loadBBBC039Mask(fullfile(P.data, 'masks', [stem '.png']));
[n, L, cents, st] = countCells(I, p);
m = matchDetections(cents, G);
rr = 121:380;  cc = 201:546;                         % zoomed region

if p.useTophat, t2 = 'Top-hat background removal'; else, t2 = 'No background step'; end
if ~strcmp(p.filterType, 'none'), t2 = [t2 ' + ' p.filterType ' filter']; end
panels = {'1  Normalized image', stretch(st.normalized); ['2-3  ' t2], stretch(st.preprocessed);
          '4  Threshold (binary mask)', st.threshold; '5  Clean-up', st.cleaned;
          '6  Watershed split', st.split};
f = newFigure(4.4);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:5
    ax = nexttile(tl);
    imshow(double(panels{k, 2}(rr, cc)), 'Parent', ax);
    title(ax, panels{k, 1}, 'FontSize', FS, 'FontWeight', 'normal');
end
ax = nexttile(tl);
In = stretch(st.normalized);
B = boundarymask(L) & L > 0;
ov = repmat(In, 1, 1, 3);
ov(cat(3, B, B, B)) = reshape(repmat([1 0.8 0], nnz(B), 1), [], 1);
imshow(ov(rr, cc, :), 'Parent', ax);  hold(ax, 'on');
in = cents(:, 2) >= rr(1) & cents(:, 2) <= rr(end) & cents(:, 1) >= cc(1) & cents(:, 1) <= cc(end);
plot(ax, cents(in, 1) - cc(1) + 1, cents(in, 2) - rr(1) + 1, '+', 'Color', [0 0.78 1], 'MarkerSize', 4);
title(ax, sprintf('7  Count: %d detected / %d annotated', n, m.nTrue), 'FontSize', FS, 'FontWeight', 'normal');
xlabel(tl, sprintf(['Validation image %s (zoomed region). Whole image: precision %.2f, ' ...
       'recall %.2f, F1 %.2f.'], stem(1:14), m.precision, m.recall, m.f1), 'FontSize', FS, 'Interpreter', 'none');
exportgraphics(f, fullfile(P.figures, 'fig2_pipeline_steps.png'), 'Resolution', 300);
close(f);
end

% -------------------------------------------------------------------------
function figResults(P, FS, C)
per  = readtable(fullfile(P.results, 'test_per_image.csv'), 'TextType', 'char');
summ = readtable(fullfile(P.results, 'test_summary.csv'), 'TextType', 'char');
tuned = per(strcmp(per.method, 'Tuned pipeline'), :);
base  = per(startsWith(per.method, 'Baseline'), :);
[~, j] = ismember(tuned.image, base.image);  base = base(j, :);

f = newFigure(5.6);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);  hold(ax, 'on');
lim = max([per.n_true; per.n_auto]) * 1.05;
plot(ax, [0 lim], [0 lim], 'k--', 'LineWidth', 0.7, 'HandleVisibility', 'off');
scatter(ax, base.n_true, base.n_auto, 12, C.grey, 'filled');
scatter(ax, tuned.n_true, tuned.n_auto, 14, C.blue, 'filled');
text(ax, lim * 0.97, lim * 0.84, 'y = x', 'FontSize', FS, 'HorizontalAlignment', 'right');
axis(ax, [0 lim 0 lim]);  axis(ax, 'square');
set(ax, 'FontSize', FS, 'Box', 'off');
xlabel(ax, 'Annotated nuclei per image');  ylabel(ax, 'Detected nuclei per image');
legend(ax, {'Baseline: Otsu + count', 'Tuned pipeline'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', FS);
title(ax, 'a  Counts on the 50 test images', 'FontSize', FS, 'FontWeight', 'normal');

ax = nexttile(tl);  hold(ax, 'on');
mc = (tuned.n_auto + tuned.n_true) / 2;  d = tuned.n_auto - tuned.n_true;
b = mean(d);  sd = std(d);
yline(ax, 0, ':', 'Color', C.grey);
yline(ax, b, 'k-', sprintf('bias %+.1f', b), 'FontSize', FS, 'LabelHorizontalAlignment', 'left');
yline(ax, b + 1.96 * sd, 'k--', sprintf('+1.96 SD = %+.1f', b + 1.96 * sd), 'FontSize', FS, 'LabelHorizontalAlignment', 'left');
yline(ax, b - 1.96 * sd, 'k--', sprintf('-1.96 SD = %+.1f', b - 1.96 * sd), 'FontSize', FS, ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
scatter(ax, mc, d, 14, C.blue, 'filled');
set(ax, 'FontSize', FS, 'Box', 'off');
xlabel(ax, 'Mean of detected and annotated count');  ylabel(ax, 'Detected - annotated (nuclei)');
title(ax, 'b  Bland-Altman, tuned pipeline', 'FontSize', FS, 'FontWeight', 'normal');

ax = nexttile(tl, [1 2]);  hold(ax, 'on');
others = summ(~strcmp(summ.method, 'Tuned pipeline'), :);
y = height(others):-1:1;
xline(ax, 0, 'k-');
for k = 1:height(others)
    v = others.delta_f1_vs_tuned(k);  pv = others.p_f1(k);
    if startsWith(others.method{k}, 'Baseline'), col = C.grey; else, col = C.orange; end
    plot(ax, [0 v], y(k) * [1 1], '-', 'Color', col, 'LineWidth', 1.4);
    plot(ax, v, y(k), 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 4.5);
    if pv < 0.001, pt = 'p < 0.001'; elseif pv < 0.01, pt = sprintf('p = %.3f', pv); else, pt = sprintf('p = %.2f', pv); end
    text(ax, 0.004, y(k), sprintf('\\DeltaF1 = %+.3f  (%s)', v, pt), 'FontSize', FS);
end
set(ax, 'YTick', fliplr(y), 'YTickLabel', flipud(strrep(others.method, ...
    'Baseline (Otsu + count only)', 'Baseline: Otsu + count only')), 'FontSize', FS, 'Box', 'off');
xlim(ax, [min(-0.14, min(others.delta_f1_vs_tuned) - 0.02) 0.075]);  ylim(ax, [0.4 height(others) + 0.6]);
xlabel(ax, 'Change in mean per-image F1 vs. the tuned pipeline (Wilcoxon signed-rank, n = 50)');
title(ax, 'c  Ablation: remove one step at a time', 'FontSize', FS, 'FontWeight', 'normal');
exportgraphics(f, fullfile(P.figures, 'fig3_test_results.png'), 'Resolution', 300);
close(f);
end

function f = newFigure(heightIn)
% 7.2 in wide, white; R2025a+ draws figures in the desktop (often dark) theme
f = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7.2 heightIn]);
if isprop(f, 'Theme'), f.Theme = 'light'; end
end

function J = stretch(I)
% contrast stretch between the 1st and 99.5th percentile, for display only
% (nearest-rank percentiles, so no Statistics Toolbox is needed)
v = sort(I(:));
q = v(max(1, round([0.01 0.995] * numel(v))));
J = min(max((I - q(1)) / (q(2) - q(1) + 1e-12), 0), 1);
end
