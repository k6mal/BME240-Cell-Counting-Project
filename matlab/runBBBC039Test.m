%% runBBBC039Test.m - evaluate the locked parameters on the 50 BBBC039 test images
%  Compares the tuned pipeline with the baseline (Otsu threshold + count only)
%  and with ablations that remove one step at a time.
%  Requires: Image Processing Toolbox; Statistics and Machine Learning Toolbox (signrank).
clear; close all;
dataDir = fullfile('..', 'data', 'BBBC039');
outDir  = fullfile('..', 'results', 'matlab');
if ~exist(outDir, 'dir'), mkdir(outDir); end

best = defaultParams();
V = struct('name', {}, 'p', {});
V(end+1) = struct('name', 'Tuned pipeline', 'p', best);
b = best; b.useTophat = false; b.filterType = 'none'; b.threshMethod = 'otsu';
b.useCleanup = false; b.useWatershed = false;
V(end+1) = struct('name', 'Baseline (Otsu + count only)', 'p', b);
if best.useTophat,                 a = best; a.useTophat = false;      V(end+1) = struct('name', 'Without background removal', 'p', a); end
if ~strcmp(best.filterType,'none'), a = best; a.filterType = 'none';   V(end+1) = struct('name', 'Without noise filter', 'p', a); end
if best.useCleanup,                a = best; a.useCleanup = false;     V(end+1) = struct('name', 'Without clean-up', 'p', a); end
if best.useWatershed,              a = best; a.useWatershed = false;   V(end+1) = struct('name', 'Without watershed', 'p', a); end

stems = sort(readSplit(dataDir, 'test'));
nI = numel(stems);  nV = numel(V);
F1 = zeros(nI, nV); NA = zeros(nI, nV); NT = zeros(nI, 1); TP = zeros(nI, nV);
for i = 1:nI
    I = imread(fullfile(dataDir, 'images', [stems{i} '.tif']));
    G = loadBBBC039Mask(fullfile(dataDir, 'masks', [stems{i} '.png']));
    for v = 1:nV
        [~, ~, cents] = countCells(I, V(v).p);
        m = matchDetections(cents, G);
        F1(i, v) = m.f1;  NA(i, v) = m.nAuto;  TP(i, v) = m.tp;  NT(i) = m.nTrue;
    end
    fprintf('%2d/%d %s\n', i, nI, stems{i});
end

absErr = abs(NA - NT);
method = {V.name}';
meanF1 = mean(F1)';  countMAE = mean(absErr)';
precision = (sum(TP) ./ sum(NA))';  recall = (sum(TP) ./ sum(NT))';
pF1 = nan(nV, 1);  pErr = nan(nV, 1);
for v = 2:nV
    pF1(v)  = signrank(F1(:, v), F1(:, 1));
    pErr(v) = signrank(absErr(:, v), absErr(:, 1));
end
T = table(method, meanF1, precision, recall, countMAE, pF1, pErr);
disp(T);
writetable(T, fullfile(outDir, 'test_summary_matlab.csv'));

% Automated vs annotated count
f = figure('Color', 'w');
plot(NT, NA(:, 2), 'o', 'Color', [0.6 0.6 0.6]); hold on
plot(NT, NA(:, 1), 'o', 'Color', [0.12 0.37 0.6], 'MarkerFaceColor', [0.12 0.37 0.6]);
lim = [0 max([NT; NA(:)]) + 10]; plot(lim, lim, 'k--');
xlabel('Annotated nuclei'); ylabel('Detected nuclei'); axis square
legend({V(2).name, V(1).name, 'y = x'}, 'Location', 'northwest');
exportgraphics(f, fullfile(outDir, 'auto_vs_true_matlab.png'), 'Resolution', 200);
