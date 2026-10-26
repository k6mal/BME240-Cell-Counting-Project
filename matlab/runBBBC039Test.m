%% runBBBC039Test.m - evaluate the locked parameters on the 50 BBBC039 test images
%  Compares the tuned pipeline with the baseline (Otsu threshold + count only)
%  and with ablations that remove one step at a time. Run once, after
%  runTuning.m has locked defaultParams.m.
%  Writes results/matlab/test_per_image.csv and results/matlab/test_summary.csv.
%  Requires the Image Processing Toolbox. Uses signrank (Statistics Toolbox)
%  if it is installed, otherwise signrankTest.m.
clear; close all;
P = projectPaths();

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

stems = sort(readSplit(P.data, 'test'));
nI = numel(stems);  nV = numel(V);
rows = cell(nI, nV);
for i = 1:nI
    I = imread(fullfile(P.data, 'images', [stems{i} '.tif']));
    G = loadBBBC039Mask(fullfile(P.data, 'masks', [stems{i} '.png']));
    for v = 1:nV
        [~, ~, cents] = countCells(I, V(v).p);
        m = matchDetections(cents, G);
        rows{i, v} = {m.nTrue, m.nAuto, m.tp, m.fp, m.fn, m.precision, m.recall, m.f1, stems{i}, V(v).name};
    end
    fprintf('%2d/%d %s\n', i, nI, stems{i});
end

% ---- per-image table (same columns as the Python reference) -------------
rows = rows';                                  % method varies fastest, as in Python
rows = vertcat(rows{:});
per = cell2table(rows, 'VariableNames', {'n_true', 'n_auto', 'tp', 'fp', 'fn', ...
                 'precision', 'recall', 'f1', 'image', 'method'});
per.error     = per.n_auto - per.n_true;
per.abs_error = abs(per.error);
per.pct_error = 100 * per.abs_error ./ per.n_true;
per.pct_error(per.n_true < 5) = NaN;           % % error only for images with >= 5 nuclei
writetable(per, fullfile(P.results, 'test_per_image.csv'));

% ---- summary ------------------------------------------------------------
haveSignrank = exist('signrank', 'file') == 2 && license('test', 'Statistics_Toolbox');
if haveSignrank, wtest = @(x, y) signrank(x, y); else, wtest = @signrankTest; end
tuned = per(strcmp(per.method, V(1).name), :);
S = cell(nV, 1);
for v = 1:nV
    d = per(strcmp(per.method, V(v).name), :);
    [~, j] = ismember(tuned.image, d.image);  d = d(j, :);   % same image order
    s.method     = V(v).name;
    s.n_images   = height(d);
    s.total_true = sum(d.n_true);
    s.total_auto = sum(d.n_auto);
    s.precision  = sum(d.tp) / sum(d.n_auto);
    s.recall     = sum(d.tp) / sum(d.n_true);
    s.mean_f1    = mean(d.f1);
    s.count_mae  = mean(d.abs_error);
    s.mape_pct   = mean(d.pct_error, 'omitnan');
    s.bias       = mean(d.error);
    s.loa_low    = s.bias - 1.96 * std(d.error);     % Bland-Altman limits of agreement
    s.loa_high   = s.bias + 1.96 * std(d.error);
    R            = corrcoef(d.n_auto, d.n_true);
    s.r2         = R(1, 2)^2;
    if v == 1
        s.delta_f1_vs_tuned = NaN;  s.p_f1 = NaN;  s.p_abs_error = NaN;
    else
        s.delta_f1_vs_tuned = mean(d.f1 - tuned.f1);
        s.p_f1        = wtest(d.f1, tuned.f1);
        s.p_abs_error = wtest(d.abs_error, tuned.abs_error);
    end
    s.pooled_f1 = 2 * s.precision * s.recall / (s.precision + s.recall);
    S{v} = s;
end
T = struct2table(vertcat(S{:}));
disp(T(:, {'method', 'precision', 'recall', 'mean_f1', 'count_mae', 'mape_pct', 'p_f1', 'p_abs_error'}));
writetable(T, fullfile(P.results, 'test_summary.csv'));
if ~haveSignrank, fprintf('p values from signrankTest.m (Statistics Toolbox not found)\n'); end
