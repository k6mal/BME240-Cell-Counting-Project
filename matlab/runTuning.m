function grid = runTuning()
%RUNTUNING  Grid search on the BBBC039 training + validation images.
%   Port of reference/python/tune_bbbc039.py. The test split is never read.
%   360 configurations: top-hat (2) x filter (5) x threshold (3) x
%   minAreaFrac (2) x watershed (6). Intermediate images are shared: each
%   image is preprocessed once per (top-hat, filter), thresholded once per
%   threshold and cleaned once per minAreaFrac.
%   Best = highest mean per-image F1, ties broken by the lower count MAE.
%   Writes results/matlab/tuning_grid.csv and results/matlab/best_params.json.

P = projectPaths();
tune = [readSplit(P.data, 'training'); readSplit(P.data, 'validation')];
r = calibrateRadius(P.data);
fprintf('calibrated radius r = %.2f px; tuning on %d images\n', r, numel(tune));

[cfgs, names] = gridConfigs(r);
nI = numel(tune);  nC = numel(cfgs);
F1 = zeros(nI, nC);  PR = F1;  RE = F1;  AE = F1;

M = 0;                                            % parfor with 0 workers = plain for loop
if canUseParfor(), pool = gcp(); M = pool.NumWorkers; end
parfor (i = 1:nI, M)
    [F1(i, :), PR(i, :), RE(i, :), AE(i, :)] = tuneImage(P.data, tune{i}, cfgs);
end

config = names;
mean_f1 = mean(F1, 1)';  mean_precision = mean(PR, 1)';
mean_recall = mean(RE, 1)';  count_mae = mean(AE, 1)';
grid = table(config, mean_f1, mean_precision, mean_recall, count_mae);
grid = sortrows(grid, {'mean_f1', 'count_mae'}, {'descend', 'ascend'});
writetable(grid, fullfile(P.results, 'tuning_grid.csv'));
disp(grid(1:10, :));

best = cfgs(strcmp(names, grid.config{1}));
fid = fopen(fullfile(P.results, 'best_params.json'), 'w');
fprintf(fid, '%s\n', jsonencode(struct('config', grid.config{1}, 'params', best)));
fclose(fid);

ref = defaultParams();
fields = {'useTophat', 'filterType', 'threshMethod', 'adaptOffset', 'minAreaFrac', ...
          'useWatershed', 'hMin', 'distSigma'};
if strcmp(best.filterType, 'gaussian'), fields{end+1} = 'sigma'; end
if strcmp(best.filterType, 'median'),   fields{end+1} = 'medianSize'; end
same = all(cellfun(@(f) isequal(best.(f), ref.(f)), fields));
if same
    fprintf('Best configuration matches defaultParams.m: %s\n', grid.config{1});
else
    warning('runTuning:differs', ['Best MATLAB configuration (%s) differs from ' ...
        'defaultParams.m. Update defaultParams.m before running runBBBC039Test.'], grid.config{1});
end
end

% -------------------------------------------------------------------------
function [f1, pr, re, ae] = tuneImage(dataDir, stem, cfgs)
% Evaluate every configuration on one image, re-using shared intermediates.
I = imread(fullfile(dataDir, 'images', [stem '.tif']));
G = loadBBBC039Mask(fullfile(dataDir, 'masks', [stem '.png']));
nC = numel(cfgs);
f1 = zeros(1, nC);  pr = f1;  re = f1;  ae = f1;
key = @(c) sprintf('%d|%s|%g|%g', c.useTophat, c.filterType, c.sigma, c.medianSize);
pre = '';  thr = '';  cln = '';
for k = 1:nC
    c = cfgs(k);
    kPre = key(c);
    kThr = sprintf('%s|%s|%g', kPre, c.threshMethod, c.adaptOffset);
    kCln = sprintf('%s|%g', kThr, c.minAreaFrac);
    if ~strcmp(kPre, pre), Ip  = preprocessImage(I, c); pre = kPre; thr = ''; end
    if ~strcmp(kThr, thr), BW0 = thresholdImage(Ip, c); thr = kThr; cln = ''; end
    if ~strcmp(kCln, cln), BW1 = cleanupMask(BW0, c);   cln = kCln; end
    [~, ~, cents] = countObjects(splitTouching(BW1, c));
    m = matchDetections(cents, G);
    f1(k) = m.f1;  pr(k) = m.precision;  re(k) = m.recall;  ae(k) = abs(m.nAuto - m.nTrue);
end
end

function [cfgs, names] = gridConfigs(r)
% Same grid and naming as tune_bbbc039.py, ordered so that configurations
% sharing preprocessing / threshold / clean-up are adjacent.
tophat  = [false true];
filters = {'none', NaN; 'gaussian', 1; 'gaussian', 2; 'median', 3; 'median', 5};
thresh  = {'otsu', NaN; 'adaptive', 0.02; 'adaptive', 0.05};
minArea = [0.10 0.25];
ws      = [NaN NaN; 0.5 2; 1 1; 1 2; 1 3; 2 2];       % (hMin, distSigma); NaN = off
base = defaultParams();
base.cellRadius = r;
cfgs = repmat(base, 0, 1);  names = {};
for t = tophat
    for f = 1:size(filters, 1)
        for h = 1:size(thresh, 1)
            for a = minArea
                for w = 1:size(ws, 1)
                    p = base;
                    p.useTophat = t;
                    p.filterType = filters{f, 1};
                    p.sigma = 1;  p.medianSize = 3;
                    if strcmp(p.filterType, 'gaussian'), p.sigma = filters{f, 2}; end
                    if strcmp(p.filterType, 'median'),   p.medianSize = filters{f, 2}; end
                    p.threshMethod = thresh{h, 1};
                    p.adaptOffset = 0.05;
                    if strcmp(p.threshMethod, 'adaptive'), p.adaptOffset = thresh{h, 2}; end
                    p.minAreaFrac = a;
                    p.useWatershed = ~isnan(ws(w, 1));
                    p.hMin = 1;  p.distSigma = 2;
                    if p.useWatershed, p.hMin = ws(w, 1);  p.distSigma = ws(w, 2); end
                    cfgs(end+1, 1) = p; %#ok<AGROW>
                    names{end+1, 1} = configName(p); %#ok<AGROW>
                end
            end
        end
    end
end
end

function s = configName(p)
% e.g. tophat1|none|adaptive0.02|minA0.1|h0.5s2 (same as tune_bbbc039.py)
switch p.filterType
    case 'none',     f = 'none';
    case 'gaussian', f = sprintf('gaussian%g', p.sigma);
    case 'median',   f = sprintf('median%g', p.medianSize);
end
if strcmp(p.threshMethod, 'otsu'), t = 'otsu'; else, t = sprintf('adaptive%g', p.adaptOffset); end
if p.useWatershed, h = sprintf('h%gs%g', p.hMin, p.distSigma); else, h = 'hoff'; end
s = sprintf('tophat%d|%s|%s|minA%g|%s', p.useTophat, f, t, p.minAreaFrac, h);
end

function tf = canUseParfor()
tf = false;
try
    tf = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
catch
end
end
