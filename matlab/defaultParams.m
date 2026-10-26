function p = defaultParams()
%DEFAULTPARAMS  Locked parameters for countCells, tuned on the BBBC039
%   training + validation images (see python/tune_bbbc039.py and
%   results/best_params.json). All sizes scale with the nucleus radius r.

p.cellRadius        = 14.13;      % r (px) = sqrt(median nucleus area / pi), training masks

% Step 2: background removal
p.useTophat         = true;
p.tophatFactor      = 3;          % top-hat disk radius = 3 r

% Step 3: noise filter ('none' | 'gaussian' | 'median')
p.filterType        = 'none';     % filtering did not help on BBBC039 (low noise)
p.sigma             = 1;
p.medianSize        = 3;

% Step 4: threshold ('otsu' | 'adaptive')
p.threshMethod      = 'adaptive'; % pixel is nucleus if I > local mean + offset
p.adaptWindowFactor = 4;          % window = 2*floor(4 r)+1 = 113 px
p.adaptOffset       = 0.02;

% Step 5: clean-up
p.useCleanup        = true;
p.openRadius        = 1;
p.minAreaFrac       = 0.10;       % remove objects < 0.10 * pi r^2 (~63 px)

% Step 6: watershed
p.useWatershed      = true;
p.hMin              = 0.5;        % h-minima depth (px)
p.distSigma         = 2;          % Gaussian smoothing of the distance map (px)
end
