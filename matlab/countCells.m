function [count, L, centroids, steps] = countCells(I, p)
%COUNTCELLS  Count fluorescently stained nuclei in a 2-D microscopy image.
%
%   [count, L, centroids, steps] = countCells(I, p)
%
%   I         : grayscale image, any numeric class (8-bit, 16-bit, double)
%   p         : parameter struct (see defaultParams.m)
%   count     : number of nuclei detected
%   L         : label matrix, one integer per detected nucleus
%   centroids : count x 2 array of [x y] = [column row] centroids
%   steps     : struct with every intermediate image (for figures)
%
%   Pipeline (each step mirrors python/cellcount.py):
%     1 normalize to [0,1]                  mat2gray
%     2 background removal                  imtophat, disk radius = tophatFactor * r
%     3 noise filter                        imgaussfilt | medfilt2
%     4 threshold                           Otsu (graythresh) | local mean + offset
%     5 clean-up                            imfill holes, imopen, bwareaopen
%     6 split touching nuclei               smoothed distance transform + h-minima + watershed
%     7 count                               bwconncomp (4-connected), regionprops
%
%   Requires the Image Processing Toolbox.

if nargin < 2, p = defaultParams(); end

% ---- 1. Normalize -------------------------------------------------------
if size(I, 3) == 3, I = rgb2gray(I); end
I = mat2gray(double(I));
steps.normalized = I;

% ---- 2. Background removal (white top-hat) ------------------------------
if p.useTophat
    I = imtophat(I, strel('disk', round(p.tophatFactor * p.cellRadius)));
    I = mat2gray(I);
end

% ---- 3. Noise filter ----------------------------------------------------
switch lower(p.filterType)
    case 'gaussian'
        I = imgaussfilt(I, p.sigma, 'Padding', 'replicate');   % kernel = 2*ceil(2*sigma)+1
    case 'median'
        I = medfilt2(I, [p.medianSize p.medianSize], 'symmetric');
end
steps.preprocessed = I;

% ---- 4. Threshold -------------------------------------------------------
switch lower(p.threshMethod)
    case 'otsu'
        BW = I > graythresh(I);
    case 'adaptive'
        w  = 2 * floor(p.adaptWindowFactor * p.cellRadius) + 1;  % odd window
        BW = I > imboxfilt(I, w, 'Padding', 'replicate') + p.adaptOffset;
    otherwise
        error('countCells:thresh', 'Unknown threshMethod "%s".', p.threshMethod);
end
steps.threshold = BW;

% ---- 5. Clean-up --------------------------------------------------------
if p.useCleanup
    BW = imfill(BW, 'holes');
    BW = imopen(BW, strel('disk', p.openRadius));
    BW = bwareaopen(BW, round(p.minAreaFrac * pi * p.cellRadius^2));
end
steps.cleaned = BW;

% ---- 6. Split touching nuclei (marker-controlled watershed) -------------
if p.useWatershed && any(BW(:))
    D = bwdist(~BW);                    % distance of each cell pixel to background
    if p.distSigma > 0
        % the raw distance map has pixel-scale ridges that create spurious
        % peaks and over-split single nuclei; light smoothing removes them
        D = imgaussfilt(D, p.distSigma, 'Padding', 'replicate');
    end
    D      = -D;                        % nucleus centres become deep minima
    marker = imextendedmin(D, p.hMin);  % keep only minima deeper than hMin px
    D      = imimposemin(D, marker);
    Lw     = watershed(D);
    BW(Lw == 0) = 0;                    % cut along the watershed lines
    if p.useCleanup                     % drop slivers left by the cuts
        BW = bwareaopen(BW, round(p.minAreaFrac * pi * p.cellRadius^2));
    end
end
steps.split = BW;

% ---- 7. Count -----------------------------------------------------------
CC    = bwconncomp(BW, 4);              % 4-conn so diagonal cuts still separate
L     = labelmatrix(CC);
count = CC.NumObjects;
if count > 0
    s = regionprops(CC, 'Centroid');
    centroids = vertcat(s.Centroid);
else
    centroids = zeros(0, 2);
end
end
