function BW = splitTouching(BW, p)
%SPLITTOUCHING  Step 6: marker-controlled watershed on the smoothed distance map.
%   Mirrors skimage watershed(-dist, markers, mask=bw, watershed_line=True):
%   flooding is restricted to the foreground, and foreground components that
%   contain no marker are dropped, because the masked skimage watershed
%   leaves them unlabelled. watershed keeps its default 8-connectivity: its
%   minima must match the 8-connected markers (skimage labels markers with
%   connectivity=2); watershed(D, 4) would split diagonally touching
%   markers into extra basins.
if ~p.useWatershed || ~any(BW(:)), return; end

D = double(bwdist(~BW));                % distance of each cell pixel to background
if p.distSigma > 0
    % the raw distance map has pixel-scale ridges that create spurious
    % peaks and over-split single nuclei; light smoothing removes them
    D = imgaussfilt(D, p.distSigma, 'Padding', 'replicate');
end
% Markers: the top of every peak at least hMin px high (skimage h_maxima).
% Not imextendedmax, which returns the whole plateau within hMin of the
% peak; without clean-up one plateau can cover several noise specks.
% The 2e-15*|D| term is h_maxima's floating-point tolerance.
rec    = imreconstruct(D - p.hMin - 2e-15 * abs(D), D);   % 8-connected
marker = (D - rec) >= p.hMin & BW;
D      = -D;                            % nucleus centres become deep minima
D(~BW) = max(D(:)) + 1;                 % background floods last: basins grow inside cells only
D      = imimposemin(D, marker);
Lw     = watershed(D);

BW = imreconstruct(marker, BW, 4);      % drop components with no marker (unlabelled in skimage)
BW = BW & Lw > 0;                       % cut along the watershed lines
if p.useCleanup                         % drop slivers left by the cuts
    BW = bwareaopen(BW, minObjectArea(p));
end
end
