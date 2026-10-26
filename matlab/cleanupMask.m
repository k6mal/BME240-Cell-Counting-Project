function BW = cleanupMask(BW, p)
%CLEANUPMASK  Step 5: fill holes, opening, remove objects < minAreaFrac * pi r^2.
if ~p.useCleanup, return; end
BW = imfill(BW, 'holes');
BW = imopen(BW, strel('disk', p.openRadius, 0));  % exact disk (radius 1 = 3x3 cross)
BW = bwareaopen(BW, minObjectArea(p));       % keeps 8-connected objects >= P px
end
