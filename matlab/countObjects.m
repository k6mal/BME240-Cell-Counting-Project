function [count, L, centroids] = countObjects(BW)
%COUNTOBJECTS  Step 7: 4-connected labeling; count = number of objects.
%   centroids are [x y] = [column row] (regionprops convention).
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
