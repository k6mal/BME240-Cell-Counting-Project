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
%   Pipeline (each step mirrors reference/python/cellcount.py):
%     1 normalize to [0,1]                  mat2gray                     (normalizeImage)
%     2 background removal                  imtophat, disk radius = 3 r  (preprocessImage)
%     3 noise filter                        imgaussfilt | medfilt2       (preprocessImage)
%     4 threshold                           Otsu | local mean + offset   (thresholdImage)
%     5 clean-up                            imfill, imopen, bwareaopen   (cleanupMask)
%     6 split touching nuclei               distance map + h-minima + watershed (splitTouching)
%     7 count                               bwconncomp (4-connected)     (countObjects)
%
%   Requires the Image Processing Toolbox.

if nargin < 2, p = defaultParams(); end

steps.normalized   = normalizeImage(I);
steps.preprocessed = preprocessImage(I, p);
steps.threshold    = thresholdImage(steps.preprocessed, p);
steps.cleaned      = cleanupMask(steps.threshold, p);
steps.split        = splitTouching(steps.cleaned, p);
[count, L, centroids] = countObjects(steps.split);
end
