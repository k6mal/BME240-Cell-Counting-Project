function [n, xy] = manualCount(imgPath)
%MANUALCOUNT  Click-to-count tool for the manual reference counts.
%   [n, xy] = manualCount(imgPath) shows the image; click once on every
%   cell, press Enter when finished. Returns the count and click positions.
%   Each team member counts the same images independently (blinded to the
%   automated result) so inter-observer variability can be reported.

I = imread(imgPath);
fig = figure('Name', 'Manual count', 'NumberTitle', 'off');
imshow(imadjust(mat2gray(im2gray(I))));
title('Click each cell, then press Enter');
hold on
xy = zeros(0, 2);
while true
    [x, y, b] = ginput(1);
    if isempty(b), break; end          % Enter pressed
    xy(end+1, :) = [x y];              %#ok<AGROW>
    plot(x, y, 'r+', 'MarkerSize', 8, 'LineWidth', 1.2);
    title(sprintf('Count = %d  (Enter to finish)', size(xy, 1)));
end
n = size(xy, 1);
close(fig);
end
