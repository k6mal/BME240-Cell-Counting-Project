function [r, areas] = calibrateRadius(dataDir)
%CALIBRATERADIUS  Typical nucleus radius r = sqrt(median nucleus area / pi).
%   Uses the training masks only, so neither the validation nor the test
%   split influences any parameter. areas lists every training nucleus (px).
if nargin < 1, P = projectPaths(); dataDir = P.data; end
stems = readSplit(dataDir, 'training');
areas = cell(numel(stems), 1);
for i = 1:numel(stems)
    G = loadBBBC039Mask(fullfile(dataDir, 'masks', [stems{i} '.png']));
    a = accumarray(double(G(G > 0)), 1);
    areas{i} = a(a > 0);
end
areas = vertcat(areas{:});
r = sqrt(median(areas) / pi);
end
