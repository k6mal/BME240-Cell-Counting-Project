%% compareStrelDisk.m - default strel('disk',R) (N = 4) vs. the exact disk (N = 0)
%  The top-hat uses a disk of radius 3r = 42 px. strel's default N = 4 is a
%  periodic-line approximation (octagon-like for large R); N = 0 is exact.
%  Runs the tuned pipeline both ways on 5 validation images and reports the
%  differences. Writes results/matlab/strel_disk_comparison.csv.
clear; close all;
P = projectPaths();
val = sort(readSplit(P.data, 'validation'));
stems = val(1:5);
p4 = defaultParams();  p4.tophatStrelN = 4;
p0 = p4;               p0.tophatStrelN = 0;
R = round(p4.tophatFactor * p4.cellRadius);
fprintf('disk radius %d px: N=4 has %d px, N=0 has %d px\n', R, ...
        nnz(getnhood(strel('disk', R, 4))), nnz(getnhood(strel('disk', R, 0))));

n = numel(stems);
image = stems;  count_n4 = zeros(n, 1);  count_n0 = count_n4;
f1_n4 = count_n4;  f1_n0 = count_n4;  mask_agreement = count_n4;  sec_n4 = count_n4;  sec_n0 = count_n4;
for i = 1:n
    I = imread(fullfile(P.data, 'images', [stems{i} '.tif']));
    G = loadBBBC039Mask(fullfile(P.data, 'masks', [stems{i} '.png']));
    tic;  [count_n4(i), ~, c4, s4] = countCells(I, p4);  sec_n4(i) = toc;
    tic;  [count_n0(i), ~, c0, s0] = countCells(I, p0);  sec_n0(i) = toc;
    m4 = matchDetections(c4, G);  m0 = matchDetections(c0, G);
    f1_n4(i) = m4.f1;  f1_n0(i) = m0.f1;
    mask_agreement(i) = mean(s4.split(:) == s0.split(:));
end
T = table(image, count_n4, count_n0, f1_n4, f1_n0, mask_agreement, sec_n4, sec_n0);
disp(T);
fprintf('mean |count difference| %.2f, mean F1 difference %+.4f, mean mask agreement %.4f\n', ...
        mean(abs(count_n4 - count_n0)), mean(f1_n4 - f1_n0), mean(mask_agreement));
writetable(T, fullfile(P.results, 'strel_disk_comparison.csv'));
