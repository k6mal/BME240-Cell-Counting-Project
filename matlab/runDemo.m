%% runDemo.m - count nuclei in one BBBC039 image and score it against the annotation
%  Run after scripts/download_data.sh has fetched the data.
clear; close all;
P       = projectPaths();
dataDir = P.data;
stems   = readSplit(dataDir, 'validation');      % any non-test image
stem    = stems{1};

I = imread(fullfile(dataDir, 'images', [stem '.tif']));
G = loadBBBC039Mask(fullfile(dataDir, 'masks', [stem '.png']));
p = defaultParams();

[count, L, cents, steps] = countCells(I, p);
m = matchDetections(cents, G);
fprintf('%s\n  detected %d | annotated %d | precision %.3f | recall %.3f | F1 %.3f\n', ...
        stem, count, m.nTrue, m.precision, m.recall, m.f1);

figure('Color', 'w', 'Position', [100 100 1400 330]);
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; imshow(imadjust(steps.normalized)); title('1 Normalized');
nexttile; imshow(imadjust(steps.preprocessed)); title('2-3 Background + filter');
nexttile; imshow(steps.threshold);  title('4 Threshold');
nexttile; imshow(steps.split);      title('5-6 Clean-up + watershed');
nexttile; imshow(labeloverlay(imadjust(steps.normalized), L, 'Transparency', 0.5)); hold on
plot(cents(:,1), cents(:,2), 'w+', 'MarkerSize', 4);
title(sprintf('7 Count: %d (true %d)', count, m.nTrue));
