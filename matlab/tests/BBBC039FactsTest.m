classdef BBBC039FactsTest < matlab.unittest.TestCase
    %BBBC039FACTSTEST  Facts about the BBBC039 annotations, from the Python reference run.
    %   Skipped (not failed) when data/BBBC039 has not been downloaded.

    properties
        dataDir
    end

    methods (TestClassSetup)
        function setup(tc)
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(fileparts(mfilename('fullpath')))));
            P = projectPaths();
            tc.dataDir = P.data;
            tc.assumeTrue(exist(fullfile(P.data, 'metadata', 'test.txt'), 'file') == 2, ...
                'BBBC039 not downloaded (run scripts/download_data.sh)');
        end
    end

    methods (Test)
        function oneMask(tc)
            G = loadBBBC039Mask(fullfile(tc.dataDir, 'masks', ...
                'IXMtest_A02_s1_w1051DAA7C-7042-435F-99F0-1E847D9B42CB.png'));
            a = accumarray(double(G(G > 0)), 1);
            tc.verifyEqual(double(max(G(:))), 110);
            tc.verifyEqual(nnz(a), 110);              % labels 1..N are all used
            tc.verifyEqual(median(a), 649);
        end

        function maskEncoding(tc)
            f = fullfile(tc.dataDir, 'masks', 'IXMtest_A02_s1_w1051DAA7C-7042-435F-99F0-1E847D9B42CB.png');
            [M, ~, alpha] = imread(f);
            tc.verifyEqual(size(M, 3), 3);            % RGB, alpha returned separately
            tc.verifyNotEmpty(alpha);
            tc.verifyTrue(all(ismember(unique(M(:, :, 1)), uint8(0:3))));
        end

        function nucleiPerSplit(tc)
            splits = {'training', 'validation', 'test'};
            n = zeros(1, 3);  nImg = zeros(1, 3);
            for s = 1:3
                stems = readSplit(tc.dataDir, splits{s});
                nImg(s) = numel(stems);
                for i = 1:numel(stems)
                    G = loadBBBC039Mask(fullfile(tc.dataDir, 'masks', [stems{i} '.png']));
                    n(s) = n(s) + double(max(G(:)));
                end
            end
            tc.verifyEqual(nImg, [100 50 50]);
            tc.verifyEqual(n(1), 12001);
            tc.verifyEqual(n(3), 5720);
            tc.verifyEqual(sum(n), 23617);
        end

        function calibratedRadius(tc)
            [r, areas] = calibrateRadius(tc.dataDir);
            tc.verifyEqual(median(areas), 627);
            tc.verifyEqual(r, sqrt(627 / pi), 'AbsTol', 1e-12);
            tc.verifyEqual(round(r, 2), 14.13);
        end
    end
end
