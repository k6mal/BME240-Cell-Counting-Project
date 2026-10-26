classdef SyntheticTest < matlab.unittest.TestCase
    %SYNTHETICTEST  matchDetections and countCells on hand-built images (no data needed).

    methods (TestClassSetup)
        function addCodePath(tc)
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(fileparts(mfilename('fullpath')))));
        end
    end

    methods (Test)
        function matchOneTpOneFpOneFn(tc)
            G = zeros(40, 40, 'uint16');
            G(5:15, 5:15)   = 1;            % matched by the first detection
            G(25:35, 25:35) = 2;            % missed (false negative)
            cents = [10 10;                 % [x y] inside nucleus 1 -> TP
                     30 8];                 % background -> FP
            m = matchDetections(cents, G);
            tc.verifyEqual([m.tp m.fp m.fn], [1 1 1]);
            tc.verifyEqual(m.precision, 0.5, 'AbsTol', 1e-12);
            tc.verifyEqual(m.recall, 0.5, 'AbsTol', 1e-12);
        end

        function matchUsesXYOrder(tc)
            % a centroid is [x y] = [column row]; swapping them must miss
            G = zeros(30, 60, 'uint16');  G(5:10, 40:50) = 1;
            tc.verifyEqual(matchDetections([45 7], G).tp, 1);
            tc.verifyEqual(matchDetections([7 45], G).tp, 0);
        end

        function matchIsOneToOne(tc)
            G = zeros(20, 20, 'uint16');  G(5:15, 5:15) = 1;
            m = matchDetections([8 8; 12 12], G);
            tc.verifyEqual([m.tp m.fp m.fn], [1 1 0]);
        end

        function threeSeparateDiscs(tc)
            [I, p] = SyntheticTest.discs([60 60; 60 160; 160 110], 15);
            tc.verifyEqual(countCells(I, p), 3);
        end

        function twoOverlappingDiscs(tc)
            [I, p] = SyntheticTest.discs([100 88; 100 112], 15);   % centres 24 px apart
            tc.verifyEqual(countCells(I, p), 2, 'watershed on');
            p.useWatershed = false;
            tc.verifyEqual(countCells(I, p), 1, 'watershed off');
        end
    end

    methods (Static)
        function [I, p] = discs(centres, radius)
            % bright discs on a dim, slightly noisy background
            [x, y] = meshgrid(1:220, 1:220);
            I = 0.1 * ones(220);
            for k = 1:size(centres, 1)
                I((y - centres(k, 1)).^2 + (x - centres(k, 2)).^2 <= radius^2) = 0.9;
            end
            rs = RandStream('mt19937ar', 'Seed', 1);
            I = I + 0.01 * randn(rs, size(I));
            p = defaultParams();
            p.cellRadius = radius;
        end
    end
end
