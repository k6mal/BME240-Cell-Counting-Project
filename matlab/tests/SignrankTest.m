classdef SignrankTest < matlab.unittest.TestCase
    %SIGNRANKTEST  signrankTest.m against p values from scipy.stats.wilcoxon 1.18.

    methods (TestClassSetup)
        function addCodePath(tc)
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(fileparts(mfilename('fullpath')))));
        end
    end

    methods (Test)
        function exactNoTies(tc)
            x = [1.83 0.50 1.62 2.48 1.68 1.88 1.55 3.06 1.30];
            y = [0.878 0.647 0.598 2.05 1.06 1.29 1.06 3.14 1.29];
            tc.verifyEqual(signrankTest(x, y), 0.0390625, 'RelTol', 1e-9);
        end

        function normalWithTiesAndZeros(tc)
            d = [3 -1 2 2 0 5 -2 4 1 0 3 -3 6 2 1];
            tc.verifyEqual(signrankTest(d, zeros(size(d))), 0.041558707296959874, 'RelTol', 1e-9);
        end

        function exactAfterDroppingZeros(tc)
            d = [0 1.5 -2.2 3.1 4.7 5.3 -0.4 6.6];
            tc.verifyEqual(signrankTest(d, zeros(size(d))), 0.109375, 'RelTol', 1e-9);
        end

        function normalLargeSample(tc)
            k = 1:60;
            d = k .* (-1).^floor(k / 3);
            tc.verifyEqual(signrankTest(d, zeros(size(d))), 0.9120730845618369, 'RelTol', 1e-9);
        end

        function allZero(tc)
            tc.verifyEqual(signrankTest([1 2 3], [1 2 3]), 1);
        end
    end
end

