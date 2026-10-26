function p = signrankTest(x, y)
%SIGNRANKTEST  Two-sided paired Wilcoxon signed-rank test (no toolbox needed).
%   p = signrankTest(x, y) tests whether the median of x - y is zero.
%   Follows scipy.stats.wilcoxon defaults, which the Python reference used:
%   zero differences are dropped; the exact null distribution is used when
%   n <= 50 and there are no tied |differences|; otherwise the normal
%   approximation with the tie correction and no continuity correction.
%   Used by runBBBC039Test when the Statistics Toolbox (signrank) is absent.
d = x(:) - y(:);
d = d(~isnan(d) & d ~= 0);
n = numel(d);
if n == 0, p = 1; return; end
[rk, ties] = tiedRanks(abs(d));
wPlus = sum(rk(d > 0));

if n <= 50 && ~ties
    % exact: number of subsets of ranks 1..n with each possible rank sum
    maxS = n * (n + 1) / 2;
    c = zeros(1, maxS + 1);  c(1) = 1;            % c(s+1) = #subsets with sum s
    for k = 1:n
        c(k+1:end) = c(k+1:end) + c(1:end-k);
    end
    prob = c / 2^n;
    lo = sum(prob(1:wPlus + 1));                  % P(W <= wPlus)
    hi = sum(prob(wPlus + 1:end));                % P(W >= wPlus)
    p  = min(1, 2 * min(lo, hi));
else
    [~, ~, g] = unique(abs(d));
    t  = accumarray(g, 1);
    mu = n * (n + 1) / 4;
    s2 = n * (n + 1) * (2 * n + 1) / 24 - sum(t.^3 - t) / 48;
    z  = (wPlus - mu) / sqrt(s2);
    p  = erfc(abs(z) / sqrt(2));                  % = 2 * (1 - normcdf(|z|))
end
end

function [r, hasTies] = tiedRanks(a)
% average ranks for ties
[s, idx] = sort(a);
n = numel(a);
r = zeros(n, 1);
i = 1;
hasTies = false;
while i <= n
    j = i;
    while j < n && s(j + 1) == s(i), j = j + 1; end
    r(idx(i:j)) = (i + j) / 2;
    hasTies = hasTies || j > i;
    i = j + 1;
end
end
