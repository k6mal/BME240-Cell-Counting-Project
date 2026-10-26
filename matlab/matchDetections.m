function m = matchDetections(centroids, G)
%MATCHDETECTIONS  Score detections against an instance label image.
%   A detection is a true positive if its centroid falls inside an annotated
%   nucleus that has not already been matched (one-to-one matching).
%   centroids : N x 2 [x y] from regionprops; G : label image (0 = background)

[H, W]  = size(G);
nTrue   = double(max(G(:)));
matched = false(nTrue, 1);
tp = 0;
for k = 1:size(centroids, 1)
    c = min(max(round(centroids(k, 1)), 1), W);
    r = min(max(round(centroids(k, 2)), 1), H);
    g = G(r, c);
    if g > 0 && ~matched(g)
        matched(g) = true;
        tp = tp + 1;
    end
end
nAuto = size(centroids, 1);
m.nTrue = nTrue;  m.nAuto = nAuto;
m.tp = tp;  m.fp = nAuto - tp;  m.fn = nTrue - tp;
m.precision = tp / max(nAuto, 1);
m.recall    = tp / max(nTrue, 1);
m.f1        = 2 * m.precision * m.recall / max(m.precision + m.recall, eps);
end
