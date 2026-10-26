function G = loadBBBC039Mask(pngPath)
%LOADBBBC039MASK  Decode a BBBC039 annotation PNG into an instance label image.
%   The red channel stores a small map colouring (values 1, 2, 3): touching
%   nuclei always get different values, so each nucleus is a 4-connected
%   component of a single value. Returns G with 0 = background and
%   1..N = individual nuclei.

M = imread(pngPath);                % RGBA PNG: M is M-by-N-by-3, alpha is a separate output
r = M(:, :, 1);
G = zeros(size(r), 'uint16');
offset = 0;
for v = reshape(unique(r(r > 0)), 1, [])
    Lv = bwlabel(r == v, 4);
    G(Lv > 0) = Lv(Lv > 0) + offset;
    offset = double(max(G(:)));
end
end
