function I = normalizeImage(I)
%NORMALIZEIMAGE  Step 1: grayscale double, min-max scaled to [0, 1].
if size(I, 3) == 3, I = rgb2gray(I); end
I = mat2gray(double(I));
end
