function BW = thresholdImage(I, p)
%THRESHOLDIMAGE  Step 4: global Otsu or local-mean adaptive threshold.
%   I must already be scaled to [0, 1] (preprocessImage).
switch lower(p.threshMethod)
    case 'otsu'
        BW = I > graythresh(I);
    case 'adaptive'
        w  = 2 * floor(p.adaptWindowFactor * p.cellRadius) + 1;  % odd window
        BW = I > imboxfilt(I, w, 'Padding', 'replicate') + p.adaptOffset;
    otherwise
        error('thresholdImage:method', 'Unknown threshMethod "%s".', p.threshMethod);
end
end
