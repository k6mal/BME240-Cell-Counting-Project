function I = preprocessImage(I, p)
%PREPROCESSIMAGE  Steps 1-3: normalize, top-hat background removal, noise filter.
I = normalizeImage(I);

% ---- 2. Background removal (white top-hat) ------------------------------
if p.useTophat
    % strel('disk', R, N): N = 4 (default) is a fast periodic-line
    % approximation, N = 0 is the exact disk (see compareStrelDisk.m)
    se = strel('disk', round(p.tophatFactor * p.cellRadius), p.tophatStrelN);
    I  = mat2gray(imtophat(I, se));
end

% ---- 3. Noise filter ----------------------------------------------------
switch lower(p.filterType)
    case 'none'
    case 'gaussian'
        I = imgaussfilt(I, p.sigma, 'Padding', 'replicate');   % kernel = 2*ceil(2*sigma)+1
    case 'median'
        I = medfilt2(I, [p.medianSize p.medianSize], 'symmetric');
    otherwise
        error('preprocessImage:filter', 'Unknown filterType "%s".', p.filterType);
end
end
