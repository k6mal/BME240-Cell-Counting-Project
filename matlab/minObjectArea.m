function P = minObjectArea(p)
%MINOBJECTAREA  Smallest object kept by the clean-up, in pixels.
P = round(p.minAreaFrac * pi * p.cellRadius^2);
end
