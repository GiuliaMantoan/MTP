function plotFanBands(ax, dates, Q)
%PLOTFANBANDS  Draw shaded quantile fan bands.
%   dates : 1×T datetime;  Q : T×7 [Q05 Q10 Q25 Q50 Q75 Q90 Q95]
    d = dates(:)';
    bands = { [1,7], [0.85 0.9 1],  '5^{th}–95^{th}'; ...
              [2,6], [0.65 0.8 1],  '10^{th}–90^{th}'; ...
              [3,5], [0.4  0.6 1],  '25^{th}–75^{th}' };
    for b = 1:3
        lo = Q(:, bands{b,1}(1))';  hi = Q(:, bands{b,1}(2))';
        fill(ax, [d, fliplr(d)], [lo, fliplr(hi)], bands{b,2}, ...
             'EdgeColor','none','FaceAlpha',1,'DisplayName',bands{b,3});
    end
end