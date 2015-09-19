function zc = zero_crossings(x, window)

% Returns zero crossings per sample
% Make sure to normalize window if you want an averaging vs. a
% summing operation.

zc(1) = 0;
zc(2:length(x)) = abs(diff(mysign(x)))/2;
zc = conv(window,zc);
% Result will have to be shifted for
% the time delay of the window function and truncated.


