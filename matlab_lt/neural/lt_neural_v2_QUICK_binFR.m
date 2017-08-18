function [Xfr_out, xtimes_out] = lt_neural_v2_QUICK_binFR(Xfr, xtimes, bininfo, TrimDown)
%% lt 8/11/17 - takes matrix of FR vectors and downsamples (nonoverlapping bins)

% Xfr, matrix of fr vectors, trial x time
% xtimes, actual times for each bin, 1 x time
% bininfo
%     if negative: numbins, number of time bins to downsample to.
%     if positive: size of bin in seconds
% TrimDown = 1, then will trim data from flanks to get a data size that
% give whole number of bins.

if ~exist('TrimDown', 'var')
    TrimDown = 0;
end
    

%%

if bininfo >0
    binsizeinds = round(bininfo/(xtimes(2)-xtimes(1)));
else
    numbins = -bininfo;
    binsizeinds = round(size(Xfr, 2)/numbins);
end


if TrimDown ==1 & mod(size(Xfr, 2), binsizeinds)~=0
    trimval = mod(size(Xfr, 2), binsizeinds);
    Xfr = Xfr(:, [1+floor(trimval/2):end-ceil(trimval/2)]);
    xtimes = xtimes([1+floor(trimval/2):end-ceil(trimval/2)]);
    
    % -- redo to get new binsizeinds
    if bininfo >0
        binsizeinds = round(bininfo/(xtimes(2)-xtimes(1)));
    else
        numbins = -bininfo;
        binsizeinds = round(size(Xfr, 2)/numbins);
    end
    
end

assert(mod(size(Xfr, 2), binsizeinds) ==0, 'problem! cannot evenly divide by numbins ...');    %     numbins = floor(size(X,2)/binsizeinds)
% 
% 
% %     assert(mod(size(Xfr, 2), numbins) ==0, 'problem! cannot evenly divide by numbins ...');
% 
% assert(binsizeinds == floor(binsizeinds), 'problem: nonintiger binsize ...');

%%

Xfr_out = reshape(Xfr, size(Xfr,1), binsizeinds, []);
Xfr_out = squeeze(mean(Xfr_out, 2));

xtimes_out = reshape(xtimes, size(xtimes,1), binsizeinds, []);
xtimes_out = squeeze(mean(xtimes_out, 2));
xtimes_out = xtimes_out';


%%  for debuging - plotting new X over old X
if (0)
    % debug
    lt_figure; hold on;
    for zz = 1:size(X,1);
        plot(xtimes, X(zz, :), '-');
        plot(xtmp, tmp(zz, :), '-ok');
        pause
        cla
    end
end

