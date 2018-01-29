function [segextract, xtimes] = lt_neural_QUICK_SpkBinned(segextract, ...
    maxdur, binsize)
%% lt 1/17/18 - takes spktime and outputs binned (1ms res)

if ~exist('binsize', 'var')
    binsize = 0.001;
end

if ~exist('maxdur', 'var')
    maxdur = []; % in sec, max dur to get bins
end

%% NOTE if TIME WARPED:
% 1) bins start at same time across trials (i.e. alignment point - motif predur). throws out spikes that occur before this
% 2) bins end at last spike across trials - not perfect, since there might be silence past that ...

%% -- run

if isempty(maxdur)
% --- what is maximum common trial dur?
stimes = {segextract.spk_Times};
maxtime = cellfun(@max, stimes);
maxtime = max(maxtime);
else
    maxtime = maxdur;
end

% --- initiate xbins
% xedges = 0:0.001:ceil(maxtime*1000)/1000;
% maxtime-mod(maxtime, binsize)
xedges = 0:binsize:maxtime;
xtimes = (xedges(2)-binsize/2):binsize:(xedges(end)-binsize/2);

%% for each trial get counts
ntrials = length(segextract);
SpkCounts = nan(ntrials, length(xedges)-1);
for t =1:ntrials
    y = histc(segextract(t).spk_Times, xedges);
    y = y(1:end-1);
    
    y = int8(y);
    
    SpkCounts(t,:) = y;
    segextract(t).spk_Binned = y';
    segextract(t).spk_Binned_x = xtimes';
end

tmp = isnan(SpkCounts);
assert(~any(isnan(tmp(:))), 'asdfas');

%%

