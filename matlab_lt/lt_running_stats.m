%% LT 3/14/15 - Running stats, but windowsize does not shrink at edges (so reduces num samples);
function RunStats=lt_running_stats(Y,BinSize, onlystd);
% Y is 1-D array of numbers
% BinSize is in samples
% RunStats is output structure
% Returns mean/median of all data if samples smaller than binsize

if ~exist('onlystd', 'var')
    onlystd=0;
end

RunStats=struct;

%% RUN
% Running average - binsize does not change at edges
NumSamps=length(Y);

% get inds
BeginInds=1:NumSamps-BinSize+1; % beginning of bins are 1:LastInd

if NumSamps<BinSize; % not enough samples, just use all sampls
    BeginInds=1;
    BinSize=NumSamps;
end

   
% Get stats
for i=1:length(BeginInds);
    Dat=Y(BeginInds(i):BeginInds(i)+BinSize-1);
    
    if onlystd==0    
    % RUNNING MEDIAN
    RunStats.Median(i)=nanmedian(Dat);
    
    % RUNNING MEAN
    RunStats.Mean(i)=nanmean(Dat);
    
    % RUNNING SEM
    RunStats.SEM(i)=nanstd(Dat)/sqrt(numel(Dat)-1);
    
    % RUNNING PERCENTILES
    RunStats.prctiles_5_30_50_70_95(i,:)=prctile(Dat, [5, 30, 50, 70, 95]);
    end
    
    % RUNNING STD
    RunStats.STD(i)=nanstd(Dat);
    
end


