%% LT 11/4/17 - rough analyses on populations


%% go thru all chans/clusts, extract data
clear all; close all;

% Batchname = 'BatchTest'; % must all have used the same batch
% ChansToGet = [9 14 17 18 21];
% ClustToGet = [1 1 1 1 1]; % should be aligned with ChansToGet
% BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA'};
% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/NOTSINGING';

Batchname = 'Batch1219to1304'; % must all have used the same batch
ChansToGet = [9 14 17 18 21];
ClustToGet = [1 1 1 1 1]; % should be aligned with ChansToGet
BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA'};
dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/NOTSINGING';

% Batchname = 'Batch1208to2205'; % must all have used the same batch
% ChansToGet = [17 18 21];
% ClustToGet = [1 1 1]; % should be aligned with ChansToGet
% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/103017_RAlearn1/';

% -- defaults
extractsound = 0;

NeurDat = struct;
numsampsAll = []; % to confirm that all data are aligned (samme num samps);
fsAll = [];
filenamesAll = {};

for i = 1:length(ChansToGet)
    cd(dirname);
    chan = ChansToGet(i);
    clust = ClustToGet(i);
    
    [SongDatTMP, NeurDatTMP, ParamsTMP] = lt_neural_ExtractDat(Batchname, chan, ...
        extractsound, clust);
    
    % - convert from ms to sec
    NeurDatTMP.spikes_cat.cluster_class(:,2) = ...
    NeurDatTMP.spikes_cat.cluster_class(:,2)./1000;
        
    
    % -
    NeurDat.metaDat = NeurDatTMP.metaDat;
    NeurDat.unit(i).spikes_cat = NeurDatTMP.spikes_cat;
    NeurDat.unit(i).chan = chan;
    NeurDat.unit(i).clust = clust;
    NeurDat.dirname = dirname;
    NeurDat.unit(i).BrainRegion = BrainRegion{i};
    
    % ====== sanity checks - all data are from same files
    numsampsAll = [numsampsAll; [NeurDatTMP.metaDat.numSamps]];
    fsAll = [fsAll; [NeurDatTMP.metaDat.fs]];
    filenamesAll = [filenamesAll; [NeurDatTMP.metaDat.filename]];
end

% === sanity checks
assert(size(unique(numsampsAll, 'rows'),1) ==1, 'problem');
assert(size(unique(fsAll, 'rows'),1)==1, 'problem');
assert(numel(unique(filenamesAll)) ==1, 'problem');

%%

NumUnits = length(NeurDat.unit);

%% PLOT ALL CHANNELS

% ============== SPIKES
lt_figure; hold on;
title('all simult. units');

chantoplot = [9]; % empty for all;

if ~isempty(chantoplot)
NeurNums = find(ChansToGet==chantoplot);
else
    NeurNums = 1:length(ChansToGet);
end

for i=NeurNums
    
    spktimes = NeurDat.unit(i).spikes_cat.cluster_class(:,2);
    spktimes = spktimes;
    if strcmp(NeurDat.unit(i).BrainRegion, 'LMAN')
        plotcol = [0.2 0.4 0.4];
    elseif strcmp(NeurDat.unit(i).BrainRegion, 'RA')
        plotcol = [0.8 0.3 0.3];
    end
    lt_neural_PLOT_rasterline(spktimes, i, plotcol);    
end


% =============== VOLTAGE (filtered)
lt_neural_PLOT_RawMultChan(NeurDat.dirname, Batchname, ChansToGet);


%% =========================== PAIRWISE CROSS-COR BETWEEN ALL PAIRS OF CHANS

% ============= 1) CROSS-CORRELATION OF BINNED SPIKES
binsize = 0.005; % in sec. shift will be identical
windowmax = 1; % in sec, extent of cross-corr (+/- windowmax)

% ----- 1) bin all spike trains
fs = NeurDat.metaDat(1).fs;
nsamps = NeurDat.metaDat(1).numSamps;
trialdur = nsamps/fs; % sec
binedges = 0:binsize:trialdur;

FRmat = nan(length(binedges)-1, length(NeurDat.unit)); % timebin x channel
for i=1:length(NeurDat.unit)
    
    x = histc(NeurDat.unit(i).spikes_cat.cluster_class(:,2), binedges); % binned spikes   
    
    FRmat(:,i) = x(1:end-1);
end

% ----- 2) calcualte all cross-corrs
lt_figure; hold on;
% [CCall, lags] = xcorr(FRmat, ceil(windowmax/binsize), 'coeff');
% [CCall, lags] = xcorr(FRmat, ceil(windowmax/binsize), 'unbiased');
[CCall, lags] = xcov(FRmat, ceil(windowmax/binsize), 'coeff');

NeurID1 = [];
NeurID2 = [];
CCall = [];
for i=1:length(NeurDat.unit)
    for ii=i:length(NeurDat.unit)
   
        [cc, lags] = xcov(FRmat(:,i), FRmat(:,ii), ceil(windowmax/binsize), 'coeff');

        CCall = [CCall cc];
        NeurID1 = [NeurID1 i];
        NeurID2 = [NeurID2 ii];
    end
end

    
% ------ 3) plot CC, separating by type of pair
figcount=1;
subplotrows=8;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ------------ a) autocorrelations
inds = find(NeurID1==NeurID2);
plotcol = 'k';

for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    
    
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2))];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
    ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    line([0 0], ylim, 'LineStyle', '--');
end



% -------------- b) same region
inds = find(strcmp(BrainRegion(NeurID1), BrainRegion(NeurID2)) & NeurID1~=NeurID2);
plotcol = 'b';

for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2)) '; ' BrainRegion{neur1} '-' BrainRegion{neur2}];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
    ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    line([0 0], ylim, 'LineStyle', '--');
end



% -------------- b) diff region
inds = find(~strcmp(BrainRegion(NeurID1), BrainRegion(NeurID2)));
plotcol = 'r';

ccall = []; % for plotting mean;
for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2)) '; ' BrainRegion{neur1} '-' BrainRegion{neur2}];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
    ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    
    % --- collect
    ccall = [ccall CCall(:,indtmp)];
    line([0 0], ylim, 'LineStyle', '--');
end
%-- plot mean
ccmean = mean(ccall,2);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean');
hsplots = [hsplots hsplot];

plot(lags.*binsize, ccmean, 'Color', plotcol);
    xlim([-windowmax windowmax]);
    ylim([-0.3 0.8]);
        line([0 0], ylim, 'LineStyle', '--');


linkaxes(hsplots, 'xy');


%% ========= calculate coherency between channels

chan1 = 9;

% load file
pj_readIntanNoGui('pu69wh78_171101_102913.rhd');

%%

indtmp = [amplifier_channels.chip_channel] == chan1;
x1 = amplifier_data(indtmp, :)';


params = struct;
params.Fs = 30000;
params.fpass = [0 1000];
% params.err = [2 0.05];
% [S, f, Serr] = mtspectrumc(x1, params);
[S, f] = mtspectrumc(x1);

lt_figure; hold on;
plot(f, S, '-ok');
% shadedErrorBar(f, S, Serr, {'Color', 'k'},1);


%% ========================================== overnight analysis
batchf = 'batchtmp';
chantype = 'ana'; % dig or ana
chan = 2; % 0, 1, ...
threshold = 1; % voltage, or rise and fall detection.
% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchf);

fline = fgetl(fid);
cumusamps = 0; 
rtAll = [];
ftAll = [];
datAll = [];
while ischar(fline)

    % ==================== extract signal
    [~,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels] = pj_readIntanNoGui(fline);

        
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==chan;
        dat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==chan;
        dat = board_adc_data(ind, :);
    end
    
    
    
    % --- detect threshold crossings
    rt = risetime(dat, 'StateLevels', [0 threshold]);
    ft = falltime(dat, 'StateLevels', [0 threshold]);
    
    % --- rise and fall times as cumu time
    rt = rt+cumusamps;
    ft = ft+cumusamps;
    
    % -- collect across all files
    rtAll = [rtAll rt];
    ftAll = [ftAll ft];
    
    if (1)
    datAll = [datAll dat];
    end
    
    % --- cumulateive time up to now
    cumusamps = cumusamps + length(dat);
    
    fline = fgetl(fid);
end



%% ############### plot batch songs, compare dig pulse with song
close all;

batchf = 'batchall';
chantype = 'ana'; % dig or ana
chan = 0; % 0, 1, ... % for pulse
audchan = 1;
threshold = 1.5; % voltage, or rise and fall detection.
% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchf);

fline = fgetl(fid);
cumusamps = 0; 
rtAll = [];
ftAll = [];
datAll = [];

figcount=1;
subplotrows=8;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


hsplots = [];
while ischar(fline)

    % ==================== extract signal
    [~,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels] = pj_readIntanNoGui(fline);

        fs = frequency_parameters.amplifier_sample_rate;
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==chan;
        dat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==chan;
        dat = board_adc_data(ind, :);
    end

    songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
    
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
hsplots = [hsplots hsplot];
title(fline);

t= [1:length(songdat)]./fs;
plot(t, songdat,'k');
plot(t, dat, 'r');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
hsplots = [hsplots hsplot];

lt_plot_spectrogram(songdat, fs, 0, 0);

%     % --- detect threshold crossings
%     rt = risetime(dat, 'StateLevels', [0 threshold]);
%     ft = falltime(dat, 'StateLevels', [0 threshold]);
%     
%     % --- rise and fall times as cumu time
%     rt = rt+cumusamps;
%     ft = ft+cumusamps;
%     
%     % -- collect across all files
%     rtAll = [rtAll rt];
%     ftAll = [ftAll ft];
%     
%     if (1)
%     datAll = [datAll dat];
%     end
%     
%     % --- cumulateive time up to now
%     cumusamps = cumusamps + length(dat);
%     
linkaxes(hsplots, 'x');
    fline = fgetl(fid);
    hsplots = [];
end
