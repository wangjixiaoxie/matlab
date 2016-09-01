%% LT 8/18/16 - takes rsults from wave clus spike sorting and aligns results to raw dat
function lt_neural_AlgnWavclus(batchf, channel_board, plotcols)

% batchf='batch_test'
% channel_board=14;
% NOTE: these have to match exactly the batch name and channel used for the
% lt_neural_concatOneChan function which did the concat before waveclus.

plotRawSound=0;

%% go to data folder
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

%% load, go back up, and plot.
% 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
load('MetaDat.mat'); % contains file names
% --- load concat neural and spikes
neural_cat=load('data.mat');
spikes_cat=load('times_data.mat');


AllSongs_old=[];
AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
cumulative_filedur=0; % keeps track as concatenating

cd ..
if isfield(metaDat, 'songDat');
    % audio was saved before
    
    AllSongs_old=[metaDat.songDat];
    AllNeural_old=neural_cat.data;
    
    for i=1:length(metaDat)
        % -- load labels, onsets, offsets (if exists)
        if exist([metaDat(i).filename '.not.mat'])
            tmp=load([metaDat(i).filename '.not.mat']);
            AllLabels=[AllLabels tmp.labels];
            
            % convert onsets to onset relative to start of entire concatenated file
            onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
            offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
            
            AllOnsets=[AllOnsets onsets_cum'];
            AllOffsets=[AllOffsets offsets_cum'];
        end
        
        % duration of this song file (sec)
        filedur=metaDat(i).numSamps/metaDat(i).fs;
        cumulative_filedur=cumulative_filedur + filedur;
    end
else
    % --- load audio and old neural (from actual files)
    
    
    for i=1:length(metaDat)
        
        % -- load original sound and neural
        [amplifier_data,~,frequency_parameters, board_adc_data, ...
            board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(i).filename);
        
        ind=find([amplifier_channels.chip_channel]==channel_board);
        
        AllSongs_old=[AllSongs_old board_adc_data(1,:)];
        AllNeural_old=[AllNeural_old amplifier_data(ind, :)];
        
        % -- load labels, onsets, offsets
        tmp=load([metaDat(i).filename '.not.mat']);
        AllLabels=[AllLabels tmp.labels];
        
        % convert onsets to onset relative to start of entire concatenated file
        onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
        offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
        
        AllOnsets=[AllOnsets onsets_cum'];
        AllOffsets=[AllOffsets offsets_cum'];
        
        % duration of this song file (sec)
        filedur=metaDat(i).numSamps/metaDat(i).fs;
        cumulative_filedur=cumulative_filedur + filedur;
    end
end
tt=[1:length(AllNeural_old)]/metaDat(1).fs;

%% == comfirm that old and new neural are the same
% assert(all(AllNeural_old==neural_cat.data), 'concat neural and old neural not identical!')
assert(length(AllNeural_old)==length(AllNeural_old), 'PROBLEM!!');

%% PLOT (segment to many figs so doesn't get too large
fs=metaDat(1).fs;
maxdur_fig=90; % sec
totaldur=length(AllSongs_old)/fs; % sec
numfigs=ceil(totaldur/maxdur_fig);

for mm=1:numfigs
    
    sampfirst=1+((mm-1)*maxdur_fig*fs);
    if mm==numfigs
        samplast=length(AllSongs_old);
    else
        samplast=((mm)*maxdur_fig*fs);
    end
    
    indsamp=sampfirst:samplast;
    
    timefirst=(mm-1)*maxdur_fig;
    timelast=(mm)*maxdur_fig;
    
    % === figure out which files are within this chunk
    startsongind=min(find(timefirst*fs<cumsum([metaDat.numSamps]))); % earliest song that ends within this window (i.e. start song)
    stargsongname=metaDat(startsongind).filename;
    
    % 2) PLOT ALL ALIGNED
    figcount=1;
    if plotRawSound==1
    subplotrows=3;
    else
    subplotrows=2;
    end
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots=[];
    
    % - a. song raw
    if plotRawSound==1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    plot(tt(indsamp)-timefirst, AllSongs_old(indsamp));
    hsplots=[hsplots hsplot];
    end
    
    % - b. song spec
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_plot_spectrogram(AllSongs_old(indsamp), metaDat(1).fs, '', 0)
    title(['first song: ' stargsongname]);
    
    hsplots=[hsplots hsplot];
    
    % onsets and labels
    inds_ons_off=AllOnsets>timefirst & AllOnsets<timelast;
    AllOnsets_sub=AllOnsets(inds_ons_off)-timefirst;
    AllOffsets_sub=AllOffsets(inds_ons_off)-timefirst;
    AllLabels_sub=AllLabels(inds_ons_off);
    for i=1:length(AllOnsets_sub)
        line([AllOnsets_sub(i) AllOffsets_sub(i)], [0 0], 'LineWidth', 2);
        lt_plot_text(AllOnsets_sub(i), 100, AllLabels_sub(i), 'r')
    end
    
    %     % - c. neural (old)
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     [datfilt] =lt_neural_filter(AllNeural_old(indsamp), frequency_parameters);
    %     plot(tt(indsamp)-timefirst, datfilt, 'r');
    %     hsplots=[hsplots hsplot];
    %     title('neural, actual dat files');
    %
    % - d. neural (cat)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    [datfilt] =lt_neural_filter(neural_cat.data(indsamp), metaDat(1).fs);
    plot(tt(indsamp)-timefirst, datfilt, 'k');
    hsplots=[hsplots hsplot];
    title('neural, concated file (extracted spikes)');
    
    % overlay spike times (cat)
    numclust=max(spikes_cat.cluster_class(:,1));
    for i=1:numclust
        inds=find([spikes_cat.cluster_class(:,1)]==i & ...
            spikes_cat.cluster_class(:,2)>timefirst*1000 ...
            & spikes_cat.cluster_class(:,2)<timelast*1000);
        
        for iii=1:length(inds)
            ii=inds(iii);
            spktime=spikes_cat.cluster_class(ii, 2); % in ms
            spktime=spktime-1000*timefirst;
            line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
        end
    end
    
    % --- link all
    linkaxes(hsplots, 'x');
end

% - e. plot waveforms for all spikes (divide into quartiles in when they
% occured)
lt_figure; hold on;
numtotspk=size(spikes_cat.spikes,1);
numchunks=10; % to divide spikes by time
binsize=numtotspk/numchunks;
hsplots2=[];
for k=1:numchunks
    hsplot2=lt_subplot(3, ceil(numchunks/3), k); hold on;
    title(['timesegment ' num2str(k)]);
    hsplots2=[hsplots2 hsplot2];
    
    minind=ceil((k-1)*binsize+1);
    maxind=floor(k*binsize);
    spkwavs=spikes_cat.spikes(minind:maxind, :);
    spkclusts=spikes_cat.cluster_class(minind:maxind, 1);
    for j=1:numclust
        
        plot(spkwavs((spkclusts==j),:)', 'color',plotcols{j},'LineStyle','--');
    end
end
linkaxes(hsplots2, 'xy');

