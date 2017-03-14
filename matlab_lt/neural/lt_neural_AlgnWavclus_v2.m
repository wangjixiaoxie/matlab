%% LT 8/18/16 - takes rsults from wave clus spike sorting and aligns results to raw dat
function lt_neural_AlgnWavclus(batchf, channel_board, plotcols, PlotSecondChan, SecondChan, maxfiguuresopen, figsstepsize)

% batchf='batch_test'
% channel_board=14;
% NOTE: these have to match exactly the batch name and channel used for the
% lt_neural_concatOneChan function which did the concat before waveclus.
% maxfiguuresopen = 1 then will close figs as they open if there are too
% many predicted figs.
plotRawSound=1;

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

AllNeural_SecondChan = [];

cd ..
% if isfield(metaDat, 'songDat') && PlotSecondChan==0
%     % audio was saved before
%
%     AllSongs_old=[metaDat.songDat];
%     AllNeural_old=neural_cat.data;
%
%     for i=1:length(metaDat)
%         % -- load labels, onsets, offsets (if exists)
%         if exist([metaDat(i).filename '.not.mat'])
%             tmp=load([metaDat(i).filename '.not.mat']);
%             AllLabels=[AllLabels tmp.labels];
%
%             % convert onsets to onset relative to start of entire concatenated file
%             onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
%             offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
%
%             AllOnsets=[AllOnsets onsets_cum'];
%             AllOffsets=[AllOffsets offsets_cum'];
%         end
%
%         % duration of this song file (sec)
%         filedur=metaDat(i).numSamps/metaDat(i).fs;
%         cumulative_filedur=cumulative_filedur + filedur;
%     end
% else
%     % --- load audio and old neural (from actual files)
%     for i=1:length(metaDat)
%
%         % -- load original sound and neural
%         [amplifier_data,~,frequency_parameters, board_adc_data, ...
%             board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(i).filename);
%
%         ind=find([amplifier_channels.chip_channel]==channel_board);
%
%         AllSongs_old=[AllSongs_old board_adc_data(1,:)];
%         AllNeural_old=[AllNeural_old amplifier_data(ind, :)];
%
%         % -- load labels, onsets, offsets
%         tmp=load([metaDat(i).filename '.not.mat']);
%         AllLabels=[AllLabels tmp.labels];
%
%         % convert onsets to onset relative to start of entire concatenated file
%         onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
%         offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
%
%         AllOnsets=[AllOnsets onsets_cum'];
%         AllOffsets=[AllOffsets offsets_cum'];
%
%         % duration of this song file (sec)
%         filedur=metaDat(i).numSamps/metaDat(i).fs;
%         cumulative_filedur=cumulative_filedur + filedur;
%
%         % =========== IF WANT SECOND CHAN, COLLECT THAT
%         if PlotSecondChan==1
%             ind=find([amplifier_channels.chip_channel]==SecondChan);
%             AllNeural_SecondChan=[AllNeural_SecondChan amplifier_data(ind, :)];
%         end
%     end
% end
% tt=[1:length(AllNeural_old)]/metaDat(1).fs;



%% ==

cumulativetime_ms = 0;
% --- load audio and old neural (from actual files)
for mm=1:length(metaDat)
    hsplots = [];
    % -- load original sound and neural
    [amplifier_data,~,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(mm).filename);
    
    indsamp=find([amplifier_channels.chip_channel]==channel_board);
    neurdat = amplifier_data(indsamp, :);
    
    % -- load labels, onsets, offsets
    tmp=load([metaDat(mm).filename '.not.mat']);
    AllLabels=[AllLabels tmp.labels];
    
    tt = [1:length(board_adc_data(1,:))]/frequency_parameters.amplifier_sample_rate;
    
    datdur_ms = 1000*length(board_adc_data(1,:))/frequency_parameters.amplifier_sample_rate;
    
    % === PLOT
    figcount=1;
    subplotrows=4;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    
    % - a. song raw
    if plotRawSound==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        plot(tt, board_adc_data(1,:));
        hsplots=[hsplots hsplot];
        
    else
        %     % - b. song spec
        %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     lt_plot_spectrogram(AllSongs_old(indsamp), metaDat(1).fs, '', 0)
        %     title(['first song: ' stargsongname]);
        %     hsplots=[hsplots hsplot];
        
    end
    
    % onsets and labels
    for i=1:length(tmp.onsets)
        line([tmp.onsets(i) tmp.offsets(i)]/1000, [1.65 1.65], 'LineWidth', 2, 'Color', 'r');
        lt_plot_text(tmp.onsets(i)/1000, 1.5, tmp.labels(i), 'r')
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
    [datfilt] =lt_neural_filter(neurdat, frequency_parameters.amplifier_sample_rate);
    plot(tt, datfilt, 'k');
    hsplots=[hsplots hsplot];
    title('neural, concated file (extracted spikes)');
    
    % overlay spike times (cat)
    numclust=max(spikes_cat.cluster_class(:,1));
    for i=1:numclust
        inds=find([spikes_cat.cluster_class(:,1)]==i & ...
            spikes_cat.cluster_class(:,2)>cumulativetime_ms ...
            & spikes_cat.cluster_class(:,2)<=cumulativetime_ms+datdur_ms);
        
        for iii=1:length(inds)
            ii=inds(iii);
            spktime=spikes_cat.cluster_class(ii, 2); % in ms
            spktime=spktime-cumulativetime_ms;
            line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
        end
    end
    
    if PlotSecondChan==1
        ind =find([amplifier_channels.chip_channel]==SecondChan);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        [datfilt] =lt_neural_filter(amplifier_data(ind, :), frequency_parameters.amplifier_sample_rate);
        plot(tt, datfilt, 'k');
        hsplots=[hsplots hsplot];
        title(['second chan (' num2str(SecondChan) ')']);
        
    end
    
    % --- link all
    linkaxes(hsplots, 'x');
    
    %         if PlotSecondChan==1
    %             ind=find([amplifier_channels.chip_channel]==SecondChan);
    %             AllNeural_SecondChan=[AllNeural_SecondChan amplifier_data(ind, :)];
    %         end
    cumulativetime_ms = cumulativetime_ms + datdur_ms;
    
    lt_subtitle(['song #' num2str(mm) '/' num2str(length(metaDat)) '(' metaDat(mm).filename ')'])
    
    tmp = input('type anything to go to next file', 's');
    
    if strcmp(tmp, 'm')
        mm = mm+4;
    end
    
    if length(metaDat) > maxfiguuresopen
        close all;
    end
end

%% - e. plot waveforms for all spikes (divide into quartiles in when they
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

