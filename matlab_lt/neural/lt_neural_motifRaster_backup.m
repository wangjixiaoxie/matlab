%% lt 8/17/16 - plot raster with syls extracted
% takes regexpr, and plots rasters aligned to onset of the first letter in
% the regexpr string. 

function [Params]=lt_neural_motifRaster(SongDat, NeurDat, Params, regexpr_str, predur, postdur)

% regexpr_str='ghh';
% predur=4; % sec
% postdur=4; % sec
% channel_board=14; neural
% batchf='BatchTest'; the file used when first concat neural data for
% waveclus. does not have to exist in folder currenrtly. 


% % -- extract labels
% datdir=['Chan' num2str(channel_board) 'amp-' batchf];
% cd(datdir)
% 
% %% 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
% load('MetaDat.mat'); % contains file names
% % --- load concat neural and spikes
% % neural_cat=load('data.mat');
% spikes_cat=load('times_data.mat');
% fs=metaDat(1).fs; % assumes all have same fs, as they should.
% 
% % --- load audio and old neural (from actual files)
% cd ..
% 
% AllSongs_old=[];
% % AllNeural_old=[];
% AllLabels=[];
% AllOnsets=[];
% AllOffsets=[];
% cumulative_filedur=0; % keeps track as concatenating
% 
% for i=1:length(metaDat)
%     
%     % -- load original sound and neural
%     [amplifier_data,~,frequency_parameters, board_adc_data, ...
%         board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(i).filename);
%    ind=find([amplifier_channels.chip_channel]==channel_board);
%    AllSongs_old=[AllSongs_old board_adc_data(1,:)];
% %    AllNeural_old=[AllNeural_old amplifier_data(ind, :)];
%     
% 
% % -- load labels, onsets, offsets
%     tmp=load([metaDat(i).filename '.not.mat']);
%     AllLabels=[AllLabels tmp.labels];
%     
%     % convert onsets to onset relative to start of entire concatenated file
%     onsets_cum=(tmp.onsets/1000)+cumulative_filedur; % sec
%     offsets_cum=(tmp.offsets/1000)+cumulative_filedur;
%         
%     AllOnsets=[AllOnsets onsets_cum'];
%     AllOffsets=[AllOffsets offsets_cum'];
%     
%     
% % duration of this song file (sec)
% filedur=metaDat(i).numSamps/metaDat(i).fs;
% cumulative_filedur=cumulative_filedur + filedur;
% end
% 
% 
% % ---- sanity check, plot onsets, offsets, cum times
% % figcount=1;
% % subplotrows=4;
% % subplotcols=1;
% % fignums_alreadyused=[];
% % hfigs=[];
% % tt=[1:length(AllSongs_old)]/frequency_parameters.amplifier_sample_rate;
% % hsplots=[];
% % 
% % % - a. song raw
% % [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% % plot(tt, AllSongs_old);
% % hsplots=[hsplots hsplot];
% % 
% % % - b. song spec
% % [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% % lt_plot_spectrogram(AllSongs_old, metaDat(1).fs, '', 0)
% % hsplots=[hsplots hsplot];
% % 
% % % onsets and labels
% % for i=1:length(AllOnsets)
% %     line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
% %     lt_plot_text(AllOnsets(i), 100, AllLabels(i), 'r')
% % end


%% ------ find motifs
[startinds, endinds, matchlabs]=regexp(AllLabels, regexpr_str, 'start', 'end', 'match');

SegmentsExtract=struct;

% -- for each match ind, extract audio + spikes
for i=1:length(startinds)
    
    % on time
    ind=startinds(i);    
    ontime=AllOnsets(ind); % sec
    ontime=ontime-predur; % - adjust on time
    onsamp=round(ontime*fs);    
    
    % off time
    ind=endinds(i);
    offtime=AllOffsets(ind);
    offtime=offtime+postdur;
    offsamp=round(offtime*fs);
    
    % - collect
    songseg=AllSongs_old(onsamp:offsamp);
    
    spkinds=(spikes_cat.cluster_class(:,2) > ontime*1000) & ...
        (spikes_cat.cluster_class(:,2) < offtime*1000);
    spk_ClustTimes = spikes_cat.cluster_class(spkinds, :); % in sec, relative to onset of the segment

    SegmentsExtract(i).songdat=songseg;
    SegmentsExtract(i).spk_Clust=spk_ClustTimes(:,1)';
    SegmentsExtract(i).spk_Times=(spk_ClustTimes(:,2)/1000)'-ontime;
    SegmentsExtract(i).global_ontime=ontime;
    SegmentsExtract(i).global_offtime=offtime;
    SegmentsExtract(i).matchlabel=matchlabs{i};
    SegmentsExtract(i).fs=fs;
end



%% =============== plot all song segments and associated spikes
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];

for i=1:length(SegmentsExtract)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    % - song spec
    datsong=SegmentsExtract(i).songdat;
    lt_plot_spectrogram(datsong, SegmentsExtract(i).fs, 1, 0)
    title(['segment starts at ' num2str(SegmentsExtract(i).global_ontime) 'sec']);
    
    Ylim=ylim;
    
    % - spike times
    numclust=max(SegmentsExtract(i).spk_Clust);
    plotcols={'k', 'r', 'b', 'm'};
    for j=1:numclust
        inds=find([SegmentsExtract(i).spk_Clust]==j);
        
        if isempty(inds)
            continue
        else
            for iii=1:length(inds)
                jj=inds(iii);
                spktime=SegmentsExtract(i).spk_Times(jj); % in ms
                line([spktime spktime], [0 -500]-j*500, 'Color', plotcols{j}, 'LineWidth',2);
            end
        end
    end
    
    ylim([Ylim(1)-(j+2)*500 Ylim(2)]);
    
end

linkaxes(hsplots, 'xy');

%% =========== Time warp songs and neural (assuming that is same song
% segments) 




%% ============ align and plot all neural and song segments
numsegs=length(SegmentsExtract);
numclusters=max(spikes_cat.cluster_class(:,1));

figcount=1;
subplotrows=4;
subplotcols=numclusters+1;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];

    % 1) song spectrogram (1st example song)
for k=1:numclusters+1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    songdat=SegmentsExtract(1).songdat;
    lt_plot_spectrogram(songdat, fs, 1, 0);
    
    % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
end


% ---- overlay sound contours for each rend
for k=1:numclusters+1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    for kk=1:length(SegmentsExtract);
    songdat=SegmentsExtract(kk).songdat;
    
    [songdat_sm, tsm]=lt_plot_AmplSm(songdat, fs);
    plot(tsm, songdat_sm);
    end
    
    % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
end


    % 2) rasters for each rendition
for k=1:numclusters+1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    for i=1:numsegs
        
        if k<numclusters+1;
            % then is still single unit
        inds=find([SegmentsExtract(i).spk_Clust]==k);
        else
            % then mult unit - combine all units
            inds = 1:length(SegmentsExtract(i).spk_Clust);
        end
        
        if isempty(inds)
            continue
        else
            for iii=1:length(inds)
                jj=inds(iii);
                spktime=SegmentsExtract(i).spk_Times(jj); % in ms
                line([spktime spktime], [i-0.4 i+0.4], 'Color', plotcols{k}, 'LineWidth',2);
            end
        end
    end
    
        % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
end

% --- mean contours
clustersall=[SegmentsExtract.spk_Clust];
spikesall=[SegmentsExtract.spk_Times];
% -- HISTOGRAM
% binsize=0.02; % sec
% for k=1:numclusters+1
%     
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title('histogram');
%     hsplots=[hsplots hsplot];
%     
%     if k<clustersall+1
%         % s.u.
%     inds=clustersall==k;
%     else
%         % multiunit
%         inds=1:length(clustersall);
%     end
%     spikes=spikesall(inds);
%     
%     % histogram
%     Xcenters=0:binsize:max(spikesall);
%     lt_plot_histogram(spikes, Xcenters, 1, 0, 1, 0, plotcols{k});
% 
%         % - plot vertical line for motif onset
%     line([predur predur], ylim, 'Color', 'y');
% end
% 
% -- KSD
% sigma=0.01; % sec
% for k=1:numclusters+1
%     
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title('ksd');
%     hsplots=[hsplots hsplot];
%     
%     if k<clustersall+1
%         % s.u.
%     inds=clustersall==k;
%     else
%         % multiunit
%         inds=1:length(clustersall);
%     end
%     spikes=spikesall(inds);
%     
%     xvals=0:0.01:max(spikesall);
%     [f, ~, bw] =ksdensity(spikes, xvals, 'bandwidth', sigma);
%     plot(xvals, f, 'k');
%         
%     % - plot vertical line for motif onset
%     line([predur predur], ylim, 'Color', 'y');
% end



% --- simply overlapping mean +/- se windows
window=0.02;
windshift=0.04;
for k=1:numclusters+1
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean(se)');
    ylabel('sp/sec');
    
    hsplots=[hsplots hsplot];
    
    % - convert this cluster's spikes into cell array
    numsegs=length(SegmentsExtract);
    Yspks={}; % one ind for each seg.
    for kk=1:numsegs
        if k<numclusters+1
            inds=find([SegmentsExtract(kk).spk_Clust]==k);
        else
            inds=1:length(SegmentsExtract(kk).spk_Clust);
        end
        spktimes=SegmentsExtract(kk).spk_Times(inds);
        Yspks{kk}=spktimes;
    end
    
    % - 
    [xbin, ymean, ysem, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
    % convert to spike rate
    
    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{k}}, 1);
    ylim([-10 150]);
    
    % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
end

linkaxes(hsplots, 'x');