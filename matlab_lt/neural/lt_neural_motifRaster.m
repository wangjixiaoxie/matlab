%% lt 8/17/16 - plot raster with syls extracted
% takes regexpr, and plots rasters aligned to onset of the first letter in
% the regexpr string. 

function [Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs)

% % regexpr_str='ghh';
% % predur=4; % sec
% % postdur=4; % sec
% % SongDat, NeurDat, Params, shoudl correspond to same dataset, extracted using lt_neural_ExtractDat
% 
% 
% %% === extract stuff from inputs
% 
% 
% AllSongs=SongDat.AllSongs;
% AllLabels=SongDat.AllLabels;
% AllOnsets=SongDat.AllOnsets;
% AllOffsets=SongDat.AllOffsets;
% spikes_cat=NeurDat.spikes_cat;
% % metaDat=NeurDat.metaDat;
% 
% fs=NeurDat.metaDat(1).fs;
% %% ------ find motifs
% [startinds, endinds, matchlabs]=regexp(AllLabels, regexpr_str, 'start', 'end', 'match');
% 
% SegmentsExtract=struct;
% 
% % -- for each match ind, extract audio + spikes
% for i=1:length(startinds)
%     
%     % on time
%     ind=startinds(i);    
%     ontime=AllOnsets(ind); % sec
%     ontime=ontime-predur; % - adjust on time
%     onsamp=round(ontime*fs);    
%     
%     % off time
%     ind=endinds(i);
%     offtime=AllOffsets(ind);
%     offtime=offtime+postdur;
%     offsamp=round(offtime*fs);
%     
%     % - collect
%     songseg=AllSongs(onsamp:offsamp);
%     
%     spkinds=(spikes_cat.cluster_class(:,2) > ontime*1000) & ...
%         (spikes_cat.cluster_class(:,2) < offtime*1000);
%     spk_ClustTimes = spikes_cat.cluster_class(spkinds, :); % in sec, relative to onset of the segment
% 
%     SegmentsExtract(i).songdat=songseg;
%     SegmentsExtract(i).spk_Clust=spk_ClustTimes(:,1)';
%     SegmentsExtract(i).spk_Times=(spk_ClustTimes(:,2)/1000)'-ontime;
%     SegmentsExtract(i).global_ontime=ontime;
%     SegmentsExtract(i).global_offtime=offtime;
%     SegmentsExtract(i).matchlabel=matchlabs{i};
%     SegmentsExtract(i).fs=fs;
% end
% 

fs=SegmentsExtract(1).fs;

predur=Params.REGEXP.predur;

    plotcols={'k', 'r', 'b', 'm'}; % each cluster.

%% =============== plot all song segments and associated spikes
if plotAllSegs==1
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
    title(['segment starts at ' num2str(SegmentsExtract(i).global_ontime_motifInclFlank) 'sec']);
    
    Ylim=ylim;
    
    % - spike times
    numclust=max(SegmentsExtract(i).spk_Clust);
    for j=1:numclust
        inds=find([SegmentsExtract(i).spk_Clust]==j);
        
        if isempty(inds)
            continue
        else
            for iii=1:length(inds)
                jj=inds(iii);
                if useRescaled==0
                spktime=SegmentsExtract(i).spk_Times(jj); 
                else
                 spktime=SegmentsExtract(i).spk_Times_scaled(jj);
                end
                line([spktime spktime], [0 -500]-j*500, 'Color', plotcols{j}, 'LineWidth',2);
            end
        end
        
        % -- overlay line of median dur if this was using scaled
        if useRescaled==1
            medianmotifdur=Params.LINTIMEWARP.motiftime_median;
            predur=Params.REGEXP.predur;
            line([predur predur+medianmotifdur], [0 0], 'Color', 'g', 'LineWidth', 2);
            lt_plot_text(predur, 100, 'rescaled', 'g');
        end
    end
    
    ylim([Ylim(1)-(j+2)*500 Ylim(2)]);
    
end

linkaxes(hsplots, 'xy');
end



%% ============ align and plot all neural and song segments
numsegs=length(SegmentsExtract);
numclusters=max([SegmentsExtract.spk_Clust]);

figcount=1;
subplotrows=4;
subplotcols=numclusters+1;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];

predur=Params.REGEXP.predur;
postdur=Params.REGEXP.postdur;


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
    
    % -- overlay line of median dur if this was using scaled
    if useRescaled==1
        medianmotifdur=Params.LINTIMEWARP.motiftime_median;
        predur=Params.REGEXP.predur;
        line([predur predur+medianmotifdur], [0.1 0.1], 'Color', 'g', 'LineWidth', 2);
        lt_plot_text(predur, 100, 'rescaled', 'g');
    end
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
                if useRescaled==0
                spktime=SegmentsExtract(i).spk_Times(jj); %
                else
                spktime=SegmentsExtract(i).spk_Times_scaled(jj); %                    
                end
                line([spktime spktime], [i-0.4 i+0.4], 'Color', plotcols{k}, 'LineWidth',2);
            end
        end
        
        % -- shade background for this raster, to indicate duration of
        % motif
        if (0) % skipped because distorts rasters
        motifdur=(SegmentsExtract(i).global_offtime_motifInclFlank - postdur) - ...
        (SegmentsExtract(i).global_ontime_motifInclFlank + predur);
        X=[predur predur+motifdur predur+motifdur predur];
        Y=[i-0.5 i-0.5 i+0.5 i+0.5];
        h=patch(X, Y, 'w', 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
        set(h, 'FaceAlpha', 0.2)
        end
    
    end
        % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
    
    % if rescaled, plot vert line for offset
    if useRescaled==1
        medianmotifdur=Params.LINTIMEWARP.motiftime_median;
        predur=Params.REGEXP.predur;
        line([predur+medianmotifdur predur+medianmotifdur], [ylim], 'Color', 'y');
    end
end

% --- mean contours
% NOTE: has not been adjusted to allow plotting lin rescaled spikes!!!!!!!

% clustersall=[SegmentsExtract.spk_Clust];
% spikesall=[SegmentsExtract.spk_Times];
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
window=0.015;
windshift=0.002;
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
        if useRescaled==0
        spktimes=SegmentsExtract(kk).spk_Times(inds);
        else
         spktimes=SegmentsExtract(kk).spk_Times_scaled(inds);           
        end
        Yspks{kk}=spktimes;
    end
    
    % - 
    [xbin, ymean, ysem, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
    % convert to spike rate
    
    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{k}}, 1);
    ylim([-10 150]);
    
    % - plot vertical line for motif onset
    line([predur predur], ylim, 'Color', 'y');
    
        % if rescaled, plot vert line for offset
    if useRescaled==1
        medianmotifdur=Params.LINTIMEWARP.motiftime_median;
        predur=Params.REGEXP.predur;
        line([predur+medianmotifdur predur+medianmotifdur], [ylim], 'Color', 'y');
    end

end

linkaxes(hsplots, 'x');