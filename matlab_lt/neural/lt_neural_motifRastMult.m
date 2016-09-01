%%  LT 8/22/16 - plots multiple motifs' rasters and overlays mean firing rate (for a given neuron

function lt_neural_motifRastMult(SongDat, NeurDat, Params, MotifList_regexp, ...
    predur, postdur, LinearWarp, suppressplots, onlyPlotMean)

% MotifList_regexp={'g(h)h', 'n(h)h', '(S)[\w-]*?E'};
% % MotifList_regexp={'g(h)h', 'n(h)h'};
%
% predur=2; % will use for all motifs
% postdur=3;
% LinearWarp=2; % 0=no; 1=yes, each motif individually; 2=yes, all motifs to global motif median dur

fs=NeurDat.metaDat(1).fs;
plotcols={'k', 'r', 'b', 'm'}; % each cluster.
numclusts=max(NeurDat.spikes_cat.cluster_class(:,1));

%%

nummotifs=length(MotifList_regexp);


%% == extract data for each motif
SegmentsALL=struct; % will hold all data

for i=1:nummotifs
    
    regexp_str=MotifList_regexp{i};
    disp(['extracting regexp str: ' regexp_str]);
    
    % == extract all data segments
    [sgmnt, prms]=lt_neural_RegExp(SongDat, NeurDat, Params, regexp_str, predur, postdur);
    
    % == linearly time warp (each segment individually)
    % IN PROGRESS
    
    SegmentsALL(i).data=sgmnt;
    SegmentsALL(i).params=prms;
end


%% === LINEARLY TIME WARP
% if do all together, then gets median duration over all motifs (tokensyl
% onset to last syl offset) and uses that as template
if LinearWarp==1
    % then each motif warp to itself
    for i=1:length(SegmentsALL)
        
        disp(['... warping: ' MotifList_regexp{i}]);
        [sgmtdata, prms] = lt_neural_LinTimeWarp(SegmentsALL(i).data, ...
            SegmentsALL(i).params, []);
        
        SegmentsALL(i).data=sgmtdata;
        SegmentsALL(i).params=prms;
    end
    
    
elseif LinearWarp==2
    % warp vs global motif median
    % -- get median over all renditions over all motifs
    OnAll=[];
    OffAll=[];
    for i=1:length(SegmentsALL)
        
        predur=SegmentsALL(i).params.REGEXP.predur;
        postdur=SegmentsALL(i).params.REGEXP.postdur;
        
        motifontimes = predur + ...
            [SegmentsALL(i).data.global_ontime_motifInclFlank];
        
        motifofftimes= [SegmentsALL(i).data.global_offtime_motifInclFlank] - ...
            postdur;
        
        % --- global
        OnAll=[OnAll motifontimes];
        OffAll=[OffAll motifofftimes];
    end
    
    MotifTimesAll=OffAll-OnAll;
    MotifTime_med=median(MotifTimesAll);
    
    % -- plot
    lt_figure; hold on;
    title('all motif durs, across motifs and rends (line=median)');
    lt_plot_histogram(MotifTimesAll);
    line([MotifTime_med MotifTime_med], ylim, 'color','r');
    
    
    
    % ======= PERFORM TIME WARP ON ALL MOTIFS
    for i=1:length(SegmentsALL)
        
        disp(['... warping: ' MotifList_regexp{i}]);
        [sgmtdata, prms] = lt_neural_LinTimeWarp(SegmentsALL(i).data, ...
            SegmentsALL(i).params, MotifTime_med);
        
        SegmentsALL(i).data=sgmtdata;
        SegmentsALL(i).params=prms;
    end
end

if suppressplots==1
    close all
end

%% ==== PLOT RASTERS FOR EACH segment separately
%  one plot for each cluster for this channel
window=0.02;
windshift=0.004;

if onlyPlotMean==0
    hsplots=[];
    
    % raster stuff
    
    for j=1:numclusts
        plotcol=plotcols{j};
        % -- Go thru all motifs
        for i=1:nummotifs
            motif=MotifList_regexp{i};
            datstruct=SegmentsALL(i).data;
            numtrials=length(datstruct);
            predur=SegmentsALL(i).params.REGEXP.predur;
            
            lt_figure; hold on;
            % == 1) spectrogram (representative)
            hsplot=lt_subplot(4,1,1); hold on;
            hsplots=[hsplots hsplot];
            
            %         title('ex. specgram');
            songdat=SegmentsALL(i).data(1).songdat;
            lt_plot_spectrogram(songdat, fs, 1, 0);
            line([predur predur], ylim, 'Color', 'y');
            
            % == 2) overlayed song amplitudes
            hsplot=lt_subplot(8,1,3); hold on;
            hsplots=[hsplots hsplot];
            
            for kk=1:numtrials;
                songdat=datstruct(kk).songdat;
                
                [songdat_sm, tsm]=lt_plot_AmplSm(songdat, fs);
                plot(tsm, songdat_sm, '--');
            end
            line([predur predur], ylim, 'Color', 'y');
            
            % == 3) rasters
            hsplot= lt_subplot(8,1, 4:6); hold on;
            hsplots=[hsplots hsplot];
            
            Yspks={};
            for kk=1:numtrials
                inds=find([datstruct(kk).spk_Clust]==j); % find this cluster's spikes
                if LinearWarp==0
                    % then no warp
                    spktimes=datstruct(kk).spk_Times(inds);
                else
                    spktimes=datstruct(kk).spk_Times_scaled(inds);
                end
                Yspks{kk}=spktimes;
            end
            [xbin, ~, ~, ymean_hz, ysem_hz]= ...
                lt_neural_plotRastMean(Yspks, window, windshift, 1, plotcol);
            line([predur predur], ylim, 'Color', 'y');
            
            % == 4) mean fr
            hsplot=lt_subplot(4, 1, 4); hold on;
            hsplots=[hsplots hsplot];
            
            shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
            ylim([-10 150]);
            line([predur predur], ylim, 'Color', 'y');
            lt_plot_zeroline
            
            lt_subtitle(['clust' num2str(j) 'motif: ' motif]);
            
            linkaxes(hsplots, 'x')
        end
    end
    
    
    % ----- distribute plots evenly
    distFig('Columns', 'Auto')
    
    
    disp('PAUSED');
    pause;
end

%% ============= OVERLAY means for diff motifs for single cluster (cell)

plotcols_motifs=lt_make_plot_colors(nummotifs, 0, 0);


for j=1:numclusts
    hsplots=[];
    lt_figure; hold on;
    
    % -- Go thru all motifs
    for i=1:nummotifs
        datstruct=SegmentsALL(i).data;
        predur=SegmentsALL(i).params.REGEXP.predur;
        numtrials=length(datstruct);
        plotcol=plotcols_motifs{i};
        hsplots=[];
        
        % === 1) plot mean amplitude contour for each
        hsplot=lt_subplot(2,1,1); hold on;
        hsplots=[hsplots hsplot];
        title('overlayed song sm ampl (NOT TIME WARPED!!!)');
        SongDatAll=[];
        for kk=1:numtrials;
            songdat=datstruct(kk).songdat;
            
            [songdat_sm, tsm]=lt_plot_AmplSm(songdat, fs);
            %             SongDatAll=[SongDatAll; songdat_sm];
            plot(tsm, songdat_sm, '-', 'Color', plotcol, 'LineWidth', 0.1);
            
        end
        %         songmean=mean(SongDatAll, 1);
        %         songstd=std(SongDatAll, 0, 1);
        %         shadedErrorBar(tsm, songmean, songstd, {'Color', plotcol}, 1);
        line([predur predur], ylim, 'Color', 'y');
        
        
        % == 4) mean fr
        hsplot=lt_subplot(2,1,2); hold on; hsplots=[hsplots hsplot];
        
        Yspks={};
        for kk=1:numtrials
            inds=find([datstruct(kk).spk_Clust]==j); % find this cluster's spikes
            if LinearWarp==0
                % then no warp
                spktimes=datstruct(kk).spk_Times(inds);
            else
                spktimes=datstruct(kk).spk_Times_scaled(inds);
            end
            Yspks{kk}=spktimes;
        end
        [xbin, ~, ~, ymean_hz, ysem_hz]= ...
            lt_neural_plotRastMean(Yspks, window, windshift, 0, plotcol);
        
        hplot=shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol, 'LineWidth', 2}, 1);
        hsplots=[hsplots hplot.mainLine];
        set(hplot.patch, 'FaceAlpha', 0.1)
        %         set(hplot.edge, 'LineStyle', 'none');
        
        line([predur predur], ylim, 'Color', plotcol);
        
        % ---
        lt_subtitle(['clust' num2str(j)]);
        linkaxes(hsplots, 'x');
        
    end
    ylim([-10 150]);
    lt_plot_zeroline;
    legend(hsplots, MotifList_regexp);
end


