function [MOTIFSTATS, MotifTime_med, spktimefield] = lt_neural_MultNeur_MotifRasters_v2(NeuronDatabase, motif_regexpr_str, motif_predur, motif_postdur, LinScaleGlobal, WHOLEBOUTS_edgedur)
%% outputs

% MotifTime_med = [], unless do global, then gives median across all.
% [currently gives [] if do warp local to a motif/neuron


%% v2 - can input multiple motifs, will align all of them and compare firing rate for each neuron across motifs
% now:
% motif_regexpr_str is cell array, each ind with one motif
suppresssomeplots = 1; % just time warping plots


%% given motif, plots rasters and smooth rate for all neurons

% LinScaleGlobal=2; % 0: NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)

%% PARAMS
NumMotifs=length(motif_regexpr_str);
NumNeurons=length(NeuronDatabase.neurons);

% smoothing window (neural)
window=0.02; windshift=0.004;

%% EXTRACT DATA
MOTIFSTATS=struct;
for i=1:NumNeurons
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
%     dirdate=NeuronDatabase.neurons(i).date;
%     tmp=dir([dirdate '*']);
%     assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
%     cd(tmp(1).name);
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(i).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);

    % - load data for this neuron
    tic
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    toc
    % --- EXTRACT DAT
    % - do this one time for each desired motif
    for j=1:NumMotifs
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        if ~exist('WHOLEBOUTS_edgedur', 'var')
            WHOLEBOUTS_edgedur = ''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        end
        keepRawSongDat = 1;
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, '', keepRawSongDat);
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
    end
end

clear SongDat
clear NeurDat

%% ========================== PLOT COMBINED RASTERS
MotifTime_med = [];
if LinScaleGlobal==1
    % =================== OPTIONAL, LINEARLY TIME WARP ALL TO SAME DURATION[DEFAULT]
    spktimefield='spk_Times_scaled';
    % -- 1) figure out median dur across all neurons (for this motif)
    OnAll=[];
    OffAll=[];
    MotifTimesAll=[];
    
    for i=1:NumNeurons
        for j=1:NumMotifs
            if ~isfield(MOTIFSTATS.neurons(i).motif(j).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            
            %         predur=MOTIFSTATS.neurons(i).Params.REGEXP.predur;
            %         postdur=MOTIFSTATS.neurons(i).Params.REGEXP.postdur;
            
            %         motifontimes = predur + ...
            %             [MOTIFSTATS.neurons(i).SegmentsExtract.global_ontime_motifInclFlank];
            %         motifofftimes=[MOTIFSTATS.neurons(i).SegmentsExtract.global_offtime_motifInclFlank] - ...
            %             postdur;
            
            %         OnAll=[OnAll motifontimes];
            %         OffAll=[OffAll motifofftimes];
            MotifTimesAll=[MotifTimesAll MOTIFSTATS.neurons(i).motif(j).SegmentsExtract.actualmotifdur];
        end
    end
    %     MotifTimesAll=OffAll-OnAll;
    lt_figure; hold on;
    lt_plot_histogram(MotifTimesAll); xlabel('all motif times (s)');
    MotifTime_med=median(MotifTimesAll);
    
    
    % -- 2) time warp
    for i=1:NumNeurons
        for j=1:NumMotifs
            if ~isfield(MOTIFSTATS.neurons(i).motif(j).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            segmentdat=MOTIFSTATS.neurons(i).motif(j).SegmentsExtract;
            prms=MOTIFSTATS.neurons(i).motif(j).Params;
            [segmentdatout, prmsout] = lt_neural_LinTimeWarp(segmentdat, ...
                prms, MotifTime_med);
            
            if suppresssomeplots==1
                close gcf
                close gcf % close 2 figs
            end
            MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=segmentdatout;
            MOTIFSTATS.neurons(i).motif(j).Params=prmsout;
        end
    end
elseif LinScaleGlobal==2
    % ======= TIME WARP WITHIN EACH NEURON AND MOTIF
    spktimefield='spk_Times_scaled';
    
    for i=1:NumNeurons
        for j=1:NumMotifs
            if ~isfield(MOTIFSTATS.neurons(i).motif(j).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            segmentdat=MOTIFSTATS.neurons(i).motif(j).SegmentsExtract;
            prms=MOTIFSTATS.neurons(i).motif(j).Params;
            [segmentdatout, prmsout] = lt_neural_LinTimeWarp(segmentdat, ...
                prms, []);
            
            if suppresssomeplots==1
                close gcf
                close gcf % close 2 figs
            end
            MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=segmentdatout;
            MOTIFSTATS.neurons(i).motif(j).Params=prmsout;
        end
    end
    
elseif LinScaleGlobal==0
    % THEN NO SCALING AT ALLL!
    spktimefield='spk_Times';
end

%% =================== PLOT ONE GLOBAL RASTER [do once for each motif]
hsplots_rasty=[];
for m=1:NumMotifs
    lt_figure; hold on;
    ylabel('trial #');
    xlabel('sec');
    yind=1; % will increment to plot over all experiments
    plotcols=lt_make_plot_colors(NumNeurons, 0, 0);
    hsplots=[];
    
    % 1) === plot spectrogram at top
    hsplot=lt_subplot(7, 1, 1); hold on
    title(motif_regexpr_str{m});
    hsplots=[hsplots hsplot];
    songdat=[];
    fs=[];
    if LinScaleGlobal==1
        % then try to find median duration trial to plot - otherwise plot
        % random song - see below.
        ylabel('median dur trial');
        % --- find the trial that has close to median motif time
        for i=1:NumNeurons
            if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            
            numtrials=length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract);
            
            for j=1:numtrials
                if abs(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur - ...
                        MotifTime_med)<0.01*MotifTime_med; % good enough, less than 1% of motif median.
                    
                    songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
                    fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
                    %                 motifdur=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur;
                    break
                end
            end
        end
    end
    if isempty(songdat)
        % then pick a random song
        i=randi(NumNeurons,1);
        j=randi(length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract),1);
        songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
        fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
    end
    
    % --- plot
    lt_plot_spectrogram(songdat, fs, 1, 0);
    line([motif_predur motif_predur], ylim, 'Color','w');
    if LinScaleGlobal==1
        line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','w'); % end
    end
    
    % 2) ===== plot rasters
    hsplot=lt_subplot(7, 1, 2:7); hold on
    hsplots=[hsplots hsplot];
    hsplots_rasty=[hsplots_rasty hsplot];
    for i=1:NumNeurons
        if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
        
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        params=MOTIFSTATS.neurons(i).motif(m).Params;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        numtrials=length(segextract);
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            
            if ~isempty(spktimes)
                for jj=1:length(spktimes)
                    spk=spktimes(jj);
                    line([spk spk], -[yind-0.4 yind+0.4], 'Color', plotcols{i}, 'LineWidth', 1);
                end
            end
            
            % -- increment y ind
            yind=yind+1;
        end
    end
    % -- line for motif onset
    line([motif_predur motif_predur], ylim, 'Color','k');
    if LinScaleGlobal==1
        line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','k'); % end
    end
    
    % --- link
    linkaxes(hsplots, 'x');
    
end

linkaxes(hsplots_rasty, 'xy');


%% === PLOT MEAN FIRING RATE SMOOTHED OVER TIME [EACH MOTIF SEPARATE PLOT]

for m=1:NumMotifs
    lt_figure; hold on;
    %     figcount=1;
    %     subplotrows=NumNeurons+1;
    %     subplotcols=1;
    %     fignums_alreadyused=[];
    %     hfigs=[];
    
    hsplots=[];
    
    % 1) === plot spectrogram at top
    hsplot=lt_subplot(NumNeurons+1, 1, 1); hold on
    title(motif_regexpr_str{m});
    hsplots=[hsplots hsplot];
    songdat=[];
    fs=[];
    if LinScaleGlobal==1
        ylabel('median dur trial');
        % --- find the trial that has close to median motif time
        for i=1:NumNeurons
            if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            
            numtrials=length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract);
            
            for j=1:numtrials
                if abs(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur - ...
                        MotifTime_med)<0.01*MotifTime_med; % good enough, less than 1% of motif median.
                    
                    songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
                    fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
                    %                 motifdur=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur;
                    break
                end
            end
        end
    end
    if isempty(songdat)
        % then pick a random song
        i=randi(NumNeurons,1);
        j=randi(length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract),1);
        songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
        fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
    end
    
    % --- plot
    lt_plot_spectrogram(songdat, fs, 1, 0);
    line([motif_predur motif_predur], ylim, 'Color','w');
    if LinScaleGlobal==1
        line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','w'); % end
    end
    
    % 2) ===== plot rasters
    for i=1:NumNeurons
        if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
        %         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplot=lt_subplot(NumNeurons+1, 1, i+1); hold on
        
        hsplots=[hsplots hsplot];
        ylabel(['neuron ' num2str(i)]);
        
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        params=MOTIFSTATS.neurons(i).motif(m).Params;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        % -- convert to smoothed rate
        if length(Yspks)>1
        [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
        
        % --- plot
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
        end
        % -- line for motif onset
        line([motif_predur motif_predur], ylim, 'Color','k');
        if LinScaleGlobal==1
            line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','k'); % end
        end
        
        ylim([0 200]);
    end
    
    % --- link
    linkaxes(hsplots, 'x');
end


%% === PLOT MEAN FIRING RATE SMOOTHED OVER TIME [OVERLAY EACH MOTIF IN ONE PLOT

PlotForIllustrator=0;
RandomSong=1;

if PlotForIllustrator==1
    figcount=1;
    subplotrows=2;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    
    
else
    figcount=1;
    subplotrows=8;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
end
hsplots=[];

plotcols_motif=lt_make_plot_colors(NumMotifs, 0, 0);

% ==== Plot spectrograms of each motif
for m=1:NumMotifs
    %     figcount=1;
    %     subplotrows=NumNeurons+1;
    %     subplotcols=1;
    %     fignums_alreadyused=[];
    %     hfigs=[];
    
    
    % 1) === plot spectrogram at top
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title(motif_regexpr_str{m});
    hsplots=[hsplots hsplot];
    songdat=[];
    fs=[];
    if LinScaleGlobal==1
        % --- find the trial that has close to median motif time
        for i=1:NumNeurons
            if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            
            numtrials=length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract);
            
            for j=1:numtrials
                if abs(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur - ...
                        MotifTime_med)<0.01*MotifTime_med; % good enough, less than 1% of motif median.
                    
                    songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
                    fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
                    %                 motifdur=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur;
                    break
                end
            end
        end
    end
    if RandomSong==1
        % make it empty so will pick random song
        songdat=[];
    end
    if isempty(songdat)
        % then pick a random song
        i=randi(NumNeurons,1);
        j=randi(length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract),1);
        songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
        fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
    end
    
    % --- plot
    lt_plot_spectrogram(songdat, fs, 1, 0);
    line([motif_predur motif_predur], ylim, 'Color','w');
    if LinScaleGlobal==1
        line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','w'); % end
    end
    if RandomSong==1
        title(['rand song: ' motif_regexpr_str{m}], 'Color', plotcols_motif{m});
    else
        title(['dur close to median: ' motif_regexpr_str{m}], 'Color', plotcols_motif{m});
    end
end


% % ===== PLOT NEURAL DATA [SHOULD BE ENTIRELY REPLICATED BELOW, BUT
% WORKING WITH LT_PLOT_MULT
% for m=1:NumMotifs
%     %     figcount=1;
%     %     subplotrows=NumNeurons+1;
%     %     subplotcols=1;
%     %     fignums_alreadyused=[];
%     %     hfigs=[];
%
%     for i=1:NumNeurons
%         if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
%             % then there is no data for this neuron
%             continue
%         end
%         if PlotForIllustrator==1
%             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         else
%         %         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         hsplot=lt_subplot(NumNeurons+NumMotifs, 1, NumMotifs+i); hold on
%         end
%         hsplots=[hsplots hsplot];
%         ylabel(['neuron ' num2str(i)]);
%
%         segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
%         params=MOTIFSTATS.neurons(i).motif(m).Params;
%         clustnum=NeuronDatabase.neurons(i).clustnum;
%
%         numtrials=length(segextract);
%         Yspks={};
%         for j=1:numtrials
%             inds=segextract(j).spk_Clust==clustnum;
%             spktimes=segextract(j).(spktimefield)(inds);
%             Yspks{j}=spktimes;
%         end
%
%         % -- convert to smoothed rate
%         [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
%
%         % --- plot
%         shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols_motif{m}}, 1);
%
%         % -- line for motif onset
%         line([motif_predur motif_predur], ylim, 'Color','k');
%         line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','k'); % end
%     end
%
% end

% ===== PLOT NEURAL DATA
%     figcount=1;
%     subplotrows=NumNeurons+1;
%     subplotcols=1;
%     fignums_alreadyused=[];
%     hfigs=[];

for i=1:NumNeurons
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for m=1:NumMotifs
        if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
        hsplots=[hsplots hsplot];
        ylabel(['neuron ' num2str(i)]);
        
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        params=MOTIFSTATS.neurons(i).motif(m).Params;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        % -- convert to smoothed rate
        if length(Yspks)>1
        [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
        
        % --- plot
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols_motif{m}}, 1);
        end
        
        % -- line for motif onset
        line([motif_predur motif_predur], ylim, 'Color','k');
        if LinScaleGlobal==1
            line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','k'); % end
        end
        ylim([0 max(ymean_hz)+20]);
    end
    
end
% --- link
linkaxes(hsplots, 'x');

