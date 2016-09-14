function lt_neural_MultNeur_MotifRasters(NeuronDatabase, motif_regexpr_str, motif_predur, motif_postdur)
%% given motif, plots rasters and smooth rate for all neurons

%  NOTE: default is to linearly scale across all neurons (motif onset to
%  offset)

%% PARAMS

NumNeurons=length(NeuronDatabase.neurons);
LinScaleGlobal=1;% scales all trials across all neurons to global median [DEFAULT, AND ALTERNATIVES NOT YET CODED]

% smoothing window (neural)
window=0.02; windshift=0.004;

%% EXTRACT DATA
MOTIFSTATS=struct;
for i=1:NumNeurons
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    % --- EXTRACT DAT
    regexpr_str=motif_regexpr_str;
    predur=motif_predur; % sec
    postdur=motif_postdur; % sec
    alignByOnset=1;
    WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
    % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
        regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
    
    % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
    MOTIFSTATS.neurons(i).SegmentsExtract=SegmentsExtract;
    MOTIFSTATS.neurons(i).Params=Params;
    
    % ------ STATISTICS
    %     motifwind=[-0.04 -0.04]; % e.g. 40ms premotif onset to 40ms pre motif offest
    %     Window_relOnset={};
    %     Window_relOffset={};
    % %     Window_relOnset{1}=[-10 -9.5]; % window rel 1st syl onset
    % %     Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
    % %     Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
    % %     Window_relOffset{2}=[9.5 10];
    %     Window_relOnset{1}=[-6 -5.5]; % window rel 1st syl onset
    %     Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
    %     Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
    %     Window_relOffset{2}=[5.5 6];
    %     ISIthreshburst=0.005; % percent of spikes that are within birsts.
    %     clustwanted=NeuronDatabase.neurons(i).clustnum; % wave_clus cluster
    %
    %     [FiringRateOut BurstFracOut] = lt_neural_CalMeanFiring(SegmentsExtract, ...
    %         Params, clustwanted, Window_relOnset, motifwind, Window_relOffset, ISIthreshburst);
    
    %     MOTIFSTATS.neurons(i).FiringRateOut=FiringRateOut;
    %     MOTIFSTATS.neurons(i).BurstFracOut=BurstFracOut;
    
end


%% ========================== PLOT COMBINED RASTERS
if LinScaleGlobal==1
    % =================== OPTIONAL, LINEARLY TIME WARP [DEFAULT]
    spktimefield='spk_Times_scaled';
    % -- 1) figure out median dur across all neurons (for this motif)
%     OnAll=[];
%     OffAll=[];
    MotifTimesAll=[];
    
    for i=1:NumNeurons
        if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
        
%         predur=MOTIFSTATS.neurons(i).Params.REGEXP.predur;
%         postdur=MOTIFSTATS.neurons(i).Params.REGEXP.postdur;
%         
%         motifontimes = predur + ...
%             [MOTIFSTATS.neurons(i).SegmentsExtract.global_ontime_motifInclFlank];
%         motifofftimes=[MOTIFSTATS.neurons(i).SegmentsExtract.global_offtime_motifInclFlank] - ...
%             postdur;
%         
%         OnAll=[OnAll motifontimes];
%         OffAll=[OffAll motifofftimes];
        
                MotifTimesAll=[MotifTimesAll MOTIFSTATS.neurons(i).SegmentsExtract.actualmotifdur];
    end
%     MotifTimesAll=OffAll-OnAll;
    lt_figure; hold on;
    lt_plot_histogram(MotifTimesAll); xlabel('all motif times (s)');
    MotifTime_med=median(MotifTimesAll);
    
    % -- 2) time warp
    for i=1:NumNeurons
        if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
        segmentdat=MOTIFSTATS.neurons(i).SegmentsExtract;
        prms=MOTIFSTATS.neurons(i).Params;
        [segmentdatout, prmsout] = lt_neural_LinTimeWarp(segmentdat, ...
            prms, MotifTime_med);
        
        MOTIFSTATS.neurons(i).SegmentsExtract=segmentdatout;
        MOTIFSTATS.neurons(i).Params=prmsout;
    end
end

%% =================== PLOT ONE GLOBAL RASTER
lt_figure; hold on;
ylabel('trial #');
xlabel('sec');
yind=1; % will increment to plot over all experiments
plotcols=lt_make_plot_colors(NumNeurons, 0, 0);
hsplots=[];

% 1) === plot spectrogram at top
hsplot=lt_subplot(8, 1, 1); hold on
ylabel('median dur trial');
hsplots=[hsplots hsplot];
% --- find the trial that has the median motif time
songdat=[];
fs=[];
tmptmp=[];
for i=1:NumNeurons
    if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
        % then there is no data for this neuron
        continue
    end
    
    numtrials=length(MOTIFSTATS.neurons(i).SegmentsExtract);
    
    for j=1:numtrials
        if abs(MOTIFSTATS.neurons(i).SegmentsExtract(j).actualmotifdur - MotifTime_med) < ...
                0.01*MotifTime_med;
            % found the median trial!!!
            disp('sdafasdf')
            songdat=MOTIFSTATS.neurons(i).SegmentsExtract(j).songdat;
            fs=MOTIFSTATS.neurons(i).SegmentsExtract(j).fs;
            motifdur=MOTIFSTATS.neurons(i).SegmentsExtract(j).actualmotifdur;
            tmptmp=[tmptmp motifdur]
            break
        end
    end
end

% --- plot
lt_plot_spectrogram(songdat, fs, 1, 0);
line([motif_predur motif_predur], ylim, 'Color','w');
line([motif_predur+motifdur motif_predur+motifdur], ylim, 'Color','w'); % end
title(motif_regexpr_str);

% 2) ===== plot rasters
hsplot=lt_subplot(8, 1, 2:8); hold on
hsplots=[hsplots hsplot];
for i=1:NumNeurons
    if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
        % then there is no data for this neuron
        continue
    end
    
    segextract=MOTIFSTATS.neurons(i).SegmentsExtract;
    params=MOTIFSTATS.neurons(i).Params;
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
line([motif_predur+motifdur motif_predur+motifdur], ylim, 'Color','k'); % end

% --- link
linkaxes(hsplots, 'x');



%% === PLOT MEAN FIRING RATE SMOOTHED OVER TIME

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

hsplots=[];

% 1) === plot spectrogram at top
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots=[hsplots hsplot];
% --- find the trial that has the median motif time
songdat=[];
fs=[];
for i=1:NumNeurons
    if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
        % then there is no data for this neuron
        continue
    end
    
    numtrials=length(MOTIFSTATS.neurons(i).SegmentsExtract);
    
    for j=1:numtrials
        if abs(MOTIFSTATS.neurons(i).SegmentsExtract(j).actualmotifdur - MotifTime_med) < ...
                0.01*MotifTime_med;
            % found the median trial!!!
            
            songdat=MOTIFSTATS.neurons(i).SegmentsExtract(j).songdat;
            fs=MOTIFSTATS.neurons(i).SegmentsExtract(j).fs;
            motifdur=MOTIFSTATS.neurons(i).SegmentsExtract(j).actualmotifdur;
            break
        end
    end
end
% --- plot
lt_plot_spectrogram(songdat, fs, 1, 0);
line([motif_predur motif_predur], ylim, 'Color','w');
line([motif_predur+motifdur motif_predur+motifdur], ylim, 'Color','w'); % end
title(motif_regexpr_str);

% 2) ===== plot rasters
for i=1:NumNeurons
    if ~isfield(MOTIFSTATS.neurons(i).SegmentsExtract, 'spk_Times')
        % then there is no data for this neuron
        continue
    end
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    ylabel(['neuron ' num2str(i)]);
    
    segextract=MOTIFSTATS.neurons(i).SegmentsExtract;
    params=MOTIFSTATS.neurons(i).Params;
    clustnum=NeuronDatabase.neurons(i).clustnum;
    
    numtrials=length(segextract);
    Yspks={};
    for j=1:numtrials
        inds=segextract(j).spk_Clust==clustnum;
        spktimes=segextract(j).(spktimefield)(inds);
        Yspks{j}=spktimes;
    end
    
    % -- convert to smoothed rate
    [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
    
    % --- plot
    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
    
    % -- line for motif onset
    line([motif_predur motif_predur], ylim, 'Color','k');
    line([motif_predur+motifdur motif_predur+motifdur], ylim, 'Color','k'); % end
end

% --- link
linkaxes(hsplots, 'x');


