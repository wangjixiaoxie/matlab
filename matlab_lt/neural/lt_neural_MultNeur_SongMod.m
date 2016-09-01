%% == analysis for song modulation (e.g. onset, offset, FR, ISI, etc)

%% === 1) Bout onset and offset modulation for each neuron
currdir=pwd;
close all;

%% ======================= FOR EACH NEURON, PLOT IT'S OWN RASTERS, AND SM FRATE

% ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
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
    regexpr_str='WHOLEBOUTS';
    predur=6; % sec
    postdur=6; % sec
    alignByOnset=1;
    WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
    % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
        regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
    
    
        % ==================== Plot individually for this neuron
        useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
        plotAllSegs=0; % then plots each trial own plot.
        [Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);
    
end



%% COLLECT AND PLOT STATS ACROSS NEURONS

% ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
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
    regexpr_str='WHOLEBOUTS';
    predur=6; % sec
    postdur=6; % sec
    alignByOnset=1;
    WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
    % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
        regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
    
    
    % ------ STATISTICS
    motifwind=[-0.04 -0.04]; % e.g. 40ms premotif onset to 40ms pre motif offest
    Window_relOnset={};
    Window_relOffset={};
%     Window_relOnset{1}=[-10 -9.5]; % window rel 1st syl onset
%     Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
%     Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
%     Window_relOffset{2}=[9.5 10];
    Window_relOnset{1}=[-6 -5.5]; % window rel 1st syl onset
    Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
    Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
    Window_relOffset{2}=[5.5 6];
    ISIthreshburst=0.005; % percent of spikes that are within birsts.
    clustwanted=NeuronDatabase.neurons(i).clustnum; % wave_clus cluster

    [FiringRateOut BurstFracOut] = lt_neural_CalMeanFiring(SegmentsExtract, ...
        Params, clustwanted, Window_relOnset, motifwind, Window_relOffset, ISIthreshburst);
    
    MOTIFSTATS.neurons(i).FiringRateOut=FiringRateOut;
    MOTIFSTATS.neurons(i).BurstFracOut=BurstFracOut;
    
end


%% =========================== PLOT SUMMARY STATS ACROSS NEURONS (E.G.
% FIRING RATE
lt_figure; hold on;
hsplots=[];
% ----- 1) using actual mean timepoints
hsplot=lt_subplot(3,2,1); hold on;
hsplots=[hsplots hsplot];
title('firing rate modulation');
xlabel('mean time of window (sec)');
ylabel('sp/sec');
Xall=[];
Yall=[];
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    
    % --- firing rate
    x=MOTIFSTATS.neurons(i).FiringRateOut.xmean_sec;
    y=MOTIFSTATS.neurons(i).FiringRateOut.Yall_rates;
    ymean=mean(y);
    ystd=std(y);
    
    lt_plot(x, ymean, {'Errors', ystd, 'LineStyle', '-', 'Color', plotcol});
        
    Xall=[Xall; x];
    Yall=[Yall; ymean];
end

% -- plot means
Xmean=mean(Xall);
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(Xmean, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});

% - line for song onset
predur=Params.REGEXP.predur;
line([predur predur], ylim, 'Color', 'b');

% ------- 3) ISI using actual mean timepoints
lt_subplot(3,2,2); hold on;
title('% sp in burst');
xlabel('mean time of window (sec)');
ylabel('%');
Xall=[];
Yall=[];
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    % --- firing rate
    x=MOTIFSTATS.neurons(i).FiringRateOut.xmean_sec;
    y=MOTIFSTATS.neurons(i).BurstFracOut.Yall_FracSpInbursts;

    lt_plot(x, y, {'LineStyle', '-', 'Color', plotcol});
        
    Xall=[Xall; x];
    Yall=[Yall; y];
end

% -- plot means
Xmean=mean(Xall);
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(Xmean, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});

% - line for song onset
predur=Params.REGEXP.predur;
line([predur predur], ylim, 'Color', 'b');
ylim([0 1]);



% ----- 2) pre (nonsinging), dur song, immediately post)
xinds=[1 3 4]; % IMPORTANT - these are the segments to plot [nonsinging, song, post];
hsplot=lt_subplot(3,2,3); hold on;
hsplots=[hsplots hsplot];
title('firing rate modulation');
xlabel('winndow number');
ylabel('sp/sec');
Yall=[];
X=1:length(xinds);
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    % --- firing rate
    y=MOTIFSTATS.neurons(i).FiringRateOut.Yall_rates;
    ymean=mean(y);
    ystd=std(y);
    
    ymean=ymean(xinds);
    ystd=ystd(xinds);
    lt_plot(X, ymean, {'Errors', ystd, 'LineStyle', '-', 'Color', plotcol});
        
    Yall=[Yall; ymean];
end

% -- plot means
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(X, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});

linkaxes(hsplots, 'y');


% ----- 4) pre (nonsinging), dur song, immediately post) [ISI]
xinds=[1 3 4]; % IMPORTANT - these are the segments to plot [nonsinging, song, post];
lt_subplot(3,2,4); hold on;
title('% sp in burst');
xlabel('window num');
ylabel('%');
Yall=[];
X=1:length(xinds);
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    % --- firing rate
    y=MOTIFSTATS.neurons(i).BurstFracOut.Yall_FracSpInbursts;

    y=y(xinds);
    lt_plot(X, y, {'LineStyle', '-', 'Color', plotcol});
        
    Yall=[Yall; y];
end

% -- plot means
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(X, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});
ylim([0 1]);




% ----- 2) pre (nonsinging), dur song, immediately post)
hsplot=lt_subplot(3,2,5); hold on;
hsplots=[hsplots hsplot];
title('firing rate modulation');
xlabel('winndow number');
ylabel('sp/sec');
Yall=[];
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    % --- firing rate
    y=MOTIFSTATS.neurons(i).FiringRateOut.Yall_rates;
    ymean=mean(y);
    ystd=std(y);
X=1:length(ymean);
    
    lt_plot(X, ymean, {'Errors', ystd, 'LineStyle', '-', 'Color', plotcol});
        
    Yall=[Yall; ymean];
end

% -- plot means
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(X, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});

linkaxes(hsplots, 'y');


% ----- 4) pre (nonsinging), dur song, immediately post) [ISI]
lt_subplot(3,2,6); hold on;
title('% sp in burst');
xlabel('window num');
ylabel('%');
Yall=[];
for i=1:NumNeurons
    
    % -- plot col (heavier if is mu);
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        plotcol=[0.6 0.6 0.6]; % mu
    else
        plotcol=[0.8 0.8 0.8]; % su
    end
    % --- firing rate
    y=MOTIFSTATS.neurons(i).BurstFracOut.Yall_FracSpInbursts;
X=1:length(ymean);

    lt_plot(X, y, {'LineStyle', '-', 'Color', plotcol});
        
    Yall=[Yall; y];
end

% -- plot means
Ymean=mean(Yall);
Ysem=lt_sem(Yall);
lt_plot(X, Ymean, {'Errors', Ysem, 'LineStyle', '-', 'Color', 'r', 'Marker', 's'});
ylim([0 1]);



