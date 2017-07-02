function lt_neural_MultNeur_SongMod(NeuronDatabase, regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, Window_relOnset, motifwind, Window_relOffset)
%% == analysis for song modulation (e.g. onset, offset, FR, ISI, etc)

% === PARAM
%     regexpr_str='WHOLEBOUTS';
%     predur=6; % sec
%     postdur=6; % sec
%     alignByOnset=1;
%     WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps

%     motifwind=[-0.04 -0.04]; % e.g. 40ms premotif onset to 40ms pre motif offest
%     Window_relOnset{1}=[-6 -5.5]; % window rel 1st syl onset
%     Window_relOnset{2}=[-0.6 -0.1]; % window rel 1st syl onset
%     Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
%     Window_relOffset{2}=[5.5 6];
    
    %% auto params
    
    ISIthreshburst=0.005; % percent of spikes that are within birsts.

    
   %% === make a plot depicting timepoint of windows
   lt_figure; hold on;
   xlabel('sec');
   song_length=10; % arbitrary, 10s
   totaldur=song_length+predur+postdur;
   
   line([0 totaldur], [0.5 0.5], 'Color','k','LineWidth', 3);
   
   
   counter=1;
   % - prewind
   for i=1:length(Window_relOnset)
       wind=Window_relOnset{i};
       
       line([predur+wind(1) predur+wind(2)], [0.5 0.5], 'Color','m', 'LineWidth', 7);
       lt_plot_text(predur+wind(1), 0.7, num2str(counter), 'm');
       counter=counter+1;
   end
   
   % - song
   line([predur predur+song_length], [0.5 0.5],'Color','b','LineWidth', 10);
   lt_plot_text(predur+1,0.7, [num2str(counter) ' (song)'], 'b');
counter=counter+1;
   % - postwin
   for i=1:length(Window_relOffset)
       wind=Window_relOffset{i};
       
       line([predur+song_length+wind(1) predur+song_length+wind(2)], [0.5 0.5], 'Color','m', 'LineWidth', 7);
        lt_plot_text(predur+song_length+wind(1), 0.7, num2str(counter), 'm');
       counter=counter+1;
  end
   
       
   
%% === 1) Bout onset and offset modulation for each neuron

NumNeurons=length(NeuronDatabase.neurons);

% %% ======================= FOR EACH NEURON, PLOT IT'S OWN RASTERS, AND SM FRATE
% 
% % ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
% for i=1:NumNeurons
%     cd(NeuronDatabase.global.basedir);
%     
%     % - find day folder
%     dirdate=NeuronDatabase.neurons(i).date;
%     tmp=dir([dirdate '*']);
%     assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
%     cd(tmp(1).name);
%     
%     % - load data for this neuron
%     batchf=NeuronDatabase.neurons(i).batchfile;
%     channel_board=NeuronDatabase.neurons(i).chan;
%     [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
%     
%     % --- EXTRACT DAT
%     regexpr_str='WHOLEBOUTS';
%     predur=6; % sec
%     postdur=6; % sec
%     alignByOnset=1;
%     WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
%     % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
%     [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
%         regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
%     
%     
%         % ==================== Plot individually for this neuron
%         useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
%         plotAllSegs=0; % then plots each trial own plot.
%         [Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);
%     
% end


%% COLLECT AND PLOT STATS ACROSS NEURONS

% ================ MOTIF STATISTICS (E.G. FIRING RATE, BURSTS, ...)
MOTIFSTATS=struct;
for i=1:NumNeurons
    try 
    cd(NeuronDatabase.global.basedir);
    catch err
        cd(NeuronDatabase.neurons(i).basedir);
    end
    
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
    clustwanted=NeuronDatabase.neurons(i).clustnum; % wave_clus cluster

    [FiringRateOut BurstFracOut] = lt_neural_CalMeanFiring(SegmentsExtract, ...
        Params, clustwanted, Window_relOnset, motifwind, Window_relOffset, ISIthreshburst);
    
    MOTIFSTATS.neurons(i).FiringRateOut=FiringRateOut;
    MOTIFSTATS.neurons(i).BurstFracOut=BurstFracOut;
    
end


%% =========================== PLOT SUMMARY STATS ACROSS NEURONS (E.G.
% FIRING RATE
if NumNeurons>1
    
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

lt_plot_text(Xmean(1), 1.2*Ymean(1), ['n=' num2str(size(Yall, 1))], 'b')

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

end


