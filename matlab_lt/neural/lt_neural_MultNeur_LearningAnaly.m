function lt_neural_MultNeur_LearningAnaly(NeuronDatabase, motif_regexpr_str, motif_predur, ...
    motif_postdur, LinScaleGlobal, FFparams, premotor_wind)
%% LT 2/7/17 - modified from lt_neural_MultNeur_MotifRasters_LearnSum
% here, looks trial by trial to analyze learning. takes FF in a premotor
% window for all syllables and analyzes in batch.


%% for each motif and each neuron, make one plot, showing progression over time of
% motif_regexpr_str={'g(h)'};
% motif_predur=0.2;
% motif_postdur=0.1;
% LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
%
%
% FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
% FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
%     % +1 is 1 after token
% FFparams.FF_sylName='h'; % Optional: what syl do you expect this to be? if incompatible will raise error
%     % not required (can set as []);
% FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
%             'v', [2450 4300]};
% % FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
% %             'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
% FFparams.cell_of_FFtimebins={'h', [0.034 0.038], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % WN on g H
% % NOTE: will also determine whether was hit or miss, based on WN sound
% % detection.
%
% % LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'};
%
% OnlyPlotNoHit=1; % then only plots trials that were not hit (WN)
% TypesToPlot={'mean', 'abs', 'std'}; % for d-prime summary metric


%% PARAMS
NumMotifs=length(motif_regexpr_str);
NumNeurons=length(NeuronDatabase.neurons);

% smoothing window (neural)
window=0.02; windshift=0.004;

plotTimesOnRaster=1; % hh:mm on raster plot

determineTrialBinBasedOnBaselineN=1; % if 1, then chooses binsize based on
% baseline num songs (i.e. baseNumSongs/DivisorBaseSongs); if 0, then = TrialBinSize
% note: assumes WNchangeDateStrings{1} is transition from baseline to WN
% on.
DivisorBaseSongs=2;
% TrialBinSize=10;

if LinScaleGlobal ==0
    spktimefield='spk_Times';
end

% for taking d-prime, relative to token onset
% premotor_wind=[-0.08 0.02]; % in sec, relative to onset of token in motif.

% UseEntireBaseline = 0; % if 1, uses entire baseline, otherwise uses 1st bin.



%% EXTRACT DATA (for each neuron x motif)
MOTIFSTATS=struct;
for i=1:NumNeurons
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
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
    % - do this one time for each desired motif
    for j=1:NumMotifs
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams);
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
    end
end

% pause;
% disp('PRESS ANYTHING TO CONTINUE');
% close all;


%% === MAKE SURE ALL TRIALS ARE IN TEMPORAL ORDER

for i=1:NumNeurons
    for m=1:NumMotifs
        
        segextract = MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        
        all_datenums=[segextract.song_datenum];
        
        [~, inds] = sort(all_datenums);
        
        if all(diff(inds))~=1
            disp('HAD TO REORDER !! ---- not a problem');
            
            MOTIFSTATS.neurons(i).motif(m).SegmentsExtract=segextract(inds);
        end
    end
end



%% CORRELATION BETWEEN CHANGE MEAN FR (FROM BASELINE) AND CHANGE IN FF
% NOTE; DIFFERENT FROM ABOVE, IN THAT HERE 1) DOES TRIAL BY TRIAL AND 2)
% DOES NOT AVERAGE OVER A MULTIPLE PREMOTOR BINS, BUT INSTEAD TAKES ONE BIN

% NOTE: ALWAYS PLOTS BOTH HITS AND MISSES, EVEN IF OnlyPlotNoHit =1
% NOTE: ALWAYS USES ENTIRE BASELINE

for i=1:NumNeurons
    
    clear WNchangeDateStrings;
    WNchangeDateStrings{1}=NeuronDatabase.neurons(i).LEARN_WNonDatestr;
    WNchangeDateStrings=[WNchangeDateStrings NeuronDatabase.neurons(i).LEARN_WNotherImportantDates];
    
    
    for m=1:NumMotifs
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        
        % ===================== collect spikes for all trials
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        
        % ===================== BASELINE
        % 1) mean and std spiking
        datestring=WNchangeDateStrings{1}; % assumes this is WN on.
        dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
        inds=[segextract.song_datenum] < dnum; %
        
        windowstart = motif_predur + premotor_wind(1);
        windowend = motif_predur + premotor_wind(2); % time of end of window (rel to data onset)
        TMPSTRUCT = lt_neural_GetStatsSingleWindow(Yspks(inds), windowstart, windowend);
        ymean_hz_base=TMPSTRUCT.FrateMean;
        ystd_hz_base=TMPSTRUCT.FrateStd;
        
        % 2) FF
        mean_FF_base=mean([segextract(inds).FF_val]);
        std_FF_base = std([segextract(inds).FF_val]);
        
        
        % ======== FOR EACH TRIAL (INCLUDING BASELINE)
        numtrials = length(segextract);
        assert(length(Yspks) == numtrials, 'asdfasdcas');
        
        Frate_All = [];
        FF_All = [];
        Is_Baseline = [];
        TimeDayvals_All = [];
        
        
        for j=1:numtrials
            % 1) FRate
            tmpstruct = lt_neural_GetStatsSingleWindow(Yspks(j), windowstart, windowend);
            Frate_All = [Frate_All tmpstruct.FrateMean];
            
            % 2) FF
            FF_All = [FF_All segextract(j).FF_val];
            
            % 3) time
            firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
            eventtime=datestr(segextract(j).song_datenum, 'ddmmmyyyy-HHMM');
            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime}); % days from start of expt
            TimeDayvals_All = [TimeDayvals_All tmp.FinalValue];
            
            % 4) is baseline?
            datestring=WNchangeDateStrings{1}; % assumes this is WN on.
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            Is_Baseline = [Is_Baseline double(segextract(j).song_datenum < dnum)];
        end
        
        %% 5) FF, deviation from local mean
        if (1) % REPLACES ABOVE VARIABLES
            localbinsize = 15;
            
            tmp = lt_running_stats(FF_All, localbinsize);
            tmp2 = FF_All(ceil(localbinsize/2):end-floor(localbinsize/2)); % get middle FF for each local bin
            FF_All_localdev = tmp2 - tmp.Mean;
            
            %         tmp = lt_running_stats(Frate_All, localbinsize);
            %         tmp2 = Frate_All(ceil(localbinsize/2):end-floor(localbinsize/2)); % get middle FF for each local bin
            %         Frate_All = tmp2 - tmp.Mean;
            
            TimeDayvals_All_localdev = TimeDayvals_All(ceil(localbinsize/2):end-floor(localbinsize/2)); % get middle FF for each local bin
            Frate_All_localdev = Frate_All(ceil(localbinsize/2):end-floor(localbinsize/2)); % get middle FF for each local bin
            Is_Baseline_localdev = Is_Baseline(ceil(localbinsize/2):end-floor(localbinsize/2)); % get middle FF for each local bin
            
        end
        
        
        %% ====================== PLOT
        lt_figure; hold on;
        
        hsplots = [];
        hsplots0 = [];
        
        % 1) FF
        hsplot = lt_subplot(6,1,1); hold on;
        hsplots0 = [hsplots0 hsplot];
        title(['neuron ' num2str(i) ', ' motif_regexpr_str{m}]);
        ylabel('FF');
        xlabel('time (days from expt start');
        plot(TimeDayvals_All, FF_All, 'or');
        inds = Is_Baseline==1;
        plot(TimeDayvals_All(inds), FF_All(inds), 'ok');
        shadedErrorBar(xlim, [mean_FF_base mean_FF_base], [std_FF_base std_FF_base], {'Color','k'}, 1)
    
        % ---- lines
    for j=1:length(WNchangeDateStrings)
    datestring=WNchangeDateStrings{j}; %
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    linetime = tmp.FinalValue;
    line([linetime linetime], ylim);        
    end
        
        
        % 2) hz
        hsplot = lt_subplot(6,1,2); hold on;
        hsplots0 = [hsplots0 hsplot];
        title(['neuron ' num2str(i) ', ' motif_regexpr_str{m}]);
        ylabel('Frate');
        xlabel('time (days from expt start');
        plot(TimeDayvals_All, Frate_All, 'o');
        inds = Is_Baseline==1;
        plot(TimeDayvals_All(inds), Frate_All(inds), 'ok');
        shadedErrorBar(xlim, [ymean_hz_base ymean_hz_base], [ystd_hz_base ystd_hz_base], {'Color','k'}, 1)
        
        linkaxes(hsplots0, 'x');
        
        % 3a) FF corr with Frate, across all data
        numcols = 3;
        counter = 0;
        
        hsplot = lt_subplot(6,numcols, counter+1+numcols*2); hold on;
        hsplots = [hsplots hsplot];
        title('all data');
        xlabel('Frate_All')
        ylabel('FF');
        lt_regress(FF_All, Frate_All, 1, 0, 1, 1, 'k')
        counter = counter+1;
        
        % 3) baseline
        hsplot = lt_subplot(6,numcols, counter+1+numcols*2); hold on;
        hsplots = [hsplots hsplot];
        title('baseline');
        inds = Is_Baseline==1;
        xlabel('Frate_All')
        ylabel('FF');
        lt_regress(FF_All(inds), Frate_All(inds), 1, 0, 1, 1, 'k')
        counter = counter+1;
        
        % --- one plot for each day
        numdays = max(floor(TimeDayvals_All));
        for k=1:numdays
            hsplot = lt_subplot(6,numcols, counter+1+numcols*2+k-1); hold on;
            hsplots = [hsplots hsplot];
            title(['[WN] day ' num2str(k)]);
            xlabel('Frate_All')
            ylabel('FF');
            
            inds = Is_Baseline == 0 & floor(TimeDayvals_All) == k;
            if any(inds)
                lt_regress(FF_All(inds), Frate_All(inds), 1, 0, 1, 1, 'r')
            end
        end
        numsubplots = counter+1+numcols*2+k-1;
        % --
        linkaxes(hsplots, 'xy');
        
        % ============ plot using local deviation of FF
        if numsubplots+numdays+1 > 6*numcols
            lt_figure; hold on;
            numsubplots = 0;
        end  
        % 1) --- baseline
        lt_subplot(6,numcols, numsubplots+1); hold on;
        title(['baseline [local FF dev, bin' num2str(localbinsize) ']']);
        inds = Is_Baseline_localdev==1;
        xlabel('Frate_All')
        ylabel('FF');
        lt_regress(FF_All_localdev(inds), Frate_All_localdev(inds), 1, 0, 1, 1, 'k')
        
        numsubplots = numsubplots + 1;
        
        % 2) WN days
        numdays = max(floor(TimeDayvals_All_localdev));
        for k=1:numdays
            lt_subplot(6,numcols,numsubplots+k); hold on;
            title(['[WN] day ' num2str(k) '[local FF dev]']);
            xlabel('Frate_All')
            ylabel('FF');
            
            inds = Is_Baseline_localdev == 0 & floor(TimeDayvals_All_localdev) == k;
            if any(inds)
                lt_regress(FF_All_localdev(inds), Frate_All_localdev(inds), 1, 0, 1, 1, 'r')
            end
        end
        
        
    end
end

%% DOES DIRECTION OF FF-NEURAL CORRELATION IN A CERTAIN WINDOW ALLOW US TO USE FRATE CHANGE
% TO PREDICT FF CHANGE?

for i=1:NumNeurons
    
    clear WNchangeDateStrings;
    WNchangeDateStrings{1}=NeuronDatabase.neurons(i).LEARN_WNonDatestr;
    WNchangeDateStrings=[WNchangeDateStrings NeuronDatabase.neurons(i).LEARN_WNotherImportantDates];
    
    TMPSTRUCT = struct;
    for m=1:NumMotifs
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        % ===================== collect spikes for all trials
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        % ======== FOR EACH TRIAL (INCLUDING BASELINE)
        Frate_All = [];
        FF_All = [];
        Is_Baseline = [];
        TimeDayvals_All = [];
        
        tmp=lt_neural_GetStatsSingleWindow(Yspks, windowstart, windowend);
        Frate_All = [tmp.NumSpksRAW]./(premotor_wind(2)-premotor_wind(1));
        
        FF_All = [segextract.FF_val];
        
        datestring=WNchangeDateStrings{1}; % assumes this is WN on.
        dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
        Is_Baseline = [segextract.song_datenum] < dnum;
        
        for j=1:numtrials
            % 3) time
            firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
            eventtime=datestr(segextract(j).song_datenum, 'ddmmmyyyy-HHMM');
            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime}); % days from start of expt
            TimeDayvals_All = [TimeDayvals_All tmp.FinalValue];
        end
        
        % === output
        TMPSTRUCT(m).motifregexp = motif_regexpr_str{m};
        TMPSTRUCT(m).Is_Baseline = Is_Baseline;
        TMPSTRUCT(m).TimeDayvals_All = TimeDayvals_All;
        TMPSTRUCT(m).FF_All = FF_All;
        TMPSTRUCT(m).Frate_All = Frate_All;
    end
    
    
    
    %% ============================= PLOTS
    plotcols = lt_make_plot_colors(NumMotifs, 0,0);
    lt_figure; hold on;
    numrows = 6;
    numcols = 1;
    trialbinsize = 10;
    hsplots = [];
    
    % ========= 1) FF (dev from base)
    hsplot = lt_subplot(numrows, numcols, 1); hold on;
    hsplots = [hsplots hsplot];
    xlabel('days');
    ylabel('FF');
    
    for m=1:NumMotifs
        x = TMPSTRUCT(m).TimeDayvals_All;
        y_basemean = mean(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        y_zscore = (TMPSTRUCT(m).FF_All - y_basemean)./y_basestd;
        
        %            plot(x, y_zscore, 'x', 'Color', plotcols{m});
        % -- running smooth
        tmp = lt_running_stats(y_zscore, trialbinsize);
        y_zscore = tmp.Median;
        y_zscore_sem = tmp.SEM;
        
        tmp = lt_running_stats(x, trialbinsize);
        x = tmp.Median;
        
        %            shadedErrorBar(x, y_zscore, y_zscore_sem, {'Marker', 'o', 'Color', plotcols{m}}, 1);
        shadedErrorBar(x, y_zscore, y_zscore_sem, {'Color', plotcols{m}}, 1);
        
        
    end
    
    % ---- lines
    for j=1:length(WNchangeDateStrings)
    datestring=WNchangeDateStrings{j}; %
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    linetime = tmp.FinalValue;
    line([linetime linetime], ylim);        
    end
    datestring=WNchangeDateStrings{1}; % assumes this is WN on.
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    BaselineEndTime = tmp.FinalValue;
    line([BaselineEndTime BaselineEndTime], ylim);
    
    lt_plot_zeroline;
    
    
    % ========= 1) FF (dev from base)
    hsplot = lt_subplot(numrows, numcols, 2); hold on;
    hsplots = [hsplots hsplot];
    xlabel('days');
    ylabel('FF');
    
    for m=1:NumMotifs
        x = TMPSTRUCT(m).TimeDayvals_All;
        y_basemean = mean(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        y_zscore = (TMPSTRUCT(m).FF_All - y_basemean)./y_basestd;
        
        %            plot(x, y_zscore, 'x', 'Color', plotcols{m});
        % -- running smooth
        tmp = lt_running_stats(y_zscore, trialbinsize);
        y_zscore = tmp.Median;
        y_zscore_sem = tmp.SEM;
        
        tmp = lt_running_stats(x, trialbinsize);
        x = tmp.Median;
        
        plot(x, y_zscore, '-', 'Color', plotcols{m}, 'LineWidth', 2);
        
        
    end
    
    % ---- lines
    datestring=WNchangeDateStrings{1}; % assumes this is WN on.
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    BaselineEndTime = tmp.FinalValue;
    line([BaselineEndTime BaselineEndTime], ylim);
    
    lt_plot_zeroline;
    
    
    % ============ 2) Frate
    hsplot =  lt_subplot(numrows, numcols, 3); hold on;
    hsplots = [hsplots hsplot];
    
    xlabel('days');
    ylabel('Frate (hz)');
    
    for m=1:NumMotifs
        x = TMPSTRUCT(m).TimeDayvals_All;
        y_basemean = mean(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_zscore = (TMPSTRUCT(m).Frate_All - y_basemean)./y_basestd;
        
        %            plot(x, y_zscore, 'x', 'Color', plotcols{m});
        % -- running smooth
        tmp = lt_running_stats(y_zscore, trialbinsize);
        y_zscore = tmp.Median;
        y_zscore_sem = tmp.SEM;
        
        tmp = lt_running_stats(x, trialbinsize);
        x = tmp.Median;
        
        %            shadedErrorBar(x, y_zscore, y_zscore_sem, {'Marker', 'o', 'Color', plotcols{m}}, 1);
        shadedErrorBar(x, y_zscore, y_zscore_sem, {'Color', plotcols{m}}, 1);
        
        
    end
    
    % ---- lines
    datestring=WNchangeDateStrings{1}; % assumes this is WN on.
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    BaselineEndTime = tmp.FinalValue;
    line([BaselineEndTime BaselineEndTime], ylim);
    
    lt_plot_zeroline;
    
    
    % ============ 2) Frate
    hsplot =  lt_subplot(numrows, numcols, 4); hold on;
    hsplots = [hsplots hsplot];
    
    xlabel('days');
    ylabel('Frate (hz)');
    
    for m=1:NumMotifs
        x = TMPSTRUCT(m).TimeDayvals_All;
        y_basemean = mean(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_zscore = (TMPSTRUCT(m).Frate_All - y_basemean)./y_basestd;
        
        %            plot(x, y_zscore, 'x', 'Color', plotcols{m});
        % -- running smooth
        tmp = lt_running_stats(y_zscore, trialbinsize);
        y_zscore = tmp.Median;
        y_zscore_sem = tmp.SEM;
        
        tmp = lt_running_stats(x, trialbinsize);
        x = tmp.Median;
        
        plot(x, y_zscore, '-', 'Color', plotcols{m}, 'LineWidth', 2);
        
    end
    
    % ---- lines
    datestring=WNchangeDateStrings{1}; % assumes this is WN on.
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    BaselineEndTime = tmp.FinalValue;
    line([BaselineEndTime BaselineEndTime], ylim);
    
    lt_plot_zeroline;
    
    
    
    % ============ 2) Frate
    hsplot =  lt_subplot(numrows, numcols, 5); hold on;
    hsplots = [hsplots hsplot];
    
    xlabel('days');
    ylabel('Frate (hz)');
    title('RAW');
    
    for m=1:NumMotifs
        x = TMPSTRUCT(m).TimeDayvals_All;
        y_basemean = mean(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_zscore = (TMPSTRUCT(m).Frate_All - y_basemean)./y_basestd;
        
        %            plot(x, y_zscore, 'x', 'Color', plotcols{m});
        
        plot(x, y_zscore, 'x', 'Color', plotcols{m});
        
    end
    
    % ---- lines
    datestring=WNchangeDateStrings{1}; % assumes this is WN on.
    dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
    tmp=lt_convert_EventTimes_to_RelTimes(firstday, dnum); % days from start of expt
    BaselineEndTime = tmp.FinalValue;
    line([BaselineEndTime BaselineEndTime], ylim);
    
    lt_plot_zeroline;
    
    
    % ============================ 3) BASELINE FF vs. FRATE CORR
    lt_subplot(numrows, 2, 1+5*2); hold on;
    title('baseline');
    ylabel('regression b (95%)');
    xlabel('motif')
    
    for m = 1:NumMotifs
        
        ff_base = TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1);
        frate_base = TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1);
        
        [b, bint] = lt_regress(ff_base, frate_base, 0, 0, 1, 1, 'k');
        
        plot(m, b(2), 'o','Color', plotcols{m});
        line([m m], [bint(2,:)], 'Color', plotcols{m}, 'LineWidth', 2);
        
    end
    set(gca, 'XTick', 1:length(motif_regexpr_str), 'XTickLabel', motif_regexpr_str);
    
    xlim([0 NumMotifs+1]);
    lt_plot_zeroline
    
    % ===========================
    linkaxes(hsplots, 'x');
    
    
    
    % ==========================================
    lt_subplot(numrows, 2, 2+5*2); hold on;
    xlabel('baseline beta (ff vs. frate)');
    ylabel('ff change (zscore)');
    zlabel('frate change (zscore)');
    
    WindowPeakLearning = input('Window for peak learning, e.g. [2.4 2.8] days ');
    BaseCorrAll = []; % base corr
    FFchangeAll = []; % Learning
    FrateChangeAll = []; % Frate change
    for m=1:NumMotifs
        % -- base corr
        ff_base = TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1);
        frate_base = TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1);
        [b] = lt_regress(ff_base, frate_base, 0, 0, 1, 1, 'k');
        
        BaseCorrAll = [BaseCorrAll b(2)];
        
        % --- Ff change
        y_basemean = mean(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).FF_All(TMPSTRUCT(m).Is_Baseline==1));
        
        inds = TMPSTRUCT(m).TimeDayvals_All > WindowPeakLearning(1) & ...
            TMPSTRUCT(m).TimeDayvals_All < WindowPeakLearning(2);
        
        ff_zscore = mean((TMPSTRUCT(m).FF_All(inds) - y_basemean)./y_basestd);
        
        FFchangeAll = [FFchangeAll ff_zscore];
        
        % --- Frate change
        y_basemean = mean(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        y_basestd = std(TMPSTRUCT(m).Frate_All(TMPSTRUCT(m).Is_Baseline==1));
        
        frate_zscore = mean((TMPSTRUCT(m).Frate_All(inds) - y_basemean)./y_basestd);
        
        FrateChangeAll = [FrateChangeAll frate_zscore];
        
    end
    
    lt_plot_stem3(BaseCorrAll, FFchangeAll, FrateChangeAll, 'r', 1);
    line([0 0], ylim);
    line(xlim, [0 0]);
    

    
    %% ==== SAVE FIGURES
    cd(NeuronDatabase.global.basedir);

    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(i).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);
    
    try cd('FIGS')
    catch err
        mkdir('FIGS')
        cd('FIGS')
    end
    
    tstamp = lt_get_timestamp(0);
    
    mkdir(['LearningAnaly_' tstamp]);
    cd(['LearningAnaly_' tstamp]);
    lt_save_all_figs;    

end

%% save all figures




%% ALL SYLS IN THE SAME PLOT: 1) FF distribution; 2)
