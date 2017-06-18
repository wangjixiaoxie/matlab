%% lt 5/22/17 - across all birds, expts, neurons, plot learning related statistics (see other code
% for visualization of raw data)
function lt_neural_v2_ANALY_LearningStats(MOTIFSTATS_Compiled, LearningMetadat, ...
    PlotTimeCourses, PlotNeuralFFDeviationCorrs, convertToZscore)
% NOTES

% DEFAULT: plots each neuron z-score relative to itself (i.e. its own
% baseline) (i.e. if multiple epochs, then does not use earliest baseline
% for a given neuron)
% TO DO: plot each neuron z-score realtive to the earliest baseline for
% itself

%%
% for smoothing fr
window_neural = 0.0075; % neural;
windshift_neural = 0.002;


premotorWind = [-0.075 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset

% === plotting similairty scores
plotRawScores=0; % if 0, just plots running avg
% convertToZscore = 1; % converts all scores (e.g neural, FF) to zscore relative to baseline
runbin = 15; % num trials to smooth over.

%%
NumBirds = length(MOTIFSTATS_Compiled.birds);



%% GET MEAN FR - go thru all experiments

for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    
    for ii=1:numexpts
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        for j=1:nummotifs
            
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            
            for nn=1:numneurons
                
                
                % ==== For each rendition, get smoothed firing rate
                clustnum = MotifStats.neurons(nn).clustnum;
                assert(clustnum == SummaryStruct.birds(1).neurons(nn).clustnum, 'asdfaf');
                
                MotifStats.neurons(nn).motif(j).SegmentsExtract = ...
                    lt_neural_SmoothFR(MotifStats.neurons(nn).motif(j).SegmentsExtract, ...
                    clustnum);
                
                
            end
        end
        
        % === stick back into main structure
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS = MotifStats;
        
    end
end


%% EXTRACT LEARNING METRIC AND OTHER THINGS

for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    for ii=1:numexpts
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        assert(length(SummaryStruct.birds)==1, 'asdfasd');
        birdname = SummaryStruct.birds(1).birdname;
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        FirstDay = []; % will fill in with first neuron/motif - other enurons use same first day.
        for j=1:nummotifs
            motifname = motiflist{j};
            targetmotifs = MotifStats.params.TargSyls;
            
            
            
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            for nn=1:numneurons
                
                segmentsextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                
                if ~isfield(segmentsextract, 'FRsmooth_xbin')
                    continue % i.e. no trials
                end
                
                % what time bins contain premotor window?
                tmp = segmentsextract(1).FRsmooth_xbin;
                premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                    & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                
                
                % what trials are baseline?
                WNonDnum = datenum(SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr, ...
                    'ddmmmyyyy-HHMM');
                baseInds = [segmentsextract.song_datenum] < WNonDnum;
                
                
                
                % =================================== Dir of learning (for
                % the one transition)
                indbird = strcmp({LearningMetadat.bird.birdname}, birdname);
                indcol = strcmp(LearningMetadat.bird(indbird).info(1,:), exptname) ...
                    & strcmp(LearningMetadat.bird(indbird).info(2,:), motifname); % column (expt and syl)
                
                transitionDates = {};
                preTranContingency = [];
                postTranContingency = []; % for all syls this is exmpty, UNLESS you are target
                
                assert(sum(indcol)<2, 'PROBBLEM !! - why more than one col');
                if sum(indcol)==1
                    transitionDates = LearningMetadat.bird(indbird).info(3:end, indcol);
                    % ---- FIND transition date corresponding to this neuron
                    for ddd = 1:length(transitionDates)
                        if strcmp(transitionDates{ddd}(1:14), ...
                                SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr)
                            % YES !!! FOUND TRANSITION INFORMATION
                            preTranContingency = transitionDates{ddd}(16:17);
                            postTranContingency = transitionDates{ddd}(19:20);
                            break
                        end
                    end
                    assert(~isempty(preTranContingency), 'DID NOT ENTER ALL TRANSITIONS INTO LEARNINGMETSTRUCT (OR made mistake)');
                end
                % note: if this is a targ, then will have
                % preTranContingency and postTranContingency. Otherwise
                % both will be empty.
                
                
                
                % -------------------------------------- trial times
                alltrialTimes = [segmentsextract.song_datenum];
                if isempty(FirstDay)
                    FirstDay = datestr(min(alltrialTimes), 'ddmmmyyyy');
                end
                tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, alltrialTimes);
                alltrialTimes = tmp.FinalValue;
                
                
                % ============== 1) FIRING RATE (NEURAL SIMILARITY TO BASELINE MEAN
                % FR)
                alltrialFR = [segmentsextract.FRsmooth_rate_CommonTrialDur]; % all trial FR (bin x trial)
                alltrialFR = alltrialFR(premotorInds, :); % extract just premotor window
                
                baseFR = alltrialFR(:, baseInds);
                baseFR_mean = mean(baseFR, 2);
                
                % ===== for each trial, calculate deviation from mean base
                % FR
                neuralsimilarity = corr(alltrialFR, baseFR_mean); % trials x 1
                
                % ==================== 2) FF
                alltrialFF = [segmentsextract.FF_val];
                
                
                % ========================== output
                tmp = num2cell(neuralsimilarity);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.neuralsimilarity] = deal(tmp{:});
                
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).Params.baseInds = baseInds;
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.preTranContingency = preTranContingency;
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.postTranContingency = postTranContingency;
                
                
            end
        end
    end
end


%% PLOT ALL EXPERIMENTS, TIMECOURSES

if PlotTimeCourses ==1
    
    for i=1:NumBirds
        numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
        for ii=1:numexpts
            
            MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
            SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
            
            motiflist = MotifStats.params.motif_regexpr_str;
            assert(length(SummaryStruct.birds)==1, 'asdfasd');
            birdname = SummaryStruct.birds(1).birdname;
            exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
            
            % ================ EACH bird, expt, syl
            nummotifs = length(motiflist);
            
            figcount=1;
            subplotcols=2;
            subplotrows=4;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            FirstDay = []; % will fill in with first neuron/motif - other enurons use same first day.
            for j=1:nummotifs
                motifname = motiflist{j};
                targetmotifs = MotifStats.params.TargSyls;
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                if any(strcmp(motifname, targetmotifs))
                    title([motifname '(targ)'], 'Color', 'r');
                else
                    title(motifname);
                end
                
                % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
                % INTERVALS FOR NEURONS
                
                % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
                numneurons = length(MotifStats.neurons);
                plotcols = lt_make_plot_colors(numneurons, 0, 0);
                for nn=1:numneurons
                    
                    segmentsextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    
                    if ~isfield(segmentsextract, 'FRsmooth_xbin')
                        continue % i.e. no trials
                    end
                    
                    % what trials are baseline?
                    baseInds = MotifStats.neurons(nn).motif(j).Params.baseInds;
                    
                    % --------------------------------- Dir of learning (for
                    % the one transition)
                    preTranContingency = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.preTranContingency;
                    postTranContingency = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.postTranContingency;
                    
                    % -------------------------------------- trial times
                    alltrialTimes = [segmentsextract.song_datenum];
                    if isempty(FirstDay)
                        FirstDay = datestr(min(alltrialTimes), 'ddmmmyyyy');
                    end
                    tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, alltrialTimes);
                    alltrialTimes = tmp.FinalValue;
                    
                    
                    % ============== 1) FIRING RATE (NEURAL SIMILARITY TO BASELINE MEAN
                    % FR)
                    neuralsimilarity = [segmentsextract.neuralsimilarity];
                    
                    % ==================== 2) FF
                    alltrialFF = [segmentsextract.FF_val];
                    
                    
                    %% ============ CONVERT TO ZSCORE (all metrics)
                    if convertToZscore==1
                        
                        % -- neural
                        basemean = mean(neuralsimilarity(baseInds));
                        basestd = std(neuralsimilarity(baseInds));
                        neuralsimilarity = (neuralsimilarity - basemean)./basestd;
                        
                        % -- learning
                        if ~all(isnan(alltrialFF))
                            basemean = mean(alltrialFF(baseInds));
                            basestd = std(alltrialFF(baseInds));
                            
                            alltrialFF = (alltrialFF - basemean)./basestd;
                        end
                    end
                    
                    
                    %% PLOT THIS NEUORN
                    if length(alltrialFF)>runbin
                        % ====== NEURAL SIMILARITY
                        if plotRawScores ==1
                            plot(alltrialTimes, neuralsimilarity, 'o');
                        end
                        Yrunning =lt_running_stats(neuralsimilarity, runbin);
                        xrunning = lt_running_stats(alltrialTimes, runbin);
                        shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', plotcols{nn}}, 1);
                        
                        % ===== FF
                        Yrunning =lt_running_stats(alltrialFF, runbin);
                        xrunning = lt_running_stats(alltrialTimes, runbin);
                        shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', 'k'}, 1);
                        
                        
                        % ---
                        lastBaseTrial = alltrialTimes(max(find(baseInds)));
                        line([lastBaseTrial lastBaseTrial], ylim, 'Color', 'k');
                        if convertToZscore==1
                            ylim([-3 3]);
                        else
                            ylim([-1 1]);
                        end
                        lt_plot_zeroline
                    end
                    
                    % ==== if is target, then indicate contingency
                    if (0) % REASON: plotting all transitions, not just the one WN on for analysis,.
                        if ~isempty(preTranContingency)
                            lt_plot_text(lastBaseTrial, 3, [preTranContingency '-' postTranContingency], ...
                                'm');
                        end
                    end
                    
                    % ------ sanity check, learning metastruct correspnd to
                    % what I think
                    assert(any(strcmp(motifname, targetmotifs)) == ~isempty(postTranContingency), 'asdafsd');
                    
                    % ================ PUT LINES FOR ALL CHANGES FOR THIS MOTIF
                    indbird = strcmp({LearningMetadat.bird.birdname}, birdname);
                    indcol = strcmp(LearningMetadat.bird(indbird).info(1,:), exptname) ...
                        & strcmp(LearningMetadat.bird(indbird).info(2,:), motifname); % column (expt and syl)
                    transitionDates = {};
                    
                    assert(sum(indcol)<2, 'PROBBLEM !! - why more than one col');
                    if sum(indcol)==1
                        transitionDates = LearningMetadat.bird(indbird).info(3:end, indcol);
                        
                        % --- for each date, plot line
                        for ddd = 1:length(transitionDates)
                            dateString = transitionDates{ddd}(1:14);
                            tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, {dateString});
                            line([tmp.FinalValue tmp.FinalValue], ylim, 'LineStyle' , '--' , 'Color', 'r');
                            
                            lt_plot_text(tmp.FinalValue, 3, [transitionDates{ddd}(16:17) '-' transitionDates{ddd}(19:20)], ...
                                'r');
                        end
                    end
                    
                end
            end
            linkaxes(hsplots, 'xy');
            lt_subtitle([birdname '-' exptname])
        end
    end
end




%% COLLECT VARIOUS STATS - PLOT CORRELATIONS ACROSS TIME BETWEEN NEURAL SIMILARITY AND FF

    
    % -- to collect across experiments (one neuron/motif one datapt)
    NeuralvsFFSlopeAll = [];
    NeuralvsFFInterceptAll = [];
    MeanFFRelBaseAll = [];
    MeanNeuralSimRelBaseAll = [];
    IsTargetAll = [];
    SameTypeAll = [];
    LearnDirAll = [];
    
    BirdnumAll = [];
    ExptnumAll = [];
    NeuronNumAll = [];
    MotifNumAll = [];
    TargLearnDirsAll=[];
    
    for i=1:NumBirds
        numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
        for ii=1:numexpts
            
            MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
            SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
            
            motiflist = MotifStats.params.motif_regexpr_str;
            sametypeSyls = MotifStats.params.SameTypeSyls;
            birdname = SummaryStruct.birds(1).birdname;
            exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
            
            % ================ EACH bird, expt, syl
            nummotifs = length(motiflist);
            
            figcount=1;
            subplotcols=3;
            subplotrows=6;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            FirstDay = []; % will fill in with first neuron/motif - other enurons use same first day.
            for j=1:nummotifs
                motifname = motiflist{j};
                targetmotifs = MotifStats.params.TargSyls;
                
                
                
                % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
                numneurons = length(MotifStats.neurons);
                plotcols = lt_make_plot_colors(numneurons, 0, 0);
                for nn=1:numneurons
                    
                    segmentsextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    
                    if ~isfield(segmentsextract, 'FRsmooth_xbin')
                        continue % i.e. no trials
                    end
                    
                    % --------------------------------- Dir of learning (for
                    % the one transition)
                    preTranContingency = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.preTranContingency;
                    postTranContingency = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.postTranContingency;
                    
                    % -------------------------------------- trial times
                    alltrialTimes = [segmentsextract.song_datenum];
                    if isempty(FirstDay)
                        FirstDay = datestr(min(alltrialTimes), 'ddmmmyyyy');
                    end
                    tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, alltrialTimes);
                    alltrialTimes = tmp.FinalValue;
                    
                    
                    % ============== 1) FIRING RATE (NEURAL SIMILARITY TO BASELINE MEAN
                    % FR)
                    neuralsimilarity = [segmentsextract.neuralsimilarity];
                    
                    % ==================== 2) FF
                    alltrialFF = [segmentsextract.FF_val];
                    
                    
                    %% ============ CONVERT TO ZSCORE (all metrics)
                    if convertToZscore==1
                        baseInds = MotifStats.neurons(nn).motif(j).Params.baseInds;
                        
                        % -- neural
                        basemean = mean(neuralsimilarity(baseInds));
                        basestd = std(neuralsimilarity(baseInds));
                        neuralsimilarity = (neuralsimilarity - basemean)./basestd;
                        
                        % -- learning
                        if ~all(isnan(alltrialFF))
                            basemean = mean(alltrialFF(baseInds));
                            basestd = std(alltrialFF(baseInds));
                            
                            alltrialFF = (alltrialFF - basemean)./basestd;
                        end
                    end
                    
                    
                    % ============= CALCULATE MEAN NEURAL SIM AND FF DURING WN (COMARED TO BASE)
                    baseInds = MotifStats.neurons(nn).motif(j).Params.baseInds;
                    neuralsim_meanchange = mean(neuralsimilarity(~baseInds)) - mean(neuralsimilarity(baseInds));
                    FF_meanchange = mean(alltrialFF(~baseInds)) - mean(alltrialFF(baseInds));
                    
                    
                    
                    
                                    % ======== targ learn dir
                targmotifinds = find(ismember(motiflist, targetmotifs));
                TargLearnDirs = [];
                for zzz=1:length(targmotifinds)
                    indtmp = targmotifinds(zzz);
                    
                    
                    disp(['----- ' MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency]);
                    
                    if strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Up')
                        TargLearnDirs = [TargLearnDirs 1];
                    elseif strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Dn')
                        TargLearnDirs = [TargLearnDirs -1];
                    elseif strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Of')
                        TargLearnDirs = [TargLearnDirs 0];
                    end
                end
                
                TargLearnDirs = unique(TargLearnDirs);
                if length(TargLearnDirs)>1
                    % diff targ do diff things, so can;t give this expt
                    % a targ dir
                    TargLearnDirs = nan;
                end

                    
                    
                    %% PLOT THIS NEUORN
                    neuralFFslope = nan;
                    neuralFFint = nan;
                    if ~all(isnan(alltrialFF))
                        plotON=0;
                        if PlotNeuralFFDeviationCorrs ==1
                            
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            if any(strcmp(motifname, targetmotifs))
                                title([motifname '(targ, ' postTranContingency ')'], 'Color', 'r');
                            else
                                title(motifname);
                            end
                            plotON = 1;
                        xlim([-3 3]);
                        ylim([-3 3]);
                        end
                        
                        
                        [~,~,~,~,~, SummaryStats] = lt_regress(neuralsimilarity, ...
                            alltrialFF, 0, 0, 0, 0, plotcols{nn}, 0);
                        
                        neuralFFint = SummaryStats.intercept;
                        neuralFFslope = SummaryStats.slope;
                    end
                    
                    %% ================ Collect stuff
                    
                    BirdnumAll = [BirdnumAll i];
                    ExptnumAll = [ExptnumAll ii];
                    NeuronNumAll = [NeuronNumAll nn];
                    MotifNumAll = [MotifNumAll j];
                    
                    if ~isempty(postTranContingency)
                        disp([postTranContingency '-' num2str(TargLearnDirs)])
                    end
                    
                    if isempty(postTranContingency)
                        LearnDirAll = [LearnDirAll nan];
                    else
                        if strcmp(postTranContingency, 'Up')
                            LearnDirAll = [LearnDirAll 1];
                        elseif strcmp(postTranContingency, 'Dn')
                            LearnDirAll = [LearnDirAll -1];
                        elseif strcmp(postTranContingency, 'Of')         
                            LearnDirAll = [LearnDirAll 0];
                            
                        end
                            
                    end
                    
                    
                    
                    if any(strcmp(motifname, sametypeSyls))
                        SameTypeAll = [SameTypeAll 1];
                    else
                        SameTypeAll = [SameTypeAll 0];
                    end
                    
                    if any(strcmp(motifname, targetmotifs))
                        IsTargetAll = [IsTargetAll 1];
                    else
                        IsTargetAll = [IsTargetAll 0];
                    end
                    
                    MeanNeuralSimRelBaseAll = [MeanNeuralSimRelBaseAll neuralsim_meanchange];
                    MeanFFRelBaseAll = [MeanFFRelBaseAll FF_meanchange];
                    
                    NeuralvsFFSlopeAll = [NeuralvsFFSlopeAll neuralFFslope];
                    NeuralvsFFInterceptAll = [NeuralvsFFInterceptAll neuralFFint];
                    
                    TargLearnDirsAll = [TargLearnDirsAll TargLearnDirs];
                    
                end
            end
            if PlotNeuralFFDeviationCorrs==1
            linkaxes(hsplots, 'xy');
            lt_subtitle([birdname '-' exptname])
            end
        end
    end


%% == troubleshooting
if (0)
i=2;
ii=1;
nn=1;
j=6

MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls
MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_regexpr_str

MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.postTranContingency
end


%% =========== PLOT SUMMARY ACROSS EXPERIMENTS

% === NEURAL VS. LEARN
lt_figure; hold on;
hsplots = [];
% -- targ
hsplot = lt_subplot(2,2,1); hold on;
hsplots = [hsplots hsplot];
title('targ');
xlabel('FF change (targdir)');
ylabel('neural change')
inds = IsTargetAll==1;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
lt_regress(Y, X, 1, 0, 1, 1, 'k', 0);
lt_plot_zeroline
lt_plot_zeroline_vert


% -- sametype
hsplot = lt_subplot(2,2,2); hold on;
hsplots = [hsplots hsplot];
title('same type');
inds = IsTargetAll==0 & SameTypeAll==1;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
if ~isempty(Y)
lt_regress(Y, X, 1, 0, 1, 1, 'b', 0);
lt_plot_zeroline
lt_plot_zeroline_vert
end
% --- diff type
hsplot = lt_subplot(2,2,3); hold on;
hsplots = [hsplots hsplot];
title('diff type');
inds = IsTargetAll==0 &   SameTypeAll==0;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
lt_regress(Y, X, 1, 0, 1, 1, 'r', 0);
lt_plot_zeroline
lt_plot_zeroline_vert

linkaxes(hsplots, 'xy');

% ======== ABS (NEURAL VS. LEARN)
lt_figure; hold on;
hsplots = [];
% -- targ
hsplot = lt_subplot(2,2,1); hold on;
hsplots = [hsplots hsplot];
title('targ');
xlabel('abs(FF change)');
ylabel('neural change')
inds = IsTargetAll==1;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
lt_regress(Y, abs(X), 1, 0, 1, 1, 'k', 0);
lt_plot_zeroline
lt_plot_zeroline_vert


% -- sametype
hsplot = lt_subplot(2,2,2); hold on;
hsplots = [hsplots hsplot];
title('same type');
inds = IsTargetAll==0 & SameTypeAll==1;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
lt_regress(Y, abs(X), 1, 0, 1, 1, 'b', 0);
lt_plot_zeroline
lt_plot_zeroline_vert

% --- diff type
hsplot = lt_subplot(2,2,3); hold on;
hsplots = [hsplots hsplot];
title('diff type');
inds = IsTargetAll==0 &   SameTypeAll==0;
learndir = TargLearnDirsAll(inds);

X = MeanFFRelBaseAll(inds) .* learndir;
Y = MeanNeuralSimRelBaseAll(inds);
lt_regress(Y, abs(X), 1, 0, 1, 1, 'r', 0);
lt_plot_zeroline
lt_plot_zeroline_vert

linkaxes(hsplots, 'xy');
