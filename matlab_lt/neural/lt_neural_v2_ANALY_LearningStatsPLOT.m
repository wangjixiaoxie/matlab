function lt_neural_v2_ANALY_LearningStatsPLOT(MOTIFSTATS_Compiled, convertToZscore, ...
    neuralmetric_toplot, PlotNeuralFFDeviationCorrs, plotOnlyTargSyl, birdtoplot, ...
    expttoplot)

%%

NumBirds = length(MOTIFSTATS_Compiled.birds);
collectLastTertile = 1; % if 1, then uses last tertile of training. if 0, then takes entire training (for calcaulting mean elarning,e tc)

%% COLLECT VARIOUS STATS - PLOT CORRELATIONS ACROSS TIME BETWEEN NEURAL SIMILARITY AND FF


% -- to collect across experiments (one neuron/motif one datapt)
NeuralvsFFSlopeAll = [];
NeuralvsFFInterceptAll = [];
MeanFFRelBaseAll = [];
MeanNeuralSimRelBaseAll = [];
IsTargetAll = [];
SameTypeAll = [];
ThisSylLearnDirAll = [];

NeuralSimBaseAll = [];
NeuralSimTrainAll = [];


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
                
                baseInds = MotifStats.neurons(nn).motif(j).Params.baseInds;
                trainInds = MotifStats.neurons(nn).motif(j).Params.trainInds;
                
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
                neuralsimilarity = [segmentsextract.(neuralmetric_toplot)];
                
                % ==================== 2) FF
                alltrialFF = [segmentsextract.FF_val];
                
             
                
                %% ============ COLLECT METRICS (all metrics)
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
                
                
                % ============= (TRAIN VS. BASE)
                if collectLastTertile==1
                    indstmp = find(trainInds);
                    trainInds_end = indstmp(2*floor(length(indstmp)/3)+1:end);
                    
                    neuralsim_base = nanmedian(neuralsimilarity(baseInds));
                    neuralsim_train = nanmedian(neuralsimilarity(trainInds_end));                
                    neuralsim_meanchange = neuralsim_train - neuralsim_base;
                    FF_meanchange = mean(alltrialFF(trainInds_end)) - mean(alltrialFF(baseInds));
                    
                else % usees bulk, all in training
                     neuralsim_base = nanmedian(neuralsimilarity(baseInds));
                    neuralsim_train = nanmedian(neuralsimilarity(trainInds));                
                   neuralsim_meanchange = neuralsim_train - neuralsim_base;
                    FF_meanchange = mean(alltrialFF(trainInds)) - mean(alltrialFF(baseInds));
                end
                
                
%                 if isnan(neuralsim_meanchange)
%                     keyboard
%                 end
                
                % ======== targ learn dir
                targmotifinds = find(ismember(motiflist, targetmotifs));
                TargLearnDirs = [];
                for zzz=1:length(targmotifinds)
                    indtmp = targmotifinds(zzz);
                    
%                     disp(['----- ' MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency]);
                    
                    if strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Up')
                        targdir = 1;
                    elseif strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Dn')
                        targdir = -1;
                    elseif strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Of') & strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.preTranContingency, ...
                            'Up')
                        targdir = -1;
                    elseif strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            'Of') & strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.preTranContingency, ...
                            'Dn')
                        targdir = 1;
                    end
                    if strcmp(MotifStats.neurons(nn).motif(indtmp).LearningInfo.postTranContingency, ...
                            MotifStats.neurons(nn).motif(indtmp).LearningInfo.preTranContingency)
                        targdir = 0;
                    end
                    TargLearnDirs = [TargLearnDirs targdir];
                end
                
                TargLearnDirs = unique(TargLearnDirs);
                if length(TargLearnDirs)>1
                    % diff targ do diff things, so can;t give this expt
                    % a targ dir
                    TargLearnDirs = nan;
                end
                
                
                % ====================== WHICH SWITCH (WITHIN THIS EXPT) IS
                % THIS NEURON ASSOCIATED WITH? (COULD BE NONE)
                
               
                
                
                
                % ====================== REGRESSION STATS
                neuralFFslope = nan;
                neuralFFint = nan;
                if ~all(isnan(alltrialFF))
                    plotON=0;
                    if PlotNeuralFFDeviationCorrs ==1
                        
                        if plotOnlyTargSyl==0 | any(strcmp(motifname, targetmotifs)) % i.e. only target if I care about that.
                            if strcmp(birdname, birdtoplot) & strcmp(exptname, expttoplot)
                                
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
                        end
                    end
                    
                    if plotON==1
                        plot(alltrialFF(baseInds), neuralsimilarity(baseInds), 'x', 'Color', [0.6 0.6 0.6]);
                    end
                    [~,~,~,~,~, SummaryStats] = lt_regress(neuralsimilarity(trainInds), ...
                        alltrialFF(trainInds), plotON, 0, 1, 1, plotcols{nn}, 0);
                    
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
                    trancode = nan;
                else
                    if strcmp(postTranContingency, 'Up')
                        trancode = 1;
                    elseif strcmp(postTranContingency, 'Dn')
                        trancode = -1;
                    elseif strcmp(postTranContingency, 'Of') & strcmp(preTranContingency, 'Up')
                        trancode = -1;
                    elseif strcmp(postTranContingency, 'Of') & strcmp(preTranContingency, 'Dn')
                        trancode = 1;
                    end
                    
                    if strcmp(preTranContingency, postTranContingency)
                        trancode = 0;
                    end
                end
                ThisSylLearnDirAll = [ThisSylLearnDirAll trancode];
                
                
                
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
                
                
                
                % =========== bulk stats (all of post vs. all of pre)N
                NeuralSimBaseAll = [NeuralSimBaseAll neuralsim_base];
                NeuralSimTrainAll = [NeuralSimTrainAll neuralsim_train];
                
                
                MeanNeuralSimRelBaseAll = [MeanNeuralSimRelBaseAll neuralsim_meanchange];
                MeanFFRelBaseAll = [MeanFFRelBaseAll FF_meanchange];
                
                NeuralvsFFSlopeAll = [NeuralvsFFSlopeAll neuralFFslope];
                NeuralvsFFInterceptAll = [NeuralvsFFInterceptAll neuralFFint];
                
                TargLearnDirsAll = [TargLearnDirsAll TargLearnDirs];
                
                
                % ============= STATS AT THE END OF THE POST PERIOD
                
                
            end
        end
        if ~isempty(hsplots)
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


%% === PLOT SUMMARY - REGRESSIONS
% for each neuron, expt, what is the relationship found between neural and
% FF change?
lt_figure; hold on;

xbins = -0.7:0.1:0.7;

% ======== targ syls
inds = IsTargetAll==1;
lt_subplot(1,3,1); hold on;
title('targ');
xlabel('slope (neural vs. ff, trial by trial');
lt_plot_annotation(1, 'blu: pos learning', 'b')
% up learning
indstmp = ThisSylLearnDirAll==1 & inds;
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'b');

% dn learning
indstmp = ThisSylLearnDirAll==-1 & inds;
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'r');


% ======== same-type syls
inds = IsTargetAll==0 & SameTypeAll==1;
lt_subplot(1,3,2); hold on;
title('same type');

% up learning
indstmp = TargLearnDirsAll==1 & inds & ~isnan(NeuralvsFFSlopeAll)
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'b');

% dn learning
indstmp = TargLearnDirsAll==-1 & inds & ~isnan(NeuralvsFFSlopeAll);
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'r');


% ======== diff syls
inds = IsTargetAll==0 & SameTypeAll==0;
lt_subplot(1,3,3); hold on;
title('diff type');

% up learning
indstmp = TargLearnDirsAll==1 & inds & ~isnan(NeuralvsFFSlopeAll);
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'b');

% dn learning
indstmp = TargLearnDirsAll==-1 & inds & ~isnan(NeuralvsFFSlopeAll);
lt_plot_histogram(NeuralvsFFSlopeAll(indstmp), xbins, 1, 0, '', 1, 'r');



%% === PLOT SUMMARY - REGRESSIONS [combined both dir learning into one plot]
% for each neuron, expt, what is the relationship found between neural and
% FF change?

lt_figure; hold on;
title('slope flipped if learn dir neg');
Xlabels={};

% ======== targ syls
inds = IsTargetAll==1;
X = 1;
Xlabels{X} = 'targ';
NeuralvsFFSlopeAll_flipped = NeuralvsFFSlopeAll;% flip if learning is down
NeuralvsFFSlopeAll_flipped(ThisSylLearnDirAll==-1) = -NeuralvsFFSlopeAll_flipped(ThisSylLearnDirAll==-1);

% - plot
Y = NeuralvsFFSlopeAll_flipped(inds);
% plot(X, Y, 'ok');
distributionPlot(Y', 'showMM', 6, 'addSpread', 1, 'xValues', X, 'histOpt', 1);
p = signrank(Y);
lt_plot_text(X, max(Y), ['p=' num2str(p)], 'r');

% ======== same type syls
inds = IsTargetAll==0 & SameTypeAll==1;
X = 2;
Xlabels{X} = 'same-type';
NeuralvsFFSlopeAll_flipped = NeuralvsFFSlopeAll;
NeuralvsFFSlopeAll_flipped(TargLearnDirsAll==-1) = -NeuralvsFFSlopeAll_flipped(TargLearnDirsAll==-1); % flip if learning is down
NeuralvsFFSlopeAll_flipped(isnan(TargLearnDirsAll)) = nan;

% - plot
Y = NeuralvsFFSlopeAll_flipped(inds);
distributionPlot(Y', 'showMM', 6, 'addSpread', 1, 'xValues', X, 'histOpt', 1)
p = signrank(Y);
lt_plot_text(X, max(Y), ['p=' num2str(p)], 'r');


% ======== diff type syls
inds = IsTargetAll==0 & SameTypeAll==0;
X = 3;
Xlabels{X} = 'diff-type';
NeuralvsFFSlopeAll_flipped = NeuralvsFFSlopeAll;
NeuralvsFFSlopeAll_flipped(TargLearnDirsAll==-1) = -NeuralvsFFSlopeAll_flipped(TargLearnDirsAll==-1); % flip if learning is down
NeuralvsFFSlopeAll_flipped(isnan(TargLearnDirsAll)) = nan;

% - plot
Y = NeuralvsFFSlopeAll_flipped(inds);
distributionPlot(Y', 'showMM', 6, 'addSpread', 1, 'xValues', X, 'histOpt', 1)
p = signrank(Y);
lt_plot_text(X, max(Y), ['p=' num2str(p)], 'r');

% ---
xlim([0 4]);

set(gca, 'XTick', 1:3, 'XTickLabel', Xlabels)
ylabel('slope (neural change vs. ff change (z))')



%% ========== [IN PROGRESS! FIRST DEFINE SWITCH] SUMMARIZE LEARNING (ONE DATAPT = ONE SWITCH/EXPT)
% CURRENTLY EACH NEURON IS ONE DATAPOINT - COULD BE MORE THAN ONE NEURON
% PER SWITCH
% === CURRENTLY ONLY FOR SINGLE CONTEXT TRAINING

exptnums = unique(ExptnumAll);
neuronnums = unique(NeuronNumAll);

LearnAll = []; % neurons x 3 (targ, same, diff)

for i=1:NumBirds
    for ii=1:max(exptnums)
        for iii=1:max(neuronnums)
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
            % -- targ
            inds = IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;
            
            if sum(inds)>1 & length(unique((ThisSylLearnDirAll(inds))))>1
                % then multiple bidir laerning
                continue
            end
            
            
            learndir = unique(ThisSylLearnDirAll(inds));
            LearnTarg = MeanFFRelBaseAll(inds)*learndir;
            LearnTarg = median(LearnTarg);
            
            % -- same type
            inds = IsTargetAll==0 & SameTypeAll==1 ...
                & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;

            LearnSameType = MeanFFRelBaseAll(inds).*learndir;
            LearnSameType = nanmedian(LearnSameType);
            
            % --- diff type
            inds = IsTargetAll==0 & SameTypeAll==0 ...
                & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;

            LearnDiffType = MeanFFRelBaseAll(inds).*learndir;
            LearnDiffType = nanmedian(LearnDiffType);
            

            % === output
            LearnAll = [LearnAll; [LearnTarg LearnSameType LearnDiffType]];
        end
    end
end

lt_figure; hold on;
lt_subplot(1,3,1); hold on;
plot([1 2 3], LearnAll', 'o-k');
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
title('[in progress - each line one neuron');
lt_plot_annotation(1, 'need to change code: one line per switch', 'r');
ylabel('learn (targ dir) (z)');
lt_plot_zeroline;
lt_plot(1.15:3.15, mean(LearnAll,1), {'Errors', lt_sem(LearnAll), 'LineStyle', '-', 'Color', 'b'});


lt_subplot(1,3,2); hold on;
[~, inds] = unique(LearnAll(:,1));
LearnAll = LearnAll(inds,:);
plot([1 2 3], LearnAll', 'o-k');
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
title('[in progress - each line one neuron (totally identical are excluded)');
lt_plot_annotation(1, 'need to change code: one line per switch', 'r');
ylabel('learn (targ dir) (z)');
lt_plot_zeroline;
lt_plot(1.15:3.15, mean(LearnAll,1), {'Errors', lt_sem(LearnAll), 'LineStyle', '-', 'Color', 'b'});


% ===================== SAME AS ABOVE, BUT ONE VALUE PER EXPERIMENT
LearnAll = []; % neurons x 3 (targ, same, diff)

for i=1:NumBirds
    for ii=1:max(exptnums)
        LearnThisExpt = []; % neuron x 3 (t, s, d)
        for iii=1:max(neuronnums)
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
            % -- targ
            inds = IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;
            
            if sum(inds)>1 & length(unique((ThisSylLearnDirAll(inds))))>1
                % then multiple bidir laerning
                continue
            end
            
            
            learndir = unique(ThisSylLearnDirAll(inds));
            LearnTarg = MeanFFRelBaseAll(inds)*learndir;
            LearnTarg = median(LearnTarg);
            
            % -- same type
            inds = IsTargetAll==0 & SameTypeAll==1 ...
                & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;

            LearnSameType = MeanFFRelBaseAll(inds).*learndir;
            LearnSameType = nanmedian(LearnSameType);
            
            % --- diff type
            inds = IsTargetAll==0 & SameTypeAll==0 ...
                & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;

            LearnDiffType = MeanFFRelBaseAll(inds).*learndir;
            LearnDiffType = nanmedian(LearnDiffType);
            
            % --- 
            LearnThisExpt = [LearnThisExpt; [LearnTarg LearnSameType LearnDiffType]];
        end
            % === output
            LearnAll = [LearnAll; median(LearnThisExpt, 1)];
    end
end

lt_subplot(1,3,3); hold on;
plot([1 2 3], LearnAll', 'o-k');
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
title('[each line one expt (median over neurons and switches)');

ylabel('learn (targ dir) (z)');
lt_plot_zeroline;
lt_plot(1.15:3.15, mean(LearnAll,1), {'Errors', lt_sem(LearnAll), 'LineStyle', '-', 'Color', 'b'});




%% ==== CHANGE IN NEURAL FOR TARG VS. SAME VS. DIFF [FOR EACH NEURON]
% ONE DATAPOINT FOR EACH NEURON/MOTIF (neuron = switch)
% takes median over same type and diff type

exptnums = unique(ExptnumAll);
neuronnums = unique(NeuronNumAll);

NeuralAll = []; % trials x 3 (targ, same, diff)
NeuralBase = [];
NeuralTrain = [];
LearningAll = [];
for i=1:NumBirds
    for ii=exptnums
        
        for iii=neuronnums
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
%             birdname = MOTIFSTATS_Compiled.birds(i).birdname;
%             exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
%             
            % -- targ
            inds = find(IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
            
            
            for j=1:length(inds)
                
            neuralbasethis = [];
            neuraltrainthis = [];
                indtmp = inds(j); % ind for the targ
                learndir = ThisSylLearnDirAll(indtmp);
                
                % -- targ
                x = MeanFFRelBaseAll(indtmp)*learndir;
                LearningAll = [LearningAll x];
                
                neuralTarg = MeanNeuralSimRelBaseAll(indtmp);
                
                neuralbasethis= [neuralbasethis NeuralSimBaseAll(indtmp)];
                neuraltrainthis= [neuraltrainthis NeuralSimTrainAll(indtmp)];
                
                
                
                % -- same
                inds_nontarg = find(IsTargetAll==0 & SameTypeAll==1 & ...
                    BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
                
%                 x = MeanFFRelBaseAll(inds_nontarg).*learndir;
                neuralSame = MeanNeuralSimRelBaseAll(inds_nontarg);
                neuralSame = nanmedian(neuralSame);
                
                neuralbasethis= [neuralbasethis nanmedian(NeuralSimBaseAll(inds_nontarg))];
                neuraltrainthis= [neuraltrainthis nanmedian(NeuralSimTrainAll(inds_nontarg))];

                % --- diff
                inds_nontarg = find(IsTargetAll==0 & SameTypeAll==0 & ...
                    BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
                
                neuralDiff = MeanNeuralSimRelBaseAll(inds_nontarg);
                neuralDiff = nanmedian(neuralDiff);

                neuralbasethis= [neuralbasethis nanmedian(NeuralSimBaseAll(inds_nontarg))];
                neuraltrainthis= [neuraltrainthis nanmedian(NeuralSimTrainAll(inds_nontarg))];

                % ==== output
                NeuralAll = [NeuralAll; [neuralTarg neuralSame neuralDiff]];
                
                NeuralBase = [NeuralBase; neuralbasethis];
                NeuralTrain = [NeuralTrain; neuraltrainthis];
            end
        end
    end
end

% ============================= DISTRIBUTISON (BASE VS. TRAIN)
lt_figure; hold on;

hsplot= lt_subplot(1,3,1); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('targ');
lt_plot_45degScatter(NeuralBase(:,1), NeuralTrain(:,1), 'k');


hsplot= lt_subplot(1,3,2); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('same');
lt_plot_45degScatter(NeuralBase(:,2), NeuralTrain(:,2), 'b');


hsplot= lt_subplot(1,3,3); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('diff');
lt_plot_45degScatter(NeuralBase(:,3), NeuralTrain(:,3), 'r');


% ============================= DISTRIBUTISON (BASE VS. TRAIN) [GOOD
% LEARNING]
learndivider = prctile(LearningAll, [66]);
indstokeep = LearningAll>=learndivider;
lt_figure; hold on;

hsplot= lt_subplot(1,3,1); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('targ (top median)');
lt_plot_45degScatter(NeuralBase(indstokeep,1), NeuralTrain(indstokeep,1), 'k');


hsplot= lt_subplot(1,3,2); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('same (top median)');
lt_plot_45degScatter(NeuralBase(indstokeep,2), NeuralTrain(indstokeep,2), 'b');


hsplot= lt_subplot(1,3,3); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('diff (top median)');
lt_plot_45degScatter(NeuralBase(indstokeep,3), NeuralTrain(indstokeep,3), 'r');


% ----- bottom median
indstokeep = LearningAll<learndivider;
lt_figure; hold on;

hsplot= lt_subplot(1,3,1); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('targ (bot median)');
lt_plot_45degScatter(NeuralBase(indstokeep,1), NeuralTrain(indstokeep,1), 'k');


hsplot= lt_subplot(1,3,2); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('same (bot median)');
lt_plot_45degScatter(NeuralBase(indstokeep,2), NeuralTrain(indstokeep,2), 'b');


hsplot= lt_subplot(1,3,3); hold on;
hsplots = [hsplots hsplot];
xlabel('corr (base)'); 
ylabel('corr (train)');

title('diff (bot median)');
lt_plot_45degScatter(NeuralBase(indstokeep,3), NeuralTrain(indstokeep,3), 'r');



% ============================== DIFFERENCES
lt_figure; hold on;

lt_subplot(1,3,1); hold on;
title('all expts');
plot([1 2 3], NeuralAll', '-ok')
lt_plot(1.15:3.15, nanmean(NeuralAll,1), {'Errors', lt_sem(NeuralAll), 'LineStyle', '-', 'Color', 'b'});
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
title('[each line one neuron/motif]');
ylabel('neural change (z)');
lt_plot_zeroline;
% anova for differences
[p, table, stats] = anova1(NeuralAll, [], 'off');
lt_plot_pvalue(p, 'anova', 1);

% ========== REPLOT, BUT THRESHOLDING ONLY FOR EXPERIMENTS WITH CLEAR
% LEARNING

% == top median
lt_subplot(1,3,2); hold on
title('top median of learning');

indstokeep = LearningAll>=learndivider;
NeuralAll_sub = NeuralAll(indstokeep, 1:3);

plot([1 2 3], NeuralAll_sub', '-ok')
lt_plot(1.15:3.15, nanmean(NeuralAll_sub,1), {'Errors', lt_sem(NeuralAll_sub), 'LineStyle', '-', 'Color', 'b'});
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
ylabel('neural change (z)');
lt_plot_zeroline;
[p, table, stats] = anova1(NeuralAll_sub, [], 'off');
lt_plot_pvalue(p, 'anova', 1);

% -- bottom median
lt_subplot(1,3,3); hold on
title('bottom median of learning');

indstokeep = LearningAll<learndivider;
NeuralAll_sub = NeuralAll(indstokeep, 1:3);

plot([1 2 3], NeuralAll_sub', '-ok')
lt_plot(1.15:3.15, nanmean(NeuralAll_sub,1), {'Errors', lt_sem(NeuralAll_sub), 'LineStyle', '-', 'Color', 'b'});
xlim([0 4]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});
ylabel('neural change (z)');
lt_plot_zeroline;
[p, table, stats] = anova1(NeuralAll_sub, [], 'off');
lt_plot_pvalue(p, 'anova', 1);



%% ==== CHANGE IN NEURAL FOR TARG VS. NONTARG [FOR EACH NEURON]
% ONE DATAPOINT FOR EACH NEURON/MOTIF (neuron = switch)
% takes median over same type and diff type

exptnums = unique(ExptnumAll);
neuronnums = unique(NeuronNumAll);

NeuralAll = []; % trials x 3 (targ, same, diff)
LearningAll = [];
for i=1:NumBirds
    for ii=exptnums
        
        for iii=neuronnums
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
%             birdname = MOTIFSTATS_Compiled.birds(i).birdname;
%             exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
%             
            % -- targ
            inds = find(IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
            
            for j=1:length(inds)
                
                indtmp = inds(j); % ind for the targ
                learndir = ThisSylLearnDirAll(indtmp);
                
                % -- targ
                x = MeanFFRelBaseAll(indtmp)*learndir;
                LearningAll = [LearningAll x];
                
                neuralTarg = MeanNeuralSimRelBaseAll(indtmp);
                
                % -- nontarg
                inds_nontarg = find(IsTargetAll==0 & ...
                    BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
                
%                 x = MeanFFRelBaseAll(inds_nontarg).*learndir;
                neuralNont = MeanNeuralSimRelBaseAll(inds_nontarg);
                neuralNont = nanmedian(neuralNont);

                
                % ==== output
                NeuralAll = [NeuralAll; [neuralTarg neuralNont]];
            end
        end
    end
end

lt_figure; hold on;

lt_subplot(1,3,1); hold on;
title('all expts');
plot([1 2], NeuralAll', '-ok')
lt_plot(1.15:2.15, nanmean(NeuralAll,1), {'Errors', lt_sem(NeuralAll), 'LineStyle', '-', 'Color', 'b'});
xlim([0 3]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'targ', 'nontarg'});
title('[each line one neuron/motif]');
ylabel('neural change (z)');
lt_plot_zeroline;


% % ========== REPLOT, BUT THRESHOLDING ONLY FOR EXPERIMENTS WITH CLEAR
% % LEARNING
% learndivider = prctile(LearningAll, [66]);
% 
% % == top median
% lt_subplot(1,3,2); hold on
% title('top median of learning');
% 
% indstokeep = LearningAll>=learndivider;
% NeuralAll_sub = NeuralAll(indstokeep, 1:3);
% 
% assert(size(NeuralAll_sub, 1)~=3, 'problem, plot cant do square matrix');
% plot([1 2], NeuralAll_sub, '-ok')
% lt_plot(1.15:2.15, nanmean(NeuralAll_sub,1), {'Errors', lt_sem(NeuralAll_sub), 'LineStyle', '-', 'Color', 'b'});
% xlim([0 3]);
% set(gca, 'XTick', 1:2, 'XTickLabel', {'targ', 'nontarg'});
% title('[each line one neuron/motif]');
% ylabel('neural change (z)');
% lt_plot_zeroline;
% 


%% =========== ONE PLOT OVER ALL EXPERIMENTS, FF CHANGE, AND NEURAL CHANGE, FOR ALL SYLS

lt_figure; hold on;
xlabel('targ (mean if mult)');
ylabel('nontarg');
% title('mean over neurons');
title('Change in neural sim (each dot one neuron)');

exptnums = unique(ExptnumAll);
neuronnums = unique(NeuronNumAll);

for i=1:NumBirds
    for ii=exptnums
        for iii=neuronnums
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
            % -- targ
            inds = IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;
            
            X = nanmean(MeanNeuralSimRelBaseAll(inds));
            
            % -- nontarg
            inds = IsTargetAll==0 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii;
            Y = MeanNeuralSimRelBaseAll(inds);
            
%             disp(X)
%             disp(Y)
%             disp('--')
            
            plot(X,nanmedian(Y), 'o');
            %         plot(X,Y, 'ok');
            
        end
    end
end
line([-3 1], [-3 1])
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% ========== WITHIN EACH EXPT, PLOT NEURAL VS. LEARN

figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


exptnums = unique(ExptnumAll);
neuronnums = unique(NeuronNumAll);
for i=1:NumBirds
    for ii=exptnums
        
        for iii=neuronnums
            
            if ~any(BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii)
                continue
            end
            
        birdname = MOTIFSTATS_Compiled.birds(i).birdname;
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
            % -- targ
            inds = find(IsTargetAll==1 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
            
            for j=1:length(inds)
                indtmp = inds(j);
                learndir = ThisSylLearnDirAll(indtmp);
                
                % -- one plot for each targ
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                targsyl = 'FILL IN';
                title([birdname '-' exptname '-' num2str(iii) '-' targsyl]);
                xlabel('ff change (targlearndir)');
                ylabel('neural change');
                
                x = MeanFFRelBaseAll(indtmp)*learndir;
                y = MeanNeuralSimRelBaseAll(indtmp);
                
                lt_plot(x,y, {'Color', 'k'});
                
                % -- one plot for each nontarg
                inds_nontarg = find(IsTargetAll==0 & BirdnumAll==i & ExptnumAll==ii & NeuronNumAll==iii);
                
                x = MeanFFRelBaseAll(inds_nontarg).*learndir;
                y = MeanNeuralSimRelBaseAll(inds_nontarg);
                
                plot(x,y, 'ob');
                
                xlim([-3 3]);
                ylim([-3 3]);
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
            end
            
        end
    end
end

        
%% =========== PLOT SUMMARY ACROSS EXPERIMENTS [COMBINE ALL SYLS, NEURAL VS. FF CHANGE]

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



%% =========== PLOT SUMMARY ACROSS EXPERIMENTS [COMBINE ALL SYLS, ABSOLUTE CHANGES, FF AND NEURAL]
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
