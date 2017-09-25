function [DATSTRUCT, BYNEURONDAT] = lt_neural_v2_ANALY_Swtch_Summary(MOTIFSTATS_Compiled, SwitchStruct, RemoveLowNumtrials, ...
    MinTrials, UseZscoreNeural, neuralmetricname, fieldname_baseneur, fieldname_trainneur, ...
    skipMultiDir, usePeakLearn, plotLearnStatsOn, learnsigalpha, OnlyKeepSigLearn, ...
    OnlyKeepWNonset, OnlyUseDatOnSwitchDay)

% usePeakLearn = 1;
% plotLearnStatsOn =1;

%% TO DO:
% LEARN FOR FOR NONNTARGET ARE STILL NOT IN CORRECT DIR - HAVE TO TAKE
% INTOA ACCOUNT MULTIPLE TARGETS

%%
% 1) corr between basemean vs. train mean ("change" is dat vs. shuffle)
% fieldname_baseneur = 'AllNeurCorrShuffMean';
% fieldname_trainneur = 'AllNeurCorrDat';

% 2) each trial corr with a base "template". change is mean of train corr vals vs.
% mean of base corr vals
% fieldname_baseneur = 'AllNeurSimBase';
% fieldname_trainneur = 'AllNeurSimTrain';

% 3) split base into two, to then get a base corr (mean of base 1 vs. mean
% of base 2) and a trining corr (mean of train vs. mean of base, averaged
% over diff pairs).
% fieldname_baseneur = 'AllNeurSplitCorrBase';
% fieldname_trainneur = 'AllNeurSplitCorrTrain';


%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);

WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

FFsmthbinsize = 10;

minrends = 5; % for both train and base

premotorWind = SwitchStruct.params.premotorWind;

numtrain = 25; % trials to get at end of training (will match base and train. will take smallest common factor)


%%
if mod(FFsmthbinsize,2)==0
    FFsmthbinsize= FFsmthbinsize+1; % conver to odd, so that median tval is at an actual datapoint.
end

%% plot, for each switch, timecourse of neural and FF for all syl types [IN PROGRESS]
% if (0)
%     for i=1:Numbirds
%         
%         numexpts = length(SwitchStruct.bird(i).exptnum);
%         birdname = SwitchStruct.bird(i).birdname;
%         
%         for ii=1:numexpts
%             exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
%             numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
%             
%             MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
%             SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
%             
%             motiflist = MotifStats.params.motif_regexpr_str;
%             targsyls = MotifStats.params.TargSyls;
%             nummotifs = length(motiflist);
%             
%             WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
%                 MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
%             
%             for iii=1:numswitches
%                 
%                 figcount=1;
%                 subplotrows=5;
%                 subplotcols=6;
%                 fignums_alreadyused=[];
%                 hfigs=[];
%                 
%                 
%                 goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] ...
%                     & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs]);
%                 
%                 if isempty(goodneurons)
%                     disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
%                     continue
%                 end
%                 
%                 swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
%                 swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
%                 swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
%                 
%                 plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
%                 
%                 
%                 % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
%                 for j=1:nummotifs
%                     
%                     for nn=goodneurons
%                         
%                         segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
%                         
%                         if ~isfield(segextract, 'spk_Times')
%                             continue
%                         end
%                         
%                         baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
%                         trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
%                         
%                         if length(baseInds)<minrends | length(trainInds) < minrends
%                             % even for good neurosn, could occur if some motifs
%                             % labeled pre but not post.
%                             continue
%                         end
%                         
%                         trialstoplot = [baseInds trainInds];
%                         clustnum = MotifStats.neurons(nn).clustnum;
%                         
%                         
%                         % ================= PLOT TIMECOURSE OF NEURAL AND FF
%                         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                         
%                         % ----- FF
%                         ffvals = [segextract.FF_val];
%                         tvals = [segextract.song_datenum];
%                         % -- convert tvals to days
%                         tmpday = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
%                         tvals = lt_convert_EventTimes_to_RelTimes(datestr(tmpday, 'ddmmmyyyy'),...
%                             tvals);
%                         tvals = tvals.FinalValue;
%                         
%                         tsmth = lt_running_stats(tvals, FFsmthbinsize);
%                         
%                         if ~all(isnan(ffvals));
%                             % ---- get zscore
%                             ffvals_basemean = mean(ffvals(baseInds));
%                             ffvalsbaseSD = std(ffvals(baseInds));
%                             
%                             ffvals = (ffvals - ffvals_basemean)./ffvalsbaseSD;
%                             
%                             ffsmth = lt_running_stats(ffvals, FFsmthbinsize);
%                             
%                             plot(tsmth.Median, ffsmth.Median, 'kx');
%                             %                       plot(tvals(trialstoplot), ffvals(trialstoplot), 'xk');
%                             %
%                         end
%                         
%                         % ---- neural (corr with baseline)
%                         neuralsim = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).(neuralmetricname);
%                         %                    if ~strcmp(neuralmetricname, 'NEURvsbase_FRcorr')
%                         % then get zscore
%                         neurbasemean = mean(neuralsim(baseInds));
%                         neurbaseSD = std(neuralsim(baseInds));
%                         neuralsim = (neuralsim - neurbasemean)./neurbaseSD;
%                         %                    end
%                         neursmth = lt_running_stats(neuralsim, FFsmthbinsize);
%                         %                     lt_plot(tsmth.Mean, neursmth.Mean, {'Errors', neursmth.SEM, ...
%                         %                         'Color', plotcols{j}});
%                         plot(tsmth.Median, neursmth.Mean, 'o', 'Color', plotcols{nn});
%                         %                    plot(tvals(trialstoplot), neuralsim(trialstoplot), 'ob');
%                         
%                         
%                         % --- stuff
%                         axis tight;
%                         ylim([-3 3]);
%                         lt_plot_zeroline;
%                         
%                         % --- line for base vs. training
%                         line([tvals(max(baseInds)) tvals(max(baseInds))], ylim, 'Color','k', 'LineWidth', 2);
%                         
%                         % --- title
%                         if any(strcmp(targsyls, motiflist{j}))
%                             title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}], 'Color', 'r');
%                         else
%                             title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}]);
%                         end
%                         
%                         % ======================== COLLECT DATA FOR PLOTTING
%                         
%                     end
%                 end
%             end
%         end
%     end
% end

%% PLOT SUMMARY ACROSS ALL SWITCHES [IN PROGRESS - JUST STARTED]

AllFFchange = [];
AllNeurSimChange = [];
AllNeurSimBase = [];
AllNeurSimTrain = [];
AllFFneurCorr = [];

AllTargLearnDir = [];
AllNumtargs = [];
AllTargssamesyls = [];
AllTargSameDir = [];

AllIsTarg = [];
AllIsSame = [];
AllIsDiff = [];

AllBirdnum = [];
AllExptnum = [];
AllSwitchnum = [];
AllSwitchnumGlobal = [];
AllNeuronnum = [];
AllMotifnum = [];
% AllSwitchCounter = [];

AllNumTrials = [];

AllBaseFFvsNeurFRCorr = [];
AllBaseFFvsNeurFRCorr_p = [];
AllBaseFFvsNeursimCorr = [];
AllBaseFFvsNeursimCorr_p = [];

AllSNRbaseEnd = [];
AllSNRtrainEnd = [];

AllNeurCorrDat = []; % corr(mean1 vs. mean2)
AllNeurCorrShuffMean = []; % shuffled ...
AllNeurCorrShuffSD = [];

AllNeurSplitCorrBase = []; % (base1 vs. base2) split into two halves
AllNeurSplitCorrTrain = []; % (train vs. base) - multiple, then take mean.
AllNeurSplitCorrTrainvsTrain = []; % train vs. train;

% --- laerning stuff
AllNeurTargLearnRate_targdir = [];
AllNeurTargLearnEnd_targdir = [];


% counter =1;
switchnum_global=1;

if plotLearnStatsOn==1
    figure(1); hold on % plots all trajectories
    figure(2); hold on; % plots population stats
end

for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        SameSyls = MotifStats.params.SameTypeSyls;
        DiffSyls = MotifStats.params.DiffTypeSyls;
        
        targsyls = MotifStats.params.TargSyls;
        nummotifs = length(motiflist);
        
        WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
            MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            %             figcount=1;
            %             subplotrows=5;
            %             subplotcols=6;
            %             fignums_alreadyused=[];
            %             hfigs=[];
            
            %             goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] ...
            %                 & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs]);
            
            
            if OnlyKeepWNonset==1
                tmptmp = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies(2:2:end);
               tmptmp = cell2mat(tmptmp');
               if ~all(tmptmp(:,1)==0)
                   disp('SKIPPED - not starting from WN off')
                   continue
               end
            elseif OnlyKeepWNonset==2
                tmptmp = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies(2:2:end);
               tmptmp = cell2mat(tmptmp');
               if any(tmptmp(:,1)==0)
                   disp('SKIPPED - starting from WN off')
                   continue
               end
            end
            
            goodneurons = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons;
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
            swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            % - things about target
            numtargs = length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs)/2;
            targssamesyl = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl;
            
            if length(unique([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}])) ==1
                TargsSameDir = 1;
            else
                TargsSameDir = 0;
            end
            
            if TargsSameDir ==0
                targlearndir = nan;
            else
                targlearndir = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2};
            end
            
            
            % --- skip if multiple targs and in different directions
            if skipMultiDir ==1
                if numtargs>1 & TargsSameDir==0
                    disp(['SKIPPING ' birdname '-' exptname '-sw' num2str(iii) ' [multidir]']);
                    continue
                end
            end
            
            % ----------- TIME OF PEAK LEARNING
            % --- will take first target
            tsyltmp = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{1};
            targdirtmp = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2};
            
            targind = strcmp(tsyltmp, ...
                {SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl.sylname});
            
            
            tvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(targind).tvals;
            ffvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(targind).ffvals;
            
            indstmp = tvals > SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            
            % - learn
            tvals = tvals(indstmp);
            ffvals = ffvals(indstmp);
            
            tvals = lt_running_stats(tvals, numtrain);
            ffvals = lt_running_stats(ffvals, numtrain);
            
            [~, indtmp] = max(targdirtmp.*ffvals.Mean);
            TimeMaxLearn = tvals.Median(indtmp);
            
            
            %% learning stats at target
            tvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(targind).tvals;
            ffvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(targind).ffvals;
            
            trainindstmp = tvals > SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            % --- how good was learning for this experiment?
            % 1) learn rate (hz/time) - can't use rends because not all
            % labeled
            x = (tvals - floor(tvals(1))) * 24; % convert from datenum to hours (relative to first day)
            [b,~,~,~,~,SummaryStats] = lt_regress(ffvals(trainindstmp), ...
                x(trainindstmp), 0, 0, 0, 0, 0, 0, learnsigalpha);
            learnrate = SummaryStats.slope;
            learnrateCI = SummaryStats.slopeCI;
            
            
            % 2) area under curve of learning (relative to base mean)
            % i.e. mean deviation from baseline (across all datapoints)
            learnMeanDevFromBase = mean(ffvals(trainindstmp) - mean(ffvals(~trainindstmp)));
            
            
            % 3) ff at end of training minus base
            ntrialstmp = min([sum(trainindstmp) numtrain]);
            learnTrainEnd = mean(ffvals(end-ntrialstmp+1:end)) - mean(ffvals(~trainindstmp));
            
            
            % --- make targ dir
            learnrate_targdir = targdirtmp * learnrate;
            learnMeanDevFromBase_targdir = targdirtmp * learnMeanDevFromBase;
            learnTrainEnd_targdir = targdirtmp * learnTrainEnd;
            
            
            % ------ IS LEARNING SIGNIFICANT?
            learnsig_regression = all(learnrateCI*targdirtmp >0); % if CI is not in 0;
            
            if targdirtmp ==1
                [~, learnsig_endmean] = ranksum(ffvals(end-ntrialstmp+1:end), ...
                    ffvals(~trainindstmp), 'tail', 'right', 'alpha', learnsigalpha);
            else
                [~, learnsig_endmean] = ranksum(ffvals(end-ntrialstmp+1:end), ...
                    ffvals(~trainindstmp), 'tail', 'left', 'alpha', learnsigalpha);
            end
            
            % ------ plot learning + stats for learning
            if plotLearnStatsOn ==1
                
                figure(1); cla;
                title([birdname '-' exptname '-sw' num2str(iii)])
                hold on
                plot(tvals, ffvals, 'ok');
                plot(tvals(trainindstmp), ffvals(trainindstmp), 'or');
                line([TimeMaxLearn TimeMaxLearn], ylim, 'Color', 'r');
                line(xlim, [mean(ffvals(~trainindstmp)) mean(ffvals(~trainindstmp))], 'Color', 'k');
                lt_plot_text(TimeMaxLearn, mean(ffvals(~trainindstmp)), ...
                    ['maxlearn, dir: ' num2str(targdirtmp)], 'r');
                
                lt_regress(ffvals(trainindstmp), tvals(trainindstmp), 1, 0, 1, 0, 'r', 0);
                
                lt_plot_text(tvals(1), max(ffvals), ['learnrate' num2str(learnrateCI, '%3.2g')], 'k')
                
                line(xlim, mean(ffvals(~trainindstmp))+ [learnMeanDevFromBase learnMeanDevFromBase], 'Color', 'b');
                lt_plot_text(tvals(end), mean(ffvals(~trainindstmp))+learnMeanDevFromBase, 'learn(mean dev)', 'b');
                
                line(xlim, mean(ffvals(~trainindstmp))+ [learnTrainEnd learnTrainEnd], 'Color', 'm');
                lt_plot_text(tvals(end),mean(ffvals(~trainindstmp))+ learnTrainEnd, 'learn(train end)', 'm');
                
                lt_plot_text(tvals(end-5), max(ffvals), ['regress sig? ' num2str(learnsig_regression)], 'b');
                lt_plot_text(tvals(end-5), max(ffvals)-20, ['trainend sig? ' num2str(learnsig_endmean)], 'b');
                
                
                
                pause
                
                
                
                figure(2)
                subplot(221); hold on;
                xlabel('learnrate');
                ylabel('learn mean dev');
                plot(learnrate_targdir, learnMeanDevFromBase_targdir, 'ok');
                
                
                subplot(222); hold on;
                xlabel('learnrate');
                ylabel('learnTrainend');
                plot(learnrate_targdir, learnTrainEnd_targdir, 'or');
                
                
                subplot(223); hold on;
                xlabel('learnTrainend');
                ylabel('learnMeandev');
                plot(learnTrainEnd_targdir, learnMeanDevFromBase_targdir, 'ob');
                
                
            end
            
            
            % -----------------------------
            if OnlyKeepSigLearn ==1
                if learnsig_regression == 0 & learnsig_endmean ==0;
                    continue
                end
            end
            
            %%
            % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
            for j=1:nummotifs
                
                sylname = motiflist{j};
                
                
                istarg = any(strcmp(sylname, targsyls));
                issame = any(strcmp(sylname, SameSyls));
                isdiff = any(strcmp(sylname, DiffSyls));
                
                %                     if i==3 & ii==1 & iii==1 & j==11
                %                         keyboard
                %                     end
                
                %                 if strcmp(birdname, 'bu77wh13') & strcmp(exptname, 'LMANlearn1') ...
                %                         & iii==2 & istarg==1
                %                     keyboard
                %                 end
                
                
                for nn=goodneurons
                    
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    
                    if ~isfield(segextract, 'spk_Times')
                        continue
                    end
                    
                    
                    if OnlyUseDatOnSwitchDay==1
                     baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds_WithinDayOfSw);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds_WithinDayOfSw);
                    else
                     baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    end
                    
                    if length(baseInds)<minrends | length(trainInds) < minrends
                        % even for good neurosn, could occur if some motifs
                        % labeled pre but not post.
                        continue
                    end
                    
                    trialstoplot = [baseInds trainInds];
                    
                    
                    % ================================================ COLLECT LEARING/NEURAL
                    % STATS
                    
                    % ----- FF (learning)
                    ffvals = [segextract.FF_val];
                    tvals = [segextract.song_datenum];
                    
                    
                    % ---- get zscore
                    if ~any(isnan(ffvals))
                        ffvals_basemean = mean(ffvals(baseInds));
                        ffvalsbaseSD = std(ffvals(baseInds));
                        
                        ffvals = (ffvals - ffvals_basemean)./ffvalsbaseSD;
                    end
                    % ---- neural (corr with baseline)
                    neuralsim = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).(neuralmetricname);
                    
                    if UseZscoreNeural==1
                        % then get zscore
                        neurbasemean = mean(neuralsim(baseInds));
                        neurbaseSD = std(neuralsim(baseInds));
                        neuralsim = (neuralsim - neurbasemean)./neurbaseSD;
                    end
                    
                    
                    % --------------------- FOR DIFFERENCES, FIGURE OUT
                    % BASE END AND TRAIN END
                    
                    % --------------------- A) N TRIALS
                    NumTrials = min([numtrain, ceil(length(trainInds)/2), length(baseInds)]);
                    %                     disp(['numtrials common to base and train (max ' numtrain '): ' num2str(NumTrials)]);
                    
                    % 1) base inds
                    baseEndInds = baseInds(end-NumTrials+1:end);
                    % 2) training inds
                    if usePeakLearn==0
                        % then take end
                        trainEndInds = trainInds(end-NumTrials+1:end);
                    else
                        % take peak learning for target
                        % find trial that is closest to peak learning time
                        indtmp = find(tvals<TimeMaxLearn, 1, 'last');
                        
                        indfirsttmp = indtmp - ceil(NumTrials/2) +1; % beginning of window
                        indfirsttmp = max([indfirsttmp max(baseInds)+1]); % make sure don't get any base trials
                        
                        trainEndInds = indfirsttmp:indfirsttmp+NumTrials-1;
                        
                        % shift back if is at the end
                        if trainEndInds(end) - trainInds(end) >0
                            % shift training Inds back
                            trainEndInds = trainInds(end-NumTrials+1:end);
                        end
                    end
                    
                    % ---------------------------- B) 2*N TRIALS - for
                    % spliiting analyses
                    NumTrials2 = min([numtrain*2, ceil(length(trainInds)/2), length(baseInds)]);
                    
                    % 1) base inds
                    baseEndInds_longer = baseInds(end-NumTrials2+1:end);
                    % 2) training inds
                    if usePeakLearn==0
                        % then take end
                        trainEndInds_longer = trainInds(end-NumTrials2+1:end);
                    else
                        % then take peak learning for target
                        % find trial that is closest to peak learning time
                        indtmp = find(tvals<TimeMaxLearn, 1, 'last');
                        
                        indfirsttmp = indtmp - ceil(NumTrials2/2) +1; % beginning of window
                        indfirsttmp = max([indfirsttmp max(baseInds)+1]); % make sure don't get any base trials
                        
                        trainEndInds_longer = indfirsttmp:indfirsttmp+NumTrials2-1;
                        
                        % shift back if is at the end
                        if trainEndInds_longer(end) - trainInds(end) >0
                            % shift training Inds back
                            trainEndInds_longer = trainInds(end-NumTrials2+1:end);
                        end
                    end
                    
                    
                    
                    
                    
                    
                    %%
                    % =================== METRICS
                    % -- ff
                    ffchange = mean(ffvals(trainEndInds)) - mean(ffvals(baseEndInds));
                    
                    
                    % -- neur sim
                    neursimbase = mean(neuralsim(baseEndInds));
                    neursimtrain = mean(neuralsim(trainEndInds));
                    
                    neursimchange = neursimtrain - neursimbase;
                    
                    if isnan(neursimtrain)
                        disp('NOTE!! there is something with nan for neural similairty - will lead to entire neuron/motif being thrown out. solve by reruning lt_neural_v2_ANALY_Swtch_Extract with interpolate');
                        
                    end
                    
                    % -- corr betwee neural and ff
                    if ~all(isnan(ffvals))
                        if all(size(ffvals) == size(neuralsim))
                            FFneurCorr = corr(ffvals(trainInds)', neuralsim(trainInds)');
                        else
                            FFneurCorr = corr(ffvals(trainInds)', neuralsim(trainInds));
                        end
                    else
                        FFneurCorr = nan;
                    end
                    
                    % --- baseline FF vs. neur corr
                    AllBaseFFvsNeurFRCorr = [AllBaseFFvsNeurFRCorr ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsMeanFR];
                    AllBaseFFvsNeurFRCorr_p = [AllBaseFFvsNeurFRCorr_p ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsMeanFR_pval];
                    AllBaseFFvsNeursimCorr = [AllBaseFFvsNeursimCorr ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsNeurSim];
                    AllBaseFFvsNeursimCorr_p = [AllBaseFFvsNeursimCorr_p ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsNeurSim_pval];
                    
                    
                    
                    %% ========= MORE METRICS
                    
                    premotorInds = MotifStats.params.premotorInds_FR;
                    alltrialFR = [segextract.FRsmooth_rate_CommonTrialDur];
                    alltrialFR = alltrialFR(premotorInds, :);
                    alltrialFR = double(alltrialFR);
                    
                    
                    % ==== 1) mean of (train) corr vs. (base) - this would
                    % be my best estimate of the correlation between signal
                    % from pre vs. post (better than taking corr on trial
                    % basis and then averaging)
                    basemeanFR = mean(alltrialFR(:, baseEndInds), 2);
                    trainmeanFR = mean(alltrialFR(:, trainEndInds), 2);
                    
                    NeurCorr_dat = corr(basemeanFR, trainmeanFR);
                    
                    % --- permutation to get mean + error
                    allFRtmp = alltrialFR(:, [baseEndInds trainEndInds]);
                    Ncycles = 15;
                    vecsize = length(baseEndInds); assert(length(baseEndInds) == length(trainEndInds), 'dasdf');
                    Z = [];
                    for mm = 1:Ncycles
                        inds = randperm(2*vecsize);
                        inds1 = inds(1:end/2);
                        inds2 = inds(end/2+1:end);
                        
                        x1 = mean(allFRtmp(:, inds1), 2);
                        x2 = mean(allFRtmp(:, inds2), 2);
                        
                        Z(mm) = corr(x1,x2);
                    end
                    NeurCorr_shuff = mean(Z);
                    NeurCorr_shuff_std = std(Z);
                    
                    AllNeurCorrDat = [AllNeurCorrDat NeurCorr_dat];
                    AllNeurCorrShuffMean = [AllNeurCorrShuffMean NeurCorr_shuff];
                    AllNeurCorrShuffSD = [AllNeurCorrShuffSD NeurCorr_shuff_std];
                    
                    %                     disp(NeurCorr_shuff); disp(NeurCorr_dat);
                    %                     disp('-');
                    
                    % ==== 2) base1 vs. base 2, then do base1 vs. train;
                    whichmethod = 1;
                    
                    % -------------- method 1 - new, shuffles
                    if whichmethod ==1
                        ntrials = floor(length(baseEndInds_longer)/2);
                        ncycles = 20;
                        
                        basecorrall = [];
                        trainvsbasecorrall = [];
                        traincorrall = [];
                        tic
                        for mm = 1:ncycles
                            
                            indtmp = randperm(length(baseEndInds_longer));
                            inds1 = indtmp(1:ntrials);
                            inds2 = indtmp(ntrials+1:ntrials*2);
                            
                            
                            % - base similarkty
                            base1 = mean(alltrialFR(:, baseEndInds_longer(inds1)),2);
                            base2 = mean(alltrialFR(:, baseEndInds_longer(inds2)), 2);
                            
                            basecorrall = [basecorrall corr(base1, base2)];
                            
                            
                            % - train vs. base
                            train2 = mean(alltrialFR(:, trainEndInds_longer(inds2)), 2);
                            
                            trainvsbasecorrall = [trainvsbasecorrall corr(base1, train2)];
                            
                            
                            % - train vs. train
                            train1 = mean(alltrialFR(:, trainEndInds_longer(inds1)), 2);
                            traincorrall = [traincorrall corr(train1, train2)];
                            
                        end
                        if (0)
                            disp(['splitcorr: ' num2str(mean(basecorrall), '%3.2g') ' - ' num2str(mean(trainvsbasecorrall), '%3.2g') ...
                            ' - ' num2str(mean(traincorrall), '%3.2g')]);
                        end
                        
                        % -- for use use mean over within base and within
                        % train. 
%                         AllNeurSplitCorrBase = [AllNeurSplitCorrBase  mean([basecorrall traincorrall])]; 
                        AllNeurSplitCorrBase = [AllNeurSplitCorrBase  mean(basecorrall)]; 
                        AllNeurSplitCorrTrain = [AllNeurSplitCorrTrain mean(trainvsbasecorrall)];
                        AllNeurSplitCorrTrainvsTrain = [AllNeurSplitCorrTrainvsTrain mean(traincorrall)];
                        
                        % ----------------- method 2 - old
                    elseif whichmethod ==2
                        ntrials = floor(length(baseEndInds_longer)/2);
                        
                        base1Inds = baseEndInds_longer(1:ntrials);
                        base2Inds = baseEndInds_longer(end-ntrials+1:end);
                        train1Inds = trainEndInds_longer(1:ntrials);
                        train2Inds = trainEndInds_longer(end-ntrials+1:end);
                        
                        base1mean = mean(alltrialFR(:, base1Inds),2);
                        base2mean = mean(alltrialFR(:, base2Inds), 2);
                        train1mean = mean(alltrialFR(:, train1Inds), 2);
                        train2mean = mean(alltrialFR(:, train2Inds), 2);
                        
                        % - base corr
                        x1 = corr(base1mean, base2mean);
                        
                        % - train corr
                        x2 = corr(base1mean, train1mean);
                        x3 = corr(base2mean, train1mean);
                        x4 = corr(base1mean, train2mean);
                        x5 = corr(base2mean, train2mean);
                        
                        
                        AllNeurSplitCorrBase = [AllNeurSplitCorrBase x1]; % (base1 vs. base2) split into two halves
                        AllNeurSplitCorrTrain = [AllNeurSplitCorrTrain mean([x2 x3 x4 x5])]; % (train vs. base) - multiple, then take mean.
                        
                    end
                    
                    %                     disp(x1); disp(x2); disp(x3); disp(x4); disp(x5);
                    %                     disp(mean([x2 x3 x4 x5]));
                    % figure
                    % hold on; plot(base1mean)
                    % hold on; plot(base2mean, 'r')
                    % hold on; plot(trainSplitmean, 'k')
                    %
                    
                    
                    
                    
                    %% diagnostics
                    
                    %                     % === 1) SNR for different trial bins
                    %                     segextract.
                    plotsubset = 0; % just when debugging. ..
                    
%                     premotorInds = MotifStats.params.premotorInds_FR;
%                     FRmat = [segextract.FRsmooth_rate_CommonTrialDur];
%                     FRmat = FRmat(premotorInds, :);
                    
                    % -- baseline
                    if plotsubset==1
                        plotOn = rand>0.90;
                    else
                        plotOn=0;
                    end
                    [SNR] = lt_neural_v2_SNR(alltrialFR(:, baseEndInds), plotOn);
                    AllSNRbaseEnd = [AllSNRbaseEnd SNR];
                    
                    
                    % -- trainEnd
                    [SNR] = lt_neural_v2_SNR(alltrialFR(:, trainEndInds), plotOn);
                    AllSNRtrainEnd = [AllSNRtrainEnd SNR];
                    
                    if plotOn==1
                        lt_plot_annotation(4, ['baseSim=' num2str(neursimbase), ...
                            ', baseFFcorr=' num2str(FFneurCorr)], 'b');
                        pause
                        close all;
                    end
                    
                    %%
                    
                    % ========================== OUTPUTS
                    AllFFchange = [AllFFchange ffchange];
                    AllNeurSimChange = [AllNeurSimChange neursimchange];
                    AllNeurSimBase = [AllNeurSimBase neursimbase];
                    AllNeurSimTrain = [AllNeurSimTrain neursimtrain];
                    
                    AllFFneurCorr = [AllFFneurCorr FFneurCorr];
                    
                    AllTargLearnDir = [AllTargLearnDir targlearndir];
                    AllNumtargs = [AllNumtargs numtargs];
                    AllTargssamesyls = [AllTargssamesyls targssamesyl];
                    AllTargSameDir = [AllTargSameDir TargsSameDir];
                    
                    AllBirdnum = [AllBirdnum i];
                    AllExptnum = [AllExptnum ii];
                    AllSwitchnum = [AllSwitchnum iii];
                    AllNeuronnum = [AllNeuronnum nn];
                    AllMotifnum = [AllMotifnum j];
                    
                    
                    AllIsTarg = [AllIsTarg istarg];
                    AllIsSame = [AllIsSame issame];
                    AllIsDiff = [AllIsDiff isdiff];
                    
                    AllNumTrials = [AllNumTrials NumTrials];
                    
                    %                     AllSwitchCounter = [AllSwitchCounter counter];
                    
                    AllSwitchnumGlobal = [AllSwitchnumGlobal switchnum_global];
                    
                    % ------ related to learning at target
                    AllNeurTargLearnRate_targdir = [AllNeurTargLearnRate_targdir learnrate_targdir];
                    AllNeurTargLearnEnd_targdir = [AllNeurTargLearnEnd_targdir learnTrainEnd_targdir];
                    
                    
                end
            end
            
            % -- increment switch counter
            switchnum_global = switchnum_global+1;
        end
    end
end


% ============= OUTPUT STRUCT
DATSTRUCT.AllFFchange = AllFFchange;

DATSTRUCT.AllNeurSimChange = AllNeurSimChange;
DATSTRUCT.AllNeurSimBase = [AllNeurSimBase ];
DATSTRUCT.AllNeurSimTrain = [AllNeurSimTrain ];

DATSTRUCT.AllFFneurCorr = [AllFFneurCorr ];

DATSTRUCT.AllTargLearnDir = [AllTargLearnDir ];
DATSTRUCT.AllNumtargs = [AllNumtargs ];
DATSTRUCT.AllTargssamesyls = [AllTargssamesyls ];
DATSTRUCT.AllTargSameDir = [AllTargSameDir ];

DATSTRUCT.AllBirdnum = [AllBirdnum ];
DATSTRUCT.AllExptnum = [AllExptnum ];
DATSTRUCT.AllSwitchnum = [AllSwitchnum];
DATSTRUCT.AllSwitchnumGlobal = [AllSwitchnumGlobal];
DATSTRUCT.AllNeuronnum = [AllNeuronnum];
DATSTRUCT.AllMotifnum = [AllMotifnum];


DATSTRUCT.AllIsTarg = [AllIsTarg ];
DATSTRUCT.AllIsSame = [AllIsSame ];
DATSTRUCT.AllIsDiff = [AllIsDiff ];

DATSTRUCT.AllNumTrials = [AllNumTrials ];

DATSTRUCT.AllBaseFFvsNeurFRCorr = AllBaseFFvsNeurFRCorr;
DATSTRUCT.AllBaseFFvsNeurFRCorr_p = [AllBaseFFvsNeurFRCorr_p];
DATSTRUCT.AllBaseFFvsNeursimCorr = [AllBaseFFvsNeursimCorr];
DATSTRUCT.AllBaseFFvsNeursimCorr_p = [AllBaseFFvsNeursimCorr_p];

DATSTRUCT.AllNeurCorrDat = AllNeurCorrDat;
DATSTRUCT.AllNeurCorrShuffMean = [AllNeurCorrShuffMean ];
DATSTRUCT.AllNeurCorrShuffSD = [AllNeurCorrShuffSD ];

DATSTRUCT.AllNeurSplitCorrBase = [AllNeurSplitCorrBase ]; % (base1 vs. base2) split into two halves
DATSTRUCT.AllNeurSplitCorrTrain = [AllNeurSplitCorrTrain ]; % (train vs. base) - multiple, then take mean.
DATSTRUCT.AllNeurSplitCorrTrainvsTrain = [AllNeurSplitCorrTrainvsTrain ]; % (train vs. base) - multiple, then take mean.

DATSTRUCT.AllSNRbaseEnd = [AllSNRbaseEnd ];
DATSTRUCT.AllSNRtrainEnd = [AllSNRtrainEnd ];

DATSTRUCT.AllNeurTargLearnRate_targdir = [AllNeurTargLearnRate_targdir ];
DATSTRUCT.AllNeurTargLearnEnd_targdir = [AllNeurTargLearnEnd_targdir];



if plotLearnStatsOn==1
    subplot(221);
    lt_plot_zeroline; lt_plot_zeroline_vert;
    
    subplot(222);
    lt_plot_zeroline; lt_plot_zeroline_vert;
    subplot(223);
    lt_plot_zeroline; lt_plot_zeroline_vert;
end

%% ===== display sample size

% numbirds = length(unique(DATSTRUCT.AllBirdnum))
% numswitches = length(unique(DATSTRUCT.AllSwitchnumGlobal))
% 
% tmp = tabulate(([num2str(DATSTRUCT.AllSwitchnumGlobal') num2str(DATSTRUCT.AllNeuronnum')]));
% size(tmp)


%% sanity check
all(sum([AllIsSame==0; AllIsDiff==0; AllIsTarg==0],1)==2); % everything is either (xor) targ, same, or diff




%% ======== remove low num trials?

if RemoveLowNumtrials==1
    %     indstoremove = AllNumTrials<MinTrials;
    %     disp(['--- REMOVED ' num2str(sum(indstoremove)) '/' num2str(length(indstoremove)) ' neurons (low N)']);
    indstokeep = find(DATSTRUCT.AllNumTrials>=MinTrials);
    disp(['--- KEPT ' num2str(length(indstokeep)) '/' num2str(length(AllNumTrials)) ' cases (fail to pass min N)']);
    
    DATSTRUCT = lt_structure_subsample_all_fields(DATSTRUCT, indstokeep);
    
end



%% == get zscore version of difference of corr of mean

DATSTRUCT.AllNeurCorrDiff_z = (DATSTRUCT.AllNeurCorrDat - DATSTRUCT.AllNeurCorrShuffMean)./DATSTRUCT.AllNeurCorrShuffSD;



%% ========== COMPARE DIFF METRICS ACROSS ALL DATA

% --- do things correlate with neursim (base)?
lt_figure; hold on;

% 1)
lt_subplot(5,2,1); hold on;
xlabel('base neur sim');
ylabel('base SNR');

plot(DATSTRUCT.AllNeurSimBase.^2, DATSTRUCT.AllSNRbaseEnd, 'x');
xlim([-0.5 1]); ylim([-0.5 1]);

% 2)
lt_subplot(5,2,2); hold on;
xlabel('train neur sim');
ylabel('SNR during training trials');
plot(DATSTRUCT.AllNeurSimTrain, DATSTRUCT.AllSNRtrainEnd, 'x');
xlim([-0.5 1]); ylim([-0.5 1]);

% 2)
lt_subplot(5,2,3); hold on;
xlabel('change in neural sim');
ylabel('change in SNR');
plot(DATSTRUCT.AllNeurSimChange, [DATSTRUCT.AllSNRtrainEnd - ...
    DATSTRUCT.AllSNRbaseEnd], 'x');
% NOTE: is interesting that see many cases without reduction in SNR, but
% clear reduction in neural sim

lt_subplot(5,2,4); hold on;
xlabel('change in neural sim');
ylabel('diff in corr of mean (dat minus shuff mean)');
plot(DATSTRUCT.AllNeurSimChange, [DATSTRUCT.AllNeurCorrDat- ...
    DATSTRUCT.AllNeurCorrShuffMean], 'xr');


lt_subplot(5,2,5); hold on;
xlabel('change in neural sim');
ylabel('diff in corr of mean (dat minus shuff mean) (zscored)');
plot(DATSTRUCT.AllNeurSimChange, DATSTRUCT.AllNeurCorrDiff_z, 'xr');


lt_subplot(5,2,6); hold on;
xlabel('neural sim');
ylabel('corr (mean vs. mean)');
title('bk: base (shuff for corr); red: train')
plot(DATSTRUCT.AllNeurSimBase, DATSTRUCT.AllNeurCorrShuffMean, 'xk');

plot(DATSTRUCT.AllNeurSimTrain, DATSTRUCT.AllNeurCorrDat, 'xr');

% -
lt_subplot(5,2,7); hold on;
xlabel('change in neural sim');
ylabel('diff in corr, using split data (train vs. base)');
plot(DATSTRUCT.AllNeurSimChange, [DATSTRUCT.AllNeurSplitCorrTrain...
    - DATSTRUCT.AllNeurSplitCorrBase] , 'xr');


%% ========= PLOT (ONE LINE FOR EACH EXPT) - NEURAL SIMILAIRTY CHANGES

% - for each neuron, go thru all motifs and compile
numbirds = max(unique([DATSTRUCT.AllBirdnum]));
numexpts = max(unique([DATSTRUCT.AllExptnum]));
numswitches = max(unique([DATSTRUCT.AllSwitchnum]));
numneurons = max(unique([DATSTRUCT.AllNeuronnum]));

figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

ByNeuron_NeuralsimPre = {}; % neuron x {[targ], [same], [diff]}
ByNeuron_NeuralsimPost = {};
ByNeuron_NeuralsimChange = {};
ByNeuron_FFchange = {};

ByNeuron_TargSameSylSameDir = []; % neuron x 1
ByNeuron_TargLearnDir = [];
ByNeuron_SwitchCounter = [];
ByNeuron_Birdname = {};
ByNeuron_Exptname = {};
ByNeuron_SWnum_real = [];

% ByNeuron_SwitchCounter = {};
swcounter = 1;
for i=1:numbirds
    for ii=1:numexpts
        for iii=1:numswitches
            
            if ~any([DATSTRUCT.AllBirdnum]==i & [DATSTRUCT.AllExptnum]==ii ...
                    & [DATSTRUCT.AllSwitchnum]==iii)
                continue
            end
            
            numuniqueneurons = length(unique(DATSTRUCT.AllNeuronnum([DATSTRUCT.AllBirdnum]==i & [DATSTRUCT.AllExptnum]==ii ...
                & [DATSTRUCT.AllSwitchnum]==iii)));
            
            birdname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            
            % === PLOT EACH NEURON OVERLAYED
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' exptname '-sw' num2str(iii)]);
            
            plotcolors = lt_make_plot_colors(numuniqueneurons, 0, 0);
            counter = 1;
            for nn=1:numneurons
                % inds to all motifs for this neuron
                inds = find([DATSTRUCT.AllBirdnum]==i & [DATSTRUCT.AllExptnum]==ii ...
                    & [DATSTRUCT.AllSwitchnum]==iii & [DATSTRUCT.AllNeuronnum]==nn);
                
                if ~any(inds)
                    continue
                end
                
                disp(['FOUND NEURON...']);
                %              Xtarg = [];
                %              Xsame = [];
                %              Xdiff = [];
                %
                %              Ytarg = [];
                %              Ysame = [];
                %              Ydiff = [];
                Xall= cell(1,3); % {targ, same, diff}
                Yall = cell(1,3);
                YminusXall = cell(1,3);
                
                FFchangeAll = cell(1,3);
                
                
                for j=1:length(inds)
                    indtmp = inds(j);
                    
                    %                     x = DATSTRUCT.AllNeurSimBase(indtmp);
                    %                     y = DATSTRUCT.AllNeurSimTrain(indtmp);
                    x = DATSTRUCT.(fieldname_baseneur)(indtmp);
                    y = DATSTRUCT.(fieldname_trainneur)(indtmp);
                    
                    istarg = DATSTRUCT.AllIsTarg(indtmp);
                    issame = DATSTRUCT.AllIsSame(indtmp);
                    isdiff = DATSTRUCT.AllIsDiff(indtmp);
                    
                    assert((istarg + issame + isdiff) ==1, 'asdasdasdf');
                    
                    ffchange = DATSTRUCT.AllFFchange(indtmp);
                    
                    % --- OUTPUT
                    Xall{logical([istarg issame isdiff])} = [Xall{logical([istarg issame isdiff])} ...
                        x];
                    Yall{logical([istarg issame isdiff])} = [Yall{logical([istarg issame isdiff])} ...
                        y];
                    
                    YminusXall{logical([istarg issame isdiff])} = [YminusXall{logical([istarg issame isdiff])} ...
                        double(y-x)];
                    
                    FFchangeAll{logical([istarg issame isdiff])} = ...
                        [FFchangeAll{logical([istarg issame isdiff])} ffchange];
                    %                 if istarg==1
                    %                     Xtarg = [Xtarg x];
                    %                     Ytarg = [];
                    %                 elseif issame==1
                    %
                    %                 elseif isdiff==1
                    %
                    %                 end
                end
                
                % ==== PLOT FOR THIS NEURON
                xvals = [1.3 2.3 3.3] - 0.6/nn;
                lt_plot_MultDist(YminusXall, xvals, 0, plotcolors{counter}, 1);
                counter = counter+1;
                
                % ====== SAVE FOR THIS NEURON
                ByNeuron_NeuralsimPre = [ByNeuron_NeuralsimPre; Xall]; % neuron x {[targ], [same], [diff]}
                ByNeuron_NeuralsimPost = [ByNeuron_NeuralsimPost; Yall];
                ByNeuron_NeuralsimChange = [ByNeuron_NeuralsimChange; YminusXall];
                
                ByNeuron_FFchange = [ByNeuron_FFchange; FFchangeAll];
                
                % ----------------- TARGET STATS
                if SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl ==1 ...
                        & length(unique([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}]))==1
                    % either only one targ, or if mult targs, then all are same
                    % syl and same dir
                    
                    ByNeuron_TargSameSylSameDir = [ByNeuron_TargSameSylSameDir 1]; % neuron x 1
                    ByNeuron_TargLearnDir = [ByNeuron_TargLearnDir ...
                        unique([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}])];
                else
                    %                     pause
                    ByNeuron_TargSameSylSameDir = [ByNeuron_TargSameSylSameDir 0]; % neuron x 1
                    ByNeuron_TargLearnDir = [ByNeuron_TargLearnDir nan];
                end
                
                % ------ other things
                ByNeuron_SwitchCounter = [ByNeuron_SwitchCounter swcounter];
                ByNeuron_Birdname = [ByNeuron_Birdname birdname];
                ByNeuron_Exptname = [ByNeuron_Exptname exptname];
                ByNeuron_SWnum_real = [ByNeuron_SWnum_real iii];
                
                
            end
            set(gca, 'XTick', [1 2 3]);
            lt_plot_zeroline;
            
            swcounter =swcounter+1;
            
        end
    end
end
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_plot_annotation(1, 'deltaNeuralSim vs. (targ, same, diff) - color=neuron; dot=motif');
linkaxes(hsplots, 'xy');


% ================ OUTPUT
BYNEURONDAT = struct;
BYNEURONDAT.ByNeuron_NeuralsimPre = ByNeuron_NeuralsimPre;
BYNEURONDAT.ByNeuron_NeuralsimPost = ByNeuron_NeuralsimPost;
BYNEURONDAT.ByNeuron_NeuralsimChange = ByNeuron_NeuralsimChange;
BYNEURONDAT.ByNeuron_FFchange = ByNeuron_FFchange;

BYNEURONDAT.ByNeuron_TargSameSylSameDir = ByNeuron_TargSameSylSameDir'; % neuron x 1
BYNEURONDAT.ByNeuron_TargLearnDir = ByNeuron_TargLearnDir';
BYNEURONDAT.ByNeuron_SwitchCounter = ByNeuron_SwitchCounter';

BYNEURONDAT.ByNeuron_Birdname = ByNeuron_Birdname';
BYNEURONDAT.ByNeuron_Exptname = ByNeuron_Exptname';
BYNEURONDAT.ByNeuron_SWnum_real = ByNeuron_SWnum_real';



%% ==========[ NEURON AS DATAPOINT]

% ===== ALL DAT
lt_neural_v2_ANALY_Swtch_Summary_c(BYNEURONDAT);



% ==== ONLY THOSE EXPT WITH GOOD LEARNING
if (0)
    [~, inds] = unique(BYNEURONDAT.ByNeuron_SwitchCounter);
    
    func = @(X)nanmean(X);
    learningall = cellfun(func, BYNEURONDAT.ByNeuron_FFchange(inds,1)).*BYNEURONDAT.ByNeuron_TargLearnDir(inds); % one value for learning  at targ for each switch
    
    threshlearn = nanmedian(learningall);
    
    % -- plot
    indstokeep = cellfun(func, BYNEURONDAT.ByNeuron_FFchange(:,1)) >= threshlearn;
    if sum(indstokeep)>0
        BYNEURONDAT_tmp = lt_structure_subsample_all_fields(BYNEURONDAT, indstokeep, 1);
        lt_neural_v2_ANALY_Swtch_Summary_c(BYNEURONDAT_tmp);
    end
end



%% ************************* NEURON/MOTIF AS DATAPOINT *******************
%% ===== PLOT DISTRIBUTIONS OF NEURAL SIMILARITY (PRE AND POST)

% --- TARG, SAME, DIFF
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT, fieldname_baseneur, fieldname_trainneur);

% ---- JUST TARG AND NONTARG
lt_neural_v2_ANALY_Swtch_Summary_b(DATSTRUCT, fieldname_baseneur, fieldname_trainneur);


%% === 1) plot distributions of baseline corrs
lt_figure; hold on

% -- targ
lt_subplot(1,3,1); hold on;
title('targ')
xlabel('ff vs. mean FR corr');
ylabel('ff vs. neural sim corr');
inds = DATSTRUCT.AllIsTarg==1;
plotcol = 'k';

X = DATSTRUCT.AllBaseFFvsNeurFRCorr(inds);
Y = DATSTRUCT.AllBaseFFvsNeursimCorr(inds);

lt_plot_45degScatter(X, Y, plotcol);

% -- same
lt_subplot(1,3,2); hold on;
title('same')
xlabel('ff vs. mean FR corr');
ylabel('ff vs. neural sim corr');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1;
plotcol = 'b';

X = DATSTRUCT.AllBaseFFvsNeurFRCorr(inds);
Y = DATSTRUCT.AllBaseFFvsNeursimCorr(inds);

lt_plot_45degScatter(X, Y, plotcol);


% -- diff
lt_subplot(1,3,3); hold on;
title('diff')
xlabel('ff vs. mean FR corr');
ylabel('ff vs. neural sim corr');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsDiff==1;
plotcol = 'r';

X = DATSTRUCT.AllBaseFFvsNeurFRCorr(inds);
Y = DATSTRUCT.AllBaseFFvsNeursimCorr(inds);

lt_plot_45degScatter(X, Y, plotcol);

%% DISTRIBUTIONS OF NEURAL SIMILARITY [SEPARATE HIGH AND LOW BASELINE PITCH CORR]
% ---- USING CORR WITH MEAN FF
% ===== SIGNIFICANT BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
inds = find(DATSTRUCT.AllBaseFFvsNeurFRCorr_p<0.05);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('significant base FF vs. neur corr (p<0.05)');

% ===== HIGH BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr));
inds = find(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr)>=medianCorr);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('high base ff vs. neurFR corr');

% ===== LOW BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
inds = find(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr)<medianCorr);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('low base ff vs. neurFR corr');



if (0)
    % ===== HIGH BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
    medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeursimCorr));
    inds = find(abs(DATSTRUCT.AllBaseFFvsNeursimCorr)>=medianCorr);
    DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
    lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
    
    % ===== LOW BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
    medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeursimCorr));
    inds = find(abs(DATSTRUCT.AllBaseFFvsNeursimCorr)<medianCorr);
    DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
    lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
end



% ======================= USING BASELINE CORR AS CRITERION
medianval = nanmedian(DATSTRUCT.AllNeurSimBase);
% -------- HIGH CORR
inds = find(DATSTRUCT.AllNeurSimBase > medianval);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('base neural metric as threshold');


% ======================= USING BASELINE SNR AS CRITERION
medianval = nanmedian(DATSTRUCT.AllSNRbaseEnd);
% -------- HIGH CORR
inds = find(DATSTRUCT.AllSNRbaseEnd > medianval);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('base SNR as threshold');


% ========================= HAVE TO HAVE HIGH SNR FOR BOTH BASE AND END
medianval = nanmedian([DATSTRUCT.AllSNRbaseEnd DATSTRUCT.AllSNRtrainEnd]);
% -------- HIGH CORR
inds = find(DATSTRUCT.AllSNRbaseEnd > medianval & DATSTRUCT.AllSNRtrainEnd > medianval);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp, fieldname_baseneur, fieldname_trainneur);
lt_subtitle('base and train SNR as threshold');


%% ============== correlation between FF and neural

lt_figure; hold on;
Yall = {};

% - targ
lt_subplot(1,4,1); hold on;
xlabel('(corr(FFchange, neur change) * learn dir)');
lt_plot_annotation(1, '(more neg is predicted if neural change with ff', 'r');
plotcol = 'k';
inds = ~isnan(DATSTRUCT.AllTargLearnDir) & ~isnan(DATSTRUCT.AllFFneurCorr) ...
    & DATSTRUCT.AllIsTarg==1;

Y = DATSTRUCT.AllTargLearnDir(inds).*DATSTRUCT.AllFFneurCorr(inds);
lt_plot_histogram(Y, '', 1, 1, '', 0, plotcol);
Yall{1} = Y;

% - same
lt_subplot(1,4,2); hold on;
plotcol = 'b';
inds = ~isnan(DATSTRUCT.AllTargLearnDir) & ~isnan(DATSTRUCT.AllFFneurCorr) ...
    & DATSTRUCT.AllIsSame==1;

Y = DATSTRUCT.AllTargLearnDir(inds).*DATSTRUCT.AllFFneurCorr(inds);
lt_plot_histogram(Y, '', 1, 1, '', 0, plotcol);
Yall{2} = Y;

% - diff
lt_subplot(1,4,3); hold on;
plotcol = 'r';
inds = ~isnan(DATSTRUCT.AllTargLearnDir) & ~isnan(DATSTRUCT.AllFFneurCorr) ...
    & DATSTRUCT.AllIsDiff==1;

Y = DATSTRUCT.AllTargLearnDir(inds).*DATSTRUCT.AllFFneurCorr(inds);
lt_plot_histogram(Y, '', 1, 1, '', 0, plotcol);
Yall{3} = Y;

% -- all 3, means
lt_subplot(1,4,4); hold on;
title('all means');
distributionPlot(Yall, 'addSpread', 1, 'showMM', 4);


%% ===============




