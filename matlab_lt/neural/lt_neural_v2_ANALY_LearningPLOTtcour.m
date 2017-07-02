function lt_neural_v2_ANALY_LearningPLOTtcour(MOTIFSTATS_Compiled, birdtoplot, ...
    expttoplot, neuralmetric_toplot, plotOnlyTargSyl)


%%
NumBirds = length(MOTIFSTATS_Compiled.birds);
runbin = 10; % num rends to smth over
plotRawScores = 0; % keep at 0, since too much data if plot
LearningMetadat = lt_neural_v2_LoadLearnMetadat;
ListOfNeuralMetrics = {'NEURrelbase_smthFrCorr', 'NEURrelbase_EuclDistance', ...
    'NEURrelbase_Norm1Distance', 'NEURrelbase_MeanFRDiff', 'NEURrelbase_EuclDistFRcentered', ...
    'NEUR_CVsmthFR', 'NEUR_SDsmthFR', 'NEUR_MeansmthFR'};

%% PLOT ALL EXPERIMENTS, TIMECOURSES
for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    for ii=1:numexpts
        
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        % --- plot this guy?
        if ~isempty(birdtoplot) & ~isempty(expttoplot)
            % then I want a specific expt
            if ~strcmp(birdtoplot, birdname) | ~strcmp(expttoplot, exptname)
                continue
            end
        end
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        assert(length(SummaryStruct.birds)==1, 'asdfasd');
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        figcount=1;
        subplotcols=2;
        subplotrows=3;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        FirstDay = []; % will fill in with first neuron/motif - other enurons use same first day.
        for j=1:nummotifs
            motifname = motiflist{j};
            targetmotifs = MotifStats.params.TargSyls;
            
            if plotOnlyTargSyl==1
                if ~any(strcmp(motifname, targetmotifs))
                    continue
                end
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            ylabel('z-score (color= neural');
            
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
                
                chan = MotifStats.neurons(nn).motif(j).Params.channel_board;
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
                neuralsimilarity = [segmentsextract.(neuralmetric_toplot)];
                
                % ==================== 2) FF
                alltrialFF = [segmentsextract.FF_val];
                
                
                %% ============ CONVERT TO ZSCORE (all metrics)
                %                     if convertToZscore==1
                
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
                %                     end
                
                
                %% PLOT THIS NEUORN
                if length(alltrialFF)>runbin
                    
                    % ===== FF
                    if plotRawScores ==1
                        plot(alltrialTimes, alltrialFF, 'x', 'Color', [0.4 0.4 0.4]);
                    end
                    Yrunning =lt_running_stats(alltrialFF, runbin);
                    xrunning = lt_running_stats(alltrialTimes, runbin);
                    %                         shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', 'k'}, 1);
                    plot(xrunning.Mean, Yrunning.Mean, 'o', 'Color', 'k');
                    
                    
                    % ====== NEURAL SIMILARITY
                    if plotRawScores ==1
                        plot(alltrialTimes, neuralsimilarity, 'o');
                    end
                    Yrunning =lt_running_stats(neuralsimilarity, runbin);
                    xrunning = lt_running_stats(alltrialTimes, runbin);
                    shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', plotcols{nn}}, 1);
                    lt_plot_text(xrunning.Mean(end), Yrunning.Mean(end), ['ch' num2str(chan)], plotcols{nn});
                    % ---
                    lastBaseTrial = alltrialTimes(max(find(baseInds)));
                    line([lastBaseTrial lastBaseTrial], ylim, 'Color', 'k');
                    ylim([-2 2]);
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


%% FOR THE TARGET SYL, COMPARE THE DIFF METRICS THAT COMPARE TO BASELINE
% ONLY PLOT TARGET MOTIFS, AS TOO MUCH DATA
for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    for ii=1:numexpts
        
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        % --- plot this guy?
        if ~isempty(birdtoplot) & ~isempty(expttoplot)
            % then I want a specific expt
            if ~strcmp(birdtoplot, birdname) | ~strcmp(expttoplot, exptname)
                continue
            end
        end
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        assert(length(SummaryStruct.birds)==1, 'asdfasd');
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        figcount=1;
        subplotcols=2;
        subplotrows=3;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        FirstDay = []; % will fill in with first neuron/motif - other enurons use same first day.
        for j=1:nummotifs
            motifname = motiflist{j};
            targetmotifs = MotifStats.params.TargSyls;
            
                if ~any(strcmp(motifname, targetmotifs))
                    continue
                end
            
            for lll = 1:length(ListOfNeuralMetrics)
                neuralmetr_temp = ListOfNeuralMetrics{lll};
                
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                ylabel([neuralmetr_temp '(z)']);
                
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
                    neuralsimilarity = [segmentsextract.(neuralmetr_temp)];
                    
                    % ==================== 2) FF
                    alltrialFF = [segmentsextract.FF_val];
                    
                    
                    %% ============ CONVERT TO ZSCORE (all metrics)
                    %                     if convertToZscore==1
                    
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
                    %                     end
                    
                    
                    %% PLOT THIS NEUORN
                    if length(alltrialFF)>runbin
                        
                        % ===== FF
                        if plotRawScores ==1
                            plot(alltrialTimes, alltrialFF, 'x', 'Color', [0.4 0.4 0.4]);
                        end
                        Yrunning =lt_running_stats(alltrialFF, runbin);
                        xrunning = lt_running_stats(alltrialTimes, runbin);
                        %                         shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', 'k'}, 1);
                        plot(xrunning.Mean, Yrunning.Mean, 'o', 'Color', 'k');
                        
                        
                        % ====== NEURAL SIMILARITY
                        if plotRawScores ==1
                            plot(alltrialTimes, neuralsimilarity, 'o');
                        end
                        Yrunning =lt_running_stats(neuralsimilarity, runbin);
                        xrunning = lt_running_stats(alltrialTimes, runbin);
                        shadedErrorBar(xrunning.Mean, Yrunning.Mean, Yrunning.SEM, {'Color', plotcols{nn}}, 1);
                        
                        % ---
                        lastBaseTrial = alltrialTimes(max(find(baseInds)));
                        line([lastBaseTrial lastBaseTrial], ylim, 'Color', 'k');
                        ylim([-2 2]);
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
        end
        linkaxes(hsplots, 'xy');
        lt_subtitle([birdname '-' exptname])
    end
end


