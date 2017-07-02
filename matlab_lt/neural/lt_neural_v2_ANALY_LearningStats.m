%% lt 5/22/17 - across all birds, expts, neurons, plot learning related statistics (see other code
% for visualization of raw data)
function [MOTIFSTATS_Compiled] = lt_neural_v2_ANALY_LearningStats(MOTIFSTATS_Compiled)
%% LT 6/29/17 - only keeps neurons that have WNdatenum entered
% NOTE: WILL REMOVE NEURON FROM ALL STRUCTS (inlcuding SummaryStruct)
RemoveIfNoWNdatenum=1;

%%

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
premotorWind = [-0.1 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset
binsize_frvec = 0.02; % for calculating FR vector on each trial (num spikes in bin) (for limits uses premotorWind)



% === plotting similairty scores
plotRawScores=0; % if 0, just plots running avg
% convertToZscore = 1; % converts all scores (e.g neural, FF) to zscore relative to baseline
runbin = 10; % num trials to smooth over.

%%
NumBirds = length(MOTIFSTATS_Compiled.birds);


LearningMetadat = lt_neural_v2_LoadLearnMetadat;

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
NeuralFRCorr_Allbase = [];
NeuralFRCorr_Allpost = [];
IsTargAll = [];
NeuralFRCorr_Allpost_lasttertrile = []; % last 33-tile of data during training.
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
        neuronstoremove = [];
        
        for j=1:nummotifs
            motifname = motiflist{j};
            targetmotifs = MotifStats.params.TargSyls;
            
            istarg = any(strcmp(targetmotifs, motifname));
            
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            for nn=1:numneurons
                
                segmentsextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                
                if ~isfield(segmentsextract, 'FRsmooth_xbin')
                    continue % i.e. no trials
                end
                
                if RemoveIfNoWNdatenum==1
                    if isempty(SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr)
                        % REMOVE THIS NEURON
                        neuronstoremove = [neuronstoremove nn];
                        continue
                    end
                end
                
                % what time bins contain premotor window?
                tmp = segmentsextract(1).FRsmooth_xbin;
                premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                    & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                
                
                % what trials are baseline?
                WNonDnum = datenum(SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr, ...
                    'ddmmmyyyy-HHMM');
                if isempty(WNonDnum)
                    baseInds = nan;
                else
                    baseInds = [segmentsextract.song_datenum] < WNonDnum;
                end
                
                
                
                % ================= get training end inds (i.e. dur
                % training, before ifirst switch that comes after training
                % onset
                if isempty(WNonDnum)
                    trainInds = [segmentsextract.song_datenum] >0;
                else
                    trainInds = [segmentsextract.song_datenum] >= WNonDnum;
                end
                
                trainOn = find(trainInds==1, 1, 'first');
                
                % ------ some expts have multiple switches - for now just
                % keep the first switch!!!
                indbird = strcmp({LearningMetadat.bird.birdname}, birdname);
                indcol = strcmp(LearningMetadat.bird(indbird).info(1,:), exptname);
                alltranstimes = LearningMetadat.bird(indbird).info(3:end, indcol);
                alltranstimes = alltranstimes(~cellfun('isempty', alltranstimes));
                alltransdnums = [];
                for ddd=1:length(alltranstimes)
                    alltransdnums = [alltransdnums ...
                        datenum(alltranstimes{ddd}(1:14), 'ddmmmyyyy-HHMM')]; % switch time
                end
                alltransdnums = unique(alltransdnums);
                
                assert(sum(alltransdnums==WNonDnum)==1, 'asasdfasd');
                
                
                trainenddnum = alltransdnums(find(alltransdnums>WNonDnum, 1, 'first')); % first trans after WN on
                if isempty(trainenddnum)
                    trainOff = find(trainInds==1, 1, 'last');
                else
                    trainOff = find(trainInds==1 & ...
                        [segmentsextract.song_datenum] <= trainenddnum, 1, 'last');
                end
                
                if ~isempty(trainenddnum)
                    if ~all([segmentsextract.song_datenum] <= trainenddnum)
                        disp(['num train days reduced: ' num2str(find(trainInds==1, 1, 'last')) ' to ' ...
                            num2str(find(trainInds==1 & ...
                            [segmentsextract.song_datenum] <= trainenddnum, 1, 'last'))])
                    end
                end
                trainInds(trainOff+1:end) = 0;
                
                
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
                
                
                %% ==== metrics relative to own baseline
                if ~isnan(baseInds)
                    % ============== 1) FIRING RATE (NEURAL SIMILARITY TO BASELINE MEAN
                    % FR)
                    alltrialFR = [segmentsextract.FRsmooth_rate_CommonTrialDur]; % all trial FR (bin x trial)
                    alltrialFR = alltrialFR(premotorInds, :); % extract just premotor window
                    
                    baseFR = alltrialFR(:, baseInds);
                    baseFR_mean = mean(baseFR, 2);
                    
                    % ---- for each trial, calculate deviation from mean base
                    % FR
                    neuralsimilarity = corr(alltrialFR, baseFR_mean); % trials x 1
                    
                    
                    
                    % ============= 2) EUCLIDIAN DISTANCE BETWEEN FR VECTOR
                    binstart = MotifStats.params.motif_predur+premotorWind(1);
                    binend = MotifStats.params.motif_predur+premotorWind(2);
                    binedges = binstart:binsize_frvec:binend;
                    
                    alltrialFRvec = []; % bin x trial
                    for k=1:length(segmentsextract)
                        inds = segmentsextract(k).spk_Clust == MotifStats.neurons(nn).clustnum;
                        spktimes = segmentsextract(k).spk_Times(inds);
                        
                        bincounts = histc(spktimes, binedges);
                        bincounts = bincounts(1:end-1);
                        alltrialFRvec = [alltrialFRvec bincounts'];
                    end
                    baseFRvec = mean(alltrialFRvec(:, baseInds),2);
                    baseFRvec_rep = repmat(baseFRvec, 1, size(alltrialFRvec,2));
                    EuclDistances = sqrt(sum((alltrialFRvec - baseFRvec_rep).^2,1));
                    
                    
                    % ==== 2.5 euclidian after subtracting mean FR (i.e. does
                    % not penalize mean FR diff, but still does penalize
                    % greater depth of modulation, etc
                    alltrialFRvec_minusmean = ...
                        alltrialFRvec - repmat(mean(alltrialFRvec, 1), size(alltrialFRvec,1), 1);
                    baseFRvec_minusmean = baseFRvec - mean(baseFRvec);
                    baseFRvec_minusmeanrep = repmat(baseFRvec_minusmean, 1, size(alltrialFRvec,2));
                    EuclDistancesFRcentered = sqrt(sum((alltrialFRvec_minusmean ...
                        - baseFRvec_minusmeanrep).^2,1));
                    
                    if (0) % troubleshooting, sanity check, plot FR vectors
                        lt_figure; hold on;
                        plot(baseFRvec_minusmean, 'ok-');
                        indstmp = randi(size(alltrialFRvec_minusmean,2), 1, 15);
                        plot(alltrialFRvec_minusmean(:,indstmp), '-');
                        
                        lt_figure; hold on;
                        plot(baseFRvec, 'ok-');
                        indstmp = randi(size(alltrialFRvec,2), 1, 15);
                        plot(alltrialFRvec(:,indstmp), '-');
                    end
                    
                    % ============ 3) MEAN FR DIFFERENCES
                    Norm1Distance =  sum(abs(alltrialFRvec - baseFRvec_rep), 1);
                    
                    % ============ 4) MEAN FR DIFF
                    MeanFRDiff = sum(alltrialFRvec,1) -  sum(baseFRvec,1);
                else
                    neuralsimilarity = nan;
                    EuclDistances = nan;
                    EuclDistancesFRcentered = nan;
                    Norm1Distance = nan;
                    MeanFRDiff = nan;
                end
                
                
                % ==================== 2) FF
                alltrialFF = [segmentsextract.FF_val];
                
                
                % ========================== output
                tmp = num2cell(neuralsimilarity);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEURrelbase_smthFrCorr] = deal(tmp{:});
                
                tmp = num2cell(EuclDistances);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEURrelbase_EuclDistance] = deal(tmp{:});
                
                tmp = num2cell(Norm1Distance);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEURrelbase_Norm1Distance] = deal(tmp{:});
                
                tmp = num2cell(MeanFRDiff);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEURrelbase_MeanFRDiff] = deal(tmp{:});
                
                tmp = num2cell(EuclDistancesFRcentered);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEURrelbase_EuclDistFRcentered] = deal(tmp{:});
                
                
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).Params.baseInds = baseInds;
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).Params.trainInds = trainInds;
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.preTranContingency = preTranContingency;
                MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).LearningInfo.postTranContingency = postTranContingency;
                
                
                NeuralFRCorr_Allbase = [NeuralFRCorr_Allbase median(neuralsimilarity(baseInds))];
                NeuralFRCorr_Allpost = [NeuralFRCorr_Allpost median(neuralsimilarity(trainInds))];
                
                indtmp = find(trainInds);
                NeuralFRCorr_Allpost_lasttertrile = [NeuralFRCorr_Allpost_lasttertrile ...
                    median(neuralsimilarity(indtmp(2*ceil(length(indtmp)/3)+1:end)))];
                
                IsTargAll = [IsTargAll istarg];
                
                %% ======= metrics agnostic to baseline (e.g. modulation ..., burstiness)
                
                % === for modulation, calculate SD of smoothed FR.
                neural_SD = std(alltrialFR, 0, 1); % bin x trial
                neural_mean = mean(alltrialFR, 1);
                neural_CV = neural_SD./neural_mean;
                
                % ========================== output
                tmp = num2cell(neural_CV);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEUR_CVsmthFR] = deal(tmp{:});
                
                tmp = num2cell(neural_SD);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEUR_SDsmthFR] = deal(tmp{:});
                
                tmp = num2cell(neural_mean);
                [MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract.NEUR_MeansmthFR] = deal(tmp{:});
                
                % ===
            end
        end
        % === remove neurons that were mission WN on time
        if ~isempty(neuronstoremove)
            neuronstoremove = unique(neuronstoremove);
            MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(neuronstoremove) = [];
            MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct.birds(1).neurons(neuronstoremove) = [];
            disp(['=========== NOTE, removed neuron (becaseu no wn on time): ' ...
                birdname '-' exptname '-' neuronstoremove]);
        end
        
    end
end

%% ======== PLOT DISTRIBUTIONS [entire base vs. train]
lt_figure; hold on;
% ===== PAIRED
lt_subplot(1,3,1); hold on;
title('corr of smth FR to base mean');
xlabel('baseline'); ylabel('dur train');
lt_plot_45degScatter(NeuralFRCorr_Allbase, NeuralFRCorr_Allpost, 'k', 1);

% targ only
lt_plot_45degScatter(NeuralFRCorr_Allbase(IsTargAll==1), ...
    NeuralFRCorr_Allpost(IsTargAll==1), 'r', 1,0);

% ==== DISTRIBUTIONS
lt_subplot(1,3,2); hold on;
title('baseline');
xlabel('corr')
xbins = min([NeuralFRCorr_Allbase NeuralFRCorr_Allpost])-0.05:0.05:max([NeuralFRCorr_Allbase NeuralFRCorr_Allpost])+0.05;
lt_plot_histogram(NeuralFRCorr_Allbase(IsTargAll==0 & ~isnan(NeuralFRCorr_Allbase)), xbins, 1, 1, '', 1, 'k'); % nontarg
lt_plot_histogram(NeuralFRCorr_Allbase(IsTargAll==1 & ~isnan(NeuralFRCorr_Allbase)), xbins, 1, 1, '', 1, 'r'); % targ
xlim([-0.6 1]);

lt_subplot(1,3,3); hold on;
title('durWN');
xlabel('corr')
lt_plot_histogram(NeuralFRCorr_Allpost(IsTargAll==0 & ~isnan(NeuralFRCorr_Allpost)), xbins, 1, 1, '', 1, 'k'); % nontarg
lt_plot_histogram(NeuralFRCorr_Allpost(IsTargAll==1 & ~isnan(NeuralFRCorr_Allpost)), xbins, 1, 1, '', 1, 'r'); % targ
xlim([-0.6 1]);

lt_plot_annotation(1, 'red = targ; bk = rest', 'r')

%% ======== PLOT DISTRIBUTIONS [same as above, but just end of train]

lt_figure; hold on;
% ===== PAIRED
lt_subplot(1,3,1); hold on;
title('corr of smth FR to base mean');
xlabel('baseline'); ylabel('dur train (last tertile of training');
lt_plot_45degScatter(NeuralFRCorr_Allbase, NeuralFRCorr_Allpost_lasttertrile, 'k', 1);

% targ only
lt_plot_45degScatter(NeuralFRCorr_Allbase(IsTargAll==1), ...
    NeuralFRCorr_Allpost_lasttertrile(IsTargAll==1), 'r', 1,0);

% ==== DISTRIBUTIONS
lt_subplot(1,3,2); hold on;
title('baseline');
xlabel('corr')
xbins = min([NeuralFRCorr_Allbase NeuralFRCorr_Allpost_lasttertrile])-0.05:0.05:max([NeuralFRCorr_Allbase NeuralFRCorr_Allpost_lasttertrile])+0.05;
lt_plot_histogram(NeuralFRCorr_Allbase(IsTargAll==0 & ~isnan(NeuralFRCorr_Allbase)), xbins, 1, 1, '', 1, 'k'); % nontarg
lt_plot_histogram(NeuralFRCorr_Allbase(IsTargAll==1 & ~isnan(NeuralFRCorr_Allbase)), xbins, 1, 1, '', 1, 'r'); % targ
xlim([-0.6 1]);

lt_subplot(1,3,3); hold on;
title('durWN');
xlabel('corr')
lt_plot_histogram(NeuralFRCorr_Allpost_lasttertrile(IsTargAll==0 & ~isnan(NeuralFRCorr_Allpost_lasttertrile)), xbins, 1, 1, '', 1, 'k'); % nontarg
lt_plot_histogram(NeuralFRCorr_Allpost_lasttertrile(IsTargAll==1 & ~isnan(NeuralFRCorr_Allpost_lasttertrile)), xbins, 1, 1, '', 1, 'r'); % targ
xlim([-0.6 1]);

lt_plot_annotation(1, 'red = targ; bk = rest', 'r')


