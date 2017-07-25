function lt_neural_v2_ANALY_Swtch_Summary(MOTIFSTATS_Compiled, SwitchStruct, RemoveLowNumtrials, ...
    MinTrials, UseZscoreNeural, neuralmetricname)
%% TO DO:
% LEARN FOR FOR NONNTARGET ARE STILL NOT IN CORRECT DIR - HAVE TO TAKE
% INTOA ACCOUNT MULTIPLE TARGETS

%%


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
if (0)
    for i=1:Numbirds
        
        numexpts = length(SwitchStruct.bird(i).exptnum);
        birdname = SwitchStruct.bird(i).birdname;
        
        for ii=1:numexpts
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
            
            MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
            SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
            
            motiflist = MotifStats.params.motif_regexpr_str;
            targsyls = MotifStats.params.TargSyls;
            nummotifs = length(motiflist);
            
            WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
                MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
            
            for iii=1:numswitches
                
                figcount=1;
                subplotrows=5;
                subplotcols=6;
                fignums_alreadyused=[];
                hfigs=[];
                
                
                goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] ...
                    & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs]);
                
                if isempty(goodneurons)
                    disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                    continue
                end
                
                swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
                swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
                swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
                
                plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
                
                
                % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
                for j=1:nummotifs
                    
                    for nn=goodneurons
                        
                        segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                        
                        if ~isfield(segextract, 'FRsmooth_xbin')
                            continue
                        end
                        
                        baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                        trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                        
                        if length(baseInds)<minrends | length(trainInds) < minrends
                            % even for good neurosn, could occur if some motifs
                            % labeled pre but not post.
                            continue
                        end
                        
                        trialstoplot = [baseInds trainInds];
                        clustnum = MotifStats.neurons(nn).clustnum;
                        
                        
                        % ================= PLOT TIMECOURSE OF NEURAL AND FF
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        
                        % ----- FF
                        ffvals = [segextract.FF_val];
                        tvals = [segextract.song_datenum];
                        % -- convert tvals to days
                        tmpday = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
                        tvals = lt_convert_EventTimes_to_RelTimes(datestr(tmpday, 'ddmmmyyyy'),...
                            tvals);
                        tvals = tvals.FinalValue;
                        
                        tsmth = lt_running_stats(tvals, FFsmthbinsize);
                        
                        if ~all(isnan(ffvals));
                            % ---- get zscore
                            ffvals_basemean = mean(ffvals(baseInds));
                            ffvalsbaseSD = std(ffvals(baseInds));
                            
                            ffvals = (ffvals - ffvals_basemean)./ffvalsbaseSD;
                            
                            ffsmth = lt_running_stats(ffvals, FFsmthbinsize);
                            
                            plot(tsmth.Median, ffsmth.Median, 'kx');
                            %                       plot(tvals(trialstoplot), ffvals(trialstoplot), 'xk');
                            %
                        end
                        
                        % ---- neural (corr with baseline)
                        neuralsim = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).(neuralmetricname);
                        %                    if ~strcmp(neuralmetricname, 'NEURvsbase_FRcorr')
                        % then get zscore
                        neurbasemean = mean(neuralsim(baseInds));
                        neurbaseSD = std(neuralsim(baseInds));
                        neuralsim = (neuralsim - neurbasemean)./neurbaseSD;
                        %                    end
                        neursmth = lt_running_stats(neuralsim, FFsmthbinsize);
                        %                     lt_plot(tsmth.Mean, neursmth.Mean, {'Errors', neursmth.SEM, ...
                        %                         'Color', plotcols{j}});
                        plot(tsmth.Median, neursmth.Mean, 'o', 'Color', plotcols{nn});
                        %                    plot(tvals(trialstoplot), neuralsim(trialstoplot), 'ob');
                        
                        
                        % --- stuff
                        axis tight;
                        ylim([-3 3]);
                        lt_plot_zeroline;
                        
                        % --- line for base vs. training
                        line([tvals(max(baseInds)) tvals(max(baseInds))], ylim, 'Color','k', 'LineWidth', 2);
                        
                        % --- title
                        if any(strcmp(targsyls, motiflist{j}))
                            title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}], 'Color', 'r');
                        else
                            title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}]);
                        end
                        
                        % ======================== COLLECT DATA FOR PLOTTING
                        
                    end
                end
            end
        end
    end
end

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

% counter =1;
switchnum_global=1;
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
                    
                    if ~isfield(segextract, 'FRsmooth_xbin')
                        continue
                    end
                    
                    baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    
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
                    NumTrials = min([numtrain, ceil(length(trainInds)/2), length(baseInds)]);
                    disp(['numtrials common to base and train (max 25): ' num2str(NumTrials)]);
                    
                    trainEndInds = trainInds(end-NumTrials+1:end);
                    baseEndInds = baseInds(end-NumTrials+1:end);
                    %                     trainEndInds = trainInds;
                    %                     baseEndInds = baseInds;
                    
                    % =================== METRICS
                    ffchange = mean(ffvals(trainEndInds)) - mean(ffvals(baseEndInds));
                    
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
                    
                    x = DATSTRUCT.AllNeurSimBase(indtmp);
                    y = DATSTRUCT.AllNeurSimTrain(indtmp);
                    
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
                        y-x];
                    
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

[~, inds] = unique(BYNEURONDAT.ByNeuron_SwitchCounter);

func = @(X)nanmean(X);
learningall = cellfun(func, BYNEURONDAT.ByNeuron_FFchange(inds,1)).*BYNEURONDAT.ByNeuron_TargLearnDir(inds); % one value for learning  at targ for each switch

threshlearn = nanmedian(learningall);


% ===== HIGH (EXPT WITH GOOD LEARNING)
indstokeep = cellfun(func, BYNEURONDAT.ByNeuron_FFchange(:,1)) >= threshlearn;
if sum(indstokeep)>0
BYNEURONDAT_tmp = lt_structure_subsample_all_fields(BYNEURONDAT, indstokeep, 1);
lt_neural_v2_ANALY_Swtch_Summary_c(BYNEURONDAT_tmp);
end



%% ===== PLOT DISTRIBUTIONS OF NEURAL SIMILARITY (PRE AND POST)

% --- TARG, SAME, DIFF
DATSTRUCT_tmp = DATSTRUCT;
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);

% ---- JUST TARG AND NONTARG
DATSTRUCT_tmp = DATSTRUCT;
lt_neural_v2_ANALY_Swtch_Summary_b(DATSTRUCT_tmp);


%% DISTRIBUTIONS OF NEURAL SIMILARITY [SEPARATE HIGH AND LOW BASELINE PITCH CORR]

% === 1) plot distributions of baseline corrs
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


% ---- USING CORR WITH MEAN FF
% ===== SIGNIFICANT BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
inds = find(DATSTRUCT.AllBaseFFvsNeurFRCorr_p<0.05);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);

% ===== HIGH BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr));
inds = find(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr)>=medianCorr);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);
lt_subtitle('high base ff vs. neurFR corr');

% ===== LOW BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
inds = find(abs(DATSTRUCT.AllBaseFFvsNeurFRCorr)<medianCorr);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);
lt_subtitle('low base ff vs. neurFR corr');



if (0)
    % ===== HIGH BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
    medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeursimCorr));
    inds = find(abs(DATSTRUCT.AllBaseFFvsNeursimCorr)>=medianCorr);
    DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
    lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);
    
    % ===== LOW BASELINE FF VS. NEURAL CORR (I.E. SITES THAT "CARE");
    medianCorr = nanmedian(abs(DATSTRUCT.AllBaseFFvsNeursimCorr));
    inds = find(abs(DATSTRUCT.AllBaseFFvsNeursimCorr)<medianCorr);
    DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
    lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);
end



% ======================= USING BASELINE CORR AS CRITERION
medianval = nanmedian(DATSTRUCT.AllNeurSimBase);
% -------- HIGH CORR
inds = find(DATSTRUCT.AllNeurSimBase > medianval);
DATSTRUCT_tmp = lt_structure_subsample_all_fields(DATSTRUCT, inds);
lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT_tmp);


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




