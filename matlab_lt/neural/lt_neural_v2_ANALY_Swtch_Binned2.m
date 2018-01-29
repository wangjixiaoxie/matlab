function lt_neural_v2_ANALY_Swtch_Binned(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs, onlySingleDir, Bregion)
%% ------------

binbysong = 1; % this first makes datapoints by song, not by rendition. this allows comparing targ and nontarg
removeOutlier =1; % if 1 then removes >3std for slow change in FF and Nspks (slopes)
useLagZeroCorr = 0; % if 0, then takes 5 trial bin centered at 0 of xcorr.
removeTrainStart = 0; % ad hoc, remove N trials from start, i.e. exclude startle


%%

winsize = 29; % for regression (to get slope of learning)

assert(onlyPlotTargNontarg~=0, 'have not coded for nontarg syls yet');
assert(onlySingleDir==1, 'have not coded for when targ multiple directions yet...');

ccmaxlagbin = 5;
ccmaxlagtrial = 25;

if any([isempty(birdname_get) isempty(exptname_get) isempty(switchnum_get)])
    plotraw = 0;
else
    plotraw =1;
end

%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);
% %
% WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

FFsmthbinsize = 10;

minrends = 4; % for both train and base

premotorWind = SwitchStruct.params.premotorWind;


% ------ for bins across renditions
minbinsize = 15; % will take this or numbase rends (divbided by 2), whichever is larger
maxbinsize = 25;
%%

% onlyPlotTargNontarg = 1; % if 1, then only targ syl, if 0, then all syls;
% if 2 then only targ


%%


if mod(FFsmthbinsize,2)==0
    FFsmthbinsize= FFsmthbinsize+1; % conver to odd, so that median tval is at an actual datapoint.
end

%%

% % -- xcorr
% AllBinned_SpkMean_vs_FFchange = [];
% AllBinned_FRsmChange_vs_FFchange = [];
% AllTrial_Spk_vs_FFslope = [];
% AllTrial_Spk_vs_FF = [];
%
% % -- corr
% AllTrial_Spk_vs_FF_RHO = [];

AllIsSame = [];
AllIsTarg = [];
AllBirdnum = [];
AllExptnum = [];
AllSwnum = [];
AllMotifnum = [];
AllNeurnum = [];

AllFFslope = [];
AllSpkSTDsongs = [];
AllFFvsSpk_corr = [];
AllFFslopeSig = [];

% AllLearnSlopeScaled = [];
% AllLearnSlopeSig = [];
%
% All_SpkMean_STDtrials = [];
% All_SpkMean_slopeScaled = [];
% %                     All_SpkMean_Change = [];
% All_SpkMean_STDtrials_notCV = [];
%
% All_SpkMean_STDtrials_base = [];
%
% All_FR_Modulation = [];
% All_StartFromWNOff = [];
% All_SingleUnit = [];



for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        % ====== sanity check
        assert(length(SummaryStruct.birds(1).neurons) == length(MotifStats.neurons), 'asdfd');
        
        
        motiflist = MotifStats.params.motif_regexpr_str;
        targsyls = MotifStats.params.TargSyls;
        samesyls = MotifStats.params.SameTypeSyls;
        nummotifs = length(motiflist);
        
        %         WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
        %             MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            goodneurons = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons;
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
            swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            %             plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            %             if onlyPlotTargNontarg==1
            %                 motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls)) ...
            %                     find(ismember(motiflist, MotifStats.params.SameTypeSyls))];
            %             elseif onlyPlotTargNontarg==2
            %                 motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls))];
            %             else
            %                 motifstoplot = 1:nummotifs;
            %             end
            
            % --- learning at target
            targlearndir = unique(cell2mat({SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}}));
            
            if onlySingleDir==1
                if length(targlearndir)>1
                    continue
                end
            end
            
            % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
            for nn=goodneurons
                
                % ============ brain region filter
                bregionThis = SummaryStruct.birds(1).neurons(nn).NOTE_Location;
                if ~isempty(Bregion)
                    if ~any(strcmp(Bregion, bregionThis))
                        continue
                    end
                end
                
                
                %% make matrices of (song x motif)
                figcount=1;
                subplotrows=5;
                subplotcols=2;
                fignums_alreadyused=[];
                hfigs=[];
                
                
                % ============= 1) get list of songs for each motif
                FnamesAll = cell(1, length(motiflist));
                NumSongs = nan(1, length(motiflist));
                IsTarg = nan(1, length(motiflist));
                IsSame = nan(1, length(motiflist));
                for j=1:length(motiflist)
                    
                    if ~isfield(MotifStats.neurons(nn).motif(j).SegmentsExtract, 'song_filename')
                        continue
                    end
                    
                    fnames = unique([MotifStats.neurons(nn).motif(j).SegmentsExtract.song_datenum]);
                    
                    FnamesAll{j} = fnames;
                    NumSongs(j) = length(fnames);
                    IsTarg(j) = any(strcmp(targsyls, motiflist{j}));
                    IsSame(j) = any(strcmp(samesyls, motiflist{j}));
                end
                
                if rand<0
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    lt_plot_histogram(NumSongs(~isnan(NumSongs)));
                end
                
                % ============ 2) throw out motifs that are rare
                n1 = min(1*NumSongs(IsTarg==1));
                n2 =NumSongs(IsTarg==1 | IsSame==1);
                N = min(n2(n2>=n1));
                %                 N = min(NumSongs(IsTarg==1 | IsSame==1));
                %                 prctile(N, []
                N = round(0.95*N); % any motifs with fewer than this songs, will throw out.
                
                motifsToKeep= find(NumSongs>=N);
                
                commonSongs = FnamesAll{motifsToKeep(1)};
                for j=2:length(motifsToKeep)
                    commonSongs = intersect(commonSongs, FnamesAll{motifsToKeep(j)});
                end
                commonSongs = sort(commonSongs);
                
                % ============= 3) prepare matrices with songs of
                % intersection
                ThisMotifOrigNums = motifsToKeep;
                ThisNSpksMat  = nan(length(commonSongs), length(motifsToKeep));
                ThisFFmat = nan(length(commonSongs), length(motifsToKeep));
                ThisTvalsVec = commonSongs;
                ThisTvalsVec = lt_convert_EventTimes_to_RelTimes(...
                    datestr(floor(min(ThisTvalsVec)), 'ddmmmyyyy'), ThisTvalsVec);
                ThisTvalsVec = ThisTvalsVec.FinalValue;
                ThisIsTargVec = IsTarg(motifsToKeep);
                ThisIsSameVec = IsSame(motifsToKeep);
                
                % ---- what are baseline/training songs?
                ThisbaseIndsSongs = find(commonSongs>=swpre & commonSongs<swthis);
                ThistrainIndsSongs = find(commonSongs>=swthis & commonSongs<swpost);
                
                % ============ 4) Collect data
                for j=1:length(motifsToKeep)
                    motifnum = motifsToKeep(j);
                    
                    segextract = MotifStats.neurons(nn).motif(motifnum).SegmentsExtract;
                    
                    % ============== GET SPIKES MATRIX
                    clustnum = MotifStats.neurons(nn).clustnum;
                    clustCell = {segextract.spk_Clust};
                    spkCell = {segextract.spk_Times};
                    % -- sanity check
                    assert(all(cell2mat(clustCell) == clustnum), 'asdfasd');
                    % --- extract numspikes for each trial (in premotor
                    % window)
                    windspk = MotifStats.params.motif_predur+premotorWind;
                    numtrials = length(spkCell);
                    Nspks = [];
                    for tt = 1:numtrials
                        spkt = spkCell{tt} >= windspk(1) & spkCell{tt}<windspk(2);
                        Nspks = [Nspks; length(spkt)];
                    end
                    
                    % ============= GET FFVALS MATRIX
                    % -- will subtract base and also flip if negative
                    % learnign
                    ffvals = [segextract.FF_val];
                    baseInds = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(motifnum).baseInds;
                    assert(length(baseInds) == length(ffvals), 'asdfasd');
                    ffvals = ffvals - mean(ffvals(baseInds)); % subtract baseline mean.
                    ffvals = ffvals * targlearndir; % positive = learning
                    
                    
                    
                    % --------------------- GO THRU ALL SONGS AND FIND
                    for jj=1:length(commonSongs)
                        songthis = commonSongs(jj);
                        
                        indthis = [segextract.song_datenum] == songthis;
                        assert(any(indthis),'asdf');
                        
                        % ------------- OUTPUT
                        ThisNSpksMat(jj, j) = mean(Nspks(indthis));
                        ThisFFmat(jj, j) = mean(ffvals(indthis));
                        
                    end
                end
                
                %% =============== GET STATS FOR EACH MOTIF
                % -- do for each motif separately
                nmotifs = size(ThisFFmat,2);
                
                % ------------------ SLOPE OF LEARNING
                learnslope = nan(1,nmotifs);
                learnsig = nan(1,nmotifs);
                for j=1:nmotifs
                    ffvals = ThisFFmat(:,j);
                    tvals = ThisTvalsVec;
                    
                    if all(isnan(ffvals))
                        learnslope(j) = nan;
                        learnsig(j) = nan;
                        continue
                    end
                    
                    [b,bint] = ...
                        lt_regress(ffvals(ThistrainIndsSongs), ...
                        tvals(ThistrainIndsSongs), 0);
                    
                    learnslope(j) = b(2)/(bint(2,2) - bint(2,1));
                    learnsig(j) = bint(2,1)*bint(2,2)>0;
                end
                
                
                % ------------------ STD OF SPIKES
                spkstd = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    nspk = ThisNSpksMat(:, j);
                    
                    spkstd(j) = std(nspk(ThistrainIndsSongs))/mean(nspk(ThistrainIndsSongs));
                end
                
                
                % ------------------ SPK vs. FF corr
                corr_FFvsSpk = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    ffvals = ThisFFmat(ThistrainIndsSongs,j);
                    nspk = ThisNSpksMat(ThistrainIndsSongs,j);
                    
                    rho = corr(nspk, ffvals);
                    
                    corr_FFvsSpk(j) = rho;
                end
                
                %% ================ OUTPUT
                AllIsTarg = [AllIsTarg; ThisIsTargVec'];
                AllIsSame = [AllIsSame; ThisIsSameVec'];
                AllBirdnum = [AllBirdnum; i*ones(size(ThisIsTargVec'))];
                AllExptnum = [AllExptnum; ii*ones(size(ThisIsTargVec'))];
                AllSwnum = [AllSwnum;  iii*ones(size(ThisIsTargVec'))];
                AllMotifnum = [AllMotifnum; ThisMotifOrigNums'];
                AllNeurnum = [AllNeurnum; nn*ones(size(ThisIsTargVec'))];
                
                % ----
                AllFFslope = [AllFFslope; learnslope'];
                AllFFslopeSig = [AllFFslopeSig; learnsig'];
                AllSpkSTDsongs = [AllSpkSTDsongs; spkstd'];
                AllFFvsSpk_corr = [AllFFvsSpk_corr; corr_FFvsSpk'];
                
                
                
                %% ======= PLOT FOR THIS NEURON
                
                if strcmp(birdname, birdname_get) & ...
                        strcmp(exptname, exptname_get) & ...
                        switchnum_get == iii
                    
                    
                    figcount=1;
                    subplotrows=6;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    hsplots = [];
                    
                    for j=1:size(ThisFFmat,2)
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(nn) '-mot' num2str(j)]);
                        ylabel('ff');
                        ff = ThisFFmat(:,j);
                        
                        % -- plotcol
                        if ThisIsTargVec(j)==1
                            pcol = 'k';
                        elseif ThisIsSameVec(j)==1
                            pcol = 'b';
                        else
                            pcol = 'r';
                        end
                        
                        if all(isnan(ff))
                            lt_plot_annotation(1, ['skipped motif' num2str(j) '(nan)'], 'r')
                        else
                            plot(ThisTvalsVec, ff, 'o', 'Color', pcol);
                        end
                        line([ThisTvalsVec(ThisbaseIndsSongs(end)) ThisTvalsVec(ThisbaseIndsSongs(end)) ], ...
                            ylim);
                    end
                    
                    axis tight
                    linkaxes(hsplots, 'xy');
                    
                    % -----------------
                    figcount=1;
                    subplotrows=6;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    hsplots = [];
                    
                    for j=1:size(ThisNSpksMat,2)
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title('nspks');
                        y = ThisNSpksMat(:,j);
                        
                        % -- plotcol
                        if ThisIsTargVec(j)==1
                            pcol = 'k';
                        elseif ThisIsSameVec(j)==1
                            pcol = 'b';
                        else
                            pcol = 'r';
                        end
                        
                        plot(ThisTvalsVec, y, 'o', 'Color', pcol);
                        line([ThisTvalsVec(ThisbaseIndsSongs(end)) ThisTvalsVec(ThisbaseIndsSongs(end)) ], ...
                            ylim);
                    end
                    axis tight
                    linkaxes(hsplots, 'xy');
                    
                end
                
                
                
            end
        end
    end
end



%% ========================= PLOT
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 0;
% onlyIfStartWNoff =0;
overlayIfSU=1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

lt_figure; hold on;
lt_subplot(3, 1, 1:2); hold on;
title('TARG (k), sametype(b), diff(r) [each neuron paired]');
xlabel('learn rate');

X = [];
Y = [];
% Birdnum = [];
TargStat = [];
SameStat = [];
% SwitchCount = [];
% swcounter = 1;
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
%                 if onlyIfStartWNoff==1
%                     if unique(All_StartFromWNOff(inds))==0
%                         continue
%                     end
%                 end
%                 
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
%                 
%                 % skip if neuron does not contribute both targ and nontarg
%                 if length(unique(AllIsTarg(inds)))==1
%                     continue
%                 end
%                 
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
%                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
%                     disp(learntarg);
                end
                
%                 if onlyKeepIfGreaterLearnNontarg==1
%                     targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
%                         & AllNeurnum==j & AllIsTarg==1));
%                     nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
%                         & AllNeurnum==j & AllIsTarg==0));
%                     
%                     if targlearn>nontarglearn
%                         continue
%                     end
%                 end


% ----------- TO COLLECT
thisLearn = [];
thisYval = [];
thisTarg = [];
thisSame = [];

                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];

                
                % ======================================= SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;

                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];


                
                
                % ======================================= DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==0;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];

          
                
                
                
                X = [X; thisLearn];
                Y = [Y; thisYval];
                TargStat = [TargStat; thisTarg];
                SameStat = [SameStat; thisSame];
                
                
                
                % #################### count switch
%                 SwitchCount = [SwitchCount swcounter];
%                 swcounter = swcounter+1;
%                 
%                 Birdnum = [Birdnum i];
                
                % --------------- plot line connecting things
%                
%                 if X(end)>X(end-1)
%                     plotcol = 'b';
%                 else
%                     plotcol = [0.7 0.7 0.7];
%                 end
%                 if Y(end)<Y(end-1)
%                     lstyle = '-';
%                 else
%                     lstyle = ':';
%                 end
%                 plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
%                     'LineStyle', lstyle);
%                 

%                 % ---- if plot SU
%                 if overlayIfSU==1
%                     if unique(All_SingleUnit(inds))==1
%                         plot(X(end-1), Y(end-1), 'sm', 'MarkerSize', 15);
%                         plot(X(end), Y(end), 'sm');
%                     end
%                 end
            end
        end
    end
end

% ================================ PLOT 1 (Y VS. LEARNING)
% lt_figure; hold on;
% lt_subplot(3, 1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');

% --------- FIRST: those that have paired data (targ vs. nontarg)
% indstmp = find(~any(isnan(X)'));
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% ----------- SECOND, plot those that are not paired (open circles)
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
%     plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    plot(X(j, 1), Y(j,1), 'ko');
    plot(X(j, 2), Y(j,2), 'bo');
    plot(X(j, 3), Y(j,3), 'ro');
end


YLIM = ylim;


% =============================================== PAIRED COMPARISON
lt_subplot(3,2,5); hold on;
xlabel('TARG -- SAMESYL -- DIFFSYL');
% ylabel('spk vs. ff corr (abs value)');
x = [1 2 3];
yy = [];

% --------- FIRST: those that have paired data (targ vs. nontarg)
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(x, Y(j,:), 'o-k');
%     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
%     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
%     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% --------- SECOND: not paired
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(x+0.1, Y(j,:), 'o', 'Color', [0.7 0.7 0.7]);
%     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
%     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
%     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

xlim([0 4]);

% --- sign rank
p = signrank(Y(:,1), Y(:,2));
lt_plot_text(1.5, 1.1*max(Y(:,1)), ['(1vs2) p=' num2str(p)]);
p = signrank(Y(:,1), Y(:,3));
lt_plot_text(2, 1.2*max(Y(:,1)), ['(1vs3) p=' num2str(p)]);



% =============================================== PAIRED COMPARISON
lt_subplot(3,2,6); hold on;
xlabel('learn diff (targ - nontarg)');
ylabel('diff in metric (targ - nontarg)');

y = Y(:,1) - Y(:,2);
x = X(:,1) - X(:,2);

lt_regress(y, x, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% ========================= PLOT [v2 -- SWITCH AS A DATPOINT]
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 0;
% onlyIfStartWNoff =0;
overlayIfSU=1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

lt_figure; hold on;
lt_subplot(3, 1, 1:2); hold on;
title('TARG (k), sametype(b), diff(r) [each neuron paired]');
xlabel('learn rate');

X = [];
Y = [];
% Birdnum = [];
TargStat = [];
SameStat = [];
% SwitchCount = [];
% swcounter = 1;
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch

            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii;
                
%                 if onlyIfStartWNoff==1
%                     if unique(All_StartFromWNOff(inds))==0
%                         continue
%                     end
%                 end
%                 
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
%                 
%                 % skip if neuron does not contribute both targ and nontarg
%                 if length(unique(AllIsTarg(inds)))==1
%                     continue
%                 end
%                 
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllIsTarg==1);
                    
                    if ~any(learnsig==1 & learntarg>0)
                        continue
                    end
%                     disp(learntarg);
                end
                
%                 if onlyKeepIfGreaterLearnNontarg==1
%                     targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
%                         & AllNeurnum==j & AllIsTarg==1));
%                     nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
%                         & AllNeurnum==j & AllIsTarg==0));
%                     
%                     if targlearn>nontarglearn
%                         continue
%                     end
%                 end


% ----------- TO COLLECT
thisLearn = [];
thisYval = [];
thisTarg = [];
thisSame = [];

                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii & AllIsTarg==1;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];

                
                % ======================================= SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==0 & AllIsSame==1;

                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];


                
                
                % ======================================= DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==0 & AllIsSame==0;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = abs(AllFFvsSpk_corr(inds));
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];

          
                
                
                
                X = [X; thisLearn];
                Y = [Y; thisYval];
                TargStat = [TargStat; thisTarg];
                SameStat = [SameStat; thisSame];
                
                
                
                % #################### count switch
%                 SwitchCount = [SwitchCount swcounter];
%                 swcounter = swcounter+1;
%                 
%                 Birdnum = [Birdnum i];
                
                % --------------- plot line connecting things
%                
%                 if X(end)>X(end-1)
%                     plotcol = 'b';
%                 else
%                     plotcol = [0.7 0.7 0.7];
%                 end
%                 if Y(end)<Y(end-1)
%                     lstyle = '-';
%                 else
%                     lstyle = ':';
%                 end
%                 plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
%                     'LineStyle', lstyle);
%                 

%                 % ---- if plot SU
%                 if overlayIfSU==1
%                     if unique(All_SingleUnit(inds))==1
%                         plot(X(end-1), Y(end-1), 'sm', 'MarkerSize', 15);
%                         plot(X(end), Y(end), 'sm');
%                     end
%                 end

        end
    end
end

% ================================ PLOT 1 (Y VS. LEARNING)
% lt_figure; hold on;
% lt_subplot(3, 1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');

% --------- FIRST: those that have paired data (targ vs. nontarg)
% indstmp = find(~any(isnan(X)'));
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% ----------- SECOND, plot those that are not paired (open circles)
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
%     plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    plot(X(j, 1), Y(j,1), 'ko');
    plot(X(j, 2), Y(j,2), 'bo');
    plot(X(j, 3), Y(j,3), 'ro');
end


YLIM = ylim;


% =============================================== PAIRED COMPARISON
lt_subplot(3,2,5); hold on;
xlabel('TARG -- SAMESYL -- DIFFSYL');
% ylabel('spk vs. ff corr (abs value)');
x = [1 2 3];
yy = [];

% --------- FIRST: those that have paired data (targ vs. nontarg)
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(x, Y(j,:), 'o-k');
%     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
%     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
%     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% --------- SECOND: not paired
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both 

for j=indstmp
    plot(x+0.1, Y(j,:), 'o', 'Color', [0.7 0.7 0.7]);
%     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
%     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
%     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

xlim([0 4]);

% --- sign rank
p = signrank(Y(:,1), Y(:,2));
lt_plot_text(1.5, 1.1*max(Y(:,1)), ['(1vs2) p=' num2str(p)]);
p = signrank(Y(:,1), Y(:,3));
lt_plot_text(2, 1.2*max(Y(:,1)), ['(1vs3) p=' num2str(p)]);



% =============================================== PAIRED COMPARISON
lt_subplot(3,2,6); hold on;
xlabel('learn diff (targ - nontarg)');
ylabel('diff in metric (targ - nontarg)');

y = Y(:,1) - Y(:,2);
x = X(:,1) - X(:,2);

lt_regress(y, x, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
