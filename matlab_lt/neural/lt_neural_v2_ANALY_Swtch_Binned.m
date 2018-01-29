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

%%
% plotneurzscore = 1; then zscore; otherwise plots difference from base
plotFFbyRend=0;
%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);

WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

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

% -- xcorr
AllBinned_SpkMean_vs_FFchange = [];
AllBinned_FRsmChange_vs_FFchange = [];
AllTrial_Spk_vs_FFslope = [];
AllTrial_Spk_vs_FF = [];

% -- corr
AllTrial_Spk_vs_FF_RHO = [];

AllIsTarg = [];
AllBirdnum = [];
AllExptnum = [];
AllSwnum = [];
AllMotifnum = [];
AllNeurnum = [];

AllLearnSlopeScaled = [];
AllLearnSlopeSig = [];

All_SpkMean_STDtrials = [];
All_SpkMean_slopeScaled = [];
%                     All_SpkMean_Change = [];
All_SpkMean_STDtrials_notCV = [];

All_SpkMean_STDtrials_base = [];

All_FR_Modulation = [];
All_StartFromWNOff = [];

All_SingleUnit = [];

for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    if ~isempty(birdname_get)
        if ~strcmp(birdname, birdname_get)
            continue
        end
    end
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        if ~isempty(exptname_get)
            if ~strcmp(exptname, exptname_get)
                continue
            end
        end
        
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        % ====== sanity check
        assert(length(SummaryStruct.birds(1).neurons) == length(MotifStats.neurons), 'asdfd');
        
        
        motiflist = MotifStats.params.motif_regexpr_str;
        targsyls = MotifStats.params.TargSyls;
        nummotifs = length(motiflist);
        
        WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
            MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            if ~isempty(switchnum_get)
                if switchnum_get ~= iii
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
            
            if onlyPlotTargNontarg==1
                motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls)) ...
                    find(ismember(motiflist, MotifStats.params.SameTypeSyls))];
            elseif onlyPlotTargNontarg==2
                motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls))];
            else
                motifstoplot = 1:nummotifs;
            end
            
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
                bregionThis = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct.birds(1).neurons(nn).NOTE_Location;
                if ~isempty(Bregion)
                    if ~any(strcmp(Bregion, bregionThis))
                        continue
                    end
                end
                
                for j=motifstoplot
                    
                    % -- direction of training;;;
                    indtmp = find(strcmp(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies, motiflist{j}));
                    if ~isempty(indtmp)
                        learnconting = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{indtmp+1};
                    else
                        learnconting = [];
                    end
                    
                    
                    %%
                    %                     lt_figure; hold on;
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    if ~isfield(segextract, 'fs')
                        continue
                    end
                    
                    
                    % =========================== FIGURE SOME THINGS OUT
                    % ABOUT BINSIZES
                    % ---- make bins relative to onset of training
                    binsize = 7; % 10 trials;
                    baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    
                    if removeTrainStart==1
                        try
                            trainInds(1:19) = [];
                        catch err
                            continue
                        end
                        
                    end
                    if length(trainInds) <winsize
                        continue
                    end
                    
                    
                    %                     binsize = max([minbinsize floor(length(baseInds)/2)]);
                    %                     binsize = min([binsize maxbinsize]);
                    
                    binOnsets1 = baseInds(end)-9:-binsize:1;
                    binOnsets1 = fliplr(binOnsets1);
                    binOnsets2 = trainInds(1):binsize:(trainInds(end)-binsize+1);
                    
                    binOnsets = [binOnsets1 binOnsets2];
                    binOffsets = binOnsets+binsize-1;
                    binIndsRelOnset = [-length(binOnsets1):1:-1 0:1:length(binOnsets2)-1]; % 0 is first training bin;
                    
                    
                    % ###################################
                    
                    % ===================================== NEURAL (spike)
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
                    
                    % -- subtract baseline
                    Nspks_orig = Nspks;
                    Nspks = Nspks - mean(Nspks(baseInds));
                    
                    
                    
                    
                    
                    % =================================== SMOOTHED FR
                    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
                    frx = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    
                    % -- get premotor window
                    indtmp = frx>=windspk(1) & frx<windspk(2);
                    frmat = frmat(indtmp, :);
                    
                    
                    
                    
                    % ======================================= FF
                    tvals = [segextract.song_datenum];
                    tvals = lt_convert_EventTimes_to_RelTimes(datestr(floor(min(tvals)), 'ddmmmyyyy'), tvals);
                    tvals = tvals.FinalValue;
                    
                    ffvals = [segextract.FF_val];
                    
                    % -- subtract baseline
                    ffvals = ffvals - mean(ffvals(baseInds));
                    % -- flip so positive is direction of learning
                    disp('FFVALS FLIPPED SO LEARNING POSITIVE');
                    ffvals = ffvals*targlearndir;
                    
                    
                    % ################################ STATS PER TRIAL
                    FRsm_modulation = std(frmat, 0, 1);
                    FRsm_modulationCV = std(frmat, 0, 1)./mean(frmat,1);
                    
                    
                    
                    
                    % ############################### CONVERT TO STATS PER BIN
                    nbins = length(binOnsets);
                    assert(length(ffvals) == ...
                        length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds), 'sadfasdf');
                    assert(length(ffvals) == length(Nspks), 'asdfga');
                    
                    SPK_mean = [];
                    SPK_change = [];
                    SPK_std = [];
                    FF_mean = [];
                    FF_std = [];
                    FF_change = [];
                    FRsmooth_corrVsPrev = [];
                    
                    frsmoothprev = [];
                    for tt = 1:nbins
                        on = binOnsets(tt);
                        off = binOffsets(tt);
                        
                        % --- spikes
                        spkmean = median(Nspks(on:off));
                        spkstd = std(Nspks(on:off));
                        
                        if tt>1
                            spkmean_last = median(Nspks(binOnsets(tt-1):binOffsets(tt-1)));
                            spkchange = spkmean - spkmean_last;
                        else
                            spkchange = nan;
                        end
                        
                        SPK_mean = [SPK_mean spkmean];
                        SPK_std = [SPK_std spkstd];
                        SPK_change = [SPK_change spkchange];
                        
                        
                        % --- FF
                        ffmean = median(ffvals(on:off));
                        ffstd = std(ffvals(on:off));
                        
                        FF_mean = [FF_mean ffmean];
                        FF_std = [FF_std ffstd];
                        
                        if tt>1
                            
                            ffchange = ffmean - median(ffvals(binOnsets(tt-1):binOffsets(tt-1)));
                        else
                            ffchange= nan;
                        end
                        
                        FF_change = [FF_change ffchange];
                        
                        % --- correlation of smoothed FR vs. previous
                        % bin
                        frsmooth = mean(frmat(:, on:off),2);
                        
                        if ~isempty(frsmoothprev)
                            % -- calc corr
                            rho = corr(frsmooth, frsmoothprev);
                            frsmoothprev = frsmooth;
                        else
                            frsmoothprev = frsmooth;
                            rho = nan;
                        end
                        FRsmooth_corrVsPrev = [FRsmooth_corrVsPrev rho];
                        
                    end
                    
                    
                    
                    % ========================= PLOT (BINNED DATA)
                    if plotraw==1
                        lt_figure; hold on;
                        
                        lt_subplot(7,1,1); hold on;
                        title('mean FF');
                        plot(binIndsRelOnset, FF_mean, 'o-k');
                        
                        lt_subplot(7,1,2); hold on;
                        title('std FF');
                        plot(binIndsRelOnset, FF_std, 'o-k');
                        
                        lt_subplot(7,1,3); hold on;
                        title('mean SPK');
                        plot(binIndsRelOnset, SPK_mean, 'o-r');
                        
                        lt_subplot(7,1,4); hold on;
                        title('std SPK');
                        plot(binIndsRelOnset, SPK_std, 'o-r');
                        
                        line([-0.5 -0.5], ylim, 'Color', 'b');
                        
                        
                        lt_subplot(7,1, 5); hold on;
                        title('FF change (rel prev)');
                        plot(binIndsRelOnset, FF_change, 'ok');
                        lt_plot_zeroline;
                        
                        lt_subplot(7,1, 6); hold on;
                        title('SPK change (rel prev)');
                        plot(binIndsRelOnset, SPK_change, 'or');
                        lt_plot_zeroline;
                        
                        lt_subplot(7,1, 7); hold on;
                        title('FRsmooth (corr vs. prev bin)');
                        plot(binIndsRelOnset, FRsmooth_corrVsPrev, 'ob');
                        lt_plot_zeroline;
                        
                    end
                    
                    
                    
                    
                    % ================= calculate local slope for all
                    % datapoints
                    % --- ffvals
                    npoints = length(ffvals);
                    FFSlopesAll = [];
                    FFSlopesScaledAll = [];
                    for tt = 1:npoints
                        if tt>floor(winsize/2) & (tt+floor(winsize/2))<=npoints
                            y = ffvals(tt-floor(winsize/2):tt+floor(winsize/2));
                            x = 1:length(y);
                            
                            [b, bint] = regress(y', [ones(size(x')) x']);
                            
                            bscale = b(2)/(bint(2,2)-bint(2,1));
                            
                            %                         fitlm(x, y, 'RobustOpts', 'on')
                            
                        else
                            b = [nan nan];
                            bscale = nan;
                        end
                        
                        FFSlopesAll = [FFSlopesAll b(2)];
                        FFSlopesScaledAll = [FFSlopesScaledAll bscale];
                    end
                    
                    
                    
                    % =============== PLOTS TRIAL BY TRIAL
                    if plotraw==1
                        
                        lt_figure; hold on;
                        
                        % --- 1) actual dat
                        lt_subplot(4,1,1); hold on;
                        plot(ffvals([baseInds trainInds]), 'ok');
                        line([length(baseInds)+0.5 length(baseInds)+0.5], ylim, 'Color', 'r');
                        lt_plot_zeroline;
                        
                        % --- 2)
                        lt_subplot(4,1,2); hold on;
                        title('regression slope (scaled in red)');
                        plot(FFSlopesAll([baseInds trainInds]), 'ok');
                        plot(FFSlopesScaledAll([baseInds trainInds]), 'or');
                        lt_plot_zeroline;
                        
                        % ---
                        lt_subplot(4,1,3); hold on;
                        title('spkcount');
                        plot(Nspks([baseInds trainInds]), 'or');
                    end
                    
                    
                    % ###################################################################################
                    % ======================= XCORR OF BINNED DAT
                    if plotraw==1
                        lt_figure; hold on;
                        
                        lt_subplot(5,1,1); hold on;
                        title('mean spk -- mean ff (binned');
                        [cc, lags] = xcov(SPK_mean(binIndsRelOnset>=0), ...
                            FF_mean(binIndsRelOnset>=0), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        lt_subplot(5,1,2); hold on;
                        title('SPK_mean -- FFchange (binned)');
                        [cc, lags] = xcov(SPK_mean(binIndsRelOnset>=0), ...
                            FF_change(binIndsRelOnset>=0), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        lt_subplot(5,1,3); hold on;
                        title('SPK_std -- FFchange (binned)');
                        [cc, lags] = xcov(SPK_std(binIndsRelOnset>=0), ...
                            FF_change(binIndsRelOnset>=0), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        
                        lt_subplot(5,1,4); hold on;
                        title('SPK_change -- FFchange (binned)');
                        [cc, lags] = xcov(SPK_change(binIndsRelOnset>=0), ...
                            FF_change(binIndsRelOnset>=0), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        
                        lt_subplot(5,1,5); hold on;
                        title('FRsmooth_corrvsprev -- FFchange (binned)');
                        [cc, lags] = xcov(FRsmooth_corrVsPrev(binIndsRelOnset>=0), ...
                            FF_change(binIndsRelOnset>=0), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        
                        
                        % ======================================= XCORR OF RAW DAT
                        lt_figure; hold on;
                        
                        % ------------- 1) SPIKES VS. FFVALS
                        lt_subplot(4,1,1); hold on;
                        title('spks -- ff (by trial)');
                        [cc, lags] = xcov(Nspks(trainInds), ffvals(trainInds), 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                        
                        lt_subplot(4,1,2); hold on;
                        title('spks -- ffslope (by trial)');
                        y = FFSlopesScaledAll(trainInds)';
                        x = Nspks(trainInds);
                        x = x(~isnan(y));
                        y = y(~isnan(y));
                        [cc, lags] = xcov(x, y, 'Unbiased');
                        plot(lags, cc, 'o-k');
                        lt_plot_zeroline;
                    end
                    
                    
                    % ================================= OUTPUT (FOR THIS
                    % NEURON/MOTIF)
                    % ------------ BINNED
                    % 1) spkmean vs. ffchange
                    [cc, lags] = xcov(SPK_mean(binIndsRelOnset>=0), ...
                        FF_change(binIndsRelOnset>=0), ccmaxlagbin, 'Coeff');
                    AllBinned_SpkMean_vs_FFchange = [AllBinned_SpkMean_vs_FFchange; ...
                        cc];
                    
                    % 2) frsmooth_change vs. ff change
                    [cc, lags] = xcov(FRsmooth_corrVsPrev(binIndsRelOnset>=0), ...
                        FF_change(binIndsRelOnset>=0), ccmaxlagbin, 'Coeff');
                    AllBinned_FRsmChange_vs_FFchange = [AllBinned_FRsmChange_vs_FFchange;
                        cc];
                    
                    % ------------ TRIAL BY TRIAL
                    % 1) spk vs. ff slope
                    y = FFSlopesScaledAll(trainInds)';
                    x = Nspks(trainInds);
                    x = x(~isnan(y));
                    y = y(~isnan(y));
                    [cc, lags] = xcov(x, y, ccmaxlagtrial,  'Coeff');
                    
                    AllTrial_Spk_vs_FFslope = [AllTrial_Spk_vs_FFslope; cc'];
                    
                    % 2) spk vs. ff
                    [cc, lags] = xcov(Nspks(trainInds), ffvals(trainInds), ccmaxlagtrial, 'Coeff');
                    AllTrial_Spk_vs_FF = [AllTrial_Spk_vs_FF; cc'];
                    
                    if useLagZeroCorr==1
                        rho = corr(Nspks(trainInds), ffvals(trainInds)');
                    else
                        indtmp = find(lags==0);
                        rho = mean(cc(indtmp-2:indtmp+2));
                    end
                    AllTrial_Spk_vs_FF_RHO = [AllTrial_Spk_vs_FF_RHO; rho];
                    %                     close all;
                    %                     figure; hold on ;
                    %                     plot(lags, cc, 'o-k');
                    %                     plot(0, tmp, 'sr');
                    
                    
                    
                    % -------------- DATA
                    istarg = any(strcmp(motiflist{j}, MotifStats.params.TargSyls));
                    AllIsTarg = [AllIsTarg;istarg];
                    AllBirdnum = [AllBirdnum; i];
                    AllExptnum = [AllExptnum; ii];
                    AllSwnum = [AllSwnum; iii];
                    AllMotifnum = [AllMotifnum; j];
                    AllNeurnum = [AllNeurnum; nn];
                    
                    
                    % ------------------- LEARNING
                    % regression to find slope of learning
                    % -- first remove extreme outliers
                    ffvalstmp = ffvals(trainInds);
                    tvalstmp = tvals(trainInds);
                    Nspkstmp = Nspks(trainInds);
                    
                    if removeOutlier==1
                        indtmp = ffvalstmp>3.5*std(ffvals);
                    else
                        indtmp = [];
                    end
                    ffvalstmp(indtmp) = [];
                    tvalstmp(indtmp) =[];
                    Nspkstmp(indtmp) = [];
                    % --- learning
                    [b,bint]=lt_regress(ffvalstmp, tvalstmp, 0);
                    learnSlopeCI = bint(2,:);
                    learnSig = (learnSlopeCI(1)*learnSlopeCI(2))>0;
                    learnSlopeScaled = b(2)/(bint(2,2) - bint(2,1));
                    
                    AllLearnSlopeScaled = [AllLearnSlopeScaled; learnSlopeScaled];
                    AllLearnSlopeSig = [AllLearnSlopeSig; learnSig];
                    
                    
                    % -- spikes
                    [b,bint]=lt_regress(Nspkstmp, tvalstmp, 0);
                    spkslopeScaled = b(2)/(bint(2,2) - bint(2,1));
                    
                    
                    
                    
                    
                    % --------------------- FR STATS
                    cvspks = std(Nspks_orig(trainInds))/mean(Nspks_orig(trainInds));
                    All_SpkMean_STDtrials = [All_SpkMean_STDtrials; cvspks];
                    All_SpkMean_STDtrials_notCV = [All_SpkMean_STDtrials_notCV; std(Nspks_orig(trainInds))];
                    
                    cvspks_base = std(Nspks_orig(baseInds))/mean(Nspks_orig(baseInds));
                    All_SpkMean_STDtrials_base = [All_SpkMean_STDtrials_base; cvspks_base];
                    All_SpkMean_slopeScaled = [All_SpkMean_slopeScaled; spkslopeScaled];
                    
                    All_FR_Modulation = [All_FR_Modulation; mean(FRsm_modulationCV(trainInds))];
                    
                    
                    
                    % ------------------- OTHER THINGS
                    tmp = cell2mat(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies(2:2:end));
                    startstate = tmp(1:2:end);
                    assert(length(unique(startstate))==1, 'safsdfasdf');
                    
                    if unique(startstate)==0
                        All_StartFromWNOff = [All_StartFromWNOff; 1];
                    else
                        All_StartFromWNOff = [All_StartFromWNOff; 0];
                    end
                    
                    su = SummaryStruct.birds(1).neurons(nn).NOTE_is_single_unit;
                    All_SingleUnit = [All_SingleUnit strcmp('su', 'yes')];
                    
                end
            end
            
            % ==================== PLOT FOR THIS SWITCH
            if plotraw==1
                %                 close all;
                lt_figure; hold on;
                
                % ---- TARG
                lt_subplot(4,2,1); hold on;
                title('targ (all neurons, motifs');
                xlabel('spkmean --- ff change');
                ylabel('Binned');
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1;
                
                ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                plot(x, ccmat, '-k');
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                % ---- NONTARG
                lt_subplot(4,2,2); hold on;
                title('nontargcontext (same syl) (all neurons, motifs');
                xlabel('spkmean --- ff change');
                ylabel('Binned');
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==0;
                
                ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                plot(x, ccmat, '-k');
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                
                % ---- TARG
                lt_subplot(4,2,3); hold on;
                title('targ (all neurons, motifs');
                xlabel('spk --- ff slope');
                ylabel('Trials');
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1;
                
                ccmat = AllTrial_Spk_vs_FFslope(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                plot(x, ccmat, '-k');
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                % ---- NONTARG
                lt_subplot(4,2,4); hold on;
                title('nontargcontext (same syl) (all neurons, motifs');
                xlabel('spk --- ff slope');
                ylabel('Trials');
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==0;
                
                ccmat = AllTrial_Spk_vs_FFslope(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                plot(x, ccmat, '-k');
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
            end
            
        end
    end
end


%% ========================== ONE PAIR OF PLOTS FOR EACH EXPT [COMBINING MOTIFS AND NEURONS]
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            % ============== targ
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==1;
            if ~any(inds)
                continue
            end
            
            birdname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            
            ccmat = AllTrial_Spk_vs_FF(inds,:);
            x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii) '[TARG]']);
            plot(x, ccmat, '-k');
            
            % -- stuff
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            ylim([-0.5 0.5]);
            
            % ============== nontarg
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii) '[NONTARG]']);
            
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0;
            if ~any(inds)
                continue
            end
            
            ccmat = AllTrial_Spk_vs_FF(inds,:);
            x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
            
            plot(x, ccmat, '-k');
            
            % stuff
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            ylim([-0.5 0.5]);
            
            
        end
    end
end

%% ========================== ONE PAIR OF PLOTS FOR EACH EXPT/NEURON
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            for j = 1:maxneuron
                
                % ============== targ
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                if ~any(inds)
                    continue
                end
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(j) '[TARG]']);
                
                birdname = SwitchStruct.bird(i).birdname;
                exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
                
                ccmat = AllTrial_Spk_vs_FF(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                
                plot(x, ccmat, '-k');
                
                % -- display learning rate
                learnRate = AllLearnSlopeScaled(inds);
                learnSig = AllLearnSlopeSig(inds);
                for jj = 1:length(learnRate)
                    lr = learnRate(jj);
                    ls = learnSig(jj);
                    if ls==1
                        lt_plot_text(x(1), 0.3+0.3*jj, ['learnrate(scaled): ' num2str(lr)], 'r');
                    else
                        lt_plot_text(x(1), 0.3+0.3*jj, ['learnrate(scaled): ' num2str(lr)], 'k');
                    end
                end
                % -- stuff
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                ylim([-0.5 0.5]);
                
                % ============== nontarg
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(j) '[NONTARG]']);
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    &  AllNeurnum==j & AllIsTarg==0;
                if ~any(inds)
                    continue
                end
                
                ccmat = AllTrial_Spk_vs_FF(inds,:);
                x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
                
                plot(x, ccmat, '-k');
                % -- display learning rate
                % -- display learning rate
                learnRate = AllLearnSlopeScaled(inds);
                learnSig = AllLearnSlopeSig(inds);
                for jj = 1:length(learnRate)
                    lr = learnRate(jj);
                    ls = learnSig(jj);
                    if ls==1
                        lt_plot_text(x(1), 0.3+0.3*jj, ['learnrate(scaled): ' num2str(lr)], 'r');
                    else
                        lt_plot_text(x(1), 0.3+0.3*jj, ['learnrate(scaled): ' num2str(lr)], 'k');
                    end
                end
                % stuff
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                ylim([-0.5 0.5]);
                
                
            end
        end
    end
end

%% ============= [SPK VS. FF correlation] vs [learning rate]

% ######################################### V1 - does not average over all
% motifs, so can have multiple motifs per neuron

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

X = [];
Y = [];
TargStat = [];
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                
                % skip if neuron does not contribute both targ and nontarg
                if length(unique(AllIsTarg(inds)))==1
                    continue
                end
                
                
                istarg = AllIsTarg(inds);
                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO(inds));
                
                X = [X; learnrate];
                Y = [Y; SpkFFcorr];
                TargStat = [TargStat; istarg];
                
            end
        end
    end
end


% ==================== PLOT
lt_figure; hold on;
hsplots = [];
% =================== TARGS
hsplot = lt_subplot(2,2,1); hold on;
hsplots = [hsplots hsplot];
title('TARG [neurons matched targ vs ntarg]');
xlabel('learn rate');
ylabel('spk vs. ff corr (abs value)');

inds = TargStat==1;
x = X(inds);
y = Y(inds);

lt_regress(y, x, 1, 0, 1, 1, 'r');


% ================= NONTARGS
hsplot = lt_subplot(2,2,2); hold on;
hsplots = [hsplots hsplot];
title('NONTARG');
xlabel('learn rate');
ylabel('spk vs. ff corr (abs value)');

inds = TargStat==0;
x = X(inds);
y = Y(inds);

lt_regress(y, x, 1, 0, 1, 1, 'k');

linkaxes(hsplots)



%% [IMPORTANT] ######################### CHANGE IN NEURAL OVER LEARNING


%% [IMPORTANT] ######################################### V2 - if multiple motifs, then
% takes average --> i.e. each neuron contributes exactly 2 points (targ,
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 1;
onlyIfStartWNoff =0;
overlayIfSU=1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);
lt_figure; hold on;
lt_subplot(3, 1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');

X = [];
Y = [];
Birdnum = [];
TargStat = [];
SwitchCount = [];
swcounter = 1;
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                if onlyIfStartWNoff==1
                    if unique(All_StartFromWNOff(inds))==0
                        continue
                    end
                end
                
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                
                % skip if neuron does not contribute both targ and nontarg
                if length(unique(AllIsTarg(inds)))==1
                    continue
                end
                
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllLearnSlopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    if any(learnsig ==0) | any(learntarg<0)
                        continue
                    end
                    disp(learntarg);
                end
                
                if onlyKeepIfGreaterLearnNontarg==1
                    targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1));
                    nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==0));
                    
                    if targlearn>nontarglearn
                        continue
                    end
                end
                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                
                istarg = AllIsTarg(inds);
                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO(inds));
                        ylabel('spk vs. ff corr (abs value)');
                    case 1
                        % std of mean spks over trials
                        SpkFFcorr = All_SpkMean_STDtrials(inds);
                        ylabel('std of mean spks over trials (CV)');
                    case 2
                        % mean FR modulation (in units of cv) over trials
                        SpkFFcorr = All_FR_Modulation(inds);
                        ylabel('mean FR modulation (in units of cv) over trials');
                    case 3
                        ylabel('corr of spk vs ff slope (trials)');
                        midind = ceil(size(AllTrial_Spk_vs_FFslope,2)/2);
                        SpkFFcorr = abs(mean(AllTrial_Spk_vs_FFslope(inds, midind-2:midind+2),2));
                    case 4
                        ylabel('corr of spk vs ff change (binned)');
                        midind = ceil(size(AllBinned_SpkMean_vs_FFchange,2)/2);
                        SpkFFcorr = abs(mean(AllBinned_SpkMean_vs_FFchange(inds, midind-2:midind+2),2));
                    case 5
                        ylabel('slope of nspk over trials (scaled)');
                        SpkFFcorr = All_SpkMean_slopeScaled(inds);
                    case 6
                        ylabel('std of mean spks base (CV)');
                        SpkFFcorr = All_SpkMean_STDtrials_base(inds);
                    case 7
                        ylabel('change in CV of mean spikes (train - base)');
                        SpkFFcorr = All_SpkMean_STDtrials(inds) - All_SpkMean_STDtrials_base(inds);
                end
                
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                istarg = unique(istarg);
                
                X = [X; learnrate];
                Y = [Y; SpkFFcorr];
                TargStat = [TargStat; istarg];
                
                
                % ======================================= NONTARG
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0;
                
                istarg = AllIsTarg(inds);
                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO(inds));
                    case 1
                        % std of mean spks over trials
                        SpkFFcorr = All_SpkMean_STDtrials(inds);
                    case 2
                        % mean FR modulation (in units of cv) over trials
                        SpkFFcorr = All_FR_Modulation(inds);
                    case 3
                        ylabel('corr of spk vs ff slope (trials)');
                        midind = ceil(size(AllTrial_Spk_vs_FFslope,2)/2);
                        SpkFFcorr = abs(mean(AllTrial_Spk_vs_FFslope(inds, midind-2:midind+2),2));
                    case 4
                        ylabel('corr of spk vs ff change (binned)');
                        midind = ceil(size(AllBinned_SpkMean_vs_FFchange,2)/2);
                        SpkFFcorr = abs(mean(AllBinned_SpkMean_vs_FFchange(inds, midind-2:midind+2),2));
                    case 5
                        ylabel('slope of nspk over trials (scaled)');
                        SpkFFcorr = All_SpkMean_slopeScaled(inds);
                    case 6
                        ylabel('std of mean spks base (CV)');
                        SpkFFcorr = All_SpkMean_STDtrials_base(inds);
                    case 7
                        ylabel('change in CV of mean spikes (train - base)');
                        SpkFFcorr = All_SpkMean_STDtrials(inds) - All_SpkMean_STDtrials_base(inds);
                end
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                istarg = unique(istarg);
                
                X = [X; learnrate];
                Y = [Y; SpkFFcorr];
                TargStat = [TargStat; istarg];
                
                % ===== count switch
                SwitchCount = [SwitchCount swcounter];
                swcounter = swcounter+1;
                
                Birdnum = [Birdnum i];
                
                % --------------- plot line connecting these two \
                if X(end)>X(end-1)
                    plotcol = 'b';
                else
                    plotcol = [0.7 0.7 0.7];
                end
                if Y(end)<Y(end-1)
                    lstyle = '-';
                else
                    lstyle = ':';
                end
                plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
                    'LineStyle', lstyle);
                
                % ---- if plot SU
                if overlayIfSU==1
                    if unique(All_SingleUnit(inds))==1
                        plot(X(end-1), Y(end-1), 'sm', 'MarkerSize', 15);
                        plot(X(end), Y(end), 'sm');
                    end
                end
            end
        end
    end
end


% =================== TARGS
inds = TargStat==1;
x = X(inds);
y = Y(inds);

lt_plot(x, y, {'Color', 'r'})

% ================= NONTARGS
inds = TargStat==0;
x = X(inds);
y = Y(inds);
lt_plot(x, y, {'Color', 'k'})

YLIM = ylim;
% =============================================== PAIRED COMPARISON
lt_subplot(3,2,5); hold on;
xlabel('TARG -- NONTARG');
% ylabel('spk vs. ff corr (abs value)');
x = [1 2];
yy = [];

% - targ
inds = TargStat==1;
y = Y(inds);
yy= [yy y];

% - nontarg
inds = TargStat==0;
y = Y(inds);
yy = [yy y];

plot(x, yy, 'o-k');
xlim([0 3]);
ylim(YLIM);

% --- sign rank
p = signrank(yy(:,1), yy(:,2));
lt_plot_pvalue(p, 'srank', 1);


% =============================================== PAIRED COMPARISON
lt_subplot(3,2,6); hold on;
xlabel('learn diff (targ - nontarg)');
ylabel('diff in metric (targ - nontarg)');

x = X(TargStat==1) - X(TargStat==0);
y = Y(TargStat==1) - Y(TargStat==0);

lt_regress(y, x, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% [IMPORTANT] [SAME AS ABOVE, BUT EACH EXPT ONE PAIR OF DATAPOINTS]
% ######################################### V2 - if multiple motifs, then
% takes average --> i.e. each neuron contributes exactly 2 points (targ,
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);
lt_figure; hold on;
lt_subplot(3,1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');
ylabel('spk vs. ff corr (abs value)');

X = [];
Y = [];
Birdnum = [];
TargStat = [];
SwitchCount = [];
swcounter = 1;
% for i=1:maxbirds
for i=1:6
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii;
            
            % skip if no dat
            if ~any(inds)
                continue
            end
            
            % skip if neuron does not contribute both targ and nontarg
            if length(unique(AllIsTarg(inds)))==1
                continue
            end
            
            % skip if no target has significant learning
            if onlySigLearning==1
                learnsig = AllLearnSlopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1);
                learntarg = AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1);
                
                if any(learnsig ==0) | any(learntarg<0)
                    continue
                end
                disp(learntarg);
            end
            
            if onlyKeepIfGreaterLearnNontarg==1
                targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1));
                nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==0));
                
                if targlearn>nontarglearn
                    continue
                end
            end
            % ======================================= TARGET
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==1;
            
            istarg = AllIsTarg(inds);
            learnrate = AllLearnSlopeScaled(inds);
            SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO(inds));
            
            learnrate = mean(learnrate);
            SpkFFcorr = mean(SpkFFcorr);
            istarg = unique(istarg);
            
            X = [X; learnrate];
            Y = [Y; SpkFFcorr];
            TargStat = [TargStat; istarg];
            
            
            % ======================================= NONTARG
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0;
            
            istarg = AllIsTarg(inds);
            learnrate = AllLearnSlopeScaled(inds);
            SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO(inds));
            
            learnrate = mean(learnrate);
            SpkFFcorr = mean(SpkFFcorr);
            istarg = unique(istarg);
            
            X = [X; learnrate];
            Y = [Y; SpkFFcorr];
            TargStat = [TargStat; istarg];
            
            % ===== count switch
            SwitchCount = [SwitchCount swcounter];
            swcounter = swcounter+1;
            
            Birdnum = [Birdnum i];
            
            % --------------- plot line connecting these two \
            if X(end)>X(end-1)
                plotcol = 'b';
            else
                plotcol = [0.7 0.7 0.7];
            end
            if Y(end)<Y(end-1)
                lstyle = '-';
            else
                lstyle = ':';
            end
            plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
                'LineStyle', lstyle);
        end
    end
end


% =================== TARGS
inds = TargStat==1;
x = X(inds);
y = Y(inds);

lt_plot(x, y, {'Color', 'r'})

% ================= NONTARGS
inds = TargStat==0;
x = X(inds);
y = Y(inds);
lt_plot(x, y, {'Color', 'k'})


% =============================================== PAIRED COMPARISON
lt_subplot(3,1,3); hold on;
xlabel('TARG -- NONTARG');
ylabel('spk vs. ff corr (abs value)');
x = [1 2];
yy = [];

% - targ
inds = TargStat==1;
y = Y(inds);
yy= [yy y];

% - nontarg
inds = TargStat==0;
y = Y(inds);
yy = [yy y];

plot(x, yy, 'o-k');
xlim([0 3]);
ylim([-0.1 0.8]);



%% ################################### ALL DATA (NOT NECESSARILY PAIRED)
lt_figure; hold on;

% ========================= TARGET
lt_subplot(3,2,1); hold on;
title('TARG (all dat)');
xlabel('learning(scaled)');
ylabel('spk vs. learning CORR');
inds = AllIsTarg==1;
x = AllLearnSlopeScaled(inds);
y = abs(AllTrial_Spk_vs_FF_RHO(inds));

plot(x, y, 'or');

% ========================= NONTARG
% lt_subplot(3,2,2); hold on;
% title('NONTARG');
% xlabel('learning(scaled)');
% ylabel('spk vs. learning CORR');
inds = AllIsTarg==0;
x = AllLearnSlopeScaled(inds);
y = abs(AllTrial_Spk_vs_FF_RHO(inds));

plot(x, y, 'ok');



%% =========================== PLOT SUMMARY ACROSS ALL EXPERIMENTS

% ==================== PLOT
lt_figure; hold on;

% ---- TARG
lt_subplot(4,2,1); hold on;
title('targ (all neurons, motifs');
xlabel('spkmean --- ff change');
ylabel('Binned');
inds = AllIsTarg==1;

ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);

plot(x, ccmat, '-k');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---- NONTARG
lt_subplot(4,2,2); hold on;
title('nontargcontext (same syl) (all neurons, motifs');
xlabel('spkmean --- ff change');
ylabel('Binned');
inds = AllIsTarg==0;

ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, ccmat, '-k');
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ---- TARG
lt_subplot(4,2,3); hold on;
title('targ (all neurons, motifs');
xlabel('spk --- ff slope');
ylabel('Trials');
inds = AllIsTarg==1;

ccmat = AllTrial_Spk_vs_FFslope(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, ccmat, '-k');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---- NONTARG
lt_subplot(4,2,4); hold on;
title('nontargcontext (same syl) (all neurons, motifs');
xlabel('spk --- ff slope');
ylabel('Trials');
inds = AllIsTarg==0;

ccmat = AllTrial_Spk_vs_FFslope(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, ccmat, '-k');
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ==================== PLOT (ABSOLUTE VALUE OF CORRELATIONS)
lt_figure; hold on;

hsplots = [];
% ---- TARG
hsplot =lt_subplot(4,2,1); hold on;
hsplots = [hsplots hsplot];
title('targ (all neurons, motifs');
xlabel('spkmean --- ff change');
ylabel('Binned');
inds = AllIsTarg==1;

ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, nanmean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---- NONTARG
hsplot = lt_subplot(4,2,2); hold on;
hsplots = [hsplots hsplot];
title('nontargcontext (same syl) (all neurons, motifs');
xlabel('spkmean --- ff change');
ylabel('Binned');
inds = AllIsTarg==0;

ccmat = AllBinned_SpkMean_vs_FFchange(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, nanmean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

linkaxes(hsplots, 'xy');

hsplots = [];
% ---- TARG
hsplot = lt_subplot(4,2,3); hold on;
hsplots = [hsplots hsplot];
title('targ (all neurons, motifs');
xlabel('spk --- ff slope');
ylabel('Trials');
inds = AllIsTarg==1;

ccmat = AllTrial_Spk_vs_FFslope(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---- NONTARG
hsplot = lt_subplot(4,2,4); hold on;
hsplots = [hsplots hsplot];
title('nontargcontext (same syl) (all neurons, motifs');
xlabel('spk --- ff slope');
ylabel('Trials');
inds = AllIsTarg==0;

ccmat = AllTrial_Spk_vs_FFslope(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

linkaxes(hsplots, 'xy')


% ==================== PLOT (SPIKES VS. FF)
lt_figure; hold on;

hsplots = [];
% ---- TARG
hsplot =lt_subplot(4,2,1); hold on;
hsplots = [hsplots hsplot];
title('targ (all neurons, motifs');
xlabel('spkmean --- ff mean');
ylabel('Trials');
inds = AllIsTarg==1;

ccmat = AllTrial_Spk_vs_FF(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, nanmean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---- NONTARG
hsplot = lt_subplot(4,2,2); hold on;
hsplots = [hsplots hsplot];
title('targ (all neurons, motifs');
xlabel('spkmean --- ff mean');
ylabel('Trials');
inds = AllIsTarg==0;

ccmat = AllTrial_Spk_vs_FF(inds,:);
x = -floor(size(ccmat,2)/2):1:floor(size(ccmat,2)/2);
plot(x, abs(ccmat), '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, nanmean(abs(ccmat)), lt_sem(abs(ccmat)), {'Color', 'r'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

linkaxes(hsplots, 'xy');


%% save Figs?
if saveFigs==1
    saveDir = '/bluejay5/lucas/analyses/neural/FIGS/ANALY_Swtch_Tcourse2';
    
    saveDir = [saveDir '/' birdname_get '_' exptname_get '_sw' num2str(switchnum_get)];
    
    try cd(saveDir)
    catch err
        mkdir(saveDir)
    end
    
    
    lt_save_figs_to_folder(saveDir,0);
    
end

