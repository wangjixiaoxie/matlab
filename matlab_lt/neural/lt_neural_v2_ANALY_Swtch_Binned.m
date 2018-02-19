function lt_neural_v2_ANALY_Swtch_Binned(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs, onlySingleDir, Bregion)
%% ------------

binbysong = 1; % this first makes datapoints by song, not by rendition. this allows comparing targ and nontarg
removeOutlier =1; % if 1 then removes >3std for slow change in FF and Nspks (slopes)
useLagZeroCorr = 0; % if 0, then takes 5 trial bin centered at 0 of xcorr.
removeTrainStart = 0; % ad hoc, remove N trials from start, i.e. exclude startle

scaleSlopesByCI = 1;
convertSDtoCV = 0;

%%

winsize = 19; % for regression (to get slope of learning)

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
AllBinned_SpkMean_vs_FFchange = {};
AllBinned_FRsmChange_vs_FFchange = {};
AllTrial_Spk_vs_FFslope = {};
AllTrial_Spk_vs_FF = {};

AllTrialXcorr_NspkStd_vs_FFmedian = {};
AllTrialXcorr_NspkMedian_vs_FFmedian = {};

% -- corr
AllTrial_Spk_vs_FF_RHO = [];
AllTrial_Spk_vs_FF_RHO_Residuals = [];

AllIsTarg = [];
AllBirdnum = [];
AllExptnum = [];
AllSwnum = [];
AllMotifnum = [];
AllNeurnum = [];
AllIsSameSyl = [];

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

% ======== save raw vectors
AllTrial_raw_ffvals = {};
AllTrial_raw_tvals = {};
AllTrial_raw_lastBaseInd = {};
AllTrial_raw_Nspks = {};

% =========== diff (end minus start)
AllDiff_ff = [];
AllDiff_ffsig = [];
AllDiff_nspk = [];
AllDiff_nspksig = [];

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
        samesyls = MotifStats.params.SameTypeSyls;
        
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
            elseif onlyPlotTargNontarg==3
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
                    binsize = 8; % 10 trials;
                    baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    
                    if isempty(baseInds) | isempty(trainInds)
                        continue
                    end
                    
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
                        if ~any(isnan(ffvals))
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
                    
                    
                    %% ================================= OUTPUT (FOR THIS
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
                    
                    
                    % ------------------- MEDIAN (END MINUS START)
                    trainIndsE = trainInds(1:floor(end/2));
                    trainIndsL = trainInds(floor(end/2)+1:end);
                    
                    ffdiff = median(ffvals(trainIndsL)) - median(ffvals(trainIndsE));
                    spkdiff = median(Nspks(trainIndsL)) - median(Nspks(trainIndsE));
                    % significance
                    if isnan(ffdiff)
                        ffdiff_sig = nan;
                    else
                    ffdiff_sig = ranksum(ffvals(trainIndsL), ffvals(trainIndsE));
                    end
                    spkdiff_sig = ranksum(Nspks(trainIndsL), Nspks(trainIndsE));
                    
                    AllDiff_ff = [AllDiff_ff; ffdiff];
                    AllDiff_ffsig = [AllDiff_ffsig; ffdiff_sig];
                    AllDiff_nspk = [AllDiff_nspk; spkdiff];
                    AllDiff_nspksig = [AllDiff_nspksig; spkdiff_sig];

                    
                    
                    % ------------ TRIAL BY TRIAL
                    % 1) spk vs. ff slope
                    y = FFSlopesScaledAll(trainInds)';
                    x = Nspks(trainInds);
                    x = x(~isnan(y));
                    y = y(~isnan(y));
                    [cc, lags] = xcov(x, y, ccmaxlagtrial,  'Coeff');
                    
                    if isempty(cc)
                        cc = nan(1, size(cc,1));
                    end
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
                    
                    
                    % 3) spk vs. ff (residuals)
                    % to check whether this is affected by slow-timescale
                    % correlations.
                    if any(isnan(ffvals))
                       rho = nan; 
                    else
                    [~, ~, r_ff]=lt_regress(ffvals(trainInds), tvals(trainInds), 0);
                    [~, ~, r_spk]=lt_regress(Nspks(trainInds), tvals(trainInds), 0);
                    
                    rho = corr(r_ff, r_spk);
                    end
                    AllTrial_Spk_vs_FF_RHO_Residuals = ...
                        [AllTrial_Spk_vs_FF_RHO_Residuals; rho];

                    
                    
                    % 4) xcorr of local variability (nspks) vs. learning
                    binrunning = 10;
                    tmp_ff = lt_running_stats(ffvals(trainInds), binrunning);
                    tmp_nspk = lt_running_stats(Nspks(trainInds), binrunning);
                    [cc, lags] = xcov(tmp_nspk.STD, tmp_ff.Median, ccmaxlagtrial, 'Coeff');
                                        
                    AllTrialXcorr_NspkStd_vs_FFmedian = [AllTrialXcorr_NspkStd_vs_FFmedian; cc];
                    
                    [cc, lags] = xcov(tmp_nspk.Median, tmp_ff.Median, ccmaxlagtrial, 'Coeff');
                    AllTrialXcorr_NspkMedian_vs_FFmedian = [AllTrialXcorr_NspkMedian_vs_FFmedian; ...
                        cc];
                    
                    % ------------------- SAVE RAW DATA 
                    AllTrial_raw_ffvals = [AllTrial_raw_ffvals; ...
                        ffvals([baseInds trainInds])];
                    AllTrial_raw_tvals = [AllTrial_raw_tvals; ...
                        tvals([baseInds trainInds])];
                    AllTrial_raw_lastBaseInd = [AllTrial_raw_lastBaseInd; ...
                        length(baseInds)];
                    AllTrial_raw_Nspks = [AllTrial_raw_Nspks; ...
                        Nspks([baseInds trainInds])];
                    
                    
                    
                    % -------------- DATA
                    istarg = any(strcmp(motiflist{j}, MotifStats.params.TargSyls));
                    issame = any(strcmp(motiflist{j}, MotifStats.params.SameTypeSyls));
                    AllIsSameSyl = [AllIsSameSyl; issame];
                    AllIsTarg = [AllIsTarg; istarg];
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
                    if any(isnan(ffvalstmp))
                        learnSlopeScaled = nan;
                        learnSig = nan;
                    else
                        [b,bint]=lt_regress(ffvalstmp, tvalstmp, 0);
                        learnSlopeCI = bint(2,:);
                        learnSig = (learnSlopeCI(1)*learnSlopeCI(2))>0;
                        if scaleSlopesByCI==1
                            learnSlopeScaled = b(2)/(bint(2,2) - bint(2,1));
                        else
                            learnSlopeScaled = b(2);
                        end
                    end
                    
                    AllLearnSlopeScaled = [AllLearnSlopeScaled; learnSlopeScaled];
                    AllLearnSlopeSig = [AllLearnSlopeSig; learnSig];
                    
                    
                    % -- spikes
                    [b,bint]=lt_regress(Nspkstmp, tvalstmp, 0);
                    if scaleSlopesByCI ==1
                    spkslopeScaled = b(2)/(bint(2,2) - bint(2,1));
                    else
                        spkslopeScaled = b(2);
                    end
                    All_SpkMean_slopeScaled = [All_SpkMean_slopeScaled; spkslopeScaled];
                    
                    
                    
                    
                    
                    % --------------------- FR STATS
                    if convertSDtoCV ==1
                        cvspks = std(Nspks_orig(trainInds))/mean(Nspks_orig(trainInds));
                    else
                        cvspks = std(Nspks_orig(trainInds));
                    end
                    All_SpkMean_STDtrials = [All_SpkMean_STDtrials; cvspks];
                    All_SpkMean_STDtrials_notCV = [All_SpkMean_STDtrials_notCV; std(Nspks_orig(trainInds))];
                    
                    if convertSDtoCV==1
                    cvspks_base = std(Nspks_orig(baseInds))/mean(Nspks_orig(baseInds));
                    else
                        cvspks_base = std(Nspks_orig(baseInds));
                    end
                    All_SpkMean_STDtrials_base = [All_SpkMean_STDtrials_base; cvspks_base];
                    
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


for j=1:length(AllTrial_Spk_vs_FFslope)
    if size(AllTrial_Spk_vs_FFslope{j},1)>1
        AllTrial_Spk_vs_FFslope{j} = AllTrial_Spk_vs_FFslope{j}';
    end
end

%% [IMPORTANT] ========================= SUMMARY - FOR EACH EXPT PLOT TIMECOURSES
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);

onlyPlotSigLearn = 1;

for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            %----------------- decide whether this has data
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==1;
            plotcol = 'k';
            if ~any(inds)
                continue
            end
            
            if onlyPlotSigLearn==1
            if ~any(AllLearnSlopeSig(inds))
                continue
            end
            end
            
            
            % -------- initiate figure
            figcount=1;
            subplotrows=6;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
            birdname = MOTIFSTATS_Compiled.birds(i).birdname;
            
            
            % --------- collect things across targs/nontargs
            RhoAll = [];
            SpksdAll = [];
            SylType = []; % 1 = targ ... 2 = samesyl ... 3 = diffsyl
            Xcorr_NspkFF = [];
            Xcorr_NspkSTD_vs_FFmedian = [];
            
            
            % ########################## targ
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==1 & ~isnan(AllTrial_Spk_vs_FF_RHO);
            
            tvals = AllTrial_raw_tvals(inds);
            ffvals = AllTrial_raw_ffvals(inds);
            lastBase = AllTrial_raw_lastBaseInd(inds);
            nspks = AllTrial_raw_Nspks(inds);
            neurnum = AllNeurnum(inds);
            learnrate = AllLearnSlopeScaled(inds);
            learnsig = AllLearnSlopeSig(inds);
            
            
            
            for j=1:length(tvals)
                
                % ==================== FF
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(neurnum(j))]);
                ylabel('FF');
                
                plot(tvals{j}, ffvals{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                ylim([-100 100]);
                lt_plot_zeroline;
                if learnsig(j)==1
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'r');
                else
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'k');
                end
                
                
                % ==================== SPIKES
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii)  '-n' num2str(neurnum(j))]);
                ylabel('spikes');
                
                plot(tvals{j}, nspks{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                lt_plot_zeroline;
                ylim([-50 50]);
                
            end
            
            % ======================= STATS
            rho = AllTrial_Spk_vs_FF_RHO(inds);
            spkSD = All_SpkMean_STDtrials(inds);
            ccmat = cell2mat(AllTrial_Spk_vs_FF(inds));
            indmid = ceil(size(ccmat,2)/2);
            
            RhoAll = [RhoAll; rho];
            SpksdAll = [SpksdAll; spkSD];
            Xcorr_NspkFF = [Xcorr_NspkFF; ccmat];
            
            ccmat = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
            Xcorr_NspkSTD_vs_FFmedian = [Xcorr_NspkSTD_vs_FFmedian; ccmat];
            
            SylType = [SylType; 1*ones(size(spkSD))]; % 1 = targ ... 2 = samesyl ... 3 = diffsyl
            

            % ######################### same syl
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0 & AllIsSameSyl==1 & ~isnan(AllTrial_Spk_vs_FF_RHO);
            plotcol = 'b';
            
            tvals = AllTrial_raw_tvals(inds);
            ffvals = AllTrial_raw_ffvals(inds);
            lastBase = AllTrial_raw_lastBaseInd(inds);
            nspks = AllTrial_raw_Nspks(inds);
            neurnum = AllNeurnum(inds);
            learnrate = AllLearnSlopeScaled(inds);
            learnsig = AllLearnSlopeSig(inds);
            
            for j=1:length(tvals)
                
                % ==================== FF
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(neurnum(j))]);
                ylabel('FF');
                
                plot(tvals{j}, ffvals{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                lt_plot_zeroline;
ylim([-100 100]);                
                if learnsig(j)==1
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'r');
                else
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'k');
                end
                
                
                % ==================== SPIKES
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(neurnum(j))]);
                ylabel('spikes');
                
                plot(tvals{j}, nspks{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                lt_plot_zeroline;
                ylim([-50 50]);
            end
            
            % ======================= STATS
            rho = AllTrial_Spk_vs_FF_RHO(inds);
            spkSD = All_SpkMean_STDtrials(inds);
            ccmat = cell2mat(AllTrial_Spk_vs_FF(inds));
            indmid = ceil(size(ccmat,2)/2);
            
            RhoAll = [RhoAll; rho];
            SpksdAll = [SpksdAll; spkSD];
            Xcorr_NspkFF = [Xcorr_NspkFF; ccmat];
            SylType = [SylType; 2*ones(size(spkSD))]; % 1 = targ ... 2 = samesyl ... 3 = diffsyl
            
                        ccmat = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
            Xcorr_NspkSTD_vs_FFmedian = [Xcorr_NspkSTD_vs_FFmedian; ccmat];

                        
            % ######################### diff syl
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0 & AllIsSameSyl==0 & ~isnan(AllTrial_Spk_vs_FF_RHO);
            plotcol = 'r';
            
            tvals = AllTrial_raw_tvals(inds);
            ffvals = AllTrial_raw_ffvals(inds);
            lastBase = AllTrial_raw_lastBaseInd(inds);
            nspks = AllTrial_raw_Nspks(inds);
            neurnum = AllNeurnum(inds);
            learnrate = AllLearnSlopeScaled(inds);
            learnsig = AllLearnSlopeSig(inds);
            
            for j=1:length(tvals)
                
                % ==================== FF
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(neurnum(j))]);
                ylabel('FF');
                
                plot(tvals{j}, ffvals{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                lt_plot_zeroline;
ylim([-100 100]);

                if learnsig(j)==1
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'r');
                else
                    lt_plot_text(tvals{j}(1), ffvals{j}(end), ...
                        ['learnrate:' num2str(learnrate(j))], 'k');
                end
                
                
                % ==================== SPIKES
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(neurnum(j))]);
                ylabel('spikes');
                
                plot(tvals{j}, nspks{j}, 'x', 'Color', plotcol);
                line([tvals{j}(lastBase{j}) tvals{j}(lastBase{j})], ylim);
                axis tight;
                lt_plot_zeroline;
                ylim([-50 50]);
            end
            
            % ======================= STATS
            rho = AllTrial_Spk_vs_FF_RHO(inds);
            spkSD = All_SpkMean_STDtrials(inds);
            ccmat = cell2mat(AllTrial_Spk_vs_FF(inds));
            indmid = ceil(size(ccmat,2)/2);
            
            RhoAll = [RhoAll; rho];
            SpksdAll = [SpksdAll; spkSD];
            Xcorr_NspkFF = [Xcorr_NspkFF; ccmat];
            SylType = [SylType; 3*ones(size(spkSD))]; % 1 = targ ... 2 = samesyl ... 3 = diffsyl
                        ccmat = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
            Xcorr_NspkSTD_vs_FFmedian = [Xcorr_NspkSTD_vs_FFmedian; ccmat];

            
            % #################################### SUMMARY
            % --- spk sd
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii)]);
            ylabel('sd of mean spike');
            
            plot(SylType, SpksdAll, 'ok');
            ymean = grpstats(SpksdAll, SylType);
            x = unique(SylType)+0.1;
            plot(x, ymean, 'or-');
            ylim([0 max(SpksdAll)]);
            xlim([0 4]);
            
            % --- rho
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii)]);
            ylabel('rho of mean spike');
            x = unique(SylType)+0.1;
            plot(x, RhoAll, 'ok');
            ymean = grpstats(RhoAll, SylType);
            plot([1.1 2.1 3.1], ymean, 'or-');
%             plot(1:length(RhoAll), RhoAll, 'ok');
            lt_plot_zeroline;
            ylim([min(RhoAll) max(RhoAll)]);
            xlim([0 4]);
            
              % --- abs(rho)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii)]);
            ylabel('abs(rho) of mean spike');
             x = unique(SylType)+0.1;
           plot(x, abs(RhoAll), 'ok');
            ymean = grpstats(abs(RhoAll), SylType);
            plot([1.1 2.1 3.1], ymean, 'or-');
%             plot(1:length(RhoAll), RhoAll, 'ok');
            lt_plot_zeroline;
            ylim([0 max(abs(RhoAll))]);
            xlim([0 4]);

            % ================== CROSS CORRELATIONS (spk vs. ff)
           for j=1:3
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii)]);
            ylabel('xcorr (spk vs. ff)');
               indtmp = SylType==j;
               
               cc = Xcorr_NspkFF(indtmp,:);
%                indmid = ceil(size(cc,2)/2);
               flanks = floor(size(cc,2)/2);
               lags = -flanks:1:flanks;
               assert(length(lags) == size(cc,2), 'asfasd');
               
               if j==1
               plotcol = 'k';
               elseif j==2
                   plotcol = 'b';
               elseif j==3
                   plotcol = 'r';
               end
                  
            plot(lags, cc, '-', 'Color', plotcol);
            plot(lags, nanmean(cc,1), '-', 'LineWidth', 2, 'Color', plotcol);
            ylim([-0.5 0.5]);
            lt_plot_zeroline;
           end
           
           % ================== CROSS CORRELATIONS (spkSTD vs. ffMEDIAN)
           for j=1:3
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' exptname '-sw' num2str(iii)]);
            ylabel('xcorr (spkstd vs. ffmedian');
               indtmp = SylType==j;
               
               cc = Xcorr_NspkSTD_vs_FFmedian(indtmp,:);
%                indmid = ceil(size(cc,2)/2);
               flanks = floor(size(cc,2)/2);
               lags = -flanks:1:flanks;
               assert(length(lags) == size(cc,2), 'asfasd');
               
               if j==1
               plotcol = 'k';
               elseif j==2
                   plotcol = 'b';
               elseif j==3
                   plotcol = 'r';
               end
                  
            plot(lags, cc, '-', 'Color', plotcol);
            plot(lags, nanmean(cc,1), '-', 'LineWidth', 2, 'Color', plotcol);
            ylim([-0.5 0.5]);
            lt_plot_zeroline;
           end
            
           
           pause
            close all;
            
            
        end
    end
end


% ===================== 1) FF

% ===================== 2) spikes

% ==================== 3) summary


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


%% [TARG-SAME-DIFF] ##################### CROSS CORRELATIONS


onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 5;
onlyIfStartWNoff =0;
onlyIfHaveAllSylTypes = 1;
lt_figure; hold on;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

LearnAll = [];
YAll = {};
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            for j = 1:maxneuron
                
                
                % --------------------------- CHECK THAT HAS DATA
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                if ~any(inds)
                    continue
                end
                
                % ------------------------------ CHECK WHETHER START WN OFF
                if onlyIfStartWNoff==1
                    if unique(All_StartFromWNOff(inds))==0
                        continue
                    end
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
                
                
                
                Ythis = {};
                % ################################ COLLECT INTO VECTOR
                
                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1 & ~isnan(AllLearnSlopeScaled);
                
                if ~any(inds)
                    y = nan;
                    learnrate = nan;
                else
                learnrate = AllLearnSlopeScaled(inds);
                y = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        y = cell2mat(AllTrial_Spk_vs_FF(inds));
                        ylabstring = 'Spk vs. FF , xcorr';
                    case 1
                        y = cell2mat(AllTrialXcorr_NspkMedian_vs_FFmedian(inds));
                        ylabstring = 'smoothed spk vs. ff , xcorr';
                    case 2
                        y = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
                        ylabstring = 'smoothed spk(STD) vs. ff(MED) , xcorr';
                    case 3
                        y = cell2mat(AllBinned_SpkMean_vs_FFchange(inds));
                        ylabstring = 'Binned spk(mean) vs. ff(change) , xcorr';
                    case 4
                        y = cell2mat(AllBinned_FRsmChange_vs_FFchange(inds));
                        ylabstring = 'Binned FRsm(change) vs. ff(change) , xcorr';
                    case 5
                        y = cell2mat(AllTrial_Spk_vs_FFslope(inds));
                        ylabstring = 'trial, Spk vs. FFslope, xcorr';
                end
               
                y = nanmean(y,1);
                end
                Ythis = [Ythis y];
                
                
               % ======================================= SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==1 & ~isnan(AllLearnSlopeScaled);
                
                if ~any(inds)
                    y = nan;
                    learnrate = nan;
                else
                learnrate = AllLearnSlopeScaled(inds);
                y = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        y = cell2mat(AllTrial_Spk_vs_FF(inds));
                    case 1
                        y = cell2mat(AllTrialXcorr_NspkMedian_vs_FFmedian(inds));
                    case 2
                        y = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
                        ylabstring = 'smoothed spk(STD) vs. ff(MED) , xcorr';
                    case 3
                        y = cell2mat(AllBinned_SpkMean_vs_FFchange(inds));
                        ylabstring = 'Binned spk(mean) vs. ff(change) , xcorr';
                    case 4
                        y = cell2mat(AllBinned_FRsmChange_vs_FFchange(inds));
                        ylabstring = 'Binned FRsm(change) vs. ff(change) , xcorr';
                    case 5
                        y = cell2mat(AllTrial_Spk_vs_FFslope(inds));
                        ylabstring = 'trial, Spk vs. FFslope, xcorr';
                end
               
                y = nanmean(y,1);
                end
                Ythis = [Ythis y];
                
               % ======================================= DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==0 & ~isnan(AllLearnSlopeScaled);
                
                if ~any(inds)
                    y = nan;
                    learnrate = nan;
                else
                learnrate = AllLearnSlopeScaled(inds);
                y = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        y = cell2mat(AllTrial_Spk_vs_FF(inds));
                     case 1
                        y = cell2mat(AllTrialXcorr_NspkMedian_vs_FFmedian(inds));
                    case 2
                        y = cell2mat(AllTrialXcorr_NspkStd_vs_FFmedian(inds));
                        ylabstring = 'smoothed spk(STD) vs. ff(MED) , xcorr';
                    case 3
                        y = cell2mat(AllBinned_SpkMean_vs_FFchange(inds));
                        ylabstring = 'Binned spk(mean) vs. ff(change) , xcorr';
                    case 4
                        y = cell2mat(AllBinned_FRsmChange_vs_FFchange(inds));
                        ylabstring = 'Binned FRsm(change) vs. ff(change) , xcorr';
                    case 5
                        y = cell2mat(AllTrial_Spk_vs_FFslope(inds));
                        ylabstring = 'trial, Spk vs. FFslope, xcorr';
               end
               
                y = nanmean(y,1);
                end
                Ythis = [Ythis y];
               
                
                % ================================== COMBINED ALL 3 !!!
                YAll = [YAll; Ythis];
            end
        end
    end
end

    
% ======== 
indstoremove = [];
for j=1:size(YAll)
    
    if any(isnan(YAll{j,1})) | any(isnan(YAll{j,2})) | any(isnan(YAll{j,3}))
        indstoremove = [indstoremove j];
    end
        
end
YAll(indstoremove, :) = [];

lt_figure; hold on;
hsplots = [];
% === targ
hsplot = lt_subplot(3,1,1); hold on;
hsplots = [hsplots hsplot];
title('targ');
y = cell2mat(YAll(:,1));
x = -floor(size(y,2)/2):1:floor(size(y,2)/2);
plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(y), lt_sem(y), {'Color', 'r'},1);
lt_plot_zeroline;
ylabel(ylabstring);

% === same
hsplot = lt_subplot(3,1,2); hold on;
hsplots = [hsplots hsplot];
title('same');
y = cell2mat(YAll(:,2));
x = -floor(size(y,2)/2):1:floor(size(y,2)/2);
plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(y), lt_sem(y), {'Color', 'r'},1);
lt_plot_zeroline;

% === diff
hsplot = lt_subplot(3,1,3); hold on;
hsplots = [hsplots hsplot];
title('diff');
y = cell2mat(YAll(:,3));
x = -floor(size(y,2)/2):1:floor(size(y,2)/2);
plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(y), lt_sem(y), {'Color', 'r'},1);
lt_plot_zeroline;

% ------------
linkaxes(hsplots, 'xy');


% ############################################# PLOT CORR AT 0 LAG
YY = [];

% ----- TARG
indtmp = 1;
y = cell2mat(YAll(:,indtmp));
tmp = ceil(size(y,2)/2);
y = y(:,tmp);

YY = [YY y];

% ----- SAME
indtmp = 2;
y = cell2mat(YAll(:,indtmp));
tmp = ceil(size(y,2)/2);
y = y(:,tmp);

YY = [YY y];

% ----- DIFF
indtmp = 3;
y = cell2mat(YAll(:,indtmp));
tmp = ceil(size(y,2)/2);
y = y(:,tmp);

YY = [YY y];

lt_figure; hold on;
plot([1 2 3], YY, '-ok');
xlim([0 4]);
lt_plot([1 2 3]+0.1, mean(YY,1), {'Errors', lt_sem(YY), 'Color','r'});


%% [TARG-SAME-DIFF] ##################### if multiple motifs, then
% takes average --> i.e. each neuron contributes exactly 2 points (targ,
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =0;
plottype = 3;
onlyIfStartWNoff =0;
onlyIfHaveAllSylTypes = 1;
CombineWithinSwitch=1;

lt_figure; hold on;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

X = [];
Y = [];
ExptCount = [];
ecounter = 0;
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
          ecounter = ecounter+1;
            for j = 1:maxneuron
                
                
                % --------------------------- CHECK THAT HAS DATA
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                if ~any(inds)
                    continue
                end
                
                % ------------------------------ CHECK WHETHER START WN OFF
                if onlyIfStartWNoff==1
                    if unique(All_StartFromWNOff(inds))==0
                        disp('skipped - not WN off start');
                        continue
                    end
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
                
                
                % ################################ COLLECT INTO VECTOR
                ythisthis = [];
                xthisthis = [];
                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1 & ~isnan(AllLearnSlopeScaled);
                
                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO(inds);
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
                        SpkFFcorr = abs(All_SpkMean_slopeScaled(inds));
                    case 6
                        ylabel('std of mean spks base (CV)');
                        SpkFFcorr = All_SpkMean_STDtrials_base(inds);
                    case 7
                        ylabel('change in CV of mean spikes (train - base)');
                        SpkFFcorr = All_SpkMean_STDtrials(inds) - All_SpkMean_STDtrials_base(inds);
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO_Residuals(inds);
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = abs(AllDiff_nspk(inds));
                        ylabel('spk diff (late minus early)');
                end
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                
                ythisthis = [ythisthis SpkFFcorr];
                xthisthis = [xthisthis learnrate];
                
                
                % ======================================= SAME SYL
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==1 & ~isnan(AllLearnSlopeScaled);
                
                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO(inds);
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
                        SpkFFcorr = abs(All_SpkMean_slopeScaled(inds));
                    case 6
                        ylabel('std of mean spks base (CV)');
                        SpkFFcorr = All_SpkMean_STDtrials_base(inds);
                    case 7
                        ylabel('change in CV of mean spikes (train - base)');
                        SpkFFcorr = All_SpkMean_STDtrials(inds) - All_SpkMean_STDtrials_base(inds);
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr =AllTrial_Spk_vs_FF_RHO_Residuals(inds);
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = abs(AllDiff_nspk(inds));
                        ylabel('spk diff (late minus early)');
                end
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                ythisthis = [ythisthis SpkFFcorr];
                xthisthis = [xthisthis learnrate];
                
                
                
                % ======================================= DIFF SYL
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==0 & ~isnan(AllLearnSlopeScaled);

                learnrate = AllLearnSlopeScaled(inds);
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO(inds);
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
                        SpkFFcorr = abs(All_SpkMean_slopeScaled(inds));
                    case 6
                        ylabel('std of mean spks base (CV)');
                        SpkFFcorr = All_SpkMean_STDtrials_base(inds);
                    case 7
                        ylabel('change in CV of mean spikes (train - base)');
                        SpkFFcorr = All_SpkMean_STDtrials(inds) - All_SpkMean_STDtrials_base(inds);
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO_Residuals(inds);
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = abs(AllDiff_nspk(inds));
                        ylabel('spk diff (late minus early)');
                end
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                ythisthis = [ythisthis SpkFFcorr];
                xthisthis = [xthisthis learnrate];

                
                
                % ====================================== APPEND TO END
                X = [X; xthisthis];
                Y = [Y; ythisthis];
                ExptCount = [ExptCount; ecounter];



            end
        end
    end
end



% ============== only keep those with all 3 syl types if so desired
if onlyIfHaveAllSylTypes==1
   indstoremove = any(isnan(Y)');
   
   X(indstoremove, :) = [];
   Y(indstoremove, :) = [];    
   ExptCount(indstoremove, :) = [];
end

% =============== if combine to get one datapt per switch
if CombineWithinSwitch==1
    Xnew = [];
    Ynew = [];
    enums = unique(ExptCount);
    for j=enums'
        indtmp = ExptCount==j;
        
        x = mean(X(indtmp,:),1);
        y = mean(Y(indtmp,:),1);
        
        Xnew = [Xnew; x];
        Ynew = [Ynew; y];
        
    end
    X = Xnew;
    Y = Ynew;
end

% =================== METRIC VS. LEARNING
xlabel('learning');

nsyls = size(X,1);
for j=1:nsyls
   x = X(j,:); % learn
   y = Y(j,:); % metric
   plot(x,y,'-', 'Color', [0.7 0.7 0.7]);
   lt_plot(x(1), y(1), {'Color', 'k'}); % targ
   lt_plot(x(2), y(2), {'Color', 'b'}); % same
   lt_plot(x(3), y(3), {'Color', 'r'}); % diff
end

YLIM = ylim;

% =============================================== PAIRED COMPARISON
lt_figure; hold on;
xlabel('TARG -- SAME -- DIFF');

x = [1 2 3];
nsyls = size(X,1);
for j=1:nsyls
   y = Y(j,:); % metric
   plot(x,y,'-', 'Color', [0.7 0.7 0.7]);
   lt_plot(x(1), y(1), {'Color', 'k'}); % targ
   lt_plot(x(2), y(2), {'Color', 'b'}); % same
   lt_plot(x(3), y(3), {'Color', 'r'}); % diff
end

xlim([0 4]);

p = signrank(Y(:,2), Y(:,1));
lt_plot_text(1.2, max(Y(:,1)), ['p(srank)=' num2str(p)]);
p = signrank(Y(:,3), Y(:,2));
lt_plot_text(2.2, max(Y(:,2)), ['p(srank)=' num2str(p)]);


% ================================================= COMPARE TARG AND SAME
% AS FUNCTION OF LEARNING
lt_figure; hold on;
xlabel('ff slope (targ minus same)');
ylabel('metric (targ minus same)');
title('targ - same');

x = X(:,1) - X(:,2);
y = Y(:,1) - Y(:,2);
lt_regress(y, x, 1, 0, 1, 1);
lt_plot(x, y);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% [IMPORTANT] ######################################### V2 - if multiple motifs, then
% takes average --> i.e. each neuron contributes exactly 2 points (targ,
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =0;
plottype = 9;
onlyIfStartWNoff =0;
overlayIfSU=0;

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
                
                disp(rand)
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
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO(inds);
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
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO_Residuals(inds));
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = AllDiff_nspk(inds);
                        ylabel('spk diff (late minus early)');
                end
                
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                istarg = unique(istarg);
                
                X = [X; learnrate];
                Y = [Y; SpkFFcorr];
                TargStat = [TargStat; istarg];
                
                
                % ======================================= SAME SYL
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==1;
                
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
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO_Residuals(inds));
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = AllDiff_nspk(inds);
                        ylabel('spk diff (late minus early)');
                end
                
                learnrate = mean(learnrate);
                SpkFFcorr = mean(SpkFFcorr);
                istarg = unique(istarg);
                
                X = [X; learnrate];
                Y = [Y; SpkFFcorr];
                TargStat = [TargStat; istarg];
                
                
                
                % ======================================= DIFF SYL
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSameSyl==0 & ~isnan(AllLearnSlopeScaled);
                if sum(inds)>2
                keyboard
                end
                
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
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO_Residuals(inds));
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                    case 9
                        SpkFFcorr = AllDiff_nspk(inds);
                        ylabel('spk diff (late minus early)');
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
plottype = 4;
onlyIfStartWNoff =0;
overlayIfSU=0;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
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
            
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii;
                
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
                SpkFFcorr = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        SpkFFcorr = AllTrial_Spk_vs_FF_RHO(inds);
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
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO_Residuals(inds));
                        ylabel('spk vs. ff corr, using residuals (abs value)');
                end
                
                
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
                    case 8
                        % trial - spk vs, ff, corr (Residuals)
                        SpkFFcorr = abs(AllTrial_Spk_vs_FF_RHO_Residuals(inds));
                        ylabel('spk vs. ff corr, using residuals (abs value)');
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

