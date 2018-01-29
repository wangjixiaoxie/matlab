function lt_neural_v2_CTXT_BRANCH_PlotByBranchID(ALLBRANCH, BrainRegions, ...
    BirdToPlot, useDprime, ExptToPlot)
%% lt 11/29/17 - plots, sepaated by unique branch points (i.e. regexpstr)

%% ======= Filter data (e.g. remove noise, poor labels, etc)

Params.LocationsToKeep = BrainRegions;
Params.birdstoexclude = {}; % performs 1st
Params.birdstokeep = BirdToPlot; % performs 2nd
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below
Params.expttokeep = ExptToPlot;

ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%%
if useDprime==1
    % DPRIME STUFF
    Niter = 3; % shuffles
    Nmin = 3; % min sample size;
    % Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;
    
    % DprimeNegVersion = 'Wohl'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
    DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.
    
    ALLBRANCH = lt_neural_v2_CTXT_AllBranchDprime(ALLBRANCH, Nmin, Niter, ...
        DprimeNegVersion);
end

%% ===== organize data by unique branches
apos =1;
DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, '', '', ...
    useDprime);

%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            
            if ~isempty(ydecode)
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % -- lines
                line([0 0], ylim);
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                plot(xdecode, ydecode_neg, '-', 'Color', [0.6 0.6 0.6]);
                plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                plot(xdecode, ydecode, '-', 'Color', 'r');
                plot(xdecode, mean(ydecode,1), '-r', 'LineWidth', 2)
                
                ymax = max(ydecode(:));
                plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
            end
            
            % ---- note down neuron/channel number
            neuronOrig = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
            for j=1:length(neuronOrig)
                
                nn = neuronOrig(j);
            
                chan = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).channel;
                lt_plot_text(xdecode(end), ydecode(j,end), ['ch' num2str(chan)], 'b', 7);
                
            end
            
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
end


%% ======================= SUMMARY PLOTS FOR EACH REGION
figcount=1;
subplotrows=4;
subplotcols=max([2 length(locationstoplot)]);
fignums_alreadyused=[];
hfigs=[];
hsplots =[];



%% =========================== BY BIRD
numbirds = length(DATSTRUCT_BYBRANCH.bird);

for loc = locationstoplot
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BY BIRD]' loc{1}]);
    
    Yall = [];
    for i=1:numbirds
        birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        Ythisbird = [];
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            Ythisbird = [Ythisbird; ydecode];
            
        end
        
        % --- plot for this bird
        Ythisbirdmean = mean(Ythisbird,1);
        if ~isempty(Ythisbirdmean)
            plot(xdecode, Ythisbirdmean, '-k');
            % ==================== COLLECT ACROSS BIRDS
            Yall = [Yall; Ythisbirdmean];
            disp(birdname);
        end
        
    end
    
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
end



%% =========================== BY BRANCH
numbirds = length(DATSTRUCT_BYBRANCH.bird);

for loc = locationstoplot
    disp(loc)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BRANCHES]' loc{1}]);
    
    Yall = [];
    Yallneg = [];
    for i=1:numbirds
        %         birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            %% take mean across all neurons for this branch
            ydecode = mean(ydecode,1);
            ydecode_neg = mean(ydecode_neg, 1);
            
            %%
            
            if ~isempty(ydecode)
                plot(xdecode, ydecode, '-', 'Color', 'k');
                
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                %                 plot(xdecode, ydecode_neg, '-', 'Color', 'k');
                %                 plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                % -- lines
                line([0 0], ylim);
                
            end
            
            
            % ==================== COLLECT
            Yall = [Yall; ydecode];
            Yallneg = [Yallneg; ydecode_neg];
        end
    end
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
    
    % -- neg
    YallmeanNEG = mean(Yallneg,1);
    YallsemNEG = lt_sem(Yallneg);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, YallmeanNEG, YallsemNEG, {'Color', 'k'},1)
    end
    
end


%% =========================== BY BRANCH/NEURON
numbirds = length(DATSTRUCT_BYBRANCH.bird);


Ymean_loc =[];
Ysem_loc = [];
for loc = locationstoplot
    disp(loc)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BRANCHES/NEURON]' loc{1}]);
    
    Yall = [];
    Yall_neg = [];
    for i=1:numbirds
        %         birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            %%
            
            if ~isempty(ydecode)
                plot(xdecode, ydecode, '-', 'Color', [0.6 0.6 0.6]);
                
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                %                 plot(xdecode, ydecode_neg, '-', 'Color', 'k');
                %                 plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                % -- lines
                line([0 0], ylim);
                
            end
            
            
            % ==================== COLLECT
            Yall = [Yall; ydecode];
            Yall_neg = [Yall_neg; ydecode_neg];
        end
    end
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
    
    % --- put neg control
    Yallneg_mean = mean(Yall_neg,1);
    Yallneg_std = std(Yall_neg,0,1);
    Yallneg_sem = lt_sem(Yall_neg);
    shadedErrorBar(xdecode, Yallneg_mean, Yallneg_sem, {'Color', 'm'},1)
    
    
    % ================== COLLECT FOR THIS LOC
    Ymean_loc = [Ymean_loc; Yallmean];
    Ysem_loc = [Ysem_loc; Yallsem];
end

% ==================== PLOT ALL LOC COMBINED
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['[COMBINED]']);

for i=1:size(Ymean_loc,1)
    shadedErrorBar(xdecode, Ymean_loc(i,:), Ysem_loc(i,:), {'Color', 'r'}, 1);
end


linkaxes(hsplots, 'xy');


if length(BrainRegions)==2
    %% ==== cross correlation between each pair of brain regions
    if (1) % OLD VERSION !!!!!!!!!!!!!!!!!
%         datwindow = [-0.1 0.05]; % data to use, rel onset
%         windowmax = 0.05; % in sec
        datwindow = []; % data to use, rel onset
        windowmax = 0.1; % in sec
        
        YallMat = [];
        binsizeholder = [];
        figure; hold on;
        for i=1:numbirds
            numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
            
            for bb=1:numbranches
                
                thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
                Yall = cell(1, length(locationstoplot));
                
                for k = 1:length(locationstoplot)
                    loc = locationstoplot{k};
                    
                    %  ================ GET brain region
                    inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc); % neurons
                    if ~any(inds)
                        continue
                    end
                    % -------------- GET DAT
                    xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
                    binsize = xdecode(2)-xdecode(1);
                    binsizeholder = [binsizeholder binsize];
                    ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                    ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                    
                    if ~isempty(datwindow)
                       indstmp = xdecode>=datwindow(1) & xdecode<=datwindow(2);
                       
                       ydecode = ydecode(:,indstmp);
                       ydecode_neg = ydecode_neg(:, indstmp);
                        
                    end
                    % ==================== COLLECT
                    Yall{k} = ydecode;
                    
                end
                
                % ==== CALCULATE ALL CROSS CORRELATIONS (between all pairs of
                % neurons)
                nneur1 = size(Yall{1},1);
                nneur2 = size(Yall{2},1);
                for k = 1:nneur1
                    for kk = 1:nneur2
                        
                        dat1 = Yall{1}(k,:);
                        dat2 = Yall{2}(kk,:);
                        
                        [cc, lags] = xcov(dat1, dat2, ceil(windowmax/binsize), 'coeff');
%                         [cc, lags] = xcorr(dat1, dat2, ceil(windowmax/binsize));
                        plot(lags*binsize, cc);
                        YallMat = [YallMat; cc];
                    end
                end
                
                
            end
        end
        
        assert(length(unique(binsizeholder))==1, 'asdfsd');
        
        plot(lags*binsize, mean(YallMat, 1), 'k');
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
        
    end
    
    %% ==== cross correlation between each pair of brain regions
    numbirds = length(DATSTRUCT_BYBRANCH.bird);
    % #################################### 1) COLLECT ALL DATA INTO ARRAYS
    YallMat = [];
    XallMat = [];
    BranchID = [];
    LocationAll = [];
    BirdID = [];
    
    binsizeholder = [];
    
    for i=1:numbirds
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            for k = 1:length(locationstoplot)
                loc = locationstoplot{k};
                
                %  ================ GET neurons
                inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc); % neurons
                if ~any(inds)
                    continue
                end
                
                % --------------- COLLECT DAT
                xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(inds,:);
                ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                
                
                % ---- collect all binsizes - will make sure all are equal
                binsize = xdecode(1,2)-xdecode(1,1);
                binsizeholder = [binsizeholder binsize];
                
                % ==================== COLLECT
                YallMat = [YallMat; ydecode];
                XallMat = [XallMat; xdecode];
                BranchID = [BranchID; bb*ones(size(ydecode,1),1)];
                LocationAll = [LocationAll; k*ones(size(ydecode,1),1)];
                BirdID = [BirdID; i*ones(size(ydecode,1),1)];
                
            end
            
            %         % ==== CALCULATE ALL CROSS CORRELATIONS (between all pairs of
            %         % neurons)
            %         nneur1 = size(Yall{1},1);
            %         nneur2 = size(Yall{2},1);
            %         for k = 1:nneur1
            %             for kk = 1:nneur2
            %
            %                 dat1 = Yall{1}(k,:);
            %                 dat2 = Yall{2}(kk,:);
            %
            %                 [cc, lags] = xcov(dat1, dat2, ceil(windowmax/binsize), 'coeff');
            %                 plot(lags*binsize, cc);
            %                 YallMat = [YallMat; cc];
            %             end
            %         end
            %
            %
        end
    end
    
    assert(length(unique(binsizeholder))==1, 'asdfsd');
    
    
    
    %% #####################################
    %% ============== CALCULATE PEAKS IN XCORR AND ALSO PLOT
    windowmax = 0.1; % sec
    
    [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll, ...
        YallMat, XallMat, locationstoplot);
    
    %% =================== SAVE LOCATION OF PEAK IN GRAND MEAN (PER BIRD)
    
    CCmean_dat = nan(size(CCall,2), numbirds);
    TimeOfPeak_dat = nan(1, numbirds);
    
    for i=1:numbirds
        
        inds = PairedBirdID==i;
    
        if ~any(inds)
            continue
        end
        
        % ------ collect mean CC
        ccmean = mean(CCall(inds,:),1);
        CCmean_dat(:, i) = ccmean;
                
        
        % ----------- find main peak of mean xcov
        [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
        timeofpeak = LagsAll(1,loctmp)*binsize;
        TimeOfPeak_dat(i) = timeofpeak;
    end
    
    
    %% ########################### NEGATIVE CONTROL - SHUFFLE LOCATION LABEL
    % ----
    nshuff = 1000;
    TimeOfPeak_Shuff = nan(nshuff, numbirds);
    CCmean_Shuff = nan(size(CCall,2), nshuff, numbirds);
    
    for nn=1:nshuff
        disp(['shuff ' num2str(nn)]);
        LocationAll_RAND = nan(size(LocationAll));
        
        % ------------ RANDOMIZE NEURON LOCATION (MAINTAINING BIRD/BRANCH ID)
        for i=1:numbirds
            numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
            for bb = 1:numbranches
%                 disp([num2str(i) '-' num2str(bb)]);
                inds = find(BirdID==i & BranchID==bb);
                
%                 if isempty(inds)
% %                     LocationAll_RAND(inds) = 97; % put some crazy number. don't want to be nan, but don't want it to affect actual data
%                     continue
%                 end
                
                indsperm = inds(randperm(length(inds)));
                
%                 assert(all(isnan(LocationAll_RAND(inds))), 'asdf');
                LocationAll_RAND(inds) = LocationAll(indsperm);
            end
        end
        
        assert(~any(isnan(LocationAll_RAND)), 'asdfds');
        
        
        % =========================== CALC COV
        if nn>2
        suppressplots = 1;
        else
            suppressplots=0;
        end
        [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll_RAND, ...
            YallMat, XallMat, locationstoplot, suppressplots);
        
        % =========================== CALCUALTE TIME OF PEAK (and collect
        % mean timecourse
        for i=1:numbirds
            
            inds = PairedBirdID==i;
            
            if ~any(inds)
                continue
            end
            
                % ------ collect mean CC
            ccmean = mean(CCall(inds,:),1);
            CCmean_Shuff(:, nn, i) = ccmean;

            % ----------- find main peak of mean xcov
            [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
            timeofpeak = LagsAll(1,loctmp)*binsize;
            
            TimeOfPeak_Shuff(nn, i) = timeofpeak;
        end
        
        
    end
   
    %% ====================== COMPARE SHUFFLES TO DATA


    % =========================1 ) time of peak
    figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:numbirds
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
    xcenters = min(TimeOfPeak_Shuff(:,i)-0.005):0.005:(max(TimeOfPeak_Shuff(:,i))+0.005);
    lt_plot_histogram(TimeOfPeak_Shuff(:,i), xcenters, 1, 1, '', 1, 'k');
    line([TimeOfPeak_dat(i) TimeOfPeak_dat(i)], ylim, 'Color', 'r', 'LineWidth', 2);
    
    % --- p value
    p = (sum(abs(TimeOfPeak_Shuff(:,i)) > abs(TimeOfPeak_dat(i))) + 1)./(length(TimeOfPeak_Shuff(:,i)) + 1);
    lt_plot_pvalue(p, '2-tailed',1);
    xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
    
    % ======================= 2) overlay data mean CC on shuffles
    figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
    for i=1:numbirds
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        % ---- overlay all shuffle dat
        ccmeanall = CCmean_Shuff(:,:,i);
        plot(LagsAll(1,:)*binsize, ccmeanall', ':', 'Color', [0.7 0.7 0.7]);
        plot(LagsAll(1,:)*binsize, CCmean_dat(:, i), 'Color', 'r');
        xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
    
    % ======================== 3) overlay mean on shuff mean + sd
    figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
    for i=1:numbirds
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        % ---- overlay all shuffle dat
        ccmeanall = CCmean_Shuff(:,:,i);
        ccmeanmean = mean(ccmeanall, 2);
        ccmeanstd = std(ccmeanall, 0, 2);
        ccmean95 = prctile(ccmeanall', [2.5 97.5]);
%         ccmean95 = prctile(ccmeanall', [2.5 97.5]);
%         ccmean95 =  ccmean95 - repmat(ccmeanmean', 2, 1);
        plot(LagsAll(1,:)*binsize, ccmeanmean, '-k', 'LineWidth', 2);
        plot(LagsAll(1,:)*binsize, ccmean95(1,:)', '-k');
        plot(LagsAll(1,:)*binsize, ccmean95(2,:)', '-k');
%                 
%         
%         shadedErrorBar(LagsAll(1,:)*binsize, ccmeanmean, ccmeanstd, {'Color', 'k'}, 1);
        plot(LagsAll(1,:)*binsize, CCmean_dat(:, i), '-r', 'LineWidth', 3);
        xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
end
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS


function [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll, ...
    YallMat, XallMat, locationstoplot, suppressplots)
if ~exist('suppressplots', 'var')
    suppressplots=0;
end

%% ############################### CALCULATE CROSS CORRELATIONS
% for a given bird/branchID, all pairwise between neurons
% windowmax = 0.1; % in sec

numbirds = max(BirdID);
numbranches = max(BranchID);
numloc = max(LocationAll);

% --- for output
CCall = [];
LagsAll = [];
PairedBirdID = [];
PairedBranchID = [];


% --- RUN
for i=1:numbirds
    for bb = 1:numbranches
        
        inds = BirdID==i & BranchID==bb;
        if ~any(inds)
            continue
        end
        
        % ======================== EXTRACT DAT FOR THIS BRANCH
        loc = LocationAll(inds);
        ymat = YallMat(inds,:);
        xmat = XallMat(inds, :);
        binsize = xmat(1,2) - xmat(1,1);
        
        if all(isnan(loc))
            continue
        end
        
        % ======================= MAKE SURE DAT IS USABLE.
        if length(unique(loc))==1
            %  then this branch only has data from one brain region ...
            continue
        end
        
        if length(unique(loc)) > 2
            % -- then this branch has data from too manyr egions ...
            disp(unique(loc));
            keyboard
            disp('=== NOTE: skipped a branch because has too many brain regions (>2)');
            continue
        end
        
        % ====================== SEPARATE INTO 2 BRAIN REGIONS
        loclist = unique(loc)';
        
        % ----------- separation matrix for each brain region
        y1 = ymat(loc==loclist(1), :);
        y2 = ymat(loc==loclist(2), :);
        
        
        for j = 1:size(y1,1)
            
            for jj=1:size(y2,1)
                
                [cc, lags] = xcov(y1(j,:), y2(jj,:), ceil(windowmax/binsize), 'coeff');
%                 [cc, lags] = xcorr(y1(j,:), y2(jj,:), ceil(windowmax/binsize), 'coeff');
                
                
                % ------------------------- OUTPUT
                CCall = [CCall; cc];
                LagsAll = [LagsAll; lags];
                
                PairedBirdID = [PairedBirdID; i];
                PairedBranchID = [PairedBranchID; bb];
                
                
            end
            
        end
        
    end
end

if suppressplots==0
    %% ======= PLOTS (BY BIRD)
    for i=1:numbirds
        figcount=1;
        subplotrows=3;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        
        % ========== all dat
        inds = PairedBirdID==i;
        if ~any(inds)
            continue
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('all neuron pairs');
        plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
        
        shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
        lt_plot_zeroline_vert;
        
        title(['bird ' num2str(i)]);
        
        % ----------- find and plot peaks
        [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
        line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
        
        for j=1:numbranches
            
            
            inds = PairedBirdID==i & PairedBranchID==j;
            
            if all(inds==0)
                continue
            end
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['branch ' num2str(j)]);
            plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
            
            shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
            xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
            lt_plot_zeroline_vert
            
            % ----------- find and plot peaks
            [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
            if ~isempty(loctmp)
            line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
            end
            
        end
    end
    
    
    %% ======= PLOTS (GRAND AVE)
    lt_figure; hold on;
    
    lt_subplot(2,2,1); hold on;
    title('all neuron pairs');
    plot(lags*binsize, CCall', 'b');
    xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
    
    lt_subplot(2,2,2); hold on;
    title('mean(sem)');
    shadedErrorBar(lags*binsize, mean(CCall, 1), lt_sem(CCall), {'Color', 'k'},1);
    xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
    
end


end

