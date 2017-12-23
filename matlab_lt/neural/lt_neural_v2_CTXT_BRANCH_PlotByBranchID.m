function lt_neural_v2_CTXT_BRANCH_PlotByBranchID(ALLBRANCH, BrainRegions, ...
    BirdToPlot)
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

ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%% ===== organize data by unique branches
apos =1;
DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH);

%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;


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


