function DATSTRUCT = lt_neural_v2_CTXT_BRANCH_DatVsShuffPLOT(analyfname)
%% lt 10/27/17 - plots CLASSES, after runninglt_neural_v2_CTXT_BRANCH_DatVsShuff
% - plots decoding of dat vs. shuffled
% - currently only works with one time bin (the first) - can easily modify
% to plot multiple time bins.

%%

decodestat = 'F1';
plotdecodenegdistr = 0;

%% load branch
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
end


%% use first time bin only (can modify)

tt = 1;
%% COLLECT

% ---- to visualize decodeneg distributions
figcount=1;
subplotrows=8;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


AllPdat = [];
AllBirdNum = [];
AllNeurNum = [];
AllBranchNum = [];
AllDecode = [];
AllDecode_z = [];

AllWindows_RelOnsActual = [];
AllWindows_RelOnsOffDesired = [];

numbirds = length(CLASSES.birds);
for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    
    for ii=1:numneurons
        
        numbranch = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for iii=1:numbranch
            
            disp(['brd' num2str(i) '-n' num2str(ii) '-br' num2str(iii)]);
            
            datstruct = CLASSES.birds(i).neurons(ii).branchnum(iii);
            
            if ~isfield(datstruct, 'SHUFFDECODE')
                disp('asdfasd'); % then lacks data
                keyboard
                continue
            end
            if isempty(datstruct.SHUFFDECODE)
                disp('asdfads');
                continue
            end
            
            
            % === other features
            AllBirdNum = [AllBirdNum; i];
            AllNeurNum = [AllNeurNum; ii];
            AllBranchNum = [AllBranchNum; iii];
            
            
            
            % === Pdats
            AllPdat = [AllPdat; datstruct.SHUFFDECODE.timebin(tt).Pdat];
            
            
            % === convert decode into z-score rel to null distribution
            cmat = datstruct.SHUFFDECODE.timebin(tt).ConfMatAll_NEG ;
            decodeneg = [];
            for j=1:length(cmat)
                sts = lt_neural_ConfMatStats(cmat{j});
                decodeneg = [decodeneg sts.(decodestat)];
            end
            
            cmat = datstruct.SHUFFDECODE.timebin(tt).ConfMatAll_DAT;
            decodedat = [];
            for j=1:length(cmat)
                sts = lt_neural_ConfMatStats(cmat{j});
                decodedat = [decodedat sts.(decodestat)];
            end
            decode = mean(decodedat);
            
            % -- out
            AllDecode = [AllDecode decode];
            
            nullmean = mean(decodeneg);
            nullstd = std(decodeneg);
            AllDecode_z = [AllDecode_z (decode - nullmean)/nullstd];
            
            
            % --------------------------- check neg decode distrubtions
            if plotdecodenegdistr ==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                lt_plot_histogram(decodeneg);
                line([decode decode], ylim, 'Color','r')
            end
            
            % ----------------- other params
            windowrelonset = datstruct.SHUFFDECODE.timebin(tt).window_relonset;
            AllWindows_RelOnsActual = [AllWindows_RelOnsActual; windowrelonset];
            
            windowrelOnsOff = CLASSES.SHUFFDECODEpar.TimeWindows_relOnsetOffset;
            AllWindows_RelOnsOffDesired = [AllWindows_RelOnsOffDesired; windowrelOnsOff];
        end
    end
end


%% plot grand histogram
for i=1:size(AllPdat,2)
    lt_figure; hold on;
    title('prob of dat rel shuff distribution');
    lt_plot_histogram(log10(AllPdat(:,i)+0.0001))
end

%% pval correlate with zscore
lt_figure; hold on;
plot(AllDecode_z, AllPdat, '.k');


%% -- separate plot for each bird, separated by neuron num
numbirds = max(AllBirdNum);
numneur = max(AllNeurNum);
numbranch = max(AllBranchNum);

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ################### PVAL
for i=1:numbirds
    
    Yvals = {};
    XNeurs = [];
    birdname = CLASSES.birds(i).birdname;
    for ii=1:numneur
        
        inds = AllBirdNum==i & AllNeurNum==ii;
        
        if ~any(inds)
            continue
        end
        
        % ==== extract all p-vals for this bird/neuron (i.e. across
        % branches)
        Yvals = [Yvals log10(AllPdat(inds)+0.0001)];
        XNeurs = [XNeurs ii];
    end
    
    % =========== plot for this bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(birdname);
    xlabel('neur num');
    ylabel('log10(prob)');
    lt_plot_MultDist(Yvals, XNeurs, 1, 'k', 1, 0);
    line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
    ylim([-5 1]);
    
    
end


% ################### DECODE (zscore)
for i=1:numbirds
    
    Yvals = {};
    XNeurs = [];
    birdname = CLASSES.birds(i).birdname;
    for ii=1:numneur
        
        inds = AllBirdNum==i & AllNeurNum==ii;
        
        if ~any(inds)
            continue
        end
        
        % ==== extract all p-vals for this bird/neuron (i.e. across
        % branches)
        Yvals = [Yvals AllDecode_z(inds)];
        XNeurs = [XNeurs ii];
    end
    
    % =========== plot for this bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(birdname);
    xlabel('neur num');
    ylabel('decode (zscore)');
    lt_plot_MultDist(Yvals, XNeurs, 1, 'k', 1, 0);
    %     line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
    %     ylim([-5 1]);
    lt_plot_zeroline;
    
end




%% -- separate plot for each bird, separated by branch num
numbirds = max(AllBirdNum);
numneur = max(AllNeurNum);
numbranch = max(AllBranchNum);

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ############################# PVAL
BranchnameAll = {};
for i=1:numbirds
    
    Yvals = {};
    Xbranch = [];
    birdname = CLASSES.birds(i).birdname;
    
    for ii=1:numbranch
        
        inds = AllBirdNum==i & AllBranchNum==ii;
        
        if ~any(inds)
            continue
        end
        
        % ==== extract all p-vals for this bird/branch (i.e. across
        % branches)
        Yvals = [Yvals log10(AllPdat(inds)+0.0001)];
        Xbranch = [Xbranch ii];
        
        % --- what is the name of this branch?
        tmp = find(inds);
        branchname = CLASSES.birds(i).neurons(AllNeurNum(tmp(1))).branchnum(ii).regexprstr;
        BranchnameAll = [BranchnameAll branchname];
    end
    % =========== plot for this bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(birdname);
    xlabel('branch ID');
    ylabel('log10(prob)');
    lt_plot_MultDist(Yvals, Xbranch, 1, 'k', 1, 0);
    line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
    ylim([-5 1]);
    
    set(gca, 'XTickLabel', BranchnameAll);
    rotateXLabels(gca, 45);
    
end


% ############################# DECODE (ZSCORE)
BranchnameAll = {};
for i=1:numbirds
    
    Yvals = {};
    Xbranch = [];
    birdname = CLASSES.birds(i).birdname;
    
    for ii=1:numbranch
        
        inds = AllBirdNum==i & AllBranchNum==ii;
        
        if ~any(inds)
            continue
        end
        
        % ==== extract all p-vals for this bird/branch (i.e. across
        % branches)
        Yvals = [Yvals AllDecode_z(inds)];
        Xbranch = [Xbranch ii];
        
        % --- what is the name of this branch?
        tmp = find(inds);
        branchname = CLASSES.birds(i).neurons(AllNeurNum(tmp(1))).branchnum(ii).regexprstr;
        BranchnameAll = [BranchnameAll branchname];
    end
    % =========== plot for this bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(birdname);
    xlabel('branch ID');
    ylabel('decode (z)');
    lt_plot_MultDist(Yvals, Xbranch, 1, 'k', 1, 0);
    lt_plot_zeroline;
    
    set(gca, 'XTickLabel', BranchnameAll);
    rotateXLabels(gca, 45);
    
end


%% -- separate plot for each bird, matrix of neuron x branch
numbirds = max(AllBirdNum);
numneur = max(AllNeurNum);
numbranch = max(AllBranchNum);

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];



for i=1:numbirds
    
    for ii=1:numbranch
        
        
    end
    
end




%% =========== for output

DATSTRUCT.AllPdat = AllPdat;
DATSTRUCT.AllBirdNum = AllBirdNum;
DATSTRUCT.AllNeurNum = AllNeurNum;
DATSTRUCT.AllBranchNum = AllBranchNum;
DATSTRUCT.AllDecode = AllDecode;
DATSTRUCT.AllDecode_z = AllDecode_z;
DATSTRUCT.AllWindows_RelOnsActual =AllWindows_RelOnsActual;
DATSTRUCT.AllWindows_RelOnsOffDesired = AllWindows_RelOnsOffDesired;







