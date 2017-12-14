function lt_neural_v2_CTXT_BranchCompareTwo(branchfname1, branchfname2)
%% lt 10/23/17 - get two save branch structs (e.g. two brain regions) and compare

% First thing is to compare the timing of decrease in context decoding
% relative to syllable onset

%% TO DO:
% go with the combined model, plot the results from that. look more closely
% at the cross validation stuff.

%% params
dattoplot = 'classperform';
LMANorX = 0; % 0, both; 1, LMAN; 2, X
birdstoexclude = {};

% durThreshOmega.syl = [0.15]; % omega2 (will only keep if lower) [leave empty to ignore]
% durThreshOmega.gappre= [];
% durThreshOmega.gappost= [0.15];
durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
durThreshOmega.gappre= [];
durThreshOmega.gappost= [];

RemoveRepeats = 0;


%% load each of the structures\
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

BRANCHES(1).dat = load([savedir '/' branchfname1]);
BRANCHES(2).dat = load([savedir '/' branchfname2]);

% === get gap duration anovas
if isempty(strfind(branchfname1, 'ALLBRANCHv2'))
    BRANCHES(1).dat.ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(BRANCHES(1).dat.ALLBRANCH);
end

if isempty(strfind(branchfname2, 'ALLBRANCHv2'))
    BRANCHES(2).dat.ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(BRANCHES(2).dat.ALLBRANCH);
end

% === get naem of struct
indstmp = strfind(branchfname1, '_');
BRANCHES(1).ID = branchfname1(indstmp(end)+1:end-4);

indstmp = strfind(branchfname2, '_');
BRANCHES(2).ID = branchfname2(indstmp(end)+1:end-4);

%% clean up

Params.LocationsToKeep = {};
Params.birdstoexclude = {};
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below

for i=1:length(BRANCHES)
    BRANCHES(i).dat.ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(...
        BRANCHES(i).dat.ALLBRANCH, Params);
end


%% convert to vectors

for i=1:length(BRANCHES)
    BRANCHES(i).DATSTRUCT = lt_neural_v2_CTXT_BranchToDatstruct(BRANCHES(i).dat.ALLBRANCH, birdstoexclude, ...
        LMANorX, RemoveRepeats, durThreshOmega, dattoplot);
end


%% ############################## METHOD 1 - COLLECTS EACH BIN FIRST

%% one plot each across all RA and all LMAN
plotScatterForEachBin = 0;


if plotScatterForEachBin==1
    Xallall = {};
    Yallall = {};
    IDallall = {};
    AlignnumAllall = [];
    YsemAllall = {};
    
    for i=1:length(BRANCHES);
        
        numalign = length(BRANCHES(i).DATSTRUCT.numalign);
        
        for aa = 1:numalign
            
            figcount=1;
            subplotrows=5;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            
            
            NN = length(BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Xcell); % num total branches/neurons
            
            
            BinMidList = [-0.1175:0.005:0.0725];
            SlopesAll = [];
            SlopesCIAll = [];
            PostMinusPreAll =[];
            PostMinusPreSEMAll = [];
            for bb = BinMidList;
                
                binEdges = [bb-0.02 bb bb+0.02];
                
                % ################################## EXTRACT ALL PRE AND POST
                % (FOR THIS BIN)
                Yall1 = [];
                Yall2= [];
                for nn = 1:NN
                    
                    xvals = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Xcell{nn};
                    yvals = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Ycell{nn};
                    assert(length(xvals) == length(yvals), 'asfd');
                    
                    % --- first half
                    inds = xvals>binEdges(1) & xvals<binEdges(2);
                    ymean1 = mean(yvals(inds));
                    
                    % --- second half
                    inds = xvals>binEdges(2) & xvals<binEdges(3);
                    ymean2 = mean(yvals(inds));
                    
                    % ---- plot
                    %             plot(ymean1, ymean2, 'ok');
                    Yall1 = [Yall1 ymean1];
                    Yall2= [Yall2 ymean2];
                end
                
                % =============================== PLOT (this bin)
                if plotScatterForEachBin ==1
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['algnSyl' num2str(aa) ',binMid' num2str(bb) ', ' BRANCHES(i).ID]);
                    xlabel('mean decode (pre)');
                    ylabel('mean decode (post)');
                    
                    [~,~,~,~,~,SummaryStats] = lt_regress(Yall2, Yall1, 1, 0, 1, 1, 'k');
                    xlim([0 1]); ylim([0 1]);
                    line([0 1], [0 1]);
                elseif plotScatterForEachBin==0
                    
                    [~,~,~,~,~,SummaryStats] = lt_regress(Yall2, Yall1, 0, 0, 0, 0, 'k');
                    
                end
                
                
                % ######################################## COLLECT
                % --------------------------------------- COLLECT SLOPES
                SlopesAll = [SlopesAll SummaryStats.slope];
                SlopesCIAll = [SlopesCIAll; SummaryStats.slopeCI];
                
                % -------------------------------------- COLLECT POST MINUS PRE
                postminuspre = nanmean(Yall2 - Yall1);
                %             postminuspre = nanmean(Yall2./Yall1);
                postminuspreSEM = lt_sem(Yall2 - Yall1);
                
                PostMinusPreAll = [PostMinusPreAll postminuspre];
                PostMinusPreSEMAll = [PostMinusPreSEMAll postminuspreSEM];
                
            end
            
            % --------------------------------------------- PLOT (slopes vs.
            % bin mid across all bins)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['algnSyl' num2str(aa) ', ' BRANCHES(i).ID]);
            xlabel('time at center (rel syl onset)');
            ylabel('slope of regression (2nd half/1st half)');
            
            %         plot(BinMidList, SlopesAll, 'ok-');
            errorbar(BinMidList, SlopesAll, SlopesAll'-SlopesCIAll(:,1), SlopesCIAll(:,2)-SlopesAll', '.k-');
            
            
            % --------------------------------------------- PLOT
            % (postminuspre)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['algnSyl' num2str(aa) ', ' BRANCHES(i).ID]);
            xlabel('time at center (rel syl onset)');
            ylabel('2nd half - 1st half');
            
            %         plot(BinMidList, SlopesAll, 'ok-');
            errorbar(BinMidList, PostMinusPreAll, PostMinusPreSEMAll, PostMinusPreSEMAll, '.k-');
            
            
            % ########################################### COLLECT ACROSS ALL
            % BRANCH
            Xallall = [Xallall BinMidList];
            Yallall = [Yallall PostMinusPreAll];
            YsemAllall = [YsemAllall PostMinusPreSEMAll];
            IDallall = [IDallall BRANCHES(i).ID];
            AlignnumAllall = [AlignnumAllall aa];
            
            
        end
    end
    
    
    %% ############################################## PLOT ALL BRANCH OVERLAID
    lt_figure; hold on;
    
    plotcols = lt_make_plot_colors(length(Xallall), 0, 0);
    
    for i=1:length(Xallall)
        
        x = Xallall{i};
        y = Yallall{i};
        ysem = YsemAllall{i};
        ID = IDallall{i};
        alignnum = AlignnumAllall(i);
        
        errorbar(x, y, ysem, '.-', 'Color', plotcols{i});
    end
    legend(IDallall)
    lt_plot_zeroline
    
end


%% ########################################### METHOD 2 - get timecourse for each neuron/branch

AllBRANCHSTRUCTnum = [];
AllAlignNum = [];
AllBirdnum = [];
AllNeurnum = [];
AllBranchID = [];
AllDifferentials = [];
BinMidList = [-0.0875:0.005:0.0625];
% BinMidList = [-0.0925:0.005:0.0625];
AllDecodeX = [];
AllDecodeY = [];
AllDecodeY_neg = [];
AllDecodeX_neg = [];
for i=1:length(BRANCHES)
    
    numalign = length(BRANCHES(i).DATSTRUCT.numalign);
    
    for aa = 1:numalign
        
        NN = length(BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Xcell); % num total branches/neurons
        
        for nn = 1:NN
            
            xvals = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Xcell{nn};
            yvals = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.Ycell{nn};
            birdnum = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.BirdnumAll(nn);
            neurnum = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.NeuronnumAll(nn);
            branchnum = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.Dat.BranchmumAll(nn);
            
            yvals_neg = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.NegContr.Ycell{nn};
            xvals_neg = BRANCHES(i).DATSTRUCT.numalign(aa).Datstruct.NegContr.Xcell{nn};
            
            % =============== evaluate differential at each time bin
            Differentials =[];
            for bb = BinMidList
                
                binEdges = [bb-0.02 bb bb+0.02];
                
                % ----------------------------- EXTRACT ALL PRE AND POST
                % (FOR THIS BIN)
                % --- first half
                inds = xvals>binEdges(1) & xvals<binEdges(2);
                ymean1 = mean(yvals(inds));
                
                % --- second half
                inds = xvals>binEdges(2) & xvals<binEdges(3);
                ymean2 = mean(yvals(inds));
                
                Differentials = [Differentials ymean2-ymean1];
            end
            
            
            % ==================================== COLLECT ALL DIFFERENTIAL
            AllBRANCHSTRUCTnum = [AllBRANCHSTRUCTnum; i];
            AllAlignNum = [AllAlignNum; aa];
            AllBirdnum = [AllBirdnum; birdnum];
            AllNeurnum = [AllNeurnum; neurnum];
            AllBranchID = [AllBranchID; branchnum];
            AllDifferentials = [AllDifferentials; Differentials];
            
            
            
            % ============================= COLLECT DECODE
            indstokeep = xvals>BinMidList(1)-0.02 & xvals<BinMidList(end)+0.02;
            
            xtokeep = xvals(indstokeep);
            ytokeep = yvals(indstokeep);
            xtokeep_neg = xvals_neg(indstokeep);
            ytokeep_neg = yvals_neg(indstokeep);

            AllDecodeX = [AllDecodeX; xtokeep];
            AllDecodeY = [AllDecodeY; ytokeep];
            
            AllDecodeY_neg = [AllDecodeY_neg; ytokeep_neg];
            AllDecodeX_neg = [AllDecodeX_neg; xtokeep_neg];

        end
    end
end


%% 
assert(max(AllAlignNum)==1, 'some plots assumes only one align num .. to solve put into extra layer of iteration over alignnums');

numbirds = max(AllBirdnum);
numanalys = max(AllBRANCHSTRUCTnum);

%% ============================ PLOTS (DECODING)
lt_figure; hold on;
hsplots = [];
% ====================== Plot MEAN OVER EVERYTHING
hsplot = lt_subplot(3,2,1); hold on;
title('grand means')
plotcols = lt_make_plot_colors(length(BRANCHES), 0,0);
hsplots =[hsplots hsplot];
for i=1:max(AllBRANCHSTRUCTnum)
    
    inds = AllBRANCHSTRUCTnum==i;
    
    x = mean(AllDecodeX(inds,:),1);
    y = mean(AllDecodeY(inds,:),1);
    ysem = lt_sem(AllDecodeY(inds, :));
    
    shadedErrorBar(x, y, ysem, {'Color',plotcols{i}},1);
    
    % ------------------ also plot negative control
    x = mean(AllDecodeX_neg(inds,:),1);
    y = mean(AllDecodeY_neg(inds,:), 1);
    shadedErrorBar(x, y, ysem, {'Color','k'},1);
    plot(x, y, '-k', 'LineWidth', 1, 'Color', plotcols{i});
end
line([0 0], ylim);


% ================================== EACH BIRD ONE LINE
for i=1:numanalys
    hsplot = lt_subplot(3,2,2+i); hold on;
    hsplots = [hsplots hsplot];
    count = 1;
    for j=1:numbirds
        
       inds = find(AllBRANCHSTRUCTnum==i & AllBirdnum==j);
       if isempty(inds) 
           continue
       end
       
       % ===== plot this bird
       ymean = mean(AllDecodeY(inds, :),1);
       x = mean(AllDecodeX(inds, :),1);
       ysem = lt_sem(AllDecodeY(inds,:));
       
       shadedErrorBar(x, ymean, ysem, {'Color', plotcols{i}},1);
       
        % ----- plot neg control
       ymean = mean(AllDecodeY_neg(inds, :),1);
       x = mean(AllDecodeX_neg(inds, :),1);
%        ysem = lt_sem(AllDecodeY(inds,:));
       
%        shadedErrorBar(x, ymean, ysem, {'Color', plotcols{i}},1);
       plot(x, ymean, '-k', 'LineWidth', 1, 'Color', 'k');
         
       % ========== count for display
       disp(['bird' num2str(count)]);
        count = count+1;
    end
    line([0 0], ylim);
    lt_plot_annotation(1, ['n=' num2str(count-1) ' birds'], 'k')
end

linkaxes(hsplots, 'xy');
axis tight


%% ============================ PLOTS (DECODING, FITTING SIGMOID MODEL)
% =========== FITS TO EITHER MEAN OR ALL RAW DATA
fittomean = 0;

lt_figure; hold on;
plotcols = lt_make_plot_colors(length(BRANCHES), 0,0);

ModelOutputs = struct;
A4fits = [];
A4fitsCI = [];

for i=1:max(AllBRANCHSTRUCTnum)
    lt_subplot(4,2, i); hold on;
    
    inds = AllBRANCHSTRUCTnum==i;
    
    x = AllDecodeX(inds,:);
    y = AllDecodeY(inds,:);
    
    if fittomean==1
        y = mean(y,1);
        x = mean(x,1);
    else
        x = x(:);
        y = y(:);
        
    end
    
    % ================= FIT MODEL
    %     modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*x + A(4)))); % logistic
    %     beta0 = [0.5 0.1 50 0.1]';
    modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*(x - A(4))))); % logistic
    beta0 = [0.5 0.1 50 0.01]';
    %     opt = statset('fitnlm');
    %     opt.RobustWgtFun = 'bisquare';
    mdl = fitnlm(x, y, modelfun, beta0);
    
    % - plot actual
    plot(x, y, 'x','Color',plotcols{i});
    
    
    xx = min(x):0.005:max(x);
    plot(xx, mdl.feval(xx), 'o-k');
    
    % -- plot mean
    x = AllDecodeX(inds,:);
    y = AllDecodeY(inds,:);
    ymean = mean(y);
    plot(x, ymean, '-k', 'LineWidth',2);
    
    % =========================== OUTPUT
    ModelOutputs(i).mdl = mdl;
    disp(mdl)
    
    tmp = mdl.coefCI;
    A4fitsCI = [A4fitsCI; tmp(4,:)];
    A4fits = [A4fits; mdl.Coefficients.Estimate(4,1)];
    
end

% === compare mdl fits (specifically A4, which is lateral shift)
lt_subplot(4,2,i+1);
errorbar(1:length(A4fits), A4fits, A4fits-A4fitsCI(:,1), A4fits-A4fitsCI(:,2))


%% ============================ PLOTS (DECODING, FITTING SIGMOID MODEL)
% ===== each neuron/bird is one datapoint


lt_figure; hold on;
plotcols = lt_make_plot_colors(length(BRANCHES), 0,0);

ModelOutputs = struct;
A4fits = [];
A4fitsCI = [];

for i=1:max(AllBRANCHSTRUCTnum)
    lt_subplot(4,2, i); hold on;
    
    inds = AllBRANCHSTRUCTnum==i;
    
    x = AllDecodeX(inds,:);
    y = AllDecodeY(inds,:);
    bnums = AllBirdnum(inds);
    bnums = repmat(bnums, 1, size(y,2)); % make sure each time bin has a bird index
    neurnums = AllNeurnum(inds);
    neurnums = repmat(neurnums, 1, size(y,2));
    
    y = y(:);
    x = x(:);
    bnums = bnums(:);
    neurnums = neurnums(:);
    
    ymean = grpstats(y, {x, bnums, neurnums}, {'mean'});
    xmean = grpstats(x, {x, bnums, neurnums}, {'mean'});
    
    
    % ================= FIT MODEL
    %     modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*x + A(4)))); % logistic
    %     beta0 = [0.5 0.1 50 0.1]';
    modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*(x - A(4))))); % logistic
    beta0 = [0.5 0.1 50 0.01]';
    %     opt = statset('fitnlm');
    %     opt.RobustWgtFun = 'bisquare';
    mdl = fitnlm(xmean, ymean, modelfun, beta0);
    
    % ============================== PLOT
    % ----------- 1) plot actual DAT
    plot(xmean, ymean, 'x','Color',plotcols{i});
    
    % -------------- 2) plot fit
    xx = min(xmean):0.005:max(xmean);
    plot(xx, mdl.feval(xx), 'o-k');
    
%     % -- plot mean
%     ymeantmp = grpstats(ymean, xmean, {'mean'});
%     xmeantmp = unique(xmean);
%     plot(xmeantmp, ymeantmp, '-k', 'LineWidth',2);
    
    % =========================== OUTPUT
    ModelOutputs(i).mdl = mdl;
    disp(mdl)
    
    tmp = mdl.coefCI;
    A4fitsCI = [A4fitsCI; tmp(4,:)];
    A4fits = [A4fits; mdl.Coefficients.Estimate(4,1)];
end

% === compare mdl fits (specifically A4, which is lateral shift)
lt_subplot(4,2,i+1);
errorbar(1:length(A4fits), A4fits, A4fits-A4fitsCI(:,1), A4fits-A4fitsCI(:,2))


    
    %% ====================== FITTING SIGMOID TO DECODE - fit both RA and LMAN at once
    % ===================== NO CROSS VALIDATION
    
    lt_figure; hold on;
    
    % ---- response var
    y = AllDecodeY;
    
    % ----- predictors
    x = AllDecodeX;
    bnums = AllBirdnum;
    bnums = repmat(bnums, 1, size(y,2)); % make sure each time bin has a bird index
    neurnums = AllNeurnum;
    neurnums = repmat(neurnums, 1, size(y,2));
    analynums = AllBRANCHSTRUCTnum;
    analynums = repmat(analynums, 1, size(y,2));
    analynums = analynums-1; % convert to 0 and 1 from 1 and 2;
    
    % ---- get into vectors
    y = y(:);
    x = x(:);
    bnums = bnums(:);
    neurnums = neurnums(:);
    analynums = analynums(:);
    
    
    % --- for each neuron/bird, get mean across branches
    ymean = grpstats(y, {x, bnums, neurnums, analynums}, {'mean'});
    xmean = grpstats(x, {x, bnums, neurnums, analynums}, {'mean'});
    analymean = grpstats(analynums, {x, bnums, neurnums, analynums}, {'mean'});
%     ymean = grpstats(y, {x, analynums}, {'mean'});
%     xmean = grpstats(x, {x, analynums}, {'mean'});
%     analymean = grpstats(analynums, {x, analynums}, {'mean'});
    
    % ---- combine predictors
    X = [xmean analymean];
    
    % ==== FIT MODEL
    % --- full
    modelfun = @(A, x) ((A(1)+A(5)*x(:,2)) + ...
        (A(2)+A(6)*x(:,2))./(1 + exp(-(A(3)+A(7)*x(:,2)).*(x(:,1) - ...
        (A(4)+A(8)*x(:,2)))))); % logistic
    
    beta0 = [0.4 0.1 50 0.01 0.02 0.02 10 0.01]';
   
    % --- reduced
%     modelfun = @(A, x) ((A(1)+A(5)*x(:,2)) + ...
%         (A(2)+A(6)*x(:,2))./(1 + exp(-(A(3)+A(7)*x(:,2)).*(x(:,1) - ...
%         (A(4)))))); % logistic
%     
%     beta0 = [0.5 0.1 50 0.01 0.01 0.01 10]';
%     
    mdl = fitnlm(X, ymean, modelfun, beta0);
    
    
    % ===== PLOT
    % ----------- 1) plot actual DAT
    for j=1:max(X(:,2)+1)
        lt_subplot(2,2,j); hold on;
        indstmp = X(:,2)+1==j;
        plot(X(indstmp,1), ymean(indstmp), 'o', 'Color', plotcols{j});
        
        % ==== also plot mean
        [ymeantmp, ystdtmp] = grpstats(ymean(indstmp), single(X(indstmp,1)), {'mean', 'std'});
        xtmp = unique(single(X(indstmp,1)));
        shadedErrorBar(xtmp, ymeantmp, ystdtmp, {'Color', 'k'},1);
    end
    
    % -------------- 2) plot fit
    for j=1:max(X(:,2)+1)
        lt_subplot(2,2,j); hold on;
        xx = min(X(:,1)):0.005:max(X(:,1));
        
        Xtmp = [xx', (j-1)*ones(size(xx'))];
        
        plot(xx, mdl.feval(Xtmp), 's-k');
    end
    
    % ------------
    disp(mdl)
   
    % ======================= OVERLAY ALL MODELS IN SAME PLOT
    lt_subplot(2,2,max(X(:,2)+1)+1); hold on;
    for j=1:max(X(:,2)+1)
        
        xx = min(X(:,1)):0.005:max(X(:,1));
        
        Xtmp = [xx', (j-1)*ones(size(xx'))];
        
        plot(xx, mdl.feval(Xtmp), 's-k');
        
        % ==== find midpoint of this model
        A = mdl.Coefficients.Estimate;
        xcenter = (A(4) + A(8)*(j-1));
        line([xcenter xcenter], ylim, 'Color', plotcols{j});
        disp(xcenter);
    end
    
    
    
    
%% ====================== FITTING SIGMOID TO DECODE - fit both RA and LMAN at once
% ==================== WITH CROSS VALIDATED
% -------------- compares full model (with lateral shift) to reduced model
% (without shift) by comparing cross validation (10x fold) residuals. 
% 
    cvfold = 5;

%  ------- full model (with lateral shift)
    modelfun1 = @(A, x) ((A(1)+A(5)*x(:,2)) + ...
        (A(2)+A(6)*x(:,2))./(1 + exp(-(A(3)+A(7)*x(:,2)).*(x(:,1) - ...
        (A(4)+A(8)*x(:,2)))))); % logistic
    
    beta1 = [0.5 0.1 50 0.01 0.01 0.01 10 0.01]';

% ----------- reduced model (no lateral shift)
    modelfun2 = @(A, x) ((A(1)+A(5)*x(:,2)) + ...
        (A(2)+A(6)*x(:,2))./(1 + exp(-(A(3)+A(7)*x(:,2)).*(x(:,1) - ...
        (A(4)))))); % logistic
    
    beta2 = [0.5 0.1 50 0.01 0.01 0.01 10]';
    
% ----------- linear model 
    modelfun3 = @(A, x) (A(1)*x(:,1)); % logistic
    
    beta3 = [0.5]';
    

ModelList = {modelfun1, modelfun2, modelfun3};
BetaList = {beta1, beta2, beta3};
ResidualsByModel = {};
for i=1:length(ModelList)
    
    modelfun = ModelList{i};
    beta0 = BetaList{i};

    lt_figure; hold on;
    
    % ---- response var
    y = AllDecodeY;
    
    % ----- predictors
    x = AllDecodeX;
    bnums = AllBirdnum;
    bnums = repmat(bnums, 1, size(y,2)); % make sure each time bin has a bird index
    neurnums = AllNeurnum;
    neurnums = repmat(neurnums, 1, size(y,2));
    analynums = AllBRANCHSTRUCTnum;
    analynums = repmat(analynums, 1, size(y,2)); 
    analynums = analynums-1; % convert to 0 and 1 from 1 and 2;

    
    % -- make unique index for each neuron
    uniqueID = [num2str(AllBirdnum) num2str(AllNeurnum) num2str(AllBRANCHSTRUCTnum)];
    [G] = grp2idx(uniqueID);
    Glevels = unique(G);
    
    % ################## PERFROM CROSS VAL BY PARTITION BY UNIQUE NEURON ID
    cpart = cvpartition(Glevels, 'KFold', cvfold);
    ResidualsAll = [];
    for j=1:cvfold
    
        % =========================== TRAINING SET
        % -- only those that have the unique ID in training set
        uniqueNeuronsToKeep = Glevels(training(cpart, j));
       
        % -- which datapoints are this unique neurons
        indstokeep = ismember(G, uniqueNeuronsToKeep);
        testinds = ~indstokeep;
        
        
        ytrain = y(indstokeep,:);
        xtrain = x(indstokeep, :);
        bnumstrain = bnums(indstokeep, :);
        neurnumstrain = neurnums(indstokeep, :);
        analynumstrain = analynums(indstokeep,:);
        
    % ---- get into vectors
    ytrain = ytrain(:);
    xtrain = xtrain(:);
    bnumstrain = bnumstrain(:);
    neurnumstrain = neurnumstrain(:);
    analynumstrain = analynumstrain(:);
    

    % --- for each neuron/bird, get mean across branches
    ymean = grpstats(ytrain, {xtrain, bnumstrain, neurnumstrain, analynumstrain}, {'mean'});
    xmean = grpstats(xtrain, {xtrain, bnumstrain, neurnumstrain, analynumstrain}, {'mean'});
    analymean = grpstats(analynumstrain, {xtrain, bnumstrain, neurnumstrain, analynumstrain}, {'mean'});
    
    % ---- combine predictors
    X = [xmean analymean];
    
    % ==== FIT MODEL

    mdl = fitnlm(X, ymean, modelfun, beta0);

    
    % =========================== GET PREDICTION RESIDUALS ON TEST SET
        ytest = y(testinds,:);
        xtest = x(testinds, :);
        analynumstest = analynums(testinds,:);
        neurnumstest = neurnums(testinds, :);
        bnumstest = bnums(testinds, :);
        
        ytest = ytest(:);
        xtest = xtest(:);
        analynumstest = analynumstest(:);
        neurnumstest = neurnumstest(:);
        bnumstest = bnumstest(:);
        
        
        % --- for each neuron/bird, get mean across branches
        ymean = grpstats(ytest, {xtest, bnumstest, neurnumstest, analynumstest}, {'mean'});
        xmean = grpstats(xtest, {xtest, bnumstest, neurnumstest, analynumstest}, {'mean'});
        analymean = grpstats(analynumstest, {xtest, bnumstest, neurnumstest, analynumstest}, {'mean'});


        Xtest = [xmean analymean];
        Ypredict = mdl.feval(Xtest);
        
        % ----- visualize prediction
        lt_subplot(5,3,j); hold on;
        title(['cvfold ' num2str(j)]);
        plot(Xtest(Xtest(:,2)==0,1), Ypredict(Xtest(:,2)==0), 'kx', ...
            Xtest(Xtest(:,2)==1,1), Ypredict(Xtest(:,2)==1), 'kx');
        
        plot(Xtest(Xtest(:,2)==0,1), ytest(Xtest(:,2)==0), 'b.', ...
            Xtest(Xtest(:,2)==1,1), ytest(Xtest(:,2)==1), 'r.');
   
        % ############## EXTRACT PREDICTION RESIDUALS
        resid = ymean - Ypredict;
        ResidualsAll = [ResidualsAll; resid];
    end
    
    lt_subplot(5,3,j+1);
    title('cv residuals');
    lt_plot_histogram(ResidualsAll);
    
    % ============== OUTPUT
    ResidualsByModel = [ResidualsByModel ResidualsAll];
end

lt_figure; hold on;
title('residuals (sq)');
ResidualsByModel_abs = ResidualsByModel;
for j=1:length(ResidualsByModel_abs)
    ResidualsByModel_abs{j} = ResidualsByModel_abs{j}.^2;
end
lt_plot_MultDist(ResidualsByModel_abs, [1 2 3]);

%% ============ PLOTS (DIFFERENTIALS)

% ================ 1) MEAN LMAN VS. RA
lt_figure; hold on;
figlegend = {};
plotcols = lt_make_plot_colors(max(AllBRANCHSTRUCTnum), 0, 0);
for i=1:max(AllBRANCHSTRUCTnum);
    for j=1:max(AllAlignNum)
        inds = AllBRANCHSTRUCTnum ==i & AllAlignNum==j;
        
        x = BinMidList;
        Ymat = AllDifferentials(inds, :);
        ymean = nanmean(Ymat, 1);
        ysem =lt_sem(Ymat);
        ID = BRANCHES(i).ID;
        figlegend = [figlegend ID];
        
        errorbar(x, ymean, ysem, '.-', 'Color', plotcols{i});
    end
end
legend(figlegend);
lt_plot_zeroline;



% ================ 1) INDIVIDUAL BIRDS (COLORED BY REGION)
lt_figure; hold on;
figlegend = {};
plotcols = lt_make_plot_colors(max(AllBRANCHSTRUCTnum), 0, 0);
for i=1:max(AllBRANCHSTRUCTnum)
    for j=1:max(AllAlignNum)
        for ii=1:max(AllBirdnum)
            inds = AllBRANCHSTRUCTnum ==i & AllBirdnum==ii ...
                & AllAlignNum==j;
            
            if ~any(inds)
                continue
            end
            
            x = BinMidList;
            Ymat = AllDifferentials(inds, :);
            ymean = nanmean(Ymat, 1);
            ysem =lt_sem(Ymat);
            ID = BRANCHES(i).ID;
            figlegend = [figlegend ID];
            
            %         errorbar(x, ymean, ysem, '.-', 'Color', plotcols{i});
            plot(x, ymean, '.-', 'Color', plotcols{i});
        end
    end
end
legend(figlegend);
lt_plot_zeroline;




% ================ 1) INDIVIDUAL BIRDS (COLORED BY REGION) (SEPARATE
% PLOTS)
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

plotcols = lt_make_plot_colors(max(AllBRANCHSTRUCTnum), 0, 0);
for i=1:max(AllBRANCHSTRUCTnum)
    for j=1:max(AllAlignNum)
        for ii=1:max(AllBirdnum)
            inds = AllBRANCHSTRUCTnum ==i & AllBirdnum==ii & AllAlignNum==j;
            
            if ~any(inds)
                continue
            end
            
            
            x = BinMidList;
            Ymat = AllDifferentials(inds, :);
            ymean = nanmean(Ymat, 1);
            ysem =lt_sem(Ymat);
            ID = BRANCHES(i).ID;
            figlegend = [figlegend ID];
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            birdname = BRANCHES(i).dat.ALLBRANCH.SummaryStruct.birds(ii).birdname;
            title([ID ', ' birdname]);
            hsplots = [hsplots hsplot];
            %         errorbar(x, ymean, ysem, '.-', 'Color', plotcols{i});
            plot(x, ymean, '.-', 'Color', plotcols{i});
            
            lt_plot_zeroline;
            xlim([-0.1 0.05]);
            %         ylim([0.5 1.2])
        end
    end
end
linkaxes(hsplots);



% ================ 1) INDIVIDUAL NEURONS (COLORED BY REGION)
TimeOfMin = [];
BRANCHSTRUCTnum = [];

lt_figure; hold on;
figlegend = {};
plotcols = lt_make_plot_colors(max(AllBRANCHSTRUCTnum), 0, 0);
for i=1:max(AllBRANCHSTRUCTnum)
    for j=1:max(AllAlignNum)
        for ii=1:max(AllBirdnum)
            for iii=1:max(AllNeurnum)
                inds = AllBRANCHSTRUCTnum ==i & AllBirdnum==ii & AllNeurnum==iii ...
                    & AllAlignNum==j;
                
                if ~any(inds)
                    continue
                end
                
                x = BinMidList;
                Ymat = AllDifferentials(inds, :);
                ymean = nanmean(Ymat, 1);
                ysem =lt_sem(Ymat);
                ID = BRANCHES(i).ID;
                figlegend = [figlegend ID];
                
                %         errorbar(x, ymean, ysem, '.-', 'Color', plotcols{i});
                plot(x, ymean, '.-', 'Color', plotcols{i});
                
                % ======================================== COLLECT MINIMUM
                [~, minind] = min(ymean);
                TimeOfMin = [TimeOfMin; x(minind)];
                BRANCHSTRUCTnum = [BRANCHSTRUCTnum; i];
            end
        end
    end
end
lt_plot_zeroline;

% ---- PLOT TIME OF MIN
lt_figure; hold on;
grpstats(TimeOfMin, BRANCHSTRUCTnum, 0.05)






















