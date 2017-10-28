function lt_neural_v2_CTXT_PlotAllBranch(ALLBRANCH, LMANorX, dattoplot, birdstoexclude, ...
    durThreshOmega, RemoveRepeats, RemovePrecededByIntro)
%% TO DO

% 1) for pos control, get syl contours. (and then can do stretching for all
% pos controls as well)
% 2) all FR bins --- fuck add everything!!!! (from previous struct... fuck
% this is anoying)


%% lt 8/27/17 - plots classifier results by branch. does time warping on classifier results, using mean contours

% LMANorX = 1; % 0, both; 1, LMAN; 2, X

CompareDprimeWohl = 1; % if 1, then overlays Sober and Wohl values

% stdbinsize = diff(ALLBRANCH.alignpos(1).ParamsFirstIter.ClassGeneral.frtimewindow(1:2))/2; % divide 2 since can be smaller
% stdbinsize = diff(ALLBRANCH.alignpos(1).ParamsFirstIter.alignOnset frtimewindow(1:2))/2; % divide 2 since can be smaller

%%  get max neurons, and branches

Maxneurons = [];
Maxbranches = [];

numalign = length(ALLBRANCH.alignpos);

for i=1:numalign
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    
    for ii = 1:numbirds
        numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        Maxbranches = max([Maxbranches numbranches]);
        
        for bb = 1:numbranches
            
            numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron);
            
            Maxneurons = max([Maxneurons numneurons]);
            
        end
        
    end
end

%% ======== dprime stuff [and mean FR and FR std]

% DPRIME STUFF
Niter = 3; % shuffles
Nmin = 3; % min sample size;
% Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;

% DprimeNegVersion = 'Wohl'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.

ALLBRANCH = lt_neural_v2_CTXT_AllBranchDprime(ALLBRANCH, Nmin, Niter, ...
    DprimeNegVersion);


%% ======== COLLECT DATA ACROSS BRANCHES [sep by alignment, stretch or no stretch

[DATSTRUCT, Diagnos] = lt_neural_v2_CTXT_BranchToDatstruct(ALLBRANCH, birdstoexclude, ...
    LMANorX, RemoveRepeats, durThreshOmega, dattoplot);


% ============ Diagnos

numbirds = length(unique(Diagnos.Allbirdnum));

numneurons = tabulate([num2str(Diagnos.Allbirdnum') num2str(Diagnos.Allneuron')]);
numneurons = size(numneurons,1);

numbranches = tabulate([num2str(Diagnos.Allbirdnum') num2str(Diagnos.Allneuron') ...
    num2str(Diagnos.Allbranchnum')]);
numbranches = size(numbranches, 1);

disp([num2str(numbirds) ' birds, ' num2str(numneurons) ' neurons,' num2str(numbranches) ' branches.'])


disp([' ============ REMOVED due to fail syl/gap dur similarity thesrhold: ' ...
    num2str(Diagnos.NumRemovedDueToThresh) '/' ...
    num2str(Diagnos.NumRemovedDueToThresh+Diagnos.NumKeptDueToThresh)]);


%% 

Allbirdnum = Diagnos.Allbirdnum;

%% ============ PLOT MEAN (STRETCHED AND UNSTRETCHED)


figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

hsplots = [];


for i=1:numalign
    
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    % ######################################################### UNSTRETCHED
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BFORE STRETCH] algnsyl' num2str(alignsyl) ', onset' num2str(alignons)]);
    
    
    % ========================================= DAT
    dattype = 'Dat';
    plotcol = 'k';
    plotMeanTraj(Datstruct, dattype, plotcol)
    
    % === NEG CONTR
    dattype = 'NegContr';
    plotcol = 'r';
    plotMeanTraj(Datstruct, dattype, plotcol)
    
    % == POS CONTR
    dattype = 'PosContr';
    plotcol = 'b';
    plotMeanTraj(Datstruct, dattype, plotcol)
    
    % -- syl contor
    dattype = 'sylcontour';
    plotMeanTraj(Datstruct, dattype)
    
    % ========= COMPARE TO SOBER, WOHLGEMUTH IF IS DPRIME
    if CompareDprimeWohl==1 & strcmp(dattoplot, 'dprime')
        
        datmean_conv = [0.228];
        datmean_div = [0.189];
        negmean = [0.148];
        posmean = [0.509];
        
        % --- place line at premotor window, lasting 70ms
        xline = [-0.025 0.05];
        line(xline, [negmean negmean], 'Color', 'r', 'LineWidth', 3);
        line(xline, [posmean posmean], 'Color', 'b', 'LineWidth', 3);
        line(xline, [datmean_conv datmean_conv], 'Color', 'k', 'LineWidth', 3);
        line(xline, [datmean_div datmean_div], 'Color', 'k', 'LineWidth', 3);
        
        ylim([0 0.7]);
        xlim([-0.15 0.15]);
    end
    
    
    
    % ##################################################### PLOT (STRETCHED)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[AFTER STRETCH] algnsyl' num2str(alignsyl) ', onset' num2str(alignons)]);
    hsplots = [hsplots hsplot];
    
    % =====================================
    % ----- Dat
    dattype = 'Dat';
    plotcol = 'k';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    % ---- Neg
    dattype = 'NegContr';
    plotcol = 'r';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    
    % ---- Pos (skip for now, need to stretch relative to itself)
    dattype = 'PosContr';
    plotcol = 'b';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    
    % -- contour
    dattype = 'sylcontour';
    plotcol = 'k';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
end
linkaxes(hsplots, 'xy')




%% ================= PLOT MEANS (SUBTRACTING CONTROLS, AND DOING STATS, PAIRWISE DIFFS)
numalign = length(ALLBRANCH.alignpos);
figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    % ====================== NO STRETCH
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[NO STRETCH] algnsyl' num2str(alignsyl) ', onset' num2str(alignons)]);
    ylabel('rel control');
    
    % ===  minus neg
    dattype = 'DatMinusNeg';
    plotcol = 'r';
    plotMeanTraj(Datstruct, dattype, plotcol)
    
    dattype = 'DatMinusPos';
    plotcol = 'b';
    plotMeanTraj(Datstruct, dattype, plotcol)
    
    % -- syl contor
    dattype = 'sylcontour';
    plotMeanTraj(Datstruct, dattype)
    
    
    % ===================== WITH STRETCH
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[STRETCH] algnsyl' num2str(alignsyl) ', onset' num2str(alignons)]);
    ylabel('rel control');
    
    % ===  minus neg
    dattype = 'DatMinusNeg';
    plotcol = 'r';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    dattype = 'DatMinusPos';
    plotcol = 'b';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    % -- syl contor
    dattype = 'sylcontour';
    warpOn=1;
    plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn);
    
    
end

linkaxes(hsplots, 'xy');



%% ======================== PLOT EACH BIRD SEPARATELY [RAW, STRETCHED]

Numbirds = max(Allbirdnum);
numalign = length(ALLBRANCH.alignpos);

figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    for ii=1:Numbirds
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        
        % ====================== NO STRETCH
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' birdname ']' num2str(alignsyl) ', onset' num2str(alignons)]);
        
        % ----- Dat
        dattype = 'Dat';
        plotcol = 'k';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        % ---- Neg
        dattype = 'NegContr';
        plotcol = 'r';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % ---- Pos (skip for now, need to stretch relative to itself)
        dattype = 'PosContr';
        plotcol = 'b';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % -- contour
        dattype = 'sylcontour';
        plotcol = 'k';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
    end
end

linkaxes(hsplots, 'xy');
lt_subtitle('stretched')

%% ======================== PLOT EACH BIRD SEPARATELY [RAW, NOT STRETCHED]

Numbirds = max(Allbirdnum);
numalign = length(ALLBRANCH.alignpos);

figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    for ii=1:Numbirds
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        
        % ====================== NO STRETCH
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' birdname ']' num2str(alignsyl) ', onset' num2str(alignons)]);
        
        % ----- Dat
        dattype = 'Dat';
        plotcol = 'k';
        warpOn=0;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        % ---- Neg
        dattype = 'NegContr';
        plotcol = 'r';
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % ---- Pos (skip for now, need to stretch relative to itself)
        dattype = 'PosContr';
        plotcol = 'b';
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % -- contour
        dattype = 'sylcontour';
        plotcol = 'k';
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
    end
end

linkaxes(hsplots, 'xy');
lt_subtitle('not strethed')


%% ======================== PLOT EACH BIRD SEPARATELY [MINUS CONTROL, STRETCHED]

Numbirds = max(Allbirdnum);
numalign = length(ALLBRANCH.alignpos);

figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    for ii=1:Numbirds
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        
        % ====================== NO STRETCH
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' birdname ']' num2str(alignsyl) ', onset' num2str(alignons)]);
        
        % ----- Dat minus neg
        dattype = 'DatMinusNeg';
        plotcol = 'r';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % ---- Dat minus pos
        dattype = 'DatMinusPos';
        plotcol = 'b';
        warpOn=1;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
        
        % -- contour
        dattype = 'sylcontour';
        plotcol = 'k';
        warpOn=0;
        plotMeanTraj(Datstruct, dattype, plotcol, 0, warpOn, ii);
        
    end
end

linkaxes(hsplots, 'xy');


%% =========================== PLOT DISTRIBUTIONS
lt_figure; hold on;

onsetbounds = [-0.001 0.001];

% --- data
dattype = 'Dat';
lt_subplot(3,1,1); hold on;
title('data');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);

tmp = min(abs(Xmat));
indtmp = Xmat>onsetbounds(1)+tmp & Xmat<onsetbounds(2)+tmp  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

% --- neg
dattype = 'NegContr';
lt_subplot(3,1,2); hold on;
title('neg contr');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);
tmp = min(abs(Xmat));
indtmp = Xmat>onsetbounds(1)+tmp & Xmat<onsetbounds(2)+tmp  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

% --- pos
dattype = 'PosContr';
lt_subplot(3,1,3); hold on;
title('pos contr');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);
tmp = min(abs(Xmat));
indtmp = Xmat>onsetbounds(1)+tmp & Xmat<onsetbounds(2)+tmp  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

lt_subtitle('performance, at bin at syl onset');


%% ============================= TIMING OF DECREASE IN CONTEXTUAL INFORMATION

% ==================== METHOD 1, fit sigmoid to mean
% ---- 1) COLLECT ALL (X,Y)
for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    Datstruct = DATSTRUCT.numalign(i).Datstruct;
    
    % ===================== collect all datapoints across all birds
    Xall = cell2mat(Datstruct.Dat.Xcell);
    Xall = Xall(:);
    Yall = cell2mat(Datstruct.Dat.Ycell);
    Yall = Yall(:);
    
    
    %     lt_figure; hold on;
    %     plot(Xall, Yall, 'x');
    
    % ----- fit sigmoid (individual points)
    %     modelfun = @(A, x) (A(1) + A(2)./(1 + A(4).*exp(-A(3)*x))); % logistic
    modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*x + A(4)))); % logistic
    beta0 = [0.5 0.1 50 0.1]';
    mdl = fitnlm(Xall, Yall, modelfun, beta0);
    
    % -- plot
    lt_figure; hold on;
    plot(Xall, Yall, 'x')
    x = min(Xall):0.001:max(Xall); plot(x, mdl.feval(x), 'o-');
    
    
    % ----- fit sigmoid (individual points, + var for neur/branch (i.e. offset))
    %     modelfun = @(A, x) (A(1)./(1 + exp(-A(2)*x))); % logistic
    %     beta0 = [0 0]';
    %     mdl = fitnlm(Xall, Yall, modelfun, beta0);
    %
    %     % -- plot
    %     lt_figure; hold on; x = -10:0.01:10; plot(x, mdl.feval(x), 'o-');
    %
    
    % --- fit sigmoid (on the mean)
    [Ymean] = grpstats(Yall, Xall, {'mean'});
    X = unique(Xall);
    %         Ymean = Ymean-mean(Ymean); %  Y vals
    
    %         modelfun = @(A, x) (A(1)./(1 + exp(-A(2)*x))); % logistic
    modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*x + A(4)))); % logistic
    %     beta0 = [0.01 0.02 1]';
    beta0 = [0.5 0.1 50 0.1]';
    %     opt = statset('fitnlm');
    %     opt.RobustWgtFun = 'bisquare';
    mdl = fitnlm(X, Ymean, modelfun, beta0);
    
    lt_figure; hold on;
    plot(X, Ymean, 'xk')
    x = min(X):0.001:max(X);
    plot(x, mdl.feval(x), 'o-');
    %     x = -1:0.001:1; plot(x, mdl.feval(x), 'o-');
    
    
    
    % ===== debug - playing around with sigmoid function
    if (0)
        lt_figure; hold on; plot(X, Ymean, 'o');
        
        %         modelfun = @(A, x) (A(1) + A(2)./(1 + A(4).*exp(-A(3)*x))); % logistic
%         modelfun = @(A, x) (A(1) + A(2)./(1 + exp((-A(3)*x) + A(4)))); % logistic
     modelfun = @(A, x) (A(1) + A(2)./(1 + exp(-A(3)*(x - A(4))))); % logistic
       
        x = -0.1:0.01:0.1;
        
        for a = [0]
            A = [0.6 -0.12 80 a];
            
            plot(x, modelfun(A, x), '-o')
        end
    end
    
    %% fitting using logistic regression
    if (0)
        %     Ymean = (Ymean - 0.42)./0.04;
        Yall = (Yall - 0.42)./0.04;
        mdl = fitglm(Xall, Yall, 'Distribution', 'normal', 'Link', 'logit');
        
        
        lt_figure; hold on;
        plot(Xall, Yall, 'xk')
        %     x = min(X):0.001:max(X);
        x =-1:0.01:1;
        plot(x, mdl.feval(x), 'o-');
    end
    
end





end








%%


function plotMeanTraj(Datstruct, dattype, plotcol, normbymaxmin, warpOn, ...
    birdnum)
%% combines all data (i.e. each datapoint grp stats, grouped by x value)

%% =========== pare down Datstruct for specific bird
% leave exmpty [] to ignore

if ~exist('birdnum', 'var')
    birdnum = [];
end

if ~isempty(birdnum)
    indstokeep = Datstruct.Dat.BirdnumAll==birdnum;
    
    datfields = fieldnames(Datstruct);
    for i=1:length(datfields)
        dfield = datfields{i};
        
        Datstruct.(dfield) = lt_structure_subsample_all_fields(Datstruct.(dfield), indstokeep);
        %         disp(indstokeep);
    end
end



%%
% --- warp
if ~exist('warpOn', 'var')
    warpOn = 0;
end

if ~exist('normbymaxmin','var')
    normbymaxmin=0;
end
CIalpha = 0.05;

if strcmp(dattype, 'sylcontour')
    % ############################################## SYL CONTOUR
    Xcontcell = Datstruct.Dat.Xcontcell;
    
    if warpOn==1
        Ycontcell = Datstruct.Dat.Ycontcell_WARP;
    else
        Ycontcell = Datstruct.Dat.Ycontcell;
    end
    
    Ycont = cell2mat(Ycontcell);
    Xcont = cell2mat(Xcontcell);
    
    % --- contour
    [Ymean, Ysem] = grpstats(Ycont, Xcont, {'mean', 'sem'});
    X = unique(Xcont);
    
    axis tight
    Ylim = ylim;
    
    plot(X, Ylim(1)+Ymean.*(Ylim(2)-Ylim(1)), '-', 'LineWidth', 1, 'Color', [0.7 0.2 0.2]);
else
    % -- MODIFY
    if warpOn==1
        Ycell = Datstruct.(dattype).Ycell_WARP;
    else
        Ycell =     Datstruct.(dattype).Ycell;
    end
    Xcell = Datstruct.(dattype).Xcell;
    
    % -- RUN
    Yall = cell2mat(Ycell);
    Xall = cell2mat(Xcell);
    
    % ---- combine all dat
    Yall = Yall(:);
    Xall = Xall(:);
    
    % ###################################################### DAT MINUS CONTROL
    if strcmp(dattype, 'DatMinusNeg')==1 | strcmp(dattype, 'DatMinusPos')==1
        % then calculate CI, for each time bin
        
        [Ymean, Ysem, yCI] = grpstats(Yall, Xall, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
        X = unique(Xall);
        
        sigvals = yCI(:,2).*yCI(:,1)>0;            % -- for each time bin as if sig
        
        % -- plot nonsig
        lt_plot(X(~sigvals), Ymean(~sigvals), {'Errors', Ysem(~sigvals), ...
            'Color', 'k'});
        
        % -- plot sig
        lt_plot(X(sigvals), Ymean(sigvals), {'Errors', Ysem(sigvals), ...
            'Color', plotcol});
        lt_plot_zeroline;
        
    else
        % ############################################## PERFORMANCE
        [Ymean, Ysem] = grpstats(Yall, Xall, {'mean', 'sem'});
        X = unique(Xall);
        
        % ------ norm by max and min?
        if normbymaxmin==1
            % -- ad hoc, for comparison to Wohlgemuth, sober
            % -- norm by diff between 50ms pre and post onset
            windowearly = [-0.06 -0.04];
            indstmp = X>= windowearly(1) & X<=windowearly(2);
            meanearly = mean(Ymean(indstmp));
            
            windowlate = [0.04 0.06];
            indstmp = X>=windowlate(1) & X<=windowlate(2);
            meanlate = mean(Ymean(indstmp));
            
            % --- scale
            Ymean = (Ymean-meanlate)./(abs(meanearly-meanlate));
            Ysem = Ysem./abs(meanearly - meanlate);
        end
        
        
        %     lt_plot(X, Ymean, {'Errors', Ysem, 'Color', plotcol});
        shadedErrorBar(X, Ymean, Ysem, {'Color', plotcol}, 1);
    end
    
    
end
end

% NOTE: moved to lt_neural_v2_CTXT_BranchToDatstruct.m
% function Datstruct = warp_stretch(Datstruct, dattype)
% %  WILL STRETCH RELATIVE TO ACTUAL DAT FOR THAT BRANCH
% 
% % -- MODIFY
% Ycell = Datstruct.(dattype).Ycell; % will be stretched
% Xcell = Datstruct.(dattype).Xcell;
% 
% % -------------- HARD PARAMS
% stretchtempl = 'Dat';
% Ycontcell = Datstruct.(stretchtempl).Ycontcell; % will determing how much to stretch
% Xcontcell = Datstruct.(stretchtempl).Xcontcell;
% 
% pkwidthtarg = 100; % 120 ms
% 
% % ----- RUN
% numsamps = length(Ycell);
% Ycell_WARP = cell(size(Ycell));
% Ycontcell_WARP = cell(size(Ycontcell));
% 
% for j=1:numsamps
%     
%     tmp = xcorr(Ycontcell{j}(~isnan(Ycontcell{j})));
%     tmp = tmp(floor(end/2):end);
%     [~, pklocs] = findpeaks(double(tmp), 'sortstr', 'descend', 'npeaks', 2, 'minpeakdistance', 20);
%     pkwidth = pklocs(2)-pklocs(1);
%     
%     % ====== stretch based on pkwidth compared to target
%     % --- performance
%     xtimes = Xcell{j};
%     yvals = Ycell{j};
%     if length(xtimes)<3
%         Ycell_WARP{j} = Ycell{j};
%         Ycontcell_WARP{j} = Ycontcell{j};
%         continue
%     end
%     xtimes_new = xtimes.*pkwidthtarg/pkwidth;
%     % resample at the original xtimes
%     yvals_new = interp1(xtimes_new, yvals, xtimes);
%     
%     % -- contour
%     xtimes_cont = Xcontcell{j};
%     ycont = Ycontcell{j};
%     xtimes_cont_new = xtimes_cont.*pkwidthtarg/pkwidth;
%     ycont_new = interp1(xtimes_cont_new, ycont, xtimes_cont);
%     
%     if (0)
%         if rand<0.1
%             lt_figure; hold on;
%             
%             subplot(311); hold on;
%             title('class performance (bk = original)');
%             plot(xtimes, yvals, '-ok');
%             plot(xtimes_new, yvals, '-r');
%             plot(xtimes, yvals_new, '-or');
%             
%             subplot(312); hold on;
%             title('syl contours');
%             plot(xtimes_cont, ycont, '-xk');
%             plot(xtimes_cont_new, ycont, '-r');
%             plot(xtimes_cont, ycont_new, '-xr');
%             
%             subplot(313);  hold on;
%             title('autocor')
%             plot(tmp);
%             line([pkwidthtarg pkwidthtarg], ylim, 'Color','k');
%             lt_plot_annotation(1, ['wdth:' num2str(pkwidth)], 'r')
%             pause; close all;
%         end
%     end
%     
%     % ==== OUTPUT (overwrite)
%     Ycell_WARP{j} = yvals_new;
%     Ycontcell_WARP{j} = ycont_new;
% end
% 
% % ---- dat
% Datstruct.(dattype).Ycell_WARP = Ycell_WARP;
% % -- template
% Datstruct.(stretchtempl).Ycontcell_WARP = Ycontcell_WARP;
% 
% end
