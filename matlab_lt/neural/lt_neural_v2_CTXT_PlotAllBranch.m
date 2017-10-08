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

stdbinsize = diff(ALLBRANCH.alignpos(1).ParamsFirstIter.ClassGeneral.frtimewindow(1:2))/2; % divide 2 since can be smaller

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


%% PLOT ACROSS BRANCHES [sep by alignment, stretch or no stretch



figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

hsplots = [];

numalign = length(ALLBRANCH.alignpos);

DATSTRUCT = struct;

Allalign = [];
Allbirdnum = [];
Allbranchnum = [];
Allneuron = [];

NumRemovedDueToThresh = 0;
NumKeptDueToThresh = 0;
for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
    
    % === one for each branch/syl
    Xcell = {};
    Ycell = {};
    Xcontcell = {};
    Ycontcell = {};
    
    Xcell_pos = {};
    Ycell_pos = {};
    
    Xcell_neg = {};
    Ycell_neg = {};
    
    BirdnumAll = [];
    
    Datstruct = struct;
    
    for ii=1:numbirds
        
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        if any(strcmp(birdstoexclude, birdname))
            disp(['skipping ' birdname]);
            continue
        end
        
        numbranch = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        for j=1:numbranch
            numneuron = length(ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron);
            
            for nn=1:numneuron
                
                datneur = ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn);
                
                
                %% if care about location
                if LMANorX==1
                    location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                    % LMAN
                    if ~strcmp(location, 'LMAN')
                        continue
                    end
                elseif LMANorX==2
                    location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                    if ~strcmp(location, 'X')
                        continue
                    end
                end
                
                %% === COLLECT FOR THIS BRANCH/ALIGNMENT
                
                if isempty(datneur.yvals)
                    continue
                end
                if length(datneur.xtimes)==1
                    continue
                end
                
               %% ==== remove anything with repeats?
               if RemoveRepeats ==1
                   isrepeat=0;
                   for k =1:length(datneur.prms_regexpstrlist)
                       thisstr = datneur.prms_regexpstrlist{k};
                       
                       tokensyl = thisstr(strfind(thisstr, '(')+1);
                       presyl = thisstr(strfind(thisstr, '(')-1);
                       
                       if tokensyl == presyl
                           isrepeat =1;
                       end
                       
                       % NOTE: ad hoc, change this to separate section.
                       if presyl =='i'
                           isrepeat=1;
                       end
                   end
                   
                   if isrepeat==1
                       continue
                   end
               end
               
               %% ===== if want to filter by syl/gap duration differences across classes (within context)
               if datneur.DurAnovas.syl_omega>durThreshOmega.syl | ...
                       datneur.DurAnovas.gappre_omega > durThreshOmega.gappre | ...
                       datneur.DurAnovas.gappost_omega > durThreshOmega.gappost
                   NumRemovedDueToThresh = NumRemovedDueToThresh+1;
                   continue
               else
                   NumKeptDueToThresh = NumKeptDueToThresh+1;
               end
                 
               
              
                %% ========================================== DATA
                
               DatN = nan(1,3);
                if strcmp(dattoplot, 'classperform')
                    Xcell = [Xcell datneur.xtimes];
                    Ycell = [Ycell datneur.yvals];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell = [Ycell nanmean(datneur.DprimeAllPairwise,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise,1))./1000;
                        Xcell = [Xcell xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell = [Ycell datneur.DprimeAllPairwise];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise,2));
                        Xcell = [Xcell xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch
                    xtimes = -motifpredur + (1:size([datneur.FR.classnum.frmean],1))./1000;
                    Xcell = [Xcell xtimes];
                    Ycell = [Ycell nanmean([datneur.FR.classnum.frmean],2)];
                    
                end
                % syl contour
                %                 plot(datneur.sylcontours_x, datneur.sylcontours_mean, '-', 'Color', [0.6 0.2 0.2]);
                Xcontcell = [Xcontcell datneur.sylcontours_x];
                Ycontcell = [Ycontcell datneur.sylcontours_mean];
                
                DatN(1) = length(Ycell{end});
                
                
                
                % =================================== POS CONTROL
                if strcmp(dattoplot, 'classperform')
                    Xcell_pos = [Xcell_pos datneur.xtimes];
                    Ycell_pos = [Ycell_pos datneur.yvals_pos];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell_pos = [Ycell_pos nanmean(datneur.DprimeAllPairwise_Pos,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Pos,1))./1000;
                        Xcell_pos = [Xcell_pos xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell_pos = [Ycell_pos datneur.DprimeAllPairwise_Pos];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Pos,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise_Pos,2));
                        Xcell_pos = [Xcell_pos xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch
                    xtimes = -motifpredur + (1:size([datneur.FR_POSCONTR.classnum.frmean],1))./1000;
                    Xcell_pos = [Xcell_pos xtimes];
                    Ycell_pos = [Ycell_pos nanmean([datneur.FR_POSCONTR.classnum.frmean],2)];
                    
                end
                DatN(2) = length(Ycell_pos{end});

                
                % ========================================= NEG CONTROL
                if strcmp(dattoplot, 'classperform')
                    Xcell_neg = [Xcell_neg datneur.xtimes];
                    Ycell_neg = [Ycell_neg datneur.yvals_neg];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell_neg = [Ycell_neg nanmean(datneur.DprimeAllPairwise_Neg,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Neg,1))./1000;
                        Xcell_neg = [Xcell_neg xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell_neg = [Ycell_neg datneur.DprimeAllPairwise_Neg];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Neg,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise_Neg,2));
                        Xcell_neg = [Xcell_neg xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch [IDENTICAL TO DATA]
                    xtimes = -motifpredur + (1:size([datneur.FR.classnum.frmean],1))./1000;
                    Xcell_neg = [Xcell_neg xtimes];
                    Ycell_neg = [Ycell_neg nanmean([datneur.FR.classnum.frmean],2)];
                end
                DatN(3) = length(Ycell_neg{end});

                
                % =================================== SAMPLE SIZE TALLY
                Allalign = [Allalign i];
                Allbirdnum = [Allbirdnum ii];
                Allbranchnum = [Allbranchnum j];
                Allneuron = [Allneuron nn];
                
                BirdnumAll = [BirdnumAll ii];
                
                % ======= confirm that is paired - i.e. each datapoint has dat and both controls
%                 disp(DatN);
                assert(length(unique(DatN))==1, 'dat and controls have diff lengths ...');
            end
        end
    end
    
    assert(length(Ycell) == length(Ycell_pos), 'asfds');
    assert(length(Ycell) == length(Ycell_neg), 'asdf');
    
    % ================== PUT INTO STRUCT
    Datstruct.Dat.Xcell = Xcell;
    Datstruct.Dat.Ycell = Ycell;
    Datstruct.Dat.Xcontcell = Xcontcell;
    Datstruct.Dat.Ycontcell = Ycontcell;
    
    Datstruct.PosContr.Xcell = Xcell_pos;
    Datstruct.PosContr.Ycell = Ycell_pos;
    
    Datstruct.NegContr.Xcell = Xcell_neg;
    Datstruct.NegContr.Ycell = Ycell_neg;
    
    Datstruct.Dat.BirdnumAll = BirdnumAll;
    
    %% ====================== SUBTRACT CONTROLS
    
    for nn = 1:length(Datstruct.Dat.Ycell);
        % ============================ VS. NEG
        % --- confirm that xvalues match
        assert(all(Datstruct.Dat.Xcell{nn} - Datstruct.NegContr.Xcell{nn} < 0.001), 'safasd');
        Datstruct.DatMinusNeg.Xcell{nn} = Datstruct.Dat.Xcell{nn};
        
        % -- get diff
        ydiff = Datstruct.Dat.Ycell{nn} - Datstruct.NegContr.Ycell{nn};
        Datstruct.DatMinusNeg.Ycell{nn} = ydiff;
        
        
        % ============================== VS POS
        % --- confirm that xvalues match
        assert(all(Datstruct.Dat.Xcell{nn} - Datstruct.PosContr.Xcell{nn} < 0.001), 'safasd');
        Datstruct.DatMinusPos.Xcell{nn} = Datstruct.Dat.Xcell{nn};
        
        % -- get diff
        ydiff = Datstruct.Dat.Ycell{nn} - Datstruct.PosContr.Ycell{nn};
        Datstruct.DatMinusPos.Ycell{nn} = ydiff;
    end
    
    %% ============= PLOT MEANS
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
    
    
    %% ========= time warping, align all neur by warping syl contours
    % ============= METHOD 1, linear stretching, based on autocorrelation
    % width between peaks
    
    % --------- Actual dat
    dattype = 'Dat';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'NegContr';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'PosContr';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'DatMinusNeg';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'DatMinusPos';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    %% ==================== PLOT (STRETCHED)
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
    
    %
    %     else
    %     % ====================================
    %     % --- MODIFY
    %     Ycell_WARP = Datstruct.Dat.Ycell_WARP;
    %     Xcell = Datstruct.Dat.Xcell;
    %     Ycontcell_WARP = Datstruct.Dat.Ycontcell_WARP;
    %     Xcontcell = Datstruct.Dat.Xcontcell;
    %
    %     % ---- RUN
    %     Yall = cell2mat(Ycell_WARP);
    %     Xall = cell2mat(Xcell);
    %
    %     % ---- combine all dat
    %     Yall = Yall(:);
    %     Xall = Xall(:);
    %
    %     Ycont = cell2mat(Ycontcell_WARP);
    %     Xcont = cell2mat(Xcontcell);
    %
    %
    %     [Ymean, Ysem] = grpstats(Yall, Xall, {'mean', 'sem'});
    %     X = unique(Xall);
    %
    %     lt_plot(X, Ymean, {'Errors', Ysem});
    %
    %     % --- contour
    %     [Ymean, Ysem] = grpstats(Ycont, Xcont, {'mean', 'sem'});
    %     X = unique(Xcont);
    %
    %     axis tight
    %     Ylim = ylim;
    %
    %     plot(X, Ylim(1)+Ymean.*(Ylim(2)-Ylim(1)), 'r-', 'LineWidth', 2);
    %     end
    
    %%
    
    DATSTRUCT.numalign(i).Datstruct = Datstruct;
    
    
end
linkaxes(hsplots, 'xy')

%%

disp([' ============ REMOVED due to fail syl/gap dur similarity thesrhold: ' ...
    num2str(NumRemovedDueToThresh) '/' num2str(NumRemovedDueToThresh+NumKeptDueToThresh)]);


%% ==== sample size, display

numbirds = length(unique(Allbirdnum));

numneurons = tabulate([num2str(Allbirdnum') num2str(Allneuron')]);
numneurons = size(numneurons,1);

numbranches = tabulate([num2str(Allbirdnum') num2str(Allneuron') num2str(Allbranchnum')]);
numbranches = size(numbranches, 1);

disp([num2str(numbirds) ' birds, ' num2str(numneurons) ' neurons,' num2str(numbranches) ' branches.'])

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

onsetbounds = [-0.02 0.01];

% --- data
dattype = 'Dat';
lt_subplot(3,1,1); hold on;
title('data');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);
indtmp = Xmat>onsetbounds(1) & Xmat<onsetbounds(2)  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

% --- neg
dattype = 'NegContr';
lt_subplot(3,1,2); hold on;
title('neg contr');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);
indtmp = Xmat>onsetbounds(1) & Xmat<onsetbounds(2)  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

% --- pos
dattype = 'PosContr';
lt_subplot(3,1,3); hold on;
title('pos contr');

Xmat = cell2mat(Datstruct.(dattype).Xcell);
Ymat = cell2mat(Datstruct.(dattype).Ycell);
indtmp = Xmat>onsetbounds(1) & Xmat<onsetbounds(2)  & ~isnan(Ymat);
lt_plot_histogram(Ymat(indtmp))
xlim([0 1]);

lt_subtitle('performance, at bin at syl onset');

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
    datfields = fieldnames(Datstruct);
    for i=1:length(datfields)
        dfield = datfields{i};
        indstokeep = Datstruct.Dat.BirdnumAll==birdnum;
        
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

function Datstruct = warp_stretch(Datstruct, dattype)
%  WILL STRETCH RELATIVE TO ACTUAL DAT FOR THAT BRANCH

% -- MODIFY
Ycell = Datstruct.(dattype).Ycell; % will be stretched
Xcell = Datstruct.(dattype).Xcell;

% -------------- HARD PARAMS
stretchtempl = 'Dat';
Ycontcell = Datstruct.(stretchtempl).Ycontcell; % will determing how much to stretch
Xcontcell = Datstruct.(stretchtempl).Xcontcell;

pkwidthtarg = 100; % 120 ms

% ----- RUN
numsamps = length(Ycell);
Ycell_WARP = cell(size(Ycell));
Ycontcell_WARP = cell(size(Ycontcell));

for j=1:numsamps
    
    tmp = xcorr(Ycontcell{j}(~isnan(Ycontcell{j})));
    tmp = tmp(floor(end/2):end);
    [~, pklocs] = findpeaks(double(tmp), 'sortstr', 'descend', 'npeaks', 2, 'minpeakdistance', 20);
    pkwidth = pklocs(2)-pklocs(1);
    
    % ====== stretch based on pkwidth compared to target
    % --- performance
    xtimes = Xcell{j};
    yvals = Ycell{j};
    if length(xtimes)<3
        Ycell_WARP{j} = Ycell{j};
        Ycontcell_WARP{j} = Ycontcell{j};
        continue
    end
    xtimes_new = xtimes.*pkwidthtarg/pkwidth;
    % resample at the original xtimes
    yvals_new = interp1(xtimes_new, yvals, xtimes);
    
    % -- contour
    xtimes_cont = Xcontcell{j};
    ycont = Ycontcell{j};
    xtimes_cont_new = xtimes_cont.*pkwidthtarg/pkwidth;
    ycont_new = interp1(xtimes_cont_new, ycont, xtimes_cont);
    
    if (0)
        if rand<0.1
            lt_figure; hold on;
            
            subplot(311); hold on;
            title('class performance (bk = original)');
            plot(xtimes, yvals, '-ok');
            plot(xtimes_new, yvals, '-r');
            plot(xtimes, yvals_new, '-or');
            
            subplot(312); hold on;
            title('syl contours');
            plot(xtimes_cont, ycont, '-xk');
            plot(xtimes_cont_new, ycont, '-r');
            plot(xtimes_cont, ycont_new, '-xr');
            
            subplot(313);  hold on;
            title('autocor')
            plot(tmp);
            line([pkwidthtarg pkwidthtarg], ylim, 'Color','k');
            lt_plot_annotation(1, ['wdth:' num2str(pkwidth)], 'r')
            pause; close all;
        end
    end
    
    % ==== OUTPUT (overwrite)
    Ycell_WARP{j} = yvals_new;
    Ycontcell_WARP{j} = ycont_new;
end

% ---- dat
Datstruct.(dattype).Ycell_WARP = Ycell_WARP;
% -- template
Datstruct.(stretchtempl).Ycontcell_WARP = Ycontcell_WARP;

end
