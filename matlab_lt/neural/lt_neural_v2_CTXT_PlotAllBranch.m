function lt_neural_v2_CTXT_PlotAllBranch(ALLBRANCH, LMANorX, dattoplot, birdstoexclude)
%% TO DO

% 1) for pos control, get syl contours. (and then can do stretching for all
% pos controls as well)
% 2) all FR bins --- fuck add everything!!!! (from previous struct... fuck
% this is anoying)


%% lt 8/27/17 - plots classifier results by branch. does time warping on classifier results, using mean contours

% LMANorX = 1; % 0, both; 1, LMAN; 2, X

Niter = 3; % shuffles
Nmin = 3; % min sample size;
% Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;

CompareDprimeWohl = 1; % if 1, then overlays Sober and Wohl values
DprimeNegVersion = 'shuff'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
% DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.

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
numalign = length(ALLBRANCH.alignpos);

for i=1:numalign
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    Nall_dat = [];
    Nall_pos = [];
    Nall_neg = [];
    for ii = 1:numbirds
        numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        Maxbranches = max([Maxbranches numbranches]);
        
        for bb = 1:numbranches
            
            numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron);
            
            for nn=1:numneurons
                
                if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).yvals)
                    continue
                end
                
                
                motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
                motifpostdur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpostdur;
                
                
                
                % ================================== CALCULATE RUNNING
                % D-PRIME AND FR
                
                % ------------- ACTUAL DAT
                dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                numclasses = length(dattmp.classnum);
                
                indstokeep = round(1:1000*(motifpredur+motifpostdur-0.005));
                %
                dprimethisbranch = [];
                
                for cc=1:numclasses
                    N1= size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                    
                    % --------------------------- GET FR
                    if N1<Nmin
                        continue
                    end
                    
                    frmean = mean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep, :),2);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).frmean = frmean;
%                     frRunningStdOfMean = lt_running_stats(frmean, round(stdbinsize*1000), 1);
%                     frRunningStdOfMean = frRunningStdOfMean.STD';
%                     numnan = length(frmean) - length(frRunningStdOfMean);
%                     frRunningStdOfMean = [nan(floor(numnan/2),1); frRunningStdOfMean ...
%                         ; nan(ceil(numnan/2),1)];
                    
%                     ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).frRunningStdOfMean ...
%                         = frRunningStdOfMean;
%                     ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).stdbinsize ...
%                         = single(stdbinsize);
                    
                    
                    
                    % -------------------------- DPRIME
                    Nall_dat = [Nall_dat N1];
                    for ccc = cc+1:numclasses
                        N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        % ===== method 2
                        alldat = [dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:) ...
                            dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur(indstokeep,:)];
                        
                        class1mean = nanmean(alldat(:, 1:N1),2);
                        class1var = nanvar(alldat(:, 1:N1)');
                        
                        class2mean = nanmean(alldat(:, N1+1:N1+N2),2);
                        class2var = nanvar(alldat(:, N1+1:N1+N2)');
                        
                        if (0)
                            lt_figure; hold on;
                            shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                            shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                        end
                        
                        
                        dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                        % from Wohlgemuth, Sober, 2010 (they took abs
                        % value)
                        
                        % ================ OUTPUT FOR THIS BRANCH
                        dprimethisbranch = [dprimethisbranch dprime];
                    end
                end
                
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise = dprimethisbranch;
                
                
                % ------------------ NEG CONTROL
                % =========== version 1 (shuffle across contexts)
                if strcmp(DprimeNegVersion, 'shuff');
                    dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    numclasses = length(dattmp.classnum);
                    
                    %                 DprimeAll = [];
                    dprimethisbranch = [];
                    
                    for cc=1:numclasses
                        N1= size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                        Nall_neg = [Nall_neg N1];
                        if N1<Nmin
                            continue
                        end
                        for ccc = cc+1:numclasses
                            N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                            
                            
                            if N1<Nmin || N2<Nmin
                                continue
                            end
                            
                            % -- actual dat
                            alldat = [dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:) ...
                                dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur(indstokeep,:)];
                            
                            % ---- shuffle
                            DprimeShuffAll = [];
                            for nshuf = 1:Niter
                                indshuff = randperm(N1+N2);
                                %                              indshuff = 1:(N1+N2);
                                alldatshuff = alldat(:, indshuff);
                                
                                class1mean = nanmean(alldatshuff(:, 1:N1),2);
                                class1var = nanvar(alldatshuff(:, 1:N1)');
                                
                                class2mean = nanmean(alldatshuff(:, N1+1:N1+N2),2);
                                class2var = nanvar(alldatshuff(:, N1+1:N1+N2)');
                                
                                %                             % --- cut slightly shorter, so no mismatch
                                %                             class1mean = class1mean(indstokeep);
                                %                             class1var = class1var(indstokeep);
                                %                             class2mean = class2mean(indstokeep);
                                %                             class2var = class2var(indstokeep);
                                %
                                if (0)
                                    lt_figure; hold on;
                                    shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                                    shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                                end
                                dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                                % from Wohlgemuth, Sober, 2010 (they took abs
                                % value)
                                
                                DprimeShuffAll = [DprimeShuffAll dprime];
                            end
                            
                            dprimethisbranch = [dprimethisbranch mean(DprimeShuffAll,2)];
                        end
                    end
                elseif strcmp(DprimeNegVersion, 'Wohl')
                    % then within each syl in stereotyped context, shuffle
                    dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    %                     dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    numclasses = length(dattmp.classnum);
                    
                    %                 DprimeAll = [];
                    dprimethisbranch = [];
                    
                    for cc=1:numclasses
                        N1= floor(size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2)./2);
                        N2 = ceil(size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2)./2);
                        
                        Nall_neg = [Nall_neg N1];
                        if N1<Nmin
                            continue
                        end
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        % -- actual dat
                        alldat = dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:);
                        
                        % ---- shuffle
                        DprimeShuffAll = [];
                        for nshuf = 1:Niter
                            indshuff = randperm(N1+N2);
                            %                              indshuff = 1:(N1+N2);
                            alldatshuff = alldat(:, indshuff);
                            
                            class1mean = nanmean(alldatshuff(:, 1:N1),2);
                            class1var = nanvar(alldatshuff(:, 1:N1)');
                            
                            class2mean = nanmean(alldatshuff(:, N1+1:N1+N2),2);
                            class2var = nanvar(alldatshuff(:, N1+1:N1+N2)');
                            
                            %                             % --- cut slightly shorter, so no mismatch
                            %                             class1mean = class1mean(indstokeep);
                            %                             class1var = class1var(indstokeep);
                            %                             class2mean = class2mean(indstokeep);
                            %                             class2var = class2var(indstokeep);
                            %
                            if (0)
                                lt_figure; hold on;
                                shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                                shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                            end
                            dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                            % from Wohlgemuth, Sober, 2010 (they took abs
                            % value)
                            
                            DprimeShuffAll = [DprimeShuffAll dprime];
                        end
                        
                        dprimethisbranch = [dprimethisbranch mean(DprimeShuffAll,2)];
                        
                    end
                    
                end
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Neg = dprimethisbranch;
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Neg_Niter = Niter;
                
                
                % ------------------- POS CONTROL
                dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR;
                numclasses = length(dattmp.classnum);
                dprimethisbranch = [];
                for cc=1:numclasses
                    N1 = size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                    Nall_pos = [Nall_pos N1];
                    if N1<Nmin
                        continue
                    end
                    
                    % --------------------------- GET FR
                    frmean = mean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep, :),2);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).frmean = frmean;
                    
                    for ccc = cc+1:numclasses
                        N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        class1mean = nanmean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                        class1var = nanvar(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur');
                        class2mean = nanmean(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        class2var = nanvar(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur');
                        
                        % --- cut slightly shorter, so no mismatch
                        class1mean = class1mean(indstokeep);
                        class1var = class1var(indstokeep);
                        class2mean = class2mean(indstokeep);
                        class2var = class2var(indstokeep);
                        
                        dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                        % from Wohlgemuth, Sober, 2010 (they took abs
                        % value)
                        
                        dprimethisbranch = [dprimethisbranch dprime];
                    end
                    
                end
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Pos = dprimethisbranch;
                
            end
            
        end
        
    end
    
    %% ==================== COMPARE SAMPLE SIZES FOR DAT AND POS
    lt_figure; hold on;
    ylabel('N');
    xlabel('dat, pos, neg');
    lt_plot_MultDist({Nall_dat, Nall_pos, Nall_neg}, [1 2 3], 1);
    
end



%% PLOT ACROSS BRANCHES [sep by alignment, stretch or no stretch



figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

hsplots = [];

numalign = length(ALLBRANCH.alignpos);

DATSTRUCT = struct;

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
    
    Datstruct = struct;
    
    for ii=1:numbirds
        

        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        if any(strcmp(birdstoexclude, birdname))
            continue
        end
            
        numbranch = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        for j=1:numbranch
            numneuron = length(ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron);
            
            for nn=1:numneuron
                
                datneur = ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn);
                
                %% if care about location
                location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                if LMANorX==1
                    % LMAN
                    if ~strcmp(location, 'LMAN')
                        continue
                    end
                elseif LMANorX==2
                    if ~strcmp(location, 'X')
                        continue
                    end
                end
                
                %% === COLLECT FOR THIS BRANCH/ALIGNMENT
                
                if isempty(datneur.yvals)
                    continue
                end
                
                % ========================================== DATA
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
                
                
                % ================================ SUBTRACT CONTROLS
                
                
            end
        end
    end
    
    
    % ================== PUT INTO STRUCT
    Datstruct.Dat.Xcell = Xcell;
    Datstruct.Dat.Ycell = Ycell;
    Datstruct.Dat.Xcontcell = Xcontcell;
    Datstruct.Dat.Ycontcell = Ycontcell;
    
    Datstruct.PosContr.Xcell = Xcell_pos;
    Datstruct.PosContr.Ycell = Ycell_pos;
    
    Datstruct.NegContr.Xcell = Xcell_neg;
    Datstruct.NegContr.Ycell = Ycell_neg;
    
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


end








%%


function plotMeanTraj(Datstruct, dattype, plotcol, normbymaxmin, warpOn)
% combines all data (i.e. each datapoint grp stats, grouped by x value)

% --- warp
if ~exist('warpOn', 'var')
    warpOn = 0;
end

if ~exist('normbymaxmin','var')
    normbymaxmin=0;
end
CIalpha = 0.01;

if strcmp(dattype, 'sylcontour')
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
        shadedErrorBar(X, Ymean, Ysem, {'Color', plotcol}, 0);
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
