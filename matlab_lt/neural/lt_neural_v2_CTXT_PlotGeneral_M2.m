function lt_neural_v2_CTXT_PlotGeneral_M2(CLASSEScompiled, plotstat)
%% lt 8/15/17 - takes output from lt_neural_v2_CTXT_PlotGeneral_M2 and plots

CIalpha = 0.01; % for CI, for each time bin

%% ---

numanalys = length(CLASSEScompiled.analynum); % actually is frwindow sizes

%%
% ---- collect params across all experiments
tmp.AllFRwindsize = [];
tmp.AllFRtwindow = [];
tmp.FRbinsize = [];
tmp.algnonset = [];
tmp.algnsyl = [];


for i = 1:numanalys % fr window size
    
    numiterations = length(CLASSEScompiled.analynum(i).iterationnum);
    frwindowsize = diff(CLASSEScompiled.analynum(i).iterationnum(1).params.ClassGeneral.frtimewindow);
    
    for ii=1:numiterations % combination of fr bin size and window location
        
        frbinsize = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frbinsize;
        frtimewindow = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frtimewindow;
        
        algnsyl = CLASSEScompiled.analynum(i).iterationnum(ii).params.alignWhichSyl;
        algnonset = CLASSEScompiled.analynum(i).iterationnum(ii).params.alignOnset;
        
        
        tmp.AllFRwindsize = [tmp.AllFRwindsize; frwindowsize];
        tmp.AllFRtwindow = [tmp.AllFRtwindow; frtimewindow];
        tmp.FRbinsize = [tmp.FRbinsize; frbinsize];
        tmp.algnonset = [tmp.algnonset; algnonset];
        tmp.algnsyl = [tmp.algnsyl; algnsyl];
        
    end
end

assert(length(unique(tmp.algnsyl)) ==1, 'problem, diff syl algn experiments included ...');
assert(length(unique(tmp.algnonset)) ==1, 'problem, diff algn onset experiments included ...');

disp(['FR window sizes used: ' num2str(unique(tmp.AllFRwindsize)')])
disp(['FR binsizes used: ' num2str(unique(tmp.FRbinsize)')])

SylAlignNum = unique(tmp.algnsyl);
SylAlignOnset = unique(tmp.algnonset);

clear tmp;

%% PLOT
AllFRwindowsize = [];
AllFRwindow_middle = [];
AllFRbinsize = [];
AllPerformance = [];
AllNumctxts = [];

% NOTE: below collects all syl contours across everything (including all
% branches and numctxts) - should change to at least split by branch type
% (or at least num ctxts);
AllSylContours = [];
AllSylContours_xtimes = [];

AllPosControlPerformance = [];
AllNegControlPerformance = [];

for i = 1:numanalys % fr window size
    
    numiterations = length(CLASSEScompiled.analynum(i).iterationnum);
    frwindowsize = diff(CLASSEScompiled.analynum(i).iterationnum(1).params.ClassGeneral.frtimewindow);
    
    for ii=1:numiterations % combination of fr bin size and window location
        
        
        % ================================== VARIOUS THINGS
        frbinsize = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frbinsize;
        frtimewindow = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frtimewindow;
        windowmid = mean(frtimewindow);
        Numctxts = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllNumCtxts];
        
            accuracy = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllAccuracy];
        % -------------- OLD VERSION ...
        if (0)
            sensitivity = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllSensitivity];
            piyactual = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllPiYActualMean];
            
            
            %         xbins = frtimewindow(1):frbinsize:frtimewindow(end);
            %         xbins = mean([xbins(1:end-1)' xbins(2:end)'], 2);
            %
            
            % === output
            if strcmp(plotstat, 'accuracy')
                AllPerformance = [AllPerformance accuracy];
            elseif strcmp(plotstat, 'sensitivity')
                AllPerformance = [AllPerformance sensitivity];
            elseif strcmp(plotstat, 'pihatmean')
                AllPerformance = [AllPerformance piyactual];
            else
                disp('problem!!!')
            end
            
        end
        
        AllNumctxts = [AllNumctxts Numctxts];
        
        % -- expand to sample size
        AllFRbinsize = [AllFRbinsize ones(size(accuracy))*frbinsize];
        AllFRwindow_middle = [AllFRwindow_middle ones(size(accuracy))*windowmid];
        AllFRwindowsize = [AllFRwindowsize ones(size(accuracy))*frwindowsize];
        
        
        
        
        % ================================================== ACTUAL DAT
        sts = lt_neural_ConfMatStats(...
            {CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllConfMat});
        AllPerformance = [AllPerformance [sts.(plotstat)]];
        
        % ================ CONTROLS
        if isfield(CLASSEScompiled.analynum(i).iterationnum(ii).allbranches, 'POSCONTR_AllConfMat')
        % --- positive
        sts = lt_neural_ConfMatStats(...
            {CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.POSCONTR_AllConfMat});
        AllPosControlPerformance = [AllPosControlPerformance [sts.(plotstat)]];
        
        % -- neg
        sts = lt_neural_ConfMatStats(...
            {CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.NEGCONTR_AllConfMat});
        AllNegControlPerformance = [AllNegControlPerformance [sts.(plotstat)]];
        end
        
        
        % =============================================== SYL CONTOURS
        % --- save mean syl traces
        % for now, just take mean over all branches
        numbranches = length(CLASSEScompiled.analynum(i).iterationnum(ii).allbranches);
        motifpredur = CLASSEScompiled.analynum(i).iterationnum(ii).params.motifpredur;
        
        SylContourstmp = [];
        
        for bb=1:numbranches
            
            sylcontours = CLASSEScompiled.analynum(i).iterationnum(ii).allbranches(bb).AllBranchSylContours;
            
            SylContourstmp = [SylContourstmp; mean(sylcontours,1)];
            
        end
        
        AllSylContours = [AllSylContours mean(SylContourstmp,1)];
        xtimes = (1:size(sylcontours,2))/1000 - motifpredur;
        AllSylContours_xtimes = [AllSylContours_xtimes xtimes];
        
        

    end
end

[SylContourMean, ystd] = grpstats(AllSylContours, AllSylContours_xtimes, {'mean', 'std'});
SylContourMean_times = unique(AllSylContours_xtimes);

%% =================== PLOT
% Separate plots for each: 0) value of numctxts 1) binsize, 2) Fr window size

figcount=1;
subplotrows=7;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

binsizelist = unique(AllFRbinsize);
numctxtslist = unique(AllNumctxts);
frwindowlist = unique(AllFRwindowsize);
hsplots = [];


for fw = frwindowlist
    for k=numctxtslist
        
        for bsize = binsizelist;
            
            inds = find(AllFRbinsize==bsize & AllFRwindowsize==fw & AllNumctxts==k);
            
            if isempty(inds)
                continue
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            xtimes = AllFRwindow_middle(inds);
            
            %  ------ 1) ACTUAL DATA
            yvals = AllPerformance(inds);
            
            plot(xtimes, yvals, 'x', 'Color', [0.6 0.6 0.6]);
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'});
            axis tight
            
            % -- line indicating null value for this number of contexts
            line(xlim, [1/k 1/k], 'Color', 'b');
            % --
            plot(SylContourMean_times, SylContourMean, '-r');
            
            
            % ------ 2) NEGATIVE CONTROL
            if ~isempty(AllNegControlPerformance)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(['[NEG] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
                
                yvals = AllNegControlPerformance(inds);
                plot(xtimes, yvals, 'x', 'Color', [0.6 0.6 0.6]);
                
                [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
                xtimes_grp = unique(xtimes);
                lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'r'});
                axis tight
                
                % -- line indicating null value for this number of contexts
                line(xlim, [1/k 1/k], 'Color', 'b');
                % --
                plot(SylContourMean_times, SylContourMean, '-r');
                
            
            % -------- 3) POSITIVE
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['[POS] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            yvals = AllPosControlPerformance(inds);
            plot(xtimes, yvals, 'x', 'Color', [0.6 0.6 0.6]);
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'b'});
            axis tight
            
            
            % -- line indicating null value for this number of contexts
            line(xlim, [1/k 1/k], 'Color', 'b');
            % --
            plot(SylContourMean_times, SylContourMean, '-r');
            end
            
            
        end
    end
end
linkaxes(hsplots, 'x');


%% ============== SINGLE PLOT FOR EACH CTXTNUM (SHOWING DAT, POS, AND NEG CONTROL)
% NOTE: only keeping branches that have data for dat, neg, and pos
% throwing out others
            if ~isempty(AllNegControlPerformance)

figcount=1;
subplotrows=7;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

binsizelist = unique(AllFRbinsize);
numctxtslist = unique(AllNumctxts);
frwindowlist = unique(AllFRwindowsize);
hsplots = [];


for fw = frwindowlist
    for k=numctxtslist
        
        for bsize = binsizelist;
            
            inds = find(AllFRbinsize==bsize & AllFRwindowsize==fw & AllNumctxts==k ...
                & ~isnan(AllPerformance.*AllNegControlPerformance.*AllPosControlPerformance));
            
            if isempty(inds)
                continue
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            xtimes = AllFRwindow_middle(inds);
            
            %  ------ 1) ACTUAL DATA
            yvals = AllPerformance(inds);
            plotcol = 'k';
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol});

            
            % ------ 2) NEGATIVE CONTROL
            yvals = AllNegControlPerformance(inds);
            plotcol = 'r';
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol});

               
                % -------- 3) POSITIVE
            yvals = AllPosControlPerformance(inds);
            plotcol = 'b';
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol});

            
           
            % ------------------- GENERAL
            % -- line indicating null value for this number of contexts
            axis tight;
            line(xlim, [1/k 1/k], 'Color', 'm');
            % --
            plot(SylContourMean_times, (1/k)*SylContourMean+0.5*1/k, '-r');
            
        end
    end
end
linkaxes(hsplots, 'x');
            end



%% ==== PAIRED COMPARISONS OF DATA, NEG, AND POS CONTROLS
% NOTE: only keeping data that have pairs (checked separately for datvs.pos
% and datvs. neg)
if   ~isempty(AllNegControlPerformance)


figcount=1;
subplotrows=7;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

binsizelist = unique(AllFRbinsize);
numctxtslist = unique(AllNumctxts);
frwindowlist = unique(AllFRwindowsize);
hsplots = [];


for fw = frwindowlist
    for k=numctxtslist
        
        for bsize = binsizelist;
            
            inds = find(AllFRbinsize==bsize & AllFRwindowsize==fw & AllNumctxts==k);
            
            if isempty(inds)
                continue
            end
            
            % ====== DAT VS. NEGATIVE CONTROL            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['[NEG] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            xtimes = AllFRwindow_middle(inds);
            
            % -- dat 
            yvals = AllPerformance(inds);
            yvals_neg = AllNegControlPerformance(inds);
            
            ydiff = yvals(~isnan(yvals_neg)) - yvals_neg(~isnan(yvals_neg));
            xtimes = xtimes(~isnan(yvals_neg));
            
            plot(xtimes, ydiff, 'x', 'Color', [0.6 0.6 0.6]);
           
            % -- mean
            [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
            xtimes_grp = unique(xtimes);
            
            sigvals = yCI(:,2).*yCI(:,1)>0;            % -- for each time bin as if sig
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
            lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', 'r'}); % SIGNIFICANT
            
            axis tight
            lt_plot_zeroline;
            plot(SylContourMean_times, SylContourMean, '-k');

            

            % ========== DAT VS POSITIVE CONTROl
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['[POS] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            xtimes = AllFRwindow_middle(inds);
            
            % -- dat 
            yvals = AllPerformance(inds);
            yvals_con = AllPosControlPerformance(inds);
            
            ydiff = yvals(~isnan(yvals_con)) - yvals_con(~isnan(yvals_con));
            xtimes = xtimes(~isnan(yvals_con));
            
            plot(xtimes, ydiff, 'x', 'Color', [0.6 0.6 0.6]);
           
            % -- mean
            [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
            xtimes_grp = unique(xtimes);
            
            sigvals = yCI(:,2).*yCI(:,1)>0;            % -- for each time bin as if sig
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
            lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', 'b'}); % SIGNIFICANT
            
            axis tight
            lt_plot_zeroline;
           
            plot(SylContourMean_times, SylContourMean, '-k');

            
            
        end
    end
end
linkaxes(hsplots, 'xy');


% =================== COMBINE ACROSS ALL CLASS SIZES (SINCE THIS IS
% DIFFERENCE FROM CONTROLS)
for fw = frwindowlist
    
    for bsize = binsizelist;
        
        inds = find(AllFRbinsize==bsize & AllFRwindowsize==fw);
        
        if isempty(inds)
            continue
        end
        
        % ====== DAT VS. NEGATIVE CONTROL
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[NEG] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', AllNumCtxts']);
        
        xtimes = AllFRwindow_middle(inds);
        
        % -- MODIFY THIS
        yvals_con = AllNegControlPerformance(inds); % CONTROL TYPE
        plotcol = 'r';
        
        % ----- LEAVE BELOW
        yvals = AllPerformance(inds);
        ydiff = yvals(~isnan(yvals_con)) - yvals_con(~isnan(yvals_con));
        xtimes = xtimes(~isnan(yvals_con));
        
        plot(xtimes, ydiff, 'x', 'Color', [0.6 0.6 0.6]);
        
        % -- mean
        [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
        xtimes_grp = unique(xtimes);
        
        sigvals = yCI(:,2).*yCI(:,1)>0;            % -- for each time bin as if sig
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        axis tight
        lt_plot_zeroline;
        
        plot(SylContourMean_times, 0.5*SylContourMean-0.25, '-k');
        % --- SAVE
        ymean_NEG = ymean;
        ysem_NEG = ysem;
        yCI_NEG = yCI;
        xtimes_grp_NEG = xtimes_grp;
        
         % ====== DAT VS. POS CONTROL
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[POS] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', AllNumCtxts']);
        
        xtimes = AllFRwindow_middle(inds);
        
        % -- MODIFY THIS
        yvals_con = AllPosControlPerformance(inds); % CONTROL TYPE
        plotcol = 'b';
        
        % ----- LEAVE BELOW
        yvals = AllPerformance(inds);
        ydiff = yvals(~isnan(yvals_con)) - yvals_con(~isnan(yvals_con));
        xtimes = xtimes(~isnan(yvals_con));
        
        plot(xtimes, ydiff, 'x', 'Color', [0.6 0.6 0.6]);
        
        % -- mean
        [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
        xtimes_grp = unique(xtimes);
        
        sigvals = yCI(:,2).*yCI(:,1)>0;            % -- for each time bin as if sig
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        axis tight
        lt_plot_zeroline;
        
        plot(SylContourMean_times, 0.5*SylContourMean-0.25, '-k');
        % --- SAVE
        ymean_POS = ymean;
        ysem_POS = ysem;
        yCI_POS = yCI;
        xtimes_grp_POS = xtimes_grp;
        
        
        % ================ COMBINE POS AND NEG CONTROLS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[NEG AND POS] wsize:' num2str(fw) ', bsize:' num2str(bsize) ', AllNumCtxts']);

        % --- neg
        plotcol = 'r';
        
        sigvals = yCI_NEG(:,2).*yCI_NEG(:,1)>0;            % -- for each time bin as if sig
        lt_plot(xtimes_grp_NEG, ymean_NEG, {'Errors', ysem_NEG, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp_NEG(sigvals), ymean_NEG(sigvals), {'Errors', ysem_NEG(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        % --- pos
        plotcol = 'b';
        
        sigvals = yCI_POS(:,2).*yCI_POS(:,1)>0;            % -- for each time bin as if sig
        lt_plot(xtimes_grp_POS, ymean_POS, {'Errors', ysem_POS, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp_POS(sigvals), ymean_POS(sigvals), {'Errors', ysem_POS(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        % --
        axis tight
        plot(SylContourMean_times, 0.3*SylContourMean-0.15, '-k');
       lt_plot_zeroline;

        
    end
end
linkaxes(hsplots, 'xy');

end




%%