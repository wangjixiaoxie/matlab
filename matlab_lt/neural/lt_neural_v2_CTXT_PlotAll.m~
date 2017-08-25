function ALLDAT = lt_neural_v2_CTXT_PlotAll(strtype, plotstat)
%% lt 8/24/17 - plot all classification results for a given string type
% NOTE: need to have run lt_neural_v2_CTXT_PlotGeneral_M for all Results
% already

CIalpha = 0.01;


%% go load all compiled data structs

savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

cd(savedir)

listoffiles = dir(['CLASSEScompiled_' strtype '_*']);

assert(~isempty(listoffiles), 'No data ...');

count =0;
ALLDAT = struct;
ALLDAT.strtype = strtype;

for i=1:length(listoffiles)
    
    tmp = load(listoffiles(i).name);
    
    uscores = strfind(listoffiles(i).name, '_');
    
    assert(uscores(3)-uscores(2)==15, 'problem, assumes single digits ...');
    
    alignsyl = str2num(listoffiles(i).name(uscores(2)+8));
    alignonset = str2num(listoffiles(i).name(uscores(3)-1));
    
    numanalys = length(tmp.CLASSEScompiled.analynum);
    for ii=1:numanalys
        
        count = count+1;
        
        ALLDAT.analynum(count).analydate = tmp.CLASSEScompiled.analynum(ii).analydate;
        ALLDAT.analynum(count).savenote = tmp.CLASSEScompiled.analynum(ii).savenote;
        ALLDAT.analynum(count).iterationnum = tmp.CLASSEScompiled.analynum(ii).iterationnum;
        
        ALLDAT.analynum(count).alignsyl = alignsyl;
        ALLDAT.analynum(count).alignonset = alignonset;
        
    end
end


%% ================ EXTRACT VECTORS OF DATA (ALL COMBINED)

% ============== METHOD 1 - NO PRESETING SIZE
numanalys = length(ALLDAT.analynum);

AllPerformance = [];
AllPerformancePOS = [];
AllPerformanceNEG = [];

AllWindowMidT = [];
AllNumCtxts = [];

AllAnalynum = [];
AllFRwindowsize = [];
AllFRbinsize = [];

AllAlignSylNum = [];
AllAlignSylOnset = [];

for i=1:numanalys
    
    frwindowsize = diff(ALLDAT.analynum(i).iterationnum(1).params.ClassGeneral.frtimewindow);
    ALLDAT.analynum(i).alignsyl;
    
    
    numiterations = length(ALLDAT.analynum(i).iterationnum);
    
    for ii=1:numiterations
        
        % ============================
        frbinsize = ALLDAT.analynum(i).iterationnum(ii).params.ClassGeneral.frbinsize;
        frtimewindow = ALLDAT.analynum(i).iterationnum(ii).params.ClassGeneral.frtimewindow;
        windowmid = mean(frtimewindow);
        numctxts = [ALLDAT.analynum(i).iterationnum(ii).allbranches.AllNumCtxts];
        
        
        % ===============================
        % ---- ACTUAL DAT
        stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.AllConfMat});
        AllPerformance = [AllPerformance [stats.(plotstat)]];
        
        
        if isfield(ALLDAT.analynum(i).iterationnum(ii).allbranches, 'NEGCONTR_AllConfMat')
            % ---- POS CONTR
            stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.POSCONTR_AllConfMat});
            AllPerformancePOS = [AllPerformancePOS [stats.(plotstat)]];
            
            % ---- NEG CONTR
            stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.NEGCONTR_AllConfMat});
            AllPerformanceNEG = [AllPerformanceNEG [stats.(plotstat)]];
            
        else
            AllPerformancePOS = [AllPerformancePOS nan(1, length(stats))];
            AllPerformanceNEG = [AllPerformanceNEG nan(1, length(stats))];
        end
        
        % ----- OTHER THINGS
        AllFRwindowsize = [AllFRwindowsize frwindowsize*ones(1,length(stats))];
        AllAnalynum = [AllAnalynum i*ones(1,length(stats))];
        AllNumCtxts = [AllNumCtxts numctxts];
        AllFRbinsize = [AllFRbinsize frbinsize*ones(1,length(stats))];
        AllWindowMidT = [AllWindowMidT windowmid*ones(1,length(stats))];
        
        AllAlignSylNum = [AllAlignSylNum ...
            ALLDAT.analynum(i).alignsyl*ones(1,length(stats))];
        
        AllAlignSylOnset = [AllAlignSylOnset ...
            ALLDAT.analynum(i).alignonset*ones(1,length(stats))];
        
        % ----------- SYL CONTOURS
        
        
    end
end

%
%
% % ============== METHOD 2
% tic
% % === quick count of datasize - to preset vectors
% numanalys = length(ALLDAT.analynum);
% veclength = 0;
% for i=1:numanalys
%     numiterations = length(ALLDAT.analynum(i).iterationnum);
%
%     for ii=1:numiterations
%
%         % ---- ACTUAL DAT
%         veclength = veclength + length(ALLDAT.analynum(i).iterationnum(ii).allbranches);
%
%     end
% end
%
%
% % ==== extract dat
% numanalys = length(ALLDAT.analynum);
%
% AllPerformance = nan(veclength,1);
% AllPerformancePOS = nan(veclength,1);
% AllPerformanceNEG = nan(veclength,1);
%
% counter = 1;
% for i=1:numanalys
%
%     frwindowsize = [];
%     frbinsize = [];
%
%     numiterations = length(ALLDAT.analynum(i).iterationnum);
%
%     for ii=1:numiterations
%
%         % ---- ACTUAL DAT
%         stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.AllConfMat});
%         AllPerformance(counter:counter+length(stats)-1) = [stats.(plotstat)];
%
%
%         if isfield(ALLDAT.analynum(i).iterationnum(ii).allbranches, 'NEGCONTR_AllConfMat')
%         % ---- POS CONTR
%         stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.POSCONTR_AllConfMat});
%         AllPerformancePOS(counter:counter+length(stats)-1) = [stats.(plotstat)];
%
%         % ---- NEG CONTR
%         stats = lt_neural_ConfMatStats({ALLDAT.analynum(i).iterationnum(ii).allbranches.NEGCONTR_AllConfMat});
%         AllPerformanceNEG(counter:counter+length(stats)-1) = [stats.(plotstat)];
%
%         end
%
%         counter = counter+length(ALLDAT.analynum(i).iterationnum(ii).allbranches);
%
%     end
% end
% toc


%% =============== PLOT 1) Dat, Pos, Neg (means)

% =============== SEPARATE PLOTS FOR EACH COMBINATION OF ALIGNSYL AND ONSET
% [DAT VS. NEG], [DAT VS POS]
lt_figure; hold on;
NumAnalys = length(unique(AllAnalynum));
numrows = NumAnalys;

for i=1:NumAnalys
    
    inds = AllAnalynum ==i;
    
    alignsylnum = unique(AllAlignSylNum(inds));
    alignsylonset = unique(AllAlignSylOnset(inds));
    windsize = unique(AllFRwindowsize(inds));
    binsize = unique(AllFRbinsize(inds));
    
    assert(numel([alignsylnum alignsylonset windsize windsize]) ==4, 'should all be unique...');
    
    lt_subplot(numrows,1,i); hold on;
    title({['alignsyl ' num2str(alignsylnum) ',onset ' num2str(alignsylonset)], ...
        ['wind ' num2str(windsize) ', bin ' num2str(binsize)]});
    ylabel([plotstat ' (rel control)']);
    
    % ---- x values
    xtimes = AllWindowMidT(inds);
    yvals = AllPerformance(inds);
    
    % ====== DAT VS. NEG
    % - modify
    yvals_con = AllPerformanceNEG(inds);
    plotcol = 'b';
    
    % -----
    ydiff = yvals - yvals_con;
    [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
    xtimes_grp = unique(xtimes);
    
    sigvals = yCI(:,2).*yCI(:,1) >0;
    lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
    lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
    
    
    % ====== DAT VS. POS
    % - modify
    yvals_con = AllPerformancePOS(inds);
    plotcol = 'r';
    
    % -----
    ydiff = yvals - yvals_con;
    [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
    xtimes_grp = unique(xtimes);
    
    sigvals = yCI(:,2).*yCI(:,1) >0;
    lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
    lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
    
    
    
    % =========== syl contour
    % - find iteration that is closest to the end
    numiterations = length(ALLDAT.analynum(i).iterationnum);
    tmpvec = [];
    for j=1:numiterations
        tmpvec = [tmpvec ALLDAT.analynum(i).iterationnum(j).params.ClassGeneral.frtimewindow(2)];
    end
    [~, indtmp] = max(tmpvec);
    
    numbranches = length(ALLDAT.analynum(i).iterationnum(indtmp).allbranches);
    sylcontours = [];
    for j=1:numbranches
        sylcontours = [sylcontours; ...
            mean([ALLDAT.analynum(i).iterationnum(indtmp).allbranches(j).AllBranchSylContours],1)]; % for each branch, get mean, then append to overal matrix
    end
    sylcontours_mean = mean(sylcontours);
    sylcontours_std = std(sylcontours, 0, 1);
    motif_predur = ALLDAT.analynum(i).iterationnum(indtmp).params.motifpredur;
    sylcontours_x = -motif_predur + (1:length(sylcontours_mean))./1000;
    
    shadedErrorBar(sylcontours_x, -0.2+sylcontours_mean./3, sylcontours_std./3, {'Color', 'k'},1);
    
    % =========== general
    axis tight
    lt_plot_zeroline;
    
end
lt_subtitle(['strtype: ' strtype])

%% ============= PLOT 2) CLASS PERFORMANCE (DIFFERENCE) - SEPARATED BY NUM CLASSES

if (0) % 

NumCtxts = unique(AllNumCtxts);

for k=NumCtxts
    
    lt_figure; hold on;
    NumAnalys = length(unique(AllAnalynum));
    numrows = NumAnalys;
    
    for i=1:NumAnalys
        
        inds = AllAnalynum ==i & AllNumCtxts==k;
        
        alignsylnum = unique(AllAlignSylNum(inds));
        alignsylonset = unique(AllAlignSylOnset(inds));
        windsize = unique(AllFRwindowsize(inds));
        binsize = unique(AllFRbinsize(inds));
        
        assert(numel([alignsylnum alignsylonset windsize windsize]) ==4, 'should all be unique...');
        
        lt_subplot(numrows,1,i); hold on;
        title({['alignsyl ' num2str(alignsylnum) ',onset ' num2str(alignsylonset)], ...
            ['wind ' num2str(windsize) ', bin ' num2str(binsize)]});
        ylabel([plotstat ' (rel control)']);
        
        % ---- x values
        xtimes = AllWindowMidT(inds);
        yvals = AllPerformance(inds);
        
        % ====== DAT VS. NEG
        % - modify
        yvals_con = AllPerformanceNEG(inds);
        plotcol = 'b';
        
        % -----
        ydiff = yvals - yvals_con;
        [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
        xtimes_grp = unique(xtimes);
        
        sigvals = yCI(:,2).*yCI(:,1) >0;
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        
        % ====== DAT VS. POS
        % - modify
        yvals_con = AllPerformancePOS(inds);
        plotcol = 'r';
        
        % -----
        ydiff = yvals - yvals_con;
        [ymean, ysem, yCI] = grpstats(ydiff, xtimes, {'mean', 'sem', 'meanci'}, 'Alpha', CIalpha);
        xtimes_grp = unique(xtimes);
        
        sigvals = yCI(:,2).*yCI(:,1) >0;
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'}); % ALL
        lt_plot(xtimes_grp(sigvals), ymean(sigvals), {'Errors', ysem(sigvals), 'Color', plotcol}); % SIGNIFICANT
        
        
        
        % =========== syl contour
        % - find iteration that is closest to the end
        numiterations = length(ALLDAT.analynum(i).iterationnum);
        tmpvec = [];
        for j=1:numiterations
            tmpvec = [tmpvec ALLDAT.analynum(i).iterationnum(j).params.ClassGeneral.frtimewindow(2)];
        end
        [~, indtmp] = max(tmpvec);
        
        numbranches = length(ALLDAT.analynum(i).iterationnum(indtmp).allbranches);
        sylcontours = [];
        for j=1:numbranches
            sylcontours = [sylcontours; ...
                mean([ALLDAT.analynum(i).iterationnum(indtmp).allbranches(j).AllBranchSylContours],1)]; % for each branch, get mean, then append to overal matrix
        end
        sylcontours_mean = mean(sylcontours);
        sylcontours_std = std(sylcontours, 0, 1);
        motif_predur = ALLDAT.analynum(i).iterationnum(indtmp).params.motifpredur;
        sylcontours_x = -motif_predur + (1:length(sylcontours_mean))./1000;
        
        shadedErrorBar(sylcontours_x, -0.2+sylcontours_mean./3, sylcontours_std./3, {'Color', 'k'},1);
        
        % =========== general
        axis tight
        lt_plot_zeroline;
        
    end
    
    lt_subtitle(['strtype: ' strtype ' numctxts ' num2str(k)])
    
    
    
end
end



%% ============= PLOT 2) CLASS PERFORMANCE (NOT DIFFERNECE) - SEPARATED BY NUM CLASSES


NumCtxts = unique(AllNumCtxts);

for k=NumCtxts
    
    lt_figure; hold on;
    NumAnalys = length(unique(AllAnalynum));
    numrows = NumAnalys;
    
    for i=1:NumAnalys
        
        inds = AllAnalynum ==i & AllNumCtxts==k;
        
        alignsylnum = unique(AllAlignSylNum(inds));
        alignsylonset = unique(AllAlignSylOnset(inds));
        windsize = unique(AllFRwindowsize(inds));
        binsize = unique(AllFRbinsize(inds));
        
        assert(numel([alignsylnum alignsylonset windsize windsize]) ==4, 'should all be unique...');
        
        lt_subplot(numrows,1,i); hold on;
        title({['alignsyl ' num2str(alignsylnum) ',onset ' num2str(alignsylonset)], ...
            ['wind ' num2str(windsize) ', bin ' num2str(binsize)]});
        ylabel(plotstat);
        
        xtimes = AllWindowMidT(inds);
        xtimes_grp = unique(xtimes);

        
        % ============== ACTUAL DAT
        % -- modify
        yvals = AllPerformance(inds);
        plotcol = 'k';
        
        [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
        
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol}); % ALL
        
        
        
              % ============== NEG CONTR
        % -- modify
        yvals = AllPerformanceNEG(inds);
        plotcol = 'r';
        
        [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
        
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol}); % ALL

        
        
              % ============== POS CONTR
        % -- modify
        yvals = AllPerformancePOS(inds);
        plotcol = 'b';
        
        [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
        
        lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', plotcol}); % ALL

        
        % =========== syl contour
        % - find iteration that is closest to the end
        numiterations = length(ALLDAT.analynum(i).iterationnum);
        tmpvec = [];
        for j=1:numiterations
            tmpvec = [tmpvec ALLDAT.analynum(i).iterationnum(j).params.ClassGeneral.frtimewindow(2)];
        end
        [~, indtmp] = max(tmpvec);
        
        numbranches = length(ALLDAT.analynum(i).iterationnum(indtmp).allbranches);
        sylcontours = [];
        for j=1:numbranches
            sylcontours = [sylcontours; ...
                mean([ALLDAT.analynum(i).iterationnum(indtmp).allbranches(j).AllBranchSylContours],1)]; % for each branch, get mean, then append to overal matrix
        end
        sylcontours_mean = mean(sylcontours);
        sylcontours_std = std(sylcontours, 0, 1);
        motif_predur = ALLDAT.analynum(i).iterationnum(indtmp).params.motifpredur;
        sylcontours_x = -motif_predur + (1:length(sylcontours_mean))./1000;
        
        shadedErrorBar(sylcontours_x, 0.4+sylcontours_mean./3, sylcontours_std./3, {'Color', 'k'},1);
        
        % =========== general
        axis tight
        lt_plot_zeroline;
        
    end
    
    lt_subtitle(['strtype: ' strtype ' numctxts ' num2str(k)])
end







