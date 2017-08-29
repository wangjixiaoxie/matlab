function [ALLBRANCH] = lt_neural_v2_CTXT_PlotAll(strtype, plotstat)
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

SummaryStructAll = [];
tmpbytes = [];
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
        
        sstruct = tmp.CLASSEScompiled.analynum(ii).SummaryStruct;
        ALLDAT.analynum(count).SummaryStruct = sstruct;
        
        dummy = whos('sstruct');
        tmpbytes = dummy.bytes;
        
        % --- original params
        ddd = [tmp.CLASSEScompiled.analynum(ii).iterationnum.params];
        assert(length(unique([ddd.motifpredur])) ==1, 'asdfas');
        assert(length(unique([ddd.motifpostdur])) ==1, 'asdfasd');
        ALLDAT.analynum(count).ParamsFirstIter = tmp.CLASSEScompiled.analynum(ii).iterationnum(1).params;
        
    end
end

assert(length(unique(tmpbytes))==1, 'rpbolem, diff symmary strucuts...');

% SummaryStructAll(1).birds(1).neurons(1).

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

AllBirdnum = [];
AllNeurNum = [];
% AllBranchNum = []; % meaningless, based on incrememtal counter
AllRegExp = {}; % actual
AllSampleSize = [];

AllIdx_analy_iter_branch = []; % to index back into ALLDAT

AllRegExpActual = {};
AllRegExpActual_PosContr = {};
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
        
        
        %% ================ branch information
        AllBirdnum = [AllBirdnum ALLDAT.analynum(i).iterationnum(ii).allbranches.AllBirdnum];
        AllNeurNum = [AllNeurNum ALLDAT.analynum(i).iterationnum(ii).allbranches.AllNeurnum];
        if (0)
            % this is arbitrary branch number, used to index into ALLDAT.
            % But is meaningless wrt to the actual syllables
            % INSTEAD, will decide branch class based on regular
            % expression string, below
            AllBranchNum = [AllBranchNum ALLDAT.analynum(i).iterationnum(ii).allbranches.AllBranchnum];
        end
        
        functmp = @(X) sum(X(:));
        AllSampleSize = [AllSampleSize ...
            cellfun(functmp, {ALLDAT.analynum(i).iterationnum(ii).allbranches.AllConfMat})];
        AllRegExp = [AllRegExp {ALLDAT.analynum(i).iterationnum(ii).allbranches.AllBranchRegexp}];
        
        AllRegExpActual = [AllRegExpActual {ALLDAT.analynum(i).iterationnum(ii).allbranches.AllBranchCtxts}];
        
        AllRegExpActual_PosContr = [AllRegExpActual_PosContr ...
            {ALLDAT.analynum(i).iterationnum(ii).allbranches.POSCONTR_regexpstrlist}];
        
        
        %% === debug
        if (0)
            FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
            FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
            % +1 is 1 after token
            FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
            collectWNhit = 0;
            LearnKeepOnlyBase = 1;
            
            branchnum = 1;
            birdnum = ALLDAT.analynum(i).iterationnum(ii).allbranches(branchnum).AllBirdnum;
            neurnum = ALLDAT.analynum(i).iterationnum(ii).allbranches(branchnum).AllNeurnum;
            
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(ALLDAT.analynum(i).SummaryStruct, birdnum, neurnum);
            
            reglistdat = ALLDAT.analynum(i).iterationnum(ii).allbranches(branchnum).AllBranchCtxts;
            for lll=1:length(reglistdat);
                regstr = reglistdat{lll};
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                    regstr, 0.1, 0.1, alignonset, '', FFparams, ...
                    0, 1, collectWNhit, 0, LearnKeepOnlyBase, 1);
                disp(length(SegmentsExtract));
            end
            
            reglistdat = ALLDAT.analynum(i).iterationnum(ii).allbranches(branchnum).POSCONTR_regexpstrlist;
            for lll=1:length(reglistdat);
                regstr = reglistdat{lll};
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                    regstr, 0.1, 0.1, alignonset, '', FFparams, ...
                    0, 1, collectWNhit, 0, LearnKeepOnlyBase, 1);
                disp(length(SegmentsExtract));
            end
        end
        
        %% ===============================
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
        
        % ----- to index back into ALLDAT
        branchnums = 1:length(ALLDAT.analynum(i).iterationnum(ii).allbranches);
        branchnums = branchnums';
        
        AllIdx_analy_iter_branch = [AllIdx_analy_iter_branch; [i*ones(size(branchnums)) ...
            ii*ones(size(branchnums)) branchnums]];
      
        
    end
end

% ================================== GIVE A UNIQUE BRANCH ID TO EACH BRANCH
% BASED ON REG EXP STRING (ACROSS ALL BIRDS, SO NEED TO DO INTERSECT WITH
% BIRD TO GET ACTUAL REAL BRANCH)
[AllBranchID] = grp2idx(AllRegExp);
AllBranchID = AllBranchID';

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



%% ============================= EACH BRANCH/NEURON AS ONE DATAPOINT
% % ======== 1) EXTRACT EACH BRANCH/NEURON
% % iterate thru each unique combination of analysis/bird/neuron/branch
% plotraw = 1;
%
% numanalys = max(AllAnalynum);
% numbirds = max(AllBirdnum);
% numneurons = max(AllNeurNum);
% numbranches = max(AllBranchID);
%
% figcount=1;
% subplotrows=6;
% subplotcols=2;
% fignums_alreadyused=[];
% hfigs=[];
%
%
% disp('alignsyl -- onset? -- binsize -- windsize -- numctxts -- rawsampsize');
%
% ALLBRANCH = struct;
%
% for i=1:numanalys
%
%
%     % - find iteration that is closest to the end [will be used to get
%     % syl contours]
%     numiterations = length(ALLDAT.analynum(i).iterationnum);
%     tmpvec = [];
%     for j=1:numiterations
%         tmpvec = [tmpvec ALLDAT.analynum(i).iterationnum(j).params.ClassGeneral.frtimewindow(2)];
%     end
%     [~, itermax] = max(tmpvec);
%
%
%     for j=1:numbirds
%         for l = 1:numbranches
%             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%
%             title(['bird' num2str(j) ', branch' num2str(l)]);
%
%             YvalsAll = [];
%             XtimesAll = [];
%             for k=1:numneurons
%
%                 inds = AllAnalynum==i & AllBirdnum==j & AllNeurNum==k ...
%                     & AllBranchID==l;
%
%                 if ~any(inds)
%                     continue
%                 end
%
%
%                 % ==== things about this neuron
%
%
%                 % ==== sanity checks, make sure classifier params are the
%                 % same for all in inds
%                 alignsyl = unique(AllAlignSylNum(inds));
%                 alignonset = unique(AllAlignSylOnset(inds));
%                 binsize = unique(AllFRbinsize(inds));
%                 windowsize = unique(AllFRwindowsize(inds));
%                 numctxts = unique(AllNumCtxts(inds));
%                 rawsampsize = unique(AllSampleSize(inds));
%
%                 assert(numel([alignsyl; alignonset; binsize; windowsize; numctxts; rawsampsize]) ==6, 'problem, should be unique///');
%                 disp(num2str([alignsyl; alignonset; binsize; windowsize; numctxts; rawsampsize]'));
%
%                 regexpstr = unique(AllRegExp(inds));
%                 assert(length(regexpstr)==1, 'should be uniuqe');
%
%                 % === for this branch and neuron, get timecourse of
%                 % classifier performance
%                 xtimes = AllWindowMidT(inds);
%                 yvals = AllPerformance(inds);
%                 yvals_neg = AllPerformanceNEG(inds);
%                 yvals_pos = AllPerformancePOS(inds);
%
%                 % -- sort
%                 [~, indstmp] = sort(xtimes);
%                 xtimes = xtimes(indstmp);
%                 yvals = yvals(indstmp);
%                 yvals_neg = yvals_neg(indstmp);
%                 yvals_pos = yvals_pos(indstmp);
%
%                 % === for this branch and neuron, get syl contours
%                 sylcontours = ALLDAT.analynum(i).iterationnum(itermax).allbranches(l).AllBranchSylContours;
%                 motif_predur = ALLDAT.analynum(i).iterationnum(itermax).params.motifpredur;
%                 sylcontours_mean = mean(sylcontours,1); % -- take mean syl contour
%                 sylcontours_x = -motif_predur + (1:length(sylcontours_mean))./1000;
%                 if (plotraw)
%                     lt_figure; hold on;
%                     subplot(311); hold on;
%                     spy(sylcontours(:,1:2:end));
%                     subplot(312); hold on;
%                     plot(sylcontours_x(1:2:end), sylcontours(:, 1:2:end), '-b');
%                     ylim([-0.1 1.1]);
%                     subplot(313); hold on;
%                     plot(sylcontours_x(1:2:end), sylcontours_mean(1:2:end), '-k');
%                 end
%
%
%
%                 % -- plot random subset
%                 % one fig for each branch, combining neurons
%                 if (plotraw)
%                     if rand<0.3;
%                         cla
%                         plot(xtimes, yvals, 'k-x');
%                         plot(sylcontours_x, sylcontours_mean, 'r');
%                         ylim([-0.1 1.1]);
%                         pause;
%                     end
%                 end
%
%                 % --- collect to plot mean
%                 try
%                     YvalsAll = [YvalsAll; yvals];
%                     XtimesAll = [XtimesAll; xtimes];
%                 catch err
%                 end
%
%                 % ========================================= OUTPUT
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).xtimes = xtimes;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).yvals = yvals;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).yvals_neg = yvals_neg;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).yvals_pos = yvals_pos;
%
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).sylcontours_mean = sylcontours_mean;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).sylcontours_x = sylcontours_x;
%
%                 if isfield(ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k), 'prm_alignsyl')
%                     assert(isempty(ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prm_alignsyl), ...
%                         'problem - multiple analyses overwriting - lilely used diff windowsize, binsize, or sliding size');
%                 end
%
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_alignsyl = alignsyl;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_alignonset = alignonset;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_binsize = binsize;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_windowsize = windowsize;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_numctxts = numctxts;
%                 ALLBRANCH.bird(j).branch(l).alignsyl(alignsyl).alignonset(alignonset).neuron(k).prms_regexpstr = regexpstr;
%
%             end
%
%             % --- plot mean across neurons
%             if size(YvalsAll,1)>1
%                 cla
%                 X = mean(XtimesAll, 1);
%                 Y = mean(YvalsAll, 1);
%                 Ysem = lt_sem(YvalsAll);
%                 shadedErrorBar(X, Y, Ysem, {'Color', 'k'},1);
%                 plot(sylcontours_x, 0.4+sylcontours_mean./3, 'r');
%                 lt_plot_annotation(1, {['alignsyl' num2str(alignsyl) ', onset' num2str(alignonset)], ...
%                     ['numctxt' num2str(numctxts) ', ' regexpstr{1}]}, 'r');
%             end
%
%
%         end
%     end
% end

%% ============================= EACH BRANCH/NEURON AS ONE DATAPOINT
% ======== 1) EXTRACT EACH BRANCH/NEURON
% iterate thru each unique combination of analysis/bird/neuron/branch
plotrawcontours = 0; % plots each neuron/branch contours
plotrawcontoursdat = 0; % plots class performance and contours
plotmeancontoursdat = 0;

numanalys = max(AllAnalynum);
numbirds = max(AllBirdnum);
numneurons = max(AllNeurNum);
numbranches = max(AllBranchID);

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


disp('alignsyl -- onset? -- binsize -- windsize -- numctxts -- rawsampsize');

ALLBRANCH = struct;

for i=1:numanalys
    
    
    % -- what unique syl and onset alignment is this?
    alignsyl = unique(AllAlignSylNum(AllAnalynum==i));
    alignonset = unique(AllAlignSylOnset(AllAnalynum==i));
    
    assert(length([alignsyl alignonset]) ==2, 'why not uniqu?');
    
    ALLBRANCH.alignpos(i).ParamsFirstIter = ALLDAT.analynum(i).ParamsFirstIter;
    
    for j=1:numbirds
        for l = 1:numbranches
            
            if plotrawcontoursdat==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['bird' num2str(j) ', branch' num2str(l)]);
            end
            YvalsAll = [];
            XtimesAll = [];
            for k=1:numneurons
                
                inds = AllAnalynum==i & AllBirdnum==j & AllNeurNum==k ...
                    & AllBranchID==l;
                
                if ~any(inds)
                    continue
                end
                
                
                % ==== things about this neuron
                
                
                % ==== sanity checks, make sure classifier params are the
                % same for all in inds
                %                 alignsyl = unique(AllAlignSylNum(inds));
                %                 alignonset = unique(AllAlignSylOnset(inds));
                binsize = unique(AllFRbinsize(inds));
                windowsize = unique(AllFRwindowsize(inds));
                numctxts = unique(AllNumCtxts(inds));
                rawsampsize = unique(AllSampleSize(inds));
                
                assert(numel([binsize; windowsize; numctxts; rawsampsize]) ==4, 'problem, should be unique///');
                disp(num2str([binsize; windowsize; numctxts; rawsampsize]'));
                
                regexpstr = unique(AllRegExp(inds));
                assert(length(regexpstr)==1, 'should be uniuqe');
                
                % ==== regexpactual
                assert(length(unique(cellfun('length', AllRegExpActual(inds))))==1, 'asdfas');
                regexpstr_actual = AllRegExpActual(inds);
                regexpstr_actual = regexpstr_actual{1};
                
                
                % ================================= stuff about positive
                % control
                % === pos control regexp
                assert(length(unique(cellfun('length', AllRegExpActual_PosContr(inds))))==1, 'adsfds');
                regexpstr_actual_pos = AllRegExpActual_PosContr(inds);
                regexpstr_actual_pos = regexpstr_actual_pos{1};
                
                indstmp = find(inds);
                kk=indstmp(1); % use first, since all should be same (i.e. each time bin. ..)
                    dattmp = ALLDAT.analynum(AllIdx_analy_iter_branch(kk, 1)).iterationnum(AllIdx_analy_iter_branch(kk,2)).allbranches(AllIdx_analy_iter_branch(kk,3));
                    POSCONTR_regexpstrlist_N = dattmp.POSCONTR_regexpstrlist_N;
                
                % === for this branch and neuron, get timecourse of
                % classifier performance
                xtimes = AllWindowMidT(inds);
                yvals = AllPerformance(inds);
                yvals_neg = AllPerformanceNEG(inds);
                yvals_pos = AllPerformancePOS(inds);
                
                % -- sort
                [~, indstmp] = sort(xtimes);
                xtimes = xtimes(indstmp);
                yvals = yvals(indstmp);
                yvals_neg = yvals_neg(indstmp);
                yvals_pos = yvals_pos(indstmp);
                
                %                 tmptmp = tabulate(round(1000*diff(xtimes)));
                %                 assert(max(tmptmp(:,3))>80, 'problem, there are numerous xslide windows ...');
                
                
                % === for this branch and neuron, get syl contours (will
                % average over all analyses (i.e. xtime bin position)
                indstmp = find(inds);
                sylcontour = [];
                for kk=indstmp
                    dattmp = ALLDAT.analynum(AllIdx_analy_iter_branch(kk, 1)).iterationnum(AllIdx_analy_iter_branch(kk,2)).allbranches(AllIdx_analy_iter_branch(kk,3));
                    assert(dattmp.AllBirdnum==j, 'asfasdf');
                    assert(dattmp.AllNeurnum==k', 'asdfasd');
                    assert(strcmp(dattmp.AllBranchRegexp, regexpstr)==1, 'asdfd');
                    
                    sylcontour = [sylcontour; mean(dattmp.AllBranchSylContours,1)];
                    motif_predur = ALLDAT.analynum(AllIdx_analy_iter_branch(kk, 1)).iterationnum(AllIdx_analy_iter_branch(kk,2)).params.motifpredur;
                    
                    if (plotrawcontours) & kk == min(indstmp)
                        lt_figure; hold on;
                        subplot(311); hold on;
                        spy(dattmp.AllBranchSylContours(:,1:2:end));
                        subplot(312); hold on;
                        plot(dattmp.AllBranchSylContours(:, 1:2:end)', '-b');
                        ylim([-0.1 1.1]);
                        subplot(313); hold on;
                        plot(mean(dattmp.AllBranchSylContours(:,1:2:end)), '-k');
                    end
                end
                
                
                % 
                sylcontours_mean = mean(sylcontour,1);
                clear sylcontour
                sylcontours_x = -motif_predur + (1:length(sylcontours_mean))./1000;
                
                
                % -- plot random subset
                % one fig for each branch, combining neurons
                if (plotrawcontoursdat)
                    if rand<0.5;
                        cla
                        plot(xtimes, yvals, 'k-x');
                        plot(sylcontours_x, sylcontours_mean, 'r');
                        ylim([-0.1 1.1]);
                        pause;
                    end
                end
                
                % --- collect to plot mean
                try
                    YvalsAll = [YvalsAll; yvals];
                    XtimesAll = [XtimesAll; xtimes];
                catch err
                end
                
                
                % ========================================= OUTPUT
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).xtimes = single(xtimes);
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).yvals = single(yvals);
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).yvals_neg = yvals_neg;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).yvals_pos = yvals_pos;
                
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).sylcontours_mean = single(sylcontours_mean);
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).sylcontours_x = single(sylcontours_x);
                
                if isfield(ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k), 'prm_alignsyl')
                    assert(isempty(ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prm_alignsyl), ...
                        'problem - multiple analyses overwriting - lilely used diff windowsize, binsize, or sliding size');
                end
                
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_alignsyl = alignsyl;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_alignonset = alignonset;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_binsize = binsize;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_windowsize = windowsize;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_numctxts = numctxts;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_regexpstr = regexpstr;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_regexpstrlist = regexpstr_actual;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_regexpstrlist_POSCONTR = regexpstr_actual_pos;
                ALLBRANCH.alignpos(i).bird(j).branch(l).neuron(k).prms_POSCONTR_regexpstrlist_N = POSCONTR_regexpstrlist_N;
                
                
            end
            
            ALLBRANCH.alignpos(i).alignsyl = alignsyl;
            ALLBRANCH.alignpos(i).alignonset = alignonset;
            
            
            
            % --- plot mean across neurons
            if size(YvalsAll,1)>1
                if plotmeancontoursdat==1
                    if  plotrawcontoursdat==0
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        
                        title(['bird' num2str(j) ', branch' num2str(l)]);
                    end
                    cla
                    X = mean(XtimesAll, 1);
                    Y = mean(YvalsAll, 1);
                    Ysem = lt_sem(YvalsAll);
                    shadedErrorBar(X, Y, Ysem, {'Color', 'k'},1);
                    plot(sylcontours_x, 0.3+sylcontours_mean./5, 'r');
                    lt_plot_annotation(1, {['alignsyl' num2str(alignsyl) ', onset' num2str(alignonset)], ...
                        ['numctxt' num2str(numctxts) ', ' regexpstr{1}]}, 'r');
                    ylim([0.2 0.7])
                end
            end
            
            
        end
    end
end

ALLBRANCH.SummaryStruct = ALLDAT.analynum(1).SummaryStruct; % already asserted that all sstructs are the same.




