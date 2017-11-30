function ALLBRANCH = lt_neural_v2_CTXT_ClassSliding(CLASSES, SummaryStruct, prms, ...
    TimeWindowDur, TimeWindowSlide, FRbinsize, savenotes, CVmethod, plotstat, ...
    saveON, LinTimeWarp, regionstowarp)
%% lt 11/29/17 - added linear time warp - TO DO, can enter which regions to use
% LinTimeWarp = 1;
%             regionstowarp = [3 4];

%% lt 10/25/17 - for each branch point/neuron, performs context decoding sliding across time.
% prms.ClassGeneral.
prms.ClassSlide.Nmin = 7;
prms.ClassSlide.frbinsize = FRbinsize;

% --- Classification params
if ~exist('CVmethod', 'var')
    CVmethod = 'Kfold'; % could also be LOO
end
CVkfoldnum = min([prms.ClassSlide.Nmin, 8]); % seems comparable to LOO at 8fold.
rebalance =1;
imbalance_thr = 0.7;
beta = 0.9;

%%

savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';
tstamp = lt_get_timestamp(0);

%%

ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)

%% initiate output

ALLBRANCH = struct;
ALLBRANCH.alignpos(1).alignonset = prms.alignOnset;
ALLBRANCH.alignpos(1).alignsyl = prms.alignWhichSyl;
ALLBRANCH.SummaryStruct = SummaryStruct;
ALLBRANCH.alignpos(1).ParamsFirstIter = prms;

%% go thru all branches and perform sliding class

numbirds = length(CLASSES.birds);

for i=1:numbirds
    
    numneur = length(CLASSES.birds(i).neurons);
    birdname = CLASSES.birds(i).birdname;
    for ii=1:numneur
        
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for iii=1:numbranches
            
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT')
                continue
            end
            
            NN = [];
            for j=1:length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum)
                n = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).SegmentsExtract);
                NN = [NN n];
            end
            
            numrows = sum(NN >= prms.ClassSlide.Nmin);
            if numrows<2
                disp(['N too small!! - Skipped, ' birdname '-n' num2str(ii) '-br' num2str(iii)]);
                continue
            else
                disp(['N fine!! - continu, ' birdname '-n' num2str(ii) '-br' num2str(iii)]);
            end
            
            % *********************** confirm that pos control has same sampel size
            NN2 = [];
            for j=1:length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum)
                n = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(j).SegmentsExtract);
                NN2 = [NN2 n];
            end
            
            NN = NN(NN>=prms.ClassSlide.Nmin);
            NN2 = NN2(NN2>=prms.ClassSlide.Nmin);
            
            tmp = tabulate(NN)==tabulate(NN2);
            assert(all(tmp(:)==1), 'pos control not same sampoel size as dat');
            % ***********************
            
            ConfMatAll = nan(numrows, numrows, size(ListOfTimeWindows,1));
            ConfMatAll_NEG = nan(numrows, numrows, size(ListOfTimeWindows,1));
            ConfMatAll_POS = nan(numrows, numrows, size(ListOfTimeWindows,1));
            
            SEGEXTRACT_DAT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT;
            
            if LinTimeWarp ==1
                % ==== 1) combine all classes into one segextract
                segextractAll = SEGEXTRACT_DAT.classnum(1).SegmentsExtract;
                originalinds = [length(SEGEXTRACT_DAT.classnum(1).SegmentsExtract)];
                for j=2:length(SEGEXTRACT_DAT.classnum)
                    segextractAll = [segextractAll ...
                        SEGEXTRACT_DAT.classnum(j).SegmentsExtract];
                    
                    originalinds = [originalinds ...
                        length(SEGEXTRACT_DAT.classnum(j).SegmentsExtract)];
                end
                
                % ==== 2) linear time warp
                segextractAll = lt_neural_LinTimeWarpSegmented(segextractAll, regionstowarp);
                
                % ===== put back into SEGEXTRACT
                count=1;
                for j=1:length(SEGEXTRACT_DAT.classnum)
%                     disp(count:(count+originalinds(j)-1));
%                     disp('--');
                    SEGEXTRACT_DAT.classnum(j).SegmentsExtract = segextractAll(count:(count+originalinds(j)-1));
                    count = (count+originalinds(j));
                end
            end
            
            % ===================== go thru all time windows and run
            % classifier
            for j=1:size(ListOfTimeWindows,1)
                
                prms.classtmp.frtimewindow = ListOfTimeWindows(j,:);
                
                % ############################# GET DATA IN FORMAT FOR CLASSIFICATION
                SEGEXTRACT = SEGEXTRACT_DAT;
                if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    clustnum = [];
                else
                    clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                end
                
                
                [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum);
                
                
                
                %% DATA
                % ############################## CLASSIFY
                if strcmp(CVmethod, 'LOO')
                    [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                        = lt_neural_v2_QUICK_classifyBACKUP(Xall, Y, 'glmnet', ...
                        rebalance, imbalance_thr, beta);
                elseif strcmp(CVmethod, 'Kfold')
                    [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                        = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
                        rebalance, imbalance_thr, beta, CVkfoldnum);
                end
                
                if isempty(ConfMat)
                    % then means failed.. (not full rank?)
                    disp('FAILED - confmat is empty!! ');
                    keyboard
                end
                
                assert(all(size(ConfMat)==[numrows numrows]), 'asdfa');
                ConfMatAll(:,:, j) = ConfMat;
                
                
                
                %% NEG CONTROL [just one iteration]
                if prms.ClassSlide.GetNegControl==1
                    
                    tmpcount=0;
                    while tmpcount==0
                        % --- shuffle Y
                        indtmp = randperm(size(Y,1));
                        Yperm = Y(indtmp);
                        
                        % --- classify
                        if strcmp(CVmethod, 'LOO')
                            [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                                = lt_neural_v2_QUICK_classifyBACKUP(Xall, Yperm, 'glmnet', ...
                                rebalance, imbalance_thr, beta);
                        elseif strcmp(CVmethod, 'Kfold')
                            [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                                = lt_neural_v2_QUICK_classify(Xall, Yperm, 'glmnet', ...
                                rebalance, imbalance_thr, beta, CVkfoldnum);
                        end
                        
                        if ~isempty(ConfMat)
                            tmpcount=1;
                        end
                    end
                    
                    ConfMatAll_NEG(:,:, j) = ConfMat;
                end
                
                
            end
            
            % ============== MAKE SURE NO NAN
            if(any(isnan(ConfMatAll(:))))
                disp('FAILED - some isnan in confmatall');
                keyboard
            end
            
            if(any(isnan(ConfMatAll_NEG(:))))
                disp('FAILED - some isnan in confmatall');
                keyboard
            end
            
            % ======================================= SAVE OUTPUT
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll = ConfMatAll;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll_NEG = ConfMatAll_NEG;
            
            
            
            %% POS CONTROL
            for j=1:size(ListOfTimeWindows,1)
                
                prms.classtmp.frtimewindow = ListOfTimeWindows(j,:);
                
                % ############################# GET DATA IN FORMAT FOR CLASSIFICATION
                SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR;
                if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    clustnum = [];
                else
                    clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                end
                [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum);
                
                
                
                % ############################## CLASSIFY
                if strcmp(CVmethod, 'LOO')
                    [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                        = lt_neural_v2_QUICK_classifyBACKUP(Xall, Y, 'glmnet', ...
                        rebalance, imbalance_thr, beta);
                elseif strcmp(CVmethod, 'Kfold')
                    [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                        = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
                        rebalance, imbalance_thr, beta, CVkfoldnum);
                end
                
                if isempty(ConfMat)
                    % then means failed.. (not full rank?)
                    disp('FAILED - confmat is empty!! ');
                    keyboard
                end
                
                assert(all(size(ConfMat)==[numrows numrows]), 'asdfa');
                ConfMatAll_POS(:,:, j) = ConfMat;
                
            end
            
            % ============== MAKE SURE NO NAN
            if(any(isnan(ConfMatAll_POS(:))))
                disp('FAILED - some isnan in confmatall');
                keyboard
            end
            
            % ======================================= SAVE OUTPUT
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll_POS = ConfMatAll_POS;
            
            
            
            
            
            %% ############ save stats about this branch/neuron in old format
            
            % === xvals (ALL)
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).xtimes = mean(ListOfTimeWindows,2)';
            
            % ===== yvals (DAT)
            confmat = ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll;
            yvals = nan(1,size(confmat,3));
            for j=1:size(confmat,3)
                sts = lt_neural_ConfMatStats(confmat(:,:,j));
                yvals(j) = sts.(plotstat);
            end
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).yvals = yvals;
            
            
            % ===== yvals (NEG)
            confmat = ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll_NEG;
            yvals = nan(1,size(confmat,3));
            for j=1:size(confmat,3)
                sts = lt_neural_ConfMatStats(confmat(:,:,j));
                yvals(j) = sts.(plotstat);
            end
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).yvals_neg = yvals;
            
            % ===== yvals (POS)
            confmat = ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).ConfMatAll_POS;
            yvals = nan(1,size(confmat,3));
            for j=1:size(confmat,3)
                sts = lt_neural_ConfMatStats(confmat(:,:,j));
                yvals(j) = sts.(plotstat);
            end
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).yvals_pos = yvals;
            
            
            
            % ##################################### OTHER INFORMATION [DAT]
            SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT;
            numclasses = length(SEGEXTRACT.classnum);
            
            maxdur = prms.motifpredur+prms.motifpostdur-0.005;
            numtrials = 10;
            SylContours = [];
            %             Ctxtnum = [];
            %             Xtimes = [];
            SylContoursByClass_means = [];
            SylContoursByClass_std = [];
            SylContoursByClass_N = [];
            SylContoursByClass_classnums = [];
            
            for j=1:numclasses
                SegmentsExtract = SEGEXTRACT.classnum(j).SegmentsExtract;
                
                if length(SegmentsExtract)<prms.ClassSlide.Nmin
                    continue
                end
                
                % ================= FR
                if isfield(ALLBRANCH.SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    clustnum = '';
                else
                    clustnum = ALLBRANCH.SummaryStruct.birds(i).neurons(ii).clustnum;
                end
                
                SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR.classnum(j).FRsmooth_rate_CommonTrialDur ...
                    = single([SegmentsExtract.FRsmooth_rate_CommonTrialDur]);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR.classnum(j).FRsmooth_xbin_CommonTrialDur ...
                    = single(SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR.classnum(j).regexpstr = ...
                    SEGEXTRACT.classnum(j).regexpstr;
                
                
                % ============== syldur
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs.classnum(j).Dur_syl = ...
                    [SegmentsExtract.Dur_syl];
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs.classnum(j).Dur_gappre = ...
                    [SegmentsExtract.Dur_gappre];
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs.classnum(j).Dur_gappost = ...
                    [SegmentsExtract.Dur_gappost];
                
                
                % =============== syl contour
                [sylcon, xtimes] = lt_neural_v2_ANALY_GetSylContours(SegmentsExtract, ...
                    numtrials, maxdur);
                
                % --- 1) ALL COMBINED
                SylContours = [SylContours; sylcon];
                %                 Ctxtnum = [Ctxtnum; j*ones(size(sylcon,1),1)];
                %                 Xtimes = [Xtimes; xtimes];
                
                % ---- 2) CLASS BY CLASS
                SylContoursByClass_means = [SylContoursByClass_means; mean(sylcon,1)];
                SylContoursByClass_std = [SylContoursByClass_std; std(sylcon, 0, 1)];
                SylContoursByClass_N = [SylContoursByClass_N; numtrials];
                SylContoursByClass_classnums = [SylContoursByClass_classnums; j];
                
            end
            
            assert(sum(isnan(SylContours(:)))./numel(SylContours) < 0.02, 'problem - >2% are nan');
            SylContours = int8(SylContours);
            sylcontours_mean = nanmean(SylContours,1);
            sylcontours_x = xtimes - prms.motifpredur;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).sylcontours_mean = sylcontours_mean;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).sylcontours_x = sylcontours_x;
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylContoursByClass_means = SylContoursByClass_means;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylContoursByClass_std = SylContoursByClass_std;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylContoursByClass_N = SylContoursByClass_N;
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylContoursByClass_classnums = SylContoursByClass_classnums;
            
            
            
            % ##################################### OTHER INFORMATION [POS CONTROL]
            SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR;
            numclasses = length(SEGEXTRACT.classnum);
            for j=1:numclasses
                SegmentsExtract = SEGEXTRACT.classnum(j).SegmentsExtract;
                
                if length(SegmentsExtract)<prms.ClassSlide.Nmin
                    continue
                end
                
                % ================= FR
                if isfield(ALLBRANCH.SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    clustnum = '';
                else
                    clustnum = ALLBRANCH.SummaryStruct.birds(i).neurons(ii).clustnum;
                end
                
                SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR_POSCONTR.classnum(j).FRsmooth_rate_CommonTrialDur ...
                    = single([SegmentsExtract.FRsmooth_rate_CommonTrialDur]);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR_POSCONTR.classnum(j).FRsmooth_xbin_CommonTrialDur ...
                    = single(SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur);
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).FR_POSCONTR.classnum(j).regexpstr = ...
                    SEGEXTRACT.classnum(j).regexpstr;
                
                
                % ============== syldur
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs_POS.classnum(j).Dur_syl = ...
                    [SegmentsExtract.Dur_syl];
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs_POS.classnum(j).Dur_gappre = ...
                    [SegmentsExtract.Dur_gappre];
                
                ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).SylGapDurs_POS.classnum(j).Dur_gappost = ...
                    [SegmentsExtract.Dur_gappost];
            end
            
            
            
            
            % ################################################### OTHER
            % INFO
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).prms_regexpstr ...
                = {CLASSES.birds(i).neurons(ii).branchnum(iii).regexprstr};
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).prms_regexpstrlist ...
                = {CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum.regexpstr};
            
            ALLBRANCH.alignpos(1).bird(i).branch(iii).neuron(ii).prms_regexpstrlist_POSCONTR ...
                = {CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum.regexpstr};
            
        end
    end
    
    %% save after each bird (overwrite)
    %% ========== save
    if saveON==1
        
        
        % --- allbranch
        fname = [savedir '/ALLBRANCHv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
            'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
        
        save(fname, 'ALLBRANCH', '-v7.3');
        
        % --- classes
        fname = [savedir '/CLASSESv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
            'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
        save(fname, 'CLASSES', '-v7.3')
        
        % --- summary struct
        fname = [savedir '/SUMMARYv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
            'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
        save(fname, 'SummaryStruct')
        
        % -- params
        fname = [savedir '/PARAMSv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
            'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
        save(fname, 'prms')
        
        
    end
    
    
end


%% ============= look at variability of gap durations

ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(ALLBRANCH);


%% ========== save
if saveON==1
    
    
    % --- allbranch
    fname = [savedir '/ALLBRANCHv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
        'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
    
    save(fname, 'ALLBRANCH', '-v7.3');
    
    % --- classes
    fname = [savedir '/CLASSESv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
        'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
    save(fname, 'CLASSES', '-v7.3')
    
    % --- summary struct
    fname = [savedir '/SUMMARYv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
        'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
    save(fname, 'SummaryStruct')
    
    % -- params
    fname = [savedir '/PARAMSv2_' prms.Extract.strtype '_Algn' num2str(prms.alignWhichSyl) ...
        'Ons' num2str(prms.alignOnset) '_' tstamp '_' savenotes '.mat'];
    save(fname, 'prms')
    
end

end

function [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum)

frtimewindow = prms.classtmp.frtimewindow; % on and off, relative to syl onset
frbinsize = prms.ClassSlide.frbinsize;
Nmin = prms.ClassSlide.Nmin;

numclasses = length(SEGEXTRACT.classnum);

% =================== COLLECT DATA TO CLASSIFY
Xall = []; % trials x bins (FR vectors)
xtimesall = []; % 1 x bins (time stamps)
Y = []; % context indicator
CtxtClasses = {};
contextcounter = 0;


% ------ methiod2
for j=1:numclasses
    
    sylname = SEGEXTRACT.classnum(j).regexpstr;
    
    % --- EXTRACT DATA
    segextract = SEGEXTRACT.classnum(j).SegmentsExtract;
    
    % -- extract FR
    segextract = lt_neural_SmoothFR(segextract, clustnum);
    
    
    if ~isfield(segextract, 'spk_Times')
        % then no data
        continue
    end
    
    if length(segextract) < Nmin
        % not neough data
        continue
    end
    
    
    % ---- EXTRAC FR VECTOR (within desired time window)
    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
    indsFR = xbin>(prms.motifpredur+frtimewindow(1)+0.0001) ...
        & xbin<=(prms.motifpredur+frtimewindow(2)+0.0001); % add/minus at ends because somtimes
    
    X = [segextract.FRsmooth_rate_CommonTrialDur];
    X = X(indsFR, :);
    X = X';
    
    xtimes = xbin(indsFR);
    xtimes = xtimes';
    
    
    % ------------ reduce Dim of FR vectors (by binning in
    % time)
    TrimDown = 1;
    [X, xtimes] = lt_neural_v2_QUICK_binFR(X, xtimes, frbinsize, TrimDown);
    
    
    % ======================== COLLECT ACROSS ALL CLASSES
    contextcounter = contextcounter+1;
    
    CtxtClasses = [CtxtClasses sylname];
    
    Xall = [Xall; X];
    xtimesall = [xtimesall; xtimes];
    Y = [Y; contextcounter*ones(size(X,1),1)];
    
end

% --- convert Y to categorical array
if version('-release')=='2013a'
    Y = nominal(Y);
else
    Y = categorical(Y);
end
end