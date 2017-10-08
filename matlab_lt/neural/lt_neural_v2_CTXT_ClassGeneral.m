function CLASSES = lt_neural_v2_CTXT_ClassGeneral(CLASSES, SummaryStruct, prms, CVmethod, CVkfoldnum)
%% lt 8/12/17 - takes output of lt_neural_v2_CTXT_GetBrnchDat and classified class based on neural FR vec
% CVmethod = 'LOO' or 'Kfold';
% CVkfoldnum = 10;


%%
frtimewindow = prms.ClassGeneral.frtimewindow; % on and off, relative to syl onset
frbinsize = prms.ClassGeneral.frbinsize;
Nmin = prms.ClassGeneral.Nmin;

rebalance =1;
imbalance_thr = 0.7;
beta = 0.9;

%%

numbirds = length(CLASSES.birds);

for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    birdname = CLASSES.birds(i).birdname;
    
    for ii=1:numneurons
        
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        exptname = SummaryStruct.birds(i).neurons(ii).exptID;
        
        for iii=1:numbranches
            
            
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT')
                continue
            end
            
            % ------------------ GET DATA IN FORMAT FOR CLASSIFICATION
            SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT;
            if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                clustnum = [];
            else
                clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
            end
            [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum);
            
            
            %% ========= MODEL
            ctxtlist = unique(Y);
            if length(ctxtlist)<2
                % then only one context, skip this
                continue
            end
            
            %             % --- if N too small, also skip
            %             if numel(Y) < size(Xall, 2)
            %                 continue
            %             end
            
            % ---------- PREDICTION
            %             rebalance=0;
            %             [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
            %                 = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
            %                 rebalance);
            %
            %             disp(tabulate(Y));
            %             disp('not rebalanced');
            %             disp(ConfMat);
            %             disp(['sensitiviy: ' num2str(sensitivity_mean)]);
            %             disp('--');
            %             rebalance=1;
            %             beta = 0.6;
            %             [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
            %                 = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
            %                 rebalance, 0.5, beta);
            %
            %             disp(['rebalanced, beta=' num2str(beta)]);
            %             disp(ConfMat);
            %             disp(['sensitiviy: ' num2str(sensitivity_mean)]);
            %             disp('--');
            %
            %             rebalance=1;
            %             beta = 1;
            %             [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
            %                 = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
            %                 rebalance, 0.5, beta);
            %
            %             disp(['rebalanced, beta=' num2str(beta)]);
            %             disp(ConfMat);
            %             disp(['sensitiviy: ' num2str(sensitivity_mean)]);
            %             disp('--');
            %             pause
            
            %             rebalance=1;
            if strcmp(CVmethod, 'LOO')
                [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                    = lt_neural_v2_QUICK_classifyBACKUP(Xall, Y, 'glmnet', ...
                    rebalance, imbalance_thr, beta);
            elseif strcmp(CVmethod, 'Kfold')
                [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
                    = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
                    rebalance, imbalance_thr, beta, CVkfoldnum);
            end
            %             [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
            %                 = lt_neural_v2_QUICK_classify(Xall, Y, 'matlab', ...
            %                 rebalance, imbalance_thr, beta);
            
            %             disp(ConfMat);
            %             lt_neural_ConfMatStats(ConfMat)
            %
            if isempty(Ypredicted)
                % then means failed.. (not full rank?)
                continue
            end
            
            %% === NEGATIVE CONTROL IF desired, get null distribution by shuffling link
            % between X-Y and classifying.
            if prms.ClassGeneral.GetNegControl==1
                
                AccuracyShuff = nan(prms.ClassGeneral.GetNegControl_N,1);
                SensitivityShuff = nan(prms.ClassGeneral.GetNegControl_N,1);
                PiYactualMeanShuff = nan(prms.ClassGeneral.GetNegControl_N,1);
                ConfMatShuff = nan(length(ctxtlist), length(ctxtlist), prms.ClassGeneral.GetNegControl_N);
                
                for nn=1:prms.ClassGeneral.GetNegControl_N
                    
                    tmpcount=0;
                    while tmpcount==0
                        % --- shuffle Y
                        indtmp = randperm(size(Y,1));
                        Yperm = Y(indtmp);
                        
                        % --- classify
                        if strcmp(CVmethod, 'LOO')
                            [~, ConfMattmp, accuracytmp, sensitivity_meantmp, PiYActualtmp] ...
                                = lt_neural_v2_QUICK_classifyBACKUP(Xall, Yperm, 'glmnet', ...
                                rebalance, imbalance_thr, beta);
                            
                        elseif strcmp(CVmethod, 'Kfold')
                            [~, ConfMattmp, accuracytmp, sensitivity_meantmp, PiYActualtmp] ...
                                = lt_neural_v2_QUICK_classify(Xall, Yperm, 'glmnet', ...
                                rebalance, imbalance_thr, beta, CVkfoldnum);
                        end
                        
                        if ~isempty(accuracytmp)
                            tmpcount=1;
                        end
                    end
                    AccuracyShuff(nn) = accuracytmp;
                    SensitivityShuff(nn) = sensitivity_meantmp;
                    PiYactualMeanShuff(nn) = mean(PiYActualtmp);
                    ConfMatShuff(:,:,nn) = ConfMattmp;
                end
                
                ShuffNeg_accuracy = mean(AccuracyShuff);
                ShuffNeg_accuracy_sem = lt_sem(AccuracyShuff);
                ShuffNeg_sensitivity = mean(SensitivityShuff);
                ShuffNeg_sensitivity_sem = lt_sem(SensitivityShuff);
                ShuffNeg_PiActualmean = mean(PiYactualMeanShuff);
                ShuffNeg_PiActualmean_sem = lt_sem(PiYactualMeanShuff);
            end
            
            
            %% ==== POS CONTROL - IF DESIRED, CLASSIFY POSITIVE CONTROLS
            if prms.ClassGeneral.GetPosControl ==1
                if isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT_POSCONTR')
                    % ================ 1) GET DATA IN CORRECT FORMAT
                    SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR;
                    [Xall_pos, xtimesall_pos, Y_pos, CtxtClasses_pos] = fn_extractClassDat(SEGEXTRACT, prms, clustnum);
                    
                    
                    % ================ 2) RUN CLASSIFIER
                    tmpcount=0;
                    while tmpcount==0
                        
                        if strcmp(CVmethod, 'LOO')
                            [~, ConfMat_pos, accuracy_pos, sensitivity_mean_pos, PiYActual_pos] ...
                                = lt_neural_v2_QUICK_classifyBACKUP(Xall_pos, Y_pos, 'glmnet', ...
                                rebalance, imbalance_thr, beta);
                        elseif strcmp(CVmethod, 'Kfold')
                            [~, ConfMat_pos, accuracy_pos, sensitivity_mean_pos, PiYActual_pos] ...
                                = lt_neural_v2_QUICK_classify(Xall_pos, Y_pos, 'glmnet', ...
                                rebalance, imbalance_thr, beta, CVkfoldnum);
                        end
                        
                        if ~isempty(accuracytmp)
                            tmpcount=1;
                        end
                    end
                    
                    % ================ OUTPUT
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_ConfMat = ConfMat_pos;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_accuracy = accuracy_pos;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_sensitivity_mean = sensitivity_mean_pos;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_PiYActual = PiYActual_pos;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_Xall = Xall_pos;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_Y = Y_pos;
                    
                    
                end
            end
            
            %% ================== OUTPUT
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ctxts_that_exist = unique(Y);
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ctxts_that_exist_name = CtxtClasses(unique(Y));
            
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_ConfMat = ConfMat;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_accuracy = accuracy;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_sensitivity_mean = sensitivity_mean;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_PiYActual = PiYActual;
            
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Dat_X = Xall;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Dat_Y = Y;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Dat_xbins = xtimesall(1,:);
            
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_accuracy = ShuffNeg_accuracy;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_accuracy_sem = ShuffNeg_accuracy_sem;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_sensitivity = ShuffNeg_sensitivity;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_sensitivity_sem =ShuffNeg_sensitivity_sem;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_PiActualmean=ShuffNeg_PiActualmean;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_PiActualmean_sem=ShuffNeg_PiActualmean_sem;
            CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_ConfMatAll=ConfMatShuff;
            
            disp(['done classify for ' birdname '-' exptname '-n' num2str(ii) '-branch' num2str(iii)]);
            
        end
    end
end


function [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum)

frtimewindow = prms.ClassGeneral.frtimewindow; % on and off, relative to syl onset
frbinsize = prms.ClassGeneral.frbinsize;
Nmin = prms.ClassGeneral.Nmin;

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
if version('-release')=='2013a';
    Y = nominal(Y);
else
    Y = categorical(Y);
end


%% TROUBLESHOOTING
if (0)
%% use this to try LOO vs. kfold CV, 
y = [];
for cvfold = 2:2:20;
    
[~, ConfMat_pos, accuracy_pos, sensitivity_mean_pos, PiYActual_pos] ...
                                = lt_neural_v2_QUICK_classify(Xall_pos, Y_pos, 'glmnet', ...
                                rebalance, imbalance_thr, beta, cvfold);
                            
                            y= [y accuracy_pos];
end


% LOO

                            [~, ConfMat_pos, accuracy_pos_LOO, sensitivity_mean_pos, PiYActual_pos] ...
                                = lt_neural_v2_QUICK_classifyBACKUP(Xall_pos, Y_pos, 'glmnet', ...
                                rebalance, imbalance_thr, beta);

%
figure; hold on;
plot(2:2:20, y, '-o');
line([2 50], [accuracy_pos_LOO accuracy_pos_LOO]);
ylim([0 1])

end