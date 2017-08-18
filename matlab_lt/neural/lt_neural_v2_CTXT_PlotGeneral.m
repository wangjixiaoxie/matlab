function lt_neural_v2_CTXT_PlotGeneral(CLASSES, SummaryStruct, prms)
%% lt 8/12/17 - plot results from classifier lt_neural_v2_CTXT_ClassGeneral

%% ==============

numbirds = length(CLASSES.birds);


%%
AllNumCtxts = [];
AllAccuracy = [];
AllSensitivity = [];
AllPihat = [];

AllF1 = [];
AllF1_neg = [];
AllF1_pos = [];

for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    birdname = CLASSES.birds(i).birdname;
    for ii=1:numneurons
        
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for iii=1:numbranches
            
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'CLASSIFIER')
                continue
            end
            
            if isempty(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER)
                continue
            end
            
            numctxts = numel(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ctxts_that_exist);
            accuracy = CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_accuracy;
            sensitivity = CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_sensitivity_mean;
            pihatactual = mean(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_PiYActual);
            
            if length(numctxts) ~= length(accuracy)
                disp('skipped - no dat since not full rank?')
                continue
            end
            
            % -------------- dat
            confmat = CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.Rslts_ConfMat;
            stats = lt_neural_ConfMatStats(confmat);
            AllF1 = [AllF1 stats.F1];
            
            % ----- pos control
            if isfield(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER, 'ContrPos_ConfMat');
            confmat = CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ContrPos_ConfMat;
            stats = lt_neural_ConfMatStats(confmat);
            AllF1_pos = [AllF1_pos stats.F1];
            end
            
             % ----- neg control
            if isfield(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER, 'ShuffNeg_ConfMatAll');
                f1tmp = [];
                for j=1:size(CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_ConfMatAll,3)
                   confmat =  CLASSES.birds(i).neurons(ii).branchnum(iii).CLASSIFIER.ShuffNeg_ConfMatAll(:,:,j);
                   stats = lt_neural_ConfMatStats(confmat);
                   f1tmp = [f1tmp stats.F1];
                   if isnan(stats.F1)
                       keyboard
                   end
                end
                AllF1_neg = [AllF1_neg mean(f1tmp)];
            end
           
            % ========== COLLECT
            AllNumCtxts = [AllNumCtxts; numctxts];
            AllAccuracy = [AllAccuracy accuracy];
            AllSensitivity = [AllSensitivity sensitivity];
            AllPihat = [AllPihat    pihatactual];
            
        end
    end
end

%% ============= PLOT ALL COMBINED
lt_figure; hold on;

lt_subplot(2,2,1); hold on;
xlabel('numctxts');
ylabel('accuracy');

plot(AllNumCtxts, AllAccuracy, 'ok');

[ymean, yN] = grpstats(AllAccuracy, AllNumCtxts, {'mean', 'numel'});

ctxtslist = unique(AllNumCtxts);

lt_plot_bar(ctxtslist, ymean, {'Color', 'b'});
lt_plot_bar(ctxtslist, 1./ctxtslist, {'Color', 'k'}); % null

xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);


% ----------------
lt_subplot(2,2,2); hold on;
xlabel('numctxts');
ylabel('sensitivity (pr(pred A if A))');

plot(AllNumCtxts, AllSensitivity, 'ok');

[ymean, yN] = grpstats(AllSensitivity, AllNumCtxts, {'mean', 'numel'});

ctxtslist = unique(AllNumCtxts);

lt_plot_bar(ctxtslist, ymean, {'Color', 'b'});
lt_plot_bar(ctxtslist, 1./ctxtslist, {'Color', 'k'}); % null

xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);

% ------------
lt_subplot(2,2,3); hold on;
xlabel('numctxts');
ylabel('prob assigned to actual outcome');

plot(AllNumCtxts, AllPihat, 'ok');

[ymean, yN] = grpstats(AllPihat, AllNumCtxts, {'mean', 'numel'});

ctxtslist = unique(AllNumCtxts);

lt_plot_bar(ctxtslist, ymean, {'Color', 'b'});
lt_plot_bar(ctxtslist, 1./ctxtslist, {'Color', 'k'}); % null

xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);

%% === plot with pos and neg contorls
lt_subplot(2,2,4); hold on;
xlabel('numctxts');
ylabel('prob assigned to actual outcome');

% - dat
plot(AllNumCtxts-0.3, AllF1, 'ok');
[ymean, yN] = grpstats(AllF1, AllNumCtxts, {'mean', 'numel'});
lt_plot_bar(ctxtslist-0.3, ymean, {'Color', 'k', 'BarWidth', 0.2});

% - pos
plot(AllNumCtxts, AllF1_pos, 'ob');
[ymean, yN] = grpstats(AllF1_pos, AllNumCtxts, {'mean', 'numel'});
lt_plot_bar(ctxtslist, ymean, {'Color', 'b', 'BarWidth', 0.2});

% -- neg
plot(AllNumCtxts+0.3, AllF1_neg, 'or');
[ymean, yN] = grpstats(AllF1_neg, AllNumCtxts, {'mean', 'numel'});
lt_plot_bar(ctxtslist+0.3, ymean, {'Color', 'r', 'BarWidth', 0.2});


xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);

