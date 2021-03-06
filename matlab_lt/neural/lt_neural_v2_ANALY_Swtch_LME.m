function lt_neural_v2_ANALY_Swtch_LME(DATSylMot)
useThreeSylTypes = 1;
%% LT 8/9/17 -

fnames = fieldnames(DATSylMot);
for i=1:length(fnames)
    DATSylMot.(fnames{i}) = DATSylMot.(fnames{i})';
end

%% modify predictor variables

tmp = cell(length(DATSylMot.AllBirdnum),1);
if useThreeSylTypes ==1
% targ same diff, as categorical

[tmp{logical(DATSylMot.AllIsTarg)}] = deal('targ');
[tmp{logical(DATSylMot.AllIsTarg==0 & DATSylMot.AllIsSame==1)}] = deal('same');
[tmp{logical(DATSylMot.AllIsTarg==0 & DATSylMot.AllIsDiff==1)}] = deal('diff');

else
 
[tmp{logical(DATSylMot.AllIsTarg)}] = deal('targ');
[tmp{logical(DATSylMot.AllIsTarg==0)}] = deal('nontarg');

end
assert(sum(logical(DATSylMot.AllIsTarg)) + sum(DATSylMot.AllIsTarg==0 & DATSylMot.AllIsSame==1) + ...
    sum(DATSylMot.AllIsTarg==0 & DATSylMot.AllIsDiff==1) ==length(DATSylMot.AllBirdnum), 'asdasfdf');

DATSylMot.SylType = tmp;

% change in 
tmp = DATSylMot.AllNeurSplitCorrTrain - DATSylMot.AllNeurSplitCorrBase;
DATSylMot.AllNeurSplitCorrChange = tmp;

% 
DATSylMot.AllBaseFFvsNeurFRCorr = double(DATSylMot.AllBaseFFvsNeurFRCorr);

%% take care of noisy data?


%% fit model
tbl = struct2table(DATSylMot);

% formula = 'AllNeurSplitCorrTrain ~ 1 + SylType + (1|AllNeurSplitCorrBase)';

% 0) most basic (using split data)
% formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + AllBaseFFvsNeurFRCorr +', ...
%     '(1|AllNeurTargLearnRate_targdir) + (1|AllBirdnum) + (1|AllBirdnum:AllSwitchnumGlobal)'];
% formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + AllBaseFFvsNeurFRCorr +', ...
%     '(1|AllNeurTargLearnRate_targdir) + (1|AllBirdnum)'];
% formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + AllBaseFFvsNeurFRCorr +', ...
%     'AllNeurTargLearnRate_targdir + (1|AllBirdnum)'];
formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + AllSNRbaseEnd +', ...
    'AllNeurTargLearnRate_targdir + (1|AllBirdnum)'];
% formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + ', ...
%     'AllNeurTargLearnRate_targdir + (1|AllBirdnum)'];
% formula = ['AllNeurSplitCorrChange ~ 1 + AllIsTarg + AllNeurSplitCorrBase', ...
%     ' + (1|AllBirdnum)'];
mdl = fitlme(tbl, formula)

% --- plot fixed effects coefficients
lt_figure; hold on;

    ymean = mdl.Coefficients.Estimate;
    ysem = mdl.Coefficients.SE;
    lt_plot_bar(1:length(ymean), ymean, {'Errors', ysem})

    % -- x labels
    set(gca, 'XTickLabel', mdl.Coefficients.Name')
    rotateXLabels(gca, 45);
    
for j=1:length(mdl.Coefficients.Estimate)
   
    % -- put p value
    p = mdl.Coefficients.pValue(j);
    lt_plot_text(j, 1.1*mdl.Coefficients.Estimate(j), ['p=' num2str(p, '%3.2g')], 'r');
    
end

% %  1) using split data
% formula = ['AllNeurSplitCorrChange ~ 1 + SylType + AllNeurSplitCorrBase + AllSNRbaseEnd +' ...
%     'AllBaseFFvsNeurFRCorr+ (1|AllBirdnum) + (1|AllBirdnum:AllSwitchnumGlobal)'];
% mdl = fitlme(tbl, formula)

% 2) using neural sim
formula = ['AllNeurSimChange ~ 1 + SylType + AllNeurSimBase + AllSNRbaseEnd +' ...
    'AllBaseFFvsNeurFRCorr + (1|AllBirdnum) + (1|AllBirdnum:AllSwitchnumGlobal)'];
% indstoexclude = find(isnan(DATSylMot.AllBaseFFvsNeurFRCorr));
mdl = fitlme(tbl, formula)