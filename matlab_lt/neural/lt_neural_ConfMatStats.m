function Stats = lt_neural_ConfMatStats(ConfMat)
%% lt 8/23/17 - modified to take in cell array, where each cell is a conf mat. can be diff dimensions in diff cells

%% lt 8/17/17 - takes in confusion matrix (contingency table) and outputs stats

% ConfMat - n x n matrix, rows = actual class; columns = predictions

%%


%%
clear Stats
if iscell(ConfMat)
    
    numcells = numel(ConfMat);
    for i=1:numcells
        
        cmat = ConfMat{i};
        
        sts = getstats(cmat);
        
        Stats(i) = sts;
    end
    
    
else
    
    Stats = getstats(ConfMat);
    
%     accuracy = trace(ConfMat)./sum(ConfMat(:));
%     
%     tmp = repmat(sum(ConfMat, 2), 1, size(ConfMat, 2));
%     sensitivity = diag(ConfMat./tmp);
%     
%     tmp = repmat(sum(ConfMat, 1), size(ConfMat,1), 1);
%     precision = diag(ConfMat./tmp);
%     
%     
%     F1 = (2*sensitivity.*precision)./(precision+sensitivity);
%     
%     % convert nan to 0
%     F1(isnan(F1)) = 0;
%     
%     % -- take means
%     sensitivity = mean(sensitivity);
%     precision = mean(precision);
%     F1 = mean(F1);
%     
%     Stats.sensitivity = sensitivity;
%     Stats.precision = precision;
%     Stats.F1 = F1;
end

end

function Stats = getstats(ConfMat)

accuracy = trace(ConfMat)./sum(ConfMat(:));

tmp = repmat(sum(ConfMat, 2), 1, size(ConfMat, 2));
sensitivity = diag(ConfMat./tmp);

tmp = repmat(sum(ConfMat, 1), size(ConfMat,1), 1);
precision = diag(ConfMat./tmp);


F1 = (2*sensitivity.*precision)./(precision+sensitivity);

% convert nan to 0
F1(isnan(F1)) = 0;

% -- take means
sensitivity = mean(sensitivity);
precision = mean(precision);
F1 = mean(F1);

% ---- one vs. others stats (allows for combining between different numbers
% of classes
numclasses = size(ConfMat,1);
OneVAllMat = nan(2,2,size(ConfMat,1));
for i=1:numclasses
   
    tmpmat = nan(2,2);
    
    tmpmat(1,1) = ConfMat(i,i); % TP
    tmpmat(1,2) = sum(ConfMat(i, [1:i-1 i+1:end])); % FN
    tmpmat(2,2) = sum(diag(ConfMat)) - ConfMat(i,i); % TN
    tmpmat(2,1) = sum(ConfMat(:)) - tmpmat(1,1) - tmpmat(1,2) - tmpmat(2,2); % FP
    
    OneVAllMat(:,:,i) = tmpmat;
    
end
OneVAllMat = sum(OneVAllMat,3);

% ------ get stats for OneVAllMat
accuracy2 = trace(OneVAllMat)./sum(OneVAllMat(:));

tmp = repmat(sum(OneVAllMat, 2), 1, size(OneVAllMat, 2));
sensitivity2 = diag(OneVAllMat./tmp);

tmp = repmat(sum(OneVAllMat, 1), size(OneVAllMat,1), 1);
precision2 = diag(OneVAllMat./tmp);

F1_2 = (2*sensitivity2.*precision2)./(precision2+sensitivity2);



% -------- OUT
Stats.accuracy = accuracy;
Stats.sensitivity = sensitivity;
Stats.precision = precision;
Stats.F1 = F1;

Stats.OneVAll_accuracy = accuracy2;
Stats.OneVAll_sensitivity = sensitivity2;
Stats.OneVAll_precision = precision2;
Stats.OneVAll_F1 = F1_2;






end