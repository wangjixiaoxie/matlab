function Stats = lt_neural_ConfMatStats(ConfMat)

%% lt 8/17/17 - takes in confusion matrix (contingency table) and outputs stats

% ConfMat - n x n matrix, rows = actual class; columns = predictions

%%


%%

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

Stats.sensitivity = sensitivity;
Stats.precision = precision;
Stats.F1 = F1;
