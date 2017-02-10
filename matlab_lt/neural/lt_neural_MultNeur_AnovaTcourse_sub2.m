function lt_neural_MultNeur_AnovaTcourse_sub2(OUTSTRUCT, transition_type, ...
    shufflechoice, predur, plotdist, timeWindRelSylOnset)

%% ------------------
NumNeurons = length(OUTSTRUCT.neuron);

if ~exist('timeWindRelSylOnset', 'var')
    timeWindRelSylOnset = [0 0];
end

%%
% +++++++++++++++++++++++ 1) COLLECT ALL PVALS
Pval_all = {};
Xval = [];
for nn=1:NumNeurons
    
    numsyls = length(OUTSTRUCT.neuron(nn).(transition_type).syl);
    % ===== EACH BRANCH POINT IS ONE DATAPOINT
    for ss = 1:numsyls
        
        % make sure this syl has data
        if isempty(OUTSTRUCT.neuron(nn).(transition_type).syl(ss).ANOVA)
            continue
        end
        
        % ===== COLLECT DATA
        % 1) pr of omega2 (actual data) relative to shuffled
        Pval = OUTSTRUCT.neuron(nn).(transition_type).syl(ss).ANOVA.(shufflechoice).pr_shuff_greaterthan_data;
        Pval_all = [Pval_all Pval];
        
        Xval = OUTSTRUCT.neuron(nn).convergent.syl(ss).ANOVA.(shufflechoice).Xvals;
        
    end
end

% ---- PLOT
xmin = min(cellfun(@length, Pval_all));
% put into matrix
Pval_all_mat = nan(length(Pval_all), xmin);
for j=1:length(Pval_all)
    Pval_all_mat(j, :) = Pval_all{j}(1:xmin);
end

if plotdist==0
    % line plot + box
    plot(Xval(1:xmin), Pval_all_mat, '-', 'Color', [0.3 0.7 0.8]);
    boxplot(Pval_all_mat, 'position', Xval(1:xmin), 'color', 'k')
    axis tight; line([predur predur], ylim, 'Color', 'b');
elseif plotdist ==1
    % heat map
    distributionPlot(Pval_all_mat,'colormap', gray, 'showMM', 5, 'variableWidth', false)
    set(gca, 'XTick', 1:4:xmin);
    set(gca, 'XTickLabel', Xval(1:4:xmin))
    ylim([0 1]);
elseif plotdist ==2
    % hist for given time window
    if timeWindRelSylOnset(1) == timeWindRelSylOnset(2)
        disp('ERROR!!!!!!!!!!! specify timeWindRelSylOnset');
    end
    % -- which inds to take
    x1 = predur + timeWindRelSylOnset(1);
    x2 = predur + timeWindRelSylOnset(2);
    inds = find(Xval>x1 & Xval<x2);
    
    % -- average over time bins
    Pval_hist = mean(Pval_all_mat(:, inds), 2);
    
    % -- plot histogram
    Xcenters = 0.0125:0.0125:0.9875;
    lt_plot_histogram(Pval_hist, Xcenters, 1, 1, 1, 1, 'r');
    xlim([0 1]);
elseif plotdist ==3
    % running histgarm
    Xcenters = 0.0125:0.0125:0.9875;
    lt_plot_mult2dhist(Pval_all_mat, Xval(1:xmin), Xcenters);
    
end
