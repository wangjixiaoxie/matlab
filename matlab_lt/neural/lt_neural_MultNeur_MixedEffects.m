function lt_neural_MultNeur_MixedEffects(SYLNEURDAT)

%% STOPPED ON: multiple model section below


%% ==== rewrite variables

numrends = length(SYLNEURDAT.data_sylrend);

y = [SYLNEURDAT.data_sylrend.SpikeRate_hz]';
y_spks = [SYLNEURDAT.data_sylrend.NumSpikes]';

y_spks_log = y_spks;
y_spks_log(y_spks_log ==0) = 1;
y_spks_log = log10(y_spks_log);

y_postfiring = [SYLNEURDAT.data_sylrend.SpikeRate_hz_postmotor]';

tmp = [SYLNEURDAT.data_sylrend.neuronID]';
x_neuron = {};
for i=1:length(tmp)
    x_neuron = [x_neuron; num2str(tmp(i))];
end
% x_neuron = num2str([SYLNEURDAT.data_sylrend.neuronID]');
x_syllable = [SYLNEURDAT.data_sylrend.Syl]';
x_syldur = [SYLNEURDAT.data_sylrend.syldur]';

x_presyl = [];
for i=1:numrends
    x_presyl = [x_presyl; SYLNEURDAT.data_sylrend(i).pre_ten_syls(end)];
end

x_postsyl = [];
for i=1:numrends
    if isempty(SYLNEURDAT.data_sylrend(i).post_ten_syls)
        x_postsyl = [x_postsyl; '-'];
    else
        x_postsyl = [x_postsyl; SYLNEURDAT.data_sylrend(i).post_ten_syls(1)];
    end
end


x_FF = [SYLNEURDAT.data_sylrend.FF]';

X = table(y, y_spks, y_postfiring, y_spks_log, x_neuron, x_syllable, x_presyl, x_FF, x_postsyl, ...
    x_syldur);

%% ==== subsample, for 

%% ===

NumNeurons = length(unique(X.x_neuron));


%% fit model

% === MODEL:
modelform = 'y ~ 1  + x_neuron + x_presyl + x_FF';
lme = fitlme(X, modelform);

% TO DO: why rank deficient? - lack data for certain transitions - what to
% do then?

% PICK OUT THE TRANSITIONS THAT DO EXIST, AND CALCULATE USING ONLY THOSE
% BRANCH POINTS.

% === IS POST SYL FIRING DEPENDENT ON POST SYL?
modelform = 'y_postfiring ~ 1 + x_syllable + (1 + x_syllable|x_postsyl) + (1 + x_syllable|x_neuron)';
lme = fitlme(X, modelform);


%% ====== CONTEXT MATTERS EVEN AFTER ACCOUNTING FOR FF?
% --- CONVERGENT POINTS ONLY
close all
branchmaster = 'v';
Xlim = [min(X.x_FF)-20 max(X.x_FF)+20];
Ylim = [min(X.y) max(X.y)];

plotSpksLog = 1; % if 0, then plots hz

for nn=1:NumNeurons
    figcount=1;
    subplotrows=3;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    
    % === 1) extract a specific branch
    % mastersyl = 'h';
    % X_onebranch = fn_extract_branch(X, x_syllable, mastersyl);
    
    filter_any = {'x_syllable', branchmaster, 'x_neuron', num2str(nn)}; % ordered pair {column_name, target};
    X_onebranch = fn_extract_anything(X, filter_any);
    
    % === 2) plot to visualize FF differences between diff contexts
    presyl_list = unique(X_onebranch.x_presyl);
    
    plotcols = lt_make_plot_colors(length(presyl_list), 0, 0);

    hsplots = [];
    % -- plot each on separate plot
    for i=1:length(presyl_list)
        presyl = presyl_list(i);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['presyl: ' presyl]);
        hsplots = [hsplots hsplot];
        
        inds = findstr(X_onebranch.x_presyl', presyl);
        
        xx_ff = X_onebranch.x_FF(inds);
        if plotSpksLog ==1
        yy_hz = X_onebranch.y_spks_log(inds);
        else
        yy_hz = X_onebranch.y(inds);
        end
        lt_regress(yy_hz, xx_ff, 1, 0, 1, 1, plotcols{i});
        
    end
    
    % -- plot each on one plot (points only);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all one plot']);
    hsplots = [hsplots hsplot];
    
    for i=1:length(presyl_list)
        presyl = presyl_list(i);
        
        inds = findstr(X_onebranch.x_presyl', presyl);
        
        xx_ff = X_onebranch.x_FF(inds);
        if plotSpksLog ==1
        yy_hz = X_onebranch.y_spks_log(inds);
        else
        yy_hz = X_onebranch.y(inds);
        end
        plot(xx_ff,yy_hz, '.', 'Color', plotcols{i});
    end
    
    
    linkaxes(hsplots, 'xy');
    xlim(Xlim); ylim(Ylim);
    
    % ---- PLOT EACH DISTRIBUTION ON ONE PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all one plot']);
    if plotSpksLog ==1
        ylabel('log10(numspikes)');
    end
    
    yvals = {};
    for i=1:length(presyl_list)
        presyl = presyl_list(i);
        
        inds = findstr(X_onebranch.x_presyl', presyl);
        
        if plotSpksLog ==1
        yy_hz = X_onebranch.y_spks_log(inds);
        else
        yy_hz = X_onebranch.y(inds);
        end
        
        yvals = [yvals yy_hz];
        if length(yy_hz)>4
        distributionPlot(yy_hz, 'xValues', i, 'color', plotcols{i});
        end
    end
%         distributionPlot(yvals, 'xValues', 1:length(presyl_list), 'color', plotcols, 'showMM', 6);

    lt_subtitle(['neuron ' num2str(nn)]);
    
end

% % === 3) fit model
% filter_any = {'x_syllable', branchmaster, 'x_neuron', '1'}; % ordered pair {column_name, target}
% X_fit = fn_extract_anything(X, filter_any);
% modelform = 'y ~ 1 + x_presyl + x_FF';

filter_any = {'x_syllable', branchmaster}; % ordered pair {column_name, target}
X_fit = fn_extract_anything(X, filter_any);
modelform = 'y ~ 1 + x_presyl + x_FF + (1+x_presyl+x_FF|x_neuron)';

lme = fitlme(X_fit, modelform);

% compare 2 models


%% ======== ACCOUNTING FOR SYL DUR

% --- CONVERGENT POINTS ONLY
branchmaster = 'h';
Xlim = [min(X.x_syldur)-0.01 max(X.x_syldur)+0.01];
Ylim = [min(X.y_spks)-2 max(X.y_spks)+2];

for nn=1:NumNeurons
    figcount=1;
    subplotrows=3;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    
    % === 1) extract a specific branch
    filter_any = {'x_syllable', branchmaster, 'x_neuron', num2str(nn)}; % ordered pair {column_name, target};
    X_onebranch = fn_extract_anything(X, filter_any);
    
    % === 2) plot to visualize FF differences between diff contexts
    presyl_list = unique(X_onebranch.x_presyl);
    plotcols = lt_make_plot_colors(length(presyl_list), 0, 0);

    hsplots = [];
    % -- plot each on separate plot
    for i=1:length(presyl_list)
        presyl = presyl_list(i);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['presyl: ' presyl]);
        hsplots = [hsplots hsplot];
        
        inds = findstr(X_onebranch.x_presyl', presyl);
        
        xx_syldur = X_onebranch.x_syldur(inds);
        yy_spks = X_onebranch.y_spks(inds);
        lt_regress(yy_spks, xx_syldur, 1, 0, 1, 1, plotcols{i});
    end
    
    % -- plot each on one plot (points only);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all one plot']);
     xlabel('syl dur (s)');
    ylabel('numspks');
   hsplots = [hsplots hsplot];
    
    for i=1:length(presyl_list)
        presyl = presyl_list(i);
        
        
        inds = findstr(X_onebranch.x_presyl', presyl);
        
        xx_syldur = X_onebranch.x_syldur(inds);
        yy_spks = X_onebranch.y_spks(inds);
        plot(xx_syldur, yy_spks, '.', 'Color', plotcols{i});
    end
    
    linkaxes(hsplots, 'xy');
    lt_subtitle(['neuron ' num2str(nn)]);
    xlim(Xlim); ylim(Ylim);
    
end



%% =========== MODEL (num spikes, controlling for FF and syl dur)

% --- IGNORE THIS WAY -- filter to get multiple syls
inds = regexp(X.x_syllable', '[hbv]');
X_fit = X(inds, :);

modelform = 'y_spks ~ 1 + (x_presyl*x_syllable|x_neuron)';
lme1 = fitlme(X_fit, modelform);

% --- filter for just one thing
branchmaster = 'h';
filter_any = {'x_syllable', branchmaster}; % ordered pair {column_name, target}
X_fit = fn_extract_anything(X, filter_any);

modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (x_presyl+x_FF+x_syldur|x_neuron)';
lme0 = fitlme(X_fit, modelform);


% 1) include all random effects
modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (1 + x_presyl + x_FF + x_syldur|x_neuron)';
lme1 = fitlme(X_fit, modelform);

% 2
modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (1 + x_presyl + x_FF|x_neuron)';
lme2 = fitlme(X_fit, modelform);

modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (1 + x_presyl + x_syldur|x_neuron)';
lme3 = fitlme(X_fit, modelform);

modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (1 + x_presyl|x_neuron)';
lme4 = fitlme(X_fit, modelform);

modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur + (1|x_neuron)';
lme5 = fitlme(X_fit, modelform);

modelform = 'y_spks ~ 1 + x_presyl + x_FF + x_syldur';
lme6 = fitlme(X_fit, modelform);

modelform = 'y_spks ~ 1 + (1 + x_FF + x_syldur|x_presyl) + (1|x_neuron)';
lme7 = fitlme(X_fit, modelform);



%% ====== DO ONE NEURON AT ATIME [IMPORTANT]

syl_list_1 = {'h', 'b', 'v'};
DATSTRUCT = struct;
p_val_matrix = nan(NumNeurons, length(syl_list_1));
p_val_matrix_FF =  nan(NumNeurons, length(syl_list_1));
p_val_matrix_syldur = nan(NumNeurons, length(syl_list_1));
for i=1:NumNeurons
    for j=1:length(syl_list_1)
        syl = syl_list_1{j};
        
        filter_any = {'x_syllable', syl, 'x_neuron', num2str(i)};
        X_fit = fn_extract_anything(X, filter_any);
        
        % --- version 1 (fixed effects)
        if (0)
            modelform = 'y_spks_log ~ 1 + x_presyl + x_FF + x_syldur';
            lme0 = fitlme(X_fit, modelform);
            
            DATSTRUCT(i,j).lme = lme0;
            p_val_matrix(i, j) = lme0.anova.pValue(2);
            p_val_matrix_FF(i,j) = lme0.anova.pValue(3);
            p_val_matrix_syldur(i,j) = lme0.anova.pValue(4);
        end
        
        % --- version 2, x_presyl random effect on intercept, FF, and
        % syldur
        if (1)
            modelform = 'y_spks_log ~ 1 + x_FF +x_syldur + (1+x_FF + x_syldur|x_presyl)';
            lme0 = fitlme(X_fit, modelform);
            
            modelform = 'y_spks_log ~ 1 + x_FF +x_syldur';
            lme0_smaller = fitlme(X_fit, modelform);
            
            % --- compare models
            results = compare(lme0_smaller, lme0); %
            
            % -- output
            DATSTRUCT(i,j).lme = lme0;
            p_val_matrix(i, j) = results.pValue(2);
            p_val_matrix_FF(i, j) = lme0.anova.pValue(2);
            p_val_matrix_syldur(i, j) = lme0.anova.pValue(3);
        end
    end
end


% ==== PLOT P-VALUE [presyl]
alpha = 0.05;
lt_figure; hold on;
for j=1:length(syl_list_1)
    lt_subplot(3,1,j);
    
    y = p_val_matrix(:,j);
    lt_plot_bar(1:length(y), log10(y));
    line(xlim, [log10(alpha) log10(alpha)])
        
    title(['branch ' syl_list_1{j}]);
    xlabel('neuron id');
    ylabel('log10(pval) (presyl0');
end


% ==== plot pval [FF]
lt_figure; hold on;
for j=1:length(syl_list_1)
    lt_subplot(3,1,j);
    
    y = p_val_matrix_FF(:,j);
    lt_plot_bar(1:length(y), log10(y));
    line(xlim, [log10(alpha) log10(alpha)])
        
    title(['branch ' syl_list_1{j}]);
    xlabel('neuron id');
    ylabel('log10(pval) (FF)');
end

% ==== plot pval [syldur]
lt_figure; hold on;
for j=1:length(syl_list_1)
    lt_subplot(3,1,j);
    
    y = p_val_matrix_syldur(:,j);
    lt_plot_bar(1:length(y), log10(y));
    line(xlim, [log10(alpha) log10(alpha)])
        
    title(['branch ' syl_list_1{j}]);
    xlabel('neuron id');
    ylabel('log10(pval) (syldur)');
end


% ===================== 


%% ====== RESTRICT TO specific branch point

% branchtype = 'conv';
% mastersyl = 'h';
%
% inds = strfind(x_syllable', mastersyl);
%
% X_onebranch = X(inds, :);
%

mastersyl = 'h';
X_onebranch = fn_extract_branch(X, x_syllable, mastersyl);

% ====== MODEL:
modelform = 'y ~ 1  + x_neuron + x_presyl + x_FF';
lme = fitlme(X_onebranch, modelform);



%% ======= BRANCH points with conv + div (firing better predicted by conv or div)

% -- find all presyls relative to a given syl
targsyl = 'h';
inds = strfind(x_syllable', targsyl);
potential_presyls = unique(x_presyl(inds));

% ==== 1) fit model to ask about presyls
X_temp = X(inds, :);

modelform = 'y ~ 1 + x_presyl + (1+x_presyl|x_neuron)';
modelform = 'y ~ 1 + x_presyl + (1|x_neuron)';
lme = fitlme(X_temp, modelform);


% ==== 2) for each presyl, ask how much its following neural activity
% is predicted by it, or by its following syl
for i=1:length(potential_presyls)
    
    presyl = potential_presyls(i);
    
    if presyl =='-'
        continue
    end
    
    inds = strfind(x_syllable', presyl);
    potentialpost_syls = unique(x_postsyl(inds));
    
    if length(potentialpost_syls)<2
        continue
    end
    
    % ----- FIT MODEL
    X_temp = X(inds, :);
    modelform = 'y_postfiring ~ 1 + x_postsyl + (1+x_postsyl|x_neuron)';
    modelform = 'y_postfiring ~ 1 + x_postsyl + (1|x_neuron)';
    lme = fitlme(X_temp, modelform);
    
    lt_figure; hold on;
    
    
end




branch_regexp = '[][]';

% ==== TO DO:
% remove '-'
% -- plot bar plot of firing rate for each presyl/postsyl - overlay
% statistics
end


%% === functions

function X_onebranch = fn_extract_branch(X, x_syllable, mastersyl)
inds = strfind(x_syllable', mastersyl);

X_onebranch = X(inds, :);
end

function X_filtered = fn_extract_anything(X, filter_any)
% like extract_branch, but makes that obsolete.

num_filters = length(filter_any)/2;

for i=1:num_filters
    col_name = filter_any{2*i-1};
    targ = filter_any{2*i};
    
    if strcmp(col_name, 'x_neuron')
        inds = strcmp(X.(col_name), targ);
    else
    inds = X.(col_name) == targ;
    end
    X = X(inds, :);
end
X_filtered = X;
end







