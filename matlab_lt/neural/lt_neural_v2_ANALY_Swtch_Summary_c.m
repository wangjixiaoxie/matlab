function lt_neural_v2_ANALY_Swtch_Summary_c(BYNEURONDAT)

% ==== EACH SWITCH SEPARATE PLOT
numswitches = max(BYNEURONDAT.ByNeuron_SwitchCounter);

figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];

hsplots = [];

% --- each switch one plot
for i=1:numswitches
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    ylabel('neural sim change');
    hsplots = [hsplots hsplot];
   inds = find(BYNEURONDAT.ByNeuron_SwitchCounter==i);
    
   if isempty(inds)
       continue
   end
   birdname = BYNEURONDAT.ByNeuron_Birdname{inds(1)};
   exptname = BYNEURONDAT.ByNeuron_Exptname{inds(1)};
   swnum = BYNEURONDAT.ByNeuron_SWnum_real(inds(1));
   title([birdname '-' exptname '-sw' num2str(swnum)]);
   
   func = @(X)nanmedian(X);
   Y = cellfun(func, BYNEURONDAT.ByNeuron_NeuralsimChange(inds,:));
%    Y = cell2mat(Y);
   
   func = @(X)lt_sem(X);
   Ysem = cellfun(func, BYNEURONDAT.ByNeuron_NeuralsimChange(inds,:));
   
   Ysem(isnan(Ysem)) = 0;
   
   plotcols = lt_make_plot_colors(length(inds),0,0);
   
   for j=1:size(Y,1);
    lt_plot(1:3, Y(j,:), {'Errors', Ysem(j,:), 'LineStyle', '-', 'Color', plotcols{j}}); 
   end
    
   xlim([0 4]);
   lt_plot_zeroline;
end

% -------------- COMBINE ACROSS ALL SWITCHES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
   plotcols = lt_make_plot_colors(numswitches,0,0);
for i=1:numswitches
   inds = find(BYNEURONDAT.ByNeuron_SwitchCounter==i);
    
   func = @(X)nanmedian(X);
   Y = cellfun(func, BYNEURONDAT.ByNeuron_NeuralsimChange(inds,:));
   
   func = @(X)lt_sem(X);
   Ysem = cellfun(func, BYNEURONDAT.ByNeuron_NeuralsimChange(inds,:));
   
   Ysem(isnan(Ysem)) = 0;
   
   for j=1:size(Y,1);
       X = [0.8:1:2.8] + 0.4*rand;
%     lt_plot(1:3, Y(j,:), {'Errors', Ysem(j,:), 'LineStyle', '-', 'Color', plotcols{i}, ...
%         'Marker', '.'}); 
    plot(X, Y(j,:), 'o-', 'Color', plotcols{i}); 
   end
    
   xlim([0 4]);
   lt_plot_zeroline;
end
% -- global average
   func = @(X)nanmedian(X);
   Y = cellfun(func, BYNEURONDAT.ByNeuron_NeuralsimChange);
   Ymean = nanmean(Y,1);
   Ysem = lt_sem(Y);
   lt_plot(1.3:3.3, Ymean, {'Errors', Ysem, 'Color', 'k'})

Yall = {};
Yall{1} = Y(:,1);
Yall{2} = Y(:,2);
Yall{3} = Y(:,3);
lt_plot_MultDist(Yall, 1:size(Yall,2), 1, '', 0, 1);
   
xlim([0 4]);
linkaxes(hsplots, 'xy');
