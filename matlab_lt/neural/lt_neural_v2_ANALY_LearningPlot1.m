function FirstDay = lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_onebird, MOTIFSTATS)
%% NOTE:
% Summarystruct must match MOTIFSTATS perfectly, i.e. must be one bird only
% 


%% Plot overview of experiment. (time, ffvals for targ and nontargs, and duration of recording for each neuron

motif_regexpr_str = MOTIFSTATS.params.motif_regexpr_str;
TargSyls = MOTIFSTATS.params.TargSyls;

NumNeurons = length(MOTIFSTATS.neurons);
NumSyls = length(motif_regexpr_str);

%%
% time range across all neurons
targnum = find(strcmp(motif_regexpr_str, TargSyls{1}));
DatenumsAll = [];
for i=1:NumNeurons
    datenumstmp = [MOTIFSTATS.neurons(i).motif(targnum).SegmentsExtract.song_datenum];
    
    DatenumsAll = [DatenumsAll datenumstmp];
    
end

DatenumsAll = sort(unique(DatenumsAll));
FirstDay = datestr(DatenumsAll(1), 'ddmmmyyyy');

% --- now plot
lt_figure; hold on;
plotcols = lt_make_plot_colors(NumNeurons, 0, 0);
for i=1:NumNeurons
    datenumstmp = [MOTIFSTATS.neurons(i).motif(targnum).SegmentsExtract.song_datenum];
    ffvals = [MOTIFSTATS.neurons(i).motif(targnum).SegmentsExtract.FF_val];
    
    tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, datenumstmp);
    tvals = tmp.FinalValue;
    
    lt_subplot(2, 1,1); hold on;
    title(TargSyls{1});
    ylabel('ff');
    xlabel('time (days)');
    lt_plot(tvals, ffvals);
    
    % plot line for dur that have data for this neuron
    lt_subplot(2,1,2); hold on;
    title('duration of data');
    ylabel('neuron #');
    line([min(tvals) max(tvals)], [i i], 'Color', plotcols{i}, 'LineWidth', 2);
    ylim([0 i+1]);
    
    % indicate which chan it is
    chan = SummaryStruct_onebird.birds(1).neurons(i).channel;
    clust = SummaryStruct_onebird.birds(1).neurons(i).clustnum;
    lt_plot_text(min(tvals), i+0.1, ['chan' num2str(chan) '-' num2str(clust)], plotcols{i});
end

% -- plot line for start of WN
lt_subplot(2,1,1);
WNon = SummaryStruct_onebird.birds(1).neurons(1).LEARN_WNonDatestr;
try
tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, {WNon});
WNon_days = tmp.FinalValue;
line([WNon_days WNon_days], ylim);
catch err
end


% --- indicate which chan it is

