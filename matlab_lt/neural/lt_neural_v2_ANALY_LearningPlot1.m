function FirstDay = lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_onebird, MOTIFSTATS)
%% NOTE:
% Summarystruct must match MOTIFSTATS perfectly, i.e. must be one bird only
% ALSO, MUST BE ONE EXPT ONLY


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
hsplots = [];
for i=1:NumNeurons
    
    hsplot = lt_subplot(3, 1,1:2); hold on;
    hsplots = [hsplots hsplot];
    plotcols_targs = lt_make_plot_colors(length(TargSyls), 0, 0);
    for ii=1:length(TargSyls)
        targind = find(strcmp(motif_regexpr_str, TargSyls{ii}));
        datenumstmp = [MOTIFSTATS.neurons(i).motif(targind).SegmentsExtract.song_datenum];
        ffvals = [MOTIFSTATS.neurons(i).motif(targind).SegmentsExtract.FF_val];
        
        % ----------
        tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, datenumstmp);
        tvals = tmp.FinalValue;
        
        %     title(TargSyls{1});
        plot(tvals, ffvals, 'o', 'Color', plotcols_targs{ii});
        if i == NumNeurons
            ylabel('ff');
            lt_plot_text(tvals(end)+0.04, ffvals(end), TargSyls{ii}, plotcols_targs{ii});
        end
    end
%     axis tight
    
    % plot line for dur that have data for this neuron
    hsplot = lt_subplot(3,1,3); hold on;
    hsplots = [hsplots hsplot];
    %     title('duration of data');
    ylabel('neuron #');
    xlabel('time (days)');
    line([min(tvals) max(tvals)], [i i], 'Color', plotcols{i}, 'LineWidth', 2);
    ylim([0 i+1]);
    
    % indicate which chan it is
    chan = SummaryStruct_onebird.birds(1).neurons(i).channel;
    clust = SummaryStruct_onebird.birds(1).neurons(i).clustnum;
    lt_plot_text(max(tvals)+0.01, i+0.1, ['chan' num2str(chan) '-' num2str(clust)], plotcols{i});
end
linkaxes(hsplots, 'x');

% -- put lines for expt
lt_subplot(3,1,1:2); hold on;
tmp = lt_neural_v2_LoadLearnMetadat;
indbird = strcmp({tmp.bird.birdname}, SummaryStruct_onebird.birds(1).birdname);
indexpt = strcmp(tmp.bird(indbird).info(1,:), SummaryStruct_onebird.birds(1).neurons(1).exptID);

for j=1:length(TargSyls)
indtmp = indexpt & strcmp(tmp.bird(indbird).info(2,:), TargSyls{j});
Transitions = tmp.bird(indbird).info(3:end,indtmp);
Transitions = Transitions(~cellfun('isempty', Transitions));

for i=1:length(Transitions)
    transtime = Transitions{i}(1:14);
    transtime = lt_convert_EventTimes_to_RelTimes(FirstDay, datenum(transtime, 'ddmmmyyyy-HHMM'));
    transtime  = transtime.FinalValue;
    line([transtime transtime], ylim, 'Color', 'k');
    precond = Transitions{i}(16:17);
    postcond = Transitions{i}(19:20);
    Ylim = ylim;
    lt_plot_text(transtime-0.01, (1-0.02*j)*Ylim(2), [precond '-' postcond],...
        plotcols_targs{j});
%     lt_plot_text(transtime+0.02, (1-0.02*j)*Ylim(2), postcond, plotcols_targs{j});
end
end
% % --- now plot
% lt_figure; hold on;
% plotcols = lt_make_plot_colors(NumNeurons, 0, 0);
%
% for i=1:NumNeurons
%     datenumstmp = [MOTIFSTATS.neurons(i).motif(targnum).SegmentsExtract.song_datenum];
%     ffvals = [MOTIFSTATS.neurons(i).motif(targnum).SegmentsExtract.FF_val];
%
%
%     % ----------
%     tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, datenumstmp);
%     tvals = tmp.FinalValue;
%
%     lt_subplot(3, 1,1:2); hold on;
%     title(TargSyls{1});
%     ylabel('ff');
%     lt_plot(tvals, ffvals);
%
%     % plot line for dur that have data for this neuron
%     lt_subplot(3,1,3); hold on;
% %     title('duration of data');
%     ylabel('neuron #');
%     xlabel('time (days)');
%     line([min(tvals) max(tvals)], [i i], 'Color', plotcols{i}, 'LineWidth', 2);
%     ylim([0 i+1]);
%
%     % indicate which chan it is
%     chan = SummaryStruct_onebird.birds(1).neurons(i).channel;
%     clust = SummaryStruct_onebird.birds(1).neurons(i).clustnum;
%     lt_plot_text(min(tvals), i+0.1, ['chan' num2str(chan) '-' num2str(clust)], plotcols{i});
% end

% subplotsqueeze(gcf, 1.1);

% -- plot line for start of WN
% lt_subplot(2,1,1);
% WNon = SummaryStruct_onebird.birds(1).neurons(1).LEARN_WNonDatestr;
% try
%     tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, {WNon});
%     WNon_days = tmp.FinalValue;
%     line([WNon_days WNon_days], ylim);
% catch err
% end


% --- indicate which chan it is

