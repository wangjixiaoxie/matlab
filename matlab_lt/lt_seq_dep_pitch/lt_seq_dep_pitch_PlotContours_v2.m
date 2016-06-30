function lt_seq_dep_pitch_PlotContours_v2(AllDays_RawDatStruct, AllDays_PlotLearning, Params, syl, onlyPlotLMANdays);
%% Plots contours - outlier not removed, actual time windows where I used for MUSC analysis - plots LMAN data automatically
% -- Does not use actual time windows of PBS, however!!!!!!!!!!!!!!!!!!
% 1) same syllable, all days individually (with one line for baseline mean) and combined
% 2) same, but for LMAN


%%

syls_unique=Params.PlotLearning.SylFields_Unique;
datafield='data_WithOutlier';

if ~exist('onlyPlotLMANdays', 'var');
    onlyPlotLMANdays=0;
end

%% Plot all days - one subplot for each day
count=1;
SubplotsPerFig=20;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


daystoplot=1:length(AllDays_RawDatStruct);
if onlyPlotLMANdays==1;
    try
        daystoplot=find(~cellfun(@isempty, Params.PlotLearning.MuscimolSchedule_ByDayInds));
    catch err
        disp('error - is this not LMAN data? then don"t ask for just LMAN days');
    end
    
end

Hsplots=[];
day1_mean=[];

for i=daystoplot;
    
    if isempty(AllDays_RawDatStruct{i})
        continue;
    end
    
    PCmat=cell2mat(AllDays_RawDatStruct{i}.(datafield).(syl)(:,2));
    
    % ===== what are boundaries?
    ind=find(strcmp(fieldnames(AllDays_RawDatStruct{i}.data_WithOutlier), syl));
    
    if isfield(Params, 'RecalculateFF'); % then means changed window, e.g. for WN songs - all days will use same
        twindow=Params.RecalculateFF.pc_time_window_list(:, ind);
        
        %         line([twindow(1) twindow(1)], ylim, 'Color','r', 'LineWidth',2);
        %         line([twindow(2) twindow(2)], ylim, 'Color','r', 'LineWidth', 2);
        %
    else
        twindow=Params.SeqFilter.pc_time_window_list{i}(:,ind);
        
    end
    
    
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT
    X=twindow(1)-100:twindow(2)+100;
    X(X<1 | X>size(PCmat,2))=[];
    
    
    % ===== DATA
    [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(['day '  num2str(i)]);
    Hsplots=[Hsplots hsplot];
    
    plot(X, PCmat(:,X)', '--','Color',[0.8 0.8 0.8]);
    
    % plot mean
    Mean1=mean(PCmat,1)';
    % --- save day 1 mean
    if isempty(day1_mean);
        day1_mean= Mean1;
    end
    
    NumRendsPBS=size(PCmat,1);
    
    
    % === PLOT LMAN DATA IF EXISTS
                NumRendsMUSC=[];
    if isfield(AllDays_RawDatStruct{i}, 'data_MUSC');
        if isfield(AllDays_RawDatStruct{i}.data_MUSC, syl);
            if ~isempty(AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{i});
                %         muscsongtimes=AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{i};
                %
                %
                %
                %         alltimes=cell2mat(AllDays_RawDatStruct{i}.data_MUSC.(syl)(:,6));
                %         inds_good_times=alltimes>min(muscsongtimes) & alltimes<max(muscsongtimes);
                
                %         muscsongtimes=AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{i};
                %
                %
                %
                % -- COllect times of all raw data for today (both musc and pbs)
                %         tmp1=cell2mat(AllDays_RawDatStruct{i}.data_MUSC.(syl)(:,6));
                %         tmp2=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,6));
                %         alltimes=[tmp1; tmp2];
                alltimes=cell2mat(AllDays_RawDatStruct{i}.data_WithOutlier.(syl)(:,6));
                [~, WithinDayVals]=lt_convert_datenum_to_hour(alltimes);
                alltimes=WithinDayVals.hours;
                
                % ---- collect all PC mat
                allPCmat=[cell2mat(AllDays_RawDatStruct{i}.data_MUSC.(syl)(:,2)); cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,2))]; % data can leak into PBS
                
                % ---- find times that are in musc window
                mintime=Params.PlotLearning.MuscimolSchedule_ByDayInds{i}.start+Params.PlotLearning.Lag_time;
                endtime=Params.PlotLearning.MuscimolSchedule_ByDayInds{i}.end + Params.PlotLearning.Dur_of_PBS_dat_that_counts_as_MUSC;
                
                inds_good_times=alltimes>mintime & alltimes<endtime;
                
                % -- plot
                if ~isempty(inds_good_times);
%                     if size(inds_good_times,1)==size(allPCmat,1)
                        
                    PCmat=allPCmat(inds_good_times,:);
                    plot(X, PCmat(:, X)', '--r');
                    
                    tmp=mean(PCmat,1)';
                    plot(X, tmp(X),'--','Color','y','LineWidth',2);
                    
                    NumRendsMUSC=size(PCmat,1);
%                     end
                end
            end
        end
    end
    
    % === data mean
    plot(X, Mean1(X),'--','Color','k','LineWidth',2);
    
    % ==== OVERLAY BASELINE MEAN
    try
        PCmean=AllDays_PlotLearning.EpochData.Baseline.(syl).pitchcontour_mean;
        plot(X, PCmean(X), 'b', 'LineWidth',2);
    catch err
        
        disp('COULD NOT PLOT BASELINE MEAN...(plotting earliest day mean instead)');
        plot(X, day1_mean(X), 'b', 'LineWidth',2);
    end
    
    %     % ==== OVERLAY BASELINE MEAN FOR MUSC
    %     AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).pitchcontour_mean
    %
    axis tight
    
    % ===== lines for time window
    line([twindow(1) twindow(1)], ylim, 'Color','k', 'LineWidth', 2);
    line([twindow(2) twindow(2)], ylim, 'Color','k', 'LineWidth', 2);
    
    % NOTE SAMPLE SIZE
    Ylim=ylim;
    Xlim=xlim;
    
    lt_plot_text(Xlim(1), Ylim(2), ['N=' num2str(NumRendsPBS)], 'm');
    if ~isempty(NumRendsMUSC) % for MUSC
        lt_plot_text(Xlim(1), Ylim(2)-100, ['N=' num2str(NumRendsMUSC)], 'k');
        
    end
    
end

lt_subtitle(syl);
linkaxes(Hsplots, 'xy');