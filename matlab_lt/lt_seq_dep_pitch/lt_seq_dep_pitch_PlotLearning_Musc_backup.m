function [Params, AllDays_RawDatStruct, DataMatrix]= lt_seq_dep_pitch_PlotLearning_Musc(Params, AllDays_RawDatStruct, DataMatrix)

%% PARAMS
Lag_time=1.7; % how many hours post MUSC start to collect data?
PBS_window=[-1.5 0]; % [-2 0] means take PBS data in window -2hrs to 0.

ExptCondition_codes={'PBS','MUSC'}; % i.e. if filename says "PBS", then is code 1 (musc=2 and so forth)
Params.PlotLearning.ExptCondition_codes=ExptCondition_codes;

%     % convert to codes
%     PBS_code=find(strcmp(ExptCondition_codes, 'PBS'));
%     MUSC_code=find(strcmp(ExptCondition_codes, 'PBS'));

plotcols=lt_make_plot_colors(length(ExptCondition_codes), 0, 0);

SylFields_Unique=Params.PlotLearning.SylFields_Unique; % all seq and syls (assume 1st day has all possible syls)

NumDays=length(AllDays_RawDatStruct);
FirstDay=Params.SeqFilter.FirstDay;
plotWNdays=Params.PlotLearning.plotWNdays;



%% RUN
for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    lt_figure; hold on; grid on;
    title(syl);
    xlabel('days'); ylabel('FF (hz)');
    
    for ii=1:NumDays;
        
        % check if day has data
        if isempty(AllDays_RawDatStruct{ii});
            continue
        end
        
%         % ==== SORT TODAY's TRIALS INTO CONDITIONS ==================
%         % get all filenames
%         filenames=AllDays_RawDatStruct{ii}.data.(syl)(:,5);
%         
%         % get conditions from filenames
%         ExptConditions=[];
%         for iii=1:length(filenames);
%             underscores=strfind(filenames{iii},'_');
%             exptcond=filenames{iii}(underscores(1)+1:underscores(2)-1);
%             
%             % convert that condition (string) to code (1, 2, ...);
%             if any(strcmp(ExptCondition_codes, exptcond));
%                 ExptConditions(iii)=find(strcmp(ExptCondition_codes, exptcond));
%             else
%                 disp('PROBLEM - a song is not PBS or MUSC based on filename');
%             end
%         end
%         
%         
%         
%         % == SAVE CONDITION TO BOTH RAW DAT AND PLOT LEARNING
%         % STRUCTURES
%         % 1) Raw Dat structure
%         AllDays_RawDatStruct{ii}.data.(syl)(:,12)=num2cell(ExptConditions');
%         
%         % update legend in params
%         Params.DayRawDat.Legend_data{1,:}={'1_FF','2_PitchContour','3_SoundDat','4_Spec',...
%             '5_FileName','6_DateNum','7_LabelStr','8_NotePos','9_Trig','10_NoteDur','11_CatchTrial','12_ExptCond'};
%         Params.DayRawDat.ExptCondition_codes=ExptCondition_codes;
%         
%         % 2) put into Data Matrix structure
%         DataMatrix.(syl).SORTED_BY_CONDITIONS.ExptConditionInds_AllData{ii}=ExptConditions;
%         
%         
        
        % === EXTRACT VARIOUS STATS - for this syl, get mean etc, for
        % each day
        
        FFvals_PBS=cell2mat(DataMatrix.(syl).FFvals{ii});
        Tvals_PBS=cell2mat(DataMatrix.(syl).Tvals{ii});
        
        FFvals_MUSC=cell2mat(DataMatrix_MUSC.(syl).FFvals{ii});
        Tvals_MUSC=cell2mat(DataMatrix_MUSC.(syl).Tvals{ii});
        
%         FFvals_all=cell2mat(DataMatrix.(syl).FFvals{ii});
%         Tvals_all=cell2mat(DataMatrix.(syl).Tvals{ii});
        
        
        % == for each condition (e.g. PBS or MUSC), get relevant data (musc
        % is lag after switch, PBS is in window before switch.
      
        for k=1:length(ExptCondition_codes);
            
            
%             ExptInds=find(ExptConditions==k);
%             
%             % -- time and vals for this condition
%             FFvals_condition=FFvals_all(ExptInds);
%             Tvals_condition=Tvals_all(ExptInds);
%             
%             
%             % -- Save those things for this condition
%             DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals{ii}=FFvals_condition;
%             DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.Tvals{ii}=Tvals_condition;
            
            
            % -- For muscimol, get data an hour after start, to account
            % for lag
                if any(strcmp(ExptCondition_codes{k}, 'MUSC'));
                    
                    if ~isempty(FFvals_MUSC)
                    
                    % When did the switch happen?
                    MUSC_start=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
                    
                    
                    % sort data by time
                    [Tvals_condition_sorted, inds] = sort(Tvals_condition);
                    FFvals_condition_sorted=FFvals_condition(inds);
                    
                    % get the data coming after the lag
                    % convert time to hours
                    [~, DataTimes] = lt_convert_datenum_to_hour(Tvals_condition_sorted);
                    
                    IndsToKeep=find(DataTimes.hours>MUSC_start+Lag_time);
                    
                    % == extract the inds that pass lag time
                    Tvals_WithinTimeWindow=Tvals_condition_sorted(IndsToKeep);
                    FFvals_WithinTimeWindow=FFvals_condition_sorted(IndsToKeep);
                    
                    
                    % == save those values
                    DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals_WithinTimeWindow{ii}=FFvals_WithinTimeWindow;
                    DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.Tvals_WithinTimeWindow{ii}=Tvals_WithinTimeWindow;
                    
                    end
                end
            end
            
            
            % -- FOR PBS, TAKE DATA IN 2 HOUR WINDOW PRECEDING SWITCH
            
            if ~isempty(ExptInds); % if this condition has data for this day
                if any(strcmp(ExptCondition_codes{k}, 'PBS'));
                    
                    % When did the switch happen?
                    if isempty(Params.PlotLearning.MuscimolSchedule_ByDayInds{ii});
                        % if no switch today, then use the median
                        % switch time
                        MUSC_start=Params.PlotLearning.MuscimolSchedule_MedianStartTime;
                    else
                        MUSC_start=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
                    end
                    
                    % sort data by time
                    [Tvals_condition_sorted, inds] = sort(Tvals_condition);
                    FFvals_condition_sorted=FFvals_condition(inds);
                    
                    % get the data coming in window 2 hrs before switch
                    % convert time to hours
                    [~, DataTimes] = lt_convert_datenum_to_hour(Tvals_condition_sorted);
                    
                    IndsToKeep=find(DataTimes.hours>MUSC_start+PBS_window(1) & DataTimes.hours<MUSC_start+PBS_window(2));
                    
                    % == extract the inds that pass lag time
                    Tvals_WithinTimeWindow=Tvals_condition_sorted(IndsToKeep);
                    FFvals_WithinTimeWindow=FFvals_condition_sorted(IndsToKeep);
                    
                    
                    % == save those values
                    DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals_WithinTimeWindow{ii}=FFvals_WithinTimeWindow;
                    DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.Tvals_WithinTimeWindow{ii}=Tvals_WithinTimeWindow;
                end
            end
            
            
            % === PLOT LEARNING, comparing muscimol to no muscimol
            if ~isempty(ExptInds);
                
                % == 1) All data each syl a plot --> mult days (for both conditions)
                Y=DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals{ii};
                X=DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.Tvals{ii};
                
                % convert time to days
                tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X);
                X=tmp.FinalValue;
                
                plot(X, Y, 'o', 'Color', plotcols{k});
                
                % plot mean
                Ymean=mean(Y);
                Ysem=lt_sem(Y);
                errorbar(ceil(max(X))-rand/5, Ymean, Ysem, 's', 'color', plotcols{k}, 'MarkerSize', 9);
                
                % -- OVERLAY with 1) data within windows for analysis, and 2)
                % lines indication swiches
                Y_windowed=DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals_WithinTimeWindow{ii};
                X_windowed=DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.Tvals_WithinTimeWindow{ii};
                
                % convert time to days
                if ~isempty(X_windowed)
                    tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X_windowed);
                    X_windowed=tmp.FinalValue;
                    
                    plot(X_windowed, Y_windowed, 'o', 'Color', plotcols{k}, 'MarkerFaceColor', plotcols{k});
                    
                    % mean
                    Ymean_windowed=mean(Y_windowed);
                    Ysem_windowed=lt_sem(Y_windowed);
                    errorbar(ceil(max(X))-rand/5, Ymean_windowed, Ysem_windowed, 's', 'color', plotcols{k}, 'MarkerFaceColor', plotcols{k}, 'MarkerSize', 9);
                end
                
                % -- Put lines indicating experiment switches
                if ~isempty(Params.PlotLearning.MuscimolSchedule_ByDayInds{ii})
                    Xline=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
                    Xline=ii+Xline/24; % convert to days;
                    line([Xline Xline], ylim);
                    
                    Xline=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.end;
                    Xline=ii+Xline/24; % convert to days;
                    line([Xline Xline], ylim);
                end
                
            end
        end
        
    end
end


% ======== PLOT LEARNING IN OTHER WAYS - MUSCIMOL
% FOR EACH SYL, PLOT JUST MEANS, SEMS, and CV

for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    lt_figure; hold on; grid on;
    xlabel('days'); ylabel('FF (hz)');
    
    for k=1:length(ExptCondition_codes);
        
        Xall=[];
        Yall.mean=[];
        Yall.sem=[];
        Yall.cv=[];
        Yall.N=[];
        
        for ii=1:NumDays;
            
            if isempty(DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals_WithinTimeWindow{ii});
                continue
            end
            
            Yvals=DataMatrix.(syl).SORTED_BY_CONDITIONS.CONDITION_NUM_data{k}.FFvals_WithinTimeWindow{ii};
            
            Yall.mean=[Yall.mean mean(Yvals)];
            Yall.sem=[Yall.sem lt_sem(Yvals)];
            Yall.N=[Yall.N length(Yvals)];
            Ystd=std(Yvals);
            Yall.cv=[Yall.cv Ystd/mean(Yvals)];
            
            Xall=[Xall ii];
            
            
            
        end
        % == Plot all days
        hsplot(1)=lt_subplot(4,1,1:2); hold on;
        
        shadedErrorBar(Xall, Yall.mean, Yall.sem, {'Color', plotcols{k}},1);
        lt_plot(Xall, Yall.mean, {'Color', plotcols{k}});
        ylabel('FF (hz) (SEM)');
        
        % == Plot N and CV
        hsplot(2)=lt_subplot(4,1,3); hold on;
        lt_plot(Xall, Yall.cv, {'Color', plotcols{k}});
        ylabel('CV');
        
        hsplot(3)=lt_subplot(4,1,4); hold on;
        lt_plot(Xall, Yall.N, {'Color', plotcols{k}});
        ylabel('N');
        
    end
    
    subplot(4,1,1:2);
    Fn_AnnotateWNLines(plotWNdays, ylim);
    
    % == annotate days that should throw out for muscmiol (e.g.
    % concentration)
    Ylim=ylim;
    plot(Params.PlotLearning.MuscimolDaysToThrowout_Inds, Ylim(2), 'xk', 'MarkerSize', 8);
    
    
    lt_subtitle(syl);
    linkaxes(hsplot,'x');
end


% === PLOT RELATIVE TO BASELINE, AND TAKING INTO ACCOUNT BASELINE
% EFFECTS OF MUSCIMOL
end

function Fn_AnnotateWNLines(plotWNdays,ylim)

global WNTimeOnInd
global WNTimeOffInd
global DaysToMarkInds

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd+0.5 WNTimeOffInd+0.5],ylim,'LineStyle','--','Color','r')
    
    for i=1:length(DaysToMarkInds);
        line([DaysToMarkInds{i} DaysToMarkInds{i}],ylim,'LineStyle','--','Color','k');
    end
end
end
