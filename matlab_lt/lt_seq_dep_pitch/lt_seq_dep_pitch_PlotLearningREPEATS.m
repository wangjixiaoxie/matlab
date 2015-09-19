function [Params, AllDays_PlotLearning]=lt_seq_dep_pitch_PlotLearningREPEATS(Params, AllDays_RawDatStruct,saveON)
% LT 3/31/15 - 
% Run after lt_seq_dep_pitch_SeqFilterCompile. Run either this function (if
% working with repeats) or lt_seq_dep_pitch_PlotLearning (if sequence
% variants, but not repeats).



%% PREPROCESSING

% Extract Params
NumDays=length(AllDays_RawDatStruct); % total days
FirstDay=Params.SeqFilter.FirstDay;
LastDay=Params.SeqFilter.LastDay;
plotWNdays=Params.PlotLearning.plotWNdays;


% Convert WNdays to index days
global WNTimeOnInd % make global, so can use in subfunction below.
global WNTimeOffInd

if Params.PlotLearning.plotWNdays==1;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeON});
    WNTimeOnInd=X.JustDays_rel;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeOFF});
    WNTimeOffInd=X.JustDays_rel;
end


% Get days to mark, if exist.
global DaysToMarkInds
DaysToMarkInds={};
X={};
isfield(Params.SeqFilter,'DaysToMark')
if isfield(Params.SeqFilter,'DaysToMark'); % if there are specific days to average over and look at.
    DaysToMark=Params.SeqFilter.DaysToMark;
    for i=1:length(DaysToMark);
        X{i}=lt_convert_EventTimes_to_RelTimes(FirstDay,{DaysToMark{i}});
        DaysToMarkInds{i}=X{i}.JustDays_rel;
    end
end


% For saving
timestampSv=lt_get_timestamp(0);
SaveDir=[Params.SeqFilter.savedir '/PlotLearning_' timestampSv];


% UPDATE PARAMS STRUCTURE
Params.PlotLearning.WNTimeOnInd=WNTimeOnInd;
Params.PlotLearning.WNTimeOffInd=WNTimeOffInd;
Params.PlotLearning.DaysToMarkInds=DaysToMarkInds;
Params.PlotLearning.savedir=SaveDir;


%% RUN, if have Params.Repeats (i.e. using old method)

if isfield(Params.SeqFilter, 'Repeats');
NumDays=length(AllDays_RawDatStruct);

figure; hold on;
PlotCols=lt_make_plot_colors(6,0);
for i=1:NumDays;
    SylFields=Params.SeqFilter.SylLists.FieldsInOrder{1};
    for ii=1:length(SylFields);
        sylfield=SylFields{ii};
        if isfield(AllDays_RawDatStruct{i}.summary_stats,sylfield);
            
         hfig(ii)=errorbar(i,AllDays_RawDatStruct{i}.summary_stats.(sylfield).meanFF,AllDays_RawDatStruct{i}.summary_stats.(sylfield).semFF,'o','Color','k','MarkerFaceColor',PlotCols{ii},'MarkerSize',8);
        end
    end
         
end

legend(hfig,SylFields);

Fn_AnnotateWNLines(plotWNdays,ylim)
end

%% FIRST, EXTRACT RELEVANT REGULAR EXPRESSION BASED DATA FROM 

%% RUN, USING REGULAR EXPRESSIONS METHOD


% === PLOT ABOVE,
for i=1:length(Params.SeqFilter.RegExpr.expressions);
    regexpr_string=Params.SeqFilter.RegExpr.expressions{i};
    
    NumRepeatClasses=length(AllDays_Repeats.RegExpr_string{i}.RepClass);
    for ii=1:NumRepeatClasses;
        
        % === One plot across days for each repeat class
        repeatclass=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.string;
        
        lt_figure; hold on;
        title([regexpr_string '; ' repeatclass])
        
        for j=1:NumDays;
            
            if isempty(AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.FFvals{j})
                continue;
            end
            
            Tvals=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.Tvals_ThisRepClass{j};
            FFvals=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.FFvals{j};
            
            % convert Tvals to days
            eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
            Tvals=eventtimes.FinalValue;
            
            plot(Tvals, FFvals, 'o');
            
            % == PLOT day mean and SEM
            FFmean=mean(FFvals);
            FFsem=lt_sem(FFvals);
            
            % remove nan
            FFmean(isnan(FFmean))=[];
            FFsem(isnan(FFsem))=[];
            
            errorbar(ceil(Tvals(1))*ones(length(FFmean),1), FFmean, FFsem, 'o');
            
            
            % ==== PLOT bar graph, to look at change relative to syl repeat
            % number
            
            
        end
        
        
        
        SylLabels=AllDays_Repeats.RegExpr_string{i}.RepClass{ii};
        
        
    end
end

end

end

%% VARIOUS SUBFUNCTIONS
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