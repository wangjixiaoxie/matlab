function [Params, AllDays_PlotLearning]=lt_seq_dep_pitch_PlotLearning(Params, AllDays_RawDatStruct,saveON);
%% LT 2/6/15 - Given raw dat struct compiled over days (from lt_seq_dep_pitch_DayRawDat) plots learning, and outputs stats struct
% Converted from lt_compile_seq_dep_pitch_data_PLOTDirUndir.
% DIFFERENCE: here plots one structure (e.g. undirected song). Old code
% also compared to directed song, so would load two structures.  Can easily
% plot to compare after this code outputs stats structure



%% PRE-PROCESSING

% Extract params
NumDays=length(AllDays_RawDatStruct); % total days
SylFieldsAll=fieldnames(AllDays_RawDatStruct{1}.data); % all seq and syls (assume 1st day has all possible syls)
FirstDay=Params.SeqFilter.FirstDay;
LastDay=Params.SeqFilter.LastDay;
plotWNdays=Params.PlotLearning.plotWNdays;


% Check other days to make sure SylFieldsAll captures all syls for all days
for i=1:NumDays;
    if ~isempty(AllDays_RawDatStruct{i})
        tmp=fieldnames(AllDays_RawDatStruct{i}.data);
        if length(tmp)>length(SylFieldsAll); % i.e. this day has more syls
            disp(['PROBLEM - day ' num2str(i) ' has more syls than baseline day 1']);
            keyboard
        elseif length(tmp)<length(SylFieldsAll);
            disp(['NOTE - day ' num2str(i) ' has fewer syls than baseline day 1']);
            disp('that day:');
            disp(tmp)
            disp('baseline:');
            disp(SylFieldsAll);
        end
    end
end


% For each syl, is is single syl or sequence?
for i=1:length(SylFieldsAll);
    X(i)=length(SylFieldsAll{i});
end
SylFieldsSingle=SylFieldsAll(X==1); % only single syls
SylFieldsSeq=SylFieldsAll(X>1); % only sequences



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
Params.PlotLearning.SylFieldsAll=SylFieldsAll;
Params.PlotLearning.SylFieldsSingle=SylFieldsSingle;
Params.PlotLearning.SylFieldsSeq=SylFieldsSeq;
Params.PlotLearning.WNTimeOnInd=WNTimeOnInd;
Params.PlotLearning.WNTimeOffInd=WNTimeOffInd;
Params.PlotLearning.DaysToMarkInds=DaysToMarkInds;
Params.PlotLearning.savedir=SaveDir;


%% PROCESS DATA - extract means, etcs

% 1) BASELINE
% although should be in "epochdata" category below, do first here becuase
% needed for day by day data.

for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    X=[];
    PCmat=[];
    for ii=Params.SeqFilter.BaselineDays;
        if ~isempty(AllDays_RawDatStruct{ii}); % check that day has data.
            X=[X cell2mat(AllDays_RawDatStruct{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
            PCmat=[PCmat; cell2mat(AllDays_RawDatStruct{ii}.data.(syl)(:,2))]; % collect PCs over all days
        end
    end

    EpochData.Baseline.(syl).pitchcontour_mean=mean(PCmat,1);
    EpochData.Baseline.(syl).pitchcontour_STD=std(PCmat,0,1);
    EpochData.Baseline.(syl).rawFF=X;
    EpochData.Baseline.(syl).meanFF=mean(X);
    EpochData.Baseline.(syl).stdFF=std(X);
    EpochData.Baseline.(syl).n=length(X);
    EpochData.Baseline.(syl).semFF=std(X)/sqrt((length(X)-1));
end


% 2) DAY-BY-DAY: Extract stats (per day) for all syls into matrix form
for jj=1:length(SylFieldsAll); % how many fields within this set?
    syl=SylFieldsAll{jj}; % actual syl name (e.g. 'a')
    
    Y={};
    Ymean=nan(NumDays,1);
    Yerr=nan(NumDays,1);
    Ystd=nan(NumDays,1);
    
    for i=1:NumDays;
        if ~isempty(AllDays_RawDatStruct{i}); % day has data?
            if isfield(AllDays_RawDatStruct{i}.summary_stats,syl); % check that day has that specific syl
                
                Yvals{i}=AllDays_RawDatStruct{i}.data.(syl)(:,1)'; % raw vals
                Ymean(i)=AllDays_RawDatStruct{i}.summary_stats.(syl).meanFF;
                Yerr(i)=AllDays_RawDatStruct{i}.summary_stats.(syl).semFF;
                Ystd(i)=AllDays_RawDatStruct{i}.summary_stats.(syl).sdFF;
            end
        end
    end
    
    DataMatrix.(syl).FFvals=Yvals;
    DataMatrix.(syl).meanFF=Ymean;
    DataMatrix.(syl).semFF=Yerr;
    DataMatrix.(syl).sdFF=Ystd;
end


% 3) Collect DEVIATION FROM BASELINE and Z-SCORE (individual days)
for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    DataMatrix.(syl).meanFF_DevFromBase=DataMatrix.(syl).meanFF-EpochData.Baseline.(syl).meanFF;
    DataMatrix.(syl).meanFF_zscore=DataMatrix.(syl).meanFF_DevFromBase./EpochData.Baseline.(syl).stdFF; % zscore = deviation from baseline mean divided by baseline std.
    DataMatrix.(syl).meanFF_PercentFromBase=DataMatrix.(syl).meanFF_DevFromBase./EpochData.Baseline.(syl).meanFF; % learning as percent of baseline FF.
end


% 4) GENERALIZATION - Transform learning scores to generalization scores (single days)
% For every day divide every syl learning by target syl learning.

% First, check how many target syls there are.
% Output: data matrix (meanFF-baseline) over all days: DataMatrix_TargDevBase
NumTargs=length(Params.SeqFilter.SylLists.TargetSyls);
if NumTargs>1; % then take average of the targets;
    X=[];
    for i=1:length(NumTargs);
        targsyl=Params.SeqFilter.SylLists.TargetSyls{i};
        X(:,i)=DataMatrix.(targsyl).meanFF_DevFromBase;
        Xz(:,i)=DataMatrix.(targsyl).meanFF_zscore;
    end
    DataMatrix_Targ.DevBaseFF=nanmean(X,2);
    DataMatrix_Targ.DevBaseZ=nanmean(Xz,2);
else % if there;s only one target
    targsyl=Params.SeqFilter.SylLists.TargetSyls{1};
    DataMatrix_Targ.DevBaseFF=DataMatrix.(targsyl).meanFF_DevFromBase;
    DataMatrix_Targ.DevBaseZ=DataMatrix.(targsyl).meanFF_zscore;
end


% Second, calculate generalizations (from all targs, and mean over targs)
for i=1:NumTargs;
    targsyl=Params.SeqFilter.SylLists.TargetSyls{i};
    
    for ii=1:length(SylFieldsAll);
        syl=SylFieldsAll{ii};
        
        DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingHz=...
            DataMatrix.(syl).meanFF_DevFromBase./DataMatrix.(targsyl).meanFF_DevFromBase;
        
        DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingZ=...
            DataMatrix.(syl).meanFF_zscore./DataMatrix.(targsyl).meanFF_zscore;
    end
end

% to average target
for ii=1:length(SylFieldsAll);
    syl=SylFieldsAll{ii};
    
    DataMatrix.(syl).GeneralizationFrom.MeanOverTargs.UsingHz=...
        DataMatrix.(syl).meanFF_DevFromBase./DataMatrix_Targ.DevBaseFF;
    
    DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingZ=...
        DataMatrix.(syl).meanFF_zscore./DataMatrix_Targ.DevBaseZ;
end




% COLLECT RUNNING AVG (multiday bins) over all days, not just WN
RunWind=Params.PlotLearning.DayBinSize;
RunWinInds=NumDays-RunWind+1; % how many windows are there, sliding, size constant (e.g. 5 days, win 3, gives 3 Inds).
DayBinsFieldname=['WindSize_' num2str(Params.PlotLearning.DayBinSize)];
for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    % Sliding window across all days - collect data and get statistics
    for nn=1:RunWinInds; % how many data bins are there (overlapping).
        windowInds=nn:nn-1+RunWind; % e.g. [5 6 7] for 1st bin, if first WN day is 5.
        
        X=[];
        for i=windowInds;
            try % missing data?
                X=[X cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,1))']; % collect FF values over all days into one matrix;
            catch err
            end
        end
        
        try % sometimes all days in bin doesn't have data.
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).rawFF{nn}=X;
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF(nn)=mean(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).stdFF(nn)=std(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).n(nn)=length(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).semFF(nn)=std(X)/sqrt((length(X)-1));
        catch err
        end
        
        % 3) DEVIATION FROM BASELINE
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline(nn)=...
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF(nn)-EpochData.Baseline.(syl).meanFF;
        
        %z-score
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).Zscore(nn)=...
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline(nn)/EpochData.Baseline.(syl).stdFF;
        
        % 4) NOTE down which windowed bins correspond to first and last day
        % of WN
        
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).FirstWNInd=WNTimeOnInd;
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).LastWNInd=WNTimeOffInd-RunWind+1;
    end
end

Params.PlotLearning.DayBinsFieldname=DayBinsFieldname;



% Convert learning scores to generalization scores (divide by learning
% at target
%
% for k=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{k}; % target;
%     for kk=1:length(SylFieldsAll);
%         syl=SylFieldsAll{kk};
%         for nn=1:RunWinInds
%
%             % UNDIR
%             % calculate
%             EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
%                 EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF_minusBaseline/EpochData.WNdaysSlidingWin{nn}.UNDIR.(targsyl).meanFF_minusBaseline;
%
%             EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
%                 EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).Zscore/EpochData.WNdaysSlidingWin{nn}.UNDIR.(targsyl).Zscore;
%

%             % DIR
%             EpochData.WNdaysSlidingWin{nn}.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
%                 EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF_minusBaseline/EpochData.WNdaysSlidingWin{nn}.DIR.(targsyl).meanFF_minusBaseline;
%
%             EpochData.WNdaysSlidingWin{nn}.DIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
%                 EpochData.WNdaysSlidingWin{nn}.DIR.(syl).Zscore/EpochData.WNdaysSlidingWin{nn}.DIR.(targsyl).Zscore;
%
%         end
%     end
% end


% CONVERT Data from sets of fields into matrices for those fields:
% e.g. for fields in order, one matrix might be (for the motif accbb, has 5 syls in order)
% columns: a, c, c, b, b,
% rows: days (start at top)

% Do this for both 3-day binned data and single-day data (here)
% 1) For 3-day binned
RunWinInds=length(EpochData.AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).meanFF); % num bins
SylCategoriesFields=fieldnames(Params.SeqFilter.SylLists); % how manyways of dividing up syls have I used?

for ll=1:length(SylCategoriesFields);
    fn=SylCategoriesFields{ll}; % e.g. fn = 'SylsSame'
    if ~isempty(Params.SeqFilter.SylLists.(fn));
        if ischar(Params.SeqFilter.SylLists.(fn){1})==1; % if this is char, then this is the end of the line for this structure, so start using these entries as the syl names. Otherwise, go one lower (below) to find syl names.
            FieldsList=Params.SeqFilter.SylLists.(fn);
            for ii=1:length(FieldsList);
                syl=FieldsList{ii};
                
                % 1) learning (hz min baseli)
                EpochData.MatrixOverDaysforSylLists.(fn).meanFF_minusBaseline(:,ii)=...
                    EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline';
                
                %             % 2) generalization
                %             for j=1:length(SylLists.TargetSyls);
                %                 targsyl=SylLists.TargetSyls{j};
                %
                %                 EpochData.MatrixOverDaysforSylLists.FieldsInOrder.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                %                     EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                %
                %             end
                
                
                % % FOR SINGLE DAYS
                DataMatrix.MatrixOverDaysforSylLists.(fn).meanFF_minusBaseline(:,ii)=...
                    DataMatrix.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                
                
                
            end
            EpochData.MatrixOverDaysforSylLists.(fn).FieldNamesInOrder=FieldsList;
            DataMatrix.MatrixOverDaysforSylLists.(fn).FieldNamesInOrder=FieldsList;
            
            
        else
            for i=1:length(Params.SeqFilter.SylLists.(fn)); % how many lists? (e.g. fields in order might have 2, but same syls is always 1)
                FieldsList=Params.SeqFilter.SylLists.(fn){i};
                for ii=1:length(FieldsList);
                    syl=FieldsList{ii};
                    
                    % 1) learning (hz min baseli)
                    EpochData.MatrixOverDaysforSylLists.(fn){i}.meanFF_minusBaseline(:,ii)=...
                        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline';
                    
                    %             % 2) generalization
                    %             for j=1:length(SylLists.TargetSyls);
                    %                 targsyl=SylLists.TargetSyls{j};
                    %
                    %                 EpochData.MatrixOverDaysforSylLists.FieldsInOrder.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                    %                     EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                    %
                    %             end
                    
                    % % FOR SINGLE DAYS
                    DataMatrix.MatrixOverDaysforSylLists.(fn){i}.meanFF_minusBaseline(:,ii)=...
                        DataMatrix.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                    
                    
                    
                end
                EpochData.MatrixOverDaysforSylLists.(fn){i}.FieldNamesInOrder=FieldsList;
                DataMatrix.MatrixOverDaysforSylLists.(fn){i}.FieldNamesInOrder=FieldsList;
                
            end
        end
    end
end



%% PLOT - Day means over learning
figure; hold on;

for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % Plot figure for this set of fields
    subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),1,j); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % Plot
        shadedErrorBar(1:length(DataMatrix.(syl).meanFF),DataMatrix.(syl).meanFF,DataMatrix.(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1)
        h(jj)=plot(DataMatrix.(syl).meanFF,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);

%         h=errorbar(DataMatrix.(syl).meanFF,DataMatrix.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%         errorbar_tick(h,200);
        
    end
    % annotate
    legend(h,FieldsList)
    title(['Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    % WN lines
    Fn_AnnotateWNLines(plotWNdays,ylim)
end


% PLOT pitch deviation and zscore - each with own axis.
figure; hold on;
hfig1=[];
hfig2=[];
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % PLOT MEAN PITCH SUBTRACTING BASELINE
    subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % compile pitch means and SEM for each day, for this syl.
        shadedErrorBar(1:length(DataMatrix.(syl).meanFF_DevFromBase),DataMatrix.(syl).meanFF_DevFromBase,DataMatrix.(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1)
        h(jj)=plot(DataMatrix.(syl).meanFF_DevFromBase,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
%         h=errorbar(DataMatrix.(syl).meanFF_DevFromBase,DataMatrix.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%         errorbar_tick(h,200);
    end
    
    legend(h,FieldsList)
    title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    % PLOT Z-score
    subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % Plot
        plot(DataMatrix.(syl).meanFF_zscore,'-','Color',plot_colors{jj},'LineWidth',1.5)
        h(jj)=plot(DataMatrix.(syl).meanFF_zscore,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);        
    end
    
    % annotate
    legend(h,FieldsList)
    title(['Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (Z-scored to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline;
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end


% PLOT pitch deviation and zscore - making yaxis similar
figure; hold on;
hfig1=[];
hfig2=[];
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % PLOT MEAN PITCH SUBTRACTING BASELINE
    hfig1(j)=subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % compile pitch means and SEM for each day, for this syl.
        shadedErrorBar(1:length(DataMatrix.(syl).meanFF_DevFromBase),DataMatrix.(syl).meanFF_DevFromBase,DataMatrix.(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1)
        h(jj)=plot(DataMatrix.(syl).meanFF_DevFromBase,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
%         h=errorbar(DataMatrix.(syl).meanFF_DevFromBase,DataMatrix.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%         errorbar_tick(h,200);
    end
    
    legend(h,FieldsList)
    title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    % PLOT Z-score
    hfig2(j)=subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % Plot
        plot(DataMatrix.(syl).meanFF_zscore,'-','Color',plot_colors{jj},'LineWidth',1.5)
        h(jj)=plot(DataMatrix.(syl).meanFF_zscore,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);        
    end
    
    % annotate
    legend(h,FieldsList)
    title(['Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (Z-scored to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline;
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end


linkaxes(hfig1,'y');
linkaxes(hfig2,'y');

% PLOT GENERALIZATION
% for kk=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{kk};
%
%     % Line plot over days.
%     % One plot for dir and one for undir
%     figure; hold on;
%     for j=1:length(SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
%         FieldsList=SylLists.FieldsToPlot{j};
%         plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%
%
%         subplot(length(SylLists.FieldsToPlot),2,2*j-1); hold on;
%
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             if strcmp(syl,targsyl)==0; % don't plot the target.
%                 % Plot
%                 plot(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%             end
%         end
%         % annotate
%         legend(FieldsList)
%         title(['UNDIR: Generalization (using Hz) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
%         ylabel('Generalization (learning (hz) divided by target syl learning)')
%         xlabel('days')
%         ylim([-1 1]);
%
%         lt_plot_zeroline;
%
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%
%
%
%         % SAME, but with z-score
%         subplot(length(SylLists.FieldsToPlot),2,2*j); hold on;
%
%         % UNDIR - Plot figure for this set of fields
%         for jj=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%
%             if strcmp(syl,targsyl)==0; % don't plot the target.
%                 % Plot
%                 plot(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ ,'o-','Color','k','MarkerFaceColor',plot_colors{jj},'MarkerSize',12);
%             end
%         end
%         % annotate
%         legend(FieldsList)
%         title(['UNDIR: Generalization (Using Z-score) from target : ' targsyl ', ' num2str(FirstDay) ' to ' num2str(LastDay)]);
%         ylabel('Generalization (learning (hz) divided by target syl learning)')
%         xlabel('days')
%         ylim([-1 1]);
%
%         lt_plot_zeroline;
%
%         Fn_AnnotateWNLines(plotWNdays,ylim)
%     end
% end
%
%



%% PLOT SYLLABLES IN SEQUENCE ORDER

   % 1st and last WN bins     
FirstWNBin=EpochData.AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).FirstWNInd;
LastWNBin=EpochData.AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).LastWNInd;

for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    figure; hold on;
    
    X=length(Params.SeqFilter.SylLists.FieldsInOrder{i}); % how many syls?
    Y=EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline;
    Ymin=min(min(Y));
    Ymax=max(max(Y));
    Ylimits{i}=[Ymin-10 Ymax+10];
    Ylimits_zoom{i}=Ylimits{i}./5;
    
    for ii=1:X;
        subplot(X,2,(-1+(X-ii+1)*2)); hold on;
        bar(Y(:,ii));
        ylim(Ylimits{i})
        ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder{ii});
        
        % annotate WN bins
        line([FirstWNBin-0.5 FirstWNBin-0.5], ylim);
        line([LastWNBin-0.5 LastWNBin-0.5], ylim);
        
        
        if ii==X;
            xlabel('Day bin #');
        end
    end
    
    
    for ii=1:X;
        subplot(X,2,((X-ii+1)*2)); hold on;
        bar(Y(:,ii));
        ylim(Ylimits_zoom{i})
        ylabel(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder{ii});
        if ii==X;
            xlabel('Day bin #');
        end
        
                % annotate WN bins
        line([FirstWNBin-0.5 FirstWNBin-0.5], ylim);
        line([LastWNBin-0.5 LastWNBin-0.5], ylim);

    end
    
    subtitle('Mean FF in 3-day bin (minus baseline, hz)');
end



% PLOT EACH DAY INDIVIDUALLY

for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{i};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    figure; hold on;
    for ii=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{ii}; % actual syl name (e.g. 'a')
        subplot(length(FieldsList),2,-1+(length(FieldsList)-ii+1)*2); hold on;
        
        % Plot pitch means and SEM for each day, for this syl.
        bar(DataMatrix.(syl).meanFF_DevFromBase);
        
        ylim(Ylimits{i})
        ylabel(syl);
        if ii==length(FieldsList);
            xlabel('3-day bin #');
        end
        
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
    
    % ZOOM IN
    for ii=1:length(FieldsList);
        syl=FieldsList{ii}; % actual syl name (e.g. 'a')
        subplot(length(FieldsList),2,(length(FieldsList)-ii+1)*2); hold on;
        bar(DataMatrix.(syl).meanFF_DevFromBase);
        
        ylim(Ylimits_zoom{i})
        ylabel(syl);
        if ii==length(FieldsList);
            xlabel('3-day bin #');
        end
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
    
    subtitle('UNDIR - (Syls in order) Mean learning per day (Hz from baseline)');
end




%% PLOT HEAT MAP

% LEARNING (Hz from baseline) - 3 day bins


figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    % WN Lines
    line(xlim,[FirstWNBin-0.5 FirstWNBin-0.5])
    line(xlim,[LastWNBin-0.5 LastWNBin-0.5])

    %     global title
    subtitle('Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end

% ZOOM IN
figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(RunWind)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    caxis(Ylimits_zoom{i})
    
    line(xlim,[FirstWNBin-0.5 FirstWNBin-0.5])
    line(xlim,[LastWNBin-0.5 LastWNBin-0.5])

    %     global title
    subtitle('ZOOM: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end





% PLOT SAME, but day by day
figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    % WN start and end line
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('Learning (mean Hz minus baseline) for each syllable.');
end

% ZOOM IN

figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    caxis(Ylimits_zoom{i})
    
    % WN start and end line
    
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('ZOOM: Learning (mean Hz minus baseline) for each syllable.');
end





%% PLOT MEAN OF LAST FEW DAYS IN ORDER OF SYLLABLES (and snapshot days, if desired)


% PLOT LEARNING AS MINUS BASELINE (HZ);
binfield=Params.PlotLearning.DayBinsFieldname;
lastday=EpochData.AllDaysSliding.(binfield).(syl).LastWNInd;

% only perform if there is data for last day
if ~isnan(EpochData.AllDaysSliding.(binfield).(syl).meanFF_minusBaseline(lastday))
    
    Y={};
    Ysem={};
    for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        % Get all syls
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        Y{ll}=[];
        Ysem{ll}=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            Y{ll}=[Y{ll} EpochData.AllDaysSliding.(binfield).(syl).meanFF_minusBaseline(lastday)];
            Ysem{ll}=[Ysem{ll} EpochData.AllDaysSliding.(binfield).(syl).semFF(lastday)];
        end
        
        Ymax(ll)=max(Y{ll});
        Ymin(ll)=min(Y{ll});
        
    end
    
    
    % Plot
    figure; hold on;
    
    for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        subplot(length(Params.SeqFilter.SylLists.FieldsInOrder),1,ll); hold on;
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
        
        bar(Y{ll})
        
        ylim([min(Ymin)-20 max(Ymax)+20]);
        
        % annotate
        Xtick=1:length(FieldsList); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',FieldsList,'FontSize',12,'FontWeight','bold')
        ylabel('FF (hz), minus baseline','FontSize',12,'FontWeight','bold')
        
    end
    [~,h]=subtitle(['Learning (Hz minus baseline) last 3 WN days']);
    set(h,'FontSize',12,'FontWeight','bold')
else
    disp('Note, not plotting summary of last WN day because lacks data!!');
end


%% TO DO: CONVERT TO FUNCTION TO PLOT ANY DESIRED DAYS
% IF WANT TO PLOT A SNAPSHOT - need to designate in params
% Params.SeqFilter.DaysForSnapshot{1}={'09Dec2014','11Dec2014'};
% Params.SeqFilter.DaysToMark= {'11Dec2014-2400'}; % will mark all plots with lines here;

% only suppports one entry for days ATM.

if isfield(Params.SeqFilter,'DaysForSnapshot');
    Params.PlotLearning.fnamedays=(['from' Params.SeqFilter.DaysForSnapshot{1}{1} 'to' Params.SeqFilter.DaysForSnapshot{1}{2}]); % field name for structure below
    
    [A]=lt_convert_EventTimes_to_RelTimes(FirstDay,Params.SeqFilter.DaysForSnapshot{1});
    
    fday=A.JustDays_rel(1);
    lday=A.JustDays_rel(2);
    
    for i=1:length(SylFieldsAll);
        syl=SylFieldsAll{i};
        
        X=[];
        for ii=fday:lday; % over the days
            if ~isempty(AllDays_RawDatStruct{ii});
            X=[X; cell2mat(AllDays_RawDatStruct{ii}.data.(syl)(:,1))];
            else % if empty, then tell user
                disp(['day ' num2str(ii) ' lacks data, but is trying to be used in Snapshot']);
            end
        end
        
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).FFvals=X;
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF=mean(X);
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).stdFF=std(X);
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).n=length(X);
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).sem=std(X)/sqrt(length(X)-1);
        
        % DEVIATION FROM BASELINE
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline=...
            EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF-EpochData.Baseline.(syl).meanFF;
        
        %z-score
        EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).Zscore=...
            EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline/EpochData.Baseline.(syl).stdFF;
        
        %             % GENERALIZATION STATS
        %         for i=1:length(SylFieldsAll);
        %             syl=SylFieldsAll{i};
        %
        %             for k=1:length(SylLists.TargetSyls);
        %                 targsyl=SylLists.TargetSyls{k}; % target;
        %
        %                 EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
        %                     EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).meanFF_minusBaseline;
        %
        %                 EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
        %                     EpochData.Snapshot.(fnamedays).UNDIR.(syl).Zscore/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).Zscore;
        %             end
        %         end
        
    end

        % NOW PLOT
        % First, gather data in plotting format for BAR
        Y={};
        Ysem={};
        for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
            % Get all syls
            FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            Y{ll}=[];
            Ysem{ll}=[];
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                Y{ll}=[Y{ll} EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline];
                Ysem{ll}=[Ysem{ll} EpochData.Snapshot.(Params.PlotLearning.fnamedays).(syl).sem];
            end
            
            Ymax(ll)=max(Y{ll});
            Ymin(ll)=min(Y{ll});
            
        end
        
        
        % Second, Plot
        figure; hold on;
        
        for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
            subplot(length(Params.SeqFilter.SylLists.FieldsInOrder),1,ll); hold on;
            FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
            
            bar(Y{ll})
            
            ylim([min(Ymin)-20 max(Ymax)+20]);
            
            % annotate
            Xtick=1:length(FieldsList); % one tick for each syl. needed.
            set(gca,'XTick',Xtick);
            set(gca,'XTickLabel',FieldsList,'FontSize',12,'FontWeight','bold')
            ylabel('FF (hz), minus baseline','FontSize',12,'FontWeight','bold')
            
        end
        [~,h]=subtitle(['Learning (Hz minus baseline) from ' Params.SeqFilter.DaysForSnapshot{1}{1} ' to ' Params.SeqFilter.DaysForSnapshot{1}{2}]);
        set(h,'FontSize',12,'FontWeight','bold')


end


%% PLOT PITCH CONTOUR OVER LEARNING
DaysToPlot=[];
% ---- Get subset of days to plot pitch contour of
% last baseline day
DaysToPlot(1)=Params.SeqFilter.BaselineDays(end);

% 4 days into learning
DaysToPlot(2)=Params.SeqFilter.BaselineDays(end)+4;

% last day of learning
DaysToPlot(3)=Params.PlotLearning.WNTimeOffInd;

% any specified day to mark
if isfield(Params.PlotLearning,'DaysToMarkInds');
    if ~isempty(Params.PlotLearning.DaysToMarkInds);
        DaysToPlot=[DaysToPlot cell2mat(Params.PlotLearning.DaysToMarkInds)];
    end
end

% reorder
DaysToPlot=sort(DaysToPlot);
% ---

plotcolors=lt_make_plot_colors(length(DaysToPlot),1,[1 0.2 0.2]);

% PLOT
for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{j};
    
    %     plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    figure; hold on;
    
    for jj=1:length(FieldsList); % how many fields within this set?
        
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        subplot(1,length(FieldsList),jj); hold on;
        title(['Syl: ' syl]);
        
        % extract baseline PC mean and std
        baselinePC_mean=EpochData.Baseline.(syl).pitchcontour_mean;
        baselinePC_std=EpochData.Baseline.(syl).pitchcontour_STD;
        
        
        for jjj=1:length(DaysToPlot);
            dayInd=DaysToPlot(jjj);
            
            if ~isempty(AllDays_RawDatStruct{dayInd})
                
                % - calcualte deviation from baseline PC
                PC_minusBase=AllDays_RawDatStruct{dayInd}.summary_stats.(syl).pitchcountour_mean...
                    -baselinePC_mean;
                PC_zscore=PC_minusBase./baselinePC_std;
                
                % -- Plot
                
                try
                    % mean PC
%                     plot(AllDays_RawDatStruct{dayInd}.summary_stats.(syl).pitchcountour_mean,'Color',plotcolors{jjj},'LineWidth',2);
%                     text(1,AllDays_RawDatStruct{dayInd}.summary_stats.(syl).pitchcountour_mean(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                    
                    % deviation of PC from baseline
%                     plot(PC_minusBase,'Color',plotcolors{jjj},'LineWidth',2);
%                     text(1,PC_minusBase(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                   
                    % zscore of PC
                    plot(PC_zscore,'Color',plotcolors{jjj},'LineWidth',2);
                    text(1,PC_zscore(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                   
                catch err
                    disp('problem with pitch contour, had to use try');
                end
                
            end
        end
        
    end
    subtitle(['Day mean pitch contour']);
end



%% PLOT HIT RATE OVER LEARNING

% collect hit statuses for all trials
figure; hold on;
title('Hit Rate each day (likely for catch songs) for all syls');
xlabel('day');
ylabel('Fraction of renditions hit');

plotcols=lt_make_plot_colors(length(SylFieldsAll),0,0);

for j=1:length(SylFieldsAll);
    syl=SylFieldsAll{j};
    
    hitrate=nan(1,NumDays); % to collect values to plot
    for i=1:NumDays;
        
        if ~isempty(AllDays_RawDatStruct{i});
            
            % collect data
            DataMatrix.(syl).HitRateStatus{i}=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,9));
            
            hitrate(i)=sum(DataMatrix.(syl).HitRateStatus{i})/length(DataMatrix.(syl).HitRateStatus{i});
        end
        
        % plot
        tmp=length(syl);
        
        if tmp==1;
            hplot(j)=plot(hitrate,'--','Color',plotcols{j},'MarkerFaceColor',plotcols{j});
        elseif tmp>1
            hplot(j)=plot(hitrate,'-o','Color',plotcols{j},'MarkerFaceColor',plotcols{j},'LineWidth',2);
        end
        
    end
end
legend(hplot,SylFieldsAll);


Fn_AnnotateWNLines(plotWNdays,ylim)


%% IS THERE CIRCADIAN FLUCTUATION THAT SHOULD BE SUBTRACTED FROM DATA?

figure; hold on;
hplot=[];
hsplot=[];
syllist=SylFieldsSingle;

plotcols=lt_make_plot_colors(length(syllist),0,0);

for j=1:length(syllist);
    syl=syllist{j};
    
    hsplot(j)=subplot(ceil(sqrt(length(syllist))),ceil(sqrt(length(syllist))),j); hold on;
    title(syl);
    
    % collect values    
    FFvals_tot=[];
    Tvals_tot=[];
    for i=Params.SeqFilter.BaselineDays;
        if ~isempty(AllDays_RawDatStruct{i});
            
            FFvals=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,1));
            tmp=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,6));
            [~, tmp] =lt_convert_datenum_to_hour(tmp);
            Tvals=tmp.hours;
            
            % collect vals across days
            FFvals_tot=[FFvals_tot; FFvals];
            Tvals_tot=[Tvals_tot; Tvals];
        end
    end
    
    % plot
    hplot(j)=plot(Tvals_tot,FFvals_tot,'.','Color',plotcols{j});
    
    % get running avg
    [Tvals_tot, ind]= sort(Tvals_tot);
    FFvals_tot=FFvals_tot(ind);
    
    Tvals_tot_sm=lt_running_stats(Tvals_tot,10);
    FFvals_tot_sm=lt_running_stats(FFvals_tot,10);
    
    plot(Tvals_tot_sm.Mean,FFvals_tot_sm.Mean,'-','Color',plotcols{j},'LineWidth',2);
end

legend(hplot,syllist);
linkaxes(hsplot,'x');
subtitle('Baseline circadian fluctuation of pitch');


%% PLOT TRIAL BY TRIAL

% TO DO - concatenate multiple days
% - plot baseline dashed line
% - overlay targ syl
% - compare across days


DaysToPlot=[16];
SmoothBin=10; % 10 renditions
targsyl=Params.SeqFilter.SylLists.TargetSyls{1};

% PLOT pitch deviation and zscore - making yaxis similar
figure; hold on;
hfig1=[];
hfig2=[];
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
    
        hsplot(jj)=subplot(ceil(length(FieldsList)/2),2,jj); hold on;
        title(syl);
        
% combine raw data across all days
        for jjj=1:length(DaysToPlot);
               dayInd=DaysToPlot(jjj);
               
               FFvals=cell2mat(AllDays_RawDatStruct{dayInd}.data.(syl)(:,1));
               % subtract beginning of day (1st 1/6 of data) from FF values
               tmp=round(length(FFvals)/6);
               baselineFF=mean(FFvals(1:tmp));
               FFvals=FFvals-baselineFF;
               
               % get times
               tmp=cell2mat(AllDays_RawDatStruct{dayInd}.data.(syl)(:,6));
               tmp=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,tmp);
               Tvals=tmp.FinalValue; % get times in form of (1.5 = day 1, noon)
               
               % smooth        
               FFvals_sm=lt_running_stats(FFvals,SmoothBin);
               Tvals_sm=lt_running_stats(Tvals,SmoothBin);
               
               % -- Plot
               shadedErrorBar(Tvals_sm.Mean,FFvals_sm.Mean,FFvals_sm.SEM,{'-','Color',plot_colors{jj}},1);
               plot(Tvals_sm.Mean,FFvals_sm.Mean,'.','Color',plot_colors{jj})
                              
               % plot line for day mean
               plot(xlim,[0 0],'--');
               Xlim=xlim;
               text(Xlim(1),0,num2str(baselineFF));
        end
    end
    linkaxes(hsplot,'xy');
end

               
               
    

%% TO DO:
% at: 
% Compile data into matrix form, for all syllable lists
% (apart from fields in order and fields of interest above

if (0)
UndirOrDirList={'UNDIR','DIR'};
lt_compile_seq_dep_pitch_data_PLOTDirUndir_ComplSylLists
end


%% OUTPUT

AllDays_PlotLearning.EpochData=EpochData;
AllDays_PlotLearning.DataMatrix=DataMatrix;
AllDays_PlotLearning.DataMatrix_Targ=DataMatrix_Targ;


%% SAVE
if saveON==1;
    mkdir(Params.PlotLearning.savedir);
    cd(Params.PlotLearning.savedir)
    
    save('Params','Params');
    save('AllDays_PlotLearning','AllDays_PlotLearning');
    lt_save_all_figs
    
end

end


%% VARIOUS SUNFUNCTIONS

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
