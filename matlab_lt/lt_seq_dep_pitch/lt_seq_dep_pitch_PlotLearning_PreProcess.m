function [EpochData, DataMatrix, DataMatrix_Targ, Params] = lt_seq_dep_pitch_PlotLearning_PreProcess(Params, AllDays_RawDatStruct, SylFields_Unique, muscimol_data)
%% called from lt_seq_dep_pitch_PlotLearning

% muscimol_data=1; then uses data_MUSC, not data.

%% PARAMS

% default - not muscimol
if ~exist('muscimol_data','var');
    muscimol_data=0;
end

% if muscimol, then use muscimol data, otherwise PBS data.
if muscimol_data==1;
    data_field='data_MUSC';
else
    data_field='data';
end
    

NumDays=length(AllDays_RawDatStruct);

% 1) BASELINE
% although should be in "epochdata" category below, do first here becuase
% needed for day by day data.

for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    X=[];
    PCmat=[];
    for ii=Params.SeqFilter.BaselineDays;
        if ~isempty(AllDays_RawDatStruct{ii}); % check that day has data.
            if isfield(AllDays_RawDatStruct{ii}.data, syl);
                X=[X cell2mat(AllDays_RawDatStruct{ii}.(data_field).(syl)(:,1))']; % collect all data points from baseline days
                
                PCmat=[PCmat; cell2mat(AllDays_RawDatStruct{ii}.(data_field).(syl)(:,2))]; % collect PCs over all days
            end
        end
    end

    EpochData.Baseline.(syl).pitchcontour_mean=nanmean(PCmat,1);
    EpochData.Baseline.(syl).pitchcontour_STD=nanstd(PCmat,0,1);
    EpochData.Baseline.(syl).rawFF=X;
    EpochData.Baseline.(syl).meanFF=nanmean(X);
    EpochData.Baseline.(syl).stdFF=nanstd(X);
    EpochData.Baseline.(syl).n=length(~isnan(X));
    EpochData.Baseline.(syl).semFF=real(nanstd(X)/sqrt((length(~isnan(X))-1)));
end


% 2) DAY-BY-DAY: Extract stats (per day) for all syls into matrix form
for jj=1:length(SylFields_Unique); % how many fields within this set?
    syl=SylFields_Unique{jj}; % actual syl name (e.g. 'a')
    
    Y={};
    Yvals={};
    Tvals={};
    CatchSongStatus={};
    Ymean=nan(NumDays,1);
    Yerr=nan(NumDays,1);
    Ystd=nan(NumDays,1);
    
    for i=1:NumDays;
        if ~isempty(AllDays_RawDatStruct{i}); % day has data?
            if isfield(AllDays_RawDatStruct{i}.(data_field),syl); % check that day has that specific syl
                
                
                Yvals{i}=AllDays_RawDatStruct{i}.(data_field).(syl)(:,1)'; % raw vals
                Tvals{i}=AllDays_RawDatStruct{i}.(data_field).(syl)(:,6)'; % datenums\
                Ymean(i)=mean(cell2mat(Yvals{i}));
                Yerr(i)=lt_sem(cell2mat(Yvals{i}));
                Ystd(i)=std(cell2mat(Yvals{i}));
                
                % get catch trial and song status, if they exist for this
                % day (previous versions of code did not have)
                if size(AllDays_RawDatStruct{i}.(data_field).(syl),2)>12;
                    % then catch song data was annotated
                    CatchSongStatus{i}=cell2mat(AllDays_RawDatStruct{i}.(data_field).(syl)(:,13));
                end
            end
        end
    end
    
    DataMatrix.(syl).FFvals=Yvals;
    DataMatrix.(syl).Tvals=Tvals;
    DataMatrix.(syl).meanFF=Ymean;
    DataMatrix.(syl).semFF=real(Yerr);
    DataMatrix.(syl).sdFF=Ystd;
    DataMatrix.(syl).CatchSongStatus=CatchSongStatus;
end


% 3) Collect DEVIATION FROM BASELINE and Z-SCORE (individual days)
for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
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
    
    for ii=1:length(SylFields_Unique);
        syl=SylFields_Unique{ii};
        
        DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingHz=...
            DataMatrix.(syl).meanFF_DevFromBase./DataMatrix.(targsyl).meanFF_DevFromBase;
        
        DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingZ=...
            DataMatrix.(syl).meanFF_zscore./DataMatrix.(targsyl).meanFF_zscore;
    end
end

% to average target
for ii=1:length(SylFields_Unique);
    syl=SylFields_Unique{ii};
    
    DataMatrix.(syl).GeneralizationFrom.MeanOverTargs.UsingHz=...
        DataMatrix.(syl).meanFF_DevFromBase./DataMatrix_Targ.DevBaseFF;
    
    DataMatrix.(syl).GeneralizationFrom.(targsyl).UsingZ=...
        DataMatrix.(syl).meanFF_zscore./DataMatrix_Targ.DevBaseZ;
end




% COLLECT RUNNING AVG (multiday bins) over all days, not just WN
RunWind=Params.PlotLearning.DayBinSize;
RunWinInds=NumDays-RunWind+1; % how many windows are there, sliding, size constant (e.g. 5 days, win 3, gives 3 Inds).
DayBinsFieldname=['WindSize_' num2str(Params.PlotLearning.DayBinSize)];
for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    % Sliding window across all days - collect data and get statistics
    for nn=1:RunWinInds; % how many data bins are there (overlapping).
        windowInds=nn:nn-1+RunWind; % e.g. [5 6 7] for 1st bin, if first WN day is 5.
        
        X=[];
        for i=windowInds;
            try % missing data?
                X=[X cell2mat(AllDays_RawDatStruct{i}.(data_field).(syl)(:,1))']; % collect FF values over all days into one matrix;
            catch err
            end
        end
        
        try % sometimes all days in bin doesn't have data.
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).rawFF{nn}=X;
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF(nn)=mean(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).stdFF(nn)=std(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).n(nn)=length(X);
            EpochData.AllDaysSliding.(DayBinsFieldname).(syl).semFF(nn)=real(std(X)/sqrt((length(X)-1)));
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
        
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).FirstWNInd=Params.PlotLearning.WNTimeOnInd;
        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).LastWNInd=Params.PlotLearning.WNTimeOffInd-RunWind+1;
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
if isfield(EpochData, 'AllDaysSliding');
RunWinInds=length(EpochData.AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).meanFF); % num bins
SylCategoriesFields=fieldnames(Params.SeqFilter.SylLists); % how manyways of dividing up syls have I used?

for ll=1:length(SylCategoriesFields);
    fn=SylCategoriesFields{ll}; % e.g. fn = 'SylsSame'
    if ~isempty(Params.SeqFilter.SylLists.(fn));
        if ischar(Params.SeqFilter.SylLists.(fn){1})==1; % if this is char, then this is the end of the line for this structure, so start using these entries as the syl names. Otherwise, go one lower (below) to find syl names.
            FieldsList=Params.SeqFilter.SylLists.(fn);
            for ii=1:length(FieldsList);
                syl=FieldsList{ii};
                
                % 1) learning (hz min baseline)
                try
                EpochData.MatrixOverDaysforSylLists.(fn).meanFF_minusBaseline(:,ii)=...
                    EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline';
                catch err
                end
                
                %             % 2) generalization
                %             for j=1:length(SylLists.TargetSyls);
                %                 targsyl=SylLists.TargetSyls{j};
                %
                %                 EpochData.MatrixOverDaysforSylLists.FieldsInOrder.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                %                     EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                %
                %             end
                
                
                % % FOR SINGLE DAYS
                try
                DataMatrix.MatrixOverDaysforSylLists.(fn).meanFF_minusBaseline(:,ii)=...
                    DataMatrix.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                catch err
                end
                
                
            end
            EpochData.MatrixOverDaysforSylLists.(fn).FieldNamesInOrder=FieldsList;
            DataMatrix.MatrixOverDaysforSylLists.(fn).FieldNamesInOrder=FieldsList;
            
            
        else
            for i=1:length(Params.SeqFilter.SylLists.(fn)); % how many lists? (e.g. fields in order might have 2, but same syls is always 1)
                FieldsList=Params.SeqFilter.SylLists.(fn){i};
                for ii=1:length(FieldsList);
                    syl=FieldsList{ii};
                    
                    % 1) learning (hz min baseli)
                    try
                    EpochData.MatrixOverDaysforSylLists.(fn){i}.meanFF_minusBaseline(:,ii)=...
                        EpochData.AllDaysSliding.(DayBinsFieldname).(syl).meanFF_minusBaseline';
                    catch err
                    end
                    
                    %             % 2) generalization
                    %             for j=1:length(SylLists.TargetSyls);
                    %                 targsyl=SylLists.TargetSyls{j};
                    %
                    %                 EpochData.MatrixOverDaysforSylLists.FieldsInOrder.FieldsInOrderNum{i}.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                    %                     EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                    %
                    %             end
                    
                    % % FOR SINGLE DAYS
                    try
                        DataMatrix.MatrixOverDaysforSylLists.(fn){i}.meanFF_minusBaseline(:,ii)=...
                        DataMatrix.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                    catch err
                    end
                        
                    
                    
                end
                EpochData.MatrixOverDaysforSylLists.(fn){i}.FieldNamesInOrder=FieldsList;
                DataMatrix.MatrixOverDaysforSylLists.(fn){i}.FieldNamesInOrder=FieldsList;
                
            end
        end
    end
end
end
