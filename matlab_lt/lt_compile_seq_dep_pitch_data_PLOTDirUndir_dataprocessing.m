%% FIRST, GET BASELINE MEANS
% although should be in "epochdata" category below, do first here becuase
% needed for day by day data.

for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    % 1) Baseline
    % UNDIR
    X=[];
    for ii=BaselineDays;
        try % in case day lacks data
            X=[X cell2mat(AllDays_compiled_DirAndUndir.UNDIR{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
        catch err
        end
    end
    
    EpochData.Baseline.UNDIR.(syl).rawFF=X;
    EpochData.Baseline.UNDIR.(syl).meanFF=mean(X);
    EpochData.Baseline.UNDIR.(syl).stdFF=std(X);
    EpochData.Baseline.UNDIR.(syl).n=length(X);
    EpochData.Baseline.UNDIR.(syl).semFF=std(X)/(length(X)-1);
    
    % DIR ------------------
    X=[];
    for ii=BaselineDays;
        try % in case day lacks data
            X=[X cell2mat(AllDays_compiled_DirAndUndir.DIR{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
        catch err
        end
    end
    
    EpochData.Baseline.DIR.(syl).rawFF=X;
    EpochData.Baseline.DIR.(syl).meanFF=mean(X);
    EpochData.Baseline.DIR.(syl).stdFF=std(X);
    EpochData.Baseline.DIR.(syl).n=length(X);
    EpochData.Baseline.DIR.(syl).semFF=std(X)/(length(X)-1);
end


%% SECOND, extract Data (day by day).

% 1) Data matrix: Extract means and SEMS for all syls into matrix form
for jj=1:length(SylFieldsAll); % how many fields within this set?
    syl=SylFieldsAll{jj}; % actual syl name (e.g. 'a')
    
    
    % compile pitch means and SEM for each day, for this syl.
    % UNDIR first
    Y=nan(NumDays,1);
    Yerr=nan(NumDays,1);
    Ystd=nan(NumDays,1);
    
    for i=1:NumDays;
        try
            Y(i)=AllDays_compiled_DirAndUndir.UNDIR{i}.summary_stats.(syl).meanFF;
            Yerr(i)=AllDays_compiled_DirAndUndir.UNDIR{i}.summary_stats.(syl).semFF;
            Ystd(i)=AllDays_compiled_DirAndUndir.UNDIR{i}.summary_stats.(syl).sdFF;
        catch err
        end
    end
    
    DataMatrix.UNDIR.(syl).meanFF=Y;
    DataMatrix.UNDIR.(syl).semFF=Yerr;
    DataMatrix.UNDIR.(syl).sdFF=Ystd;
    
    % DIR second
    Y=nan(NumDays,1);
    Yerr=nan(NumDays,1);
    Ystd=nan(NumDays,1);
    
    for i=1:NumDays;
        try
            Y(i)=AllDays_compiled_DirAndUndir.DIR{i}.summary_stats.(syl).meanFF;
            Yerr(i)=AllDays_compiled_DirAndUndir.DIR{i}.summary_stats.(syl).semFF;
            Ystd(i)=AllDays_compiled_DirAndUndir.DIR{i}.summary_stats.(syl).sdFF;
        catch err
        end
    end
    
    DataMatrix.DIR.(syl).meanFF=Y;
    DataMatrix.DIR.(syl).semFF=Yerr;
    DataMatrix.DIR.(syl).sdFF=Ystd;
    
end


% 3) Collect DEVIATION FROM BASELINE and Z-SCORE (individual days)
for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    %     UNDIR
    DataMatrix.UNDIR.(syl).meanFF_DevFromBase=DataMatrix.UNDIR.(syl).meanFF-EpochData.Baseline.UNDIR.(syl).meanFF;
    DataMatrix.UNDIR.(syl).meanFF_zscore=DataMatrix.UNDIR.(syl).meanFF_DevFromBase./EpochData.Baseline.UNDIR.(syl).stdFF; % zscore = deviation from baseline mean divided by baseline std.
    DataMatrix.UNDIR.(syl).meanFF_PercentFromBase=DataMatrix.UNDIR.(syl).meanFF_DevFromBase./EpochData.Baseline.UNDIR.(syl).meanFF; % learning as percent of baseline FF.
    
    
    %     DIR
    DataMatrix.DIR.(syl).meanFF_DevFromBase=DataMatrix.DIR.(syl).meanFF-EpochData.Baseline.DIR.(syl).meanFF;
    DataMatrix.DIR.(syl).meanFF_zscore=DataMatrix.DIR.(syl).meanFF_DevFromBase./EpochData.Baseline.DIR.(syl).stdFF;
    DataMatrix.DIR.(syl).meanFF_PercentFromBase=DataMatrix.DIR.(syl).meanFF_DevFromBase./EpochData.Baseline.DIR.(syl).meanFF; % learning as percent of baseline FF.
    
    
end

% 5) Transform learning scores to generalization scores (single days)
% For every day divide every syl learning by target syl learning.
for i=1:length(SylLists.TargetSyls);
    targsyl=SylLists.TargetSyls{i};
    for ii=1:length(SylFieldsAll);
        syl=SylFieldsAll{ii};
        
        % UNDIR
        DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
            DataMatrix.UNDIR.(syl).meanFF_DevFromBase./DataMatrix.UNDIR.(targsyl).meanFF_DevFromBase;
        
        DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
            DataMatrix.UNDIR.(syl).meanFF_zscore./DataMatrix.UNDIR.(targsyl).meanFF_zscore;
        
        
        % DIR - do same thing
        DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
            DataMatrix.DIR.(syl).meanFF_DevFromBase./DataMatrix.DIR.(targsyl).meanFF_DevFromBase;

        DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
            DataMatrix.DIR.(syl).meanFF_zscore./DataMatrix.DIR.(targsyl).meanFF_zscore;
        
    end
end



%% THIRD, Will collect three-day (default) running averages across learning. (POINT-BY POINT AVERAGE).
RunWind=3;
NumWNdays=WNTimeOffInd-WNTimeOnInd;
RunWinInds=NumWNdays-RunWind+1; % how many windows are there, sliding, size constant (e.g. 5 days, win 3, gives 3 Inds).
for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    % 1) Baseline
    % UNDIR
    X=[];
    for ii=BaselineDays;
        try % in case day lacks data
            X=[X cell2mat(AllDays_compiled_DirAndUndir.UNDIR{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
        catch err
        end
    end
    
    EpochData.Baseline.UNDIR.(syl).rawFF=X;
    EpochData.Baseline.UNDIR.(syl).meanFF=mean(X);
    EpochData.Baseline.UNDIR.(syl).stdFF=std(X);
    EpochData.Baseline.UNDIR.(syl).n=length(X);
    EpochData.Baseline.UNDIR.(syl).semFF=std(X)/(length(X)-1);
    
    %     DIR ------------------
    X=[];
    for ii=BaselineDays;
        try % in case day lacks data
            X=[X cell2mat(AllDays_compiled_DirAndUndir.DIR{ii}.data.(syl)(:,1))']; % collect all data points from baseline days
        catch err
        end
    end
    
    EpochData.Baseline.DIR.(syl).rawFF=X;
    EpochData.Baseline.DIR.(syl).meanFF=mean(X);
    EpochData.Baseline.DIR.(syl).stdFF=std(X);
    EpochData.Baseline.DIR.(syl).n=length(X);
    EpochData.Baseline.DIR.(syl).semFF=std(X)/(length(X)-1);
    
    
    % SECOND, Sliding window across WN days - collect data and get statistics
    % UNDIR
    for nn=1:RunWinInds; % how many data bins are there (overlapping).
        windowInds=WNTimeOnInd+nn-1:WNTimeOnInd+nn+1; % e.g. [5 6 7] for 1st bin, if first WN day is 5.
        
        
        X=[];
        for i=windowInds;
            try % missing data?
                X=[X cell2mat(AllDays_compiled_DirAndUndir.UNDIR{i}.data.(syl)(:,1))']; % collect FF values over all days into one matrix;
            catch err
            end
        end
        
        try
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).rawFF=X;
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF=mean(X);
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).stdFF=std(X);
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).n=length(X);
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).semFF=std(X)/(length(X)-1);
            
        catch err
        end
        
        % 3) DEVIATION FROM BASELINE
        EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF_minusBaseline=...
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF-EpochData.Baseline.UNDIR.(syl).meanFF;
        
        %z-score
        EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).Zscore=...
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF_minusBaseline/EpochData.Baseline.UNDIR.(syl).stdFF;
        
        
        % 4) Note down the actual indices of these days
        EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).DaysIncluded=windowInds;
        
        % DIR
        X=[];
        for i=windowInds;
            try % possibly missing data
                X=[X cell2mat(AllDays_compiled_DirAndUndir.DIR{i}.data.(syl)(:,1))']; % collect FF values over all days into one matrix;
            catch err
            end
            
        end
        try
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).rawFF=X;
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF=mean(X);
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).stdFF=std(X);
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).n=length(X);
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).semFF=std(X)/(length(X)-1);
        catch err
        end
        
        
        % 3) DEVIATION from baseline
        EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF_minusBaseline=...
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF-EpochData.Baseline.DIR.(syl).meanFF;
        
        %z-score
        EpochData.WNdaysSlidingWin{nn}.DIR.(syl).Zscore=...
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF_minusBaseline/EpochData.Baseline.DIR.(syl).stdFF;
        
        
        % 4) note down the actual indices of these days
        EpochData.WNdaysSlidingWin{nn}.DIR.(syl).DaysIncluded=windowInds;
        
    end
end


% 2b) Convert learning scores to generalization scores (divide by learning
% at target

for k=1:length(SylLists.TargetSyls);
    targsyl=SylLists.TargetSyls{k}; % target;
    for kk=1:length(SylFieldsAll);
        syl=SylFieldsAll{kk};
        for nn=1:RunWinInds
            
            % UNDIR
            % calculate
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
                EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).meanFF_minusBaseline/EpochData.WNdaysSlidingWin{nn}.UNDIR.(targsyl).meanFF_minusBaseline;
            
            EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
                EpochData.WNdaysSlidingWin{nn}.UNDIR.(syl).Zscore/EpochData.WNdaysSlidingWin{nn}.UNDIR.(targsyl).Zscore;
            
            
            % DIR
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
                EpochData.WNdaysSlidingWin{nn}.DIR.(syl).meanFF_minusBaseline/EpochData.WNdaysSlidingWin{nn}.DIR.(targsyl).meanFF_minusBaseline;
        
            EpochData.WNdaysSlidingWin{nn}.DIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
                EpochData.WNdaysSlidingWin{nn}.DIR.(syl).Zscore/EpochData.WNdaysSlidingWin{nn}.DIR.(targsyl).Zscore;

        end
    end
end


%% FOURTH, COllect data comparing Dir and Undir (e.g. difference) (individual
% days)

for j=1:length(SylFieldsAll);
    syl=SylFieldsAll{j};
    
    % collect data
    DataMatrix.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR=DataMatrix.DIR.(syl).meanFF-DataMatrix.UNDIR.(syl).meanFF;
    
    % BASELINE DIFF
    EpochData.Baseline.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR=EpochData.Baseline.DIR.(syl).meanFF-EpochData.Baseline.UNDIR.(syl).meanFF;
end





%% DON'T USE BELOW- averaging using day as datapoints can strongly affect
% result (e.g. low learning biases gen score to be high for one day.

% 6) Get generalization at beginnign, middle, and end of WN, by binning days
% into those three epochs. (averaging day values, not individaul data
% points) ( also will take days even if no song).
% for i=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{i};
%     for ii=1:length(SylFieldsAll);
%         syl=SylFieldsAll{ii};
%
%         % UNDIR
%         LearningEpochs.UNDIR.(syl).FirstThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz(FirstThreeDays));
%
%         LearningEpochs.UNDIR.(syl).MiddleDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz(MiddleDays));
%
%         LearningEpochs.UNDIR.(syl).LastThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz(LastThreeDays));
%
%
%
%         % DIR
%         LearningEpochs.DIR.(syl).FirstThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz(FirstThreeDays));
%
%         LearningEpochs.DIR.(syl).MiddleDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz(MiddleDays));
%
%         LearningEpochs.DIR.(syl).LastThreeDays.GeneralizationFrom.(targsyl).UsingHz.meanFF=...
%             nanmean(DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz(LastThreeDays));
%     end
% end

