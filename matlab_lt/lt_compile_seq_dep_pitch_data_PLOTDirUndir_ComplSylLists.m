% FIRST COMPILE DATA BASED ON SIMILAR/DIFFERENT SYLS - convert to matrix
% form
% First convert that data into matrix form
% Do this for both 3-day binned data and single-day data

AllSylListFields=fieldnames(SylLists);

% Remove fields already analyzed previously or not wanted
IndsToRemove=[];
for i=1:length(AllSylListFields);
    if strcmp(AllSylListFields{i},'FieldsToPlot')==1;
        IndsToRemove=[IndsToRemove i];
    elseif strcmp(AllSylListFields{i},'FieldsInOrder')==1;
        IndsToRemove=[IndsToRemove i];
    elseif strcmp(AllSylListFields{i},'TargetSyls')==1;
        IndsToRemove=[IndsToRemove i];
    end
end

AllSylListFields(IndsToRemove)=[];


% COMPILE
for jj=1:length(AllSylListFields);
    for jjj=1:2;
        UndirOrDir=UndirOrDirList{jjj};
        
        SylListsField=AllSylListFields{jj};
        
        % 1) For 3-day binned
        for ii=1:length(SylLists.(SylListsField));
            syl=SylLists.(SylListsField){ii};
            for iii=1:RunWinInds; % how many day bins are there?
                
                % (UndirOrDir)
                % 1) learning (hz min baseli)
                EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline(iii,ii)=...
                    EpochData.WNdaysSlidingWin{iii}.(UndirOrDir).(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                
                % zscore
                EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore(iii,ii)=...
                    EpochData.WNdaysSlidingWin{iii}.(UndirOrDir).(syl).Zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                
                % 2) generalization
                for j=1:length(SylLists.TargetSyls);
                    targsyl=SylLists.TargetSyls{j};
                    
                    EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                        EpochData.WNdaysSlidingWin{iii}.(UndirOrDir).(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
                end
            end
        end
        
        
        % 2) For single-days
        for ii=1:length(SylLists.(SylListsField));
            syl=SylLists.(SylListsField){ii};
            
            % (UndirOrDir)
            % 1) learning (hz min baseli)
            DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline(:,ii)=...
                DataMatrix.(UndirOrDir).(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
            
            DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore(:,ii)=...
                DataMatrix.(UndirOrDir).(syl).meanFF_zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
            
            
            % 2) generalization
            for j=1:length(SylLists.TargetSyls);
                targsyl=SylLists.TargetSyls{j};
                
                DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
                    DataMatrix.(UndirOrDir).(syl).GeneralizationFrom.(targsyl).UsingHz;
            end
        end
        
        EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).FieldNames=SylLists.(SylListsField);
        EpochData.MatrixOverDaysforSylLists.(SylListsField).DIR.FieldNames=SylLists.(SylListsField);
        
        
        % GET MEANS ACROSS SYLLABLES
        % Binned days
        EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.FF_minusBase.mean=...
            mean(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline,2);
        EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.FF_minusBase.SD=...
            std(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline,0,2);
        
        
        EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.Zscore.mean=...
            mean(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore,2);
        EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.Zscore.SD=...
            std(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore,0,2);
        
        for i=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{i};
            
            EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.mean=...
                mean(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz,2);
            EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.SD=...
                std(EpochData.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz,0,2);
        end
        
        
        % Single days
        DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.FF_MinusBase.mean=...
            mean(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline,2);
        DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.FF_MinusBase.SD=...
            std(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).meanFF_minusBaseline,0,2);
        
        
        DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.Zscore.mean=...
            mean(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore,2);
        DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.Zscore.SD=...
            std(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).Zscore,0,2);
        
        
        for i=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{i};
            
            DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.mean=...
                mean(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz,2);
            DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.SD=...
                std(DataMatrix.MatrixOverDaysforSylLists.(SylListsField).(UndirOrDir).GeneralizationFrom.(targsyl).UsingHz,0,2);
        end
    end
end
