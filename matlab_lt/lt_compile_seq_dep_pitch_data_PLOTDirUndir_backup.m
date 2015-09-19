function lt_compile_seq_dep_pitch_data_PLOTDirUndir(DirFilename,UndirFilename,plotDIR,BaselineDays,plotWNdays,WNTimeON,WNTimeOFF,SylLists);
%% TO ADD:
% 1) plot generalization for each syllable
% 2) get avearge generalziation for DIR vs. UNDIR - I bet that it is much larger for DIR
% 3) Plot bar graph of summary genearlization over an image of the sequence itself (e.g. abb...)
% 4) Plot similarity



%% LT 11/18/14 - TO plot over days, pitch data.
% After compiled data with lt_compile_seq_dep_pitch_data_LOADSAVEFILES, run this, entering the filenames
% of teh data structures, and it plots data over days. It plots DIR and UNDIR, and makes comparison plots
% If only UNDIR is desired, then enter the same filename for ddir and undir and have plotDIR=0;
% Example:
%
% plotDIR=0;
% DirFilename='/bluejay3/lucas/birds/pu53wh88/compile_seq_dep_pitch_data_SeqDepPitchShiftDIR/SEQFILTER/AllDays_Compiled/AllDays_Compiled_27Oct2014_to_11Nov2014.mat';
% UndirFilename='/bluejay3/lucas/birds/pu11wh87/compile_seq_dep_pitch_data_SeqDepPitchShift/SEQFILTER/AllDays_Compiled/AllDays_Compiled_26Oct2014_to_17Nov2014';
%
% BaselineDays=1:4; % INCLUSIVE
%
% PLOT LINES TO DESIGNATE WN DAYS
% plotWNdays=1; % 1 is on, 0 is off.
% WNTimeON='30Oct2014-0000'; % Time WN turned on
% WNTimeOFF= '06Nov2014-2400'; % Time WN turned off ( use 0000 and 2400 if only plotting days)
%
% SylLists is a structure that lists syls categorized in ways useful for
% various plots. specifically:
% (EACH ELEM IN THIS CELL ARRAY IS A SET OF THINGS THAT YOU WOULD LIKE TO
% PLOT TOGETHER. EACH SET WILL GET ITS OWN SET OF PLOTS)
% SylLists.FieldsToPlot{1}={'aB','abB','bccB','bccbB','dccB','dccbB'};
% SylLists.FieldsToPlot{2}={'bC','bcC','dC','dcC'};

% SylLists.FieldsInOrder{1}={'a','aB','abB','bC','bcC','bccB','bccbB'};
% SylLists.FieldsInOrder{2}={'dC','dcC','dccB','dccbB'};
%
% SylLists.TargetSyls={'dccB'};
%
% SylLists.SylsSame={'aB','abB','bccB','bccbB','dccbB'};
% SylLists.SylsDifferent={'a','bC','bcC','dC','dcC'};


%% LOAD DIR and UNDIR

if plotDIR==0;
    % If not DIR data, then quick fix - load the same UNDIR data and put
    % that into both the DIR and UNDIR structures. So will plot everything
    % twice, or overlayed.
    
    load(UndirFilename);
    AllDays_compiled_DirAndUndir.DIR=AllDays_compiled_seqdep_pitch;
    
    load(UndirFilename);
    AllDays_compiled_DirAndUndir.UNDIR=AllDays_compiled_seqdep_pitch;
    
elseif plotDIR==1;
    load(DirFilename);
    AllDays_compiled_DirAndUndir.DIR=AllDays_compiled_seqdep_pitch;
    
    load(UndirFilename);
    AllDays_compiled_DirAndUndir.UNDIR=AllDays_compiled_seqdep_pitch;
    
end
%% EXTRACT parameters

CurrDir=pwd;
NumDays=length(AllDays_compiled_DirAndUndir.DIR); % total days
SylFieldsAll=fieldnames(AllDays_compiled_DirAndUndir.UNDIR{1}.data); % all seq and syls

% single syl or sequence?
for i=1:length(SylFieldsAll);
    X(i)=length(SylFieldsAll{i});
end

SylFieldsSingle=SylFieldsAll(X==1); % only single syls
SylFieldsSeq=SylFieldsAll(X>1); % only sequences


% first and last days
try % first day might not have data. but either DIR or UNDIR has to have. (index 1 in both strucutres are the same real day)
    FirstDay=AllDays_compiled_DirAndUndir.UNDIR{1}.PARAMETERS.date;
catch err
    FirstDay=AllDays_compiled_DirAndUndir.DIR{1}.PARAMETERS.date;
end

% last day is that of undir.
LastDay=AllDays_compiled_DirAndUndir.UNDIR{end}.PARAMETERS.date;


% Convert WNdays to index days
global WNTimeOnInd
global WNTimeOffInd

if plotWNdays==1;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{WNTimeON});
    WNTimeOnInd=X.FinalValue;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{WNTimeOFF});
    WNTimeOffInd=X.FinalValue;
end


% What are the early and late epochs of learning?  divide into 1st 3 days, last 3 days, and middle.
LastThreeDays=WNTimeOffInd-3:WNTimeOffInd-1;
FirstThreeDays=WNTimeOnInd:WNTimeOnInd+2;
MiddleDays=WNTimeOnInd+3:WNTimeOffInd-4;

% what days overlap between DIR and UNDIR?


% For saving
bluejaynum=AllDays_compiled_DirAndUndir.UNDIR{1}.PARAMETERS.bluejaynum;
birdname=AllDays_compiled_DirAndUndir.UNDIR{1}.PARAMETERS.birdname;
phrase=AllDays_compiled_DirAndUndir.UNDIR{1}.PARAMETERS.phrase;

timestampSv=lt_get_timestamp(0);
SaveDir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/compile_seq_dep_pitch_data_' phrase '_PLOTS'];


% Save to params structure
params.NumDays=NumDays;
params.SylFieldsAll=SylFieldsAll;
params.SylFieldsSingle=SylFieldsSingle;
params.SylFieldsSeq=SylFieldsSeq;
params.FirstDay=FirstDay;
params.LastDay=LastDay;
params.WNTimeOnInd=WNTimeOnInd;
params.WNTimeOffInd=WNTimeOffInd;
params.birdname=birdname;
params.phrase=phrase;

params.ARGUMENTS.DirFilename=DirFilename;
params.ARGUMENTS.UndirFilename=UndirFilename;
params.ARGUMENTS.plotDIR=plotDIR;
params.ARGUMENTS.BaselineDays=BaselineDays;
params.ARGUMENTS.plotWNdays=plotWNdays;
params.ARGUMENTS.WNTimeON=WNTimeON;
params.ARGUMENTS.WNTimeOFF=WNTimeOFF;
params.ARGUMENTS.SylLists=SylLists;


%% Check that a given real day is the same index for both structures (dir and undir)
% Allows easy plotting on same plots
% This check: there has to be at least one index where both structures have data. Check
% to make sure that day is the same.
check=0;
for i=1:length(AllDays_compiled_DirAndUndir.UNDIR); % num of total days
    try % in case index does not have data in both structures.
        if AllDays_compiled_DirAndUndir.UNDIR{i}.PARAMETERS.date==...
                AllDays_compiled_DirAndUndir.DIR{i}.PARAMETERS.date;
            check=check+1;
        end
    catch err
        continue
    end
end

if check==0;
    disp('ERROR, First days for Dir and Undir do not match up; Aborting');
end


%% DATA PROCESSING --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------

lt_compile_seq_dep_pitch_data_PLOTDirUndir_dataprocessing


%% PLOTTING --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------

%% FIRST, Plot day means over learning.

% PLOT ABSOLUTE VALUES
lt_compile_seq_dep_pitch_data_PLOTDirUndir_MeanPitchDays;

% PLOT RELATIVE TO BASELINE (diff and Z-score)
lt_compile_seq_dep_pitch_data_PLOTDirUndir_LearningOverDays;

% PLOT GENERALIZATION SCORE
lt_compile_seq_dep_pitch_data_PLOTDirUndir_GeneralizationOverDays;


%% PLOT GENERALIZATION
%% PLOT SYLLABLES IN ACTUAL SEQUENCE ORDER - 3 day bins
% PLOT IN 3D bar graph, SO CAN LOOK AT CHANGES OVER DAYS

% First convert that data into matrix form
% Do this for both 3-day binned data and single-day data

lt_compile_seq_dep_pitch_data_PLOTDirUndir_SylInOrder;


%% PLOTTING IN HEAT MAP - as above

lt_compile_seq_dep_pitch_data_PLOTDirUndir_PlotHeatMap;



%% PLOT scatter of Learning and Generalization for same and different syllables, and comparing early to late learning.

% FIRST COMPILE DATA BASED ON SIMILAR/DIFFERENT SYLS - convert to matrix
% form
% First convert that data into matrix form
% Do this for both 3-day binned data and single-day data

% SAME SYLS
% 1) For 3-day binned
for ii=1:length(SylLists.SylsSame);
    syl=SylLists.SylsSame{ii};
    for iii=1:RunWinInds; % how many day bins are there?
        
        % UNDIR
        % 1) learning (hz min baseli)
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        % zscore
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).Zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        end
        
        
        % DIR
        % 1) learning (hz min baseli)
        EpochData.MatrixOverDaysforSylLists.SylsSame.DIR.meanFF_minusBaseline(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.DIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        EpochData.MatrixOverDaysforSylLists.SylsSame.DIR.Zscore(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.DIR.(syl).Zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            EpochData.MatrixOverDaysforSylLists.SylsSame.DIR.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        end
    end
end


% 2) For single-days
for ii=1:length(SylLists.SylsSame);
    syl=SylLists.SylsSame{ii};
    
    % UNDIR
    % 1) learning (hz min baseli)
    DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(:,ii)=...
        DataMatrix.UNDIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(:,ii)=...
        DataMatrix.UNDIR.(syl).meanFF_zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    
    % 2) generalization
    for j=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{j};
        
        DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
            DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
    end
    
    
    % DIR
    % 1) learning (hz min baseli)
    DataMatrix.MatrixOverDaysforSylLists.SylsSame.DIR.meanFF_minusBaseline(:,ii)=...
        DataMatrix.DIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    DataMatrix.MatrixOverDaysforSylLists.SylsSame.DIR.Zscore(:,ii)=...
        DataMatrix.DIR.(syl).meanFF_zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    % 2) generalization
    for j=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{j};
        
        DataMatrix.MatrixOverDaysforSylLists.SylsSame.DIR.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
            DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
    end
end

EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.FieldNames=SylLists.SylsSame;
EpochData.MatrixOverDaysforSylLists.SylsSame.DIR.FieldNames=SylLists.SylsSame;


% DIFF SYLS
for ii=1:length(SylLists.SylsDifferent);
    syl=SylLists.SylsDifferent{ii};
    for iii=1:RunWinInds; % how many day bins are there?
        
        % UNDIR
        % 1) learning (hz min baseli)
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).Zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        end
        
        
        % DIR
        % 1) learning (hz min baseli)
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.DIR.meanFF_minusBaseline(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.DIR.(syl).meanFF_minusBaseline; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.DIR.Zscore(iii,ii)=...
            EpochData.WNdaysSlidingWin{iii}.DIR.(syl).Zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        
        % 2) generalization
        for j=1:length(SylLists.TargetSyls);
            targsyl=SylLists.TargetSyls{j};
            
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.DIR.GeneralizationFrom.(targsyl).UsingHz(iii,ii)=...
                EpochData.WNdaysSlidingWin{iii}.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
        end
    end
end



% 2) For single-days
for ii=1:length(SylLists.SylsDifferent);
    syl=SylLists.SylsDifferent{ii};
    
    % UNDIR
    % 1) learning (hz min baseli)
    DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(:,ii)=...
        DataMatrix.UNDIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(:,ii)=...
        DataMatrix.UNDIR.(syl).meanFF_zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    
    % 2) generalization
    for j=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{j};
        
        DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
            DataMatrix.UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
    end
    
    
    % DIR
    % 1) learning (hz min baseli)
    DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.DIR.meanFF_minusBaseline(:,ii)=...
        DataMatrix.DIR.(syl).meanFF_DevFromBase; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    
    DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.DIR.Zscore(:,ii)=...
        DataMatrix.DIR.(syl).meanFF_zscore; % format: columns are fields in order, and rows are mean data, each being one 3-day bin.
    
    % 2) generalization
    for j=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{j};
        
        DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.DIR.GeneralizationFrom.(targsyl).UsingHz(:,ii)=...
            DataMatrix.DIR.(syl).GeneralizationFrom.(targsyl).UsingHz;
    end
end

EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.FieldNames=SylLists.SylsDifferent;
EpochData.MatrixOverDaysforSylLists.SylsDifferent.DIR.FieldNames=SylLists.SylsDifferent;


% GET MEANS ACROSS SYLLABLES
LL={'SylsSame','SylsDifferent'};
for i=1:2;
    field=LL{i};
    
    % Binned days
    EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.FF_minusBase.mean=...
        mean(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.meanFF_minusBaseline,2);
    EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.FF_minusBase.SD=...
        std(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.meanFF_minusBaseline,0,2);
    
    
    EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.Zscore.mean=...
        mean(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.Zscore,2);
    EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.Zscore.SD=...
        std(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.Zscore,0,2);
    
    for i=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{i};
        
        EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.mean=...
            mean(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.GeneralizationFrom.(targsyl).UsingHz,2);
        EpochData.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.SD=...
            std(EpochData.MatrixOverDaysforSylLists.(field).UNDIR.GeneralizationFrom.(targsyl).UsingHz,0,2);
    end
    
    
    % Single days
    DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.FF_MinusBase.mean=...
        mean(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.meanFF_minusBaseline,2);
    DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.FF_MinusBase.SD=...
        std(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.meanFF_minusBaseline,0,2);
    
    
    DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.Zscore.mean=...
        mean(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.Zscore,2);
    DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.Zscore.SD=...
        std(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.Zscore,0,2);
    
    
    for i=1:length(SylLists.TargetSyls);
        targsyl=SylLists.TargetSyls{i};
        
        DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.mean=...
            mean(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.GeneralizationFrom.(targsyl).UsingHz,2);
        DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.StatsAcrossSyls.GeneralizationFrom.(targsyl).UsingHz.SD=...
            std(DataMatrix.MatrixOverDaysforSylLists.(field).UNDIR.GeneralizationFrom.(targsyl).UsingHz,0,2);
    end
    
end





% PLOT LEARNING in hz and zscore.
figure; hold on;
% 1) learning as hz
subplot(1,2,1); hold on;

% PLOT SIMILAR SYLS
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii)];
    end
    
    plot(Y,'-ob');
    
    Xlims=xlim;
    text(Xlims(2)+0.1,Y(2),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','b');

% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);




% PLOT DIFF SYLS
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii)];
    end
    
    plot(Y,'-or');
    Xlims=xlim;
    text(Xlims(2)+0.15,Y(2),syl,'Color','r'); % label the datapoint
    
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_minusBase.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','r');

% Format plot
xlim([Xlims(1)-0.2 Xlims(2)+0.4])

ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);
title('FF (minus baseline): Same (RED) and Different (BLUE) syls');
lt_plot_zeroline;



% PLOT SAME, BUT Z-score
subplot(1,2,2); hold on;

% PLOT SIMILAR SYLS
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.Zscore(end,ii)];
    end
    
    plot(Y,'-ob');
    
    Xlims=xlim;
    text(Xlims(2)+0.1,Y(2),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','b');

% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);




% PLOT DIFF SYLS
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;
    if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
        MidBinInd=ceil(RunWinInds/2); % middle 3-day bin index
        X=[1 2 3];
        Xtl={'Early','Middle','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(MidBinInd,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(end,ii)];
        
    else
        X=[1 2];
        Xtl={'Early','Late'};
        Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(1,ii),...
            EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.Zscore(end,ii)];
    end
    
    plot(Y,'-or');
    Xlims=xlim;
    text(Xlims(2)+0.15,Y(2),syl,'Color','r'); % label the datapoint
    
end

% plot mean and SEM
if NumWNdays>=3*RunWind; % i.e. then have at least 3 independent bins
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(MidBinInd),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
else
    Y=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.mean(end)];
    
    Ysd=[EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(1),...
        EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.Zscore.SD(end)];
    
    Ysem=Ysd./sqrt(length(FieldsList)-1);
    
end

errorbar(X-0.1,Y,Ysem,'-ok','MarkerFaceColor','r');

% Format plot
xlim([Xlims(1)-0.2 Xlims(2)+0.4])

ylabel('Pitch (z-score)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);
title('Z-score: Same (RED) and Different (BLUE) syls');
lt_plot_zeroline;


subtitle('Mean values from binned days - first and last days of WN');








% WILL NOT PLOT ABSOLUTE VALUES

% % PLOT AS ABOVE, but absolute values.
% figure; hold on;
% % 1) learning as hz
% subplot(1,2,1); hold on;
% 
% % PLOT SIMILAR SYLS
% FieldsList=SylLists.SylsSame;
% YY=[];
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).meanFF_minusBaseline),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline)];
%     
%     plot(X,Y,'-ob');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
% 
% title('Similar Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % PLOT DIFF SYLS
% FieldsList=SylLists.SylsDifferent;
% YY=[];
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).meanFF_minusBaseline),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline)];
%     
%     plot(X,Y,'-or');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
% 
% 
% title('Different Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % zero line
% line(xlim,[0 0],'Color','k','LineStyle','--')
% 
% % global title
% title('ABSOLUTE LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
% 
% 
% 
% % PLOT SAME, BUT Z-score
% subplot(1,2,2); hold on;
% 
% % PLOT SIMILAR SYLS
% FieldsList=SylLists.SylsSame;
% YY=[];
% 
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).Zscore),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore)];
%     
%     plot(X,Y,'-ob');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
% 
% 
% title('Similar Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% 
% % PLOT DIFF SYLS
% FieldsList=SylLists.SylsDifferent;
% YY=[];
% 
% for ii=1:length(FieldsList); % how many fields within this set?
%     syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%     
%     % Plot LEARNING (dev from baseline)
%     % put early and late into 2 elements of a matrix;
%     X=[1 2];
%     Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).Zscore),...
%         abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore)];
%     
%     plot(X,Y,'-or');
%     text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).Zscore),syl); % label the datapoint
%     xlim([0 3]);
%     
%     % Compile syls data to get mean
%     YY=[YY; Y];
%     
% end
% 
% % PLOT MEAN
% Ymean=mean(YY,1);
% Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
% 
% errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
% 
% 
% title('Different Syllables');
% ylabel('Pitch (minus baseline) (hz)');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Early','Late'});
% 
% line(xlim,[0 0],'Color','k','LineStyle','--')
% 
% % global title
% title('ABSOLUTE LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
% 



% WILL NOT PLOT GENERALIZATION - only meaningful when the target is not
% changing pitch
% PLOT SAME, BUT Generalization
% figure; hold on;
% 
% for i=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{i};
%     subplot(length(SylLists.TargetSyls),1,i); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%             EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%             EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     subtitle(['GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
% end



% PLOT SAME, BUT ABSOLUTE Generalization
% figure; hold on;
% 
% for i=1:length(SylLists.TargetSyls);
%     targsyl=SylLists.TargetSyls{i};
%     subplot(length(SylLists.TargetSyls),1,i); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl)),...
%             abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl))];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.UNDIR.(syl).GeneralizationFrom.(targsyl)),...
%             abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl))];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     subtitle(['Absolute GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
% end


%% EXAACTLY LIKE ABOVE, BUT DIR
% NEED TO DO, below is not latest code

% if plotDIR==1;
%     % PLOT LEARNING in hz and zscore.
%     figure; hold on;
%     % 1) learning as hz
%     subplot(1,2,1); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % zero line
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     % PLOT SAME, BUT Z-score
%     subplot(1,2,2); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore];
%         
%         plot(X,Y,'-ob');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore,...
%             EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore];
%         
%         plot(X,Y,'-or');
%         text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore,syl); % label the datapoint
%         xlim([0 3]);
%         
%     end
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     % PLOT AS ABOVE, but absolute values.
%     figure; hold on;
%     % 1) learning as hz
%     subplot(1,2,1); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).meanFF_minusBaseline),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline)];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % zero line
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('ABSOLUTE LEARNING (hz minus baseline) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     
%     
%     
%     
%     
%     % PLOT SAME, BUT Z-score
%     subplot(1,2,2); hold on;
%     
%     % PLOT SIMILAR SYLS
%     FieldsList=SylLists.SylsSame;
%     YY=[];
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore)];
%         
%         plot(X,Y,'-ob');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%     
%     
%     title('Similar Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     
%     % PLOT DIFF SYLS
%     FieldsList=SylLists.SylsDifferent;
%     YY=[];
%     
%     for ii=1:length(FieldsList); % how many fields within this set?
%         syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%         
%         % Plot LEARNING (dev from baseline)
%         % put early and late into 2 elements of a matrix;
%         X=[1 2];
%         Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).Zscore),...
%             abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore)];
%         
%         plot(X,Y,'-or');
%         text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).Zscore),syl); % label the datapoint
%         xlim([0 3]);
%         
%         % Compile syls data to get mean
%         YY=[YY; Y];
%         
%     end
%     
%     % PLOT MEAN
%     Ymean=mean(YY,1);
%     Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%     
%     errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%     
%     
%     title('Different Syllables');
%     ylabel('Pitch (minus baseline) (hz)');
%     set(gca,'XTick',[1 2]);
%     set(gca,'XTickLabel',{'Early','Late'});
%     
%     line(xlim,[0 0],'Color','k','LineStyle','--')
%     
%     % global title
%     title('ABSOLUTE LEARNING (Z-scored) - categorizing into 1) early vs. late, and 2) similar vs. different syls');
%     
%     
%     
%     
%     
%     % PLOT SAME, BUT Generalization
%     figure; hold on;
%     
%     for i=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{i};
%         subplot(length(SylLists.TargetSyls),1,i); hold on;
%         
%         % PLOT SIMILAR SYLS
%         FieldsList=SylLists.SylsSame;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                 EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)];
%             
%             plot(X,Y,'-ob');
%             text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%         
%         
%         title('Similar Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         % PLOT DIFF SYLS
%         FieldsList=SylLists.SylsDifferent;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                 EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)];
%             
%             plot(X,Y,'-or');
%             text(2.1,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%         
%         title('Different Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--')
%         
%         % global title
%         subtitle(['GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
%     end
%     
%     
%     
%     % PLOT SAME, BUT ABSOLUTE Generalization
%     figure; hold on;
%     
%     for i=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{i};
%         subplot(length(SylLists.TargetSyls),1,i); hold on;
%         
%         % PLOT SIMILAR SYLS
%         FieldsList=SylLists.SylsSame;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl)),...
%                 abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl))];
%             
%             plot(X,Y,'-ob');
%             text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');
%         
%         
%         title('Similar Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         % PLOT DIFF SYLS
%         FieldsList=SylLists.SylsDifferent;
%         YY=[];
%         for ii=1:length(FieldsList); % how many fields within this set?
%             syl=FieldsList{ii}; % actual syl name (e.g. 'a')
%             
%             % Plot LEARNING (dev from baseline)
%             % put early and late into 2 elements of a matrix;
%             X=[1 2];
%             Y=[abs(EpochData.WNdaysSlidingWin{1}.DIR.(syl).GeneralizationFrom.(targsyl)),...
%                 abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl))];
%             
%             plot(X,Y,'-or');
%             text(2.1,abs(EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl)),syl); % label the datapoint
%             xlim([0 3]);
%             
%             % Compile syls data to get mean
%             YY=[YY; Y];
%             
%             
%         end
%         
%         % PLOT MEAN
%         Ymean=mean(YY,1);
%         Ysem=std(YY,0,1)/sqrt(size(YY,1)-1);
%         
%         errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');
%         
%         title('Different Syllables');
%         ylabel('Pitch (minus baseline) (hz)');
%         set(gca,'XTick',[1 2]);
%         set(gca,'XTickLabel',{'Early','Late'});
%         
%         
%         
%         line(xlim,[0 0],'Color','k','LineStyle','--')
%         
%         % global title
%         subtitle(['Absolute GENERALIZATION to ' targsyl '(based on mean FF minus baseline)']);
%     end
%     
%     
% end
% 

%% SAME AS ABOVE, BUT ALL DAYS, NOT BINNED DAYS


% PLOT LEARNING in hz and zscore.
figure; hold on;


% PLOT SIMILAR SYLS
subplot(1,2,1); hold on;
FieldsList=SylLists.SylsSame;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;

    
    plot(DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(:,ii),'-ob')
    text(NumDays+0.4,DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii),syl,'Color','b'); % label the datapoint
end

% plot mean and SEM
Y=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.mean;
Ysd=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.SD;
Ysem=DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.StatsAcrossSyls.FF_MinusBase.SD/...
    sqrt(length(FieldsList)-1);

errorbar(Y,Ysem,'-ok','MarkerFaceColor','b');

% Plot the target syl in diff color
targsyl=SylLists.TargetSyls{1};
    plot(DataMatrix.UNDIR.(targsyl).meanFF_DevFromBase,'-og')
    text(NumDays+0.4,DataMatrix.UNDIR.(targsyl).meanFF_DevFromBase(end),targsyl,'Color','k'); % label the datapoint



% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
xlabel('days');

Ylimits=ylim;
Fn_AnnotateWNLines(plotWNdays,ylim);
lt_plot_zeroline;



% PLOT DIFF SYLS
subplot(1,2,2); hold on;
FieldsList=SylLists.SylsDifferent;

for ii=1:length(FieldsList); % how many fields within this set?
    syl=FieldsList{ii}; % actual syl name (e.g. 'a')
    
    % Plot LEARNING (dev from baseline)
    % put early late into 2 elements of a matrix;

    
    plot(DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(:,ii),'-or')
    text(NumDays+0.4,DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii),syl,'Color','r'); % label the datapoint
end

% plot mean and SEM
Y=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.mean;
Ysd=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.SD;
Ysem=DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.StatsAcrossSyls.FF_MinusBase.SD/...
    sqrt(length(FieldsList)-1);

errorbar(Y,Ysem,'-ok','MarkerFaceColor','r');
ylim(Ylimits)

% Format plot
title('Different Syllables');
ylabel('Pitch (minus baseline) (hz)');
xlabel('days');

Fn_AnnotateWNLines(plotWNdays,ylim);
lt_plot_zeroline;


%% PLOT MEAN OF LAST FEW DAYS IN ORDER OF SYLLABLES

    % PLOT LEARNING AS MINUS BASELINE (HZ);
        % UNDIR - Plot figure for this set of fields

        figure; hold on;
        
        for ll=1:length(SylLists.FieldsInOrder);
            subplot(length(SylLists.FieldsInOrder),1,ll); hold on;
            % Get all syls
            FieldsList=SylLists.FieldsInOrder{ll};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            
            Y=[];
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                Y=[Y EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).meanFF_minusBaseline];
                
            end
            
            Ymax=max(Y);
            Ymin=min(Y);
            
            % Plot
            bar(Y);
            
            
            
            % annotate
            Xtick=1:length(FieldsList); % one tick for each syl. needed.
            set(gca,'XTick',Xtick);
            set(gca,'XTickLabel',FieldsList)
            ylabel('FF (hz), minus baseline')
            
        end
        subtitle(['UNDIR: Learning (Hz minus baseline) last 3 WN days']);
        
        
        
        
        
        % DIR - Plot figure for this set of fields
        if plotDIR==1;
        figure; hold on;
        for ll=1:length(SylLists.FieldsInOrder);
            subplot(length(SylLists.FieldsInOrder),1,ll); hold on;
            % Get all syls
            FieldsList=SylLists.FieldsInOrder{ll};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            
            Y=[];
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                Y=[Y EpochData.WNdaysSlidingWin{end}.DIR.(syl).meanFF_minusBaseline];
                
            end
            
            Ymax=max(Y);
            Ymin=min(Y);
            
            % Plot
            bar(Y);
            
            
            
            % annotate
            Xtick=1:length(FieldsList); % one tick for each syl. needed.
            set(gca,'XTick',Xtick);
            set(gca,'XTickLabel',FieldsList)
            ylabel('FF (hz), minus baseline')
            
        end
        subtitle(['DIR: Learning (Hz minus baseline) last 3 WN days']);
        end
    
        
        
    % PLOT LEARNING AS GENERALIZATION COEFFICIENT - IGNORE, just scales
    % learning.  not useful.
%     for kk=1:length(SylLists.TargetSyls);
%         targsyl=SylLists.TargetSyls{kk};
%         
%         for ll=1:length(SylLists.FieldsInOrder);
%             figure; hold on;
%             
%             % Get all syls
%             FieldsList=SylLists.FieldsInOrder{ll};
%             plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
%             
%             
%             % UNDIR - Plot figure for this set of fields
%             subplot(2,1,1); hold on;
%             for jj=1:length(FieldsList); % how many fields within this set?
%                 syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%                 
%                 % Plot
%                 plot(jj,EpochData.WNdaysSlidingWin{end}.UNDIR.(syl).GeneralizationFrom.(targsyl),...
%                     'o-','Color','k','MarkerFaceColor','k','MarkerSize',10);
%                 
%                 
%             end
%             line(xlim,[0 0],'LineStyle','--')
%             
%             % annotate
%             Xtick=1:length(FieldsList); % one tick for each syl. needed.
%             set(gca,'XTick',Xtick);
%             set(gca,'XTickLabel',FieldsList)
%             title(['UNDIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%             ylabel('Generalization (learning (hz) divided by target syl learning)')
%             ylim([-0.5 0.5]);
%             
%             
%             % DIR - Plot figure for this set of fields
%             subplot(2,1,2); hold on;
%             
%             for jj=1:length(FieldsList); % how many fields within this set?
%                 syl=FieldsList{jj}; % actual syl name (e.g. 'a')
%                 
%                 
%                 plot(jj,EpochData.WNdaysSlidingWin{end}.DIR.(syl).GeneralizationFrom.(targsyl),...
%                     's-','Color','k','MarkerFaceColor','k','MarkerSize',12);
%             end
%             
%             line(xlim,[0 0],'LineStyle','--')
%             
%             %         annotate
%             Xtick=1:length(FieldsList); % one tick for each syl. needed.
%             set(gca,'XTick',Xtick);
%             set(gca,'XTickLabel',FieldsList)
%             title(['DIR: Generalization from ' targsyl ', average (day = datapoint) over last 3 WN days']);
%             ylabel('Generalization (learning (hz) divided by target syl learning)')
%             ylim([-0.5 0.5]);
%         end
%     end


%% ONE PLOT FOR EACH SYLLABLE
% Combines DIR and UNDIR
% Plots: single rendition and means, deviations, and "AFP COMPONENT (UNDIR minus DIR)
% Plots single-rendition data points as well.

% PLOT DEVIATION OF MEAN FROM BASELINE
plot_colors=lt_make_plot_colors(length(SylFieldsAll),0); % initiate colors for plot
for j=1:length(SylFieldsAll); % how many fields (i.e. syls)?
    syl=SylFieldsAll{j};
    
    figure; hold on;
    
    % 1) PLOT individual renditions overlayed with mean PITCH (absolute)
    subplot(3,1,1); hold on;
    
    for i=1:NumDays;
        
        % UNDIR
        try % in case no data.
            % get times
            time=cell2mat(AllDays_compiled_DirAndUndir.UNDIR{i}.data.(syl)(:,6));
            [~, B]=lt_convert_datenum_to_hour(time);
            time=B.days+-0.5+i; % so that noon of 1st day is on x=1.
            % get FF
            FF=cell2mat(AllDays_compiled_DirAndUndir.UNDIR{i}.data.(syl)(:,1));
            % Plot
            plot(time,FF,'.','Color','k');
        catch err
        end
        
        % DIR
        try % in case no data.
            % get times
            time=cell2mat(AllDays_compiled_DirAndUndir.DIR{i}.data.(syl)(:,6));
            [~, B]=lt_convert_datenum_to_hour(time);
            time=B.days+-0.5+i; % so that noon of 1st day is on x=1.
            % get FF
            FF=cell2mat(AllDays_compiled_DirAndUndir.DIR{i}.data.(syl)(:,1));
            % Plot
            plot(time,FF,'.','Color','b');
        catch err
        end
        
        
    end
    % Overlay with mean values.
    h=errorbar(DataMatrix.UNDIR.(syl).meanFF,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',8);
    errorbar_tick(h,200);
    
    h=errorbar(DataMatrix.DIR.(syl).meanFF,DataMatrix.UNDIR.(syl).semFF,'s--','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',8);
    errorbar_tick(h,200);
    
    
    %annotate
    legend(SylFieldsAll{j});
    title(['UNDIR (circle) and DIR (square): Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz)')
    xlabel('days')
    xlim([1 NumDays]);
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    % 2) PLOT MEAN PITCH (rel to baseline)
    subplot(3,1,2); hold on;
    
    h=errorbar(DataMatrix.UNDIR.(syl).meanFF_DevFromBase ,DataMatrix.UNDIR.(syl).semFF,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',12);
    errorbar_tick(h,200);
    
    h=errorbar(DataMatrix.DIR.(syl).meanFF_DevFromBase ,DataMatrix.UNDIR.(syl).semFF,'s--','Color',plot_colors{j},'MarkerSize',12);
    errorbar_tick(h,200);
    
    
    %annotate
    %     legend(SylFieldsAll{j});
    title(['UNDIR (circle) and DIR (square): Mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)')
    xlabel('days')
    xlim([1 NumDays]);
    
    % Line for 0
    line(xlim, [0 0],'Color','k','LineStyle','--')
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    % 3) % PLOT DIFFERENCE BETWEEN DIR AND UNDIR
    subplot(3,1,3); hold on;
    DirMinusUndir=DataMatrix.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR;
    BaselineDirMinusUndir=EpochData.Baseline.DIRvsUNDIR.(syl).meanFF_DIRminusUNDIR;
    
    PutAFPComp=-(DirMinusUndir-BaselineDirMinusUndir); % i.e. UNDIR minus DIR, taking into account DIR positive pitch effect at baseline.
    
    plot(PutAFPComp,'o-','Color','k','MarkerFaceColor',plot_colors{j},'MarkerSize',12)
    xlim([1 NumDays]);
    
    % Line for 0
    line(xlim, [0 0],'Color','k','LineStyle','--')
    
    %annotate
    %     legend(SylFieldsAll{j});
    title(['PUTATIVE AFP COMPONENT (UNDIR minus DIR, baseline subtracted): ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('hz')
    xlabel('days')
    
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end





%% SAVE
try
    cd(SaveDir);
catch err
    mkdir(SaveDir);
    cd(SaveDir);
end

savedir2=[FirstDay 'to' LastDay '_made' timestampSv]; % dir for this specific run of code.
mkdir(savedir2);
cd(savedir2);

savemultfigs;

% TO ADD - Save data structure
FinalStruct.AllDays_compiled_DirAndUndir=AllDays_compiled_DirAndUndir;
FinalStruct.DataMatrix=DataMatrix;
FinalStruct.EpochData=EpochData;
FinalStruct.params=params;

save('SeqDepPitch', 'FinalStruct','-v7.3')


end

%% VARIOUS SUNFUNCTIONS

function Fn_AnnotateWNLines(plotWNdays,ylim)

global WNTimeOnInd
global WNTimeOffInd

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd-0.5 WNTimeOffInd-0.5],ylim,'LineStyle','--')
end
end
