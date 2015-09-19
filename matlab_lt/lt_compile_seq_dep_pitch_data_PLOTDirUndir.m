function lt_compile_seq_dep_pitch_data_PLOTDirUndir(DirFilename,UndirFilename,plotDIR,BaselineDays,plotWNdays,WNTimeON,WNTimeOFF,SylLists,DaysForSnapshot,DaysToMark);
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

% DaysForSnapshot{1}={'19Dec2014','21Dec2014'}; then will average from 19
% to 21st of Dec, and plot syl in order with average learning.  For
% multiple epochs, have multiple indices.



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


% Get days to mark, if exist.
global DaysToMarkInds
DaysToMarkInds={};
X={};
if exist('DaysToMark','var');
    for i=1:length(DaysToMark);
        X{i}=lt_convert_EventTimes_to_RelTimes(FirstDay,{DaysToMark{i}});
        DaysToMarkInds{i}=X{i}.FinalValue;
    end
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
lt_compile_seq_dep_pitch_data_PLOTDirUndir_GenOverDays;


%% PLOT SYLLABLES IN ACTUAL SEQUENCE ORDER - 3 day bins
% PLOT IN 3D bar graph, SO CAN LOOK AT CHANGES OVER DAYS

% First convert that data into matrix form
% Do this for both 3-day binned data and single-day data
lt_compile_seq_dep_pitch_data_PLOTDirUndir_SylInOrder;


% PLOTTING IN HEAT MAP - as above
lt_compile_seq_dep_pitch_data_PLOTDirUndir_PlotHeatMap;


% PLOT average over days, or single days (i.e. snapshot)
lt_compile_seq_dep_pitch_data_PLOTDirUndir_Snapshot



%% Compile data into matrix form, for all syllable lists
% (apart from fields in order and fields of interest above

UndirOrDirList={'UNDIR','DIR'};
lt_compile_seq_dep_pitch_data_PLOTDirUndir_ComplSylLists


%% PLOTS of same and different syllables, and comparing early to late learning.

if isfield(SylLists,'SylsSame');
    
    lt_compile_seq_dep_pitch_data_PLOTDirUndir_SameDiff
end


%% PLOT DAY BY DAY CV of pitch





%% PLOTS OF PITCH CONTOURS

% 1) PC for each day, separate windows and overlayed

% 2) PC for select epochs

% 3) shape of PC for each syl compared


%% ONE PLOT FOR EACH SYLLABLE
% Combines DIR and UNDIR
% Plots: single rendition and means, deviations, and "AFP COMPONENT (UNDIR minus DIR)
% Plots single-rendition data points as well.

if (0)
    lt_compile_seq_dep_pitch_data_PLOTDirUndir_EachSyl;
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
global DaysToMarkInds

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd-0.5 WNTimeOffInd-0.5],ylim,'LineStyle','--','Color','r')
    
        for i=1:length(DaysToMarkInds);
            line([DaysToMarkInds{i}-0.5 DaysToMarkInds{i}-0.5],ylim,'LineStyle','--','Color','k');
        end
end
end
