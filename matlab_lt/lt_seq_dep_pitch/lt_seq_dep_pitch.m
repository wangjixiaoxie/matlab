%% LT 1/31/15 - rewriting lt_compile_seq_dep..., this time much more modular, allowing for variuos analyses to be separately performed on raw data 


%% 1) Gather raw data for this day
clear all; close all;

Params.DayRawDat.fs=32000;
Params.DayRawDat.pc_harms=1; % harmonics to take weighted avg over. 1 or 2 is good.
Params.DayRawDat.batch='batch.labeled.all';
Params.DayRawDat.syllables={'a','b','c'};
Params.DayRawDat.frequency_range={[1300 2200], [2800 3950],[2150 3150]};
Params.DayRawDat.pc_dur=[0.12,0.11,0.10];
% Params.DayRawDat.pc_time_window={[375 525],[60 220],[55 320]};
Params.DayRawDat.pc_time_window={[375 525],[30 50],[55 320]}; % WN over b
Params.DayRawDat.pc_sigma=1;

% --- trying to get all syllables - look at code to check
% Params.DayRawDat.syllables={'a','ab','cb','bb','cc','bc','dc'};
% Params.DayRawDat.frequency_range={[1300 2200], [2800 3950],[2800 3950],[2800 3950],[2150 3150],[2150 3150],[2150 3150]};
% Params.DayRawDat.pc_dur=[0.12,0.09,0.09,0.09,0.11,0.11,0.11];
% Params.DayRawDat.pc_time_window={[375 525],[60 220],[60 220],[60 220],[55 320],[55 320],[55 320]};

% plot and save?
plotON=1;
saveON=1;

% Related to LMAN inactivation
plotLMANinact=1;
Params.DayRawDat.Musc_On_Time='1153'; % time given muscimol - will plot data with temporal lag after this.
Params.DayRawDat.Musc_Off_Time='1700';

[Params, RawDatStruct]=lt_seq_dep_pitch_DayRawDat(Params,1,1,'',plotLMANinact);


%% Script to change name of all song files in a day to stick "PBS" or "MUSC" right after bird name

StringToAdd='PBS';

FilesInFolder=dir('*'); % get all cbins, cbinnotmat, and rec

% copy all stuff to backup folder
mkdir OldSongFiles
!cp * OldSongFiles;

% continue
for i=1:length(FilesInFolder);
    fn=FilesInFolder(i).name;
    
    if any(strfind(fn,'.cbin')) || any(strfind(fn,'.rec')) || any(strfind(fn,'.not.mat'));
    
        fn_new=[fn(1:9) StringToAdd '_' fn(10:end)];
        
        eval(['!mv ' fn ' ' fn_new]);
        
    end
end


%% TO DO OVER ALL DAYS
clear all; close all
phrase = 'SeqDepPitchShift2';
first_day= '23Nov2014';
last_day= '20Dec2014';
% first_day= '19Nov2014';
% last_day= '22Nov2014';
save_results=1;

% functions to run (SAME FOR ALL MOTIFS)
FcnAll={'seq_dep_pitch_2'};

% Parameters for functions within
Params.DayRawDat.fs=32000;
Params.DayRawDat.pc_harms=1; % harmonics to take weighted avg over. 1 or 2 is good.
Params.DayRawDat.batch='batch.labeled.catch';
Params.DayRawDat.syllables={'a','b','c'};
Params.DayRawDat.frequency_range={[1300 2200], [2800 3950],[2150 3150]};
Params.DayRawDat.pc_dur=[0.12,0.11,0.10];

Params.DayRawDat.pc_time_window={[375 525],[60 220],[55 320]};
Params.DayRawDat.pc_sigma=1;


plotON=0;
saveON=1;

WithinParams={'ParamsSDP',Params,'plotON_SDP',plotON,'saveON_SDP',saveON};

[filename_save all_days_various]=lt_all_days_various_calculations_FUNCTION(phrase,first_day,last_day,FcnAll,WithinParams,save_results);



%% TROUBLESHOOTING - to look at ampl contour of specific syl to check if got enough data
% load RawDatStruct, then do:


figure; hold on;

trials_to_plot=[5 10 15 20];

for i=1:length(trials_to_plot);
    trial=trials_to_plot(i);
plot(smooth(log(RawDatStruct.data.c{trial,3}.^2),40));

end




%% 2) Seq filter, remove outliers, and compile raw data, and enter experiment info into params

clear all; close all;

% 0) keep?
Params.SeqFilter.AmplThr=2200;

% 1) Seq filter and remove outliers and compile into one struct
Params.SeqFilter.all_daysON=1; % If 1, then doesn't matter what I enter for days argumemtns.
Params.SeqFilter.FirstDay='';
Params.SeqFilter.LastDay='';

Params.SeqFilter.SeqPreList={'a','ab','c','cb','bcc','bccb','dcc','dccb','b','bc','d','dc'}; % To skip seq filter, keep blank. (i.e. {})
Params.SeqFilter.SylTargList={'b','b','b','b','b','b','b','b','c','c','c','c'};
Params.SeqFilter.SeqPostList={'','','','','','','','','','','',''};

% repeats?
Params.SeqFilter.repeats={'jBc','jBd'}; % to skip repeat filter, don't even define the field "repeats". Will find all things that are 1) category: j[BBB...]c (j and c are optional), 2) category j[BBBB]

% Regular expressions
Params.SeqFilter.RegExpr.expressions={'acb+g', 'acb+'}; 


% 2) experiment info
Params.SeqFilter.WNTimeON='23Nov2014-0000'; % Time WN turned on (1st WN day)
Params.SeqFilter.WNTimeOFF= '20Dec2014-2400'; % Time WN turned off (last WN day) ( use 0000 and 2400 if only plotting days)
Params.SeqFilter.BaselineDays=1:4;

Params.SeqFilter.SylLists.FieldsInOrder{1}={'a','aB','abB','bC','bcC','bccB','bccbB'};
Params.SeqFilter.SylLists.FieldsInOrder{2}={'dC','dcC','dccB','dccbB'};
Params.SeqFilter.SylLists.TargetSyls={'dccB'};
Params.SeqFilter.SylLists.SylsSame={'aB','abB','bccB','bccbB','dccbB'};
Params.SeqFilter.SylLists.SylsDifferent={'a','bC','bcC','dC','dcC'};

Params.SeqFilter.DaysForSnapshot{1}={'09Dec2014','11Dec2014'};
Params.SeqFilter.DaysToMark= {'11Dec2014-2400'}; % will mark all plots with lines here;


Params.SeqFilter.SylLists.FieldsToPlot{1}=[Params.SeqFilter.SylLists.TargetSyls Params.SeqFilter.SylLists.SylsSame];
Params.SeqFilter.SylLists.FieldsToPlot{2}=Params.SeqFilter.SylLists.SylsDifferent;

% 3) RUN
plotON=1;
[Params, AllDays_RawDatStruct]=lt_seq_dep_pitch_SeqFilterCompile(Params,plotON);



%% (OPTIONAL) - Change time windows for syllables in specific contexts.
% replaces FF value with new value, taken over average over pitch contour

% indices refer to location in fieldnames
Params.RecalculateFF.pc_time_window_list=Params.SeqFilter.pc_time_window_list{1}; % list is the same as specified for seq filter, except one change:
Params.RecalculateFF.pc_time_window_list(:,8)=[30 50]; % bccb should have different time window.
plotON=1;

[Params, AllDays_RawDatStruct] = lt_seq_dep_pitch_RecalculateFF(Params, AllDays_RawDatStruct, plotON);



%% 3) Perform various analyses on that data structure

% Extra params, for muscimol only.  do not specific this field if don't care about muscimol.  
% below, times of muscimol in and out (2 times per day)
Params.PlotLearning.MuscimolSchedule={...
    {'03May2015', '1337', '1745'}, ...
    {'04May2015', '1259', '1700'}, ...
    {'05May2015', '1401', '2036'}, ...
    {'07May2015', '1125', '1719'}, ...
    {'08May2015', '1155', '1636'}, ...
    {'09May2015', '1404', '1915'}, ...
    {'10May2015', '1238', '1831'}, ...
    {'11May2015', '1215', '1555'}, ...
    {'18May2015', '1150', '1643'}};


Params.PlotLearning.MuscimolDaysToThrowout={'08May2015', '10May2015'};

% params
Params.PlotLearning.plotWNdays=1; % if 1, then plots WN lines, if 0, then no.
Params.PlotLearning.DayBinSize=3; % 3 day running avg.
saveON=1;


[Params, AllDays_PlotLearning]=lt_seq_dep_pitch_PlotLearning(Params, AllDays_RawDatStruct,saveON);



%% 4) PLOT - looking at effects of LMAN inactivation
    Params.PlotLearning.Lag_time=1.7; % time from switch to musc
    Params.PlotLearning.PBS_window=[-1.5 0]; % time before switch for PBS
    
Params.PlotLearning.timeline.consolid_start='04Jul2015';
Params.PlotLearning.timeline.consolid_end='08Jul2015';
Params.PlotLearning.timeline.bidir_start='09Jul2015';
Params.PlotLearning.timeline.bidir_end='12Jul2015';


[Params, AllDays_PlotLearning]= lt_seq_dep_pitch_PlotLearning_Musc(Params, AllDays_PlotLearning);



%% 3) Extract structure statistics

% to extract data
[Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct);

% to look at syllable similarity - being replaced by "Correlations" below
% [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats_SylSimilarity(Params, AllDays_StructStatsStruct);

% to plot and perform PCA on syl structure
Params.PCA.epoch='baseline'; % will look at baseline
% Params.PCA.epoch=[21 24]; % days 21 to 24

[Params, PCA_Struct]=lt_seq_dep_pitch_StructureStats_featurePCA(Params, AllDays_StructStatsStruct,1);



%% 4) Look at correlations between syllables
% work in progress
DaysWanted='baseline'; % either baseline (astring) or array
lt_seq_dep_pitch_Correlations(Params, AllDays_StructStatsStruct,DaysWanted);


%% ANOTHER WAY TO DO CORRELATIONS, USING RAW DAT STRUCT - e.g. parsing by motif and so on.





%% 5) Across Birds/Expts
% script:

lt_seq_dep_pitch_ACROSSBIRDS;



%% ========================== Use regular expressions to sort data from Raw data
% Params.RegExpr.expressions={'acb+g', 'acb+'};
Params.RegExpr.expressions={'abbccbb', 'dccbb'};
DoLMAN=1;
saveON=1;
[Params, AllDays_RegExpr] = lt_seq_dep_pitch_RegExpr(Params, AllDays_RawDatStruct, saveON, DoLMAN, AllDays_PlotLearning);


% ========================== PLOT REG EXPR DATA
Params.PlotRegExpr.plotWNdays=1;
saveON=1;
LMANon=0;
[Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon);


