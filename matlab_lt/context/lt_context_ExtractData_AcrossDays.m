function [Dirs_Missing_notegroups]=lt_context_ExtractData_AcrossDays(Params)
%% LT 9/2/15 - Goes into each day and performs extraction of data - saves to birdname/context_data/ExperimentID
% INSTRUCTIONS:
% run in bird folder.

% example:
% % metadata params
% Params.metadata.experiment = 'LightDepPitch';
% Params.metadata.condition='';
% Params.metadata.notes='';
% Params.metadata.date_range={'05May2015','25May2015'};
% Params.metadata.only_labeled_dirs=0;
% 
% % analysis params
% Params.analysis.batch='batch.keep';
% % last version config (template used at very end of SeqDepPitch)
% Params.analysis.config = '/bluejay4/lucas/birds/rd23gr89/050715_LightDepPitch_durWN_day1/config.evconfig2';
% Params.analysis.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.
% Params.analysis.Make_BatchKeep=1; % will make a new batch.keep for each day

%%

currdir=pwd;


%% EXTRACT dir names that will go into

MetadataStruct=lt_metadata_collect;
ListOfDirs=lt_metadata_find_dirs(MetadataStruct, Params.metadata.experiment, Params.metadata.condition, Params.metadata.notes, Params.metadata.date_range, Params.metadata.only_labeled_dirs);


%% GO THROUGH all days and perform analysis
Dirs_Missing_notegroups={};

for i=1:length(ListOfDirs)
    
    cd(ListOfDirs{i});
    
    
    % ===== FIRST, make a new batch.keep?
    if Params.analysis.Make_BatchKeep==1;
        lt_make_batch(1);
    end
    
    
    % ===== SECOND, perform data extraction
    Params_DayExtraction=Params.analysis;
    saveON=1;
    OnlyContinueIfNoteGroups=0;
    [~, error_NoNoteGroup]=lt_context_ExtractDayData(Params_DayExtraction, saveON, OnlyContinueIfNoteGroups);
    
    % ----- COLLECT ALL THE DAYS THAT HAD INCOMPLETE NOTE GROUPS
    if error_NoNoteGroup==1;
        Dirs_Missing_notegroups=[Dirs_Missing_notegroups ListOfDirs{i}];
    end
    
    cd(currdir)
end

%% DONE

disp(' ');
disp('--------------------');
disp('lt_context_ExtractData_AcrossDays DONE!');
disp('Days with incomplete note groups:');
if isempty(Dirs_Missing_notegroups);
    disp('NONE');
else
    disp(Dirs_Missing_notegroups);
end

