%% lt 12/5/15 - delete noise files from multiple days
function lt_batchsong_delete_noisefiles_MULTDAYS(day1, day2, quick_mode)
% Run in bird folder
% day1='12Sep2015';
% day2='16Sep2015';
% quick_mode (See single day code)

%% collect dirs within dates
MetadataStruct=lt_metadata_collect;

experiment = '';
condition='';
notes='';
date_range={day1,day2};
only_labeled_dirs=0;

ListOfDirs=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs, 2);

%% go through all days and run

for i=1:length(ListOfDirs)
    close all;
    cd(ListOfDirs(i).dirname)
    
    disp(ListOfDirs(i).dirname);
    
%     quick_mode=1;
    rerun_cleandir=1;
    lt_batchsong_delete_noisefiles(quick_mode, rerun_cleandir)
    
    cd ..
    close all;
end
    
    