%% ========================== 11/18/15 - wh25 context autolabel
clear all; close all;
batch = 'batch.keep';
% config='/bluejay4/lucas/birds/wh25pk77/111515_CtxtDepPitch_away_durWN/CONFIG_away.evconfig2';
config='/bluejay4/lucas/birds/wh25pk77/config_112015.evconfig2';


syl.targ='b';
syl.pre='';
syl.post=''; 
NoteNum=0; 

ampThresh=4000;
min_dur=20;
min_int=4;

overwrite_notmat=1;

% ---- RUN
MetadataStruct=lt_metadata_collect;

experiment = 'CtxtDepPitch';
condition='';
notes='';
date_range={'06Dec2015', '14Dec2015'};
only_labeled_dirs=0;
ListOfDirs2=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs);

for i=1:length(ListOfDirs2);
   
    cd(ListOfDirs2{i});
   
   lt_make_batch(1);
    
[fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);

cd ..
end


