%% or1 autolabel PART 1

clear all; close all;

ListOfDirs = {...
    '/bluejay5/lucas/birds/or1/051917_Association1_WNoff_STIMon', ...
    '/bluejay5/lucas/birds/or1/052017_Association1_WNoff_STIMon', ...
    '/bluejay5/lucas/birds/or1/052117_Association1_WNoff_STIMon'};


%% PART 2

for i=1:length(ListOfDirs)
    
    cd(ListOfDirs{i});
    
%     lt_make_batch(4);
%     lt_sort_batch_by_labeled('batch.rec_FB');
%     batch = 'batch.rec_FB.UNLABELED';

    lt_make_batch(1);
    lt_sort_batch_by_labeled('batch.keep');
    batch = 'batch.keep.UNLABELED';
    
    config='/bluejay5/lucas/birds/or1/config_autolabel.evconfig2';
    
    config = '/bluejay5/lucas/birds/or1/config_autolabel.evconfig2';

    syl.targ='c';
    syl.pre='kkmhdef';
    syl.post=''; 
    NoteNum=0; 

    ampThresh=9000;
    min_dur=20;
    min_int=5;

    overwrite_notmat=1; % will always make backup folder
    
    [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);
    
    
end

if (0)
    %% ======= POST STUFF
    
    lt_v2_db_transfer_calls(2);
    
    batch = 'batch.labeled.all';
    syl.targ='c';
    syl.pre='kkmhdef';
    syl.post=''; 
    
    [fnames, sylnum]=lt_jc_chcklbl(batch, syl.targ, 0.025,0.025,'','','');
    [vlsorfn vlsorind]=jc_vlsorfn(batch, syl.targ,'','');
    
    %%
    evsonganaly
    
    
    %% ============ 3) Replace hand-checekd mislabeld syls
    lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)
    
    
end




%% ========================== 12/6/15 - pu11 context autolabel
clear all; close all;
batch = 'batch.keep';
% config='/bluejay4/lucas/birds/pu11wh87/config_111715.evconfig2'; 
config='/bluejay4/lucas/birds/pu11wh87/config_112015.evconfig2'; 

syl.targ='b';
syl.pre='cc';
syl.post=''; 
NoteNum=0; 

ampThresh=19000;
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



% 
% % ==== additional
% cd /bluejay4/lucas/birds/pu11wh87/111315_CtxtDepPitch_away_preWN
%         lt_make_batch(1);
%         
%         [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);
% 
% 
% cd /bluejay4/lucas/birds/pu11wh87/111415_CtxtDepPitch_away_preWN
%         lt_make_batch(1);
%         
%         [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);
%       
        
%% =====
[fnames, sylnum]=lt_jc_chcklbl(batch, syl.targ, 0.025,0.025,'','','');
[vlsorfn vlsorind]=jc_vlsorfn(batch, syl.targ,'','');

lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)




%% ================================================
clear all; close all
config= '/bluejay4/lucas/birds/pk32/config_081515.evconfig2'; 

syl.targ='c';
syl.pre='d';
syl.post=''; 
NoteNum=0; 

ampThresh=53000;
min_dur=20;
min_int=4;

overwrite_notmat=1;

% ====== all syls
cd /bluejay4/lucas/birds/pk32/091715_Consolidation2_TargB_STIMon_day7
lt_make_batch(4);
batch = 'batch.rec_FB';
lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)

cd /bluejay4/lucas/birds/pk32/091815_Consolidation2_TargB_STIMon_day8
lt_make_batch(4);
batch = 'batch.rec_FB';
lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)

cd /bluejay4/lucas/birds/pk32/091915_Consolidation2_TargB_STIMon_day9
lt_make_batch(4);
batch = 'batch.rec_FB';
lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)


% ==== subset syls
p=0.5;
cd /bluejay4/lucas/birds/pk32/092015_Consolidation2_TargB_STIMoff
lt_make_batch(4);
randsamp('batch.rec_FB',p)
batch = 'batch.rec_FB.rand';
lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)



%% 5/13 - 

%% wh73 autolabel

clear all; close all;

% done:

DirList={
 };

curdir=pwd;

for i=1:length(DirList);

    cd(DirList{i});
% Collect all songs
lt_make_batch(4);

% Autolabel
if any(strfind(DirList{i},'StimON'))
    Params.batch='batch.rec_FB';
else
    randsamp('batch.rec_FB',0.4);
    Params.batch='batch.rec_FB.rand';
end

Params.ampThresh=21000;
Params.min_dur=13;
Params.min_int=1;

Params.syl.pre='';
Params.syl.post='';
Params.syl.targ='b';

Params.overwrite_notmat=1;


% TEMPLATE SETTINGS
Params.TEMPLATE.templatefile = 'autolabel_templ_b1_v2.dat';

Params.TEMPLATE.cntrng(1).MIN=1;
Params.TEMPLATE.cntrng(1).MAX=3;
Params.TEMPLATE.cntrng(1).NOT=0;
Params.TEMPLATE.cntrng(1).MODE=1;
Params.TEMPLATE.cntrng(1).TH=1;
Params.TEMPLATE.cntrng(1).AND=0;
Params.TEMPLATE.cntrng(1).BTMIN=0;

Params.TEMPLATE.cntrng(2).MIN=1;
Params.TEMPLATE.cntrng(2).MAX=3;
Params.TEMPLATE.cntrng(2).NOT=0;
Params.TEMPLATE.cntrng(2).MODE=1;
Params.TEMPLATE.cntrng(2).TH=2.2;
Params.TEMPLATE.cntrng(2).AND=0;
Params.TEMPLATE.cntrng(2).BTMIN=0;

Params.TEMPLATE.cntrng(3).MIN=1;
Params.TEMPLATE.cntrng(3).MAX=3;
Params.TEMPLATE.cntrng(3).NOT=0;
Params.TEMPLATE.cntrng(3).MODE=1;
Params.TEMPLATE.cntrng(3).TH=2.2;
Params.TEMPLATE.cntrng(3).AND=1;
Params.TEMPLATE.cntrng(3).BTMIN=0;

Params.TEMPLATE.cntrng(4).MIN=0;
Params.TEMPLATE.cntrng(4).MAX=10;
Params.TEMPLATE.cntrng(4).NOT=0;
Params.TEMPLATE.cntrng(4).MODE=0;
Params.TEMPLATE.cntrng(4).TH=3;
Params.TEMPLATE.cntrng(4).AND=0;
Params.TEMPLATE.cntrng(4).BTMIN=0;

Params.TEMPLATE.refract=0.2;

% RUN
[fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_function(Params);

end
 



%% 5/11/15 -



%% pu35 
clear all; close all;
cd /bluejay4/lucas/birds/pu35wh17/051015_LightDepPitch_durWN_Probe_day4

% -- 1) what batch and config file to use to run evtafsim.  if needed, make a batch.
Params.batch='batch.keep';
Params.config='/bluejay4/lucas/birds/pu35wh17/050715_LightDepPitch_durWN_day1/config.evconfig2'; % LightDepPitch (WN ON) 5/7 +
Params.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.

% -- Run evtafsim and save information
AllData=lt_context_ExtractDayData(Params);



%% rd23


clear all; close all

cd /bluejay4/lucas/birds/rd23gr89/051015_LightDepPitch_durWN_Probe_day4;

lt_make_batch(1);

Params.batch='batch.keep';
Params.config='/bluejay4/lucas/birds/rd23gr89/050715_LightDepPitch_durWN_day1/config.evconfig2'; % for LightDepPitch
Params.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.

% -- Run evtafsim and save information
AllData=lt_context_ExtractDayData(Params);





%% PREVIOUS TO 5/11 ================================================



%% BEFORE or INCLUDING 4/27/15
%% Stuff to run opto analysis on
% 
% clear all; close all
% BaseDir='/bluejay3/lucas/birds/wh73pk61/';
% 
% ListOfDirs=...
%     {'012115_HVCChR2_StimAndLearn_durStim_WNoff_day2',...
%     '012215_HVCChR2_StimAndLearn_durStim_WNoff_day3',...
%     '012315_HVCChR2_StimAndLearn_durStim_WNoff_day4',...
%     '012415_HVCChR2_StimAndLearn_durStim_WNoff_day5',...
%     '012515_HVCChR2_StimAndLearn_durStim_WNoff_day6',...
%     '012615_HVCChR2_StimAndLearn_durStim_WNoff_day7'};
% 
% ListOfBatchNames={'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.all'};
% 
% ListOfTrialTypes=[0 1 0 0 0 0 1 1 1];



%% AUTOLABEL STUFF
if (0)
    clear all; close all;

DirList={
'/bluejay3/lucas/birds/wh73pk61/012115_HVCChR2_StimAndLearn_durStim_WNoff_day2',...
'/bluejay3/lucas/birds/wh73pk61/012215_HVCChR2_StimAndLearn_durStim_WNoff_day3',...
'/bluejay3/lucas/birds/wh73pk61/012315_HVCChR2_StimAndLearn_durStim_WNoff_day4',...
'/bluejay3/lucas/birds/wh73pk61/012415_HVCChR2_StimAndLearn_durStim_WNoff_day5',...
'/bluejay3/lucas/birds/wh73pk61/012515_HVCChR2_StimAndLearn_durStim_WNoff_day6',...
'/bluejay3/lucas/birds/wh73pk61/012615_HVCChR2_StimAndLearn_durStim_WNoff_day7'};

curdir=pwd;

for i=1:length(DirList);
    cd(curdir);
    cd(DirList{i});
% Collect all songs
% lt_make_batch(1)

% Autolabel
ampThresh=21000;
min_dur=13;
min_int=1;
batch='batch.rec_FB';

templatefile = 'autolabel_templ_b1_v2.dat';
syl.pre='';
syl.post='';
syl.targ='b';

overwrite_notmat=1;

% TEMPLATE SETTINGS
% TEMPLATE SETTINGS
cntrng(1).MIN=1;
cntrng(1).MAX=3;
cntrng(1).NOT=0;
cntrng(1).MODE=1;
cntrng(1).TH=1;
cntrng(1).AND=0;
cntrng(1).BTMIN=0;

cntrng(2).MIN=1;
cntrng(2).MAX=3;
cntrng(2).NOT=0;
cntrng(2).MODE=1;
cntrng(2).TH=2.2;
cntrng(2).AND=0;
cntrng(2).BTMIN=0;

cntrng(3).MIN=1;
cntrng(3).MAX=3;
cntrng(3).NOT=0;
cntrng(3).MODE=1;
cntrng(3).TH=2.2;
cntrng(3).AND=1;
cntrng(3).BTMIN=0;

cntrng(4).MIN=0;
cntrng(4).MAX=10;
cntrng(4).NOT=0;
cntrng(4).MODE=0;
cntrng(4).TH=3;
cntrng(4).AND=0;
cntrng(4).BTMIN=0;


refract=0.2;

% 1) RUN
lt_jc_autoLabel(batch,templatefile,cntrng,syl,refract,ampThresh,min_dur,min_int,overwrite_notmat)

end
end

%%
lt_make_batch(4);

lt_read_ltrec2('note_group_log.ltrec2');

% check note group 2 nothing missed? (since only WN at trigger) - there
% shouldn't be any songs here:
cleandirAuto('batch.NoteGroup_2.rec_noFB',1000,4,4);



% clean note group 3
cleandirAuto('batch.NoteGroup_3',1000,4,4);


%% Checking if discarded songs have any real songs

cleandirAuto('batch.dcrd',1000,4,4)

lt_disp_random_songs_syls('batch.dcrd.keep',0,1,50)