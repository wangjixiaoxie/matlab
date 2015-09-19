%% 7/3/15 - Notes on all autolabel methods up to today

%% Evtaf amp sim version
% =============== Using modified version of Joanne's, which uses evtaf_amp simulation
% This autolabels, and does post-analysis, looking at all syls and
% replacing any mistakes.

Params.batch='batch.rec_FB.rand';

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


%% === Evtafv4 sim version
% After this, can run stuff same as for evtaf amp version to replace false
% positives.

batch = 'batch.rand.keep';
config= '/bluejay4/lucas/birds/pk32/config.evconfig2';

syl.targ='c';
syl.pre='d';
syl.post='cb'; 
NoteNum=0; 

ampThresh=53000;
min_dur=20;
min_int=4;

overwrite_notmat=1;

lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)



%% using feature vectors

lt_autolabel
