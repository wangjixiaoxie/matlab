%% LT 4/13/15 - Compiles data and plots across experiments/birds.  Need to specific dir of saved structures
%% TO DO:
% 1) remove repeat experiments - i.e. average over them.


%% INPUT PARAMS
clear all; close all;

PARAMS.global.use_zscore_learning=1; % not important. keep as 1.
PARAMS.global.ConsolBinDaySize=2;  % days at start and end of consolid
PARAMS.global.learning_metric='zscore';
% options: 'zscore' (zscore in epoch first pass threshold, defined
% below)
%               'ff_consolstart' (mean ff, at consol start, defined in cell)
PARAMS.global.remove_bad_syls=1; % e.g. those hard to quantify
PARAMS.global.SylsToRemove=...
    {'gr41gr90','SeqDepPitchShift',{'d', 'Jbba','g'}, ...
    'gr41gr90','SeqDepPitchShift2',{'d', 'Jbba','g'}, ...
    'gr41gr90','SeqDepPitchLMAN',{'Jbba','d', 'g'}, ...
    'gr41gr90','SeqDepPitchLMAN2',{'d','Dbba','dBba','dbBa','dbbA', 'dbbaC','dbbacB','dbbacbB', 'Jbba', 'g'}, ... 
    'rd12pu6','SeqDepPitch',{'Jjb','jJb','Ja','h'}, ...
    'rd12pu6','SeqDepPitch2',{'Jjb','jJb','Ja','h'}, ...
    'rd12pu6','SeqDepPitch3',{'Jjb','jJb','Ja','h'}, ...
    'pu35wh17','SeqDepPitch',{'Jb'}, ...
    'pu35wh17','SeqDepPitch2',{'Jb'}, ...
    'pu37wh20', 'SeqDepPitchShift2', {'Qc', 'qC', 'qcC', 'qccB', 'qccbB', 'qccbbG'}, ...
    'rd23gr89','SeqDepPitch',{'k'},...
    'rd23gr89','SeqDepPitchLMAN',{'k','g'}, ...
    'rd23gr89','SeqDepPitchLMAN2',{'Nah', 'Nak', 'Nd', 'h', 'hJ', 'kJ'}, ...
    'pu64bk13','SeqDepPitchShift', {'d'}, ...
    'pu11wh87','SeqDepPitchLMAN',{'bccbB'}, ...
    'pu11wh87','SeqDepPitchLMAN2', {'bC', 'bcC'}, ...
    'rd28pu64', 'SeqDepPitch', {'Jjb', 'jJb', 'Ja', 'kJ', 'jjbB','dG', 'kjbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'rd28pu64', 'SeqDepPitchLMAN', {'Jjb', 'jJb', 'Ja', 'kJ', 'jjbB','dG', 'kjbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'rd28pu64', 'SeqDepPitchLMAN2', {'Jjb', 'jJb', 'Ja', 'kJ', 'jjbB','dG', 'kjbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bu11or68', 'RepDepPitchShift', {'g'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch2', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch3', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitchLMAN', {'acbbbbB', 'acbbbbbB', 'acbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitchLMAN2', {'g'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or73pu40', 'SeqDepPitch', {'acbbbbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or73pu40', 'SeqDepPitch2', {'acbbbbbbbbbB', 'acbbbbbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'gr87bu18', 'SeqDepPitch2', {'k', 'hbbbbbB', 'abbbbbB'}, ... % NOTE: k contour should be altered, potentailyl could be fine.
    'wh4wh77', 'SeqDepPitchLMAN', {'dbbG', 'cbbG', 'dbbgH', 'cbbgH'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'wh25pk77', 'SeqDepPitchLMAN', {'Jk', 'jK','aJ', 'n'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bk34bk68', 'SeqDepPitchLMAN', {'n', 'nJ', 'njJ', 'njjJ', 'jjbbG', 'k','l','lJ','ljbbG'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bk34bk68', 'SeqDepPitchLMAN3', {'n', 'nJ', 'njJ', 'njjJ', 'jjbbG', 'k','l','lJ','ljbbG'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
%     'pu53wh88','SeqDepPitchLMAN', {'abbB'}, ...



% GR41...
% PARAMS.global.SylsToRemove=...
%     {'gr41gr90','SeqDepPitchShift',{'d'}, ...
%     'gr41gr90','SeqDepPitchShift2',{'d'}, ...
%     'gr41gr90','SeqDepPitchLMAN',{'d'}, ...
%     'gr41gr90','SeqDepPitchLMAN2',{'d','dBba','dbBa'}, ... 
%     'rd12pu6','SeqDepPitch',{'h'}, ...
%     'rd12pu6','SeqDepPitch2',{'h'}, ...
%     'rd23gr89','SeqDepPitch',{'n', 'k'},...
%     'rd23gr89','SeqDepPitchLMAN',{'n','k','g'}, ...
%     'pu64bk13','SeqDepPitchShift', {'d'}, ...
%     'pu11wh87','SeqDepPitchLMAN',{'bccbB'}, ...3263.5
%     'pu11wh87','SeqDepPitchLMAN2', {'bC', 'bcC'}, ...
%     'pu53wh88','SeqDepPitchLMAN', {'abbB'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});


PARAMS.global.timestamp_analy=lt_get_timestamp(0);

%% ============= SAME SOBER - ONLY KEEP UP TO 4 RENDITIONS FOR REPEATS

PARAMS.global.SylsToRemove=...
    {'gr41gr90','SeqDepPitchShift',{'d', 'Jbba','g'}, ...
    'gr41gr90','SeqDepPitchShift2',{'d', 'Jbba','g'}, ...
    'gr41gr90','SeqDepPitchLMAN',{'Jbba','d', 'g'}, ...
    'gr41gr90','SeqDepPitchLMAN2',{'d','Dbba','dBba','dbBa','dbbA', 'dbbaC','dbbacB','dbbacbB', 'Jbba', 'g'}, ... 
    'rd12pu6','SeqDepPitch',{'Jjb','jJb','Ja','h'}, ...
    'rd12pu6','SeqDepPitch2',{'Jjb','jJb','Ja','h'}, ...
    'rd12pu6','SeqDepPitch3',{'Jjb','jJb','Ja','h'}, ...
    'pu35wh17','SeqDepPitch',{'Jb'}, ...
    'pu35wh17','SeqDepPitch2',{'Jb'}, ...
    'pu37wh20', 'SeqDepPitchShift2', {'Qc', 'qC', 'qcC', 'qccB', 'qccbB', 'qccbbG'}, ...
    'rd23gr89','SeqDepPitch',{'k'},...
    'rd23gr89','SeqDepPitchLMAN',{'k','g'}, ...
    'rd23gr89','SeqDepPitchLMAN2',{'Nah', 'Nak', 'Nd', 'h', 'hJ', 'kJ'}, ...
    'pu64bk13','SeqDepPitchShift', {'d'}, ...
    'pu11wh87','SeqDepPitchLMAN',{'bccbB'}, ...
    'pu11wh87','SeqDepPitchLMAN2', {'bC', 'bcC'}, ...
    'rd28pu64', 'SeqDepPitch', {'Jjb', 'jJb', 'Ja', 'kJ', 'dG'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'rd28pu64', 'SeqDepPitchLMAN', {'Jjb', 'jJb', 'Ja', 'kJ', 'dG'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'rd28pu64', 'SeqDepPitchLMAN2', {'Jjb', 'jJb', 'Ja', 'kJ','dG'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bu11or68', 'RepDepPitchShift', {'g', 'jbbbbB', 'jbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch2', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitch3', {'acbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitchLMAN', {'acbbbbB', 'acbbbbbB', 'acbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or100pu10', 'SeqDepPitchLMAN2', {'g'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or73pu40', 'SeqDepPitch', {'acbbbbB', 'acbbbbbB', 'acbbbbbbB', 'acbbbbbbbB', 'acbbbbbbbbB', 'acbbbbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or73pu40', 'SeqDepPitch2', {'acbbbbB', 'acbbbbbB', 'acbbbbbbB', 'acbbbbbbbB', 'acbbbbbbbbB', 'acbbbbbbbbbB', 'acbbbbbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'or73pu40', 'SeqDepPitch3', {'acbbbbB', 'acbbbbbB', 'acbbbbbbB', 'acbbbbbbbB', 'acbbbbbbbbB', 'acbbbbbbbbbB'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'gr87bu18', 'SeqDepPitch2', {'k', 'h', 'hbbbbB','hbbbbbB', 'abbbbB', 'abbbbbB'}, ... % NOTE: k contour should be altered, potentailyl could be fine.
    'gr87bu18', 'RepDepPitch', {'k', 'h', 'hbbbbB','hbbbbbB', 'abbbbbB','abbbbB'}, ... % NOTE: k contour should be altered, potentailyl could be fine.
    'wh4wh77', 'SeqDepPitchLMAN', {'dbbG', 'cbbG', 'dbbgH', 'cbbgH'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'wh25pk77', 'SeqDepPitchLMAN', {'Jk', 'jK','aJ', 'n'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bk34bk68', 'SeqDepPitchLMAN', {'n', 'nJ', 'njJ', 'njjJ', 'jjbbG', 'k','l','lJ','ljbbG'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'bk34bk68', 'SeqDepPitchLMAN3', {'n', 'nJ', 'njJ', 'njjJ', 'jjbbG', 'k','l','lJ','ljbbG'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
%     'pu53wh88','SeqDepPitchLMAN', {'abbB'}, ...



%% consolidation - these are days relative to WN day 1. (i.e. 1 = WN day 1)

% PARAMS.TimeCourse.ConsolidationList={...
%     'pu53wh88', 'SeqDepPitchShift', [8 13], ...
%     'pu53wh88', 'SeqDepPitchShift2', [10 16], ...
%     'pu53wh88', 'SeqDepPitchShift3', [8 13], ...
%     'pu53wh88', 'SeqDepPitchLMAN', [6 10], ... %
%     'pu11wh87', 'SyntaxDepPitchShift_abUP', [5 8], ... %
%     'pu11wh87', 'SeqDepPitchShift2', [11 19], ...
%     'pu11wh87', 'SeqDepPitchShift3', [12 20], ... %
%     'pu11wh87', 'SyntaxDepPitchShift_cbDOWN', [5 11], ...
%     'pu11wh87', 'SeqDepPitchLMAN', [15 28], ...
%     'pu11wh87', 'SeqDepPitchLMAN2', [8 16], ...
%     'pu11wh87', 'SeqDepPitchLMAN6', [4 5], ...
%     'pu37wh20', 'SeqDepPitchShift', [4 6], ...
%     'pu37wh20', 'SeqDepPitchShift2', [4 5], ...
%     'pu64bk13', 'SeqDepPitchShift', [4 8], ...
%     'pu64bk13', 'SeqDepPitchShift2', [5 12], ...
%     'gr41gr90', 'SeqDepPitchShift', [3 12], ...
%     'gr41gr90', 'SeqDepPitchShift2', [10 17], ...
%     'gr41gr90', 'SeqDepPitchLMAN', [4 8], ...
%     'gr41gr90', 'SeqDepPitchLMAN2', [3 6], ...
%     'rd23gr89', 'SeqDepPitch', [7 12], ...
%     'rd23gr89', 'SeqDepPitchLMAN', [4 5], ...
%     'pu35wh17', 'SeqDepPitch', [7 11], ...
%     'rd12pu6', 'SeqDepPitch', [5 9], ...
%     'rd12pu6', 'SeqDepPitch2', [3 8], ...
%     'rd12pu6', 'SeqDepPitch3', [4 12], ...
%     'rd28pu64', 'SeqDepPitch', [5 11], ...
%     'rd28pu64', 'SeqDepPitchLMAN', [6 6], ...
%     'bk34bk68', 'SeqDepPitchLMAN', [5 7]};


%% INPUT DATA STRUCTURES
% from PlotLearning: AllDays_PlotLearning and Params

% Cell holding experiments {birdname, experient ID, NumTargets(see below), directory of plotlearning, ...
% ConsolidationStart, ConsolidationEnd, RegularExpressionsCell, {One day before start bidir, day of last bidir, bidir syl 1, bidir syl 2, multidir syl 3...}, ...
% LMANinactivation, RestOfMultiTargExpts, DayOfLastThreshAdjust};

% Note:
% NumTargets: 
    % 1 (just one dir)
    % 2 (start off two targs
    % 0 - first 1, then 2 (diff dir) 
    % 3 - first 1, then 2 (same dir)
    % 4 - first 1, then 2 (diff dir), then 2 (same dir)
    % 5 - first 1, then 2 (same), then 2 (diff)
    % 9 - no learning, just baseline [in PROGRESS]
    
    
% "bidir" is actually multidir, just enter additional syls if you want
% multidir
% "consolidation" dates are mainly important for determing what direction
% the syl should be learning.
% RestOfMultiTargExpts is cell containing rest of multi target epochs
    % (either same dir or diff dir), only use if there is more than one
    % multidir epoch. If only one, then use the first cell that is inputed
    % earlier. If there are multiple, use that first cell for the first one,
    % and this cell for later ones.
    % e.g. RestOfMultiTargExpts = {one day before start epoch, last day of
    % epoch, {syl1, syl2, ...} ...} i.e. 3 entries for each epoch.
    
% DayOfLastThreshAdjust=[a b]; a: excluding minor adjust (can use); col 2:
% actual, even with minor adjust; note: days are starting from 1st day with
% usable data.
    

% for example:
% RegularExpressionsCell={'abbccbb', 'dccbb'}; OPTIONAL, only if will have
% to perform reg expr analysis.
% LMANinactivation= 1 or 0; (if does not exist, then is 0)

ExperimentList{1}={'pu53wh88','SyntaxDepPitchShift_abDowncbUp',2,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift/SeqFilterCompile',...
    '16Jul2014','21Jul2014', {'abbb', 'accbb'}, {'07Jul2014','16Jul2014', 'aB', 'cB'}, 0, {}, []};

% ExperimentList{22}={'pu53wh88','SyntaxDepPitchShift_cbUp',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift_cbUp/SeqFilterCompile',...
%     '04Aug2014','06Aug2014', {'abbb', 'accbb'}, {}, 0}; % THROW OUT AS
%     MAX LEARNING ONLY 50HZ/ONLY 10 days since previous learning

ExperimentList{2}={'pu53wh88','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '06Nov2014','11Nov2014', {'abbb', 'accbb'}, {}, 0, {}, [2 3]};

ExperimentList{3}={'pu53wh88','SeqDepPitchShift2',1','/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '04Dec2014','11Dec2014', {'abbb', 'accbb'}, {}, 0, {}, [4 4]};

ExperimentList{4}={'pu53wh88','SeqDepPitchShift3',0,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '21Feb2015','26Feb2015', {'abbb', 'accbb'}, {'01Mar2015', '07Mar2015', 'abB', 'abbB'}, 0, {}, [2 2]};


ExperimentList{5}={'pu11wh87','SyntaxDepPitchShift_abUP',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_abUP/SeqFilterCompile',...
    '12Jul2014','15Jul2014', {'abbccbb', 'dccbb'}, {}, 0, {}, [4 4]};

ExperimentList{6}={'pu11wh87','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '03Dec2014','11Dec2014', {'abbccbb', 'dccbb'}, {'11Dec2014', '20Dec2014','dccB', 'bccB'}, 0, {}, [4 4]};

ExperimentList{7}={'pu11wh87','SeqDepPitchShift3',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '26Feb2015','06Mar2015', {'abbccbb', 'dccbb'}, {}, 0, {}, [8 8]}; % note, no bidir, actually switched target after end of unidir learning.

ExperimentList{8}={'pu11wh87','SyntaxDepPitchShift_cbDOWN',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_cbDOWN/SeqFilterCompile',...
    '28Jul2014','03Aug2014', {'abbccbb', 'dccbb'}, {}, 0, {}, [3 3]}; % note, only analyzing abbccbb motif (ignore dccbb data entirely)

ExperimentList{9}={'pu11wh87','SeqDepPitchShift',2,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '02Nov2014','06Nov2014', {'abbccbb', 'dccbb'},{'29Oct2014', '02Nov2014', 'aB', 'bccB'}, 0, {}, []}; % note, only analyzing abbccbb motif, ignore dccbb motif entirely (they look the same)



ExperimentList{10}={'pu37wh20','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu37wh20/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '05Nov2014','06Nov2014', {'abbbg', 'abccbbg','qdh'}, {}, 0, {}, [5 5]};

ExperimentList{11}={'pu37wh20','SeqDepPitchShift2',1,'/bluejay3/lucas/birds/pu37wh20/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '25Nov2014','27Nov2014', {'abbbg', 'qccbbg','abccbbg','qdh'}, {}, 0, {}, [5 5]};

ExperimentList{12}={'pu64bk13','SeqDepPitchShift',0,'/bluejay4/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '17Dec2014','21Dec2014', {'dbbb'}, {'21Dec2014', '28Dec2014', 'dbB', 'dbbB'}, 0, {}, [3 3]};

ExperimentList{13}={'pu64bk13','SeqDepPitchShift2',0,'/bluejay4/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '17Feb2015','22Feb2015', {'dbbb'}, {'28Feb2015', '07Mar2015', 'dB', 'dbB'}, 0, {}, [3 3]};


ExperimentList{14}={'gr41gr90','SeqDepPitchShift',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '15Feb2015','24Feb2015', {'jbbacbb'}, {'28Feb2015', '07Mar2015', 'jbbacbB', 'jbBa'}, 0, {}, [2 2]}; % NOTE - need to separate all gr41 expts to j motif only

ExperimentList{15}={'gr41gr90','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '18Apr2015','23Apr2015', {'jbbacbb'}, {'10May2015', '19May2015', 'jbbacB', 'jbbacbB', 'jBba', 'jbBa'}, 0,{},  [2 2]};


ExperimentList{16}={'rd23gr89','SeqDepPitch',1,'/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '09Jun2015','14Jun2015', {'nakj','nahj', 'ndbbgcb'}, {}, 0, {}, [4 6]};


ExperimentList{17}={'pu35wh17','SeqDepPitch',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '10Jun2015','14Jun2015', {'abbb', 'jbbb'}, {}, 0, {}, [4 4]};

ExperimentList{18}={'pu35wh17','SeqDepPitch2',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '29Jun2015','30Jun2015', {'abbb', 'jbbb'}, {}, 0, {}, [3 3]};

ExperimentList{19}={'pu35wh17','SeqDepPitch3',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
    '12Jul2015','13Jul2015', {'abbb', 'jbbb'}, {}, 0, {}, [11 11]};


ExperimentList{20}={'rd12pu6','SeqDepPitch',0,'/bluejay4/lucas/birds/rd12pu6/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '07Aug2015','11Aug2015', {'jabbbhk', 'jjbdaabdgg'}, {'11Aug2015','16Aug2015','aaB','jB'}, 0, {}, [3 3]};

ExperimentList{37}={'rd12pu6','SeqDepPitch2',0,'/bluejay4/lucas/birds/rd12pu6/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '21Sep2015','29Sep2015', {'jabbbhk', 'jjbdaabdgg'}, {'29Sep2015','09Oct2015','aaB','aBb'}, 0, {}, [3 3]};

ExperimentList{43}={'rd12pu6','SeqDepPitch3',0,'/bluejay4/lucas/birds/rd12pu6/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
    '14Nov2015','22Nov2015', {'jabbbhk', 'jjbdaabdgg'}, {'22Nov2015','28Nov2015','aBb','aaB'}, 0, {}, [2 2]};


ExperimentList{38}={'rd28pu64','SeqDepPitch',0,'/bluejay4/lucas/birds/rd28pu64/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '17Nov2015','22Nov2015', {'jabdgkjbb', 'jjbb'}, {'22Nov2015','28Nov2015','kjB','jjB'}, 0, {}, [5 5]}; % used to be {'jabdgkjb+g', 'jjb+g'}
 



% ======================================================================
% ---------------------- LMAN EXPERIMENTS
ExperimentList{21}={'pu11wh87','SeqDepPitchLMAN', 0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '27May2015','08Jun2015', {'abbccbb', 'dccbb'}, {'08Jun2015', '16Jun2015', 'bccB', 'dccB'}, 1, {}, [9 9]}; % 
    % NEED TO USE EARLIER DAY THAN what I called consolidation start (based
    % on stability, about 8 days post 1st musc day (which is about 8 days post WN start), although they look
    % similar
    % - need to clean up labeling, as, likely because of WN, I get some
    % outliers on each day
    % - NEED TO LABEL CATCH SONG for bidir (that's when started catch song
    % method.) ONce do that, turn off skipping of post syls.
    % NOT CATCH

ExperimentList{22}={'pu11wh87','SeqDepPitchLMAN2', 0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '07Jul2015','12Jul2015', {'abbccbb', 'dccbb'}, {'14Jul2015', '25Jul2015', 'aB', 'abB'}, 1, {}, [4 4]}; % LABEL/MUSC GOOD - but need to account for WN during single dir.  
    % remove post syl (how many? abB is target. I will start just removing bC). SOME DAYS could use more labels, but important days
    % are relatively well labeled. Consol start is earliest MUSC day
    % rise: 8 days
    % hold: 8
    % bidir rise: 3
    % bidir hold: 8
    % NOT CATCH (SINGLE DIR)
    
ExperimentList{23}={'pu11wh87','SeqDepPitchLMAN3', 1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN3/SeqFilterCompile',...
    '17Aug2015','24Aug2015', {'abbccbb', 'dccbb'}, {}, 1, {}, [7 7]}; 

ExperimentList{24}={'pu53wh88','SeqDepPitchLMAN', 1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '23Jun2015','24Jun2015', {'abbb', 'accbb'}, {}, 1, {}, [12 12]}; % CONFIRM GOOD LABELING/INACTIVATION [to do: day 19 and 21 look potentially good.  need more labels late (after 3.5 hr) for those days.
%     the late days look good earlier (e.g. 2hrs post) (with higher musc).
%     Either way bird learned slowly and inactivation days are ~wk and more
%     during learning, so hard to compare to other experiments.

ExperimentList{25}={'gr41gr90','SeqDepPitchLMAN', 1,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '23Jul2015','24Jul2015', {'jbbacbb'}, {}, 1, {}, [3 3]}; % NOTE - LABEL/MUSC checked good. PROBLEM: FP on couple other B was high.

ExperimentList{26}={'gr41gr90','SeqDepPitchLMAN2', 0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '15Aug2015','17Aug2015', {'jbbacbb', 'dbbacbb'}, {'17Aug2015', '24Aug2015', 'jbBa', 'jbbacbB'}, 1, {}, [3 3]}; 
    % NOTE: bidir, syl 2 was actually any cbB (incluyding in motif after d
    % or g. need to limit bidir analysis to just one motif.

ExperimentList{35}={'rd23gr89','SeqDepPitchLMAN', 0, '/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '05Oct2015','06Oct2015', {'nakj','nahj', 'ndbbgcb'}, {'06Oct2015', '14Oct2015', 'dB', 'dbB'}, 1, {}, [2 2]};    
   

ExperimentList{36}={'pu11wh87','SeqDepPitchLMAN6', 0, '/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN6/SeqFilterCompile',...
    '05Oct2015','07Oct2015', {'abbccbb', 'dccbb'}, {'07Oct2015','14Oct2015', 'aB','abB'}, 1,{},  [2 2]};    
    

ExperimentList{44}={'bk34bk68','SeqDepPitchLMAN', 0, '/bluejay4/lucas/birds/bk34bk68/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '18Dec2015','20Dec2015', {'njjjbbga', 'kljbbga'}, {'21Dec2015','31Dec2015', 'ljbB','jjbB'}, 1, {}, [5 6]};    

ExperimentList{45}={'rd28pu64','SeqDepPitchLMAN', 0, '/bluejay4/lucas/birds/rd28pu64/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '18Dec2015','20Dec2015', {'jabdgkjbb', 'jjbb'}, {'21Dec2015','31Dec2015', 'kjB','jjB'}, 1, {}, [4 5]}; % used to be {'jabdgkjb+g', 'jjb+g'}


% ------ IN PROGRESS BELOW
% ExperimentList{46}={'wh4wh77','SeqDepPitchLMAN', 4, '/bluejay4/lucas/birds/wh4wh77/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
%     '14Jan2016','16Jan2016', {'jdbbghn', 'jakccbbghn'}, {'22Jan2016', '03Feb2016', 'cbB','cB'}, 1, {'03Feb2016', '16Feb2016', {'cbB','cB'}}}; %
 ExperimentList{46}={'wh4wh77','SeqDepPitchLMAN', 4, '/bluejay4/lucas/birds/wh4wh77/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '14Jan2016','16Jan2016', {'jdbbghn', 'jakccbb'}, {'22Jan2016', '03Feb2016', 'cbB','cB'}, 1, {'03Feb2016', '16Feb2016', {'cbB','cB'}}, [5 6]}; %

% ExperimentList{47}={'rd28pu64','SeqDepPitchLMAN2', 1, '/bluejay4/lucas/birds/rd28pu64/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
%     '21Jan2016','26Jan2016', {'jabdgkjbb', 'jjbb'}, {}, 1}; % CHANGE CONSOL DATES/ ADD BIDIR
ExperimentList{47}={'rd28pu64','SeqDepPitchLMAN2', 3, '/bluejay4/lucas/birds/rd28pu64/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '21Jan2016','26Jan2016', {'jabdgkjbb', 'jjbb'}, {'12Feb2016', '25Feb2016', 'jjB','kjB'}, 1, {}, [5 5]}; %

% ExperimentList{48}={'rd23gr89','SeqDepPitchLMAN2', 1, '/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
%     '23Jan2016','26Jan2016', {'ak','ah', 'dbbgcb'}, {}, 1}; % CHANGE CONSOL DATES/ ADD BIDIR
ExperimentList{48}={'rd23gr89','SeqDepPitchLMAN2', 3, '/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '23Jan2016','24Jan2016', {'nakj','nahj', 'ndbbgcb'}, {'03Feb2016','14Feb2016' , 'cB','dB', 'dbB'}, 1, {}, [7 7]}; % 

% ExperimentList{49}={'wh25pk77','SeqDepPitchLMAN', 1, '/bluejay4/lucas/birds/wh25pk77/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
%     '23Jan2016','26Jan2016', {'jjkbbg', 'ajnhbbg'}, {'25Jan2016', '31Jan2016', 'hbB', 'kbB'}, 1}; % CHANGE CONSOL DATES/ ADD BIDIR
ExperimentList{49}={'wh25pk77','SeqDepPitchLMAN', 4, '/bluejay4/lucas/birds/wh25pk77/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '23Jan2016','26Jan2016', {'jjkbbg', 'ajnhbbg'}, {'25Jan2016', '04Feb2016', 'hbB', 'kbB'}, 1, {'04Feb2016', '14Feb2016', {'hbB','kbB'}}, [4 4]}; % CHANGE CONSOL DATES/ ADD BIDIR

% ExperimentList{50}={'bk34bk68','SeqDepPitchLMAN3', 1, '/bluejay4/lucas/birds/bk34bk68/seq_dep_pitch_SeqDepPitchLMAN3/SeqFilterCompile',...
%     '26Jan2016','27Jan2016', {'jjbbga', 'ljbbga'}, {}, 1}; % CHANGE CONSOL DATES/ ADD BIDIR
ExperimentList{50}={'bk34bk68','SeqDepPitchLMAN3', 5, '/bluejay4/lucas/birds/bk34bk68/seq_dep_pitch_SeqDepPitchLMAN3/SeqFilterCompile',...
    '26Jan2016','27Jan2016', {'njjjbbga', 'kljbbga'}, {'01Feb2016', '13Feb2016', 'jjbB','ljbB'}, 1, {'13Feb2016', '29Feb2016', {'jjbB','ljbB'}}, [4 4]}; % CHANGE CONSOL DATES/ ADD BIDIR

%% ===================== ABORTED EXPERIMENTS, FOR HAMISH ANALYSIS (AFP BIAS PREDICT FIRST DAY LAERNING?)

ExperimentList{51}={'pu11wh87','SeqDepPitchLMAN4', 1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN4/SeqFilterCompile',...
    '08Sep2015','08Sep2015', {'abbccbb', 'dccbb'}, {}, 1, {}, [1 1]}; % NOTE - LABEL/MUSC checked good. PROBLEM: FP on couple other B was high.

ExperimentList{52}={'pu11wh87','SeqDepPitchLMAN5', 1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN5/SeqFilterCompile',...
    '19Sep2015','19Sep2015', {'abbccbb', 'dccbb'}, {}, 1, {}, [1 1]}; % NOTE - LABEL/MUSC checked good. PROBLEM: FP on couple other B was high.

ExperimentList{53}={'bk34bk68','SeqDepPitchLMAN2', 1,'/bluejay4/lucas/birds/bk34bk68/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '12Jan2016','12Jan2016', {'njjjbbga', 'kljbbga'}, {}, 1, {}, [1 1]}; % NOTE - LABEL/MUSC checked good. PROBLEM: FP on couple other B was high.


%% ======================================================================
%% REPEATS - only select these if you want to analyze repeats
ExperimentList{27}={'or100pu10','SeqDepPitch',0,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '10Jun2015','15Jun2015', {'acb+'}, {'15Jun2015','17Jun2015', 'acbB', 'acbbB'}, 0, {}, []};

ExperimentList{28}={'or100pu10','SeqDepPitch2',1,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '01Jul2015','02Jul2015', {'acb+'}, {}, 0, {}, []};

ExperimentList{29}={'or100pu10','SeqDepPitch3',0,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
    '24Jul2015','25Jul2015', {'acb+'}, {'25Jul2015','01Aug2015','acbbB','acbB'}, 0, {}, []};

% ExperimentList{30}={'or100pu10','SeqDepPitch3_duplicate',0,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
%     '24Jul2015','25Jul2015', {'acb+'}, {'03Aug2015','08Aug2015','acbbB','acB'}, 0}; % COPY OF ABOVE, BUT SECOND BIDIR

ExperimentList{31}={'or73pu40','SeqDepPitch',1,'/bluejay4/lucas/birds/or73pu40/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '12Jun2015','15Jun2015', {'acb+'}, {}, 0, {}, []};

ExperimentList{32}={'or73pu40','SeqDepPitch2',1,'/bluejay4/lucas/birds/or73pu40/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '29Jun2015','08Jul2015', {'acb+'}, {}, 0, {}, []}; % NOTE, SAME DIR MULTI SYL AFTER CONSOLID END

ExperimentList{33}={'or73pu40','SeqDepPitch3',1,'/bluejay4/lucas/birds/or73pu40/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
    '01Aug2015','05Aug2015', {'acb+'}, {}, 0, {}, []}; % NOTE, SAME DIR MULTI SYL AFTER CONSOLID END

ExperimentList{34}={'bu11or68','RepDepPitchShift',1,'/bluejay4/lucas/birds/bu11or68/seq_dep_pitch_RepDepPitchShift/SeqFilterCompile',...
    '04Dec2014','08Dec2014', {'jb+'}, {}, 0, {}, []}; % NOTE, SAME DIR MULTI SYL AFTER CONSOLID END

ExperimentList{41}={'gr87bu18','SeqDepPitch2',1,'/bluejay1/lucas/birds/gr87bu18/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '04Oct2015','07Oct2015', {'kb', 'ghb+', 'dccc', 'ab+'}, {}, 0, {}, []}; % 

ExperimentList{42}={'gr87bu18','RepDepPitch',1,'/bluejay1/lucas/birds/gr87bu18/seq_dep_pitch_RepDepPitch/SeqFilterCompile',...
    '11Nov2015','12Nov2015', {'kb', 'ghb+', 'dccc', 'ab+'}, {}, 0, {}, []}; % 


% ======== LMAN REPEAT EXPERIMENTS
ExperimentList{39}={'or100pu10','SeqDepPitchLMAN',0,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '05Oct2015','09Oct2015', {'acb+'}, {'09Oct2015','14Oct2015','acbB','acbbB'}, 1, {}, []};

ExperimentList{40}={'or100pu10','SeqDepPitchLMAN2',1,'/bluejay4/lucas/birds/or100pu10/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '16Nov2015','17Nov2015', {'acb+'}, {}, 1, {}, []};



%% =========== RAW DATA ANALYSIS
%% ====== Recalculate acoustic features, to get things that don't have well-defined pitch

lt_seq_dep_pitch_ACROSSBIRDS_RAW_Acoustic

%% ============ RECALC ACOUSTIC (USING POWER INSTEAD OF AMPLITUDE FOR ENTROPY ANALYSES)

lt_seq_dep_pitch_ACROSSBIRDS_RAW_Acoustic2(ExperimentList)
% NOTE: have done for 


%% ============

%% PREPROCESS

PARAMS.global.ExperimentList=ExperimentList;

close all;
% -- acoustic stuff
PARAMS.preprocess.plotPCA_stuff=0; % default: 0
PARAMS.preprocess.Zscore_use_rends=1; % default 1 - takes zscore over all rends (if 0, then over all syl means) % NOTE: it makes almost no difference whether use rends or means, so just use rends.
PARAMS.preprocess.load_old_vector_space=1; % default, 0, recreates vector space on each run. if 1, then loads previous version.
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PREPROCESS(ExperimentList, PARAMS);


% =============== FILTER OUT EXPERIMENTS THAT START WITH 2-DIR 
% SeqDepPitch_AcrossBirds_ORIG=SeqDepPitch_AcrossBirds;
% Remove expts starting with 2 targs;
filter='RemoveStartTwoTargs';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


%% ============= KEEP/REMOVE REPEATS - DON'T RUN THIS IF DON'T CARE ONE WAY OR ANOTHER [USING HAND LABELED]
% TO DO: REMOVE EXPERIMENTS IF ALL SYLS HAPPEN TO BE REMOVED
remove_repeats=9; % removes all repeats in targ motif.
extract_repeats=0;
% lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);


% TO DO: REMOVE EXPERIMENTS IF ALL SYLS HAPPEN TO BE REMOVED
remove_repeats=10; % if targ is in repeat, remove syls in that repeat coming after targ, but keep all other syls
extract_repeats=0;
% lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);


% TO DO: REMOVE EXPERIMENTS IF ALL SYLS HAPPEN TO BE REMOVED
remove_repeats=3; % [syl pos 3 +] a) throws out all experiments if target was syl3 or 
% higher in repeat, and b) For all experiments, throw out syls that are at repeat position 3+ [EVEN IF IS DIFF TYPE SYL] [CAN END UP WITH EXPTS WITH NO/FEW NONTARGS]
extract_repeats=0;
% lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);




%% +++++++++++++++++++++++ SHOULD DO
%% =============== [OPTIONAL - EITHER DO HERE OR AFTER REMOVE REPEATS - REDEFINE SYL SIMILARITY AND SEQUENCE
% SIMILARITY USING ACOUSTIC DISTANCE]
% ==== ACOUSTIC SIMILARITY STUFF
close all;
OverwriteSylIDs=1; % replaces syl and presyl similar (1 or 0) - saves old in new field
rd28_AdHoc=0; % then doesn't have acoustic distance yet, use hand labeled (only for presyl, syl distance is still automatic)
% dimension=7;
dimension=8;
TestIfDistanceSignDiff=0;
% DIST: 1.43
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_AcousticAnaly(SeqDepPitch_AcrossBirds, PARAMS, OverwriteSylIDs, rd28_AdHoc, dimension, TestIfDistanceSignDiff);

% === compile same-type and diff type unique syls
SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_CompSyls(SeqDepPitch_AcrossBirds);

%% == RECALCULATE LEARNING AND DRIFT FROM 2 BASELINE DAYS 
close all
OverwriteLearningMetric=1; % OVERWRITES LEARNING METRIC
ConvertToHz=1; % 1, uses FF diff, not zscore; 0, uses zscore [default];
SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_RecalcBaseline(SeqDepPitch_AcrossBirds, PARAMS, OverwriteLearningMetric, ConvertToHz);
% lt_seq_dep_pitch_ACROSSBIRDS_RecalcBaseline(SeqDepPitch_AcrossBirds, PARAMS,OverwriteLearningMetric);

%% ==== PLOT ALL LEARNING DISTRIBUTIONS IN RASTER PLOT
% also outputs percentiles of baseline drift, and distribution

close all
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LearnRaster(SeqDepPitch_AcrossBirds, PARAMS);

%% ==== AD HOC THINGS SHOULD DO. [ONLY MATTERS FOR REPEATS]
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_AdHocChanges(SeqDepPitch_AcrossBirds, PARAMS);


%% +++++++++++++++++++++++++++++++++++++++++++++++++==


%% ========== plots distribition of all pairwise acoustic distances across all experiments
close all
IncludeThrownOutSyls=0; % i.e. those that were specificed above, e.g. noisy PC
threshold=1.43; % based on minimum of distribution
lt_seq_dep_pitch_ACROSSBIRDS_AllPairAcDist(SeqDepPitch_AcrossBirds, PARAMS, IncludeThrownOutSyls, threshold);



%% VARIOUS SANITY CHECKS

% open and run what you need
lt_seq_dep_pitch_ACROSSBIRDS_SanityChecks


%% PLOT all learning trajectories
close all
lt_seq_dep_pitch_ACROSSBIRDS_PlotAllTraj(SeqDepPitch_AcrossBirds, PARAMS);
close all;
lt_seq_dep_pitch_ACROSSBIRDS_PlotAllTraj(SeqDepPitch_AcrossBirds_filtered, PARAMS);

%% === plot each expt in syl acoustic space (pca space)



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% PLOTS
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ===== PLOT SPECTROGRAMS OF MOTIFS, SYLS, FOR EACH BIRD
close all;
plotForIllustrator=1; % then greater resolution

% == to plot specific or all birds/expts [for all, keep both entries as '')
PARAMS.PlotSpec.BirdToPlot='bk34bk68'; % make this empty to plot all birds/expts
PARAMS.PlotSpec.ExptToPlot='SeqDepPitchLMAN';
PARAMS.PlotSpec.PlotAcousticDist=1;
PARAMS.PlotSpec.PlotGeneralization=0;
% --- syls
PARAMS.PlotSpec.NumPreNotes=3;
PARAMS.PlotSpec.NumPostNotes=0;
% --- motifs
PARAMS.PlotSpec.NumRendsToPlot=4; % make 0 to skip plotting motifs
lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec2(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator)

%% PLOT RAW (FFVALS) ACROSS DAYS FOR GIVEN EXPERIMENT (OR ALL)
close all
% BirdToPlot='rd28pu64'; % bird or 'all'
% ExptToPlot='SeqDepPitchLMAN';
% % SylsToPlot={'aBb','aaB','jB','daA','Abb', 'gG'}; % either cell array, or 'all';
% SylsToPlot={'aB'}; % either cell array, or 'all';
BirdToPlot='rd12pu6'; % bird or 'all'
ExptToPlot='SeqDepPitch3';
% SylsToPlot={'aBb','aaB','jB','daA','Abb', 'gG'}; % either cell array, or 'all';
SylsToPlot={'aBb','aaB', 'jB', 'dG'}; % either cell array, or 'all';
use_rand_color=0; % if 1, then random color for syls, if 0, then use color scheme (targ, etc)

UseRecalcBaseline=1; % then limits to days used to calcualte deviation of learning

ThreshChangeTime={'rd12pu6', 'SeqDepPitch3', 'aBb', ... % quadruplets {birdname, experiment, syl, {{11Nov2015-0000, 3620}, {11Nov2015-1322, 3680}}} 
    {{'11Nov2015-0000', 3620}, {'11Nov2015-1322', 3680},...
    {'11Nov2015-2400', 3725}, {'12Nov2015-1355', 3770}, ...
    {'12Nov2015-2400', 3820}}, ...
    'rd12pu6', 'SeqDepPitch3', 'aaB', ...
    {{'22Nov2015-2400', 3815}, {'23Nov2015-1003', 3650}, ...
    {'24Nov2015-1258', 3625}, {'24Nov2015-2400', 3640}, ...
    {'25Nov2015-1418', 3600}, {'26Nov2015-1400', 3560}, ...
    {'26Nov2015-1421', 3565}}, ...
    };

lt_seq_dep_pitch_ACROSSBIRDS_PlotRawTraj(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, use_rand_color, UseRecalcBaseline, ThreshChangeTime)


%% === [MUSC AND PBS] PLOT RAW FFVALS, INCLUDING BOTH PBS AND MUSC, FOR AS MANY SYLLABLES AS DESIRED
close all;
% BirdToPlot='gr41gr90';
% ExptToPlot='SeqDepPitchLMAN2';
% SylsToPlot={'jbBa','jbbacbB'};
% BirdToPlot='rd23gr89';
% ExptToPlot='SeqDepPitchLMAN';
% SylsToPlot={'dB','dbB', 'cB'};
% BirdToPlot='pu11wh87';
% ExptToPlot='SeqDepPitchLMAN3';
% SylsToPlot={'bccB', 'dccB'};
% BirdToPlot='rd23gr89';
% ExptToPlot='SeqDepPitchLMAN2';
% SylsToPlot={'dB', 'cB'};
BirdToPlot='pu11wh87';
ExptToPlot='SeqDepPitchLMAN';
SylsToPlot={'bccB', 'dccB'};
BirdToPlot='bk34bk68';
ExptToPlot='SeqDepPitchLMAN3';
SylsToPlot={'jjbB', 'ljbB'};
BirdToPlot='rd23gr89';
ExptToPlot='SeqDepPitchLMAN';
SylsToPlot={'dB'};

BirdToPlot='rd23gr89';
ExptToPlot='SeqDepPitchLMAN2';
SylsToPlot={'Ah', 'Ak', 'k' , 'd' ,'dB' ,'dbB' ,'g' ,'c' ,'cB'};

BirdToPlot='pu53wh88';
ExptToPlot='SeqDepPitchLMAN';
SylsToPlot={'dB'};

% -- use below to plot all LMAN experiments
BirdToPlot='';
ExptToPlot='';
SylsToPlot={'same'};
% -- use below to plot all LMAN experiments
overlayMeans=1;
plotRawFF=0; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
UseSylColors=0; % 0 is black;
flipsign=1; % plot pos, if neg
use_std=0; % applies to mean plots only!! (std, instead of sem)

% ---- NO LMAN STATS
OverlayLMANStats=0;
lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats)

% -- OVERLAY LMAN LEARNING STATS? - can only do this if using
% SeqDepPitch_AcrossBirds_LMAN (after do filter below)
OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
OverlayMUSC_days=[];
OverlayMUSC_days=[20 21]; % plots mean shift (MUSC and PBS) on user-specified days (if empty, then uses learning window).
    % (if not empty, then takes precedence over OverlayLMANStats);
    % OverlayLMANStats needs to be 1 for this to work. format: [start end]
    % (day inds)
lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds_LMAN, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, OverlayMUSC_days)

% --- TO PLOT WITH CV AND ALL RAW LMAN STUFF
plotRawFF=1; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
OverlayMUSC_days=[];
lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds_LMAN, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, OverlayMUSC_days)



%% PLOT BARS OVER MOTIFS ACROSS DAYS - GOOD FOR REPEATS
close all;
DaysToPlot=[3 4];
% DaysToPlot=[];
plot_hit_rate=0;
equal_y_scale=1;
AnnotateSylPositions=1;

lt_seq_dep_pitch_ACROSSBIRDS_MotifInputDays(SeqDepPitch_AcrossBirds, PARAMS, DaysToPlot, plot_hit_rate, equal_y_scale, AnnotateSylPositions)
lt_seq_dep_pitch_ACROSSBIRDS_MotifInputDays(SeqDepPitch_AcrossBirds_REPEATS, PARAMS, DaysToPlot, plot_hit_rate, equal_y_scale, AnnotateSylPositions)


%% PATTERNS OF LEARNING [TAILORED FOR REPEATS!!]
% - plots bunch of things, like PCA space, learnnig, acoustic. Q: for
% repeats, is pattern explained by acoustic dist?

close all;
DaysToPlot=[3 4];
plot_hit_rate=0;
equal_y_scale=1;
AnnotateSylPositions=1;
plotSizeForIllustrator=1; % then 2 x 2 subplots (not 4 x 2)

remove_repeats=0;
extract_repeats=1;
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);

lt_seq_dep_pitch_ACROSSBIRDS_MotifInputDays_v2(SeqDepPitch_AcrossBirds, PARAMS, DaysToPlot, plot_hit_rate, equal_y_scale, AnnotateSylPositions, plotSizeForIllustrator)


%% ==== TRIAL TO TRIAL LEARNING/ DAY VS. OVERNIGHT?
% NOTE: FLIP SIGN OF SLOPE TO COMPARE ACROSS EXPERIMENTS

close all;
BirdToPlot='pu37wh20';
ExptToPlot='SeqDepPitchShift2';
OnlyWellLearned=1;

lt_seq_dep_pitch_ACROSSBIRDS_TrialByTrial(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, OnlyWellLearned);


%% ====== [SAM SOBER] PLOT TRIAL-BY-TRIAL LEARNING, COMPARE TRAJECTORY RELATIVE TO LEARN MAGNITUDE, ETC (A LA SOBER)
% NOTE: in progress. see keyboard section within code. Interpolating.
close all;
HourBins=[9:2:19]; % in 24 hour clock.
lt_seq_dep_pitch_ACROSSBIRDS_Trial2(SeqDepPitch_AcrossBirds, PARAMS, HourBins)

% ------------------------- PARAMS TO PLOT LIKE SOBER
% ------------ PLOT ONLY LOW LEARNING EXPERIMENTS
% -- USE -2 AND 2 FOR Z-SCORE RANGE
filter='learning_range';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
filter='did_not_sing2';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);
% ------------------------

close all;
HourBins=[9 19]; % in 24 hour clock.
plotCents=1;
lt_seq_dep_pitch_ACROSSBIRDS_Trial2(SeqDepPitch_AcrossBirds_filtered, PARAMS, HourBins, plotCents)


%% =============== PLOT MEAN TRAJECTORY DAY TO DAY


% ------------------------------------ LOW LEARNERS
% ---- use -2 to 1.15; 1.15 to 2.4; 2.4 to 5; [REPEATS]
% ---- [-0.15 1.22 2 2.5 4.2]; [NO REPEATS]
% ---- [-100 to 66.56], [66.56 to 134.38], [134.38 400] ----  [HZ, NO REPEATS (MAIN DATASET)]
filter='learning_range';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
% --- 
filter='did_not_sing2';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);


close all;
OnlyExptsWithNoStartDelay=0; % removes expts with that start with delay from WN to learning.
TakeIntoAccountStartDelay=1; % then will start WN1 from 1st day with data
plotCents=1;
plotZ=0;
DayWindow=[-3 4]; % [baseDays, WNdays];
throwOutIfAnyDayEmpty=0; % then if any day is nan, throws out just that syl (but will likely apply to all syls in an expt)
recalcBaseMean=0;
useHandLabSame=1; % default = 0; if 1, then uses hand lab. if 0, comp labeled.

lt_seq_dep_pitch_ACROSSBIRDS_MeanTraj(SeqDepPitch_AcrossBirds_filtered, PARAMS, ...
    OnlyExptsWithNoStartDelay, TakeIntoAccountStartDelay, plotCents, DayWindow, ...
    throwOutIfAnyDayEmpty, plotZ, recalcBaseMean, useHandLabSame)


%% ++++++++++++++++++++ SINGLE DIR LEARNING ANALYSIS
%% =====================

% Remove experiments where learning metric is "nan";
filter='learning_metric';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
% 40.34hz
% 0.8z

% Remove experiments where late learning is poor; (e.g. for
% antigeneralizatioj analysis, looking at last window of single targ
% training)
filter='AdHocExptRemove';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

% Only keep experiments where learning was within a certain range
filter='learning_range';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

% --------------- ONLY KEEP EXPERIMENTS THAT ARE WITHIN LEARNING RANGE, AND
% ALSO DON'T HAVE ANY MISSING DAYS FROM START OF LEARNING
filter='learning_range';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
filter='DelayExptStart';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);
% ------------------------



% ===== OPTIONS - removing experiments before analysis [I DON'T USE]
% remove if learning less than threshold - don't use
    filter='good_learning';
    [SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);
    % Remove LMAN experiments
    filter='notLMAN';
    [SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);
    % Remove experiments that are missing data in beginning of learning (i.e.
    % did not sing)
    filter='did_not_sing';
    [SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);
    % Remove experiments that are missing data at least 2 days in beginning of learning (i.e.
    % did not sing)
    filter='did_not_sing2';
    [SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_filtered, filter);


% ====== LEARNING ANALYSIS [NOTED INSIDE WHAT IS USING CONSOLD START VS.
% LEARNING METRIC]
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds_filtered, PARAMS);

% OPTIONAL, do without filtering out poor learners
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds, PARAMS);



% Option 2) ==== RUN THIS IF DON'T WANT TO FILTER
% ====== LEARNING ANALYSIS
% [SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds, PARAMS);


%% ===== linear model of generalization
close all
onlyLMANexpts=0; % only uses LMAN expts that pass learning filter
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_LinModel(SeqDepPitch_AcrossBirds_filtered, PARAMS, onlyLMANexpts);
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_LinModel(SeqDepPitch_AcrossBirds, PARAMS, onlyLMANexpts);


%% ==== SINGLE DIR timecourse (from start of WN, and during manually chosen "consolidation" period);
% TO DO: BE ABLE TO START FROM WN DAY 1, USING CONSOLID END AT END, AND
% TAKING INTO ACCOUNT EMPTY DAYS AT START

% IN PROGRESS [NOT USED. IS INCORPORATED INTO CONSOLID PERIOD BELOW]
DriveMoreList={'pu64bk13', 'SeqDepPitchShift2', 12, ...
    'gr41gr90', 'SeqDepPitchShift', 13, ...
    'gr41gr90', 'SeqDepPitchShift2', 17, ...
    'pu11wh87','SeqDepPitchShift2', 11,...
    'pu11wh87','SeqDepPitchLMAN2', 1, ...
    'pu53wh88', 'SeqDepPitchShift3',13}; % i.e. value triplets {birdname, exptname, LastDayBeforeDriveMore_WNdayInds}

% DONE:
% PARAMS.TimeCourse.ConsolidationList={...
%     'pu53wh88', 'SeqDepPitchShift', [8 13], ...
%     'pu53wh88', 'SeqDepPitchShift2', [10 16], ...
%     'pu53wh88', 'SeqDepPitchShift3', [8 13], ...
%     'pu53wh88', 'SeqDepPitchLMAN', [6 10], ... %
%     'pu11wh87', 'SyntaxDepPitchShift_abUP', [5 8], ... %
%     'pu11wh87', 'SeqDepPitchShift2', [11 19], ...
%     'pu11wh87', 'SeqDepPitchShift3', [12 20], ... %
%     'pu11wh87', 'SyntaxDepPitchShift_cbDOWN', [5 11], ...
%     'pu11wh87', 'SeqDepPitchLMAN', [15 28], ...
%     'pu11wh87', 'SeqDepPitchLMAN2', [8 16], ...
%     'pu11wh87', 'SeqDepPitchLMAN6', [4 5], ...
%     'pu37wh20', 'SeqDepPitchShift', [4 6], ...
%     'pu37wh20', 'SeqDepPitchShift2', [4 5], ...
%     'pu64bk13', 'SeqDepPitchShift', [4 8], ...
%     'pu64bk13', 'SeqDepPitchShift2', [5 12], ...
%     'gr41gr90', 'SeqDepPitchShift', [3 12], ...
%     'gr41gr90', 'SeqDepPitchShift2', [10 17], ...
%     'gr41gr90', 'SeqDepPitchLMAN', [4 8], ...
%     'gr41gr90', 'SeqDepPitchLMAN2', [3 6], ...
%     'rd23gr89', 'SeqDepPitch', [7 12], ...
%     'rd23gr89', 'SeqDepPitchLMAN', [4 5], ...
%     'pu35wh17', 'SeqDepPitch', [7 11], ...
%     'rd12pu6', 'SeqDepPitch', [5 9], ...
%     'rd12pu6', 'SeqDepPitch2', [3 8], ...
%     'rd12pu6', 'SeqDepPitch3', [4 12], ...
%     'rd28pu64', 'SeqDepPitch', [5 11], ...
%     'rd28pu64', 'SeqDepPitchLMAN', [6 6], ...
%     'bk34bk68', 'SeqDepPitchLMAN', [5 7]}; 
% OLD above

PARAMS.TimeCourse.ConsolidationList={...
    'pu53wh88', 'SeqDepPitchShift', [8 13], ...
    'pu53wh88', 'SeqDepPitchShift2', [8 16], ...
    'pu53wh88', 'SeqDepPitchShift3', [5 12], ...
    'pu53wh88', 'SeqDepPitchLMAN', [14 21], ... 
    'pu11wh87', 'SyntaxDepPitchShift_abUP', [5 8], ... 
    'pu11wh87', 'SeqDepPitchShift2', [11 19], ...
    'pu11wh87', 'SeqDepPitchShift3', [12 17], ... 
    'pu11wh87', 'SyntaxDepPitchShift_cbDOWN', [5 11], ...
    'pu11wh87', 'SeqDepPitchLMAN', [13 28], ...
    'pu11wh87', 'SeqDepPitchLMAN2', [5 15], ... % if want to avoid gap, use [8 15]
    'pu11wh87', 'SeqDepPitchLMAN3', [8 15], ... 
    'pu11wh87', 'SeqDepPitchLMAN6', [4 5], ...
    'pu37wh20', 'SeqDepPitchShift', [5 6], ...
    'pu37wh20', 'SeqDepPitchShift2', [4 5], ...
    'pu64bk13', 'SeqDepPitchShift', [4 8], ...
    'pu64bk13', 'SeqDepPitchShift2', [4 12], ...
    'gr41gr90', 'SeqDepPitchShift', [3 12], ...
    'gr41gr90', 'SeqDepPitchShift2', [3 14], ...
    'gr41gr90', 'SeqDepPitchLMAN', [5 8], ...
    'gr41gr90', 'SeqDepPitchLMAN2', [4 6], ...
    'rd23gr89', 'SeqDepPitch', [6 12], ...
    'rd23gr89', 'SeqDepPitchLMAN', [4 5], ...
    'rd23gr89', 'SeqDepPitchLMAN2', [8 14], ...
    'pu35wh17', 'SeqDepPitch', [6 11], ...
    'pu35wh17', 'SeqDepPitch2', [4 4], ... 
    'pu35wh17', 'SeqDepPitch3', [5 5], ... 
    'rd12pu6', 'SeqDepPitch', [4 9], ...
    'rd12pu6', 'SeqDepPitch2', [4 12], ...
    'rd12pu6', 'SeqDepPitch3', [4 12], ...
    'rd28pu64', 'SeqDepPitch', [6 11], ...
    'rd28pu64', 'SeqDepPitchLMAN', [6 6], ...
    'rd28pu64', 'SeqDepPitchLMAN2', [7 11], ... 
    'bk34bk68', 'SeqDepPitchLMAN', [7 7], ... 
    'bk34bk68', 'SeqDepPitchLMAN3', [5 8], ... 
    'wh4wh77', 'SeqDepPitchLMAN', [6 12], ...
    'wh25pk77', 'SeqDepPitchLMAN', [5 7]};
% triplets {birdname, exptname, [day1 daylast] (relative to WN day1 (not
% modified to remove gaps - i.e. WN day 1 = 1);
% UP TO DATE (for all expt, even if fail to pass learning thresh)



% AD HOC - STARTING FROM WN DAY 1 (AND TAKING INTO ACCOUNT EARLY WN NO DATA
% DAYS)
PARAMS.TimeCourse.ConsolidationList={...
    'pu53wh88', 'SeqDepPitchShift', [4 13], ...
    'pu53wh88', 'SeqDepPitchShift2', [3 16], ...
    'pu53wh88', 'SeqDepPitchShift3', [3 12], ...
    'pu53wh88', 'SeqDepPitchLMAN', [2 21], ... 
    'pu11wh87', 'SyntaxDepPitchShift_abUP', [1 8], ... 
    'pu11wh87', 'SeqDepPitchShift2', [1 19], ...
    'pu11wh87', 'SeqDepPitchShift3', [1 17], ... 
    'pu11wh87', 'SyntaxDepPitchShift_cbDOWN', [1 11], ...
    'pu11wh87', 'SeqDepPitchLMAN', [1 28], ...
    'pu11wh87', 'SeqDepPitchLMAN2', [1 15], ... % if want to avoid gap, use [8 15]
    'pu11wh87', 'SeqDepPitchLMAN3', [1 15], ... 
    'pu11wh87', 'SeqDepPitchLMAN6', [1 5], ...
    'pu37wh20', 'SeqDepPitchShift', [1 6], ...
    'pu37wh20', 'SeqDepPitchShift2', [1 5], ...
    'pu64bk13', 'SeqDepPitchShift', [1 8], ...
    'pu64bk13', 'SeqDepPitchShift2', [1 12], ...
    'gr41gr90', 'SeqDepPitchShift', [1 12], ...
    'gr41gr90', 'SeqDepPitchShift2', [1 14], ...
    'gr41gr90', 'SeqDepPitchLMAN', [2 8], ...
    'gr41gr90', 'SeqDepPitchLMAN2', [1 6], ...
    'rd23gr89', 'SeqDepPitch', [1 12], ...
    'rd23gr89', 'SeqDepPitchLMAN', [2 5], ...
    'rd23gr89', 'SeqDepPitchLMAN2', [1 14], ...
    'pu35wh17', 'SeqDepPitch', [1 11], ...
    'pu35wh17', 'SeqDepPitch2', [1 4], ... 
    'pu35wh17', 'SeqDepPitch3', [1 5], ... 
    'rd12pu6', 'SeqDepPitch', [1 9], ...
    'rd12pu6', 'SeqDepPitch2', [1 12], ...
    'rd12pu6', 'SeqDepPitch3', [1 12], ...
    'rd28pu64', 'SeqDepPitch', [1 11], ...
    'rd28pu64', 'SeqDepPitchLMAN', [1 6], ...
    'rd28pu64', 'SeqDepPitchLMAN2', [1 11], ... 
    'bk34bk68', 'SeqDepPitchLMAN', [1 7], ... 
    'bk34bk68', 'SeqDepPitchLMAN3', [1 8], ... 
    'wh4wh77', 'SeqDepPitchLMAN', [1 12], ...
    'wh25pk77', 'SeqDepPitchLMAN', [1 7]};



close all;
% DayToLockAllTo='WNday1';
DayToLockAllTo='ConsolidDay1';

% --- for plotting locked
PARAMS.TimeCourse.DaysToPlot=8; % num days past lock (leave [] if don't care. has to be defined in PlotOnlyGaplessExpt=1)
PARAMS.TimeCourse.PlotOnlyGaplessExpt=1; % only plots if expt has no holes in data for that period
PARAMS.TimeCourse.AccountForEmptyDays_WNon=1; % account for start of experiment
PARAMS.TimeCourse.NumDaysInBin=2; % for unpaired analysis.
PARAMS.TimeCourse.ThrowOutPoorLearner=0; % based on z score learning;


% DayToLockAllTo='SingleDirConsolid'; % locks all to start of consolid, and plots only things within time window
% PARAMS.TimeCourse.SingleDirConsolid= ...
% {'pu53wh88', 'SeqDepPitchLMAN', [19 21], ...
% 'rd23gr89', 'SeqDepPitchLMAN', [4 5]}; % triplets {bird, expt, [day1 last day] (days rel to WN on)

lt_seq_dep_pitch_ACROSSBIRDS_TimeCourse(SeqDepPitch_AcrossBirds, PARAMS, DayToLockAllTo)


%% PLOT LEARNING VS. DISTANCE FROM TARGET 
% Modified from singledir analysis, here  is simple, good for repeats.
close all;
norm_all_to_targ=1;
plot_error_bars_on_raw=0;
equal_y_axis=1;
lt_seq_dep_pitch_ACROSSBIRDS_LearnVsPosition(SeqDepPitch_AcrossBirds, PARAMS, norm_all_to_targ, plot_error_bars_on_raw, equal_y_axis);
lt_seq_dep_pitch_ACROSSBIRDS_LearnVsPosition(SeqDepPitch_AcrossBirds_filtered, PARAMS, norm_all_to_targ, plot_error_bars_on_raw, equal_y_axis);


%% ==== REPEATS (good code) 
close all;
% ========= SINGLE DIR
remove_repeats=0;
extract_repeats=1;
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats);

% RUN
lt_seq_dep_pitch_ACROSSBIRDS_PlotRepeats(SeqDepPitch_AcrossBirds, PARAMS);

% =========== LMAN
% SEE BELOW


% ======== MULTIDIR
close all;
PARAMS.global.MULTIDIR.DayBinSize=2; % num days to take at start and end. (also = num baseline days)
PARAMS.global.MULTIDIR.WNdaynum=[3 4]; % days to use, post bidir WN on
RepeatsOnly=1;
OnlyUseSylsInSylsUnique=1;

[SeqDepPitch_AcrossBirds_MULTIDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, PARAMS, RepeatsOnly, OnlyUseSylsInSylsUnique);




%% MULTIDIR LEARNING ANALYSES [HAVE NOT IMPLEMENTED LEARNING METRIC]
close all;
PARAMS.global.MULTIDIR.DayBinSize=2; % num days to take at start and end. (also = num baseline days)
PARAMS.global.MULTIDIR.WNdaynum=[3 4]; % days to use, post bidir WN on
RepeatsOnly=0;
OnlyUseSylsInSylsUnique=1; % 

[SeqDepPitch_AcrossBirds_MULTIDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, PARAMS, RepeatsOnly,OnlyUseSylsInSylsUnique);
 

% ==== PLOT ACROSS DAYS - CONSTRAINED TO SHIFT TOGETHER?
close all;
RepeatsOnly=0;
OnlyUseSylsInSylsUnique=1; % 
ExcludeSeqLearning=0;
ExcludeNotFullyLabeled=1; % ad hoc, expt with unresolvable holes. (e.g. reprobing, mistakes, see code)
ExcludeIfHasSameDirBeforeDiffDir=1;
ExcludeIfFirstTargDriveMore=1; % ad hoc, remove expt where first targ drove more.
DaysToPlot=[3 5]; % [num days before bidir day 1, num days during bidir (inclusize)]
[SeqDepPitch_AcrossBirds_MULTIDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2(SeqDepPitch_AcrossBirds, PARAMS, RepeatsOnly,OnlyUseSylsInSylsUnique, DaysToPlot, ExcludeSeqLearning, ExcludeNotFullyLabeled, ExcludeIfHasSameDirBeforeDiffDir, ExcludeIfFirstTargDriveMore);




%% LMAN LEARNING ANALYSIS
SeqDepPitch_AcrossBirds_LMAN=SeqDepPitch_AcrossBirds;


%% [HAMISH] ============= BASELINE LMAN BIAS PREDICT LEARNING?

close all;
plotExptRawDat=0;
lt_seq_dep_pitch_ACROSSBIRDS_Hamish(SeqDepPitch_AcrossBirds, PARAMS, plotExptRawDat)

%% ===================== OPTIONAL
% ++++++++++++++++++++++++++++++++++++++ BASELINE EFFECT OF INACTIVATION
% NOTE: this shows that baseline effect of musc is to bring same types
% closer together, with no effect on diff-types. don't see effect if do
% control analysis using early and late day.
% TO DO: perform within experiment controls - label days fully for those
% expts.
close all;
NullControl_splitday=0; 
same_type_thr=1.43;
[~, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbase(SeqDepPitch_AcrossBirds_LMAN, PARAMS, same_type_thr, NullControl_splitday);


% ++++++++++++++++++++++++++++++++++++ LEARNING
% 1) =============== PREPROCESS
% --- reextract MUSc data without normalizing for baseline musc effect?
[SeqDepPitch_AcrossBirds_LMAN, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANrenorm(SeqDepPitch_AcrossBirds_LMAN, PARAMS, musc_day_window_WNdayinds, musc_day_window_Bidirdayinds_Bidir, pu53_use_later_days, debugON);

%% MAIN LMAN STUFF
% --- what day windows to use for analysis?
close all;
debugON=0; % plots 
musc_day_window_WNdayinds=[4 10]; % from WN day 3 to 10, collect all for analysis
musc_day_window_Bidirdayinds_Bidir=[4 10];
% musc_day_window_WNdayinds=[9 12]; % from WN day 3 to 10, collect all for analysis
pu53_use_later_days=2; % if 1, then goes to late days that had inactivaiotn; if 2, then throws out pu53
single_targ_UseMaintainedEarly=1; % then use initial learning to day 4 of maintained shift [max is day 10]
% ONLY REQUIRE BELOW IF USING 1 FOR single_targ_UseMaintainedEarly [ I.E.
% these are consolidation windows for single-targ consolidation code below]
PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={...
    'pu11wh87', 'SeqDepPitchLMAN',[11 28], ... % avoid early as still adjusting
    'pu11wh87', 'SeqDepPitchLMAN2',[5 15], ... % one day skipped labeling, but is fine.
    'gr41gr90', 'SeqDepPitchLMAN2', [2 6], ...
    'rd23gr89',  'SeqDepPitchLMAN2', [7 14], ...
    'rd28pu64',  'SeqDepPitchLMAN2', [7 11], ... % gap in day
    'bk34bk68',  'SeqDepPitchLMAN3', [5 8], ...
    'bk34bk68',  'SeqDepPitchLMAN', [3 7], ...
    'wh4wh77',  'SeqDepPitchLMAN', [4 12], ...
    'wh25pk77',  'SeqDepPitchLMAN', [3 7]};
[SeqDepPitch_AcrossBirds_LMAN, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANpreprocess(SeqDepPitch_AcrossBirds_LMAN, PARAMS, musc_day_window_WNdayinds, musc_day_window_Bidirdayinds_Bidir, pu53_use_later_days, debugON, single_targ_UseMaintainedEarly);


% === THROW OUT BAD SYLS? i.e. those immediately post WN since labeled
% not-catch songs. (will keep for plots summarizing experiments, but will
% throw out for actual analyses);
remove_bad_syls=1;
% SylsToRemove_SingleDir= ...
%     {'pu11wh87','SeqDepPitchLMAN',{'bccbB'}, ...
%     'pu11wh87','SeqDepPitchLMAN2', {'bC', 'bcC'}, ...
%     'pu53wh88','SeqDepPitchLMAN', {'abbB'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
% SylsToRemove_BiDir= ...
%     {'pu11wh87','SeqDepPitchLMAN',{'bccbB', 'dccbB'}};

SylsToRemove_SingleDir= ...
    {'pu53wh88','SeqDepPitchLMAN', {'abbB'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
SylsToRemove_BiDir= ...
    {};
PARAMS.global.LMAN.SylsToRemove_SingleDir=SylsToRemove_SingleDir;
PARAMS.global.LMAN.SylsToRemove_BiDir=SylsToRemove_BiDir;



% 2) RUN VARIOUS ANALYSIS
% ==== BASELINE CORRELATIONS
close all;
[SeqDepPitch_AcrossBirds_LMAN, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMAN(SeqDepPitch_AcrossBirds_LMAN, PARAMS);





% ++++++++++++++++++++++++++++++++++++ SINGLE DIR
% ===== SINGLE DIR LEARNING [INLUDING AFP/AFP, ETC]
close all;
norm_by_targsyl=0; % normalize within each experiment
epochfield_input='days_consolid_early';
epochfield_input='final_extracted_window';
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl, epochfield_input);

% ====== SINGLE DIR - AFP BIAS LOCALIZED?
close all;
norm_by_targsyl=1; % normalize within each experiment
epochfield_input='days_consolid_early';
epochfield_input='final_extracted_window';
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANLocal(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl, epochfield_input);

% ===== SINGLE DIR - AFP BIAS, CLEAN PLOTS
close all;
norm_by_targsyl=0; % normalize within each experiment
% epochfield_input='days_consolid_early';
epochfield_input='final_extracted_window';
UseBaselineForCV=0; % then uses baseline data for CV reduction analysis
DispEachSylCVpval=0; % if 1, lists p vals
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANLocalClean(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl, epochfield_input, UseBaselineForCV, DispEachSylCVpval);


% ====== METRIC FOR SINGLE DIR LEARNING (METRIC OF SPECIFICITY)
close all;
norm_by_targsyl=0; % normalize within each experiment
% epochfield_input='days_consolid_early';
epochfield_input='final_extracted_window';
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANmetric(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl, epochfield_input);



% ===== SINGLE DIR - CONSOLIDATION DRIVES GENERALIZATION?
close all;
norm_by_targsyl=0; % normalize within each experiment
% epochfield_input='days_consolid_early';
epochfield_input='final_extracted_window';
similar_only=1; % 0=only diff; 1=similar only; 2 = all;

[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANconsolGen(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl, epochfield_input, similar_only);


% +++++++++++++++++++++++++++++++ 

% ++++++++++++++++++++++++++= PLOT FOR BIDIR LEARNING
close all;
use_final_extracted_windows=1; % if 0, then uses end of consolid and end of bidir (early if no). if 1, then uses consolid end + final extracted bidir window (above)
NormToTarg=0;
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds_LMAN, PARAMS, use_final_extracted_windows,NormToTarg);


% ==== BIDIR, SEPARATION METRIC
close all;
use_final_extracted_windows=1; % if 0, then uses end of consolid and end of bidir (early if no). if 1, then uses consolid end + final extracted bidir window (above)
NormToTarg=0;
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANmetric2(SeqDepPitch_AcrossBirds_LMAN, PARAMS, use_final_extracted_windows,NormToTarg);




% ++++++++++++++++++++++++++= TRAJECTORY OVER TIME  [SINGLE DIR AND BIDIR]
close all
% === MULTIDIR
OnlyConsolPeriod_bidir=1; %  still locked to start of bidir, but only plots days that fall in this window
BinRelConsolDay1=1; % if 0, then bins rel to bidir day 1
TakeAverageWithinExpt=1; % if 1, then average within each expt if multiple datapts in one time bin
NumDaysInBin=2;
GetAutoConsolWindows=1; % replace manual consol windows with auto?
% NOTE: ACTUAL WINDOWS ARE AUTO, THEY ARE CALCULATED ON EACH RUN, CODE
% OUTPUTS. THESE WINDOWS BELOW ARE OLD HAND CHOSEN.
PARAMS.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart={...
    'pu11wh87', 'SeqDepPitchLMAN',[4 8], ... %     'pu11wh87', 'SeqDepPitchLMAN2',[6 10], ... THIS IS DIFF SYLS
      'gr41gr90', 'SeqDepPitchLMAN2',[3 4], ... %      'pu11wh87', 'SeqDepPitchLMAN6',[3 7], ... this is diff syl
      'rd23gr89', 'SeqDepPitchLMAN',[5 8], ...
       'rd28pu64', 'SeqDepPitchLMAN',[3 8], ...
      'bk34bk68', 'SeqDepPitchLMAN3',[11 16], ...
      'wh4wh77', 'SeqDepPitchLMAN',[4 12], ...
      'wh25pk77', 'SeqDepPitchLMAN',[4 10], ...
      'bk34bk68', 'SeqDepPitchLMAN',[3 8]}
  
% PARAMS.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart={...
%     'pu11wh87', 'SeqDepPitchLMAN',[4 8], ... %     'pu11wh87', 'SeqDepPitchLMAN2',[6 10], ... THIS IS DIFF SYLS
%      'pu11wh87', 'SeqDepPitchLMAN6',[3 7], ...
%       'gr41gr90', 'SeqDepPitchLMAN2',[3 4], ...
%       'rd23gr89', 'SeqDepPitchLMAN',[4 8], ...
%        'rd28pu64', 'SeqDepPitchLMAN',[4 10], ...
%       'wh4wh77', 'SeqDepPitchLMAN',[4 12], ...
%       'wh25pk77', 'SeqDepPitchLMAN',[4 10], ...
%       'bk34bk68', 'SeqDepPitchLMAN',[4 10]} % only keeping expt where is first phase
   

[OUTPUT_multidir, DATSTRUCT_multidir] =lt_seq_dep_pitch_ACROSSBIRDS_LMANTimeCourse(SeqDepPitch_AcrossBirds_LMAN, PARAMS, OnlyConsolPeriod_bidir, BinRelConsolDay1, TakeAverageWithinExpt, NumDaysInBin, GetAutoConsolWindows);


% === SAME DIR
close all;
OnlyConsolPeriod_bidir=1; %  still locked to start of bidir, but only plots days that fall in this window
BinRelConsolDay1=1; % if 0, then bins rel to bidir day 1
TakeAverageWithinExpt=1; % if 1, then average within each expt if multiple datapts in one time bin
NumDaysInBin=2;
GetAutoConsolWindows=1;
% NOTE: ACTUAL WINDOWS ARE AUTO, THEY ARE CALCULATED ON EACH RUN, CODE
% OUTPUTS. THESE WINDOWS BELOW ARE OLD HAND CHOSEN.
PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={...
      'rd23gr89', 'SeqDepPitchLMAN2',[5 10], ...
      'bk34bk68', 'SeqDepPitchLMAN3',[2 11], ...
      'wh4wh77', 'SeqDepPitchLMAN',[6 13], ...
      'wh25pk77', 'SeqDepPitchLMAN',[4 9], ...
      'rd28pu64', 'SeqDepPitchLMAN2',[5 12]};
     
% PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={...
%       'rd23gr89', 'SeqDepPitchLMAN2',[5 10], ...
%       'bk34bk68', 'SeqDepPitchLMAN3',[2 11]}; % only experiments that are first phase
% 

[OUTPUT_samedir, DATSTRUCT_samedir]=lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_SamDir(SeqDepPitch_AcrossBirds_LMAN, PARAMS, OnlyConsolPeriod_bidir, BinRelConsolDay1, TakeAverageWithinExpt, NumDaysInBin, GetAutoConsolWindows);


% ++++++++++++++++++++ SINGLE TARG
close all;
OnlyConsolPeriod_bidir=1; %  still locked to start of bidir, but only plots days that fall in this window
BinRelConsolDay1=1; % if 0, then bins rel to bidir day 1
TakeAverageWithinExpt=1; % if 1, then average within each expt if multiple datapts in one time bin
NumDaysInBin=2;
GetAutoConsolWindows=0; % KEEP THIS 0 for single targ (see below) but 1 for samedir and diff dir
% PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={... % ONLY THOSE PAIRED WITH TWO-TARG SAME DIR
%     'rd23gr89',  'SeqDepPitchLMAN2', [10 14], ...
%     'rd28pu64',  'SeqDepPitchLMAN2', [7 11], ...
%     'bk34bk68',  'SeqDepPitchLMAN3', [3 8], ...
%     'wh4wh77',  'SeqDepPitchLMAN', [4 12]};
% 
% % EXPERIMENTS THAT HAVE DATA FOR ALL DAY BINS
% PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={...
%     'pu11wh87', 'SeqDepPitchLMAN',[13 28], ...
%     'pu11wh87', 'SeqDepPitchLMAN2',[5 15], ...
%     'rd23gr89',  'SeqDepPitchLMAN2', [10 14], ...
%     'rd28pu64',  'SeqDepPitchLMAN2', [7 11], ...
%     'bk34bk68',  'SeqDepPitchLMAN3', [3 8], ...
%     'wh4wh77',  'SeqDepPitchLMAN', [4 12]};


% DEFAULT: ALL EXPERIMENTS (NOTE: IS UP TO DATE AND CONSISTENT WITH
% MAINTAINED SHIFT RULES) [NOTE, not exactly same as computer decided,
% because some cases WN is still being updated. or because computer does
% not take into account start of bidir or samedir.
PARAMS.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart={...
    'pu11wh87', 'SeqDepPitchLMAN',[11 28], ... % avoid early as still adjusting
    'pu11wh87', 'SeqDepPitchLMAN2',[5 15], ... % one day skipped labeling, but is fine.
    'gr41gr90', 'SeqDepPitchLMAN2', [2 6], ...
    'rd23gr89',  'SeqDepPitchLMAN2', [7 14], ...
    'rd28pu64',  'SeqDepPitchLMAN2', [7 11], ... % gap in day
    'bk34bk68',  'SeqDepPitchLMAN3', [5 8], ...
    'bk34bk68',  'SeqDepPitchLMAN', [3 7], ...
    'wh4wh77',  'SeqDepPitchLMAN', [4 12], ...
    'wh25pk77',  'SeqDepPitchLMAN', [3 7]};

 
[OUTPUT_singleTarg, DATSTRUCT_singleTarg]=lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_SinTarg(SeqDepPitch_AcrossBirds_LMAN, PARAMS, OnlyConsolPeriod_bidir, BinRelConsolDay1, TakeAverageWithinExpt, NumDaysInBin, GetAutoConsolWindows);


% ===== STATS COMPARE TWO DIR TO SAME DIR
lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_stats; % hand entered values
lt_seq_dep_pitch_ACROSSBIRDS_LMANhitrate_1; % ad hoc, hand entered values 


% note: for below, make sure the same dir and multidir codes used the same
% day bin sizes (qwithin the code). also need to make sure they are using
% the same syllables
close all;
DayBinToUse=3;
lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_stats2(OUTPUT_multidir, OUTPUT_samedir, DATSTRUCT_multidir, DATSTRUCT_samedir,DayBinToUse); % using output structures


% =============================================== COMPILE ALL CONSOLIDATION
% (SINGLE TARG, TWO TARG BIDIR, TWO TARG SAME DIR

% note: for below, make sure the same dir and multidir codes used the same
% day bin sizes (qwithin the code). also need to make sure they are using
% the same syllables
close all;
DayBinToUse=3;
AdHocExpt={};
% --- if using MUSC baseline (default)
if (0)
    AdHocExpt={'wh25pk77', 'SeqDepPitchLMAN', 'multiDir2', ...
    {'ff_combined_PBS', 308.4 , 'ff_combined_MUSC', 191.3, 'ff_first_PBS', 215.4, 'ff_first_MUSC', 206.7, ...
    'ff_second_PBS', -93, 'ff_second_MUSC', 15.4, 'phase_day1_relWN', 29},...
    'bin 2 of 3-day binds'}; % 5-tuple cell arrays
% ---- if using PBS baseline for both musc and pbs data 
    AdHocExpt={'wh25pk77', 'SeqDepPitchLMAN', 'multiDir2', ...
    {'ff_combined_PBS', 308.45, 'ff_combined_MUSC', 230.3, 'ff_first_PBS', 215.35, 'ff_first_MUSC', 206.19, ...
    'ff_second_PBS', -93.1 , 'ff_second_MUSC', -24.1 , 'phase_day1_relWN', 29},...
    'bin 2 of 3-day binds'}; % 5-tuple cell arrays
end

% BELOW: used 1st learning epoch, wrong vals.
% % --- if using MUSC baseline (default)
% AdHocExpt={'wh25pk77', 'SeqDepPitchLMAN', 'multiDir2', ...
%     {'ff_combined_PBS', 242.5, 'ff_combined_MUSC', 180.8, 'ff_first_PBS', 217.5, 'ff_first_MUSC', 169.7, ...
%     'ff_second_PBS', -25.1, 'ff_second_MUSC', -10.93, 'phase_day1_relWN', 29},...
%     'bin 2 of 3-day binds'}; % 5-tuple cell arrays
% % ---- if using PBS baseline for both musc and pbs data 
% AdHocExpt={'wh25pk77', 'SeqDepPitchLMAN', 'multiDir2', ...
%     {'ff_combined_PBS', 199.7, 'ff_combined_MUSC', 186.9, 'ff_first_PBS', 215.35, 'ff_first_MUSC', 206.19, ...
%     'ff_second_PBS', 15.62 , 'ff_second_MUSC', 19.25 , 'phase_day1_relWN', 29},...
%     'bin 2 of 3-day binds'}; % 5-tuple cell arrays

TakeAvg=0;
ExpToAvg={};
ExpToAvg{1}={'wh25pk77', 'SeqDepPitchLMAN', 'multiDir', 'wh25pk77', 'SeqDepPitchLMAN',...
    'multiDir2'}; % {bird1, expt1, phase1, bird2, expt2, phase2};
ExpToAvg{1}={'bk34bk68', 'SeqDepPitchLMAN', 'multiDir', 'bk34bk68', 'SeqDepPitchLMAN3', ...
    'multiDir'};


lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_stats4(SeqDepPitch_AcrossBirds_LMAN, OUTPUT_multidir, ...
    OUTPUT_samedir, OUTPUT_singleTarg, DATSTRUCT_multidir, DATSTRUCT_samedir, ...
    DATSTRUCT_singleTarg,  DayBinToUse, AdHocExpt, TakeAvg, ExpToAvg); % using output structures




%% LMAN, learning state space day to day + SINGLE --> BIDIR EXPERIMENTS


lt_seq_dep_pitch_ACROSSBIRDS_LMAN_State(SeqDepPitch_AcrossBirds_LMAN, PARAMS);



%% --- [OLD!!! - does not combine multiple datapioints within bins, etc, use above] v2, get timecourse locked to arbitrary day, and bin. like above, but
% above was locked to bidir start
% -- [USEFUL FOR PLOTTING ENTIRE TRAJECTORY WITH LMAN DATA]
close all;
DayToLockAllTo='BidirDay1';
DayToLockAllTo='WNday1';
DayToLockAllTo='BaseDay1';
DayToLockAllTo='SingleDirConsolid'; % locks all to start of consolid, and plots only things within time window
PARAMS.LMANTimeCourse_v2.SingleDirConsolid= ...
    {'pu53wh88', 'SeqDepPitchLMAN', [19 21], ...
    'pu11wh87', 'SeqDepPitchLMAN', [15 28], ...
    'pu11wh87', 'SeqDepPitchLMAN2', [6 15], ...
    'pu11wh87', 'SeqDepPitchLMAN3', [8 15], ...
    'pu11wh87', 'SeqDepPitchLMAN6', [4 5], ...
    'gr41gr90', 'SeqDepPitchLMAN', [4 8], ...
    'gr41gr90', 'SeqDepPitchLMAN2', [3 6], ...
    'rd23gr89', 'SeqDepPitchLMAN', [4 5], ...
    'rd23gr89',  'SeqDepPitchLMAN2', [10 14], ...
    'rd28pu64',  'SeqDepPitchLMAN', [4 6], ...
    'rd28pu64',  'SeqDepPitchLMAN2', [7 11], ...
    'bk34bk68',  'SeqDepPitchLMAN', [5 7], ...
    'bk34bk68',  'SeqDepPitchLMAN3', [3 8], ...
    'wh4wh77',  'SeqDepPitchLMAN', [4 12], ...
    'wh25pk77',  'SeqDepPitchLMAN', [5 7]}; % triplets {bird, expt, [day1 last day] (days rel to WN on day 1)

% OnlyConsolPeriod_bidir=0;
lt_seq_dep_pitch_ACROSSBIRDS_LMANTimeCourse_v2(SeqDepPitch_AcrossBirds_LMAN, PARAMS, DayToLockAllTo);




% +++++++++++++++++++++++++++++++++++++++++++++++++

% --- Separation - divide MUSC days into 2 bins, early and late
% [MODIFY THIS SO I CAN INPUT THE RANGE OF DAY BINS - I.E. CONSOLIDATION
% BINS]
lt_seq_dep_pitch_ACROSSBIRDS_LMANseparation(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl)



%% +++++++++++++++++++++ SEQUENCE LEARNING - EFFECT OF MUSC
close all; 
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANsquence(SeqDepPitch_AcrossBirds_LMAN, PARAMS);

% -- using user entered transitions
close all; 

PARAMS.sequence_v2.birdname_exptname_Motifs= ...
    {'pu53wh88', 'SeqDepPitchLMAN', {'ab', 'ac'}, ...
    'pu11wh87', 'SeqDepPitchLMAN', {'ab','dc'}, ...
    'pu11wh87', 'SeqDepPitchLMAN2', {'ab','dc'}, ...
    'pu11wh87', 'SeqDepPitchLMAN3', {'ab','dc'}, ...
    'pu11wh87', 'SeqDepPitchLMAN6', {'ab', 'dc'}, ...
    'gr41gr90', 'SeqDepPitchLMAN', {'jbba', 'dbba', 'gbba'}, ...
    'gr41gr90', 'SeqDepPitchLMAN2', {'jbba', 'dbba', 'gbba'}, ...
    'rd23gr89', 'SeqDepPitchLMAN', {'db', 'ah', 'ak'}, ...
    'rd28pu64', 'SeqDepPitchLMAN', {'jjb', 'kjb', 'jab'}, ...
    'bk34bk68', 'SeqDepPitchLMAN', {'jjb', 'ljb'}};



[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANsquence_v2(SeqDepPitch_AcrossBirds_LMAN, PARAMS);






%% ++++++++++++++++++++++++++++ REPEATS
% ==== x
epoch_field='final_extracted_window';
lt_seq_dep_pitch_ACROSSBIRDS_LMANRepeats(SeqDepPitch_AcrossBirds_LMAN, PARAMS, epoch_field);



% +++++++++++++++++++++++ [DIRTY!!!!] PLOT STUFF RELATIVE TO MOTIF DISTANCE
close all;
plotLMAN_musc=1; % if 0, then plots LMAN experiments, but PBS data
% FIRST, change the name of the field so the MUSC data will be where normal
% data are
SeqDepPitch_AcrossBirds_TEMP=SeqDepPitch_AcrossBirds_LMAN;
if plotLMAN_musc==1;
    for i=1:length(SeqDepPitch_AcrossBirds_TEMP.birds);
        numexperiments=length(SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment);
        
        for ii=1:numexperiments;
            SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment{ii}.Syl_ID_Dimensions=SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN;
        end
    end
end
% Things like correlations, learning etc.
[SeqDepPitch_AcrossBirds_TEMP, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

% Rel Target correaltions/acoustic vs. motif pos
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_reltarg(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

% All pairs corre/acoustic vs. motif post
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_allpairs(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% SAVING
savedir='/bluejay4/lucas/across_birds/seq_dep_pitch';
currdir=pwd;

try
    cd(savedir);
catch err
    mkdir(savedir);
    cd(savedir);
end

save('SeqDepPitch_AcrossBirds', 'SeqDepPitch_AcrossBirds')
try
    save('SeqDepPitch_AcrossBirds_ONETARG', 'SeqDepPitch_AcrossBirds_ONETARG')
catch err
end

cd(currdir)

