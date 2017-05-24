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
