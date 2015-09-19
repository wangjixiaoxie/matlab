%% LT 6/3 - took from gr66 analysis file and made into its own program

% LT 6/2/14 - modifying this for context experiment 2, whcih will be done on bluejay, not gullserver

% NOTE: THIS is the latest version (5/26/14) that is being edited on gullserver.

% LT - 4/18/14 - Context dependent learning

%% NOTE TO SELF
% 1) make the data acquisition iterate over all contexts.
% 2) make the plots bars or lines
% 3) also plot total length of song (not just renditions of a) --> a halting effect?
%

%% Goals:
% 1) plot in bins (# rendition) transition probabilties for targets
% 2) plot rate of song halting
% 3) plot sining rate
% 4) put all on same graph that shows across day sand context switches


%% FIRST, going into each folder and getting transition data - do this after labeling and before compiling data into structure below
% run the below for each day folder (run either C1 or C2)
% if C1, then first separate b.c.k into b.c.k.early and b.c.k.late. Then
% run functions.  For C2, just use b.c.k. and run function.
clear all;close all

% INPUTS
folder_phrase='ContextSeq'; % marking the save folder
[nameofbird, bluejaynum, date, ~]=lt_get_birdname_date_from_dir(1);
datestring=date{2};

% RUN DIFFERENT DEPENDING ON CONTEXT

% BASELINE - no switching yet
% context A
lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep',nameofbird,folder_phrase,'contextA_early',datestring); % just call it contextA_early instead of contextA (for ease of code)

% context A
edit batch.catch.keep
lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep.early',nameofbird,folder_phrase,'contextA_early',datestring);
lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep.late',nameofbird,folder_phrase,'contextA_late',datestring);

%  context B
lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep',nameofbird,folder_phrase,'contextB',datestring);



%% SECOND, COMPILE DATA for all contexts for a single day into one data structure.
% INSTRUCTIONS: 1) first run lt_get_all_transition_probabilities. 2) Then enter
% your parameters below and run to compile transition info over multiple
% contexts into one structure. 3) Run this in the bird folder!

clear all; close all
% INPUTS
first_day='31May2014';
last_day='09Jun2014';
folder_phrase='ContextSeq'; %this marks the all_transitions folder (e.g. all_days_transition_matrix_ContextSeq)
motifs={'aa','ab','aj','an'}; % BEST TO TAKE ALL TRANSITIONS AT A BRANCH POINT
% motifs={'aa','ab','aj','an'}; % BEST TO TAKE ALL TRANSITIONS AT A BRANCH POINT
context_label_set={'contextA_early','contextB','contextA_late'}; % enter the contexts that are stored in that structure (exactly as they are written in the saved structures (e.g. rd66gr93_01Jun2014_contextA_early.mat)
curr_dir=pwd;
save_phrase='div_a';

lt_all_days_ContextSeq_COMPILE(first_day,last_day,folder_phrase,motifs,context_label_set,save_phrase)

cd(curr_dir);


% It will ask you for context switch times.  Keep track of them below:
% -------------------------------------------------------------
% LIST OF CONTEXT SWITCH TIMES
% rd66 trial1:
% {'31May2014-1205','31May2014-1658'}
% {'01Jun2014-1305','01Jun2014-1730','01Jun2014-2030'}
% {'02Jun2014-1210','02Jun2014-1707','02Jun2014-2030'}
% {'03Jun2014-1225','03Jun2014-1723','03Jun2014-2030'}
% {'04Jun2014-1205','04Jun2014-1213','04Jun2014-1703','04Jun2014-1757','04Jun2014-2035'}
% {'05Jun2014-1133','05Jun2014-1135','05Jun2014-1208','05Jun2014-1716'}
% {'06Jun2014-1201','06Jun2014-1714','06Jun2014-2019'}
% 
% {'07Jun2014-1218','07Jun2014-1724','07Jun2014-2032'}
% {'08Jun2014-1238','08Jun2014-1722','08Jun2014-2005'}
% {'09Jun2014-1317','09Jun2014-1712','09Jun2014-2100'}

%% PLOT
% Uses the same parameters in "SECOND" above.  
% GO INTO PROGRAM AND CHANGE PARAMETERS

% Uses the same parameters in "SECOND" above.  

% PLOT PARAMETERS
BaselineDays=1:4;
WNDays=5:10;
ProbeDays=[];

% PLOT PARAMETERS
% plotting across days
trn_to_plot='ab'; % for summary cross day data
num_renditions=6; % EDGES: how many songs to take at the edges
% plotting individual days
motifs_to_plot={trn_to_plot}; % for plotting each individual day
bin_size=3; % to change this you have to go back to making transitions, or incorporate that part here to make afresh

% ADDITIONAL INPUTS FOR PLOTTING
[birdname, bluejay_num]=lt_get_birdname_date_from_dir(0);
bird_dir=['/bluejay' num2str(bluejay_num) '/lucas/birds/' birdname];
phrase='div_a';

% RUN
lt_all_days_ContextSeq_PLOT



