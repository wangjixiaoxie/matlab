function [filename_save, all_days_various]=lt_all_days_various_calculations_FUNCTION_MOTIF(phrase, first_day, last_day, which_analy, within_analy_params, save_results, MOTIF_params)
% All arguments same as lt_all_days_various_calculations_FUNCTION except
% MOTIF_params: {syl_targ, motif_pre, motif_post, syl_replace};
% e.g. if you want to only analyze c in context of ab[c]de, have:
% syl_targ='c';
% motif_pre='ab';
% motif_post='de';
% syl_replace='Z' (i.e. you want to use Z as a placeholder for c in that context)
    

%% LT 7/7/14 - Perform any analysis on syllables in speicifc sequential contexts (e.g. b only after agj.).  Do this as follows:
% 1) go through all label files and change that motif to use a unique label
% (e.g. capital letter)
% 2) run analysis as normal
% 3) go back and change that letter back to normal.  

% ADV specifi

%% First, go through each day and change desired syl in motif to an unused label (only does if the replace label is different from the target)

curr_dir=pwd; % this should be bird dir;

% INPUTS - specific for doing replacement
syl_targ = MOTIF_params{1};
motif_pre= MOTIF_params{2};
motif_post= MOTIF_params{3};
syl_replace=MOTIF_params{4};

% AUTOMATIC STUFF
MotifOrig=[motif_pre syl_targ motif_post];
MotifNew=[motif_pre syl_replace motif_post];

if ~strcmp(syl_replace,syl_targ);
% CHANGE those labels (error-proof: will automatically go thru days, and give
% error and user control if it finds that syl already used that day)
lt_all_days_various_calculations_FUNCTION(phrase,first_day,last_day,{'changesyl'},{'old_motif',MotifOrig,...
    'new_motif',MotifNew,'syl_that_should_not_already_be_used',syl_replace},0);
end

%% Second, run analyses as usual, targetting the syl_replace

cd(curr_dir);
[filename_save, all_days_various]=lt_all_days_various_calculations_FUNCTION(phrase, first_day, last_day, which_analy, within_analy_params,save_results);

%% Third, change the syl label back to normal

if ~strcmp(syl_replace,syl_targ);
lt_all_days_various_calculations_FUNCTION(phrase,first_day,last_day,{'changesyl'},{'old_motif',MotifNew,...
    'new_motif',MotifOrig,'syl_that_should_not_already_be_used',''},0);
end

%% Fourth, save a note to know what the motif "replacement" is (i.e. this is like a code that shoudl always be linked to this motif:

log_note=['Actual: ' MotifOrig ' | Code: ' MotifNew ' | ADV structure: ' filename_save ' | Date: ' all_days_various.parameters.time_of_analysis];

try
cd([curr_dir '/all_days_various_calculations']);
catch err
    mkdir([curr_dir '/all_days_various_calculations']);
    cd([curr_dir '/all_days_various_calculations']);
end
fid=fopen('MOTIF_CODE_all_days_various_calculations_FUNCTION_MOTIF.txt','a');
fprintf(fid,'%s \n \n',log_note);
fclose(fid);
cd(curr_dir);

