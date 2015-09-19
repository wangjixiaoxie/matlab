function [filename_save, all_days_various]=lt_all_days_various_calculations_FUNCTION(phrase, first_day, last_day, which_analy, within_analy_params,save_results)
%% LT 7/4/14 - started converting to function. 7/7/ complete
% varargin{1}= 1 (save) or 0 (no save);
% varargin{2}= 1 (keep plots) or 0 (close plots);


% Instructions:
% 1) phrase= phase indexing the experiments? i.e. % will only go through folders with this in the name
% 2) first_day e.g. 06Jun2012
% 3) last_day e.g. 06Jun2012
% 4) list of analyses to perform: the possible analyses names are listed in
% the function. e.g. to run 3 specific analyses:
% which_analy={'dbseqf','rendamnt','allvocal'};
% 5) then you want to enter the parameters within each analysis.  use
% within_analy_params . every 2 entries designates one paramter.  odd entries are paramter names, even entries are the values
% e.g. % within_analy_params={'syllables',{'d'},'syllables_pre',{'f'},'syllables_post',{''},'btaf',''};
% 6) save (1 = yes, 0 = no)


%% LT started 11/1/13 - calculating entropy, amplitude, sung fraction, using db_seq_func_day_save code
% 12/2013 - because of pretime and posttime for calcualting entropy (rel
% to onset and offset).  change for each bird and syllable. (see line 25)

% this goes through all day folders and takes labeled notes and does
% calcualtiones (which you determine), and then you ahve the option of
% plotting all days for each calculation.

% taken from the end of lt_db_contour...v4 on 11/18.

% 1/4/14 - renamed v2 from "EDIT" version, use this one instead of v1 (i.e.
% no version).
% 1/14/14 - added HIT rate calc.

% clear all

%% input parameters
% curr_directory=pwd;
% cd ..;
% one_dir_up=pwd; cd(curr_directory);
%

% if strcmp(input('do you want to use previous parameters?','s'),'y')==1;
%     disp('go to all_days_various_calculations folder and load the parameters, then type return to return to script')
%     keyboard
%
%     % re-assign the variables which were stored in a structure, but now need
%     % to be in single-variable form (not structure).
%     phrase=all_days_various.parameters.phrase;
%     nameofbird=all_days_various.parameters.nameofbird;
%
%     days = all_days_various.parameters.days;
%
%     omit = 0;
%     syllables=all_days_various.parameters.syllables;
%
%     % ask what analysis you want to do.
%     seq_func_yes_or_no=all_days_various.parameters.questions.rendition_amount_yes_or_no;
%     rendition_amount_yes_or_no=all_days_various.parameters.questions.seq_func_yes_or_no;
%     plot_all_days_yes_or_no=all_days_various.parameters.questions.plot_all_days_yes_or_no;
%     get_all_vocalization_stats_yes_or_no=all_days_various.parameters.questions.get_all_vocalization_stats_yes_or_no;
%     mult_seq_yes_or_no=all_days_various.parameters.questions.mult_seq_yes_or_no;
%     findwnote2tw_yes_or_no=all_days_various.parameters.questions.findwnote2tw_yes_or_no;
%
%     if strcmp(findwnote2tw_yes_or_no,'y')==1;
%         findwnote_syl=all_days_various.findwnote2tw.findwnote_syl;
%     end
%
%
%     if strcmp(mult_seq_yes_or_no,'y')==1;
%         num_tp=all_days_various.mult_seq_analy.num_tp;
%         tpofinterest=all_days_various.mult_seq_analy.tpofinterest;
%         con_or_div=all_days_various.mult_seq_analy.con_or_div;
%     end
%
%     if seq_func_yes_or_no=='y'||rendition_amount_yes_or_no=='y';
%         syllables=all_days_various.parameters.syllables;
%     end
%
% %     if seq_func_yes_or_no=='y';
% %         entropy_pretime=all_days_various.seq_func_day_save.entropy_pretime;
% %         entropy_posttime=all_days_various.seq_func_day_save.entropy_posttime;
% %     end
%
%     if plot_all_days_yes_or_no=='y';
%         metrics=all_days_various.seq_func_day_save.metrics;
%     end
%
%
% else



%% FUNCTION ONLY - Initiate parameters


keep_plots=0; % 0 to close, 1 to keep


% initiate as "no":
seq_func_yes_or_no='n';
rendition_amount_yes_or_no='n';
get_all_vocalization_stats_yes_or_no='n';
mult_seq_yes_or_no='n';
triglabel_yes_or_no='n';
change_syl_yes_or_no='n';
get_labels_yes_or_no='n';
get_transitions_yes_or_no='n';
lt_db_check_yes_or_no='n';
lt_calc_day_pitch='n';
check_HTF='n';
seq_dep_pitch='n';

num_args=size(which_analy,2);
for i=1:num_args;
    switch which_analy{i};
        case 'dbseqf';
            seq_func_yes_or_no='y'; % 'analyze spectral features (e.g. entropy)  (y or n? ', 's');
        case 'rendamnt';
            rendition_amount_yes_or_no='y'; %input('calculate # of labeled songs and syl renditions for each day? ', 's');
        case 'allvocal';
            get_all_vocalization_stats_yes_or_no='y'; %input('perform lt_get_all_vocalization_stats? (i.e. get stats, such as duration, for each vocalization in batch (even unlabeled) ','s');
        case 'multseq';
            mult_seq_yes_or_no='y'; %input('perform mult_seq_analysis on transitions? ','s');
        case 'triglabel'
            triglabel_yes_or_no='y'; %input('perform triglabel to get WN hits, misses, for specific syls? ','s');
        case 'changesyl'
            change_syl_yes_or_no='y'; %input('change syl labels using db_change_syl_... ? ','s');
        case 'getlabels'
            get_labels_yes_or_no='y'; %input('perform lt_db_get_labels to get labels, onsets, offsets, time, for all syl? ','s');
        case 'ltgettrans'
            get_transitions_yes_or_no='y'; %input('perform lt_get_all_transition... on each day? ','s');
        case 'checktemplfreq';
            lt_db_check_yes_or_no='y'; %input('perform lt_db_check_..._BIRDTAF to get hit rates, frequencies, for pitch training? ','s');
            %         case 'optionalstuff' % add this if run the optional things.
        case 'calcdaypitch';
            lt_calc_day_pitch='y';
        case 'check_HTF'; % lt_check_hit_templ_freq
            check_HTF='y';
        case 'seq_dep_pitch';
            seq_dep_pitch='y'; % lt_compile_seq_dep_pitch_data
        case 'seq_dep_pitch_2'; % lt_seq_dep_pitch_DayRawDat
            seq_dep_pitch_2='y';
            
        otherwise
            disp(['the varargin ' which_analy{i} ' was not understood']); % give error saywing that that input does not match
    end
end

% GETTING WITHIN-ANALYSIS PARAMETERS:
try
    num_params=size(within_analy_params,2)/2;
catch err
    disp('mistake in entering arguments - within_analy_params should have even # of entries');
end

for i=1:num_params;
%     params.(within_analy_params{i*2-1})=within_analy_params{i*2};
    eval([within_analy_params{i*2-1} '= within_analy_params{i*2}']);
%     eval([within_analy_params{i*2-1} '= params.(within_analy_params{i*2-1})']) 
%     assignin('caller',within_analy_params{i*2-1},within_analy_params{i*2}); % assign param variables
%     eval([within_analy_params{i*2-1} '=' sprintf(within_analy_params{i*2})])
end


% CONVERTING from function parameters to parameters originally used in script.
days{1} = first_day;
days{2} = last_day;
days{3} = datenum(days{1});
days{4} = datenum(days{2});

%% OTHER PARAMETERS - DISCONTINUED here (function ver)
% phrase=input('what is your phase indexing the experiments? ','s');
curr_directory=pwd;
cd ..;
one_dir_up=pwd; cd(curr_directory);

[nameofbird]=lt_get_birdname_date_from_dir(0);
pwd;
omit = 0;

% WITHIN ANALYSIS PARAMS:  Some are just listed, and not "input" (for function ver)
if (0) % only use this for string.  in function will have to enter variables as arguments
    if strcmp(lt_db_check_yes_or_no,'y')==1;
        if strcmp(input('lt_db_check: use batch.catch.keep? ','s'),'y')==1;
            syl_check_batch='batch.catch.keep';
        else
            syl_check_batch=input('what batch to use? ','s');
        end
        syllables=input('what syllable was targeted (code only supports one)? (e.g. {"d"}) ');
        syllables_pre=input('what syllable comes before ? e.g. {"f"} or {''}) ');
        syllables_post=input('what syllable comes after? e.g. {"g"}) ');
        btaf=input('template used birdtaf? (y or n) ','s');
    end
    
    if strcmp(change_syl_yes_or_no,'y')==1;
        old_motif=input('old motif you want to change (e.g. ab). ', 's');
        new_motif=input('new motif to change it into (e.g. ac), ','s');
        syl_that_should_not_already_be_used=input('syl_that_should_not_already_be_used? (enter "" if don"t care)','s');
    end
    
    if strcmp(triglabel_yes_or_no,'y')==1;
        triglabel_birdtaf=input('triglabel: did template use birdtaf ? (yes: type y, otherwise type n) ','s');
        if triglabel_birdtaf=='y';
            trig_syl=input('enter the syllables to check for WN hit rate (e.g. abc )','s');
            trig_syl_pre=input('enter syls that come before (in same order; e.g. jkm [blank = none]) ','s');
            %             trig_syl_post=input('enter syls that come after (in same order) ','s');
        elseif triglabel_birdtaf=='n';
            trig_syl=input('enter the syllables to check for WN hit rate (e.g. abc )','s');
        end
    end
    
        
    if strcmp(mult_seq_yes_or_no,'y')==1;
        num_tp=input('How many transition points?  ');
        for i = 1:num_tp;
            %asks for the name of the transition points
            tpofinterest{i} = input(['What is the name of transition #' num2str(i) '?  '], 's');
        end
        con_or_div=input('con or div transitions? (i.e. con or div)','s');
        
    end
    
    
    if rendition_amount_yes_or_no=='y';
        syllables_RendAmnt = input('what are the syllables for rendition amount analyses)? (e.g. ab) ','s');
    end
    
    
    if seq_func_yes_or_no=='y';
        %         for i=1:length(syllables);
        %             entropy_pretime.(syllables(i))= input(['pretime for entropy calculation for syllable ' ...
        %                 syllables(i) '? (i.e. 0.2 to start 0.2s before onset. negative # to start after onset) ']);
        %             entropy_posttime.(syllables(i))=input(['posttime for entropy calculation for syllable ' ...
        %                 syllables(i) '? (i.e. 0.2 to end 0.2s after onset.)']);
        %         end
        %         plot_all_days_yes_or_no=input('do you want to plot results for all days after calculations finish? (y or n) ', 's'); % do you want to plot all days after this is done?
        syllables_seq_func = input('what are the syllables for seq_func? (e.g. ab) ','s');
        if plot_all_days_yes_or_no=='y';
            metrics=input('what metrics do you want to look at? (e.g. for multiple, {"entropy", "amplitude", "fraction", "timing"}. otherwise enter 1 to use those 4 ');
            if metrics==1;
                metrics={'entropy', 'amplitude', 'fraction', 'timing'};
            end
        end
    end
    
    if lt_calc_day_pitch=='y';
        disp('lt_calc_day_pitch: enter parameters');
        syl_targetDP=input('targets for day pitch (e.g. "ab" (perform for both a and b)','s');
        syl_preDP=input('syl_pre: e.g. {"c","d"} (preceding syllables. leave as {"",""} if not needed)');
        phraseDP=input('what is phrase (e.g. CPseq (indexing both data folders and the future save folder)','s');
        freq_rangeDP= input('freq_range: e.g. {[2000 2500],[1000 1500]} - one for each syllable - keep this as tight as possible');
        pc_time_windowDP=input('pc_time_window. e.g. {[200 400], [100 300]}');
        plot_resultDP=('want to plot each day? (1 for yes, 0 for no)');
        pc_windowDP=('pc window (e.g. 3000)');
    end
    
%     Parameters for check_HTF
%     Already set:
%     batchHTF='batch.catch.keep';
%     syl_refractHTF=0.2
%     You enter:     
%     sylHTF={'a','b',}
%     syl_preHTF=same
%     syl_postHTF=same
%     freq_rangeHTF{}
%     evtaf_verHTF
%     get_WN_hitsHTF
%     get_offline_matchHTF
%     get_FFHTF

% For seq_dep_pitch;
% batchSDP='batch.labeled.all';
% syllables_SDP={'b','c'};
% frequency_range_SDP ={[3000 3900],[2000 3000]}; % for findwnote
% pc_window_SDP =[0.06,0.1]; % size of syl data window (in sec); (how much data to get for each rend (relative to onset), in sec)
% pc_time_window ={[25 120],[45 220]}; %for pitch contour (time bins to avg).
% 
% SeqPreList_SDP={'ac','ab'}; % format: 1st elem of all three lists should combine
% SylTargList_SDP={'c','b'}; % these must already have raw data compiled above
% SeqPostList_SDP={'b',''};


end


%% FIRST, set up loop to go through every day with songs
%finds the folder names to use
song_folders = lt_db_list_song_folders(phrase, omit,...
    curr_directory);
fid = fopen([curr_directory ...
    '/song_folders_' phrase '.txt'],'r');
folders = textscan(fid,'%s');
fclose(fid);

% %picks out only the folders within the dates specified
% folders{1} = folders{1}(find(datenum(song_folders) ==
% days{3}):find(datenum(song_folders) == days{4})); LT added, but no
% necessary
%

%% SECOND, go through each day and calculate whatever you want

for j = days{3}:days{4};
    if keep_plots==0; 
        close all; % suppress continuous figure output.
    end
    for k=1:length(song_folders);
        if j == datenum(song_folders(k));
            try
                cd([curr_directory '/' folders{1}{k}])
            catch err
                display([ datestr(j, 'ddmmmyyyy') ' is missing'])
                cd(curr_directory)
                break
            end
            
            %-----------------------------------------------------
            % INSERT ANALYSIS STUFF FOR EACH DAY HERE
            
            % ONE, for each day analyze spectral features:
            if seq_func_yes_or_no=='y';
                for i=1:length(syllables_seq_func);
                    lt_db_seq_func_day_save_rd3gr35('batch.catch.keep','',nameofbird, datestr(j,'ddmmmyyyy'), 1, syllables_seq_func(i), phrase)
                    disp([ datestr(j, 'ddmmmyyyy') ' for ' syllables_seq_func(i) ' done']);
                end
            end
            
            % TWO, INSERT ANALYSIS HERE (PITCH WAVERINESS?)
            
            % THREE, STORE NUMBER OF RENDITIONS FOR EACH SYLLABLE
            if rendition_amount_yes_or_no=='y';
                for i=1:length(syllables_RendAmnt);
                    jj=(j-days{3})+1; % make 1st day index 1, not the datenum.
                    [song_amount(jj), syl_renditions(i,jj)]=lt_calc_syl_rendition_amount('batch.catch.keep',syllables_RendAmnt(i));
                end
            end
            
            % FOURTH, calculate the durations for each (even unlabeled) vocalization in the labeled songs.
            if strcmp(get_all_vocalization_stats_yes_or_no,'y')==1;
                jj=(j-days{3})+1; % make 1st day index 1, not the datenum.
                [durations_all_songs{jj}, durations_per_song{jj}]=lt_get_all_vocalization_stats('batch.catch.keep');
            end
            
            %             % FIFTH, perform findwnote2tw
            %             if strcmp(findwnote2tw_yes_or_no,'y')==1;
            %                 jj=(j-days{3})+1;
            %                 fvalsstr{jj}=findwnote2tw('all_cbin_not_mat',findwnote_syl,'',-0.016,[1000 3000],1024,0,'obs0');
            %             end
            %
            % SIXTH, calculate transition stats
            if strcmp(mult_seq_yes_or_no,'y')==1;
                try
                    lt_db_mult_seq_analy_FUNCTION(one_dir_up,nameofbird,datestr(j, 'ddmmmyyyy'),phrase,'batch.catch.keep',num_tp,tpofinterest,con_or_div);
                    disp([ datestr(j, 'ddmmmyyyy') ' for mult_seq_analy done'])
                    close all; % close figures
                catch err
                    continue
                end
            end
            
            % 7/23 - commented out since lt_check_hit_templ_freq_FUNCTION
            % is better
            % SEVENTH, calculates WN hit rate
%             if strcmp(triglabel_yes_or_no,'y')==1;
%                 for ii=1:length(trig_syl);
%                     jj=(j-days{3})+1;
%                     if strcmp(triglabel_birdtaf,'n')==1;
%                         try
%                             [vals{j}{ii}, trigs{j}{ii}]=triglabel('batch.catch.keep',trig_syl(ii),1,0,0,0);
%                             trig_vals.(trig_syl(ii)){jj}=vals{j}{ii};
%                             trig_trigs.(trig_syl(ii)){jj}=trigs{j}{ii};
%                         catch err
%                             disp([trig_syl(ii) 'had error for day ' datestr(j, 'ddmmmyyyy') '; skipped syl']);
%                             continue
%                         end
%                     elseif strcmp(triglabel_birdtaf,'y')==1;
%                         try
%                             [vals{j}{ii}, trigs{j}{ii}]=triglabel2('batch.catch.keep',trig_syl(ii),trig_syl_pre(ii),...
%                                 '',1,0,0,0);
%                             trig_vals.(trig_syl(ii)){jj}=vals{j}{ii};
%                             trig_trigs.(trig_syl(ii)){jj}=trigs{j}{ii};
%                         catch err
%                             disp([trig_syl(ii) 'had error for day ' datestr(j, 'ddmmmyyyy') '; skipped syl']);
%                             continue
%                         end
%                     end
%                 end
%             end
            
            % EIGHTH, changes syl labels
            if strcmp(change_syl_yes_or_no,'y')==1;
                % Check whether that replacement syl used before
                syl_list_used=lt_get_syls_used_as_labels('batch.catch.keep');
                if sum(cell2mat(strfind(syl_list_used,syl_that_should_not_already_be_used)))>0; % ie.. replacement syl has been used today
                    disp(['problem: ' syl_that_should_not_already_be_used ' has been used today'])
                    if strcmp(input('do you still want to continue? (y continues, n terminates)','s'),'y');
                        %                         try
                        db_change_syllable_in_batchfile('batch.catch.keep',old_motif,new_motif);
                        disp(['changing syl: ' datestr(j, 'ddmmmyyyy') ' done.'])
                        %                         catch err
                        %                             disp(['changing syl: error for day ' datestr(j, 'ddmmmyyyy') '; skipped day']);
                        %                             continue
                        %                         end
                    else
                        dafasdfv; % use this to initiate error to halt program
                    end
                else
                    db_change_syllable_in_batchfile('batch.catch.keep',old_motif,new_motif);
                    disp(['changing syl: ' datestr(j, 'ddmmmyyyy') ' done.'])
                end
            end
            
            % NINTH, collect labels, onsets, offsets, and absolute times,
            % for all syllables.
            if strcmp(get_labels_yes_or_no,'y')==1;
                jj=(j-days{3})+1;
                disp([datestr(j,'ddmmmyyyy') '. lt_db_get_labels: gethering labels, onsets, etc, for all syl']);
                [syl{jj} time_syl{jj} filenames{jj} onsets_syl{jj} offsets_syl{jj}]=lt_db_get_labels('batch.catch.keep','n');
            end
            
            % TENTH, perform lt_get_all_transition...
            if strcmp(get_transitions_yes_or_no,'y')==1;
                jj=(j-days{3})+1;
                disp([datestr(j,'ddmmmyyyy') '. performing lt_get_all_transition...']);
                lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep',nameofbird,phrase,'',datestr(j,'ddmmmyyyy'));
                close all;
            end
            
            % ELEVENTH, perform % lt_db_check_template_timing_and_hit_rate_BIRDTAF
            % use this to get frequency at trigger timepoints (online trig)
            % using fft.
            
            if strcmp(lt_db_check_yes_or_no,'y')==1;
                jj=(j-days{3})+1;
                disp([datestr(j,'ddmmmyyyy') '. performing lt_db_check_templ_freq...']);
                all_days_various.lt_db_check_templ_freq.data{jj}=...
                    lt_db_check_template_timing_and_hit_rate_BIRDTAF_FUNCTION(syl_check_batch,1,syllables,...
                    syllables_pre,syllables_post,btaf);
                close all
            end
            
            % TWELFTH - lt_calc_day_pitch_v2 - this calcs pitch using picth
            % contour, and saves a stcuture for each day and also
            % incorporates into lt_all_days_various.
            
            if strcmp(lt_calc_day_pitch,'y')==1;
                jj=(j-days{3})+1;
                % check to make sure syl exists today
                syls_usedDP=lt_get_syls_used_as_labels('batch.catch.keep');
                if length(intersect(strjoin(syls_usedDP),syl_targetDP))>0;
                    disp([datestr(j,'ddmmmyyyy') '. performing lt_calc_day_pitch...']);
                    if exist(batchDP);
                    all_days_various.lt_calc_day_pitch{jj}=...
                        lt_calc_day_pitch_v2_FUNCTION(syl_targetDP, syl_preDP, phraseDP, freq_rangeDP, pc_time_windowDP, plot_resultDP,pc_windowDP,batchDP); % NOTE: automatically uses batch.catch.keep. to do otherwise, add batch as last argument.
                    else % batch chosen automatically, see code within.
                    all_days_various.lt_calc_day_pitch{jj}=...
                        lt_calc_day_pitch_v2_FUNCTION(syl_targetDP, syl_preDP, phraseDP, freq_rangeDP, pc_time_windowDP, plot_resultDP,pc_windowDP); % NOTE: automatically uses batch.catch.keep. to do otherwise, add batch as last argument.
                    end
                    else
                    disp(['calc_day_pitch: skipped ' datestr(j,'ddmmmyyyy') ' for syl ' syl_targetDP '. Did not occur today']);
                end
            end
            
            % THIRTEENTH, lt_check_hit_templ_freq - for each day gets hits,
            % uses template info to get matches, and get freq at
            % (triggered+labeled). 
            
            if strcmp(check_HTF,'y')==1;
                jj=(j-days{3})+1;
                disp([datestr(j,'ddmmmyyyy') '. performing lt_check_hit_templ_freq...']);
                batchHTF='batch.catch.keep';
                syl_refractHTF=0.2;
                for ii=1:length(sylHTF);
                    syl_labelHTF=[syl_preHTF{ii} sylHTF{ii} syl_postHTF{ii}];
                    if get_offline_matchHTF==1;
                        all_days_various.lt_check_hit_templ_freq.(syl_labelHTF){jj}=...
                            lt_check_hit_templ_freq_FUNCTION(batchHTF, sylHTF{ii}, syl_preHTF{ii}, syl_postHTF{ii}, syl_refractHTF,...
                            freq_rangeHTF{ii}, evtaf_verHTF,get_WN_hitsHTF,get_offline_matchHTF,get_FFHTF,template_name,cntrng_values,col_logic);
                    else
                        all_days_various.lt_check_hit_templ_freq.(syl_labelHTF){jj}=...
                            lt_check_hit_templ_freq_FUNCTION(batchHTF, sylHTF{ii}, syl_preHTF{ii}, syl_postHTF{ii}, syl_refractHTF,...
                            freq_rangeHTF{ii}, evtaf_verHTF,get_WN_hitsHTF,get_offline_matchHTF,get_FFHTF);
                    end
                end
            end
            
            
            % FOURTEENTH, lt_compile_seq_dep_pitch_data and
            % lt_compile_seq_dep_pitch_data_SEQFILTER; compiles pitch info
            % and filters sequence of the syl
            if strcmp(seq_dep_pitch,'y')==1;
                jj=(j-days{3})+1;
                disp([datestr(j,'ddmmmyyyy') '. performing lt_compile_seq_dep_pitch_data and lt_compile_seq_dep_pitch_data_SEQFILTER...']);
                plotON=0;
                clear compiled_seqdep_pitch;
                
                compiled_seqdep_pitch=lt_compile_seq_dep_pitch_data(batchSDP, syllables_SDP, frequency_range_SDP, pc_window_SDP, ...
                    pc_time_window, plotON); % compile data for single syls (e.g. a, b, c)
                
                compiled_seqdep_pitch=lt_compile_seq_dep_pitch_data_SEQFILTER(compiled_seqdep_pitch,...
                    SeqPreList_SDP,SylTargList_SDP,SeqPostList_SDP); % take that combiled data, and filter out sequences (e.g. b after a)
                
                all_days_various.compiled_seqdep_pitch{jj}=compiled_seqdep_pitch;
                
            end
  
            
            % FIFTEENTH, lt_seq_dep_pitch_DayRawDat (replaces
            % lt_compile_seq_dep_pitch_data above)
            
            if strcmp(seq_dep_pitch_2,'y')==1;
                
                % only continue if the batch desired actually contains
                % songs
                fid_temp=fopen(ParamsSDP.DayRawDat.batch,'r');
                line=fgetl(fid_temp);
                if line~=-1; % then this has string, which should be a song.
                    jj=(j-days{3})+1;
                    disp([datestr(j,'ddmmmyyyy') '. performing lt_seq_dep_pitch_DayRawDat']);
                    
                    [~, ~]=lt_seq_dep_pitch_DayRawDat(ParamsSDP,plotON_SDP,saveON_SDP,phrase);
                else
                    disp([datestr(j,'ddmmmyyyy') '. no songs in batch, not performing lt_seq_dep_pitch_DayRawDat']);
                end
                
            end
            
            
            
                % --------------------------------------------------------
                % OPTIONAL: ADD THINGS HERE THAT WILL NOT SAVE, BUT RANDOM
            % THINGS
            
            % FOR rd66gr93, context seq, going through folders and getting
            % transition data. 2 kinds of folders indexed by phrase
            if (0)
                folder_phrase='ContextSeq'; % marking the save folder
                [nameofbird, bluejaynum, date, ~]=lt_get_birdname_date_from_dir(1);
                datestring=date{2};
                
                if strcmp(phrase, 'ContextSeq_ContextA')==1;
                    lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep.early',nameofbird,folder_phrase,'contextA_early',datestring);
                    try
                        lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep.late',nameofbird,folder_phrase,'contextA_late',datestring);
                    catch err
                        disp(['late batch error for day ' datestr(j,'ddmmmyyyy') ])
                    end
                elseif strcmp(phrase, 'ContextSeq_ContextB')==1;
                    lt_get_all_transition_probabilities_FUNCTION('batch.catch.keep',nameofbird,folder_phrase,'contextB',datestring);
                end
            end
            
            % AFTER DONE ITERATE ONE STEP IN LOOP-------------------------
            cd(curr_directory)
            break
        else
            cd(curr_directory)
            continue
        end
    end
end

cd(curr_directory)
fprintf('done!\n')

%% saving variables

% FIRST, assign variables fields in one structure
all_days_various.parameters.phrase=phrase;
all_days_various.parameters.curr_directory=curr_directory;
all_days_various.parameters.nameofbird=nameofbird;
all_days_various.parameters.days=days;
all_days_various.parameters.time_of_analysis=lt_get_timestamp;

all_days_various.parameters.questions.seq_func_yes_or_no=seq_func_yes_or_no;
all_days_various.parameters.questions.rendition_amount_yes_or_no=rendition_amount_yes_or_no;
all_days_various.parameters.questions.get_all_vocalization_stats_yes_or_no=get_all_vocalization_stats_yes_or_no;
all_days_various.parameters.questions.mult_seq_yes_or_no=mult_seq_yes_or_no;
% all_days_various.parameters.questions.findwnote2tw_yes_or_no=findwnote2tw_yes_or_no;
all_days_various.parameters.questions.triglabel_yes_or_no=triglabel_yes_or_no;
all_days_various.parameters.questions.change_syl_yes_or_no=change_syl_yes_or_no;
all_days_various.parameters.questions.lt_db_get_labels=get_labels_yes_or_no;
all_days_various.parameters.questions.get_transitions_yes_or_no = get_transitions_yes_or_no;
all_days_various.parameters.questions.lt_db_check_yes_or_no = lt_db_check_yes_or_no;
all_days_various.parameters.questions.lt_calc_day_pitch = lt_calc_day_pitch;
all_days_various.parameters.questions.lt_check_hit_templ_freq = check_HTF;
all_days_various.parameters.questions.seq_dep_pitch = seq_dep_pitch;




things_analyzed=[];
things_analyzed_deep=['|']; % more information to log in text file

if strcmp(lt_db_check_yes_or_no,'y')==1;
    all_days_various.lt_db_check_templ_freq.parameters.syl_check_batch=syl_check_batch;
    all_days_various.lt_db_check_templ_freq.parameters.syllables=syllables;
    all_days_various.lt_db_check_templ_freq.parameters.syllables_pre=syllables_pre;
    all_days_various.lt_db_check_templ_freq.parameters.syllables_post=syllables_post;
    all_days_various.lt_db_check_templ_freq.parameters.btaf=btaf;
    things_analyzed=[things_analyzed '_LtDbCheckTemplFreq'];
    things_analyzed_deep=[things_analyzed_deep '_LtDbCheckTemplFreq_' syllables_pre syllables syllables_post '_|'];
end

if strcmp(get_transitions_yes_or_no,'y')==1;
    things_analyzed=[things_analyzed '_GetTransitions'];
    things_analyzed_deep=[things_analyzed_deep '_GetTransitions_|'];
end

if strcmp(get_labels_yes_or_no,'y')==1;
    all_days_various.lt_db_get_labels.labels=syl;
    all_days_various.lt_db_get_labels.abs_time_syl=time_syl;
    all_days_various.lt_db_get_labels.filenames=filenames;
    all_days_various.lt_db_get_labels.onsets=onsets_syl;
    all_days_various.lt_db_get_labels.offsets=offsets_syl;
    things_analyzed=[things_analyzed '_GetLabels'];
    things_analyzed_deep=[things_analyzed_deep '_GetLabels_|'];
end


if strcmp(change_syl_yes_or_no,'y')==1;
    all_days_various.change_syl.old_motif=old_motif;
    all_days_various.change_syl.new_motif=new_motif;
    all_days_various.change_syl.syl_that_should_not_already_be_used=syl_that_should_not_already_be_used;
    things_analyzed=[things_analyzed '_ChangedSyl'];
    things_analyzed_deep=[things_analyzed_deep '_ChangedSyl_' old_motif 'to' new_motif '_|'];
    
end


if strcmp(triglabel_yes_or_no,'y')==1;
    all_days_various.trig_label.vals=trig_vals;
    all_days_various.trig_label.trigs=trig_trigs;
    all_days_various.trig_label.trig_syl=trig_syl;
    try
        all_days_various.trig_label.trig_syl_pre=trig_syl_pre;
        all_days_various.trig_label.trig_syl_post=trig_syl_post;
    catch err
    end
    things_analyzed=[things_analyzed '_TrigLabel'];
    things_analyzed_deep=[things_analyzed_deep '_TrigLabel_' trig_syl '_|'];
end


% if strcmp(findwnote2tw_yes_or_no,'y')==1;
%     all_days_various.findwnote2tw.findwnote_syl=findwnote_syl;
%     all_days_various.findwnote2tw.fvalsstr=fvalsstr;
%         things_analyzed=[things_analyzed '_findwnote2tw'];
%     things_analyzed_deep=[things_analyzed_deep '_findwnote2tw_' findwnote_syl '_|'];
%
% end

if rendition_amount_yes_or_no=='y'
    all_days_various.syl_rendition_amount.song_amount=song_amount;
    all_days_various.syl_rendition_amount.syl_renditions=syl_renditions;
    all_days_various.syl_rendition_amount.syllables=syllables_RendAmnt;
    things_analyzed=[things_analyzed '_RenditionAmnt'];
    things_analyzed_deep=[things_analyzed_deep '_RenditionAmnt_' syllables_RendAmnt '_|'];
    
end

if seq_func_yes_or_no=='y'
    all_days_various.seq_func_day_save.syllables=syllables_seq_func;
    %     all_days_various.seq_func_day_save.metrics=metrics;
%     all_days_various.parameters.questions.plot_all_days_yes_or_no=plot_all_days_yes_or_no;
    
    things_analyzed=[things_analyzed '_SeqFunc'];
    things_analyzed_deep=[things_analyzed_deep '_SeqFunc_' syllables_seq_func '_|'];
    %     for i=1:length(metrics);
    %         things_analyzed=[things_analyzed metrics{i}];
    %         things_analyzed_deep=[things_analyzed_deep metrics{i}];
    %     end
    %     things_analyzed_deep=[things_analyzed_deep '_' syllables_seq_func];
end


if strcmp(get_all_vocalization_stats_yes_or_no,'y')==1;
    all_days_various.get_all_vocalization_stats.durations_all_songs=durations_all_songs;
    all_days_various.get_all_vocalization_stats.durations_per_song=durations_per_song;
    things_analyzed=[things_analyzed '_AllVocalStats'];
    things_analyzed_deep=[things_analyzed_deep '_AllVocalizationStats_|'];
    
end

if strcmp(mult_seq_yes_or_no,'y')==1;
    all_days_various.mult_seq_analy.num_tp=num_tp;
    all_days_various.mult_seq_analy.tpofinterest=tpofinterest;
    all_days_various.mult_seq_analy.con_or_div=con_or_div;
    things_analyzed=[things_analyzed '_MultSeq'];
    things_analyzed_deep=[things_analyzed_deep '_MultSeq_'];
    for i=1:length(tpofinterest);
        if i==length(tpofinterest);
            things_analyzed_deep=[things_analyzed_deep tpofinterest{i} '_|'];
        else
            things_analyzed_deep=[things_analyzed_deep tpofinterest{i} '_'];
        end
    end
end

if strcmp(lt_calc_day_pitch,'y')==1;
    things_analyzed=[things_analyzed '_CalcDayPitch'];
    things_analyzed_deep=[things_analyzed_deep '_CalcDayPitch_'];
    for i=1:length(syl_targetDP);
        if i==length(syl_targetDP);
            things_analyzed_deep=[things_analyzed_deep syl_targetDP '_|'];
        else
            things_analyzed_deep=[things_analyzed_deep syl_targetDP '_'];
        end
    end
end

if strcmp(check_HTF,'y')==1;
    things_analyzed=[things_analyzed '_CheckHitTemplFreq'];
    things_analyzed_deep=[things_analyzed_deep '_CheckHitTemplFreq'];
    all_days_various.lt_check_hit_templ_freq.parameters.sylHTF=sylHTF;
    all_days_various.lt_check_hit_templ_freq.parameters.syl_preHTF=syl_preHTF;
    all_days_various.lt_check_hit_templ_freq.parameters.syl_postHTF=syl_postHTF;
    all_days_various.lt_check_hit_templ_freq.parameters.evtaf_verHTF=evtaf_verHTF;  
end

if strcmp(seq_dep_pitch,'y')==1;
    things_analyzed=[things_analyzed '_seqdeppitch'];
    things_analyzed_deep=[things_analyzed_deep '_seqdeppitch'];
end



all_days_various.parameters.things_analyzed=things_analyzed;
all_days_various.parameters.things_analyzed_deep=things_analyzed_deep;


%% SAVING

if save_results==1;
try cd('all_days_various_calculations')
catch error
    mkdir('all_days_various_calculations')
    cd('all_days_various_calculations')
end

% save(['all_days_various_' all_days_various.parameters.time_of_analysis '_' things_analyzed],'all_days_various')
filename_save=['all_days_various_' all_days_various.parameters.time_of_analysis];
save(filename_save,'all_days_various','-v7.3');

% MAKES a note in a text file in the folder of 1) analyses compiled, 2) syl/motifs analyzed, 3) dates analyzed, 4) date of analysis (structure identifier)
log_note=[all_days_various.parameters.time_of_analysis '__|' things_analyzed_deep '|__' all_days_various.parameters.days{1} '_to_' all_days_various.parameters.days{2}];
fid=fopen('log.txt','a');
fprintf(fid,'%s \n \n',log_note);
fclose(fid);

cd(curr_directory)

else filename_save='did_not_save';
end
%% Below is to compile all days' data and plot. DO NOT USE. USE lt_all_days_all_analysis instead.
    if seq_func_yes_or_no=='y'
        if plot_all_days_yes_or_no=='y';
            lt_plot_all_days_spectral_features(metrics,all_days_various.parameters.curr_directory,all_days_various.parameters.phrase,...
                all_days_various.parameters.syllables)
        end
    end
end
