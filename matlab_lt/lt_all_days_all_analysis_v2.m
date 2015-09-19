%% LT 5/6/14 - MODIFIED in large way - the day argument number now goes on the end.  Modified everything that this code links out to and also the _PLOT program.
% Kept the old version as lt_all_days_all_analysis, and all related
% functions are not version "2"

%% USE THIS

%% LT 12/6/13 - use this to make one overarching structure for all days for all analyses for this bird.
% structure will be: day --> all songs --> different analyses --> all
% syllalbles;  and day--> songs separate --> analyses.
% to compile results from lt_db_seq_func_day_save_rd3gr35
% analysis of entropy, timing, amplitude, fraction (of syllables)
% and pitch data db_contour_and_FF_analysis_over_time_v3.


% THIS IS MODIFIED FROM lt_plot_all_days_spectral_features, mainly to store
% into structure instead of plotting. Start off in the main bird directory.

%% entering variables


% load_previous_data_y_or_n=input('load previous saved data structure and add new stuff to it? (y or n) ','s');
% if strcmp(load_previous_data_y_or_n, 'y')==1;
%     disp('do this manually: go to all_days_all_analysis folder and load the structure you want, then type "return"')
%     keyboard
%
%     % take saved variables and reassign into the simpler names used in this
%     % program.
%
%
%
% else
clear all;
[name_of_bird bluejay_num]=lt_get_birdname_date_from_dir(0);
curr_dir=pwd;

% asking which data you would like to compile
seq_func_day_question=input('compile syllable structure data from lt_db_seq_func_day_save_rd3gr35? ','s');
pitch_question=input('compile pitch data from db_contour_and_FF_analysis_over_time_v3 ','s');
syl_rendition_question=input('compile data on syl and song rendition amounts? (i.e. from all_days_various_calculations?) ','s');
get_all_vocalization_stats_question=input('compile data on stats of all vocalizations in labeled songs? (not just labeled) ','s');
gap_duration_question=input('compile gap_duration stats? ','s');
WN_trig_question=input('compile WN trig rates?' ,'s');
transitions_lt_question=input('compile transitions (lt function)? ','s');

% LOADING all_days_various_structures for those analyses that need them. --
if strcmp(WN_trig_question,'y')==1;
    cd('all_days_various_calculations'); ls
    open('log.txt')
    date_and_time_of_various_struct_WNtrig=input('enter date and time of struct you want to load for WN trig analysis (e.g. 06Jan2014_1737) ','s');
    cd(curr_dir);
end

%COMMENTED OUT.  See lt_all_days_various for reason. use data not from
%findwnote, but from lt_db_get_labels.
% if strcmp(gap_duration_question,'y')==1;
%     transitions_gd=input('what transitions to analyze for gap durations? (e.g. {"ab" "cd"})');
%     cd('all_days_various_calculations'); ls
%     open('log.txt')
%     date_and_time_of_various_struct_GapDur=input('enter date and time of struct you want to load for gap duration analysis (i.e. find2note2tw) (e.g. 06Jan2014_1737 ','s');
%     cd(curr_dir);
% end

if strcmp(transitions_lt_question,'y')==1;
transitions_phrase=input('phrase for all_days_transition_matrix? ','s');
end


if strcmp(gap_duration_question,'y')==1;
    transitions_gd=input('what transitions to analyze for gap durations? (e.g. {"ab" "cd"})');
    cd('all_days_various_calculations'); ls
    open('log.txt')
    date_and_time_of_various_struct_GapDur=input('enter date and time of struct you want to load for gap duration analysis (i.e. lt_db_get_labels) (e.g. 06Jan2014_1737 ','s');
    cd(curr_dir);
end

if strcmp(get_all_vocalization_stats_question,'y')==1;
    cd('all_days_various_calculations'); ls
    open('log.txt')
    date_and_time_of_various_struct_AllVocalStats=input('enter date and time of struct you want to load for AllVocalStats analysis (e.g. 06Jan2014_1737) ','s');
    cd(curr_dir);
end

if strcmp(syl_rendition_question,'y')==1;
    cd('all_days_various_calculations'); ls
    open('log.txt')
    date_and_time_of_various_struct_SylRend=input('enter date and time of struct you want to load for SylRendition analysis (e.g. 06Jan2014_1737) ','s');
    cd(curr_dir);
end

% if strcmp(input('automatically find all syllables that were used today? (y or n) ','s'),'y')==1;
%     % ENTER FUNCTION TO PULL OUT SYLLABLES HERE
%     were_syllables_chosen_automatically=1;
% else

% INPUTTING other variables ---------------------
if strcmp(seq_func_day_question,'y')==1;
syllables_seq_func = input('compile data for which syllables? (for seq_func)','s');
% were_syllables_chosen_automatically=0;
end

if strcmp(pitch_question,'y')==1;
syllables_pitch = input('compile data for which syllables? (for pitch analyses)','s');
% were_syllables_chosen_automatically=0;
end


if strcmp(seq_func_day_question,'y')==1;
    if strcmp(input('for syllable structure analysis, compile data for all 4 metrics (y or n)? (i.e. entropy, amplitude, timing, fraction) ','s'),'y')==1;
        metrics={'entropy','amplitude','timing','fraction'};
    else
        metrics=input('enter the metrics to look at (e.g. {"entropy","amplitude"} ');
    end
end
% 
% if strcmp(syl_rendition_question,'y')==1;
%     disp('will use all_days_various_calculations previously calculated syl. rendition data. MANUALLY go to folder and load latest structure, then type "return"');
%     keyboard
%     cd(curr_dir);
% end

first_day = input('first day to analyze? ','s');
last_day=input('last day to analyze? ','s');
phrase=input('phrase identifying this experiment? ','s');
things_analyzed=cell(0); % this will be a high level field that you can glance at quickly to see what you analyzed this session.
things_analyzed_deep='|'; % this will write to a log file saved in the folder.

all_days_all_analysis=struct('data',[]); % declare the structure where you wills ave things.
% end

%% Adding data related to the experiment (e.g. lesion dates, WN dates)
if strcmp(input('do you want to enter lesion dates? (y or n) ','s'),'y')==1;
    all_days_all_analysis.summary_of_experiment.inter_lesion_days=input('inter-lesion days (i.e. [A B C...] where A=days prelesion, B=days between lesion 1 and 2...');
end

if strcmp(input('do you want to input WN training days? (y or n) ','s'),'y')==1;
    all_days_all_analysis.summary_of_experiment.WN_baseline_days=input('how many baseline (pre-WN) days?');
    all_days_all_analysis.summary_of_experiment.WN_driving_days=input('how many days driving (pre-consolidation) with WN?');
    all_days_all_analysis.summary_of_experiment.WN_consolidation_days=input('how many days of consolidation');
end



%% compiling data for lt_get_all_transitions...
if strcmp(transitions_lt_question,'y')==1;
    disp('lt_get_all_transition...: compiling data')
    cd(['/' curr_dir '/all_days_transition_matrix_' transitions_phrase])
    try
        [batch]=lt_write_all_folder_contents_to_batch(0,'x'); % put all the day .mat files into a batch file
    catch
        err
        disp(['cannot find .mat analysis files lt_get_all_transitions.'])        
    end
    % to plot above data
    %     batch=input('what is name of batch? ','s');
    %     syllable=input('what is the name of the syllable? ','s');
    
    all_days_all_analysis=lt_save_data_over_experiment_to_structure_TRANSITIONS_v2(batch, first_day, all_days_all_analysis);
    
    %     saveas(figure(gcf),[metrics{i} '_over_days_of_' num2str(syllables(j))],'fig')
    
    things_analyzed=[things_analyzed, '_AllTrans'];
    things_analyzed_deep=[things_analyzed_deep '_AllTrans_|'];
    cd(curr_dir);
end

%% compileing data from lt_db_seq_func_day_save_..., syllable metrics:
% Going thru each metric folder, each syllable, and each day's data and compiling into one structure.

if strcmp(seq_func_day_question,'y')==1;
    disp('Syllable metrics: compiling data from lt_db_seq_func_day_save_rd3gr35 analysis')
    i_max=length(metrics);
    for i=1:i_max; % for each metric
        for j=1:length(syllables_seq_func);
            
            cd(['/' curr_dir '/all_days_' metrics{i} '_' phrase])
            try
                [batch]=lt_write_all_folder_contents_to_batch(syllables_seq_func(j)); % put all the day .mat files into a batch file
            catch
                err
                disp(['cannot find .mat analysis files for syllable: ' syllables_seq_func(j) '. skipping to next syllable.'])
                continue
            end
            % to plot above data
            %     batch=input('what is name of batch? ','s');
            %     syllable=input('what is the name of the syllable? ','s');
            
            all_days_all_analysis=lt_save_data_over_experiment_to_structure_SEQFUNC_v2(batch,syllables_seq_func(j),metrics{i},0,first_day,all_days_all_analysis);
            
            %         saveas(figure(gcf),[metrics{i} '_over_days_of_' num2str(syllables(j))],'fig')
        end
    end
    
    things_analyzed=[things_analyzed, '_' metrics];
    for i=1:length(metrics);
        things_analyzed_deep=[things_analyzed_deep '_' metrics{i}];
    end
    things_analyzed_deep=[things_analyzed_deep '_' syllables_seq_func '_|'];
    cd(curr_dir);
end

%% pitch

% %% saving and entering parameters used in this function
% all_days_all_analysis.date_of_analysis
% all_days_all_analysis.parameters. % e.g. dates looked over,
% sylallbles (depend on if collated or entered)
%
% also put in things such as baseline days, days of lesion, etc., things related to plotting it.

if strcmp(pitch_question,'y')
    disp('Pitch: compiling data from db_contour... analysis')
    cd(['all_days_pitch_' phrase]);
    load(name_of_bird);
    if strcmp(multiple_pitch.parameters.days{1},first_day)==1;
        number_of_days=multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1;
        for i=1:length(syllables_pitch);
            for j=1:number_of_days;
                all_days_all_analysis.data.all_songs.(syllables_pitch(i)).pitch.FF{j}=...
                    multiple_pitch.FF.(syllables_pitch(i)).time_and_FF{j};
                all_days_all_analysis.data.summary_of_day.(syllables_pitch(i)).pitch.mean{j}=multiple_pitch.FF.(syllables_pitch(i)).mean_FF{j};
                all_days_all_analysis.data.summary_of_day.(syllables_pitch(i)).pitch.sd{j} = multiple_pitch.FF.(syllables_pitch(i)).sd_FF{j};
            end
            all_days_all_analysis.summary_of_experiment.(syllables_pitch(i)).pitch.baseline_mean=multiple_pitch.FF.(syllables_pitch(i)).baseline_mean;
        end
        
        all_days_all_analysis.parameters.pitch=multiple_pitch.parameters;
    else
        disp('days previously used in db_contour... to get pitch data do not match the dates you just entered. will not compile pitch data.')
    end
    
    
    things_analyzed=[things_analyzed, '_Pitch'];
    things_analyzed_deep=[things_analyzed_deep '_Pitch_' syllables_pitch '_|'];

    cd(curr_dir)
end

%% Syllable renditions
% ANAlysis computed using lt_all_days_various_calculations, and saved in
% that structure within the all_days_various... folder in bird folder.
% Here I extract that info.
% Data: for each syllable, how many 

if strcmp(syl_rendition_question,'y')==1;
    disp('Syl Rendition Amount');
    [all_days_various]=lt_load_all_days_various_structure(bluejay_num,name_of_bird,date_and_time_of_various_struct_SylRend);
    for i=1:length(all_days_various.syl_rendition_amount.syllables); % use set of syllables specific to this analysis, as there are syls for which I rendition amnt is calculated, but structure analysis was not.
        if strcmp(all_days_various.parameters.days{1},first_day)==1;
            number_of_days=all_days_various.parameters.days{4}-all_days_various.parameters.days{3}+1;
            for j=1:number_of_days;
                all_days_all_analysis.data.summary_of_day.(all_days_various.syl_rendition_amount.syllables(i)).syl_rendition_amount{j}= ...
                    all_days_various.syl_rendition_amount.syl_renditions(i,j);
                all_days_all_analysis.data.summary_of_day.song_amount_labeled{j}= ...
                    all_days_various.syl_rendition_amount.song_amount(j);
            end
        else
            disp('starting day from old analysis in all_days_various_calculations does not match the first day in current all_days_all_analysis program. will not compile syl_rendition data');
        end
    end
    things_analyzed=[things_analyzed '_RendAmnt'];
    things_analyzed_deep=[things_analyzed_deep '_RendAmnt_' all_days_various.syl_rendition_amount.syllables '_|'];

end

%% Compile data for durations of all vocalizations (not just labeld ones) in labeled songs.

if strcmp(get_all_vocalization_stats_question,'y')==1;
    disp('Getting All Vocalization Stats');
    [all_days_various]=lt_load_all_days_various_structure(bluejay_num,name_of_bird,date_and_time_of_various_struct_AllVocalStats);
    if strcmp(all_days_various.parameters.days{1},first_day)==1;
        number_of_days=all_days_various.parameters.days{4}-all_days_various.parameters.days{3}+1;
        for j=1:number_of_days;
            all_days_all_analysis.data.all_songs.all_vocalization_stats.durations_all_songs{j}= all_days_various.get_all_vocalization_stats.durations_all_songs{j};
            all_days_all_analysis.data.per_song.all_vocalization_stats.durations_per_song{j}= all_days_various.get_all_vocalization_stats.durations_per_song{j};
        end
    else
        disp('starting day from old analysis in all_days_various_calculations does not match the first day in current all_days_all_analysis program. will not compile vocalization stats data');
    end
    things_analyzed=[things_analyzed '_VocalStats'];
        things_analyzed_deep=[things_analyzed_deep '_VocalStats_|'];

end

    
%% Compile data for transition probabilities.  IN PROGRESS

%% Compile data for gap durations
% first do lt_all_days_all_various_v2, which uses findwnote2tw to get
% labels, onset, and offsets, for all days (saves to all_days_all_var...
% structure, then use this to extract gap durations for all days.

if strcmp(gap_duration_question,'y')==1; disp('getting gap durations');
    [days_gd gap_durations_over_day]=lt_all_days_get_gap_durations_v2(transitions_gd, date_and_time_of_various_struct_GapDur, curr_dir);
    if strcmp(days_gd{1},first_day)==1;
        num_days=days_gd{4}-days_gd{3}+1;
        
        % put output structure into correct structural format
        for tr=1:size(transitions_gd,2);
            for j=1:num_days;
                try % USE TRY becuase sometimes some days are empty, and would give error
                    all_days_all_analysis.data.all_songs.(transitions_gd{tr}).gap_duration{j}=gap_durations_over_day{j}.(transitions_gd{tr})';
                catch err
                end
            end
        end
        
        things_analyzed=[things_analyzed '_GapDur'];
        things_analyzed_deep=[things_analyzed_deep '_GapDur_'];

        for i=1:length(transitions_gd);
            if i==length(transitions_gd);
                things_analyzed_deep=[things_analyzed_deep transitions_gd{i} '_|'];
            else
                things_analyzed_deep=[things_analyzed_deep transitions_gd{i} '_'];
            end
        end
            else
                disp('days previously used in ...all_various... to get gap duration data do not match the dates you just entered. will not compile data.')
        end
        cd(curr_dir);
end

%% Compile data for WN hit rates
% first do lt_all_days_various...v2. data will be stores in the
% all_days_various structure

if strcmp(WN_trig_question,'y')==1;
    disp('compiling WN trig rate data')
    [all_days_various]=lt_load_all_days_various_structure(bluejay_num,name_of_bird,date_and_time_of_various_struct_WNtrig);
    if strcmp(all_days_various.parameters.days{1},first_day)==1;
        number_of_days=all_days_various.parameters.days{4}-all_days_various.parameters.days{3}+1;
        trig_syl=all_days_various.trig_label.trig_syl;
        for j=1:number_of_days;
            for ii=1:length(trig_syl);
                all_days_all_analysis.data.WN_hit_rate.per_song.(trig_syl(ii)).hits_labels_allhits{j}=all_days_various.trig_label.vals.(trig_syl(ii)){j};
                all_days_all_analysis.data.WN_hit_rate.per_song.(trig_syl(ii)).all_values{j}=all_days_various.trig_label.trigs.(trig_syl(ii)){j};
                % compile over day
                all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}=sum(all_days_various.trig_label.vals.(trig_syl(ii)){j});
                % get stats (hit rate, false positive over hits)
                try
                    all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_divide_labels{j}=all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}(1)/all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}(2); %hit/labels
                    all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).FalsePos{j}=all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}(3)-all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}(1); % amnt of f.p.
                    all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).FalsePos_divide_TotalWN{j}=all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).FalsePos{j}/all_days_all_analysis.data.WN_hit_rate.all_songs.(trig_syl(ii)).hits_labels_allhits{j}(3);  % f.p./AllWN
                catch err
                end
            end
        end
        things_analyzed=[things_analyzed '_WNtrigLabel'];
        things_analyzed_deep=[things_analyzed_deep '_WNtrigLabel_' trig_syl '_|'];

    else
        disp('WN trig rate compilation error: starting day from old analysis in all_days_various_calculations does not match the first day in current all_days_all_analysis program. will not compile syl_rendition data');
    end
    cd(curr_dir);
end

%% SAVING

% FIRST compile all variables into structure
timestamp=lt_get_timestamp(0);

all_days_all_analysis.things_analyzed=things_analyzed;

all_days_all_analysis.parameters.date_of_analysis=timestamp;
all_days_all_analysis.parameters.first_day=first_day;
all_days_all_analysis.parameters.last_day=last_day;
all_days_all_analysis.parameters.phrase=phrase;
all_days_all_analysis.parameters.name_of_bird=name_of_bird;

all_days_all_analysis.parameters.seq_func_day_question=seq_func_day_question;
all_days_all_analysis.parameters.pitch_question=pitch_question;
all_days_all_analysis.parameters.syl_rendition_question=syl_rendition_question;
all_days_all_analysis.parameters.get_all_vocalization_stats_question=get_all_vocalization_stats_question;
all_days_all_analysis.parameters.gap_duration_question=gap_duration_question;
all_days_all_analysis.parameters.WN_trig_question=WN_trig_question;
all_days_all_analysis.parameters.transitions_lt_question=transitions_lt_question;

if strcmp(transitions_lt_question,'y')==1;
    all_days_all_analysis.parameters.all_trans.transitions_phrase=transitions_phrase;
end

if strcmp(WN_trig_question,'y')==1;
    all_days_all_analysis.parameters.WN_trig.date_and_time_of_various_struct=date_and_time_of_various_struct_WNtrig;
    all_days_all_analysis.parameters.WN_trig.trig_syl=trig_syl;
end    

if strcmp(seq_func_day_question,'y')==1;
    all_days_all_analysis.parameters.syllables_seq_func=syllables_seq_func;
%     all_days_all_analysis.parameters.were_syllables_chosen_automatically=were_syllables_chosen_automatically;
end

if strcmp(pitch_question,'y')==1;
    all_days_all_analysis.parameters.syllables=syllables_pitch;
%     all_days_all_analysis.parameters.were_syllables_chosen_automatically=were_syllables_chosen_automatically;
end


if strcmp(gap_duration_question,'y')==1;
    all_days_all_analysis.parameters.gap_duration.transitions=transitions_gd;
    all_days_all_analysis.parameters.gap_duration.date_and_time_of_various_struct=date_and_time_of_various_struct_GapDur;
    all_days_all_analysis.parameters.gap_duration.dates=days_gd;
end

if strcmp(seq_func_day_question,'y')==1;
    all_days_all_analysis.parameters.seq_func_day.metrics=metrics;
end

if strcmp(get_all_vocalization_stats_question,'y')==1;
    all_days_all_analysis.parameters.gap_duration.date_and_time_of_various_struct=date_and_time_of_various_struct_AllVocalStats;
end

if strcmp(syl_rendition_question,'y')==1;
    all_days_all_analysis.parameters.gap_duration.date_and_time_of_various_struct=date_and_time_of_various_struct_SylRend;
end

%% specific to v2, change name of structure

all_days_all_analysis_v2=all_days_all_analysis; % only save v2.

%% SECOND, save
try cd([curr_dir '/all_days_all_analysis_v2'])
catch error
    mkdir([curr_dir '/all_days_all_analysis_v2'])
    cd([curr_dir '/all_days_all_analysis_v2'])
end

%         if ~exist('syllables','var')==1; % use this just for naming structure
%             syllables='';
%         end
%
% if ~exist('transitions_gd','var')==1; % use this just for naming structure
%     transitions_string='';
% else
%     transitions_string=[];
%     for ii=1:size(transitions_gd,2);
%         transitions_string=[transitions_string '_' transitions_gd{ii}];
%     end
% end

% save(['all_days_all_analysis_' timestamp '_' phrase things_analyzed{1:end} '_' syllables transitions_string],'all_days_all_analysis')
save(['all_days_all_analysis_v2_' timestamp '_' phrase things_analyzed{1:end}],'all_days_all_analysis_v2')



% LOGGING ANALYSIS INTO log.txt ------------------------------------
log_note=[timestamp '__|' things_analyzed_deep '|__' first_day '_to_' last_day];
fid=fopen('log.txt','a');
fprintf(fid,'%s \n \n',log_note);
fclose(fid);

cd(curr_dir)

%% To plot stuff, use lt_all_days_all_analysis_PLOT.

cd([curr_dir '/all_days_all_analysis_v2'])
if isdir('PLOT')==0;
    mkdir('PLOT')
end

if strcmp(input('do you want to use lt_all_days_all_analysis_PLOT to plot data? (y or n) ','s'),'y')==1;
    open('lt_all_days_all_analysis_PLOT')
end

