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
name_of_bird=input('name of bird? ','s');
curr_dir=pwd;

% asking which data you would like to compile
seq_func_day_question=input('compile syllable structure data from lt_db_seq_func_day_save_rd3gr35? ','s');
pitch_question=input('compile pitch data from db_contour_and_FF_analysis_over_time_v3 ','s');
syl_rendition_question=input('compile data on syl and song rendition amounts? (i.e. from all_days_various_calculations?) ','s');
get_all_vocalization_stats_question=input('compile data on stats of all vocalizations in labeled songs? (not just labeled) ','s');
gap_duration_question=input('compile gap_duration stats? ','s');

if strcmp(gap_duration_question,'y')==1;
    transitions_gd=input('what transitions to analyze for gap durations? (e.g. {"ab" "cd"})');
    cd('all_days_various_calculations'); ls
    date_and_time_of_various_struct=input('enter date and time of struct you want to load for gap duration analysis (e.g. {"06Jan2014","1737"}');
    cd(curr_dir);
end

% if strcmp(input('automatically find all syllables that were used today? (y or n) ','s'),'y')==1;
%     % ENTER FUNCTION TO PULL OUT SYLLABLES HERE
%     were_syllables_chosen_automatically=1;
% else

if strcmp(seq_func_day_question,'y')==1||strcmp(pitch_question,'y')==1;
syllables = input('compile data for which syllables? (for seq_func and pitch analyses)','s');
were_syllables_chosen_automatically=0;
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

if strcmp(get_all_vocalization_stats_question,'y')==1 || strcmp(syl_rendition_question,'y')==1|| strcmp(syl_rendition_question,'y')==1;
    disp('will use all_days_various_calculations previously calculated data. MANUALLY go to folder and load latest structure, then type "return"');
    disp('warning, if doing gap duration as well, be careful.  GD code can overwrite the structure you will now load')
    keyboard
    cd(curr_dir);
end


first_day = input('first day of analysis?','s');
phrase=input('phrase identifying this experiment? ','s');
things_analyzed=cell(0); % this will be a high level field that you can glance at quickly to see what you analyzed this session.

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


%% compileing data from lt_db_seq_func_day_save_..., syllable metrics:
% Going thru each metric folder, each syllable, and each day's data and compiling into one structure.

if strcmp(seq_func_day_question,'y')==1;
    disp('Syllable metrics: compiling data from lt_db_seq_func_day_save_rd3gr35 analysis')
    
    i_max=length(metrics);
    
    for i=1:i_max; % for each metric
        for j=1:length(syllables);
            
            cd(['/' curr_dir '/all_days_' metrics{i} '_' phrase])
            try
                [batch]=lt_write_all_folder_contents_to_batch(syllables(j)); % put all the day .mat files into a batch file
            catch
                err
                disp(['cannot find .mat analysis files for syllable: ' syllables(j) '. skipping to next syllable.'])
                continue
            end
            % to plot above data
            %     batch=input('what is name of batch? ','s');
            %     syllable=input('what is the name of the syllable? ','s');
            
            all_days_all_analysis=lt_save_data_over_experiment_to_structure(batch,syllables(j),metrics{i},0,first_day,all_days_all_analysis);
            
            %         saveas(figure(gcf),[metrics{i} '_over_days_of_' num2str(syllables(j))],'fig')
        end
    end
    
    things_analyzed=[things_analyzed, '_' metrics];
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
        for i=1:length(syllables);
            for j=1:number_of_days;
                all_days_all_analysis.data{j}.all_songs.(syllables(i)).pitch.FF=...
                    multiple_pitch.FF.(syllables(i)).time_and_FF{j};
                all_days_all_analysis.data{j}.summary_of_day.(syllables(i)).pitch.mean=multiple_pitch.FF.(syllables(i)).mean_FF{j};
                all_days_all_analysis.data{j}.summary_of_day.(syllables(i)).pitch.sd = multiple_pitch.FF.(syllables(i)).sd_FF{j};
            end
            all_days_all_analysis.summary_of_experiment.(syllables(i)).pitch.baseline_mean=multiple_pitch.FF.(syllables(i)).baseline_mean;
        end
        
        all_days_all_analysis.parameters.pitch=multiple_pitch.parameters;
    else
        disp('days previously used in db_contour... to get pitch data do not match the dates you just entered. will not compile pitch data.')
    end
    
    
    things_analyzed=[things_analyzed, '_Pitch'];
    cd(curr_dir)
end

%% Syllable renditions
% ANAlysis computed using lt_all_days_various_calculations, and saved in
% that structure within the all_days_various... folder in bird folder.
% Here I extract that info.
% Data: for each syllable, how many 

if strcmp(syl_rendition_question,'y')==1;
    disp('Syl Rendition Amount: using previously loaded all_days_various_calculations structure to compile syl. rendition data...');
    for i=1:length(all_days_various.parameters.syllables); % use set of syllables specific to this analysis, as there are syls for which I rendition amnt is calculated, but structure analysis was not.
        if strcmp(all_days_various.parameters.days{1},first_day)==1;
            number_of_days=all_days_various.parameters.days{4}-all_days_various.parameters.days{3}+1;
            for j=1:number_of_days;
                all_days_all_analysis.data{j}.summary_of_day.(all_days_various.parameters.syllables(i)).syl_rendition_amount= ...
                    all_days_various.syl_rendition_amount.syl_renditions(i,j);
                all_days_all_analysis.data{j}.summary_of_day.song_amount_labeled= ...
                    all_days_various.syl_rendition_amount.song_amount(j);
            end
        else
            disp('starting day from old analysis in all_days_various_calculations does not match the first day in current all_days_all_analysis program. will not compile syl_rendition data');
        end
    end
    things_analyzed=[things_analyzed '_RendAmnt'];
end

%% Compile data for durations of all vocalizations (not just labeld ones) in labeled songs.

if strcmp(get_all_vocalization_stats_question,'y')==1;
    disp('Vocalization durations: compiling stats from lt_get_all_vocalization_stats.');
    if strcmp(all_days_various.parameters.days{1},first_day)==1;
        number_of_days=all_days_various.parameters.days{4}-all_days_various.parameters.days{3}+1;
        for j=1:number_of_days;
            all_days_all_analysis.data{j}.all_songs.all_vocalization_stats.durations_all_songs= all_days_various.get_all_vocalization_stats.durations_all_songs{j};
            all_days_all_analysis.data{j}.per_song.all_vocalization_stats.durations_per_song= all_days_various.get_all_vocalization_stats.durations_per_song{j};
        end
    else
        disp('starting day from old analysis in all_days_various_calculations does not match the first day in current all_days_all_analysis program. will not compile vocalization stats data');
    end
    things_analyzed=[things_analyzed '_VocalStats'];
end

    
%% Compile data for transition probabilities.  IN PROGRESS

%% Compile data for gap durations
% first do lt_all_days_all_various_v2, which uses findwnote2tw to get
% labels, onset, and offsets, for all days (saves to all_days_all_var...
% structure, then use this to extract gap durations for all days.

[days_gd gap_durations_over_day]=lt_all_days_get_gap_durations(transitions_gd, date_and_time_of_various_struct, curr_dir);
if strcmp(days_gd{1},first_day)==1;
    num_days=days_gd{4}-days_gd{3}+1;
    
    % put output structure into correct structural format
    for tr=1:size(transitions_gd,2);
        for j=1:num_days;
            try % USE TRY becuase sometimes some days are empty, and would give error
                all_days_all_analysis.data{j}.all_songs.(transitions_gd{tr}).gap_duration=gap_durations_over_day{j}.(transitions_gd{tr})';
            catch err
            end
        end
    end
    
    things_analyzed=[things_analyzed '_GapDur'];
else
    disp('days previously used in ...all_various... to get gap duration data do not match the dates you just entered. will not compile data.')
end
cd(curr_dir);

%% SAVING

% FIRST compile all variables into structure
timestamp=lt_get_timestamp(0);

all_days_all_analysis.things_analyzed=things_analyzed;

all_days_all_analysis.parameters.date_of_analysis=timestamp;
all_days_all_analysis.parameters.first_day=first_day;
all_days_all_analysis.parameters.phrase=phrase;
all_days_all_analysis.parameters.name_of_bird=name_of_bird;

all_days_all_analysis.parameters.seq_func_day_question=seq_func_day_question;
all_days_all_analysis.parameters.pitch_question=pitch_question;
all_days_all_analysis.parameters.syl_rendition_question=syl_rendition_question;
all_days_all_analysis.parameters.get_all_vocalization_stats_question=get_all_vocalization_stats_question;
all_days_all_analysis.parameters.gap_duration_question=gap_duration_question;

if strcmp(seq_func_day_question,'y')==1||strcmp(pitch_question,'y')==1;
    all_days_all_analysis.parameters.syllables=syllables;
    all_days_all_analysis.parameters.were_syllables_chosen_automatically=were_syllables_chosen_automatically;
end

if strcmp(gap_duration_question,'y')==1;
    all_days_all_analysis.parameters.gap_duration.transitions=transitions_gd;
    all_days_all_analysis.parameters.gap_duration.date_and_time_of_various_struct=date_and_time_of_various_struct;
    all_days_all_analysis.parameters.gap_duration.dates=days_gd;
end

if strcmp(seq_func_day_question,'y')==1;
    all_days_all_analysis.parameters.seq_func_day.metrics=metrics;
end

%% SECOND, save
try cd([curr_dir '/all_days_all_analysis'])
catch error
    mkdir([curr_dir '/all_days_all_analysis'])
    cd([curr_dir '/all_days_all_analysis'])
end

if ~exist('syllables','var')==1; % use this just for naming structure
    syllables=[];
end

if ~exist('transitions_gd','var')==1; % use this just for naming structure
    transitions_gd=[];
end

save(['all_days_all_analysis_' phrase '_' syllables '_' transitions_gd things_analyzed{1:end} '_' timestamp],'all_days_all_analysis')

%% To plot stuff, use lt_all_days_all_analysis_PLOT.

if isdir('PLOT')==0;
    mkdir('PLOT')
end

if strcmp(input('do you want to use lt_all_days_all_analysis_PLOT to plot data? (y or n) ','s'),'y')==1;
    lt_all_days_all_analysis_PLOT
end

