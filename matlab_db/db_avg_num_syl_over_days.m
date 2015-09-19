%Calculates the average occurences of a syllable per song. Used to see if
%100% white noise reduces wiener entropy and causes syllable to decrease in
%frequency
%
%Written by DM Brady 11/2012


%% Making an all average number of syllable per song folder on your computer

% write the name of your computer and/or directory here, it will make an
% all_calls_analysis folder that will save your three large structure
% variables and folders with figures and workspace variables for each bird.

which_computer = input('Bluejay number? ');
avg_num_syl.parameters.computer = ['/bluejay' num2str(which_computer) '/dbrady'];
clear which_computer

if exist([avg_num_syl.parameters.computer '/all_avg_num_syl'], 'dir') ~= 7
    mkdir([avg_num_syl.parameters.computer '/all_avg_num_syl'])
else
end

%% Entering or loading parameters for your bird

% asks for name of bird. your batch file should be named birdname_date.
avg_num_syl.parameters.nameofbird = input('What is the name of your bird? \n(format should be name of batch file excluding date)  ', 's');

%asks for experiment type (either pitch, seq, or entropy)
avg_num_syl.parameters.experiment_type = input('What type of experiment did you do? \n(pitch, sequence, entropy)  ', 's');

% makes a directory for a particular bird if it does not exist yet.
if exist([avg_num_syl.parameters.computer '/all_avg_num_syl/' avg_num_syl.parameters.experiment_type...
        '/' avg_num_syl.parameters.nameofbird], 'dir') ~= 7;
    mkdir([avg_num_syl.parameters.computer '/all_avg_num_syl/' avg_num_syl.parameters.experiment_type...
        '/' avg_num_syl.parameters.nameofbird])
else
end

% tries to load structure with data from previous runs
try
    load([ avg_num_syl.parameters.computer '/all_avg_num_syl/' avg_num_syl.parameters.experiment_type '/'...
        avg_num_syl.parameters.nameofbird '/' avg_num_syl.parameters.nameofbird '.mat'])
catch err
    display('Cannot find files, hope this is your first time running this program on this bird')
    display(' ')
    
    
    % asks for the syllables for the particular bird
    avg_num_syl.parameters.syllables =...
        input(['What are the syllables for ' avg_num_syl.parameters.nameofbird '?\n(no spaces between letters)  '], 's');
    
    avg_num_syl.parameters.days{1} =...
        input(['What is the first day for ' avg_num_syl.parameters.nameofbird '?\n(format should be day month year. i.e. 06Jun2012)    '], 's');
    avg_num_syl.parameters.days{3} = datenum(avg_num_syl.parameters.days{1});
    
    % Asks for length before WN delivery, pitch shift, consolidation, etc.
    avg_num_syl.parameters.duration.begin = input('How many days before WN delivery?  ');
    avg_num_syl.parameters.duration.shift = input('How many days of driving learning?  ');
    avg_num_syl.parameters.duration.consolidation = input('How many days of consolidation?  ');
end

% puts a time stamp on all figures and variables when you did the analysis
avg_num_syl.parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM');
avg_num_syl.parameters.today_date = avg_num_syl.parameters.today_date(avg_num_syl.parameters.today_date ~= ' ');

% adds the last date  for the bird. you must have already created batch
% files for the dates of interest.
avg_num_syl.parameters.days{2} =...
    input(['What is the last day for ' avg_num_syl.parameters.nameofbird '?\n(format should be day month year. i.e. 31Dec2013)   '], 's');
avg_num_syl.parameters.days{4} = datenum(avg_num_syl.parameters.days{2});

%checks to make sure you have at least two days
if avg_num_syl.parameters.days{3} == avg_num_syl.parameters.days{4}
    display('Sorry, need to have at least two days before you can run this program')
    return
else
end


%% Making structures for results of avg num syl/song analysis

fprintf('\nGetting number of syllables/song ...')

for i = 1:length(avg_num_syl.parameters.syllables)
    for j = avg_num_syl.parameters.days{3}:avg_num_syl.parameters.days{4}
        
        %if ~isempty checks to see if cell has been made
        %for this bird, syllable, and day
        
        try
            iscell(avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1});
            iscell(avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1});
            iscell(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1});
            iscell(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1});
        catch err
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1} = {};
        end
            
    end
end

%% Getting number of syllable per song


curr_directory = pwd;


for i = 1:length(avg_num_syl.parameters.syllables)
    for j = avg_num_syl.parameters.days{3}:avg_num_syl.parameters.days{4}
        
        %if ~isempty checks to see if findnote2tw or jc_pitchcontour
        %has been performed on this bird, syllable, and day
        try
            cd([avg_num_syl.parameters.computer '/' avg_num_syl.parameters.nameofbird '/all_calls/' datestr(j, 'ddmmmyyyy')])
            
            try
                if ~isempty(avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1}) &&...
                        ~isempty(avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1}) &&...
                        ~isempty(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1}) &&...
                        ~isempty(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1})
                    continue
                else
                    try
                        %gives the number of the that syllable in each song
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_get_num_syl_batch([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')],...
                            avg_num_syl.parameters.syllables(i));
                        
                        %says when the song occurred
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_get_timing_from_batch([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')]);
                        
                        %gets the length of the song
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_song_length([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')]);
                        
                        %calculates the rate of syllable singing
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1} = ...
                            avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1}./...
                            avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1};
                        
                    catch err
                        db_batch_convert([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')], 'temp_cbin');
                        
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_get_num_syl_batch('temp_cbin',...
                            avg_num_syl.parameters.syllables(i));
                        
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_get_timing_from_batch([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')]);
                        
                        %gets the length of the song
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1} = ...
                            db_song_length([avg_num_syl.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')])./1000;
                        
                        %calculates the rate of syllable singing
                        avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1} = ...
                            avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song{j-avg_num_syl.parameters.days{3}+1}./...
                            avg_num_syl.(avg_num_syl.parameters.syllables(i)).length_song{j-avg_num_syl.parameters.days{3}+1};
                        
                    end
                end
            catch err
                continue
            end
            
            cd([avg_num_syl.parameters.computer '/' avg_num_syl.parameters.nameofbird '/all_calls'])
        catch err
            display([ datestr(j, 'ddmmmyyyy') ' is missing'])
            continue
        end
    end
end

cd(curr_directory)
clear curr_directory

fprintf('done!\n')

%% Calculating the mean, sd, median, and iqr for each day

fprintf('\nCalculating descriptive statistics of syllable rate per song...')

for i = 1:length(avg_num_syl.parameters.syllables)
    for j = 1:length(avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song)
        try
            %mean number of syllables per song per day
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.mean{j} =...
                mean(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j});
            
            %sd of number of syllables per song per day
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.sd{j} =...
                std(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j});
            
            %median number of syllables per song per day
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.median{j} =...
                median(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j});
            
            %iqr of number of syllables per song per day
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.iqr{j} =...
                iqr(avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j});
        catch err
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.mean{j} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.sd{j} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.iqr{j} = {};
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.median{j} = {};
        end
    end
end

fprintf('done!\n')

 %% Calculating baseline syl/song

fprintf('\nCalculating baseline syllable rate...')
for i = 1:length(avg_num_syl.parameters.syllables)
    %calcualtes baseline syllables/song
    try
        for j = 1:avg_num_syl.parameters.duration.begin
            temp_base.mean(j) = avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.mean{j};
        end
        avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.baseline_mean = mean(temp_base.mean);
    catch err
        continue
    end
    
    clear('temp_base')
    
end
fprintf('done!\n')

%% Figure of syl/song values over time

%Determining the color for the various conditions (grey = preWN, red =
%shift, blue = consolidation, black = postWN)
fprintf('\nDetermining color of each day according to epoch during paradigm...')
for i = 1:length(avg_num_syl.parameters.syllables)
    for j = 1:avg_num_syl.parameters.days{4}-avg_num_syl.parameters.days{3}+1
        if j <= avg_num_syl.parameters.duration.begin
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).color{j} = [.6 .6 .6];
        elseif j > avg_num_syl.parameters.duration.begin && j <= avg_num_syl.parameters.duration.begin+avg_num_syl.parameters.duration.shift
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).color{j} = [1 0 0];
        elseif j > avg_num_syl.parameters.duration.begin+avg_num_syl.parameters.duration.shift &&...
                j <= avg_num_syl.parameters.duration.begin+avg_num_syl.parameters.duration.shift+avg_num_syl.parameters.duration.consolidation
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).color{j} = [0 0 1];
        elseif j > avg_num_syl.parameters.duration.begin+avg_num_syl.parameters.duration.shift+avg_num_syl.parameters.duration.consolidation
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).color{j} = [.3 .3 .3];
        end
    end
end
fprintf('done!\n')
        
fprintf('\nCreating figures of syl/song over time...')
for i = 1:length(avg_num_syl.parameters.syllables)
    figure(), hold on
    title(['Syllable rate for ' avg_num_syl.parameters.nameofbird '   duration: '...
        num2str(avg_num_syl.parameters.days{4}-avg_num_syl.parameters.days{3}+1)...
        ' days   syllable: ' avg_num_syl.parameters.syllables(i)])
    ylabel('Syllable/sec')
    xlabel('Time (days)')
    xlim([-1-avg_num_syl.parameters.duration.begin...
        length(avg_num_syl.(avg_num_syl.parameters.syllables(i)).per_song)+1])
    for j = avg_num_syl.parameters.days{3}:avg_num_syl.parameters.days{4}
        try

            %plots syl/song for time range
            plot(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1}-...
                min(floor(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{1}))-...
                avg_num_syl.parameters.duration.begin,...
                avg_num_syl.(avg_num_syl.parameters.syllables(i)).rate{j-avg_num_syl.parameters.days{3}+1},...
                '.',...
                'Color', avg_num_syl.(avg_num_syl.parameters.syllables(i)).color{j-avg_num_syl.parameters.days{3}+1})
            
            %plots syl/song mean per day
            plot(median(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1})-...
                min(floor(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{1}))-...
                avg_num_syl.parameters.duration.begin,...
                avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.mean{j-avg_num_syl.parameters.days{3}+1},...
                'ko',...
                'MarkerFaceColor','k',...
                'LineWidth',2)
            
            %plots entropy sd per day
            errorbar(median(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{j-avg_num_syl.parameters.days{3}+1})-...
                min(floor(avg_num_syl.(avg_num_syl.parameters.syllables(i)).time_of_day{1}))-...
                avg_num_syl.parameters.duration.begin,...
                avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.mean{j-avg_num_syl.parameters.days{3}+1},...
                avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.sd{j-avg_num_syl.parameters.days{3}+1},...
                'ko',...
                'LineWidth',2)
            

        catch err
            continue
        end
    end
    
    ylim(get(gca, 'YLim'));
    
    %plots lines for WN begin, consolidation begin, and WN end
    %line for WN begin
    try
        line([0 0],...
            [min(get(gca, 'YLim')) max(get(gca, 'YLim'))],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for consolidation begin
    try
        line([avg_num_syl.parameters.duration.shift...
            avg_num_syl.parameters.duration.shift],...
            [min(get(gca, 'YLim')) max(get(gca, 'YLim'))],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for WN end
    try
        line([avg_num_syl.parameters.duration.shift+avg_num_syl.parameters.duration.consolidation...
            avg_num_syl.parameters.duration.shift+avg_num_syl.parameters.duration.consolidation],...
            [min(get(gca, 'YLim')) max(get(gca, 'YLim'))],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for baseline entropy
    try
        line(get(gca, 'XLim'),...
            [avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.baseline_mean...
            avg_num_syl.(avg_num_syl.parameters.syllables(i)).describe.baseline_mean],...
            'Color', [.3 .3 .3],...
            'LineStyle', '--')
    catch err
        continue
    end
    
    
    
    saveas(figure(gcf), [avg_num_syl.parameters.computer '/all_avg_num_syl/' avg_num_syl.parameters.experiment_type '/'...
        avg_num_syl.parameters.nameofbird...
        '/Avg_syl_over_time_' avg_num_syl.parameters.nameofbird '_' 'duration_'...
         num2str(avg_num_syl.parameters.days{4}-avg_num_syl.parameters.days{3}+1)...
        'days' '_' avg_num_syl.parameters.syllables(i)], 'fig')
end

fprintf('done!\n')
        

%% Save data

clear i
clear j

fprintf('\nSaving data...')

save([avg_num_syl.parameters.computer '/all_avg_num_syl/' avg_num_syl.parameters.experiment_type '/'...
        avg_num_syl.parameters.nameofbird '/' avg_num_syl.parameters.nameofbird '.mat'], 'avg_num_syl', '-v7.3')

fprintf('done!\n')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
