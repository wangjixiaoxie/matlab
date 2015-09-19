% Modified by L. Tian on 7/24/13, to work with my folders.
%%
% Measures wiener entropy over an experiment
%
% Written by DM Brady 10/2012
%
%Stores information about the power spectrum for each syllable of
%interest. Uses this information to calculate the wiener entropy
%(geomean(powerspectrum)./mean(powerspectrum). db_contour already includes
%this analysis. This is for experiments where FF is not a concern. Folder
%needs to be organized like db_contour (an all calls folder in your bird's
%folder with each date organized in a ddmmmyyyy label). Batch files can be
%.cbin or .cbin.not.mat, but must have a bird_date format.

%% Making an all wiener entropy folder on your computer

% write the name of your computer and/or directory here, it will make an
% all_calls_analysis folder that will save your three large structure
% variables and folders with figures and workspace variables for each bird.

%which_computer = input('Bluejay number? ');
wiener_entropy.parameters.computer = '/home/lucas/data/song';
%clear which_computer

if exist([wiener_entropy.parameters.computer '/all_wiener_entropy'], 'dir') ~= 7
    mkdir([wiener_entropy.parameters.computer '/all_wiener_entropy'])
else
end

%% Entering or loading parameters for your bird

% asks for name of bird. your batch file should be named birdname_date.
wiener_entropy.parameters.nameofbird = input('What is the name of your bird? \n(format should be name of batch file excluding date)  ', 's');

%asks for experiment type (either pitch, seq, or entropy)
wiener_entropy.parameters.experiment_type = input('What type of experiment did you do? \n(pitch, sequence, entropy)  ', 's');

% makes a directory for a particular bird if it does not exist yet.
if exist([wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type...
        '/' wiener_entropy.parameters.nameofbird], 'dir') ~= 7;
    mkdir([wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type...
        '/' wiener_entropy.parameters.nameofbird])
else
end

% tries to load structure with data from previous runs
try
    load([ wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type '/'...
        wiener_entropy.parameters.nameofbird '/' wiener_entropy.parameters.nameofbird '.mat'])
catch err
    display('Cannot find files, hope this is your first time running this program on this bird')
    display(' ')
    
    
    % asks for the syllables for the particular bird
    wiener_entropy.parameters.syllables =...
        input(['What are the syllables for ' wiener_entropy.parameters.nameofbird '?\n(no spaces between letters)  '], 's');
    
    wiener_entropy.parameters.days{1} =...
        input(['What is the first day for ' wiener_entropy.parameters.nameofbird '?\n(format should be day month year. i.e. 06Jun2012)    '], 's');
    wiener_entropy.parameters.days{3} = datenum(wiener_entropy.parameters.days{1});
    
    % Asks for length before WN delivery, pitch shift, consolidation, etc.
    wiener_entropy.parameters.duration.begin = input('How many days before WN delivery?  ');
    wiener_entropy.parameters.duration.shift = input('How many days of driving learning?  ');
    wiener_entropy.parameters.duration.consolidation = input('How many days of consolidation?  ');
end

% puts a time stamp on all figures and variables when you did the analysis
wiener_entropy.parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM');
wiener_entropy.parameters.today_date = wiener_entropy.parameters.today_date(wiener_entropy.parameters.today_date ~= ' ');

% adds the last date  for the bird. you must have already created batch
% files for the dates of interest.
wiener_entropy.parameters.days{2} =...
    input(['What is the last day for ' wiener_entropy.parameters.nameofbird '?\n(format should be day month year. i.e. 31Dec2013)   '], 's');
wiener_entropy.parameters.days{4} = datenum(wiener_entropy.parameters.days{2});

%checks to make sure you have at least two days
if wiener_entropy.parameters.days{3} == wiener_entropy.parameters.days{4}
    display('Sorry, need to have at least two days before you can run this program')
    return
else
end


%% Making structures for results of entropy analysis

fprintf('\nGetting spectrograms ...')

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = wiener_entropy.parameters.days{3}:wiener_entropy.parameters.days{4}
        
        %if ~isempty checks to see if cell has been made
        %for this bird, syllable, and day
        
        try
            iscell(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j-wiener_entropy.parameters.days{3}+1});
            iscell(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1});
            iscell(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).freq{j-wiener_entropy.parameters.days{3}+1});
            iscell(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1});
        catch err
            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j-wiener_entropy.parameters.days{3}+1} = {};
            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1} = {};
            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).freq{j-wiener_entropy.parameters.days{3}+1} = {};
            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1} = {};
            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1} = {};
        end
            
    end
end

%% Getting spectrograms


curr_directory = pwd;


for i = 1:length(wiener_entropy.parameters.syllables)
    for j = wiener_entropy.parameters.days{3}:wiener_entropy.parameters.days{4}
        
        %if ~isempty checks to see if findnote2tw or jc_pitchcontour
        %has been performed on this bird, syllable, and day
        try
            cd([wiener_entropy.parameters.computer '/' wiener_entropy.parameters.nameofbird '/all_calls/' datestr(j, 'ddmmmyyyy')])
            
            try
                if ~isempty(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j-wiener_entropy.parameters.days{3}+1}) &&...
                     ~isempty(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1}) && ...
                     ~isempty(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).freq{j-wiener_entropy.parameters.days{3}+1}) &&...
                     ~isempty(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).filename{j-wiener_entropy.parameters.days{3}+1}) &&...
                     ~isempty(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1})
                    continue
                else
                    try
                        [wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).freq{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).filename{j-wiener_entropy.parameters.days{3}+1}]= ...
                            db_spec_syllable([wiener_entropy.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')],...
                            wiener_entropy.parameters.syllables(i), 0.20, 0.20);
                        wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1} = ...
                            db_glt_et_timing(...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).filename{j-wiener_entropy.parameters.days{3}+1});
                    catch err
                        db_batch_convert([wiener_entropy.parameters.nameofbird '_' datestr(j, 'ddmmmyyyy')], 'temp_cbin');
                        [wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).freq{j-wiener_entropy.parameters.days{3}+1}...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).filename{j-wiener_entropy.parameters.days{3}+1}]= ...
                            db_spec_syllable('temp_cbin', wiener_entropy.parameters.syllables(i), 0.20, 0.20);
                        delete('temp_cbin');
                        wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j-wiener_entropy.parameters.days{3}+1} = ...
                            db_get_timing(...
                            wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).filename{j-wiener_entropy.parameters.days{3}+1});
                    end
                end
            catch err
                continue
            end
            
            cd([wiener_entropy.parameters.computer '/' wiener_entropy.parameters.nameofbird '/all_calls_analysis_pitch'])
        catch err
            display([ datestr(j, 'ddmmmyyyy') ' is missing'])
            continue
        end
    end
end

cd(curr_directory)
clear curr_directory

fprintf('done!\n')



%% Calculating entropy for all time (.2 before to .2 after syllable onset)

fprintf('\nCalculating entropy all time...')

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec)
        try
            %calcualtes entropy along entire timeline (0.2 before and 0.2
            %seconds after start of syllable
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).all_time{j} =...
                db_calc_wiener_entropy(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec{j});
            
        catch err
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).all_time{j} = {};
        end
    end
end
    
 fprintf('done!\n')


%% Figure of entropy contours over time
 
 fprintf('\nPreparing entropy contour figures over time...')
 
 for i = 1:length(wiener_entropy.parameters.syllables)
    figure(), hold on
    title(['Wiener entropy contour for ' wiener_entropy.parameters.nameofbird '   duration: '...
        num2str(wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}+1) ' days'...
        '   syllable: ' wiener_entropy.parameters.syllables(i)])
    xlabel('Time (msec)')
    ylabel('Wiener Entropy')
    for j = wiener_entropy.parameters.days{3}:wiener_entropy.parameters.days{4}
        %plots the mean
        try
            plot(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j-wiener_entropy.parameters.days{3}+1},...
                mean(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).all_time{j-wiener_entropy.parameters.days{3}+1},1),...
                'Linewidth',2,...
                'Color',[1-(j-wiener_entropy.parameters.days{3})./(wiener_entropy.parameters.days{4} - wiener_entropy.parameters.days{3})...
                0 0+(j-wiener_entropy.parameters.days{3})./(wiener_entropy.parameters.days{4} - wiener_entropy.parameters.days{3})])
            %             %plots mean + sd
            %             plot(mean(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}')...
            %                 +std(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}'),...
            %                 'Color',[1-(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})...
            %                 0 0+(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})])
            %             %plots mean -sd
            %              plot(mean(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}')...
            %                 -std(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}'),...
            %                 'Color',[1-(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})...
            %                 0 0+(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})])
        catch err
            continue
        end
    end
    saveas(figure(gcf), [wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type '/'...
        wiener_entropy.parameters.nameofbird ...
        '/WE_contour_' wiener_entropy.parameters.nameofbird '_' 'duration_'...
        num2str(wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}+1)...
        'days' '_' wiener_entropy.parameters.syllables(i)], 'fig')
 end

fprintf('done!\n')

%% Variable for time range for calculating Entropy

% need to enter time range to calculate entropy, use the entropy contour figure to
% determine the best time range.

for i = 1:length(wiener_entropy.parameters.syllables)
    if isfield(wiener_entropy.parameters, 'time_range') == 1
        if isfield(wiener_entropy.parameters.time_range, wiener_entropy.parameters.syllables(i)) == 1
            if strcmpi(input(['Do you want to select a new time range for ' wiener_entropy.parameters.syllables(i) '?\n'...
                      'Current range:  [' num2str(wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i))) ']   ' ],'s'),'y')
                wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i)) =...
                    input(['Start and stop time for entropy plotting of ' wiener_entropy.parameters.nameofbird...
                    ' ' wiener_entropy.parameters.syllables(i) '?\n(format [start stop])  ']);

                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)) =...
                    rmfield(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)),'window');
                
            else
            end
        else
           wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i)) =...
            input(['Start and stop time for entropy plotting of ' wiener_entropy.parameters.nameofbird...
            ' ' wiener_entropy.parameters.syllables(i) '?\n(format [start stop])  ']);
        end
    else
        wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i)) =...
            input(['Start and stop time for entropy plotting of ' wiener_entropy.parameters.nameofbird...
            ' ' wiener_entropy.parameters.syllables(i) '?\n(format [start stop])  ']);
    end
end

%% Calculating entropy during window

fprintf('\nCalculating entropy during window...')

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec)
        try
            
            %entropy only during window of entropy calculation
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).window{j} = ...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).all_time{j}(:,...
                wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j}>...
                min(wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i))) &...
                wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).time{j}<...
                max(wiener_entropy.parameters.time_range.(wiener_entropy.parameters.syllables(i))));
           
            
            %mean entropy during window (column 2) and time of day (column
            %1)
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j} = [...
                wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).date{j}...
                mean(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).window{j},2)];
            
            %mean entropy per day
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean_day{j} = mean(...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j}(:,2),1);

            %sd entropy per day
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).sd_day{j} = std(...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j}(:,2),1);
            
            
            
        catch err
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).window{j} = {};
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j} = {};
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean_day{j} = {};
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).sd_day{j} = {};
        end
    end
end
    
 fprintf('done!\n')
 
 %% Calculating baseline entropy

fprintf('\nCalculating baseline entropy...')
for i = 1:length(wiener_entropy.parameters.syllables)
    %calcualtes baseline frequency (kHz), sd (kHz), and cv
    try
        for j = 1:wiener_entropy.parameters.duration.begin
            temp_base.mean{j} = wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean_day{j};
        end
        temp_base.mean = temp_base.mean(~cellfun('isempty',temp_base.mean));
        temp_base.mean = cell2mat(temp_base.mean);
        wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).baseline_mean = mean(temp_base.mean);
    catch err
        continue
    end
    
    clear('temp_base')
    
end
fprintf('done!\n')

%% Figure of entropy values over time

%Determining the color for the various conditions (grey = preWN, red =
%shift, blue = consolidation, black = postWN)
fprintf('\nDetermining color of each day according to epoch during paradigm...')
for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}+1
        if j <= wiener_entropy.parameters.duration.begin
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j} = [.6 .6 .6];
        elseif j > wiener_entropy.parameters.duration.begin && j <= wiener_entropy.parameters.duration.begin+wiener_entropy.parameters.duration.shift
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j} = [1 0 0];
        elseif j > wiener_entropy.parameters.duration.begin+wiener_entropy.parameters.duration.shift &&...
                j <= wiener_entropy.parameters.duration.begin+wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j} = [0 0 1];
        elseif j > wiener_entropy.parameters.duration.begin+wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j} = [.3 .3 .3];
        end
    end
end
fprintf('done!\n')

fprintf('\nCreating figures of entropy over time...')
for i = 1:length(wiener_entropy.parameters.syllables)
    figure(), hold on
    title(['Wiener entropy for ' wiener_entropy.parameters.nameofbird '   duration: '...
        num2str(wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}+1)...
        ' days   syllable: ' wiener_entropy.parameters.syllables(i)])
    ylabel('Wiener Entropy')
    xlabel('Time (days)')
    xlim([-1-wiener_entropy.parameters.duration.begin...
        length(wiener_entropy.spectrogram.(wiener_entropy.parameters.syllables(i)).spec)+1])
    for j = wiener_entropy.parameters.days{3}:wiener_entropy.parameters.days{4}
        try

            %plots WE for time range
            plot(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j-wiener_entropy.parameters.days{3}+1}(:,1)-...
                min(floor(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{1}(:,1)))-...
                wiener_entropy.parameters.duration.begin,...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j-wiener_entropy.parameters.days{3}+1}(:,2),...
                '.',...
                'Color', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j-wiener_entropy.parameters.days{3}+1})
            
            %plots WE mean per day
            plot(median(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j-wiener_entropy.parameters.days{3}+1}(:,1))-...
                min(floor(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{1}(:,1)))-...
                wiener_entropy.parameters.duration.begin,...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean_day{j-wiener_entropy.parameters.days{3}+1},...
                'ko',...
                'MarkerFaceColor','k',...
                'LineWidth',2)
            
            %plots entropy sd per day
            errorbar(median(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j-wiener_entropy.parameters.days{3}+1}(:,1))-...
                min(floor(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{1}(:,1)))-...
                wiener_entropy.parameters.duration.begin,...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean_day{j-wiener_entropy.parameters.days{3}+1},...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).sd_day{j-wiener_entropy.parameters.days{3}+1},...
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
            [min(get(gca, 'YLim')) 0],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for consolidation begin
    try
        line([wiener_entropy.parameters.duration.shift...
            wiener_entropy.parameters.duration.shift],...
            [min(get(gca, 'YLim')) 0],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for WN end
    try
        line([wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation...
            wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation],...
            [min(get(gca, 'YLim')) 0],...
            'Color', [0 0 0],...
            'LineStyle', '--')
    catch err
        continue
    end
    %line for baseline entropy
    try
        line(get(gca, 'XLim'),...
            [wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).baseline_mean...
            wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).baseline_mean],...
            'Color', [.3 .3 .3],...
            'LineStyle', '--')
    catch err
        continue
    end
    
    
    
    saveas(figure(gcf), [wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type '/'...
        wiener_entropy.parameters.nameofbird...
        '/WE_over_time_' wiener_entropy.parameters.nameofbird '_' 'duration_'...
         num2str(wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}+1)...
        'days' '_' wiener_entropy.parameters.syllables(i)], 'fig')
end

fprintf('done!\n')

%% Save data

clear i
clear j

fprintf('\nSaving data...')

save([wiener_entropy.parameters.computer '/all_wiener_entropy/' wiener_entropy.parameters.experiment_type '/'...
        wiener_entropy.parameters.nameofbird '/' wiener_entropy.parameters.nameofbird '.mat'], 'wiener_entropy', '-v7.3')

fprintf('done!\n')

%% Do want to run db_WE_bootstrap

run_WE = input('Do you want to run WE_mult_bootstrap?  (y or n)  ', 's');

if strcmpi(run_WE,'y') == 1
    close all
    clear all
    db_WE_mult_bootstrap
else
    clear run_WE
    return
end
