%% modified by LT from db_contour_and_FF_analysis_over_time_v3
% 11/1/13, in order to automatically go through each song folder and
% perform various analyses, similarly to how pitch was done). - this is
% obsolete: use lt_all_days_various_calculations instead.

% also modified line 56 to add lesion dates to parameters

%% FF and entropy analysis over several days

%need to run db_transfer_calls(2) first


%% Entering or loading parameters for your bird

%gets current directory
multiple_pitch.parameters.curr_directory = pwd;
slashes = strfind(multiple_pitch.parameters.curr_directory,'/');

multiple_pitch.parameters.nameofbird = multiple_pitch.parameters.curr_directory(slashes(end)+1:end);

multiple_pitch.parameters.phrase = input('What is your phrase?  ','s');

if exist(['all_days_pitch_' multiple_pitch.parameters.phrase],'dir') == 0
    mkdir(['all_days_pitch_' multiple_pitch.parameters.phrase])
end

% tries to load structure with data from previous runs
try
    load([ 'all_days_pitch_' multiple_pitch.parameters.phrase '/' multiple_pitch.parameters.nameofbird '.mat'])
catch err
    display('Cannot find files, hope this is your first time running this program on this bird')
    display(' ')
    
    
    % asks for the syllables for the particular bird
    multiple_pitch.parameters.syllables =...
        input(['What are the syllables for ' multiple_pitch.parameters.nameofbird '?\n(no spaces between letters)  '], 's');
    
    multiple_pitch.parameters.days{1} =...
        input(['What is the first day for ' multiple_pitch.parameters.nameofbird '?\n(format should be day month year. i.e. 06Jun2012)    '], 's');
    multiple_pitch.parameters.days{3} = datenum(multiple_pitch.parameters.days{1});
    
    % Asks for the frequency range and refractory period for each syllable
    for i = 1:length(multiple_pitch.parameters.syllables)
        multiple_pitch.parameters.frequency_range.(multiple_pitch.parameters.syllables(i)) =...
            input(['Frequency range for ' multiple_pitch.parameters.syllables(i) ' for ' multiple_pitch.parameters.nameofbird '?\n(format [start stop])  ']);
        multiple_pitch.parameters.refractory_period.(multiple_pitch.parameters.syllables(i)) =...
            input(['Refractory period for ' multiple_pitch.parameters.syllables(i) ' for ' multiple_pitch.parameters.nameofbird '?  ']);
    end
    
    % Asks for length before WN delivery, pitch shift, consolidation, etc.
    multiple_pitch.parameters.WN_question=input('is this a WN experiment? (y or n)','s');
    if multiple_pitch.parameters.WN_question=='y';
        multiple_pitch.parameters.duration.begin = input('How many days before WN delivery?  ');
        multiple_pitch.parameters.duration.shift = input('How many days of driving learning?  ');
        multiple_pitch.parameters.duration.consolidation = input('How many days of consolidation?  ');
    end
    
    multiple_pitch.parameters.lesion_question=input('is this a lesion experiment? (y or n)','s');
    if multiple_pitch.parameters.lesion_question=='y';
        multiple_pitch.parameters.lesion_amount=input('how many lesions?');
        for i=1:multiple_pitch.parameters.lesion_amount;
            multiple_pitch.parameters.lesion_dates{i}=input(['what is the date of lesion # ' num2str(i) '? (e.g. 10Oct2013)'], 's');
            multiple_pitch.parameters.lesion_dates_datenum{i}=datenum(multiple_pitch.parameters.lesion_dates{i});
        end
        multiple_pitch.parameters.lesion_begin=1+multiple_pitch.parameters.lesion_dates_datenum{1}-multiple_pitch.parameters.days{3}; %plus 1 becuase last dates are inclusive
    end
end
    
% puts a time stamp on all figures and variables when you did the analysis
multiple_pitch.parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM');
multiple_pitch.parameters.today_date = multiple_pitch.parameters.today_date(multiple_pitch.parameters.today_date ~= ' ');


% adds the last date  for the bird. you must have already created batch
% files for the dates of interest.
multiple_pitch.parameters.days{2} =...
    input(['What is the last day for ' multiple_pitch.parameters.nameofbird '?\n(format should be day month year. i.e. 31Dec2013)   '], 's');
multiple_pitch.parameters.days{4} = datenum(multiple_pitch.parameters.days{2});

%checks to make sure you have at least two dayssong_folders = db_list_song_folders(multiple_pitch.parameters.phrase, omit);
if multiple_pitch.parameters.days{3} == multiple_pitch.parameters.days{4}
    display('Sorry, need to have at least two days before you can run this program')
    return
else
end




%% Making structures for results of pitch contour and FF analysis

fprintf('\nCalculating pitch contours...')

for i = 1:length(multiple_pitch.parameters.syllables)
    for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
        
        %if ~isempty checks to see if cell has been made
        %for this bird, syllable, and day
        try
            iscell(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1});
        catch err
            multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1} = {};
        end
        
        try
            iscell(multiple_pitch.fvalsstr_forpc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1});
        catch err
            multiple_pitch.fvalsstr_forpc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1} = {};
        end
        
        try
            iscell(multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1});
        catch err
            multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1} = {};
        end
        
    end
end


%% Calculating pitch contours


% this part takes a long time to process (specifically, the jc_pitchcontour
% analysis)



omit = input('Number of folders to omit (from end)?  ');

%finds the folder names to use
song_folders = db_list_song_folders(multiple_pitch.parameters.phrase, omit,...
    multiple_pitch.parameters.curr_directory);
fid = fopen([multiple_pitch.parameters.curr_directory ...
    '/song_folders_' multiple_pitch.parameters.phrase '.txt'],'r');
folders = textscan(fid,'%s');
fclose(fid);

%picks out only the folders within the dates specified
folders{1} = folders{1}(find(datenum(song_folders) == multiple_pitch.parameters.days{3}):find(datenum(song_folders) == multiple_pitch.parameters.days{4}));
    

for i = 1:length(multiple_pitch.parameters.syllables)
    k = 1;
    for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
        
        %if ~isempty checks to see if findnote2tw or jc_pitchcontour
        %has been performed on this bird, syllable, and day
        try
            if j == datenum(song_folders(k))
                cd([multiple_pitch.parameters.curr_directory '/' folders{1}{k}])
            end
            try
                if ~isempty(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1})
                    k = k+1;
                    continue
                else
                    [multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1}] =...
                        findwnote2tw('all_cbin_not_mat',...
                        multiple_pitch.parameters.syllables(i),'',...
                        multiple_pitch.parameters.refractory_period.(multiple_pitch.parameters.syllables(i)),...
                        multiple_pitch.parameters.frequency_range.(multiple_pitch.parameters.syllables(i)),1024,0,'obs0');
                    k = k+1;
                end
            catch err
                continue
            end
            
            try
                if ~isempty(multiple_pitch.fvalsstr_forpc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1})
                    continue
                else
                    [multiple_pitch.fvalsstr_forpc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1}] =...
                        findwnote2tw('all_cbin_not_mat',...
                        multiple_pitch.parameters.syllables(i),'',-0.016,...
                        multiple_pitch.parameters.frequency_range.(multiple_pitch.parameters.syllables(i)),8000,0,'obs0');
                end
            catch err
                continue
            end
            
            try
                if ~isempty(multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1})
                    continue
                else
                    multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1} = ...
                        jc_pitchcontourFV(multiple_pitch.fvalsstr_forpc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1},...
                        1024,1020,1, min(multiple_pitch.parameters.frequency_range.(multiple_pitch.parameters.syllables(i))),...
                        max(multiple_pitch.parameters.frequency_range.(multiple_pitch.parameters.syllables(i))),[1 2],'obs0');
                end
            catch err
                continue
            end
            
            cd(multiple_pitch.parameters.curr_directory)
        catch err
            display([ datestr(j, 'ddmmmyyyy') ' is missing'])
            continue
        end
    end
end

cd(multiple_pitch.parameters.curr_directory)

fprintf('done!\n')



%% Figure of pitch contour comparing mean (and sd) over time

fprintf('\nCreating figures of pitch contours over time...')

for i = 1:length(multiple_pitch.parameters.syllables)
    figure(i), hold on
    title([multiple_pitch.parameters.nameofbird '   duration: '...
        num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1) ' days'...
        '   syllable: ' multiple_pitch.parameters.syllables(i)])
    xlabel('Time (10^-4 sec)')
    ylabel('Frequency (Hz)')
    for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
        %plots the mean
        try
            plot(mean(multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1}'),...
                'Linewidth',2,...
                'Color',[1-(j-multiple_pitch.parameters.days{3})./(multiple_pitch.parameters.days{4} - multiple_pitch.parameters.days{3})...
                0 0+(j-multiple_pitch.parameters.days{3})./(multiple_pitch.parameters.days{4} - multiple_pitch.parameters.days{3})])
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
    saveas(figure(i), ['all_days_pitch_' multiple_pitch.parameters.phrase...
        '/PC_over_time_' multiple_pitch.parameters.nameofbird '_' 'duration_'...
        num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
        'days' '_' multiple_pitch.parameters.syllables(i)], 'fig')
end

fprintf('done!\n')

%% Variable for time range for calculating FF

% need to enter time range to calculate entropy, use the entropy contour figure to
% determine the best time range.

for i = 1:length(multiple_pitch.parameters.syllables)
    if isfield(multiple_pitch.parameters, 'time_range') == 1
        if isfield(multiple_pitch.parameters.time_range, multiple_pitch.parameters.syllables(i)) == 1
            if strcmpi(input(['Do you want to select a new time range for ' multiple_pitch.parameters.syllables(i) '?\n'...
                      'Current range:  [' num2str(multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i))) ']  ' ],'s'),'y')
                multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i)) =...
                    input(['Start and stop time for FF plotting of ' multiple_pitch.parameters.nameofbird...
                    ' ' multiple_pitch.parameters.syllables(i) '?\n(format [start stop])  ']);

                multiple_pitch.FF = rmfield(multiple_pitch.FF,multiple_pitch.parameters.syllables(i));
                
            else
            end
        else
           multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i)) =...
            input(['Start and stop time for FF plotting of ' multiple_pitch.parameters.nameofbird...
            ' ' multiple_pitch.parameters.syllables(i) '?\n(format [start stop])  ']);
        end
    else
        multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i)) =...
            input(['Start and stop time for FF plotting of ' multiple_pitch.parameters.nameofbird...
            ' ' multiple_pitch.parameters.syllables(i) '?\n(format [start stop])  ']);
    end
end

%% FF Plots over several days



% Calculates FF value, time during day, and gets rid of outliers
fprintf('\nCalculating FF over time and getting rid of outliers...')
for i = 1:length(multiple_pitch.parameters.syllables)
    for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
        try
            %calculates time of note during the day - DOES diff things if
            %this is WN versus lesion experiment. (but all have same
            %outcome of getting last baseline day to equal "0"
            if multiple_pitch.parameters.WN_question=='y';
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1) = ...
                    (db_timing4(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1})-...
                    floor(db_timing4(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){1}(1)))-...
                    multiple_pitch.parameters.duration.begin);
                
            elseif multiple_pitch.parameters.lesion_question=='y'
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1) = ...
                    (db_timing4(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1})-...
                    floor(db_timing4(multiple_pitch.fvalsstr_all.(multiple_pitch.parameters.syllables(i)){1}(1)))-...
                    multiple_pitch.parameters.lesion_begin);
            end
            %calculates FF of each note during day
            multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,2) = ...
                mean(multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){j-multiple_pitch.parameters.days{3}+1} ...
                (min(multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i))):...
                max(multiple_pitch.parameters.time_range.(multiple_pitch.parameters.syllables(i))),:))';
            
            %gets rid of outliers (Tukey's method)
            [multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}, ...
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).outliers_high{j-multiple_pitch.parameters.days{3}+1},...
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).outliers_low{j-multiple_pitch.parameters.days{3}+1}] = ...
                db_tukey_outlier(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1},2);
            
            %calculates mean for day
            multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1} = ...
                mean(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,2));
            
            %calculates sd for day
            multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).sd_FF{j-multiple_pitch.parameters.days{3}+1} = ...
                std(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,2));
            
        catch err
            % print(num2str(j))
            continue
        end
    end
end
fprintf('done!\n')

%% PLOTTING

% IF WN experiment:
if multiple_pitch.parameters.WN_question=='y';
    %Calculates baseline pitch
    fprintf('\nCalculating baseline FF...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        %calcualtes baseline frequency (kHz), sd (kHz), and cv
        try
            for j = 1:multiple_pitch.parameters.duration.begin
                temp_base.mean{j} = multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j};
            end
            temp_base.mean = temp_base.mean(~cellfun('isempty',temp_base.mean));
            temp_base.mean = cell2mat(temp_base.mean);
            multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean = mean(temp_base.mean);
        catch err
            continue
        end
        
        clear('temp_base')
        
    end
    
    fprintf('done!\n')%Determining the color for the various conditions (grey = preWN, red =
    %shift, blue = consolidation, black = postWN)
    fprintf('\nDetermining color of each day according to epoch during paradigm...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        for j = 1:multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1
            if j <= multiple_pitch.parameters.duration.begin
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j} = [.6 .6 .6];
            elseif j > multiple_pitch.parameters.duration.begin && j <= multiple_pitch.parameters.duration.begin+multiple_pitch.parameters.duration.shift
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j} = [1 0 0];
            elseif j > multiple_pitch.parameters.duration.begin+multiple_pitch.parameters.duration.shift &&...
                    j <= multiple_pitch.parameters.duration.begin+multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j} = [0 0 1];
            elseif j > multiple_pitch.parameters.duration.begin+multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j} = [.3 .3 .3];
            end
        end
    end
    fprintf('done!\n')
    
    % Makes a figure of FF for each day during learning
    fprintf('\nCreating figures of FF over time...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        figure(), hold on
        title([multiple_pitch.parameters.nameofbird '   duration: ' num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
            ' days   syllable: ' multiple_pitch.parameters.syllables(i)])
        ylabel('Frequency (Hz)')
        xlabel('Time (days)')
        if isempty(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}) == 0
            xlim([-1-multiple_pitch.parameters.duration.begin...
                ceil(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}(end,1))+1])
        else
            xlim([-1-multiple_pitch.parameters.duration.begin ...
                -1-multiple_pitch.parameters.duration.begin + length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)])
        end
        for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
            try
                
                %plots FF for time range
                plot(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,2),...
                    '.',...
                    'Color', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j-multiple_pitch.parameters.days{3}+1})
                
                %plots FF mean per day
                plot(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
                    'ko',...
                    'MarkerFaceColor','k',...
                    'LineWidth',2)
                
                %plots FF sd per day
                errorbar(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).sd_FF{j-multiple_pitch.parameters.days{3}+1},...
                    'ko',...
                    'LineWidth',2)
                
                
            catch err
                continue
            end
        end
        
        
        %line for end of baseline period
        try
            line([0-.25 0-.25],...
                get(gca, 'YLim'),...
                'Color', [0 0 0],...
                'LineStyle', '--')
        catch err
            continue
        end
        
        %line for consolidation begin
        try
            line([multiple_pitch.parameters.duration.shift-.25...
                multiple_pitch.parameters.duration.shift-.25],...
                get(gca, 'YLim'),...
                'Color', [0 0 0],...
                'LineStyle', '--')
        catch err
            continue
        end
        %line for WN end
        try
            line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.25...
                multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.25],...
                get(gca, 'YLim'),...
                'Color', [0 0 0],...
                'LineStyle', '--')
        catch err
            continue
        end
        %line for baseline pitch
        try
            line(get(gca, 'XLim'),...
                [multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean...
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean],...
                'Color', [.3 .3 .3],...
                'LineStyle', '--')
        catch err
            continue
        end
        
        %SAVE, each fig for each syllable
        saveas(figure(gcf), ['all_days_pitch_' multiple_pitch.parameters.phrase...
            '/FF_over_time_' multiple_pitch.parameters.nameofbird '_' 'duration_'...
            num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
            'days' '_' multiple_pitch.parameters.syllables(i)], 'fig')
    end
    
    
    % IF lesion experiment, then do the following to plot
elseif multiple_pitch.parameters.lesion_question=='y';
    
    %Calculates baseline pitch
    fprintf('\nCalculating baseline FF...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        %calcualtes baseline frequency (kHz), sd (kHz), and cv
        try
            for j = 1:multiple_pitch.parameters.lesion_begin
                temp_base.mean{j} = multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j};
            end
            temp_base.mean = temp_base.mean(~cellfun('isempty',temp_base.mean));
            temp_base.mean = cell2mat(temp_base.mean);
            multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean = mean(temp_base.mean);
        catch err
            continue
        end
        
        clear('temp_base')
        
    end
    fprintf('done!\n')
    
    % Makes a figure of FF for each day during learning
    fprintf('\nCreating figures of FF over time...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        figure(), hold on
        title([multiple_pitch.parameters.nameofbird '   duration: ' num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
            ' days   syllable: ' multiple_pitch.parameters.syllables(i)])
        ylabel('Frequency (Hz)')
        xlabel('Time (days)')
        if isempty(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}) == 0
            xlim([-1-multiple_pitch.parameters.lesion_begin...
                ceil(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}(end,1))+1])
        else
            xlim([-1-multiple_pitch.parameters.lesion_begin ...
                -1-multiple_pitch.parameters.lesion_begin + length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)])
        end
        for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
            try
                
                %plots FF for time range
                plot(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,2),...
                    '.')
                
                %plots FF mean per day
                plot(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
                    'ko',...
                    'MarkerFaceColor','k',...
                    'LineWidth',2)
                
                %plots FF sd per day
                errorbar(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
                    multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).sd_FF{j-multiple_pitch.parameters.days{3}+1},...
                    'ko',...
                    'LineWidth',2)
                
                
            catch err
                continue
            end
        end
        % line for baseline pitch
        try
            line(get(gca, 'XLim'),...
                [multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean...
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean],...
                'Color', [.3 .3 .3],...
                'LineStyle', '--')
        catch err
            continue
        end
        %line for end of baseline period
        try
            line([0-.25 0-.25],...
                get(gca, 'YLim'),...
                'Color', [0 0 0],...
                'LineStyle', '--')
        catch err
            continue
        end
        
        % LINES demarcating lesion days (all are on end of day because
        % generlaly have pre-lesion recordings for that same day, but not
        % post)
        for l=1:multiple_pitch.parameters.lesion_amount-1;
            line_index(l)=multiple_pitch.parameters.lesion_dates_datenum{l+1}-multiple_pitch.parameters.lesion_dates_datenum{1};
        end
        for i=1:length(line_index);
            try
                line([line_index(i)+0.25...
                    line_index(i)]+0.25,...
                    get(gca, 'YLim'),...
                    'Color', [0 0 0],...
                    'LineStyle', '--')
            catch err
                continue
            end
        end
        %SAVE, each fig for each syllable
        saveas(figure(gcf), ['all_days_pitch_' multiple_pitch.parameters.phrase...
            '/FF_over_time_' multiple_pitch.parameters.nameofbird '_' 'duration_'...
            num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
            'days' '_' multiple_pitch.parameters.syllables(i)], 'fig')
        
    end
    
    
end

fprintf('done!\n')

%% Save workspace

clear i
clear j

fprintf('\nSaving data...')

%saves just the FF data in your bird folder with the song folders
save(['all_days_pitch_' multiple_pitch.parameters.phrase '/' multiple_pitch.parameters.nameofbird '_' ...
    multiple_pitch.parameters.phrase '_FF_data.mat'], '-struct', 'multiple_pitch', 'FF')

save(['all_days_pitch_' multiple_pitch.parameters.phrase '/' multiple_pitch.parameters.nameofbird '.mat'], 'multiple_pitch', '-v7.3')

fprintf('done!\n')

%% Do want to run db_FF_mult_bootstrap

run_db_FF = input('Do you want to run db_FF_mult_bootstrap?  (y or n)  ', 's');

if strcmpi(run_db_FF,'y') == 1
    close all
    clear all
    lt_db_FF_mult_bootstrap_v2
else
    clear run_db_FF
    return
end


%% LT added 11/1/13 - calculating entropy, amplitude, sung fraction, using db_seq_func_day_save code



if input('do you want to run lt_db_seq_func_day_save_rd3gr35? (y or n) ', 's') =='y';
    
    plot_all_days_answer=input('do you want to plot results for all days? (y or n) ', 's'); % do you want to plot all days after this is done?

    % FIRST, set up loop to go through every day with song
    %finds the folder names to use
    omit = input('Number of folders to omit (from end)?  ');
    song_folders = db_list_song_folders(multiple_pitch.parameters.phrase, omit,...
        multiple_pitch.parameters.curr_directory);
    fid = fopen([multiple_pitch.parameters.curr_directory ...
        '/song_folders_' multiple_pitch.parameters.phrase '.txt'],'r');
    folders = textscan(fid,'%s');
    fclose(fid);
    
    %picks out only the folders within the dates specified
    folders{1} = folders{1}(find(datenum(song_folders) == multiple_pitch.parameters.days{3}):find(datenum(song_folders) == multiple_pitch.parameters.days{4}));
    
    
    for i = 1:length(multiple_pitch.parameters.syllables)
        k = 1;
        for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
            
            try
                if j == datenum(song_folders(k))
                    cd([multiple_pitch.parameters.curr_directory '/' folders{1}{k}])
                end
                % INSERT ANALYSIS STUFF FOR EACH DAY HERE
                
                % for each day analyze spectral features:
                lt_db_seq_func_day_save_rd3gr35('batch.catch.keep','',multiple_pitch.parameters.nameofbird, datestr(j,'ddmmmyyyy'), 1, multiple_pitch.parameters.syllables(i), multiple_pitch.parameters.phrase)
                cd(multiple_pitch.parameters.curr_directory)
                k=k+1;
                disp([ datestr(j, 'ddmmmyyyy') ' for ' multiple_pitch.parameters.syllables(i) ' done']);
            catch err
                display([ datestr(j, 'ddmmmyyyy') ' is missing'])
                continue
                % end spectral features
                
                
            end
        end
    end
    
    cd(multiple_pitch.parameters.curr_directory)
    
    fprintf('done!\n')
    
    if plot_all_days_answer=='y';
        display('plotting over all days')
        cd(['/' computer '/all_days_timing_' phrase])
        
        lt_write_all_folder_contents_to_batch;
        
        % to plot above data
        batch=input('what is name of batch? ','s');
        syllable=input('what is the name of the syllable? ','s');
        metric=input('what metric is this? (e.g. entropy, pitch, amplitude; all lower case)', 's');
        
        
        
        lt_db_plot_over_experiment(batch,'r', 1, 1,1,0) % with removal of tukey outliers
        
        title([metric ' vs. day; tukey outliers removed. Syllable: ' num2str(syllable)])
        % figure; lt_db_plot_over_experiment(batch,'g', 1, 0,1,0) % without removal of tukey outliers
        % title('entropy vs. day; tukey outliers not removed')
        
        
        saveas(figure(gcf),[metric '_vs_day_all_points_of_' num2str(syllable),'fig'])
    end
    
end

