%% To check sequence/entropy/etc. template matching, hit rate, and trigger timing

% This just looks at trigger timing and hit rate (both triggered (offline) 
% and all labeled (online)) with your current settings. Use for any
% paradigm that is not FF (use db_check_template_timing instead).
%
% Written by DM Brady 20/06/2012

%% Example of everything you need to do before you run this program

% transfer files from acquistion computer

% make a batch file of all cbins
% !ls *bk49_27*cbin > batch

% run db_catch_and_clean (input batch file name)
% will create a batch file of 'catch' trials and get rid of suspected noise
% files

% label syllables using evsonganaly
% evsonganaly
% (use filtfilt and load your batch.catch.keep file)

%% Asks if you loaded your template
% will stop program if you have not loaded your template

which_computer = input('Bluejay number? ');
computer = ['/bluejay' num2str(which_computer) '/dbrady/'];
clear which_computer


bird_name = input('What is the name of your bird? ', 's');
date_of_exp = input('what is the date of your experiment?\n(ex: 03Aug2012)  ', 's');

% A date stamp that will label all your saved variables and figures with the
% time that you did the analysis
today_date = datestr(now, 'ddmmmyyyy_HHMMAM');
today_date = today_date(today_date ~= ' ');

% Make a directory to store figures and variables
mkdir(['check_template_timing_' today_date])

% which batch file do want to analyze (probably batch.catch.keep)
batchfile = input('what is the name of your batch file?  ', 's');

%% Load previous parameters or not

load_previous_parameters = input('Would you like to load your last set of parameters? (y or n)  ', 's');
if strcmpi(load_previous_parameters, 'y')
    num_syl = input('how many syllables or transition points?  ');
    tpofinterest = cell(1,num_syl);
    for i = 1:num_syl
        tpofinterest{i} = input(['what is the name of syllable/transition point #' num2str(i) '?  '], 's');
        load([computer bird_name '/current_parameters_' tpofinterest{i} '.mat'])
        syllables.(tpofinterest{i}).parameters = current_parameters.(tpofinterest{i});
    end
    
elseif strcmpi(load_previous_parameters, 'n')
    did_load_template = input('Did you already load your template(s)? (y or n)  ','s');
    if strcmpi(did_load_template,'y')
    else
        num_temps = input('How many templates do you want to load?  ');
        template_name = cell(1, num_temps);
        for i = 1:num_temps
            template_name{i} = input(['What is the name of template #' num2str(i) '?  '], 's');
            try
                load([computer bird_name '/' template_name{i} '.dat'])
            catch err
                template_name{i} = input(['Sorry, what is the name of template #' num2str(i) ' again?  '], 's');
                load([computer bird_name '/' template_name{i} '.dat'])
            end
        end
    end
    
    %% Setting up batch file, template, and syllables
    
    
    
    % Input variables
    numberoftranisitionpoints = input('how many branch points do you have? ');
    
    
    tpofinterest = cell(1,numberoftranisitionpoints);
    syllables = struct;
    
    for i = 1:numberoftranisitionpoints
        tpofinterest{i} = input(['what is syllable/transition #' num2str(i) '?  '], 's');
        syllables.(tpofinterest{i}) = struct;
        if strcmpi(did_load_template,'y')
            syllables.(tpofinterest{i}).parameters.template = input(['what is the name of your template for ' tpofinterest{i} '?  ']);
        else
            syllables.(tpofinterest{i}).parameters.template = eval(template_name{i});
        end
        syllables.(tpofinterest{i}).parameters.numberoftemps = size(syllables.(tpofinterest{i}).parameters.template,2);
        syllables.(tpofinterest{i}).parameters.refractory_time = input(['what is the refractory time for ' tpofinterest{i} '?  ']);
        % syllables.(tpofinterest{i}).parameters.num_diff_cntrng = input(['how many different sets of cntrng values for ' tpofinterest{i} '?  ']);
        syllables.(tpofinterest{i}).parameters.same_cntrng = {};
        for m = 1:syllables.(tpofinterest{i}).parameters.numberoftemps
            if m == 1
                
                syllables.(tpofinterest{i}).parameters.cntrng_values{m} =...
                    input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
                    tpofinterest{i} ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
                
            elseif m > 1
                
                syllables.(tpofinterest{i}).parameters.same_cntrng{m-1} =...
                    input(['are the values for ' tpofinterest{i} ' in column ' num2str(m) ' the same as ' num2str(m-1) '? (y or n)  '],'s');
                
                if strcmpi(syllables.(tpofinterest{i}).parameters.same_cntrng{m-1}, 'y')
                    syllables.(tpofinterest{i}).parameters.cntrng_values{m} = syllables.(tpofinterest{i}).parameters.cntrng_values{m-1};
                    
                elseif strcmpi(syllables.(tpofinterest{i}).parameters.same_cntrng{m-1}, 'n')
                    syllables.(tpofinterest{i}).parameters.cntrng_values{m} =...
                        input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
                        tpofinterest{i} ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
                end
            end
        end
    end
end

%% saves current parameters

for i = 1:length(tpofinterest)
    current_parameters.(tpofinterest{i}) = syllables.(tpofinterest{i}).parameters;
    save([computer bird_name '/current_parameters_' tpofinterest{i}], 'current_parameters') 
end

%% Make X.tmp files and set cntrng

for i = 1:length(tpofinterest)
    %Makes X.tmp files
    mk_tempf(batchfile, syllables.(tpofinterest{i}).parameters.template ,2,'obs0');


%     %asks if you are using birdtaf
%     %syllables.(tpofinterest{i}).parameters.run_bt = input(['Birdtaf for ' tpofinterest{i} '? (y or n)  '], 's');
%     syllables.(tpofinterest{i}).parameters.run_bt = 'y'; %automatically asks for birdtaf templates
%     
%     %if you are not using birdtaf, the cntrng values are all the same,
%     %comment this section out and use the ones below if this is not true.
%     if strcmpi(syllables.(tpofinterest{i}).parameters.run_bt,'n')
        for j = 1:syllables.(tpofinterest{i}).parameters.numberoftemps
            cntrng(j).MIN=syllables.(tpofinterest{i}).parameters.cntrng_values{j}{1}(1);
            cntrng(j).MAX=syllables.(tpofinterest{i}).parameters.cntrng_values{j}{1}(2);
            cntrng(j).TH=syllables.(tpofinterest{i}).parameters.cntrng_values{j}{1}(3);
            
            if strcmpi(syllables.(tpofinterest{i}).parameters.cntrng_values{j}{2}, 'and')
                cntrng(j).AND=1;
            else
                cntrng(j).AND=0;
            end
            
            if strcmpi(syllables.(tpofinterest{i}).parameters.cntrng_values{j}{3}, 'y')
                cntrng(j).NOT=1;
            else
                cntrng(j).NOT=0;
            end
            
            if strcmpi(syllables.(tpofinterest{i}).parameters.cntrng_values{j}{4}, 'y')
                cntrng(j).MODE=0;
            else
                cntrng(j).MODE=1;
            end
            
            cntrng(j).BTMIN=0;
        end
    
%     %if you are using birdtaf, it will ask how many birdtaf columns are you using    
%     elseif strcmpi(syllables.(tpofinterest{i}).parameters.run_bt,'y')
%         syllables.(tpofinterest{i}).parameters.num_birdtaf = input(['How many Birdtaf templates for ' tpofinterest{i} '?  ']);
%         
%         %if you are using more than one birdtaf template, is it 'and' or
%         %'or' logic
%         if syllables.(tpofinterest{i}).parameters.num_birdtaf > 1
%             and_or = input('"And" or "Or" for birdtaf templates?', 's');
%         else
%             and_or = 'and';
%         end
%         
%         syllables.(tpofinterest{i}).parameters.birdtaf_values = cell(1,syllables.(tpofinterest{i}).parameters.numberoftemps); %for storing birdtaf cntrng values
%         
%         for kk = 1:syllables.(tpofinterest{i}).parameters.numberoftemps
%             
%             %I am assuming that you are not using 'not' logiv and BTMIN is
%             %0, please change if this is not true.
%             cntrng(kk).NOT=0;
%             cntrng(kk).BTMIN=0;
%             
%             %sets the cntrng values for your non-birdtaf templates
%             if kk <= syllables.(tpofinterest{i}).parameters.numberoftemps-syllables.(tpofinterest{i}).parameters.num_birdtaf
%                 cntrng(kk).MODE=1;
%                 cntrng(kk).MIN=syllables.(tpofinterest{i}).parameters.cntrng_values{kk}(1);
%                 cntrng(kk).MAX=syllables.(tpofinterest{i}).parameters.cntrng_values{kk}(2);
%                 cntrng(kk).TH=syllables.(tpofinterest{i}).parameters.cntrng_values{kk}(3);
%                 
%                 %sets last non-birdtaf template to 'and' logic, all others
%                 %are set to 'or' logic
%                 if kk == syllables.(tpofinterest{i}).parameters.numberoftemps-syllables.(tpofinterest{i}).parameters.num_birdtaf
%                     cntrng(kk).AND=1;
%                 else
%                     cntrng(kk).AND=0;
%                 end
%             
%             %sets the cntrng values for birdtaf templates
%             elseif kk > syllables.(tpofinterest{i}).parameters.numberoftemps-syllables.(tpofinterest{i}).parameters.num_birdtaf
%                 cntrng(kk).MODE=0; %this means it is birdtaf mode
%                 
%                 %asks for each birdtaf template your min, max, and
%                 %threshold values, then sets them.
%                 syllables.(tpofinterest{i}).parameters.birdtaf_values{kk} = input(['what are the min, max, and threshold for birdtaf in column ' num2str(kk) '?\n ([min max threshold])  ']);
%                 cntrng(kk).MIN=syllables.(tpofinterest{i}).parameters.birdtaf_values{kk}(1);
%                 cntrng(kk).MAX=syllables.(tpofinterest{i}).parameters.birdtaf_values{kk}(2);
%                 cntrng(kk).TH=syllables.(tpofinterest{i}).parameters.birdtaf_values{kk}(3);
%                 
%                 %sets birdtaf template to 'and' or 'or'
%                 if strcmpi(and_or,'and')
%                     cntrng(kk).AND=1;
%                 elseif strcmpi(and_or,'or')
%                     cntrng(kk).AND=0;
%                 end
%             end
%         end
%     end
            

            
            

            
    % otherwise use below
    % cntrng(1).MIN=3;
    % cntrng(1).MAX=4;
    % cntrng(1).NOT=0;
    % cntrng(1).MODE=1;
    % cntrng(1).TH=2.7;
    % cntrng(1).AND=0;
    % cntrng(1).BTMIN=0
    %
    % cntrng(2).MIN=3;
    % cntrng(2).MAX=4;
    % cntrng(2).NOT=0;
    % cntrng(2).MODE=1;
    % cntrng(2).TH=2.7;
    % cntrng(2).AND=0;
    % cntrng(2).BTMIN=0
    %
    % cntrng(3).MIN=3;
    % cntrng(3).MAX=4;
    % cntrng(3).NOT=0;
    % cntrng(3).MODE=1;
    % cntrng(3).TH=2.7;
    % cntrng(3).AND=0;
    % cntrng(3).BTMIN=0
    
%% Get triggers
    
    get_trigt2(batchfile,cntrng,syllables.(tpofinterest{i}).parameters.refractory_time,128,1,1)

%% Calculating hit rate (using ADDX = 0), offline
    
    %Finds triggers of interest for syllables that set off WN)
    [syllables.(tpofinterest{i}).hit.values syllables.(tpofinterest{i}).hit.trigger] = triglabel(batchfile,tpofinterest{i}(end),1,1,0,0);
    
    %Calculates hit rate
    syllables.(tpofinterest{i}).hit.sum = sum(syllables.(tpofinterest{i}).hit.values);
    syllables.(tpofinterest{i}).hit.hit_rate = [num2str(syllables.(tpofinterest{i}).hit.sum(1)...
        ./syllables.(tpofinterest{i}).hit.sum(2)*100) '%'];
    
    %Calculates trigger timing
    syllables.(tpofinterest{i}).hit.toff = [];
    for ii = 1:length(syllables.(tpofinterest{i}).hit.trigger)
        syllables.(tpofinterest{i}).hit.toff = [syllables.(tpofinterest{i}).hit.toff; syllables.(tpofinterest{i}).hit.trigger(ii).toffset];
    end
        
    syllables.(tpofinterest{i}).hit.toff_mean = mean(syllables.(tpofinterest{i}).hit.toff);
    syllables.(tpofinterest{i}).hit.toff_sd = std(syllables.(tpofinterest{i}).hit.toff);
    syllables.(tpofinterest{i}).hit.toff_median = median(syllables.(tpofinterest{i}).hit.toff);
    syllables.(tpofinterest{i}).hit.toff_iqr = iqr(syllables.(tpofinterest{i}).hit.toff);
    
    %To see if trigger timing is normally distributed. Stops program if we
    %catch err for some reason
%     try
%         if kstest(syllables.(sofinterest{i}).hit.toff) == 1;
%             syllables.(sofinterest{i}).hit.toff_normal = 'no';
%         elseif kstest(syllables.(sofinterest{i}).hit.toff) == 0;
%             syllables.(sofinterest{i}).hit.toff_normal = 'yes';
%         end
%     catch err
%         continue
%     end

    %Plots trigger timing vs. trigger number so you can see outliers    
    figure()
    plot(syllables.(tpofinterest{i}).hit.toff,...
        'o',...
        'Color', [1-(length(tpofinterest)-i)/length(tpofinterest) 0 0+(length(tpofinterest)-i)/length(tpofinterest)]);
    title(['Distribution of timing offsets for ' tpofinterest{i} ' (triggered)'])
    xlabel('Syllable number')
    ylabel('Timing offset (ms)')
    saveas(figure(gcf), ['check_template_timing_' today_date '/Toff_distro_' tpofinterest{i} '_(triggered)_' date_of_exp], 'fig')
    
    syllables.(tpofinterest{i}).hit.date = datestr(now);

    
    
%% Calculating template matching and timing (using ADDX = 1), online
    
    %Finds triggers for all labeled syllables 
    [syllables.(tpofinterest{i}).match.values syllables.(tpofinterest{i}).match.trigger] = triglabel(batchfile,tpofinterest{i}(end),1,1,0,1);
    
    %Calculates match rate
    syllables.(tpofinterest{i}).match.sum = sum(syllables.(tpofinterest{i}).match.values);
    syllables.(tpofinterest{i}).match.match_rate = [num2str(syllables.(tpofinterest{i}).match.sum(1)...
        ./syllables.(tpofinterest{i}).match.sum(2)*100) '%'];
    
    %Calculates timing of template matching
    syllables.(tpofinterest{i}).match.toff = [];
    for jj = 1:length(syllables.(tpofinterest{i}).match.trigger)
        syllables.(tpofinterest{i}).match.toff = [syllables.(tpofinterest{i}).match.toff; syllables.(tpofinterest{i}).match.trigger(jj).toffset];
    end
        
    syllables.(tpofinterest{i}).match.toff_mean = mean(syllables.(tpofinterest{i}).match.toff);
    syllables.(tpofinterest{i}).match.toff_sd = std(syllables.(tpofinterest{i}).match.toff);
    syllables.(tpofinterest{i}).match.toff_median = median(syllables.(tpofinterest{i}).match.toff);
    syllables.(tpofinterest{i}).match.toff_iqr = iqr(syllables.(tpofinterest{i}).match.toff);
    
    %To see if trigger timing is normally distributed. Stops program if we
    %catch err for some reason
%     try
%         if kstest(syllables.(sofinterest{i}).match.toff) == 1;
%             syllables.(sofinterest{i}).match.toff_normal = 'no';
%         elseif kstest(syllables.(sofinterest{i}).match.toff) == 0;
%             syllables.(sofinterest{i}).match.toff_normal = 'yes';
%         end
%     catch err
%         continue
%     end
    
    %Plots trigger timing vs. trigger number so you can see outliers
    figure()
    plot(syllables.(tpofinterest{i}).match.toff,...
        'o',...
        'Color', [1-(length(tpofinterest)-i)/length(tpofinterest) 0 0+(length(tpofinterest)-i)/length(tpofinterest)]);
    title(['Distribution of timing offsets for ' tpofinterest{i} ' (all)'])
    xlabel('Syllable number')
    ylabel('Timing offset (ms)')
    saveas(figure(gcf), ['check_template_timing_' today_date '/Toff_distro_' tpofinterest{i} '_(all)_' date_of_exp], 'fig')
    
    syllables.(tpofinterest{i}).match.date = datestr(now);

    
    clear cntrng
end


%% Displaying information

%Summary information about hit rate
for i = 1:length(tpofinterest)
    display(' ')
    display(['syllable/branch point: ' tpofinterest{i} ' hit rate information'])
    display(syllables.(tpofinterest{i}).hit)
    display(' ')
end

%Summary information about match rate
for i = 1:length(tpofinterest)
    display(' ')
    display(['syllable/branch point: ' tpofinterest{i} ' match rate information'])
    display(syllables.(tpofinterest{i}).match)
    display(' ')
end


%% Save Data




save(['check_template_timing_' today_date '/check_template_timing_' today_date])

%move batch files to new folder, make sure the batch file name is not too
%similar to the name of your bird
if strcmpi(input(['Move batch files to check_template_timing_' today_date ' folder? (y or n)  '], 's'),'y')
    movefile([batchfile(1:4) '*'],['check_template_timing_' today_date])
else
end

