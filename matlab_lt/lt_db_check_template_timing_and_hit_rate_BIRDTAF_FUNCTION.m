%% LT - 7/22 - do not use this. the non-function is most updated, then I rewrote in lt_check_hit_templ_freq

function [syllables] = lt_db_check_template_timing_and_hit_rate_BIRDTAF_FUNCTION(batch, num_syl, syllable, syllable_pre, syllable_post,birdtaf)
%% modified by LT 3/24 to become function
% changes: need to save parameters by running the script version of this function once on a
% specific day using the desired template. Then run this function at will.
% batch: 'batch.catch.keep'
% num_syl: 1
% syllable:{'d'}
% syllable_pre and _post: {'f'}
% birdtaf="y" (yes) or "n" (no)
%     (only handling one birdtaf template at the moment)

% output: syllables - structure

%% To check template matching, hit rate, and trigger timing

% This just looks at trigger timing and hit rate (both triggered (offline) 
% and all labeled (online)) with your current settings, 
% then calculates online FF to give you your new threshold.
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

%% Setting up batch file, template, and syllables

[check_pitch.parameters.bird_name, which_computer, ~] = lt_get_birdname_date_from_dir(1);
% which_computer = input('Bluejay number? ','s');
check_pitch.parameters.computer = ['/bluejay' which_computer '/lucas/birds/'];
%clear which_computer

% check_pitch.parameters.bird_name = input('What is the name of your bird? ', 's');

% A date stamp that will label all your saved variables and figures with the
% time that you did the analysis
check_pitch.parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM');
check_pitch.parameters.today_date = check_pitch.parameters.today_date(check_pitch.parameters.today_date ~= ' ');

% Make a directory to store figures and variables
mkdir(['check_template_' check_pitch.parameters.today_date])

% which batch file do want to analyze (probably batch.catch.keep)
% if strcmp(input('use batch.catch.keep? ','s'),'y')==1;
%     check_pitch.parameters.batchfile='batch.catch.keep';
% else
%     check_pitch.parameters.batchfile = input('what is the name of your batch file?  ', 's');
% end
check_pitch.parameters.batchfile = batch;


%% Load previous db_check_template_timing_and_hit_rate parameters

% load_previous_parameters = input('Would you like to load your last set of parameters? (y or n)  ', 's');
% if strcmpi(load_previous_parameters, 'y')
% num_syl = input('how many syllables?  ');
check_pitch.parameters.sofinterest = cell(1,num_syl);
for i = 1:num_syl
%     check_pitch.parameters.sofinterest{i} = input(['what is the name of syllable #' num2str(i) '?  '], 's');
%     check_pitch.parameters.sofinterest_pre{i} = input('what is the preceding syllable? - RETURN if dont care (for birdtaf) ','s');
%     check_pitch.parameters.sofinterest_post{i} = input('what is the subsequent syllable? - RETURN if dont care (for birdtaf) ','s');
%     load([check_pitch.parameters.computer check_pitch.parameters.bird_name '/current_parameters_' check_pitch.parameters.sofinterest{i} '.mat'])
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters = current_parameters.(check_pitch.parameters.sofinterest{i});
    
    check_pitch.parameters.sofinterest{i} = syllable{i};
    check_pitch.parameters.sofinterest_pre{i} = syllable_pre{i};
    check_pitch.parameters.sofinterest_post{i} = syllable_post{i};
    
    load([check_pitch.parameters.computer check_pitch.parameters.bird_name '/current_parameters_' check_pitch.parameters.sofinterest{i} '.mat'])
    syllables.(check_pitch.parameters.sofinterest{i}).parameters = current_parameters.(check_pitch.parameters.sofinterest{i});

end
% elseif strcmpi(load_previous_parameters, 'n')
%% Asks if you loaded your template
% will stop program if you have not loaded your template


%     did_load_template = input('Did you already load your template(s)? (y or n)  ', 's');
%     if strcmpi(did_load_template,'y')
%     else
%         num_temps = input('How many templates do you want to load?  ');
%         template_name = cell(1, num_temps);
%         for i = 1:num_temps
%             template_name{i} = input(['What is the name of template #' num2str(i) '?  '], 's');
%             load([check_pitch.parameters.computer check_pitch.parameters.bird_name '/' template_name{i} '.dat'])
%         end
%     end



%% Input variables
% numberofsyllables = size(syllable,2);

%     %note that this will only work if you want to load the same template. if
%     %you need to load an unusual combination of templates, act as if they are
%     %all different
%     try
%         if numberofsyllables > length(template_name)
%             for i = 2:numberofsyllables
%                 template_name{i} = template_name{1};
%             end
%         else
%         end
%     catch err
%         continue
%     end

% check_pitch.parameters.sofinterest = cell(1,numberofsyllables);

% for i = 1:numberofsyllables
%     check_pitch.parameters.sofinterest{i} = input(['what is syllable #' num2str(i) '?  '], 's');
%     syllables.(check_pitch.parameters.sofinterest{i}) = struct;
%     if strcmpi(did_load_template,'y')
%         syllables.(check_pitch.parameters.sofinterest{i}).parameters.template = input(['what is the name of your template for ' check_pitch.parameters.sofinterest{i} '?  ']);
%     else
%         syllables.(check_pitch.parameters.sofinterest{i}).parameters.template = eval(template_name{i});
%     end
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters.numberoftemps = size(syllables.(check_pitch.parameters.sofinterest{i}).parameters.template,2);
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters.refractory_time = input(['what is the refractory time for ' check_pitch.parameters.sofinterest{i} '?  ']);
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters.freq_range = input(['What is the frequency range for ' check_pitch.parameters.sofinterest{i} '?\n([min max])  ']);
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters.num_diff_cntrng = size(syllables.(check_pitch.parameters.sofinterest{i}).parameters.template,2);
%     syllables.(check_pitch.parameters.sofinterest{i}).parameters.same_cntrng = {};
%     
%     yes_or_no=input('use "template_parameters" structure for cntr parameters? (y or n): ', 's');
%     
%     if yes_or_no=='n';
%         for m = 1:syllables.(check_pitch.parameters.sofinterest{i}).parameters.num_diff_cntrng
%             if m == 1
%                 
%                 syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{m} =...
%                     input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
%                     check_pitch.parameters.sofinterest{i} ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
%                 
%             elseif m > 1
%                 
%                 syllables.(check_pitch.parameters.sofinterest{i}).parameters.same_cntrng{m-1} =...
%                     input(['are the values for ' check_pitch.parameters.sofinterest{i} ' in column ' num2str(m) ' the same as ' num2str(m-1) '? (y or n)  '],'s');
%                 
%                 if strcmpi(syllables.(check_pitch.parameters.sofinterest{i}).parameters.same_cntrng{m-1}, 'y')
%                     syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{m} = syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{m-1};
%                     
%                 elseif strcmpi(syllables.(check_pitch.parameters.sofinterest{i}).parameters.same_cntrng{m-1}, 'n')
%                     syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{m} =...
%                         input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
%                         check_pitch.parameters.sofinterest{i} ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
%                 end
%             end
%         end
%         
%     elseif yes_or_no=='y';
%         
%         load('/home/lucas/data/song/pu13bk43/template_parameters.mat')
%         syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values=eval(['template_parameters.' template_name{i}]);
%     end
%     
% end
% end



%% makes a variable to load recently used parameters
for i = 1:length(check_pitch.parameters.sofinterest)
    current_parameters.(check_pitch.parameters.sofinterest{i}) = syllables.(check_pitch.parameters.sofinterest{i}).parameters;
    save([check_pitch.parameters.computer check_pitch.parameters.bird_name '/current_parameters_' check_pitch.parameters.sofinterest{i}], 'current_parameters') 
end




%% Make X.tmp files and set cntrng
sofinterest=check_pitch.parameters.sofinterest; % just renaming.

for i = 1:length(check_pitch.parameters.sofinterest)
   %Makes X.tmp files
    mk_tempf(check_pitch.parameters.batchfile, syllables.(check_pitch.parameters.sofinterest{i}).parameters.template ,2,'obs0');


    %asks if you are using birdtaf
%     syllables.(sofinterest{i}).parameters.run_bt = input(['Birdtaf for syllable ' sofinterest{i} '? (y or n)  '], 's');
    syllables.(sofinterest{i}).parameters.run_bt = birdtaf;

    %if you are not using birdtaf, the cntrng values are all the same,
    %comment this section out and use the ones below if this is not true.
    if strcmpi(syllables.(sofinterest{i}).parameters.run_bt,'n')
        for j = 1:syllables.(check_pitch.parameters.sofinterest{i}).parameters.numberoftemps
            cntrng(j).MIN=syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{1}(1);
            cntrng(j).MAX=syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{1}(2);
            cntrng(j).TH=syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{1}(3);
            
            if strcmpi(syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{2}, 'and')
                cntrng(j).AND=1;
            else
                cntrng(j).AND=0;
            end
            
            if strcmpi(syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{3}, 'y')
                cntrng(j).NOT=1;
            else
                cntrng(j).NOT=0;
            end
            
            if strcmpi(syllables.(check_pitch.parameters.sofinterest{i}).parameters.cntrng_values{j}{4}, 'y')
                cntrng(j).MODE=0;
            else
                cntrng(j).MODE=1;
            end
            
            cntrng(j).BTMIN=0;
        end
    
%     %if you are using birdtaf, it will ask how many birdtaf columns are you using    
    elseif strcmpi(syllables.(sofinterest{i}).parameters.run_bt,'y')
%         syllables.(sofinterest{i}).parameters.num_birdtaf = input(['How many Birdtaf templates for syllable ' sofinterest{i} '?  ']);
         syllables.(sofinterest{i}).parameters.num_birdtaf = 1;
       
        %if you are using more than one birdtaf template, is it 'and' or
        %'or' logic
        if syllables.(sofinterest{i}).parameters.num_birdtaf > 1
            and_or = input('"And" or "Or" for birdtaf templates?', 's');
        else
            and_or = 'and';
        end
        
        syllables.(sofinterest{i}).parameters.birdtaf_values = cell(1,syllables.(sofinterest{i}).parameters.numberoftemps); %for storing birdtaf cntrng values
        
        for kk = 1:syllables.(sofinterest{i}).parameters.numberoftemps
            
            %I am assuming that you are not using 'not' logiv and BTMIN is
            %0, please change if this is not true.
            cntrng(kk).NOT=0;
            cntrng(kk).BTMIN=0;
            
            %sets the cntrng values for your non-birdtaf templates
            if kk <= syllables.(sofinterest{i}).parameters.numberoftemps-syllables.(sofinterest{i}).parameters.num_birdtaf
                cntrng(kk).MODE=1;
                cntrng(kk).MIN=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(1);
                cntrng(kk).MAX=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(2);
                cntrng(kk).TH=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(3);
                
                %sets last non-birdtaf template to 'and' logic, all others
                %are set to 'or' logic
                if kk == syllables.(sofinterest{i}).parameters.numberoftemps-syllables.(sofinterest{i}).parameters.num_birdtaf
                    cntrng(kk).AND=1;
                else
                    cntrng(kk).AND=0;
                end
            
            %sets the cntrng values for birdtaf templates
            elseif kk > syllables.(sofinterest{i}).parameters.numberoftemps-syllables.(sofinterest{i}).parameters.num_birdtaf
                cntrng(kk).MODE=0; %this means it is birdtaf mode
                
                %asks for each birdtaf template your min, max, and
                %threshold values, then sets them.
%                 syllables.(sofinterest{i}).parameters.birdtaf_values{kk} = input(['what are the min, max, and threshold for birdtaf in column ' num2str(kk) '?\n ([min max threshold])  ']);
                cntrng(kk).MIN=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(1);
                cntrng(kk).MAX=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(2);
                cntrng(kk).TH=syllables.(sofinterest{i}).parameters.cntrng_values{kk}{1}(3);
                
                %sets birdtaf template to 'and' or 'or'
                if strcmpi(and_or,'and')
                    cntrng(kk).AND=1;
                elseif strcmpi(and_or,'or')
                    cntrng(kk).AND=0;
                end
            end
        end
    end
            

            
            

            
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
    
    get_trigt2(check_pitch.parameters.batchfile,cntrng,syllables.(check_pitch.parameters.sofinterest{i}).parameters.refractory_time,128,1,1)

%% Calculating hit rate (using ADDX = 0), offline
%     
%     %Finds triggers of interest for syllables that set off WN)
%     [syllables.(check_pitch.parameters.sofinterest{i}).hit.values syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger] = triglabel(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},1,0,0,0);
%     
%     %Calculates hit rate
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.sum = sum(syllables.(check_pitch.parameters.sofinterest{i}).hit.values);
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.hit_rate = [num2str(syllables.(check_pitch.parameters.sofinterest{i}).hit.sum(1)...
%         ./syllables.(check_pitch.parameters.sofinterest{i}).hit.sum(2)*100) '%'];
%     
%     %Calculates trigger timing
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.toff = [];
%     for ii = 1:length(syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger)
%         syllables.(check_pitch.parameters.sofinterest{i}).hit.toff = [syllables.(check_pitch.parameters.sofinterest{i}).hit.toff; syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger(ii).toffset];
%     end
%         
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_mean = mean(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_sd = std(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_median = median(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_iqr = iqr(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
%     triglabel(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},1,1,0,0);
%     %To see if trigger timing is normally distributed. Stops program if we
%     %catch err for some reason
% %     try
% %         if kstest(syllables.(sofinterest{i}).hit.toff) == 1;
% %             syllables.(sofinterest{i}).hit.toff_normal = 'no';
% %         elseif kstest(syllables.(sofinterest{i}).hit.toff) == 0;
% %             syllables.(sofinterest{i}).hit.toff_normal = 'yes';
% %         end
% %     catch err
% %         continue
% %     end
% 
%     %Plots trigger timing vs. trigger number so you can see outliers    
%     figure()
%     plot(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff,...
%         'o',...
%         'Color', [1-(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest) 0 0+(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest)]);
%     title(['Distribution of timing offsets for ' check_pitch.parameters.sofinterest{i} ' (triggered)'])
%     xlabel('Syllable number')
%     ylabel('Timing offset (ms)')
%     saveas(figure(gcf), ['check_template_' check_pitch.parameters.today_date '/Toff_distro_' check_pitch.parameters.sofinterest{i} '_(triggered)_' check_pitch.parameters.today_date], 'fig')
%     
%     syllables.(check_pitch.parameters.sofinterest{i}).hit.date = datestr(now);

%% Calculating hit rate (using ADDX = 0), offline - LT MODIFIED 3/24/14 - from commented out section above.  Goal: used triglabel2 to get hitrate even with birdtaf.
    
    %Finds triggers of interest for syllables that set off WN)
    [syllables.(check_pitch.parameters.sofinterest{i}).hit.values, syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger] = triglabel2(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},...
        check_pitch.parameters.sofinterest_pre{i},check_pitch.parameters.sofinterest_post{i},1,0,0,0);
    
    %Calculates hit rate
    syllables.(check_pitch.parameters.sofinterest{i}).hit.sum = sum(syllables.(check_pitch.parameters.sofinterest{i}).hit.values);
    syllables.(check_pitch.parameters.sofinterest{i}).hit.hit_rate = [num2str(syllables.(check_pitch.parameters.sofinterest{i}).hit.sum(1)...
        ./syllables.(check_pitch.parameters.sofinterest{i}).hit.sum(2)*100) '%'];
    
    %Calculates trigger timing
    syllables.(check_pitch.parameters.sofinterest{i}).hit.toff = [];
    for ii = 1:length(syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger)
        syllables.(check_pitch.parameters.sofinterest{i}).hit.toff = [syllables.(check_pitch.parameters.sofinterest{i}).hit.toff; syllables.(check_pitch.parameters.sofinterest{i}).hit.trigger(ii).toffset];
    end
        
    syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_mean = mean(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
    syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_sd = std(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
    syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_median = median(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
    syllables.(check_pitch.parameters.sofinterest{i}).hit.toff_iqr = iqr(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff);
    triglabel(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},1,1,0,0);
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
    plot(syllables.(check_pitch.parameters.sofinterest{i}).hit.toff,...
        'o',...
        'Color', [1-(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest) 0 0+(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest)]);
    title(['Distribution of timing offsets for ' check_pitch.parameters.sofinterest{i} ' (triggered)'])
    xlabel('Syllable number')
    ylabel('Timing offset (ms)')
    saveas(figure(gcf), ['check_template_' check_pitch.parameters.today_date '/Toff_distro_' check_pitch.parameters.sofinterest{i} '_(triggered)_' check_pitch.parameters.today_date], 'fig')
    
    syllables.(check_pitch.parameters.sofinterest{i}).hit.date = datestr(now);
    
    
%% Calculating template matching and timing (using ADDX = 1), online
%     %Finds triggers for all labeled syllables 
%     [syllables.(check_pitch.parameters.sofinterest{i}).match.values, syllables.(check_pitch.parameters.sofinterest{i}).match.trigger] = ...
%         triglabel(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},1,0,0,1);
%     
%     %Calculates match rate
%     syllables.(check_pitch.parameters.sofinterest{i}).match.sum = sum(syllables.(check_pitch.parameters.sofinterest{i}).match.values);
%     syllables.(check_pitch.parameters.sofinterest{i}).match.match_rate = [num2str(syllables.(check_pitch.parameters.sofinterest{i}).match.sum(1)...
%         ./syllables.(check_pitch.parameters.sofinterest{i}).match.sum(2)*100) '%'];
%     
%     %Calculates timing of template matching
%     syllables.(check_pitch.parameters.sofinterest{i}).match.toff = [];
%     for jj = 1:length(syllables.(check_pitch.parameters.sofinterest{i}).match.trigger)
%         syllables.(check_pitch.parameters.sofinterest{i}).match.toff = [syllables.(check_pitch.parameters.sofinterest{i}).match.toff; syllables.(check_pitch.parameters.sofinterest{i}).match.trigger(jj).toffset];
%     end
%         
%     syllables.(check_pitch.parameters.sofinterest{i}).match.toff_mean = mean(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).match.toff_sd = std(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).match.toff_median = median(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
%     syllables.(check_pitch.parameters.sofinterest{i}).match.toff_iqr = iqr(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
%     
%     %To see if trigger timing is normally distributed. Stops program if we
%     %catch err for some reason
% %     try
% %         if kstest(syllables.(sofinterest{i}).match.toff) == 1;
% %             syllables.(sofinterest{i}).match.toff_normal = 'no';
% %         elseif kstest(syllables.(sofinterest{i}).match.toff) == 0;
% %             syllables.(sofinterest{i}).match.toff_normal = 'yes';
% %         end
% %     catch err
% %         continue
% %     end
%     
%     %Plots trigger timing vs. trigger number so you can see outliers
%     figure()
%     plot(syllables.(check_pitch.parameters.sofinterest{i}).match.toff,...
%         'o',...
%         'Color', [1-(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest) 0 0+(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest)]);
%     title(['Distribution of timing offsets for ' check_pitch.parameters.sofinterest{i} ' (all)'])
%     xlabel('Syllable number')
%     ylabel('Timing offset (ms)')
%     saveas(figure(gcf), ['check_template_' check_pitch.parameters.today_date '/Toff_distro_' check_pitch.parameters.sofinterest{i} '_(all)_' check_pitch.parameters.today_date], 'fig')
%     
%     syllables.(check_pitch.parameters.sofinterest{i}).match.date = datestr(now);

%% Calculating template matching and timing (using ADDX = 1), online - LT MODIFIED 3/24 to use triglabel2 (in the case of birdtaf)
%Finds triggers for all labeled syllables
[syllables.(check_pitch.parameters.sofinterest{i}).match.values, syllables.(check_pitch.parameters.sofinterest{i}).match.trigger] = ...
    triglabel2(check_pitch.parameters.batchfile,check_pitch.parameters.sofinterest{i},...
        check_pitch.parameters.sofinterest_pre{i},check_pitch.parameters.sofinterest_post{i},1,0,0,1);

%Calculates match rate
syllables.(check_pitch.parameters.sofinterest{i}).match.sum = sum(syllables.(check_pitch.parameters.sofinterest{i}).match.values);
syllables.(check_pitch.parameters.sofinterest{i}).match.match_rate = [num2str(syllables.(check_pitch.parameters.sofinterest{i}).match.sum(1)...
    ./syllables.(check_pitch.parameters.sofinterest{i}).match.sum(2)*100) '%'];

%Calculates timing of template matching
syllables.(check_pitch.parameters.sofinterest{i}).match.toff = [];
for jj = 1:length(syllables.(check_pitch.parameters.sofinterest{i}).match.trigger)
    syllables.(check_pitch.parameters.sofinterest{i}).match.toff = [syllables.(check_pitch.parameters.sofinterest{i}).match.toff; syllables.(check_pitch.parameters.sofinterest{i}).match.trigger(jj).toffset];
end

syllables.(check_pitch.parameters.sofinterest{i}).match.toff_mean = mean(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
syllables.(check_pitch.parameters.sofinterest{i}).match.toff_sd = std(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
syllables.(check_pitch.parameters.sofinterest{i}).match.toff_median = median(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);
syllables.(check_pitch.parameters.sofinterest{i}).match.toff_iqr = iqr(syllables.(check_pitch.parameters.sofinterest{i}).match.toff);

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
plot(syllables.(check_pitch.parameters.sofinterest{i}).match.toff,...
    'o',...
    'Color', [1-(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest) 0 0+(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest)]);
title(['Distribution of timing offsets for ' check_pitch.parameters.sofinterest{i} ' (all)'])
xlabel('Syllable number')
ylabel('Timing offset (ms)')
saveas(figure(gcf), ['check_template_' check_pitch.parameters.today_date '/Toff_distro_' check_pitch.parameters.sofinterest{i} '_(all)_' check_pitch.parameters.today_date], 'fig')

syllables.(check_pitch.parameters.sofinterest{i}).match.date = datestr(now);

    

%% Online FF calculation (ADDX = 1)
    
    threshold = [30 50 70 90]; %Percentiles for new threshold
    
    %Calculates FF for all labeled syllables
    syllables.(check_pitch.parameters.sofinterest{i}).freq.vals = evtaf_freq(check_pitch.parameters.batchfile, syllables.(check_pitch.parameters.sofinterest{i}).parameters.freq_range, check_pitch.parameters.sofinterest{i}, 128, 'obs0', 1);
    
    %Calculates running FF for all labeled syllables
    try
        running_avg_window = 10;
        syllables.(check_pitch.parameters.sofinterest{i}).freq.running_avg = db_runningaverage(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2),running_avg_window);
        figure
        grid on
        title('Running avg of FF')
        ylabel('FF')
        xlabel('Syllable number')
        hold on
        plot(syllables.(check_pitch.parameters.sofinterest{i}).freq.running_avg,'-','LineWidth',2,'MarkerSize',4);
        hold off
        saveas(figure(gcf),['check_template_' check_pitch.parameters.today_date '/Running_avg_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date], 'fig')
    catch err
        continue
    end
    
    %calculates cumulative FF for all labeled syllables
    for jj = 1:length(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2))
        syllables.(check_pitch.parameters.sofinterest{i}).freq.cumulative_avg(jj) =...
            mean(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(1:jj,2));
    end
    figure
    grid on
    title('Cumulative avg of FF')
    ylabel('FF')
    xlabel('Syllable number')
    hold on
    plot(syllables.(check_pitch.parameters.sofinterest{i}).freq.cumulative_avg,'-','LineWidth',2,'MarkerSize',4);
    hold off
    saveas(figure(gcf),['check_template_' check_pitch.parameters.today_date '/Cumulative_avg_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date], 'fig')
    
    %Summary statistics of FF for labeled syllables
    syllables.(check_pitch.parameters.sofinterest{i}).freq.mean = mean(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2));
    syllables.(check_pitch.parameters.sofinterest{i}).freq.sd = std(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2));
    syllables.(check_pitch.parameters.sofinterest{i}).freq.median = median(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2));
    syllables.(check_pitch.parameters.sofinterest{i}).freq.iqr = iqr(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2));
    
    syllables.(check_pitch.parameters.sofinterest{i}).freq.prctiles = num2str(threshold);
    syllables.(check_pitch.parameters.sofinterest{i}).freq.threshold = prctile(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2), threshold);
    
    
    %To see if FF  is normally distributed
    try
        if kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 1;
            syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'no';
        elseif kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 0;
            syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'yes';
        end
    catch err
        continue
    end
    
    %Histogram of FF distrobution
    figure()
    [f,x] = hist(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2),30);
    bar(x,f/trapz(x,f))
    title(['FF Probability density function for ' check_pitch.parameters.sofinterest{i}])
    xlabel('Frequency (Hz)')
    ylabel('Density')
    saveas(figure(gcf), ['check_template_' check_pitch.parameters.today_date '/FF_PDF_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date], 'fig')
        
    %Plots FF of labeled syllables
    figure(), hold on
    plot(24.*(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,1)-floor(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,1))),...
        syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2),...
        'o',...
        'Color', [1-(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest) 0 0+(length(check_pitch.parameters.sofinterest)-i)/length(check_pitch.parameters.sofinterest)]);
    title(['FF of ' check_pitch.parameters.sofinterest(i)])
    xlabel('Time (hours)')
    ylabel('Frequency (Hz)')
    hold off
    saveas(figure(gcf),['check_template_' check_pitch.parameters.today_date '/check_pitch_shift_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date],'fig')
        
        
    syllables.(check_pitch.parameters.sofinterest{i}).freq.date = datestr(now);
    
    clear cntrng

end



%% Displaying information

%To know the time range that you analyzed
timerange{1} = datestr(min(syllables.(check_pitch.parameters.sofinterest{1}).freq.vals(:,1)));
timerange{2} = datestr(max(syllables.(check_pitch.parameters.sofinterest{1}).freq.vals(:,1)));

display(timerange)

%Summary information about hit rate
for i = 1:length(check_pitch.parameters.sofinterest)
    display(' ')
    display(['syllable: ' check_pitch.parameters.sofinterest{i} ' hit rate information'])
    display(syllables.(check_pitch.parameters.sofinterest{i}).hit)
    display(' ')
end

%Summary information about match rate
for i = 1:length(check_pitch.parameters.sofinterest)
    display(' ')
    display(['syllable: ' check_pitch.parameters.sofinterest{i} ' match rate information'])
    display(syllables.(check_pitch.parameters.sofinterest{i}).match)
    display(' ')
end

%Summary information about FF
for i = 1:length(check_pitch.parameters.sofinterest)
    display(' ')
    display(['syllable: ' check_pitch.parameters.sofinterest{i} ' FF information'])
    display(syllables.(check_pitch.parameters.sofinterest{i}).freq)
    display(' ')
end

% QUICK SUMMARY
for i = 1:length(check_pitch.parameters.sofinterest)
    display(' ')
    display(['syllable: ' check_pitch.parameters.sofinterest{i} ' QUICK SUMMARY'])
    disp('# labeled (birdtaf)')
    disp(syllables.(check_pitch.parameters.sofinterest{i}).match.sum(2));
    disp('# matched with template (regardless whether hit)')
    disp(syllables.(check_pitch.parameters.sofinterest{i}).match.sum(1));
    disp('# hit with WN')
    disp(syllables.(check_pitch.parameters.sofinterest{i}).hit.sum(1));
    display(' ')
end


%% Save Data

save(['check_template_' check_pitch.parameters.today_date '/check_pitch_shift_' check_pitch.parameters.today_date])

%move batch files to new folder, make sure the batch file name is not too
%similar to the name of your bird

% if strcmpi(input(['Move batch files to check_template_timing_' check_pitch.parameters.today_date ' folder? (y or n)  '], 's'),'y')
%     movefile([check_pitch.parameters.batchfile(1:4) '*'],['check_template_' check_pitch.parameters.today_date])
% else
% end


















