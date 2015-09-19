%% LT 7/22 - directly modified from lt_db_check_template_timing_and_hit_rate_BIRDTAF.
% USE THIS NOW
% NOTE: might not work if try to do multiple syllables at once - run
% separates for each analysis.

% What this does:
% - Looks at template timing and hit rate (both online and offline), and gets
% a tonne of information on that.
% - Also calculates pitch at the location of the hit (time bin directly
% rpeceding match), offline,
% - Accomodates birdtaf templates, and isolates labels consistent with
% birdtaf.
% - Accomodates EvtafAmp and v4.
% - Does not support multiple bracketed logic (e.g. (a+b+c)*(d+e+f) for
% evtaf_v4 templates. Can modify code to do that see uievtaf4sim,
% findtrigtimes, and CounterSetup.

% OPTIONS: to analyze just:
% just hits - then don't need template info
% hits and matches
% hits and matches with freq

%% Before running:
% 1) Label songs
% 2) run lt_db_transfer_calls(2)
% 3) load template if doing offline checks.

% NOTES ON COL LOGIC:
% 1) needs to have brackets around everything - e.g. even if simple a+b+c, needs to be (a+b+c);
% 2) Might not handle recursiveness well (since i am not sure how actual evtaf handles it) --> e.g. (a+(a*b))*(c)
% 3) goes from left to right within brackets (i assume that is how real evtaf does it)
% 4) related to 3--> everything within set of brackets shoudl be same logic


%% PARAMETERS

% INPUTS:
evtaf_ver=input('which evtaf ver? (amp or v4) ?','s');
get_offline_match=1; % do offline matching using template? (ADDX=1)
get_FF=1; % Analyze FF using offline matching?
% which batch file do want to analyze (probably batch.catch.keep)
if strcmp(input('use batch.catch.keep? ','s'),'y')==1;
    parameters.batchfile='batch.catch.keep';
else
    parameters.batchfile = input('what is the name of your batch file?  ', 's');
end
[parameters.bird_name, which_computer, data_date, parameters.phrase] = lt_get_birdname_date_from_dir(1);
parameters.computer = ['/bluejay' which_computer '/lucas/birds/'];

% Syllable to analyze
parameters.sofinterest='b';
parameters.sofinterest_pre='a';
parameters.sofinterest_post='';
parameters.refractory_time = input(['what is the refractory time for ' parameters.sofinterest '?  ']);
parameters.freq_range = input(['What is the frequency range for ' parameters.sofinterest '?\n([min max])  ']);

% Template: (has to be saved in the bird folder)
if get_offline_match==1; % then need template information
    template_name = input(['What is the name of template ?  '], 's');
    load([parameters.computer parameters.bird_name '/' template_name '.dat'])
    % cntrng parameters
    parameters.template=eval(template_name);
    parameters.numberoftemps = size(parameters.template,2);
    parameters.num_diff_cntrng = size(parameters.template,2);
    
    %     for m = 1:parameters.num_diff_cntrng
    %         if m == 1
    %             parameters.cntrng_values{m} =...
    %                 input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
    %                 parameters.sofinterest ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
    %         elseif m > 1
    %             parameters.same_cntrng{m-1} =...
    %                 input(['are the values for ' parameters.sofinterest ' in column ' num2str(m) ' the same as ' num2str(m-1) '? (y or n)  '],'s');
    %             if strcmpi(parameters.same_cntrng{m-1}, 'y')
    %                 parameters.cntrng_values{m} = parameters.cntrng_values{m-1};
    %             elseif strcmpi(parameters.same_cntrng{m-1}, 'n')
    %                 parameters.cntrng_values{m} =...
    %                     input(['what are the cntrng min, max, threshold, "and/or",\n "not", and "birdtaf" for '...
    %                     parameters.sofinterest ' in column ' num2str(m) '?\n ex: {[min max threshold] "and" "y" "y"})  ']);
    %             end
    %         end
    %     end
    % end
    parameters.cntrng_values{1}={[5 14 1] 'or' 'n' 'y'};
    parameters.cntrng_values{2}={[5 14 0.1] 'or' 'n' 'y'};
    parameters.cntrng_values{3}={[5 14 1.1] 'and' 'n' 'y'};
    parameters.cntrng_values{4}={[1 4 1.5] 'or' 'n' 'n'};
    parameters.cntrng_values{5}={[1 4 1.5] 'or' 'n' 'n'};
    parameters.cntrng_values{6}={[1 4 1.5] 'or' 'n' 'n'};
    
    if strcmp(evtaf_ver,'v4');
        parameters.col_logic='(a+b+c)*d+e+f';
    end
    
end



% AUTOMATIC:
parameters.data_date=data_date{2};
parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM'); % time-stamp for analysis
parameters.today_date = parameters.today_date(parameters.today_date ~= ' ');

% Save dir
SaveDir1=[parameters.computer parameters.bird_name '/check_hit_templ_freq/'];
SaveDir2=[SaveDir1 parameters.phrase '_' parameters.data_date];
SaveDir3=[SaveDir2 '/Analysis_' parameters.today_date];
mkdir(SaveDir3);


%% FIRST, Calculating hit rate (using ADDX = 0), (using actual WN trig records)
sofinterest=parameters.sofinterest; % just renaming for ease of use

%Finds triggers of interest for syllables that set off WN)
[check_stuff.(sofinterest).hit.values, check_stuff.(sofinterest).hit.trigger] = triglabel2_LT(parameters.batchfile,parameters.sofinterest,...
    parameters.sofinterest_pre,parameters.sofinterest_post,1,0,0,0);

%Calculates hit rate
check_stuff.(sofinterest).hit.sum = sum(check_stuff.(sofinterest).hit.values);
check_stuff.(sofinterest).hit.hit_rate = [num2str(check_stuff.(sofinterest).hit.sum(1)...
    ./check_stuff.(sofinterest).hit.sum(2)*100) '%'];

%Calculates trigger timing
check_stuff.(sofinterest).hit.toff = [];
for ii = 1:length(check_stuff.(sofinterest).hit.trigger)
    check_stuff.(sofinterest).hit.toff = [check_stuff.(sofinterest).hit.toff; check_stuff.(sofinterest).hit.trigger(ii).toffset];
end

check_stuff.(sofinterest).hit.toff_mean = mean(check_stuff.(sofinterest).hit.toff);
check_stuff.(sofinterest).hit.toff_sd = std(check_stuff.(sofinterest).hit.toff);
check_stuff.(sofinterest).hit.toff_median = median(check_stuff.(sofinterest).hit.toff);
check_stuff.(sofinterest).hit.toff_iqr = iqr(check_stuff.(sofinterest).hit.toff);


%Plots trigger timing vs. trigger number so you can see outliers
figure()
plot(check_stuff.(sofinterest).hit.toff,'o');
title(['Distribution of timing offsets for ' parameters.sofinterest ' (triggered)'])
xlabel('Syllable number')
ylabel('Timing offset (ms)')
saveas(figure(gcf), [SaveDir3 '/Toff_distro_' parameters.sofinterest '_(correct triggered)_' parameters.today_date], 'fig')


%% SECOND(a) - OFFLINE TEMPLATE MATCHING (using ADDX = 1) - Make X.tmp files and set cntrng (i.e. use template to get offline matches)

if get_offline_match==1;
    %Makes X.tmp files
    mk_tempf(parameters.batchfile, parameters.template ,2,'obs0');
    
    % Converts cntrng value into format for get_trigs
    for j = 1:parameters.numberoftemps;
        cntrng(j).MIN=parameters.cntrng_values{j}{1}(1);
        cntrng(j).MAX=parameters.cntrng_values{j}{1}(2);
        cntrng(j).TH=parameters.cntrng_values{j}{1}(3);
        
        if strcmpi(parameters.cntrng_values{j}{2}, 'and')
            cntrng(j).AND=1;
        else
            cntrng(j).AND=0;
        end
        
        if strcmpi(parameters.cntrng_values{j}{3}, 'y')
            cntrng(j).NOT=1;
        else
            cntrng(j).NOT=0;
        end
        
        if strcmpi(parameters.cntrng_values{j}{4}, 'y')
            cntrng(j).MODE=0;
        else
            cntrng(j).MODE=1;
        end
        
        cntrng(j).BTMIN=0;
    end
    
    parameters.cntrng=cntrng;
    
    %% SECOND(b) Using the X.tmp file above, gets offline triggers
    
    if strcmp(evtaf_ver,'v4');
        get_trigt2_v4_LT(parameters.batchfile,cntrng,parameters.refractory_time,128,1,1,parameters.col_logic)
    elseif strcmp(evtaf_ver,'amp');
        get_trigt2(parameters.batchfile,cntrng,parameters.refractory_time,128,1,1)
    else
        disp('problem - what version of get_trig2 to run?');
    end
    
    
    %% SECOND(c), Calculating template matching
    %Finds triggers for all labeled syllables
    [check_stuff.(sofinterest).match.values, check_stuff.(sofinterest).match.trigger] = ...
        triglabel2_LT(parameters.batchfile,parameters.sofinterest,...
        parameters.sofinterest_pre,parameters.sofinterest_post,1,0,0,1);
    
    %Calculates match rate
    check_stuff.(sofinterest).match.sum = sum(check_stuff.(sofinterest).match.values);
    check_stuff.(sofinterest).match.match_rate = [num2str(check_stuff.(sofinterest).match.sum(1)...
        ./check_stuff.(sofinterest).match.sum(2)*100) '%'];
    
    %Calculates timing of template matching
    check_stuff.(sofinterest).match.toff = [];
    for jj = 1:length(check_stuff.(sofinterest).match.trigger)
        check_stuff.(sofinterest).match.toff = [check_stuff.(sofinterest).match.toff; check_stuff.(sofinterest).match.trigger(jj).toffset];
    end
    
    check_stuff.(sofinterest).match.toff_mean = mean(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_sd = std(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_median = median(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_iqr = iqr(check_stuff.(sofinterest).match.toff);
    
    
    %Plots trigger timing vs. trigger number so you can see outliers
    figure()
    plot(check_stuff.(sofinterest).match.toff,'o');
    title(['Distribution of timing offsets for ' parameters.sofinterest ' (Offline match; correct label)'])
    xlabel('Syllable number')
    ylabel('Timing offset (ms)')
    saveas(figure(gcf), [SaveDir3 '/Toff_distro_' parameters.sofinterest '_(all)_' parameters.today_date], 'fig')
    
    
    %% THIRD, FF calculation (Only using offline, in order to get even those that failed WN contingency - ADDX = 1)
    
    if get_FF==1;
        threshold = [30 50 70 90]; %Percentiles for new threshold
        
        %Calculates FF for all labeled syllables
        check_stuff.(sofinterest).freq.vals = evtaf_freq2_LT(parameters.batchfile, parameters.freq_range, parameters.sofinterest,...
            parameters.sofinterest_pre,parameters.sofinterest_post, 128, 'obs0', 1,evtaf_ver);
        
        %Calculates running FF for all labeled syllables
        try
            running_avg_window = 10;
            check_stuff.(sofinterest).freq.running_avg = db_runningaverage(check_stuff.(sofinterest).freq.vals(:,2),running_avg_window);
            figure
            grid on
            title('Running avg of FF')
            ylabel('FF')
            xlabel('Syllable number')
            hold on
            plot(check_stuff.(sofinterest).freq.running_avg,'-','LineWidth',2,'MarkerSize',4);
            saveas(figure(gcf),['check_template_' parameters.today_date '/Running_avg_' parameters.sofinterest{i} '_' parameters.today_date], 'fig')
        catch err
        end
        
        %calculates cumulative FF for all labeled syllables
        for jj = 1:length(check_stuff.(sofinterest).freq.vals(:,2))
            check_stuff.(sofinterest).freq.cumulative_avg(jj) = mean(check_stuff.(sofinterest).freq.vals(1:jj,2));
        end
        figure
        grid on
        title('Cumulative avg of FF')
        ylabel('FF')
        xlabel('Syllable number')
        hold on
        plot(check_stuff.(sofinterest).freq.cumulative_avg,'-','LineWidth',2,'MarkerSize',4);
        saveas(figure(gcf),[SaveDir3 '/Cumulative_avg_' parameters.sofinterest '_' parameters.today_date], 'fig')
        
        %Summary statistics of FF for labeled syllables
        check_stuff.(sofinterest).freq.mean = mean(check_stuff.(sofinterest).freq.vals(:,2));
        check_stuff.(sofinterest).freq.sd = std(check_stuff.(sofinterest).freq.vals(:,2));
        check_stuff.(sofinterest).freq.median = median(check_stuff.(sofinterest).freq.vals(:,2));
        check_stuff.(sofinterest).freq.iqr = iqr(check_stuff.(sofinterest).freq.vals(:,2));
        
        check_stuff.(sofinterest).freq.prctiles = num2str(threshold);
        check_stuff.(sofinterest).freq.threshold = prctile(check_stuff.(sofinterest).freq.vals(:,2), threshold);
        
        
        %To see if FF  is normally distributed - LT commented out, not needed
        %     try
        %         if kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 1;
        %             syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'no';
        %         elseif kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 0;
        %             syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'yes';
        %         end
        %     catch err
        %         continue
        %     end
        
        %Histogram of FF distrobution
        figure()
        [f,x] = hist(check_stuff.(sofinterest).freq.vals(:,2),30);
        bar(x,f/trapz(x,f))
        title(['FF Probability density function for ' parameters.sofinterest])
        xlabel('Frequency (Hz)')
        ylabel('Density')
        saveas(figure(gcf), [SaveDir3 '/FF_PDF_' parameters.sofinterest '_' parameters.today_date], 'fig')
        
        %Plots FF of labeled syllables
        figure(), hold on
        plot(24.*(check_stuff.(sofinterest).freq.vals(:,1)-floor(check_stuff.(sofinterest).freq.vals(:,1))),...
            check_stuff.(sofinterest).freq.vals(:,2), 'o');
        title(['FF of ' parameters.sofinterest])
        xlabel('Time (hours)')
        ylabel('Frequency (Hz)')
        hold off
        saveas(figure(gcf),[SaveDir3 '/check_pitch_shift_' parameters.sofinterest '_' parameters.today_date],'fig')
        
    end
end

%% FOURTH, Displaying information

%To know the time range that you analyzed
try
    timerange{1} = datestr(min(check_stuff.(sofinterest).freq.vals(:,1)));
    timerange{2} = datestr(max(check_stuff.(sofinterest).freq.vals(:,1)));
    
    display(timerange)
catch err
end

%Summary information about hit rate
display(' ')
display(['syllable: ' parameters.sofinterest ' hit rate information'])
display(check_stuff.(sofinterest).hit)
display(' ')

%Summary information about match rate
try % maybe did not analyze
    display(' ')
    display(['syllable: ' parameters.sofinterest ' match rate information'])
    display(check_stuff.(sofinterest).match)
    display(' ')
catch err
end

%Summary information about FF
try % i.e. maybe did not analyze FF
    display(' ')
    display(['syllable: ' parameters.sofinterest ' FF information'])
    display(check_stuff.(sofinterest).freq)
    display(' ')
catch err
end

% QUICK SUMMARY
display(' ')
display(['syllable: ' parameters.sofinterest ' QUICK SUMMARY'])
disp('# Labeled')
disp(check_stuff.(sofinterest).match.sum(2));
disp('# WN hit (on target)')
disp(check_stuff.(sofinterest).hit.sum(1));
disp('# WN hit (all)')
disp(check_stuff.(sofinterest).hit.sum(3));
disp('# Offline match (on target, ignoring pitch contingency)')
disp(check_stuff.(sofinterest).match.sum(1));
disp('# Offline match (all)')
disp(check_stuff.(sofinterest).match.sum(3));

% - TO DO: give summary of pitch of hits, escapes, and misses, and the types of things hit that are false positive (syllables.b.hit.trigger.nontrignt);

%% FIFTH, Save Data

check_stuff.parameters=parameters;
save([SaveDir3 '/check_stuff'],'check_stuff')

