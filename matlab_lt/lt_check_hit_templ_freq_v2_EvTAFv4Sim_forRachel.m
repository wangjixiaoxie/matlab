function check_stuff=lt_check_hit_templ_freq_v2_EvTAFv4Sim(batchf, syl, syl_pre, syl_post, get_WN_hits, get_offline_match, get_FF, config, NoteNum)

%% LT 12/3/14 - This is directly modifed from lt_check_hit_templ_freq_FUNCTION
% Only difference being that here the offline triggers (i.e. what are the
% times within songs when WN would be delivered in an actual evtaf
% experiemnts (ignoreing pitch/ampl contingencies), are detected using
% EvTAFv4Sim_LT, and not get_trigt2. The difference is that here the code
% faithfully simulates EvTAF and thus takes as input simply the .evconfig2
% config file.

% Another minor point is that need to specify NoteNum as input - i.e. 0, 1,
% 2,... This indicates which note to use as the template, in case you have
% multiple notes (only applies to Evtafv4).  Specify 0 if only one note.


% THEREFORE: this code is preferred if you are dealing with repeats and
% delays, which are not included in the simulation in get_trigt2.  A
% drawback is that this code requires an actual config file, and not
% template parameters that replicate the actual config.

% THE codes are otherwise identical.  The outputs, e.g. hit rate, match rate, etc. are identical for experiemnts
% not using repeats or delays.  Any changes can be done in the
% scripts that both of these codes call.

% ONE last minor difference is that this code outputs repeat number (in
% progress).

% RUN FROM day folder.

% TO DO:
% 1) At the moment only supports having one note (i.e. evtafv4 can have
% multiple templates). Modify to do for each note in config.  Is OK now
% since can just run multiple times.
% 2) outputs repeat number

%% HOW TO RUN
% Run from day folder

% Arguments (some examples):
% batchf= 'batch.labeled.all';
% get_WN_hits=0;
% get_offline_match=1; % do offline matching using template? (ADDX=1)
% get_FF=1; % Analyze FF using offline matching?
% syl = 'b';
% syl_pre = 'c';
% syl_post = '';
% config= '/bluejay4/lucas/birds/or73pu40/060415_SeqDepPitch_preWN/config.evconfig2';
% NoteNum = 0; 
% 
% lt_check_hit_templ_freq_v2_EvTAFv4Sim(batchf, syl, syl_pre, syl_post, get_WN_hits, get_offline_match, get_FF, config, NoteNum)


% What this does:
% - Looks at template timing and hit rate (both online and offline), and gets
% a tonne of information on that.
% - Also calculates pitch at the location of the hit (time bin directly
% rpeceding match), offline,

% Before running:
% 1) Label songs
% 2) run lt_v2_db_transfer_calls(2)


%% PARAMETERS

if ~exist('NoteNum');
    NoteNum=0;
end


% Read evtafv4 config for params
ND=ReadEvTAFv4ConfigFile(config);


% INPUTS:
parameters.batchfile=batchf;
parameters.sofinterest=syl;
parameters.sofinterest_pre=syl_pre;
parameters.sofinterest_post=syl_post;
parameters.config=config;
parameters.freq_range=[ND(NoteNum+1).FreqRng.MinFreq ND(NoteNum+1).FreqRng.MaxFreq];


% [parameters.bird_name, which_computer, data_date, parameters.phrase] = lt_get_birdname_date_from_dir(1);
% parameters.computer = ['/bluejay' which_computer '/lucas/birds/'];

currdir = pwd;

slashes = findstr(currdir, '/');
uscores = findstr(currdir, '_');

tmp = find(uscores>slashes(end));
if length(tmp)>1
parameters.phrase = currdir(slashes(end)+8:uscores(tmp(2))-1);
else
  parameters.phrase =  currdir(slashes(end)+8:end);
end

% parameters.bird_name=input('Name of bird? ','s');
parameters.bird_name=currdir(slashes(end-1)+1:slashes(end)-1);
% data_date=input('what is date? (e.g. 26Jun2015) ','s');
% parameters.phrase=input('phrase that identified experiment? (e.g. RAimplant) ','s');

% OTHERS:
% parameters.data_date=data_date{2};
% parameters.data_date=input('what is date? (e.g. 26Jun2015) ','s');

parameters.data_date=datestr(datenum(currdir(slashes(end)+1:slashes(end)+6), 'mmddyy'), 'ddmmmyyyy');

parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM'); % time-stamp for analysis
parameters.today_date = parameters.today_date(parameters.today_date ~= ' ');

% Save dir
% SaveDir1=[parameters.computer parameters.bird_name '/check_hit_templ_freq/'];

SaveDir1=[pwd '/check_hit_templ_freq/'];
SaveDir2=[SaveDir1 parameters.phrase '_' parameters.data_date];
SaveDir3=[SaveDir2 '/Analysis_' parameters.today_date '_' parameters.sofinterest_pre parameters.sofinterest parameters.sofinterest_post];
mkdir(SaveDir3);

sofinterest=parameters.sofinterest; % just renaming for ease of use

%% FIRST, Calculating hit rate (using ADDX = 0), (using actual WN trig records)
if get_WN_hits==1;
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
end

%% SECOND(a) - OFFLINE TEMPLATE MATCHING (using ADDX = 1) - Make X.tmp files and set cntrng (i.e. use template to get offline matches)

if get_offline_match==1;
    
    % Using EvTAFv4Sim_LT, calculates triggers offline.
    EvTAFv4Sim_LT(parameters.batchfile,parameters.config,'obs0');
    
    
    % SECOND(c), Finds triggers for all labeled syllables
    [check_stuff.(sofinterest).match.values, check_stuff.(sofinterest).match.trigger] = ...
        triglabel2_LT_Evtafv4_v2(parameters.batchfile,parameters.sofinterest,...
        parameters.sofinterest_pre,parameters.sofinterest_post,1,0,0,1,NoteNum);

    % Calculate match and timing statistics.
    lt_check_hit_templ_freq_FUNCTION_TemplMatch
    
    
    % THIRD, FF calculation (Only using offline, in order to get even those that failed WN contingency - ADDX = 1)
    evtaf_ver='v4';
    if get_FF==1 ;
        lt_check_hit_templ_freq_FUNCTION_FFcheck
    end
end

%% FOURTH, Displaying summary information

lt_check_hit_templ_freq_FUNCTION_summary


%% FIFTH, Save Data

check_stuff.parameters=parameters;
save([SaveDir3 '/check_stuff'],'check_stuff')

end

