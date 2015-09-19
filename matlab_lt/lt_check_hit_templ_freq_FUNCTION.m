function check_stuff=lt_check_hit_templ_freq_FUNCTION(batchf, syl, syl_pre, syl_post, syl_refract, freq_range, evtaf_ver,get_WN_hits,get_offline_match,get_FF,varargin);


%% LT 12/2/14 - NOTE ON doing repeats
% TO simulate evtaf when repeats are an important contingency, use lt_check_hit_templ_freq_v2_EvTAFv4Sim
% It is modified from this code. exactly the same except with how it
% simulates EvtafV4. see its notes for details.


%% LT 7/22 - directly modified from lt_db_check_template_timing_and_hit_rate_BIRDTAF.
% USE THIS version
% Run from day folder

% EXAMPLE SYNTAX:
% batchf='batch.labeled.all';
% syl='c';
% syl_pre='-';
% syl_post='c';
% syl_refract=0.15;
% freq_range=[2800 3800];
% evtaf_ver='v4';
% get_WN_hits=1;
% get_offline_match=1;
% get_FF=0;
% 
% 
% template_name='wh73pk61_c1_v1_2';
% cntrng_values{1}={[3 4 2.2] 'or' 'n' 'n'};
% cntrng_values{2}={[3 4 2.2] 'and' 'n' 'n'};
% col_logic='(a+b)';
% 
% check_stuff=lt_check_hit_templ_freq_FUNCTION(batchf, syl, syl_pre, syl_post, syl_refract, freq_range, ...
%     evtaf_ver,get_WN_hits,get_offline_match,get_FF,template_name,cntrng_values,col_logic);


% ARGUMENTS (some examples):
% batchf= 'batch.catch.keep';
% evtaf_ver= 'amp' or 'v4' - important for getting date info and get trigs.
% get_WN_hits=0 or 1 (get actual WN hit info?)
% get_offline_match=1; % do offline matching using template? (ADDX=1)
% get_FF=1; % Analyze FF using offline matching?
% syl = 'b'
% syl_pre = 'a'
% syl_post = ''
% syl_refract = 0.2 % seconds
% freq_range = [2800 3800];
% varargin          % contains template information
%     varargin{1}= template_name % e.g. pu11wh87ab_SyntaxDepPitchShift_v1_6; (string, without dat)
%     varargin{2}= cntrng_values (format: {[min max threshold] "and"
%     "btaf"}
%     e.g. cntrng_values{1}={[5 14 1] 'or' 'n' 'y'}, one cell entry
%     for each column
%     varargin{3} = col_logic % only necessary for evtafv4 - e.g. (a+b+c)*(d+e+f). See notes below for important thigns.

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
parameters.batchfile=batchf;
parameters.sofinterest=syl;
parameters.sofinterest_pre=syl_pre;
parameters.sofinterest_post=syl_post;
parameters.refractory_time=syl_refract; % in seconds.
parameters.freq_range=freq_range;


[parameters.bird_name, which_computer, data_date, parameters.phrase] = lt_get_birdname_date_from_dir(1);
parameters.computer = ['/bluejay' which_computer '/lucas/birds/'];

% TEMPLATE:
% (it has to be saved in the bird folder)
if get_offline_match==1; % then need template info
    template_name=varargin{1};
    load([parameters.computer parameters.bird_name '/' template_name '.dat'])
    
    % cntrng parameters
    parameters.template=eval(template_name);
    parameters.numberoftemps = size(parameters.template,2);
    parameters.num_diff_cntrng = size(parameters.template,2);
    parameters.cntrng_values=varargin{2};
    
    if strcmp(evtaf_ver,'v4');
        parameters.col_logic=varargin{3};
    end
end


% OTHERS:
parameters.data_date=data_date{2};
parameters.today_date = datestr(now, 'ddmmmyyyy_HHMMAM'); % time-stamp for analysis
parameters.today_date = parameters.today_date(parameters.today_date ~= ' ');

% Save dir
SaveDir1=[parameters.computer parameters.bird_name '/check_hit_templ_freq/'];
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
    
    
    % SECOND(b) Using the X.tmp file above, gets offline triggers
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

    % Calculate template matching and timing information
    lt_check_hit_templ_freq_FUNCTION_TemplMatch
    

    %% THIRD, FF calculation (Only using offline, in order to get even those that failed WN contingency - ADDX = 1)
    
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