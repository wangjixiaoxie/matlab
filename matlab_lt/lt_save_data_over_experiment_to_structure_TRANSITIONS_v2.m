function [all_days_all_analysis] = lt_save_data_over_experiment_to_structure_v2(batch, first_day,structure)
%% LT 5/6/14 - modified to go along with lt_all_days_all_analysis_v2 (see that code).
%% 3/14 - Specifically for lt_get_all_transitions...

%% LT 12/6/13 - 
% modified from db_plot_over_experiment.  Latter takes individual day data
% saved in single folder (e.g. of entropy), and plots over days.  Here I
% modify to save all days into one structure, instead of plotting over
% days.  



% structure: input structure.  this stores all data.  want to have as an
% argument because repeatedly running this function "tacks" on data, and
% that would erase the old structure if it is declared anew each iteration.
% generally call it "all_days_all_analysis" (or same thing as output
% structure)
%% to do: separate by song (have to do before making the say structures that are opened by this function

% FIRST, 

% PUT the tukey stuff in.

% %% variables
% batch
% tukey=0;
% metric='entropy'
% day (this is entered later)
% first_day=input('what is the first day? (e.g. 12Nov2013)','s');
% first_day='04Oct2013'
% syllable = 'a';

all_days_all_analysis=structure;

%%
% open and cycle thru batch of data (of each day)
fid = fopen(batch,'r');

tline = fgetl(fid);
i=1;
while ischar(tline)
    load(tline)
    
    rel_day=datenum(all_trans.parameters.date_of_exp)-datenum(first_day)+1;
    all_days_all_analysis.data.all_songs.all_trans{rel_day}=all_trans;
    
    
%                     if isfield(all_days_all_analysis.data.date{rel_day},'var')==0;
%                     all_days_all_analysis.data.date{rel_day}=all_trans.parameters.date_of_exp;
%                 else
%                     if strcmp(all_days_all_analysis.data.date{rel_day},all_trans.parameters.date_of_exp)==0;
%                         disp(['all_trans error: the (index-wise) same all_days_all_analysis structure day was ' all_days_all_analysis.data.date{rel_day}...
%                             ', but the equivalent day for all_trans structure was ' all_trans.parameters.date_of_exp]);
%                     end
%                     end
try % first, the field will not exist (will get err), then field exists (so need to not pass if - if pass "if," then make sure date is correct)
    if ~isempty(all_days_all_analysis.data.date{rel_day})==1; 
        if strcmp(all_days_all_analysis.data.date{rel_day},all_trans.parameters.date_of_exp)==0;
        disp(['all_trans error: the (index-wise) same all_days_all_analysis structure day was ' all_days_all_analysis.data.date{rel_day}...
            ', but the equivalent day for all_trans structure was ' all_trans.parameters.date_of_exp]);
        end
    else
        all_days_all_analysis.data.date{rel_day}=all_trans.parameters.date_of_exp;
    end
catch err % if that date does not exist yet in the all_days structure.
    all_days_all_analysis.data.date{rel_day}=all_trans.parameters.date_of_exp;
end


    
    %
    %     if i == 1
    %         if isempty(time) == 1
    %             find_dashes = strfind(tline,'_');
    %             start_date = tline(find_dashes(1)+1:find_dashes(2)-1);
    %             start_date = datenum(start_date)-1;
    %         end
    %     end
    
    
    
    
    
    
    
%     if strcmp(metric,'fraction')==0;
%         if isempty(data) == 1;
%         else
%             if tukey == 1;
%                 %gets rid of tukey outliers (3.5*IQR)
%                 [tukey_data, high, low] = db_tukey_outlier(data);
%                 tukey_time = removerows(time, 'ind', [high;low]);
%                 
%                 
%                 
%                 
%             elseif tukey ~= 1;
%                 
%                 day=floor(time(1));
%                 relative_day=day-datenum(first_day)+1; % converts date number into 1,2,3...
%                 
%                 all_days_all_analysis.data{relative_day}.all_songs.(syllable).(metric)(:,1)=time;% compiles data into strucutre. column 1: time.
%                 all_days_all_analysis.data{relative_day}.all_songs.(syllable).(metric)(:,2)=data;% data
%                 if isfield(all_days_all_analysis.data{relative_day},'date')==0;
%                     all_days_all_analysis.data{relative_day}.date=datestr(day,'ddmmmyyyy');
%                 else
%                     if strcmp(all_days_all_analysis.data{relative_day}.date,datestr(day,'ddmmmyyyy'))==0;
%                         disp(['error: in structure the day is previously called ' all_days_all_analysis.data{relative_day}.date...
%                             ', but now the equivalent day is ' datestr(day,'ddmmmyyyy')]);
%                     end
%                 end
%                 
%                 
%             end
%         end
%     end
%     
%     if strcmp(metric,'fraction')==1;
%         underscores=findstr(tline,'_');
%         day=datenum(tline(underscores(1)+1:underscores(2)-1));
%         relative_day=day-datenum(first_day)+1;
%         all_days_all_analysis.data{relative_day}.per_song.(syllable).(metric)=frac_per_song;
%         all_days_all_analysis.data{relative_day}.summary_of_day.(syllable).(metric){1}=percentiles;
%         all_days_all_analysis.data{relative_day}.summary_of_day.(syllable).(metric){2}='percentiles = 2.5 50 97.5';
%         
%         if isfield(all_days_all_analysis.data{relative_day},'date')==0;
%             all_days_all_analysis.data{relative_day}.date=datestr(day,'ddmmmyyyy');
%         else
%             if strcmp(all_days_all_analysis.data{relative_day}.date,datestr(day,'ddmmmyyyy'))==0;
%                 disp(['error: in structure the day is previously called ' all_days_all_analysis.data{relative_day}.date...
%                     ', but now the equivalent day is ' datestr(day,'ddmmmyyyy')]);
%             end
%         end
%         
%         
%         
%     end
    tline = fgetl(fid);
    i = i+1;
    
end