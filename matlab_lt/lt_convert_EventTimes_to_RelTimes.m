%% LT - given cell array of event times {'20Mar2014-1624',...}, and a date of first day in analysis (e.g. '15Mar2014'), convert event times to format of days (e.g. output: 1.5 is day1, 12th hr)

% POSSIBLE WAYS TO CALL:
% 1) event times include hour info
% Example:
% event times: {'24Dec2014-1240','27Dec2014-1714'}
% first day: '22Dec2014';
% Output:
%         EventTimes.Input % the original input event times
%         EventTimes.JustDays_rel % index - only looks at the day info (i.e. if ddmmmyyyy-HHMM, then will only use days (so 2400 will not roll over to next day)
%         EventTimes.JustDays_actual % same as above, but string form
%         EventTimes.FinalValue=day_plus_hour; % adds up day and time info (i.e. if 2430, then rolls half hour into next day)


% 2) event times only have day info.
% Example:
% event times: {'24Dec2014','27Dec2014'}
% first day: '22Dec2014';

% 3) event times are matrix of datenums
% Example:
% event times: (735925.639456019, 735925.639456019)
% first day: '22Dec2014';
% output: as above, but justdaysrel will be the same day as that of
% finalvalue.




function [EventTimes]= lt_convert_EventTimes_to_RelTimes(first_day,events_Actual);

% Global things
first_day_datenum=datenum(first_day,'ddmmmyyyy');


% Run
if iscell(events_Actual) % DO this if events are in cell array
    
    if length(events_Actual{1})>9; % i.e. this has time and date info
        
        % get day + time info (i.e. 2400 rolls into next day)
        events_datenum=datenum(events_Actual,'ddmmmyyyy-HHMM'); % takes hours into acct (i.e. if 2400, then is next day)
        
        [DayVals, WithinDayVals]=lt_convert_datenum_to_hour(events_datenum);
        
        event_RelDays=DayVals.datenum-first_day_datenum+1;
        day_plus_hour=event_RelDays+WithinDayVals.days;
        
        % just get the day (ignore time info)
        events_datenum_DayOnly=datenum(events_Actual,'ddmmmyyyy'); % only looks at the ddmmmyyyy part (so 2400 is still the same day).
        event_RelDays_DayOnly=events_datenum_DayOnly-first_day_datenum+1;
        
        %output
        EventTimes.Input=events_Actual;
        EventTimes.JustDays_rel=event_RelDays_DayOnly; % only looks at the day info
        EventTimes.JustDays_actual=datestr(events_datenum_DayOnly,'ddmmmyyyy');
        EventTimes.FinalValue=day_plus_hour; % adds up day and time info (i.e. if 2400, then rolls into next day)
        
    else % i.e. only has date ddmmmyyyy
        events_datenum=datenum(events_Actual,'ddmmmyyyy');
        
        
        event_RelDays=events_datenum-first_day_datenum+1;
        
        EventTimes.Input=events_Actual;
        EventTimes.JustDays_rel=event_RelDays;
        EventTimes.JustDays_actual=events_Actual;
        EventTimes.FinalValue=event_RelDays;
        
    end
elseif isnumeric(events_Actual); % DO THIS if events are in matrix array
    
    
    events_datenum=events_Actual;
    [DayVals WithinDayVals]=lt_convert_datenum_to_hour(events_datenum);
    
    event_RelDays=DayVals.datenum-first_day_datenum+1;
    
    if size(event_RelDays,1)~=size(WithinDayVals.days,1);
        event_RelDays=event_RelDays';
    end
    
    day_plus_hour=event_RelDays+WithinDayVals.days;
    
    %output
    EventTimes.Input=events_Actual;
    EventTimes.JustDays_rel=event_RelDays;
    EventTimes.JustDays_actual=DayVals.ddmmmyyyy;
    EventTimes.FinalValue=day_plus_hour;
    
end
end
%
% for i=1:length(events_Actual);
%     day_rel(i)=events_Actual{i}
%
%     lt_convert_datenum_to_hour
