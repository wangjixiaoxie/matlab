%% LT - given datenums that incorporate real time info, convert to within day values in units of hour or day (e.g. may 2nd 12pm --> 12.0hr or 0.5days)
% that info ignores day. to get info on day, the output day is a datenum
% that incorporates only day info.
% e.g. input_date=datenum('31May2014-1237','ddmmmyyyy-HHMM')
%   [a b]=lt_convert_datenum_to_hour(input_date);
%   a = 
%     ddmmmyyyy: '31May2014'
%       datenum: 735750
%   b = 
%     hours: 12.6167
%      days: 0.5257

function [DaysVals WithinDayVals]=lt_convert_datenum_to_hour(values_datenum);

[~,~,~,hours,minutes,seconds]=datevec(values_datenum); 
values_hours=hours+minutes.*(1/60)+seconds.*(1/3600);
values_days=values_hours/24;
WithinDayVals.hours=values_hours;
WithinDayVals.days=values_days;


days_ddmmmyyyy=datestr(values_datenum,'ddmmmyyyy');
days_datenum=datenum(days_ddmmmyyyy,'ddmmmyyyy');
DaysVals.ddmmmyyyy=days_ddmmmyyyy;
DaysVals.datenum=days_datenum;

end

