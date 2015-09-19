function timestamp=lt_get_timestamp(input)
% if input==0; % get date_time
% elseif input==1; % get date
% elseif input==2; % get time


if exist('input','var')==0;
    input=0;
end

if input==0; % get date_time
timestamp = datestr(now, 'ddmmmyyyy_HHMM');
timestamp = timestamp(timestamp ~= ' ');

elseif input==1; % get date
timestamp = datestr(now, 'ddmmmyyyy');
timestamp = timestamp(timestamp ~= ' ');

elseif input==2; % get time
timestamp = datestr(now, 'HHMM');
timestamp = timestamp(timestamp ~= ' ');
   
end
