function [ time, year, month, days, hours, minutes, seconds ] = db_timing4_file( filename )
%db_timing4 Gets the timing for a specific cbin file
%   Reads a filename from evtaf and gives you time (absolute matlab time)


underscore_locator = find(filename=='_');
year = str2double(filename(underscore_locator(1)+5:underscore_locator(1)+6))+2000;

month = str2double(filename(underscore_locator(1)+3:underscore_locator(1)+4));

days = str2double(filename(underscore_locator(1)+1:underscore_locator(1)+2));

hours = str2double(filename(underscore_locator(2)+1:underscore_locator(2)+2));

minutes = str2double(filename(underscore_locator(2)+3:underscore_locator(2)+4));

if strcmpi(filename(underscore_locator(2)+5:underscore_locator(2)+6),'.-')
    seconds = str2double([filename(underscore_locator(2)+7) filename(underscore_locator(2)+8)]);
else
    seconds = str2double(filename(underscore_locator(2)+6:underscore_locator(2)+7));
end

time = datenum(year, month, days, hours, minutes, seconds);


