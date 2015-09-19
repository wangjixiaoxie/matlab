function [ timing ] = db_get_timing( filenames )
%db_get_timing Reads a list of filenames (in a cell) of songs to get timing info


%% Making a matrix of dates

date_matrix = cell(1,length(filenames));
timing = cell(1,length(filenames));
for i = 1:length(filenames)
    hyphen_locator = find(filenames{i} == '_');
    %locates the year
    date_matrix{i}{1} = str2double(['20' filenames{i}(hyphen_locator(1)+5:hyphen_locator(1)+6)]);
    %locates the month
    date_matrix{i}{2} = str2double(filenames{i}(hyphen_locator(1)+3:hyphen_locator(1)+4));
    %locates the day
    date_matrix{i}{3} = str2double(filenames{i}(hyphen_locator(1)+1:hyphen_locator(1)+2));
    %locates the hour
    date_matrix{i}{4} = str2double(filenames{i}(hyphen_locator(2)+1:hyphen_locator(2)+2));
    %locates the minutes
    date_matrix{i}{5} = str2double(filenames{i}(hyphen_locator(2)+3:hyphen_locator(2)+4));
    %sets seconds to 0
    date_matrix{i}{6} = 0;
    
    timing{i} = datenum(date_matrix{i}{1},date_matrix{i}{2},date_matrix{i}{3},date_matrix{i}{4},date_matrix{i}{5},date_matrix{i}{6});
end
timing = timing';

%converts to a vector
timing = cell2mat(timing);

