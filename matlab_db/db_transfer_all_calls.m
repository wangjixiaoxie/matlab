%After you are done with labeling for the day, this will create a
%*.cbin.not.mat file and move all songs to a folder within all_calls that
%is from this day.



%% Gets info on bird and date
%gets current working directory
folder = pwd;
%finds dashes to pull out date and bird name later
dashes = strfind(folder,'/');

%bird name
birdname = folder(dashes(3)+1:dashes(4)-1);

%Calculating the date
month = folder(dashes(4)+1:dashes(4)+2);
day = folder(dashes(4)+3:dashes(4)+4);
year = ['20' folder(dashes(4)+5:dashes(4)+6)];

date_folder = datenum(str2double(year), str2double(month), str2double(day));
date_folder = datestr(date_folder,'ddmmmyyyy');


%% Writes a *cbin.not.mat file named birdname_date
fid = fopen([birdname '_' date_folder], 'w');
dir_contents = dir('*cbin.not.mat');
for i = 1:length(dir_contents);
    fprintf(fid, '%s\n', dir_contents(i).name);
end
fclose(fid);

%% Makes a folder with the date in all calls
if ~exist([folder(1:dashes(4)) 'all_calls/'],'dir')
    mkdir([folder(1:dashes(4)) 'all_calls/'])
    mkdir([folder(1:dashes(4)) 'all_calls/' date_folder])
elseif ~exist([folder(1:dashes(4)) 'all_calls/' date_folder],'dir')
    mkdir([folder(1:dashes(4)) 'all_calls/' date_folder])
else
end

%% Moves file to date folder in all calls
movefile([birdname '_' day '*'], [folder(1:dashes(4)) 'all_calls/' date_folder]);

%% Clear variables

clear all

