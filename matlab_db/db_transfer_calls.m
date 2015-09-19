function [] = db_transfer_calls( copy_or_move )
%db_transfer_calls Copies (0) or moves (1) all song files to an all_calls folder
    %After you are done with labeling for the day, this will create a
%*.cbin.not.mat file and copy/move all songs to a folder within all_calls that
%is from this day. It will overwrite any file called 'batch.catch.keep'

%% Gets info on bird and date

fprintf('\nGetting directory info...')

%gets current working directory
folder = pwd;
%finds dashes to pull out date and bird name later
dashes = strfind(folder,'/');

%bird name
birdname = folder(dashes(end-1)+1:dashes(end)-1);

%Calculating the date
month = folder(dashes(end)+1:dashes(end)+2);
day = folder(dashes(end)+3:dashes(end)+4);
year = ['20' folder(dashes(end)+5:dashes(end)+6)];

date_folder = datenum(str2double(year), str2double(month), str2double(day));
date_folder = datestr(date_folder,'ddmmmyyyy');

fprintf('done!\n')

%% Makes a folder with the date in all calls
fprintf('\nMaking directory in all_calls...')

if copy_or_move == 0 || copy_or_move == 1
    if ~exist([folder(1:dashes(4)) 'all_calls/'],'dir')
        mkdir([folder(1:dashes(4)) 'all_calls/'])
        mkdir([folder(1:dashes(4)) 'all_calls/' date_folder])
    elseif ~exist([folder(1:dashes(4)) 'all_calls/' date_folder],'dir')
        mkdir([folder(1:dashes(4)) 'all_calls/' date_folder])
    else
    end
end

fprintf('done!\n')

%% Writes a *cbin.not.mat file named birdname_date
fprintf('\nWriting all_cbin_not_mat...')

fid_0 = fopen('all_cbin_not_mat', 'w');
dir_contents = dir('*cbin.not.mat');
for i = 1:length(dir_contents);
    fprintf(fid_0, '%s\n', dir_contents(i).name);
end
fclose(fid_0);

fprintf('done!\n')

%% Goes through all_cbin_not_mat and saves only the songs to two files: 
% one for the current file called 'batch.catch.keep', and one for
% the all calls folder called bird_name_date
fprintf('\nFiguring out which files have labeled songs...')

fid = fopen('all_cbin_not_mat', 'r');

%The song list for the all calls folder
fid_2 = fopen([birdname '_' date_folder], 'w');

%The song list for the current folder (cbin)
fid_3 = fopen('batch.catch.keep', 'w');



next_line = fgetl(fid);
% cbin_line = next_line(1:end-8);

while ischar(next_line)
    load(next_line)
    if sum(labels ~= '-') >= 1
        fprintf(fid_2, '%s\n', next_line);
        fprintf(fid_3, '%s\n', next_line(1:end-8));
        disp(next_line)
        if copy_or_move == 0
            copyfile([next_line(1:end-12) '*'], [folder(1:dashes(4)) 'all_calls/' date_folder]);
        elseif copy_or_move == 1
            movefile([next_line(1:end-12) '*'], [folder(1:dashes(4)) 'all_calls/' date_folder]);
        else
        end
    else
    end
    next_line = fgetl(fid);
end
fclose(fid);
fclose(fid_2);
fclose(fid_3);

%The song list for the current folder (cbin.not.mat)
if copy_or_move == 0 || copy_or_move == 1
    copyfile([birdname '_' date_folder], 'all_song_cbin_not_mat')
    
    movefile([birdname '_' date_folder], [folder(1:dashes(4)) 'all_calls/' date_folder])
end

if copy_or_move == 0
    fid_4 = fopen('1_COPIED_SONGS','w');
    fclose(fid_4);
elseif copy_or_move == 1
    fid_4 = fopen('1_TRANSFERRED_SONGS','w');
    fclose(fid_4);
else
    fid_4 = fopen('1_RAN_ALL_CALLS','w');
    fclose(fid_4);
end

fprintf('done\n')
    
%% Creates another file converting all_cbin_not_mat to just a list of cbin files

db_batch_convert('all_cbin_not_mat','all_cbin');

%% Clear variables

clear all

end

