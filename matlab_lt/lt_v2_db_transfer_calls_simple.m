%% LT - use 2
%% 10/30/14 - LT modified:
% 1) saves song files to : a) catch compiled, notcatch compiled, and all compiled (what was formerly batch.catch.keep)
% 2) also have option to only look at songs within a specified batch file (and then divide those songs into the three compiled versions above)
% NOTE: the other docs that are saved (all_cbin, all_cbin_not_mat, and
% 3) Also made the total labeled batch to be called batch.labeled.All (instead of batch.catch.keep).
% 4) Still makes batch.catch.keep (= batch.labeled.all)
% [birdname]_[date]) will have all songs.

%% 1/28/14 - LT modified to work with all dir structure, not just [name]/birds/[birdname]/...
function [] = lt_db_transfer_calls
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
% birdname = folder(dashes(4)+1:dashes(5)-1);
birdname = folder(dashes(end-1)+1:dashes(end)-1);




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

% Batches:
%The song list for the current folder (cbin) (ALL SONGS)
fid_3 = fopen('batch.labeled.all', 'w');

% Catch songs
fid_5=fopen('batch.labeled.catch','w');

% Not-catch songs only
fid_6=fopen('batch.labeled.notcatch','w');

fid_7=fopen('batch.catch.keep','w');

next_line = fgetl(fid);
% cbin_line = next_line(1:end-8);

count=0;
while ischar(next_line)
    load(next_line)
    if sum(labels ~= '-') >= 1
        count=count+1;
        % This is real song - now, what batch to save it in?
        % Is it a catch song?
        
        fnr=[next_line(1:end-13) '.rec'];
        if (~exist(fnr,'file'))
            disp(['Rec file does not exist = ',fnr '; DID NOT KEEP SONG!']);
        end
        
        frec=fopen(fnr,'r');
        while (1)
            ln = fgetl(frec);
            if (~ischar(ln))
                break;
            end
            pp = findstr(ln,'Catch');
            if (length(pp)==1)
                pp = findstr(ln,'=');
                iscatch = str2num(ln((pp(end)+1):end));
                
                if (iscatch==1)
                    fprintf(fid_5,'%s\n',next_line(1:end-8));
                else
                    fprintf(fid_6,'%s\n',next_line(1:end-8));
                end
                
            end
        end
        
        fclose(frec);
        
        fprintf(fid_3, '%s\n', next_line(1:end-8));
        fprintf(fid_7, '%s\n', next_line(1:end-8));
        disp(next_line)
        %         if copy_or_move == 0
        %             copyfile([next_line(1:end-12) '*'], [folder(1:dashes(4)) 'all_calls/' date_folder]);
        %         elseif copy_or_move == 1
        %             movefile([next_line(1:end-12) '*'], [folder(1:dashes(4)) 'all_calls/' date_folder]);

    else
    end
    next_line = fgetl(fid);
end
fclose(fid);
fclose(fid_3);
fclose(fid_5);
fclose(fid_6);
fclose(fid_7);

fprintf(['Number of songs: ' num2str(count)]);

%% Creates another file converting all_cbin_not_mat to just a list of cbin files

db_batch_convert('all_cbin_not_mat','all_cbin');

%% Clear variables

clear all

end

