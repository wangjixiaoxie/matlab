function [ songlength ] = db_song_length( batchfile )
%db_song_length Takes a batchfile and finds the song lengths (in msecs)
%   

%Opens the batchfile
fid = fopen(batchfile);

%creates relevant variables
songlength = {};
%gets first line of batchfile
tline = fgetl(fid);
i = 1;

while ischar(tline)
    %removes the .cbin or .cbin.not.mat part of the line
    rec_file = tline(1:strfind(tline,'.cbin'));
    
    %opens the relevant rec file
    fid2 = fopen([rec_file 'rec']);
    if fid2 ~= -1
        %skips the first four lines of rec file
        for j = 1:4
            fgetl(fid2);
        end
        %gets line with song length
        songlength{i} = fgetl(fid2);
        %gets rid of all characters before song length
        songlength{i} = songlength{i}(strfind(songlength{i},'=')+2:end);
        %gets rid of all characters after song length
        songlength{i} = str2num(songlength{i}(1:strfind(songlength{i},' ')-1));
    else
    end
    
    %goes to next line of batchfile
    tline = fgetl(fid);
    %goes to next cell of songlength
    i = i+1;
end

%orients songlength to be one column
songlength = songlength';
%converts songlength to a vector
songlength = cell2mat(songlength);
        
            


end

