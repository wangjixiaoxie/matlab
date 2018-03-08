function lt_batchsong_ChangeSongFname(StringToAdd, StringToRemove)

% StringToAdd='PBS'; % will add after birdname and before date string
% -- assumes that name is e.g. "wh44wh39_[datestring].[extension]

% StringToRemove = 'PBS'; assumes is structure e.g. "wh44wh39_[StringToRemove]_...

% can leave either empty.

% will change if extension is cbin, rec, notmat, or rhd

if ~exist('StringToRemove', 'var')
    StringToRemove = '';
end
%% Script to change name of all song files in a day to stick "PBS" or "MUSC" right after bird name


FilesInFolder=dir('*'); % get all cbins, cbinnotmat, and rec

% copy all stuff to backup folder
mkdir OldSongFiles
!cp * OldSongFiles;

% continue
for i=1:length(FilesInFolder);
    fn=FilesInFolder(i).name;
    
    if ~isempty(StringToAdd)
        if any(strfind(fn,'.cbin')) || any(strfind(fn,'.rec')) || ...
                any(strfind(fn,'.not.mat')) || any(strfind(fn,'.rhd'));
            
            fn_new=[fn(1:9) StringToAdd '_' fn(10:end)];
            
            eval(['!mv ' fn ' ' fn_new]);
        end
    end
    
    if ~isempty(StringToRemove)
        
        if any(strfind(fn,'.cbin')) || any(strfind(fn,'.rec')) || ...
                any(strfind(fn,'.not.mat')) || any(strfind(fn,'.rhd'));
            
            indstoremove = strfind(fn, StringToRemove);
            indstoremove = [indstoremove indstoremove+length(StringToRemove)]
            
            fn_new=[fn(1:indstoremove(1)-1) fn(indstoremove(2)+1:end)];
            
            eval(['!mv ' fn ' ' fn_new]);
        end
        
        
    end
end
