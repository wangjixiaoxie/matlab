%% LT 11/13/15 - changes, add, or remove, suffix of bird name of all song files (i.e. must be uscore separated pu11wh87_SUFFIX_ .....cbin)
% CHANGES CBIN, .REC AND NOTMAT FILES
function lt_change_birdname_suffix(StringToRemove, StringToAdd, saveoldfiles);
% RUN IN DAY FOLDER

% e.g.
% StringToRemove='MUSC'; If StringToRemove='', then there is no suffix, and
% will add suffix.
% StringToAdd='PBS'; If '', then will remove, and not add.
% saveoldfiles = 1; then backups in subfolder


%% Script to replace PBS with MUSC (or vice versa) in file names of all songs in day

FilesInFolder=dir('*'); % get all cbins, cbinnotmat, and rec

% copy all stuff to backup folder
if saveoldfiles==1;
    mkdir OldSongFiles
!cp * OldSongFiles;
end

% continue
for i=1:length(FilesInFolder);
    fn=FilesInFolder(i).name;
    
    if any(strfind(fn,'.cbin')) || any(strfind(fn,'.rec')) || any(strfind(fn,'.not.mat'));
    
        length_remove=length(StringToRemove);
        inds=findstr(fn, '_');
        
        if ~isempty(StringToRemove) & ~isempty(StringToAdd)
            % then replace something
        fn_new=[fn(1:inds(1)) StringToAdd fn(inds(1)+1+length_remove:end)];        
        elseif isempty(StringToRemove) & ~isempty(StringToAdd)
            % then no suffix present, add new
        fn_new=[fn(1:inds(1)) StringToAdd '_' fn(inds(1)+1+length_remove:end)];     
        elseif ~isempty(StringToRemove) & isempty(StringToAdd)
            % then remove something, don't add anything
        fn_new=[fn(1:inds(1)) StringToAdd fn(inds(1)+2+length_remove:end)];     
        end
        eval(['!mv ' fn ' ' fn_new]);
    end
    disp(num2str(i));
end

disp('DONE!');

