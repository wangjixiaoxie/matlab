function [] = db_change_filename(old_string, new_string )
%db_change_filename Replaces the substring old_string with the substring
%new_string for all matching files within a folder.

%Finds all the files with old_string
files = dir(['*' old_string '*']);

%If substring does not match any files, will say so
if isempty(files) == 1
   display('No files with that substring')
   return
end

%loops through each filename, replaces it with new_string, then 'moves' the
%file in the directory (and therefore replaces it)
for i = 1:length(files)
    newname = files(i).name;
    newname = strrep(newname, old_string, new_string);
    
    movefile(files(i).name, newname);
    display(['Changing ' files(i).name ' to ' newname])
end



end

