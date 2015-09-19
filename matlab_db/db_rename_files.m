function [] = db_rename_files(original_name,new_name)
%db_rename_files To rename files in a folder
%   Small program to rename files, please enter your 'original_name',
%   'new_name' 



files = dir(['*' original_name '*']);
for i = 1:length(files)
    nameOrig = files(i).name;
    nameNew = strrep(nameOrig,original_name,new_name);
    movefile(nameOrig,nameNew);
end



end

