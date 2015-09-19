function [ names ] = db_list_folder_contents( directory )
%db_list_folder_contents Makes a cell with the filenames in a specified
%folder
%   Detailed explanation goes here

all_info = dir(directory);

for i = 1:max(size(all_info))
    names{i} = all_info(i).name;
end

names = names';

end

