function [ordered_filenames] = db_day_order(birdname, phrase, batch)
%db_day_order Puts the .mat file in order by date
%   Files should be labeled birdname_date_datatype. It compiles a list of
%   all the files with birdname and datatype, then organizes them by date.
%   If batch == 1, it will create a separate batchfile called
%   'Ordered_filenames_datatype'.

directory_contents = dir([birdname '*' phrase '.mat']);

filenames = cell(size(directory_contents,1),1);
dates = [];

for i = 1:size(directory_contents,1)
    filenames{i} = directory_contents(i).name;
    dates(i,:) = directory_contents(i).name(length(birdname)+2:length(birdname)+10);
end

dates = char(dates);

dates_num = datenum(dates);
dates_num = sort(dates_num);
dates_num = datestr(dates_num,'ddmmmyyyy');

ordered_filenames = {};
for i = 1:size(directory_contents,1)
    cell_location = find(~cellfun('isempty',regexp(filenames,dates_num(i,:))) == 1);
    for j = 1:length(cell_location)
        ordered_filenames{i}{j} = filenames{cell_location(j)};
    end
end

if batch == 1
    fid = fopen(phrase,'w');
    for i = 1:size(ordered_filenames,2)
        for j = 1:size(ordered_filenames{i},2)
            fprintf(fid, '%s\n', ordered_filenames{i}{j});
            display(ordered_filenames{i}{j})
        end
    end
    fclose(fid);
end

end

