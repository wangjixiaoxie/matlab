

calls_folders = dir('all_calls/');

j=1;
for i = 1:length(calls_folders)
    if ~isempty(regexp(calls_folders(i).name,'[1-9]*', 'once'))
        calls{j} = calls_folders(i).name;
        j = j+1;
    end
end


folders = dir('.');
j=1;
for i = 1:length(folders)
    if ~isempty(regexp(folders(i).name,'[1-9]*', 'once'))
        folds{j} = folders(i).name;
        j = j+1;
    end
end

for i = 1:length(calls)
    for j = 1:length(folds)
        if datenum(calls{i}) == datenum(str2double(['20' folds{j}(5:6)]), str2double(folds{j}(1:2)), str2double(folds{j}(3:4)))
            try
                movefile(['all_calls/' calls{i} '/*'], [folds{j} '/'])
                display(['moving ' calls{i} ' to ' folds{j}])
            catch err
                display(['error in movefile for ' calls{i} ' to ' folds{j}])
                continue
            end
        end
    end
end